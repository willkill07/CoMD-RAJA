/// \file
/// Communicate halo data such as "ghost" atoms with neighboring tasks.

#ifndef __HALO_EXCHANGE_
#define __HALO_EXCHANGE_

#include "Atoms.hpp"
#include "MyTypes.hpp"
#include "HaloTypes.hpp"
#include "Domain.hpp"
#include "LinkCell.hpp"
#include "Parallel.hpp"
#include "Timers.hpp"

#include <stdlib.h>

struct Domain;

/// A polymorphic structure to store information about a halo exchange.
/// This structure can be thought of as an abstract base class that
/// specifies the interface and implements the communication patterns of
/// a halo exchange.  Concrete sub-classes supply actual implementations
/// of the loadBuffer, unloadBuffer, and destroy functions, that are
/// specific to the actual data being exchanged.  If the subclass needs
/// additional data members, these can be stored in a structure that is
/// pointed to by parms.
///
/// Designing the structure this way allows us to re-use the
/// communication code for both atom data and partial force data.
///
/// \see eamForce
/// \see redistributeAtoms
template <typename Child>
struct HaloExchange {
  /// The MPI ranks of the six face neighbors of the local domain.
  /// Ranks are stored in the order specified in HaloFaceOrder.
  int nbrRank[6] = {0,0,0,0,0,0};
  /// The maximum send/recv buffer size (in bytes) that will be needed
  /// for this halo exchange.
  int bufCapacity = 0;

  HaloExchange() = default;

  HaloExchange(Domain const& domain) {
    nbrRank[HALO_X_MINUS] = domain.processorNum(-1, 0, 0);
    nbrRank[HALO_X_PLUS]  = domain.processorNum(+1, 0, 0);
    nbrRank[HALO_Y_MINUS] = domain.processorNum(0, -1, 0);
    nbrRank[HALO_Y_PLUS]  = domain.processorNum(0, +1, 0);
    nbrRank[HALO_Z_MINUS] = domain.processorNum(0, 0, -1);
    nbrRank[HALO_Z_PLUS]  = domain.processorNum(0, 0, +1);
  }

  template <typename T>
  void
  exchange(T& data) {
    for (int iAxis = 0; iAxis < 3; ++iAxis)
      exchangeData(data, iAxis);
  }

private:
  template <typename T, typename Buffer>
  int
  load(T& data, HaloFaceOrder face, Buffer* buffer) {
    static_cast<Child*>(this)->loadBuffer(data, face, buffer);
  }

  template <typename T, typename Buffer>
  void
  unload(T& data, HaloFaceOrder face, int bufferSize,
         Buffer* buffer) {
    static_cast<Child*>(this)->unloadBuffer(data, face, bufferSize, buffer);
  }

  template <typename T>
  void
  exchangeData(T& data, int iAxis) {
    Child* self = static_cast<Child*>(this);
    using MsgType = typename Child::MsgType;

    enum HaloFaceOrder faceM = static_cast<HaloFaceOrder>(2 * iAxis);
    enum HaloFaceOrder faceP = static_cast<HaloFaceOrder>(faceM + 1);
    int nbrRankM             = nbrRank[faceM];
    int nbrRankP             = nbrRank[faceP];
    char* block              = (char*)::malloc(4 * bufCapacity);
    MsgType* sendBufM           = reinterpret_cast<MsgType*>(block);
    MsgType* sendBufP           = reinterpret_cast<MsgType*>(block + bufCapacity);
    MsgType* recvBufM           = reinterpret_cast<MsgType*>(block + 2 * bufCapacity);
    MsgType* recvBufP           = reinterpret_cast<MsgType*>(block + 3 * bufCapacity);
    int nSendM               = self->loadBuffer(data, faceM, sendBufM);
    int nSendP               = self->loadBuffer(data, faceP, sendBufP);
    startTimer(commHaloTimer);
    int nRecvP = Parallel::sendReceive(sendBufM, nSendM, nbrRankM, recvBufP,
                                       bufCapacity, nbrRankP);
    int nRecvM = Parallel::sendReceive(sendBufP, nSendP, nbrRankP, recvBufM,
                                       bufCapacity, nbrRankM);
    stopTimer(commHaloTimer);
    self->unloadBuffer(data, faceM, nRecvM, recvBufM);
    self->unloadBuffer(data, faceP, nRecvP, recvBufP);
    ::free(block);
  }
};

struct HaloForceExchange : public HaloExchange<HaloForceExchange> {
  using MsgType       = ForceMsg;
  using ExchangeParms = ForceExchangeParms;

  ExchangeParms parms;

  HaloForceExchange();
  HaloForceExchange(HaloForceExchange const&) = default;
  HaloForceExchange(HaloForceExchange&&)      = default;
  HaloForceExchange(Domain const& domain, LinkCell const& boxes);
  ~HaloForceExchange();
  HaloForceExchange&
  operator=(HaloForceExchange const&) = default;

  int
  loadBuffer(ForceExchangeData& data, HaloFaceOrder face, ForceMsg* buf);
  void
  unloadBuffer(ForceExchangeData& data, HaloFaceOrder face, int bufSize,
               ForceMsg* buf);
};

struct HaloAtomExchange : public HaloExchange<HaloAtomExchange> {
  using MsgType       = AtomMsg;
  using ExchangeParms = AtomExchangeParms;

  ExchangeParms parms;

  HaloAtomExchange();
  HaloAtomExchange(HaloAtomExchange const&) = default;
  HaloAtomExchange(HaloAtomExchange&&)      = default;
  HaloAtomExchange(Domain const& domain, LinkCell const& boxes);
  ~HaloAtomExchange();
  HaloAtomExchange&
  operator=(HaloAtomExchange const&) = default;

  template <typename Simulation>
  int
  loadBuffer(Simulation& s, HaloFaceOrder face, AtomMsg* buf) {
    real_t* pbcFactor = parms.pbcFactor[face];
    real3 shift;
    shift[0] = pbcFactor[0] * s.domain.globalExtent[0];
    shift[1] = pbcFactor[1] * s.domain.globalExtent[1];
    shift[2] = pbcFactor[2] * s.domain.globalExtent[2];

    int nCells    = parms.nCells[face];
    int* cellList = parms.cellList[face];
    int nBuf      = 0;
    for (int iCell = 0; iCell < nCells; ++iCell) {
      int iBox = cellList[iCell];
      int iOff = iBox * MAXATOMS;
      for (int ii = iOff; ii < iOff + s.boxes.nAtoms[iBox]; ++ii) {
        buf[nBuf].gid  = s.atoms.gid[ii];
        buf[nBuf].type = s.atoms.iSpecies[ii];
        buf[nBuf].rx   = s.atoms.r[ii][0] + shift[0];
        buf[nBuf].ry   = s.atoms.r[ii][1] + shift[1];
        buf[nBuf].rz   = s.atoms.r[ii][2] + shift[2];
        buf[nBuf].px   = s.atoms.p[ii][0];
        buf[nBuf].py   = s.atoms.p[ii][1];
        buf[nBuf].pz   = s.atoms.p[ii][2];
        ++nBuf;
      }
    }
    return nBuf * sizeof(AtomMsg);
  }

  template <typename Simulation>
  void
  unloadBuffer(Simulation& s, HaloFaceOrder /*face*/, int bufSize,
               AtomMsg* buf) {
    int nBuf = bufSize / sizeof(AtomMsg);
    for (int ii = 0; ii < nBuf; ++ii) {
      int gid   = buf[ii].gid;
      int type  = buf[ii].type;
      real_t rx = buf[ii].rx;
      real_t ry = buf[ii].ry;
      real_t rz = buf[ii].rz;
      real_t px = buf[ii].px;
      real_t py = buf[ii].py;
      real_t pz = buf[ii].pz;
      s.boxes.putAtomInBox(s.atoms, gid, type, rx, ry, rz, px, py, pz);
    }
  }
};

/// Sort the atoms by gid in the specified link cell.
void
sortAtomsInCell(Atoms& atoms, LinkCell const& boxes, int iBox);

#endif
