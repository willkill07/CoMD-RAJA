#include "HaloExchange.hpp"
#include "CoMDTypes.hpp"
#include "Domain.hpp"
#include "EAMPotential.hpp"
#include "Memory.hpp"
#include "Parallel.hpp"
#include "Timers.hpp"

#include <assert.h>

#include <algorithm>

int *
mkForceSendCellList(LinkCell const &boxes, HaloFaceOrder face, int nCells) {
  int *list = comdMalloc<int>(nCells);
  int xBegin, xEnd, yBegin, yEnd, zBegin, zEnd;
  int nx = boxes.gridSize[0];
  int ny = boxes.gridSize[1];
  int nz = boxes.gridSize[2];
  switch (face) {
  case HALO_X_MINUS:
    xBegin = 0;
    xEnd   = 1;
    yBegin = 0;
    yEnd   = ny;
    zBegin = 0;
    zEnd   = nz;
    break;
  case HALO_X_PLUS:
    xBegin = nx - 1;
    xEnd   = nx;
    yBegin = 0;
    yEnd   = ny;
    zBegin = 0;
    zEnd   = nz;
    break;
  case HALO_Y_MINUS:
    xBegin = -1;
    xEnd   = nx + 1;
    yBegin = 0;
    yEnd   = 1;
    zBegin = 0;
    zEnd   = nz;
    break;
  case HALO_Y_PLUS:
    xBegin = -1;
    xEnd   = nx + 1;
    yBegin = ny - 1;
    yEnd   = ny;
    zBegin = 0;
    zEnd   = nz;
    break;
  case HALO_Z_MINUS:
    xBegin = -1;
    xEnd   = nx + 1;
    yBegin = -1;
    yEnd   = ny + 1;
    zBegin = 0;
    zEnd   = 1;
    break;
  case HALO_Z_PLUS:
    xBegin = -1;
    xEnd   = nx + 1;
    yBegin = -1;
    yEnd   = ny + 1;
    zBegin = nz - 1;
    zEnd   = nz;
    break;
  default: exit(-1);
  }

  int count = 0;
  for (int ix = xBegin; ix < xEnd; ++ix)
    for (int iy = yBegin; iy < yEnd; ++iy)
      for (int iz     = zBegin; iz < zEnd; ++iz)
        list[count++] = boxes.getBoxFromTuple(ix, iy, iz);
  assert(count == nCells);
  return list;
}

int *
mkForceRecvCellList(LinkCell const &boxes, int face, int nCells) {
  int *list = comdMalloc<int>(nCells);
  int xBegin, xEnd, yBegin, yEnd, zBegin, zEnd;
  int nx = boxes.gridSize[0];
  int ny = boxes.gridSize[1];
  int nz = boxes.gridSize[2];
  switch (face) {
  case HALO_X_MINUS:
    xBegin = -1;
    xEnd   = 0;
    yBegin = 0;
    yEnd   = ny;
    zBegin = 0;
    zEnd   = nz;
    break;
  case HALO_X_PLUS:
    xBegin = nx;
    xEnd   = nx + 1;
    yBegin = 0;
    yEnd   = ny;
    zBegin = 0;
    zEnd   = nz;
    break;
  case HALO_Y_MINUS:
    xBegin = -1;
    xEnd   = nx + 1;
    yBegin = -1;
    yEnd   = 0;
    zBegin = 0;
    zEnd   = nz;
    break;
  case HALO_Y_PLUS:
    xBegin = -1;
    xEnd   = nx + 1;
    yBegin = ny;
    yEnd   = ny + 1;
    zBegin = 0;
    zEnd   = nz;
    break;
  case HALO_Z_MINUS:
    xBegin = -1;
    xEnd   = nx + 1;
    yBegin = -1;
    yEnd   = ny + 1;
    zBegin = -1;
    zEnd   = 0;
    break;
  case HALO_Z_PLUS:
    xBegin = -1;
    xEnd   = nx + 1;
    yBegin = -1;
    yEnd   = ny + 1;
    zBegin = nz;
    zEnd   = nz + 1;
    break;
  default: assert(1 == 0);
  }

  int count = 0;
  for (int ix = xBegin; ix < xEnd; ++ix)
    for (int iy = yBegin; iy < yEnd; ++iy)
      for (int iz     = zBegin; iz < zEnd; ++iz)
        list[count++] = boxes.getBoxFromTuple(ix, iy, iz);

  assert(count == nCells);
  return list;
}

int *
mkAtomCellList(LinkCell const &boxes, HaloFaceOrder iFace, const int nCells) {
  int *list  = comdMalloc<int>(nCells);
  int xBegin = -1;
  int xEnd   = boxes.gridSize[0] + 1;
  int yBegin = -1;
  int yEnd   = boxes.gridSize[1] + 1;
  int zBegin = -1;
  int zEnd   = boxes.gridSize[2] + 1;

  if (iFace == HALO_X_MINUS)
    xEnd = xBegin + 2;
  if (iFace == HALO_X_PLUS)
    xBegin = xEnd - 2;
  if (iFace == HALO_Y_MINUS)
    yEnd = yBegin + 2;
  if (iFace == HALO_Y_PLUS)
    yBegin = yEnd - 2;
  if (iFace == HALO_Z_MINUS)
    zEnd = zBegin + 2;
  if (iFace == HALO_Z_PLUS)
    zBegin = zEnd - 2;

  int count = 0;
  for (int ix = xBegin; ix < xEnd; ++ix)
    for (int iy = yBegin; iy < yEnd; ++iy)
      for (int iz     = zBegin; iz < zEnd; ++iz)
        list[count++] = boxes.getBoxFromTuple(ix, iy, iz);
  assert(count == nCells);
  return list;
}

HaloForceExchange::HaloForceExchange() : HaloExchange<HaloForceExchange>() {
  bufCapacity = 0;
  for (int i = 0; i < 6; ++i) {
    parms.nCells[i]    = 0;
    parms.sendCells[i] = nullptr;
    parms.recvCells[i] = nullptr;
  }
}

HaloForceExchange::HaloForceExchange(Domain const &domain,
                                     LinkCell const &boxes)
    : HaloExchange<HaloForceExchange>(domain) {
  int size0   = (boxes.gridSize[1]) * (boxes.gridSize[2]);
  int size1   = (boxes.gridSize[0] + 2) * (boxes.gridSize[2]);
  int size2   = (boxes.gridSize[0] + 2) * (boxes.gridSize[1] + 2);
  int maxSize = std::max({size0, size1, size2});
  bufCapacity = (maxSize)*MAXATOMS * sizeof(ForceMsg);
  parms.nCells[HALO_X_MINUS] = size0;
  parms.nCells[HALO_Y_MINUS] = size1;
  parms.nCells[HALO_Z_MINUS] = size2;
  parms.nCells[HALO_X_PLUS]  = parms.nCells[HALO_X_MINUS];
  parms.nCells[HALO_Y_PLUS]  = parms.nCells[HALO_Y_MINUS];
  parms.nCells[HALO_Z_PLUS]  = parms.nCells[HALO_Z_MINUS];
  for (int ii = 0; ii < 6; ++ii) {
    parms.sendCells[ii] = mkForceSendCellList(
        boxes, static_cast<HaloFaceOrder>(ii), parms.nCells[ii]);
    parms.recvCells[ii] = mkForceRecvCellList(
        boxes, static_cast<HaloFaceOrder>(ii), parms.nCells[ii]);
  }
}

HaloForceExchange::~HaloForceExchange() {
  for (int ii = 0; ii < 6; ++ii) {
    comdFree(parms.sendCells[ii]);
    comdFree(parms.recvCells[ii]);
  }
}

HaloAtomExchange::HaloAtomExchange() : HaloExchange<HaloAtomExchange>() {
  bufCapacity = 0;
  for (int i = 0; i < 6; ++i) {
    parms.nCells[i]   = 0;
    parms.cellList[i] = nullptr;
  }
}

HaloAtomExchange::HaloAtomExchange(Domain const &domain, LinkCell const &boxes)
    : HaloExchange<HaloAtomExchange>(domain) {
  int size0   = (boxes.gridSize[1] + 2) * (boxes.gridSize[2] + 2);
  int size1   = (boxes.gridSize[0] + 2) * (boxes.gridSize[2] + 2);
  int size2   = (boxes.gridSize[0] + 2) * (boxes.gridSize[1] + 2);
  int maxSize = std::max({size0, size1, size2});
  bufCapacity = maxSize * 2 * MAXATOMS * sizeof(AtomMsg);
  parms.nCells[HALO_X_MINUS] =
      2 * (boxes.gridSize[1] + 2) * (boxes.gridSize[2] + 2);
  parms.nCells[HALO_Y_MINUS] =
      2 * (boxes.gridSize[0] + 2) * (boxes.gridSize[2] + 2);
  parms.nCells[HALO_Z_MINUS] =
      2 * (boxes.gridSize[0] + 2) * (boxes.gridSize[1] + 2);
  parms.nCells[HALO_X_PLUS] = parms.nCells[HALO_X_MINUS];
  parms.nCells[HALO_Y_PLUS] = parms.nCells[HALO_Y_MINUS];
  parms.nCells[HALO_Z_PLUS] = parms.nCells[HALO_Z_MINUS];

  for (int ii = 0; ii < 6; ++ii)
    parms.cellList[ii] =
        mkAtomCellList(boxes, static_cast<HaloFaceOrder>(ii), parms.nCells[ii]);

  for (int ii = 0; ii < 6; ++ii)
    for (int jj               = 0; jj < 3; ++jj)
      parms.pbcFactor[ii][jj] = 0.0;

  if (domain.procCoord[HALO_X_AXIS] == 0)
    parms.pbcFactor[HALO_X_MINUS][HALO_X_AXIS] = +1.0;
  if (domain.procCoord[HALO_X_AXIS] == domain.procGrid[HALO_X_AXIS] - 1)
    parms.pbcFactor[HALO_X_PLUS][HALO_X_AXIS] = -1.0;
  if (domain.procCoord[HALO_Y_AXIS] == 0)
    parms.pbcFactor[HALO_Y_MINUS][HALO_Y_AXIS] = +1.0;
  if (domain.procCoord[HALO_Y_AXIS] == domain.procGrid[HALO_Y_AXIS] - 1)
    parms.pbcFactor[HALO_Y_PLUS][HALO_Y_AXIS] = -1.0;
  if (domain.procCoord[HALO_Z_AXIS] == 0)
    parms.pbcFactor[HALO_Z_MINUS][HALO_Z_AXIS] = +1.0;
  if (domain.procCoord[HALO_Z_AXIS] == domain.procGrid[HALO_Z_AXIS] - 1)
    parms.pbcFactor[HALO_Z_PLUS][HALO_Z_AXIS] = -1.0;
}

HaloAtomExchange::~HaloAtomExchange() {
  for (int ii = 0; ii < 6; ++ii) {
    comdFree(parms.cellList[ii]);
  }
}

int
HaloForceExchange::loadBuffer(ForceExchangeData &data, HaloFaceOrder face,
                              ForceMsg *buf) {
  int nCells    = parms.nCells[face];
  int *cellList = parms.recvCells[face];
  int iBuf      = 0;
  for (int iCell = 0; iCell < nCells; ++iCell) {
    int iBox = cellList[iCell];
    int iOff = iBox * MAXATOMS;
    for (int ii = iOff; ii < iOff + data.boxes->nAtoms[iBox]; ++ii) {
      data.dfEmbed[ii] = buf[iBuf].dfEmbed;
      ++iBuf;
    }
  }
  return iBuf * sizeof(ForceMsg);
}

void
HaloForceExchange::unloadBuffer(ForceExchangeData &data, HaloFaceOrder face,
                                int bufSize, ForceMsg *buf) {
  assert(bufSize % sizeof(ForceMsg) == 0);

  int nCells    = parms.nCells[face];
  int *cellList = parms.recvCells[face];
  int iBuf      = 0;
  for (int iCell = 0; iCell < nCells; ++iCell) {
    int iBox = cellList[iCell];
    int iOff = iBox * MAXATOMS;
    for (int ii = iOff; ii < iOff + data.boxes->nAtoms[iBox]; ++ii) {
      data.dfEmbed[ii] = buf[iBuf].dfEmbed;
      ++iBuf;
    }
  }
  assert(iBuf == bufSize / sizeof(ForceMsg));
}

void
sortAtomsInCell(Atoms &atoms, LinkCell const &boxes, int iBox) {
  const int nAtoms = boxes.nAtoms[iBox];

  AtomMsg *tmp = comdMalloc<AtomMsg>(nAtoms);

  int begin = iBox * MAXATOMS;
  int end   = begin + nAtoms;
  for (int ii = begin, iTmp = 0; ii < end; ++ii, ++iTmp) {
    tmp[iTmp].gid  = atoms.gid[ii];
    tmp[iTmp].type = atoms.iSpecies[ii];
    tmp[iTmp].rx   = atoms.r[ii][0];
    tmp[iTmp].ry   = atoms.r[ii][1];
    tmp[iTmp].rz   = atoms.r[ii][2];
    tmp[iTmp].px   = atoms.p[ii][0];
    tmp[iTmp].py   = atoms.p[ii][1];
    tmp[iTmp].pz   = atoms.p[ii][2];
  }

  std::sort(tmp, tmp + nAtoms,
            [](AtomMsg const &a, AtomMsg const &b) { return (a.gid < b.gid); });

  // for (int i = 1; i < nAtoms; ++i)
  //  assert(tmp[i - 1].gid != tmp[i].gid);

  for (int ii = begin, iTmp = 0; ii < end; ++ii, ++iTmp) {
    atoms.gid[ii]      = tmp[iTmp].gid;
    atoms.iSpecies[ii] = tmp[iTmp].type;
    atoms.r[ii][0]     = tmp[iTmp].rx;
    atoms.r[ii][1]     = tmp[iTmp].ry;
    atoms.r[ii][2]     = tmp[iTmp].rz;
    atoms.p[ii][0]     = tmp[iTmp].px;
    atoms.p[ii][1]     = tmp[iTmp].py;
    atoms.p[ii][2]     = tmp[iTmp].pz;
  }

  comdFree(tmp);
}
