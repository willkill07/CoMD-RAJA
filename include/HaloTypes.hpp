#ifndef HALO_TYPES_HPP_
#define HALO_TYPES_HPP_

#include "LinkCell.hpp"
#include "MyTypes.hpp"

struct AtomMsg {
  int gid;
  int type;
  real_t rx, ry, rz;
  real_t px, py, pz;
};

struct AtomExchangeParms {
  int nCells[6];    //!< Number of cells in cellList for each face.
  int* cellList[6]; //!< List of link cells from which to load data for each
                    //! face.
  real_t pbcFactor[6][3]; //!< Whether this face is a periodic boundary.
};

struct ForceMsg {
  real_t dfEmbed;
};

struct ForceExchangeData {
  real_t* dfEmbed; //<! derivative of embedding energy
  LinkCell* boxes;
};

struct ForceExchangeParms {
  int nCells[6];     //!< Number of cells to send/recv for each face.
  int* sendCells[6]; //!< List of link cells to send for each face.
  int* recvCells[6]; //!< List of link cells to recv for each face.
};

enum HaloFaceOrder {
  HALO_X_MINUS,
  HALO_X_PLUS,
  HALO_Y_MINUS,
  HALO_Y_PLUS,
  HALO_Z_MINUS,
  HALO_Z_PLUS
};

enum HaloAxisOrder { HALO_X_AXIS, HALO_Y_AXIS, HALO_Z_AXIS };

#endif
