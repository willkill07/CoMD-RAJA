#ifndef HALO_TYPES_HPP_
#define HALO_TYPES_HPP_

#include "MyTypes.hpp"
#include "LinkCell.hpp"

/// A structure to package data for a single atom to pack into a
/// send/recv buffer.  Also used for sorting atoms within link cells.
struct AtomMsg {
  int gid;
  int type;
  real_t rx, ry, rz;
  real_t px, py, pz;
};

/// Extra data members that are needed for the exchange of atom data.
/// For an atom exchange, the HaloExchangeSt::parms will point to a
/// structure of this type.
struct AtomExchangeParms {
  int nCells[6];    //!< Number of cells in cellList for each face.
  int* cellList[6]; //!< List of link cells from which to load data for each
                    //! face.
  real_t pbcFactor[6][3]; //!< Whether this face is a periodic boundary.
};

/// Package data for the force exchange.
struct ForceMsg {
  real_t dfEmbed;
};

struct ForceExchangeData {
  real_t* dfEmbed; //<! derivative of embedding energy
  LinkCell* boxes;
};

/// Extra data members that are needed for the exchange of force data.
/// For an force exchange, the HaloExchangeSt::parms will point to a
/// structure of this type.
struct ForceExchangeParms {
  int nCells[6];     //!< Number of cells to send/recv for each face.
  int* sendCells[6]; //!< List of link cells to send for each face.
  int* recvCells[6]; //!< List of link cells to recv for each face.
};

/// Don't change the order of the faces in this enum.
enum HaloFaceOrder {
  HALO_X_MINUS,
  HALO_X_PLUS,
  HALO_Y_MINUS,
  HALO_Y_PLUS,
  HALO_Z_MINUS,
  HALO_Z_PLUS
};

/// Don't change the order of the axes in this enum.
enum HaloAxisOrder { HALO_X_AXIS, HALO_Y_AXIS, HALO_Z_AXIS };

#endif
