/// \file
/// Functions to maintain link cell structures for fast pair finding.

#ifndef __LINK_CELLS_H_
#define __LINK_CELLS_H_

#include "MyTypes.hpp"

/// The maximum number of atoms that can be stored in a link cell.
#define MAXATOMS 64

struct Domain;
struct Atoms;

/// Link cell data.  For convenience, we keep a copy of the localMin and
/// localMax coordinates that are also found in the DomainsSt.
struct LinkCell {
    int gridSize[3];     //!< number of boxes in each dimension on processor
    int nLocalBoxes;     //!< total number of local boxes on processor
    int nHaloBoxes;      //!< total number of remote halo/ghost boxes on processor
    int nTotalBoxes;     //!< total number of boxes on processor
    //!< nLocalBoxes + nHaloBoxes
    real3 localMin;      //!< minimum local bounds on processor
    real3 localMax;      //!< maximum local bounds on processor
    real3 boxSize;       //!< size of box in each dimension
    real3 invBoxSize;    //!< inverse size of box in each dimension

    int *nAtoms;         //!< total number of atoms in each box

    using Neighbor = int[27];
    Neighbor* nbrBoxes;      //!< neighbor boxes for each box

    LinkCell(Domain const &domain, real_t cutoff);

    ~LinkCell();

    void putAtomInBox(Atoms &atoms, const int gid, const int iType,
                      const real_t x, const real_t y, const real_t z,
                      const real_t px, const real_t py, const real_t pz);

    int getNeighborBoxes(int iBox, int* nbrBoxes) const;

/// Calculates the link cell index from the grid coords.  The valid
/// coordinate range in direction ii is [-1, gridSize[ii]].  Any
/// coordinate that involves a -1 or gridSize[ii] is a halo link cell.
/// Because of the order in which the local and halo link cells are
/// stored the indices of the halo cells are special cases.
/// \see initLinkCells for an explanation of storage order.
    int getBoxFromTuple(int ix, int iy, int iz) const;

/// Move an atom from one link cell to another.
/// \param iId [in]  The index with box iBox of the atom to be moved.
/// \param iBox [in] The index of the link cell the particle is moving from.
/// \param jBox [in] The index of the link cell the particle is moving to.
    void moveAtom(Atoms& atoms, int iId, int iBox, int jBox);

/// \details
/// This is the first step in returning data structures to a consistent
/// state after the atoms move each time step.  First we discard all
/// atoms in the halo link cells.  These are all atoms that are
/// currently stored on other ranks and so any information we have about
/// them is stale.  Next, we move any atoms that have crossed link cell
/// boundaries into their new link cells.  It is likely that some atoms
/// will be moved into halo link cells.  Since we have deleted halo
/// atoms from other tasks, it is clear that any atoms that are in halo
/// cells at the end of this routine have just transitioned from local
/// to halo atoms.  Such atom must be sent to other tasks by a halo
/// exchange to avoid being lost.
/// \see redistributeAtoms
    void updateLinkCells(Atoms & atoms);

/// \return The largest number of atoms in any link cell.
    int maxOccupancy() const;

/// Get the index of the link cell that contains the specified
/// coordinate.  This can be either a halo or a local link cell.
///
/// Because the rank ownership of an atom is strictly determined by the
/// atom's position, we need to take care that all ranks will agree which
/// rank owns an atom.  The conditionals at the end of this function are
/// special care to ensure that all ranks make compatible link cell
/// assignments for atoms that are near a link cell boundaries.  If no
/// ranks claim an atom in a local cell it will be lost.  If multiple
/// ranks claim an atom it will be duplicated.
    int getBoxFromCoord(real3 rr) const;

/// Set the number of atoms to zero in all halo link cells.
    void emptyHaloCells();

/// Get the grid coordinates of the link cell with index iBox.  Local
/// cells are easy as they use a standard 1D->3D mapping.  Halo cell are
/// special cases.
/// \see initLinkCells for information on link cell order.
/// \param [in]  iBox Index to link cell for which tuple is needed.
/// \param [out] ixp  x grid coord of link cell.
/// \param [out] iyp  y grid coord of link cell.
/// \param [out] izp  z grid coord of link cell.
    void getTuple(int iBox, int& ixp, int& iyp, int& izp) const;

};


#endif
