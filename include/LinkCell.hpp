#ifndef LINK_CELL_HPP_
#define LINK_CELL_HPP_

#include "MyTypes.hpp"

static constexpr const int MAXATOMS = 64;

struct Domain;
struct Atoms;

struct LinkCell {
  int gridSize[3];
  int nLocalBoxes;
  int nHaloBoxes;
  int nTotalBoxes;

  real3 localMin;
  real3 localMax;
  real3 boxSize;
  real3 invBoxSize;

  int* nAtoms;

  using Neighbor = int[27];
  Neighbor* nbrBoxes;

  LinkCell(Domain const& domain, real_t cutoff);

  ~LinkCell();

  void
  putAtomInBox(Atoms& atoms, const int gid, const int iType, const real_t x,
               const real_t y, const real_t z, const real_t px, const real_t py,
               const real_t pz);

  int
  getNeighborBoxes(int iBox, int* nbrBoxes) const;

  int
  getBoxFromTuple(int ix, int iy, int iz) const;

  void
  moveAtom(Atoms& atoms, int iId, int iBox, int jBox);

  void
  updateLinkCells(Atoms& atoms);

  int
  maxOccupancy() const;

  int
  getBoxFromCoord(real3 rr) const;

  void
  emptyHaloCells();

  void
  getTuple(int iBox, int& ixp, int& iyp, int& izp) const;
};

#endif
