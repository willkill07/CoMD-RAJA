#include "LinkCell.hpp"
#include "CoMDTypes.hpp"
#include "Domain.hpp"
#include "Memory.hpp"
#include "Parallel.hpp"
#include "Timers.hpp"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>

void
copyAtom(Atoms &atoms, int iAtom, int iBox, int jAtom, int jBox);

LinkCell::LinkCell(Domain const &domain, real_t cutoff) {
  for (int i = 0; i < 3; i++) {
    localMin[i]   = domain.localMin[i];
    localMax[i]   = domain.localMax[i];
    gridSize[i]   = domain.localExtent[i] / cutoff; // local number of boxes
    boxSize[i]    = domain.localExtent[i] / ((real_t)gridSize[i]);
    invBoxSize[i] = 1.0 / boxSize[i];
  }

  nLocalBoxes = gridSize[0] * gridSize[1] * gridSize[2];

  nHaloBoxes = 2 * ((gridSize[0] + 2) * (gridSize[1] + gridSize[2] + 2) +
                    (gridSize[1] * gridSize[2]));

  nTotalBoxes = nLocalBoxes + nHaloBoxes;

  nAtoms = comdMalloc<int>(nTotalBoxes);
  for (int iBox  = 0; iBox < nTotalBoxes; ++iBox)
    nAtoms[iBox] = 0;

  assert((gridSize[0] >= 2) && (gridSize[1] >= 2) && (gridSize[2] >= 2));

  // Added creating neighbors once
  nbrBoxes = comdMalloc<Neighbor>(nTotalBoxes);
}

LinkCell::~LinkCell() {
  comdFree(nbrBoxes);
  comdFree(nAtoms);
}

int
LinkCell::getNeighborBoxes(int iBox, int *nbrBoxes) const {
  int ix, iy, iz;
  getTuple(iBox, ix, iy, iz);
  int count = 0;
  for (int i = ix - 1; i <= ix + 1; i++)
    for (int j = iy - 1; j <= iy + 1; j++)
      for (int k          = iz - 1; k <= iz + 1; k++)
        nbrBoxes[count++] = getBoxFromTuple(i, j, k);
  return count;
}

void
LinkCell::putAtomInBox(Atoms &atoms, const int gid, const int iType,
                       const real_t x, const real_t y, const real_t z,
                       const real_t px, const real_t py, const real_t pz) {
  real_t xyz[3] = {x, y, z};

  // Find correct box.
  int iBox = getBoxFromCoord(xyz);
  int iOff = iBox * MAXATOMS;
  iOff += nAtoms[iBox];

  // assign values to array elements
  if (iBox < nLocalBoxes)
    atoms.nLocal++;
  nAtoms[iBox]++;
  atoms.gid[iOff]      = gid;
  atoms.iSpecies[iOff] = iType;

  atoms.r[iOff][0] = x;
  atoms.r[iOff][1] = y;
  atoms.r[iOff][2] = z;

  atoms.p[iOff][0] = px;
  atoms.p[iOff][1] = py;
  atoms.p[iOff][2] = pz;
}

int
LinkCell::getBoxFromTuple(int ix, int iy, int iz) const {
  int iBox = 0;

  // Halo in Z+
  if (iz == gridSize[2]) {
    iBox = nLocalBoxes + 2 * gridSize[2] * gridSize[1] +
           2 * gridSize[2] * (gridSize[0] + 2) +
           (gridSize[0] + 2) * (gridSize[1] + 2) +
           (gridSize[0] + 2) * (iy + 1) + (ix + 1);
  }
  // Halo in Z-
  else if (iz == -1) {
    iBox = nLocalBoxes + 2 * gridSize[2] * gridSize[1] +
           2 * gridSize[2] * (gridSize[0] + 2) + (gridSize[0] + 2) * (iy + 1) +
           (ix + 1);
  }
  // Halo in Y+
  else if (iy == gridSize[1]) {
    iBox = nLocalBoxes + 2 * gridSize[2] * gridSize[1] +
           gridSize[2] * (gridSize[0] + 2) + (gridSize[0] + 2) * iz + (ix + 1);
  }
  // Halo in Y-
  else if (iy == -1) {
    iBox = nLocalBoxes + 2 * gridSize[2] * gridSize[1] +
           iz * (gridSize[0] + 2) + (ix + 1);
  }
  // Halo in X+
  else if (ix == gridSize[0]) {
    iBox = nLocalBoxes + gridSize[1] * gridSize[2] + iz * gridSize[1] + iy;
  }
  // Halo in X-
  else if (ix == -1) {
    iBox = nLocalBoxes + iz * gridSize[1] + iy;
  }
  // local link celll.
  else {
    iBox = ix + gridSize[0] * iy + gridSize[0] * gridSize[1] * iz;
  }
  assert(iBox >= 0);
  assert(iBox < nTotalBoxes);

  return iBox;
}

void
LinkCell::moveAtom(Atoms &atoms, int iId, int iBox, int jBox) {
  int nj = nAtoms[jBox];
  copyAtom(atoms, iId, iBox, nj, jBox);
  nAtoms[jBox]++;

  assert(nAtoms[jBox] < MAXATOMS);

  nAtoms[iBox]--;
  int ni = nAtoms[iBox];
  if (ni)
    copyAtom(atoms, ni, iBox, iId, iBox);

  if (jBox > nLocalBoxes)
    --atoms.nLocal;

  return;
}

void
LinkCell::updateLinkCells(Atoms &atoms) {
  emptyHaloCells();
  for (int iBox = 0; iBox < nLocalBoxes; ++iBox) {
    int iOff = iBox * MAXATOMS;
    int ii   = 0;
    while (ii < nAtoms[iBox]) {
      int jBox = getBoxFromCoord(atoms.r[iOff + ii]);
      if (jBox != iBox)
        moveAtom(atoms, ii, iBox, jBox);
      else
        ++ii;
    }
  }
}

int
LinkCell::maxOccupancy() const {
  int localMax = 0;
  for (int ii = 0; ii < nLocalBoxes; ++ii)
    localMax = std::max(localMax, nAtoms[ii]);

  int globalMax;

  startTimer(commReduceTimer);
  Parallel::max(&localMax, &globalMax, 1);
  stopTimer(commReduceTimer);

  return globalMax;
}

int
LinkCell::getBoxFromCoord(real3 rr) const {
  int ix = (int)(floor((rr[0] - localMin[0]) * invBoxSize[0]));
  int iy = (int)(floor((rr[1] - localMin[1]) * invBoxSize[1]));
  int iz = (int)(floor((rr[2] - localMin[2]) * invBoxSize[2]));
  // For each axis, if we are inside the local domain, make sure we get
  // a local link cell.  Otherwise, make sure we get a halo link cell.
  if (rr[0] < localMax[0]) {
    if (ix == gridSize[0])
      ix = gridSize[0] - 1;
  } else
    ix = gridSize[0]; // assign to halo cell
  if (rr[1] < localMax[1]) {
    if (iy == gridSize[1])
      iy = gridSize[1] - 1;
  } else
    iy = gridSize[1];
  if (rr[2] < localMax[2]) {
    if (iz == gridSize[2])
      iz = gridSize[2] - 1;
  } else
    iz = gridSize[2];

  return getBoxFromTuple(ix, iy, iz);
}

/// Set the number of atoms to zero in all halo link cells.
void
LinkCell::emptyHaloCells() {
  for (int ii  = nLocalBoxes; ii < nTotalBoxes; ++ii)
    nAtoms[ii] = 0;
}

void
LinkCell::getTuple(int iBox, int &ixp, int &iyp, int &izp) const {
  int ix, iy, iz;
  // If a local box
  if (iBox < nLocalBoxes) {
    ix = iBox % gridSize[0];
    iBox /= gridSize[0];
    iy = iBox % gridSize[1];
    iz = iBox / gridSize[1];
  }
  // It's a halo box
  else {
    int ink;
    ink = iBox - nLocalBoxes;
    if (ink < 2 * gridSize[1] * gridSize[2]) {
      if (ink < gridSize[1] * gridSize[2]) {
        ix = 0;
      } else {
        ink -= gridSize[1] * gridSize[2];
        ix = gridSize[0] + 1;
      }
      iy = 1 + ink % gridSize[1];
      iz = 1 + ink / gridSize[1];
    } else if (ink < (2 * gridSize[2] * (gridSize[1] + gridSize[0] + 2))) {
      ink -= 2 * gridSize[2] * gridSize[1];
      if (ink < ((gridSize[0] + 2) * gridSize[2])) {
        iy = 0;
      } else {
        ink -= (gridSize[0] + 2) * gridSize[2];
        iy = gridSize[1] + 1;
      }
      ix = ink % (gridSize[0] + 2);
      iz = 1 + ink / (gridSize[0] + 2);
    } else {
      ink -= 2 * gridSize[2] * (gridSize[1] + gridSize[0] + 2);
      if (ink < ((gridSize[0] + 2) * (gridSize[1] + 2))) {
        iz = 0;
      } else {
        ink -= (gridSize[0] + 2) * (gridSize[1] + 2);
        iz = gridSize[2] + 1;
      }
      ix = ink % (gridSize[0] + 2);
      iy = ink / (gridSize[0] + 2);
    }

    // Calculated as off by 1
    ix--;
    iy--;
    iz--;
  }
  ixp = ix;
  iyp = iy;
  izp = iz;
}

void
copyAtom(Atoms &atoms, int iAtom, int iBox, int jAtom, int jBox) {
  const int iOff       = MAXATOMS * iBox + iAtom;
  const int jOff       = MAXATOMS * jBox + jAtom;
  atoms.gid[jOff]      = atoms.gid[iOff];
  atoms.iSpecies[jOff] = atoms.iSpecies[iOff];
  memcpy(atoms.r[jOff], atoms.r[iOff], sizeof(real3));
  memcpy(atoms.p[jOff], atoms.p[iOff], sizeof(real3));
  memcpy(atoms.f[jOff], atoms.f[iOff], sizeof(real3));
  memcpy(atoms.U + jOff, atoms.U + iOff, sizeof(real_t));
}
