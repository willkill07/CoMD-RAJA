#include "LJPotential.hpp"
#include "Constants.hpp"

#include <assert.h>

constexpr const real_t POT_SHIFT = 1.0;

LJPotential::LJPotential()
    : BasePotential<LJPotential>(2.5 * 2.315, 63.55 * amuToInternalMass, 3.615,
                                 "FCC", "Cu", 29),
      sigma(2.315),
      epsilon(0.167) {}

void
LJPotential::prepare(SimFlat<LJPotential> &) {}

void
LJPotential::force_impl(SimFlat<LJPotential> &s) {
  real_t rCut2 = cutoff * cutoff;

  // zero forces and energy
  real_t ePot  = 0.0;
  s.ePotential = 0.0;
  int fSize    = s.boxes.nTotalBoxes * MAXATOMS;
#pragma omp parallel for
  for (int ii = 0; ii < fSize; ++ii) {
    zeroArray(s.atoms.f[ii]);
    s.atoms.U[ii] = 0.;
  }

  real_t s6 = sigma * sigma * sigma * sigma * sigma * sigma;

  real_t rCut6  = s6 / (rCut2 * rCut2 * rCut2);
  real_t eShift = POT_SHIFT * rCut6 * (rCut6 - 1.0);

  int nNbrBoxes = 27;

// loop over local boxes
#pragma omp parallel for reduction(+ : ePot)
  for (int iBox = 0; iBox < s.boxes.nLocalBoxes; iBox++) {
    int nIBox = s.boxes.nAtoms[iBox];

    // loop over neighbors of iBox
    for (int jTmp = 0; jTmp < nNbrBoxes; jTmp++) {
      int jBox = s.boxes.nbrBoxes[iBox][jTmp];

      assert(jBox >= 0);

      int nJBox = s.boxes.nAtoms[jBox];

      // loop over atoms in iBox
      for (int iOff = MAXATOMS * iBox; iOff < (iBox * MAXATOMS + nIBox);
           iOff++) {
        // loop over atoms in jBox
        for (int jOff = jBox * MAXATOMS; jOff < (jBox * MAXATOMS + nJBox);
             jOff++) {
          real3 dr;
          real_t r2 = 0.0;
          for (int m = 0; m < 3; m++) {
            dr[m] = s.atoms.r[iOff][m] - s.atoms.r[jOff][m];
            r2 += dr[m] * dr[m];
          }

          if (r2 <= rCut2 && r2 > 0.0) {
            // Important note:
            // from this point on r actually refers to 1.0/r
            r2            = 1.0 / r2;
            real_t r6     = s6 * (r2 * r2 * r2);
            real_t eLocal = r6 * (r6 - 1.0) - eShift;
            s.atoms.U[iOff] += 0.5 * eLocal;
            ePot += 0.5 * eLocal;

            // different formulation to avoid sqrt computation
            real_t fr = -4.0 * epsilon * r6 * r2 * (12.0 * r6 - 6.0);
            for (int m = 0; m < 3; m++) {
              s.atoms.f[iOff][m] -= dr[m] * fr;
            }
          }
        } // loop over atoms in jBox
      }   // loop over atoms in iBox
    }     // loop over neighbor boxes
  }       // loop over local boxes in system

  ePot         = ePot * 4.0 * epsilon;
  s.ePotential = ePot;
}

void
LJPotential::print_impl(FILE *file) {
  fprintf(file, "  Potential type   : Lennard-Jones\n");
  fprintf(file, "  Species name     : %s\n", name);
  fprintf(file, "  Atomic number    : %d\n", atomicNo);
  fprintf(file, "  Mass             : " FMT1 " amu\n",
          mass / amuToInternalMass);
  fprintf(file, "  Lattice Type     : %s\n", latticeType);
  fprintf(file, "  Lattice spacing  : " FMT1 " Angstroms\n", lat);
  fprintf(file, "  Cutoff           : " FMT1 " Angstroms\n", cutoff);
  fprintf(file, "  Epsilon          : " FMT1 " eV\n", epsilon);
  fprintf(file, "  Sigma            : " FMT1 " Angstroms\n", sigma);
}

template <>
LJPotential
make_potential<LJPotential>(Command const &) {
  return {};
}
