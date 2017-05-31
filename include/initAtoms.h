/// \file
/// Initialize the atom configuration.

#ifndef __INIT_ATOMS_H
#define __INIT_ATOMS_H

#include "mytype.h"

struct SimFlat;
struct LinkCell;

/// Atom data
struct Atoms {
   // atom-specific data
   int nLocal;    //!< total number of atoms on this processor
   int nGlobal;   //!< total number of atoms in simulation

   int* gid;      //!< A globally unique id for each atom
   int* iSpecies; //!< the species index of the atom

   real3*  r;     //!< positions
   real3*  p;     //!< momenta of atoms
   real3*  f;     //!< forces
   real_t* U;     //!< potential energy per atom
};


/// Allocates memory to store atom data.
Atoms* initAtoms(LinkCell* boxes);
void destroyAtoms(Atoms* atoms);

void createFccLattice(int nx, int ny, int nz, real_t lat, SimFlat* s);

void setVcm(SimFlat* s, real_t vcm[3]);
void setTemperature(SimFlat* s, real_t temperature);
void randomDisplacements(SimFlat* s, real_t delta);
#endif
