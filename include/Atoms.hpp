/// \file
/// Initialize the atom configuration.

#ifndef __INIT_ATOMS_H
#define __INIT_ATOMS_H

#include "MyTypes.hpp"

struct LinkCell;

/// Atom data
struct Atoms {
    int nLocal;    //!< total number of atoms on this processor
    int nGlobal;   //!< total number of atoms in simulation

    int *gid;      //!< A globally unique id for each atom
    int *iSpecies; //!< the species index of the atom
    real3 *r;     //!< positions
    real3 *p;     //!< momenta of atoms
    real3 *f;     //!< forces
    real_t *U;     //!< potential energy per atom

    Atoms() = delete;
    Atoms(Atoms const&) = default;
    Atoms(Atoms&&) = default;
    Atoms& operator=(Atoms const&) = default;

    Atoms(LinkCell const &boxes);

    ~Atoms();

};

#endif
