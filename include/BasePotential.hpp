#ifndef BASE_POTENTIAL_HPP
#define BASE_POTENTIAL_HPP

#include "HaloExchange.hpp"
#include "MyTypes.hpp"
#include "Parallel.hpp"

#include <stdio.h>
#include <string.h>

template <typename PotentialType>
struct BasePotential {
  real_t cutoff{0.0}; //!< potential cutoff distance in Angstroms
  real_t mass{0.0};   //!< mass of atoms in intenal units
  real_t lat{0.0};    //!< lattice spacing (angs) of unit cell
  char latticeType[8];
  char name[3];
  int atomicNo{-1}; //!< atomic number

  BasePotential() = default;

  BasePotential(real_t cutoff_, real_t mass_, real_t lat_,
                char const latticeType_[8], char const name_[3], int atomicNo_)
      : cutoff(cutoff_), mass(mass_), lat(lat_), atomicNo(atomicNo_) {
    strncpy(latticeType, latticeType_, 8);
    strncpy(name, name_, 3);
  }

  template <typename Simulation>
  int
  force(Simulation &s) {
    static_cast<PotentialType &>(*this).force_impl(s);
  }

  void
  print(FILE *file) const {
    auto asConstChild = static_cast<PotentialType const *>(this);
    auto asChild      = const_cast<PotentialType *>(asConstChild);
    asChild->print_impl(file);
  }
};

struct Command;

template <typename T>
T
make_potential(Command const &cmd);

#endif
