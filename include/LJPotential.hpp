#ifndef LJ_POTENTIAL_HPP_
#define LJ_POTENTIAL_HPP_

#include "CoMDTypes.hpp"

template <typename T>
struct SimFlat;

struct LJPotential : public BasePotential<LJPotential> {
  real_t sigma;
  real_t epsilon;

  LJPotential();
  LJPotential(LJPotential const &) = default;
  LJPotential(LJPotential &&)      = default;
  LJPotential &
  operator=(LJPotential const &) = default;

  void
  prepare(SimFlat<LJPotential> &s);
  void
  force_impl(SimFlat<LJPotential> &s);
  void
  print_impl(FILE *file);
};

template <>
LJPotential
make_potential<LJPotential>(Command const &cmd);

#endif
