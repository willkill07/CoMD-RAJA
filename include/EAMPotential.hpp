/// \file
/// Compute forces for the Embedded Atom Model (EAM).

#ifndef __EAM_H
#define __EAM_H

#include "CoMDTypes.hpp"
#include "HaloTypes.hpp"
#include "HaloExchange.hpp"

/// Handles interpolation of tabular data.
///
/// \see initInterpolationObject
/// \see interpolate
struct InterpolationObject {
    real_t *values; //!< the abscissa values
    real_t x0;      //!< the starting ordinate range
    real_t invDx;   //!< the inverse of the table spacing
    int n;          //!< the number of values in the table

    InterpolationObject() = default;
    InterpolationObject(InterpolationObject const &) = default;
    InterpolationObject(InterpolationObject&&) = default;
    InterpolationObject& operator=(InterpolationObject const &) = default;

    InterpolationObject(int n, real_t x0, real_t dx, real_t *data);

    ~InterpolationObject();

    void interpolate(real_t r, real_t *f, real_t *df);

    void broadcast();

    void print(const char *fileName);
};

template <typename T>
struct SimFlat;

/// Derived struct for an EAM potential.
/// Uses table lookups for function evaluation.
/// Polymorphic with BasePotential.
/// \see BasePotential
struct EAMPotential : public BasePotential<EAMPotential> {
    InterpolationObject phi;             //!< Pair energy
    InterpolationObject rho;             //!< Electron Density
    InterpolationObject f;               //!< Embedding Energy
    HaloForceExchange forceExchange;
    ForceExchangeData forceExchangeData;
    real_t* rhobar;                       //!< per atom storage for rhobar
    real_t* dfEmbed; //!< per atom storage for derivative of Embedding

    EAMPotential(const char *dir, const char *file, const char *type);

    void prepare(SimFlat<EAMPotential> &s);
    int force_impl(SimFlat<EAMPotential> & s);
    void print_impl(FILE* file);

    friend EAMPotential make_eam_setfl(const char*, const char*);
    friend EAMPotential make_eam_funcfl(const char*, const char*);

protected:
    EAMPotential() = default;

private:
    void broadcast();
};

template <>
EAMPotential
make_potential<EAMPotential>(Command const &cmd);

#endif
