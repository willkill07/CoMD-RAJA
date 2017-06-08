#ifndef COMD_TYPES_HPP_
#define COMD_TYPES_HPP_

#include "MyTypes.hpp"
#include "Command.hpp"
#include "Domain.hpp"
#include "Atoms.hpp"
#include "HaloExchange.hpp"
#include "LinkCell.hpp"
#include "BasePotential.hpp"
#include "EAMPotential.hpp"
#include "LJPotential.hpp"
#include "Timers.hpp"
#include "YAML.hpp"
#include "Constants.hpp"
#include "Random.hpp"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

struct SpeciesData {
  char name[3]; //!< element name
  int atomicNo; //!< atomic number
  real_t mass;  //!< mass in internal units

  template <typename Potential>
  SpeciesData(Potential const &pot) : atomicNo(pot.atomicNo), mass(pot.mass) {
    strncpy(name, pot.name, 3);
  }
};

struct Validate {
  double eTot0; //<! Initial total energy
  int nAtoms0;  //<! Initial global number of atoms

  template <typename Simulation>
  Validate(Simulation &sim) {
    sim.sumAtoms();
    eTot0   = (sim.ePotential + sim.eKinetic) / sim.atoms.nGlobal;
    nAtoms0 = sim.atoms.nGlobal;

    if (Parallel::printRank()) {
      fprintf(screenOut, "\n");
      printSeparator(screenOut);
      fprintf(screenOut, "Initial energy : %14.12f, atom count : %d \n", eTot0,
              nAtoms0);
      fprintf(screenOut, "\n");
    }
  }
};

template <typename PotentialType>
struct SimFlat {
  real_t ePotential{0.0}; //!< the total potential energy of the system
  real_t eKinetic{0.0};   //!< the total kinetic energy of the system
  int nSteps;             //<! number of time steps to run
  int printRate;          //<! number of steps between output
  double dt;              //<! time step

  PotentialType pot; //!< the potential
  real_t latticeConstant;
  SpeciesData species;   //<! species data (per species, not per atom)
  Domain domain;        //<! domain decomposition data
  LinkCell boxes;       //<! link-cell data
  Atoms atoms;          //<! atom data (positions, momenta, ...)

  HaloAtomExchange atomExchange;

  SimFlat(Command const &cmd)
      : nSteps(cmd.nSteps),
        printRate(cmd.printRate),
        dt(cmd.dt),
        pot(make_potential<PotentialType>(cmd)),
        latticeConstant((cmd.lat < 0.0) ? pot.lat : cmd.lat),
        species(pot),
        domain(make_domain(cmd, latticeConstant)),
        boxes{domain, pot.cutoff},
        atoms{boxes} {

    // create lattice with desired temperature and displacement.
    createFccLattice(cmd.nx, cmd.ny, cmd.nz, latticeConstant);
    setTemperature(cmd.temperature);
    randomDisplacements(cmd.initialDelta);

    atomExchange = HaloAtomExchange(domain, boxes);
    pot.prepare(*this);

    // Forces must be computed before we call the time stepper.
    startTimer(redistributeTimer);
    redistributeAtoms();
    stopTimer(redistributeTimer);

    startTimer(computeForceTimer);
    computeForce();
    stopTimer(computeForceTimer);

    kineticEnergy();
  }

  double
  timestep(int /*n*/, real_t dt) {
    for (int ii = 0; ii < nSteps; ++ii) {
      startTimer(velocityTimer);
      advanceVelocity(boxes.nLocalBoxes, 0.5 * dt);
      stopTimer(velocityTimer);

      startTimer(positionTimer);
      advancePosition(boxes.nLocalBoxes, dt);
      stopTimer(positionTimer);

      startTimer(redistributeTimer);
      redistributeAtoms();
      stopTimer(redistributeTimer);

      startTimer(computeForceTimer);
      computeForce();
      stopTimer(computeForceTimer);

      startTimer(velocityTimer);
      advanceVelocity(boxes.nLocalBoxes, 0.5 * dt);
      stopTimer(velocityTimer);
    }
    kineticEnergy();
    return ePotential;
  }

  void
  advanceVelocity(int nBoxes, real_t dt) {
    auto nAtoms = boxes.nAtoms;
    auto p      = atoms.p;
    auto f      = atoms.f;

#pragma omp parallel for
    for (int iBox = 0; iBox < nBoxes; iBox++) {
      for (int ii = 0; ii < nAtoms[iBox]; ii++) {
        int iOff = MAXATOMS * iBox + ii;
        p[iOff][0] += dt * f[iOff][0];
        p[iOff][1] += dt * f[iOff][1];
        p[iOff][2] += dt * f[iOff][2];
      }
    }
  }

  void
  advancePosition(int nBoxes, real_t dt) {
    auto nAtoms   = boxes.nAtoms;
    auto p        = this->atoms.p;
    auto r        = this->atoms.r;
    auto species  = this->species;
    auto iSpecies = this->atoms.iSpecies;

    #pragma omp parallel for
    for (int iBox = 0; iBox < nBoxes; iBox++) {
      for (int ii = 0; ii < nAtoms[iBox]; ++ii) {
        int iOff       = MAXATOMS * iBox + ii;
        real_t invMass = 1.0 / (&species)[iSpecies[iOff]].mass;
        r[iOff][0] += dt * p[iOff][0] * invMass;
        r[iOff][1] += dt * p[iOff][1] * invMass;
        r[iOff][2] += dt * p[iOff][2] * invMass;
      }
    }
  }

  void
  computeForce() {
    pot.force(*this);
  }

  void
  kineticEnergy() {
    auto nLocalBoxes = this->boxes.nLocalBoxes;
    auto nAtoms      = this->boxes.nAtoms;
    auto iSpecies    = this->atoms.iSpecies;
    auto species     = this->species;
    auto p           = this->atoms.p;

    real_t kenergy   = 0.0;

    #pragma omp parallel for reduction(+ : kenergy)
    for (int iBox = 0; iBox < nLocalBoxes; iBox++) {
      for (int ii = 0; ii < nAtoms[iBox]; ii++) {
        int iOff       = MAXATOMS * iBox + ii;
        real_t invMass = 0.5 / (&species)[iSpecies[iOff]].mass;
        kenergy += (p[iOff][0] * p[iOff][0] + p[iOff][1] * p[iOff][1] +
                    p[iOff][2] * p[iOff][2]) *
                   invMass;
      }
    }

    real_t eLocal[2] = {ePotential, kenergy};
    real_t eSum[2];
    startTimer(commReduceTimer);
    Parallel::add(eLocal, eSum, 2);
    stopTimer(commReduceTimer);

    ePotential = eSum[0];
    eKinetic   = eSum[1];
  }

  /// Update local and remote link cells after atoms have moved.
  void
  redistributeAtoms() {
    boxes.updateLinkCells(atoms);

    startTimer(atomHaloTimer);
    atomExchange.exchange(*this);
    stopTimer(atomHaloTimer);

#pragma omp parallel for
    for (int ii = 0; ii < boxes.nTotalBoxes; ++ii)
      sortAtomsInCell(atoms, boxes, ii);
  }

  void
  sumAtoms() {
    // sum atoms across all processers
    atoms.nLocal = 0;
    for (int i = 0; i < boxes.nLocalBoxes; i++) {
      atoms.nLocal += boxes.nAtoms[i];
    }

    startTimer(commReduceTimer);
    Parallel::add(&atoms.nLocal, &atoms.nGlobal, 1);
    stopTimer(commReduceTimer);
  }

  void
  printThings(int iStep, double elapsedTime) {
    static int iStepPrev = -1;
    static int firstCall = 1;
    int nEval = iStep - iStepPrev; // gives nEval = 1 for zeroth step.
    iStepPrev = iStep;
    if (!Parallel::printRank())
      return;
    if (firstCall) {
      firstCall = 0;
      fprintf(
          screenOut,
          "#                                                                                         Performance\n"
          "#  Loop   Time(fs)       Total Energy   Potential Energy     Kinetic Energy  Temperature   (us/atom)     # Atoms\n");
      fflush(screenOut);
    }
    real_t time        = iStep * dt;
    real_t eTotal      = (ePotential + eKinetic) / atoms.nGlobal;
    real_t eK          = eKinetic / atoms.nGlobal;
    real_t eU          = ePotential / atoms.nGlobal;
    real_t Temp        = (eKinetic / atoms.nGlobal) / (kB_eV * 1.5);
    double timePerAtom = 1.0e6 * elapsedTime / (double)(nEval * atoms.nLocal);
    fprintf(screenOut,
            " %6d %10.2f %18.12f %18.12f %18.12f %12.4f %10.4f %12d\n", iStep,
            time, eTotal, eU, eK, Temp, timePerAtom, atoms.nGlobal);
  }

  void
  printYaml(FILE *file) {
    // All ranks get maxOccupancy
    int maxOcc = boxes.maxOccupancy();

    // Only rank 0 prints
    if (!Parallel::printRank())
      return;

    fprintf(file, "Simulation data: \n");
    fprintf(file, "  Total atoms        : %d\n", atoms.nGlobal);
    fprintf(file, "  Min global bounds  : [ %14.10f, %14.10f, %14.10f ]\n",
            domain.globalMin[0], domain.globalMin[1], domain.globalMin[2]);
    fprintf(file, "  Max global bounds  : [ %14.10f, %14.10f, %14.10f ]\n",
            domain.globalMax[0], domain.globalMax[1], domain.globalMax[2]);
    printSeparator(file);
    fprintf(file, "Decomposition data: \n");
    fprintf(file, "  Processors         : %6d,%6d,%6d\n", domain.procGrid[0],
            domain.procGrid[1], domain.procGrid[2]);
    fprintf(file, "  Local boxes        : %6d,%6d,%6d = %8d\n",
            boxes.gridSize[0], boxes.gridSize[1], boxes.gridSize[2],
            boxes.gridSize[0] * boxes.gridSize[1] * boxes.gridSize[2]);
    fprintf(file, "  Box ze           : [ %14.10f, %14.10f, %14.10f ]\n",
            boxes.boxSize[0], boxes.boxSize[1], boxes.boxSize[2]);
    fprintf(file, "  Box factor         : [ %14.10f, %14.10f, %14.10f ] \n",
            boxes.boxSize[0] / pot.cutoff, boxes.boxSize[1] / pot.cutoff,
            boxes.boxSize[2] / pot.cutoff);
    fprintf(file, "  Max Link Cell Occupancy: %d of %d\n", maxOcc, MAXATOMS);
    printSeparator(file);
    fprintf(file, "Potential data: \n");
    pot.print(file);
    fflush(file);
  }

  void
  createFccLattice(int nx, int ny, int nz, real_t lat) {
    const real_t *localMin = domain.localMin; // alias
    const real_t *localMax = domain.localMax; // alias

    int nb         = 4; // number of atoms in the basis
    real3 basis[4] = {{0.25, 0.25, 0.25},
                      {0.25, 0.75, 0.75},
                      {0.75, 0.25, 0.75},
                      {0.75, 0.75, 0.25}};

    // create and place atoms
    int begin[3];
    int end[3];
    for (int ii = 0; ii < 3; ++ii) {
      begin[ii] = floor(localMin[ii] / lat);
      end[ii]   = ceil(localMax[ii] / lat);
    }

    real_t px, py, pz;
    px = py = pz = 0.0;
    for (int ix = begin[0]; ix < end[0]; ++ix)
      for (int iy = begin[1]; iy < end[1]; ++iy)
        for (int iz = begin[2]; iz < end[2]; ++iz)
          for (int ib = 0; ib < nb; ++ib) {
            real_t rx = (ix + basis[ib][0]) * lat;
            real_t ry = (iy + basis[ib][1]) * lat;
            real_t rz = (iz + basis[ib][2]) * lat;
            if (rx < localMin[0] || rx >= localMax[0])
              continue;
            if (ry < localMin[1] || ry >= localMax[1])
              continue;
            if (rz < localMin[2] || rz >= localMax[2])
              continue;
            int id = ib + nb * (iz + nz * (iy + ny * (ix)));
            boxes.putAtomInBox(atoms, id, 0, rx, ry, rz, px, py, pz);
          }

    // set total atoms in simulation
    startTimer(commReduceTimer);
    Parallel::add(&atoms.nLocal, &atoms.nGlobal, 1);
    stopTimer(commReduceTimer);

    assert(atoms.nGlobal == nb * nx * ny * nz);
  }

  /// Computes the center of mass velocity of the system.
  void
  computeVcm(real_t vcm[3]) {
    real_t vcmLocal[4] = {0., 0., 0., 0.};
    real_t vcmSum[4]   = {0., 0., 0., 0.};
    real_t v0          = 0.0;
    real_t v1          = 0.0;
    real_t v2          = 0.0;
    real_t v3          = 0.0;

// sum the momenta and particle masses
#pragma omp parallel for reduction(+ : v0) reduction(+ : v1) reduction( \
    + : v2) reduction(+ : v3)
    for (int iBox = 0; iBox < boxes.nLocalBoxes; ++iBox) {
      for (int ii = 0; ii < boxes.nAtoms[iBox]; ++ii) {
        int iOff     = MAXATOMS * iBox + ii;
        int iSpecies = atoms.iSpecies[iOff];
        v0 += atoms.p[iOff][0];
        v1 += atoms.p[iOff][1];
        v2 += atoms.p[iOff][2];
        v3 += (&species)[iSpecies].mass;
      }
    }

    vcmLocal[0] = v0;
    vcmLocal[1] = v1;
    vcmLocal[2] = v2;
    vcmLocal[3] = v3;

    startTimer(commReduceTimer);
    Parallel::add(vcmLocal, vcmSum, 4);
    stopTimer(commReduceTimer);

    real_t totalMass = vcmSum[3];
    vcm[0]           = vcmSum[0] / totalMass;
    vcm[1]           = vcmSum[1] / totalMass;
    vcm[2]           = vcmSum[2] / totalMass;
  }

  void
  setVcm(real_t newVcm[3]) {
    real_t oldVcm[3];
    computeVcm(oldVcm);

    real_t vShift[3] = {(newVcm[0] - oldVcm[0]), (newVcm[1] - oldVcm[1]),
                        (newVcm[2] - oldVcm[2])};

#pragma omp parallel for
    for (int iBox = 0; iBox < boxes.nLocalBoxes; ++iBox) {
      for (int ii = 0; ii < boxes.nAtoms[iBox]; ++ii) {
        int iOff     = MAXATOMS * iBox + ii;
        int iSpecies = atoms.iSpecies[iOff];
        real_t mass  = (&species)[iSpecies].mass;
        atoms.p[iOff][0] += mass * vShift[0];
        atoms.p[iOff][1] += mass * vShift[1];
        atoms.p[iOff][2] += mass * vShift[2];
      }
    }
  }

  void
  setTemperature(real_t temperature) {
// set initial velocities for the distribution
#pragma omp parallel for
    for (int iBox = 0; iBox < boxes.nLocalBoxes; ++iBox) {
      for (int ii = 0; ii < boxes.nAtoms[iBox]; ++ii) {
        int iOff         = MAXATOMS * iBox + ii;
        int iType        = atoms.iSpecies[iOff];
        real_t mass      = (&species)[iType].mass;
        real_t sigma     = sqrt(kB_eV * temperature / mass);
        uint64_t seed    = mkSeed(atoms.gid[iOff], 123);
        atoms.p[iOff][0] = mass * sigma * gasdev(&seed);
        atoms.p[iOff][1] = mass * sigma * gasdev(&seed);
        atoms.p[iOff][2] = mass * sigma * gasdev(&seed);
      }
    }
    // compute the resulting temperature
    // kinetic energy  = 3/2 kB * Temperature
    if (temperature == 0.0)
      return;
    real_t vZero[3];
    setVcm(vZero);
    kineticEnergy();
    real_t temp = (eKinetic / atoms.nGlobal) / kB_eV / 1.5;
    // scale the velocities to achieve the target temperature
    real_t scaleFactor = sqrt(temperature / temp);
#pragma omp parallel for
    for (int iBox = 0; iBox < boxes.nLocalBoxes; ++iBox) {
      for (int ii = 0; ii < boxes.nAtoms[iBox]; ++ii) {
        int iOff = MAXATOMS * iBox + ii;
        atoms.p[iOff][0] *= scaleFactor;
        atoms.p[iOff][1] *= scaleFactor;
        atoms.p[iOff][2] *= scaleFactor;
      }
    }
    kineticEnergy();
    temp = eKinetic / atoms.nGlobal / kB_eV / 1.5;
  }

  void
  randomDisplacements(real_t delta) {
#pragma omp parallel for
    for (int iBox = 0; iBox < boxes.nLocalBoxes; ++iBox) {
      for (int ii = 0; ii < boxes.nAtoms[iBox]; ++ii) {
        int iOff      = MAXATOMS * iBox + ii;
        uint64_t seed = mkSeed(atoms.gid[iOff], 457);
        atoms.r[iOff][0] += (2.0 * lcg61(&seed) - 1.0) * delta;
        atoms.r[iOff][1] += (2.0 * lcg61(&seed) - 1.0) * delta;
        atoms.r[iOff][2] += (2.0 * lcg61(&seed) - 1.0) * delta;
      }
    }
  }
};

#endif
