/// \file
/// CoMD data structures.

#ifndef __COMDTYPES_H_
#define __COMDTYPES_H_

#include <cstdio>

#include "mytype.hpp"
#include "haloExchange.hpp"
#include "linkCells.hpp"
#include "decomposition.hpp"
#include "initAtoms.hpp"

struct SimFlat;

/// The base struct from which all potentials derive.  Think of this as an
/// abstract base class.
///
/// CoMD uses the following units:
///  - distance is in Angstroms
///  - energy is in eV
///  - time is in fs
///  - force in in eV/Angstrom
///
///  The choice of distance, energy, and time units means that the unit
///  of mass is eV*fs^2/Angstrom^2.  Hence, we must convert masses that
///  are input in AMU (atomic mass units) into internal mass units.
struct BasePotential
{
   real_t cutoff;          //!< potential cutoff distance in Angstroms
   real_t mass;            //!< mass of atoms in intenal units
   real_t lat;             //!< lattice spacing (angs) of unit cell
   char latticeType[8];    //!< lattice type, e.g. FCC, BCC, etc.
   char  name[3];          //!< element name
   int   atomicNo;         //!< atomic number
   int  (*force)(SimFlat* s); //!< function pointer to force routine
   void (*print)(FILE* file, BasePotential* pot);
   void (*destroy)(BasePotential** pot); //!< destruction of the potential
};


/// species data: chosen to match the data found in the setfl/funcfl files
struct SpeciesData {
   char  name[3];   //!< element name
   int   atomicNo;  //!< atomic number
   real_t mass;     //!< mass in internal units
};

/// Simple struct to store the initial energy and number of atoms.
/// Used to check energy conservation and atom conservation.
struct Validate {
   double eTot0; //<! Initial total energy
   int nAtoms0;  //<! Initial global number of atoms
};

///
/// The fundamental simulation data structure with MAXATOMS in every
/// link cell.
///
struct SimFlat {
   int nSteps;            //<! number of time steps to run
   int printRate;         //<! number of steps between output
   double dt;             //<! time step

   Domain* domain;        //<! domain decomposition data

   LinkCell* boxes;       //<! link-cell data

   Atoms* atoms;          //<! atom data (positions, momenta, ...)

   SpeciesData* species;  //<! species data (per species, not per atom)

   real_t ePotential;     //!< the total potential energy of the system
   real_t eKinetic;       //!< the total kinetic energy of the system

   BasePotential *pot;    //!< the potential

   HaloExchange* atomExchange;

};

#endif
