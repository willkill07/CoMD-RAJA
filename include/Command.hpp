/// \file
/// Handle command line arguments.

#ifndef MYCOMMAND_H
#define MYCOMMAND_H

#include <cstdio>
#include <string>

/// A structure to hold the value of every run-time parameter that can
/// be read from the command line.
struct Command {
    char potDir[1024]{"pots"};  //!< the directory where EAM potentials reside
    char potName[1024]{""}; //!< the name of the potential
    char potType[1024]{"funcfl"}; //!< the type of the potential (funcfl or setfl)
    int doeam{0};          //!< a flag to determine whether we're running EAM potentials
    int nx{20};             //!< number of unit cells in x
    int ny{20};             //!< number of unit cells in y
    int nz{20};             //!< number of unit cells in z
    int xproc{1};          //!< number of processors in x direction
    int yproc{1};          //!< number of processors in y direction
    int zproc{1};          //!< number of processors in z direction
    int nSteps{100};         //!< number of time steps to run
    int printRate{10};      //!< number of steps between output
    double dt{1.0};          //!< time step (in femtoseconds)
    double lat{-1.0};         //!< lattice constant (in Angstroms)
    double temperature{600.0}; //!< simulation initial temperature (in Kelvin)
    double initialDelta{0.0}; //!< magnitude of initial displacement from lattice (in Angstroms)

    Command(int argc, char *argv[]);

    Command() = delete;

    Command(Command const &) = default;

    Command(Command &&) = default;

    Command &operator=(Command const &) = default;

    void printYAML(FILE *);

};

#endif
