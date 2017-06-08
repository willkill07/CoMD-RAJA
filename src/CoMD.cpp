#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include "CoMDTypes.hpp"
#include "Command.hpp"
#include "EAMPotential.hpp"
#include "LJPotential.hpp"
#include "Parallel.hpp"
#include "YAML.hpp"

#define REDIRECT_OUTPUT 0

void
initSubsystems();

void
finalizeSubsystems();

template <typename T>
void
validateResult(const SimFlat<T> &sim, const Validate &val);

void
sanityChecks(Command const &cmd, double cutoff, double latticeConst,
             char const latticeType[8]);

template <typename PotentialType>
void
doSimulation(Command const &cmd);

int
main(int argc, char **argv) {
  // Prolog
  Parallel::init(&argc, &argv);
  profileStart(totalTimer);
  initSubsystems();
  Parallel::timestampBarrier("Starting Initialization\n");

  yamlAppInfo(yamlFile);
  yamlAppInfo(screenOut);

  Command cmd(argc, argv);
  cmd.printYAML(yamlFile);
  cmd.printYAML(screenOut);

  if (cmd.doeam) {
    doSimulation<EAMPotential>(cmd);
  } else {
    doSimulation<LJPotential>(cmd);
  }

  finalizeSubsystems();
  profileStop(totalTimer);
  Parallel::timestampBarrier("CoMD Ending\n");
  Parallel::destroy();
  return 0;
}

void
initSubsystems(void) {
#if REDIRECT_OUTPUT
  freopen("testOut.txt", "w", screenOut);
#endif

  yamlBegin();
}

void
finalizeSubsystems(void) {
#if REDIRECT_OUTPUT
  fclose(screenOut);
#endif
  yamlEnd();
}

template <typename PotentialType>
void
doSimulation(Command const &cmd) {
  SimFlat<PotentialType> sim{cmd};
  sanityChecks(cmd, sim.pot.cutoff, sim.latticeConstant, sim.pot.latticeType);
  sim.printYaml(yamlFile);
  sim.printYaml(screenOut);

  Validate validate(sim);

  Parallel::timestampBarrier("Initialization Finished\n");
  Parallel::timestampBarrier("Starting simulation\n");

  const int nSteps    = sim.nSteps;
  const int printRate = sim.printRate;
  int iStep           = 0;

  profileStart(loopTimer);
  while (iStep < nSteps) {
    startTimer(commReduceTimer);
    sim.sumAtoms();
    stopTimer(commReduceTimer);

    sim.printThings(iStep, getElapsedTime(timestepTimer));

    startTimer(timestepTimer);
    sim.timestep(printRate, sim.dt);
    stopTimer(timestepTimer);

    iStep += printRate;
  }
  profileStop(loopTimer);

  sim.sumAtoms();
  sim.printThings(iStep, getElapsedTime(timestepTimer));
  Parallel::timestampBarrier("Ending simulation\n");

  validateResult(sim, validate);
  printPerformanceResults(sim.atoms.nGlobal, sim.printRate);
  printPerformanceResultsYaml(yamlFile);
}

template <typename T>
void
validateResult(const SimFlat<T> &sim, const Validate &val) {
  if (Parallel::printRank()) {
    real_t eFinal = (sim.ePotential + sim.eKinetic) / sim.atoms.nGlobal;

    int nAtomsDelta = (sim.atoms.nGlobal - val.nAtoms0);

    fprintf(screenOut, "\n");
    fprintf(screenOut, "\n");
    fprintf(screenOut, "Simulation Validation:\n");

    fprintf(screenOut, "  Initial energy  : %14.12f\n", val.eTot0);
    fprintf(screenOut, "  Final energy    : %14.12f\n", eFinal);
    fprintf(screenOut, "  eFinal/eInitial : %f\n", eFinal / val.eTot0);
    if (nAtomsDelta == 0) {
      fprintf(screenOut, "  Final atom count : %d, no atoms lost\n",
              sim.atoms.nGlobal);
    } else {
      fprintf(screenOut, "#############################\n");
      fprintf(screenOut, "# WARNING: %6d atoms lost #\n", nAtomsDelta);
      fprintf(screenOut, "#############################\n");
    }
  }
}

void
sanityChecks(Command const &cmd, double cutoff, double latticeConst,
             char const latticeType[8]) {
  int failCode = 0;
  int nProcs = cmd.xproc * cmd.yproc * cmd.zproc;
  if (nProcs != Parallel::totalRanks()) {
    failCode |= 1;
    if (Parallel::printRank())
      fprintf(screenOut,
              "\nNumber of MPI ranks must match xproc * yproc * zproc\n");
  }

  double minx  = 2 * cutoff * cmd.xproc;
  double miny  = 2 * cutoff * cmd.yproc;
  double minz  = 2 * cutoff * cmd.zproc;
  double sizex = cmd.nx * latticeConst;
  double sizey = cmd.ny * latticeConst;
  double sizez = cmd.nz * latticeConst;

  if (sizex < minx || sizey < miny || sizez < minz) {
    failCode |= 2;
    if (Parallel::printRank())
      fprintf(screenOut,
              "\nSimulation too small.\n"
              "  Increase the number of unit cells to make the simulation\n"
              "  at least (%3.2f, %3.2f. %3.2f) Ansgstroms in size\n",
              minx, miny, minz);
  }

  if (strcasecmp(latticeType, "FCC") != 0) {
    failCode |= 4;
    if (Parallel::printRank())
      fprintf(screenOut,
              "\nOnly FCC Lattice type supported, not %s. Fatal Error.\n",
              latticeType);
  }
  int checkCode = failCode;
  Parallel::bcast(&checkCode, sizeof(int));
  assert(checkCode == failCode);
  if (failCode != 0)
    exit(failCode);
}
