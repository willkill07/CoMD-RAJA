#include "Domain.hpp"
#include "Command.hpp"
#include "Parallel.hpp"

#include <assert.h>

Domain::Domain(int xproc, int yproc, int zproc, real3 globalExtent) {
  assert(xproc * yproc * zproc == Parallel::totalRanks());
  procGrid[0] = xproc;
  procGrid[1] = yproc;
  procGrid[2] = zproc;

  // calculate grid coordinates i,j,k for this processor
  int myRank   = Parallel::myRank();
  procCoord[0] = myRank % procGrid[0];
  myRank /= procGrid[0];
  procCoord[1] = myRank % procGrid[1];
  procCoord[2] = myRank / procGrid[1];
  // initialialize global bounds
  for (int i = 0; i < 3; i++) {
    globalMin[i]    = 0;
    globalMax[i]    = globalExtent[i];
    globalExtent[i] = globalMax[i] - globalMin[i];
  }

  // initialize local bounds on this processor
  for (int i = 0; i < 3; i++) {
    localExtent[i] = globalExtent[i] / procGrid[i];
    localMin[i]    = globalMin[i] + procCoord[i] * localExtent[i];
    localMax[i]    = globalMin[i] + (procCoord[i] + 1) * localExtent[i];
  }
}

int
Domain::processorNum(int dix, int diy, int diz) const {
  int ix = (procCoord[0] + dix + procGrid[0]) % procGrid[0];
  int iy = (procCoord[1] + diy + procGrid[1]) % procGrid[1];
  int iz = (procCoord[2] + diz + procGrid[2]) % procGrid[2];
  return ix + procGrid[0] * (iy + procGrid[1] * iz);
}

Domain
make_domain(Command const& cmd, real_t latticeConstant) {
  real3 globalExtent = {cmd.nx * latticeConstant, cmd.ny * latticeConstant,
                        cmd.nz * latticeConstant};

  return {cmd.xproc, cmd.yproc, cmd.zproc, globalExtent};
}
