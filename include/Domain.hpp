/// \file
/// Parallel domain decomposition.

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "MyTypes.hpp"

/// Domain decomposition information.
struct Domain {
    // process-layout data
    int procGrid[3];        //!< number of processors in each dimension
    int procCoord[3];       //!< i,j,k for this processor

    // global bounds data
    real3 globalMin;        //!< minimum global coordinate (angstroms)
    real3 globalMax;        //!< maximum global coordinate (angstroms)
    real3 globalExtent;     //!< global size: globalMax - globalMin

    // local bounds data
    real3 localMin;         //!< minimum coordinate on local processor
    real3 localMax;         //!< maximum coordinate on local processor
    real3 localExtent;      //!< localMax - localMin

    Domain() = delete;

    Domain(Domain const &) = default;

    Domain(Domain &&) = default;

    Domain &operator=(Domain const &) = default;

    Domain(int xproc, int yproc, int zproc, real3 globalExtent);

    int processorNum(int dix, int diy, int diz) const;
};

struct Command;

Domain make_domain(Command const &, real_t);

#endif
