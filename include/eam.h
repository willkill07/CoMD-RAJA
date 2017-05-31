/// \file
/// Compute forces for the Embedded Atom Model (EAM).

#ifndef __EAM_H
#define __EAM_H

#include "mytype.h"

struct BasePotential;
struct LinkCell;

/// Pointers to the data that is needed in the load and unload functions
/// for the force halo exchange.
/// \see loadForceBuffer
/// \see unloadForceBuffer
struct ForceExchangeData
{
   real_t* dfEmbed; //<! derivative of embedding energy
   LinkCell* boxes;
};

BasePotential* initEamPot(const char* dir, const char* file, const char* type);
#endif
