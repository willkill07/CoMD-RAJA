#include "Random.hpp"

#include <stdint.h>
#include <stddef.h>
#include <math.h>

real_t gasdev(uint64_t* seed)
{
   real_t rsq,v1,v2;
   do
   {
      v1 = 2.0*lcg61(seed)-1.0;
      v2 = 2.0*lcg61(seed)-1.0;
      rsq = v1*v1+v2*v2;
   } while (rsq >= 1.0 || rsq == 0.0);

   return v2 * sqrt(-2.0*log(rsq)/rsq);
}

double lcg61(uint64_t* seed)
{
   static const double convertToDouble = 1.0/UINT64_C(2305843009213693951);

   *seed *= UINT64_C(437799614237992725);
   *seed %= UINT64_C(2305843009213693951);

   return *seed*convertToDouble;
}

uint64_t mkSeed(uint32_t id, uint32_t callSite)
{
   uint32_t s1 = id * UINT32_C(2654435761);
   uint32_t s2 = (id+callSite) * UINT32_C(2654435761);

   uint64_t iSeed = (UINT64_C(0x100000000) * s1) + s2;
   for (unsigned jj=0; jj<10; ++jj)
      lcg61(&iSeed);

   return iSeed;
}
