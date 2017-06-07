/// \file
/// Frequently needed typedefs.

#ifndef __MYTYPE_H_
#define __MYTYPE_H_

/// \def SINGLE determines whether single or double precision is built
#ifdef SINGLE
using real_t = float;  //!< define native type for CoMD as single precision
  #define FMT1 "%g"    //!< /def format argument for floats
  #define EMT1 "%e"    //!< /def format argument for eng floats
#else
using real_t = double; //!< define native type for CoMD as double precision
  #define FMT1 "%lg"   //!< \def format argument for doubles
  #define EMT1 "%le"   //!< \def format argument for eng doubles
#endif

using real3 = real_t[3]; //!< a convenience vector with three real_t

template <typename T, unsigned int N>
void zeroArray(T (&a)[N])
{
  constexpr const T ZERO(0);
  for (unsigned int i = 0; i < N; ++i)
    a[i] = ZERO;
}

#define screenOut stdout

#endif
