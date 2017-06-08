#ifndef MYTYPE_HPP_
#define MYTYPE_HPP_

#include <stdio.h>

#ifdef SINGLE
using real_t = float;
#define FMT1 "%g"
#define EMT1 "%e"
#else
using real_t = double;
#define FMT1 "%lg"
#define EMT1 "%le"
#endif

using real3 = real_t[3];

template <typename T, unsigned int N>
void
zeroArray(T (&a)[N]) {
  constexpr const T ZERO(0);
  for (unsigned int i = 0; i < N; ++i)
    a[i]              = ZERO;
}

static FILE* const screenOut = stdout;

#endif
