/// \file
/// Wrappers for memory allocation.

#ifndef _MEMUTILS_H_
#define _MEMUTILS_H_

#include <stdlib.h>

#define freeMe(s, element)                                                     \
  {                                                                            \
    if (s->element)                                                            \
      comdFree(s->element);                                                    \
    s->element = NULL;                                                         \
  }

template <typename T> inline T *comdMalloc(size_t iSize) {
  return (T *)malloc(iSize * sizeof(T));
}

template <typename T> inline T *comdCalloc(size_t num, size_t iSize) {
  return (T *)calloc(num, sizeof(T) * iSize);
}

template <typename T> inline T *comdRealloc(void *ptr, size_t iSize) {
  return (T *)realloc(ptr, sizeof(T) * iSize);
}

template <typename T> inline void comdFree(T *ptr) { if (ptr) free((void *)ptr); }

#endif
