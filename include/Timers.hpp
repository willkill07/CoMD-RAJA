#ifndef PERFORMANCE_TIMERS_HPP_
#define PERFORMANCE_TIMERS_HPP_

#include <stdio.h>

enum TimerHandle {
  totalTimer,
  loopTimer,
  timestepTimer,
  positionTimer,
  velocityTimer,
  redistributeTimer,
  atomHaloTimer,
  computeForceTimer,
  eamHaloTimer,
  commHaloTimer,
  commReduceTimer,
  numberOfTimers
};

#ifndef NTIMING
#define startTimer(handle) \
  do {                     \
    profileStart(handle);  \
  } while (0)
#define stopTimer(handle) \
  do {                    \
    profileStop(handle);  \
  } while (0)
#else
#define startTimer(handle)
#define stopTimer(handle)
#endif

void
profileStart(const TimerHandle handle);

void
profileStop(const TimerHandle handle);

double
getElapsedTime(const TimerHandle handle);

void
printPerformanceResults(int nGlobalAtoms, int printRate);

void
printPerformanceResultsYaml(FILE* file);
#endif
