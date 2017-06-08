#ifndef PARALLEL_HPP_
#define PARALLEL_HPP_

#include "MyTypes.hpp"

struct RankReduceData {
  double val;
  int rank;
};

#ifdef DO_MPI
template <typename T>
struct MPIType;

template <>
struct MPIType<int> : std::integral_constant<int, MPI_INT> {};
template <>
struct MPIType<float> : std::integral_constant<int, MPI_FLOAT> {};
template <>
struct MPIType<double> : std::integral_constant<int, MPI_DOUBLE> {};

#endif

struct Parallel {
  static int
  totalRanks();

  static int
  myRank();

  static int
  printRank();

  static void
  timestampBarrier(const char *msg);

  static void
  init(int *argc, char ***argv);

  static void
  destroy();

  static void
  barrier();

  static int
  sendReceive(void *sendBuf, int sendLen, int dest, void *recvBuf, int recvLen,
              int source);

  template <typename T>
  static void
  add(T *sendBuf, T *recvBuf, int count) {
#ifdef DO_MPI
    MPI_Allreduce(sendBuf, recvBuf, count, MPIType<T>::value, MPI_SUM,
                  MPI_COMM_WORLD);
#else
    for (int ii   = 0; ii < count; ++ii)
      recvBuf[ii] = sendBuf[ii];
#endif
  }

  template <typename T>
  static void
  max(T *sendBuf, T *recvBuf, int count) {
#ifdef DO_MPI
    MPI_Allreduce(sendBuf, recvBuf, count, MPIType<T>::value, MPI_MAX,
                  MPI_COMM_WORLD);
#else
    for (int ii   = 0; ii < count; ++ii)
      recvBuf[ii] = sendBuf[ii];
#endif
  }

  static void
  minRankDouble(RankReduceData *sendBuf, RankReduceData *recvBuf, int count);

  static void
  maxRankDouble(RankReduceData *sendBuf, RankReduceData *recvBuf, int count);

  static void
  bcast(void *buf, int len);

  static bool
  builtWithMpi();
};

#endif
