/// \file
/// Wrappers for MPI functions.

#ifndef _PARALLEL_H_
#define _PARALLEL_H_

#include "MyTypes.hpp"

/// Structure for use with MPI_MINLOC and MPI_MAXLOC operations.
struct RankReduceData {
    double val;
    int rank;
};

#ifdef DO_MPI
template <typename T>
struct MPIType;

template <> struct MPIType<int> : std::integral_constant<int, MPI_INT> {};
template <> struct MPIType<float> : std::integral_constant<int, MPI_FLOAT> {};
template <> struct MPIType<double> : std::integral_constant<int, MPI_DOUBLE> {};

#endif

struct Parallel {

/// Return total number of processors.
    static int totalRanks();

/// Return local rank.
    static int myRank();

/// Return non-zero if printing occurs from this rank.
    static int printRank();

/// Print a timestamp and message when all tasks arrive.
    static void timestampBarrier(const char *msg);

/// Wrapper for MPI_Init.
  static void init(int *argc, char*** argv);

/// Wrapper for MPI_Finalize.
    static void destroy();

/// Wrapper for MPI_Barrier(MPI_COMM_WORLD).
    static void barrier();

/// Wrapper for MPI_Sendrecv.
    static int sendReceive(void *sendBuf, int sendLen, int dest, void *recvBuf, int recvLen, int source);

    template <typename T>
    static void add(T* sendBuf, T* recvBuf, int count) {
#ifdef DO_MPI
        MPI_Allreduce(sendBuf, recvBuf, count, MPIType<T>::value, MPI_SUM, MPI_COMM_WORLD);
#else
        for (int ii=0; ii<count; ++ii)
            recvBuf[ii] = sendBuf[ii];
#endif
    }

    template <typename T>
    static void max(T* sendBuf, T* recvBuf, int count) {
#ifdef DO_MPI
        MPI_Allreduce(sendBuf, recvBuf, count, MPIType<T>::value, MPI_MAX, MPI_COMM_WORLD);
#else
        for (int ii=0; ii<count; ++ii)
            recvBuf[ii] = sendBuf[ii];
#endif
    }

/// Wrapper for MPI_Allreduce double min with rank.
    static void minRankDouble(RankReduceData *sendBuf, RankReduceData *recvBuf, int count);

/// Wrapper for MPI_Allreduce double max with rank.
    static void maxRankDouble(RankReduceData *sendBuf, RankReduceData *recvBuf, int count);

/// Wrapper for MPI_Bcast
    static void bcast(void *buf, int len);

///  Return non-zero if code was built with MPI active.
    static bool builtWithMpi();

};

#endif
