/// \file
/// Wrappers for MPI functions.  This should be the only compilation
/// unit in the code that directly calls MPI functions.  To build a pure
/// serial version of the code with no MPI, do not define DO_MPI.  If
/// DO_MPI is not defined then all MPI functionality is replaced with
/// equivalent single task behavior.

#include "Parallel.hpp"

#ifdef DO_MPI
#include <mpi.h>
#endif

#include <cstdio>
#include <ctime>
#include <cstring>
#include <cassert>

static int rank = 0;
static int nRanks = 1;

#ifdef DO_MPI
#ifdef SINGLE
#define REAL_MPI_TYPE MPI_FLOAT
#else
#define REAL_MPI_TYPE MPI_DOUBLE
#endif

#endif

int Parallel::totalRanks()
{
   return nRanks;
}

int Parallel::myRank()
{
   return rank;
}

/// \details
/// For now this is just a check for rank 0 but in principle it could be
/// more complex.  It is also possible to suppress practically all
/// output by causing this function to return 0 for all ranks.
int Parallel::printRank()
{
   if (rank == 0) return 1;
   return 0;
}

void Parallel::timestampBarrier(const char* msg)
{
   barrier();
   if (! printRank())
      return;
   time_t t= time(NULL);
   char* timeString = ctime(&t);
   timeString[24] = '\0'; // clobber newline
   fprintf(screenOut, "%s: %s\n", timeString, msg);
   fflush(screenOut);
}

void Parallel::init(int *argc, char*** argv)
{
#ifdef DO_MPI
   MPI_Init(argc, argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
#endif
}

void Parallel::destroy()
{
#ifdef DO_MPI
   MPI_Finalize();
#endif
}

void Parallel::barrier()
{
#ifdef DO_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/// \param [in]  sendBuf Data to send.
/// \param [in]  sendLen Number of bytes to send.
/// \param [in]  dest    Rank in MPI_COMM_WORLD where data will be sent.
/// \param [out] recvBuf Received data.
/// \param [in]  recvLen Maximum number of bytes to receive.
/// \param [in]  source  Rank in MPI_COMM_WORLD from which to receive.
/// \return Number of bytes received.
int Parallel::sendReceive(void* sendBuf, int sendLen, int dest, void* recvBuf, int recvLen, int source)
{
#ifdef DO_MPI
   int bytesReceived;
   MPI_Status status;
   MPI_Sendrecv(sendBuf, sendLen, MPI_BYTE, dest,   0,
                recvBuf, recvLen, MPI_BYTE, source, 0,
                MPI_COMM_WORLD, &status);
   MPI_Get_count(&status, MPI_BYTE, &bytesReceived);

   return bytesReceived;
#else
   assert(source == dest);
   memcpy(recvBuf, sendBuf, sendLen);

   return sendLen;
#endif
}

void Parallel::minRankDouble(RankReduceData* sendBuf, RankReduceData* recvBuf, int count)
{
#ifdef DO_MPI
   MPI_Allreduce(sendBuf, recvBuf, count, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
#else
   for (int ii=0; ii<count; ++ii)
   {
      recvBuf[ii].val = sendBuf[ii].val;
      recvBuf[ii].rank = sendBuf[ii].rank;
   }
#endif
}

void Parallel::maxRankDouble(RankReduceData* sendBuf, RankReduceData* recvBuf, int count)
{
#ifdef DO_MPI
   MPI_Allreduce(sendBuf, recvBuf, count, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
#else
   for (int ii=0; ii<count; ++ii)
   {
      recvBuf[ii].val = sendBuf[ii].val;
      recvBuf[ii].rank = sendBuf[ii].rank;
   }
#endif
}

/// \param [in] count Length of buf in bytes.
void Parallel::bcast(void *buf, int count)
{
#ifdef DO_MPI
   MPI_Bcast(buf, count, MPI_BYTE, root, MPI_COMM_WORLD);
#endif
}

bool Parallel::builtWithMpi()
{
#ifdef DO_MPI
   return 1;
#else
   return 0;
#endif
}
