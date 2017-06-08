#include "Command.hpp"
#include "Memory.hpp"
#include "Parallel.hpp"

#include <getopt.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>

int
addArg(const char *longOption, const char shortOption, bool has_arg,
       const char type, void *dataPtr, int dataSize, const char *help);

void
processArgs(int argc, char **argv);

void
printArgs(void);

void
freeArgs(void);

Command::Command(int argc, char *argv[]) {
  bool help = 0;
  // add arguments for processing.  Please update the html documentation too!
  addArg("help", 'h', 0, 'i', &(help), 0, "print this message");
  addArg("potDir", 'd', 1, 's', potDir, sizeof(potDir), "potential directory");
  addArg("potName", 'p', 1, 's', potName, sizeof(potName), "potential name");
  addArg("potType", 't', 1, 's', potType, sizeof(potType),
         "potential type (funcfl or setfl)");
  addArg("doeam", 'e', 0, 'i', &(doeam), 0, "compute eam potentials");
  addArg("nx", 'x', 1, 'i', &(nx), 0, "number of unit cells in x");
  addArg("ny", 'y', 1, 'i', &(ny), 0, "number of unit cells in y");
  addArg("nz", 'z', 1, 'i', &(nz), 0, "number of unit cells in z");
  addArg("xproc", 'i', 1, 'i', &(xproc), 0, "processors in x direction");
  addArg("yproc", 'j', 1, 'i', &(yproc), 0, "processors in y direction");
  addArg("zproc", 'k', 1, 'i', &(zproc), 0, "processors in z direction");
  addArg("nSteps", 'N', 1, 'i', &(nSteps), 0, "number of time steps");
  addArg("printRate", 'n', 1, 'i', &(printRate), 0,
         "number of steps between output");
  addArg("dt", 'D', 1, 'd', &(dt), 0, "time step (in fs)");
  addArg("lat", 'l', 1, 'd', &(lat), 0, "lattice parameter (Angstroms)");
  addArg("temp", 'T', 1, 'd', &(temperature), 0, "initial temperature (K)");
  addArg("delta", 'r', 1, 'd', &(initialDelta), 0, "initial delta (Angstroms)");

  processArgs(argc, argv);

  // If user didn't set potName, set type dependent default.

  if (strlen(potName) == 0) {
    if (strcmp(potType, "setfl") == 0)
      strncpy(potName, "Cu01.eam.alloy", 1024);
    if (strcmp(potType, "funcfl") == 0)
      strncpy(potName, "Cu_u6.eam", 1024);
  }

  if (help) {
    printArgs();
    freeArgs();
    exit(2);
  }
  freeArgs();
}

void
Command::printYAML(FILE *file) {
  if (!Parallel::printRank())
    return;
  fprintf(file, "Command Line Parameters:\n"
                "  doeam: %d\n"
                "  potDir: %s\n"
                "  potName: %s\n"
                "  potType: %s\n"
                "  nx: %d\n"
                "  ny: %d\n"
                "  nz: %d\n"
                "  xproc: %d\n"
                "  yproc: %d\n"
                "  zproc: %d\n"
                "  Lattice constant: %g Angstroms\n"
                "  nSteps: %d\n"
                "  printRate: %d\n"
                "  Time step: %g fs\n"
                "  Initial Temperature: %g K\n"
                "  Initial Delta: %g Angstroms\n"
                "\n",
          doeam, potDir, potName, potType, nx, ny, nz, xproc, yproc, zproc, lat,
          nSteps, printRate, dt, temperature, initialDelta);
  fflush(file);
}

struct MyOption {
  char *help;
  char *longArg;
  unsigned char shortArg[2];
  int argFlag;
  char type;
  int sz;
  void *ptr;
  void *next;
};

inline MyOption *
nextOption(MyOption *o) {
  return (MyOption *)o->next;
}

static size_t longest   = 1;
static MyOption *myargs = NULL;

static MyOption *
myOptionAlloc(const char *longOption, const char shortOption, int has_arg,
              const char type, void *dataPtr, int dataSize, const char *help) {
  static int iBase = 129;
  MyOption *o      = comdCalloc<MyOption>(1, 1);
  o->help          = strdup(help);
  o->longArg       = strdup(longOption);
  if (shortOption)
    o->shortArg[0] = (unsigned char)shortOption;
  else {
    o->shortArg[0] = iBase;
    iBase++;
  }
  o->argFlag = has_arg;
  o->type    = type;
  o->ptr     = dataPtr;
  o->sz      = dataSize;
  if (longOption)
    longest = std::max(longest, strlen(longOption));
  return o;
}

static MyOption *
myOptionFree(MyOption *o) {
  MyOption *r;
  if (!o)
    return NULL;
  r = nextOption(o);
  if (o->longArg)
    free(o->longArg);
  if (o->help)
    free(o->help);
  free(o);
  return r;
}

static MyOption *
lastOption(MyOption *o) {
  if (!o)
    return o;
  while (nextOption(o))
    o = nextOption(o);
  return o;
}

static MyOption *
findOption(MyOption *o, unsigned char shortArg) {
  while (o) {
    if (o->shortArg[0] == shortArg)
      return o;
    o = nextOption(o);
  }
  return o;
}

int
addArg(const char *longOption, const char shortOption, bool has_arg,
       const char type, void *dataPtr, int dataSize, const char *help) {
  MyOption *o;
  MyOption *p;
  o = myOptionAlloc(longOption, shortOption, has_arg, type, dataPtr, dataSize,
                    help);
  if (!o)
    return 1;
  if (!myargs)
    myargs = o;
  else {
    p       = lastOption(myargs);
    p->next = (void *)o;
  }
  return 0;
}

void
freeArgs() {
  while (myargs) {
    myargs = myOptionFree(myargs);
  }
  return;
}

void
printArgs() {
  MyOption *o = myargs;
  char s[4096];
  unsigned char *shortArg;
  fprintf(screenOut, "\n  Arguments are: \n");
  sprintf(s, "   --%%-%lus", longest);
  while (o) {
    if (o->shortArg[0] < 0xFF)
      shortArg = o->shortArg;
    else
      shortArg = (unsigned char *)"---";
    fprintf(screenOut, s, o->longArg);
    fprintf(screenOut, " -%c  arg=%1d type=%c  %s\n", shortArg[0], o->argFlag,
            o->type, o->help);
    o = nextOption(o);
  }
  fprintf(screenOut, "\n\n");
  return;
}

void
processArgs(int argc, char **argv) {
  MyOption *o;
  int n = 0;
  int i;
  struct option *opts;
  char *sArgs;
  int c;

  if (!myargs)
    return;
  o = myargs;
  while (o) {
    n++, o = nextOption(o);
  }

  o     = myargs;
  sArgs = comdCalloc<char>(2 * (n + 2), 1);
  opts  = comdCalloc<option>(n, 1);
  for (i = 0; i < n; i++) {
    opts[i].name    = o->longArg;
    opts[i].has_arg = o->argFlag;
    opts[i].flag    = 0;
    opts[i].val     = o->shortArg[0];

    strcat(sArgs, (char *)o->shortArg);
    if (o->argFlag)
      strcat(sArgs, ":");
    o = nextOption(o);
  }

  while (1) {
    int option_index = 0;

    c = getopt_long(argc, argv, sArgs, opts, &option_index);
    if (c == -1)
      break;
    o = findOption(myargs, c);
    if (!o) {
      fprintf(screenOut, "\n\n"
                         "    invalid switch : -%c in getopt()\n"
                         "\n\n",
              c);
      break;
    }
    if (!o->argFlag) {
      int *ii = static_cast<int *>(o->ptr);
      *ii     = 1;
    } else {
      switch (o->type) {
      case 'i': sscanf(optarg, "%d", (int *)o->ptr); break;
      case 'f': sscanf(optarg, "%f", (float *)o->ptr); break;
      case 'd': sscanf(optarg, "%lf", (double *)o->ptr); break;
      case 's':
        strncpy((char *)o->ptr, optarg, o->sz);
        ((char *)o->ptr)[o->sz - 1] = '\0';
        break;
      case 'c': sscanf(optarg, "%c", (char *)o->ptr); break;
      default:
        fprintf(screenOut,
                "\n\n"
                "    invalid type : %c in getopt()\n"
                "    valid values are 'e', 'z'. 'i','d','f','s', and 'c'\n"
                "\n\n",
                c);
      }
    }
  }

  free(opts);
  free(sArgs);

  return;
}
