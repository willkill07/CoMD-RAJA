/// \file
/// Handle command line arguments.



#include <cstring>
#include <getopt.h>

#include "Command.hpp"
#include "Memory.hpp"
#include "Parallel.hpp"

#include <algorithm>

/// \page pg_running_comd Running CoMD
///
/// \section sec_command_line_options Command Line Options
///
/// CoMD accepts a number of command line options to set the parameters
/// of the simulation.  Every option has both a long form and a short
/// form.  The long and short form of the arguments are entirely
/// interchangeable and may be mixed. All the arguments are independent
/// with the exception of the \--potDir, \--potName, and \--potType,
/// (short forms -d, -n, and -t) arguments which are only relevant when
/// used in conjunction with \--doeam, (-e).
///
/// Supported options are:
///
/// | Long  Form    | Short Form  | Default Value | Description
/// | :------------ | :---------: | :-----------: | :----------
/// | \--help       | -h          | N/A           | print this message
/// | \--potDir     | -d          | pots          | potential directory
/// | \--potName    | -p          | Cu_u6.eam     | potential name
/// | \--potType    | -t          | funcfl        | potential type (funcfl or setfl)
/// | \--doeam      | -e          | N/A           | compute eam potentials (default is LJ)
/// | \--nx         | -x          | 20            | number of unit cells in x
/// | \--ny         | -y          | 20            | number of unit cells in y
/// | \--nz         | -z          | 20            | number of unit cells in z
/// | \--xproc      | -i          | 1             | number of ranks in x direction
/// | \--yproc      | -j          | 1             | number of ranks in y direction
/// | \--zproc      | -k          | 1             | number of ranks in z direction
/// | \--nSteps     | -N          | 100           | total number of time steps
/// | \--printRate  | -n          | 10            | number of steps between output
/// | \--dt         | -D          | 1             | time step (in fs)
/// | \--lat        | -l          | -1            | lattice parameter (Angstroms)
/// | \--temp       | -T          | 600           | initial temperature (K)
/// | \--delta      | -r          | 0             | initial delta (Angstroms)
///
/// Notes:
///
/// The negative value for the lattice parameter (such as the default
/// value, -1) is interpreted as a flag to indicate that the lattice
/// parameter should be set from the potential. All supplied potentials
/// are for copper and have a lattice constant of 3.615
/// Angstroms. Setting the lattice parameter to any positive value will
/// override the values provided in the potential files.
///
/// The default potential name for the funcfl potential type is
/// Cu_u6.eam (Adams potential).  For the setfl type the default
/// potential name is Cu01.eam.alloy (Mishin potential).  Although these
/// will yield similar dynamics, the table have a very different number
/// of entries (500 vs. 10,000 points, respectively) This may give very
/// different performance, depending on the hardware.
///
/// The default temperature is 600K.  However, when using a perfect
/// lattice the system will rapidly cool to 300K due to equipartition of
/// energy.
///
///
/// \subsection ssec_example_command_lines Examples
///
/// All of the examples below assume:
/// - The current working directory contains a copy of the pots dir (or
///   a link to it).
/// - The CoMD bin directory is located in ../bin
///
/// Running in the examples directory will satisfy these requirements.
///
/// ------------------------------
///
/// The canonical base simulation, is
///
///     $ mpirun -np 1 ../bin/CoMD-mpi
///
/// Or, if the code was built without MPI:
///
///     $ ../bin/CoMD-serial
///
/// ------------------------------
///
/// \subsubsection cmd_examples_potential Changing Potentials
///
/// To run with the default (Adams) EAM potential, specify -e:
///
///     $ ../bin/CoMD-mpi -e
///
/// ------------------------------
///
/// To run using the Mishin EAM potential contained in the setfl file
/// Cu01.eam.alloy. This potential uses much larger tables (10,000
/// entries vs. 500 for the Adams potential).
///
///     $ ../bin/CoMD-mpi -e -t setfl
///
/// ------------------------------
///
/// Selecting the name of a setfl file without setting the appropriate
/// potential type
///
///     $ ../bin/CoMD-mpi -e -p Cu01.eam.alloy
///
/// will result in an error message:
///
/// Only FCC Lattice type supported, not . Fatal Error.
///
/// Instead use:
///
///     $ ../bin/CoMD-mpi -e -t setfl -p Cu01.eam.alloy
///
/// ------------------------------
///
/// \subsubsection cmd_example_struct Initial Structure Modifications
///
/// To change the lattice constant and run with an expanded or
/// compressed lattice:
///
///     $ ../bin/CoMD-mpi -l 3.5
///
/// This can be useful to test that the potential is being correctly
/// evaluated as a function of interatomic spacing (the cold
/// curve). However, due to the high degree of symmetry of a perfect
/// lattice, this type of test is unlikely to detect errors in the force
/// computation.
///
/// ------------------------------
///
/// Initialize with zero temperature (zero instantaneous particle
/// velocity) but with a random displacements of the atoms (in this
/// case the maximum displacement is 0.1 Angstrom along each axis).
///
///      $ ../bin/CoMD-mpi --delta 0.1 -T 0
///
/// Typical values of delta are in the range of 0.1 to 0.5 Angstroms.
/// Larger values of delta correspond to higher initial potential energy
/// which in turn produce higer temperatures as the structure
/// equilibrates.
///
/// ------------------------------
///
///
/// \subsubsection cmd_examples_scaling Scaling Examples
///
/// Simple shell scripts that demonstrate weak and strong scaling
/// studies are provided in the examples directory.
///
/// ------------------------------
///
/// Run the default global simulation size (32,000 atoms) distributed
/// over 8 cubic subdomains, an example of strong scaling.  If the
/// number of processors does not equal (i*j*k) the run will abort.
/// Notice that spaces are optional between short form options and their
/// arguments.
///
///     $ mpirun -np 8 ../bin/CoMD-mpi -i2 -j2 -k2
///
/// ------------------------------
///
/// Run a weak scaling example: the simulation is doubled in each
/// dimension from the default 20 x 20 x 20 and the number of subdomains
/// in each direction is also doubled.
///
///     $ mpirun -np 8 ../bin/CoMD-mpi -i2 -j2 -k2 -x 40 -y 40 -z 40
///
/// ------------------------------
///
/// The same weak scaling run, but for 10,000 timesteps, with output
/// only every 100 steps.
///
///     $ mpirun -np 8 ../bin/CoMD-mpi -i2 -j2 -k2 -x 40 -y 40 -z 40 -N 10000 -n 100
///

/// Specifies a command line argument that should be accepted by the program.
/// \param [in]  longOption  The long name of option i.e., --optionname
/// \param [in]  shortOption The short name of option i.e., -o
/// \param [in]  has_arg  Whether this option has an argument i.e., -o value.
///                       If has_arg is 0, then dataPtr must be an integer
///                       pointer.
/// \param [in]  type  The type of the argument. Valid values are:
///                    -  i   integer
///                    -  f   float
///                    -  d   double
///                    -  s   string
///                    -  c   character
///
/// \param [in]  dataPtr  A pointer to where the value will be stored.
/// \param [in]  dataSize The length of dataPtr, only useful for character
///                       strings.
/// \param [in]  help     A short help string, preferably a single line or
///                       less.
int addArg(const char *longOption, const char shortOption,
           bool has_arg, const char type, void *dataPtr, int dataSize,
           const char *help);

/// Call this to process your arguments.
void processArgs(int argc, char **argv);

/// Prints the arguments to the stdout stream.
void printArgs(void);

void freeArgs(void);

/// \details Initialize a Command structure with default values, then
/// parse any command line arguments that were supplied to overwrite
/// defaults.
///
/// \param [in] argc the number of command line arguments
/// \param [in] argv the command line arguments array
Command::Command(int argc, char *argv[]) {
    bool help = 0;
    // add arguments for processing.  Please update the html documentation too!
    addArg("help", 'h', 0, 'i', &(help), 0, "print this message");
    addArg("potDir", 'd', 1, 's', potDir, sizeof(potDir), "potential directory");
    addArg("potName", 'p', 1, 's', potName, sizeof(potName), "potential name");
    addArg("potType", 't', 1, 's', potType, sizeof(potType), "potential type (funcfl or setfl)");
    addArg("doeam", 'e', 0, 'i', &(doeam), 0, "compute eam potentials");
    addArg("nx", 'x', 1, 'i', &(nx), 0, "number of unit cells in x");
    addArg("ny", 'y', 1, 'i', &(ny), 0, "number of unit cells in y");
    addArg("nz", 'z', 1, 'i', &(nz), 0, "number of unit cells in z");
    addArg("xproc", 'i', 1, 'i', &(xproc), 0, "processors in x direction");
    addArg("yproc", 'j', 1, 'i', &(yproc), 0, "processors in y direction");
    addArg("zproc", 'k', 1, 'i', &(zproc), 0, "processors in z direction");
    addArg("nSteps", 'N', 1, 'i', &(nSteps), 0, "number of time steps");
    addArg("printRate", 'n', 1, 'i', &(printRate), 0, "number of steps between output");
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

void Command::printYAML(FILE *file) {
    if (!Parallel::printRank())
        return;
    fprintf(file,
            "Command Line Parameters:\n"
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
            doeam,
            potDir,
            potName,
            potType,
            nx, ny, nz,
            xproc, yproc, zproc,
            lat,
            nSteps,
            printRate,
            dt,
            temperature,
            initialDelta
    );
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

inline MyOption *nextOption(MyOption *o) { return (MyOption *) o->next; }

static size_t longest = 1;
static MyOption *myargs = NULL;

static MyOption *myOptionAlloc(const char *longOption, const char shortOption,
                               int has_arg, const char type, void *dataPtr,
                               int dataSize, const char *help) {
    static int iBase = 129;
    MyOption *o = comdCalloc<MyOption>(1, 1);
    o->help = strdup(help);
    o->longArg = strdup(longOption);
    if (shortOption)
        o->shortArg[0] = (unsigned char) shortOption;
    else {
        o->shortArg[0] = iBase;
        iBase++;
    }
    o->argFlag = has_arg;
    o->type = type;
    o->ptr = dataPtr;
    o->sz = dataSize;
    if (longOption)
        longest = std::max(longest, strlen(longOption));
    return o;
}

static MyOption *myOptionFree(MyOption *o) {
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

static MyOption *lastOption(MyOption *o) {
    if (!o)
        return o;
    while (nextOption(o))
        o = nextOption(o);
    return o;
}

static MyOption *findOption(MyOption *o, unsigned char shortArg) {
    while (o) {
        if (o->shortArg[0] == shortArg)
            return o;
        o = nextOption(o);
    }
    return o;
}

int addArg(const char *longOption, const char shortOption, bool has_arg,
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
        p = lastOption(myargs);
        p->next = (void *) o;
    }
    return 0;
}

void freeArgs() {
    while (myargs) {
        myargs = myOptionFree(myargs);
    }
    return;
}

void printArgs() {
    MyOption *o = myargs;
    char s[4096];
    unsigned char *shortArg;
    fprintf(screenOut, "\n"
            "  Arguments are: \n");
    sprintf(s, "   --%%-%lus", longest);
    while (o) {
        if (o->shortArg[0] < 0xFF)
            shortArg = o->shortArg;
        else
            shortArg = (unsigned char *) "---";
        fprintf(screenOut, s, o->longArg);
        fprintf(screenOut, " -%c  arg=%1d type=%c  %s\n", shortArg[0], o->argFlag,
                o->type, o->help);
        o = nextOption(o);
    }
    fprintf(screenOut, "\n\n");
    return;
}

void processArgs(int argc, char **argv) {
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

    o = myargs;
    sArgs = comdCalloc<char>(2 * (n + 2), 1);
    opts = comdCalloc<option>(n, 1);
    for (i = 0; i < n; i++) {
        opts[i].name = o->longArg;
        opts[i].has_arg = o->argFlag;
        opts[i].flag = 0;
        opts[i].val = o->shortArg[0];

        strcat(sArgs, (char *) o->shortArg);
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
            *ii = 1;
        } else {
            switch (o->type) {
                case 'i':
                    sscanf(optarg, "%d", (int *) o->ptr);
                    break;
                case 'f':
                    sscanf(optarg, "%f", (float *) o->ptr);
                    break;
                case 'd':
                    sscanf(optarg, "%lf", (double *) o->ptr);
                    break;
                case 's':
                    strncpy((char *) o->ptr, optarg, o->sz);
                    ((char *) o->ptr)[o->sz - 1] = '\0';
                    break;
                case 'c':
                    sscanf(optarg, "%c", (char *) o->ptr);
                    break;
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
