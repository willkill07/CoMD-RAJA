#include "EAMPotential.hpp"
#include "CoMDTypes.hpp"
#include "Constants.hpp"
#include "Memory.hpp"
#include "Parallel.hpp"

#include <algorithm>

#include <math.h>
#include <string.h>

void
fileNotFound(const char *callSite, const char *filename);

void
notAlloyReady(const char *callSite);

void
typeNotSupported(const char *callSite, const char *type);

EAMPotential
make_eam_setfl(const char *dir, const char *potName);
EAMPotential
make_eam_funcfl(const char *dir, const char *potName);

EAMPotential::EAMPotential(const char *dir, const char *file,
                           const char *type) {
  if (Parallel::myRank() == 0) {
    if (strcmp(type, "setfl") == 0)
      *this = make_eam_setfl(dir, file);
    else if (strcmp(type, "funcfl") == 0)
      *this = make_eam_funcfl(dir, file);
    else
      typeNotSupported("initEamPot", type);
  }
  broadcast();
}

void
EAMPotential::prepare(SimFlat<EAMPotential> &s) {
  int maxTotalAtoms         = MAXATOMS * s.boxes.nTotalBoxes;
  dfEmbed                   = comdMalloc<real_t>(maxTotalAtoms);
  rhobar                    = comdMalloc<real_t>(maxTotalAtoms);
  forceExchange             = HaloForceExchange(s.domain, s.boxes);
  forceExchangeData.dfEmbed = dfEmbed;
  forceExchangeData.boxes   = &(s.boxes);
}

int
EAMPotential::force_impl(SimFlat<EAMPotential> &s) {
  real_t rCut2 = cutoff * cutoff;
  real_t etot  = 0.;

  int fsize = s.boxes.nTotalBoxes * MAXATOMS;
#pragma omp parallel for
  for (int ii = 0; ii < fsize; ii++) {
    s.atoms.f[ii][0] = 0.;
    s.atoms.f[ii][1] = 0.;
    s.atoms.f[ii][2] = 0.;
    s.atoms.U[ii]    = 0.;
    dfEmbed[ii]      = 0.;
    rhobar[ii]       = 0.;
  }

  constexpr const int nNbrBoxes = 27;

/*
RAJA::forallN<
  RAJA::NestedPolicy<
    RAJA::ExecPolicy<
      RAJA::omp_for_collapse_nowait_exec,
      RAJA::omp_for_collapse_nowait_exec>,
    RAJA::OMP_Parallel>>(
      RAJA::RangeSegment(0, s.boxes.nLocalBoxes),
      RAJA::RangeSegment(0, nNbrBoxes), [=](int iBox, int jTmp) {
        int jBox = s.boxes.nbrBoxes[iBox][jTmp];
        int nIBox = s.boxes.nAtoms[iBox];
        int nJBox = s.boxes.nAtoms[jBox];
        RAJA::forallN<
          RAJA::NestedPolicy<
            RAJA::ExecPolicy<
              RAJA::simd_exec,
              RAJA::simd_exec>>>(
                RAJA::RangeSegment(MAXATOMS * iBox, MAXATOMS * iBox + nIBox),
                RAJA::RangeSegment(MAXATOMS * jBox, MAXATOMS * jBox + nJBox),
                [=](int iOff, int jOff) {
                  real3 dr = {
                    s.atoms.r[iOff][0] - s.atoms.r[jOff][0],
                    s.atoms.r[iOff][1] - s.atoms.r[jOff][1],
                    s.atoms.r[iOff][2] - s.atoms.r[jOff][2]
                  };
                  real_t r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                  if (r2 <= rCut2 && r2 > 0.0) {
                    real_t r = sqrt(r2);
                    real_t phiTmp, dPhi, rhoTmp, dRho;
                    phi->interpolate(r, &phiTmp, &dPhi);
                    rho->interpolate(r, &rhoTmp, &dRho);
                    s.atoms.f[iOff][0] -= dPhi * dr[0] / r;
                    s.atoms.f[iOff][1] -= dPhi * dr[1] / r;
                    s.atoms.f[iOff][2] -= dPhi * dr[2] / r;
                    s.atoms.U[iOff] += 0.5 * phiTmp;
                    etot += 0.5 * phiTmp;
                    rhobar[iOff] += rhoTmp;
                  }
                });
      });
*/
#pragma omp parallel for reduction(+ : etot)
  for (int iBox = 0; iBox < s.boxes.nLocalBoxes; iBox++) {
    // loop over neighbor boxes of iBox (some may be halo boxes)
    for (int jTmp = 0; jTmp < nNbrBoxes; jTmp++) {
      int jBox  = s.boxes.nbrBoxes[iBox][jTmp];
      int nIBox = s.boxes.nAtoms[iBox];
      int nJBox = s.boxes.nAtoms[jBox];

      // loop over atoms in iBox
      for (int iOff = MAXATOMS * iBox; iOff < (iBox * MAXATOMS + nIBox);
           iOff++) {
        // loop over atoms in jBox
        for (int jOff = MAXATOMS * jBox; jOff < (jBox * MAXATOMS + nJBox);
             jOff++) {
          real3 dr  = {0, 0, 0};
          real_t r2 = 0.0;
          for (int k = 0; k < 3; k++) {
            dr[k] = s.atoms.r[iOff][k] - s.atoms.r[jOff][k];
            r2 += dr[k] * dr[k];
          }

          if (r2 <= rCut2 && r2 > 0.0) {
            real_t r = sqrt(r2);

            real_t phiTmp, dPhi, rhoTmp, dRho;
            phi.interpolate(r, &phiTmp, &dPhi);
            rho.interpolate(r, &rhoTmp, &dRho);
            s.atoms.f[iOff][0] -= dPhi * dr[0] / r;
            s.atoms.f[iOff][1] -= dPhi * dr[1] / r;
            s.atoms.f[iOff][2] -= dPhi * dr[2] / r;
            s.atoms.U[iOff] += 0.5 * phiTmp;
            etot += 0.5 * phiTmp;
            rhobar[iOff] += rhoTmp;
          }
        } // loop over atoms in jBox
      }   // loop over atoms in iBox
    }     // loop over neighbor boxes
  }       // loop over local boxes

// Compute Embedding Energy
// loop over all local boxes
/*{
  RAJA::ReduceSum<real_t, ReducePol> etotal;
  RAJA::forall<OuterExecPol>(0, s.boxes.nLocalBoxes, [=](int iBox) {
    RAJA::forall<InnerExecPol>(0, s.boxes.nAtoms[iBox], [=](int i) {
      int iOff = (MAXATOMS * iBox) + i;
      real_t fEmbed, dfEmbed;
      f->interpolate(rhobar[iOff], &fEmbed, &dfEmbed);
      dfEmbed[iOff] = dfEmbed; // save derivative for halo exchange
      s.atoms.U[iOff] += fEmbed;
      etotal += fEmbed;
    });
  });
  etot = etotal.get();
}*/

#pragma omp parallel for reduction(+ : etot)
  for (int iBox = 0; iBox < s.boxes.nLocalBoxes; iBox++) {
    int nIBox = s.boxes.nAtoms[iBox];

    // loop over atoms in iBox
    for (int iOff = MAXATOMS * iBox; iOff < (MAXATOMS * iBox + nIBox); iOff++) {
      real_t fEmbed, dfembed;
      f.interpolate(rhobar[iOff], &fEmbed, &dfembed);
      // interpolate(f, rhobar[iOff], &fEmbed, &dfEmbed);
      dfEmbed[iOff] = dfembed; // save derivative for halo exchange
      s.atoms.U[iOff] += fEmbed;
      etot += fEmbed;
    }
  }

  // exchange derivative of the embedding energy with repsect to rhobar
  startTimer(eamHaloTimer);
  forceExchange.exchange(forceExchangeData);
  stopTimer(eamHaloTimer);

// third pass
// loop over local boxes
#pragma omp parallel for
  for (int iBox = 0; iBox < s.boxes.nLocalBoxes; iBox++) {
    int nIBox = s.boxes.nAtoms[iBox];

    // loop over neighbor boxes of iBox (some may be halo boxes)
    for (int jTmp = 0; jTmp < nNbrBoxes; jTmp++) {
      int jBox  = s.boxes.nbrBoxes[iBox][jTmp];
      int nJBox = s.boxes.nAtoms[jBox];

      // loop over atoms in iBox
      for (int iOff = MAXATOMS * iBox; iOff < (MAXATOMS * iBox + nIBox);
           iOff++) {
        // loop over atoms in jBox
        for (int jOff = MAXATOMS * jBox; jOff < (MAXATOMS * jBox + nJBox);
             jOff++) {
          real_t r2 = 0.0;
          real3 dr  = {0, 0, 0};
          for (int k = 0; k < 3; k++) {
            dr[k] = s.atoms.r[iOff][k] - s.atoms.r[jOff][k];
            r2 += dr[k] * dr[k];
          }

          if (r2 <= rCut2 && r2 > 0.0) {
            real_t r = sqrt(r2);

            real_t rhoTmp, dRho;
            rho.interpolate(r, &rhoTmp, &dRho);

            for (int k = 0; k < 3; k++) {
              s.atoms.f[iOff][k] -=
                  (dfEmbed[iOff] + dfEmbed[jOff]) * dRho * dr[k] / r;
            }
          }

        } // loop over atoms in jBox
      }   // loop over atoms in iBox
    }     // loop over neighbor boxes
  }       // loop over local boxes

  s.ePotential = (real_t)etot;

  return 0;
}

void
EAMPotential::print_impl(FILE *file) {
  fprintf(file, "  Potential type  : EAM\n");
  fprintf(file, "  Species name    : %s\n", name);
  fprintf(file, "  Atomic number   : %d\n", atomicNo);
  fprintf(file, "  Mass            : " FMT1 " amu\n",
          mass / amuToInternalMass); // print in amu
  fprintf(file, "  Lattice type    : %s\n", latticeType);
  fprintf(file, "  Lattice spacing : " FMT1 " Angstroms\n", lat);
  fprintf(file, "  Cutoff          : " FMT1 " Angstroms\n", cutoff);
}

void
EAMPotential::broadcast() {
  struct Buffer {
    real_t cutoff, mass, lat;
    int atomicNo;
    char latticeType[8];
    char name[3];
  };

  Buffer buf;
  if (Parallel::myRank() == 0) {
    buf.cutoff   = cutoff;
    buf.mass     = mass;
    buf.lat      = lat;
    buf.atomicNo = atomicNo;
    strcpy(buf.latticeType, latticeType);
    strcpy(buf.name, name);
  }
  Parallel::bcast(&buf, sizeof(buf));
  cutoff   = buf.cutoff;
  mass     = buf.mass;
  lat      = buf.lat;
  atomicNo = buf.atomicNo;
  strcpy(latticeType, buf.latticeType);
  strcpy(name, buf.name);

  phi.broadcast();
  rho.broadcast();
  f.broadcast();
}

void
InterpolationObject::broadcast() {
  struct Buffer {
    real_t x0, invDx;
    int n;
  };
  Buffer buf;

  if (Parallel::myRank() == 0) {
    buf.n     = n;
    buf.x0    = x0;
    buf.invDx = invDx;
  }
  Parallel::bcast(&buf, sizeof(Buffer));

  if (Parallel::myRank() != 0) {
    n      = buf.n;
    x0     = buf.x0;
    invDx  = buf.invDx;
    values = comdMalloc<real_t>(buf.n + 3);
    values++;
  }

  int valuesSize = sizeof(real_t) * (n + 3);
  Parallel::bcast(values - 1, valuesSize);
}

EAMPotential
make_eam_setfl(const char *dir, const char *potName) {
  EAMPotential pot;
  char tmp[4096];
  sprintf(tmp, "%s/%s", dir, potName);

  FILE *potFile = fopen(tmp, "r");
  if (potFile == nullptr)
    fileNotFound("eamReadSetfl", tmp);

  // read the first 3 lines (comments)
  fgets(tmp, sizeof(tmp), potFile);
  fgets(tmp, sizeof(tmp), potFile);
  fgets(tmp, sizeof(tmp), potFile);

  // line 4
  fgets(tmp, sizeof(tmp), potFile);
  int nElems;
  sscanf(tmp, "%d", &nElems);
  if (nElems != 1)
    notAlloyReady("eamReadSetfl");

  // line 5
  int nRho, nR;
  double dRho, dR, cutoff;
  //  The same cutoff is used by all alloys, NB: cutoff = nR * dR is redundant
  fgets(tmp, sizeof(tmp), potFile);
  sscanf(tmp, "%d %le %d %le %le", &nRho, &dRho, &nR, &dR, &cutoff);
  pot.cutoff = cutoff;

  // **** THIS CODE IS RESTRICTED TO ONE ELEMENT
  // Per-atom header
  fgets(tmp, sizeof(tmp), potFile);
  int nAtomic;
  double mass, lat;
  char latticeType[8];
  sscanf(tmp, "%d %le %le %s", &nAtomic, &mass, &lat, latticeType);
  pot.atomicNo = nAtomic;
  pot.lat      = lat;
  pot.mass     = mass * amuToInternalMass; // file has mass in AMU.
  strncpy(pot.latticeType, latticeType, 3);

  // allocate read buffer
  int bufSize = std::max(nRho, nR);
  real_t *buf = comdMalloc<real_t>(bufSize);
  real_t x0   = 0.0;

  // Read embedding energy F(rhobar)
  for (int ii = 0; ii < nRho; ++ii)
    fscanf(potFile, FMT1, buf + ii);
  pot.f = InterpolationObject(nRho, x0, dRho, buf);

  // Read electron density rho(r)
  for (int ii = 0; ii < nR; ++ii)
    fscanf(potFile, FMT1, buf + ii);
  pot.rho = InterpolationObject(nR, x0, dR, buf);

  // Read phi(r)*r and convert to phi(r)
  for (int ii = 0; ii < nR; ++ii)
    fscanf(potFile, FMT1, buf + ii);
  for (int ii = 1; ii < nR; ++ii) {
    real_t r = x0 + ii * dR;
    buf[ii] /= r;
  }

  buf[0]  = buf[1] + (buf[1] - buf[2]); // Linear interpolation to get phi[0].
  pot.phi = InterpolationObject(nR, x0, dR, buf);
  comdFree(buf);
  return pot;
}

EAMPotential
make_eam_funcfl(const char *dir, const char *potName) {
  char tmp[4096];
  EAMPotential pot;

  sprintf(tmp, "%s/%s", dir, potName);
  FILE *potFile = fopen(tmp, "r");
  if (potFile == NULL)
    fileNotFound("eamReadFuncfl", tmp);

  // line 1
  fgets(tmp, sizeof(tmp), potFile);
  char name[3];
  sscanf(tmp, "%s", name);
  strcpy(pot.name, name);

  // line 2
  int nAtomic;
  double mass, lat;
  char latticeType[8];
  fgets(tmp, sizeof(tmp), potFile);
  sscanf(tmp, "%d %le %le %s", &nAtomic, &mass, &lat, latticeType);
  pot.atomicNo = nAtomic;
  pot.lat      = lat;
  pot.mass     = mass * amuToInternalMass; // file has mass in AMU.
  strcpy(pot.latticeType, latticeType);

  // line 3
  int nRho, nR;
  double dRho, dR, cutoff;
  fgets(tmp, sizeof(tmp), potFile);
  sscanf(tmp, "%d %le %d %le %le", &nRho, &dRho, &nR, &dR, &cutoff);
  pot.cutoff = cutoff;
  real_t x0  = 0.0; // tables start at zero.

  // allocate read buffer
  int bufSize = std::max(nRho, nR);
  real_t *buf = comdMalloc<real_t>(bufSize);

  // read embedding energy
  for (int ii = 0; ii < nRho; ++ii)
    fscanf(potFile, FMT1, buf + ii);
  pot.f = InterpolationObject(nRho, x0, dRho, buf);

  // read Z(r) and convert to phi(r)
  for (int ii = 0; ii < nR; ++ii)
    fscanf(potFile, FMT1, buf + ii);
  for (int ii = 1; ii < nR; ++ii) {
    real_t r = x0 + ii * dR;
    buf[ii] *= buf[ii] / r;
    buf[ii] *= hartreeToEv * bohrToAngs; // convert to eV
  }
  buf[0]  = buf[1] + (buf[1] - buf[2]); // linear interpolation to get phi[0].
  pot.phi = InterpolationObject(nR, x0, dR, buf);

  // read electron density rho
  for (int ii = 0; ii < nR; ++ii)
    fscanf(potFile, FMT1, buf + ii);
  pot.rho = InterpolationObject(nR, x0, dR, buf);

  comdFree(buf);
}

void
fileNotFound(const char *callSite, const char *filename) {
  fprintf(screenOut, "%s: Can't open file %s.  Fatal Error.\n", callSite,
          filename);
  exit(-1);
}

void
notAlloyReady(const char *callSite) {
  fprintf(screenOut,
          "%s: CoMD 1.1 does not support alloys and cannot\n"
          "   read setfl files with multiple species.  Fatal Error.\n",
          callSite);
  exit(-1);
}

void
typeNotSupported(const char *callSite, const char *type) {
  fprintf(screenOut, "%s: Potential type %s not supported. Fatal Error.\n",
          callSite, type);
  exit(-1);
}

InterpolationObject::InterpolationObject(int n, real_t x0, real_t dx,
                                         real_t *data)
    : values(comdCalloc<real_t>(1, n + 3)), x0(x0), invDx(1.0 / dx), n(n) {
  values[0] = data[0];

  ++values;
  for (int i  = 0; i < n; ++i)
    values[i] = data[i];

  real_t last   = data[n - 1];
  values[n]     = last;
  values[n + 1] = last;
}

InterpolationObject::~InterpolationObject() {
  --values;
  comdFree(values);
}

void
InterpolationObject::interpolate(real_t r, real_t *f, real_t *df) {
  const real_t *tt = values; // alias
  if (r < x0)
    r    = x0;
  r      = (r - x0) * (invDx);
  int ii = (int)floor(r);
  if (ii > n) {
    ii = n;
    r  = n / invDx;
  }
  // reset r to fractional distance
  r         = r - floor(r);
  real_t g1 = tt[ii + 1] - tt[ii - 1];
  real_t g2 = tt[ii + 2] - tt[ii];
  *f  = tt[ii] + 0.5 * r * (g1 + r * (tt[ii + 1] + tt[ii - 1] - 2.0 * tt[ii]));
  *df = 0.5 * (g1 + r * (g2 - g1)) * invDx;
}

void
InterpolationObject::print(const char *fileName) {
  if (!Parallel::printRank())
    return;
  FILE *potData;
  potData   = fopen(fileName, "w");
  real_t dR = 1.0 / invDx;
  for (int i = 0; i < n; i++) {
    real_t r = x0 + i * dR;
    fprintf(potData, "%d %e %e\n", i, r, values[i]);
  }
  fclose(potData);
}

template <>
EAMPotential
make_potential<EAMPotential>(Command const &cmd) {
  return {cmd.potDir, cmd.potName, cmd.potType};
}
