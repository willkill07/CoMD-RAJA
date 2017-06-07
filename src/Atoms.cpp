#include "Atoms.hpp"
#include "LinkCell.hpp"

#include "Memory.hpp"

Atoms::Atoms(LinkCell const &boxes) {
    int maxTotalAtoms = MAXATOMS * boxes.nTotalBoxes;
    gid = comdMalloc<int>(maxTotalAtoms);
    iSpecies = comdMalloc<int>(maxTotalAtoms);
    r = comdMalloc<real3>(maxTotalAtoms);
    p = comdMalloc<real3>(maxTotalAtoms);
    f = comdMalloc<real3>(maxTotalAtoms);
    U = comdMalloc<real_t>(maxTotalAtoms);

    nLocal = 0;
    nGlobal = 0;

    for (int iOff = 0; iOff < maxTotalAtoms; iOff++) {
        gid[iOff] = 0;
        iSpecies[iOff] = 0;
        zeroArray(r[iOff]);
        zeroArray(p[iOff]);
        zeroArray(f[iOff]);
        U[iOff] = 0.;
    }
}

Atoms::~Atoms() {
    comdFree(gid);
    comdFree(iSpecies);
    comdFree(r);
    comdFree(p);
    comdFree(f);
    comdFree(U);
}