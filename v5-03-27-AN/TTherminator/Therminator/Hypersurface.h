/***************************************************************************
 *                             L H Y Q U I D                               *
 *              reLativistic HYdrodynamics of Quark glUon fluID            *
 *                                                                         *
 *                         T H E R M I N A T O R                           *
 *                      THERMal heavy-IoN generATOR                        *
 *                                                                         *
 *                              INTERFACE                                  *
 *                                                                         *
 * Author of LHYQUID THERMINATOR INTERFACE code:                           *
 *   Miko�aj Chojnacki, e-mail: Mikolaj.Chojnacki@ifj.edu.pl               *
 *                                                                         *
 * About the interface:                                                    *
 *   Code in the files "Hypersurface.h" and "Hypersurface.cxx" enables the *
 * usage of a paramatrized hypersurface from LHYQUID (written in           *
 * Mathematica).                                                           *
 * Set the option FreezeOutModel = Lhyquid in the file "therminator.in".   *
 * Modifications in the original files (Integrator.cxx, Integrator.h,      *
 * therm_events, Makefile, therminator.in are marked by "MCH" comment      *
 *                                                                         *
 ***************************************************************************/


#ifndef _MCH_HYPERSURFACE_
#define _MCH_HYPERSURFACE_

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;

class Hypersurface {
  public:
    // constructor
    Hypersurface(const char *dirname);
    Hypersurface(void);
    Hypersurface(const Hypersurface &aSurf);
    Hypersurface& operator=(const Hypersurface &aSurf);
    // destructor
    ~Hypersurface(void);
    // variables
    double TFO;
    double tau0;
    // functions
    double fahs  (double p, double z);
    double fvhs  (double p, double z);
    double fdhs  (double p, double z);
    double fDpdhs(double p, double z);
    double fDzdhs(double p, double z);
  private:
    unsigned int i, j;
    unsigned int Np, Nz;
    double   ip, fp, dp, iz, fz, dz;
    double   **aArr, **vArr, **dArr, **DpdArr, **DzdArr;
    char     *FName;
    char     buff[100];
    ifstream *HSFile;
};


#endif /* _MCH_HYPERSURFACE_ */
