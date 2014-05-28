#ifndef _f_FilHep_h_included_
#define _f_FilHep_h_included_

/**
 * This file contains an interface between TAUOLA FORTRAN routines 
 * and the C++ event. These methods, are defined so they can be called
 * by tauola (in tauola.f). They will never be called by the C++ code. 
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "TauolaParticle.h"
#include "DecayList.h"

using namespace std;

namespace Tauolapp
{

/** Fill a particle into the TauolaEvent. This relies heavily
    on the static data structure DecayList */
extern "C" void filhep_(int * n, int * status, int * pdg_id,
                        int * mother_first, int * mother_last, 
                        int * daughter_first, int * daughter_last, 
                        float p4[4], float * p_inv_mass, bool * photos_flag);

/** This function defines lorentz transformationfrom
    first (kto=1) or second (kto=2) tau to laboratory frame.
    It's heavily used in Fortran code. */
extern "C" void tralo4_(float * kto, float p[4], float q[4], float * ams);

} // namespace Tauolapp
#endif
