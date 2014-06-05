#ifndef _f_Decay_h_included_
#define _f_Decay_h_included_

/**
 * This file contains an interface between the C++ code and TAUOLA
 * FORTRAN routines for decaying taus. TauolaDecay() should be used
 * by C++ code. This call the dexay_ or dekay_ routines defined in
 * tauola.f
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 */

#include <vector>
#include <iostream>

namespace Tauolapp
{

extern "C" {

  extern struct{//positions of taus in the LUND common block
    int npa;
    int npb;
  } taupos_;

  //extern void dexay_(int *state, double pol[4]);
  extern void dekay_(int *state, double pol[4]);

  extern void taupi0_(double pp[4],int *k);
  extern void tauk0s_(double pp[4],int *k);
  extern void taueta_(double pp[4],int *k);
}

/** Invokes DEKAY with "1" or "2" to get the polarization information. */
void TauolaDecay(int sign_type, double *polx, double *poly,
                 double *polz, double *poln);

/** Invokes DEKAY with "11" or "12" to produce the decay. */
void TauolaWriteDecayToEventRecord(int sign_type);

} // namespace Tauolapp
#endif
