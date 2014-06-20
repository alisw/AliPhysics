#ifndef _f_Init_h_included_
#define _f_Init_h_included_

/**
 * This file contains an interface between the C++ code and TAUOLA
 * FORTRAN routines for tauola initalization.  
 * f_interface_tauolaInitialize() should be used
 * by C++ code. This call the initiphy_ or inimas_ routines defined in
 * tauola.f and tauola_extras.f
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 */

#include "f_Variables.h"
#include "f_Decay.h"

namespace Tauolapp
{

extern "C" {

  extern struct {
    int idff; //tau pdg id
  } idfc_;

  extern void inietc_(float jak1=0,float jak2=0,float itdkrc=1,float ifphot=1);
  extern void inimas_();
  extern void iniphx_(float *i);
  extern void initdk_();
  extern void iniphy_(float *i);
}

void f_interface_tauolaInitialize(int pdg_id, int firstDecayMode, 
                                  int secondDecayMode, bool rad,
                                  double rad_cut_off, double iniphy);

/** DEPRECATED: Use 'f_interface_tauolaInitialize' instead. */
void f_interface_tauolaInitialise(int pdg_id, int firstDecayMode, 
                                  int secondDecayMode, bool rad,
                                  double rad_cut_off, double iniphy);

double f_getTauMass();

} // namespace Tauolapp
#endif
