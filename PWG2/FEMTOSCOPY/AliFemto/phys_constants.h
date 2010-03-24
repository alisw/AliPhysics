// malisa - there was somewhere a "phys_constants.h" in the STAR hierarchy
// the Strangeness HBT guys used very little of it, so I just put that stuff
// in here...

#ifndef _StHbt_phys_constants_h_
#define _StHbt_phys_constants_h_

#include "PhysicalConstants.h"  // from StarClassLibrary

// the strangeness guys used a different naming system (of course)
static const double M_LAMBDA        = kLambdaMassC2;
static const double M_KAON_0_SHORT  = kKaon0ShortMassC2;
static const double M_PROTON        = kProtonMassC2;
static const double M_PION_PLUS     = kPionPlusMassC2;
static const double M_PION_MINUS    = kPionMinusMassC2;
static const double M_XI_MINUS      = kXiMinusMassC2;
static const double M_OMEGA_MINUS   = 1672.45;
static const double M_KAON_MINUS    = 493.677;
static const double M_KAON_PLUS     = 493.677;

#endif
