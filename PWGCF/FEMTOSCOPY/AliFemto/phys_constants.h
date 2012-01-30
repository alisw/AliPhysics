// malisa - there was somewhere a "phys_constants.h" in the STAR hierarchy
// the Strangeness HBT guys used very little of it, so I just put that stuff
// in here...

#ifndef _StHbt_phys_constants_h_
#define _StHbt_phys_constants_h_

#include "PhysicalConstants.h"  // from StarClassLibrary

// the strangeness guys used a different naming system (of course)
static const double kMLAMBDA        = kLambdaMassC2;
static const double kMKAON0SHORT  = kKaon0ShortMassC2;
static const double kMPROTON        = kProtonMassC2;
static const double kMPIONPLUS     = kPionPlusMassC2;
static const double kMPIONMINUS    = kPionMinusMassC2;
static const double kMXIMINUS      = kXiMinusMassC2;
static const double kMOMEGAMINUS   = 1672.45;
static const double kMKAONMINUS    = 493.677;
static const double kMKAONPLUS     = 493.677;

#endif
