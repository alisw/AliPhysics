#ifndef AliHLTMUONPRINT_H
#define AliHLTMUONPRINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
// Print routines to display internal dHLT data on the console.
//
////////////////////////////////////////////////////////////////////////////////

#include <ostream>

#include "AliHLTMUONCoreEventID.h"
#include "AliHLTMUONCorePoint.h"
#include "AliHLTMUONCoreTriggerRecord.h"
#include "AliHLTMUONCoreRegionOfInterest.h"

std::ostream& operator << (std::ostream& os, const AliHLTMUONCoreEventID& id);
std::ostream& operator << (std::ostream& os, const AliHLTMUONCorePoint& p);
std::ostream& operator << (std::ostream& os, const AliHLTMUONCoreParticleSign s);
std::ostream& operator << (std::ostream& os, const AliHLTMUONCoreTriggerRecord& rec);
std::ostream& operator << (std::ostream& os, const AliHLTMUONCoreChamberID chamber);

#endif // AliHLTMUONPRINT_H
