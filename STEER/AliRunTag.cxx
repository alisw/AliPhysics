/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------
//           Implementation of the RunTag class
//   This is the class to deal with the tags in the run level
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include <stdlib.h>
#include <Riostream.h>

#include "AliRunTag.h"
#include "AliLHCTag.h"
#include "AliDetectorTag.h"
#include "AliEventTag.h"

ClassImp(AliRunTag)

TClonesArray *AliRunTag::fgEvents = 0;
TClonesArray *AliRunTag::fgDetectors = 0;

//______________________________________________________________________________
AliRunTag::AliRunTag()
{
  if (!fgEvents) fgEvents = new TClonesArray("AliEventTag", 1000);
  fEventTag = fgEvents;
  fNumEvents = 0;
  
  if (!fgDetectors) fgDetectors = new TClonesArray("AliDetectorTag", 1000);
  fDetectorTag = fgDetectors;
  fNumDetectors = 0;
  
  fAliceMagneticField = 0.0;
  fAliceRunStartTime = 0;
  fAliceRunStopTime = 0;
  fAliceReconstructionVersion = 0;
  fAliceRunQuality = 0;
  fAliceBeamEnergy = 0.0;
  fAliceCalibrationVersion = 0;	
  fAliceDataType = 0;
}

//______________________________________________________________________________
AliRunTag::~AliRunTag()
{
}

//______________________________________________________________________________
void AliRunTag::SetLHCTag(Float_t lumin, char *type)
{
  fLHCTag.SetLHCTag(lumin,type);
}

//______________________________________________________________________________
void AliRunTag::SetDetectorTag(AliDetectorTag *DetTag)
{
  TClonesArray &detectors = *fDetectorTag;
  new(detectors[fNumDetectors++]) AliDetectorTag(DetTag);
}

//______________________________________________________________________________
void AliRunTag::AddEventTag(AliEventTag *EvTag)
{
  TClonesArray &events = *fEventTag;
  new(events[fNumEvents++]) AliEventTag(EvTag);
}

//______________________________________________________________________________
void AliRunTag::Clear(const char *)
{
  fNumEvents = 0;
  fNumDetectors = 0;
}
