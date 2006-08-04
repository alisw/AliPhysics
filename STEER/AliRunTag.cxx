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

#include "AliRunTag.h"
#include "AliDetectorTag.h"
#include "AliEventTag.h"

ClassImp(AliRunTag)

//___________________________________________________________________________
  AliRunTag::AliRunTag() :
    TObject(),
    fAliceRunId(-1),
    fAliceMagneticField(0.0),
    fAliceRunStartTime(0),
    fAliceRunStopTime(0),
    fAlirootVersion(0),
    fRootVersion(0),
    fGeant3Version(0),
    fAliceRunQuality(0),
    fAliceBeamEnergy(0.0),
    fAliceBeamType(0),
    fAliceCalibrationVersion(0),
    fAliceDataType(0),
    fNumEvents(0),
    fNumDetectors(0),
    fEventTag("AliEventTag", 1000),
    fDetectorTag("AliDetectorTag", 1000),
    fLHCTag()
{
  //Default constructor
}

//___________________________________________________________________________
AliRunTag::~AliRunTag() {
  //Destructor
  fEventTag.Delete();
  fDetectorTag.Delete();
}

//___________________________________________________________________________
void AliRunTag::SetLHCTag(Float_t lumin, char *type) {
  //Setter for the LHC tags
  fLHCTag.SetLHCTag(lumin,type);
}

//___________________________________________________________________________
void AliRunTag::SetDetectorTag(const AliDetectorTag &DetTag) {
  //Setter for the detector tags
  new(fDetectorTag[fNumDetectors++]) AliDetectorTag(DetTag);
}

//___________________________________________________________________________
void AliRunTag::AddEventTag(const AliEventTag & EvTag) {
  //Adds an entry to the event tag TClonesArray
  new(fEventTag[fNumEvents++]) AliEventTag(EvTag);
}

//___________________________________________________________________________
void AliRunTag::Clear(const char *) {
  //Resets the number of events and detectors
  fEventTag.Clear();
  fNumEvents = 0;
  fDetectorTag.Clear();
  fNumDetectors = 0;
}
