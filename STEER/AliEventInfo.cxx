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

//////////////////////////////////////////////////////////////////////////////
//                          Class AliEventInfo                              //
//   Container class for all the information related to LHCstate, run and   //
//   event types, trigger mask and trigger clusters.                        //
//   It is used in order to provide the detector's AliRecoParam objects with//
//   the necessary information so that they can decide which instance of    //
//   AliDetectorRecoParam to use in reconstruction one particular event.    //
//                                                                          //
//   cvetan.cheshkov@cern.ch 12/06/2008                                     //
//////////////////////////////////////////////////////////////////////////////

#include "AliEventInfo.h"

ClassImp(AliEventInfo)

//______________________________________________________________________________
AliEventInfo::AliEventInfo():
  TObject(),
  fLHCState("UNKNOWN"),
  fRunType("UNKNOWN"),
  fActiveDetectors(""),
  fEventType(0),
  fTriggerClasses(""),
  fTriggerMask(0),
  fTriggerCluster("")
{
  // default constructor
  // ...
}

//______________________________________________________________________________
AliEventInfo::AliEventInfo(const char *lhcState, const char *runType, const char *activeDetectors):
  TObject(),
  fLHCState(lhcState),
  fRunType(runType),
  fActiveDetectors(activeDetectors),  
  fEventType(0),
  fTriggerClasses(""),
  fTriggerMask(0),
  fTriggerCluster("")
{
  // constructor
  // ...
}

//______________________________________________________________________________
AliEventInfo::AliEventInfo(const AliEventInfo &evInfo):
  TObject(evInfo),
  fLHCState(evInfo.fLHCState),
  fRunType(evInfo.fRunType),
  fActiveDetectors(evInfo.fActiveDetectors),
  fEventType(evInfo.fEventType),
  fTriggerClasses(evInfo.fTriggerClasses),
  fTriggerMask(evInfo.fTriggerMask),
  fTriggerCluster(evInfo.fTriggerCluster)
{
  // Copy constructor
  // ...
}

//_____________________________________________________________________________
AliEventInfo &AliEventInfo::operator =(const AliEventInfo& evInfo)
{
  // assignment operator
  // ...
  if(this==&evInfo) return *this;
  ((TObject *)this)->operator=(evInfo);

  fLHCState = evInfo.fLHCState;
  fRunType = evInfo.fRunType;
  fActiveDetectors = evInfo.fActiveDetectors;
  fEventType = evInfo.fEventType;
  fTriggerClasses = evInfo.fTriggerClasses;
  fTriggerMask = evInfo.fTriggerMask; 
  fTriggerCluster = evInfo.fTriggerCluster;

  return *this;
}

//______________________________________________________________________________
void AliEventInfo::Reset()
{
  // Reset the contents
  // ...
  fLHCState = "UNKNOWN";
  fRunType = "UNKNOWN";
  fActiveDetectors = "";
  fEventType = 0;
  fTriggerClasses = "";
  fTriggerMask = 0;
  fTriggerCluster = "";
}
