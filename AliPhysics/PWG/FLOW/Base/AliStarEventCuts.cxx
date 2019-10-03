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

// AliStarEventCuts:
// An event cut class for the AliStarEvent
//
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#include <limits.h>
#include <float.h>
#include "TNamed.h"
#include "AliStarEvent.h"
#include "AliStarEventCuts.h"

ClassImp(AliStarEventCuts)

//-----------------------------------------------------------------------
AliStarEventCuts::AliStarEventCuts():
  TNamed(),
  fCutRunID(kFALSE),
  fRunIDMax(INT_MAX),
  fRunIDMin(INT_MIN),
  fCutEventNumber(kFALSE),
  fEventNumberMax(INT_MAX),
  fEventNumberMin(INT_MIN),
  fCutVtxX(kFALSE),
  fVtxXMax(FLT_MAX),
  fVtxXMin(-FLT_MAX),
  fCutVtxY(kFALSE),
  fVtxYMax(FLT_MAX),
  fVtxYMin(-FLT_MAX),
  fCutVtxZ(kFALSE),
  fVtxZMax(FLT_MAX),
  fVtxZMin(-FLT_MAX),
  fCutBField(kFALSE),
  fBFieldMax(FLT_MAX),
  fBFieldMin(-FLT_MAX),
  fCutRefMult(kFALSE),
  fRefMultMax(INT_MAX),
  fRefMultMin(INT_MIN),
  fCutCentralityID(kFALSE),
  fCentralityIDMax(INT_MAX),
  fCentralityIDMin(INT_MIN),
  fCutNumberOfPrimaryTracks(kFALSE),
  fNumberOfPrimaryTracksMax(INT_MAX),
  fNumberOfPrimaryTracksMin(INT_MIN),
  fCutNumberOfTracks(kFALSE),
  fNumberOfTracksMax(INT_MAX),
  fNumberOfTracksMin(INT_MIN)
{
  //constructor 
}

////-----------------------------------------------------------------------
//AliStarEventCuts::AliStarEventCuts(const AliStarEventCuts& someCuts):
//  TNamed(),
//  fCutID(that.fCutID),
//  fIDMax(that.fIDMax),
//  fIDMin(that.fIDMin),
//{
//  //copy constructor 
//}
//
////-----------------------------------------------------------------------
//AliStarEventCuts& AliStarEventCuts::operator=(const AliStarEventCuts& someCuts)
//{
//  //assignment
//  fCutID=that.fCutID;
//  fIDMax=that.fIDMax;
//  fIDMin=that.fIDMin;
//
//  return *this;
//}

//----------------------------------------------------------------------- 
Bool_t AliStarEventCuts::PassesCuts(const AliStarEvent *event) const
{
  ///check if event passes cuts
  if(fCutRunID) {if (event->GetRunID() < fRunIDMin || event->GetRunID() > fRunIDMax ) return kFALSE;} //integer value: non inclusive bounds!
  if(fCutEventNumber) {if (event->GetEventNumber() < fEventNumberMin || event->GetEventNumber() > fEventNumberMax ) return kFALSE;}
  if(fCutVtxX) {if (event->GetVtxX() < fVtxXMin || event->GetVtxX() >= fVtxXMax ) return kFALSE;}
  if(fCutVtxY) {if (event->GetVtxY() < fVtxYMin || event->GetVtxY() >= fVtxYMax ) return kFALSE;}
  if(fCutVtxZ) {if (event->GetVtxZ() < fVtxZMin || event->GetVtxZ() >= fVtxZMax ) return kFALSE;}
  if(fCutBField) {if (event->GetBField() < fBFieldMin || event->GetBField() >= fBFieldMax ) return kFALSE;}
  if(fCutRefMult) {if (event->GetRefMult() < fRefMultMin || event->GetRefMult() > fRefMultMax ) return kFALSE;}
  if(fCutCentralityID) {if (event->GetCentralityID() < fCentralityIDMin || event->GetCentralityID() > fCentralityIDMax ) return kFALSE;}
  if(fCutNumberOfPrimaryTracks) {if (event->GetNumberOfPrimaryTracks() < fNumberOfPrimaryTracksMin || event->GetNumberOfPrimaryTracks() > fNumberOfPrimaryTracksMax ) return kFALSE;}
  if(fCutNumberOfTracks) {if (event->GetNumberOfTracks() < fNumberOfTracksMin || event->GetNumberOfTracks() > fNumberOfTracksMax ) return kFALSE;}
  return kTRUE;
}

//----------------------------------------------------------------------- 
AliStarEventCuts* AliStarEventCuts::StandardCuts()
{
  //make a set of standard event cuts, caller becomes owner
  AliStarEventCuts* cuts = new AliStarEventCuts();
  cuts->SetVtxXMin(-1.0);
  cuts->SetVtxXMax(1.0);
  cuts->SetVtxYMin(-1.0);
  cuts->SetVtxYMax(1.0);
  cuts->SetVtxZMin(-30.0);
  cuts->SetVtxZMax(30.0);
  cuts->SetRefMultMin(10);
  cuts->SetRefMultMax(1000);
  cuts->SetNumberOfPrimaryTracksMin(0);
  cuts->SetNumberOfPrimaryTracksMax(3000);
  return cuts;
}
