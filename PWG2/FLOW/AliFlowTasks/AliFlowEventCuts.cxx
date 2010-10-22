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

// AliFlowEventCuts:
// An event cut class for the flow framework
//
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#include <limits.h>
#include <float.h>
#include "TNamed.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliMCEvent.h"
#include "AliFlowEventCuts.h"
#include "AliFlowTrackCuts.h"

ClassImp(AliFlowEventCuts)

//-----------------------------------------------------------------------
AliFlowEventCuts::AliFlowEventCuts():
  TNamed(),
  fCutNumberOfTracks(kFALSE),
  fNumberOfTracksMax(INT_MAX),
  fNumberOfTracksMin(INT_MIN),
  fCutRefMult(kFALSE),
  fRefMultMethod(kTPConly),
  fRefMultMax(INT_MAX),
  fRefMultMin(INT_MIN),
  fRefMult(-1)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowEventCuts::AliFlowEventCuts(const char* name, const char* title):
  TNamed(name, title),
  fCutNumberOfTracks(kFALSE),
  fNumberOfTracksMax(INT_MAX),
  fNumberOfTracksMin(INT_MIN),
  fCutRefMult(kFALSE),
  fRefMultMethod(kTPConly),
  fRefMultMax(INT_MAX),
  fRefMultMin(INT_MIN),
  fRefMult(-1)
{
  //constructor 
}

////-----------------------------------------------------------------------
//AliFlowEventCuts::AliFlowEventCuts(const AliFlowEventCuts& someCuts):
//  TNamed(someCuts),
//  fCutNumberOfTracks(that.fCutNumberOfTracks),
//  fNumberOfTracksMax(that.fNumberOfTracksMax),
//  fNumberOfTracksMin(that.fNumberOfTracksMin),
//{
//  //copy constructor 
//}
//
////-----------------------------------------------------------------------
//AliFlowEventCuts& AliFlowEventCuts::operator=(const AliFlowEventCuts& someCuts)
//{
//  //assignment
//  fCutNumberOfTracks=that.fCutNumberOfTracks;
//  fNumberOfTracksMax=that.fNumberOfTracksMax;
//  fNumberOfTracksMin=that.fNumberOfTracksMin;
//
//  return *this;
//}

//----------------------------------------------------------------------- 
Bool_t AliFlowEventCuts::IsSelected(const TObject* obj)
{
  //check cuts
  const AliVEvent* vevent = dynamic_cast<const AliVEvent*>(obj);
  if (vevent) return PassesCuts(vevent);
  return kFALSE;  //when passed wrong type of object
}
//----------------------------------------------------------------------- 
Bool_t AliFlowEventCuts::PassesCuts(const AliVEvent *event)
{
  ///check if event passes cuts
  if(fCutNumberOfTracks) {if (event->GetNumberOfTracks() < fNumberOfTracksMin || event->GetNumberOfTracks() >= fNumberOfTracksMax ) return kFALSE;}
  if(fCutRefMult)
  {
    //reference multiplicity still to be defined
    fRefMult = RefMult(event);
    if (fRefMult < fRefMultMin || fRefMult >= fRefMultMax )
      return kFALSE;
  }
  return kTRUE;
}

//----------------------------------------------------------------------- 
AliFlowEventCuts* AliFlowEventCuts::StandardCuts()
{
  //make a set of standard event cuts, caller becomes owner
  AliFlowEventCuts* cuts = new AliFlowEventCuts();
  return cuts;
}

//----------------------------------------------------------------------- 
Int_t AliFlowEventCuts::RefMult(const AliVEvent* event, AliFlowTrackCuts* cuts)
{
  //calculate the reference multiplicity, if all fails return 0
  Int_t refmult=0;

  //in the case of an ESD
  const AliESDEvent* esd=dynamic_cast<const AliESDEvent*>(event);
  if (esd) 
  {
    AliMultiplicity* tracklets=NULL;
    switch (fRefMultMethod)
    {
      case kTPConly:
        if (!cuts)
        {
          //if not explicitly passed, make default cuts
          cuts = AliFlowTrackCuts::GetStandardTPCOnlyTrackCuts();
          cuts->SetEtaRange(-0.8,0.8);
          cuts->SetPtMin(0.15);
        }
        for (Int_t i=0; i<esd->GetNumberOfTracks();i++)
        {
          AliESDtrack* track = esd->GetTrack(i);
          if (cuts->IsSelected(track)) refmult++;
        }
        break;
      case kSPDtracklets:
        if (!cuts)
        {
          //if not explicitly passed, make default cuts
          cuts = new AliFlowTrackCuts();
          cuts->SetEtaRange(-0.8,0.8);
        }
        tracklets = const_cast<AliMultiplicity*>(esd->GetMultiplicity());
        for (Int_t i=0; i<tracklets->GetNumberOfTracklets(); i++)
        {
          if (cuts->IsSelected(tracklets,i)) refmult++;
        }
        break;
      default:
        return 0;
    }
    return refmult;
  }

  //in case of MC event
  AliMCEvent* mc=const_cast<AliMCEvent*>(dynamic_cast<const AliMCEvent*>(event));
  if (mc) return mc->GetNumberOfPrimaries();

  return event->GetNumberOfTracks(); //default, at least returns some number
}
