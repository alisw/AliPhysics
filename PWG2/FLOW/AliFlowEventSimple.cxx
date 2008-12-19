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

#include "Riostream.h"
#include "TObjArray.h"
#include "TMath.h"
#include "AliFlowVector.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowEventSimple.h"

// AliFlowEventSimple:
// A simple event for flow analysis
//
//
// authors: N. van der Kolk (kolk@nikhef.nl), A. Bilandzic (anteb@nikhef.nl)


ClassImp(AliFlowEventSimple)

//-----------------------------------------------------------------------

  AliFlowEventSimple::AliFlowEventSimple(Int_t aLenght):
    fTrackCollection(NULL),
    fNumberOfTracks(0),
    fEventNSelTracksIntFlow(0)
{
  //constructor 
  fTrackCollection =  new TObjArray(aLenght) ;
}

//-----------------------------------------------------------------------

AliFlowEventSimple::AliFlowEventSimple(const AliFlowEventSimple& anEvent):
  TObject(),
  fTrackCollection(anEvent.fTrackCollection),
  fNumberOfTracks(anEvent.fNumberOfTracks),
  fEventNSelTracksIntFlow(anEvent.fEventNSelTracksIntFlow)
{
  //copy constructor 
}

//-----------------------------------------------------------------------

AliFlowEventSimple& AliFlowEventSimple::operator=(const AliFlowEventSimple& anEvent)
{
  *fTrackCollection =  *anEvent.fTrackCollection ;
  fNumberOfTracks = anEvent.fNumberOfTracks;
  fEventNSelTracksIntFlow = anEvent.fEventNSelTracksIntFlow;

  return *this;

}


//----------------------------------------------------------------------- 

AliFlowEventSimple::~AliFlowEventSimple()
{
  //destructor
  fTrackCollection->Delete() ; delete fTrackCollection ;
}

//----------------------------------------------------------------------- 

AliFlowTrackSimple* AliFlowEventSimple::GetTrack(Int_t i)
{
  //get track i from collection
  AliFlowTrackSimple* pTrack = (AliFlowTrackSimple*)TrackCollection()->At(i) ;
  return pTrack;
}

//-----------------------------------------------------------------------   
 AliFlowVector AliFlowEventSimple::GetQ(Int_t n) 
{
  //calculate Q. 
  
  Double_t dQX = 0.;
  Double_t dQY = 0.;
  AliFlowVector vQ;
  vQ.Set(0.,0.);
  
  Int_t iOrder = n;
  Int_t iUsedTracks = 0;

  for (Int_t i=0;i<fNumberOfTracks;i++)                  
    {
      AliFlowTrackSimple* pTrack = (AliFlowTrackSimple*)TrackCollection()->At(i) ; 
      if (pTrack){
	if (pTrack->UseForIntegratedFlow()) {
	  Double_t dPhi = pTrack->Phi();
	  dQX += TMath::Cos(iOrder*dPhi);
	  dQY += TMath::Sin(iOrder*dPhi);
	  iUsedTracks++;
	}
      } //if particle
      else {cerr << "no particle!!!"<<endl;}
    }//loop over particles

  vQ.Set(dQX,dQY);
  vQ.SetMult(iUsedTracks);
   
  return vQ;
  
}


