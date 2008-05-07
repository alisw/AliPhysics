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

  AliFlowEventSimple::AliFlowEventSimple(Int_t lenght):
    fTrackCollection(0),
    fTrack(0),
    fNumberOfTracks(0),
    fEventNSelTracksIntFlow(0)
{
  //constructor 
  fTrackCollection =  new TObjArray(lenght) ;


}

//-----------------------------------------------------------------------

AliFlowEventSimple::AliFlowEventSimple(const AliFlowEventSimple& event):
  TObject(),
  fTrackCollection(event.fTrackCollection),
  fTrack(event.fTrack),
  fNumberOfTracks(event.fNumberOfTracks),
  fEventNSelTracksIntFlow(event.fEventNSelTracksIntFlow)
{
  //copy constructor 
  //  *fTrack = *event.fTrack;
  //  *fTrackCollection =  *event.fTrackCollection ;

}

//-----------------------------------------------------------------------

AliFlowEventSimple& AliFlowEventSimple::operator=(const AliFlowEventSimple& event)
{
  *fTrack = *event.fTrack;
  *fTrackCollection =  *event.fTrackCollection ;
  fNumberOfTracks = event.fNumberOfTracks;
  fEventNSelTracksIntFlow = event.fEventNSelTracksIntFlow;

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
  fTrack = (AliFlowTrackSimple*)TrackCollection()->At(i) ;
  return fTrack;
}

//-----------------------------------------------------------------------   
 AliFlowVector AliFlowEventSimple::GetQ() 
{
  //calculate Q. 
  
  Double_t fQX = 0.;
  Double_t fQY = 0.;
  AliFlowVector fQ;
  fQ.Set(0.,0.);
  Double_t fOrder = 2.;
  Int_t fUsedTracks = 0;

  for (Int_t i=0;i<fNumberOfTracks;i++)                  
    {
      fTrack = (AliFlowTrackSimple*)TrackCollection()->At(i) ; 
      if (fTrack){
	if (fTrack->UseForIntegratedFlow()) {
	  Double_t fPhi = fTrack->Phi();
	  fQX += TMath::Cos(fOrder*fPhi);
	  fQY += TMath::Sin(fOrder*fPhi);
	  fUsedTracks++;
	}
      } //if particle
      else {cerr << "no particle!!!"<<endl;}
    }//loop over particles

  fQ.Set(fQX,fQY);
  fQ.SetMult(fUsedTracks);
   
  return fQ;
  
}
