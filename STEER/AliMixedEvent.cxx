/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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


//-------------------------------------------------------------------------
//                          Class AliMixedEvent
// VEvent which is the container of several VEvents 
// Use Case: Event Mixing     
// Origin: Andreas Morsch, CERN, Andreas.Morsch@cern.ch 
//-------------------------------------------------------------------------


#include "AliMixedEvent.h"
#include "AliExternalTrackParam.h"
#include "TVector3.h"
#include "AliVVertex.h"
#include <TMath.h>
#include <TMatrix.h>
#include <TMatrixD.h>

ClassImp(AliMixedEvent)


AliMixedEvent::AliMixedEvent() :
    AliVEvent(),
    fEventList(),
    fNEvents(0),       
    fNumberOfTracks(0),
    fNTracksCumul(0),
    fMeanVertex(0)
{
    // Default constructor
}

AliMixedEvent::AliMixedEvent(const AliMixedEvent& Evnt) :
    AliVEvent(Evnt),
    fEventList(),
    fNEvents(0),
    fNumberOfTracks(0),
    fNTracksCumul(0),
    fMeanVertex(0)
{ } // Copy constructor

AliMixedEvent& AliMixedEvent::operator=(const AliMixedEvent& vEvnt)
{ if (this!=&vEvnt) { 
    AliVEvent::operator=(vEvnt); 
  }
  
  return *this; 
}


void AliMixedEvent::AddEvent(AliVEvent* evt)
{
    // Add a new event to the list
    fEventList.AddLast(evt);
}


void AliMixedEvent::Init()
{
    // Initialize meta information
    fNEvents = fEventList.GetEntries();
    fNTracksCumul = new Int_t[fNEvents];
    fNumberOfTracks = 0;
    TIter next(&fEventList);
    AliVEvent* event;
    Int_t iev = 0;
    
    while((event = (AliVEvent*)next())) {
	fNTracksCumul[iev++] = fNumberOfTracks;
	fNumberOfTracks += (event->GetNumberOfTracks());
    }
}


AliVParticle* AliMixedEvent::GetTrack(Int_t i) const
{
    // Return track # i
    Int_t iEv  = TMath::BinarySearch(fNEvents, fNTracksCumul, i);
    while((iEv < (fNEvents - 1)) && (fNTracksCumul[iEv] == fNTracksCumul[iEv+1])) {iEv++;}

    Int_t irel = i - fNTracksCumul[iEv];
    AliVEvent* evt = (AliVEvent*) (fEventList.At(iEv));
    return (evt->GetTrack(irel));
}

const AliVVertex* AliMixedEvent::GetEventVertex(Int_t i) const
{
    // Return track # i
    Int_t iEv  = TMath::BinarySearch(fNEvents, fNTracksCumul, i);
    while((iEv < (fNEvents - 1)) && (fNTracksCumul[iEv] == fNTracksCumul[iEv+1])) {iEv++;}
    AliVEvent* evt = (AliVEvent*) (fEventList.At(iEv));
    return (evt->GetPrimaryVertex());
}

void AliMixedEvent::Reset()
{
    // Reset the event
    fEventList.Clear();
    fNEvents = 0;
    fNumberOfTracks = 0;
    if (fNTracksCumul) {
	delete[]  fNTracksCumul;
	fNTracksCumul = 0;
    }
}

Int_t AliMixedEvent::EventIndex(Int_t itrack)
{
  // Return the event index for track #itrack
  return  TMath::BinarySearch(fNEvents, fNTracksCumul, itrack);
}

Bool_t AliMixedEvent::ComputeVtx(TObjArray *vertices, Double_t *pos,Double_t *sig,Int_t *nContributors){
//
// Calculate the mean vertex psoitions from events in the buffer
 
    Int_t nentries = vertices->GetEntriesFast();
    Double_t sum[3]={0.,0.,0.};
    Double_t sumsigma[6]={0.,0.,0.,0.,0.,0.};

    
    for(Int_t ivtx = 0; ivtx < nentries; ivtx++){
	AliVVertex *vtx=(AliVVertex*)vertices->UncheckedAt(ivtx);
	Double_t covariance[6];
	vtx->GetCovarianceMatrix(covariance);
	Double_t vtxPos[3];
	vtx->GetXYZ(vtxPos);
	if(TMath::Abs(covariance[0])<1.e-13) {
	return kFALSE;
	}else{
	sum[0]+=vtxPos[0]*(1./covariance[0]);
	sumsigma[0]+=(1./covariance[0]);
	}
	if(TMath::Abs(covariance[2])<1.e-13) {
	return kFALSE;
	}else{
	sum[1]+=vtxPos[1]*(1./covariance[2]);
	sumsigma[2]+=(1./covariance[2]);
	}
	if(TMath::Abs(covariance[5])<1.e-13) {
	return kFALSE;
	}else{
	sum[2]+=vtxPos[2]*(1./covariance[5]);
	sumsigma[5]+=(1./covariance[5]);
	}
	if(TMath::Abs(covariance[1])<1.e-13) {
         sumsigma[1]+=0.;
	}else{
	sumsigma[1]+=(1./covariance[1]);
	}
	if(TMath::Abs(covariance[3])<1.e-13) {
	sumsigma[3]+=0.;
	}else{
	sumsigma[3]+=(1./covariance[3]);
	}
	if(TMath::Abs(covariance[4])<1.e-13) {
	sumsigma[4]+=0.;
	}else{
	sumsigma[4]+=(1./covariance[4]);
	}

     nContributors[0]=nContributors[0]+vtx->GetNContributors();
    }
    
    for(Int_t i=0;i<3;i++){
	if(TMath::Abs(sumsigma[i])<1.e-13) continue;
	pos[i]=sum[i]/sumsigma[i];
    }
    for(Int_t i2=0;i2<6;i2++){
	if(TMath::Abs(sumsigma[i2])<1.e-13) {sig[i2]=0.; continue;}
	sig[i2]=1./sumsigma[i2];
    }
    return kTRUE;
}


Double_t AliMixedEvent::GetMagneticField() const
{
    // Return magnetic field of the first event in the list
    if (fEventList.GetEntries() == 0) return -999.;
    
    AliVEvent* evt = (AliVEvent*) (fEventList.At(0));
    return evt->GetMagneticField();
}
