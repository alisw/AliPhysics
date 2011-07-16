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
#include "AliLog.h"
#include "AliVCaloCells.h"


ClassImp(AliMixedEvent)


AliMixedEvent::AliMixedEvent() :
  AliVEvent(),
  fEventList(),
  fNEvents(0),       
  fNumberOfTracks(0),
  fNumberOfCaloClusters(0), 
  fNumberOfPHOSCells(0), 
  fNumberOfEMCALCells(0),
  fNTracksCumul(0),
  fNCaloClustersCumul(0),
  fNPHOSCellsCumul(0), 
  fNEMCALCellsCumul(0), 
  fPHOSCells(NULL), 
  fEMCALCells(NULL), 
  fMeanVertex(0)
{
    // Default constructor
}

AliMixedEvent::AliMixedEvent(const AliMixedEvent& Evnt) :
  AliVEvent(Evnt),
  fEventList(),
  fNEvents(0),
  fNumberOfTracks(0),
  fNumberOfCaloClusters(0), 
  fNumberOfPHOSCells(0), 
  fNumberOfEMCALCells(0),
  fNTracksCumul(0),
  fNCaloClustersCumul(0),
  fNPHOSCellsCumul(0), 
  fNEMCALCellsCumul(0), 
  fPHOSCells(NULL), 
  fEMCALCells(NULL), 
  fMeanVertex(0)
{ } // Copy constructor

AliMixedEvent& AliMixedEvent::operator=(const AliMixedEvent& vEvnt)
{ 
// Assignment operator
    if (this!=&vEvnt) { 
    AliVEvent::operator=(vEvnt); 
}
  
  return *this; 
}

AliMixedEvent::~AliMixedEvent() 
{
    // dtor
  Reset();
  delete fPHOSCells ; 
  delete fEMCALCells ; 
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
  fNCaloClustersCumul = new Int_t[fNEvents];
  fNumberOfCaloClusters = 0;
  fNumberOfPHOSCells    = 0;  
  fNumberOfEMCALCells   = 0; 
  fNPHOSCellsCumul  = new Int_t[fNEvents];
  fNEMCALCellsCumul = new Int_t[fNEvents];

  TIter next(&fEventList);
  AliVEvent* event;
  Int_t iev = 0;
    
  while((event = (AliVEvent*)next())) {
    fNTracksCumul[iev] = fNumberOfTracks;
    fNumberOfTracks += (event->GetNumberOfTracks());
    fNCaloClustersCumul[iev] = fNumberOfCaloClusters;
    fNumberOfCaloClusters += event->GetNumberOfCaloClusters(); 
    fNPHOSCellsCumul[iev] = fNumberOfPHOSCells;
    if (event->GetPHOSCells()) 
      fNumberOfPHOSCells += event->GetPHOSCells()->GetNumberOfCells(); 
    fNEMCALCellsCumul[iev] = fNumberOfEMCALCells;
    if (event->GetEMCALCells()) 
      fNumberOfEMCALCells += event->GetEMCALCells()->GetNumberOfCells(); 
    iev++ ;  
  }

  next.Reset() ; 
  Short_t phosPos = 0, emcalPos = 0; 
  Int_t firstPHOSEvent  = kTRUE;
  Int_t firstEMCALEvent = kTRUE;
  
  while((event = (AliVEvent*)next())) {
    AliVCaloCells * phosCells = event->GetPHOSCells() ; 
    if (phosCells) {
      
      //Create the container
      if(firstPHOSEvent)
      {
        if(!fPHOSCells) fPHOSCells = phosCells->CopyCaloCells(kFALSE) ;// Just recover the first event type:  ESD/AOD
        else fPHOSCells->DeleteContainer(); //delete the previous container 
        //Now create a new container with the adequate size
        fPHOSCells->SetType(AliVCaloCells::kPHOSCell) ; 
        fPHOSCells->CreateContainer(fNumberOfPHOSCells) ;
        firstPHOSEvent=kFALSE;

      }//First event

      Int_t ncells = event->GetPHOSCells()->GetNumberOfCells() ;
      for (Int_t icell = 0; icell < ncells; icell++) {
          fPHOSCells->SetCell(phosPos++, phosCells->GetCellNumber(icell), phosCells->GetAmplitude(icell), phosCells->GetTime(icell)) ; 
      }
     
    }// phos cells
    
    AliVCaloCells * emcalCells = event->GetEMCALCells() ; 
    if (emcalCells) {
      
      //Create the container
      if(firstEMCALEvent)
      {
        if(!fEMCALCells)fEMCALCells = emcalCells->CopyCaloCells(kFALSE) ; // Just recover the first event type:  ESD/AOD
        else fEMCALCells->DeleteContainer();       // delete the previous container
        //Now create a new container with the adequate size
        fEMCALCells->SetType(AliVCaloCells::kEMCALCell) ; 
        fEMCALCells->CreateContainer(fNumberOfEMCALCells) ;
        firstEMCALEvent=kFALSE;
      }//First event
      
      Int_t ncells = emcalCells->GetNumberOfCells() ;
      for (Int_t icell = 0; icell < ncells; icell++) {
          fEMCALCells->SetCell(emcalPos++, emcalCells->GetCellNumber(icell), emcalCells->GetAmplitude(icell), emcalCells->GetTime(icell)) ; 
      }
    }//EMCAL cells
  }//while event
  
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

AliVCluster* AliMixedEvent::GetCaloCluster(Int_t i) const
{
    // Return calo cluster # i
  Int_t iEv  = TMath::BinarySearch(fNEvents, fNCaloClustersCumul, i);
  while((iEv < (fNEvents - 1)) && (fNCaloClustersCumul[iEv] == fNCaloClustersCumul[iEv+1])) {iEv++;}
  
  Int_t irel = i - fNCaloClustersCumul[iEv];
  AliVEvent* evt = (AliVEvent*) (fEventList.At(iEv));
  return (evt->GetCaloCluster(irel));
}

const AliVVertex* AliMixedEvent::GetEventVertex(Int_t i) const
{
    // Return vertex of track # i
    Int_t iEv  = TMath::BinarySearch(fNEvents, fNTracksCumul, i);
    while((iEv < (fNEvents - 1)) && (fNTracksCumul[iEv] == fNTracksCumul[iEv+1])) {iEv++;}
    AliVEvent* evt = (AliVEvent*) (fEventList.At(iEv));
    return (evt->GetPrimaryVertex());
}

const AliVVertex* AliMixedEvent::GetVertexOfEvent(Int_t i) const
{
    // Return vertex of event # i
  if (i > fNEvents)
    AliFatal(Form("%d events in buffer, event %d requested", fNEvents, i)) ;  
  AliVEvent* evt = (AliVEvent*) (fEventList.At(i));
  return (evt->GetPrimaryVertex());
}

void AliMixedEvent::Reset()
{
    // Reset the event
  fEventList.Clear();
  fNEvents = 0;
  fNumberOfTracks = 0;
  fNumberOfCaloClusters = 0;
  fNumberOfPHOSCells = 0;
  fNumberOfEMCALCells = 0;
  if (fNTracksCumul) {
    delete[]  fNTracksCumul;
    fNTracksCumul = 0;
  }
  if (fNCaloClustersCumul) {
    delete[]  fNCaloClustersCumul;
    fNCaloClustersCumul = 0;
  }
  if (fNPHOSCellsCumul) {
    delete[]  fNPHOSCellsCumul;
    fNPHOSCellsCumul = 0;
  }
  if (fNEMCALCellsCumul) {
    delete[]  fNEMCALCellsCumul;
    fNEMCALCellsCumul = 0;
  }
  
  if (fPHOSCells) {	 
    fPHOSCells->DeleteContainer();	 
  }	 
  if (fEMCALCells) {	 
    fEMCALCells->DeleteContainer();	 
  }
  
}

Int_t AliMixedEvent::EventIndex(Int_t itrack) const
{
  // Return the event index for track #itrack
  return  TMath::BinarySearch(fNEvents, fNTracksCumul, itrack);
}

Int_t AliMixedEvent::EventIndexForCaloCluster(Int_t icluster) const
{
    // Return the event index for track #itrack
  return  TMath::BinarySearch(fNEvents, fNCaloClustersCumul, icluster);
}

Int_t AliMixedEvent::EventIndexForPHOSCell(Int_t icell) const
{
    // Return the event index for track #itrack
  return  TMath::BinarySearch(fNEvents, fNPHOSCellsCumul, icell);
}

Int_t AliMixedEvent::EventIndexForEMCALCell(Int_t icell) const
{
    // Return the event index for track #itrack
  return  TMath::BinarySearch(fNEvents, fNEMCALCellsCumul, icell);
}

Bool_t AliMixedEvent::ComputeVtx(const TObjArray *vertices, Double_t *pos,Double_t *sig,Int_t *nContributors)  {
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
    for(Int_t i2=0;i2<3;i2++){
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
