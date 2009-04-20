
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
/* $Id: $ */

//_________________________________________________________________________
// Class for reading data (AODs) in order to do prompt gamma
//  or other particle identification and correlations
// 
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
//#include "Riostream.h"

//---- ANALYSIS system ----
#include "AliCaloTrackAODReader.h" 
#include "AliAODEvent.h"
#include "AliAODCaloCluster.h"
#include "AliAODTrack.h"
#include "AliFidutialCut.h"

ClassImp(AliCaloTrackAODReader)

//____________________________________________________________________________
AliCaloTrackAODReader::AliCaloTrackAODReader() : 
  AliCaloTrackReader()
{
  //Default Ctor
  
  //Initialize parameters
  fDataType=kAOD;
  
}

//____________________________________________________________________________
AliCaloTrackAODReader::AliCaloTrackAODReader(const AliCaloTrackAODReader & g) :   
  AliCaloTrackReader(g)
{
  // cpy ctor
}

//_________________________________________________________________________
//AliCaloTrackAODReader & AliCaloTrackAODReader::operator = (const AliCaloTrackAODReader & source)
//{
//  // assignment operator
//
//  if(&source == this) return *this;
//
//  return *this;
//
//}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputCTS() {
  //Return array with CTS tracks

  Int_t nTracks   = fInputEvent->GetNumberOfTracks() ;
  Bool_t first = kTRUE;
  Double_t p[3];
  
  for (Int_t itrack =  0; itrack <  nTracks; itrack++) {////////////// track loop
    AliAODTrack * track = ((AliAODEvent*)fInputEvent)->GetTrack(itrack) ; // retrieve track from esd
    
    //     //We want tracks fitted in the detectors:
    //     ULong_t status=AliAODTrack::kTPCrefit;
    //     status|=AliAODTrack::kITSrefit;
    
    //We want tracks whose PID bit is set:
    //     ULong_t status =AliAODTrack::kITSpid;
    //     status|=AliAODTrack::kTPCpid;
    
    //   if ( (track->GetStatus() & status) == status) {//Check if the bits we want are set
    
    track->GetPxPyPz(p) ;
    TLorentzVector momentum(p[0],p[1],p[2],0);
    
    if(fCTSPtMin < momentum.Pt() && fFidutialCut->IsInFidutialCut(momentum,"CTS")){
      
      if(fDebug > 2 && momentum.Pt() > 0.1)printf("FillInputCTS():: Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						  momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
      
      if(first){
	new (fAODCTS) TRefArray(TProcessID::GetProcessWithUID(track)); 
	first=kFALSE;
      }     
      fAODCTS->Add(track);
      
    }//Pt and Fidutial cut passed. 
    //}// track status
  }// track loop
  if(fDebug > 1) printf("FillInputCTS():: aod entries %d\n", fAODCTS->GetEntriesFast());
}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputEMCAL() {
  //Return array with EMCAL clusters in aod format

  Bool_t first = kTRUE;
  
  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  GetVertex(v);
  
  //Loop to select clusters in fidutial cut and fill container with aodClusters
  //Int_t naod = 0;
  for (Int_t iclus =  0; iclus < ((AliAODEvent*)fInputEvent)->GetNCaloClusters(); iclus++) {
    AliAODCaloCluster * clus = 0;
    if ( (clus = ((AliAODEvent*)fInputEvent)->GetCaloCluster(iclus)) ) {
      if (clus->IsEMCALCluster()){
	TLorentzVector momentum ;
	clus->GetMomentum(momentum, v);      
	
	if(fEMCALPtMin < momentum.Pt() && fFidutialCut->IsInFidutialCut(momentum,"EMCAL")){
    
	  if(fDebug > 2 && momentum.E() > 0.1)printf("FillInputEMCAL():: Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						     momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	  
	  if(first){
	    new (fAODEMCAL) TRefArray(TProcessID::GetProcessWithUID(clus)); 
	    first=kFALSE;
	  }
	  fAODEMCAL->Add(clus);	
	}//Pt and Fidutial cut passed.
      }//EMCAL cluster
    }// cluster exists
  }//esd cluster loop
  
  if(fDebug > 1) printf("FillInputEMCAL():: aod entries %d\n", fAODEMCAL->GetEntriesFast());

}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputPHOS() {
  //Return array with PHOS clusters in aod format

   Bool_t first = kTRUE;
  
  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  GetVertex(v);

  //Loop to select clusters in fidutial cut and fill container with aodClusters
  //Int_t naod = 0;

  for (Int_t iclus =  0; iclus < ((AliAODEvent*)fInputEvent)->GetNCaloClusters(); iclus++) {
    AliAODCaloCluster * clus = 0;
    if ( (clus = ((AliAODEvent*)fInputEvent)->GetCaloCluster(iclus)) ) {
      if (clus->IsPHOSCluster()){
	TLorentzVector momentum ;
	clus->GetMomentum(momentum, v);      
	
	if(fPHOSPtMin < momentum.Pt() && fFidutialCut->IsInFidutialCut(momentum,"PHOS")){
	  
	  if(fDebug > 2 && momentum.E() > 0.1)printf("FillInputPHOS():: Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						     momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	  
	  if(first){
	    new (fAODPHOS) TRefArray(TProcessID::GetProcessWithUID(clus)); 
	    first=kFALSE;
	  }     
	  fAODPHOS->Add(clus);	
	}//Pt and Fidutial cut passed.
      }//PHOS cluster
    }//cluster exists
  }//esd cluster loop
  
  if(fDebug > 1) printf("FillInputPHOS():: aod entries %d\n", fAODPHOS->GetEntriesFast());


}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputEMCALCells() {
  //Return array with EMCAL cells in aod format

  fEMCALCells = (TNamed*) ((AliAODEvent*)fInputEvent)->GetEMCALCells(); 

}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputPHOSCells() {
  //Return array with PHOS cells in aod format

  fPHOSCells = (TNamed*) ((AliAODEvent*)fInputEvent)->GetPHOSCells(); 

}

//____________________________________________________________________________
void AliCaloTrackAODReader::GetVertex(Double_t  v[3]) const {
  //Return vertex position

  v[0] = ((AliAODEvent*)fInputEvent)->GetVertex(0)->GetX() ;//CHECK!!!
  v[1] = ((AliAODEvent*)fInputEvent)->GetVertex(0)->GetY() ;//CHECK!!!
  v[2] = ((AliAODEvent*)fInputEvent)->GetVertex(0)->GetZ() ;//CHECK!!!
}


//____________________________________________________________________________
void AliCaloTrackAODReader::SetInputOutputMCEvent(AliVEvent* input, AliAODEvent* aod, AliMCEvent* mc) {
  // Connect the data pointers
  // If input is AOD, do analysis with input, if not, do analysis with the output aod.
  if(!strcmp(input->GetName(),"AliESDEvent"))   {
    SetInputEvent(aod);
    SetOutputEvent(aod);
  }
  else if(!strcmp(input->GetName(),"AliAODEvent")){
    SetInputEvent(input);
    SetOutputEvent(aod);
  }
  else{ 
    printf("AliCaloTrackAODReader::SetInputOutputMCEvent() - ABORT::Unknown data format: %s",input->GetName());
    abort();
  }
  
  SetMC(mc);
  
}
