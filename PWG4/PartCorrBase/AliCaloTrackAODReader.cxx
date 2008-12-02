
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
#include "AliAODVertex.h"
#include "AliAODCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliLog.h"

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
  fAODCTS = new TClonesArray("AliAODTrack",0);

  Int_t nTracks   = fAOD->GetNumberOfTracks() ;
  Int_t naod = 0;
  Double_t p[3];
  
  for (Int_t itrack =  0; itrack <  nTracks; itrack++) {////////////// track loop
    AliAODTrack * track = fAOD->GetTrack(itrack) ; // retrieve track from esd
    
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
      
      new((*fAODCTS)[naod++])  AliAODTrack(*track);
      
    }//Pt and Fidutial cut passed. 
    //}// track status
  }// track loop
  if(fDebug > 1) printf("FillInputCTS():: aod entries %d\n", fAODCTS->GetEntriesFast());
}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputEMCAL() {
  //Return array with EMCAL clusters in aod format

   fAODEMCAL = new TClonesArray("AliAODCaloCluster",0);
   TRefArray * caloClusters = new TRefArray();
   fAOD->GetEMCALClusters(caloClusters);

  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  GetVertex(v);

  //Loop to select clusters in fidutial cut and fill container with aodClusters
  Int_t naod = 0;
  for (Int_t iclus =  0; iclus <  caloClusters->GetEntriesFast(); iclus++) {
    AliAODCaloCluster * clus = (AliAODCaloCluster *) caloClusters->At(iclus) ;
    TLorentzVector momentum ;
    clus->GetMomentum(momentum, v);      
    
    if(fEMCALPtMin < momentum.Pt() && fFidutialCut->IsInFidutialCut(momentum,"EMCAL")){
    
      if(fDebug > 2 && momentum.E() > 0.1)printf("FillInputEMCAL():: Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						 momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());

      new((*fAODEMCAL)[naod++]) AliAODCaloCluster(*clus);

    }//Pt and Fidutial cut passed.
  }//esd cluster loop
  
  if(fDebug > 1) printf("FillInputEMCAL():: aod entries %d\n", fAODEMCAL->GetEntriesFast());

}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputPHOS() {
  //Return array with PHOS clusters in aod format

  fAODPHOS = new TClonesArray("AliAODCaloCluster",0);
  TRefArray * caloClusters = new TRefArray();
  fAOD->GetPHOSClusters(caloClusters);

  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  GetVertex(v);

  //Loop to select clusters in fidutial cut and fill container with aodClusters
  Int_t naod = 0;

  for (Int_t iclus =  0; iclus <  caloClusters->GetEntriesFast(); iclus++) {
    AliAODCaloCluster * clus = (AliAODCaloCluster *) caloClusters->At(iclus) ;
    TLorentzVector momentum ;
    clus->GetMomentum(momentum, v);      
    
    if(fPHOSPtMin < momentum.Pt() && fFidutialCut->IsInFidutialCut(momentum,"PHOS")){
      
      if(fDebug > 2 && momentum.E() > 0.1)printf("FillInputPHOS():: Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						 momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());

      new((*fAODPHOS)[naod++])  AliAODCaloCluster(*clus);
      
    }//Pt and Fidutial cut passed.
  }//esd cluster loop
  
  if(fDebug > 1) printf("FillInputPHOS():: aod entries %d\n", fAODPHOS->GetEntriesFast());


}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputEMCALCells() {
  //Return array with EMCAL cells in aod format

  fEMCALCells = (TNamed*) fAOD->GetEMCALCells(); 

}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputPHOSCells() {
  //Return array with PHOS cells in aod format

  fPHOSCells = (TNamed*) fAOD->GetPHOSCells(); 

}

//____________________________________________________________________________
void AliCaloTrackAODReader::GetVertex(Double_t  v[3]) const {
  //Return vertex position

  v[0]=fAOD->GetVertex(0)->GetX() ;//CHECK!!!
  v[1]=fAOD->GetVertex(0)->GetY() ;//CHECK!!!
  v[2]=fAOD->GetVertex(0)->GetZ() ;//CHECK!!!
}


//____________________________________________________________________________
void AliCaloTrackAODReader::SetInputEvent(TObject* input, TObject* aod, TObject* mc) {
  // Connect the data pointers

  //If input is AOD, do analysis with input, if not, do analysis with the output aod.
  if(!strcmp(input->GetName(),"AliESDEvent"))   SetAOD((AliAODEvent*) aod);
  else if(!strcmp(input->GetName(),"AliAODEvent")) SetAOD((AliAODEvent*) input);
  else AliFatal(Form("Unknown data format: %s",input->GetName()));
  
  SetMC((AliMCEvent*) mc);

}
