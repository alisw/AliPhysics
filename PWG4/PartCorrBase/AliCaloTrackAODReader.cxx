
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
//  or other particle identification and correlations.
// Mixing analysis can be done, input AOD with events
// is opened in the AliCaloTrackReader::Init()
// 
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
//#include "Riostream.h"

//---- ANALYSIS system ----
#include "AliCaloTrackAODReader.h" 
#include "AliAODCaloCluster.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliFiducialCut.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

ClassImp(AliCaloTrackAODReader)

//____________________________________________________________________________
AliCaloTrackAODReader::AliCaloTrackAODReader() : 
  AliCaloTrackReader(), fWriteOutputStdAOD(kFALSE)
{
  //Default Ctor
  
  //Initialize parameters
  fDataType=kAOD;
  fReadStack          = kTRUE;
  fReadAODMCParticles = kFALSE;
  //We want tracks fitted in the detectors:
  fTrackStatus=AliESDtrack::kTPCrefit;
  fTrackStatus|=AliESDtrack::kITSrefit;

}

//____________________________________________________________________________
AliCaloTrackAODReader::AliCaloTrackAODReader(const AliCaloTrackAODReader & aodr) :   
  AliCaloTrackReader(aodr), fWriteOutputStdAOD(aodr.fWriteOutputStdAOD)
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
  //Return array with Central Tracking System (CTS) tracks
  if(fDebug > 2 ) printf("AliCaloTrackAODReader::FillInputCTS()\n");
  Int_t nTracks   = fInputEvent->GetNumberOfTracks() ;
  Int_t naod = 0;
  Double_t p[3];
  for (Int_t itrack =  0; itrack <  nTracks; itrack++) {////////////// track loop
    AliAODTrack * track = ((AliAODEvent*)fInputEvent)->GetTrack(itrack) ; // retrieve track from esd
    	  
    //Select tracks under certain conditions, TPCrefit, ITSrefit ... check the set bits
	if (fTrackStatus && !((track->GetStatus() & fTrackStatus) == fTrackStatus)) continue ;
    
    track->GetPxPyPz(p) ;
    TLorentzVector momentum(p[0],p[1],p[2],0);
    
    if(fCTSPtMin < momentum.Pt() && fFiducialCut->IsInFiducialCut(momentum,"CTS")){
      
      if(fDebug > 2 && momentum.Pt() > 0.1) printf("AliCaloTrackAODReader::FillInputCTS() - Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						  momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	  
	  if(fWriteOutputStdAOD){
		   AliAODTrack* newtrack = new((*(fOutputEvent->GetTracks()))[naod++]) AliAODTrack(*track);
		   fAODCTS->Add(newtrack); //Use AOD stored in output for references.
	  }
	  else fAODCTS->Add(track);
    }//Pt and Fiducial cut passed. 
  }// track loop
	
  fAODCTSNormalInputEntries = fAODCTS->GetEntriesFast();
  if(fDebug > 1) printf("AliCaloTrackAODReader::FillInputCTS()   - aod entries %d\n", fAODCTSNormalInputEntries);

  //If second input event available, add the clusters.
  if(fSecondInputAODTree && fSecondInputAODEvent){
	  nTracks   = fSecondInputAODEvent->GetNumberOfTracks() ;
	  if(fDebug > 1) printf("AliCaloTrackAODReader::FillInputCTS()   - Add second input tracks, entries %d\n", nTracks) ;
	  for (Int_t itrack =  0; itrack <  nTracks; itrack++) {////////////// track loop
		  AliAODTrack * track = ((AliAODEvent*)fSecondInputAODEvent)->GetTrack(itrack) ; // retrieve track from esd
		  
		  //Select tracks under certain conditions, TPCrefit, ITSrefit ... check the set bits
		  if (fTrackStatus && !((track->GetStatus() & fTrackStatus) == fTrackStatus)) continue ;
		  
		  track->GetPxPyPz(p) ;
		  TLorentzVector momentum(p[0],p[1],p[2],0);
		  
		  if(fCTSPtMin < momentum.Pt() && fFiducialCut->IsInFiducialCut(momentum,"CTS")){
			  
			  if(fDebug > 2 && momentum.Pt() > 0.1) printf("AliCaloTrackAODReader::FillInputCTS() - Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
														   momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
			  
			  if(fWriteOutputStdAOD){
				  AliAODTrack* newtrack = new((*(fOutputEvent->GetTracks()))[naod++]) AliAODTrack(*track);
				  fAODCTS->Add(newtrack); //Use AOD stored in output for references.
			  }
			  else fAODCTS->Add(track);
			  
		  }//Pt and Fiducial cut passed. 
	  }// track loop
	  
	  if(fDebug > 1) printf("AliCaloTrackAODReader::FillInputCTS()   - aod normal entries %d, after second input %d\n", fAODCTSNormalInputEntries, fAODCTS->GetEntriesFast());
  }	//second input loop
	
}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputEMCAL() {
  //Return array with EMCAL clusters in aod format
  if(fDebug > 2 ) printf("AliCaloTrackAODReader::FillInputEMCAL()\n");
	
  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  GetVertex(v);

  Int_t naod =  (fOutputEvent->GetCaloClusters())->GetEntriesFast();
  //Loop to select clusters in fiducial cut and fill container with aodClusters
  Int_t nclusters = ((AliAODEvent*)fInputEvent)->GetNCaloClusters();
  for (Int_t iclus =  0; iclus <  nclusters; iclus++) {
    AliAODCaloCluster * clus = 0;
    if ( (clus = ((AliAODEvent*)fInputEvent)->GetCaloCluster(iclus)) ) {
      if (clus->IsEMCALCluster()){
		  
	//Check if the cluster contains any bad channel and if close to calorimeter borders
	if(GetCaloUtils()->ClusterContainsBadChannel("EMCAL",clus->GetCellsAbsId(), clus->GetNCells())) continue;
	if(!GetCaloUtils()->CheckCellFiducialRegion(clus,((AliAODEvent*)fInputEvent)->GetEMCALCells())) continue;

	TLorentzVector momentum ;
	clus->GetMomentum(momentum, v);      
	
	if(fEMCALPtMin < momentum.Pt() && fFiducialCut->IsInFiducialCut(momentum,"EMCAL")){
    
	  if(fDebug > 2 && momentum.E() > 0.1) printf("AliCaloTrackAODReader::FillInputEMCAL() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						     momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	  
	  if(fWriteOutputStdAOD){
		AliAODCaloCluster * newclus = new((*(fOutputEvent->GetCaloClusters()))[naod++])AliAODCaloCluster(*clus);
		fAODEMCAL->Add(newclus);	
	  }
	  else fAODEMCAL->Add(clus);	
	}//Pt and Fiducial cut passed.
      }//EMCAL cluster
    }// cluster exists
  }// cluster loop
	
  fAODEMCALNormalInputEntries = fAODEMCAL->GetEntriesFast();
  if(fDebug > 1) printf("AliCaloTrackAODReader::FillInputEMCAL() - aod entries %d\n", fAODEMCALNormalInputEntries);

  //If second input event available, add the clusters.
  if(fSecondInputAODTree && fSecondInputAODEvent){
	  GetSecondInputAODVertex(v);
	  nclusters = ((AliAODEvent*)fSecondInputAODEvent)->GetNCaloClusters();
	  if(fDebug > 1) printf("AliCaloTrackAODReader::FillInputEMCAL() - Add second input clusters, entries %d\n", nclusters) ;
		for (Int_t iclus =  0; iclus < nclusters; iclus++) {
			AliAODCaloCluster * clus = 0;
			if ( (clus = ((AliAODEvent*)fSecondInputAODEvent)->GetCaloCluster(iclus)) ) {
				if (clus->IsEMCALCluster()){
					TLorentzVector momentum ;
					clus->GetMomentum(momentum, v);      
					
					if(fEMCALPtMin < momentum.Pt() && fFiducialCut->IsInFiducialCut(momentum,"EMCAL")){
						
						if(fDebug > 2 && momentum.E() > 0.1) printf("AliCaloTrackAODReader::FillInputEMCAL() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
																	momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
						if(fWriteOutputStdAOD){
							AliAODCaloCluster * newclus = new((*(fOutputEvent->GetCaloClusters()))[naod++])AliAODCaloCluster(*clus);
						    fAODEMCAL->Add(newclus);	
						}
						else fAODEMCAL->Add(clus);	
					}//Pt and Fiducial cut passed.
				}//EMCAL cluster
			}// cluster exists
		}// cluster loop
		
	  if(fDebug > 1) printf("AliCaloTrackAODReader::FillInputEMCAL() - aod normal entries %d, after second input %d\n", fAODEMCALNormalInputEntries, fAODEMCAL->GetEntriesFast());

	} //second input loop
}

//____________________________________________________________________________
void AliCaloTrackAODReader::FillInputPHOS() {
  //Return array with PHOS clusters in aod format
  if(fDebug > 2 ) printf("AliCaloTrackAODReader::FillInputPHOS()\n");
	
  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  GetVertex(v);
	
  Int_t naod =  (fOutputEvent->GetCaloClusters())->GetEntriesFast();
  //Loop to select clusters in fiducial cut and fill container with aodClusters
  Int_t nclusters = ((AliAODEvent*)fInputEvent)->GetNCaloClusters();
  for (Int_t iclus =  0; iclus < nclusters; iclus++) {
    AliAODCaloCluster * clus = 0;
    if ( (clus = ((AliAODEvent*)fInputEvent)->GetCaloCluster(iclus)) ) {
      if (clus->IsPHOSCluster()){
		  
	//Check if the cluster contains any bad channel and if close to calorimeter borders
	if( GetCaloUtils()->ClusterContainsBadChannel("PHOS",clus->GetCellsAbsId(), clus->GetNCells())) continue;
	if(!GetCaloUtils()->CheckCellFiducialRegion(clus, ((AliAODEvent*)fInputEvent)->GetPHOSCells())) continue;

	TLorentzVector momentum ;
	clus->GetMomentum(momentum, v);      
	
	if(fPHOSPtMin < momentum.Pt() && fFiducialCut->IsInFiducialCut(momentum,"PHOS")){
	  
	  if(fDebug > 2 && momentum.E() > 0.1) printf("AliCaloTrackAODReader::FillInputPHOS() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						     momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());

		if(fWriteOutputStdAOD){
			AliAODCaloCluster * newclus = new((*(fOutputEvent->GetCaloClusters()))[naod++])AliAODCaloCluster(*clus);
			fAODPHOS->Add(newclus);	
		}
		else fAODPHOS->Add(clus);	
	}//Pt and Fiducial cut passed.
      }//PHOS cluster
    }//cluster exists
  }//esd cluster loop
	
  fAODPHOSNormalInputEntries = fAODPHOS->GetEntriesFast() ;
  if(fDebug > 1) printf("AliCaloTrackAODReader::FillInputPHOS()  - aod entries %d\n", fAODPHOSNormalInputEntries);

  //If second input event available, add the clusters.
  if(fSecondInputAODTree && fSecondInputAODEvent){  
	  GetSecondInputAODVertex(v);
	  nclusters = ((AliAODEvent*)fSecondInputAODEvent)->GetNCaloClusters();
	  if(fDebug > 1) printf("AliCaloTrackAODReader::FillInputPHOS()  - Add second input clusters, entries %d\n", nclusters);
		for (Int_t iclus =  0; iclus < nclusters; iclus++) {
			AliAODCaloCluster * clus = 0;
			if ( (clus = ((AliAODEvent*)fSecondInputAODEvent)->GetCaloCluster(iclus)) ) {
				if (clus->IsPHOSCluster()){
					TLorentzVector momentum ;
					clus->GetMomentum(momentum, v);      
					
					if(fPHOSPtMin < momentum.Pt() && fFiducialCut->IsInFiducialCut(momentum,"PHOS")){
						
						if(fDebug > 2 && momentum.E() > 0.1) printf("AliCaloTrackAODReader::FillInputPHOS() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
																	momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
						if(fWriteOutputStdAOD){
							AliAODCaloCluster * newclus = new((*(fOutputEvent->GetCaloClusters()))[naod++])AliAODCaloCluster(*clus);
						    fAODPHOS->Add(newclus);	
						}
						else fAODPHOS->Add(clus);	
					}//Pt and Fiducial cut passed.
				}//PHOS cluster
			}// cluster exists
		}// cluster loop
		if(fDebug > 1) printf("AliCaloTrackAODReader::FillInputPHOS()  - aod normal entries %d, after second input %d\n", fAODPHOSNormalInputEntries, fAODPHOS->GetEntriesFast());
  }	//second input loop

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

 ((AliAODEvent*)fInputEvent)->GetPrimaryVertex()->GetXYZ(v);
	
}

//____________________________________________________________________________
void AliCaloTrackAODReader::GetSecondInputAODVertex(Double_t  v[3]) const {
	//Return vertex position of second AOD input
	
	fSecondInputAODEvent->GetPrimaryVertex()->GetXYZ(v);

}

//____________________________________________________________________________
Double_t AliCaloTrackAODReader::GetBField() const {
  //Return magnetic field

  Double_t bfield =  ((AliAODEvent*)fInputEvent)->GetMagneticField();

  return bfield;

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
	  
	  AliAODInputHandler* aodIH = static_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
	  //printf("AODInputHandler %p, MergeEvents %d \n",aodIH, aodIH->GetMergeEvents());
	  if (aodIH && aodIH->GetMergeEvents()) {
		  //Merged events, use output AOD.
		  SetInputEvent(aod);
		  SetOutputEvent(aod);
	  }
	  else{
		  SetInputEvent(input);
		  SetOutputEvent(aod);
	  }
  }
  else{ 
    printf("AliCaloTrackAODReader::SetInputOutputMCEvent() - STOP : Wrong data format: %s\n",input->GetName());
    abort();
  }
  
  SetMC(mc);
  
}

//________________________________________________________________
void AliCaloTrackAODReader::Print(const Option_t * opt) const
{
	
	//Print some relevant parameters set for the analysis
	AliCaloTrackReader::Print(opt);
	
	printf("Write std AOD       =     %d\n", fWriteOutputStdAOD) ;
	
	printf("    \n") ;
} 

