
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
/* $Id:  $ */

//_________________________________________________________________________
// Class for reading data (ESDs) in order to do prompt gamma 
// or other particle identification and correlations
//
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
//#include "Riostream.h"

//---- ANALYSIS system ----
#include "AliCaloTrackESDReader.h" 
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliFidutialCut.h"

ClassImp(AliCaloTrackESDReader)

//____________________________________________________________________________
AliCaloTrackESDReader::AliCaloTrackESDReader() : 
AliCaloTrackReader()
{
	//Default Ctor
	
	//Initialize parameters
	fDataType=kESD;
	
}

//____________________________________________________________________________
AliCaloTrackESDReader::AliCaloTrackESDReader(const AliCaloTrackESDReader & g) :   
  AliCaloTrackReader(g)
{
  // cpy ctor
}

//_________________________________________________________________________
//AliCaloTrackESDReader & AliCaloTrackESDReader::operator = (const AliCaloTrackESDReader & source)
//{
//	// assignment operator
//	
//	if(&source == this) return *this;
//	
//	return *this;
//	
//}

//____________________________________________________________________________
void AliCaloTrackESDReader::FillInputCTS() {
  //Return array with CTS tracks
  //TRefArray * fAODCTS = new TRefArray();
  Int_t nTracks   = fInputEvent->GetNumberOfTracks() ;
  Int_t naod = 0;
  Double_t pos[3];
  Double_t p[3];
  Double_t covTr[21];
  Double_t pid[10];
  for (Int_t itrack =  0; itrack <  nTracks; itrack++) {////////////// track loop
    AliESDtrack * track = (AliESDtrack*) ((AliESDEvent*)fInputEvent)->GetTrack(itrack) ; // retrieve track from esd
    
    //We want tracks fitted in the detectors:
    ULong_t status=AliESDtrack::kTPCrefit;
    status|=AliESDtrack::kITSrefit;
    
    if ( (track->GetStatus() & status) == status) {//Check if the bits we want are set
      
      track->GetPxPyPz(p) ;
      TLorentzVector momentum(p[0],p[1],p[2],0);
      
      if(fCTSPtMin < momentum.Pt() &&fFidutialCut->IsInFidutialCut(momentum,"CTS")){
	
	if(fDebug > 3 && momentum.Pt() > 0.2) printf("AliCaloTrackESDReader::FillInputCTS() - Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						    momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	
	track->GetXYZ(pos);
	track->GetCovarianceXYZPxPyPz(covTr);
	track->GetESDpid(pid);
	
	Float_t impactXY, impactZ;
	
	track->GetImpactParameters(impactXY,impactZ);
	
	if (impactXY<3) {
	  // track inside the beam pipe
	  //Put new aod object in file in AOD tracks array
	  AliAODTrack *aodTrack = new((*(fOutputEvent->GetTracks()))[naod++]) 
	    AliAODTrack(track->GetID(), track->GetLabel(), p, kTRUE, pos, kFALSE,covTr, (Short_t)track->GetSign(), track->GetITSClusterMap(), 
			pid,
			0x0,//primary,
			kTRUE, // check if this is right
			kTRUE, // check if this is right
			AliAODTrack::kPrimary, 
			0);
	  
	  aodTrack->SetFlags(track->GetStatus());
	  aodTrack->ConvertAliPIDtoAODPID();
	  
	  fAODCTS->Add(aodTrack);					}
	else continue;   // outside the beam pipe: orphan track	
      }//Pt and Fidutial cut passed. 
    }// track status
  }// track loop
  
  //Put references to selected tracks in array
  for(Int_t itrack = 0; itrack < (fOutputEvent->GetTracks())->GetEntriesFast(); itrack++){
    AliAODTrack * track =  (AliAODTrack*) (fOutputEvent->GetTracks())->At(itrack);	
    fAODCTS->Add(track);				
  }	
  
  if(fDebug > 1) printf("AliCaloTrackESDReader::FillInputCTS() - aod entries %d\n", fAODCTS->GetEntriesFast());
}

//____________________________________________________________________________
void AliCaloTrackESDReader::FillInputEMCAL() {
  //Return array with EMCAL clusters in aod format
  
  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  GetVertex(v);
  
  Float_t pos[3] ;
  Int_t naod      = (fOutputEvent->GetCaloClusters())->GetEntriesFast();
  Int_t nTracks   = fInputEvent->GetNumberOfTracks() ;
  
  //Loop to select clusters in fidutial cut and fill container with aodClusters
  for (Int_t iclus = 0; iclus < ((AliESDEvent*)fInputEvent)->GetNumberOfCaloClusters(); iclus++) {
    AliESDCaloCluster * clus = 0;
    if ( (clus = ((AliESDEvent*)fInputEvent)->GetCaloCluster(iclus)) ) {
      if (clus->IsEMCAL()){
	TLorentzVector momentum ;
	clus->GetMomentum(momentum, v);      
	if(fDebug > 3 && momentum.E() > 0.2) printf("AliCaloTrackESDReader::FillInputEMCAL() - all clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						   momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta()); 
	if(fEMCALPtMin < momentum.Pt() && fFidutialCut->IsInFidutialCut(momentum,"EMCAL")){
	  
	  if(fDebug > 2 && momentum.E() > 0.2) printf("AliCaloTrackESDReader::FillInputEMCAL() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						     momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	  
	  clus->GetPosition(pos) ;
	  Int_t id = clus->GetID();
	  Int_t nLabel = clus->GetNLabels();
	  Int_t *labels=0x0;
	  if(clus->GetLabels()) labels =  (clus->GetLabels())->GetArray();
	  
	  Float_t energy = clus->E();
	  Char_t ttype= AliAODCluster::kEMCALClusterv1;
	  
	  //Put new aod object in file in AOD calo clusters array
	  AliAODCaloCluster *caloCluster = new((*(fOutputEvent->GetCaloClusters()))[naod++]) 
	    AliAODCaloCluster(id,nLabel,labels,energy, pos, NULL,ttype,0);
	  
	  //       printf("Reader PID ESD: ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f,pi %0.2f, k %0.2f, p %0.2f, k0 %0.2f, n %0.2f, mu %0.2f \n",
	  // 	     pid[AliPID::kPhoton],pid[AliPID::kPi0],pid[AliPID::kElectron],pid[AliPID::kEleCon],pid[AliPID::kPion],
	  // 	     pid[AliPID::kKaon],pid[AliPID::kProton], pid[AliPID::kKaon0],pid[AliPID::kNeutron], pid[AliPID::kMuon]);
	  caloCluster->SetPIDFromESD(clus->GetPid());
	  caloCluster->SetCaloCluster(clus->GetDistanceToBadChannel(), clus->GetClusterDisp(), 
				      clus->GetM20(), clus->GetM02(),   
				      clus->GetEmcCpvDistance(),  clus->GetNExMax(), clus->GetTOF()) ;
	  
	  
	  if(fDebug > 3 && momentum.E() > 0.2)
	    printf("AliCaloTrackESDReader::FillInputEMCAL() - Selected clusters Distance BC %2.2f, dispersion %2.2f, M20 %f, M02 %3.2f, NexMax %d, TOF %e\n",
		   clus->GetDistanceToBadChannel(), clus->GetClusterDisp(),clus->GetM20(), clus->GetM02(),
		   clus->GetNExMax(), clus->GetTOF());
	  
	  caloCluster->SetNCells(clus->GetNCells());
	  caloCluster->SetCellsAbsId(clus->GetCellsAbsId());
	  caloCluster->SetCellsAmplitudeFraction(clus->GetCellsAmplitudeFraction());
	  
	  TArrayI* matchedT = 	clus->GetTracksMatched();
	  if (nTracks > 0 && matchedT && clus->GetTrackMatched() >= 0) {	
	    for (Int_t im = 0; im < matchedT->GetSize(); im++) {
	      Int_t iESDtrack = matchedT->At(im);
	      if(iESDtrack < nTracks && iESDtrack > -1)
		caloCluster->AddTrackMatched((fInputEvent->GetTrack(iESDtrack)));
	    }
	  }
	  //Fill reference array
	}//Pt and Fidutial cut passed.
      }//EMCAL cluster
    }//cluster exists
  }//esd cluster loop
  
  //Put references to selected clusters in array
  for(Int_t iclus = 0; iclus < (fOutputEvent->GetCaloClusters())->GetEntriesFast(); iclus++){
    AliAODCaloCluster * clus =  (AliAODCaloCluster*) (fOutputEvent->GetCaloClusters())->At(iclus);	
    fAODEMCAL->Add(clus);				
  }
  if(fDebug > 1) printf("AliCaloTrackESDReader::FillInputEMCAL() - aod entries %d\n", fAODEMCAL->GetEntriesFast());
  
}

//____________________________________________________________________________
void AliCaloTrackESDReader::FillInputPHOS() {
  //Return array with PHOS clusters in aod format
  Int_t nEMCAL = (fOutputEvent->GetCaloClusters())->GetEntriesFast();
  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  GetVertex(v);
  
  Float_t pos[3] ;
  Double_t * pid = new Double_t[AliPID::kSPECIESN];
  Int_t naod      = (fOutputEvent->GetCaloClusters())->GetEntriesFast();
  Int_t nTracks   = fInputEvent->GetNumberOfTracks() ;

  //Loop to select clusters in fidutial cut and fill container with aodClusters
  for (Int_t iclus = 0; iclus < ((AliESDEvent*)fInputEvent)->GetNumberOfCaloClusters(); iclus++) {
    AliESDCaloCluster * clus = 0;
    if ( (clus = ((AliESDEvent*)fInputEvent)->GetCaloCluster(iclus)) ) {
      if (clus->IsPHOS()){
	TLorentzVector momentum ;
	clus->GetMomentum(momentum, v);      
	if(fDebug > 3 && momentum.E() > 0.1)printf("AliCaloTrackESDReader::FillInputPHOS() - all clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						   momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	if(fPHOSPtMin < momentum.Pt() && fFidutialCut->IsInFidutialCut(momentum,"PHOS")){
	  
	  if(fDebug > 2 && momentum.E() > 0.1)printf("AliCaloTrackESDReader::FillInputPHOS() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						     momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	  
	  clus->GetPosition(pos) ;
	  Int_t id = clus->GetID();
	  Int_t nLabel = clus->GetNLabels();
	  Int_t *labels=0x0;
	  if(clus->GetLabels()) labels =  (clus->GetLabels())->GetArray();
	  Float_t energy = clus->E();
	  
	  //Phos cluster type
	  pid = clus->GetPid();
	  // printf("Reader PID ESD: ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f,pi %0.2f, k %0.2f, p %0.2f, k0 %0.2f, n %0.2f, mu %0.2f \n",
	  // 	 pid[AliPID::kPhoton],pid[AliPID::kPi0],pid[AliPID::kElectron],pid[AliPID::kEleCon],pid[AliPID::kPion],
	  // 	 pid[AliPID::kKaon],pid[AliPID::kProton], pid[AliPID::kKaon0],pid[AliPID::kNeutron], pid[AliPID::kMuon]);
	  Char_t ttype= AliAODCluster::kPHOSNeutral;
	  Float_t wNeutral = pid[AliPID::kNeutron]+ pid[AliPID::kKaon0]+pid[AliPID::kPhoton]+pid[AliPID::kPi0];
	  Float_t wCharged = pid[AliPID::kMuon] + pid[AliPID::kElectron] + pid[AliPID::kEleCon]+ 
	    pid[AliPID::kProton]+pid[AliPID::kKaon]+pid[AliPID::kPion];
	  if( wCharged > wNeutral)  ttype= AliAODCluster::kPHOSCharged;
	  
	  //Put new aod object in file in AOD calo clusters array
	  AliAODCaloCluster *caloCluster = new((*(fOutputEvent->GetCaloClusters()))[naod++]) 
	    AliAODCaloCluster(id,nLabel,labels,energy, pos, NULL, ttype, 0);


	  caloCluster->SetPIDFromESD(pid);   
	  caloCluster->SetCaloCluster(clus->GetDistanceToBadChannel(), clus->GetClusterDisp(), 
				      clus->GetM20(), clus->GetM02(),  
				      clus->GetEmcCpvDistance(),  clus->GetNExMax()) ;
	  
	  if(fDebug > 3 && momentum.E() > 0.2)
	    printf("AliCaloTrackESDReader::FillInputPHOS() - Selected clusters Distance BC %2.2f, dispersion %2.2f, M20 %f, M02 %3.2f, EmcCpvDist %3.3f, NexMax %d, TOF %e\n",
		   clus->GetDistanceToBadChannel(), clus->GetClusterDisp(),clus->GetM20(), clus->GetM02(),
		   clus->GetEmcCpvDistance(),  clus->GetNExMax(), clus->GetTOF());
	  
	  caloCluster->SetNCells(clus->GetNCells());
	  caloCluster->SetCellsAbsId(clus->GetCellsAbsId());
	  caloCluster->SetCellsAmplitudeFraction(clus->GetCellsAmplitudeFraction());

	  TArrayI* matchedT = 	clus->GetTracksMatched();
	  if (nTracks > 0 && matchedT && clus->GetTrackMatched() >= 0) {	
	    for (Int_t im = 0; im < matchedT->GetSize(); im++) {
	      Int_t iESDtrack = matchedT->At(im);
	      if(iESDtrack < nTracks && iESDtrack > -1)
		caloCluster->AddTrackMatched((fInputEvent->GetTrack(iESDtrack)));
	    }
	  }
	}//Pt and Fidutial cut passed.
      }//PHOS cluster
    }//cluster exists
  }//esd cluster loop
  
  //Put references to selected clusters in array
  for(Int_t iclus = nEMCAL; iclus < (fOutputEvent->GetCaloClusters())->GetEntriesFast(); iclus++){
    AliAODCaloCluster * clus =  (AliAODCaloCluster*) (fOutputEvent->GetCaloClusters())->At(iclus);	
    fAODPHOS->Add(clus);				
  }	
  if(fDebug > 1) printf("FillInputPHOS():: aod entries %d\n", fAODPHOS->GetEntriesFast());
  
}

//____________________________________________________________________________
void AliCaloTrackESDReader::FillInputEMCALCells() {
  //Return array with EMCAL cells in esd format
  
  fEMCALCells = (TNamed*) ((AliESDEvent*)fInputEvent)->GetEMCALCells(); 
  
}

//____________________________________________________________________________
void AliCaloTrackESDReader::FillInputPHOSCells() {
  //Return array with PHOS cells in esd format
  
  fPHOSCells = (TNamed*) ((AliESDEvent*)fInputEvent)->GetPHOSCells(); 
  
}

//____________________________________________________________________________
void AliCaloTrackESDReader::GetVertex(Double_t  v[3]) const {
  //Return vertex position
  
  ((AliESDEvent*)fInputEvent)->GetVertex()->GetXYZ(v) ;
  
}


//____________________________________________________________________________
void AliCaloTrackESDReader::SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) {
  // Connect the data pointers
  
  if(strcmp(esd->GetName(),"AliESDEvent")){
    printf("AliCaloTrackESDReader::SetInputOutputMCEvent() - ABORT::Wrong reader, here only ESDs. Input name: %s != AliESDEvent \n",esd->GetName());
    abort();
  }
  
  SetInputEvent(esd);
  SetOutputEvent(aod);
  SetMC(mc);
  
}
