
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
// It is a filtering class, transforms ESD tracks or CaloClusters
// into AOD tracks and calocluters, which are the basic input of the analysis
// classes in this frame.
// It is recommended to use the official filter AliAnalysisTaskESDfilter, and 
// then the reader for AODs AliCaloTrackAODReader.
//
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
//#include "Riostream.h"

//---- ANALYSIS system ----
#include "AliCaloTrackESDReader.h" 
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliAODEvent.h"
#include "AliFiducialCut.h"

ClassImp(AliCaloTrackESDReader)

//____________________________________________________________________________
AliCaloTrackESDReader::AliCaloTrackESDReader() : 
AliCaloTrackReader()
{
	//Default Ctor
	
	//Initialize parameters
	fDataType=kESD;
	fReadStack          = kTRUE;
	fReadAODMCParticles = kFALSE;
	//We want tracks fitted in the detectors:
	fTrackStatus=AliESDtrack::kTPCrefit;
    fTrackStatus|=AliESDtrack::kITSrefit;

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
  if(fDebug > 2 ) printf("AliCaloTrackESDReader::FillInputCTS()\n");

  //TObjArray * fAODCTS = new TObjArray();
  Int_t nTracks   = fInputEvent->GetNumberOfTracks() ;
  Int_t naod = 0;
  Double_t pos[3];
  Double_t p[3];
  Double_t covTr[21];
  Double_t pid[10];
  Double_t bfield = ((AliESDEvent*)fInputEvent)->GetMagneticField();

  Double_t      timezero        = 0;   //TO BE FIXED

  //To be replaced by call to AliEMCALGeoUtils when the class becomes available
  Double_t radius = 441.0; //[cm] EMCAL radius +13cm
  if(fDebug > 1) printf("AliCaloTrackESDReader::FillInputCTS() - org entries %d\n", nTracks);
  for (Int_t itrack =  0; itrack <  nTracks; itrack++) {////////////// track loop
    AliESDtrack * track = (AliESDtrack*) ((AliESDEvent*)fInputEvent)->GetTrack(itrack) ; // retrieve track from esd
    
    //We want tracks fitted in the detectors: TPCrefit, ITSrefit ... check the set bits.
	if (fTrackStatus && !((track->GetStatus() & fTrackStatus) == fTrackStatus)) continue ;
      
	track->GetPxPyPz(p) ;
	TLorentzVector momentum(p[0],p[1],p[2],0);
      
	if(fCTSPtMin < momentum.Pt() &&fFiducialCut->IsInFiducialCut(momentum,"CTS")){
	
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

		  //fill needed AliAODPid information, including
		  //Extrapolate track to EMCAL surface for AOD-level track-cluster matching
		  AliAODPid *aodpid = new AliAODPid;
		  aodpid->SetITSsignal(track->GetITSsignal());
		  aodpid->SetTPCsignal(track->GetTPCsignal());
		  //n TRD planes = 6
		  Int_t nslices = track->GetNumberOfTRDslices()*6;
		  Double_t *trdslices = new Double_t[nslices];
		  for(Int_t iSl =0; iSl < track->GetNumberOfTRDslices(); iSl++) {
		    for(Int_t iPl =0; iPl<6; iPl++) trdslices[iPl*track->GetNumberOfTRDslices()+iSl] = track->GetTRDslice(iPl,iSl);
		  }
		  aodpid->SetTRDsignal(track->GetNumberOfTRDslices()*6,trdslices);
		  Double_t times[AliAODPid::kSPECIES]; track->GetIntegratedTimes(times);
		  aodpid->SetIntegratedTimes(times);

		  aodpid->SetTOFsignal(track->GetTOFsignal()-timezero); // to be fixed 
		  aodpid->SetHMPIDsignal(track->GetHMPIDsignal());

		  Double_t emcpos[3] = {0.,0.,0.};
		  Double_t emcmom[3] = {0.,0.,0.};
		  aodpid->SetEMCALPosition(emcpos);
		  aodpid->SetEMCALMomentum(emcmom);
		  
		  AliExternalTrackParam *outerparam = (AliExternalTrackParam*)track->GetOuterParam();
		  if(!outerparam) continue;
		  
		  Bool_t okpos = outerparam->GetXYZAt(radius,bfield,emcpos);
		  Bool_t okmom = outerparam->GetPxPyPzAt(radius,bfield,emcmom);
		  if(!(okpos && okmom)) continue;
		  
		  aodpid->SetEMCALPosition(emcpos);
		  aodpid->SetEMCALMomentum(emcmom);
		  
		  aodTrack->SetDetPID(aodpid);
		}
		else continue;   // outside the beam pipe: orphan track	
	}//Pt and Fiducial cut passed. 
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
  if(fDebug > 2 ) printf("AliCaloTrackESDReader::FillInputEMCAL()\n");
	
  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  GetVertex(v);
  
  Float_t pos[3] ;
  Int_t naod      = (fOutputEvent->GetCaloClusters())->GetEntriesFast();
  Int_t nTracks   = fInputEvent->GetNumberOfTracks() ;
  
  //Loop to select clusters in fiducial cut and fill container with aodClusters
  for (Int_t iclus = 0; iclus < ((AliESDEvent*)fInputEvent)->GetNumberOfCaloClusters(); iclus++) {
    AliESDCaloCluster * clus = 0;
    if ( (clus = ((AliESDEvent*)fInputEvent)->GetCaloCluster(iclus)) ) {
      if (clus->IsEMCAL()){
		  
	//Check if the cluster contains any bad channel and if close to calorimeter borders
	if( GetCaloUtils()->ClusterContainsBadChannel("EMCAL",clus->GetCellsAbsId(), clus->GetNCells())) continue;
	if(!GetCaloUtils()->CheckCellFiducialRegion(clus, ((AliESDEvent*)fInputEvent)->GetEMCALCells())) continue;

	TLorentzVector momentum ;
	clus->GetMomentum(momentum, v);      
	if(fDebug > 3 && momentum.E() > 0.1) printf("AliCaloTrackESDReader::FillInputEMCAL() - all clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						   momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta()); 
	if(fEMCALPtMin < momentum.Pt() && fFiducialCut->IsInFiducialCut(momentum,"EMCAL")){
	  
	  if(fDebug > 2 && momentum.E() > 0.1) printf("AliCaloTrackESDReader::FillInputEMCAL() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						     momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	  
	  clus->GetPosition(pos) ;
	  Int_t id = clus->GetID();
	  Int_t nLabel = clus->GetNLabels();
	  Int_t *labels=0x0;
	  if(clus->GetLabels()) labels =  (clus->GetLabels())->GetArray();
	  
	  Float_t energy = clus->E();
	  Char_t ttype= AliAODCluster::kEMCALClusterv1;

	  //Recalibrate the cluster energy 
		if(GetCaloUtils()->IsRecalibrationOn())	energy = GetCaloUtils()->RecalibrateClusterEnergy(clus, (AliESDCaloCells*)GetEMCALCells());
		
		
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
	  
	  
	  if(fDebug > 3 && momentum.E() > 0.1)
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
	}//Pt and Fiducial cut passed.
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
  if(fDebug > 2 ) printf("AliCaloTrackESDReader::FillInputPHOS()\n");

  Int_t nEMCAL = (fOutputEvent->GetCaloClusters())->GetEntriesFast();
  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  GetVertex(v);
  
  Float_t pos[3] ;
  Double_t * pid = new Double_t[AliPID::kSPECIESN];
  Int_t naod      = (fOutputEvent->GetCaloClusters())->GetEntriesFast();
  Int_t nTracks   = fInputEvent->GetNumberOfTracks() ;

  //Loop to select clusters in fiducial cut and fill container with aodClusters
  for (Int_t iclus = 0; iclus < ((AliESDEvent*)fInputEvent)->GetNumberOfCaloClusters(); iclus++) {
    AliESDCaloCluster * clus = 0;
    if ( (clus = ((AliESDEvent*)fInputEvent)->GetCaloCluster(iclus)) ) {
      if (clus->IsPHOS()){
		  
	//Check if the cluster contains any bad channel and if close to calorimeter borders
	if( GetCaloUtils()->ClusterContainsBadChannel("PHOS",clus->GetCellsAbsId(), clus->GetNCells())) continue;
	if(!GetCaloUtils()->CheckCellFiducialRegion(clus, ((AliESDEvent*)fInputEvent)->GetPHOSCells())) continue;

	TLorentzVector momentum ;
	clus->GetMomentum(momentum, v);      
	if(fDebug > 3 && momentum.E() > 0.1)printf("AliCaloTrackESDReader::FillInputPHOS() - all clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
						   momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	if(fPHOSPtMin < momentum.Pt() && fFiducialCut->IsInFiducialCut(momentum,"PHOS")){
	  
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
	  
	  //Recalibrate the cluster energy 
	  if(GetCaloUtils()->IsRecalibrationOn()) energy = GetCaloUtils()->RecalibrateClusterEnergy(clus, (AliESDCaloCells*)GetPHOSCells());
		
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
	}//Pt and Fiducial cut passed.
      }//PHOS cluster
    }//cluster exists
  }//esd cluster loop
  
  //Put references to selected clusters in array
  for(Int_t iclus = nEMCAL; iclus < (fOutputEvent->GetCaloClusters())->GetEntriesFast(); iclus++){
    AliAODCaloCluster * clus =  (AliAODCaloCluster*) (fOutputEvent->GetCaloClusters())->At(iclus);	
    fAODPHOS->Add(clus);				
  }	
  if(fDebug > 1) printf("AliCaloTrackESDReader::FillInputPHOS() - aod entries %d\n", fAODPHOS->GetEntriesFast());
  
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
  
  //((AliESDEvent*)fInputEvent)->GetVertex()->GetXYZ(v) ;//SPD
  ((AliESDEvent*)fInputEvent)->GetPrimaryVertex()->GetXYZ(v);
	
}


//____________________________________________________________________________
Double_t AliCaloTrackESDReader::GetBField() const {
  //Return magnetic field

  Double_t bfield = ((AliESDEvent*)fInputEvent)->GetMagneticField();

  return bfield;
}


//____________________________________________________________________________
void AliCaloTrackESDReader::SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) {
  // Connect the data pointers
  
  if(strcmp(esd->GetName(),"AliESDEvent")){
    printf("AliCaloTrackESDReader::SetInputOutputMCEvent() - STOP ::Wrong reader, here only ESDs. Input name: %s != AliESDEvent \n",esd->GetName());
    abort();
  }
  
  SetInputEvent(esd);
  SetOutputEvent(aod);
  SetMC(mc);
  
}
