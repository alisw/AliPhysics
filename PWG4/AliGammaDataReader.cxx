
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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class for reading data (ESDs) in order to do prompt gamma correlations
//  Class created from old AliPHOSGammaJet 
//  (see AliRoot versions previous Release 4-09)
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TParticle.h>
#include <TFormula.h>

//---- ANALYSIS system ----
#include "AliGammaDataReader.h" 
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDCaloCluster.h"
#include "AliLog.h"

ClassImp(AliGammaDataReader)

//____________________________________________________________________________
AliGammaDataReader::AliGammaDataReader() : 
  AliGammaReader()
{
  //Default Ctor
  
  //Initialize parameters
  fDataType=kData;
  
}

//____________________________________________________________________________
AliGammaDataReader::AliGammaDataReader(const AliGammaDataReader & g) :   
  AliGammaReader(g)
{
  // cpy ctor
}

//_________________________________________________________________________
AliGammaDataReader & AliGammaDataReader::operator = (const AliGammaDataReader & source)
{
  // assignment operator

  if(&source == this) return *this;

  return *this;

}

//____________________________________________________________________________
void AliGammaDataReader::CreateParticleList(TObject * data, TObject *,
					    TClonesArray * plCTS, 
					    TClonesArray * plEMCAL,  
					    TClonesArray * plPHOS, 
					    TClonesArray * , 
					    TClonesArray * ,  
					    TClonesArray * ){
  
  //Create a list of particles from the ESD. These particles have been measured 
  //by the Central Tracking system (TPC+ITS+...), PHOS and EMCAL 

  AliESDEvent* esd = (AliESDEvent*) data;

  Int_t npar  = 0 ;
  Double_t *pid = new Double_t[AliPID::kSPECIESN];  
   AliDebug(3,"Fill particle lists");
  
  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  esd->GetVertex()->GetXYZ(v) ; 
  
  //########### CALORIMETERS ##############
    
  Int_t nCaloCluster = esd->GetNumberOfCaloClusters() ;  
  Int_t indexPH = plPHOS->GetEntries() ;
  Int_t indexEM = plEMCAL->GetEntries() ;
  

  for (npar =  0; npar <  nCaloCluster; npar++) {//////////////CaloCluster loop
    AliESDCaloCluster * clus = esd->GetCaloCluster(npar) ; // retrieve cluster from esd
    Int_t type = clus->GetClusterType();
    
    //########### PHOS ##############
    if(fSwitchOnPHOS && type ==  AliESDCaloCluster::kPHOSCluster){
      AliDebug(4,Form("PHOS clusters: E %f, match %d", clus->E(),clus->GetTrackMatched()));
      
      if(clus->GetTrackMatched()==-1){
	TLorentzVector momentum ;
	clus->GetMomentum(momentum, v);      
	Double_t phi = momentum.Phi();
	if(phi<0) phi+=TMath::TwoPi() ;
	if(momentum.Pt() > fNeutralPtCut &&  TMath::Abs(momentum.Eta()) < fPHOSEtaCut &&
	   phi > fPhiPHOSCut[0] && phi < fPhiPHOSCut[1] ) {
	  
	  pid=clus->GetPid();	
	  Int_t pdg = 22;
	    
	  if(IsPHOSPIDOn()){
	    AliDebug(5,Form("E %1.2f; PID: ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f,pi %0.2f, k %0.2f, p %0.2f, k0 %0.2f, n %0.2f, mu %0.2f ",
			    momentum.E(),pid[AliPID::kPhoton],pid[AliPID::kPi0],pid[AliPID::kElectron],pid[AliPID::kEleCon],pid[AliPID::kPion],
			    pid[AliPID::kKaon],pid[AliPID::kProton], pid[AliPID::kKaon0],pid[AliPID::kNeutron], pid[AliPID::kMuon]));
	    
	    Float_t wPhoton =  fPHOSPhotonWeight;
	    Float_t wPi0 =  fPHOSPi0Weight;
	    
	    if(fPHOSWeightFormula){
	      wPhoton = fPHOSPhotonWeightFormula->Eval(momentum.E()) ;
	      wPi0 =    fPHOSPi0WeightFormula->Eval(momentum.E());
	    }
	    
	    if(pid[AliPID::kPhoton] > wPhoton) 
	      pdg = kPhoton ;
	    else if(pid[AliPID::kPi0] > wPi0) 
	      pdg = kPi0 ; 
	    else if(pid[AliPID::kElectron] > fPHOSElectronWeight)  
	      pdg = kElectron ;
	    else if(pid[AliPID::kEleCon] > fPHOSElectronWeight) 
	      pdg = kEleCon ;
	    else if(pid[AliPID::kPion]+pid[AliPID::kKaon]+pid[AliPID::kProton] > fPHOSChargeWeight) 
	      pdg = kChargedHadron ;  
	    else if(pid[AliPID::kKaon0]+pid[AliPID::kNeutron] > fPHOSNeutralWeight) 
	      pdg = kNeutralHadron ; 
	    
	    else if(pid[AliPID::kElectron]+pid[AliPID::kEleCon]+pid[AliPID::kPion]+pid[AliPID::kKaon]+pid[AliPID::kProton]  >  
		    pid[AliPID::kPhoton] + pid[AliPID::kPi0]+pid[AliPID::kKaon0]+pid[AliPID::kNeutron]) 
	      pdg = kChargedUnknown  ; 
	    else 
	      pdg = kNeutralUnknown ; 
	    //neutral cluster, unidentifed.
	  }
	  
	  if(pdg != kElectron && pdg != kEleCon && pdg !=kChargedHadron && pdg !=kChargedUnknown ){//keep only neutral particles in the array
	    TParticle * particle = new TParticle(pdg, 1, -1, -1, -1, -1, 
						 momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E(), v[0], v[1], v[2], 0);
	    
	    AliDebug(4,Form("PHOS added: pdg %d, pt %f, phi %f, eta %f", pdg, particle->Pt(),particle->Phi(),particle->Eta()));
	    
	    new((*plPHOS)[indexPH++])   TParticle(*particle) ;
	  }
	  else AliDebug(4,Form("PHOS charged cluster NOT added: pdg %d, pt %f, phi %f, eta %f\n", 
			       pdg, momentum.Pt(),momentum.Phi(),momentum.Eta()));	
	  
	}//pt, eta, phi cut
	else 	AliDebug(4,"Particle not added");
      }//track-match?
    }//PHOS cluster

    //################ EMCAL ##############
    else if(fSwitchOnEMCAL &&  type ==  AliESDCaloCluster::kEMCALClusterv1){
      AliDebug(4,Form("EMCAL clusters: E %f, match %d", clus->E(),clus->GetTrackMatched()));
      
      if(clus->GetTrackMatched()==-1 ){
	TLorentzVector momentum ;
	clus->GetMomentum(momentum, v); 
	Double_t phi = momentum.Phi();
	if(phi<0) phi+=TMath::TwoPi() ;
	if(momentum.Pt() > fNeutralPtCut &&  TMath::Abs(momentum.Eta()) < fEMCALEtaCut &&
	   phi > fPhiEMCALCut[0] && phi < fPhiEMCALCut[1] ) {
	  
	  pid=clus->GetPid();	
	  Int_t pdg = 22;
	  
	  if(IsEMCALPIDOn()){
	    AliDebug(5,Form("E %1.2f; PID: ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f,pi %0.2f, k %0.2f, p %0.2f, k0 %0.2f, n %0.2f, mu %0.2f ",
			    momentum.E(),pid[AliPID::kPhoton],pid[AliPID::kPi0],pid[AliPID::kElectron],pid[AliPID::kEleCon],pid[AliPID::kPion],
			    pid[AliPID::kKaon],pid[AliPID::kProton], pid[AliPID::kKaon0],pid[AliPID::kNeutron], pid[AliPID::kMuon]));
	    
	    if(pid[AliPID::kPhoton] > fEMCALPhotonWeight) 
	      pdg = kPhoton ;
	    else if(pid[AliPID::kPi0] > fEMCALPi0Weight) 
	      pdg = kPi0 ; 
	    else if(pid[AliPID::kElectron] > fEMCALElectronWeight)  
	      pdg = kElectron ;
	    else if(pid[AliPID::kEleCon] > fEMCALElectronWeight) 
	      pdg = kEleCon ;
	    else if(pid[AliPID::kPion]+pid[AliPID::kKaon]+pid[AliPID::kProton] > fEMCALChargeWeight) 
	      pdg = kChargedHadron ;  
	    else if(pid[AliPID::kKaon0]+pid[AliPID::kNeutron] > fEMCALNeutralWeight) 
	      pdg = kNeutralHadron ; 
	    else if(pid[AliPID::kElectron]+pid[AliPID::kEleCon]+pid[AliPID::kPion]+pid[AliPID::kKaon]+pid[AliPID::kProton]  >  
		    pid[AliPID::kPhoton] + pid[AliPID::kPi0]+pid[AliPID::kKaon0]+pid[AliPID::kNeutron]) 
	      pdg = kChargedUnknown ; 
	    else 
	      pdg = kNeutralUnknown ;
	  }
	  
	  if(pdg != kElectron && pdg != kEleCon && pdg !=kChargedHadron && pdg !=kChargedUnknown){//keep only neutral particles in the array
	    
	    TParticle * particle = new TParticle(pdg, 1, -1, -1, -1, -1, 
						 momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E(), v[0], v[1], v[2], 0);
	    AliDebug(4,Form("EMCAL cluster added: pdg %f, pt %f, phi %f, eta %f", pdg, particle->Pt(),particle->Phi(),particle->Eta()));
	    
	    new((*plEMCAL)[indexEM++])   TParticle(*particle) ;
	  }
	  else AliDebug(4,Form("EMCAL charged cluster NOT added: pdg %d, pt %f, phi %f, eta %f", 
			       pdg, momentum.Pt(),momentum.Phi(),momentum.Eta()));
	  
	}//pt, phi, eta cut
	else 	AliDebug(4,"Particle not added");
      }//track-matched
    }//EMCAL cluster

  }//cluster loop


  //########### CTS (TPC+ITS) #####################
  Int_t nTracks   = esd->GetNumberOfTracks() ;
  Int_t indexCh  = plCTS->GetEntries() ;
  
  if(fSwitchOnCTS){
    AliDebug(3,Form("Number of tracks %d",nTracks));
    
    for (npar =  0; npar <  nTracks; npar++) {////////////// track loop
      AliESDtrack * track = esd->GetTrack(npar) ; // retrieve track from esd
      
      //We want tracks fitted in the detectors:
      ULong_t status=AliESDtrack::kTPCrefit;
      status|=AliESDtrack::kITSrefit;
      
      //We want tracks whose PID bit is set:
      //     ULong_t status =AliESDtrack::kITSpid;
      //     status|=AliESDtrack::kTPCpid;
      
      if ( (track->GetStatus() & status) == status) {//Check if the bits we want are set
	// Do something with the tracks which were successfully
	// re-fitted 
	Double_t en = 0; //track ->GetTPCsignal() ;
	Double_t mom[3];
	track->GetPxPyPz(mom) ;
	Double_t px = mom[0];
	Double_t py = mom[1];
	Double_t pz = mom[2]; //Check with TPC people if this is correct.
	Int_t pdg = 11; //Give any charged PDG code, in this case electron.
	//I just want to tag the particle as charged
	TParticle * particle = new TParticle(pdg, 1, -1, -1, -1, -1, 
					     px, py, pz, en, v[0], v[1], v[2], 0);
	
	//TParticle * particle = new TParticle() ;
	//particle->SetMomentum(px,py,pz,en) ;
	if(particle->Pt() > fChargedPtCut && TMath::Abs(particle->Eta())<fCTSEtaCut)
	  new((*plCTS)[indexCh++])       TParticle(*particle) ;    
      }// select track from refit
    }//track loop
  }//CTS
  
  AliDebug(3,Form("Particle lists filled, tracks  %d , clusters: EMCAL %d, PHOS %d", indexCh,indexEM,indexPH));
  
}

