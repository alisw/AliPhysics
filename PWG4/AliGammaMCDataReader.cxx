
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
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class for reading data (Kinematics and ESDs) in order to do prompt gamma correlations
//  Class created from old AliPHOSGammaJet 
//  (see AliRoot versions previous Release 4-09)
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "Riostream.h"
#include <TParticle.h>

//---- ANALYSIS system ----
#include "AliGammaMCDataReader.h" 
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDCaloCluster.h"
#include "AliStack.h"

ClassImp(AliGammaMCDataReader)

//____________________________________________________________________________
  AliGammaMCDataReader::AliGammaMCDataReader() : 
    AliGammaReader(),  
    fEMCALPID(0),fPHOSPID(0),
    fEMCALPhotonWeight(0.), fEMCALPi0Weight(0.), fPHOSPhotonWeight(0.),  fPHOSPi0Weight(0.) 
{
  //Default Ctor
  
  //Initialize parameters
  fDataType=kMCData;
  InitParameters();
  
}

//____________________________________________________________________________
AliGammaMCDataReader::AliGammaMCDataReader(const AliGammaMCDataReader & g) :   
  AliGammaReader(g),
  fEMCALPID(g.fEMCALPID), 
  fPHOSPID(g.fPHOSPID),
  fEMCALPhotonWeight(g.fEMCALPhotonWeight), 
  fEMCALPi0Weight(g.fEMCALPi0Weight), 
  fPHOSPhotonWeight(g.fPHOSPhotonWeight),
  fPHOSPi0Weight(g.fPHOSPi0Weight)
{
  // cpy ctor
}

//_________________________________________________________________________
AliGammaMCDataReader & AliGammaMCDataReader::operator = (const AliGammaMCDataReader & source)
{
  // assignment operator
  
  if(&source == this) return *this;
  
  fEMCALPID = source.fEMCALPID ;
  fPHOSPID = source.fPHOSPID ;
  fEMCALPhotonWeight = source. fEMCALPhotonWeight ;
  fEMCALPi0Weight = source.fEMCALPi0Weight ;
  fPHOSPhotonWeight = source.fPHOSPhotonWeight ;
  fPHOSPi0Weight = source.fPHOSPi0Weight ;
  
  return *this;
  
}

//____________________________________________________________________________
void AliGammaMCDataReader::CreateParticleList(TObject * data, TObject * kine,
					      TClonesArray * plCTS, 
					      TClonesArray * plEMCAL,  
					      TClonesArray * plPHOS, TClonesArray *){
  
  //Create a list of particles from the ESD. These particles have been measured 
  //by the Central Tracking system (TPC+ITS), PHOS and EMCAL 

  AliESDEvent* esd = (AliESDEvent*) data;
  AliStack* stack = (AliStack*) kine;
  
  Int_t npar  = 0 ;
  Float_t *pid = new Float_t[AliPID::kSPECIESN];  
  AliDebug(3,"Fill particle lists");
  
  //Get vertex for momentum calculation  
  Double_t v[3] ; //vertex ;
  esd->GetVertex()->GetXYZ(v) ; 
  
  //########### PHOS ##############
  
  Int_t begphos = esd->GetFirstPHOSCluster();  
  Int_t endphos = esd->GetFirstPHOSCluster() + 
    esd->GetNumberOfPHOSClusters() ;  
  Int_t indexPH = plPHOS->GetEntries() ;
  AliDebug(3,Form("First PHOS particle %d, last particle %d", begphos,endphos));
  
  for (npar =  begphos; npar <  endphos; npar++) {//////////////PHOS track loop
    AliESDCaloCluster * clus = esd->GetCaloCluster(npar) ; // retrieve track from esd
    if(clus->GetTrackMatched()==-1 ){
      //Create a TParticle to fill the particle list
      TLorentzVector momentum ;
      clus->GetMomentum(momentum, v);
      pid=clus->GetPid();	
      Int_t pdg = 22;

      if(fPHOSPID){
	if( pid[AliPID::kPhoton] > fPHOSPhotonWeight) pdg=22;
	if( pid[AliPID::kPi0] > fPHOSPi0Weight) pdg=111;
	else pdg = 0 ;
      }
      
      TParticle * particle = new TParticle(pdg, 1, -1, -1, -1, -1, 
					   momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E(), v[0], v[1], v[2], 0);

      AliDebug(4,Form("PHOS clusters: pt %f, phi %f, eta %f", particle->Pt(),particle->Phi(),particle->Eta()));
      
      //Look  if parent is prompt photon
      Int_t label = clus->GetLabel();
      if(label>=0){
	TParticle * pmother = stack->Particle(label);
	Int_t imother = pmother->GetFirstMother(); 
	if(imother == 6 || imother == 7)
	  pmother->SetFirstMother(22);
      }
      
      new((*plPHOS)[indexPH++])   TParticle(*particle) ;
    
    }//not charged
  }//cluster loop

  //########### CTS (TPC+ITS) #####################
  Int_t begtpc   = 0 ;  
  Int_t endtpc   = esd->GetNumberOfTracks() ;
  Int_t indexCh  = plCTS->GetEntries() ;
  AliDebug(3,Form("First CTS particle %d, last particle %d", begtpc,endtpc));
  
  for (npar =  begtpc; npar <  endtpc; npar++) {////////////// track loop
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
 
      new((*plCTS)[indexCh++])       TParticle(*particle) ;    
    }
  }
  
  //################ EMCAL ##############
  
  Int_t indexEM  = plEMCAL->GetEntries() ; 
  Int_t begem = esd->GetFirstEMCALCluster();  
  Int_t endem = esd->GetFirstEMCALCluster() + 
    esd->GetNumberOfEMCALClusters() ;  
  
  AliDebug(3,Form("First EMCAL particle %d, last particle %d",begem,endem));
  
  for (npar =  begem; npar <  endem; npar++) {//////////////EMCAL track loop
    AliESDCaloCluster * clus = esd->GetCaloCluster(npar) ; // retrieve track from esd
    Int_t clustertype= clus->GetClusterType();
    if(clustertype == AliESDCaloCluster::kClusterv1 && clus->GetTrackMatched()==-1 ){
      
      TLorentzVector momentum ;
      clus->GetMomentum(momentum, v);            
      pid=clus->GetPid();	
      Int_t pdg = 22;

      if(fEMCALPID){
	if( pid[AliPID::kPhoton] > fEMCALPhotonWeight) pdg=22;
	else if( pid[AliPID::kPi0] > fEMCALPi0Weight) pdg=111;
	else pdg = 0;
      }
      
      TParticle * particle = new TParticle(pdg, 1, -1, -1, -1, -1, 
					   momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E(), v[0], v[1], v[2], 0);
      AliDebug(4,Form("EMCAL clusters: pt %f, phi %f, eta %f", particle->Pt(),particle->Phi(),particle->Eta()));

      //Look  if parent is prompt photon
      Int_t label = clus->GetLabel();
      if(label>=0){
	TParticle * pmother = stack->Particle(label);
	Int_t imother = pmother->GetFirstMother(); 
	if(imother == 6 || imother == 7)
	  pmother->SetFirstMother(22);
      }
      
      new((*plEMCAL)[indexEM++])   TParticle(*particle) ;
      
    }//not charged, not pseudocluster
  }//Cluster loop
  
  AliDebug(3,"Particle lists filled");
  
}

  //____________________________________________________________________________
void AliGammaMCDataReader::InitParameters()
{
 
  //Initialize the parameters of the analysis.

  //Fill particle lists when PID is ok
  fEMCALPID = kFALSE;
  fPHOSPID = kFALSE;
  fEMCALPhotonWeight = 0.5 ;
  fEMCALPi0Weight = 0.5 ;
  fPHOSPhotonWeight = 0.8 ;
  fPHOSPi0Weight = 0.5 ;

}


void AliGammaMCDataReader::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  printf("PHOS PID on?               =     %d\n",  fPHOSPID) ; 
  printf("EMCAL PID  on?         =     %d\n",  fEMCALPID) ;
  printf("PHOS PID weight , photon  %f, pi0 %f\n\n",  fPHOSPhotonWeight,  fPHOSPi0Weight) ; 
  printf("EMCAL PID weight, photon %f, pi0 %f\n",   fEMCALPhotonWeight,  fEMCALPi0Weight) ; 

}
