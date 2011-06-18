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
/* $Id: AliCaloPID.cxx 21839 2007-10-29 13:49:42Z gustavo $ */

//_________________________________________________________________________
// Class for PID selection with calorimeters
// The Output of the 2 main methods GetPdg is a PDG number identifying the cluster, 
// being kPhoton, kElectron, kPi0 ... as defined in the header file
//   - GetPdg(const TString calo, const Double_t * pid, const Float_t energy)
//      Reads the PID weights array of the ESDs and depending on its magnitude identifies the particle
//   - GetPdg(const TString calo,const TLorentzVector mom, const AliVCluster * cluster)
//      Recalcultes PID, the bayesian or any new one to be implemented in the future
//      Right now only the possibility to recalculate EMCAL with bayesian and simple PID.
//      In order to recalculate Bayesian, it is necessary to load the EMCALUtils library
//      and do SwitchOnBayesianRecalculation().
//      To change the PID parameters from Low to High like the ones by default, use the constructor 
//      AliCaloPID(flux)
//      where flux is AliCaloPID::kLow or AliCaloPID::kHigh
//      If it is necessary to change the parameters use the constructor 
//      AliCaloPID(AliEMCALPIDUtils *utils) and set the parameters before.
//   - SetPIDBits: Simple PID, depending on the thresholds fDispCut fTOFCut and even the
//     result of the PID bayesian a different PID bit is set. 
//
//  All these methods can be called in the analysis you are interested.
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TMath.h>
#include <TString.h>
//#include <TFormula.h>

//---- ANALYSIS system ----
#include "AliCaloPID.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliAODPWG4Particle.h"
#include "AliEMCALPIDUtils.h"
#include "AliCalorimeterUtils.h"

ClassImp(AliCaloPID)


//________________________________________________
AliCaloPID::AliCaloPID() : 
TObject(), fEMCALPhotonWeight(0.), fEMCALPi0Weight(0.),  
fEMCALElectronWeight(0.),  fEMCALChargeWeight(0.),
fEMCALNeutralWeight(0.),
fPHOSPhotonWeight(0.), fPHOSPi0Weight(0.),  
fPHOSElectronWeight(0.), fPHOSChargeWeight(0.) , 
fPHOSNeutralWeight(0.), //fPHOSWeightFormula(0), 
//fPHOSPhotonWeightFormula(0x0), fPHOSPi0WeightFormula(0x0),
fDispCut(0.),fTOFCut(0.), fDebug(-1), 
fRecalculateBayesian(kFALSE), fParticleFlux(kLow), 
fEMCALPIDUtils(new AliEMCALPIDUtils),
fHistoNEBins(100),   fHistoEMax(100.),   fHistoEMin(0.),
fHistoNDEtaBins(100),fHistoDEtaMax(0.15),fHistoDEtaMin(-0.15),
fHistoNDPhiBins(100),fHistoDPhiMax(0.15),fHistoDPhiMin(-0.15),
fhTrackMatchedDEta(0x0),fhTrackMatchedDPhi(0x0),fhTrackMatchedDEtaDPhi(0x0)

{
  //Ctor
  
  //Initialize parameters
  InitParameters();
}

//________________________________________________
AliCaloPID::AliCaloPID(const Int_t flux) : 
TObject(), fEMCALPhotonWeight(0.), fEMCALPi0Weight(0.),  
fEMCALElectronWeight(0.),  fEMCALChargeWeight(0.),
fEMCALNeutralWeight(0.),
fPHOSPhotonWeight(0.), fPHOSPi0Weight(0.),  
fPHOSElectronWeight(0.), fPHOSChargeWeight(0.) , 
fPHOSNeutralWeight(0.), //fPHOSWeightFormula(0), 
//fPHOSPhotonWeightFormula(0x0), fPHOSPi0WeightFormula(0x0),
fDispCut(0.),fTOFCut(0.), fDebug(-1), 
fRecalculateBayesian(kTRUE), fParticleFlux(flux), 
fEMCALPIDUtils(new AliEMCALPIDUtils),
fHistoNEBins(100),   fHistoEMax(100.),   fHistoEMin(0.),
fHistoNDEtaBins(100),fHistoDEtaMax(0.15),fHistoDEtaMin(-0.15),
fHistoNDPhiBins(100),fHistoDPhiMax(0.15),fHistoDPhiMin(-0.15),
fhTrackMatchedDEta(),fhTrackMatchedDPhi(),fhTrackMatchedDEtaDPhi()

{
	//Ctor
	
	//Initialize parameters
	InitParameters();
}

//________________________________________________
AliCaloPID::AliCaloPID(const TTask * emcalpid) : 
TObject(), fEMCALPhotonWeight(0.), fEMCALPi0Weight(0.),  
fEMCALElectronWeight(0.),  fEMCALChargeWeight(0.),
fEMCALNeutralWeight(0.),
fPHOSPhotonWeight(0.), fPHOSPi0Weight(0.),  
fPHOSElectronWeight(0.), fPHOSChargeWeight(0.) , 
fPHOSNeutralWeight(0.), //fPHOSWeightFormula(0), 
//fPHOSPhotonWeightFormula(0x0), fPHOSPi0WeightFormula(0x0),
fDispCut(0.),fTOFCut(0.), fDebug(-1), 
fRecalculateBayesian(kTRUE), fParticleFlux(-1), 
fEMCALPIDUtils( (AliEMCALPIDUtils*) emcalpid),
fHistoNEBins(100),   fHistoEMax(100.),   fHistoEMin(0.),
fHistoNDEtaBins(100),fHistoDEtaMax(0.15),fHistoDEtaMin(-0.15),
fHistoNDPhiBins(100),fHistoDPhiMax(0.15),fHistoDPhiMin(-0.15),
fhTrackMatchedDEta(),fhTrackMatchedDPhi(),fhTrackMatchedDEtaDPhi()

{
	//Ctor
	
	//Initialize parameters
	InitParameters();
}

//_________________________________
AliCaloPID::~AliCaloPID() {
  //Dtor
  
//  if(fPHOSPhotonWeightFormula) delete fPHOSPhotonWeightFormula ;
//  if(fPHOSPi0WeightFormula)    delete fPHOSPi0WeightFormula ;
  if(fEMCALPIDUtils)           delete fEMCALPIDUtils ;
}


//________________________________________________________________________
TList *  AliCaloPID::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer of the analysis class that calls this class.
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("CaloPIDHistos") ; 
  
  outputContainer->SetOwner(kFALSE);
  
  fhTrackMatchedDEta  = new TH2F
  ("TrackMatchedDEta",
   "d#eta of cluster-track vs cluster energy",fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNDEtaBins,fHistoDEtaMin,fHistoDEtaMax); 
  fhTrackMatchedDEta->SetYTitle("d#eta");
  fhTrackMatchedDEta->SetXTitle("E_{cluster} (GeV)");
  
  fhTrackMatchedDPhi  = new TH2F
  ("TrackMatchedDPhi",
   "d#phi of cluster-track vs cluster energy"
   ,fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNDPhiBins,fHistoDPhiMin,fHistoDPhiMax); 
  fhTrackMatchedDPhi->SetYTitle("d#phi (rad)");
  fhTrackMatchedDPhi->SetXTitle("E_{cluster} (GeV)");
  
  fhTrackMatchedDEtaDPhi  = new TH2F
  ("TrackMatchedDEtaDPhi",
   "d#eta vs d#phi of cluster-track vs cluster energy"
   ,fHistoNDEtaBins,fHistoDEtaMin,fHistoDEtaMax,fHistoNDPhiBins,fHistoDPhiMin,fHistoDPhiMax); 
  fhTrackMatchedDEtaDPhi->SetYTitle("d#phi (rad)");
  fhTrackMatchedDEtaDPhi->SetXTitle("d#eta");   
  
  outputContainer->Add(fhTrackMatchedDEta) ; 
  outputContainer->Add(fhTrackMatchedDPhi) ;
  outputContainer->Add(fhTrackMatchedDEtaDPhi) ; 
  
  return outputContainer;
}



//_______________________________________________________________
void AliCaloPID::InitParameters()
{
  //Initialize the parameters of the PID.
  
  fEMCALPhotonWeight   = 0.6 ;
  fEMCALPi0Weight      = 0.6 ;
  fEMCALElectronWeight = 0.6 ;
  fEMCALChargeWeight   = 0.6 ;
  fEMCALNeutralWeight  = 0.6 ;
  
  fPHOSPhotonWeight    = 0.6 ;
  fPHOSPi0Weight       = 0.6 ;
  fPHOSElectronWeight  = 0.6 ;
  fPHOSChargeWeight    = 0.6 ;
  fPHOSNeutralWeight   = 0.6 ;
  
  //Formula to set the PID weight threshold for photon or pi0
  //fPHOSWeightFormula = kTRUE;
  //if(fPHOSPhotonWeightFormula) delete fPHOSPhotonWeightFormula;
  //if(fPHOSPi0WeightFormula)    delete fPHOSPi0WeightFormula;
  //fPHOSPhotonWeightFormula = 
    //new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
  //fPHOSPi0WeightFormula = 
    //new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
  
  fDispCut  = 0.25;
  fTOFCut   = 5.e-9;
  fDebug = -1;
	
  if(fRecalculateBayesian){
	if(fParticleFlux == kLow){
		printf("AliCaloPID::Init() - SetLOWFluxParam\n");
		fEMCALPIDUtils->SetLowFluxParam() ;
	}
	else if (fParticleFlux == kHigh){
		printf("AliCaloPID::Init() - SetHIGHFluxParam\n");
		fEMCALPIDUtils->SetHighFluxParam() ;
	}
  }
}

//_______________________________________________________________
Int_t AliCaloPID::GetPdg(const TString calo, const Double_t * pid, const Float_t energy) const {
  //Return most probable identity of the particle.
  
  if(!pid){ 
    printf("AliCaloPID::GetPdg() - pid pointer not initialized!!!\n");
    abort();
  }
  
  Float_t wPh  =  fPHOSPhotonWeight ;
  Float_t wPi0 =  fPHOSPi0Weight ;
  Float_t wE   =  fPHOSElectronWeight ;
  Float_t wCh  =  fPHOSChargeWeight ;
  Float_t wNe  =  fPHOSNeutralWeight ;
  
  
//  if(calo == "PHOS" && fPHOSWeightFormula){
//    wPh  = fPHOSPhotonWeightFormula->Eval(energy) ;
//    wPi0 = fPHOSPi0WeightFormula->Eval(energy);
//  }
  
  if(calo == "EMCAL"){
    
    wPh  =  fEMCALPhotonWeight ;
    wPi0 =  fEMCALPi0Weight ;
    wE   =  fEMCALElectronWeight ;
    wCh  =  fEMCALChargeWeight ;
    wNe  =  fEMCALNeutralWeight ;
    
  }
  
  if(fDebug > 0)  printf("AliCaloPID::GetPdg: calo %s, ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, hadrons: pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
			 calo.Data(),pid[AliVCluster::kPhoton], pid[AliVCluster::kPi0],
			 pid[AliVCluster::kElectron], pid[AliVCluster::kEleCon],
			 pid[AliVCluster::kPion], pid[AliVCluster::kKaon], pid[AliVCluster::kProton],
			 pid[AliVCluster::kNeutron], pid[AliVCluster::kKaon0]);
  
  Int_t pdg = kNeutralUnknown ;
  Float_t chargedHadronWeight = pid[AliVCluster::kProton]+pid[AliVCluster::kKaon]+
    pid[AliVCluster::kPion]+pid[AliVCluster::kMuon];
  Float_t neutralHadronWeight = pid[AliVCluster::kNeutron]+pid[AliVCluster::kKaon0];
  Float_t allChargedWeight    = pid[AliVCluster::kElectron]+pid[AliVCluster::kEleCon]+ chargedHadronWeight;
  Float_t allNeutralWeight    = pid[AliVCluster::kPhoton]+pid[AliVCluster::kPi0]+ neutralHadronWeight;
  
  //Select most probable ID
  if(calo=="PHOS"){
    if(pid[AliVCluster::kPhoton] > wPh) pdg = kPhoton ;
    else if(pid[AliVCluster::kPi0] > wPi0) pdg = kPi0 ; 
    else if(pid[AliVCluster::kElectron] > wE)  pdg = kElectron ;
    else if(pid[AliVCluster::kEleCon] >  wE) pdg = kEleCon ;
    else if(chargedHadronWeight > wCh) pdg = kChargedHadron ;  
    else if(neutralHadronWeight > wNe) pdg = kNeutralHadron ; 
    else if(allChargedWeight >  allNeutralWeight)
      pdg = kChargedUnknown ; 
    else 
      pdg = kNeutralUnknown ;
  }
  else{//EMCAL
    
    if(pid[AliVCluster::kPhoton]  > wPh)                     pdg = kPhoton ;
    else if(pid[AliVCluster::kElectron]  > wE)               pdg = kElectron ;
    else if(pid[AliVCluster::kPhoton]+pid[AliVCluster::kElectron]  > wPh) pdg = kPhoton ; //temporal sollution until track matching for electrons is considered
    else if(pid[AliVCluster::kPi0] > wPi0)                   pdg = kPi0 ; 
    else if(chargedHadronWeight + neutralHadronWeight > wCh) pdg = kChargedHadron ;  
    else if(neutralHadronWeight + chargedHadronWeight > wNe) pdg = kNeutralHadron ; 
    else                                                     pdg = kNeutralUnknown ;
  }
  
  if(fDebug > 0)printf("AliCaloPID::GetPdg:Final Pdg: %d, cluster energy %2.2f \n", pdg,energy);
   
  return pdg ;
  
}

//_______________________________________________________________
Int_t AliCaloPID::GetPdg(const TString calo,const TLorentzVector mom, const AliVCluster * cluster) const {
  //Recalculated PID with all parameters
  Float_t lambda0 = cluster->GetM02();
  Float_t energy = mom.E();	

  if(fDebug > 0) printf("AliCaloPID::GetPdg: Calorimeter %s, E %3.2f, l0 %3.2f, l1 %3.2f, disp %3.2f, tof %1.11f, distCPV %3.2f, distToBC %1.1f, NMax %d\n",
						calo.Data(),energy,lambda0,cluster->GetM20(),cluster->GetDispersion(),cluster->GetTOF(), 
						cluster->GetEmcCpvDistance(), cluster->GetDistanceToBadChannel(),cluster->GetNExMax());

  if(calo == "EMCAL") {
	  //Recalculate Bayesian
	  if(fRecalculateBayesian){	  
		  if(fDebug > 0)  {
			  const Double_t  *pid0 = cluster->GetPID();
			  printf("AliCaloPID::GetPdg: BEFORE calo %s, ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, hadrons: pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
								 calo.Data(),pid0[AliVCluster::kPhoton], pid0[AliVCluster::kPi0],
								 pid0[AliVCluster::kElectron], pid0[AliVCluster::kEleCon],
								 pid0[AliVCluster::kPion], pid0[AliVCluster::kKaon], pid0[AliVCluster::kProton],
								 pid0[AliVCluster::kNeutron], pid0[AliVCluster::kKaon0]);
		  }
		  
		 fEMCALPIDUtils->ComputePID(energy, lambda0);
		 Double_t pid[AliPID::kSPECIESN];
		 for(Int_t i = 0; i < AliPID::kSPECIESN; i++) pid[i] = fEMCALPIDUtils->GetPIDFinal(i);
		 return GetPdg(calo, pid, energy);
		  
	  
	}
	  
	  // If no use of bayesian, simple PID  
	  if(lambda0 < 0.25) return kPhoton ;
	  //else return  kNeutralHadron ; 
	  else return  kPi0 ;
  }
  
  //   if(calo == "PHOS") {
  //    if(cluster->GetM02()< 0.25) return kPhoton ;
  //    else return  kNeutralHadron ; 
  //  }
  
  return  kNeutralHadron ; 
  
}

//__________________________________________________
TString  AliCaloPID::GetPIDParametersList()  {
  //Put data member values in string to keep in output container
  
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  snprintf(onePar,buffersize,"--- AliCaloPID ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"fDispCut =%2.2f (Cut on dispersion, used in PID evaluation) \n",fDispCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fTOFCut  =%e (Cut on TOF, used in PID evaluation) \n",fTOFCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fEMCALPhotonWeight =%2.2f (EMCAL bayesian weight for photons)\n",fEMCALPhotonWeight) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fEMCALPi0Weight =%2.2f (EMCAL bayesian weight for pi0)\n",fEMCALPi0Weight) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fEMCALElectronWeight =%2.2f(EMCAL bayesian weight for electrons)\n",fEMCALElectronWeight) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fEMCALChargeWeight =%2.2f (EMCAL bayesian weight for charged hadrons)\n",fEMCALChargeWeight) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fEMCALNeutralWeight =%2.2f (EMCAL bayesian weight for neutral hadrons)\n",fEMCALNeutralWeight) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPHOSPhotonWeight =%2.2f (PHOS bayesian weight for photons)\n",fPHOSPhotonWeight) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPHOSPi0Weight =%2.2f (PHOS bayesian weight for pi0)\n",fPHOSPi0Weight) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPHOSElectronWeight =%2.2f(PHOS bayesian weight for electrons)\n",fPHOSElectronWeight) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPHOSChargeWeight =%2.2f (PHOS bayesian weight for charged hadrons)\n",fPHOSChargeWeight) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPHOSNeutralWeight =%2.2f (PHOS bayesian weight for neutral hadrons)\n",fPHOSNeutralWeight) ;
  parList+=onePar ;
  
//  if(fPHOSWeightFormula){
//	sprintf(onePar,buffersize,"PHOS Photon Weight Formula: %s\n",(fPHOSPhotonWeightFormula->GetExpFormula("p")).Data()) ;
//    parList+=onePar;
//	sprintf(onePar,buffersize,"PHOS Pi0    Weight Formula: %s\n",(fPHOSPi0WeightFormula->GetExpFormula("p")).Data()) ;
//	parList+=onePar;	  
//  }
  
  return parList; 
  
}

//________________________________________________________________
void AliCaloPID::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
  
  printf("PHOS PID weight , photon %0.2f, pi0 %0.2f, e %0.2f, charge %0.2f, neutral %0.2f \n",  
	 fPHOSPhotonWeight,  fPHOSPi0Weight, 
	 fPHOSElectronWeight,  fPHOSChargeWeight,   fPHOSNeutralWeight) ; 
  printf("EMCAL PID weight, photon %0.2f, pi0 %0.2f, e %0.2f, charge %0.2f, neutral %0.2f\n",   
	 fEMCALPhotonWeight,  fEMCALPi0Weight, 
	 fEMCALElectronWeight,  fEMCALChargeWeight,  fEMCALNeutralWeight) ; 
  
//  printf("PHOS Parametrized weight on?  =     %d\n",  fPHOSWeightFormula) ; 
//  if(fPHOSWeightFormula){
//    printf("Photon weight formula = %s\n", (fPHOSPhotonWeightFormula->GetExpFormula("p")).Data());
//    printf("Pi0    weight formula = %s\n", (fPHOSPhotonWeightFormula->GetExpFormula("p")).Data());
//  }
  
  printf("TOF cut        = %e\n",fTOFCut);
  printf("Dispersion cut = %2.2f\n",fDispCut);
  printf("Debug level    = %d\n",fDebug);
  printf("Recalculate Bayesian?    = %d\n",fRecalculateBayesian);
 if(fRecalculateBayesian) printf("Particle Flux?    = %d\n",fParticleFlux);
  printf(" \n");
  
} 

//_______________________________________________________________
void AliCaloPID::SetPIDBits(const TString calo, const AliVCluster * cluster, AliAODPWG4Particle * ph, const AliCalorimeterUtils* cu) {
  //Set Bits for PID selection
  
  //Dispersion/lambdas
  //Double_t disp= cluster->GetDispersion()  ;
  Double_t l1  = cluster->GetM20() ;
  Double_t l0  = cluster->GetM02() ; 
  Bool_t isDispOK = kTRUE ;
  if(cluster->IsPHOS()){
    
    if(l1>= 0   && l0>= 0   && l1 < 0.1 && l0 < 0.1) isDispOK=kFALSE ;
    if(l1>= 0   && l0 > 0.5 && l1 < 0.1 && l0 < 1.5) isDispOK=kTRUE ;
    if(l1>= 0   && l0 > 2.0 && l1 < 0.1 && l0 < 2.7) isDispOK=kFALSE ;
    if(l1>= 0   && l0 > 2.7 && l1 < 0.1 && l0 < 4.0) isDispOK=kFALSE ;
    if(l1 > 0.1 && l1 < 0.7 && l0 > 0.7 && l0 < 2.1) isDispOK=kTRUE ;
    if(l1 > 0.1 && l1 < 0.3 && l0 > 3.0 && l0 < 5.0) isDispOK=kFALSE  ;
    if(l1 > 0.3 && l1 < 0.7 && l0 > 2.5 && l0 < 4.0) isDispOK=kFALSE ;
    if(l1 > 0.7 && l1 < 1.3 && l0 > 1.0 && l0 < 1.6) isDispOK=kTRUE ;
    if(l1 > 0.7 && l1 < 1.3 && l0 > 1.6 && l0 < 3.5) isDispOK=kTRUE ;
    if(l1 > 1.3 && l1 < 3.5 && l0 > 1.3 && l0 < 3.5) isDispOK=kTRUE ;    
    
  }
  else{//EMCAL
    
    if(l0 > fDispCut || l0 < 0.) isDispOK = kFALSE;
    
  }
  
  ph->SetDispBit(isDispOK) ;
  
  //TOF
  Double_t tof=cluster->GetTOF()  ;
  ph->SetTOFBit(TMath::Abs(tof)<fTOFCut) ; 
  
  //Charged veto  
  //Bool_t isNeutral = kTRUE ;
  //if(cluster->IsPHOS())  isNeutral = cluster->GetEmcCpvDistance() > 5. ;
  //else 
  Bool_t isNeutral = IsTrackMatched(cluster,cu);
  
  ph->SetChargedBit(isNeutral);
  
  //Set PID pdg
  ph->SetPdg(GetPdg(calo,cluster->GetPID(),ph->E()));
  
  if(fDebug > 0){ 
    printf("AliCaloPID::SetPIDBits: TOF %e, Lambda0 %2.2f, Lambda1 %2.2f\n",tof , l0, l1); 	
    printf("AliCaloPID::SetPIDBits: pdg %d, bits: TOF %d, Dispersion %d, Charge %d\n",
	     ph->GetPdg(), ph->GetTOFBit() , ph->GetDispBit() , ph->GetChargedBit()); 
  }
}

//__________________________________________________________________________
Bool_t AliCaloPID::IsTrackMatched(const AliVCluster* cluster, const AliCalorimeterUtils * cu) const {
  //Check if there is any track attached to this cluster
  
  Int_t nMatches = cluster->GetNTracksMatched();
//  if(nMatches>0){
//    printf("N matches %d, first match (ESD) %d or (AOD) %d\n",nMatches,cluster->GetTrackMatchedIndex(), cluster->GetTrackMatched(0));
//    if     (cluster->GetTrackMatched(0))        printf("\t matched track id %d\n",((AliVTrack*)cluster->GetTrackMatched(0)) ->GetID() ) ;
//  }
//  else {
//    printf("Not Matched");
//  }

  //If EMCAL track matching needs to be recalculated
  if(cluster->IsEMCAL() && cu && cu->IsRecalculationOfClusterTrackMatchingOn()){
    Float_t dR = 999., dZ = 999.;
    cu->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dR,dZ);

    if(dR < 999) {     
      
      if(fhTrackMatchedDEta){
        fhTrackMatchedDEta->Fill(cluster->E(),dZ);
        fhTrackMatchedDPhi->Fill(cluster->E(),dR);
        if(cluster->E() > 0.5)fhTrackMatchedDEtaDPhi->Fill(dZ,dR);
      }
      //printf("dR %f, dZ %f \n",dR,dZ);
      return kTRUE;
    }
    else         
      return kFALSE;
  }//EMCAL cluster and recalculation of matching on
  
  if(fhTrackMatchedDEta){
    fhTrackMatchedDEta->Fill(cluster->GetTrackDz(),cluster->E());
    fhTrackMatchedDEta->Fill(cluster->GetTrackDx(),cluster->E());
    if(cluster->E() > 0.5)fhTrackMatchedDEtaDPhi->Fill(cluster->GetTrackDz(),cluster->GetTrackDx());
  }
  
  if(!strcmp("AliESDCaloCluster",Form("%s",cluster->ClassName())))
  {    
    if (nMatches > 0) {
      if (nMatches == 1 ) {
        Int_t iESDtrack = cluster->GetTrackMatchedIndex();
        //printf("Track Matched index %d\n",iESDtrack);
        if(iESDtrack==-1) return kFALSE ;// Default value of array, there is no match
        else              return kTRUE;
      }//Just one, check
      else return kTRUE ;//More than one, there is a match.
    }// > 0
    else return kFALSE; //It does not happen, but in case
    
  }//ESDs
  else
  {
//    //Patched AODs
//    if(cu && cu->IsRecalculationOfClusterTrackMatchingOn()){
//      Float_t dR = TMath::Abs(cluster->GetEmcCpvDistance()) ; //FIXME
//      if(dR > cu->GetCutR()) return kFALSE;
//      else                   return kTRUE;
//    }
//    else {
//      if(nMatches > 0) return kTRUE; //There is at least one match.
//      else             return kFALSE;
    
//    }
    
    if(nMatches > 0) return kTRUE; //There is at least one match.
    else             return kFALSE;
    
  }//AODs or MC (copy into AOD)
  
  return kFALSE;
  
}


