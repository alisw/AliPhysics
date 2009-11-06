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
//   - GetPdg(const TString calo,const TLorentzVector mom, const AliAODCaloCluster * cluster)
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
#include <TFormula.h>

//---- ANALYSIS system ----
#include "AliCaloPID.h"
#include "AliAODCaloCluster.h"
#include "AliAODPWG4Particle.h"
#include "AliEMCALPIDUtils.h"
ClassImp(AliCaloPID)


//________________________________________________
AliCaloPID::AliCaloPID() : 
TObject(), fEMCALPhotonWeight(0.), fEMCALPi0Weight(0.),  
fEMCALElectronWeight(0.),  fEMCALChargeWeight(0.),
fEMCALNeutralWeight(0.),
fPHOSPhotonWeight(0.), fPHOSPi0Weight(0.),  
fPHOSElectronWeight(0.), fPHOSChargeWeight(0.) , 
fPHOSNeutralWeight(0.), fPHOSWeightFormula(0), 
fPHOSPhotonWeightFormula(0x0), fPHOSPi0WeightFormula(0x0),
fDispCut(0.),fTOFCut(0.), fDebug(-1), 
fRecalculateBayesian(kFALSE), fParticleFlux(kLow), fEMCALPIDUtils(new AliEMCALPIDUtils)
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
fPHOSNeutralWeight(0.), fPHOSWeightFormula(0), 
fPHOSPhotonWeightFormula(0x0), fPHOSPi0WeightFormula(0x0),
fDispCut(0.),fTOFCut(0.), fDebug(-1), 
fRecalculateBayesian(kTRUE), fParticleFlux(flux), fEMCALPIDUtils(new AliEMCALPIDUtils)
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
fPHOSNeutralWeight(0.), fPHOSWeightFormula(0), 
fPHOSPhotonWeightFormula(0x0), fPHOSPi0WeightFormula(0x0),
fDispCut(0.),fTOFCut(0.), fDebug(-1), 
fRecalculateBayesian(kTRUE), fParticleFlux(-1), fEMCALPIDUtils( (AliEMCALPIDUtils*) emcalpid)
{
	//Ctor
	
	//Initialize parameters
	InitParameters();
}

//____________________________________________________________________________
AliCaloPID::AliCaloPID(const AliCaloPID & pid) :   
TObject(pid), fEMCALPhotonWeight(pid.fEMCALPhotonWeight), 
fEMCALPi0Weight(pid.fEMCALPi0Weight), 
fEMCALElectronWeight(pid.fEMCALElectronWeight), 
fEMCALChargeWeight(pid.fEMCALChargeWeight), 
fEMCALNeutralWeight(pid.fEMCALNeutralWeight), 
fPHOSPhotonWeight(pid.fPHOSPhotonWeight),
fPHOSPi0Weight(pid.fPHOSPi0Weight),
fPHOSElectronWeight(pid.fPHOSElectronWeight), 
fPHOSChargeWeight(pid.fPHOSChargeWeight),
fPHOSNeutralWeight(pid.fPHOSNeutralWeight),
fPHOSWeightFormula(pid.fPHOSWeightFormula), 
fPHOSPhotonWeightFormula(pid.fPHOSPhotonWeightFormula), 
fPHOSPi0WeightFormula(pid.fPHOSPi0WeightFormula), 
fDispCut(pid.fDispCut),fTOFCut(pid.fTOFCut),
fDebug(pid.fDebug), fRecalculateBayesian(pid.fRecalculateBayesian),
fParticleFlux(pid.fParticleFlux), fEMCALPIDUtils(pid.fEMCALPIDUtils)
{
  // cpy ctor
	
}

//_________________________________________________________________________
AliCaloPID & AliCaloPID::operator = (const AliCaloPID & pid)
{
  // assignment operator
  
  if(&pid == this) return *this;
  
  fEMCALPhotonWeight = pid. fEMCALPhotonWeight ;
  fEMCALPi0Weight = pid.fEMCALPi0Weight ;
  fEMCALElectronWeight = pid.fEMCALElectronWeight; 
  fEMCALChargeWeight = pid.fEMCALChargeWeight;
  fEMCALNeutralWeight = pid.fEMCALNeutralWeight;
  
  fPHOSPhotonWeight = pid.fPHOSPhotonWeight ;
  fPHOSPi0Weight = pid.fPHOSPi0Weight ;
  fPHOSElectronWeight = pid.fPHOSElectronWeight; 
  fPHOSChargeWeight = pid.fPHOSChargeWeight;
  fPHOSNeutralWeight = pid.fPHOSNeutralWeight;
  
  fPHOSWeightFormula       = pid.fPHOSWeightFormula; 
  fPHOSPhotonWeightFormula = pid.fPHOSPhotonWeightFormula; 
  fPHOSPi0WeightFormula    = pid.fPHOSPi0WeightFormula;
  
  fDispCut  = pid.fDispCut;
  fTOFCut   = pid.fTOFCut;
  fDebug    = pid.fDebug;
  
  fRecalculateBayesian = pid.fRecalculateBayesian;
  fParticleFlux        = pid.fParticleFlux;
  fEMCALPIDUtils       = pid.fEMCALPIDUtils;
		
  return *this;
  
}

//_________________________________
AliCaloPID::~AliCaloPID() {
  //Dtor
  
  if(fPHOSPhotonWeightFormula) delete  fPHOSPhotonWeightFormula ;
  if(fPHOSPi0WeightFormula) delete  fPHOSPi0WeightFormula ;
  if(fEMCALPIDUtils) delete fEMCALPIDUtils ;
}


//_______________________________________________________________
void AliCaloPID::InitParameters()
{
  //Initialize the parameters of the PID.
  
  fEMCALPhotonWeight   = 0.5 ;
  fEMCALPi0Weight      = 0.5 ;
  fEMCALElectronWeight = 0.5 ;
  fEMCALChargeWeight   = 0.5 ;
  fEMCALNeutralWeight  = 0.5 ;
  
  fPHOSPhotonWeight    = 0.5 ;
  fPHOSPi0Weight       = 0.5 ;
  fPHOSElectronWeight  = 0.5 ;
  fPHOSChargeWeight    = 0.5 ;
  fPHOSNeutralWeight   = 0.5 ;
  
  //Formula to set the PID weight threshold for photon or pi0
  fPHOSWeightFormula = kTRUE;
  fPHOSPhotonWeightFormula = 
    new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
  fPHOSPi0WeightFormula = 
    new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
  
  fDispCut  = 1.5;
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
  
  
  if(calo == "PHOS" && fPHOSWeightFormula){
    wPh  = fPHOSPhotonWeightFormula->Eval(energy) ;
    wPi0 = fPHOSPi0WeightFormula->Eval(energy);
  }
  
  if(calo == "EMCAL"){
    
    wPh  =  fEMCALPhotonWeight ;
    wPi0 =  fEMCALPi0Weight ;
    wE   =  fEMCALElectronWeight ;
    wCh  =  fEMCALChargeWeight ;
    wNe  =  fEMCALNeutralWeight ;
    
  }
  
  if(fDebug > 0)  printf("AliCaloPID::GetPdg: calo %s, ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, hadrons: pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
			 calo.Data(),pid[AliAODCluster::kPhoton], pid[AliAODCluster::kPi0],
			 pid[AliAODCluster::kElectron], pid[AliAODCluster::kEleCon],
			 pid[AliAODCluster::kPion], pid[AliAODCluster::kKaon], pid[AliAODCluster::kProton],
			 pid[AliAODCluster::kNeutron], pid[AliAODCluster::kKaon0]);
  
  Int_t pdg = kNeutralUnknown ;
  Float_t chargedHadronWeight = pid[AliAODCluster::kProton]+pid[AliAODCluster::kKaon]+
    pid[AliAODCluster::kPion]+pid[AliAODCluster::kMuon];
  Float_t neutralHadronWeight = pid[AliAODCluster::kNeutron]+pid[AliAODCluster::kKaon0];
  Float_t allChargedWeight    = pid[AliAODCluster::kElectron]+pid[AliAODCluster::kEleCon]+ chargedHadronWeight;
  Float_t allNeutralWeight    = pid[AliAODCluster::kPhoton]+pid[AliAODCluster::kPi0]+ neutralHadronWeight;
  
  //Select most probable ID
  if(calo=="PHOS"){
    if(pid[AliAODCluster::kPhoton] > wPh) pdg = kPhoton ;
    else if(pid[AliAODCluster::kPi0] > wPi0) pdg = kPi0 ; 
    else if(pid[AliAODCluster::kElectron] > wE)  pdg = kElectron ;
    else if(pid[AliAODCluster::kEleCon] >  wE) pdg = kEleCon ;
    else if(chargedHadronWeight > wCh) pdg = kChargedHadron ;  
    else if(neutralHadronWeight > wNe) pdg = kNeutralHadron ; 
    else if(allChargedWeight >  allNeutralWeight)
      pdg = kChargedUnknown ; 
    else 
      pdg = kNeutralUnknown ;
  }
  else{//EMCAL

    if(pid[AliAODCluster::kPhoton]+pid[AliAODCluster::kElectron]  > wPh) pdg = kPhoton ; //temporal sollution until track matching for electrons is considered
    //if(pid[AliAODCluster::kPhoton]  > wPh) pdg = kPhoton ;
    else if(pid[AliAODCluster::kPi0] > wPi0) pdg = kPi0 ; 
    //else if(pid[AliAODCluster::kElectron]  > wE) pdg = kElectron ;
    else if(chargedHadronWeight + neutralHadronWeight > wCh) pdg = kChargedHadron ;  
    else if(neutralHadronWeight + chargedHadronWeight > wNe) pdg = kNeutralHadron ; 
    else pdg =  kNeutralUnknown ;
    
  }
  
  if(fDebug > 0)printf("AliCaloPID::GetPdg:Final Pdg: %d \n", pdg);
   
  return pdg ;
  
}

//_______________________________________________________________
Int_t AliCaloPID::GetPdg(const TString calo,const TLorentzVector mom, const AliAODCaloCluster * cluster) const {
  //Recalculated PID with all parameters
  Float_t lambda0 = cluster->GetM02();
  Float_t energy = mom.E();	

  if(fDebug > 0) printf("AliCaloPID::GetPdg: Calorimeter %s, E %3.2f, l0 %3.2f, l1 %3.2f, disp %3.2f, tof %1.11f, distCPV %3.2f, distToBC %1.1f, NMax %d\n",
						calo.Data(),energy,lambda0,cluster->GetM20(),cluster->GetDispersion(),cluster->GetTOF(), 
						cluster->GetEmcCpvDistance(), cluster->GetDistToBadChannel(),cluster->GetNExMax());

  if(calo == "EMCAL") {
	  //Recalculate Bayesian
	  if(fRecalculateBayesian){	  
		  if(fDebug > 0)  {
			  const Double_t  *pid0 = cluster->PID();
			  printf("AliCaloPID::GetPdg: BEFORE calo %s, ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, hadrons: pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
								 calo.Data(),pid0[AliAODCluster::kPhoton], pid0[AliAODCluster::kPi0],
								 pid0[AliAODCluster::kElectron], pid0[AliAODCluster::kEleCon],
								 pid0[AliAODCluster::kPion], pid0[AliAODCluster::kKaon], pid0[AliAODCluster::kProton],
								 pid0[AliAODCluster::kNeutron], pid0[AliAODCluster::kKaon0]);
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
  char onePar[255] ;
  sprintf(onePar,"--- AliCaloPID ---\n") ;
  parList+=onePar ;	
  sprintf(onePar,"fDispCut =%2.2f (Cut on dispersion, used in PID evaluation) \n",fDispCut) ;
  parList+=onePar ;
  sprintf(onePar,"fTOFCut  =%e (Cut on TOF, used in PID evaluation) \n",fTOFCut) ;
  parList+=onePar ;
  sprintf(onePar,"fEMCALPhotonWeight =%2.2f (EMCAL bayesian weight for photons)\n",fEMCALPhotonWeight) ;
  parList+=onePar ;
  sprintf(onePar,"fEMCALPi0Weight =%2.2f (EMCAL bayesian weight for pi0)\n",fEMCALPi0Weight) ;
  parList+=onePar ;
  sprintf(onePar,"fEMCALElectronWeight =%2.2f(EMCAL bayesian weight for electrons)\n",fEMCALElectronWeight) ;
  parList+=onePar ;
  sprintf(onePar,"fEMCALChargeWeight =%2.2f (EMCAL bayesian weight for charged hadrons)\n",fEMCALChargeWeight) ;
  parList+=onePar ;
  sprintf(onePar,"fEMCALNeutralWeight =%2.2f (EMCAL bayesian weight for neutral hadrons)\n",fEMCALNeutralWeight) ;
  parList+=onePar ;
  sprintf(onePar,"fPHOSPhotonWeight =%2.2f (PHOS bayesian weight for photons)\n",fPHOSPhotonWeight) ;
  parList+=onePar ;
  sprintf(onePar,"fPHOSPi0Weight =%2.2f (PHOS bayesian weight for pi0)\n",fPHOSPi0Weight) ;
  parList+=onePar ;
  sprintf(onePar,"fPHOSElectronWeight =%2.2f(PHOS bayesian weight for electrons)\n",fPHOSElectronWeight) ;
  parList+=onePar ;
  sprintf(onePar,"fPHOSChargeWeight =%2.2f (PHOS bayesian weight for charged hadrons)\n",fPHOSChargeWeight) ;
  parList+=onePar ;
  sprintf(onePar,"fPHOSNeutralWeight =%2.2f (PHOS bayesian weight for neutral hadrons)\n",fPHOSNeutralWeight) ;
  parList+=onePar ;
  
  if(fPHOSWeightFormula){
	sprintf(onePar,"PHOS Photon Weight Formula: %s\n",(fPHOSPhotonWeightFormula->GetExpFormula("p")).Data()) ;
    parList+=onePar;
	sprintf(onePar,"PHOS Pi0    Weight Formula: %s\n",(fPHOSPi0WeightFormula->GetExpFormula("p")).Data()) ;
	parList+=onePar;	  
  }
  
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
  
  printf("PHOS Parametrized weight on?  =     %d\n",  fPHOSWeightFormula) ; 
  if(fPHOSWeightFormula){
    printf("Photon weight formula = %s\n", (fPHOSPhotonWeightFormula->GetExpFormula("p")).Data());
    printf("Pi0    weight formula = %s\n", (fPHOSPhotonWeightFormula->GetExpFormula("p")).Data());
  }
  
  printf("TOF cut        = %e\n",fTOFCut);
  printf("Dispersion cut = %2.2f\n",fDispCut);
  printf("Debug level    = %d\n",fDebug);
  printf("Recalculate Bayesian?    = %d\n",fRecalculateBayesian);
 if(fRecalculateBayesian) printf("Particle Flux?    = %d\n",fParticleFlux);
  printf(" \n");
  
} 

//_______________________________________________________________
void AliCaloPID::SetPIDBits(const TString calo, const AliAODCaloCluster * cluster, AliAODPWG4Particle * ph) {
  //Set Bits for PID selection
  
  //Dispersion/lambdas
  Double_t disp=cluster->GetDispersion()  ;
	//    Double_t m20=calo->GetM20() ;
	//    Double_t m02=calo->GetM02() ; 
  ph->SetDispBit(disp<fDispCut) ;  
  
  //TOF
  Double_t tof=cluster->GetTOF()  ;
  ph->SetTOFBit(TMath::Abs(tof)<fTOFCut) ; 
  
  //Charged veto
  //    Double_t cpvR=calo->GetEmcCpvDistance() ; 
  Int_t ntr=cluster->GetNTracksMatched();  //number of track matched
  ph->SetChargedBit(ntr>0) ;  //Temporary cut, should we evaluate distance?
  
  //Set PID pdg
  ph->SetPdg(GetPdg(calo,cluster->PID(),ph->E()));
  
  if(fDebug > 0){ 
    printf("AliCaloPID::SetPIDBits: TOF %e, Dispersion %2.2f, NTracks %d\n",tof , disp, ntr); 	
    printf("AliCaloPID::SetPIDBits: pdg %d, bits: TOF %d, Dispersion %d, Charge %d\n",
	     ph->GetPdg(), ph->GetTOFBit() , ph->GetDispBit() , ph->GetChargedBit()); 
  }
}


