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
// Class for track/cluster acceptance selection
// Selection in Central barrel, EMCAL and PHOS
//                
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TMath.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TFormula.h>
  
//---- ANALYSIS system ----
#include "AliLog.h"
#include "AliCaloPID.h"
#include "AliAODCluster.h"

ClassImp(AliCaloPID)


//________________________________________________
AliCaloPID::AliCaloPID() : 
  TObject(), fEMCALPhotonWeight(0.), fEMCALPi0Weight(0.),  
  fEMCALElectronWeight(0.),  fEMCALChargeWeight(0.),
  fEMCALNeutralWeight(0.),
  fPHOSPhotonWeight(0.), fPHOSPi0Weight(0.),  
  fPHOSElectronWeight(0.), fPHOSChargeWeight(0.) , 
  fPHOSNeutralWeight(0.), fPHOSWeightFormula(0), 
  fPHOSPhotonWeightFormula(0x0), fPHOSPi0WeightFormula(0x0) 
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
  fPHOSPi0WeightFormula(pid.fPHOSPi0WeightFormula) 


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

  return *this;

}

//_________________________________
AliCaloPID::~AliCaloPID() {
  //Dtor

  if(fPHOSPhotonWeightFormula) delete  fPHOSPhotonWeightFormula ;
  if(fPHOSPi0WeightFormula) delete  fPHOSPi0WeightFormula ;

}


//_______________________________________________________________
void AliCaloPID::InitParameters()
{
  //Initialize the parameters of the PID.

  fEMCALPhotonWeight = 0.8 ;
  fEMCALPi0Weight = 0.5 ;
  fEMCALElectronWeight = 0.8 ;
  fEMCALChargeWeight = 0.5 ;
  fEMCALNeutralWeight = 0.5 ;

  fPHOSPhotonWeight = 0.75 ;
  fPHOSPi0Weight = 0.8 ;
  fPHOSElectronWeight = 0.5 ;
  fPHOSChargeWeight = 0.5 ;
  fPHOSNeutralWeight = 0.5 ;

  //Formula to set the PID weight threshold for photon or pi0
  fPHOSWeightFormula = kTRUE;
  fPHOSPhotonWeightFormula = 
    new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
  fPHOSPi0WeightFormula = 
    new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
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
//   if(fPHOSWeightFormula){
//     printf(">>>>>>>>>>> Photon weight formula<<<<<<<<<<<<\n");
//     fPHOSPhotonWeightFormula->Print();
//     printf(">>>>>>>>>>> Pi0    weight formula<<<<<<<<<<<<\n");
//     fPHOSPhotonWeightFormula->Print();
//   }
  printf(" \n");

} 

//_______________________________________________________________
Int_t AliCaloPID::GetPdg(const TString calo, const Double_t * pid, const Float_t energy) const {
  //Return most probable identity of the particle.
  
  if(!pid) AliFatal("pid pointer not initialized!!!");

  Float_t wPh =  fPHOSPhotonWeight ;
  Float_t wPi0 =  fPHOSPi0Weight ;
  Float_t wE =  fPHOSElectronWeight ;
  Float_t wCh =  fPHOSChargeWeight ;
  Float_t wNe =  fPHOSNeutralWeight ;
  
  
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
  
//   printf("PID: calo %s, ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, hadrons: pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
// 	 calo.Data(),pid[AliAODCluster::kPhoton], pid[AliAODCluster::kPi0],
// 	 pid[AliAODCluster::kElectron], pid[AliAODCluster::kEleCon],
// 	 pid[AliAODCluster::kPion], pid[AliAODCluster::kKaon], pid[AliAODCluster::kProton],
// 	 pid[AliAODCluster::kNeutron], pid[AliAODCluster::kKaon0]);

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
    //Temporal solution, electrons and photons not differenciated
    if(pid[AliAODCluster::kPhoton] + pid[AliAODCluster::kElectron]  > wPh) pdg = kPhoton ;
    else if(pid[AliAODCluster::kPi0] > wPi0) pdg = kPi0 ; 
    else if(chargedHadronWeight + neutralHadronWeight > wCh) pdg = kChargedHadron ;  
    else if(neutralHadronWeight + chargedHadronWeight > wNe) pdg = kNeutralHadron ; 
    else pdg =  kNeutralUnknown ;

  }


  //printf("Final Pdg: %d \n", pdg);
  


  return pdg ;
  
}

//_______________________________________________________________
Int_t AliCaloPID::GetPdg(const TString calo, const TLorentzVector mom, const Double_t l0, 
			 const Double_t l1, const Double_t disp, const Double_t tof, 
			 const Double_t distCPV) const {
  //Recalculated PID with all parameters
  AliDebug(2,Form("Calorimeter %s, E %3.2f, l0 %3.2f, l1 %3.2f, disp %3.2f, tof %3.2f, distCPV %3.2f",
		  calo.Data(),mom.E(),l0,l1,disp,tof,distCPV));

  if(calo == "EMCAL") {
    if(l0 < 0.25) return kPhoton ;
    else return  kNeutralHadron ; 
  }
  
  return  kNeutralHadron ; 
  
}

