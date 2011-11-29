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
// The Output of the main method GetIdentifiedParticleType is a PDG number identifying the cluster, 
// being kPhoton, kElectron, kPi0 ... as defined in the header file
//   - GetIdentifiedParticleType(const TString calo, const TLorentzVector mom, const AliVCluster * cluster) 
//      Assignes a PID tag to the cluster, right now there is the possibility to : use bayesian weights from reco, 
//      recalculate them (EMCAL) or use other procedures not used in reco.
//      In order to recalculate Bayesian, it is necessary to load the EMCALUtils library
//      and do SwitchOnBayesianRecalculation().
//      To change the PID parameters from Low to High like the ones by default, use the constructor 
//      AliCaloPID(flux)
//      where flux is AliCaloPID::kLow or AliCaloPID::kHigh
//      If it is necessary to change the parameters use the constructor 
//      AliCaloPID(AliEMCALPIDUtils *utils) and set the parameters before.

//   - GetGetIdentifiedParticleTypeFromBayesian(const TString calo, const Double_t * pid, const Float_t energy)
//      Reads the PID weights array of the ESDs and depending on its magnitude identifies the particle, 
//      executed when bayesian is ON by GetIdentifiedParticleType(const TString calo, const TLorentzVector mom, const AliVCluster * cluster) 
//   - SetPIDBits: Simple PID, depending on the thresholds fLOCut fTOFCut and even the
//     result of the PID bayesian a different PID bit is set. 
//
//   - IsTrackMatched(): Independent method that needs to be combined with GetIdentifiedParticleType to know if the cluster was matched
//
//*-- Author: Gustavo Conesa (INFN-LNF)
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TMath.h>
#include <TString.h>
#include <TList.h>

// ---- ANALYSIS system ----
#include "AliCaloPID.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliAODPWG4Particle.h"
#include "AliCalorimeterUtils.h"
#include "AliVEvent.h"

// ---- Detector ----
#include "AliEMCALPIDUtils.h"

ClassImp(AliCaloPID)


//________________________
AliCaloPID::AliCaloPID() : 
TObject(),                fDebug(-1),                  fParticleFlux(kLow),
//Bayesian
fEMCALPIDUtils(),         fUseBayesianWeights(kFALSE), fRecalculateBayesian(kFALSE),
fEMCALPhotonWeight(0.),   fEMCALPi0Weight(0.),  
fEMCALElectronWeight(0.), fEMCALChargeWeight(0.),      fEMCALNeutralWeight(0.),
fPHOSPhotonWeight(0.),    fPHOSPi0Weight(0.),  
fPHOSElectronWeight(0.),  fPHOSChargeWeight(0.) ,      fPHOSNeutralWeight(0.), 
fPHOSWeightFormula(0),    fPHOSPhotonWeightFormula(0), fPHOSPi0WeightFormula(0),
fPHOSPhotonWeightFormulaExpression(""), 
fPHOSPi0WeightFormulaExpression(""),
//PID calculation
fEMCALL0CutMax(100.),     fEMCALL0CutMin(0),           
fEMCALDEtaCut(2000.),     fEMCALDPhiCut(2000.),
fTOFCut(0.), 
fPHOSDispersionCut(1000), fPHOSRCut(1000),                 
// Histogram settings
fHistoNEBins(100),        fHistoEMax(100.),            fHistoEMin(0.),
fHistoNDEtaBins(100),     fHistoDEtaMax(0.15),         fHistoDEtaMin(-0.15),
fHistoNDPhiBins(100),     fHistoDPhiMax(0.15),         fHistoDPhiMin(-0.15),
fhTrackMatchedDEta(0x0),  fhTrackMatchedDPhi(0x0),     fhTrackMatchedDEtaDPhi(0x0)
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
}

//________________________________________
AliCaloPID::AliCaloPID(const Int_t flux) : 
TObject(),                fDebug(-1),                  fParticleFlux(flux),
//Bayesian
fEMCALPIDUtils(),         fUseBayesianWeights(kFALSE), fRecalculateBayesian(kFALSE),
fEMCALPhotonWeight(0.),   fEMCALPi0Weight(0.),  
fEMCALElectronWeight(0.), fEMCALChargeWeight(0.),      fEMCALNeutralWeight(0.),
fPHOSPhotonWeight(0.),    fPHOSPi0Weight(0.),  
fPHOSElectronWeight(0.),  fPHOSChargeWeight(0.) ,      fPHOSNeutralWeight(0.), 
fPHOSWeightFormula(0),    fPHOSPhotonWeightFormula(0), fPHOSPi0WeightFormula(0),
fPHOSPhotonWeightFormulaExpression(""), 
fPHOSPi0WeightFormulaExpression(""),
//PID calculation
fEMCALL0CutMax(100.),     fEMCALL0CutMin(0),           
fEMCALDEtaCut(2000.),     fEMCALDPhiCut(2000.),
fTOFCut(0.), 
fPHOSDispersionCut(1000), fPHOSRCut(1000),                 
// Histogram settings
fHistoNEBins(100),        fHistoEMax(100.),            fHistoEMin(0.),
fHistoNDEtaBins(100),     fHistoDEtaMax(0.15),         fHistoDEtaMin(-0.15),
fHistoNDPhiBins(100),     fHistoDPhiMax(0.15),         fHistoDPhiMin(-0.15),
fhTrackMatchedDEta(0x0),  fhTrackMatchedDPhi(0x0),     fhTrackMatchedDEtaDPhi(0x0)
{
  //Ctor
	
  //Initialize parameters
  InitParameters();
  
}

//_______________________________________________
AliCaloPID::AliCaloPID(const TNamed * emcalpid) : 
TObject(),                   fDebug(-1),                  fParticleFlux(kLow),
//Bayesian
fEMCALPIDUtils((AliEMCALPIDUtils*)emcalpid),         
fUseBayesianWeights(kFALSE), fRecalculateBayesian(kFALSE),
fEMCALPhotonWeight(0.),      fEMCALPi0Weight(0.),  
fEMCALElectronWeight(0.),    fEMCALChargeWeight(0.),      fEMCALNeutralWeight(0.),
fPHOSPhotonWeight(0.),       fPHOSPi0Weight(0.),  
fPHOSElectronWeight(0.),     fPHOSChargeWeight(0.) ,      fPHOSNeutralWeight(0.), 
fPHOSWeightFormula(0),       fPHOSPhotonWeightFormula(0), fPHOSPi0WeightFormula(0),
fPHOSPhotonWeightFormulaExpression(""), 
fPHOSPi0WeightFormulaExpression(""),
//PID calculation
fEMCALL0CutMax(100.),        fEMCALL0CutMin(0),   
fEMCALDEtaCut(2000.),        fEMCALDPhiCut(2000.),
fTOFCut(0.), 
fPHOSDispersionCut(1000),    fPHOSRCut(1000),                 
// Histogram settings
fHistoNEBins(100),           fHistoEMax(100.),            fHistoEMin(0.),
fHistoNDEtaBins(100),        fHistoDEtaMax(0.15),         fHistoDEtaMin(-0.15),
fHistoNDPhiBins(100),        fHistoDPhiMax(0.15),         fHistoDPhiMin(-0.15),
fhTrackMatchedDEta(0x0),     fhTrackMatchedDPhi(0x0),     fhTrackMatchedDEtaDPhi(0x0)
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
}

//_______________________
AliCaloPID::~AliCaloPID() 
{
  //Dtor
  
  delete fPHOSPhotonWeightFormula ;
  delete fPHOSPi0WeightFormula ;
  delete fEMCALPIDUtils ;
  
}

//___________________________________________
TList *  AliCaloPID::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer of the analysis class that calls this class.
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("CaloPIDHistos") ; 
  
  outputContainer->SetOwner(kFALSE);
  
  fhTrackMatchedDEta  = new TH2F
  ("TrackMatchedDEta",
   "d#eta of cluster-track vs cluster energy",
   fHistoNEBins,fHistoEMin,fHistoEMax,fHistoNDEtaBins,fHistoDEtaMin,fHistoDEtaMax); 
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



//_______________________________
void AliCaloPID::InitParameters()
{
  //Initialize the parameters of the PID.
  
  // Bayesian
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
  fPHOSWeightFormula       = kFALSE;
  fPHOSPhotonWeightFormulaExpression = "0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))";
  fPHOSPi0WeightFormulaExpression    = "0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))"   ;
  
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
  
  //PID recalculation, not bayesian
  
  //EMCAL
  fEMCALL0CutMax = 0.3 ;
  fEMCALL0CutMin = 0.01;
  
  fEMCALDPhiCut  = 0.05; // Same cut as in AliEMCALRecoUtils
  fEMCALDEtaCut  = 0.025;// Same cut as in AliEMCALRecoUtils

  // PHOS / EMCAL, not used
  fTOFCut        = 1.e-6;
  
  //PHOS
  fPHOSRCut          = 2. ; 
  fPHOSDispersionCut = 2.5;
  
}

//______________________________________________
AliEMCALPIDUtils *AliCaloPID::GetEMCALPIDUtils() 
{
  // return pointer to AliEMCALPIDUtils, create it if needed
  
  if(!fEMCALPIDUtils) fEMCALPIDUtils = new AliEMCALPIDUtils ; 
  return fEMCALPIDUtils ; 
  
}


//______________________________________________________________________
Int_t AliCaloPID::GetIdentifiedParticleType(const TString calo,
                                            const TLorentzVector mom, 
                                            const AliVCluster * cluster) 
{
  // Returns a PDG number corresponding to the likely ID of the cluster
  
  Float_t energy  = mom.E();	
  Float_t lambda0 = cluster->GetM02();
  Float_t lambda1 = cluster->GetM20();
  
  // ---------------------
  // Use bayesian approach
  // ---------------------
  
  if(fUseBayesianWeights){
    
    Double_t weights[AliPID::kSPECIESN];
    
    if(calo == "EMCAL"&& fRecalculateBayesian){	        
      fEMCALPIDUtils->ComputePID(energy, lambda0);
      for(Int_t i = 0; i < AliPID::kSPECIESN; i++) weights[i] = fEMCALPIDUtils->GetPIDFinal(i);
    }
    else {
      for(Int_t i = 0; i < AliPID::kSPECIESN; i++) weights[i] = cluster->GetPID()[i];
    }

    if(fDebug > 0)  {
      printf("AliCaloPID::GetIdentifiedParticleType: BEFORE calo %s, ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, hadrons: pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
             calo.Data(),
             weights[AliVCluster::kPhoton],    weights[AliVCluster::kPi0],
             weights[AliVCluster::kElectron],  weights[AliVCluster::kEleCon],
             weights[AliVCluster::kPion],      weights[AliVCluster::kKaon], 
             weights[AliVCluster::kProton],
             weights[AliVCluster::kNeutron],   weights[AliVCluster::kKaon0]);
    }
    
    return GetIdentifiedParticleTypeFromBayesWeights(calo, weights, energy);
  }
  
  // -------------------------------------------------------
  // Calculate PID SS from data, do not use bayesian weights
  // -------------------------------------------------------
  
  if(fDebug > 0) printf("AliCaloPID::GetIdentifiedParticleType: Calorimeter %s, E %3.2f, l0 %3.2f, l1 %3.2f, disp %3.2f, tof %1.11f, distCPV %3.2f, distToBC %1.1f, NMax %d\n",
                        calo.Data(),energy,lambda0,cluster->GetM20(),cluster->GetDispersion(),cluster->GetTOF(), 
                        cluster->GetEmcCpvDistance(), cluster->GetDistanceToBadChannel(),cluster->GetNExMax());
  
  if(cluster->IsEMCAL()){
    
    if(fDebug > 0) printf("AliCaloPID::GetIdentifiedParticleType() - EMCAL SS %f <%f < %f?\n",fEMCALL0CutMin, lambda0, fEMCALL0CutMax);
    
    if(lambda0 < fEMCALL0CutMax && lambda0 > fEMCALL0CutMin) return kPhoton ;
    else                                                     return kNeutralUnknown ; 
  }//EMCAL
  else {//PHOS
    if(TestPHOSDispersion(mom.Pt(),lambda0,lambda1) < fPHOSDispersionCut) return kPhoton;
    else                                                                  return kNeutralUnknown;
  }
  
}

//_______________________________________________________________________________
Int_t AliCaloPID::GetIdentifiedParticleTypeFromBayesWeights(const TString calo, 
                                                            const Double_t * pid, 
                                                            const Float_t energy) 
{
  //Return most probable identity of the particle after bayesian weights calculated in reconstruction
  
  if(!pid){ 
    printf("AliCaloPID::GetIdentifiedParticleType() - pid pointer not initialized!!!\n");
    abort();
  }
  
  Float_t wPh  =  fPHOSPhotonWeight ;
  Float_t wPi0 =  fPHOSPi0Weight ;
  Float_t wE   =  fPHOSElectronWeight ;
  Float_t wCh  =  fPHOSChargeWeight ;
  Float_t wNe  =  fPHOSNeutralWeight ;
  
  if(calo == "PHOS" && fPHOSWeightFormula){
    wPh  = GetPHOSPhotonWeightFormula()->Eval(energy) ;
    wPi0 = GetPHOSPi0WeightFormula()   ->Eval(energy);
  }
  
  if(calo == "EMCAL"){
    
    wPh  =  fEMCALPhotonWeight ;
    wPi0 =  fEMCALPi0Weight ;
    wE   =  fEMCALElectronWeight ;
    wCh  =  fEMCALChargeWeight ;
    wNe  =  fEMCALNeutralWeight ;
    
  }
  
  if(fDebug > 0)  printf("AliCaloPID::GetIdentifiedParticleType: calo %s, ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, hadrons: pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
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
    if(pid[AliVCluster::kPhoton] > wPh)        pdg = kPhoton ;
    else if(pid[AliVCluster::kPi0] > wPi0)     pdg = kPi0 ; 
    else if(pid[AliVCluster::kElectron] > wE)  pdg = kElectron ;
    else if(pid[AliVCluster::kEleCon] >  wE)   pdg = kEleCon ;
    else if(chargedHadronWeight > wCh)         pdg = kChargedHadron ;  
    else if(neutralHadronWeight > wNe)         pdg = kNeutralHadron ; 
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
  
  if(fDebug > 0)printf("AliCaloPID::GetIdentifiedParticleType:Final Pdg: %d, cluster energy %2.2f \n", pdg,energy);

  return pdg ;
  
}

//_________________________________________
TString  AliCaloPID::GetPIDParametersList()  
{
  //Put data member values in string to keep in output container
  
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  snprintf(onePar,buffersize,"--- AliCaloPID ---\n") ;
  parList+=onePar ;	
  if(fUseBayesianWeights){
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
    
    if(fPHOSWeightFormula){
      snprintf(onePar,buffersize,"PHOS Photon Weight Formula: %s\n",fPHOSPhotonWeightFormulaExpression.Data() ) ;
      parList+=onePar;
      snprintf(onePar,buffersize,"PHOS Pi0    Weight Formula: %s\n",fPHOSPi0WeightFormulaExpression.Data()    ) ;
      parList+=onePar;	  
    }
  }
  else {
    snprintf(onePar,buffersize,"EMCAL: fEMCALL0CutMin =%2.2f, fEMCALL0CutMax =%2.2f  (Cut on Shower Shape) \n",fEMCALL0CutMin, fEMCALL0CutMax) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"EMCAL: fEMCALDEtaCut =%2.2f, fEMCALDPhiCut =%2.2f  (Cut on track matching) \n",fEMCALDEtaCut, fEMCALDPhiCut) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fTOFCut  =%e (Cut on TOF, used in PID evaluation) \n",fTOFCut) ;
    parList+=onePar ;	
    snprintf(onePar,buffersize,"fPHOSRCut =%2.2f, fPHOSDispersionCut =%2.2f  (Cut on Shower Shape and CPV) \n",fPHOSRCut,fPHOSDispersionCut) ;
    parList+=onePar ;
    
  }
  
  return parList; 
  
}

//________________________________________________
void AliCaloPID::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
  
  if(fUseBayesianWeights){
    printf("PHOS PID weight , photon %0.2f, pi0 %0.2f, e %0.2f, charge %0.2f, neutral %0.2f \n",  
           fPHOSPhotonWeight,  fPHOSPi0Weight, 
           fPHOSElectronWeight,  fPHOSChargeWeight,   fPHOSNeutralWeight) ; 
    printf("EMCAL PID weight, photon %0.2f, pi0 %0.2f, e %0.2f, charge %0.2f, neutral %0.2f\n",   
           fEMCALPhotonWeight,  fEMCALPi0Weight, 
           fEMCALElectronWeight,  fEMCALChargeWeight,  fEMCALNeutralWeight) ; 
    
    printf("PHOS Parametrized weight on?  =     %d\n",  fPHOSWeightFormula) ; 
    if(fPHOSWeightFormula){
      printf("Photon weight formula = %s\n", fPHOSPhotonWeightFormulaExpression.Data());
      printf("Pi0    weight formula = %s\n", fPHOSPi0WeightFormulaExpression   .Data());
    }
    if(fRecalculateBayesian) printf(" Recalculate bayesian with Particle Flux?    = %d\n",fParticleFlux);
  }
  else {
    printf("TOF cut        = %e\n",fTOFCut);
    printf("EMCAL Lambda0 cut min = %2.2f; max = %2.2f\n",fEMCALL0CutMin, fEMCALL0CutMax);
    printf("EMCAL cluster-track dEta < %2.3f; dPhi < %2.3f\n",fEMCALDEtaCut, fEMCALDPhiCut);
    printf("PHOS Treac matching cut =%2.2f, Dispersion Cut =%2.2f \n",fPHOSRCut,fPHOSDispersionCut) ;
    
  }
  
  printf(" \n");
  
} 

//___________________________________________________________________________
void AliCaloPID::SetPIDBits(const TString calo, AliVCluster * cluster, 
                            AliAODPWG4Particle * ph, AliCalorimeterUtils* cu, 
                            AliVEvent* event) 
{
  //Set Bits for PID selection
  
  //Dispersion/lambdas
  //Double_t disp= cluster->GetDispersion()  ;
  Double_t l1  = cluster->GetM20() ;
  Double_t l0  = cluster->GetM02() ; 
  Bool_t isDispOK = kTRUE ;
  if(cluster->IsPHOS()){ 
    if(TestPHOSDispersion(ph->Pt(),l0,l1) < fPHOSDispersionCut) isDispOK = kTRUE;
    else                                                        isDispOK = kFALSE; 
  }
  else{//EMCAL
    
    if(l0 > fEMCALL0CutMin && l0 < fEMCALL0CutMax) isDispOK = kTRUE;

  }
  
  ph->SetDispBit(isDispOK) ;
  
  //TOF
  Double_t tof=cluster->GetTOF()  ;
  ph->SetTOFBit(TMath::Abs(tof)<fTOFCut) ; 
  
  //Charged 
  Bool_t isNeutral = IsTrackMatched(cluster,cu,event);
  
  ph->SetChargedBit(isNeutral);
  
  //Set PID pdg
  ph->SetIdentifiedParticleType(GetIdentifiedParticleType(calo,*ph->Momentum(),cluster));
  
  if(fDebug > 0){ 
    printf("AliCaloPID::SetPIDBits: TOF %e, Lambda0 %2.2f, Lambda1 %2.2f\n",tof , l0, l1); 	
    printf("AliCaloPID::SetPIDBits: pdg %d, bits: TOF %d, Dispersion %d, Charge %d\n",
           ph->GetIdentifiedParticleType(), ph->GetTOFBit() , ph->GetDispBit() , ph->GetChargedBit()); 
  }
}

//_____________________________________________________________________
Bool_t AliCaloPID::IsTrackMatched(AliVCluster* cluster,
                                  AliCalorimeterUtils * cu, 
                                  AliVEvent* event) const 
{
  //Check if there is any track attached to this cluster
  
  Int_t nMatches = cluster->GetNTracksMatched();
  AliVTrack * track = 0;
  Double_t p[3];

  if(nMatches > 0){
    
    //In case of ESDs, by default without match one entry with negative index, no match, reject.
    if(!strcmp("AliESDCaloCluster",Form("%s",cluster->ClassName())))
    {    
      Int_t iESDtrack = cluster->GetTrackMatchedIndex();
      if(iESDtrack >= 0) track = dynamic_cast<AliVTrack*> (event->GetTrack(iESDtrack));
      else return kFALSE;
      
      if (!track){
        printf("AliCaloPID::IsTrackMatched() - Null matched track in ESD when index is OK!\n");
        return kFALSE;
      }
    }      
    else { // AOD
      track = dynamic_cast<AliVTrack*> (cluster->GetTrackMatched(0));
      if (!track){
        printf("AliCaloPID::IsTrackMatched() - Null matched track in AOD!\n");
        return kFALSE;
      }
    }
    
    Float_t clE = cluster->E();
    Float_t dZ  = cluster->GetTrackDz();
    Float_t dR  = cluster->GetTrackDx();
    
    // if track matching was recalculated
    if(cluster->IsEMCAL() &&cu && cu->IsRecalculationOfClusterTrackMatchingOn()){
      dR = 2000., dZ = 2000.;
      cu->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dR,dZ);
    }
    
    // Fill control histograms
    if(fhTrackMatchedDEta && TMath::Abs(dR) < 999){
      fhTrackMatchedDEta->Fill(clE,dZ);
      fhTrackMatchedDPhi->Fill(clE,dR);
      if(clE > 0.5)fhTrackMatchedDEtaDPhi->Fill(dZ,dR);
      //printf("AliCaloPID::IsTrackMatched - %d dR %f , dZ %f \n",cluster->IsEMCAL(),dR, dZ);
    }  
    
    if(cluster->IsPHOS()) {
      
      track->GetPxPyPz(p) ;
      TLorentzVector trackmom(p[0],p[1],p[2],0);
      Int_t charge = track->Charge();
      Double_t mf  = event->GetMagneticField();
      if(TestPHOSChargedVeto(dR, dZ, trackmom.Pt(), charge, mf ) < fPHOSRCut) return kTRUE;
      else                                                                    return kFALSE;
      
    }//PHOS
    else {//EMCAL
      
      if(fDebug > 0) 
        printf("AliCaloPID::IsTrackMatched - EMCAL dR %f < %f, dZ %f < %f \n",dR, fEMCALDPhiCut, dZ, fEMCALDEtaCut);
      
      if(TMath::Abs(dR) < fEMCALDPhiCut && 
         TMath::Abs(dZ) < fEMCALDEtaCut)   return kTRUE;
      else                                 return kFALSE;
      
    }//EMCAL cluster 
    
    
  } // more than 1 match, at least one track in array
  else return kFALSE;
    
}

//___________________________________________________________________________________________________
Float_t AliCaloPID::TestPHOSDispersion(const Double_t pt, const Double_t l1, const Double_t l2) const 
{
  //Check if cluster photon-like. Uses photon cluster parameterization in real pp data 
  //Returns distance in sigmas. Recommended cut 2.5
  
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c       =-0.35-0.550*TMath::Exp(-0.390730*pt) ;
  Double_t R2      = 0.5*  (l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
                     0.5*  (l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
                     0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  
  if(fDebug > 0) printf("AliCaloPID::TestPHOSDispersion() - PHOS SS R %f < %f?\n", TMath::Sqrt(R2), fPHOSDispersionCut);
  
  return TMath::Sqrt(R2) ; 
  
  
}

//_______________________________________________________________________________________________
Float_t AliCaloPID::TestPHOSChargedVeto(const Double_t dx,  const Double_t dz, const Double_t pt, 
                                        const Int_t charge, const Double_t mf) const 
{
  //Checks distance to the closest track. Takes into account 
  //non-perpendicular incidence of tracks.
  //returns distance in sigmas. Recommended cut: 2.
  //Requires (sign) of magnetic filed. onc can find it for example as following
  //  Double_t mf=0. ;
  //  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  //  if(event)
  //    mf = event->GetMagneticField(); //Positive for ++ and negative for --
  
  
  Double_t meanX = 0.;
  Double_t meanZ = 0.;
  Double_t sx = TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
                           6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+
                           1.59219);
  Double_t sz = TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+
                           1.60) ;
  
  if(mf<0.){ //field --
    meanZ = -0.468318 ;
    if(charge>0)
      meanX = TMath::Min(7.3, 3.89994*1.20679 *1.20679 /(pt*pt+1.20679*1.20679)+  
                         0.249029+2.49088e+07*TMath::Exp(-pt*3.33650e+01)) ;
    else
      meanX =-TMath::Min(7.7, 3.86040*0.912499*0.912499/(pt*pt+0.912499*0.912499)+
                         1.23114 +4.48277e+05*TMath::Exp(-pt*2.57070e+01)) ;
  }
  else{ //Field ++
    meanZ = -0.468318;
    if(charge>0)
      meanX =-TMath::Min(8.0,3.86040*1.31357*1.31357/(pt*pt+1.31357*1.31357)+
                         0.880579+7.56199e+06*TMath::Exp(-pt*3.08451e+01)) ;
    else
      meanX = TMath::Min(6.85, 3.89994*1.16240*1.16240/(pt*pt+1.16240*1.16240)-
                         0.120787+2.20275e+05*TMath::Exp(-pt*2.40913e+01)) ;     
  }
  
  Double_t rz = (dz-meanZ)/sz ;
  Double_t rx = (dx-meanX)/sx ;
  
  if(fDebug > 0) 
    printf("AliCaloPID::TestPHOSDispersion() - PHOS Matching R %f < %f\n",TMath::Sqrt(rx*rx+rz*rz), fPHOSRCut);
  
  return TMath::Sqrt(rx*rx+rz*rz) ;
  
}
