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

//_________________________________________________________________________
// Class for PID selection with calorimeters
// The Output of the main method GetIdentifiedParticleType is a PDG number identifying the cluster, 
// being kPhoton, kElectron, kPi0 ... as defined in the header file
//   - GetIdentifiedParticleType(const AliVCluster * cluster) 
//      Assignes a PID tag to the cluster, right now there is the possibility to : use bayesian weights from reco, 
//      recalculate them (EMCAL) or use other procedures not used in reco.
//      In order to recalculate Bayesian, it is necessary to load the EMCALUtils library
//      and do SwitchOnBayesianRecalculation().
//      To change the PID parameters from Low to High like the ones by default, use the constructor 
//      AliCaloPID(flux)
//      where flux is AliCaloPID::kLow or AliCaloPID::kHigh
//      If it is necessary to change the parameters use the constructor 
//      AliCaloPID(AliEMCALPIDUtils *utils) and set the parameters before.

//   - GetGetIdentifiedParticleTypeFromBayesian(const Double_t * pid, const Float_t energy)
//      Reads the PID weights array of the ESDs and depending on its magnitude identifies the particle, 
//      executed when bayesian is ON by GetIdentifiedParticleType(const AliVCluster * cluster) 
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
#include "AliAODCaloCluster.h"
#include "AliVCaloCells.h"
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
//Split
fDoClusterSplitting(kFALSE),
fUseSimpleMassCut(kFALSE),
fUseSimpleM02Cut(kFALSE),
fUseSplitAsyCut(kFALSE),
fSplitM02MaxCut(0),       fSplitM02MinCut(0),          fSplitMinNCells(0),
fMassEtaMin(0),           fMassEtaMax(0),
fMassPi0Min(0),           fMassPi0Max(0),
fMassPhoMin(0),           fMassPhoMax(0),
fSplitWidthSigma(0)
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
//Split
fDoClusterSplitting(kFALSE),
fUseSimpleMassCut(kFALSE),
fUseSimpleM02Cut(kFALSE),
fUseSplitAsyCut(kFALSE),
fSplitM02MaxCut(0),       fSplitM02MinCut(0),          fSplitMinNCells(0),
fMassEtaMin(0),           fMassEtaMax(0),
fMassPi0Min(0),           fMassPi0Max(0),
fMassPhoMin(0),           fMassPhoMax(0),
fSplitWidthSigma(0)
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
//Split
fDoClusterSplitting(kFALSE),
fUseSimpleMassCut(kFALSE),
fUseSimpleM02Cut(kFALSE),
fUseSplitAsyCut(kFALSE),
fSplitM02MaxCut(0),       fSplitM02MinCut(0),          fSplitMinNCells(0),
fMassEtaMin(0),           fMassEtaMax(0),
fMassPi0Min(0),           fMassPi0Max(0),
fMassPhoMin(0),           fMassPhoMax(0),
fSplitWidthSigma(0)

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
  
  // Cluster splitting
  
  fSplitM02MinCut = 0.3 ;
  fSplitM02MaxCut = 5   ;
  fSplitMinNCells = 4   ;

  fMassEtaMin  = 0.4;
  fMassEtaMax  = 0.6;
  
  fMassPi0Min  = 0.11;
  fMassPi0Max  = 0.18;
  
  fMassPhoMin  = 0.0;
  fMassPhoMax  = 0.08;
  
  fMassWidthPi0Param[0] = 0.100;  // Absolute Low mass cut for NLM=1 and E < 10 GeV
  fMassWidthPi0Param[1] = 0.050;  // Absolute Low mass cut for NLM=2 and E < 10 GeV
  fMassWidthPi0Param[2] = 0.009;  // constant width for E < 8 GeV, 9 MeV
  fMassWidthPi0Param[3] = 0.0023; // pol1 param0 of width for 8 < E < 16 GeV
  fMassWidthPi0Param[4] = 0.0008; // pol1 param1 of width for 8 < E < 16 GeV
  fMassWidthPi0Param[5] = 0.130;  // Mean mass value for NLM=1, pp (121 PbPb)
  fMassWidthPi0Param[6] = 0.134;  // Mean mass value for NLM=2, pp (125 PbPb)
  
  
  fM02MinParam[0][0] = 5.762   ; // pol3 param0 for NLM=1 , E < 20 GeV, pp/PbPb
  fM02MinParam[0][1] =-0.880   ; // pol3 param1 for NLM=1 , E < 20 GeV, pp/PbPb
  fM02MinParam[0][2] = 0.0487  ; // pol3 param2 for NLM=1 , E < 20 GeV, pp/PbPb
  fM02MinParam[0][3] =-0.000913; // pol3 param2 for NLM=1 , E < 20 GeV, pp/PbPb
  fM02MinParam[0][4] = 0.3     ; // cut for E > 20 GeV, pp/PbPb
  fM02MinParam[0][5] = 20.     ; // E cut change

  fM02MinParam[1][0] = 8.297  ;  // pol3 param0 for NLM>2 , E < 14 GeV, pp/PbPb
  fM02MinParam[1][1] =-1.485  ;  // pol3 param1 for NLM>2 , E < 14 GeV, pp/PbPb
  fM02MinParam[1][2] = 0.0949 ;  // pol3 param2 for NLM>2 , E < 14 GeV, pp/PbPb
  fM02MinParam[1][3] =-0.00200;  // pol3 param2 for NLM>2 , E < 14 GeV, pp/PbPb
  fM02MinParam[1][4] = 0.6;     // cut for E > 14 GeV, pp/PbPb
  fM02MinParam[1][5] = 14.;     // E cut change

  fM02MaxParam[0][0] = 4.826   ; // pol3 param0 for NLM=1 , E < 20 GeV, pp/PbPb
  fM02MaxParam[0][1] =-0.395   ; // pol3 param1 for NLM=1 , E < 20 GeV, pp/PbPb
  fM02MaxParam[0][2] = 0.0123  ; // pol3 param2 for NLM=1 , E < 20 GeV, pp/PbPb
  fM02MaxParam[0][3] =-0.000121; // pol3 param2 for NLM=1 , E < 20 GeV, pp/PbPb
  fM02MaxParam[0][4] = 0.75    ; // cut for E > 23 GeV, pp/PbPb
  fM02MaxParam[0][5] = 25.     ; // E cut change
  
  fM02MaxParam[1][0] =11.4    ;  // pol3 param0 for NLM>2 , E < 14 GeV, pp/PbPb
  fM02MaxParam[1][1] =-1.466  ;  // pol3 param1 for NLM>2 , E < 14 GeV, pp/PbPb
  fM02MaxParam[1][2] = 0.0726 ;  // pol3 param2 for NLM>2 , E < 14 GeV, pp/PbPb
  fM02MaxParam[1][3] =-0.00121;  // pol3 param2 for NLM>2 , E < 14 GeV, pp/PbPb
  fM02MaxParam[1][4] = 1.45  ;  // cut for E > 21 GeV, pp/PbPb
  fM02MaxParam[1][5] = 21.   ;  // E cut change
  
  fAsyMinParam[0][0] =-0.245  ;  // pol3 param0 for NLM=1 , E < 25 GeV, pp
  fAsyMinParam[0][1] = 0.0873 ;  // pol3 param1 for NLM=1 , E < 25 GeV, pp
  fAsyMinParam[0][2] =-0.00173;  // pol3 param2 for NLM=1 , E < 25 GeV, pp
  fAsyMinParam[0][3] = 0     ;  // pol3 param2 for NLM=1 , E < 25 GeV, pp
  fAsyMinParam[0][4] = 1.0   ;  // cut for NLM=1 , E > 25 GeV, pp/PbPb
  fAsyMinParam[0][5] = 25.   ;  // E cut change

  fAsyMinParam[1][0] =-1.31  ;   // pol3 param0 for NLM>2 , E < 18 GeV, pp
  fAsyMinParam[1][1] = 0.387 ;   // pol3 param1 for NLM>2 , E < 18 GeV, pp
  fAsyMinParam[1][2] =-0.0217;   // pol3 param2 for NLM>2 , E < 18 GeV, pp
  fAsyMinParam[1][3] = 0.000409; // pol3 param2 for NLM>2 , E < 18 GeV, pp
  fAsyMinParam[1][4] = 1.0  ;   // cut for NLM>2 , E > 18 GeV, pp/PbPb
  fAsyMinParam[1][5] = 18.  ;   // E cut change
  
//  fAsyMinParam[0][0] =-0.663 ;  // pol3 param0 for NLM=1 , E < 25 GeV, PbPb
//  fAsyMinParam[0][1] = 0.131 ;  // pol3 param1 for NLM=1 , E < 25 GeV, PbPb
//  fAsyMinParam[0][2] =-0.00278;  // pol3 param2 for NLM=1 , E < 25 GeV, PbPb
//  fAsyMinParam[0][3] = 0     ;  // pol3 param2 for NLM=1 , E < 25 GeV, PbPb

//  fAsyMinParam[1][0] =-1.308 ;   // pol3 param0 for NLM>2 , E < 25 GeV, PbPb
//  fAsyMinParam[1][1] = 0.257 ;   // pol3 param1 for NLM>2 , E < 25 GeV, PbPb
//  fAsyMinParam[1][2] =-0.00788;  // pol3 param2 for NLM>2 , E < 25 GeV, PbPb
//  fAsyMinParam[1][3] = 0.0000377; // pol3 param2 for NLM>2 , E < 25 GeV, PbPb

  
  fSplitEFracMin[0]   = 0.0 ; // 0.96
  fSplitEFracMin[1]   = 0.0 ; // 0.96
  fSplitEFracMin[2]   = 0.0 ; // 0.7

  fSplitWidthSigma = 3. ;
  
}


//_____________________________________________________________________________________________________
Bool_t AliCaloPID::IsInPi0SplitAsymmetryRange(const Float_t energy, const Float_t asy, const Int_t nlm)
{
  // Select the appropriate mass range for pi0 selection in splitting method
  // No used yet in splitting ID decision
  
  Float_t abasy = TMath::Abs(asy);

  Int_t inlm = nlm-1;
  if(nlm > 2) inlm=1; // only 2 cases defined nlm=1 and nlm>=2
  
  // Get the parametrized min cut of asymmetry for NLM=2 up to 11 GeV
  Float_t cut = fAsyMinParam[inlm][0]+
                fAsyMinParam[inlm][1]*energy+
                fAsyMinParam[inlm][2]*energy*energy+
                fAsyMinParam[inlm][3]*energy*energy*energy;

  // In any case and beyond validity energy range of the function,
  // the parameter cannot be smaller than 1
  if(cut > fAsyMinParam[inlm][4] || energy > fAsyMinParam[inlm][5] ) cut = fAsyMinParam[inlm][4];
  
  //printf("energy %2.2f - nlm: %d (%d)- p0 %f, p1 %f, p2 %f, p3 %f ; cut: %2.2f\n",energy,nlm,inlm,
  //       fAsyMinParam[inlm][0],fAsyMinParam[inlm][1],fAsyMinParam[inlm][2],fAsyMinParam[inlm][3],cut);
  
  if(abasy < cut) return kTRUE;
  else            return kFALSE;
  
}

//_________________________________________________________________________________________________
Bool_t AliCaloPID::IsInPi0SplitMassRange(const Float_t energy, const Float_t mass, const Int_t nlm)
{
  // Select the appropriate mass range for pi0 selection in splitting method
  
  if(fUseSimpleMassCut)
  {
    if(mass < fMassPi0Max && mass > fMassPi0Min) return kTRUE;
    else                                         return kFALSE;
  }
  
  // Get the selected mean value as reference for the mass
  Float_t       meanMass = fMassWidthPi0Param[6];
  if(nlm == 1)  meanMass = fMassWidthPi0Param[5];
  
  if(energy > 16) meanMass+= (energy*fMassWidthPi0Param[4]-0.013); // Increase mean mass with energy from 16 GeV, same slope as for width
    
  // Get the parametrized width of the mass
  Float_t width   = 0.009;
  if      (energy < 8 ) width = fMassWidthPi0Param[2];
  else if (energy < 16) width = fMassWidthPi0Param[3]+energy*fMassWidthPi0Param[4];
  else                  width = fMassWidthPi0Param[3]+16    *fMassWidthPi0Param[4]; // Fixed value at 16 GeV
  
  // Calculate the 2 sigma cut
  Float_t minMass = meanMass-fSplitWidthSigma*width;
  Float_t maxMass = meanMass+fSplitWidthSigma*width;

  // In case of low energy, hard cut to avoid conversions
  if(energy < 10  && nlm == 1 && minMass < fMassWidthPi0Param[0] ) minMass = fMassWidthPi0Param[0];
  if(energy < 10  && nlm == 2 && minMass < fMassWidthPi0Param[1] ) minMass = fMassWidthPi0Param[1];
  
  //printf("\t \t sigma %1.1f width %3.1f, mean Mass %3.0f, minMass %3.0f, maxMass %3.0f\n ", 
  //       fSplitWidthSigma, width*1000, meanMass*1000,minMass*1000,maxMass*1000);
  
  if(mass < maxMass && mass > minMass) return kTRUE;
  else                                 return kFALSE;
  
}

//_____________________________________________________________________________________________
Bool_t AliCaloPID::IsInPi0M02Range(const Float_t energy, const Float_t m02,  const Int_t nlm)
{
  // Select the appropriate m02 range in splitting method for pi0
  
  Float_t minCut = fSplitM02MinCut;
  Float_t maxCut = fSplitM02MaxCut;
  
  if(!fUseSimpleM02Cut)
  {
    Int_t inlm = nlm-1;
    if(nlm > 2) inlm=1; // only 2 cases defined nlm=1 and nlm>=2
    
    minCut = fM02MinParam[inlm][0]+
             fM02MinParam[inlm][1]*energy+
             fM02MinParam[inlm][2]*energy*energy+
             fM02MinParam[inlm][3]*energy*energy*energy;
    
    maxCut = fM02MaxParam[inlm][0]+
             fM02MaxParam[inlm][1]*energy+
             fM02MaxParam[inlm][2]*energy*energy+
             fM02MaxParam[inlm][3]*energy*energy*energy;

    
    // In any case and beyond validity energy range of the function,
    // the parameter cannot be smaller than this (0.3 for nlm=1 and 0.6 for the rest)
    if( minCut < fM02MinParam[inlm][4] || energy > fM02MinParam[inlm][5] ) minCut = fM02MinParam[inlm][4];
    if( maxCut < fM02MaxParam[inlm][4] || energy > fM02MaxParam[inlm][5] ) maxCut = fM02MaxParam[inlm][4];
    if( nlm>2 ) maxCut+=0.75;
  }
  
  //if(energy > 7) printf("\t \t E %2.2f, nlm %d, m02 %2.2f, minM02 %2.2f, maxM02 %2.2f\n",energy, nlm, m02,minCut,maxCut);
  
  if(m02 < maxCut && m02 > minCut) return kTRUE;
  else                             return kFALSE;

}


//_____________________________________________________________________________________________
Bool_t AliCaloPID::IsInEtaM02Range(const Float_t energy, const Float_t m02,  const Int_t nlm)
{
  // Select the appropriate m02 range in splitting method to select eta's
  // Use same parametrization as pi0, just shift the distributions (to be tuned)
  
  Float_t minCut = fSplitM02MinCut;
  Float_t maxCut = fSplitM02MaxCut;
  
  if(!fUseSimpleM02Cut)
  {
    Int_t inlm = nlm-1;
    if(nlm > 2) inlm=1; // only 2 cases defined nlm=1 and nlm>=2
    
    Float_t shiftE = energy-20; // to be tuned
    if(nlm==1) shiftE=energy-28;
    
    minCut = fM02MinParam[inlm][0]+
    fM02MinParam[inlm][1]*shiftE+
    fM02MinParam[inlm][2]*shiftE*shiftE+
    fM02MinParam[inlm][3]*shiftE*shiftE*shiftE;
    
    // In any case and beyond validity energy range of the function,
    // the parameter cannot be smaller than this (0.3 for nlm=1 and 0.6 for the rest)
    Float_t minE = 50;
    if(nlm>1) minE = 35;
    if( minCut < fM02MinParam[inlm][4]+1 || shiftE > minE ) minCut = fM02MinParam[inlm][4]+1;

    shiftE = energy+20; // to be tuned

    maxCut = 1+fM02MaxParam[inlm][0]+
    fM02MaxParam[inlm][1]*shiftE+
    fM02MaxParam[inlm][2]*shiftE*shiftE+
    fM02MaxParam[inlm][3]*shiftE*shiftE*shiftE;
    
    Float_t maxE = 50;
    if(nlm>1) maxE = 40;
    if( maxCut < fM02MaxParam[inlm][4]+1 || shiftE > maxE ) maxCut = fM02MaxParam[inlm][4]+1;
    if( nlm>2 ) maxCut+=0.75;
  }
  
  //if(energy>6)printf("\t \t E %2.2f, nlm %d, m02 %2.2f, minM02 %2.2f, maxM02 %2.2f\n",energy, nlm, m02,minCut,maxCut);
  
  if(m02 < maxCut && m02 > minCut) return kTRUE;
  else                             return kFALSE;
  
}

//_____________________________________________________________________________________________
Bool_t AliCaloPID::IsInConM02Range(const Float_t energy, const Float_t m02,  const Int_t nlm)
{
  // Select the appropriate m02 range in splitting method for converted photons
  // Just min limit for pi0s is max for conversion.
  
  Float_t minCut = 0.1;
  Float_t maxCut = 0.3;
  
  if(!fUseSimpleM02Cut)
  {
    Int_t inlm = nlm-1;
    if(nlm > 2) inlm=1; // only 2 cases defined nlm=1 and nlm>=2
    
    maxCut = fM02MinParam[inlm][0]+
    fM02MinParam[inlm][1]*energy+
    fM02MinParam[inlm][2]*energy*energy+
    fM02MinParam[inlm][3]*energy*energy*energy;
    
    // In any case and beyond validity energy range of the function,
    // the parameter cannot be smaller than this (0.3 for nlm=1 and 0.6 for the rest)
    if( maxCut < fM02MaxParam[inlm][4] || energy > fM02MaxParam[inlm][5] ) maxCut = fM02MaxParam[inlm][4];
  }
  
  
  if(m02 < maxCut && m02 > minCut) return kTRUE;
  else                             return kFALSE;
  
}

//______________________________________________
AliEMCALPIDUtils *AliCaloPID::GetEMCALPIDUtils() 
{
  // return pointer to AliEMCALPIDUtils, create it if needed
  
  if(!fEMCALPIDUtils) fEMCALPIDUtils = new AliEMCALPIDUtils ; 
  return fEMCALPIDUtils ; 
  
}


//______________________________________________________________________
Int_t AliCaloPID::GetIdentifiedParticleType(const AliVCluster * cluster) 
{
  // Returns a PDG number corresponding to the likely ID of the cluster
  
  Float_t energy  = cluster->E();	
  Float_t lambda0 = cluster->GetM02();
  Float_t lambda1 = cluster->GetM20();
  
  // ---------------------
  // Use bayesian approach
  // ---------------------
  
  if(fUseBayesianWeights)
  {
    Double_t weights[AliPID::kSPECIESCN];
    
    if(cluster->IsEMCAL() && fRecalculateBayesian)
    {	        
      fEMCALPIDUtils->ComputePID(energy, lambda0);
      for(Int_t i = 0; i < AliPID::kSPECIESCN; i++) weights[i] = fEMCALPIDUtils->GetPIDFinal(i);
    }
    else 
    {
      for(Int_t i = 0; i < AliPID::kSPECIESCN; i++) weights[i] = cluster->GetPID()[i];
    }

    if(fDebug > 0)  PrintClusterPIDWeights(weights);
    
    return GetIdentifiedParticleTypeFromBayesWeights(cluster->IsEMCAL(), weights, energy);
  }
  
  // -------------------------------------------------------
  // Calculate PID SS from data, do not use bayesian weights
  // -------------------------------------------------------
  
  if(fDebug > 0) printf("AliCaloPID::GetIdentifiedParticleType: EMCAL %d?, E %3.2f, l0 %3.2f, l1 %3.2f, disp %3.2f, tof %1.11f, distCPV %3.2f, distToBC %1.1f, NMax %d\n",
                        cluster->IsEMCAL(),energy,lambda0,cluster->GetM20(),cluster->GetDispersion(),cluster->GetTOF(), 
                        cluster->GetEmcCpvDistance(), cluster->GetDistanceToBadChannel(),cluster->GetNExMax());
  
  if(cluster->IsEMCAL())
  {
    if(fDebug > 0) printf("AliCaloPID::GetIdentifiedParticleType() - EMCAL SS %f <%f < %f?\n",fEMCALL0CutMin, lambda0, fEMCALL0CutMax);
    
    if(lambda0 < fEMCALL0CutMax && lambda0 > fEMCALL0CutMin) return kPhoton ;
    else                                                     return kNeutralUnknown ; 
  }    // EMCAL
  else // PHOS
  {    
    if(TestPHOSDispersion(energy,lambda0,lambda1) < fPHOSDispersionCut) return kPhoton;
    else                                                                return kNeutralUnknown;
  }
  
}

//_______________________________________________________________________________
Int_t AliCaloPID::GetIdentifiedParticleTypeFromBayesWeights(const Bool_t isEMCAL, 
                                                            const Double_t * pid, 
                                                            const Float_t energy) 
{
  //Return most probable identity of the particle after bayesian weights calculated in reconstruction
  
  if(!pid)
  { 
    printf("AliCaloPID::GetIdentifiedParticleType() - pid pointer not initialized!!!\n");
    abort();
  }
  
  Float_t wPh  =  fPHOSPhotonWeight ;
  Float_t wPi0 =  fPHOSPi0Weight ;
  Float_t wE   =  fPHOSElectronWeight ;
  Float_t wCh  =  fPHOSChargeWeight ;
  Float_t wNe  =  fPHOSNeutralWeight ;
  
  if(!isEMCAL && fPHOSWeightFormula){
    wPh  = GetPHOSPhotonWeightFormula()->Eval(energy) ;
    wPi0 = GetPHOSPi0WeightFormula()   ->Eval(energy);
  }
  else
  {
    wPh  =  fEMCALPhotonWeight ;
    wPi0 =  fEMCALPi0Weight ;
    wE   =  fEMCALElectronWeight ;
    wCh  =  fEMCALChargeWeight ;
    wNe  =  fEMCALNeutralWeight ;
  }
  
  if(fDebug > 0) PrintClusterPIDWeights(pid);
    
  Int_t pdg = kNeutralUnknown ;
  Float_t chargedHadronWeight = pid[AliVCluster::kProton]+pid[AliVCluster::kKaon]+
  pid[AliVCluster::kPion]+pid[AliVCluster::kMuon];
  Float_t neutralHadronWeight = pid[AliVCluster::kNeutron]+pid[AliVCluster::kKaon0];
  Float_t allChargedWeight    = pid[AliVCluster::kElectron]+pid[AliVCluster::kEleCon]+ chargedHadronWeight;
  Float_t allNeutralWeight    = pid[AliVCluster::kPhoton]+pid[AliVCluster::kPi0]+ neutralHadronWeight;
  
  //Select most probable ID
  if(!isEMCAL) // PHOS
  {
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
  else //EMCAL
  {
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

//____________________________________________________________________________________________________
Int_t AliCaloPID::GetIdentifiedParticleTypeFromClusterSplitting(AliVCluster* cluster, 
                                                                AliVCaloCells* cells,
                                                                AliCalorimeterUtils * caloutils,
                                                                Double_t   vertex[3],
                                                                Int_t    & nMax,
                                                                Double_t & mass,   Double_t & angle,
                                                                Double_t & e1  ,   Double_t & e2,
                                                                Int_t    & absId1, Int_t    & absId2 )
{
  // Split the cluster in 2, do invariant mass, get the mass and decide 
  // if this is a photon, pi0, eta, ...
  
  Float_t eClus  = cluster->E();
  Float_t m02    = cluster->GetM02();
  const Int_t nc = cluster->GetNCells();
  Int_t   absIdList[nc]; 
  Float_t maxEList [nc];
  
  mass  = -1.;
  angle = -1.;
  
  //If too low number of cells, skip it
  if ( nc < fSplitMinNCells)  return kNeutralUnknown ; 
  
  if(fDebug > 0) printf("\t pass nCells cut\n");
  
  // Get Number of local maxima
  nMax  = caloutils->GetNumberOfLocalMaxima(cluster, cells, absIdList, maxEList) ; 
  
  if(fDebug > 0) printf("AliCaloPID::GetIdentifiedParticleTypeFromClusterSplitting() - Cluster : E %1.1f, M02 %1.2f, NLM %d, N Cells %d\n",
                        eClus,m02,nMax,nc);

  //---------------------------------------------------------------------
  // Get the 2 max indeces and do inv mass
  //---------------------------------------------------------------------
  
  TString  calorimeter = "EMCAL";
  if(cluster->IsPHOS()) calorimeter = "PHOS";

  if     ( nMax == 2 )
  {
    absId1 = absIdList[0];
    absId2 = absIdList[1];
    
    //Order in energy
    Float_t en1 = cells->GetCellAmplitude(absId1);
    caloutils->RecalibrateCellAmplitude(en1,calorimeter,absId1);
    Float_t en2 = cells->GetCellAmplitude(absId2);
    caloutils->RecalibrateCellAmplitude(en2,calorimeter,absId2);
    if(en1 < en2)
    {
      absId2 = absIdList[0];
      absId1 = absIdList[1];
    }
  }
  else if( nMax == 1 )
  {
    
    absId1 = absIdList[0];
    
    //Find second highest energy cell
    
    Float_t enmax = 0 ;
    for(Int_t iDigit = 0 ; iDigit < cluster->GetNCells() ; iDigit++)
    {
      Int_t absId = cluster->GetCellsAbsId()[iDigit];
      if( absId == absId1 ) continue ; 
      Float_t endig = cells->GetCellAmplitude(absId);
      caloutils->RecalibrateCellAmplitude(endig,calorimeter,absId); 
      if(endig > enmax) 
      {
        enmax  = endig ;
        absId2 = absId ;
      }
    }// cell loop
  }// 1 maxima 
  else
  {  // n max > 2
    // loop on maxima, find 2 highest
    
    // First max
    Float_t enmax = 0 ;
    for(Int_t iDigit = 0 ; iDigit < nMax ; iDigit++)
    {
      Float_t endig = maxEList[iDigit];
      if(endig > enmax) 
      {
        enmax  = endig ;
        absId1 = absIdList[iDigit];
      }
    }// first maxima loop
    
    // Second max 
    Float_t enmax2 = 0;
    for(Int_t iDigit = 0 ; iDigit < nMax ; iDigit++)
    {
      if(absIdList[iDigit]==absId1) continue;
      Float_t endig = maxEList[iDigit];
      if(endig > enmax2) 
      {
        enmax2  = endig ;
        absId2 = absIdList[iDigit];
      }
    }// second maxima loop
    
  } // n local maxima > 2
  
  if(absId2<0 || absId1<0) 
  {
    if(fDebug > 0) printf("AliCaloPID::GetIdentifiedParticleTypeFromClusterSplitting() - Bad index for local maxima : N max %d, i1 %d, i2 %d, cluster E %2.2f, ncells %d, m02 %2.2f\n",
                          nMax,absId1,absId2,eClus,nc,m02);
    return kNeutralUnknown ; 
  }
  
  //---------------------------------------------------------------------
  // Split the cluster energy in 2, around the highest 2 local maxima
  //---------------------------------------------------------------------  
  
  AliAODCaloCluster cluster1(0, 0,NULL,0.,NULL,NULL,1,0);
  AliAODCaloCluster cluster2(1, 0,NULL,0.,NULL,NULL,1,0);
  
  caloutils->SplitEnergy(absId1,absId2,cluster, cells, &cluster1, &cluster2,nMax); /*absIdList, maxEList,*/
  
  TLorentzVector cellMom1; 
  TLorentzVector cellMom2;  
  
  cluster1.GetMomentum(cellMom1,vertex);
  cluster2.GetMomentum(cellMom2,vertex);
  
  mass  = (cellMom1+cellMom2).M();
  angle = cellMom2.Angle(cellMom1.Vect());
  e1    = cluster1.E();
  e2    = cluster2.E();
    
  // Consider clusters with splitted energy not too different to original cluster energy
  Float_t splitFracCut = 0;
  if(nMax < 3)  splitFracCut = fSplitEFracMin[nMax-1];
  else          splitFracCut = fSplitEFracMin[2];
  if((e1+e2)/eClus < splitFracCut) return kNeutralUnknown ;
  
  if(fDebug > 0) printf("\t pass Split E frac cut\n");
    
  // Asymmetry of cluster
  Float_t asy =-10;
  if(e1+e2 > 0) asy = (e1-e2) / (e1+e2);
  if( fUseSplitAsyCut &&  !IsInPi0SplitAsymmetryRange(eClus,asy,nMax) ) return kNeutralUnknown ;
  
  if (fDebug>0) printf("\t pass asymmetry cut\n");
  
  Bool_t pi0OK = kFALSE;
  Bool_t etaOK = kFALSE;
  Bool_t conOK = kFALSE;
  
  //If too small or big M02, skip it
  if     (IsInPi0M02Range(eClus,m02,nMax))  pi0OK = kTRUE;
  else if(IsInEtaM02Range(eClus,m02,nMax))  etaOK = kTRUE;
  else if(IsInConM02Range(eClus,m02,nMax))  conOK = kTRUE;
  
  // Check the mass, and set an ID to the splitted cluster
  if     ( conOK && mass < fMassPhoMax && mass > fMassPhoMin     ) { if(fDebug > 0) printf("\t Split Conv \n"); return kPhoton ; }
  else if( etaOK && mass < fMassEtaMax && mass > fMassEtaMin     ) { if(fDebug > 0) printf("\t Split Eta \n");  return kEta    ; }
  else if( pi0OK && IsInPi0SplitMassRange(cluster->E(),mass,nMax)) { if(fDebug > 0) printf("\t Split Pi0 \n");  return kPi0    ; }
  else                                                                                                          return kNeutralUnknown ;
  
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
  
  if(fDoClusterSplitting)
  {
    snprintf(onePar,buffersize,"%2.2f< M02 < %2.2f \n",    fSplitM02MinCut, fSplitM02MaxCut) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fMinNCells =%d \n",        fSplitMinNCells) ;
    parList+=onePar ;    
    snprintf(onePar,buffersize,"pi0 : %2.1f < m <%2.1f\n", fMassPi0Min,fMassPi0Max);
    parList+=onePar ;
    snprintf(onePar,buffersize,"eta : %2.1f < m <%2.1f\n", fMassEtaMin,fMassEtaMax);
    parList+=onePar ;
    snprintf(onePar,buffersize,"conv: %2.1f < m <%2.1f\n", fMassPhoMin,fMassPhoMax);
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
  
  if(fUseBayesianWeights)
  {
    printf("PHOS PID weight , photon %0.2f, pi0 %0.2f, e %0.2f, charge %0.2f, neutral %0.2f \n",  
           fPHOSPhotonWeight,   fPHOSPi0Weight, 
           fPHOSElectronWeight, fPHOSChargeWeight, fPHOSNeutralWeight) ; 
    printf("EMCAL PID weight, photon %0.2f, pi0 %0.2f, e %0.2f, charge %0.2f, neutral %0.2f\n",   
           fEMCALPhotonWeight,   fEMCALPi0Weight, 
           fEMCALElectronWeight, fEMCALChargeWeight, fEMCALNeutralWeight) ; 
    
    printf("PHOS Parametrized weight on?  =     %d\n",  fPHOSWeightFormula) ; 
    if(fPHOSWeightFormula)
    {
      printf("Photon weight formula = %s\n", fPHOSPhotonWeightFormulaExpression.Data());
      printf("Pi0    weight formula = %s\n", fPHOSPi0WeightFormulaExpression   .Data());
    }
    if(fRecalculateBayesian) printf(" Recalculate bayesian with Particle Flux?    = %d\n",fParticleFlux);
  }
  else 
  {
    printf("TOF cut        = %e\n",                                   fTOFCut);
    printf("EMCAL Lambda0 cut min = %2.2f; max = %2.2f\n",            fEMCALL0CutMin,fEMCALL0CutMax);
    printf("EMCAL cluster-track dEta < %2.3f; dPhi < %2.3f\n",        fEMCALDEtaCut, fEMCALDPhiCut);
    printf("PHOS Treac matching cut =%2.2f, Dispersion Cut =%2.2f \n",fPHOSRCut,     fPHOSDispersionCut) ;
    
  }
  
  if(fDoClusterSplitting)
  {
    printf("Min. N Cells =%d \n",         fSplitMinNCells) ;
    printf("%2.2f < lambda_0^2 <%2.2f \n",fSplitM02MinCut,fSplitM02MaxCut);
    printf("pi0 : %2.2f<m<%2.2f \n",      fMassPi0Min,fMassPi0Max);
    printf("eta : %2.2f<m<%2.2f \n",      fMassEtaMin,fMassEtaMax);
    printf("phot: %2.2f<m<%2.2f \n",      fMassPhoMin,fMassPhoMax);
  }
  
  printf(" \n");
  
} 

//_________________________________________________________________
void AliCaloPID::PrintClusterPIDWeights(const Double_t * pid) const
{
  // print PID of cluster, (AliVCluster*)cluster->GetPID()
  
  printf("AliCaloPID::PrintClusterPIDWeights() \n \t ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, \n \t \
         pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
         pid[AliVCluster::kPhoton],    pid[AliVCluster::kPi0],
         pid[AliVCluster::kElectron],  pid[AliVCluster::kEleCon],
         pid[AliVCluster::kPion],      pid[AliVCluster::kKaon], 
         pid[AliVCluster::kProton],
         pid[AliVCluster::kNeutron],   pid[AliVCluster::kKaon0]);
  
}

//___________________________________________________________________________
void AliCaloPID::SetPIDBits(AliVCluster * cluster, 
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
  ph->SetIdentifiedParticleType(GetIdentifiedParticleType(cluster));
  
  if(fDebug > 0)
  { 
    printf("AliCaloPID::SetPIDBits: TOF %e, Lambda0 %2.2f, Lambda1 %2.2f\n",tof , l0, l1); 	
    printf("AliCaloPID::SetPIDBits: pdg %d, bits: TOF %d, Dispersion %d, Charge %d\n",
           ph->GetIdentifiedParticleType(), ph->GetTOFBit() , ph->GetDispBit() , ph->GetChargedBit()); 
  }
}

//_________________________________________________________
Bool_t AliCaloPID::IsTrackMatched(AliVCluster* cluster,
                                  AliCalorimeterUtils * cu, 
                                  AliVEvent* event) const 
{
  //Check if there is any track attached to this cluster
  
  Int_t nMatches = cluster->GetNTracksMatched();
  AliVTrack * track = 0;
  Double_t p[3];

  if(nMatches > 0)
  {
    //In case of ESDs, by default without match one entry with negative index, no match, reject.
    if(!strcmp("AliESDCaloCluster",Form("%s",cluster->ClassName())))
    {    
      Int_t iESDtrack = cluster->GetTrackMatchedIndex();
      if(iESDtrack >= 0) track = dynamic_cast<AliVTrack*> (event->GetTrack(iESDtrack));
      else return kFALSE;
      
      if (!track)
      {
        if(fDebug > 0) printf("AliCaloPID::IsTrackMatched() - Null matched track in ESD when index is OK!\n");
        return kFALSE;
      }
    }      
    else { // AOD
      track = dynamic_cast<AliVTrack*> (cluster->GetTrackMatched(0));
      if (!track)
      {
        if(fDebug > 0) printf("AliCaloPID::IsTrackMatched() - Null matched track in AOD!\n");
        return kFALSE;
      }
    }
    
    Float_t dZ  = cluster->GetTrackDz();
    Float_t dR  = cluster->GetTrackDx();
    
    // if track matching was recalculated
    if(cluster->IsEMCAL() && cu && cu->IsRecalculationOfClusterTrackMatchingOn())
    {
      dR = 2000., dZ = 2000.;
      cu->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
    }
        
    if(cluster->IsPHOS()) 
    {
      track->GetPxPyPz(p) ;
      TLorentzVector trackmom(p[0],p[1],p[2],0);
      Int_t charge = track->Charge();
      Double_t mf  = event->GetMagneticField();
      if(TestPHOSChargedVeto(dR, dZ, trackmom.Pt(), charge, mf ) < fPHOSRCut) return kTRUE;
      else                                                                    return kFALSE;
      
    }//PHOS
    else //EMCAL
    {
    if(fDebug > 1) 
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
  Double_t r2      = 0.5*  (l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
                     0.5*  (l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
                     0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  
  if(fDebug > 0) printf("AliCaloPID::TestPHOSDispersion() - PHOS SS R %f < %f?\n", TMath::Sqrt(r2), fPHOSDispersionCut);
  
  return TMath::Sqrt(r2) ; 
  
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
