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

// --- ROOT system ---
#include <TMath.h>
#include <TString.h>
#include <TList.h>

// ---- ANALYSIS system ----
#include "AliCaloPID.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVCaloCells.h"
#include "AliVTrack.h"
#include "AliAODPWG4Particle.h"
#include "AliCalorimeterUtils.h"
#include "AliFiducialCut.h" // detector enum definition
#include "AliVEvent.h"
#include "AliLog.h"

// ---- Detector ----
#include "AliEMCALPIDUtils.h"

/// \cond CLASSIMP
ClassImp(AliCaloPID) ;
/// \endcond

//________________________
/// Default constructor.
/// Initialize parameters.
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
fEMCALUseTrackPtDepMatchingCut(0), 
fEMCALFuncTrackPtDepDEta(0), fEMCALFuncTrackPtDepDPhi(0),
fEMCALFuncTrackPtDepDEtaString(""), fEMCALFuncTrackPtDepDPhiString(""), 
fEMCALFuncTrackPtDepDEtaNParam(0) , fEMCALFuncTrackPtDepDPhiNParam(0),
fTOFCut(0.), 
fPHOSDispersionCut(1000), fPHOSRCut(1000),
//Split
fUseSimpleMassCut(kFALSE),
fUseSimpleM02Cut(kFALSE),
fUseSplitAsyCut(kFALSE),
fUseSplitSSCut(kTRUE),
fSplitM02MaxCut(0),       fSplitM02MinCut(0),          fSplitMinNCells(0),
fMassEtaMin(0),           fMassEtaMax(0),
fMassPi0Min(0),           fMassPi0Max(0),
fMassPhoMin(0),           fMassPhoMax(0),
fM02MaxParamShiftNLMN(0),
fSplitWidthSigma(0),      fMassShiftHighECell(0)
{
  InitParameters();
}

//________________________________________
/// Constructor.
/// \param flux: high or low particle environment. 
/// To be used when recalculating bayesian PID. Not used currently.
/// Initialize parameters.
//________________________________________
AliCaloPID::AliCaloPID(Int_t flux) :
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
fEMCALUseTrackPtDepMatchingCut(0), 
fEMCALFuncTrackPtDepDEta(0), fEMCALFuncTrackPtDepDPhi(0),
fEMCALFuncTrackPtDepDEtaString(""), fEMCALFuncTrackPtDepDPhiString(""), 
fTOFCut(0.), 
fPHOSDispersionCut(1000), fPHOSRCut(1000),
//Split
fUseSimpleMassCut(kFALSE),
fUseSimpleM02Cut(kFALSE),
fUseSplitAsyCut(kFALSE),
fUseSplitSSCut(kTRUE),
fSplitM02MaxCut(0),       fSplitM02MinCut(0),          fSplitMinNCells(0),
fMassEtaMin(0),           fMassEtaMax(0),
fMassPi0Min(0),           fMassPi0Max(0),
fMassPhoMin(0),           fMassPhoMax(0),
fM02MaxParamShiftNLMN(0),
fSplitWidthSigma(0),      fMassShiftHighECell(0)
{
  InitParameters();
}

//_______________________________________________
/// Constructor.
/// \param emcalpid: pointer to EMCal class to recalculate PID weights. 
/// To be used when recalculating bayesian PID and need different parameters. Not used currently.
/// Initialize parameters.
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
fEMCALUseTrackPtDepMatchingCut(0), 
fEMCALFuncTrackPtDepDEta(0), fEMCALFuncTrackPtDepDPhi(0),
fEMCALFuncTrackPtDepDEtaString(""), fEMCALFuncTrackPtDepDPhiString(""), 
fTOFCut(0.), 
fPHOSDispersionCut(1000),    fPHOSRCut(1000),
//Split
fUseSimpleMassCut(kFALSE),
fUseSimpleM02Cut(kFALSE),
fUseSplitAsyCut(kFALSE),
fUseSplitSSCut(kTRUE),
fSplitM02MaxCut(0),       fSplitM02MinCut(0),          fSplitMinNCells(0),
fMassEtaMin(0),           fMassEtaMax(0),
fMassPi0Min(0),           fMassPi0Max(0),
fMassPhoMin(0),           fMassPhoMax(0),
fM02MaxParamShiftNLMN(0),
fSplitWidthSigma(0),      fMassShiftHighECell(0)
{
  InitParameters();
}

//_______________________
// Destructor.
//_______________________
AliCaloPID::~AliCaloPID() 
{  
  delete fPHOSPhotonWeightFormula ;
  delete fPHOSPi0WeightFormula ;
  delete fEMCALPIDUtils ;
  delete fEMCALFuncTrackPtDepDEta;
  delete fEMCALFuncTrackPtDepDPhi;
  delete [] fEMCALFuncTrackPtDepDEtaParam;
  delete [] fEMCALFuncTrackPtDepDPhiParam;
}

//_______________________________
// Initialize the parameters of the PID.
//_______________________________
void AliCaloPID::InitParameters()
{  
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
  
  if(fRecalculateBayesian)
  {
    if(fParticleFlux == kLow)
    {
      AliInfo("SetLOWFluxParam");
      fEMCALPIDUtils->SetLowFluxParam() ;
    }
    else if (fParticleFlux == kHigh)
    {
      AliInfo("SetHighFluxParam");
      fEMCALPIDUtils->SetHighFluxParam() ;
    }
  }
  
  //PID recalculation, not bayesian
  
  //EMCAL  
  fEMCALL0CutMax = 0.3 ;
  fEMCALL0CutMin = 0.01;
  
  // Fix Track Matching
  fEMCALDPhiCut  = 0.05; // Same cut as in AliEMCALRecoUtils
  fEMCALDEtaCut  = 0.025;// Same cut as in AliEMCALRecoUtils

  // Pt dependent track matching
  // In case we change the default setting to true
  if(fEMCALUseTrackPtDepMatchingCut) 
    InitParamTrackMatchPtDependent();
  
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
  
  fMassPi0Param[0][0] = 0     ; // Constant term on mass dependence
  fMassPi0Param[0][1] = 0     ; // slope term on mass dependence
  fMassPi0Param[0][2] = 0     ; // E function change
  fMassPi0Param[0][3] = 0.044 ; // constant term on mass dependence
  fMassPi0Param[0][4] = 0.0049; // slope term on mass dependence
  fMassPi0Param[0][5] = 0.070 ; // Absolute low mass cut
  
  fMassPi0Param[1][0] = 0.115 ; // Constant term below 21 GeV
  fMassPi0Param[1][1] = 0.00096; // slope term below 21 GeV
  fMassPi0Param[1][2] = 21    ; // E function change
  fMassPi0Param[1][3] = 0.10  ; // constant term on mass dependence
  fMassPi0Param[1][4] = 0.0017; // slope term on mass dependence
  fMassPi0Param[1][5] = 0.070 ; // Absolute low mass cut
  
  fWidthPi0Param[0][0] = 0.012 ; // Constant term on width dependence
  fWidthPi0Param[0][1] = 0.0   ; // Slope term on width dependence
  fWidthPi0Param[0][2] = 19    ; // E function change
  fWidthPi0Param[0][3] = 0.0012; // Constant term on width dependence
  fWidthPi0Param[0][4] = 0.0006; // Slope term on width dependence
  fWidthPi0Param[0][5] = 0.0   ; // xx term

  fWidthPi0Param[1][0] = 0.009 ; // Constant term on width dependence
  fWidthPi0Param[1][1] = 0.000 ; // Slope term on width dependence
  fWidthPi0Param[1][2] = 10    ; // E function change
  fWidthPi0Param[1][3] = 0.0023 ; // Constant term on width dependence
  fWidthPi0Param[1][4] = 0.00067; // Slope term on width dependence
  fWidthPi0Param[1][5] = 0.000 ;// xx term

  fMassShiftHighECell = 0; // Shift of cuts in case of higher energy threshold in cells, 5 MeV when Ecell>150 MeV
  
  //TF1 *lM02MinNLM1 = new TF1("M02MinNLM1","exp(2.135-0.245*x)",6,13.6);
  fM02MinParam[0][0] = 2.135  ; 
  fM02MinParam[0][1] =-0.245  ;
  fM02MinParam[0][2] = 0.0    ; 
  fM02MinParam[0][3] = 0.0    ;
  fM02MinParam[0][4] = 0.0    ;

  // Same as NLM=1 for NLM=2
  fM02MinParam[1][0] = 2.135  ;
  fM02MinParam[1][1] =-0.245  ;
  fM02MinParam[1][2] = 0.0    ;
  fM02MinParam[1][3] = 0.0    ;
  fM02MinParam[1][4] = 0.0    ;

  //TF1 *lM02MaxNLM1 = new TF1("M02MaxNLM1","exp(0.0662-0.0201*x)-0.0955+0.00186*x[0]+9.91/x[0]",6,100);
  fM02MaxParam[0][0] = 0.0662 ;
  fM02MaxParam[0][1] =-0.0201 ;
  fM02MaxParam[0][2] =-0.0955 ;
  fM02MaxParam[0][3] = 0.00186;
  fM02MaxParam[0][4] = 9.91   ;
  
  //TF1 *lM02MaxNLM2 = new TF1("M02MaxNLM2","exp(0.353-0.0264*x)-0.524+0.00559*x[0]+21.9/x[0]",6,100);
  fM02MaxParam[1][0] = 0.353  ;  
  fM02MaxParam[1][1] =-0.0264 ;  
  fM02MaxParam[1][2] =-0.524  ; 
  fM02MaxParam[1][3] = 0.00559;
  fM02MaxParam[1][4] = 21.9   ;
  
  fM02MaxParamShiftNLMN = 0.75;
  
  //TF1 *lAsyNLM1 = new TF1("lAsyNLM1","0.96-879/(x*x*x)",5,100);
  fAsyMinParam[0][0] = 0.96 ;
  fAsyMinParam[0][1] = 0    ;
  fAsyMinParam[0][2] =-879  ;
  fAsyMinParam[0][3] = 0.96 ; // Absolute max

  //TF1 *lAsyNLM2 = new TF1("lAsyNLM2","0.95+0.0015*x-233/(x*x*x)",5,100);
  fAsyMinParam[1][0] = 0.95  ;
  fAsyMinParam[1][1] = 0.0015;
  fAsyMinParam[1][2] =-233   ;
  fAsyMinParam[1][3] = 1.0   ; // Absolute max

  fSplitEFracMin[0]   = 0.0 ; // 0.96
  fSplitEFracMin[1]   = 0.0 ; // 0.96
  fSplitEFracMin[2]   = 0.0 ; // 0.7

  fSubClusterEMin[0]  = 0.0; // 3 GeV
  fSubClusterEMin[1]  = 0.0; // 1 GeV
  fSubClusterEMin[2]  = 0.0; // 1 GeV
  
  fSplitWidthSigma = 3. ;
}


//_______________________________
/// Initialize the default parameters of the pT dep track matching.
///
/// Borrowed from PWGGA/GammaCong/AliCaloPhotonCuts::SetTrackMatchingCut() case 7
/// Used in neutral mesons analysis at 2.76 and 8 TeV (F. Bock and D. Mulheim)
///
/// Called in InitParameters() and SwitchOnEMCTrackPtDepReaMatching()
//_______________________________
void AliCaloPID::InitParamTrackMatchPtDependent()
{
  fEMCALFuncTrackPtDepDEtaString = "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])" ;
  fEMCALFuncTrackPtDepDPhiString = "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])" ;
  
  fEMCALFuncTrackPtDepDEtaNParam = 3;
  fEMCALFuncTrackPtDepDPhiNParam = 3;
  
  fEMCALFuncTrackPtDepDEtaParam = new Float_t[fEMCALFuncTrackPtDepDEtaNParam];
  fEMCALFuncTrackPtDepDPhiParam = new Float_t[fEMCALFuncTrackPtDepDPhiNParam];
  
  fEMCALFuncTrackPtDepDEtaParam[0] = 0.04;
  fEMCALFuncTrackPtDepDPhiParam[0] = 0.09;
  fEMCALFuncTrackPtDepDEtaParam[1] = 0.010;
  fEMCALFuncTrackPtDepDPhiParam[1] = 0.015;
  fEMCALFuncTrackPtDepDEtaParam[2] = 2.5;
  fEMCALFuncTrackPtDepDPhiParam[2] = 2.;
}

//_________________________________________________________________________________________
/// Select the appropriate asymmetry range for pi0 selection in splitting method.
/// \param energy: of the cluster.
/// \param asy: asymmetry after splitting.
/// \param nlm: number of local maxima of the cluster.
/// \return kTRUE if asy is within the defined range.
//_________________________________________________________________________________________
Bool_t AliCaloPID::IsInPi0SplitAsymmetryRange(Float_t energy, Float_t asy, Int_t nlm) const
{
  if(!fUseSplitAsyCut) return kTRUE ;
  
  Float_t abasy = TMath::Abs(asy);

  Int_t inlm = nlm-1;
  if(nlm > 2) inlm=1; // only 2 cases defined nlm=1 and nlm>=2
  
  // Get the parametrized min cut of asymmetry for NLM=2 up to 11 GeV

  Float_t cut = fAsyMinParam[inlm][0] + fAsyMinParam[inlm][1]*energy + fAsyMinParam[inlm][2]/energy/energy/energy ;
  
  // In any case and beyond validity energy range of the function,
  // the parameter cannot be smaller than 1
  if( cut > fAsyMinParam[inlm][3] ) cut = fAsyMinParam[inlm][3];
  
  //printf("energy %2.2f - nlm: %d (%d)- p0 %f, p1 %f, p2 %f, p3 %f ; cut: %2.2f\n",energy,nlm,inlm,
  //       fAsyMinParam[inlm][0],fAsyMinParam[inlm][1],fAsyMinParam[inlm][2],fAsyMinParam[inlm][3],cut);
  
  if(abasy < cut) return kTRUE;
  else            return kFALSE;
}

//______________________________________________________________________________________
/// Select the appropriate mass range for pi0 selection in splitting method
/// \param energy: of the cluster.
/// \param mass: mass after splitting.
/// \param nlm: number of local maxima of the cluster.
/// \return kTRUE if mass is within the defined range.
//______________________________________________________________________________________
Bool_t AliCaloPID::IsInPi0SplitMassRange(Float_t energy, Float_t mass, Int_t nlm) const
{
  if(fUseSimpleMassCut)
  {
    if(mass < fMassPi0Max && mass > fMassPi0Min) return kTRUE;
    else                                         return kFALSE;
  }
  
  // Get the selected mean value as reference for the mass
  Int_t inlm = nlm-1;
  if(nlm > 2) inlm=1; // only 2 cases defined nlm=1 and nlm>=2
    
  Float_t meanMass = energy * fMassPi0Param[inlm][1] + fMassPi0Param[inlm][0];
  if(energy > fMassPi0Param[inlm][2]) meanMass = energy * fMassPi0Param[inlm][4] + fMassPi0Param[inlm][3];
  
  // In case of higher energy cell cut than 50 MeV, smaller mean mass 0-5 MeV, not really necessary
  meanMass -= fMassShiftHighECell;
  
  // Get the parametrized width of the mass
  Float_t width   = 0.009;
  if     (energy > 8 && energy < fWidthPi0Param[inlm][2])
    width = energy * fWidthPi0Param[inlm][1] + fWidthPi0Param[inlm][0];
  else if(              energy > fWidthPi0Param[inlm][2])
    width = energy * energy * fWidthPi0Param[inlm][5] + energy * fWidthPi0Param[inlm][4] + fWidthPi0Param[inlm][3];

  // Calculate the 2 sigma cut
  Float_t minMass = meanMass-fSplitWidthSigma*width;
  Float_t maxMass = meanMass+fSplitWidthSigma*width;

  // In case of low energy, hard cut to avoid conversions
  if(energy < 10  && minMass < fMassPi0Param[inlm][5] ) minMass = fMassPi0Param[inlm][5];
  
  //printf("E %2.2f, mass %1.1f, nlm %d: sigma %1.1f width %3.1f, mean Mass %3.0f, minMass %3.0f, maxMass %3.0f\n ",
  //       energy,mass *1000, inlm, fSplitWidthSigma, width*1000, meanMass*1000,minMass*1000,maxMass*1000);
  
  if(mass < maxMass && mass > minMass) return kTRUE;
  else                                 return kFALSE;
}

//________________________________________________
/// Select the appropriate shower shape range, simple fix range, not energy dependent.
/// \param m02: shower shape main axis of the cluster.
/// \return kTRUE if m02 is within the defined range
//________________________________________________
Bool_t AliCaloPID::IsInM02Range(Float_t m02) const
{    
  Float_t minCut = fSplitM02MinCut;
  Float_t maxCut = fSplitM02MaxCut;

  if(m02 < maxCut && m02 > minCut) return kTRUE;
  else                             return kFALSE;
}

//_______________________________________________________________________________
/// Select the appropriate shower shape range in splitting method for pi0
/// \param energy: of the cluster.
/// \param m02: shower shape main axis of the cluster.
/// \param nlm: number of local maxima of the cluster.
/// \return kTRUE if m02 is within the defined range.
//_______________________________________________________________________________
Bool_t AliCaloPID::IsInPi0M02Range(Float_t energy, Float_t m02,  Int_t nlm) const
{
  if(!fUseSplitSSCut) return kTRUE ;

  //First check the absolute minimum and maximum
  if(!IsInM02Range(m02)) return kFALSE ;
  
  //If requested, check the E dependent cuts
  else if(!fUseSimpleM02Cut)
  {
    Int_t inlm = nlm-1;
    if(nlm > 2) inlm=1; // only 2 cases defined nlm=1 and nlm>=2

    Float_t minCut = fSplitM02MinCut;
    Float_t maxCut = fSplitM02MaxCut;
    
    //e^{a+bx} + c + dx + e/x
    if(energy > 1) minCut = TMath::Exp( fM02MinParam[inlm][0] + fM02MinParam[inlm][1]*energy ) +
                            fM02MinParam[inlm][2] + fM02MinParam[inlm][3]*energy + fM02MinParam[inlm][4]/energy;
    
    if(energy > 1) maxCut = TMath::Exp( fM02MaxParam[inlm][0] + fM02MaxParam[inlm][1]*energy ) +
                            fM02MaxParam[inlm][2] + fM02MaxParam[inlm][3]*energy + fM02MaxParam[inlm][4]/energy;
    
    // In any case and beyond validity energy range of the function,
    // the parameter cannot be smaller than 0.3 or larger than 4-5
    if( minCut < fSplitM02MinCut) minCut = fSplitM02MinCut;
    if( maxCut > fSplitM02MaxCut) maxCut = fSplitM02MaxCut;
    if( nlm > 2 ) maxCut+=fM02MaxParamShiftNLMN;
    
    //if(energy > 7) printf("\t \t E %2.2f, nlm %d, m02 %2.2f, minM02 %2.2f, maxM02 %2.2f\n",energy, nlm, m02,minCut,maxCut);
    
    if(m02 < maxCut && m02 > minCut) return kTRUE;
    else                             return kFALSE;
    
  }
  
  else return kTRUE;
}


//______________________________________________________________________________
// Select the appropriate shower shape range in splitting method to select eta's
// Use same parametrization as pi0, just shift the distributions (to be tuned)
/// \param energy: of the cluster.
/// \param m02: shower shape main axis of the cluster.
/// \param nlm: number of local maxima of the cluster.
/// \return kTRUE if m02 is within the defined range.
//______________________________________________________________________________
Bool_t AliCaloPID::IsInEtaM02Range(Float_t energy, Float_t m02, Int_t nlm) const
{  
  if(!fUseSplitSSCut) return kTRUE ;
  
  //First check the absolute minimum and maximum
  if(!IsInM02Range(m02)) return kFALSE ;
  
  //DO NOT USE, study parametrization
  
  //If requested, check the E dependent cuts
  else if(!fUseSimpleM02Cut)
  {
    Int_t inlm = nlm-1;
    if(nlm > 2) inlm=1; // only 2 cases defined nlm=1 and nlm>=2
    
    Float_t minCut = fSplitM02MinCut;
    Float_t maxCut = fSplitM02MaxCut;
    
    Float_t shiftE = energy-20; // to be tuned
    if(nlm==1) shiftE=energy-28;
    
    //e^{a+bx} + c + dx + e/x
    if(shiftE > 1) minCut = TMath::Exp( fM02MinParam[inlm][0] + fM02MinParam[inlm][1]*shiftE ) +
                  fM02MinParam[inlm][2] + fM02MinParam[inlm][3]*shiftE + fM02MinParam[inlm][4]/shiftE;
    
    // In any case the parameter cannot be smaller than 0.3
    if( minCut < fSplitM02MinCut) minCut = fSplitM02MinCut;
    
    shiftE = energy+20; // to be tuned
    
    if(shiftE > 1)  maxCut = 1 + TMath::Exp( fM02MaxParam[inlm][0] + fM02MaxParam[inlm][1]*shiftE ) +
                             fM02MaxParam[inlm][2] + fM02MaxParam[inlm][3]*shiftE + fM02MaxParam[inlm][4]/shiftE;
    
    // In any case the parameter cannot be smaller than 4-5
    if( maxCut > fSplitM02MaxCut) maxCut = fSplitM02MaxCut;
    if( nlm > 2 ) maxCut+=fM02MaxParamShiftNLMN;
    
    //if(energy>6)printf("\t \t E %2.2f, nlm %d, m02 %2.2f, minM02 %2.2f, maxM02 %2.2f\n",energy, nlm, m02,minCut,maxCut);
    
    if(m02 < maxCut && m02 > minCut) return kTRUE;
    else                             return kFALSE;
    
  }
  
  else return kTRUE;
}

//______________________________________________________________________________
/// Select the appropriate shower shape range in splitting method for converted photons
/// Just minimum limit for pi0s is max for conversion.
/// \param energy: of the cluster.
/// \param m02: shower shape main axis of the cluster.
/// \param nlm: number of local maxima of the cluster.
/// \return kTRUE if m02 is within the defined range.
//______________________________________________________________________________
Bool_t AliCaloPID::IsInConM02Range(Float_t energy, Float_t m02, Int_t nlm) const
{
  if(!fUseSplitSSCut) return kTRUE ;
  
  Float_t minCut = 0.1;
  Float_t maxCut = fSplitM02MinCut;
  
  if(!fUseSimpleM02Cut)
  {
    Int_t inlm = nlm-1;
    if(nlm > 2) inlm=1; // only 2 cases defined nlm=1 and nlm>=2
    
    //e^{a+bx} + c + dx + e/x
    if(energy > 1) maxCut = TMath::Exp( fM02MinParam[inlm][0] + fM02MinParam[inlm][1]*energy ) +
                            fM02MinParam[inlm][2] + fM02MinParam[inlm][3]*energy + fM02MinParam[inlm][4]/energy;
    
    if( maxCut < fSplitM02MinCut) maxCut = fSplitM02MinCut;
  }
  
  if(m02 < maxCut && m02 > minCut) return kTRUE;
  else                             return kFALSE;
  
}

//______________________________________________
/// \return pointer to AliEMCALPIDUtils, create it if needed.
//______________________________________________
AliEMCALPIDUtils *AliCaloPID::GetEMCALPIDUtils() 
{  
  if(!fEMCALPIDUtils) fEMCALPIDUtils = new AliEMCALPIDUtils ; 
 
  return fEMCALPIDUtils ; 
}


//________________________________________________________________ 
/// \return PDG number corresponding to the likely ID of the cluster
/// \param cluster: input cluster pointer
//________________________________________________________________ 
Int_t AliCaloPID::GetIdentifiedParticleType(AliVCluster * cluster)
{  
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
  
  AliDebug(1,Form("EMCAL %d?, E %3.2f, l0 %3.2f, l1 %3.2f, disp %3.2f, tof %1.11f, distCPV %3.2f, distToBC %1.1f, NMax %d",
                  cluster->IsEMCAL(),energy,lambda0,cluster->GetM20(),cluster->GetDispersion(),cluster->GetTOF(),
                  cluster->GetEmcCpvDistance(), cluster->GetDistanceToBadChannel(),cluster->GetNExMax()));
  
  if(cluster->IsEMCAL())
  {
    AliDebug(1,Form("EMCAL SS %f <%f < %f?",fEMCALL0CutMin, lambda0, fEMCALL0CutMax));
    
    if(lambda0 < fEMCALL0CutMax && lambda0 > fEMCALL0CutMin) return kPhoton ;
    else                                                     return kNeutralUnknown ; 
  }    // EMCAL
  else // PHOS
  {    
    if(TestPHOSDispersion(energy,lambda0,lambda1) < fPHOSDispersionCut) return kPhoton;
    else                                                                return kNeutralUnknown;
  }
}

//_________________________________________________________________________________________________________
/// \return most probable identity (PDG) of the particle after bayesian weights calculated in reconstruction.
/// \param isEMCAL: which calo.
/// \param pid: array with bayesian weights
/// \param energy: of the cluster
//_________________________________________________________________________________________________________
Int_t AliCaloPID::GetIdentifiedParticleTypeFromBayesWeights(Bool_t isEMCAL, Double_t * pid, Float_t energy)
{  
  if(!pid)
  { 
    AliFatal("pid pointer not initialized!!!");
    return kNeutralUnknown; // not needed, added to content coverity
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
  
  AliDebug(1,Form("Final Pdg: %d, cluster energy %2.2f", pdg,energy));

  return pdg ;
  
}

//____________________________________________________________________________________________________________
/// Split the cluster in 2, do invariant mass, get the mass and decide 
/// if this is a photon, pi0, eta, .... Only implemented for EMCal, for now.
/// \param cluster: input cluster pointer.
/// \param cells: list of calorimeter cells.
/// \param caloutils: pointer to AliCalorimeterUtils, the split is implemented in this class.
/// \param vertex: Event vertex array.
/// \param nMax: number of local maxima of the input cluster, output.
/// \param mass: splitted cluster mass, output.  
/// \param angle: opening angle of sub-clusters, output.
/// \param l1: sub-cluster kinematics, output.
/// \param l2: sub-cluster kinematics, output.
/// \param absId1: sub-cluster main cell index, output.
/// \param absId2: sub-cluster main cell index, output.
/// \param distbad1: sub-cluster distance to bad channel, output. 
/// \param distbad2: sub-cluster distance to bad channel, output.
/// \param fidcut1: sub-cluster close to border, output.  
/// \param fidcut2: sub-cluster close to border, output. 
/// \return kTRUE if cluster type: kPi0, kEta, kPhoton, kNeutralUnknown
//____________________________________________________________________________________________________________
Int_t AliCaloPID::GetIdentifiedParticleTypeFromClusterSplitting(AliVCluster* cluster, 
                                                                AliVCaloCells* cells,
                                                                AliCalorimeterUtils * caloutils,
                                                                Double_t   vertex[3],
                                                                Int_t    & nMax,
                                                                Double_t & mass,   Double_t & angle,
                                                                TLorentzVector & l1, TLorentzVector & l2,
                                                                Int_t   & absId1,   Int_t   & absId2,
                                                                Float_t & distbad1, Float_t & distbad2,
                                                                Bool_t  & fidcut1,  Bool_t  & fidcut2  ) const
{  
  Float_t eClus  = cluster->E();
  Float_t m02    = cluster->GetM02();
  const Int_t nc = cluster->GetNCells();
  Int_t   absIdList[nc]; 
  Float_t maxEList [nc];
  
  mass  = -1.;
  angle = -1.;
  
  //If too low number of cells, skip it
  if ( nc < fSplitMinNCells)  return kNeutralUnknown ; 
  
  AliDebug(2,"\t pass nCells cut");
  
  // Get Number of local maxima
  nMax  = caloutils->GetNumberOfLocalMaxima(cluster, cells, absIdList, maxEList) ; 
  
  AliDebug(1,Form("Cluster : E %1.1f, M02 %1.2f, NLM %d, N Cells %d",eClus,m02,nMax,nc));

  //---------------------------------------------------------------------
  // Get the 2 max indeces and do inv mass
  //---------------------------------------------------------------------
  
  Int_t  calorimeter = AliFiducialCut::kEMCAL;
  if(cluster->IsPHOS()) calorimeter = AliFiducialCut::kPHOS;

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
    AliDebug(1,Form("Bad index for local maxima : N max %d, i1 %d, i2 %d, cluster E %2.2f, ncells %d, m02 %2.2f",
                    nMax,absId1,absId2,eClus,nc,m02));
    return kNeutralUnknown ;
  }
  
  //---------------------------------------------------------------------
  // Split the cluster energy in 2, around the highest 2 local maxima
  //---------------------------------------------------------------------  
  
  AliAODCaloCluster cluster1(0, 0,NULL,0.,NULL,NULL,1,0);
  AliAODCaloCluster cluster2(1, 0,NULL,0.,NULL,NULL,1,0);
  
  caloutils->SplitEnergy(absId1,absId2,cluster, cells, &cluster1, &cluster2,nMax); /*absIdList, maxEList,*/
  
  fidcut1 = caloutils->GetEMCALRecoUtils()->CheckCellFiducialRegion(caloutils->GetEMCALGeometry(), &cluster1,cells);
  fidcut2 = caloutils->GetEMCALRecoUtils()->CheckCellFiducialRegion(caloutils->GetEMCALGeometry(), &cluster2,cells);

  caloutils->GetEMCALRecoUtils()->RecalculateClusterDistanceToBadChannel(caloutils->GetEMCALGeometry(),cells,&cluster1);
  caloutils->GetEMCALRecoUtils()->RecalculateClusterDistanceToBadChannel(caloutils->GetEMCALGeometry(),cells,&cluster2);

  distbad1 = cluster1.GetDistanceToBadChannel();
  distbad2 = cluster2.GetDistanceToBadChannel();
//  if(!fidcut2 || !fidcut1 || distbad1 < 2 || distbad2 < 2)
//    printf("*** Dist to bad channel cl %f, cl1 %f, cl2 %f; fid cut cl %d, cl1 %d, cl2 %d \n",
//           cluster->GetDistanceToBadChannel(),distbad1,distbad2,
//           caloutils->GetEMCALRecoUtils()->CheckCellFiducialRegion(caloutils->GetEMCALGeometry(), cluster,cells),fidcut1,fidcut2);
  
  cluster1.GetMomentum(l1,vertex);
  cluster2.GetMomentum(l2,vertex);
  
  mass  = (l1+l2).M();
  angle = l2.Angle(l1.Vect());
  Float_t e1 = cluster1.E();
  Float_t e2 = cluster2.E();
  
  // Consider clusters with splitted energy not too different to original cluster energy
  Float_t splitFracCut = 0;
  if(nMax < 3)  splitFracCut = fSplitEFracMin[nMax-1];
  else          splitFracCut = fSplitEFracMin[2];
  if((e1+e2)/eClus < splitFracCut) return kNeutralUnknown ;

  AliDebug(2,"\t pass Split E frac cut");
  
  // Consider sub-clusters with minimum energy
  Float_t minECut = fSubClusterEMin[2];
  if     (nMax == 2)  minECut = fSubClusterEMin[1];
  else if(nMax == 1)  minECut = fSubClusterEMin[0];
  if(e1 < minECut || e2 < minECut)
  {
    //printf("Reject: e1 %2.1f, e2 %2.1f, cut %2.1f\n",e1,e2,minECut);
    return kNeutralUnknown ;
  }

  AliDebug(2,"\t pass min sub-cluster E cut");
  
  // Asymmetry of cluster
  Float_t asy =-10;
  if(e1+e2 > 0) asy = (e1-e2) / (e1+e2);

  if( !IsInPi0SplitAsymmetryRange(eClus,asy,nMax) ) return kNeutralUnknown ;
  
  
  AliDebug(2,"\t pass asymmetry cut");
  
  Bool_t pi0OK = kFALSE;
  Bool_t etaOK = kFALSE;
  Bool_t conOK = kFALSE;
  
  //If too small or big M02, skip it
  if     (IsInPi0M02Range(eClus,m02,nMax))  pi0OK = kTRUE;
  else if(IsInEtaM02Range(eClus,m02,nMax))  etaOK = kTRUE;
  else if(IsInConM02Range(eClus,m02,nMax))  conOK = kTRUE;
  
  Float_t energy = eClus;
  if(nMax > 2) energy = e1+e2; // In case of NLM>2 use mass cut for NLM=2 but for the split sum not the cluster energy that is not the pi0 E.
  
  // Check the mass, and set an ID to the splitted cluster
  if     ( conOK && mass < fMassPhoMax && mass > fMassPhoMin ) { AliDebug(2,"\t Split Conv"); return kPhoton ; }
  else if( etaOK && mass < fMassEtaMax && mass > fMassEtaMin ) { AliDebug(2,"\t Split Eta" ); return kEta    ; }
  else if( pi0OK && IsInPi0SplitMassRange(energy,mass,nMax)  ) { AliDebug(2,"\t Split Pi0" ); return kPi0    ; }
  else                                                                                        return kNeutralUnknown ;
  
}

//_________________________________________
/// Put data member values in string to keep in output container.
//_________________________________________
TString  AliCaloPID::GetPIDParametersList()  
{  
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  snprintf(onePar,buffersize,"--- AliCaloPID ---") ;
  parList+=onePar ;	
  if(fUseBayesianWeights)
  {
    snprintf(onePar,buffersize,"fEMCALPhotonWeight =%2.2f (EMCAL bayesian weight for photons)",fEMCALPhotonWeight) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fEMCALPi0Weight =%2.2f (EMCAL bayesian weight for pi0)",fEMCALPi0Weight) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fEMCALElectronWeight =%2.2f(EMCAL bayesian weight for electrons)",fEMCALElectronWeight) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fEMCALChargeWeight =%2.2f (EMCAL bayesian weight for charged hadrons)",fEMCALChargeWeight) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fEMCALNeutralWeight =%2.2f (EMCAL bayesian weight for neutral hadrons)",fEMCALNeutralWeight) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fPHOSPhotonWeight =%2.2f (PHOS bayesian weight for photons)",fPHOSPhotonWeight) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fPHOSPi0Weight =%2.2f (PHOS bayesian weight for pi0)",fPHOSPi0Weight) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fPHOSElectronWeight =%2.2f(PHOS bayesian weight for electrons)",fPHOSElectronWeight) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fPHOSChargeWeight =%2.2f (PHOS bayesian weight for charged hadrons)",fPHOSChargeWeight) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fPHOSNeutralWeight =%2.2f (PHOS bayesian weight for neutral hadrons)",fPHOSNeutralWeight) ;
    parList+=onePar ;
    
    if(fPHOSWeightFormula)
    {
      snprintf(onePar,buffersize,"PHOS Photon Weight Formula: %s",fPHOSPhotonWeightFormulaExpression.Data() ) ;
      parList+=onePar;
      snprintf(onePar,buffersize,"PHOS Pi0    Weight Formula: %s",fPHOSPi0WeightFormulaExpression.Data()    ) ;
      parList+=onePar;	  
    }
  }
  else
  {
    snprintf(onePar,buffersize,"EMCAL: fEMCALL0CutMin =%2.2f, fEMCALL0CutMax =%2.2f  (Cut on Shower Shape)",fEMCALL0CutMin, fEMCALL0CutMax) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"EMCAL: fEMCALDEtaCut =%2.2f, fEMCALDPhiCut =%2.2f  (Cut on track matching)",fEMCALDEtaCut, fEMCALDPhiCut) ;
    parList+=onePar ;
    snprintf(onePar,buffersize,"fTOFCut  =%e (Cut on TOF, used in PID evaluation)",fTOFCut) ;
    parList+=onePar ;	
    snprintf(onePar,buffersize,"fPHOSRCut =%2.2f, fPHOSDispersionCut =%2.2f  (Cut on Shower Shape and CPV)",fPHOSRCut,fPHOSDispersionCut) ;
    parList+=onePar ;
    
  }
  
  if(fUseSimpleM02Cut)
  {
    snprintf(onePar,buffersize,"%2.2f< M02 < %2.2f",    fSplitM02MinCut, fSplitM02MaxCut) ;
    parList+=onePar ;
  }
  snprintf(onePar,buffersize,"fMinNCells =%d",        fSplitMinNCells) ;
  parList+=onePar ;
  if(fUseSimpleMassCut)
  {
    snprintf(onePar,buffersize,"pi0 : %2.1f < m <%2.1f", fMassPi0Min,fMassPi0Max);
    parList+=onePar ;
  }
  snprintf(onePar,buffersize,"eta : %2.1f < m <%2.1f", fMassEtaMin,fMassEtaMax);
  parList+=onePar ;
  snprintf(onePar,buffersize,"conv: %2.1f < m <%2.1f", fMassPhoMin,fMassPhoMax);
  parList+=onePar ;
  
  
  return parList;
}

//________________________________________________
/// Print some relevant parameters set for the analysis
//________________________________________________
void AliCaloPID::Print(const Option_t * opt) const
{
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
  
  printf("Min. N Cells =%d \n",         fSplitMinNCells) ;
  if(fUseSimpleM02Cut) printf("%2.2f < lambda_0^2 <%2.2f \n",fSplitM02MinCut,fSplitM02MaxCut);
  if(fUseSimpleMassCut)printf("pi0 : %2.2f<m<%2.2f \n",      fMassPi0Min,fMassPi0Max);
  printf("eta : %2.2f<m<%2.2f \n",      fMassEtaMin,fMassEtaMax);
  printf("phot: %2.2f<m<%2.2f \n",      fMassPhoMin,fMassPhoMax);
  
  printf(" \n");
} 

//_________________________________________________________________
// Print PID of cluster, (AliVCluster*)cluster->GetPID()
//_________________________________________________________________
void AliCaloPID::PrintClusterPIDWeights(const Double_t * pid) const
{  
  printf("AliCaloPID::PrintClusterPIDWeights() \n \t ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, \n \t \
         pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
         pid[AliVCluster::kPhoton],    pid[AliVCluster::kPi0],
         pid[AliVCluster::kElectron],  pid[AliVCluster::kEleCon],
         pid[AliVCluster::kPion],      pid[AliVCluster::kKaon], 
         pid[AliVCluster::kProton],
         pid[AliVCluster::kNeutron],   pid[AliVCluster::kKaon0]);
}

//___________________________________________________________________________
/// Set Bits for PID selection
//___________________________________________________________________________
void AliCaloPID::SetPIDBits(AliVCluster * cluster, 
                            AliAODPWG4Particle * ph, AliCalorimeterUtils* cu, 
                            AliVEvent* event) 
{  
  // Dispersion/lambdas
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
 
  AliDebug(1,Form("TOF %e, Lambda0 %2.2f, Lambda1 %2.2f",tof , l0, l1));
  AliDebug(1,Form("pdg %d, bits: TOF %d, Dispersion %d, Charge %d",
           ph->GetIdentifiedParticleType(), ph->GetTOFBit() , ph->GetDispBit() , ph->GetChargedBit()));
}

//_________________________________________________________
/// \return function with eta residual matching cut.
///
/// Consider as charged clusters those with a eta residual larger 
/// than the value of the function for a given track pT.
///
/// Init the function here the first time it is called.
//_________________________________________________________
TF1 * AliCaloPID::GetEMCALFuncTrackPtDepDEta()      
{ 
  if( !fEMCALFuncTrackPtDepDEta ) 
  {
    fEMCALFuncTrackPtDepDEta = new TF1("emc_dEta",fEMCALFuncTrackPtDepDEtaString) ;
    
    for(Int_t iparam = 0; iparam < fEMCALFuncTrackPtDepDEtaNParam; iparam++)
      fEMCALFuncTrackPtDepDEta->SetParameter(iparam,fEMCALFuncTrackPtDepDEtaParam[iparam]);
  }
  
  return fEMCALFuncTrackPtDepDEta ; 
} 

//_________________________________________________________
/// \return function with phi residual matching cut.
///
/// Consider as charged clusters those with a phi residual larger 
/// than the value of the function for a given track pT.
///
/// Init the function here the first time it is called.
//_________________________________________________________
TF1 * AliCaloPID::GetEMCALFuncTrackPtDepDPhi()      
{ 
  if ( !fEMCALFuncTrackPtDepDPhi ) 
  {
    fEMCALFuncTrackPtDepDPhi = new TF1("emc_dPhi",fEMCALFuncTrackPtDepDPhiString) ;
    
    for(Int_t iparam = 0; iparam < fEMCALFuncTrackPtDepDPhiNParam; iparam++)
      fEMCALFuncTrackPtDepDPhi->SetParameter(iparam,fEMCALFuncTrackPtDepDPhiParam[iparam]);
  }
  
  return fEMCALFuncTrackPtDepDPhi; 
} 

//_________________________________________________________
/// Check if there is any track attached to this cluster.
/// \param cluster: pointer to calorimeter cluster.
/// \param cu: pointer to AliCalorimeterUtils, needed if track matching is recalculated in the fly
/// \param event: AliVEvent pointer. Needed to get the tracks or the magnetic field.
/// \return kTRUE if cluster is matched by a track.
//_________________________________________________________
Bool_t AliCaloPID::IsTrackMatched(AliVCluster* cluster,
                                  AliCalorimeterUtils * cu,
                                  AliVEvent* event) 
{  
  Int_t nMatches = cluster->GetNTracksMatched();
  AliVTrack * track = 0;
  
  // At least one match
  //
  if(nMatches <= 0) return kFALSE;
  
  // Select the track, depending on ESD or AODs
  //
  //In case of ESDs, 
  //by default without match one entry with negative index, no match, reject.
  //
  if(!strcmp("AliESDCaloCluster",Form("%s",cluster->ClassName())))
  {
    Int_t iESDtrack = ((AliESDCaloCluster*)cluster)->GetTracksMatched()->At(0); //cluster->GetTrackMatchedIndex();
    
    if(iESDtrack >= 0) track = dynamic_cast<AliVTrack*> (event->GetTrack(iESDtrack));
    else return kFALSE;
    
    if (!track)
    {
      AliWarning(Form("Null matched track in ESD for index %d",iESDtrack));
      return kFALSE;
    }
  }    // ESDs
  else
  {    // AODs
    track = dynamic_cast<AliVTrack*> (cluster->GetTrackMatched(0));
    if (!track)
    {
      AliWarning("Null matched track in AOD!");
      return kFALSE;
    }
  }   // AODs
  
  Float_t dEta  = cluster->GetTrackDz();
  Float_t dPhi  = cluster->GetTrackDx();
  
  // Comment out, new value already set in AliCalorimeterUtils::RecalculateClusterTrackMatching()
  // when executed in the reader.
  //    // if track matching was recalculated
  //    if(cluster->IsEMCAL() && cu && cu->IsRecalculationOfClusterTrackMatchingOn())
  //    {
  //      dR = 2000., dZ = 2000.;
  //      cu->GetEMCALRecoUtils()->GetMatchedResiduals(cluster->GetID(),dZ,dR);
  //      //AliDebug(2,"Residuals, (Old, New): z (%2.4f,%2.4f), x (%2.4f,%2.4f)\n", cluster->GetTrackDz(),dZ,cluster->GetTrackDx(),dR));
  //    }
  
  if(cluster->IsPHOS())
  {
    Int_t charge = track->Charge();
    Double_t mf  = event->GetMagneticField();
    if(TestPHOSChargedVeto(dPhi, dEta, track->Pt(), charge, mf ) < fPHOSRCut) 
      return kTRUE;
    else                                                                  
      return kFALSE;
    
  }    // PHOS
  else // EMCAL
  {
    AliDebug(1,Form("EMCAL dPhi %f < %f, dEta %f < %f ",dPhi, fEMCALDPhiCut, dEta, fEMCALDEtaCut));
    
    if(!fEMCALUseTrackPtDepMatchingCut)
    {
      if(TMath::Abs(dPhi) < fEMCALDPhiCut &&
         TMath::Abs(dEta) < fEMCALDEtaCut)   return kTRUE;
      else                                   return kFALSE;
    }
    else
    {
      Float_t trackPt = track->Pt();

      Bool_t matchDEta = kFALSE;
      if( TMath::Abs(dEta) < GetEMCALFuncTrackPtDepDEta()->Eval(trackPt)) 
        matchDEta = kTRUE;
      else 
        matchDEta = kFALSE;
      
      Bool_t matchDPhi = kFALSE;
      if( TMath::Abs(dPhi) < GetEMCALFuncTrackPtDepDPhi()->Eval(trackPt)) 
        matchDPhi = kTRUE;
      else 
        matchDPhi = kFALSE;

//      printf("Cluster E %2.2f, track pT %2.2f, dEta %2.2f, dPhi %2.2f, cut eta %2.2f, cut phi %2.2f, match eta %d, match phi %d\n",
//             cluster->E(),trackPt,dEta,dPhi,
//             GetEMCALFuncTrackPtDepDEta()->Eval(trackPt), GetEMCALFuncTrackPtDepDPhi()->Eval(trackPt),
//             matchDEta, matchDPhi);
      
      if(matchDPhi && matchDEta) return kTRUE ;
      else                       return kFALSE;
      
    }
  }// EMCAL cluster
}

//___________________________________________________________________________________________________
/// Check if PHOS cluster is photon-like. Uses photon cluster parameterization in real pp data.
/// \param pt: cluster momentum
/// \param l1: first ellipse axis
/// \param l2: sencond ellipse axis
/// \return distance in sigmas. 
///
/// Recommended cut 2.5
//___________________________________________________________________________________________________
Float_t AliCaloPID::TestPHOSDispersion(Double_t pt, Double_t l1, Double_t l2) const
{
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c       =-0.35-0.550*TMath::Exp(-0.390730*pt) ;
  Double_t r2      = 0.5*  (l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma +
  0.5*  (l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
  0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  
  AliDebug(1,Form("PHOS SS R %f < %f?", TMath::Sqrt(r2), fPHOSDispersionCut));
  
  return TMath::Sqrt(r2) ;
}

//_______________________________________________________________________________________________
/// Checks distance to the closest track. Takes into account 
/// non-perpendicular incidence of tracks.
/// \param dx: track-cluster residual x.
/// \param dz: track-cluster residual z.
/// \param pt: cluster (or track?) pt.
/// \param charge: sign of charge.
/// \param mf: magnetic field sign.
/// \return distance in sigmas. Recommended cut: 2.
///
/// Requires (sign) of magnetic filed. onc can find it for example as following
///   Double_t mf=0. ;
///   AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
///   if(event)
///     mf = event->GetMagneticField(); //Positive for ++ and negative for --
//_______________________________________________________________________________________________
Float_t AliCaloPID::TestPHOSChargedVeto(Double_t dx,  Double_t dz, Double_t pt, 
                                        Int_t charge, Double_t mf) const 
{
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
  
  AliDebug(1,Form("PHOS Matching R %f < %f",TMath::Sqrt(rx*rx+rz*rz), fPHOSRCut));
  
  return TMath::Sqrt(rx*rx+rz*rz) ;
}

