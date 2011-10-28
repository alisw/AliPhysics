/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/********************************** 
 * flow analysis with Q-cumulants * 
 *                                * 
 * author: Ante Bilandzic         * 
 *        (abilandzic@gmail.com)  *
 *********************************/ 

#define AliFlowAnalysisWithQCumulants_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "TChain.h"

#include "TFile.h"
#include "TList.h"
#include "TGraph.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TProfile2D.h" 
#include "TMath.h"
#include "TArrow.h"
#include "TPaveLabel.h"
#include "TCanvas.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithQCumulants.h"
#include "TArrayD.h"
#include "TRandom.h"
#include "TF1.h"

class TH1;
class TH2;
class TGraph;
class TPave;
class TLatex;
class TMarker;
class TRandom3;
class TObjArray;
class TList;
class TCanvas;
class TSystem;
class TROOT;
class AliFlowVector;
class TVector;

//================================================================================================================

ClassImp(AliFlowAnalysisWithQCumulants)

AliFlowAnalysisWithQCumulants::AliFlowAnalysisWithQCumulants(): 
 // 0.) base:
 fHistList(NULL),
 // 1.) common:
 fBookOnlyBasicCCH(kFALSE),
 fCommonHists(NULL),
 fCommonHists2nd(NULL), 
 fCommonHists4th(NULL),
 fCommonHists6th(NULL),
 fCommonHists8th(NULL),
 fCommonHistsResults2nd(NULL),
 fCommonHistsResults4th(NULL),
 fCommonHistsResults6th(NULL),
 fCommonHistsResults8th(NULL),
 fnBinsPhi(0),
 fPhiMin(0),
 fPhiMax(0),
 fPhiBinWidth(0),
 fnBinsPt(0),
 fPtMin(0),
 fPtMax(0),
 fPtBinWidth(0),
 fnBinsEta(0),
 fEtaMin(0),
 fEtaMax(0),
 fEtaBinWidth(0),
 fCommonConstants(NULL),
 fFillMultipleControlHistograms(kFALSE),
 fHarmonic(2),
 fAnalysisLabel(NULL),
 // 2a.) particle weights:
 fWeightsList(NULL),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fUseTrackWeights(kFALSE),
 fUseParticleWeights(NULL),
 fPhiWeights(NULL),
 fPtWeights(NULL),
 fEtaWeights(NULL),
 // 2b.) event weights:
 fMultiplicityWeight(NULL),
 // 3.) integrated flow:
 fIntFlowList(NULL), 
 fIntFlowProfiles(NULL),
 fIntFlowResults(NULL),
 fIntFlowAllCorrelationsVsM(NULL),
 fIntFlowFlags(NULL),
 fApplyCorrectionForNUA(kFALSE),  
 fApplyCorrectionForNUAVsM(kFALSE),
 fnBinsMult(10000),
 fMinMult(0.),  
 fMaxMult(10000.), 
 fPropagateErrorAlsoFromNIT(kFALSE), 
 fCalculateCumulantsVsM(kFALSE),
 fCalculateAllCorrelationsVsM(kFALSE), 
 fMinimumBiasReferenceFlow(kTRUE), 
 fForgetAboutCovariances(kFALSE), 
 fStorePhiDistributionForOneEvent(kFALSE),
 fReQ(NULL),
 fImQ(NULL),
 fSpk(NULL),
 fIntFlowCorrelationsEBE(NULL),
 fIntFlowEventWeightsForCorrelationsEBE(NULL),
 fIntFlowCorrelationsAllEBE(NULL),
 fReferenceMultiplicityEBE(0.),  
 fAvMultiplicity(NULL),
 fIntFlowCorrelationsPro(NULL),
 fIntFlowSquaredCorrelationsPro(NULL),
 fIntFlowCorrelationsAllPro(NULL),
 fIntFlowExtraCorrelationsPro(NULL),
 fIntFlowProductOfCorrelationsPro(NULL),
 fIntFlowProductOfCorrectionTermsForNUAPro(NULL),
 fIntFlowCorrelationsHist(NULL),
 fIntFlowCorrelationsAllHist(NULL),
 fIntFlowCovariances(NULL),
 fIntFlowSumOfProductOfEventWeights(NULL),
 fIntFlowCovariancesNUA(NULL),
 fIntFlowSumOfProductOfEventWeightsNUA(NULL),
 fIntFlowQcumulants(NULL),
 fIntFlowQcumulantsRebinnedInM(NULL), 
 fIntFlowQcumulantsErrorSquaredRatio(NULL), 
 fIntFlow(NULL),
 fIntFlowRebinnedInM(NULL),
 fIntFlowDetectorBias(NULL),
 // 4.) differential flow:
 fDiffFlowList(NULL),
 fDiffFlowProfiles(NULL),
 fDiffFlowResults(NULL),
 fDiffFlow2D(NULL),
 fDiffFlowFlags(NULL),
 fCalculateDiffFlow(kTRUE),
 fCalculate2DDiffFlow(kFALSE),
 fCalculateDiffFlowVsEta(kTRUE),
 // 5.) other differential correlators:
 fOtherDiffCorrelatorsList(NULL),
 // 6.) distributions:
 fDistributionsList(NULL),
 fDistributionsFlags(NULL),
 fStoreDistributions(kFALSE),
 // 7.) various:
 fVariousList(NULL),
 fPhiDistributionForOneEvent(NULL),
 // x.) debugging and cross-checking:
 fNestedLoopsList(NULL),
 fEvaluateIntFlowNestedLoops(kFALSE),
 fEvaluateDiffFlowNestedLoops(kFALSE),
 fMaxAllowedMultiplicity(10),
 fEvaluateNestedLoops(NULL),
 fIntFlowDirectCorrelations(NULL),
 fIntFlowExtraDirectCorrelations(NULL),
 fCrossCheckInPtBinNo(10),
 fCrossCheckInEtaBinNo(20),
 fNoOfParticlesInBin(NULL)
 {
  // constructor  
  
  // base list to hold all output objects:
  fHistList = new TList();
  fHistList->SetName("cobjQC");
  fHistList->SetOwner(kTRUE);
  
  // list to hold histograms with phi, pt and eta weights:      
  fWeightsList = new TList();
  
  // multiplicity weight:
  fMultiplicityWeight = new TString("combinations");
    
  // analysis label;
  fAnalysisLabel = new TString();
      
  // initialize all arrays:  
  this->InitializeArraysForIntFlow();
  this->InitializeArraysForDiffFlow();
  this->InitializeArraysForDistributions();
  this->InitializeArraysForVarious();
  this->InitializeArraysForNestedLoops();
  
 } // end of constructor
 
//================================================================================================================  

AliFlowAnalysisWithQCumulants::~AliFlowAnalysisWithQCumulants()
{
 // destructor
 
 delete fHistList;

} // end of AliFlowAnalysisWithQCumulants::~AliFlowAnalysisWithQCumulants()

//================================================================================================================

void AliFlowAnalysisWithQCumulants::Init()
{
 // a) Cross check if the settings make sense before starting the QC adventure;
 // b) Access all common constants;
 // c) Book all objects;
 // d) Store flags for integrated and differential flow;
 // e) Store flags for distributions of corelations;
 // f) Store harmonic which will be estimated.
  
 //save old value and prevent histograms from being added to directory
 //to avoid name clashes in case multiple analaysis objects are used
 //in an analysis
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
 TH1::AddDirectory(kFALSE);
 
 // a) Cross check if the settings make sense before starting the QC adventure; 
 this->CrossCheckSettings();
 // b) Access all common constants and book a profile to hold them:
 this->CommonConstants("Init");
 // c) Book all objects:
 this->BookAndFillWeightsHistograms(); 
 this->BookAndNestAllLists();
 this->BookCommonHistograms();
 this->BookEverythingForIntegratedFlow(); 
 this->BookEverythingForDifferentialFlow(); 
 this->BookEverythingFor2DDifferentialFlow(); 
 this->BookEverythingForDistributions();
 this->BookEverythingForVarious();
 this->BookEverythingForNestedLoops();
 // d) Store flags for integrated and differential flow:
 this->StoreIntFlowFlags();
 this->StoreDiffFlowFlags();
 // e) Store flags for distributions of corelations:
 this->StoreFlagsForDistributions();
 // f) Store harmonic which will be estimated:
 this->StoreHarmonic();
 
 TH1::AddDirectory(oldHistAddStatus);
} // end of void AliFlowAnalysisWithQCumulants::Init()

//================================================================================================================

void AliFlowAnalysisWithQCumulants::Make(AliFlowEventSimple* anEvent)
{
 // Running over data only in this method.
 
 // a) Check all pointers used in this method;
 // b) Define local variables;
 // c) Fill the common control histograms and call the method to fill fAvMultiplicity;
 // d) Loop over data and calculate e-b-e quantities Q_{n,k}, S_{p,k} and s_{p,k};
 // e) Calculate the final expressions for S_{p,k} and s_{p,k} (important !!!!); 
 // f) Call the methods which calculate correlations for reference flow;
 // g) Call the methods which calculate correlations for differential flow;
 // h) Call the methods which calculate correlations for 2D differential flow;
 // i) Call the methods which calculate other differential correlators;
 // j) Distributions of correlations;
 // k) Store phi distribution for one event to illustrate flow;
 // l) Cross-check with nested loops correlators for reference flow;
 // m) Cross-check with nested loops correlators for differential flow;
 // n) Reset all event-by-event quantities (very important !!!!). 
 
 // a) Check all pointers used in this method:
 this->CheckPointersUsedInMake();
 
 // b) Define local variables:
 Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
 Double_t dPt  = 0.; // transverse momentum
 Double_t dEta = 0.; // pseudorapidity
 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight
 Double_t wTrack = 1.; // track weight
 Int_t nRP = anEvent->GetEventNSelTracksRP(); // number of RPs (i.e. number of reference particles)
 fReferenceMultiplicityEBE = anEvent->GetReferenceMultiplicity(); // reference multiplicity for current event
 Double_t ptEta[2] = {0.,0.}; // 0 = dPt, 1 = dEta
  
 // c) Fill the common control histograms and call the method to fill fAvMultiplicity:
 this->FillCommonControlHistograms(anEvent);                                                               
 this->FillAverageMultiplicities(nRP);                                                                  
                                                                                                                                                                                                                                                                                        
 // d) Loop over data and calculate e-b-e quantities Q_{n,k}, S_{p,k} and s_{p,k}:
 Int_t nPrim = anEvent->NumberOfTracks();  // nPrim = total number of primary tracks, i.e. nPrim = nRP + nPOI where:
                                           //  nRP   = # of reference particles;
                                           //  nPOI  = # of particles of interest.
 AliFlowTrackSimple *aftsTrack = NULL;
 Int_t n = fHarmonic; // shortcut for the harmonic 
 for(Int_t i=0;i<nPrim;i++) 
 { 
  aftsTrack=anEvent->GetTrack(i);
  if(aftsTrack)
  {
   if(!(aftsTrack->InRPSelection() || aftsTrack->InPOISelection())){continue;} // safety measure: consider only tracks which are RPs or POIs
   if(aftsTrack->InRPSelection()) // RP condition:
   {    
    dPhi = aftsTrack->Phi();
    dPt  = aftsTrack->Pt();
    dEta = aftsTrack->Eta();
    if(fUsePhiWeights && fPhiWeights && fnBinsPhi) // determine phi weight for this particle:
    {
     wPhi = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*fnBinsPhi/TMath::TwoPi())));
    }
    if(fUsePtWeights && fPtWeights && fnBinsPt) // determine pt weight for this particle:
    {
     wPt = fPtWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/fPtBinWidth))); 
    }              
    if(fUseEtaWeights && fEtaWeights && fEtaBinWidth) // determine eta weight for this particle: 
    {
     wEta = fEtaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/fEtaBinWidth))); 
    }      
    // Access track weight:
    if(fUseTrackWeights)
    {
     wTrack = aftsTrack->Weight(); 
    }
    // Calculate Re[Q_{m*n,k}] and Im[Q_{m*n,k}] for this event (m = 1,2,...,6, k = 0,1,...,8):
    for(Int_t m=0;m<6;m++) // to be improved - hardwired 6 
    {
     for(Int_t k=0;k<9;k++) // to be improved - hardwired 9
     {
      (*fReQ)(m,k)+=pow(wPhi*wPt*wEta*wTrack,k)*TMath::Cos((m+1)*n*dPhi); 
      (*fImQ)(m,k)+=pow(wPhi*wPt*wEta*wTrack,k)*TMath::Sin((m+1)*n*dPhi); 
     } 
    }
    // Calculate S_{p,k} for this event (Remark: final calculation of S_{p,k} follows after the loop over data bellow):
    for(Int_t p=0;p<8;p++)
    {
     for(Int_t k=0;k<9;k++)
     {     
      (*fSpk)(p,k)+=pow(wPhi*wPt*wEta*wTrack,k);
     }
    } 
    // Differential flow:
    if(fCalculateDiffFlow || fCalculate2DDiffFlow)
    {
     ptEta[0] = dPt; 
     ptEta[1] = dEta; 
     // Calculate r_{m*n,k} and s_{p,k} (r_{m,k} is 'p-vector' for RPs): 
     for(Int_t k=0;k<9;k++) // to be improved - hardwired 9
     {
      for(Int_t m=0;m<4;m++) // to be improved - hardwired 4
      {
       if(fCalculateDiffFlow)
       {
        for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
        {
         fReRPQ1dEBE[0][pe][m][k]->Fill(ptEta[pe],pow(wPhi*wPt*wEta*wTrack,k)*TMath::Cos((m+1.)*n*dPhi),1.);
         fImRPQ1dEBE[0][pe][m][k]->Fill(ptEta[pe],pow(wPhi*wPt*wEta*wTrack,k)*TMath::Sin((m+1.)*n*dPhi),1.);          
         if(m==0) // s_{p,k} does not depend on index m
         {
          fs1dEBE[0][pe][k]->Fill(ptEta[pe],pow(wPhi*wPt*wEta*wTrack,k),1.);
         } // end of if(m==0) // s_{p,k} does not depend on index m
        } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
       } // end of if(fCalculateDiffFlow) 
       if(fCalculate2DDiffFlow)
       {
        fReRPQ2dEBE[0][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta*wTrack,k)*TMath::Cos((m+1.)*n*dPhi),1.);
        fImRPQ2dEBE[0][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta*wTrack,k)*TMath::Sin((m+1.)*n*dPhi),1.);      
        if(m==0) // s_{p,k} does not depend on index m
        {
         fs2dEBE[0][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta*wTrack,k),1.);
        } // end of if(m==0) // s_{p,k} does not depend on index m
       } // end of if(fCalculate2DDiffFlow)
      } // end of for(Int_t m=0;m<4;m++) // to be improved - hardwired 4
     } // end of for(Int_t k=0;k<9;k++) // to be improved - hardwired 9
     // Checking if RP particle is also POI particle:      
     if(aftsTrack->InPOISelection())
     {
      // Calculate q_{m*n,k} and s_{p,k} ('q-vector' and 's' for RPs && POIs): 
      for(Int_t k=0;k<9;k++) // to be improved - hardwired 9
      {
       for(Int_t m=0;m<4;m++) // to be improved - hardwired 4
       {
        if(fCalculateDiffFlow)
        {
         for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
         {
          fReRPQ1dEBE[2][pe][m][k]->Fill(ptEta[pe],pow(wPhi*wPt*wEta*wTrack,k)*TMath::Cos((m+1.)*n*dPhi),1.);
          fImRPQ1dEBE[2][pe][m][k]->Fill(ptEta[pe],pow(wPhi*wPt*wEta*wTrack,k)*TMath::Sin((m+1.)*n*dPhi),1.);          
          if(m==0) // s_{p,k} does not depend on index m
          {
           fs1dEBE[2][pe][k]->Fill(ptEta[pe],pow(wPhi*wPt*wEta*wTrack,k),1.);
          } // end of if(m==0) // s_{p,k} does not depend on index m
         } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
        } // end of if(fCalculateDiffFlow) 
        if(fCalculate2DDiffFlow)
        {
         fReRPQ2dEBE[2][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta*wTrack,k)*TMath::Cos((m+1.)*n*dPhi),1.);
         fImRPQ2dEBE[2][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta*wTrack,k)*TMath::Sin((m+1.)*n*dPhi),1.);      
         if(m==0) // s_{p,k} does not depend on index m
         {
          fs2dEBE[2][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta*wTrack,k),1.);
         } // end of if(m==0) // s_{p,k} does not depend on index m
        } // end of if(fCalculate2DDiffFlow)
       } // end of for(Int_t m=0;m<4;m++) // to be improved - hardwired 4
      } // end of for(Int_t k=0;k<9;k++) // to be improved - hardwired 9    
     } // end of if(aftsTrack->InPOISelection())  
    } // end of if(fCalculateDiffFlow || fCalculate2DDiffFlow)         
   } // end of if(pTrack->InRPSelection())
   if(aftsTrack->InPOISelection())
   {
    dPhi = aftsTrack->Phi();
    dPt  = aftsTrack->Pt();
    dEta = aftsTrack->Eta();
    wPhi = 1.;
    wPt  = 1.;
    wEta = 1.;
    wTrack = 1.;
    if(fUsePhiWeights && fPhiWeights && fnBinsPhi && aftsTrack->InRPSelection()) // determine phi weight for POI && RP particle:
    {
     wPhi = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*fnBinsPhi/TMath::TwoPi())));
    }
    if(fUsePtWeights && fPtWeights && fnBinsPt && aftsTrack->InRPSelection()) // determine pt weight for POI && RP particle:
    {
     wPt = fPtWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/fPtBinWidth))); 
    }              
    if(fUseEtaWeights && fEtaWeights && fEtaBinWidth && aftsTrack->InRPSelection()) // determine eta weight for POI && RP particle: 
    {
     wEta = fEtaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/fEtaBinWidth))); 
    }      
    // Access track weight for POI && RP particle:
    if(aftsTrack->InRPSelection() && fUseTrackWeights)
    {
     wTrack = aftsTrack->Weight(); 
    }
    ptEta[0] = dPt;
    ptEta[1] = dEta;
    // Calculate p_{m*n,k} ('p-vector' for POIs): 
    for(Int_t k=0;k<9;k++) // to be improved - hardwired 9
    {
     for(Int_t m=0;m<4;m++) // to be improved - hardwired 4
     {
      if(fCalculateDiffFlow)
      {
       for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
       {
        fReRPQ1dEBE[1][pe][m][k]->Fill(ptEta[pe],pow(wPhi*wPt*wEta*wTrack,k)*TMath::Cos((m+1.)*n*dPhi),1.);
        fImRPQ1dEBE[1][pe][m][k]->Fill(ptEta[pe],pow(wPhi*wPt*wEta*wTrack,k)*TMath::Sin((m+1.)*n*dPhi),1.);          
       } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
      } // end of if(fCalculateDiffFlow) 
      if(fCalculate2DDiffFlow)
      {
       fReRPQ2dEBE[1][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta*wTrack,k)*TMath::Cos((m+1.)*n*dPhi),1.);
       fImRPQ2dEBE[1][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta*wTrack,k)*TMath::Sin((m+1.)*n*dPhi),1.);      
      } // end of if(fCalculate2DDiffFlow)
     } // end of for(Int_t m=0;m<4;m++) // to be improved - hardwired 4
    } // end of for(Int_t k=0;k<9;k++) // to be improved - hardwired 9    
   } // end of if(pTrack->InPOISelection())    
  } else // to if(aftsTrack)
    {
     printf("\n WARNING (QC): No particle (i.e. aftsTrack is a NULL pointer in AFAWQC::Make())!!!!\n\n");
    }
 } // end of for(Int_t i=0;i<nPrim;i++) 

 // e) Calculate the final expressions for S_{p,k} and s_{p,k} (important !!!!):
 for(Int_t p=0;p<8;p++)
 {
  for(Int_t k=0;k<9;k++)
  {
   (*fSpk)(p,k)=pow((*fSpk)(p,k),p+1);
   // ... for the time being s_{p,k} dosn't need higher powers, so no need to finalize it here ...
  } // end of for(Int_t k=0;k<9;k++)  
 } // end of for(Int_t p=0;p<8;p++)
 
 // f) Call the methods which calculate correlations for reference flow:
 if(!fEvaluateIntFlowNestedLoops)
 {
  if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
  {
   if(nRP>1){this->CalculateIntFlowCorrelations();} // without using particle weights
  } else // to if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
    {
     if(nRP>1){this->CalculateIntFlowCorrelationsUsingParticleWeights();} // with using particle weights   
    }        
  // Whether or not using particle weights the following is calculated in the same way:  
  if(nRP>3){this->CalculateIntFlowProductOfCorrelations();}
  if(nRP>1){this->CalculateIntFlowSumOfEventWeights();}
  if(nRP>1){this->CalculateIntFlowSumOfProductOfEventWeights();}  
  // Non-isotropic terms:
  if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
  {
   if(nRP>0){this->CalculateIntFlowCorrectionsForNUASinTerms();}
   if(nRP>0){this->CalculateIntFlowCorrectionsForNUACosTerms();}
  } else // to if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
    {
     if(nRP>0){this->CalculateIntFlowCorrectionsForNUASinTermsUsingParticleWeights();}
     if(nRP>0){this->CalculateIntFlowCorrectionsForNUACosTermsUsingParticleWeights();}     
    }      
  // Whether or not using particle weights the following is calculated in the same way:  
  if(nRP>0){this->CalculateIntFlowProductOfCorrectionTermsForNUA();}     
  if(nRP>0){this->CalculateIntFlowSumOfEventWeightsNUA();}     
  if(nRP>0){this->CalculateIntFlowSumOfProductOfEventWeightsNUA();}      
 } // end of if(!fEvaluateIntFlowNestedLoops)

 // g) Call the methods which calculate correlations for differential flow:
 if(!fEvaluateDiffFlowNestedLoops && fCalculateDiffFlow)
 {
  if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
  {
   // Without using particle weights:
   this->CalculateDiffFlowCorrelations("RP","Pt"); 
   if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrelations("RP","Eta");}
   this->CalculateDiffFlowCorrelations("POI","Pt");
   if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrelations("POI","Eta");}
   // Non-isotropic terms:
   this->CalculateDiffFlowCorrectionsForNUASinTerms("RP","Pt");
   if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrectionsForNUASinTerms("RP","Eta");}
   this->CalculateDiffFlowCorrectionsForNUASinTerms("POI","Pt");
   if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrectionsForNUASinTerms("POI","Eta");}
   this->CalculateDiffFlowCorrectionsForNUACosTerms("RP","Pt");
   if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrectionsForNUACosTerms("RP","Eta");}
   this->CalculateDiffFlowCorrectionsForNUACosTerms("POI","Pt");
   if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrectionsForNUACosTerms("POI","Eta");}   
  } else // to if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
    {
     // With using particle weights:   
     this->CalculateDiffFlowCorrelationsUsingParticleWeights("RP","Pt"); 
     if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrelationsUsingParticleWeights("RP","Eta");} 
     this->CalculateDiffFlowCorrelationsUsingParticleWeights("POI","Pt"); 
     if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrelationsUsingParticleWeights("POI","Eta");} 
     // Non-isotropic terms:
     this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("RP","Pt");
     if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("RP","Eta");}
     this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("POI","Pt");
     if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("POI","Eta");}
     this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("RP","Pt");
     if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("RP","Eta");}
     this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("POI","Pt");
     if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("POI","Eta");}   
    }     
  // Whether or not using particle weights the following is calculated in the same way:  
  this->CalculateDiffFlowProductOfCorrelations("RP","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowProductOfCorrelations("RP","Eta");}
  this->CalculateDiffFlowProductOfCorrelations("POI","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowProductOfCorrelations("POI","Eta");}
  this->CalculateDiffFlowSumOfEventWeights("RP","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowSumOfEventWeights("RP","Eta");}
  this->CalculateDiffFlowSumOfEventWeights("POI","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowSumOfEventWeights("POI","Eta");}
  this->CalculateDiffFlowSumOfProductOfEventWeights("RP","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowSumOfProductOfEventWeights("RP","Eta");}
  this->CalculateDiffFlowSumOfProductOfEventWeights("POI","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowSumOfProductOfEventWeights("POI","Eta");}   
 } // end of if(!fEvaluateDiffFlowNestedLoops && fCalculateDiffFlow)

 // h) Call the methods which calculate correlations for 2D differential flow:
 if(!fEvaluateDiffFlowNestedLoops && fCalculate2DDiffFlow)
 {
  if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
  {
   // Without using particle weights:
   this->Calculate2DDiffFlowCorrelations("RP"); 
   this->Calculate2DDiffFlowCorrelations("POI");
   // Non-isotropic terms:
   // ... to be ctd ...
  } else // to if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
    {
     // With using particle weights:   
     // ... to be ctd ...  
     // Non-isotropic terms:
     // ... to be ctd ...
    }     
  // Whether or not using particle weights the following is calculated in the same way:  
  // ... to be ctd ...   
 } // end of if(!fEvaluateDiffFlowNestedLoops && fCalculate2DDiffFlow)
 
 // i) Call the methods which calculate other differential correlators:
 if(!fEvaluateDiffFlowNestedLoops && fCalculateDiffFlow)
 {
  if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
  {
   // Without using particle weights:
   this->CalculateOtherDiffCorrelators("RP","Pt"); 
   if(fCalculateDiffFlowVsEta){this->CalculateOtherDiffCorrelators("RP","Eta");}
   this->CalculateOtherDiffCorrelators("POI","Pt"); 
   if(fCalculateDiffFlowVsEta){this->CalculateOtherDiffCorrelators("POI","Eta");}     
  } else // to if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
    {
     // With using particle weights:   
     // ... to be ctd ...  
    }     
  // Whether or not using particle weights the following is calculated in the same way:  
  // ... to be ctd ...   
 } // end of if(!fEvaluateDiffFlowNestedLoops)
 
 // j) Distributions of correlations:
 if(fStoreDistributions){this->StoreDistributionsOfCorrelations();}
 
 // k) Store phi distribution for one event to illustrate flow: 
 if(fStorePhiDistributionForOneEvent){this->StorePhiDistributionForOneEvent(anEvent);}
   
 // l) Cross-check with nested loops correlators for reference flow:
 if(fEvaluateIntFlowNestedLoops){this->EvaluateIntFlowNestedLoops(anEvent);} 

 // m) Cross-check with nested loops correlators for differential flow:
 if(fEvaluateDiffFlowNestedLoops){this->EvaluateDiffFlowNestedLoops(anEvent);} 
 
 // n) Reset all event-by-event quantities (very important !!!!):
 this->ResetEventByEventQuantities();
 
} // end of AliFlowAnalysisWithQCumulants::Make(AliFlowEventSimple* anEvent)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::Finish()
{
 // Calculate the final results.
 
 // a) Check all pointers used in this method;
 // b) Acces the constants;
 // c) Access the flags;
 // d) Calculate reference cumulants (not corrected for detector effects);
 // e) Correct reference cumulants for detector effects;
 // f) Calculate reference flow;
 // g) Store results for reference flow in AliFlowCommonHistResults and print them on the screen;
 // h) Calculate the final results for differential flow (without/with weights);
 // i) Correct the results for differential flow (without/with weights) for effects of non-uniform acceptance (NUA);
 // j) Calculate the final results for integrated flow (RP/POI) and store in AliFlowCommonHistResults;
 // k) Store results for differential flow in AliFlowCommonHistResults;
 // l) Print the final results for integrated flow (RP/POI) on the screen; 
 // m) Cross-checking: Results from Q-vectors vs results from nested loops.
 
 // a) Check all pointers used in this method:
 this->CheckPointersUsedInFinish();
  
 // b) Acces the constants:
 this->CommonConstants("Finish");          
 
 if(fCommonHists && fCommonHists->GetHarmonic()) // to be improved (moved somewhere else)
 {
  fHarmonic = (Int_t)(fCommonHists->GetHarmonic())->GetBinContent(1);
 } 
 
 // c) Access the flags: // to be improved (implement a method for this? should I store again the flags becose they can get modified with redoFinish?)
 fUsePhiWeights = (Bool_t)fUseParticleWeights->GetBinContent(1); 
 fUsePtWeights = (Bool_t)fUseParticleWeights->GetBinContent(2); 
 fUseEtaWeights = (Bool_t)fUseParticleWeights->GetBinContent(3);  
 fUseTrackWeights = (Bool_t)fUseParticleWeights->GetBinContent(4);  
 fApplyCorrectionForNUA = (Bool_t)fIntFlowFlags->GetBinContent(3); 
 fPrintFinalResults[0] = (Bool_t)fIntFlowFlags->GetBinContent(4);
 fPrintFinalResults[1] = (Bool_t)fIntFlowFlags->GetBinContent(5);
 fPrintFinalResults[2] = (Bool_t)fIntFlowFlags->GetBinContent(6);
 fPrintFinalResults[3] = (Bool_t)fIntFlowFlags->GetBinContent(7);
 fApplyCorrectionForNUAVsM = (Bool_t)fIntFlowFlags->GetBinContent(8);  
 fPropagateErrorAlsoFromNIT = (Bool_t)fIntFlowFlags->GetBinContent(9);  
 fCalculateCumulantsVsM = (Bool_t)fIntFlowFlags->GetBinContent(10); 
 fMinimumBiasReferenceFlow = (Bool_t)fIntFlowFlags->GetBinContent(11); 
 fForgetAboutCovariances = (Bool_t)fIntFlowFlags->GetBinContent(12);
 fStorePhiDistributionForOneEvent = (Bool_t)fIntFlowFlags->GetBinContent(13);
 fFillMultipleControlHistograms = (Bool_t)fIntFlowFlags->GetBinContent(14); 
 fCalculateAllCorrelationsVsM = (Bool_t)fIntFlowFlags->GetBinContent(15);
 fEvaluateIntFlowNestedLoops = (Bool_t)fEvaluateNestedLoops->GetBinContent(1);
 fEvaluateDiffFlowNestedLoops = (Bool_t)fEvaluateNestedLoops->GetBinContent(2); 
 fCrossCheckInPtBinNo = (Int_t)fEvaluateNestedLoops->GetBinContent(3);
 fCrossCheckInEtaBinNo = (Int_t)fEvaluateNestedLoops->GetBinContent(4); 
          
 // d) Calculate reference cumulants (not corrected for detector effects):
 this->FinalizeCorrelationsIntFlow();
 this->CalculateCovariancesIntFlow();
 this->CalculateCumulantsIntFlow();

 // e) Correct reference cumulants for detector effects:
 this->FinalizeCorrectionTermsForNUAIntFlow();
 this->CalculateCovariancesNUAIntFlow(); 
 this->CalculateQcumulantsCorrectedForNUAIntFlow();  

 // f) Calculate reference flow:
 this->CalculateReferenceFlow(); 
  
 // g) Store results for reference flow in AliFlowCommonHistResults and print them on the screen:
 this->FillCommonHistResultsIntFlow();  
 if(fPrintFinalResults[0]){this->PrintFinalResultsForIntegratedFlow("RF");}
 if(fPrintFinalResults[3] && fCalculateCumulantsVsM){this->PrintFinalResultsForIntegratedFlow("RF, rebinned in M");}
 
 // h) Calculate the final results for differential flow (without/with weights):
 if(fCalculateDiffFlow)
 {
  this->FinalizeReducedCorrelations("RP","Pt"); 
  if(fCalculateDiffFlowVsEta){this->FinalizeReducedCorrelations("RP","Eta");} 
  this->FinalizeReducedCorrelations("POI","Pt"); 
  if(fCalculateDiffFlowVsEta){this->FinalizeReducedCorrelations("POI","Eta");}
  this->CalculateDiffFlowCovariances("RP","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCovariances("RP","Eta");}
  this->CalculateDiffFlowCovariances("POI","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCovariances("POI","Eta");}
  this->CalculateDiffFlowCumulants("RP","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCumulants("RP","Eta");}
  this->CalculateDiffFlowCumulants("POI","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCumulants("POI","Eta");}
  this->CalculateDiffFlow("RP","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlow("RP","Eta");}
  this->CalculateDiffFlow("POI","Pt");
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlow("POI","Eta");}
 } // if(fCalculateDiffFlow)
 
 // i) Correct the results for differential flow (without/with weights) for effects of non-uniform acceptance (NUA):
 if(fCalculateDiffFlow)
 {
  this->FinalizeCorrectionTermsForNUADiffFlow("RP","Pt");
  if(fCalculateDiffFlowVsEta){this->FinalizeCorrectionTermsForNUADiffFlow("RP","Eta");}
  this->FinalizeCorrectionTermsForNUADiffFlow("POI","Pt");
  if(fCalculateDiffFlowVsEta){this->FinalizeCorrectionTermsForNUADiffFlow("POI","Eta");}      
  this->CalculateDiffFlowCumulantsCorrectedForNUA("RP","Pt");   
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCumulantsCorrectedForNUA("RP","Eta");}   
  this->CalculateDiffFlowCumulantsCorrectedForNUA("POI","Pt");   
  if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCumulantsCorrectedForNUA("POI","Eta");}  
  if(fApplyCorrectionForNUA)
  {
   this->CalculateDiffFlowCorrectedForNUA("RP","Pt"); 
   if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrectedForNUA("RP","Eta");} 
   this->CalculateDiffFlowCorrectedForNUA("POI","Pt"); 
   if(fCalculateDiffFlowVsEta){this->CalculateDiffFlowCorrectedForNUA("POI","Eta");} 
  }
 } // end of if(fCalculateDiffFlow && fApplyCorrectionForNUA)
 
 // i) Calcualate final results for 2D differential flow: 
 if(fCalculate2DDiffFlow)
 {
  this->Calculate2DDiffFlowCumulants("RP");
  this->Calculate2DDiffFlowCumulants("POI");
  this->Calculate2DDiffFlow("RP");  
  this->Calculate2DDiffFlow("POI");  
 } // end of if(fCalculate2DDiffFlow)
    
 // j) Calculate the final results for integrated flow (RP/POI) and store in AliFlowCommonHistResults:
 if(fCalculateDiffFlow)
 {
  this->CalculateFinalResultsForRPandPOIIntegratedFlow("RP");
  this->CalculateFinalResultsForRPandPOIIntegratedFlow("POI");
 }
 
 // k) Store results for differential flow in AliFlowCommonHistResults:
 if(fCalculateDiffFlow)
 {
  this->FillCommonHistResultsDiffFlow("RP");
  this->FillCommonHistResultsDiffFlow("POI");
 }
 
 // l) Print the final results for integrated flow (RP/POI) on the screen:
 if(fPrintFinalResults[1] && fCalculateDiffFlow){this->PrintFinalResultsForIntegratedFlow("RP");} 
 if(fPrintFinalResults[2] && fCalculateDiffFlow){this->PrintFinalResultsForIntegratedFlow("POI");}
    
 // m) Cross-checking: Results from Q-vectors vs results from nested loops:
 //  m1) Reference flow:
 if(fEvaluateIntFlowNestedLoops)
 {
  this->CrossCheckIntFlowCorrelations();
  this->CrossCheckIntFlowCorrectionTermsForNUA(); 
  if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights){this->CrossCheckIntFlowExtraCorrelations();}     
 } // end of if(fEvaluateIntFlowNestedLoops)  
 //  m2) Differential flow: 
 if(fEvaluateDiffFlowNestedLoops && fCalculateDiffFlow) 
 {
  // Correlations:
  this->PrintNumberOfParticlesInSelectedBin();
  this->CrossCheckDiffFlowCorrelations("RP","Pt");  
  if(fCalculateDiffFlowVsEta){this->CrossCheckDiffFlowCorrelations("RP","Eta");} 
  this->CrossCheckDiffFlowCorrelations("POI","Pt");  
  if(fCalculateDiffFlowVsEta){this->CrossCheckDiffFlowCorrelations("POI","Eta");}
  // Correction terms for non-uniform acceptance:
  this->CrossCheckDiffFlowCorrectionTermsForNUA("RP","Pt");      
  if(fCalculateDiffFlowVsEta){this->CrossCheckDiffFlowCorrectionTermsForNUA("RP","Eta");}       
  this->CrossCheckDiffFlowCorrectionTermsForNUA("POI","Pt");      
  if(fCalculateDiffFlowVsEta){this->CrossCheckDiffFlowCorrectionTermsForNUA("POI","Eta");}
  // Other differential correlators:       
  this->CrossCheckOtherDiffCorrelators("RP","Pt");  
  if(fCalculateDiffFlowVsEta){this->CrossCheckOtherDiffCorrelators("RP","Eta");} 
  this->CrossCheckOtherDiffCorrelators("POI","Pt");  
  if(fCalculateDiffFlowVsEta){this->CrossCheckOtherDiffCorrelators("POI","Eta");}
 } // end of if(fEvaluateDiffFlowNestedLoops)
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithQCumulants::Finish()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::EvaluateIntFlowNestedLoops(AliFlowEventSimple* anEvent)
{
 // Evalauted all correlators for reference flow with nested loops.
 
 Int_t nPrim = anEvent->NumberOfTracks(); // nPrim = nRP + nPOI 
 if(nPrim>0 && nPrim<=fMaxAllowedMultiplicity) // by default fMaxAllowedMultiplicity = 10 
 {
  // Without using particle weights:
  if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
  {
   // Correlations:
   this->CalculateIntFlowCorrelations(); // from Q-vectors
   this->EvaluateIntFlowCorrelationsWithNestedLoops(anEvent); // from nested loops (to be improved: do I have to pass here anEvent or not?)
   // Correction for non-uniform acceptance:
   this->CalculateIntFlowCorrectionsForNUASinTerms(); // from Q-vectors (sin terms)
   this->CalculateIntFlowCorrectionsForNUACosTerms(); // from Q-vectors (cos terms)
   this->EvaluateIntFlowCorrectionsForNUAWithNestedLoops(anEvent); // from nested loops (both sin and cos terms)
  }
  // Using particle weights:
  if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)
  {
   // Correlations
   this->CalculateIntFlowCorrelationsUsingParticleWeights(); // from Q-vectors
   this->EvaluateIntFlowCorrelationsWithNestedLoopsUsingParticleWeights(anEvent); // from nested loops (to be improved: do I have to pass here anEvent or not?)
   // Correction for non-uniform acceptance:
   this->CalculateIntFlowCorrectionsForNUASinTermsUsingParticleWeights(); // from Q-vectors (sin terms)
   this->CalculateIntFlowCorrectionsForNUACosTermsUsingParticleWeights(); // from Q-vectors (cos terms)
   this->EvaluateIntFlowCorrectionsForNUAWithNestedLoopsUsingParticleWeights(anEvent); // from nested loops (both sin and cos terms)   
  }
 } else if(nPrim>fMaxAllowedMultiplicity) // to if(nPrim>0 && nPrim<=fMaxAllowedMultiplicity)
   {
    cout<<endl;
    cout<<"Skipping the event because multiplicity is "<<nPrim<<". Too high to evaluate nested loops!"<<endl;
   } else
     {
      cout<<endl;
      cout<<"Skipping the event because multiplicity is "<<nPrim<<"."<<endl;      
     } 

} // end of void AliFlowAnalysisWithQCumulants::EvaluateIntFlowNestedLoops(AliFlowEventSimple* anEvent)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::EvaluateDiffFlowNestedLoops(AliFlowEventSimple* anEvent)
{
 // Evalauted all correlators for differential flow with nested loops.

 if(!fCalculateDiffFlow){return;}

 Int_t nPrim = anEvent->NumberOfTracks(); // nPrim = nRP + nPOI 
 if(nPrim>0 && nPrim<=fMaxAllowedMultiplicity) // by default fMaxAllowedMultiplicity = 10
 {
  // Without using particle weights:
  if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
  {
   // 1.) Reduced correlations:
   //  Q-vectors:
   this->CalculateDiffFlowCorrelations("RP","Pt");
   this->CalculateDiffFlowCorrelations("RP","Eta");
   this->CalculateDiffFlowCorrelations("POI","Pt");
   this->CalculateDiffFlowCorrelations("POI","Eta");
   //  Nested loops:
   this->EvaluateDiffFlowCorrelationsWithNestedLoops(anEvent,"RP","Pt"); 
   this->EvaluateDiffFlowCorrelationsWithNestedLoops(anEvent,"RP","Eta"); 
   this->EvaluateDiffFlowCorrelationsWithNestedLoops(anEvent,"POI","Pt"); 
   this->EvaluateDiffFlowCorrelationsWithNestedLoops(anEvent,"POI","Eta"); 
   // 2.) Reduced corrections for non-uniform acceptance:
   //  Q-vectors:
   this->CalculateDiffFlowCorrectionsForNUASinTerms("RP","Pt");
   this->CalculateDiffFlowCorrectionsForNUASinTerms("RP","Eta");
   this->CalculateDiffFlowCorrectionsForNUASinTerms("POI","Pt");
   this->CalculateDiffFlowCorrectionsForNUASinTerms("POI","Eta");
   this->CalculateDiffFlowCorrectionsForNUACosTerms("RP","Pt");
   this->CalculateDiffFlowCorrectionsForNUACosTerms("RP","Eta");
   this->CalculateDiffFlowCorrectionsForNUACosTerms("POI","Pt");
   this->CalculateDiffFlowCorrectionsForNUACosTerms("POI","Eta");
   //  Nested loops:
   this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(anEvent,"RP","Pt");
   this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(anEvent,"RP","Eta");
   this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(anEvent,"POI","Pt"); 
   this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(anEvent,"POI","Eta"); 
   // 3.) Other differential correlators:
   //  Q-vectors:
   this->CalculateOtherDiffCorrelators("RP","Pt");
   this->CalculateOtherDiffCorrelators("RP","Eta");
   this->CalculateOtherDiffCorrelators("POI","Pt");
   this->CalculateOtherDiffCorrelators("POI","Eta");   
   //  Nested loops:
   this->EvaluateOtherDiffCorrelatorsWithNestedLoops(anEvent,"RP","Pt");
   this->EvaluateOtherDiffCorrelatorsWithNestedLoops(anEvent,"RP","Eta");
   this->EvaluateOtherDiffCorrelatorsWithNestedLoops(anEvent,"POI","Pt");
   this->EvaluateOtherDiffCorrelatorsWithNestedLoops(anEvent,"POI","Eta");   
  } // end of if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
  // Using particle weights:
  if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)
  {
   this->CalculateDiffFlowCorrelationsUsingParticleWeights("RP","Pt"); 
   this->CalculateDiffFlowCorrelationsUsingParticleWeights("RP","Eta"); 
   this->CalculateDiffFlowCorrelationsUsingParticleWeights("POI","Pt"); 
   this->CalculateDiffFlowCorrelationsUsingParticleWeights("POI","Eta"); 
   this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("RP","Pt");
   this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("RP","Eta");
   this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("POI","Pt");
   this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("POI","Eta");
   this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("RP","Pt");
   this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("RP","Eta");
   this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("POI","Pt");
   this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("POI","Eta");
   this->EvaluateDiffFlowCorrelationsWithNestedLoopsUsingParticleWeights(anEvent,"RP","Pt"); 
   this->EvaluateDiffFlowCorrelationsWithNestedLoopsUsingParticleWeights(anEvent,"RP","Eta");
   this->EvaluateDiffFlowCorrelationsWithNestedLoopsUsingParticleWeights(anEvent,"POI","Pt"); 
   this->EvaluateDiffFlowCorrelationsWithNestedLoopsUsingParticleWeights(anEvent,"POI","Eta");   
   this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoopsUsingParticleWeights(anEvent,"RP","Pt"); 
   this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoopsUsingParticleWeights(anEvent,"RP","Eta"); 
   this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoopsUsingParticleWeights(anEvent,"POI","Pt"); 
   this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoopsUsingParticleWeights(anEvent,"POI","Eta"); 
  } // end of if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)
 } // end of if(nPrim>0 && nPrim<=fMaxAllowedMultiplicity) // by default fMaxAllowedMultiplicity = 10

} // end of void AliFlowAnalysisWithQCumulants::EvaluateDiffFlowNestedLoops(AliFlowEventSimple* anEvent)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUACosTerms()
{
 // Calculate correction terms for non-uniform acceptance of the detector for reference flow (cos terms).
 
 // multiplicity:
 Double_t dMult = (*fSpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 //Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 //Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);
        
 //                                  *************************************************************
 //                                  **** corrections for non-uniform acceptance (cos terms): ****
 //                                  *************************************************************
 //
 // Remark 1: corrections for non-uniform acceptance (cos terms) calculated with non-weighted Q-vectors 
 //           are stored in 1D profile fQCorrectionsCos.
 // Remark 2: binning of fIntFlowCorrectionTermsForNUAPro[1] is organized as follows:
 // --------------------------------------------------------------------------------------------------------------------
 // 1st bin: <<cos(n*(phi1))>> = cosP1n
 // 2nd bin: <<cos(n*(phi1+phi2))>> = cosP1nP1n
 // 3rd bin: <<cos(n*(phi1-phi2-phi3))>> = cosP1nM1nM1n
 // 4th bin: <<cos(n*(2phi1-phi2))>> = cosP2nM1n
 // --------------------------------------------------------------------------------------------------------------------
  
 // 1-particle:
 Double_t cosP1n = 0.; // <<cos(n*(phi1))>>
   
 if(dMult>0)
 {
  cosP1n = dReQ1n/dMult; 
  
  // average non-weighted 1-particle correction (cos terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[1]->SetBinContent(1,cosP1n);
  // event weights for NUA terms:
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->SetBinContent(1,dMult);
  
  // final average non-weighted 1-particle correction (cos terms) for non-uniform acceptance for all events:
  fIntFlowCorrectionTermsForNUAPro[1]->Fill(0.5,cosP1n,dMult);  
  if(fCalculateCumulantsVsM){fIntFlowCorrectionTermsForNUAVsMPro[1][0]->Fill(dMult+0.5,cosP1n,dMult);}    
 } 
 
 // 2-particle:
 Double_t cosP1nP1n = 0.; // <<cos(n*(phi1+phi2))>>
 Double_t cosP2nM1n = 0.; // <<cos(n*(2phi1-phi2))>>
 
 if(dMult>1)
 {
  cosP1nP1n = (pow(dReQ1n,2)-pow(dImQ1n,2)-dReQ2n)/(dMult*(dMult-1)); 
  cosP2nM1n = (dReQ2n*dReQ1n+dImQ2n*dImQ1n-dReQ1n)/(dMult*(dMult-1)); 
  
  // average non-weighted 2-particle correction (cos terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[1]->SetBinContent(2,cosP1nP1n);
  fIntFlowCorrectionTermsForNUAEBE[1]->SetBinContent(4,cosP2nM1n);
  // event weights for NUA terms:
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->SetBinContent(2,dMult*(dMult-1));
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->SetBinContent(4,dMult*(dMult-1));
      
  // final average non-weighted 2-particle correction (cos terms) for non-uniform acceptance for all events:
  fIntFlowCorrectionTermsForNUAPro[1]->Fill(1.5,cosP1nP1n,dMult*(dMult-1));  
  fIntFlowCorrectionTermsForNUAPro[1]->Fill(3.5,cosP2nM1n,dMult*(dMult-1));
  if(fCalculateCumulantsVsM)
  {
   fIntFlowCorrectionTermsForNUAVsMPro[1][1]->Fill(dMult+0.5,cosP1nP1n,dMult*(dMult-1));  
   fIntFlowCorrectionTermsForNUAVsMPro[1][3]->Fill(dMult+0.5,cosP2nM1n,dMult*(dMult-1));
  }
 } 
 
 // 3-particle:
 Double_t cosP1nM1nM1n = 0.; // <<cos(n*(phi1-phi2-phi3))>>
 
 if(dMult>2)
 {
  cosP1nM1nM1n = (dReQ1n*(pow(dReQ1n,2)+pow(dImQ1n,2))-dReQ1n*dReQ2n-dImQ1n*dImQ2n-2.*(dMult-1)*dReQ1n)
               / (dMult*(dMult-1)*(dMult-2)); 
  
  // average non-weighted 3-particle correction (cos terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[1]->SetBinContent(3,cosP1nM1nM1n);
  // event weights for NUA terms:
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->SetBinContent(3,dMult*(dMult-1)*(dMult-2));
  
  // final average non-weighted 3-particle correction (cos terms) for non-uniform acceptance for all events:
  fIntFlowCorrectionTermsForNUAPro[1]->Fill(2.5,cosP1nM1nM1n,dMult*(dMult-1)*(dMult-2));
  if(fCalculateCumulantsVsM){fIntFlowCorrectionTermsForNUAVsMPro[1][2]->Fill(dMult+0.5,cosP1nM1nM1n,dMult*(dMult-1)*(dMult-2));}  
 } 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUACosTerms()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUASinTerms()
{
 // calculate corrections for non-uniform acceptance of the detector for no-name integrated flow (sin terms)
 
 // multiplicity:
 Double_t dMult = (*fSpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 //Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 //Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);
        
 //                                  *************************************************************
 //                                  **** corrections for non-uniform acceptance (sin terms): ****
 //                                  *************************************************************
 //
 // Remark 1: corrections for non-uniform acceptance (sin terms) calculated with non-weighted Q-vectors 
 //           are stored in 1D profile fQCorrectionsSin.
 // Remark 2: binning of fIntFlowCorrectionTermsForNUAPro[0] is organized as follows:
 // --------------------------------------------------------------------------------------------------------------------
 // 1st bin: <<sin(n*(phi1))>> = sinP1n
 // 2nd bin: <<sin(n*(phi1+phi2))>> = sinP1nP1n
 // 3rd bin: <<sin(n*(phi1-phi2-phi3))>> = sinP1nM1nM1n
 // 4th bin: <<sin(n*(2phi1-phi2))>> = sinP2nM1n
 // --------------------------------------------------------------------------------------------------------------------
 
 // 1-particle:
 Double_t sinP1n = 0.; // <sin(n*(phi1))>
 
 if(dMult>0)
 {
  sinP1n = dImQ1n/dMult; 
     
  // average non-weighted 1-particle correction (sin terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[0]->SetBinContent(1,sinP1n);  
  // event weights for NUA terms:
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->SetBinContent(1,dMult);
  
  // final average non-weighted 1-particle correction (sin terms) for non-uniform acceptance for all events:   
  fIntFlowCorrectionTermsForNUAPro[0]->Fill(0.5,sinP1n,dMult);  
  if(fCalculateCumulantsVsM){fIntFlowCorrectionTermsForNUAVsMPro[0][0]->Fill(dMult+0.5,sinP1n,dMult);} 
 } 
 
 // 2-particle:
 Double_t sinP1nP1n = 0.; // <<sin(n*(phi1+phi2))>>
 Double_t sinP2nM1n = 0.; // <<sin(n*(2phi1-phi2))>>
 if(dMult>1)
 {
  sinP1nP1n = (2.*dReQ1n*dImQ1n-dImQ2n)/(dMult*(dMult-1)); 
  sinP2nM1n = (dImQ2n*dReQ1n-dReQ2n*dImQ1n-dImQ1n)/(dMult*(dMult-1)); 
     
  // average non-weighted 2-particle correction (sin terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[0]->SetBinContent(2,sinP1nP1n);
  fIntFlowCorrectionTermsForNUAEBE[0]->SetBinContent(4,sinP2nM1n);
  // event weights for NUA terms:
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->SetBinContent(2,dMult*(dMult-1));
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->SetBinContent(4,dMult*(dMult-1));
  
  // final average non-weighted 1-particle correction (sin terms) for non-uniform acceptance for all events:      
  fIntFlowCorrectionTermsForNUAPro[0]->Fill(1.5,sinP1nP1n,dMult*(dMult-1));  
  fIntFlowCorrectionTermsForNUAPro[0]->Fill(3.5,sinP2nM1n,dMult*(dMult-1));  
  if(fCalculateCumulantsVsM)
  {
   fIntFlowCorrectionTermsForNUAVsMPro[0][1]->Fill(dMult+0.5,sinP1nP1n,dMult*(dMult-1));  
   fIntFlowCorrectionTermsForNUAVsMPro[0][3]->Fill(dMult+0.5,sinP2nM1n,dMult*(dMult-1));    
  }
 } 
 
 // 3-particle:
 Double_t sinP1nM1nM1n = 0.; // <<sin(n*(phi1-phi2-phi3))>>
 
 if(dMult>2)
 {
  sinP1nM1nM1n = (-dImQ1n*(pow(dReQ1n,2)+pow(dImQ1n,2))+dReQ1n*dImQ2n-dImQ1n*dReQ2n+2.*(dMult-1)*dImQ1n)
               / (dMult*(dMult-1)*(dMult-2)); 
  
  // average non-weighted 3-particle correction (sin terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[0]->SetBinContent(3,sinP1nM1nM1n);
  // event weights for NUA terms:
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->SetBinContent(3,dMult*(dMult-1)*(dMult-2));
  
  // final average non-weighted 3-particle correction (sin terms) for non-uniform acceptance for all events:  
  fIntFlowCorrectionTermsForNUAPro[0]->Fill(2.5,sinP1nM1nM1n,dMult*(dMult-1)*(dMult-2));
  if(fCalculateCumulantsVsM){fIntFlowCorrectionTermsForNUAVsMPro[0][2]->Fill(dMult+0.5,sinP1nM1nM1n,dMult*(dMult-1)*(dMult-2));}  
 } 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUASinTerms()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::GetOutputHistograms(TList *outputListHistos)
{
 // a) Get pointers for common control and common result histograms;
 // b) Get pointers for histograms holding particle weights;
 // c) Get pointers for reference flow histograms;
 // d) Get pointers for differential flow histograms;
 // e) Get pointers for 2D differential flow histograms;
 // f) Get pointers for other differential correlators;
 // g) Get pointers for nested loops' histograms.
 
 if(outputListHistos)
 {	
  this->SetHistList(outputListHistos);
  if(!fHistList)
  {
   printf("\n WARNING (QC): fHistList is NULL in AFAWQC::GOH() !!!!\n\n");
   exit(0);
  }
  this->GetPointersForCommonHistograms(); 
  this->GetPointersForParticleWeightsHistograms(); 
  this->GetPointersForIntFlowHistograms();
  this->GetPointersForDiffFlowHistograms(); 
  this->GetPointersFor2DDiffFlowHistograms(); 
  this->GetPointersForOtherDiffCorrelators();  
  this->GetPointersForNestedLoopsHistograms(); 
 } else 
   {
    printf("\n WARNING (QC): outputListHistos is NULL in AFAWQC::GOH() !!!!\n\n");
    exit(0);
   }
   
} // end of void AliFlowAnalysisWithQCumulants::GetOutputHistograms(TList *outputListHistos)

//================================================================================================================================

TProfile* AliFlowAnalysisWithQCumulants::MakePtProjection(TProfile2D *profilePtEta) const
{
 // project 2D profile onto pt axis to get 1D profile
 
 Int_t nBinsPt   = profilePtEta->GetNbinsX();
 Double_t dPtMin = (profilePtEta->GetXaxis())->GetXmin();
 Double_t dPtMax = (profilePtEta->GetXaxis())->GetXmax();
 
 Int_t nBinsEta   = profilePtEta->GetNbinsY();
 
 TProfile *profilePt = new TProfile("","",nBinsPt,dPtMin,dPtMax); 
 
 for(Int_t p=1;p<=nBinsPt;p++)
 {
  Double_t contentPt = 0.;
  Double_t entryPt = 0.;
  Double_t spreadPt = 0.;
  Double_t sum1 = 0.;
  Double_t sum2 = 0.;
  Double_t sum3 = 0.;
  for(Int_t e=1;e<=nBinsEta;e++)
  {
   contentPt += (profilePtEta->GetBinContent(profilePtEta->GetBin(p,e)))
              * (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)));
   entryPt   += (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)));
   
   sum1 += (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)))
         * (pow(profilePtEta->GetBinError(profilePtEta->GetBin(p,e)),2.)
            + pow(profilePtEta->GetBinContent(profilePtEta->GetBin(p,e)),2.)); 
   sum2 += (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)));
   sum3 += (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)))
         * (profilePtEta->GetBinContent(profilePtEta->GetBin(p,e)));            
  }
  if(sum2>0. && sum1/sum2-pow(sum3/sum2,2.) > 0.)
  {
   spreadPt = pow(sum1/sum2-pow(sum3/sum2,2.),0.5);
  }
  profilePt->SetBinContent(p,contentPt);
  profilePt->SetBinEntries(p,entryPt);
  {
   profilePt->SetBinError(p,spreadPt);
  }
  
 }
 
 return profilePt;
 
} // end of TProfile* AliFlowAnalysisWithQCumulants::MakePtProjection(TProfile2D *profilePtEta)


//================================================================================================================================


TProfile* AliFlowAnalysisWithQCumulants::MakeEtaProjection(TProfile2D *profilePtEta) const
{
 // project 2D profile onto eta axis to get 1D profile
 
 Int_t nBinsEta   = profilePtEta->GetNbinsY();
 Double_t dEtaMin = (profilePtEta->GetYaxis())->GetXmin();
 Double_t dEtaMax = (profilePtEta->GetYaxis())->GetXmax();
 
 Int_t nBinsPt = profilePtEta->GetNbinsX();
 
 TProfile *profileEta = new TProfile("","",nBinsEta,dEtaMin,dEtaMax); 
 
 for(Int_t e=1;e<=nBinsEta;e++)
 {
  Double_t contentEta = 0.;
  Double_t entryEta = 0.;
  for(Int_t p=1;p<=nBinsPt;p++)
  {
   contentEta += (profilePtEta->GetBinContent(profilePtEta->GetBin(p,e)))
              * (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)));
   entryEta   += (profilePtEta->GetBinEntries(profilePtEta->GetBin(p,e)));
  }
  profileEta->SetBinContent(e,contentEta);
  profileEta->SetBinEntries(e,entryEta);
 }
 
 return profileEta;
 
} // end of TProfile* AliFlowAnalysisWithQCumulants::MakeEtaProjection(TProfile2D *profilePtEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::PrintFinalResultsForIntegratedFlow(TString type)
{
 // Printing on the screen the final results for integrated flow (RF, POI and RP). 
 
 Int_t n = fHarmonic; 
 
 Double_t dVn[4] = {0.}; // array to hold Vn{2}, Vn{4}, Vn{6} and Vn{8}   
 Double_t dVnErr[4] = {0.}; // array to hold errors of Vn{2}, Vn{4}, Vn{6} and Vn{8}   
 
 if(type == "RF")
 {
  for(Int_t b=0;b<4;b++)
  {
   dVn[0] = (fCommonHistsResults2nd->GetHistIntFlow())->GetBinContent(1); 
   dVnErr[0] = (fCommonHistsResults2nd->GetHistIntFlow())->GetBinError(1); 
   dVn[1] = (fCommonHistsResults4th->GetHistIntFlow())->GetBinContent(1); 
   dVnErr[1] = (fCommonHistsResults4th->GetHistIntFlow())->GetBinError(1); 
   dVn[2] = (fCommonHistsResults6th->GetHistIntFlow())->GetBinContent(1); 
   dVnErr[2] = (fCommonHistsResults6th->GetHistIntFlow())->GetBinError(1); 
   dVn[3] = (fCommonHistsResults8th->GetHistIntFlow())->GetBinContent(1); 
   dVnErr[3] = (fCommonHistsResults8th->GetHistIntFlow())->GetBinError(1);    
  }  
 } else if(type == "RP")
   {
    dVn[0] = (fCommonHistsResults2nd->GetHistIntFlowRP())->GetBinContent(1); 
    dVnErr[0] = (fCommonHistsResults2nd->GetHistIntFlowRP())->GetBinError(1); 
    dVn[1] = (fCommonHistsResults4th->GetHistIntFlowRP())->GetBinContent(1); 
    dVnErr[1] = (fCommonHistsResults4th->GetHistIntFlowRP())->GetBinError(1); 
    dVn[2] = (fCommonHistsResults6th->GetHistIntFlowRP())->GetBinContent(1); 
    dVnErr[2] = (fCommonHistsResults6th->GetHistIntFlowRP())->GetBinError(1); 
    dVn[3] = (fCommonHistsResults8th->GetHistIntFlowRP())->GetBinContent(1); 
    dVnErr[3] = (fCommonHistsResults8th->GetHistIntFlowRP())->GetBinError(1); 
   } else if(type == "POI")
     {
      dVn[0] = (fCommonHistsResults2nd->GetHistIntFlowPOI())->GetBinContent(1); 
      dVnErr[0] = (fCommonHistsResults2nd->GetHistIntFlowPOI())->GetBinError(1); 
      dVn[1] = (fCommonHistsResults4th->GetHistIntFlowPOI())->GetBinContent(1); 
      dVnErr[1] = (fCommonHistsResults4th->GetHistIntFlowPOI())->GetBinError(1); 
      dVn[2] = (fCommonHistsResults6th->GetHistIntFlowPOI())->GetBinContent(1); 
      dVnErr[2] = (fCommonHistsResults6th->GetHistIntFlowPOI())->GetBinError(1); 
      dVn[3] = (fCommonHistsResults8th->GetHistIntFlowPOI())->GetBinContent(1); 
      dVnErr[3] = (fCommonHistsResults8th->GetHistIntFlowPOI())->GetBinError(1); 
     } else if(type == "RF, rebinned in M" && fCalculateCumulantsVsM)
       {
        for(Int_t b=0;b<4;b++)
        {
         dVn[b] = fIntFlowRebinnedInM->GetBinContent(b+1); 
         dVnErr[b] = fIntFlowRebinnedInM->GetBinError(b+1);
        }  
       }
 
 TString title = " flow estimates from Q-cumulants"; 
 TString subtitle = "    ("; 
 TString subtitle2 = "       (rebinned in M)"; 
 
 if(type != "RF, rebinned in M")
 {
  if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
  {
   subtitle.Append(type);
   subtitle.Append(", without weights)");
  } else  
    {
     subtitle.Append(type);
     subtitle.Append(", with weights)");
    }
 } else
   {
    if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
    {
     subtitle.Append("RF");
     subtitle.Append(", without weights)");
    } else  
      {
       subtitle.Append("RF");
       subtitle.Append(", with weights)");      
      }
   } 
   
 cout<<endl;
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<title.Data()<<endl; 
 cout<<subtitle.Data()<<endl; 
 if(type == "RF, rebinned in M"){cout<<subtitle2.Data()<<endl;}
 cout<<endl;
  
 for(Int_t i=0;i<4;i++)
 {
  cout<<"  v_"<<n<<"{"<<2*(i+1)<<"} = "<<dVn[i]<<" +/- "<<dVnErr[i]<<endl;
 }
 
 cout<<endl;
 if(type == "RF")
 {
  if(fApplyCorrectionForNUA)
  {
   cout<<" detector bias (corrected for): "<<endl;
  } else
    {
     cout<<" detector bias (not corrected for):"<<endl;  
    }
  cout<<"  to QC{2}: "<<fIntFlowDetectorBias->GetBinContent(1)<<" +/- "<<fIntFlowDetectorBias->GetBinError(1)<<endl;
  cout<<"  to QC{4}: "<<fIntFlowDetectorBias->GetBinContent(2)<<" +/- "<<fIntFlowDetectorBias->GetBinError(2)<<endl;
  cout<<endl;
 }
 if(type == "RF" || type == "RF, rebinned in M")
 {
  cout<<"     nEvts = "<<(Int_t)fCommonHists->GetHistMultRP()->GetEntries()<<", <M> = "<<(Double_t)fCommonHists->GetHistMultRP()->GetMean()<<endl; 
 }
 else if (type == "RP")
 {
  cout<<"     nEvts = "<<(Int_t)fCommonHists->GetHistMultRP()->GetEntries()<<", <M> = "<<(Double_t)fCommonHists->GetHistMultRP()->GetMean()<<endl;  
 } 
 else if (type == "POI")
 {
  cout<<"     nEvts = "<<(Int_t)fCommonHists->GetHistMultPOI()->GetEntries()<<", <M> = "<<(Double_t)fCommonHists->GetHistMultPOI()->GetMean()<<endl;
 }  
 
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<endl; 
  
}// end of AliFlowAnalysisWithQCumulants::PrintFinalResultsForIntegratedFlow(TString type="RF");

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 //output->WriteObject(fHistList, "cobjQC","SingleKey");
 fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
 delete output;
}


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::WriteHistograms(TDirectoryFile *outputFileName)
{
 //store the final results in output .root file
 fHistList->SetName("cobjQC");
 fHistList->SetOwner(kTRUE);
 outputFileName->Add(fHistList);
 outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);
}

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::BookCommonHistograms()
{
 // Book common control histograms and common histograms for final results.
 //  a) Book common control histograms;
 //  b) Book common result histograms.
 
 // a) Book common control histograms: 
 //  Common control histograms (all events):
 TString commonHistsName = "AliFlowCommonHistQC";
 commonHistsName += fAnalysisLabel->Data();
 fCommonHists = new AliFlowCommonHist(commonHistsName.Data(),commonHistsName.Data(),fBookOnlyBasicCCH);
 fHistList->Add(fCommonHists);  
 //  Common control histograms (selected events):
 if(fFillMultipleControlHistograms)
 {
  // Common control histogram filled for events with 2 and more reference particles:
  TString commonHists2ndOrderName = "AliFlowCommonHist2ndOrderQC";
  commonHists2ndOrderName += fAnalysisLabel->Data();
  fCommonHists2nd = new AliFlowCommonHist(commonHists2ndOrderName.Data(),commonHists2ndOrderName.Data(),fBookOnlyBasicCCH);
  fHistList->Add(fCommonHists2nd);  
  // Common control histogram filled for events with 2 and more reference particles:
  TString commonHists4thOrderName = "AliFlowCommonHist4thOrderQC";
  commonHists4thOrderName += fAnalysisLabel->Data();
  fCommonHists4th = new AliFlowCommonHist(commonHists4thOrderName.Data(),commonHists4thOrderName.Data(),fBookOnlyBasicCCH);
  fHistList->Add(fCommonHists4th);  
  // Common control histogram filled for events with 6 and more reference particles:
  TString commonHists6thOrderName = "AliFlowCommonHist6thOrderQC";
  commonHists6thOrderName += fAnalysisLabel->Data();
  fCommonHists6th = new AliFlowCommonHist(commonHists6thOrderName.Data(),commonHists6thOrderName.Data(),fBookOnlyBasicCCH);
  fHistList->Add(fCommonHists6th);  
  // Common control histogram filled for events with 8 and more reference particles:
  TString commonHists8thOrderName = "AliFlowCommonHist8thOrderQC";
  commonHists8thOrderName += fAnalysisLabel->Data();
  fCommonHists8th = new AliFlowCommonHist(commonHists8thOrderName.Data(),commonHists8thOrderName.Data(),fBookOnlyBasicCCH);
  fHistList->Add(fCommonHists8th);    
 } // end of if(fFillMultipleControlHistograms)
 
 // b) Book common result histograms: 
 //  Common result histograms for QC{2}:
 TString commonHistResults2ndOrderName = "AliFlowCommonHistResults2ndOrderQC";
 commonHistResults2ndOrderName += fAnalysisLabel->Data();
 fCommonHistsResults2nd = new AliFlowCommonHistResults(commonHistResults2ndOrderName.Data(),"",fHarmonic);
 fHistList->Add(fCommonHistsResults2nd);  
 //  Common result histograms for QC{4}:
 TString commonHistResults4thOrderName = "AliFlowCommonHistResults4thOrderQC";
 commonHistResults4thOrderName += fAnalysisLabel->Data();
 fCommonHistsResults4th = new AliFlowCommonHistResults(commonHistResults4thOrderName.Data(),"",fHarmonic);
 fHistList->Add(fCommonHistsResults4th); 
 //  Common result histograms for QC{6}:
 TString commonHistResults6thOrderName = "AliFlowCommonHistResults6thOrderQC";
 commonHistResults6thOrderName += fAnalysisLabel->Data();
 fCommonHistsResults6th = new AliFlowCommonHistResults(commonHistResults6thOrderName.Data(),"",fHarmonic);
 fHistList->Add(fCommonHistsResults6th);  
 //  Common result histograms for QC{8}:
 TString commonHistResults8thOrderName = "AliFlowCommonHistResults8thOrderQC";
 commonHistResults8thOrderName += fAnalysisLabel->Data();
 fCommonHistsResults8th = new AliFlowCommonHistResults(commonHistResults8thOrderName.Data(),"",fHarmonic);
 fHistList->Add(fCommonHistsResults8th); 
 
} // end of void AliFlowAnalysisWithQCumulants::BookCommonHistograms()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::BookAndFillWeightsHistograms()
{
 // Book and fill histograms which hold phi, pt and eta weights.

 if(!fWeightsList)
 {
  printf("\n WARNING (QC): fWeightsList is NULL in AFAWQC::BAFWH() !!!! \n\n");
  exit(0);  
 }
    
 TString fUseParticleWeightsName = "fUseParticleWeightsQC";
 fUseParticleWeightsName += fAnalysisLabel->Data();
 fUseParticleWeights = new TProfile(fUseParticleWeightsName.Data(),"0 = particle weight not used, 1 = particle weight used ",4,0,4);
 fUseParticleWeights->SetLabelSize(0.06);
 (fUseParticleWeights->GetXaxis())->SetBinLabel(1,"w_{#phi}");
 (fUseParticleWeights->GetXaxis())->SetBinLabel(2,"w_{p_{T}}");
 (fUseParticleWeights->GetXaxis())->SetBinLabel(3,"w_{#eta}");
 (fUseParticleWeights->GetXaxis())->SetBinLabel(4,"w_{track}");
 fUseParticleWeights->Fill(0.5,(Int_t)fUsePhiWeights);
 fUseParticleWeights->Fill(1.5,(Int_t)fUsePtWeights);
 fUseParticleWeights->Fill(2.5,(Int_t)fUseEtaWeights);
 fUseParticleWeights->Fill(3.5,(Int_t)fUseTrackWeights);
 fWeightsList->Add(fUseParticleWeights); 
  
 if(fUsePhiWeights)
 {
  if(fWeightsList->FindObject("phi_weights"))
  {
   fPhiWeights = dynamic_cast<TH1F*>(fWeightsList->FindObject("phi_weights"));
   if(!fPhiWeights)
   {
    printf("\n WARNING (QC): fPhiWeights is NULL in AFAWQC::BAFWH() !!!!\n\n");
    exit(0);
   }
   if(TMath::Abs(fPhiWeights->GetBinWidth(1)-fPhiBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (QC): Inconsistent binning in histograms for phi-weights throughout the code."<<endl;
    cout<<endl;
    //exit(0);
   }
  } else 
    {
     cout<<"WARNING: fWeightsList->FindObject(\"phi_weights\") is NULL in AFAWQC::BAFWH() !!!!"<<endl;
     exit(0);
    }
 } // end of if(fUsePhiWeights)
 
 if(fUsePtWeights) 
 {
  if(fWeightsList->FindObject("pt_weights"))
  {
   fPtWeights = dynamic_cast<TH1D*>(fWeightsList->FindObject("pt_weights"));
   if(!fPtWeights)
   {
    printf("\n WARNING (QC): fPtWeights is NULL in AFAWQC::BAFWH() !!!!\n\n");
    exit(0);
   }
   if(TMath::Abs(fPtWeights->GetBinWidth(1)-fPtBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (QC): Inconsistent binning in histograms for pt-weights throughout the code."<<endl;
    cout<<endl;
    //exit(0);
   }
  } else 
    {
     cout<<"WARNING: fWeightsList->FindObject(\"pt_weights\") is NULL in AFAWQC::BAFWH() !!!!"<<endl;
     exit(0);
    }
 } // end of if(fUsePtWeights)    

 if(fUseEtaWeights) 
 {
  if(fWeightsList->FindObject("eta_weights"))
  {
   fEtaWeights = dynamic_cast<TH1D*>(fWeightsList->FindObject("eta_weights"));
   if(!fEtaWeights)
   {
    printf("\n WARNING (QC): fEtaWeights is NULL in AFAWQC::BAFWH() !!!!\n\n");
    exit(0);
   }
   if(TMath::Abs(fEtaWeights->GetBinWidth(1)-fEtaBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (QC): Inconsistent binning in histograms for eta-weights throughout the code."<<endl;
    cout<<endl;
    //exit(0);
   }
  } else 
    {
     cout<<"WARNING: fUseEtaWeights && fWeightsList->FindObject(\"eta_weights\") is NULL in AFAWQC::BAFWH() !!!!"<<endl;
     exit(0);
    }
 } // end of if(fUseEtaWeights)
 
} // end of AliFlowAnalysisWithQCumulants::BookAndFillWeightsHistograms()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::BookEverythingForIntegratedFlow()
{
 // Book all objects for integrated flow:
 //  a) Book profile to hold all flags for integrated flow;
 //  b) Book event-by-event quantities;
 //  c) Book profiles; // to be improved (comment)
 //  d) Book histograms holding the final results.
 
 TString sinCosFlag[2] = {"sin","cos"}; // to be improved (should I promote this to data members?)
 TString powerFlag[2] = {"linear","quadratic"}; // to be improved (should I promote this to data members?)
 
 // a) Book profile to hold all flags for integrated flow:
 TString intFlowFlagsName = "fIntFlowFlags";
 intFlowFlagsName += fAnalysisLabel->Data();
 fIntFlowFlags = new TProfile(intFlowFlagsName.Data(),"Flags for Integrated Flow",15,0,15);
 fIntFlowFlags->SetTickLength(-0.01,"Y");
 fIntFlowFlags->SetMarkerStyle(25);
 fIntFlowFlags->SetLabelSize(0.04);
 fIntFlowFlags->SetLabelOffset(0.02,"Y");
 fIntFlowFlags->GetXaxis()->SetBinLabel(1,"Particle Weights");
 fIntFlowFlags->GetXaxis()->SetBinLabel(2,"Event Weights");
 fIntFlowFlags->GetXaxis()->SetBinLabel(3,"Corrected for NUA?");
 fIntFlowFlags->GetXaxis()->SetBinLabel(4,"Print RF results");
 fIntFlowFlags->GetXaxis()->SetBinLabel(5,"Print RP results");
 fIntFlowFlags->GetXaxis()->SetBinLabel(6,"Print POI results");
 fIntFlowFlags->GetXaxis()->SetBinLabel(7,"Print RF (rebinned in M) results");
 fIntFlowFlags->GetXaxis()->SetBinLabel(8,"Corrected for NUA vs M?");
 fIntFlowFlags->GetXaxis()->SetBinLabel(9,"Propagate errors to v_{n} from correlations?");
 fIntFlowFlags->GetXaxis()->SetBinLabel(10,"Calculate cumulants vs M");
 fIntFlowFlags->GetXaxis()->SetBinLabel(11,"fMinimumBiasReferenceFlow");
 fIntFlowFlags->GetXaxis()->SetBinLabel(12,"fForgetAboutCovariances");
 fIntFlowFlags->GetXaxis()->SetBinLabel(13,"fStorePhiDistributionForOneEvent");
 fIntFlowFlags->GetXaxis()->SetBinLabel(14,"fFillMultipleControlHistograms");
 fIntFlowFlags->GetXaxis()->SetBinLabel(15,"Calculate all correlations vs M");
 fIntFlowList->Add(fIntFlowFlags);

 // b) Book event-by-event quantities:
 // Re[Q_{m*n,k}], Im[Q_{m*n,k}] and S_{p,k}^M: 
 fReQ  = new TMatrixD(6,9);
 fImQ  = new TMatrixD(6,9);
 fSpk = new TMatrixD(8,9);
 // average correlations <2>, <4>, <6> and <8> for single event (bining is the same as in fIntFlowCorrelationsPro and fIntFlowCorrelationsHist):
 TString intFlowCorrelationsEBEName = "fIntFlowCorrelationsEBE";
 intFlowCorrelationsEBEName += fAnalysisLabel->Data();
 fIntFlowCorrelationsEBE = new TH1D(intFlowCorrelationsEBEName.Data(),intFlowCorrelationsEBEName.Data(),4,0,4);
 // weights for average correlations <2>, <4>, <6> and <8> for single event:
 TString intFlowEventWeightsForCorrelationsEBEName = "fIntFlowEventWeightsForCorrelationsEBE";
 intFlowEventWeightsForCorrelationsEBEName += fAnalysisLabel->Data();
 fIntFlowEventWeightsForCorrelationsEBE = new TH1D(intFlowEventWeightsForCorrelationsEBEName.Data(),intFlowEventWeightsForCorrelationsEBEName.Data(),4,0,4);
 // average all correlations for single event (bining is the same as in fIntFlowCorrelationsAllPro and fIntFlowCorrelationsAllHist):
 TString intFlowCorrelationsAllEBEName = "fIntFlowCorrelationsAllEBE";
 intFlowCorrelationsAllEBEName += fAnalysisLabel->Data();
 fIntFlowCorrelationsAllEBE = new TH1D(intFlowCorrelationsAllEBEName.Data(),intFlowCorrelationsAllEBEName.Data(),64,0,64);
 // average correction terms for non-uniform acceptance for single event 
 // (binning is the same as in fIntFlowCorrectionTermsForNUAPro[2] and fIntFlowCorrectionTermsForNUAHist[2]):
 TString fIntFlowCorrectionTermsForNUAEBEName = "fIntFlowCorrectionTermsForNUAEBE";
 fIntFlowCorrectionTermsForNUAEBEName += fAnalysisLabel->Data();
 for(Int_t sc=0;sc<2;sc++) // sin or cos terms
 {
  fIntFlowCorrectionTermsForNUAEBE[sc] = new TH1D(Form("%s: %s terms",fIntFlowCorrectionTermsForNUAEBEName.Data(),sinCosFlag[sc].Data()),Form("Correction terms for non-uniform acceptance (%s terms)",sinCosFlag[sc].Data()),4,0,4);  
 }
 // event weights for terms for non-uniform acceptance: 
 TString fIntFlowEventWeightForCorrectionTermsForNUAEBEName = "fIntFlowEventWeightForCorrectionTermsForNUAEBE";
 fIntFlowEventWeightForCorrectionTermsForNUAEBEName += fAnalysisLabel->Data();
 for(Int_t sc=0;sc<2;sc++) // sin or cos terms
 {
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[sc] = new TH1D(Form("%s: %s terms",fIntFlowEventWeightForCorrectionTermsForNUAEBEName.Data(),sinCosFlag[sc].Data()),Form("Event weights for terms for non-uniform acceptance (%s terms)",sinCosFlag[sc].Data()),4,0,4); // to be improved - 4  
 }
 // c) Book profiles: // to be improved (comment)
 // profile to hold average multiplicities and number of events for events with nRP>=0, nRP>=1, ... , and nRP>=8:
 TString avMultiplicityName = "fAvMultiplicity";
 avMultiplicityName += fAnalysisLabel->Data();
 fAvMultiplicity = new TProfile(avMultiplicityName.Data(),"Average multiplicities of reference particles (RPs)",9,0,9);
 fAvMultiplicity->SetTickLength(-0.01,"Y");
 fAvMultiplicity->SetMarkerStyle(25);
 fAvMultiplicity->SetLabelSize(0.05);
 fAvMultiplicity->SetLabelOffset(0.02,"Y");
 fAvMultiplicity->SetYTitle("Average multiplicity");
 (fAvMultiplicity->GetXaxis())->SetBinLabel(1,"all evts");
 (fAvMultiplicity->GetXaxis())->SetBinLabel(2,"n_{RP} #geq 1");
 (fAvMultiplicity->GetXaxis())->SetBinLabel(3,"n_{RP} #geq 2");
 (fAvMultiplicity->GetXaxis())->SetBinLabel(4,"n_{RP} #geq 3");
 (fAvMultiplicity->GetXaxis())->SetBinLabel(5,"n_{RP} #geq 4");
 (fAvMultiplicity->GetXaxis())->SetBinLabel(6,"n_{RP} #geq 5");
 (fAvMultiplicity->GetXaxis())->SetBinLabel(7,"n_{RP} #geq 6");
 (fAvMultiplicity->GetXaxis())->SetBinLabel(8,"n_{RP} #geq 7");
 (fAvMultiplicity->GetXaxis())->SetBinLabel(9,"n_{RP} #geq 8");
 fIntFlowProfiles->Add(fAvMultiplicity);
 // Average correlations <<2>>, <<4>>, <<6>> and <<8>> for all events (with wrong errors!):
 TString correlationFlag[4] = {"#LT#LT2#GT#GT","#LT#LT4#GT#GT","#LT#LT6#GT#GT","#LT#LT8#GT#GT"};
 TString intFlowCorrelationsProName = "fIntFlowCorrelationsPro";
 intFlowCorrelationsProName += fAnalysisLabel->Data();
 fIntFlowCorrelationsPro = new TProfile(intFlowCorrelationsProName.Data(),"Average correlations for all events",4,0,4,"s");
 fIntFlowCorrelationsPro->Sumw2();
 fIntFlowCorrelationsPro->SetTickLength(-0.01,"Y");
 fIntFlowCorrelationsPro->SetMarkerStyle(25);
 fIntFlowCorrelationsPro->SetLabelSize(0.06);
 fIntFlowCorrelationsPro->SetLabelOffset(0.01,"Y");
 for(Int_t b=0;b<4;b++)
 {
  (fIntFlowCorrelationsPro->GetXaxis())->SetBinLabel(b+1,correlationFlag[b].Data());
 }
 fIntFlowProfiles->Add(fIntFlowCorrelationsPro);
 // Average correlations squared <<2>^2>, <<4>^2>, <<6>^2> and <<8>^2> for all events:
 TString squaredCorrelationFlag[4] = {"#LT#LT2#GT^{2}#GT","#LT#LT4#GT^{2}#GT","#LT#LT6#GT^{2}#GT","#LT#LT8#GT^{2}#GT"};
 TString intFlowSquaredCorrelationsProName = "fIntFlowSquaredCorrelationsPro";
 intFlowSquaredCorrelationsProName += fAnalysisLabel->Data();
 fIntFlowSquaredCorrelationsPro = new TProfile(intFlowSquaredCorrelationsProName.Data(),"Average squared correlations for all events",4,0,4,"s");
 fIntFlowSquaredCorrelationsPro->Sumw2();
 fIntFlowSquaredCorrelationsPro->SetTickLength(-0.01,"Y");
 fIntFlowSquaredCorrelationsPro->SetMarkerStyle(25);
 fIntFlowSquaredCorrelationsPro->SetLabelSize(0.06);
 fIntFlowSquaredCorrelationsPro->SetLabelOffset(0.01,"Y");
 for(Int_t b=0;b<4;b++)
 {
  (fIntFlowSquaredCorrelationsPro->GetXaxis())->SetBinLabel(b+1,squaredCorrelationFlag[b].Data());
 }
 fIntFlowProfiles->Add(fIntFlowSquaredCorrelationsPro);
 if(fCalculateCumulantsVsM)
 {
  for(Int_t ci=0;ci<4;ci++) // correlation index
  {
   // average correlations <<2>>, <<4>>, <<6>> and <<8>> versus multiplicity for all events (with wrong errors):
   TString intFlowCorrelationsVsMProName = "fIntFlowCorrelationsVsMPro";
   intFlowCorrelationsVsMProName += fAnalysisLabel->Data();
   fIntFlowCorrelationsVsMPro[ci] = new TProfile(Form("%s, %s",intFlowCorrelationsVsMProName.Data(),correlationFlag[ci].Data()),
                                                 Form("%s vs multiplicity",correlationFlag[ci].Data()),
                                                 fnBinsMult,fMinMult,fMaxMult,"s");   
   fIntFlowCorrelationsVsMPro[ci]->Sumw2();                                                                                       
   fIntFlowCorrelationsVsMPro[ci]->GetYaxis()->SetTitle(correlationFlag[ci].Data());
   fIntFlowCorrelationsVsMPro[ci]->GetXaxis()->SetTitle("M");
   fIntFlowProfiles->Add(fIntFlowCorrelationsVsMPro[ci]);
   // average squared correlations <<2>^2>, <<4>^2>, <<6>^2> and <<8>^2> versus multiplicity for all events:  
   TString intFlowSquaredCorrelationsVsMProName = "fIntFlowSquaredCorrelationsVsMPro";
   intFlowSquaredCorrelationsVsMProName += fAnalysisLabel->Data();
   fIntFlowSquaredCorrelationsVsMPro[ci] = new TProfile(Form("%s, %s",intFlowSquaredCorrelationsVsMProName.Data(),squaredCorrelationFlag[ci].Data()),
                                                        Form("%s vs multiplicity",squaredCorrelationFlag[ci].Data()),
                                                        fnBinsMult,fMinMult,fMaxMult,"s");   
   fIntFlowSquaredCorrelationsVsMPro[ci]->Sumw2();                                                                                              
   fIntFlowSquaredCorrelationsVsMPro[ci]->GetYaxis()->SetTitle(squaredCorrelationFlag[ci].Data());
   fIntFlowSquaredCorrelationsVsMPro[ci]->GetXaxis()->SetTitle("M");
   fIntFlowProfiles->Add(fIntFlowSquaredCorrelationsVsMPro[ci]);
  } // end of for(Int_t ci=0;ci<4;ci++) // correlation index  
 } // end of if(fCalculateCumulantsVsM)
 // averaged all correlations for all events (with wrong errors!):
 TString intFlowCorrelationsAllProName = "fIntFlowCorrelationsAllPro";
 intFlowCorrelationsAllProName += fAnalysisLabel->Data();
 fIntFlowCorrelationsAllPro = new TProfile(intFlowCorrelationsAllProName.Data(),"Average all correlations for all events",64,0,64);
 fIntFlowCorrelationsAllPro->Sumw2();
 fIntFlowCorrelationsAllPro->SetTickLength(-0.01,"Y");
 fIntFlowCorrelationsAllPro->SetMarkerStyle(25);
 fIntFlowCorrelationsAllPro->SetLabelSize(0.03);
 fIntFlowCorrelationsAllPro->SetLabelOffset(0.01,"Y");
 // 2-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(1,"#LT#LT2#GT#GT_{n|n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(2,"#LT#LT2#GT#GT_{2n|2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(3,"#LT#LT2#GT#GT_{3n|3n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(4,"#LT#LT2#GT#GT_{4n|4n}");
 // 3-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(6,"#LT#LT3#GT#GT_{2n|n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(7,"#LT#LT3#GT#GT_{3n|2n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(8,"#LT#LT3#GT#GT_{4n|2n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(9,"#LT#LT3#GT#GT_{4n|3n,n}");
 // 4-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(11,"#LT#LT4#GT#GT_{n,n|n,n}"); 
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(12,"#LT#LT4#GT#GT_{2n,n|2n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(13,"#LT#LT4#GT#GT_{2n,2n|2n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(14,"#LT#LT4#GT#GT_{3n|n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(15,"#LT#LT4#GT#GT_{3n,n|3n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(16,"#LT#LT4#GT#GT_{3n,n|2n,2n}"); 
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(17,"#LT#LT4#GT#GT_{4n|2n,n,n}");
 // 5-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(19,"#LT#LT5#GT#GT_{2n,n|n,n,n}"); 
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(20,"#LT#LT5#GT#GT_{2n,2n|2n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(21,"#LT#LT5#GT#GT_{3n,n|2n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(22,"#LT#LT5#GT#GT_{4n|n,n,n,n}");
 // 6-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(24,"#LT#LT6#GT#GT_{n,n,n|n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(25,"#LT#LT6#GT#GT_{2n,n,n|2n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(26,"#LT#LT6#GT#GT_{2n,2n|n,n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(27,"#LT#LT6#GT#GT_{3n,n|n,n,n,n}");
 // 7-p correlations:  
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(29,"#LT#LT7#GT#GT_{2n,n,n|n,n,n,n}");
 // 8-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(31,"#LT#LT8#GT#GT_{n,n,n,n|n,n,n,n}");
 //  EXTRA correlations for v3{5} study:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(33,"#LT#LT4#GT#GT_{4n,2n|3n,3n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(34,"#LT#LT5#GT#GT_{3n,3n|2n,2n,2n}");
 //  EXTRA correlations for Teaney-Yan study:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(35,"#LT#LT2#GT#GT_{5n|5n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(36,"#LT#LT2#GT#GT_{6n|6n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(37,"#LT#LT3#GT#GT_{5n|3n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(38,"#LT#LT3#GT#GT_{5n|4n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(39,"#LT#LT3#GT#GT_{6n|3n,3n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(40,"#LT#LT3#GT#GT_{6n|4n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(41,"#LT#LT3#GT#GT_{6n|5n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(42,"#LT#LT4#GT#GT_{6n|3n,2n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(43,"#LT#LT4#GT#GT_{3n,2n|3n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(44,"#LT#LT4#GT#GT_{4n,1n|3n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(45,"#LT#LT4#GT#GT_{3n,3n|3n,3n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(46,"#LT#LT4#GT#GT_{4n,2n|3n,3n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(47,"#LT#LT4#GT#GT_{5n,1n|3n,3n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(48,"#LT#LT4#GT#GT_{4n,2n|4n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(49,"#LT#LT4#GT#GT_{5n,1n|4n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(50,"#LT#LT4#GT#GT_{5n|3n,1n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(51,"#LT#LT4#GT#GT_{5n|2n,2n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(52,"#LT#LT4#GT#GT_{5n,1n|5n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(53,"#LT#LT5#GT#GT_{3n,3n|3n,2n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(54,"#LT#LT5#GT#GT_{4n,2n|3n,2n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(55,"#LT#LT5#GT#GT_{3n,2n|3n,1n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(56,"#LT#LT5#GT#GT_{3n,2n|2n,2n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(57,"#LT#LT5#GT#GT_{5n,1n|3n,2n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(58,"#LT#LT6#GT#GT_{3n,2n,1n|3n,2n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(59,"#LT#LT4#GT#GT_{6n|4n,1n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(60,"#LT#LT4#GT#GT_{6n|2n,2n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(61,"#LT#LT5#GT#GT_{6n|2n,2n,1n,1n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(62,"#LT#LT5#GT#GT_{4n,1n,1n|3n,3n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(63,"#LT#LT6#GT#GT_{3n,3n|2n,2n,1n,1n}");
 fIntFlowProfiles->Add(fIntFlowCorrelationsAllPro);
 // average all correlations versus multiplicity (errors via Sumw2 - to be improved):
 if(fCalculateAllCorrelationsVsM)
 {
  // 2-p correlations vs M:  
  fIntFlowCorrelationsAllVsMPro[0] = new TProfile("two1n1n","#LT#LT2#GT#GT_{n|n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[0]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[0]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[0]);  
  fIntFlowCorrelationsAllVsMPro[1] = new TProfile("two2n2n","#LT#LT2#GT#GT_{2n|2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[1]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[1]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[1]);
  fIntFlowCorrelationsAllVsMPro[2] = new TProfile("two3n3n","#LT#LT2#GT#GT_{3n|3n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[2]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[2]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[2]);
  fIntFlowCorrelationsAllVsMPro[3] = new TProfile("two4n4n","#LT#LT2#GT#GT_{4n|4n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[3]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[3]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[3]);
  // 3-p correlations vs M:
  fIntFlowCorrelationsAllVsMPro[5] = new TProfile("three2n1n1n","#LT#LT3#GT#GT_{2n|n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[5]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[5]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[5]);
  fIntFlowCorrelationsAllVsMPro[6] = new TProfile("three3n2n1n","#LT#LT3#GT#GT_{3n|2n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[6]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[6]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[6]);
  fIntFlowCorrelationsAllVsMPro[7] = new TProfile("three4n2n2n","#LT#LT3#GT#GT_{4n|2n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[7]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[7]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[7]);
  fIntFlowCorrelationsAllVsMPro[8] = new TProfile("three4n3n1n","#LT#LT3#GT#GT_{4n|3n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[8]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[8]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[8]);
  // 4-p correlations vs M:
  fIntFlowCorrelationsAllVsMPro[10] = new TProfile("four1n1n1n1n","#LT#LT4#GT#GT_{n,n|n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[10]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[10]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[10]);
  fIntFlowCorrelationsAllVsMPro[11] = new TProfile("four2n1n2n1n","#LT#LT4#GT#GT_{2n,n|2n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[11]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[11]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[11]);
  fIntFlowCorrelationsAllVsMPro[12] = new TProfile("four2n2n2n2n","#LT#LT4#GT#GT_{2n,2n|2n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[12]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[12]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[12]);
  fIntFlowCorrelationsAllVsMPro[13] = new TProfile("four3n1n1n1n","#LT#LT4#GT#GT_{3n|n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[13]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[13]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[13]);
  fIntFlowCorrelationsAllVsMPro[14] = new TProfile("four3n1n3n1n","#LT#LT4#GT#GT_{3n,n|3n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[14]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[14]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[14]);
  fIntFlowCorrelationsAllVsMPro[15] = new TProfile("four3n1n2n2n","#LT#LT4#GT#GT_{3n,n|2n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[15]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[15]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[15]);
  fIntFlowCorrelationsAllVsMPro[16] = new TProfile("four4n2n1n1n","#LT#LT4#GT#GT_{4n|2n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[16]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[16]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[16]);
  // 5-p correlations vs M:
  fIntFlowCorrelationsAllVsMPro[18] = new TProfile("five2n1n1n1n1n","#LT#LT5#GT#GT_{2n,n|n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[18]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[18]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[18]);
  fIntFlowCorrelationsAllVsMPro[19] = new TProfile("five2n2n2n1n1n","#LT#LT5#GT#GT_{2n,2n|2n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[19]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[19]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[19]);
  fIntFlowCorrelationsAllVsMPro[20] = new TProfile("five3n1n2n1n1n","#LT#LT5#GT#GT_{3n,n|2n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[20]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[20]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[20]);
  fIntFlowCorrelationsAllVsMPro[21] = new TProfile("five4n1n1n1n1n","#LT#LT5#GT#GT_{4n|n,n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[21]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[21]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[21]);
  // 6-p correlations vs M:
  fIntFlowCorrelationsAllVsMPro[23] = new TProfile("six1n1n1n1n1n1n","#LT#LT6#GT#GT_{n,n,n|n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[23]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[23]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[23]);
  fIntFlowCorrelationsAllVsMPro[24] = new TProfile("six2n1n1n2n1n1n","#LT#LT6#GT#GT_{2n,n,n|2n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[24]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[24]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[24]);
  fIntFlowCorrelationsAllVsMPro[25] = new TProfile("six2n2n1n1n1n1n","#LT#LT6#GT#GT_{2n,2n|n,n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[25]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[25]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[25]);
  fIntFlowCorrelationsAllVsMPro[26] = new TProfile("six3n1n1n1n1n1n","#LT#LT6#GT#GT_{3n,n|n,n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[26]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[26]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[26]);
  // 7-p correlations vs M:
  fIntFlowCorrelationsAllVsMPro[28] = new TProfile("seven2n1n1n1n1n1n1n","#LT#LT7#GT#GT_{2n,n,n|n,n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[28]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[28]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[28]);
  // 8-p correlations vs M:
  fIntFlowCorrelationsAllVsMPro[30] = new TProfile("eight1n1n1n1n1n1n1n1n","#LT#LT8#GT#GT_{n,n,n,n|n,n,n,n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[30]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[30]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[30]);
  // EXTRA correlations vs M for v3{5} study (to be improved - put them in a right order somewhere):
  fIntFlowCorrelationsAllVsMPro[32] = new TProfile("four4n2n3n3n","#LT#LT4#GT#GT_{4n,2n|3n,3n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[32]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[32]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[32]);
  fIntFlowCorrelationsAllVsMPro[33] = new TProfile("five3n3n2n2n2n","#LT#LT5#GT#GT_{3n,3n|2n,2n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[33]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[33]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[33]);
  // EXTRA correlations vs M for Teaney-Yan study (to be improved - put them in a right order somewhere):
  fIntFlowCorrelationsAllVsMPro[34] = new TProfile("two5n5n","#LT#LT2#GT#GT_{5n|5n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[34]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[34]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[34]);  
  fIntFlowCorrelationsAllVsMPro[35] = new TProfile("two6n6n","#LT#LT2#GT#GT_{6n|6n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[35]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[35]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[35]);  
  fIntFlowCorrelationsAllVsMPro[36] = new TProfile("three5n3n2n","#LT#LT3#GT#GT_{5n|3n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[36]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[36]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[36]); 
  fIntFlowCorrelationsAllVsMPro[37] = new TProfile("three5n4n1n","#LT#LT3#GT#GT_{5n|4n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[37]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[37]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[37]);
  fIntFlowCorrelationsAllVsMPro[38] = new TProfile("three6n3n3n","#LT#LT3#GT#GT_{6n|3n,3n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[38]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[38]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[38]);  
  fIntFlowCorrelationsAllVsMPro[39] = new TProfile("three6n4n2n","#LT#LT3#GT#GT_{6n|4n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[39]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[39]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[39]);
  fIntFlowCorrelationsAllVsMPro[40] = new TProfile("three6n5n1n","#LT#LT3#GT#GT_{6n|5n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[40]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[40]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[40]);
  fIntFlowCorrelationsAllVsMPro[41] = new TProfile("four6n3n2n1n","#LT#LT4#GT#GT_{6n|3n,2n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[41]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[41]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[41]);
  fIntFlowCorrelationsAllVsMPro[42] = new TProfile("four3n2n3n2n","#LT#LT4#GT#GT_{3n,2n|3n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[42]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[42]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[42]);
  fIntFlowCorrelationsAllVsMPro[43] = new TProfile("four4n1n3n2n","#LT#LT4#GT#GT_{4n,1n|3n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[43]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[43]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[43]);
  fIntFlowCorrelationsAllVsMPro[44] = new TProfile("four3n3n3n3n","#LT#LT4#GT#GT_{3n,3n|3n,3n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[44]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[44]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[44]);
  fIntFlowCorrelationsAllVsMPro[45] = new TProfile("four4n2n3n3n","#LT#LT4#GT#GT_{4n,2n|3n,3n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[45]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[45]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[45]);
  fIntFlowCorrelationsAllVsMPro[46] = new TProfile("four5n1n3n3n","#LT#LT4#GT#GT_{5n,1n|3n,3n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[46]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[46]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[46]);
  fIntFlowCorrelationsAllVsMPro[47] = new TProfile("four4n2n4n2n","#LT#LT4#GT#GT_{4n,2n|4n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[47]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[47]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[47]);
  fIntFlowCorrelationsAllVsMPro[48] = new TProfile("four5n1n4n2n","#LT#LT4#GT#GT_{5n,1n|4n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[48]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[48]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[48]);
  fIntFlowCorrelationsAllVsMPro[49] = new TProfile("four5n3n1n1n","#LT#LT4#GT#GT_{5n|3n,1n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[49]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[49]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[49]);
  fIntFlowCorrelationsAllVsMPro[50] = new TProfile("four5n2n2n1n","#LT#LT4#GT#GT_{5n|2n,2n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[50]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[50]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[50]);
  fIntFlowCorrelationsAllVsMPro[51] = new TProfile("four5n1n5n1n","#LT#LT4#GT#GT_{5n,1n|5n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[51]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[51]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[51]);
  fIntFlowCorrelationsAllVsMPro[52] = new TProfile("five3n3n3n2n1n","#LT#LT5#GT#GT_{3n,3n|3n,2n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[52]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[52]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[52]);
  fIntFlowCorrelationsAllVsMPro[53] = new TProfile("five4n2n3n2n1n","#LT#LT5#GT#GT_{4n,2n|3n,2n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[53]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[53]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[53]);
  fIntFlowCorrelationsAllVsMPro[54] = new TProfile("five3n2n3n1n1n","#LT#LT5#GT#GT_{3n,2n|3n,1n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[54]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[54]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[54]);
  fIntFlowCorrelationsAllVsMPro[55] = new TProfile("five3n2n2n2n1n","#LT#LT5#GT#GT_{3n,2n|2n,2n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[55]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[55]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[55]);
  fIntFlowCorrelationsAllVsMPro[56] = new TProfile("five5n1n3n2n1n","#LT#LT5#GT#GT_{5n,1n|3n,2n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[56]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[56]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[56]);
  fIntFlowCorrelationsAllVsMPro[57] = new TProfile("six3n2n1n3n2n1n","#LT#LT6#GT#GT_{3n,2n,1n|3n,2n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[57]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[57]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[57]);  
  fIntFlowCorrelationsAllVsMPro[58] = new TProfile("four6n4n1n1n","#LT#LT4#GT#GT_{6n|4n,1n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[58]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[58]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[58]);
  fIntFlowCorrelationsAllVsMPro[59] = new TProfile("four6n2n2n2n","#LT#LT4#GT#GT_{6n|2n,2n,2n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[59]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[59]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[59]);
  fIntFlowCorrelationsAllVsMPro[60] = new TProfile("five6n2n2n1n1n","#LT#LT5#GT#GT_{6n|2n,2n,1n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[60]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[60]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[60]);
  fIntFlowCorrelationsAllVsMPro[61] = new TProfile("five4n1n1n3n3n","#LT#LT5#GT#GT_{4n,1n,1n|3n,3n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[61]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[61]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[61]);
  fIntFlowCorrelationsAllVsMPro[62] = new TProfile("six3n3n2n2n1n1n","#LT#LT6#GT#GT_{3n,3n|2n,2n,1n,1n}",fnBinsMult,fMinMult,fMaxMult);
  fIntFlowCorrelationsAllVsMPro[62]->Sumw2();
  fIntFlowCorrelationsAllVsMPro[62]->GetXaxis()->SetTitle("M");  
  fIntFlowAllCorrelationsVsM->Add(fIntFlowCorrelationsAllVsMPro[62]);
 } // end of if(fCalculateAllCorrelationsVsM)
 // when particle weights are used some extra correlations appear:
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights) 
 {
  TString intFlowExtraCorrelationsProName = "fIntFlowExtraCorrelationsPro";
  intFlowExtraCorrelationsProName += fAnalysisLabel->Data();
  fIntFlowExtraCorrelationsPro = new TProfile(intFlowExtraCorrelationsProName.Data(),"Average extra correlations for all events",100,0,100,"s");
  fIntFlowExtraCorrelationsPro->SetTickLength(-0.01,"Y");
  fIntFlowExtraCorrelationsPro->SetMarkerStyle(25);
  fIntFlowExtraCorrelationsPro->SetLabelSize(0.03);
  fIntFlowExtraCorrelationsPro->SetLabelOffset(0.01,"Y");
  // extra 2-p correlations:
  (fIntFlowExtraCorrelationsPro->GetXaxis())->SetBinLabel(1,"<<w1^3 w2 cos(n*(phi1-phi2))>>");
  (fIntFlowExtraCorrelationsPro->GetXaxis())->SetBinLabel(2,"<<w1 w2 w3^2 cos(n*(phi1-phi2))>>");
  fIntFlowProfiles->Add(fIntFlowExtraCorrelationsPro);
 } // end of if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)
 // average product of correlations <2>, <4>, <6> and <8>:  
 TString productFlag[6] = {"#LT#LT2#GT#LT4#GT#GT","#LT#LT2#GT#LT6#GT#GT","#LT#LT2#GT#LT8#GT#GT",
                           "#LT#LT4#GT#LT6#GT#GT","#LT#LT4#GT#LT8#GT#GT","#LT#LT6#GT#LT8#GT#GT"};
 TString intFlowProductOfCorrelationsProName = "fIntFlowProductOfCorrelationsPro";
 intFlowProductOfCorrelationsProName += fAnalysisLabel->Data();
 fIntFlowProductOfCorrelationsPro = new TProfile(intFlowProductOfCorrelationsProName.Data(),"Average products of correlations",6,0,6);
 fIntFlowProductOfCorrelationsPro->SetTickLength(-0.01,"Y");
 fIntFlowProductOfCorrelationsPro->SetMarkerStyle(25); 
 fIntFlowProductOfCorrelationsPro->SetLabelSize(0.05);
 fIntFlowProductOfCorrelationsPro->SetLabelOffset(0.01,"Y");
 for(Int_t b=0;b<6;b++)
 {
  (fIntFlowProductOfCorrelationsPro->GetXaxis())->SetBinLabel(b+1,productFlag[b].Data());
 }
 fIntFlowProfiles->Add(fIntFlowProductOfCorrelationsPro); 
 // average product of correlations <2>, <4>, <6> and <8> versus multiplicity
 // [0=<<2><4>>,1=<<2><6>>,2=<<2><8>>,3=<<4><6>>,4=<<4><8>>,5=<<6><8>>]  
 if(fCalculateCumulantsVsM)
 {
  TString intFlowProductOfCorrelationsVsMProName = "fIntFlowProductOfCorrelationsVsMPro";
  intFlowProductOfCorrelationsVsMProName += fAnalysisLabel->Data();
  for(Int_t pi=0;pi<6;pi++)
  { 
   fIntFlowProductOfCorrelationsVsMPro[pi] = new TProfile(Form("%s, %s",intFlowProductOfCorrelationsVsMProName.Data(),productFlag[pi].Data()),
                                                          Form("%s versus multiplicity",productFlag[pi].Data()),
                                                          fnBinsMult,fMinMult,fMaxMult);             
   fIntFlowProductOfCorrelationsVsMPro[pi]->GetXaxis()->SetTitle("M");
   fIntFlowProfiles->Add(fIntFlowProductOfCorrelationsVsMPro[pi]);
  } // end of for(Int_t pi=0;pi<6;pi++)
 } // end of if(fCalculateCumulantsVsM) 
 // average product of correction terms for NUA:  
 TString intFlowProductOfCorrectionTermsForNUAProName = "fIntFlowProductOfCorrectionTermsForNUAPro";
 intFlowProductOfCorrectionTermsForNUAProName += fAnalysisLabel->Data();
 fIntFlowProductOfCorrectionTermsForNUAPro = new TProfile(intFlowProductOfCorrectionTermsForNUAProName.Data(),"Average products of correction terms for NUA",27,0,27);
 fIntFlowProductOfCorrectionTermsForNUAPro->SetTickLength(-0.01,"Y");
 fIntFlowProductOfCorrectionTermsForNUAPro->SetMarkerStyle(25); 
 fIntFlowProductOfCorrectionTermsForNUAPro->SetLabelSize(0.03);
 fIntFlowProductOfCorrectionTermsForNUAPro->SetLabelOffset(0.01,"Y");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(1,"<<2><cos(#phi)>>");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(2,"<<2><sin(#phi)>>");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(3,"<<cos(#phi)><sin(#phi)>>");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(4,"Cov(<2>,<cos(#phi_{1}+#phi_{2})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(5,"Cov(<2>,<sin(#phi_{1}+#phi_{2})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(6,"Cov(<2>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(7,"Cov(<2>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(8,"Cov(<4>,<cos(#phi)>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(9,"Cov(<4>,<sin(#phi)>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(10,"Cov(<4>,<cos(#phi_{1}+#phi_{2})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(11,"Cov(<4>,<sin(#phi_{1}+#phi_{2})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(12,"Cov(<4>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(13,"Cov(<4>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(14,"Cov(<cos(#phi)>,<cos(#phi_{1}+#phi_{2})>)"); 
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(15,"Cov(<cos(#phi)>,<sin(#phi_{1}+#phi_{2})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(16,"Cov(<cos(#phi)>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(17,"Cov(<cos(#phi)>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(18,"Cov(<sin(#phi)>,<cos(#phi_{1}+#phi_{2})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(19,"Cov(<sin(#phi)>,<sin(#phi_{1}+#phi_{2})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(20,"Cov(<sin(#phi)>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(21,"Cov(<sin(#phi)>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(22,"Cov(<cos(#phi_{1}+#phi_{2})>,<sin(#phi_{1}+#phi_{2})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(23,"Cov(<cos(#phi_{1}+#phi_{2})>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(24,"Cov(<cos(#phi_{1}+#phi_{2})>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(25,"Cov(<sin(#phi_{1}+#phi_{2})>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(26,"Cov(<sin(#phi_{1}+#phi_{2})>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowProductOfCorrectionTermsForNUAPro->GetXaxis())->SetBinLabel(27,"Cov(<cos(#phi_{1}-#phi_{2}-#phi_{3}>,<sin(#phi_{1}-#phi_{2}-#phi_{3}>)");
 fIntFlowProfiles->Add(fIntFlowProductOfCorrectionTermsForNUAPro);
 // average correction terms for non-uniform acceptance (with wrong errors!):
 for(Int_t sc=0;sc<2;sc++) // sin or cos terms
 {
  TString intFlowCorrectionTermsForNUAProName = "fIntFlowCorrectionTermsForNUAPro";
  intFlowCorrectionTermsForNUAProName += fAnalysisLabel->Data();
  fIntFlowCorrectionTermsForNUAPro[sc] = new TProfile(Form("%s: %s terms",intFlowCorrectionTermsForNUAProName.Data(),sinCosFlag[sc].Data()),Form("Correction terms for non-uniform acceptance (%s terms)",sinCosFlag[sc].Data()),4,0,4,"s");
  fIntFlowCorrectionTermsForNUAPro[sc]->SetTickLength(-0.01,"Y");
  fIntFlowCorrectionTermsForNUAPro[sc]->SetMarkerStyle(25);
  fIntFlowCorrectionTermsForNUAPro[sc]->SetLabelSize(0.05);
  fIntFlowCorrectionTermsForNUAPro[sc]->SetLabelOffset(0.01,"Y");
  (fIntFlowCorrectionTermsForNUAPro[sc]->GetXaxis())->SetBinLabel(1,Form("#LT#LT%s(n(#phi_{1}))#GT#GT",sinCosFlag[sc].Data()));
  (fIntFlowCorrectionTermsForNUAPro[sc]->GetXaxis())->SetBinLabel(2,Form("#LT#LT%s(n(#phi_{1}+#phi_{2}))#GT#GT",sinCosFlag[sc].Data()));  
  (fIntFlowCorrectionTermsForNUAPro[sc]->GetXaxis())->SetBinLabel(3,Form("#LT#LT%s(n(#phi_{1}-#phi_{2}-#phi_{3}))#GT#GT",sinCosFlag[sc].Data()));  
  (fIntFlowCorrectionTermsForNUAPro[sc]->GetXaxis())->SetBinLabel(4,Form("#LT#LT%s(n(2#phi_{1}-#phi_{2}))#GT#GT",sinCosFlag[sc].Data()));  
  fIntFlowProfiles->Add(fIntFlowCorrectionTermsForNUAPro[sc]);
  // versus multiplicity:
  if(fCalculateCumulantsVsM)
  {
   TString correctionTermFlag[4] = {"(n(phi1))","(n(phi1+phi2))","(n(phi1-phi2-phi3))","(n(2phi1-phi2))"}; // to be improved - hardwired 4
   for(Int_t ci=0;ci<4;ci++) // correction term index (to be improved - hardwired 4)
   {
    TString intFlowCorrectionTermsForNUAVsMProName = "fIntFlowCorrectionTermsForNUAVsMPro";
    intFlowCorrectionTermsForNUAVsMProName += fAnalysisLabel->Data();
    fIntFlowCorrectionTermsForNUAVsMPro[sc][ci] = new TProfile(Form("%s: #LT#LT%s%s#GT#GT",intFlowCorrectionTermsForNUAVsMProName.Data(),sinCosFlag[sc].Data(),correctionTermFlag[ci].Data()),Form("#LT#LT%s%s#GT#GT vs M",sinCosFlag[sc].Data(),correctionTermFlag[ci].Data()),fnBinsMult,fMinMult,fMaxMult,"s");
    fIntFlowProfiles->Add(fIntFlowCorrectionTermsForNUAVsMPro[sc][ci]);
   }
  } // end of if(fCalculateCumulantsVsM)
 } // end of for(Int_t sc=0;sc<2;sc++) 
 
 // d) Book histograms holding the final results:
 // average correlations <<2>>, <<4>>, <<6>> and <<8>> for all events (with correct errors!):
 TString intFlowCorrelationsHistName = "fIntFlowCorrelationsHist";
 intFlowCorrelationsHistName += fAnalysisLabel->Data();
 fIntFlowCorrelationsHist = new TH1D(intFlowCorrelationsHistName.Data(),"Average correlations for all events",4,0,4);
 fIntFlowCorrelationsHist->SetTickLength(-0.01,"Y");
 fIntFlowCorrelationsHist->SetMarkerStyle(25);
 fIntFlowCorrelationsHist->SetLabelSize(0.06);
 fIntFlowCorrelationsHist->SetLabelOffset(0.01,"Y");
 (fIntFlowCorrelationsHist->GetXaxis())->SetBinLabel(1,"#LT#LT2#GT#GT");
 (fIntFlowCorrelationsHist->GetXaxis())->SetBinLabel(2,"#LT#LT4#GT#GT");
 (fIntFlowCorrelationsHist->GetXaxis())->SetBinLabel(3,"#LT#LT6#GT#GT");
 (fIntFlowCorrelationsHist->GetXaxis())->SetBinLabel(4,"#LT#LT8#GT#GT");
 fIntFlowResults->Add(fIntFlowCorrelationsHist);
 // average correlations <<2>>, <<4>>, <<6>> and <<8>> for all events (with correct errors!) vs M:
 if(fCalculateCumulantsVsM)
 {
  for(Int_t ci=0;ci<4;ci++) // correlation index
  {
   TString intFlowCorrelationsVsMHistName = "fIntFlowCorrelationsVsMHist";
   intFlowCorrelationsVsMHistName += fAnalysisLabel->Data();
   fIntFlowCorrelationsVsMHist[ci] = new TH1D(Form("%s, %s",intFlowCorrelationsVsMHistName.Data(),correlationFlag[ci].Data()),
                                              Form("%s vs multiplicity",correlationFlag[ci].Data()),
                                              fnBinsMult,fMinMult,fMaxMult);                                            
   fIntFlowCorrelationsVsMHist[ci]->GetYaxis()->SetTitle(correlationFlag[ci].Data());
   fIntFlowCorrelationsVsMHist[ci]->GetXaxis()->SetTitle("M");
   fIntFlowResults->Add(fIntFlowCorrelationsVsMHist[ci]);
  } // end of for(Int_t ci=0;ci<4;ci++) // correlation index   
 } // end of if(fCalculateCumulantsVsM) 
 // average all correlations for all events (with correct errors!):
 TString intFlowCorrelationsAllHistName = "fIntFlowCorrelationsAllHist";
 intFlowCorrelationsAllHistName += fAnalysisLabel->Data();
 fIntFlowCorrelationsAllHist = new TH1D(intFlowCorrelationsAllHistName.Data(),"Average correlations for all events",34,0,34);
 fIntFlowCorrelationsAllHist->SetTickLength(-0.01,"Y");
 fIntFlowCorrelationsAllHist->SetMarkerStyle(25);
 fIntFlowCorrelationsAllHist->SetLabelSize(0.03);
 fIntFlowCorrelationsAllHist->SetLabelOffset(0.01,"Y");
 // 2-p correlations:
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(1,"<<2>>_{n|n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(2,"<<2>>_{2n|2n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(3,"<<2>>_{3n|3n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(4,"<<2>>_{4n|4n}");
 // 3-p correlations:
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(6,"<<3>>_{2n|n,n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(7,"<<3>>_{3n|2n,n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(8,"<<3>>_{4n|2n,2n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(9,"<<3>>_{4n|3n,n}");
 // 4-p correlations:
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(11,"<<4>>_{n,n|n,n}"); 
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(12,"<<4>>_{2n,n|2n,n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(13,"<<4>>_{2n,2n|2n,2n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(14,"<<4>>_{3n|n,n,n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(15,"<<4>>_{3n,n|3n,n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(16,"<<4>>_{3n,n|2n,2n}"); 
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(17,"<<4>>_{4n|2n,n,n}");
 // 5-p correlations:
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(19,"<<5>>_{2n|n,n,n,n}"); 
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(20,"<<5>>_{2n,2n|2n,n,n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(21,"<<5>>_{3n,n|2n,n,n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(22,"<<5>>_{4n|n,n,n,n}");
 // 6-p correlations:
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(24,"<<6>>_{n,n,n|n,n,n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(25,"<<6>>_{2n,n,n|2n,n,n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(26,"<<6>>_{2n,2n|n,n,n,n}");
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(27,"<<6>>_{3n,n|n,n,n,n}");
 // 7-p correlations:  
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(29,"<<7>>_{2n,n,n|n,n,n,n}");
 // 8-p correlations:
 (fIntFlowCorrelationsAllHist->GetXaxis())->SetBinLabel(31,"<<8>>_{n,n,n,n|n,n,n,n}");
 fIntFlowResults->Add(fIntFlowCorrelationsAllHist);
 // average correction terms for non-uniform acceptance (with correct errors!):
 for(Int_t sc=0;sc<2;sc++) // sin or cos terms
 {
  TString intFlowCorrectionTermsForNUAHistName = "fIntFlowCorrectionTermsForNUAHist";
  intFlowCorrectionTermsForNUAHistName += fAnalysisLabel->Data();
  fIntFlowCorrectionTermsForNUAHist[sc] = new TH1D(Form("%s: %s terms",intFlowCorrectionTermsForNUAHistName.Data(),sinCosFlag[sc].Data()),Form("Correction terms for non-uniform acceptance (%s terms)",sinCosFlag[sc].Data()),4,0,4);
  fIntFlowCorrectionTermsForNUAHist[sc]->SetTickLength(-0.01,"Y");
  fIntFlowCorrectionTermsForNUAHist[sc]->SetMarkerStyle(25);
  fIntFlowCorrectionTermsForNUAHist[sc]->SetLabelSize(0.05);
  fIntFlowCorrectionTermsForNUAHist[sc]->SetLabelOffset(0.01,"Y");
  (fIntFlowCorrectionTermsForNUAHist[sc]->GetXaxis())->SetBinLabel(1,Form("#LT#LT%s(n(#phi_{1}))#GT#GT",sinCosFlag[sc].Data()));
  (fIntFlowCorrectionTermsForNUAHist[sc]->GetXaxis())->SetBinLabel(2,Form("#LT#LT%s(n(#phi_{1}+#phi_{2}))#GT#GT",sinCosFlag[sc].Data()));  
  (fIntFlowCorrectionTermsForNUAHist[sc]->GetXaxis())->SetBinLabel(3,Form("#LT#LT%s(n(#phi_{1}-#phi_{2}-#phi_{3}))#GT#GT",sinCosFlag[sc].Data()));  
  (fIntFlowCorrectionTermsForNUAHist[sc]->GetXaxis())->SetBinLabel(4,Form("#LT#LT%s(n(2#phi_{1}-#phi_{2}))#GT#GT",sinCosFlag[sc].Data()));   
  fIntFlowResults->Add(fIntFlowCorrectionTermsForNUAHist[sc]);
 } // end of for(Int_t sc=0;sc<2;sc++) 
 // covariances (multiplied with weight dependent prefactor):
 TString intFlowCovariancesName = "fIntFlowCovariances";
 intFlowCovariancesName += fAnalysisLabel->Data();
 fIntFlowCovariances = new TH1D(intFlowCovariancesName.Data(),"Covariances (multiplied with weight dependent prefactor)",6,0,6);
 fIntFlowCovariances->SetLabelSize(0.04);
 fIntFlowCovariances->SetMarkerStyle(25);
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(1,"Cov(#LT2#GT,#LT4#GT)");
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(2,"Cov(#LT2#GT,#LT6#GT)");
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(3,"Cov(#LT2#GT,#LT8#GT)");
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(4,"Cov(#LT4#GT,#LT6#GT)");
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(5,"Cov(#LT4#GT,#LT8#GT)");
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(6,"Cov(#LT6#GT,#LT8#GT)");  
 fIntFlowResults->Add(fIntFlowCovariances);
 // sum of linear and quadratic event weights for <2>, <4>, <6> and <8>:
 TString intFlowSumOfEventWeightsName = "fIntFlowSumOfEventWeights";
 intFlowSumOfEventWeightsName += fAnalysisLabel->Data();
 for(Int_t power=0;power<2;power++)
 {
  fIntFlowSumOfEventWeights[power] = new TH1D(Form("%s: %s",intFlowSumOfEventWeightsName.Data(),powerFlag[power].Data()),Form("Sum of %s event weights for correlations",powerFlag[power].Data()),4,0,4);
  fIntFlowSumOfEventWeights[power]->SetLabelSize(0.04);
  fIntFlowSumOfEventWeights[power]->SetMarkerStyle(25);
  if(power == 0)
  {
   (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(1,"#sum_{i=1}^{N} w_{#LT2#GT}");
   (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(2,"#sum_{i=1}^{N} w_{#LT4#GT}");
   (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(3,"#sum_{i=1}^{N} w_{#LT6#GT}");
   (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(4,"#sum_{i=1}^{N} w_{#LT8#GT}");
  } else if (power == 1) 
    {
     (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(1,"#sum_{i=1}^{N} w_{#LT2#GT}^{2}");
     (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(2,"#sum_{i=1}^{N} w_{#LT4#GT}^{2}");
     (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(3,"#sum_{i=1}^{N} w_{#LT6#GT}^{2}");
     (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(4,"#sum_{i=1}^{N} w_{#LT8#GT}^{2}");
    }
  fIntFlowResults->Add(fIntFlowSumOfEventWeights[power]);
 } 
 // sum of products of event weights for correlations <2>, <4>, <6> and <8>:  
 TString intFlowSumOfProductOfEventWeightsName = "fIntFlowSumOfProductOfEventWeights";
 intFlowSumOfProductOfEventWeightsName += fAnalysisLabel->Data();
 fIntFlowSumOfProductOfEventWeights = new TH1D(intFlowSumOfProductOfEventWeightsName.Data(),"Sum of product of event weights for correlations",6,0,6);
 fIntFlowSumOfProductOfEventWeights->SetLabelSize(0.04);
 fIntFlowSumOfProductOfEventWeights->SetMarkerStyle(25);
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(1,"#sum_{i=1}^{N} w_{#LT2#GT} w_{#LT4#GT}");
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(2,"#sum_{i=1}^{N} w_{#LT2#GT} w_{#LT6#GT}");
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(3,"#sum_{i=1}^{N} w_{#LT2#GT} w_{#LT8#GT}");
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(4,"#sum_{i=1}^{N} w_{#LT4#GT} w_{#LT6#GT}");
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(5,"#sum_{i=1}^{N} w_{#LT4#GT} w_{#LT8#GT}");
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(6,"#sum_{i=1}^{N} w_{#LT6#GT} w_{#LT8#GT}");
 fIntFlowResults->Add(fIntFlowSumOfProductOfEventWeights);
 // final result for covariances of correlations (multiplied with weight dependent prefactor) versus M
 // [0=Cov(2,4),1=Cov(2,6),2=Cov(2,8),3=Cov(4,6),4=Cov(4,8),5=Cov(6,8)]:
 if(fCalculateCumulantsVsM)
 {
  TString intFlowCovariancesVsMName = "fIntFlowCovariancesVsM";
  intFlowCovariancesVsMName += fAnalysisLabel->Data();
  TString covarianceFlag[6] = {"Cov(<2>,<4>)","Cov(<2>,<6>)","Cov(<2>,<8>)","Cov(<4>,<6>)","Cov(<4>,<8>)","Cov(<6>,<8>)"};
  for(Int_t ci=0;ci<6;ci++)
  {
   fIntFlowCovariancesVsM[ci] = new TH1D(Form("%s, %s",intFlowCovariancesVsMName.Data(),covarianceFlag[ci].Data()),
                                         Form("%s vs multiplicity",covarianceFlag[ci].Data()),
                                         fnBinsMult,fMinMult,fMaxMult);
   fIntFlowCovariancesVsM[ci]->GetYaxis()->SetTitle(covarianceFlag[ci].Data());
   fIntFlowCovariancesVsM[ci]->GetXaxis()->SetTitle("M");
   fIntFlowResults->Add(fIntFlowCovariancesVsM[ci]);
  }
 } // end of if(fCalculateCumulantsVsM) 
 // sum of linear and quadratic event weights for <2>, <4>, <6> and <8> versus multiplicity
 // [0=sum{w_{<2>}},1=sum{w_{<4>}},2=sum{w_{<6>}},3=sum{w_{<8>}}][0=linear 1,1=quadratic]:
 if(fCalculateCumulantsVsM)
 {
  TString intFlowSumOfEventWeightsVsMName = "fIntFlowSumOfEventWeightsVsM";
  intFlowSumOfEventWeightsVsMName += fAnalysisLabel->Data();
  TString sumFlag[2][4] = {{"#sum_{i=1}^{N} w_{<2>}","#sum_{i=1}^{N} w_{<4>}","#sum_{i=1}^{N} w_{<6>}","#sum_{i=1}^{N} w_{<8>}"},
                           {"#sum_{i=1}^{N} w_{<2>}^{2}","#sum_{i=1}^{N} w_{<4>}^{2}","#sum_{i=1}^{N} w_{<6>}^{2}","#sum_{i=1}^{N} w_{<8>}^{2}"}};
  for(Int_t si=0;si<4;si++)
  {
   for(Int_t power=0;power<2;power++)
   {
    fIntFlowSumOfEventWeightsVsM[si][power] = new TH1D(Form("%s, %s",intFlowSumOfEventWeightsVsMName.Data(),sumFlag[power][si].Data()),
                                                       Form("%s vs multiplicity",sumFlag[power][si].Data()),
                                                       fnBinsMult,fMinMult,fMaxMult);    
    fIntFlowSumOfEventWeightsVsM[si][power]->GetYaxis()->SetTitle(sumFlag[power][si].Data());  
    fIntFlowSumOfEventWeightsVsM[si][power]->GetXaxis()->SetTitle("M");  
    fIntFlowResults->Add(fIntFlowSumOfEventWeightsVsM[si][power]);
   } // end of for(Int_t power=0;power<2;power++)
  } // end of for(Int_t si=0;si<4;si++)   
 } // end of if(fCalculateCumulantsVsM)
 // sum of products of event weights for correlations <2>, <4>, <6> and <8> vs M
 // [0=sum{w_{<2>}w_{<4>}},1=sum{w_{<2>}w_{<6>}},2=sum{w_{<2>}w_{<8>}},
 //  3=sum{w_{<4>}w_{<6>}},4=sum{w_{<4>}w_{<8>}},5=sum{w_{<6>}w_{<8>}}]:  
 if(fCalculateCumulantsVsM)
 {
  TString intFlowSumOfProductOfEventWeightsVsMName = "fIntFlowSumOfProductOfEventWeightsVsM";
  intFlowSumOfProductOfEventWeightsVsMName += fAnalysisLabel->Data();
  TString sopowFlag[6] = {"#sum_{i=1}^{N} w_{<2>} w_{<4>}","#sum_{i=1}^{N} w_{<2>} w_{<6>}","#sum_{i=1}^{N} w_{<2>} w_{<8>}",
                          "#sum_{i=1}^{N} w_{<4>} w_{<6>}","#sum_{i=1}^{N} w_{<4>} w_{<8>}","#sum_{i=1}^{N} w_{<6>} w_{<8>}"}; 
  for(Int_t pi=0;pi<6;pi++)
  {
   fIntFlowSumOfProductOfEventWeightsVsM[pi] = new TH1D(Form("%s, %s",intFlowSumOfProductOfEventWeightsVsMName.Data(),sopowFlag[pi].Data()),
                                                        Form("%s versus multiplicity",sopowFlag[pi].Data()),
                                                        fnBinsMult,fMinMult,fMaxMult); 
   fIntFlowSumOfProductOfEventWeightsVsM[pi]->GetXaxis()->SetTitle("M");
   fIntFlowSumOfProductOfEventWeightsVsM[pi]->GetYaxis()->SetTitle(sopowFlag[pi].Data()); 
   fIntFlowResults->Add(fIntFlowSumOfProductOfEventWeightsVsM[pi]);
  } // end of for(Int_t pi=0;pi<6;pi++) 
 } // end of if(fCalculateCumulantsVsM)
 // covariances of NUA terms (multiplied with weight dependent prefactor):
 TString intFlowCovariancesNUAName = "fIntFlowCovariancesNUA";
 intFlowCovariancesNUAName += fAnalysisLabel->Data();
 fIntFlowCovariancesNUA = new TH1D(intFlowCovariancesNUAName.Data(),"Covariances for NUA (multiplied with weight dependent prefactor)",27,0,27);
 fIntFlowCovariancesNUA->SetLabelSize(0.04);
 fIntFlowCovariancesNUA->SetMarkerStyle(25);
 fIntFlowCovariancesNUA->GetXaxis()->SetLabelSize(0.02);
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(1,"Cov(<2>,<cos(#phi)>");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(2,"Cov(<2>,<sin(#phi)>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(3,"Cov(<cos(#phi)>,<sin(#phi)>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(4,"Cov(<2>,<cos(#phi_{1}+#phi_{2})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(5,"Cov(<2>,<sin(#phi_{1}+#phi_{2})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(6,"Cov(<2>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(7,"Cov(<2>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(8,"Cov(<4>,<cos(#phi)>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(9,"Cov(<4>,<sin(#phi)>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(10,"Cov(<4>,<cos(#phi_{1}+#phi_{2})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(11,"Cov(<4>,<sin(#phi_{1}+#phi_{2})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(12,"Cov(<4>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(13,"Cov(<4>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(14,"Cov(<cos(#phi)>,<cos(#phi_{1}+#phi_{2})>)"); 
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(15,"Cov(<cos(#phi)>,<sin(#phi_{1}+#phi_{2})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(16,"Cov(<cos(#phi)>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(17,"Cov(<cos(#phi)>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(18,"Cov(<sin(#phi)>,<cos(#phi_{1}+#phi_{2})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(19,"Cov(<sin(#phi)>,<sin(#phi_{1}+#phi_{2})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(20,"Cov(<sin(#phi)>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(21,"Cov(<sin(#phi)>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(22,"Cov(<cos(#phi_{1}+#phi_{2})>,<sin(#phi_{1}+#phi_{2})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(23,"Cov(<cos(#phi_{1}+#phi_{2})>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(24,"Cov(<cos(#phi_{1}+#phi_{2})>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(25,"Cov(<sin(#phi_{1}+#phi_{2})>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(26,"Cov(<sin(#phi_{1}+#phi_{2})>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)");
 (fIntFlowCovariancesNUA->GetXaxis())->SetBinLabel(27,"Cov(<cos(#phi_{1}-#phi_{2}-#phi_{3}>,<sin(#phi_{1}-#phi_{2}-#phi_{3}>)");
 fIntFlowResults->Add(fIntFlowCovariancesNUA);
 // sum of linear and quadratic event weights for NUA terms:
 TString intFlowSumOfEventWeightsNUAName = "fIntFlowSumOfEventWeightsNUA";
 intFlowSumOfEventWeightsNUAName += fAnalysisLabel->Data();
 for(Int_t sc=0;sc<2;sc++)
 {
  for(Int_t power=0;power<2;power++)
  {
   fIntFlowSumOfEventWeightsNUA[sc][power] = new TH1D(Form("%s: %s, %s",intFlowSumOfEventWeightsNUAName.Data(),powerFlag[power].Data(),sinCosFlag[sc].Data()),Form("Sum of %s event weights for NUA %s terms",powerFlag[power].Data(),sinCosFlag[sc].Data()),4,0,4); // to be improved - 4
   fIntFlowSumOfEventWeightsNUA[sc][power]->SetLabelSize(0.05);
   fIntFlowSumOfEventWeightsNUA[sc][power]->SetMarkerStyle(25);
   if(power == 0)
   {
    (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(1,Form("#sum_{i=1}^{N} w_{<%s(#phi)>}",sinCosFlag[sc].Data()));
    (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(2,Form("#sum_{i=1}^{N} w_{<%s(#phi_{1}+#phi_{2})>}",sinCosFlag[sc].Data()));
    (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(3,Form("#sum_{i=1}^{N} w_{<%s(#phi_{1}-#phi_{2}-#phi_{3})>}",sinCosFlag[sc].Data()));   
    (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(4,Form("#sum_{i=1}^{N} w_{<%s(2#phi_{1}-#phi_{2})>}",sinCosFlag[sc].Data()));
   } else if(power == 1) 
     {
      (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(1,Form("#sum_{i=1}^{N} w_{<%s(#phi)>}^{2}",sinCosFlag[sc].Data()));
      (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(2,Form("#sum_{i=1}^{N} w_{<%s(#phi_{1}+#phi_{2})>}^{2}",sinCosFlag[sc].Data()));
      (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(3,Form("#sum_{i=1}^{N} w_{<%s(#phi_{1}-#phi_{2}-#phi_{3})>}^{2}",sinCosFlag[sc].Data()));
      (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(4,Form("#sum_{i=1}^{N} w_{<%s(2#phi_{1}-#phi_{2})>}^{2}",sinCosFlag[sc].Data()));
     }
   fIntFlowResults->Add(fIntFlowSumOfEventWeightsNUA[sc][power]);
  }
 }  
 // sum of products of event weights for NUA terms:  
 TString intFlowSumOfProductOfEventWeightsNUAName = "fIntFlowSumOfProductOfEventWeightsNUA";
 intFlowSumOfProductOfEventWeightsNUAName += fAnalysisLabel->Data();
 fIntFlowSumOfProductOfEventWeightsNUA = new TH1D(intFlowSumOfProductOfEventWeightsNUAName.Data(),"Sum of product of event weights for NUA terms",27,0,27);
 fIntFlowSumOfProductOfEventWeightsNUA->SetLabelSize(0.02);
 fIntFlowSumOfProductOfEventWeightsNUA->SetMarkerStyle(25);
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(1,"#sum_{i=1}^{N} w_{#LT2#GT} w_{#LTcos(#phi)#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(2,"#sum_{i=1}^{N} w_{#LT2#GT} w_{#LTsin(#phi)#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(3,"#sum_{i=1}^{N} w_{#LTcos(#phi)#GT} w_{#LTsin(#phi)#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(4,"#sum_{i=1}^{N} w_{#LT2#GT} w_{#LTcos(#phi_{1}+#phi_{2})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(5,"#sum_{i=1}^{N} w_{#LT2#GT} w_{#LTsin(#phi_{1}+#phi_{2})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(6,"#sum_{i=1}^{N} w_{#LT2#GT} w_{#LTcos(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(7,"#sum_{i=1}^{N} w_{#LT2#GT} w_{#LTsin(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(8,"#sum_{i=1}^{N} w_{#LT4#GT} w_{#LTcos(#phi)#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(9,"#sum_{i=1}^{N} w_{#LT4#GT} w_{#LTsin(#phi)#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(10,"#sum_{i=1}^{N} w_{#LT4#GT} w_{#LTcos(#phi_{1}+#phi_{2})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(11,"#sum_{i=1}^{N} w_{#LT4#GT} w_{#LTsin(#phi_{1}+#phi_{2})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(12,"#sum_{i=1}^{N} w_{#LT4#GT} w_{#LTcos(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(13,"#sum_{i=1}^{N} w_{#LT4#GT} w_{#LTsin(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(14,"#sum_{i=1}^{N} w_{#LTcos(#phi)#GT} w_{#LTcos(#phi_{1}+#phi_{2})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(15,"#sum_{i=1}^{N} w_{#LTcos(#phi)#GT} w_{#LTsin(#phi_{1}+#phi_{2})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(16,"#sum_{i=1}^{N} w_{#LTcos(#phi)#GT} w_{#LTcos(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(17,"#sum_{i=1}^{N} w_{#LTcos(#phi)#GT} w_{#LTsin(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(18,"#sum_{i=1}^{N} w_{#LTsin(#phi)#GT} w_{#LTcos(#phi_{1}+#phi_{2})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(19,"#sum_{i=1}^{N} w_{#LTsin(#phi)#GT} w_{#LTsin(#phi_{1}+#phi_{2})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(20,"#sum_{i=1}^{N} w_{#LTsin(#phi)#GT} w_{#LTcos(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(21,"#sum_{i=1}^{N} w_{#LTsin(#phi)#GT} w_{#LTsin(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(22,"#sum_{i=1}^{N} w_{#LTcos(#phi_{1}+#phi_{2})#GT} w_{#LTsin(#phi_{1}+#phi_{2})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(23,"#sum_{i=1}^{N} w_{#LTcos(#phi_{1}+#phi_{2})#GT} w_{#LTcos(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(24,"#sum_{i=1}^{N} w_{#LTcos(#phi_{1}+#phi_{2})#GT} w_{#LTsin(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(25,"#sum_{i=1}^{N} w_{#LTsin(#phi_{1}+#phi_{2})#GT} w_{#LTcos(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(26,"#sum_{i=1}^{N} w_{#LTsin(#phi_{1}+#phi_{2})#GT} w_{#LTsin(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(27,"#sum_{i=1}^{N} w_{#LTcos(#phi_{1}-#phi_{2}-#phi_{3})#GT} w_{#LTsin(#phi_{1}-#phi_{2}-#phi_{3})#GT}");
 fIntFlowResults->Add(fIntFlowSumOfProductOfEventWeightsNUA);
 // Final results for reference Q-cumulants:
 TString cumulantFlag[4] = {"QC{2}","QC{4}","QC{6}","QC{8}"};
 TString intFlowQcumulantsName = "fIntFlowQcumulants";
 intFlowQcumulantsName += fAnalysisLabel->Data();
 fIntFlowQcumulants = new TH1D(intFlowQcumulantsName.Data(),"Reference Q-cumulants",4,0,4);
 if(fPropagateErrorAlsoFromNIT)
 {
  fIntFlowQcumulants->SetTitle("Reference Q-cumulants (error from non-isotropic terms also propagated)");
 }
 fIntFlowQcumulants->SetLabelSize(0.05);
 fIntFlowQcumulants->SetMarkerStyle(25);
 for(Int_t b=0;b<4;b++)
 {
  (fIntFlowQcumulants->GetXaxis())->SetBinLabel(b+1,cumulantFlag[b].Data());
 } 
 fIntFlowResults->Add(fIntFlowQcumulants);
 // Final results for reference Q-cumulants rebinned in M: 
 if(fCalculateCumulantsVsM)
 {
  TString intFlowQcumulantsRebinnedInMName = "fIntFlowQcumulantsRebinnedInM";
  intFlowQcumulantsRebinnedInMName += fAnalysisLabel->Data();
  fIntFlowQcumulantsRebinnedInM = new TH1D(intFlowQcumulantsRebinnedInMName.Data(),"Reference Q-cumulants rebinned in M",4,0,4);
  fIntFlowQcumulantsRebinnedInM->SetLabelSize(0.05);
  fIntFlowQcumulantsRebinnedInM->SetMarkerStyle(25);
  for(Int_t b=0;b<4;b++)
  {
   (fIntFlowQcumulantsRebinnedInM->GetXaxis())->SetBinLabel(b+1,cumulantFlag[b].Data());
  } 
  fIntFlowResults->Add(fIntFlowQcumulantsRebinnedInM);
 } // end of if(fCalculateCumulantsVsM) 
 // Ratio between error squared: with/without non-isotropic terms:
 TString intFlowQcumulantsErrorSquaredRatioName = "fIntFlowQcumulantsErrorSquaredRatio";
 intFlowQcumulantsErrorSquaredRatioName += fAnalysisLabel->Data();
 fIntFlowQcumulantsErrorSquaredRatio = new TH1D(intFlowQcumulantsErrorSquaredRatioName.Data(),"Error squared of reference Q-cumulants: #frac{with NUA terms}{without NUA terms}",4,0,4);
 fIntFlowQcumulantsErrorSquaredRatio->SetLabelSize(0.05);
 fIntFlowQcumulantsErrorSquaredRatio->SetMarkerStyle(25);
 for(Int_t b=0;b<4;b++)
 {
  (fIntFlowQcumulantsErrorSquaredRatio->GetXaxis())->SetBinLabel(b+1,cumulantFlag[b].Data());
 } 
 fIntFlowResults->Add(fIntFlowQcumulantsErrorSquaredRatio);
 // final results for integrated Q-cumulants versus multiplicity:
 if(fCalculateCumulantsVsM)
 {
  TString intFlowQcumulantsVsMName = "fIntFlowQcumulantsVsM";
  intFlowQcumulantsVsMName += fAnalysisLabel->Data();
  for(Int_t co=0;co<4;co++) // cumulant order
  {
   fIntFlowQcumulantsVsM[co] = new TH1D(Form("%s, %s",intFlowQcumulantsVsMName.Data(),cumulantFlag[co].Data()),
                                        Form("%s vs multipicity",cumulantFlag[co].Data()),
                                        fnBinsMult,fMinMult,fMaxMult);
   fIntFlowQcumulantsVsM[co]->GetXaxis()->SetTitle("M");                                     
   fIntFlowQcumulantsVsM[co]->GetYaxis()->SetTitle(cumulantFlag[co].Data());  
   fIntFlowResults->Add(fIntFlowQcumulantsVsM[co]);                                    
  } // end of for(Int_t co=0;co<4;co++) // cumulant order
 } // end of if(fCalculateCumulantsVsM)
 // final integrated flow estimates from Q-cumulants:
 TString flowFlag[4] = {Form("v_{%d}{2,QC}",fHarmonic),Form("v_{%d}{4,QC}",fHarmonic),Form("v_{%d}{6,QC}",fHarmonic),Form("v_{%d}{8,QC}",fHarmonic)};
 TString intFlowName = "fIntFlow";
 intFlowName += fAnalysisLabel->Data();  
 // integrated flow from Q-cumulants:
 fIntFlow = new TH1D(intFlowName.Data(),"Reference flow estimates from Q-cumulants",4,0,4);
 fIntFlow->SetLabelSize(0.05);
 fIntFlow->SetMarkerStyle(25);
 for(Int_t b=0;b<4;b++)
 {
  (fIntFlow->GetXaxis())->SetBinLabel(b+1,flowFlag[b].Data()); 
 }
 fIntFlowResults->Add(fIntFlow); 
 // Reference flow vs M rebinned in one huge bin:
 if(fCalculateCumulantsVsM)
 { 
  TString intFlowRebinnedInMName = "fIntFlowRebinnedInM";
  intFlowRebinnedInMName += fAnalysisLabel->Data();  
  fIntFlowRebinnedInM = new TH1D(intFlowRebinnedInMName.Data(),"Reference flow estimates from Q-cumulants (rebinned in M)",4,0,4);
  fIntFlowRebinnedInM->SetLabelSize(0.05);
  fIntFlowRebinnedInM->SetMarkerStyle(25);
  for(Int_t b=0;b<4;b++)
  {
   (fIntFlowRebinnedInM->GetXaxis())->SetBinLabel(b+1,flowFlag[b].Data()); 
  }
  fIntFlowResults->Add(fIntFlowRebinnedInM); 
 } 
 // integrated flow from Q-cumulants: versus multiplicity:
 if(fCalculateCumulantsVsM)
 {
  TString intFlowVsMName = "fIntFlowVsM";
  intFlowVsMName += fAnalysisLabel->Data();
  for(Int_t co=0;co<4;co++) // cumulant order
  {
   fIntFlowVsM[co] = new TH1D(Form("%s, %s",intFlowVsMName.Data(),flowFlag[co].Data()),
                              Form("%s vs multipicity",flowFlag[co].Data()),
                              fnBinsMult,fMinMult,fMaxMult);
   fIntFlowVsM[co]->GetXaxis()->SetTitle("M");                                     
   fIntFlowVsM[co]->GetYaxis()->SetTitle(flowFlag[co].Data());  
   fIntFlowResults->Add(fIntFlowVsM[co]);                                    
  } // end of for(Int_t co=0;co<4;co++) // cumulant order
 } // end of if(fCalculateCumulantsVsM)
 // quantifying detector effects effects to correlations:
 TString intFlowDetectorBiasName = "fIntFlowDetectorBias";
 intFlowDetectorBiasName += fAnalysisLabel->Data();  
 fIntFlowDetectorBias = new TH1D(intFlowDetectorBiasName.Data(),"Quantifying detector bias",4,0,4);
 fIntFlowDetectorBias->SetLabelSize(0.05);
 fIntFlowDetectorBias->SetMarkerStyle(25);
 for(Int_t ci=0;ci<4;ci++)
 {  
  (fIntFlowDetectorBias->GetXaxis())->SetBinLabel(ci+1,Form("#frac{corrected}{measured} %s",cumulantFlag[ci].Data()));
 }
 fIntFlowResults->Add(fIntFlowDetectorBias); 
 // quantifying detector effects to correlations versus multiplicity:
 if(fCalculateCumulantsVsM)
 {
  TString intFlowDetectorBiasVsMName = "fIntFlowDetectorBiasVsM";
  intFlowDetectorBiasVsMName += fAnalysisLabel->Data();
  for(Int_t ci=0;ci<4;ci++) // correlation index
  {
   fIntFlowDetectorBiasVsM[ci] = new TH1D(Form("%s for %s",intFlowDetectorBiasVsMName.Data(),cumulantFlag[ci].Data()),
                                          Form("Quantifying detector bias for %s vs multipicity",cumulantFlag[ci].Data()),
                                          fnBinsMult,fMinMult,fMaxMult);
   fIntFlowDetectorBiasVsM[ci]->GetXaxis()->SetTitle("M");                                     
   fIntFlowDetectorBiasVsM[ci]->GetYaxis()->SetTitle("#frac{corrected}{measured}");  
   fIntFlowResults->Add(fIntFlowDetectorBiasVsM[ci]);                                    
  } // end of for(Int_t co=0;co<4;co++) // cumulant order
 } // end of if(fCalculateCumulantsVsM)
   
} // end of AliFlowAnalysisWithQCumulants::BookEverythingForIntegratedFlow()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::InitializeArraysForNestedLoops()
{
 // Initialize arrays of all objects relevant for calculations with nested loops.
 
 // integrated flow:
 for(Int_t sc=0;sc<2;sc++) // sin or cos terms
 {
  fIntFlowDirectCorrectionTermsForNUA[sc] = NULL;
 } 

 // differential flow:  
 // correlations:
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  {
   for(Int_t ci=0;ci<4;ci++) // correlation index
   {
    fDiffFlowDirectCorrelations[t][pe][ci] = NULL;
   } // end of for(Int_t ci=0;ci<4;ci++) // correlation index  
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
 // correction terms for non-uniform acceptance:
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  {
   for(Int_t sc=0;sc<2;sc++) // sin or cos terms
   {
    for(Int_t cti=0;cti<9;cti++) // correction term index
    {
     fDiffFlowDirectCorrectionTermsForNUA[t][pe][sc][cti] = NULL;
    }   
   }
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI

 // other differential correlators: 
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  {
   for(Int_t sc=0;sc<2;sc++) // sin or cos terms
   {
    for(Int_t ci=0;ci<1;ci++) // correlator index
    {
     fOtherDirectDiffCorrelators[t][pe][sc][ci] = NULL;
    }   
   }
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI

} // end of void AliFlowAnalysisWithQCumulants::InitializeArraysForNestedLoops()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::BookEverythingForNestedLoops()
{
 // Book all objects relevant for calculations with nested loops.
 
 TString sinCosFlag[2] = {"sin","cos"}; // to be improved (should I promote this to data members?)
 TString typeFlag[2] = {"RP","POI"}; // to be improved (should I promote this to data members?)
 TString ptEtaFlag[2] = {"p_{T}","#eta"}; // to be improved (should I promote this to data members?)
 TString reducedCorrelationIndex[4] = {"<2'>","<4'>","<6'>","<8'>"}; // to be improved (should I promote this to data members?)
 Double_t lowerPtEtaEdge[2] = {fPtMin+(fCrossCheckInPtBinNo-1)*fPtBinWidth,fEtaMin+(fCrossCheckInEtaBinNo-1)*fEtaBinWidth};
 Double_t upperPtEtaEdge[2] = {fPtMin+fCrossCheckInPtBinNo*fPtBinWidth,fEtaMin+fCrossCheckInEtaBinNo*fEtaBinWidth};

 TString evaluateNestedLoopsName = "fEvaluateNestedLoops";
 evaluateNestedLoopsName += fAnalysisLabel->Data();
 fEvaluateNestedLoops = new TProfile(evaluateNestedLoopsName.Data(),"Flags for nested loops",4,0,4);
 fEvaluateNestedLoops->SetLabelSize(0.03);
 (fEvaluateNestedLoops->GetXaxis())->SetBinLabel(1,"fEvaluateIntFlowNestedLoops");
 (fEvaluateNestedLoops->GetXaxis())->SetBinLabel(2,"fEvaluateDiffFlowNestedLoops");
 (fEvaluateNestedLoops->GetXaxis())->SetBinLabel(3,"fCrossCheckInPtBinNo");
 (fEvaluateNestedLoops->GetXaxis())->SetBinLabel(4,"fCrossCheckInEtaBinNo");
 fEvaluateNestedLoops->Fill(0.5,(Int_t)fEvaluateIntFlowNestedLoops);
 fEvaluateNestedLoops->Fill(1.5,(Int_t)fEvaluateDiffFlowNestedLoops);
 fEvaluateNestedLoops->Fill(2.5,fCrossCheckInPtBinNo);
 fEvaluateNestedLoops->Fill(3.5,fCrossCheckInEtaBinNo);
 fNestedLoopsList->Add(fEvaluateNestedLoops);
 // nested loops for integrated flow:
 if(fEvaluateIntFlowNestedLoops)
 {
  // correlations:
  TString intFlowDirectCorrelationsName = "fIntFlowDirectCorrelations";
  intFlowDirectCorrelationsName += fAnalysisLabel->Data();
  fIntFlowDirectCorrelations = new TProfile(intFlowDirectCorrelationsName.Data(),"Multiparticle correlations calculated with nested loops (for int. flow)",64,0,64,"s");
  fNestedLoopsList->Add(fIntFlowDirectCorrelations);
  if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)
  {
   TString intFlowExtraDirectCorrelationsName = "fIntFlowExtraDirectCorrelations";
   intFlowExtraDirectCorrelationsName += fAnalysisLabel->Data();
   fIntFlowExtraDirectCorrelations = new TProfile(intFlowExtraDirectCorrelationsName.Data(),"Extra multiparticle correlations calculated with nested loops (for int. flow)",100,0,100,"s");
   fNestedLoopsList->Add(fIntFlowExtraDirectCorrelations);  
  } // end of if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)
  // correction terms for non-uniform acceptance:
  for(Int_t sc=0;sc<2;sc++) // sin or cos terms
  {
   TString intFlowDirectCorrectionTermsForNUAName = "fIntFlowDirectCorrectionTermsForNUA";
   intFlowDirectCorrectionTermsForNUAName += fAnalysisLabel->Data();
   fIntFlowDirectCorrectionTermsForNUA[sc] = new TProfile(Form("%s: %s terms",intFlowDirectCorrectionTermsForNUAName.Data(),sinCosFlag[sc].Data()),Form("Correction terms for non-uniform acceptance (%s terms)",sinCosFlag[sc].Data()),10,0,10,"s");
   fNestedLoopsList->Add(fIntFlowDirectCorrectionTermsForNUA[sc]);
  } // end of for(Int_t sc=0;sc<2;sc++) 
 } // end of if(fEvaluateIntFlowNestedLoops)
 
 // nested loops for differential flow: 
 if(fEvaluateDiffFlowNestedLoops)
 {
  // reduced correlations:
  TString diffFlowDirectCorrelationsName = "fDiffFlowDirectCorrelations";
  diffFlowDirectCorrelationsName += fAnalysisLabel->Data();
  for(Int_t t=0;t<2;t++) // type: RP or POI
  { 
   for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
   {
    for(Int_t rci=0;rci<4;rci++) // reduced correlation index
    {
     // reduced correlations:
     fDiffFlowDirectCorrelations[t][pe][rci] = new TProfile(Form("%s, %s, %s, %s",diffFlowDirectCorrelationsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[rci].Data()),Form("%s, %s, %s, %s",diffFlowDirectCorrelationsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[rci].Data()),1,lowerPtEtaEdge[pe],upperPtEtaEdge[pe],"s");
     fDiffFlowDirectCorrelations[t][pe][rci]->SetXTitle(ptEtaFlag[pe].Data());
     fNestedLoopsList->Add(fDiffFlowDirectCorrelations[t][pe][rci]); // to be improved (add dedicated list to hold reduced correlations)
    } // end of for(Int_t rci=0;rci<4;rci++) // correlation index
   } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta 
  } // end of for(Int_t t=0;t<2;t++) // type: RP or POI 
  
  
  // correction terms for non-uniform acceptance:
  TString diffFlowDirectCorrectionTermsForNUAName = "fDiffFlowDirectCorrectionTermsForNUA";
  diffFlowDirectCorrectionTermsForNUAName += fAnalysisLabel->Data();
  for(Int_t t=0;t<2;t++) // typeFlag (0 = RP, 1 = POI)
  { 
   for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
   {
    for(Int_t sc=0;sc<2;sc++) // sin or cos
    {
     for(Int_t cti=0;cti<9;cti++) // correction term index
     {
      fDiffFlowDirectCorrectionTermsForNUA[t][pe][sc][cti] = new TProfile(Form("%s, %s, %s, %s, cti = %d",diffFlowDirectCorrectionTermsForNUAName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),cti+1),Form("%s, %s, %s, %s, cti = %d",diffFlowDirectCorrectionTermsForNUAName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),cti+1),1,lowerPtEtaEdge[pe],upperPtEtaEdge[pe],"s"); 
      fNestedLoopsList->Add(fDiffFlowDirectCorrectionTermsForNUA[t][pe][sc][cti]);
     }
    }
   }
  }
  // other differential correlators: 
  TString otherDirectDiffCorrelatorsName = "fOtherDirectDiffCorrelators";
  otherDirectDiffCorrelatorsName += fAnalysisLabel->Data();
  for(Int_t t=0;t<2;t++) // typeFlag (0 = RP, 1 = POI)
  { 
   for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
   {
    for(Int_t sc=0;sc<2;sc++) // sin or cos
    {
     for(Int_t ci=0;ci<1;ci++) // correlator index
     {
      fOtherDirectDiffCorrelators[t][pe][sc][ci] = new TProfile(Form("%s, %s, %s, %s, ci = %d",otherDirectDiffCorrelatorsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),ci+1),Form("%s, %s, %s, %s, ci = %d",otherDirectDiffCorrelatorsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),ci+1),1,lowerPtEtaEdge[pe],upperPtEtaEdge[pe]); 
      fNestedLoopsList->Add(fOtherDirectDiffCorrelators[t][pe][sc][ci]);
     }
    }
   }
  }
  // number of RPs and POIs in selected pt and eta bins for cross-checkings:
  TString noOfParticlesInBinName = "fNoOfParticlesInBin";
  fNoOfParticlesInBin = new TH1D(noOfParticlesInBinName.Data(),"Number of RPs and POIs in selected p_{T} and #eta bin",4,0,4);
  fNoOfParticlesInBin->GetXaxis()->SetBinLabel(1,"# of RPs in p_{T} bin");
  fNoOfParticlesInBin->GetXaxis()->SetBinLabel(2,"# of RPs in #eta bin");
  fNoOfParticlesInBin->GetXaxis()->SetBinLabel(3,"# of POIs in p_{T} bin");
  fNoOfParticlesInBin->GetXaxis()->SetBinLabel(4,"# of POIs in #eta bin");
  fNestedLoopsList->Add(fNoOfParticlesInBin);
 } // end of if(fEvaluateDiffFlowNestedLoops)

} // end of AliFlowAnalysisWithQCumulants::BookEverythingForNestedLoops()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrelations()
{
 // Calculate in this method all multiparticle azimuthal correlations.
 //
 // Remark 1: All multiparticle correlations are stored in TProfile fIntFlowCorrelationsAllPro;
 // Remark 2: There is a special TProfile fIntFlowCorrelationsPro holding results 
 //           only for same harmonic's correlations <<2>>, <<4>>, <<6>> and <<8>>;  
 // Remark 3: Binning of fIntFlowCorrelationsAllPro is organized as follows:
 // --------------------------------------------------------------------------------------------------------------------
 //  1st bin: <2>_{1n|1n} = two1n1n = cos(1n(phi1-phi2))>
 //  2nd bin: <2>_{2n|2n} = two2n2n = cos(2n(phi1-phi2))>
 //  3rd bin: <2>_{3n|3n} = two3n3n = cos(3n(phi1-phi2))> 
 //  4th bin: <2>_{4n|4n} = two4n4n = cos(4n(phi1-phi2))>
 //  5th bin:           ----  EMPTY ----
 //  6th bin: <3>_{2n|1n,1n} = three2n1n1n = <cos(n(2*phi1-phi2-phi3))>
 //  7th bin: <3>_{3n|2n,1n} = three3n2n1n = <cos(n(3*phi1-2*phi2-phi3))>
 //  8th bin: <3>_{4n|2n,2n} = three4n2n2n = <cos(n(4*phi1-2*phi2-2*phi3))>
 //  9th bin: <3>_{4n|3n,1n} = three4n3n1n = <cos(n(4*phi1-3*phi2-phi3))>
 // 10th bin:           ----  EMPTY ----
 // 11th bin: <4>_{1n,1n|1n,1n} = four1n1n1n1n = <cos(n(phi1+phi2-phi3-phi4))>
 // 12th bin: <4>_{2n,1n|2n,1n} = four2n1n2n1n = <cos(n(2*phi1+phi2-2*phi3-phi4))>
 // 13th bin: <4>_{2n,2n|2n,2n} = four2n2n2n2n = <cos(2n(phi1+phi2-phi3-phi4))>
 // 14th bin: <4>_{3n|1n,1n,1n} = four3n1n1n1n = <cos(n(3*phi1-phi2-phi3-phi4))> 
 // 15th bin: <4>_{3n,1n|3n,1n} = four3n1n3n1n = <cos(n(3*phi1+phi2-3*phi3-phi4))>
 // 16th bin: <4>_{3n,1n|2n,2n} = four3n1n2n2n = <cos(n(3*phi1+phi2-2*phi3-2*phi4))>
 // 17th bin: <4>_{4n|2n,1n,1n} = four4n2n1n1n = <cos(n(4*phi1-2*phi2-phi3-phi4))>
 // 18th bin:           ----  EMPTY ----
 // 19th bin: <5>_{2n,1n|1n,1n,1n} = five2n1n1n1n1n = <cos(n(2*phi1+phi2-phi3-phi4-phi5))>
 // 20th bin: <5>_{2n,2n|2n,1n,1n} = five2n2n2n1n1n = <cos(n(2*phi1+2*phi2-2*phi3-phi4-phi5))>
 // 21st bin: <5>_{3n,1n|2n,1n,1n} = five3n1n2n1n1n = <cos(n(3*phi1+phi2-2*phi3-phi4-phi5))>
 // 22nd bin: <5>_{4n|1n,1n,1n,1n} = five4n1n1n1n1n = <cos(n(4*phi1-phi2-phi3-phi4-phi5))>
 // 23rd bin:           ----  EMPTY ----
 // 24th bin: <6>_{1n,1n,1n|1n,1n,1n} = six1n1n1n1n1n1n = <cos(n(phi1+phi2+phi3-phi4-phi5-phi6))>
 // 25th bin: <6>_{2n,1n,1n|2n,1n,1n} = six2n1n1n2n1n1n = <cos(n(2*phi1+phi2+phi3-2*phi4-phi5-phi6))>
 // 26th bin: <6>_{2n,2n|1n,1n,1n,1n} = six2n2n1n1n1n1n = <cos(n(2*phi1+2*phi2-phi3-phi4-phi5-phi6))>
 // 27th bin: <6>_{3n,1n|1n,1n,1n,1n} = six3n1n1n1n1n1n = <cos(n(3*phi1+phi2-phi3-phi4-phi5-phi6))> 
 // 28th bin:           ----  EMPTY ----
 // 29th bin: <7>_{2n,1n,1n|1n,1n,1n,1n} = seven2n1n1n1n1n1n1n =  <cos(n(2*phi1+phi2+phi3-phi4-phi5-phi6-phi7))>
 // 30th bin:           ----  EMPTY ----
 // 31st bin: <8>_{1n,1n,1n,1n|1n,1n,1n,1n} = eight1n1n1n1n1n1n1n1n = <cos(n(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8))>
 // 32nd bin:           ----  EMPTY ----
 //  Extra correlations for v3{5} study: 
 // 33rd bin: <4>_{4n,2n|3n,3n} = four4n2n3n3n = <cos(n(4*phi1+2*phi2-3*phi3-3*phi4))>
 // 34th bin: <5>_{3n,3n|2n,2n,2n} = five3n3n2n2n2n = <cos(n(3*phi1+3*phi2-2*phi3-2*phi4-2*phi5))> 
 //  Extra correlations for Teaney-Yan study: 
 // 35th bin: <2>_{5n|5n} = two5n5n = <cos(5n(phi1-phi2)> 
 // 36th bin: <2>_{6n|6n} = two6n6n = <cos(6n(phi1-phi2)> 
 // 37th bin: <3>_{5n|3n,2n} = three5n3n2n = <cos(n(5*phi1-3*phi2-2*phi3)> 
 // 38th bin: <3>_{5n|4n,1n} = three5n4n1n = <cos(n(5*phi1-4*phi2-1*phi3)> 
 // 39th bin: <3>_{6n|3n,3n} = three6n3n3n = <cos(n(6*phi1-3*phi2-3*phi3)> 
 // 40th bin: <3>_{6n|4n,2n} = three6n4n2n = <cos(n(6*phi1-4*phi2-2*phi3)> 
 // 41st bin: <3>_{6n|5n,1n} = three6n5n1n = <cos(n(6*phi1-5*phi2-1*phi3)>
 // 42nd bin: <4>_{6n|3n,2n,1n} = four6n3n2n1n = <cos(n(6*phi1-3*phi2-2*phi3-1*phi4)>
 // 43rd bin: <4>_{3n,2n|3n,2n} = four3n2n3n2n = <cos(n(3*phi1+2*phi2-3*phi3-2*phi4)>
 // 44th bin: <4>_{4n,1n|3n,2n} = four4n1n3n2n = <cos(n(4*phi1+1*phi2-3*phi3-2*phi4)>
 // 45th bin: <4>_{3n,3n|3n,3n} = four3n3n3n3n = <cos(3n*(phi1+phi2-phi3-phi4))> 
 // 46th bin: <4>_{4n,2n|3n,3n} = four4n2n3n3n = <cos(n(4*phi1+2*phi2-3*phi3-3*phi4)>
 // 47th bin: <4>_{5n,1n|3n,3n} = four5n1n3n3n = <cos(n(5*phi1+1*phi2-3*phi3-3*phi4)>
 // 48th bin: <4>_{4n,2n|4n,2n} = four4n2n4n2n = <cos(n(4*phi1+2*phi2-4*phi3-2*phi4)> 
 // 49th bin: <4>_{5n,1n|4n,2n} = four5n1n4n2n = <cos(n(5*phi1+1*phi2-4*phi3-2*phi4)>
 // 50th bin: <4>_{5n|3n,1n,1n} = four5n3n1n1n = <cos(n(5*phi1-3*phi2-1*phi3-1*phi4)>
 // 51st bin: <4>_{5n|2n,2n,1n} = four5n2n2n1n = <cos(n(5*phi1-2*phi2-2*phi3-1*phi4)>
 // 52nd bin: <4>_{5n,1n|5n,1n} = four5n1n5n1n = <cos(n(5*phi1+1*phi2-5*phi3-1*phi4)>
 // 53rd bin: <5>_{3n,3n|3n,2n,1n} = five3n3n3n2n1n = <cos(n(3*phi1+3*phi2-3*phi3-2*phi4-1*phi5)>
 // 54th bin: <5>_{4n,2n|3n,2n,1n} = five4n2n3n2n1n = <cos(n(4*phi1+2*phi2-3*phi3-2*phi4-1*phi5)>
 // 55th bin: <5>_{3n,2n|3n,1n,1n} = five3n2n3n1n1n = <cos(n(3*phi1+2*phi2-3*phi3-1*phi4-1*phi5)>
 // 56th bin: <5>_{3n,2n|2n,2n,1n} = five3n2n2n2n1n = <cos(n(3*phi1+2*phi2-2*phi3-2*phi4-1*phi5)>
 // 57th bin: <5>_{5n,1n|3n,2n,1n} = five5n1n3n2n1n = <cos(n(5*phi1+1*phi2-3*phi3-2*phi4-1*phi5)>
 // 58th bin: <6>_{3n,2n,1n|3n,2n,1n} = six3n2n1n3n2n1n = <cos(n(3*phi1+2*phi2+1*phi3-3*phi4-2*phi5-1*phi6)>
 //  Extra correlations for Teaney-Yan study (B): 
 // 59th bin: <4>_{6n|4n,1n,1n} = four6n4n1n1n = <cos(n(6*phi1-4*phi2-1*phi3-1*phi4)>
 // 60th bin: <4>_{6n|2n,2n,2n} = four6n2n2n2n = <cos(n(6*phi1-2*phi2-2*phi3-2*phi4)>
 // 61st bin: <5>_{6n|2n,2n,1n,1n} = five6n2n2n1n1n = <cos(n(6*phi1-2*phi2-2*phi3-1*phi4-1*phi5)>
 // 62nd bin: <5>_{4n,1n,1n|3n,3n} = five4n1n1n3n3n = <cos(n(4*phi1+1*phi2+1*phi3-3*phi4-3*phi5)>
 // 63rd bin: <6>_{3n,3n|2n,2n,1n,1n} = six3n3n2n2n1n1n = <cos(n(3*phi1+3*phi2-2*phi3-2*phi4-1*phi5-1*phi6)>
 // --------------------------------------------------------------------------------------------------------------------

 // Multiplicity of an event: 
 Double_t dMult = (*fSpk)(0,0);
 // Real parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n, 4n, 5n and 6n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 Double_t dReQ3n = (*fReQ)(2,0);
 Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dReQ5n = (*fReQ)(4,0); 
 Double_t dReQ6n = (*fReQ)(5,0);
 // Imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n, 4n, 5n and 6n:
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 Double_t dImQ3n = (*fImQ)(2,0);
 Double_t dImQ4n = (*fImQ)(3,0);
 Double_t dImQ5n = (*fImQ)(4,0); 
 Double_t dImQ6n = (*fImQ)(5,0);
  
 // Real parts of expressions involving various combinations of Q-vectors which appears
 // simultaneously in several equations for multiparticle correlations bellow: 
 // Re[Q_{2n}Q_{n}^*Q_{n}^*]
 Double_t reQ2nQ1nstarQ1nstar = pow(dReQ1n,2.)*dReQ2n+2.*dReQ1n*dImQ1n*dImQ2n-pow(dImQ1n,2.)*dReQ2n; 
 // Re[Q_{6n}Q_{3n}^*Q_{3n}^*]
 Double_t reQ6nQ3nstarQ3nstar = pow(dReQ3n,2.)*dReQ6n+2.*dReQ3n*dImQ3n*dImQ6n-pow(dImQ3n,2.)*dReQ6n;  
 // Re[Q_{4n}Q_{2n}^*Q_{2n}^*]
 Double_t reQ4nQ2nstarQ2nstar = pow(dReQ2n,2.)*dReQ4n+2.*dReQ2n*dImQ2n*dImQ4n-pow(dImQ2n,2.)*dReQ4n;
 // Re[Q_{4n}Q_{3n}^*Q_{n}^*]
 Double_t reQ4nQ3nstarQ1nstar = dReQ4n*(dReQ3n*dReQ1n-dImQ3n*dImQ1n)+dImQ4n*(dReQ3n*dImQ1n+dImQ3n*dReQ1n);
 // Re[Q_{3n}Q_{2n}^*Q_{n}^*]
 Double_t reQ3nQ2nstarQ1nstar = dReQ3n*dReQ2n*dReQ1n-dReQ3n*dImQ2n*dImQ1n+dImQ3n*dReQ2n*dImQ1n
                              + dImQ3n*dImQ2n*dReQ1n; 
 // Re[Q_{5n}Q_{3n}^*Q_{2n}^*]
 Double_t reQ5nQ3nstarQ2nstar = dReQ5n*dReQ2n*dReQ3n-dReQ5n*dImQ2n*dImQ3n+dImQ5n*dReQ2n*dImQ3n
                              + dImQ5n*dImQ2n*dReQ3n;                             
 // Re[Q_{5n}Q_{4n}^*Q_{1n}^*]
 Double_t reQ5nQ4nstarQ1nstar = dReQ5n*dReQ4n*dReQ1n-dReQ5n*dImQ4n*dImQ1n+dImQ5n*dReQ4n*dImQ1n
                              + dImQ5n*dImQ4n*dReQ1n;                              
 // Re[Q_{6n}Q_{5n}^*Q_{1n}^*]                              
 Double_t reQ6nQ5nstarQ1nstar = dReQ6n*dReQ5n*dReQ1n-dReQ6n*dImQ5n*dImQ1n+dImQ6n*dReQ5n*dImQ1n
                              + dImQ6n*dImQ5n*dReQ1n;
 // Re[Q_{6n}Q_{4n}^*Q_{2n}^*]                              
 Double_t reQ6nQ4nstarQ2nstar = dReQ6n*dReQ4n*dReQ2n-dReQ6n*dImQ4n*dImQ2n+dImQ6n*dReQ4n*dImQ2n
                              + dImQ6n*dImQ4n*dReQ2n;
 // Re[Q_{3n}Q_{n}Q_{2n}^*Q_{2n}^*]
 Double_t reQ3nQ1nQ2nstarQ2nstar = (pow(dReQ2n,2.)-pow(dImQ2n,2.))*(dReQ3n*dReQ1n-dImQ3n*dImQ1n) 
                                 + 2.*dReQ2n*dImQ2n*(dReQ3n*dImQ1n+dImQ3n*dReQ1n);                                
 // Re[Q_{3n}Q_{n}^*Q_{n}^*Q_{n}^*]
 Double_t reQ3nQ1nstarQ1nstarQ1nstar = dReQ3n*pow(dReQ1n,3)-3.*dReQ1n*dReQ3n*pow(dImQ1n,2)
                                     + 3.*dImQ1n*dImQ3n*pow(dReQ1n,2)-dImQ3n*pow(dImQ1n,3);
 // Re[Q_{6n}Q_{2n}^*Q_{2n}^*Q_{2n}^*]
 Double_t reQ6nQ2nstarQ2nstarQ2nstar = dReQ6n*pow(dReQ2n,3)-3.*dReQ2n*dReQ6n*pow(dImQ2n,2)
                                     + 3.*dImQ2n*dImQ6n*pow(dReQ2n,2)-dImQ6n*pow(dImQ2n,3);
 // Re[Q_{4n}Q_{2n}^*Q_{n}^*Q_{n}^*]
 Double_t reQ4nQ2nstarQ1nstarQ1nstar = (dReQ4n*dReQ2n+dImQ4n*dImQ2n)*(pow(dReQ1n,2)-pow(dImQ1n,2)) 
                                     + 2.*dReQ1n*dImQ1n*(dImQ4n*dReQ2n-dReQ4n*dImQ2n);  
 // Re[Q_{4n}Q_{2n}^*Q_{3n}^*Q_{3n}^*]
 Double_t reQ4nQ2nQ3nstarQ3nstar = (dReQ4n*dReQ2n-dImQ4n*dImQ2n)*(dReQ3n*dReQ3n-dImQ3n*dImQ3n)
                                 + 2.*(dReQ4n*dImQ2n+dImQ4n*dReQ2n)*dReQ3n*dImQ3n;                    
 // Re[Q_{4n}Q_{n}Q_{3n}^*Q_{2n}^*]
 Double_t reQ4nQ1nQ3nstarQ2nstar = dImQ1n*dImQ2n*dImQ3n*dImQ4n+dImQ3n*dImQ4n*dReQ1n*dReQ2n 
                                 + dImQ2n*dImQ4n*dReQ1n*dReQ3n-dImQ1n*dImQ4n*dReQ2n*dReQ3n
                                 - dImQ2n*dImQ3n*dReQ1n*dReQ4n+dImQ1n*dImQ3n*dReQ2n*dReQ4n 
                                 + dImQ1n*dImQ2n*dReQ3n*dReQ4n+dReQ1n*dReQ2n*dReQ3n*dReQ4n;
 // Re[Q_{5n}Q_{n}Q_{4n}^*Q_{2n}^*]
 Double_t reQ5nQ1nQ4nstarQ2nstar = dImQ1n*dImQ2n*dImQ4n*dImQ5n+dImQ4n*dImQ5n*dReQ1n*dReQ2n 
                                 + dImQ2n*dImQ5n*dReQ1n*dReQ4n-dImQ1n*dImQ5n*dReQ2n*dReQ4n
                                 - dImQ2n*dImQ4n*dReQ1n*dReQ5n+dImQ1n*dImQ4n*dReQ2n*dReQ5n 
                                 + dImQ1n*dImQ2n*dReQ4n*dReQ5n+dReQ1n*dReQ2n*dReQ4n*dReQ5n;                                  
 // Re[Q_{5n}Q_{n}Q_{3n}^*Q_{3n}^*]                                  
 Double_t reQ5nQ1nQ3nstarQ3nstar = dImQ1n*pow(dImQ3n,2.)*dImQ5n+2.*dImQ3n*dImQ5n*dReQ1n*dReQ3n
                                 - dImQ1n*dImQ5n*pow(dReQ3n,2.)-pow(dImQ3n,2.)*dReQ1n*dReQ5n 
                                 + 2.*dImQ1n*dImQ3n*dReQ3n*dReQ5n+dReQ1n*pow(dReQ3n,2.)*dReQ5n;
 // Re[Q_{5n}Q_{3n}^*Q_{n}^*Q_{n}^*]                                  
 Double_t reQ5nQ3nstarQ1nstarQ1nstar = -pow(dImQ1n,2.)*dImQ3n*dImQ5n+dImQ3n*dImQ5n*pow(dReQ1n,2.)
                                     + 2.*dImQ1n*dImQ5n*dReQ1n*dReQ3n-2.*dImQ1n*dImQ3n*dReQ1n*dReQ5n 
                                     - pow(dImQ1n,2.)*dReQ3n*dReQ5n+pow(dReQ1n,2.)*dReQ3n*dReQ5n;                     
 // Re[Q_{5n}Q_{2n}^*Q_{2n}^*Q_{n}^*]                                  
 Double_t reQ5nQ2nstarQ2nstarQ1nstar = -pow(dImQ2n,2.)*dImQ1n*dImQ5n+dImQ1n*dImQ5n*pow(dReQ2n,2.)
                                     + 2.*dImQ2n*dImQ5n*dReQ2n*dReQ1n-2.*dImQ2n*dImQ1n*dReQ2n*dReQ5n 
                                     - pow(dImQ2n,2.)*dReQ1n*dReQ5n+pow(dReQ2n,2.)*dReQ1n*dReQ5n;                                     
 // Re[Q_{6n}Q_{4n}^*Q_{n}^*Q_{n}^*]                                  
 Double_t reQ6nQ4nstarQ1nstarQ1nstar = -pow(dImQ1n,2.)*dImQ4n*dImQ6n+dImQ4n*dImQ6n*pow(dReQ1n,2.) 
                                     +  2.*dImQ1n*dImQ6n*dReQ1n*dReQ4n-2.*dImQ1n*dImQ4n*dReQ1n*dReQ6n 
                                     -  pow(dImQ1n,2.)*dReQ4n*dReQ6n+pow(dReQ1n,2.)*dReQ4n*dReQ6n;
 // |Q_{2n}|^2 |Q_{n}|^2
 Double_t dQ2nQ1nQ2nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.));
 // |Q_{4n}|^2 |Q_{2n}|^2
 Double_t dQ4nQ2nQ4nstarQ2nstar = (pow(dReQ4n,2.)+pow(dImQ4n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.));
 // |Q_{3n}|^2 |Q_{2n}|^2
 Double_t dQ3nQ2nQ3nstarQ2nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ3n,2.)+pow(dImQ3n,2.));
 // |Q_{5n}|^2 |Q_{n}|^2
 Double_t dQ5nQ1nQ5nstarQ1nstar = (pow(dReQ5n,2.)+pow(dImQ5n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.));
 // Re[Q_{2n}Q_{n}Q_{n}^*Q_{n}^*Q_{n}^*]
 Double_t reQ2nQ1nQ1nstarQ1nstarQ1nstar = (dReQ2n*dReQ1n-dImQ2n*dImQ1n)*(pow(dReQ1n,3)-3.*dReQ1n*pow(dImQ1n,2))
                                        + (dReQ2n*dImQ1n+dReQ1n*dImQ2n)*(3.*dImQ1n*pow(dReQ1n,2)-pow(dImQ1n,3)); 
 // Re[Q_{2n}Q_{2n}Q_{2n}^*Q_{n}^*Q_{n}^*]
 Double_t reQ2nQ2nQ2nstarQ1nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))
                                        * (dReQ2n*(pow(dReQ1n,2.)-pow(dImQ1n,2.)) + 2.*dImQ2n*dReQ1n*dImQ1n);
 // Re[Q_{4n}Q_{n}^*Q_{n}^*Q_{n}^*Q_{n}^*]
 Double_t reQ4nQ1nstarQ1nstarQ1nstarQ1nstar = pow(dReQ1n,4.)*dReQ4n-6.*pow(dReQ1n,2.)*dReQ4n*pow(dImQ1n,2.)
                                            + pow(dImQ1n,4.)*dReQ4n+4.*pow(dReQ1n,3.)*dImQ1n*dImQ4n
                                            - 4.*pow(dImQ1n,3.)*dReQ1n*dImQ4n;
 // Re[Q_{3n}Q_{n}Q_{2n}^*Q_{n}^*Q_{n}^*]
 Double_t reQ3nQ1nQ2nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)+pow(dImQ1n,2.))
                                        * (dReQ1n*dReQ2n*dReQ3n-dReQ3n*dImQ1n*dImQ2n
                                        + dReQ2n*dImQ1n*dImQ3n+dReQ1n*dImQ2n*dImQ3n);
 // Re[Q_{6n}Q_{n}Q_{3n}^*Q_{2n}^*Q_{n}^*]
 Double_t reQ6nQ3nstarQ2nstarQ1nstar = dReQ1n*dReQ2n*dReQ3n*dReQ6n-dReQ3n*dReQ6n*dImQ1n*dImQ2n
                                     - dReQ2n*dReQ6n*dImQ1n*dImQ3n-dReQ1n*dReQ6n*dImQ2n*dImQ3n
                                     + dReQ2n*dReQ3n*dImQ1n*dImQ6n+dReQ1n*dReQ3n*dImQ2n*dImQ6n 
                                     + dReQ1n*dReQ2n*dImQ3n*dImQ6n-dImQ1n*dImQ2n*dImQ3n*dImQ6n;
 // Re[Q_{3n}Q_{3n}Q_{3n}^*Q_{2n}^*Q_{n}^*]
 Double_t reQ3nQ3nQ3nstarQ2nstarQ1nstar = (pow(dImQ3n,2.)+pow(dReQ3n,2.))
                                        * (dImQ2n*dImQ3n*dReQ1n+dImQ1n*dImQ3n*dReQ2n
                                        - dImQ1n*dImQ2n*dReQ3n+dReQ1n*dReQ2n*dReQ3n);   
 // Re[Q_{3n}Q_{3n}Q_{2n}^*Q_{2n}^*Q_{2n}^*]
 Double_t reQ3nQ3nQ2nstarQ2nstarQ2nstar = pow(dReQ2n,3.)*pow(dReQ3n,2.) 
                                        - 3.*dReQ2n*pow(dReQ3n,2.)*pow(dImQ2n,2.)
                                        + 6.*pow(dReQ2n,2.)*dReQ3n*dImQ2n*dImQ3n 
                                        - 2.*dReQ3n*pow(dImQ2n,3.)*dImQ3n-pow(dReQ2n,3.)*pow(dImQ3n,2.) 
                                        + 3.*dReQ2n*pow(dImQ2n,2.)*pow(dImQ3n,2.);
 // Re[Q_{4n}Q_{2n}Q_{3n}^*Q_{2n}^*Q_{n}^*]
 Double_t reQ4nQ2nQ3nstarQ2nstarQ1nstar = (pow(dImQ2n,2.)+pow(dReQ2n,2.))
                                        * (dImQ3n*dImQ4n*dReQ1n+dImQ1n*dImQ4n*dReQ3n 
                                        - dImQ1n*dImQ3n*dReQ4n+dReQ1n*dReQ3n*dReQ4n);
 // Re[Q_{3n}Q_{2n}Q_{3n}^*Q_{n}^*Q_{n}^*]
 Double_t reQ3nQ2nQ3nstarQ1nstarQ1nstar = -(pow(dImQ3n,2.)+pow(dReQ3n,2.))
                                        * (-2.*dImQ1n*dImQ2n*dReQ1n+pow(dImQ1n,2.)*dReQ2n-pow(dReQ1n,2.)*dReQ2n);                              
 // Re[Q_{3n}Q_{2n}Q_{2n}^*Q_{2n}^*Q_{n}^*]
 Double_t reQ3nQ2nQ2nstarQ2nstarQ1nstar = (pow(dImQ2n,2.)+pow(dReQ2n,2.))
                                        * (dImQ2n*dImQ3n*dReQ1n+dImQ1n*dImQ3n*dReQ2n 
                                        - dImQ1n*dImQ2n*dReQ3n+dReQ1n*dReQ2n*dReQ3n);
 // Re[Q_{5n}Q_{n}Q_{3n}^*Q_{2n}^*Q_{n}^*]
 Double_t reQ5nQ1nQ3nstarQ2nstarQ1nstar = (pow(dImQ1n,2.)+pow(dReQ1n,2.))
                                        * (dImQ3n*dImQ5n*dReQ2n+dImQ2n*dImQ5n*dReQ3n 
                                        - dImQ2n*dImQ3n*dReQ5n+dReQ2n*dReQ3n*dReQ5n);   
 // Re[Q_{2n}Q_{2n}Q_{n}^*Q_{n}^*Q_{n}^*Q_{n}^*]
 Double_t reQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)*dReQ2n-2.*dReQ1n*dReQ2n*dImQ1n-dReQ2n*pow(dImQ1n,2.)
                                               + dImQ2n*pow(dReQ1n,2.)+2.*dReQ1n*dImQ1n*dImQ2n-pow(dImQ1n,2.)*dImQ2n)
                                               * (pow(dReQ1n,2.)*dReQ2n+2.*dReQ1n*dReQ2n*dImQ1n-dReQ2n*pow(dImQ1n,2.)
                                               - dImQ2n*pow(dReQ1n,2.)+2.*dReQ1n*dImQ1n*dImQ2n+pow(dImQ1n,2.)*dImQ2n); 
 // Re[Q_{3n}Q_{n}Q_{n}^*Q_{n}^*Q_{n}^*Q_{n}^*]
 Double_t reQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)+pow(dImQ1n,2.))
                                               * (pow(dReQ1n,3.)*dReQ3n-3.*dReQ1n*dReQ3n*pow(dImQ1n,2.)
                                               + 3.*pow(dReQ1n,2.)*dImQ1n*dImQ3n-pow(dImQ1n,3.)*dImQ3n);
 // |Q_{2n}|^2 |Q_{n}|^4
 Double_t dQ2nQ1nQ1nQ2nstarQ1nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))*pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.); 
 // |Q_{3n}|^2 |Q_{2n}|^2 |Q_{n}|^2
 Double_t dQ3nQ2nQ1nQ3nstarQ2nstarQ1nstar = (pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                                          * (pow(dReQ1n,2.)+pow(dImQ1n,2.));
 // Re[Q_{2n}Q_{n}Q_{n}Q_{n}^*Q_{n}^*Q_{n}^*Q_{n}^*]
 Double_t reQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                                                  * (pow(dReQ1n,2.)*dReQ2n-dReQ2n*pow(dImQ1n,2.)
                                                  + 2.*dReQ1n*dImQ1n*dImQ2n);
    
 // Results for multiparticle azimuthal correlations:
 // 2-particle:
 Double_t two1n1n = 0.; // <cos(n(phi1-phi2))>
 Double_t two2n2n = 0.; // <cos(2n(phi1-phi2))>
 Double_t two3n3n = 0.; // <cos(3n(phi1-phi2))>
 Double_t two4n4n = 0.; // <cos(4n(phi1-phi2))>
 if(dMult>1)
 {
  two1n1n = (pow(dReQ1n,2.)+pow(dImQ1n,2.)-dMult)/(dMult*(dMult-1.)); 
  two2n2n = (pow(dReQ2n,2.)+pow(dImQ2n,2.)-dMult)/(dMult*(dMult-1.)); 
  two3n3n = (pow(dReQ3n,2.)+pow(dImQ3n,2.)-dMult)/(dMult*(dMult-1.)); 
  two4n4n = (pow(dReQ4n,2.)+pow(dImQ4n,2.)-dMult)/(dMult*(dMult-1.)); 
  // Average 2-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(1,two1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(2,two2n2n);
  fIntFlowCorrelationsAllEBE->SetBinContent(3,two3n3n);
  fIntFlowCorrelationsAllEBE->SetBinContent(4,two4n4n);         
  // Average 2-particle correlations for all events:      
  fIntFlowCorrelationsAllPro->Fill(0.5,two1n1n,dMult*(dMult-1.));
  fIntFlowCorrelationsAllPro->Fill(1.5,two2n2n,dMult*(dMult-1.)); 
  fIntFlowCorrelationsAllPro->Fill(2.5,two3n3n,dMult*(dMult-1.)); 
  fIntFlowCorrelationsAllPro->Fill(3.5,two4n4n,dMult*(dMult-1.)); 
  // Store separetately <2>:
  fIntFlowCorrelationsEBE->SetBinContent(1,two1n1n); // <2>  
  // Testing other multiplicity weights:
  Double_t mWeight2p = 0.;
  if(!strcmp(fMultiplicityWeight->Data(),"combinations"))
  {
   mWeight2p = dMult*(dMult-1.);
  } else if(!strcmp(fMultiplicityWeight->Data(),"unit"))
    {
     mWeight2p = 1.;    
    } else if(!strcmp(fMultiplicityWeight->Data(),"multiplicity"))
      {
       mWeight2p = dMult;           
      }          
  fIntFlowEventWeightsForCorrelationsEBE->SetBinContent(1,mWeight2p); // eW_<2>
  fIntFlowCorrelationsPro->Fill(0.5,two1n1n,mWeight2p);
  fIntFlowSquaredCorrelationsPro->Fill(0.5,two1n1n*two1n1n,mWeight2p);
  if(fCalculateCumulantsVsM)
  {
   fIntFlowCorrelationsVsMPro[0]->Fill(dMult+0.5,two1n1n,mWeight2p);
   fIntFlowSquaredCorrelationsVsMPro[0]->Fill(dMult+0.5,two1n1n*two1n1n,mWeight2p);
  } 
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[0]->Fill(dMult+0.5,two1n1n,mWeight2p);
   fIntFlowCorrelationsAllVsMPro[1]->Fill(dMult+0.5,two2n2n,mWeight2p);
   fIntFlowCorrelationsAllVsMPro[2]->Fill(dMult+0.5,two3n3n,mWeight2p);
   fIntFlowCorrelationsAllVsMPro[3]->Fill(dMult+0.5,two4n4n,mWeight2p);
  }  
 } // end of if(dMult>1)
 
 // 3-particle:
 Double_t three2n1n1n = 0.; // <cos(n(2*phi1-phi2-phi3))>
 Double_t three3n2n1n = 0.; // <cos(n(3*phi1-2*phi2-phi3))>
 Double_t three4n2n2n = 0.; // <cos(n(4*phi1-2*phi2-2*phi3))>
 Double_t three4n3n1n = 0.; // <cos(n(4*phi1-3*phi2-phi3))> 
 if(dMult>2)
 {
  three2n1n1n = (reQ2nQ1nstarQ1nstar-2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
              - (pow(dReQ2n,2.)+pow(dImQ2n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.));                     
  three3n2n1n = (reQ3nQ2nstarQ1nstar-(pow(dReQ3n,2.)+pow(dImQ3n,2.))
              - (pow(dReQ2n,2.)+pow(dImQ2n,2.))
              - (pow(dReQ1n,2.)+pow(dImQ1n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.));
  three4n2n2n = (reQ4nQ2nstarQ2nstar-2.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
              - (pow(dReQ4n,2.)+pow(dImQ4n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.)); 
  three4n3n1n = (reQ4nQ3nstarQ1nstar-(pow(dReQ4n,2.)+pow(dImQ4n,2.))
              - (pow(dReQ3n,2.)+pow(dImQ3n,2.))
              - (pow(dReQ1n,2.)+pow(dImQ1n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.));              
  // Average 3-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(6,three2n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(7,three3n2n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(8,three4n2n2n);
  fIntFlowCorrelationsAllEBE->SetBinContent(9,three4n3n1n);
  // Average 3-particle correlations for all events:                
  fIntFlowCorrelationsAllPro->Fill(5.5,three2n1n1n,dMult*(dMult-1.)*(dMult-2.)); 
  fIntFlowCorrelationsAllPro->Fill(6.5,three3n2n1n,dMult*(dMult-1.)*(dMult-2.));
  fIntFlowCorrelationsAllPro->Fill(7.5,three4n2n2n,dMult*(dMult-1.)*(dMult-2.)); 
  fIntFlowCorrelationsAllPro->Fill(8.5,three4n3n1n,dMult*(dMult-1.)*(dMult-2.));  
  // Average 3-particle correlations vs M for all events:                
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[5]->Fill(dMult+0.5,three2n1n1n,dMult*(dMult-1.)*(dMult-2.));
   fIntFlowCorrelationsAllVsMPro[6]->Fill(dMult+0.5,three3n2n1n,dMult*(dMult-1.)*(dMult-2.));
   fIntFlowCorrelationsAllVsMPro[7]->Fill(dMult+0.5,three4n2n2n,dMult*(dMult-1.)*(dMult-2.));
   fIntFlowCorrelationsAllVsMPro[8]->Fill(dMult+0.5,three4n3n1n,dMult*(dMult-1.)*(dMult-2.));
  }    
 } // end of if(dMult>2)
 
 // 4-particle:
 Double_t four1n1n1n1n = 0.; // <cos(n(phi1+phi2-phi3-phi4))>
 Double_t four2n2n2n2n = 0.; // <cos(2n(phi1+phi2-phi3-phi4))>
 Double_t four2n1n2n1n = 0.; // <cos(n(2*phi1+phi2-2*phi3-phi4))> 
 Double_t four3n1n1n1n = 0.; // <cos(n(3*phi1-phi2-phi3-phi4))> 
 Double_t four4n2n1n1n = 0.; // <cos(n(4*phi1-2*phi2-phi3-phi4))> 
 Double_t four3n1n2n2n = 0.; // <cos(n(3*phi1+phi2-2*phi3-2*phi4))> 
 Double_t four3n1n3n1n = 0.; // <cos(n(3*phi1+phi2-3*phi3-phi4))>    
 if(dMult>3)
 {
  four1n1n1n1n = (2.*dMult*(dMult-3.)+pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)-4.*(dMult-2.)*(pow(dReQ1n,2.)
               + pow(dImQ1n,2.))-2.*reQ2nQ1nstarQ1nstar+(pow(dReQ2n,2.)+pow(dImQ2n,2.)))
               / (dMult*(dMult-1)*(dMult-2.)*(dMult-3.));     
  four2n2n2n2n = (2.*dMult*(dMult-3.)+pow((pow(dReQ2n,2.)+pow(dImQ2n,2.)),2.)-4.*(dMult-2.)*(pow(dReQ2n,2.)
               + pow(dImQ2n,2.))-2.*reQ4nQ2nstarQ2nstar+(pow(dReQ4n,2.)+pow(dImQ4n,2.)))
               / (dMult*(dMult-1)*(dMult-2.)*(dMult-3.));
  four2n1n2n1n = (dQ2nQ1nQ2nstarQ1nstar-2.*reQ3nQ2nstarQ1nstar-2.*reQ2nQ1nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - ((dMult-5.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
               + (dMult-4.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-(pow(dReQ3n,2.)+pow(dImQ3n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               + (dMult-6.)/((dMult-1.)*(dMult-2.)*(dMult-3.));
  four3n1n1n1n = (reQ3nQ1nstarQ1nstarQ1nstar-3.*reQ3nQ2nstarQ1nstar-3.*reQ2nQ1nstarQ1nstar
               + 2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+3.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
               + 6.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-6.*dMult)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  four4n2n1n1n = (reQ4nQ2nstarQ1nstarQ1nstar-2.*reQ4nQ3nstarQ1nstar-reQ4nQ2nstarQ2nstar-2.*reQ3nQ2nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - (reQ2nQ1nstarQ1nstar-2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))-2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
               - 3.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - 6./((dMult-1.)*(dMult-2.)*(dMult-3.));
  four3n1n2n2n = (reQ3nQ1nQ2nstarQ2nstar-reQ4nQ2nstarQ2nstar-reQ4nQ3nstarQ1nstar-2.*reQ3nQ2nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - (2.*reQ2nQ1nstarQ1nstar-(pow(dReQ4n,2.)+pow(dImQ4n,2.))-2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
               - 4.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - 6./((dMult-1.)*(dMult-2.)*(dMult-3.)); 
  four3n1n3n1n = ((pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
               - 2.*reQ4nQ3nstarQ1nstar-2.*reQ3nQ2nstarQ1nstar
               + pow(dReQ4n,2.)+pow(dImQ4n,2.)-(dMult-4.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
               + pow(dReQ2n,2.)+pow(dImQ2n,2.)-(dMult-4.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
               + dMult*(dMult-6.))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));           
  // Average 4-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(11,four1n1n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(12,four2n1n2n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(13,four2n2n2n2n);
  fIntFlowCorrelationsAllEBE->SetBinContent(14,four3n1n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(15,four3n1n3n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(16,four3n1n2n2n);
  fIntFlowCorrelationsAllEBE->SetBinContent(17,four4n2n1n1n);       
  // Average 4-particle correlations for all events:                
  fIntFlowCorrelationsAllPro->Fill(10.5,four1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(11.5,four2n1n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(12.5,four2n2n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(13.5,four3n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(14.5,four3n1n3n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(15.5,four3n1n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));  
  fIntFlowCorrelationsAllPro->Fill(16.5,four4n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));  
  // Average 4-particle correlations vs M for all events:                
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[10]->Fill(dMult+0.5,four1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[11]->Fill(dMult+0.5,four2n1n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[12]->Fill(dMult+0.5,four2n2n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[13]->Fill(dMult+0.5,four3n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[14]->Fill(dMult+0.5,four3n1n3n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[15]->Fill(dMult+0.5,four3n1n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[16]->Fill(dMult+0.5,four4n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  }       
  // Store separetately <4>:
  fIntFlowCorrelationsEBE->SetBinContent(2,four1n1n1n1n); // <4>
  // Testing other multiplicity weights:
  Double_t mWeight4p = 0.;
  if(!strcmp(fMultiplicityWeight->Data(),"combinations"))
  {
   mWeight4p = dMult*(dMult-1.)*(dMult-2.)*(dMult-3.);
  } else if(!strcmp(fMultiplicityWeight->Data(),"unit"))
    {
     mWeight4p = 1.;    
    } else if(!strcmp(fMultiplicityWeight->Data(),"multiplicity"))
      {
       mWeight4p = dMult;           
      }      
  fIntFlowEventWeightsForCorrelationsEBE->SetBinContent(2,mWeight4p); // eW_<4>
  fIntFlowCorrelationsPro->Fill(1.5,four1n1n1n1n,mWeight4p);
  fIntFlowSquaredCorrelationsPro->Fill(1.5,four1n1n1n1n*four1n1n1n1n,mWeight4p);
  if(fCalculateCumulantsVsM)
  {
   fIntFlowCorrelationsVsMPro[1]->Fill(dMult+0.5,four1n1n1n1n,mWeight4p);
   fIntFlowSquaredCorrelationsVsMPro[1]->Fill(dMult+0.5,four1n1n1n1n*four1n1n1n1n,mWeight4p);
  }   
 } // end of if(dMult>3)

 // 5-particle:
 Double_t five2n1n1n1n1n = 0.; // <cos(n(2*phi1+phi2-phi3-phi4-phi5))>
 Double_t five2n2n2n1n1n = 0.; // <cos(n(2*phi1+2*phi2-2*phi3-phi4-phi5))>
 Double_t five3n1n2n1n1n = 0.; // <cos(n(3*phi1+phi2-2*phi3-phi4-phi5))>
 Double_t five4n1n1n1n1n = 0.; // <cos(n(4*phi1-phi2-phi3-phi4-phi5))>
 if(dMult>4)
 {            
  five2n1n1n1n1n = (reQ2nQ1nQ1nstarQ1nstarQ1nstar-reQ3nQ1nstarQ1nstarQ1nstar+5.*reQ3nQ2nstarQ1nstar
                 - 3.*(dMult-5.)*reQ2nQ1nstarQ1nstar-2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 - 3.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))     
                 + 3.*(dMult-4.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 - 3.*pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                 + 6.*(2.*dMult-5.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-6.*dMult*(dMult-4.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  five2n2n2n1n1n = (reQ2nQ2nQ2nstarQ1nstarQ1nstar-reQ4nQ2nstarQ1nstarQ1nstar-2.*reQ3nQ1nQ2nstarQ2nstar
                 + 3.*reQ4nQ2nstarQ2nstar+8.*reQ3nQ2nstarQ1nstar+2.*reQ4nQ3nstarQ1nstar
                 - 2.*(dMult-6.)*reQ2nQ1nstarQ1nstar
                 - 2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))-4.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 - pow((pow(dReQ2n,2.)+pow(dImQ2n,2.)),2.)
                 + 2.*(3.*dMult-10.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 - 4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 + 4.*(dMult-5.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-4.*dMult*(dMult-6.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  five4n1n1n1n1n = (reQ4nQ1nstarQ1nstarQ1nstarQ1nstar-6.*reQ4nQ2nstarQ1nstarQ1nstar-4.*reQ3nQ1nstarQ1nstarQ1nstar
                 + 8.*reQ4nQ3nstarQ1nstar+3.*reQ4nQ2nstarQ2nstar+12.*reQ3nQ2nstarQ1nstar+12.*reQ2nQ1nstarQ1nstar
                 - 6.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))-8.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 - 12.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-24.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))+24.*dMult)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  five3n1n2n1n1n = (reQ3nQ1nQ2nstarQ1nstarQ1nstar-reQ4nQ2nstarQ1nstarQ1nstar-reQ3nQ1nstarQ1nstarQ1nstar
                 - reQ3nQ1nQ2nstarQ2nstar+4.*reQ4nQ3nstarQ1nstar+reQ4nQ2nstarQ2nstar
                 - (2.*dMult-13.)*reQ3nQ2nstarQ1nstar+7.*reQ2nQ1nstarQ1nstar
                 - 2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
                 + 2.*(dMult-5.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 - 2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 + 2.*(dMult-6.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 - 2.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                 + 2.*(3.*dMult-11.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-4.*dMult*(dMult-6.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));                 
  // Average 5-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(19,five2n1n1n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(20,five2n2n2n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(21,five3n1n2n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(22,five4n1n1n1n1n);        
  // Average 5-particle correlations for all events:                         
  fIntFlowCorrelationsAllPro->Fill(18.5,five2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)); 
  fIntFlowCorrelationsAllPro->Fill(19.5,five2n2n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(20.5,five3n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(21.5,five4n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)); 
  // Average 5-particle correlations vs M for all events:                
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[18]->Fill(dMult+0.5,five2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
   fIntFlowCorrelationsAllVsMPro[19]->Fill(dMult+0.5,five2n2n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
   fIntFlowCorrelationsAllVsMPro[20]->Fill(dMult+0.5,five3n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
   fIntFlowCorrelationsAllVsMPro[21]->Fill(dMult+0.5,five4n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  }    
 } // end of if(dMult>4)
    
 // 6-particle:
 Double_t six1n1n1n1n1n1n = 0.; // <cos(n(phi1+phi2+phi3-phi4-phi5-phi6))>
 Double_t six2n2n1n1n1n1n = 0.; // <cos(n(2*phi1+2*phi2-phi3-phi4-phi5-phi6))>
 Double_t six3n1n1n1n1n1n = 0.; // <cos(n(3*phi1+phi2-phi3-phi4-phi5-phi6))>
 Double_t six2n1n1n2n1n1n = 0.; // <cos(n(2*phi1+phi2+phi3-2*phi4-phi5-phi6))>
 if(dMult>5)
 {
  six1n1n1n1n1n1n = (pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),3.)-6.*reQ2nQ1nQ1nstarQ1nstarQ1nstar
                  + 4.*reQ3nQ1nstarQ1nstarQ1nstar-12.*reQ3nQ2nstarQ1nstar+18.*(dMult-4.)*reQ2nQ1nstarQ1nstar
                  + 9.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                  + 4.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))-9.*(dMult-4.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                  - 9.*(dMult-4.)*pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.) 
                  + 18.*(dMult*dMult-7.*dMult+10.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                  - 6.*dMult*(dMult*dMult-9.*dMult+20.))
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  six2n1n1n2n1n1n = (dQ2nQ1nQ1nQ2nstarQ1nstarQ1nstar-4.*reQ3nQ1nQ2nstarQ1nstarQ1nstar
                  - 4.*reQ2nQ1nQ1nstarQ1nstarQ1nstar-2.*reQ2nQ2nQ2nstarQ1nstarQ1nstar
                  + 4.*reQ4nQ2nstarQ1nstarQ1nstar+4.*reQ3nQ1nQ2nstarQ2nstar+4.*reQ3nQ1nstarQ1nstarQ1nstar
                  - 8.*reQ4nQ3nstarQ1nstar-4.*reQ4nQ2nstarQ2nstar+4.*(2.*dMult-13.)*reQ3nQ2nstarQ1nstar
                  + 2.*(7.*dMult-34.)*reQ2nQ1nstarQ1nstar+4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                  - 4.*(dMult-7.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                  + 4.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))-4.*(dMult-6.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                  + pow((pow(dReQ2n,2.)+pow(dImQ2n,2.)),2.)+(2.*dMult*dMult-27.*dMult+76.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                  - (dMult-12.)*pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                  + 4.*(dMult*dMult-15.*dMult+34.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                  - 2.*dMult*(dMult*dMult-17.*dMult+60.))
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  six2n2n1n1n1n1n = (reQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar-6.*reQ2nQ2nQ2nstarQ1nstarQ1nstar-reQ4nQ1nstarQ1nstarQ1nstarQ1nstar
                  - 8.*reQ2nQ1nQ1nstarQ1nstarQ1nstar+8.*reQ3nQ1nstarQ1nstarQ1nstar+6.*reQ4nQ2nstarQ1nstarQ1nstar
                  + 8.*reQ3nQ1nQ2nstarQ2nstar-40.*reQ3nQ2nstarQ1nstar-8.*reQ4nQ3nstarQ1nstar-9.*reQ4nQ2nstarQ2nstar
                  + 24.*(dMult-4.)*reQ2nQ1nstarQ1nstar+24.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                  + 6.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+16.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                  + 3.*pow((pow(dReQ2n,2.)+pow(dImQ2n,2.)),2.)-12.*(2.*dMult-7.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                  + 12.*pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)-48.*(dMult-3.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                  + 24.*dMult*(dMult-5.))
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  six3n1n1n1n1n1n = (reQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar-6.*reQ3nQ1nQ2nstarQ1nstarQ1nstar+6.*reQ4nQ2nstarQ1nstarQ1nstar
                  - reQ4nQ1nstarQ1nstarQ1nstarQ1nstar-4.*reQ2nQ1nQ1nstarQ1nstarQ1nstar+3.*reQ3nQ1nQ2nstarQ2nstar
                  - 4.*(dMult-5.)*reQ3nQ1nstarQ1nstarQ1nstar-14.*reQ4nQ3nstarQ1nstar
                  - 3.*reQ4nQ2nstarQ2nstar+4.*(3.*dMult-17.)*reQ3nQ2nstarQ1nstar+12.*(dMult-6.)*reQ2nQ1nstarQ1nstar
                  + 12.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))  
                  + 8.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ3n,2.)+pow(dImQ3n,2.))  
                  + 6.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))-8.*(dMult-5.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.)) 
                  - 12.*(dMult-5.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-48.*(dMult-3.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                  + 12.*pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)+24.*dMult*(dMult-5.))
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  // Average 6-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(24,six1n1n1n1n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(25,six2n1n1n2n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(26,six2n2n1n1n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(27,six3n1n1n1n1n1n);
  // Average 6-particle correlations for all events:         
  fIntFlowCorrelationsAllPro->Fill(23.5,six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  fIntFlowCorrelationsAllPro->Fill(24.5,six2n1n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  fIntFlowCorrelationsAllPro->Fill(25.5,six2n2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  fIntFlowCorrelationsAllPro->Fill(26.5,six3n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  // Average 6-particle correlations vs M for all events:                
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[23]->Fill(dMult+0.5,six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
   fIntFlowCorrelationsAllVsMPro[24]->Fill(dMult+0.5,six2n1n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
   fIntFlowCorrelationsAllVsMPro[25]->Fill(dMult+0.5,six2n2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
   fIntFlowCorrelationsAllVsMPro[26]->Fill(dMult+0.5,six3n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  }    
  // Store separetately <6>:
  fIntFlowCorrelationsEBE->SetBinContent(3,six1n1n1n1n1n1n); // <6>
  // Testing other multiplicity weights:
  Double_t mWeight6p = 0.;
  if(!strcmp(fMultiplicityWeight->Data(),"combinations"))
  {
   mWeight6p = dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.);
  } else if(!strcmp(fMultiplicityWeight->Data(),"unit"))
    {
     mWeight6p = 1.;    
    } else if(!strcmp(fMultiplicityWeight->Data(),"multiplicity"))
      {
       mWeight6p = dMult;           
      }
  fIntFlowEventWeightsForCorrelationsEBE->SetBinContent(3,mWeight6p); // eW_<6>
  fIntFlowCorrelationsPro->Fill(2.5,six1n1n1n1n1n1n,mWeight6p);
  fIntFlowSquaredCorrelationsPro->Fill(2.5,six1n1n1n1n1n1n*six1n1n1n1n1n1n,mWeight6p);
  if(fCalculateCumulantsVsM)
  {
   fIntFlowCorrelationsVsMPro[2]->Fill(dMult+0.5,six1n1n1n1n1n1n,mWeight6p);
   fIntFlowSquaredCorrelationsVsMPro[2]->Fill(dMult+0.5,six1n1n1n1n1n1n*six1n1n1n1n1n1n,mWeight6p);
  }    
 } // end of if(dMult>5)
 
 // 7-particle:
 Double_t seven2n1n1n1n1n1n1n = 0.; // <cos(n(2*phi1+phi2+phi3-phi4-phi5-phi6-phi7))>
 if(dMult>6)
 {
  seven2n1n1n1n1n1n1n = (reQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar-4.*pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),3.)
                      - reQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar-2.*reQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar
                      + 9.*reQ2nQ2nQ2nstarQ1nstarQ1nstar+20.*reQ3nQ1nQ2nstarQ1nstarQ1nstar 
                      + 2.*reQ4nQ1nstarQ1nstarQ1nstarQ1nstar-8.*(dMult-8.)*reQ2nQ1nQ1nstarQ1nstarQ1nstar
                      - 18.*reQ4nQ2nstarQ1nstarQ1nstar-14.*reQ3nQ1nQ2nstarQ2nstar
                      + 8.*(dMult-7.)*reQ3nQ1nstarQ1nstarQ1nstar+28.*reQ4nQ3nstarQ1nstar
                      + 12.*reQ4nQ2nstarQ2nstar-8.*(5.*dMult-31.)*reQ3nQ2nstarQ1nstar      
                      + 12.*(dMult*dMult-15.*dMult+46.)*reQ2nQ1nstarQ1nstar
                      - 16.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                      - 6.*pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),2.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                      - 3.*pow(pow(dReQ2n,2.)+pow(dImQ2n,2.),2.)
                      + 12.*(2.*dMult-13.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                      - 12.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+16.*(dMult-6.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                      - 12.*(dMult-8.)*(dMult-4.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                      + 12.*(3.*dMult-14.)*pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),2.)
                      - 24.*(3.*dMult-7.)*(dMult-6.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                      + 24.*dMult*(dMult-5.)*(dMult-6.))
                      / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.));   
  // Average 7-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(29,seven2n1n1n1n1n1n1n);       
  // Average 7-particle correlations for all events:                      
  fIntFlowCorrelationsAllPro->Fill(28.5,seven2n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)
                                                                 *(dMult-4.)*(dMult-5.)*(dMult-6.));
  // Average 7-particle correlations vs M for all events:                
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[28]->Fill(dMult+0.5,seven2n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)
                                                                              *(dMult-4.)*(dMult-5.)*(dMult-6.));
  }    
 } // end of if(dMult>6)
 
 // 8-particle:
 Double_t eight1n1n1n1n1n1n1n1n = 0.; // <cos(n(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8))>
 if(dMult>7)
 {  
  eight1n1n1n1n1n1n1n1n = (pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),4.)-12.*reQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar
                        + 16.*reQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar+6.*reQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar
                        - 12.*reQ4nQ1nstarQ1nstarQ1nstarQ1nstar-36.*reQ2nQ2nQ2nstarQ1nstarQ1nstar
                        - 96.*reQ3nQ1nQ2nstarQ1nstarQ1nstar
                        + 72.*reQ4nQ2nstarQ1nstarQ1nstar+48.*reQ3nQ1nQ2nstarQ2nstar
                        - 64.*(dMult-6.)*reQ3nQ1nstarQ1nstarQ1nstar
                        + 96.*(dMult-6.)*reQ2nQ1nQ1nstarQ1nstarQ1nstar
                        - 96.*reQ4nQ3nstarQ1nstar-36.*reQ4nQ2nstarQ2nstar
                        + 192.*(dMult-6.)*reQ3nQ2nstarQ1nstar
                        - 144.*(dMult-7.)*(dMult-4.)*reQ2nQ1nstarQ1nstar
                        + 64.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                        - 144.*(dMult-6.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                        + 72.*(dMult-7.)*(dMult-4.)*(pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),2.)+pow(dReQ2n,2.)+pow(dImQ2n,2.))
                        - 96.*(dMult-7.)*(dMult-6.)*(dMult-2.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                        + 36.*pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),2.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                        + 9.*pow(pow(dReQ2n,2.)+pow(dImQ2n,2.),2.)
                        - 64.*(dMult-6.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                        + 36.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
                        - 16.*(dMult-6.)*pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),3.)
                        + 24.*dMult*(dMult-7.)*(dMult-6.)*(dMult-5.))
                        / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));   
  // Average 8-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(31,eight1n1n1n1n1n1n1n1n);      
  // Average 8-particle correlations for all events:                       
  fIntFlowCorrelationsAllPro->Fill(30.5,eight1n1n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)
                                                                   *(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
  // Average 8-particle correlations vs M for all events:                
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[30]->Fill(dMult+0.5,eight1n1n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)
                                                                                *(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
  }     
  // Store separetately <8>:
  fIntFlowCorrelationsEBE->SetBinContent(4,eight1n1n1n1n1n1n1n1n); // <8>
  // Testing other multiplicity weights:
  Double_t mWeight8p = 0.;
  if(!strcmp(fMultiplicityWeight->Data(),"combinations"))
  {
   mWeight8p = dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.);
  } else if(!strcmp(fMultiplicityWeight->Data(),"unit"))
    {
     mWeight8p = 1.;    
    } else if(!strcmp(fMultiplicityWeight->Data(),"multiplicity"))
      {
       mWeight8p = dMult;           
      }        
  fIntFlowEventWeightsForCorrelationsEBE->SetBinContent(4,mWeight8p); // eW_<8>
  fIntFlowCorrelationsPro->Fill(3.5,eight1n1n1n1n1n1n1n1n,mWeight8p);
  fIntFlowSquaredCorrelationsPro->Fill(3.5,eight1n1n1n1n1n1n1n1n*eight1n1n1n1n1n1n1n1n,mWeight8p);  
  if(fCalculateCumulantsVsM)
  {
   fIntFlowCorrelationsVsMPro[3]->Fill(dMult+0.5,eight1n1n1n1n1n1n1n1n,mWeight8p);
   fIntFlowSquaredCorrelationsVsMPro[3]->Fill(dMult+0.5,eight1n1n1n1n1n1n1n1n*eight1n1n1n1n1n1n1n1n,mWeight8p);
  }    
 } // end of if(dMult>7) 
 
 // EXTRA correlations for v3{5} study:
 // 4-particle:
 Double_t four4n2n3n3n = 0.; // <cos(n(4*phi1+2*phi2-3*phi3-3*phi4))>
 if(dMult>3.)
 {
  four4n2n3n3n = (reQ4nQ2nQ3nstarQ3nstar-reQ6nQ4nstarQ2nstar-reQ6nQ3nstarQ3nstar
               - 2.*reQ4nQ3nstarQ1nstar-2.*reQ3nQ2nstarQ1nstar
               + (pow(dReQ6n,2.)+pow(dImQ6n,2.))+2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
               + 2.*(2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+(pow(dReQ2n,2.)+pow(dImQ2n,2.))
               + (pow(dReQ1n,2.)+pow(dImQ1n,2.))-3.*dMult))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));               
  fIntFlowCorrelationsAllPro->Fill(32.5,four4n2n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  // Average 4-particle correlations vs M for all events:                
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[32]->Fill(dMult+0.5,four4n2n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  }    
 } // end of if(dMult>3.)
 
 // 5-particle:
 Double_t five3n3n2n2n2n = 0.; // <cos(n(3*phi1+3*phi2-2*phi3-2*phi4-2*phi5))>                                    
 if(dMult>4.)
 {
  five3n3n2n2n2n = (reQ3nQ3nQ2nstarQ2nstarQ2nstar-reQ6nQ2nstarQ2nstarQ2nstar-3.*reQ4nQ2nQ3nstarQ3nstar 
                 - 6.*reQ3nQ1nQ2nstarQ2nstar+2.*reQ6nQ3nstarQ3nstar+3.*reQ6nQ4nstarQ2nstar
                 + 6.*reQ4nQ3nstarQ1nstar+6.*reQ4nQ2nstarQ2nstar
                 + 12.*reQ3nQ2nstarQ1nstar+6.*reQ2nQ1nstarQ1nstar
                 - 2.*((pow(dReQ6n,2.)+pow(dImQ6n,2.)) 
                 + 3.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
                 + 6.*(pow(dReQ3n,2.)+pow(dImQ3n,2.)) 
                 + 9.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 + 6.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-12.*dMult))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(33.5,five3n3n2n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[33]->Fill(dMult+0.5,five3n3n2n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  }     
 } // end of if(dMult>4.)
 
 // EXTRA correlations for Teaney-Yan study: 
 // 2-particle:
 Double_t two5n5n = 0.; // <cos(5n(phi1-phi2))>
 Double_t two6n6n = 0.; // <cos(6n(phi1-phi2))> 
 if(dMult>1)
 {
  two5n5n = (pow(dReQ5n,2.)+pow(dImQ5n,2.)-dMult)/(dMult*(dMult-1.)); 
  two6n6n = (pow(dReQ6n,2.)+pow(dImQ6n,2.)-dMult)/(dMult*(dMult-1.));        
  // Average 2-particle correlations for all events:      
  fIntFlowCorrelationsAllPro->Fill(34.5,two5n5n,dMult*(dMult-1.));
  fIntFlowCorrelationsAllPro->Fill(35.5,two6n6n,dMult*(dMult-1.)); 
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[34]->Fill(dMult+0.5,two5n5n,dMult*(dMult-1.));
   fIntFlowCorrelationsAllVsMPro[35]->Fill(dMult+0.5,two6n6n,dMult*(dMult-1.));
  }       
 } // end of if(dMult>1)
 
 // 3-particle:
 Double_t three5n3n2n = 0.; // <cos(n(5*phi1-3*phi2-2*phi3)> 
 Double_t three5n4n1n = 0.; // <cos(n(5*phi1-4*phi2-1*phi3)> 
 Double_t three6n3n3n = 0.; // <cos(n(6*phi1-3*phi2-3*phi3)> 
 Double_t three6n4n2n = 0.; // <cos(n(6*phi1-4*phi2-2*phi3)> 
 Double_t three6n5n1n = 0.; // <cos(n(6*phi1-5*phi2-1*phi3)> 
 if(dMult>2)
 {
  three5n3n2n = (reQ5nQ3nstarQ2nstar-(pow(dReQ5n,2.)+pow(dImQ5n,2.))
              - (pow(dReQ3n,2.)+pow(dImQ3n,2.))
              - (pow(dReQ2n,2.)+pow(dImQ2n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.));           
  three5n4n1n = (reQ5nQ4nstarQ1nstar-(pow(dReQ5n,2.)+pow(dImQ5n,2.))
              - (pow(dReQ4n,2.)+pow(dImQ4n,2.))
              - (pow(dReQ1n,2.)+pow(dImQ1n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.));                          
  three6n3n3n = (reQ6nQ3nstarQ3nstar-2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
              - (pow(dReQ6n,2.)+pow(dImQ6n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.)); 
  three6n4n2n = (reQ6nQ4nstarQ2nstar-(pow(dReQ6n,2.)+pow(dImQ6n,2.))
              - (pow(dReQ4n,2.)+pow(dImQ4n,2.))
              - (pow(dReQ2n,2.)+pow(dImQ2n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.));
  three6n5n1n = (reQ6nQ5nstarQ1nstar-(pow(dReQ6n,2.)+pow(dImQ6n,2.))
              - (pow(dReQ5n,2.)+pow(dImQ5n,2.))
              - (pow(dReQ1n,2.)+pow(dImQ1n,2.))+2.*dMult)
              / (dMult*(dMult-1.)*(dMult-2.));            
  // Average 3-particle correlations for all events:      
  fIntFlowCorrelationsAllPro->Fill(36.5,three5n3n2n,dMult*(dMult-1.)*(dMult-2.)); // <<cos(n(5*phi1-3*phi2-2*phi3)>> 
  fIntFlowCorrelationsAllPro->Fill(37.5,three5n4n1n,dMult*(dMult-1.)*(dMult-2.)); // <<cos(n(5*phi1-4*phi2-1*phi3)>> 
  fIntFlowCorrelationsAllPro->Fill(38.5,three6n3n3n,dMult*(dMult-1.)*(dMult-2.)); // <<cos(n(6*phi1-3*phi2-3*phi3)>> 
  fIntFlowCorrelationsAllPro->Fill(39.5,three6n4n2n,dMult*(dMult-1.)*(dMult-2.)); // <<cos(n(6*phi1-4*phi2-2*phi3)>>
  fIntFlowCorrelationsAllPro->Fill(40.5,three6n5n1n,dMult*(dMult-1.)*(dMult-2.)); // <<cos(n(6*phi1-5*phi2-1*phi3)>>
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[36]->Fill(dMult+0.5,three5n3n2n,dMult*(dMult-1.)*(dMult-2.));
   fIntFlowCorrelationsAllVsMPro[37]->Fill(dMult+0.5,three5n4n1n,dMult*(dMult-1.)*(dMult-2.));
   fIntFlowCorrelationsAllVsMPro[38]->Fill(dMult+0.5,three6n3n3n,dMult*(dMult-1.)*(dMult-2.));
   fIntFlowCorrelationsAllVsMPro[39]->Fill(dMult+0.5,three6n4n2n,dMult*(dMult-1.)*(dMult-2.));
   fIntFlowCorrelationsAllVsMPro[40]->Fill(dMult+0.5,three6n5n1n,dMult*(dMult-1.)*(dMult-2.));
  }       
 } // end of if(dMult>2)
 
 // 4-particle:
 Double_t four6n3n2n1n = 0.; // <cos(n(6*phi1-3*phi2-2*phi3-1*phi4)>
 Double_t four3n2n3n2n = 0.; // <cos(n(3*phi1+2*phi2-3*phi3-2*phi4)>
 Double_t four4n1n3n2n = 0.; // <cos(n(4*phi1+1*phi2-3*phi3-2*phi4)> 
 Double_t four3n3n3n3n = 0.; // <cos(3n(phi1+phi2-phi3-phi4))> 
 //Double_t four4n2n3n3n = 0.; // <cos(n(4*phi1+2*phi2-3*phi3-3*phi4)> // I already have this one above
 Double_t four5n1n3n3n = 0.; // <cos(n(5*phi1+1*phi2-3*phi3-3*phi4)>
 Double_t four4n2n4n2n = 0.; // <cos(n(4*phi1+2*phi2-4*phi3-2*phi4)> 
 Double_t four5n1n4n2n = 0.; // <cos(n(5*phi1+1*phi2-4*phi3-2*phi4)> 
 Double_t four5n3n1n1n = 0.; // <cos(n(5*phi1-3*phi2-1*phi3-1*phi4)> 
 Double_t four5n2n2n1n = 0.; // <cos(n(5*phi1-2*phi2-2*phi3-1*phi4)>
 Double_t four5n1n5n1n = 0.; // <cos(n(5*phi1+1*phi2-5*phi3-1*phi4)>
 Double_t four6n4n1n1n = 0.; // <cos(n(6*phi1-4*phi2-1*phi3-1*phi4)>
 Double_t four6n2n2n2n = 0.; // <cos(n(6*phi1-2*phi2-2*phi3-2*phi4)>
 if(dMult>3)
 {
  four6n3n2n1n = (reQ6nQ3nstarQ2nstarQ1nstar-reQ6nQ4nstarQ2nstar-reQ6nQ3nstarQ3nstar-reQ6nQ5nstarQ1nstar
               - reQ5nQ3nstarQ2nstar-reQ4nQ3nstarQ1nstar-reQ3nQ2nstarQ1nstar
               + 2.*(pow(dReQ6n,2.)+pow(dImQ6n,2.))+pow(dReQ5n,2.)+pow(dImQ5n,2.)
               + pow(dReQ4n,2.)+pow(dImQ4n,2.)+3.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
               + 2.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))+2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-6.*dMult)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  four3n2n3n2n = (dQ3nQ2nQ3nstarQ2nstar-2.*reQ5nQ3nstarQ2nstar-2.*reQ3nQ2nstarQ1nstar
               + pow(dReQ5n,2.)+pow(dImQ5n,2.)+pow(dReQ1n,2.)+pow(dImQ1n,2.)               
               -(dMult-4.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.)+pow(dReQ2n,2.)+pow(dImQ2n,2.))
               + dMult*(dMult-6.))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));     
  four4n1n3n2n = (reQ4nQ1nQ3nstarQ2nstar-reQ5nQ3nstarQ2nstar-reQ5nQ4nstarQ1nstar-reQ4nQ3nstarQ1nstar
               - reQ4nQ2nstarQ2nstar-reQ3nQ2nstarQ1nstar-reQ2nQ1nstarQ1nstar
               + pow(dReQ5n,2.)+pow(dImQ5n,2.)+2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
               + 2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+3.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
               + 3.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-6.*dMult) 
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));                             
  four3n3n3n3n = (2.*dMult*(dMult-3.)+pow((pow(dReQ3n,2.)+pow(dImQ3n,2.)),2.)-4.*(dMult-2.)*(pow(dReQ3n,2.)
               + pow(dImQ3n,2.))-2.*reQ6nQ3nstarQ3nstar+(pow(dReQ6n,2.)+pow(dImQ6n,2.)))
               / (dMult*(dMult-1)*(dMult-2.)*(dMult-3.));   
  //four4n2n3n3n = ; // I already have this one above
  four5n1n3n3n = (reQ5nQ1nQ3nstarQ3nstar-reQ6nQ5nstarQ1nstar-reQ6nQ3nstarQ3nstar-2.*reQ5nQ3nstarQ2nstar
               - 2.*reQ3nQ2nstarQ1nstar+pow(dReQ6n,2.)+pow(dImQ6n,2.)+2.*(pow(dReQ5n,2.)+pow(dImQ5n,2.))
               + 4.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+2.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
               + 2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-6.*dMult)                                  
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));   
  four4n2n4n2n = (dQ4nQ2nQ4nstarQ2nstar-2.*reQ6nQ4nstarQ2nstar-2.*reQ4nQ2nstarQ2nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - ((dMult-5.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
               + (dMult-4.)*(pow(dReQ4n,2.)+pow(dImQ4n,2.))-(pow(dReQ6n,2.)+pow(dImQ6n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               + (dMult-6.)/((dMult-1.)*(dMult-2.)*(dMult-3.)); 
  four5n1n4n2n = (reQ5nQ1nQ4nstarQ2nstar-reQ6nQ5nstarQ1nstar-reQ6nQ4nstarQ2nstar-reQ5nQ4nstarQ1nstar
               - reQ5nQ3nstarQ2nstar-reQ4nQ3nstarQ1nstar-reQ2nQ1nstarQ1nstar+pow(dReQ6n,2.)+pow(dImQ6n,2.)
               + 2.*(pow(dReQ5n,2.)+pow(dImQ5n,2.))+2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
               + pow(dReQ3n,2.)+pow(dImQ3n,2.)+2.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
               + 3.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-6.*dMult)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));     
  four5n3n1n1n = (reQ5nQ3nstarQ1nstarQ1nstar-2.*reQ5nQ4nstarQ1nstar-reQ5nQ3nstarQ2nstar-2.*reQ4nQ3nstarQ1nstar
               - reQ2nQ1nstarQ1nstar+2.*(pow(dReQ5n,2.)+pow(dImQ5n,2.))+2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
               + 2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+pow(dReQ2n,2.)+pow(dImQ2n,2.)
               + 4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-6.*dMult) 
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));                    
  four5n2n2n1n = (reQ5nQ2nstarQ2nstarQ1nstar-reQ5nQ4nstarQ1nstar-2.*reQ5nQ3nstarQ2nstar-reQ4nQ2nstarQ2nstar
               - 2.*reQ3nQ2nstarQ1nstar+2.*(pow(dReQ5n,2.)+pow(dImQ5n,2.))+pow(dReQ4n,2.)+pow(dImQ4n,2.)
               + 2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+4.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
               + 2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-6.*dMult)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));   
  four5n1n5n1n = (dQ5nQ1nQ5nstarQ1nstar-2.*reQ6nQ5nstarQ1nstar-2.*reQ5nQ4nstarQ1nstar
               + pow(dReQ6n,2.)+pow(dImQ6n,2.)-(dMult-4.)*(pow(dReQ5n,2.)+pow(dImQ5n,2.))
               + pow(dReQ4n,2.)+pow(dImQ4n,2.)-(dMult-4.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))+dMult*(dMult-6.))  
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));                  
  
  four6n4n1n1n = (reQ6nQ4nstarQ1nstarQ1nstar
               - dMult*(dMult-1.)*(dMult-2.)*(three2n1n1n+2.*three5n4n1n+2.*three6n5n1n+three6n4n2n)
               - dMult*(dMult-1.)*(2.*two1n1n+1.*two4n4n+1.*two6n6n+1.*two2n2n+2.*two5n5n)
               - 1.*dMult)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
 
  four6n2n2n2n = (reQ6nQ2nstarQ2nstarQ2nstar-3.*reQ6nQ4nstarQ2nstar-3.*reQ4nQ2nstarQ2nstar
               + 2.*(pow(dReQ6n,2.)+pow(dImQ6n,2.))+3.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
               + 6.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-6.*dMult)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  // Average 4-particle correlations for all events:      
  fIntFlowCorrelationsAllPro->Fill(41.5,four6n3n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(42.5,four3n2n3n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(43.5,four4n1n3n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(44.5,four3n3n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  //fIntFlowCorrelationsAllPro->Fill(45.5,four4n2n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)); // I already have this one above
  fIntFlowCorrelationsAllPro->Fill(46.5,four5n1n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(47.5,four4n2n4n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(48.5,four5n1n4n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(49.5,four5n3n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(50.5,four5n2n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(51.5,four5n1n5n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(58.5,four6n4n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(59.5,four6n2n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[41]->Fill(dMult+0.5,four6n3n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[42]->Fill(dMult+0.5,four3n2n3n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[43]->Fill(dMult+0.5,four4n1n3n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[44]->Fill(dMult+0.5,four3n3n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   //fIntFlowCorrelationsAllVsMPro[45]->Fill(dMult+0.5,four4n2n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[46]->Fill(dMult+0.5,four5n1n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[47]->Fill(dMult+0.5,four4n2n4n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[48]->Fill(dMult+0.5,four5n1n4n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[49]->Fill(dMult+0.5,four5n3n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[50]->Fill(dMult+0.5,four5n2n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[51]->Fill(dMult+0.5,four5n1n5n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[58]->Fill(dMult+0.5,four6n4n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
   fIntFlowCorrelationsAllVsMPro[59]->Fill(dMult+0.5,four6n2n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  }       
 } // end of if(dMult>3)

 // 5-particle:
 Double_t five3n3n3n2n1n = 0.; // <cos(n(3*phi1+3*phi2-3*phi3-2*phi4-1*phi5)>
 Double_t five4n2n3n2n1n = 0.; // <cos(n(4*phi1+2*phi2-3*phi3-2*phi4-1*phi5)>
 Double_t five3n2n3n1n1n = 0.; // <cos(n(3*phi1+2*phi2-3*phi3-1*phi4-1*phi5)>
 Double_t five3n2n2n2n1n = 0.; // <cos(n(3*phi1+2*phi2-2*phi3-2*phi4-1*phi5)>
 Double_t five5n1n3n2n1n = 0.; // <cos(n(5*phi1+1*phi2-3*phi3-2*phi4-1*phi5)>
 Double_t five6n2n2n1n1n = 0.; // <cos(n(6*phi1-2*phi2-2*phi3-1*phi4-1*phi5)>
 Double_t five4n1n1n3n3n = 0.; // <cos(n(4*phi1+1*phi2+1*phi3-3*phi4-3*phi5)>
 if(dMult>4)
 { 
  five3n3n3n2n1n = (reQ3nQ3nQ3nstarQ2nstarQ1nstar-reQ6nQ3nstarQ2nstarQ1nstar-reQ5nQ1nQ3nstarQ3nstar-reQ4nQ2nQ3nstarQ3nstar
                 + reQ6nQ5nstarQ1nstar+reQ6nQ4nstarQ2nstar+3.*reQ6nQ3nstarQ3nstar+4.*reQ5nQ3nstarQ2nstar+4.*reQ4nQ3nstarQ1nstar
                 - 2.*(dMult-6.)*reQ3nQ2nstarQ1nstar-2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.)) 
                 - 2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - 2.*(pow(dReQ6n,2.)+pow(dImQ6n,2.))-2.*(pow(dReQ5n,2.)+pow(dImQ5n,2.))
                 - 2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+2.*(3.*dMult-10.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 - pow((pow(dReQ3n,2.)+pow(dImQ3n,2.)),2.)+2.*(dMult-5.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 + 2.*(dMult-5.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-4.*dMult*(dMult-6.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));                                  
  five4n2n3n2n1n = (reQ4nQ2nQ3nstarQ2nstarQ1nstar-reQ6nQ3nstarQ2nstarQ1nstar-reQ5nQ1nQ4nstarQ2nstar
                 - reQ4nQ2nQ3nstarQ3nstar-reQ4nQ1nQ3nstarQ2nstar-reQ4nQ2nstarQ1nstarQ1nstar
                 - reQ3nQ1nQ2nstarQ2nstar-(pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 + 3.*reQ6nQ4nstarQ2nstar+reQ6nQ5nstarQ1nstar+reQ6nQ3nstarQ3nstar+reQ5nQ4nstarQ1nstar
                 + 3.*reQ5nQ3nstarQ2nstar-(dMult-7.)*reQ4nQ3nstarQ1nstar+3.*reQ4nQ2nstarQ2nstar+7.*reQ3nQ2nstarQ1nstar
                 + 4.*reQ2nQ1nstarQ1nstar-(pow(dReQ4n,2.)+pow(dImQ4n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 - (pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - 2.*(pow(dReQ6n,2.)+pow(dImQ6n,2.))-2.*(pow(dReQ5n,2.)+pow(dImQ5n,2.))
                 + (dMult-8.)*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+(dMult-10.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 + 2.*(dMult-7.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))+(dMult-12.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - 2.*dMult*(dMult-12.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));   
  five3n2n3n1n1n = (reQ3nQ2nQ3nstarQ1nstarQ1nstar-reQ5nQ3nstarQ1nstarQ1nstar-2.*reQ4nQ1nQ3nstarQ2nstar-reQ3nQ1nstarQ1nstarQ1nstar
                 - 2.*reQ3nQ1nQ2nstarQ2nstar+2.*reQ5nQ4nstarQ1nstar+3.*reQ5nQ3nstarQ2nstar+6.*reQ4nQ3nstarQ1nstar
                 + 2.*reQ4nQ2nstarQ2nstar+9.*reQ3nQ2nstarQ1nstar-(dMult-8.)*reQ2nQ1nstarQ1nstar
                 - (pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.)) 
                 - 2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.)) 
                 - 2.*(pow(dReQ5n,2.)+pow(dImQ5n,2.))-4.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
                 + 2.*(dMult-6.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+(dMult-12.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 + 2.*(dMult-9.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-2.*dMult*(dMult-12.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  five3n2n2n2n1n = (reQ3nQ2nQ2nstarQ2nstarQ1nstar-reQ5nQ2nstarQ2nstarQ1nstar-reQ4nQ1nQ3nstarQ2nstar-reQ3nQ1nQ2nstarQ2nstar
                 - 2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))+reQ5nQ4nstarQ1nstar
                 + 4.*reQ5nQ3nstarQ2nstar+reQ4nQ3nstarQ1nstar+3.*reQ4nQ2nstarQ2nstar-2.*(dMult-6.)*reQ3nQ2nstarQ1nstar
                 + 4.*reQ2nQ1nstarQ1nstar-2.*(pow(dReQ5n,2.)+pow(dImQ5n,2.))
                 - 2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+2.*(dMult-5.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 - 2.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - pow((pow(dReQ2n,2.)+pow(dImQ2n,2.)),2.)+2.*(3.*dMult-10.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 + 2.*(dMult-6.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-4.*dMult*(dMult-6.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)); 
  five5n1n3n2n1n = (reQ5nQ1nQ3nstarQ2nstarQ1nstar-reQ6nQ3nstarQ2nstarQ1nstar-reQ5nQ1nQ4nstarQ2nstar-reQ5nQ1nQ3nstarQ3nstar
                 - reQ4nQ1nQ3nstarQ2nstar-reQ5nQ3nstarQ1nstarQ1nstar-reQ5nQ2nstarQ2nstarQ1nstar 
                 + 3.*reQ6nQ5nstarQ1nstar+reQ6nQ4nstarQ2nstar+reQ6nQ3nstarQ3nstar+4.*reQ5nQ4nstarQ1nstar
                 - (dMult-7.)*reQ5nQ3nstarQ2nstar+4.*reQ4nQ3nstarQ1nstar+reQ4nQ2nstarQ2nstar+6.*reQ3nQ2nstarQ1nstar
                 + 3.*reQ2nQ1nstarQ1nstar-(pow(dReQ5n,2.)+pow(dImQ5n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - (pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - (pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - 2.*(pow(dReQ6n,2.)+pow(dImQ6n,2.))+(dMult-8.)*(pow(dReQ5n,2.)+pow(dImQ5n,2.))  
                 - 4.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+(dMult-10.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 + (dMult-10.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))+2.*(dMult-7.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - 2.*dMult*(dMult-12.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));      
  
  // five6n2n2n1n1n = ;
  // five4n1n1n3n3n = ; 
  
  // Average 5-particle correlations for all events:      
  fIntFlowCorrelationsAllPro->Fill(52.5,five3n3n3n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(53.5,five4n2n3n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(54.5,five3n2n3n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(55.5,five3n2n2n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(56.5,five5n1n3n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(60.5,five6n2n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(61.5,five4n1n1n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[52]->Fill(dMult+0.5,five3n3n3n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
   fIntFlowCorrelationsAllVsMPro[53]->Fill(dMult+0.5,five4n2n3n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
   fIntFlowCorrelationsAllVsMPro[54]->Fill(dMult+0.5,five3n2n3n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
   fIntFlowCorrelationsAllVsMPro[55]->Fill(dMult+0.5,five3n2n2n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
   fIntFlowCorrelationsAllVsMPro[56]->Fill(dMult+0.5,five5n1n3n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
   fIntFlowCorrelationsAllVsMPro[60]->Fill(dMult+0.5,five6n2n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
   fIntFlowCorrelationsAllVsMPro[61]->Fill(dMult+0.5,five4n1n1n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  }         
 } // end of if(dMult>4)

 // 6-particle:
 Double_t six3n2n1n3n2n1n = 0.; // <cos(n(3*phi1+2*phi2+1*phi3-3*phi4-2*phi5-1*phi6)>
 Double_t six3n3n2n2n1n1n = 0.; // <cos(n(3*phi1+3*phi2-2*phi3-2*phi4-1*phi5-1*phi6)>   
 if(dMult>5.)
 { 
  six3n2n1n3n2n1n = (dQ3nQ2nQ1nQ3nstarQ2nstarQ1nstar-2.*reQ3nQ3nQ3nstarQ2nstarQ1nstar
                  - 2.*reQ3nQ2nQ2nstarQ2nstarQ1nstar-2.*reQ3nQ1nQ2nstarQ1nstarQ1nstar
                  - 2.*reQ3nQ2nQ3nstarQ1nstarQ1nstar-2.*reQ4nQ2nQ3nstarQ2nstarQ1nstar
                  - 2.*reQ5nQ1nQ3nstarQ2nstarQ1nstar+4.*reQ6nQ3nstarQ2nstarQ1nstar
                  + 2.*reQ5nQ1nQ4nstarQ2nstar+2.*reQ5nQ1nQ3nstarQ3nstar
                  + 2.*reQ4nQ2nQ3nstarQ3nstar+6.*reQ4nQ1nQ3nstarQ2nstar
                  + 2.*reQ5nQ3nstarQ1nstarQ1nstar+2.*reQ5nQ2nstarQ2nstarQ1nstar
                  + 6.*reQ3nQ1nQ2nstarQ2nstar+2.*reQ4nQ2nstarQ1nstarQ1nstar
                  - 4.*reQ6nQ5nstarQ1nstar-4.*reQ6nQ4nstarQ2nstar-6.*reQ5nQ4nstarQ1nstar
                  - 4.*reQ6nQ3nstarQ3nstar+2.*(dMult-11.)*reQ5nQ3nstarQ2nstar
                  + 2.*(dMult-13.)*reQ4nQ3nstarQ1nstar-8.*reQ4nQ2nstarQ2nstar
                  + 2.*(5.*dMult-32.)*reQ3nQ2nstarQ1nstar+2.*reQ3nQ1nstarQ1nstarQ1nstar
                  + 2.*(dMult-13.)*reQ2nQ1nstarQ1nstar
                  - (dMult-10.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                  + (pow(dReQ5n,2.)+pow(dImQ5n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                  + (pow(dReQ4n,2.)+pow(dImQ4n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                  - (dMult-11.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                  - (dMult-10.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                  + 4.*(pow(dReQ6n,2.)+pow(dImQ6n,2.))-(dMult-12.)*(pow(dReQ5n,2.)+pow(dImQ5n,2.))
                  - (dMult-16.)*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+pow((pow(dReQ3n,2.)+pow(dImQ3n,2.)),2.)
                  + (dMult*dMult-19.*dMult+68.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                  + (dMult*dMult-19.*dMult+72.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                  + pow((pow(dReQ2n,2.)+pow(dImQ2n,2.)),2.)+pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                  + (dMult*dMult-20.*dMult+80.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                  - dMult*(dMult-12.)*(dMult-10.))
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));   

  // six3n3n2n2n1n1n = ;

  // Average 6-particle correlations for all events:      
  fIntFlowCorrelationsAllPro->Fill(57.5,six3n2n1n3n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  fIntFlowCorrelationsAllPro->Fill(62.5,six3n3n2n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  if(fCalculateAllCorrelationsVsM)
  {
   fIntFlowCorrelationsAllVsMPro[57]->Fill(dMult+0.5,six3n2n1n3n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
   fIntFlowCorrelationsAllVsMPro[62]->Fill(dMult+0.5,six3n3n2n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  }          
 } // end of if(dMult>5.)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrelations()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::StorePhiDistributionForOneEvent(AliFlowEventSimple *anEvent)
{
 // Store phi distribution for one event to illustrate flow.
 
 if(fPhiDistributionForOneEvent->GetEntries()>0){return;} // store only phi distribution for one event
 
 Double_t vMin = fPhiDistributionForOneEventSettings[0]; 
 Double_t vMax = fPhiDistributionForOneEventSettings[1]; 
 Double_t refMultMin = fPhiDistributionForOneEventSettings[2]; 
 Double_t refMultMax = fPhiDistributionForOneEventSettings[3]; 
 
 Double_t vEBE = 0.;
 Double_t cumulant4thEBE = fIntFlowCorrelationsEBE->GetBinContent(2)-2.*pow(fIntFlowCorrelationsEBE->GetBinContent(1),2.);
 if(cumulant4thEBE<0.)
 {
  vEBE = pow(-1.*cumulant4thEBE,0.25);
  if((vEBE>vMin && vEBE<vMax) && (fReferenceMultiplicityEBE>refMultMin && fReferenceMultiplicityEBE<refMultMax))
  {
   fPhiDistributionForOneEvent->SetTitle(Form("v_{%i} = %f",fHarmonic,vEBE));
   for(Int_t p=0;p<anEvent->NumberOfTracks();p++)
   {
    if(anEvent->GetTrack(p)->InRPSelection())
    {
     fPhiDistributionForOneEvent->Fill(anEvent->GetTrack(p)->Phi());
    }
   } // end of for(Int_t p=0;p<anEvent->NumberOfTracks();p++)
  } else
    {
     fPhiDistributionForOneEvent->SetTitle(Form("v_{%i} = %f, out of specified boundaries",fHarmonic,vEBE));  
    } 
   
 } // end of if(cumulant4thEBE<0.)
 
} // end of void AliFlowAnalysisWithQCumulants::StorePhiDistributionForOneEvent(AliFlowEventSimple *anEvent)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateIntFlowProductOfCorrelations()
{
 // Calculate averages of products of correlations for integrated flow.
 
 // multiplicity:
 Double_t dMult = (*fSpk)(0,0);
 
 Int_t counter = 0;
 
 for(Int_t ci1=1;ci1<4;ci1++)
 {
  for(Int_t ci2=ci1+1;ci2<=4;ci2++)
  {
   fIntFlowProductOfCorrelationsPro->Fill(0.5+counter,
                                          fIntFlowCorrelationsEBE->GetBinContent(ci1)*
                                          fIntFlowCorrelationsEBE->GetBinContent(ci2),
                                          fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci1)*
                                          fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci2));
   // products versus multiplicity:  // [0=<<2><4>>,1=<<2><6>>,2=<<2><8>>,3=<<4><6>>,4=<<4><8>>,5=<<6><8>>]
   if(fCalculateCumulantsVsM)
   {
    fIntFlowProductOfCorrelationsVsMPro[counter]->Fill(dMult+0.5, // to be improved: dMult => sum of weights ?
                                                       fIntFlowCorrelationsEBE->GetBinContent(ci1)*
                                                       fIntFlowCorrelationsEBE->GetBinContent(ci2),
                                                       fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci1)*
                                                       fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci2));
   } // end of if(fCalculateCumulantsVsM)
   counter++;                                                                                                                        
  }
 }
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowProductOfCorrelations()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateIntFlowProductOfCorrectionTermsForNUA()
{
 // Calculate averages of products of correction terms for NUA.
 
 // a) Binning of fIntFlowProductOfCorrectionTermsForNUAPro is organized as follows:
 //     1st bin: <<2><cos(phi)>> 
 //     2nd bin: <<2><sin(phi)>>
 //     3rd bin: <<cos(phi)><sin(phi)>>
 //     4th bin: <<2><cos(phi1+phi2)>> 
 //     5th bin: <<2><sin(phi1+phi2)>>
 //     6th bin: <<2><cos(phi1-phi2-phi3)>> 
 //     7th bin: <<2><sin(phi1-phi2-phi3)>>
 //     8th bin: <<4><cos(phi1)>>
 //     9th bin: <<4><sin(phi1)>>
 //    10th bin: <<4><cos(phi1+phi2)>>
 //    11th bin: <<4><sin(phi1+phi2)>>
 //    12th bin: <<4><cos(phi1-phi2-phi3)>>
 //    13th bin: <<4><sin(phi1-phi2-phi3)>>
 //    14th bin: <<cos(phi1)><cos(phi1+phi2)>>
 //    15th bin: <<cos(phi1)><sin(phi1+phi2)>> 
 //    16th bin: <<cos(phi1)><cos(phi1-phi2-phi3)>>
 //    17th bin: <<cos(phi1)><sin(phi1-phi2-phi3)>> 
 //    18th bin: <<sin(phi1)><cos(phi1+phi2)>>
 //    19th bin: <<sin(phi1)><sin(phi1+phi2)>> 
 //    20th bin: <<sin(phi1)><cos(phi1-phi2-phi3)>>
 //    21st bin: <<sin(phi1)><sin(phi1-phi2-phi3)>>
 //    22nd bin: <<cos(phi1+phi2)><sin(phi1+phi2)>>
 //    23rd bin: <<cos(phi1+phi2)><cos(phi1-phi2-phi3)>>
 //    24th bin: <<cos(phi1+phi2)><sin(phi1-phi2-phi3)>>
 //    25th bin: <<sin(phi1+phi2)><cos(phi1-phi2-phi3)>>
 //    26th bin: <<sin(phi1+phi2)><sin(phi1-phi2-phi3)>>
 //    27th bin: <<cos(phi1-phi2-phi3)><sin(phi1-phi2-phi3)>>
 
 // <<2><cos(phi)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(0.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(1),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1));
 // <<2><sin(phi)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(1.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(1),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1));
 // <<cos(phi)><sin(phi)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(2.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(1),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1));
 // <<2><cos(phi1+phi2)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(3.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(2),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)); 
 // <<2><sin(phi1+phi2)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(4.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(2),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2)); 
 // <<2><cos(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(5.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(3),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3)); 
 // <<2><sin(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(6.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(3),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3)); 
 // <<4><cos(phi1)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(7.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(1),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1));
 // <<4><sin(phi1)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(8.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(1),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1));
 // <<4><cos(phi1+phi2)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(9.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(2),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)); 
 // <<4><sin(phi1+phi2)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(10.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(2),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2));
 // <<4><cos(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(11.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(3),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3)); 
 // <<4><sin(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(12.5,
                                                 fIntFlowCorrelationsEBE->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(3),
                                                 fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3));   
 // <<cos(phi1)><cos(phi1+phi2)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(13.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(2),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)); 
 // <<cos(phi1)><sin(phi1+phi2)>>: 
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(14.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(2),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2)); 
 // <<cos(phi1)><cos(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(15.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(3),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3)); 
 // <<cos(phi1)><sin(phi1-phi2-phi3)>>: 
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(16.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(3),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3));  
 // <<sin(phi1)><cos(phi1+phi2)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(17.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(2),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2));  
 // <<sin(phi1)><sin(phi1+phi2)>>: 
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(18.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(2),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2));  
 // <<sin(phi1)><cos(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(19.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(3),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3)); 
 // <<sin(phi1)><sin(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(20.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(1)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(3),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3)); 
 // <<cos(phi1+phi2)><sin(phi1+phi2)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(21.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(2),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2)); 
 // <<cos(phi1+phi2)><cos(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(22.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(3),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3));   
 // <<cos(phi1+phi2)><sin(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(23.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(3),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3));    
 // <<sin(phi1+phi2)><cos(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(24.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(3),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3));    
 // <<sin(phi1+phi2)><sin(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(25.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(2)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(3),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3));    
 // <<cos(phi1-phi2-phi3)><sin(phi1-phi2-phi3)>>:
 fIntFlowProductOfCorrectionTermsForNUAPro->Fill(26.5,
                                                 fIntFlowCorrectionTermsForNUAEBE[1]->GetBinContent(3)*fIntFlowCorrectionTermsForNUAEBE[0]->GetBinContent(3),
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3)
                                                 *fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3));    

} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowProductOfCorrectionTermsForNUA()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateCovariancesIntFlow()
{
 // a) Calculate unbiased estimators Cov(<2>,<4>), Cov(<2>,<6>), Cov(<2>,<8>), Cov(<4>,<6>), Cov(<4>,<8>) and Cov(<6>,<8>)
 //    for covariances V_(<2>,<4>), V_(<2>,<6>), V_(<2>,<8>), V_(<4>,<6>), V_(<4>,<8>) and V_(<6>,<8>).
 // b) Store in histogram fIntFlowCovariances for instance the following: 
 //
 //             Cov(<2>,<4>) * (sum_{i=1}^{N} w_{<2>}_i w_{<4>}_i )/[(sum_{i=1}^{N} w_{<2>}_i) * (sum_{j=1}^{N} w_{<4>}_j)]
 // 
 //    where N is the number of events, w_{<2>} is event weight for <2> and w_{<4>} is event weight for <4>.
 // c) Binning of fIntFlowCovariances is organized as follows:
 // 
 //     1st bin: Cov(<2>,<4>) * (sum_{i=1}^{N} w_{<2>}_i w_{<4>}_i )/[(sum_{i=1}^{N} w_{<2>}_i) * (sum_{j=1}^{N} w_{<4>}_j)] 
 //     2nd bin: Cov(<2>,<6>) * (sum_{i=1}^{N} w_{<2>}_i w_{<6>}_i )/[(sum_{i=1}^{N} w_{<2>}_i) * (sum_{j=1}^{N} w_{<6>}_j)]
 //     3rd bin: Cov(<2>,<8>) * (sum_{i=1}^{N} w_{<2>}_i w_{<8>}_i )/[(sum_{i=1}^{N} w_{<2>}_i) * (sum_{j=1}^{N} w_{<8>}_j)]
 //     4th bin: Cov(<4>,<6>) * (sum_{i=1}^{N} w_{<4>}_i w_{<6>}_i )/[(sum_{i=1}^{N} w_{<4>}_i) * (sum_{j=1}^{N} w_{<6>}_j)]
 //     5th bin: Cov(<4>,<8>) * (sum_{i=1}^{N} w_{<4>}_i w_{<8>}_i )/[(sum_{i=1}^{N} w_{<4>}_i) * (sum_{j=1}^{N} w_{<8>}_j)]
 //     6th bin: Cov(<6>,<8>) * (sum_{i=1}^{N} w_{<6>}_i w_{<8>}_i )/[(sum_{i=1}^{N} w_{<6>}_i) * (sum_{j=1}^{N} w_{<8>}_j)]
 //
    
 // Average 2-, 4-, 6- and 8-particle correlations for all events:
 Double_t correlation[4] = {0.};
 for(Int_t ci=0;ci<4;ci++)
 {
  correlation[ci] = fIntFlowCorrelationsPro->GetBinContent(ci+1);
 } 
 // Average products of 2-, 4-, 6- and 8-particle correlations: 
 Double_t productOfCorrelations[4][4] = {{0.}};
 Int_t productOfCorrelationsLabel = 1;
 // Denominators in the expressions for the unbiased estimator for covariance:
 Double_t denominator[4][4] = {{0.}};
 Int_t sumOfProductOfEventWeightsLabel1 = 1;
 // Weight dependent prefactor which multiply unbiased estimators for covariances:
 Double_t wPrefactor[4][4] = {{0.}}; 
 Int_t sumOfProductOfEventWeightsLabel2 = 1;
 for(Int_t c1=0;c1<4;c1++)
 {
  for(Int_t c2=c1+1;c2<4;c2++)
  {
   productOfCorrelations[c1][c2] = fIntFlowProductOfCorrelationsPro->GetBinContent(productOfCorrelationsLabel);
   if(TMath::Abs(fIntFlowSumOfEventWeights[0]->GetBinContent(c1+1)) > 1.e-44 && TMath::Abs(fIntFlowSumOfEventWeights[0]->GetBinContent(c2+1)) > 1.e-44)
   {
    denominator[c1][c2] = 1.-(fIntFlowSumOfProductOfEventWeights->GetBinContent(sumOfProductOfEventWeightsLabel1))
                        / (fIntFlowSumOfEventWeights[0]->GetBinContent(c1+1) 
                        * fIntFlowSumOfEventWeights[0]->GetBinContent(c2+1));                              
    wPrefactor[c1][c2] = fIntFlowSumOfProductOfEventWeights->GetBinContent(sumOfProductOfEventWeightsLabel2)
                       / (fIntFlowSumOfEventWeights[0]->GetBinContent(c1+1)
                       * fIntFlowSumOfEventWeights[0]->GetBinContent(c2+1));                                                       
   }
   productOfCorrelationsLabel++; // to be improved - do I need here all 3 counters?
   sumOfProductOfEventWeightsLabel1++;
   sumOfProductOfEventWeightsLabel2++;  
  } // end of for(Int_t c2=c1+1;c2<4;c2++)
 } // end of for(Int_t c1=0;c1<4;c1++)
 
 Int_t covarianceLabel = 1;
 for(Int_t c1=0;c1<4;c1++)
 {
  for(Int_t c2=c1+1;c2<4;c2++)
  {
   if(TMath::Abs(denominator[c1][c2]) > 1.e-44)
   {
    // Covariances:
    Double_t cov = (productOfCorrelations[c1][c2]-correlation[c1]*correlation[c2])/denominator[c1][c2]; 
    // Covariances multiplied with weight dependent prefactor:
    Double_t wCov = cov * wPrefactor[c1][c2];
    fIntFlowCovariances->SetBinContent(covarianceLabel,wCov);
   }
   covarianceLabel++;
  } // end of for(Int_t c2=c1+1;c2<4;c2++) 
 } // end of for(Int_t c1=0;c1<4;c1++)
 
 // Versus multiplicity: 
 if(!fCalculateCumulantsVsM){return;}
 Int_t nBins = fIntFlowCorrelationsVsMPro[0]->GetNbinsX(); // to be improved (hardwired 0) 
 for(Int_t b=1;b<=nBins;b++)
 {
  // Average 2-, 4-, 6- and 8-particle correlations for all events:
  Double_t correlationVsM[4] = {0.};
  for(Int_t ci=0;ci<4;ci++)
  {
   correlationVsM[ci] = fIntFlowCorrelationsVsMPro[ci]->GetBinContent(b);
  } // end of for(Int_t ci=0;ci<4;ci++)
  // Average products of 2-, 4-, 6- and 8-particle correlations: 
  Double_t productOfCorrelationsVsM[4][4] = {{0.}};
  Int_t productOfCorrelationsLabelVsM = 1;
  // Denominators in the expressions for the unbiased estimator for covariance:
  Double_t denominatorVsM[4][4] = {{0.}};
  Int_t sumOfProductOfEventWeightsLabel1VsM = 1;
  // Weight dependent prefactor which multiply unbiased estimators for covariances:
  Double_t wPrefactorVsM[4][4] = {{0.}}; 
  Int_t sumOfProductOfEventWeightsLabel2VsM = 1;
  for(Int_t c1=0;c1<4;c1++)
  {
   for(Int_t c2=c1+1;c2<4;c2++)
   {
    productOfCorrelationsVsM[c1][c2] = fIntFlowProductOfCorrelationsVsMPro[productOfCorrelationsLabelVsM-1]->GetBinContent(b);
    if(TMath::Abs(fIntFlowSumOfEventWeightsVsM[c1][0]->GetBinContent(b)) > 1.e-44 && TMath::Abs(fIntFlowSumOfEventWeightsVsM[c2][0]->GetBinContent(b)) > 1.e-44)
    {
     denominatorVsM[c1][c2] = 1.-(fIntFlowSumOfProductOfEventWeightsVsM[sumOfProductOfEventWeightsLabel1VsM-1]->GetBinContent(b))
                            / (fIntFlowSumOfEventWeightsVsM[c1][0]->GetBinContent(b) 
                            * fIntFlowSumOfEventWeightsVsM[c2][0]->GetBinContent(b));                              
     wPrefactorVsM[c1][c2] = fIntFlowSumOfProductOfEventWeightsVsM[sumOfProductOfEventWeightsLabel2VsM-1]->GetBinContent(b)
                           / (fIntFlowSumOfEventWeightsVsM[c1][0]->GetBinContent(b)
                           * fIntFlowSumOfEventWeightsVsM[c2][0]->GetBinContent(b));                                                       
    }
    productOfCorrelationsLabelVsM++;
    sumOfProductOfEventWeightsLabel1VsM++;
    sumOfProductOfEventWeightsLabel2VsM++;  
   } // end of for(Int_t c1=0;c1<4;c1++) 
  } // end of for(Int_t c2=c1+1;c2<4;c2++)
 
  Int_t covarianceLabelVsM = 1;
  for(Int_t c1=0;c1<4;c1++)
  {
   for(Int_t c2=c1+1;c2<4;c2++)
   {
    if(TMath::Abs(denominatorVsM[c1][c2]) > 1.e-44)
    {
     // Covariances:
     Double_t covVsM = (productOfCorrelationsVsM[c1][c2]-correlationVsM[c1]*correlationVsM[c2])/denominatorVsM[c1][c2]; 
     // Covariances multiplied with weight dependent prefactor:
     Double_t wCovVsM = covVsM * wPrefactorVsM[c1][c2];
     fIntFlowCovariancesVsM[covarianceLabelVsM-1]->SetBinContent(b,wCovVsM);
    }
    covarianceLabelVsM++;
   } // end of for(Int_t c2=c1+1;c2<4;c2++)
  } // end of for(Int_t c1=0;c1<4;c1++)
 } // end of for(Int_t b=1;b<=nBins;b++)
  
} // end of AliFlowAnalysisWithQCumulants::CalculateCovariancesIntFlow()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateCovariancesNUAIntFlow()
{
 // a) Calculate unbiased estimators Cov(*,*) for true covariances V_(*,*) for NUA terms.
 // b) Store in histogram fIntFlowCovariancesNUA for instance the following: 
 //
 //             Cov(<2>,<cos(phi)>) * (sum_{i=1}^{N} w_{<2>}_i w_{<cos(phi)>}_i )/[(sum_{i=1}^{N} w_{<2>}_i) * (sum_{j=1}^{N} w_{<cos(phi)>}_j)]
 // 
 //    where N is the number of events, w_{<2>} is event weight for <2> and w_{<cos(phi)>} is event weight for <cos(phi)>.
 // c) Binning of fIntFlowCovariancesNUA is organized as follows:
 // 
 //     1st bin: Cov(<2>,<cos(phi)>) * (sum_{i=1}^{N} w_{<2>}_i w_{<cos(phi)>}_i )/[(sum_{i=1}^{N} w_{<2>}_i) * (sum_{j=1}^{N} w_{<cos(phi)>}_j)] 
 //     2nd bin: Cov(<2>,<sin(phi)>) * (sum_{i=1}^{N} w_{<2>}_i w_{<sin(phi)>}_i )/[(sum_{i=1}^{N} w_{<2>}_i) * (sum_{j=1}^{N} w_{<sin(phi)>}_j)]
 //     3rd bin: Cov(<cos(phi)>,<sin(phi)>) * (sum_{i=1}^{N} w_{<cos(phi)>}_i w_{<sin(phi)>}_i )/[(sum_{i=1}^{N} w_{<cos(phi)>}_i) * (sum_{j=1}^{N} w_{<sin(phi)>}_j)]
 // ...
      
 // Cov(<2>,<cos(phi)>):
 Double_t product1 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(1); // <<2><cos(phi)>> 
 Double_t term1st1 = fIntFlowCorrelationsPro->GetBinContent(1); // <<2>>
 Double_t term2nd1 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(1); // <<cos(phi)>>
 Double_t sumOfW1st1 = fIntFlowSumOfEventWeights[0]->GetBinContent(1); // W_{<2>}
 Double_t sumOfW2nd1 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(1); // W_{<cos(phi)>}
 Double_t sumOfWW1 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(1); // W_{<2>} * W_{<cos(phi)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator1 = product1 - term1st1*term2nd1; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator1 = 0.;
 if(TMath::Abs(sumOfW1st1*sumOfW2nd1)>0.)
 {
  denominator1 = 1.-sumOfWW1/(sumOfW1st1*sumOfW2nd1);
  if(TMath::Abs(denominator1)>0.)
  {
   // covariance:
   Double_t covariance1 = numerator1/denominator1;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor1 = sumOfWW1/(sumOfW1st1*sumOfW2nd1);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(1,wPrefactor1*covariance1);
  } // end of if(TMath::Abs(denominator)>0.)
 } // end of if(TMath::Abs(sumOfW1st1*sumOfW2nd1)>0.)
 
 // Cov(<2>,<sin(phi)>):
 Double_t product2 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(2); // <<2><sin(phi)>> 
 Double_t term1st2 = fIntFlowCorrelationsPro->GetBinContent(1); // <<2>>
 Double_t term2nd2 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(1); // <<sin(phi)>>
 Double_t sumOfW1st2 = fIntFlowSumOfEventWeights[0]->GetBinContent(1); // W_{<2>}
 Double_t sumOfW2nd2 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(1); // W_{<sin(phi)>}
 Double_t sumOfWW2 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(2); // W_{<2>} * W_{<sin(phi)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator2 = product2 - term1st2*term2nd2;
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator2 = 0.;
 if(TMath::Abs(sumOfW1st2*sumOfW2nd2)>0.)
 {  
  denominator2 = 1.-sumOfWW2/(sumOfW1st2*sumOfW2nd2);
  if(TMath::Abs(denominator2)>0.)
  {
   // covariance:
   Double_t covariance2 = numerator2/denominator2;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor2 = sumOfWW2/(sumOfW1st2*sumOfW2nd2);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(2,wPrefactor2*covariance2);
  } // end of if(TMath::Abs(denominator2)>0.)
 } // end of if(TMath::Abs(sumOfW1st2*sumOfW2nd2)>0.)
 
 // Cov(<cos(phi)>,<sin(phi)>):
 Double_t product3 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(3); // <<cos(phi)><sin(phi)>> 
 Double_t term1st3 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(1); // <<cos(phi)>>
 Double_t term2nd3 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(1); // <<sin(phi)>>
 Double_t sumOfW1st3 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(1); // W_{<cos(phi)>}
 Double_t sumOfW2nd3 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(1); // W_{<sin(phi)>}
 Double_t sumOfWW3 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(3); // W_{<cos(phi)>} * W_{<sin(phi)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator3 = product3 - term1st3*term2nd3; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator3 = 0;
 if(TMath::Abs(sumOfW1st3*sumOfW2nd3)>0.)
 { 
  denominator3 = 1.-sumOfWW3/(sumOfW1st3*sumOfW2nd3);
  if(TMath::Abs(denominator3)>0.)
  {
   // covariance:
   Double_t covariance3 = numerator3/denominator3;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor3 = sumOfWW3/(sumOfW1st3*sumOfW2nd3);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(3,wPrefactor3*covariance3);
  } // end of if(TMath::Abs(denominator3)>0.)
 } // end of if(TMath::Abs(sumOfW1st3*sumOfW2nd3)>0.)
 
 // Cov(<2>,<cos(phi1+phi2)>):
 Double_t product4 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(4); // <<2><cos(phi1+phi2)>> 
 Double_t term1st4 = fIntFlowCorrelationsPro->GetBinContent(1); // <<2>>
 Double_t term2nd4 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(2); // <<cos(phi1+phi2)>>
 Double_t sumOfW1st4 = fIntFlowSumOfEventWeights[0]->GetBinContent(1); // W_{<2>}
 Double_t sumOfW2nd4 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(2); // W_{<cos(phi1+phi2)>}
 Double_t sumOfWW4 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(4); // W_{<2>} * W_{<cos(phi1+phi2)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator4 = product4 - term1st4*term2nd4; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator4 = 0.;
 if(TMath::Abs(sumOfW1st4*sumOfW2nd4)>0.)
 { 
  denominator4 = 1.-sumOfWW4/(sumOfW1st4*sumOfW2nd4);
  if(TMath::Abs(denominator4)>0.)
  {  
   // covariance:
   Double_t covariance4 = numerator4/denominator4;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor4 = sumOfWW4/(sumOfW1st4*sumOfW2nd4);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(4,wPrefactor4*covariance4);
  } // end of if(TMath::Abs(denominator4)>0.)
 } // end of if(TMath::Abs(sumOfW1st4*sumOfW2nd4)>0.)
 
 // Cov(<2>,<sin(phi1+phi2)>):
 Double_t product5 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(5); // <<2><sin(phi1+phi2)>> 
 Double_t term1st5 = fIntFlowCorrelationsPro->GetBinContent(1); // <<2>>
 Double_t term2nd5 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(2); // <<sin(phi1+phi2)>>
 Double_t sumOfW1st5 = fIntFlowSumOfEventWeights[0]->GetBinContent(1); // W_{<2>}
 Double_t sumOfW2nd5 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(2); // W_{<sin(phi1+phi2)>}
 Double_t sumOfWW5 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(5); // W_{<2>} * W_{<sin(phi1+phi2)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator5 = product5 - term1st5*term2nd5; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator5 = 0.;
 if(TMath::Abs(sumOfW1st5*sumOfW2nd5)>0.)
 {  
  denominator5 = 1.-sumOfWW5/(sumOfW1st5*sumOfW2nd5);
  if(TMath::Abs(denominator5)>0.)
  {  
   // covariance:
   Double_t covariance5 = numerator5/denominator5;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor5 = sumOfWW5/(sumOfW1st5*sumOfW2nd5);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(5,wPrefactor5*covariance5);
  } // end of if(TMath::Abs(denominator5)>0.)
 } // end of if(TMath::Abs(sumOfW1st5*sumOfW2nd5)>0.)
 
 // Cov(<2>,<cos(phi1-phi2-phi3)>):
 Double_t product6 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(6); // <<2><cos(phi1-phi2-phi3)>> 
 Double_t term1st6 = fIntFlowCorrelationsPro->GetBinContent(1); // <<2>>
 Double_t term2nd6 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(3); // <<cos(phi1-phi2-phi3)>>
 Double_t sumOfW1st6 = fIntFlowSumOfEventWeights[0]->GetBinContent(1); // W_{<2>}
 Double_t sumOfW2nd6 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(3); // W_{<cos(phi1-phi2-phi3)>}
 Double_t sumOfWW6 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(6); // W_{<2>} * W_{<cos(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator6 = product6 - term1st6*term2nd6; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator6 = 0.;
 if(TMath::Abs(sumOfW1st6*sumOfW2nd6)>0.)
 {  
  denominator6 = 1.-sumOfWW6/(sumOfW1st6*sumOfW2nd6);
  if(TMath::Abs(denominator6)>0.)
  {  
   // covariance:
   Double_t covariance6 = numerator6/denominator6;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor6 = sumOfWW6/(sumOfW1st6*sumOfW2nd6);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(6,wPrefactor6*covariance6);
  } // end of if(TMath::Abs(denominator6)>0.)
 } // end of if(TMath::Abs(sumOfW1st6*sumOfW2nd6)>0.)
 
 // Cov(<2>,<sin(phi1-phi2-phi3)>):
 Double_t product7 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(7); // <<2><sin(phi1-phi2-phi3)>> 
 Double_t term1st7 = fIntFlowCorrelationsPro->GetBinContent(1); // <<2>>
 Double_t term2nd7 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(3); // <<sin(phi1-phi2-phi3)>>
 Double_t sumOfW1st7 = fIntFlowSumOfEventWeights[0]->GetBinContent(1); // W_{<2>}
 Double_t sumOfW2nd7 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(3); // W_{<sin(phi1-phi2-phi3)>}
 Double_t sumOfWW7 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(7); // W_{<2>} * W_{<sin(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator7 = product7 - term1st7*term2nd7; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator7 = 0.;
 if(TMath::Abs(sumOfW1st7*sumOfW2nd7)>0.)
 {  
  denominator7 = 1.-sumOfWW7/(sumOfW1st7*sumOfW2nd7);
  if(TMath::Abs(denominator7)>0.)
  {   
   // covariance:
   Double_t covariance7 = numerator7/denominator7;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor7 = sumOfWW7/(sumOfW1st7*sumOfW2nd7);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(7,wPrefactor7*covariance7);
  } // end of if(TMath::Abs(denominator7)>0.)
 } // end of if(TMath::Abs(sumOfW1st7*sumOfW2nd7)>0.)
 
 // Cov(<4>,<cos(phi1>):
 Double_t product8 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(8); // <<4><cos(phi1)>> 
 Double_t term1st8 = fIntFlowCorrelationsPro->GetBinContent(2); // <<4>>
 Double_t term2nd8 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(1); // <<cos(phi1)>>
 Double_t sumOfW1st8 = fIntFlowSumOfEventWeights[0]->GetBinContent(2); // W_{<4>}
 Double_t sumOfW2nd8 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(1); // W_{<cos(phi1)>}
 Double_t sumOfWW8 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(8); // W_{<4>} * W_{<cos(phi1)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator8 = product8 - term1st8*term2nd8; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator8 = 0.;
 if(TMath::Abs(sumOfW1st8*sumOfW2nd8)>0.)
 { 
  denominator8 = 1.-sumOfWW8/(sumOfW1st8*sumOfW2nd8);
  if(TMath::Abs(denominator8)>0.)
  {     
   // covariance:
   Double_t covariance8 = numerator8/denominator8;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor8 = sumOfWW8/(sumOfW1st8*sumOfW2nd8);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(8,wPrefactor8*covariance8);
  } // end of if(TMath::Abs(denominator8)>0.)
 } // end of if(TMath::Abs(sumOfW1st8*sumOfW2nd8)>0.)
 
 // Cov(<4>,<sin(phi1)>):
 Double_t product9 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(9); // <<4><sin(phi1)>> 
 Double_t term1st9 = fIntFlowCorrelationsPro->GetBinContent(2); // <<4>>
 Double_t term2nd9 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(1); // <<sin(phi1)>>
 Double_t sumOfW1st9 = fIntFlowSumOfEventWeights[0]->GetBinContent(2); // W_{<4>}
 Double_t sumOfW2nd9 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(1); // W_{<sin(phi1)>}
 Double_t sumOfWW9 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(9); // W_{<4>} * W_{<sin(phi1)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator9 = product9 - term1st9*term2nd9; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator9 = 0.;
 if(TMath::Abs(sumOfW1st9*sumOfW2nd9)>0.) 
 {
  denominator9 = 1.-sumOfWW9/(sumOfW1st9*sumOfW2nd9);
  if(TMath::Abs(denominator9)>0.)
  {     
   // covariance:
   Double_t covariance9 = numerator9/denominator9;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor9 = sumOfWW9/(sumOfW1st9*sumOfW2nd9);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(9,wPrefactor9*covariance9); 
  }
 } // end of if(TMath::Abs(sumOfW1st9*sumOfW2nd9)>0.) 
 
 // Cov(<4>,<cos(phi1+phi2)>):
 Double_t product10 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(10); // <<4><cos(phi1+phi2)>> 
 Double_t term1st10 = fIntFlowCorrelationsPro->GetBinContent(2); // <<4>>
 Double_t term2nd10 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(2); // <<cos(phi1+phi2)>>
 Double_t sumOfW1st10 = fIntFlowSumOfEventWeights[0]->GetBinContent(2); // W_{<4>}
 Double_t sumOfW2nd10 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(2); // W_{<cos(phi1+phi2)>}
 Double_t sumOfWW10 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(10); // W_{<4>} * W_{<cos(phi1+phi2)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator10 = product10 - term1st10*term2nd10; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator10 = 0.;
 if(TMath::Abs(sumOfW1st10*sumOfW2nd10)>0.) 
 { 
  denominator10 = 1.-sumOfWW10/(sumOfW1st10*sumOfW2nd10);
  if(TMath::Abs(denominator10)>0.) 
  { 
   // covariance:
   Double_t covariance10 = numerator10/denominator10;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor10 = sumOfWW10/(sumOfW1st10*sumOfW2nd10);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(10,wPrefactor10*covariance10);
  } // end of if(TMath::Abs(denominator10)>0.) 
 } // end of if(TMath::Abs(sumOfW1st10*sumOfW2nd10)>0.) 
 
 // Cov(<4>,<sin(phi1+phi2)>):
 Double_t product11 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(11); // <<4><sin(phi1+phi2)>> 
 Double_t term1st11 = fIntFlowCorrelationsPro->GetBinContent(2); // <<4>>
 Double_t term2nd11 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(2); // <<sin(phi1+phi2)>>
 Double_t sumOfW1st11 = fIntFlowSumOfEventWeights[0]->GetBinContent(2); // W_{<4>}
 Double_t sumOfW2nd11 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(2); // W_{<sin(phi1+phi2)>}
 Double_t sumOfWW11 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(11); // W_{<4>} * W_{<sin(phi1+phi2)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator11 = product11 - term1st11*term2nd11; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator11 = 0.;
 if(TMath::Abs(sumOfW1st11*sumOfW2nd11)>0.) 
 {  
  denominator11 = 1.-sumOfWW11/(sumOfW1st11*sumOfW2nd11);
  if(TMath::Abs(denominator11)>0.) 
  { 
   // covariance:
   Double_t covariance11 = numerator11/denominator11;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor11 = sumOfWW11/(sumOfW1st11*sumOfW2nd11);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(11,wPrefactor11*covariance11);
  } // end of if(TMath::Abs(denominator11)>0.) 
 } // end of if(TMath::Abs(sumOfW1st11*sumOfW2nd11)>0.) 

 // Cov(<4>,<cos(phi1-phi2-phi3)>):
 Double_t product12 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(12); // <<4><cos(phi1-phi2-phi3)>> 
 Double_t term1st12 = fIntFlowCorrelationsPro->GetBinContent(2); // <<4>>
 Double_t term2nd12 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(3); // <<cos(phi1-phi2-phi3)>>
 Double_t sumOfW1st12 = fIntFlowSumOfEventWeights[0]->GetBinContent(2); // W_{<4>}
 Double_t sumOfW2nd12 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(3); // W_{<cos(phi1-phi2-phi3)>}
 Double_t sumOfWW12 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(12); // W_{<4>} * W_{<cos(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator12 = product12 - term1st12*term2nd12; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator12 = 0.;
 if(TMath::Abs(sumOfW1st12*sumOfW2nd12)>0.) 
 {   
  denominator12 = 1.-sumOfWW12/(sumOfW1st12*sumOfW2nd12);
  if(TMath::Abs(denominator12)>0.) 
  { 
   // covariance:
   Double_t covariance12 = numerator12/denominator12;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor12 = sumOfWW12/(sumOfW1st12*sumOfW2nd12);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(12,wPrefactor12*covariance12);
  } // end of if(TMath::Abs(denominator12)>0.)
 } // end of if(TMath::Abs(sumOfW1st12*sumOfW2nd12)>0.)  

 // Cov(<4>,<sin(phi1-phi2-phi3)>):
 Double_t product13 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(13); // <<4><sin(phi1-phi2-phi3)>> 
 Double_t term1st13 = fIntFlowCorrelationsPro->GetBinContent(2); // <<4>>
 Double_t term2nd13 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(3); // <<sin(phi1-phi2-phi3)>>
 Double_t sumOfW1st13 = fIntFlowSumOfEventWeights[0]->GetBinContent(2); // W_{<4>}
 Double_t sumOfW2nd13 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(3); // W_{<sin(phi1-phi2-phi3)>}
 Double_t sumOfWW13 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(13); // W_{<4>} * W_{<sin(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator13 = product13 - term1st13*term2nd13; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator13 = 0.;
 if(TMath::Abs(sumOfW1st13*sumOfW2nd13)>0.) 
 {   
  denominator13 = 1.-sumOfWW13/(sumOfW1st13*sumOfW2nd13);
  if(TMath::Abs(denominator13)>0.) 
  { 
   // covariance:
   Double_t covariance13 = numerator13/denominator13;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor13 = sumOfWW13/(sumOfW1st13*sumOfW2nd13);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(13,wPrefactor13*covariance13);
  } // end of if(TMath::Abs(denominator13)>0.) 
 } // end of if(TMath::Abs(sumOfW1st13*sumOfW2nd13)>0.) 

 // Cov(<cos(phi1)>,<cos(phi1+phi2)>):
 Double_t product14 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(14); // <<cos(phi1)><cos(phi1+phi2)>> 
 Double_t term1st14 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(1); // <<cos(phi1)>>
 Double_t term2nd14 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(2); // <<cos(phi1+phi2)>>
 Double_t sumOfW1st14 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(1); // W_{<cos(phi1)>}
 Double_t sumOfW2nd14 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(2); // W_{<cos(phi1+phi2)>}
 Double_t sumOfWW14 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(14); // W_{<cos(phi1)>} * W_{<cos(phi1+phi2)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator14 = product14 - term1st14*term2nd14; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator14 = 0.;
 if(TMath::Abs(sumOfW1st14*sumOfW2nd14)>0.) 
 {  
  denominator14 = 1.-sumOfWW14/(sumOfW1st14*sumOfW2nd14);
  if(TMath::Abs(denominator14)>0.) 
  { 
   // covariance:
   Double_t covariance14 = numerator14/denominator14;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor14 = sumOfWW14/(sumOfW1st14*sumOfW2nd14);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(14,wPrefactor14*covariance14);
  } // end of if(TMath::Abs(denominator14)>0.) 
 } // end of if(TMath::Abs(sumOfW1st14*sumOfW2nd14)>0.) 

 // Cov(<cos(phi1)>,<sin(phi1+phi2)>):
 Double_t product15 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(15); // <<cos(phi1)><sin(phi1+phi2)>> 
 Double_t term1st15 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(1); // <<cos(phi1)>>
 Double_t term2nd15 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(2); // <<sin(phi1+phi2)>>
 Double_t sumOfW1st15 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(1); // W_{<cos(phi1)>}
 Double_t sumOfW2nd15 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(2); // W_{<sin(phi1+phi2)>}
 Double_t sumOfWW15 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(15); // W_{<cos(phi1)>} * W_{<sin(phi1+phi2)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator15 = product15 - term1st15*term2nd15; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator15 = 0.;
 if(TMath::Abs(sumOfW1st15*sumOfW2nd15)>0.) 
 {  
  denominator15 = 1.-sumOfWW15/(sumOfW1st15*sumOfW2nd15);
  if(TMath::Abs(denominator15)>0.) 
  { 
   // covariance:
   Double_t covariance15 = numerator15/denominator15;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor15 = sumOfWW15/(sumOfW1st15*sumOfW2nd15);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(15,wPrefactor15*covariance15);
  } // end of if(TMath::Abs(denominator15)>0.)  
 } // end of if(TMath::Abs(sumOfW1st15*sumOfW2nd15)>0.)  
 
 // Cov(<cos(phi1)>,<cos(phi1-phi2-phi3)>):
 Double_t product16 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(16); // <<cos(phi1)><cos(phi1-phi2-phi3)>> 
 Double_t term1st16 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(1); // <<cos(phi1)>>
 Double_t term2nd16 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(3); // <<cos(phi1-phi2-phi3)>>
 Double_t sumOfW1st16 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(1); // W_{<cos(phi1)>}
 Double_t sumOfW2nd16 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(3); // W_{<cos(phi1-phi2-phi3)>}
 Double_t sumOfWW16 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(16); // W_{<cos(phi1)>} * W_{<cos(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator16 = product16 - term1st16*term2nd16; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator16 = 0.;
 if(TMath::Abs(sumOfW1st16*sumOfW2nd16)>0.) 
 {   
  denominator16 = 1.-sumOfWW16/(sumOfW1st16*sumOfW2nd16);
  if(TMath::Abs(denominator16)>0.) 
  {   
   // covariance:
   Double_t covariance16 = numerator16/denominator16;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor16 = sumOfWW16/(sumOfW1st16*sumOfW2nd16);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(16,wPrefactor16*covariance16);
  } // end of if(TMath::Abs(denominator16)>0.)
 } // end ofif(TMath::Abs(sumOfW1st16*sumOfW2nd16)>0.)  
 
 // Cov(<cos(phi1)>,<sin(phi1-phi2-phi3)>):
 Double_t product17 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(17); // <<cos(phi1)><sin(phi1-phi2-phi3)>> 
 Double_t term1st17 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(1); // <<cos(phi1)>>
 Double_t term2nd17 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(3); // <<sin(phi1-phi2-phi3)>>
 Double_t sumOfW1st17 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(1); // W_{<cos(phi1)>}
 Double_t sumOfW2nd17 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(3); // W_{<sin(phi1-phi2-phi3)>}
 Double_t sumOfWW17 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(17); // W_{<cos(phi1)>} * W_{<sin(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator17 = product17 - term1st17*term2nd17; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator17 = 0.;
 if(TMath::Abs(sumOfW1st17*sumOfW2nd17)>0.) 
 {
  denominator17 = 1.-sumOfWW17/(sumOfW1st17*sumOfW2nd17);
  if(TMath::Abs(denominator17)>0.) 
  {   
   // covariance:
   Double_t covariance17 = numerator17/denominator17;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor17 = sumOfWW17/(sumOfW1st17*sumOfW2nd17);
    // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(17,wPrefactor17*covariance17);
  } // end of if(TMath::Abs(denominator17)>0.) 
 } // end of if(TMath::Abs(sumOfW1st17*sumOfW2nd17)>0.) 

 // Cov(<sin(phi1)>,<cos(phi1+phi2)>):
 Double_t product18 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(18); // <<sin(phi1)><cos(phi1+phi2)>> 
 Double_t term1st18 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(1); // <<sin(phi1)>>
 Double_t term2nd18 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(2); // <<cos(phi1+phi2)>>
 Double_t sumOfW1st18 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(1); // W_{<sin(phi1)>}
 Double_t sumOfW2nd18 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(2); // W_{<cos(phi1+phi2)>}
 Double_t sumOfWW18 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(18); // W_{<sin(phi1)>} * W_{<cos(phi1+phi2)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator18 = product18 - term1st18*term2nd18; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator18 = 0.;
 if(TMath::Abs(sumOfW1st18*sumOfW2nd18)>0.) 
 { 
  denominator18 = 1.-sumOfWW18/(sumOfW1st18*sumOfW2nd18);
  if(TMath::Abs(denominator18)>0.) 
  {   
   // covariance:
   Double_t covariance18 = numerator18/denominator18;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor18 = sumOfWW18/(sumOfW1st18*sumOfW2nd18);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(18,wPrefactor18*covariance18); 
  } // end of if(TMath::Abs(denominator18)>0.) 
 } // end of if(TMath::Abs(sumOfW1st18*sumOfW2nd18)>0.) 

 // Cov(<sin(phi1)>,<sin(phi1+phi2)>):
 Double_t product19 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(19); // <<sin(phi1)><sin(phi1+phi2)>> 
 Double_t term1st19 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(1); // <<sin(phi1)>>
 Double_t term2nd19 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(2); // <<sin(phi1+phi2)>>
 Double_t sumOfW1st19 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(1); // W_{<sin(phi1)>}
 Double_t sumOfW2nd19 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(2); // W_{<sin(phi1+phi2)>}
 Double_t sumOfWW19 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(19); // W_{<sin(phi1)>} * W_{<sin(phi1+phi2)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator19 = product19 - term1st19*term2nd19; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator19 = 0.;
 if(TMath::Abs(sumOfW1st19*sumOfW2nd19)>0.) 
 { 
  denominator19 = 1.-sumOfWW19/(sumOfW1st19*sumOfW2nd19);
  if(TMath::Abs(denominator19)>0.) 
  {   
   // covariance:
   Double_t covariance19 = numerator19/denominator19;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor19 = sumOfWW19/(sumOfW1st19*sumOfW2nd19);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(19,wPrefactor19*covariance19);
  } // end of if(TMath::Abs(denominator19)>0.)
 } // end of if(TMath::Abs(sumOfW1st19*sumOfW2nd19)>0.)
 
 // Cov(<sin(phi1)>,<cos(phi1-phi2-phi3)>):
 Double_t product20 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(20); // <<sin(phi1)><cos(phi1-phi2-phi3)>> 
 Double_t term1st20 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(1); // <<sin(phi1)>>
 Double_t term2nd20 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(3); // <<cos(phi1-phi2-phi3)>>
 Double_t sumOfW1st20 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(1); // W_{<sin(phi1)>}
 Double_t sumOfW2nd20 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(3); // W_{<cos(phi1-phi2-phi3)>}
 Double_t sumOfWW20 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(20); // W_{<sin(phi1)>} * W_{<cos(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator20 = product20 - term1st20*term2nd20; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator20 = 0.;
 if(TMath::Abs(sumOfW1st20*sumOfW2nd20)>0.)
 { 
  denominator20 = 1.-sumOfWW20/(sumOfW1st20*sumOfW2nd20);
  if(TMath::Abs(denominator20)>0.) 
  { 
   // covariance:
   Double_t covariance20 = numerator20/denominator20;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor20 = sumOfWW20/(sumOfW1st20*sumOfW2nd20);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(20,wPrefactor20*covariance20);
  } // end of if(TMath::Abs(denominator20)>0.) 
 } // end of if(TMath::Abs(sumOfW1st20*sumOfW2nd20)>0.)

 // Cov(<sin(phi1)>,<sin(phi1-phi2-phi3)>):
 Double_t product21 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(21); // <<sin(phi1)><sin(phi1-phi2-phi3)>> 
 Double_t term1st21 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(1); // <<sin(phi1)>>
 Double_t term2nd21 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(3); // <<sin(phi1-phi2-phi3)>>
 Double_t sumOfW1st21 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(1); // W_{<sin(phi1)>}
 Double_t sumOfW2nd21 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(3); // W_{<sin(phi1-phi2-phi3)>}
 Double_t sumOfWW21 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(21); // W_{<sin(phi1)>} * W_{<sin(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator21 = product21 - term1st21*term2nd21; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator21 = 0.;
 if(TMath::Abs(sumOfW1st21*sumOfW2nd21)>0.)
 { 
  denominator21 = 1.-sumOfWW21/(sumOfW1st21*sumOfW2nd21);
  if(TMath::Abs(denominator21)>0.) 
  {   
   // covariance:
   Double_t covariance21 = numerator21/denominator21;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor21 = sumOfWW21/(sumOfW1st21*sumOfW2nd21);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(21,wPrefactor21*covariance21);
  } // end of if(TMath::Abs(denominator21)>0.)
 } // end of if(TMath::Abs(sumOfW1st21*sumOfW2nd21)>0.)

 // Cov(<cos(phi1+phi2)>,<sin(phi1+phi2)>):
 Double_t product22 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(22); // <<cos(phi1+phi2)><sin(phi1+phi2)>> 
 Double_t term1st22 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(2); // <<cos(phi1+phi2)>>
 Double_t term2nd22 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(2); // <<sin(phi1+phi2)>>
 Double_t sumOfW1st22 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(2); // W_{<cos(phi1+phi2)>}
 Double_t sumOfW2nd22 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(2); // W_{<sin(phi1+phi2)>}
 Double_t sumOfWW22 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(22); // W_{<cos(phi1+phi2)>} * W_{<sin(phi1+phi2)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator22 = product22 - term1st22*term2nd22; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator22 = 0.;
 if(TMath::Abs(sumOfW1st22*sumOfW2nd22)>0.)
 { 
  denominator22 = 1.-sumOfWW22/(sumOfW1st22*sumOfW2nd22);
  if(TMath::Abs(denominator22)>0.) 
  {   
   // covariance:
   Double_t covariance22 = numerator22/denominator22;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor22 = sumOfWW22/(sumOfW1st22*sumOfW2nd22);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(22,wPrefactor22*covariance22);
  } // end of if(TMath::Abs(denominator22)>0.) 
 } // end of if(TMath::Abs(sumOfW1st22*sumOfW2nd22)>0.) 

 // Cov(<cos(phi1+phi2)>,<cos(phi1-phi2-phi3)>):
 Double_t product23 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(23); // <<cos(phi1+phi2)><cos(phi1-phi2-phi3)>> 
 Double_t term1st23 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(2); // <<cos(phi1+phi2)>>
 Double_t term2nd23 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(3); // <<cos(phi1-phi2-phi3)>>
 Double_t sumOfW1st23 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(2); // W_{<cos(phi1+phi2)>}
 Double_t sumOfW2nd23 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(3); // W_{<cos(phi1-phi2-phi3)>}
 Double_t sumOfWW23 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(23); // W_{<cos(phi1+phi2)>} * W_{<cos(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator23 = product23 - term1st23*term2nd23; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator23 = 0.;
 if(TMath::Abs(sumOfW1st23*sumOfW2nd23)>0.)
 {  
  denominator23 = 1.-sumOfWW23/(sumOfW1st23*sumOfW2nd23);
  if(TMath::Abs(denominator23)>0.) 
  {   
   // covariance:
   Double_t covariance23 = numerator23/denominator23;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor23 = sumOfWW23/(sumOfW1st23*sumOfW2nd23);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(23,wPrefactor23*covariance23);
  } // end of if(TMath::Abs(denominator23)>0.) 
 } // end of if(TMath::Abs(sumOfW1st23*sumOfW2nd23)>0.)
 
 // Cov(<cos(phi1+phi2)>,<sin(phi1-phi2-phi3)>):
 Double_t product24 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(24); // <<cos(phi1+phi2)><sin(phi1-phi2-phi3)>> 
 Double_t term1st24 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(2); // <<cos(phi1+phi2)>>
 Double_t term2nd24 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(3); // <<sin(phi1-phi2-phi3)>>
 Double_t sumOfW1st24 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(2); // W_{<cos(phi1+phi2)>}
 Double_t sumOfW2nd24 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(3); // W_{<sin(phi1-phi2-phi3)>}
 Double_t sumOfWW24 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(24); // W_{<cos(phi1+phi2)>} * W_{<sin(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator24 = product24 - term1st24*term2nd24; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator24 = 0.;
 if(TMath::Abs(sumOfW1st24*sumOfW2nd24)>0.)
 {   
  denominator24 = 1.-sumOfWW24/(sumOfW1st24*sumOfW2nd24);
  if(TMath::Abs(denominator24)>0.) 
  {   
   // covariance:
   Double_t covariance24 = numerator24/denominator24;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor24 = sumOfWW24/(sumOfW1st24*sumOfW2nd24);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(24,wPrefactor24*covariance24);
  } // end of if(TMath::Abs(denominator24)>0.)  
 } // end of if(TMath::Abs(sumOfW1st24*sumOfW2nd24)>0.)

 // Cov(<sin(phi1+phi2)>,<cos(phi1-phi2-phi3)>):
 Double_t product25 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(25); // <<sin(phi1+phi2)><cos(phi1-phi2-phi3)>> 
 Double_t term1st25 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(2); // <<sin(phi1+phi2)>>
 Double_t term2nd25 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(3); // <<cos(phi1-phi2-phi3)>>
 Double_t sumOfW1st25 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(2); // W_{<sin(phi1+phi2)>}
 Double_t sumOfW2nd25 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(3); // W_{<cos(phi1-phi2-phi3)>}
 Double_t sumOfWW25 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(25); // W_{<sin(phi1+phi2)>} * W_{<cos(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator25 = product25 - term1st25*term2nd25; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator25 = 0.;
 if(TMath::Abs(sumOfW1st25*sumOfW2nd25)>0.)
 { 
  denominator25 = 1.-sumOfWW25/(sumOfW1st25*sumOfW2nd25);
  if(TMath::Abs(denominator25)>0.) 
  { 
   // covariance:
   Double_t covariance25 = numerator25/denominator25;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor25 = sumOfWW25/(sumOfW1st25*sumOfW2nd25);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(25,wPrefactor25*covariance25);
  } // end of if(TMath::Abs(denominator25)>0.)
 } // end of if(TMath::Abs(sumOfW1st25*sumOfW2nd25)>0.)
 
 // Cov(<sin(phi1+phi2)>,<sin(phi1-phi2-phi3)>):
 Double_t product26 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(26); // <<sin(phi1+phi2)><sin(phi1-phi2-phi3)>> 
 Double_t term1st26 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(2); // <<sin(phi1+phi2)>>
 Double_t term2nd26 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(3); // <<sin(phi1-phi2-phi3)>>
 Double_t sumOfW1st26 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(2); // W_{<sin(phi1+phi2)>}
 Double_t sumOfW2nd26 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(3); // W_{<sin(phi1-phi2-phi3)>}
 Double_t sumOfWW26 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(26); // W_{<sin(phi1+phi2)>} * W_{<sin(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator26 = product26 - term1st26*term2nd26; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator26 = 0.;
 if(TMath::Abs(sumOfW1st26*sumOfW2nd26)>0.)
 { 
  denominator26 = 1.-sumOfWW26/(sumOfW1st26*sumOfW2nd26);
  if(TMath::Abs(denominator26)>0.) 
  { 
   // covariance:
   Double_t covariance26 = numerator26/denominator26;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor26 = sumOfWW26/(sumOfW1st26*sumOfW2nd26);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(26,wPrefactor26*covariance26);
  } // end of if(TMath::Abs(denominator26)>0.) 
 } // end of if(TMath::Abs(sumOfW1st26*sumOfW2nd26)>0.)
 
 // Cov(<cos(phi1-phi2-phi3)>,<sin(phi1-phi2-phi3)>):
 Double_t product27 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(27); // <<cos(phi1-phi2-phi3)><sin(phi1-phi2-phi3)>> 
 Double_t term1st27 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(3); // <<cos(phi1-phi2-phi3)>>
 Double_t term2nd27 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(3); // <<sin(phi1-phi2-phi3)>>
 Double_t sumOfW1st27 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(3); // W_{<cos(phi1-phi2-phi3)>}
 Double_t sumOfW2nd27 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(3); // W_{<sin(phi1-phi2-phi3)>}
 Double_t sumOfWW27 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(27); // W_{<cos(phi1-phi2-phi3)>} * W_{<sin(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator27 = product27 - term1st27*term2nd27; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator27 = 0.;
 if(TMath::Abs(sumOfW1st27*sumOfW2nd27)>0.)
 { 
  denominator27 = 1.-sumOfWW27/(sumOfW1st27*sumOfW2nd27);
  if(TMath::Abs(denominator27)>0.) 
  { 
   // covariance:
   Double_t covariance27 = numerator27/denominator27;
   // weight dependent prefactor for covariance:
   Double_t wPrefactor27 = sumOfWW27/(sumOfW1st27*sumOfW2nd27);
   // finally, store "weighted" covariance:
   fIntFlowCovariancesNUA->SetBinContent(27,wPrefactor27*covariance27);
  } // end of if(TMath::Abs(denominator27)>0.) 
 } // end of if(TMath::Abs(sumOfW1st27*sumOfW2nd27)>0.)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateCovariancesNUAIntFlow()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::FinalizeCorrelationsIntFlow() 
{
 // From profile fIntFlowCorrelationsPro access measured correlations and spread, 
 // correctly calculate the statistical errors and store the final results and 
 // statistical errors for correlations in histogram fIntFlowCorrelationsHist.
 //
 // Remark: Statistical error of correlation is calculated as:
 //
 //          statistical error = termA * spread * termB:
 //          termA = sqrt{sum_{i=1}^{N} w^2}/(sum_{i=1}^{N} w)
 //          termB = 1/sqrt(1-termA^2)   
 //
   
 for(Int_t ci=1;ci<=4;ci++) // correlation index
 {
  if(fIntFlowCorrelationsPro->GetBinEffectiveEntries(ci) < 2 || fIntFlowSquaredCorrelationsPro->GetBinEffectiveEntries(ci) < 2)
  {
   fIntFlowCorrelationsPro->SetBinError(ci,0.);
   fIntFlowSquaredCorrelationsPro->SetBinError(ci,0.);
   continue;
  } 
  Double_t correlation = fIntFlowCorrelationsPro->GetBinContent(ci);
  Double_t squaredCorrelation = fIntFlowSquaredCorrelationsPro->GetBinContent(ci);
  Double_t spread = 0.;
  if(squaredCorrelation-correlation*correlation >= 0.)
  {
   spread = pow(squaredCorrelation-correlation*correlation,0.5);
  } else
    {
     cout<<endl;
     cout<<Form(" WARNING: Imaginary 'spread' for %d-particle correlation!!!! ",2*ci)<<endl;
     cout<<endl;
    }
  Double_t sumOfLinearEventWeights = fIntFlowSumOfEventWeights[0]->GetBinContent(ci);
  Double_t sumOfQuadraticEventWeights = fIntFlowSumOfEventWeights[1]->GetBinContent(ci);
  Double_t termA = 0.;
  Double_t termB = 0.;
  if(TMath::Abs(sumOfLinearEventWeights) > 0.) // to be improved - shall I omitt here Abs() ?
  {
   termA = pow(sumOfQuadraticEventWeights,0.5)/sumOfLinearEventWeights;
  } else
    {
     cout<<endl;
     cout<<" WARNING (QC): sumOfLinearEventWeights == 0 in method FinalizeCorrelationsIntFlow() !!!!"<<endl;
     cout<<"               (for "<<2*ci<<"-particle correlation)"<<endl;
     cout<<endl;
    }
  if(1.-pow(termA,2.) > 0.)
  {
   termB = 1./pow(1-pow(termA,2.),0.5);
  } else
    {
     cout<<endl;
     cout<<" WARNING (QC): 1.-pow(termA,2.) <= 0 in method FinalizeCorrelationsIntFlow() !!!!"<<endl;   
     cout<<"               (for "<<2*ci<<"-particle correlation)"<<endl;
     cout<<endl;
    }     
  Double_t statisticalError = termA * spread * termB;
  fIntFlowCorrelationsHist->SetBinContent(ci,correlation);
  fIntFlowCorrelationsHist->SetBinError(ci,statisticalError);
 } // end of for(Int_t ci=1;ci<=4;ci++) // correlation index     
 
 // Versus multiplicity: 
 if(!fCalculateCumulantsVsM){return;}
 for(Int_t ci=0;ci<=3;ci++) // correlation index
 {
  Int_t nBins = fIntFlowCorrelationsVsMPro[ci]->GetNbinsX(); 
  for(Int_t b=1;b<=nBins;b++) // looping over multiplicity bins
  {
   if(fIntFlowCorrelationsVsMPro[ci]->GetBinEffectiveEntries(b) < 2 || fIntFlowSquaredCorrelationsVsMPro[ci]->GetBinEffectiveEntries(b) < 2)
   {
    fIntFlowCorrelationsVsMPro[ci]->SetBinError(b,0.);
    fIntFlowSquaredCorrelationsVsMPro[ci]->SetBinError(b,0.);
    continue;
   } 
   Double_t correlationVsM = fIntFlowCorrelationsVsMPro[ci]->GetBinContent(b);
   Double_t squaredCorrelationVsM = fIntFlowSquaredCorrelationsVsMPro[ci]->GetBinContent(b);
   Double_t spreadVsM = 0.;
   if(squaredCorrelationVsM-correlationVsM*correlationVsM >= 0.)
   {
    spreadVsM = pow(squaredCorrelationVsM-correlationVsM*correlationVsM,0.5);
   } else
     {
      cout<<endl;
      cout<<Form(" WARNING (QC): Imaginary 'spreadVsM' for ci = %d, bin = %d, entries = %f !!!!",
                 ci,b,fIntFlowCorrelationsVsMPro[ci]->GetBinEffectiveEntries(b))<<endl; 
      cout<<endl;
     }     
   Double_t sumOfLinearEventWeightsVsM = fIntFlowSumOfEventWeightsVsM[ci][0]->GetBinContent(b);
   Double_t sumOfQuadraticEventWeightsVsM = fIntFlowSumOfEventWeightsVsM[ci][1]->GetBinContent(b);
   Double_t termAVsM = 0.;
   Double_t termBVsM = 0.;
   if(sumOfLinearEventWeightsVsM > 0.) 
   {
    termAVsM = pow(sumOfQuadraticEventWeightsVsM,0.5)/sumOfLinearEventWeightsVsM;
   }
   if(1.-pow(termAVsM,2.) > 0.)
   {
    termBVsM = 1./pow(1-pow(termAVsM,2.),0.5);
   }     
   Double_t statisticalErrorVsM = termAVsM * spreadVsM * termBVsM;
   fIntFlowCorrelationsVsMHist[ci]->SetBinContent(b,correlationVsM);
   fIntFlowCorrelationsVsMHist[ci]->SetBinError(b,statisticalErrorVsM);  
  } // end of for(Int_t b=1;b<=nBins;b++)
 } // end of for(Int_t ci=1;ci<=4;ci++) // correlation index                                                        
                                                                                                                           
} // end of AliFlowAnalysisWithQCumulants::FinalizeCorrelationsIntFlow()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::FillAverageMultiplicities(Int_t nRP)
{
 // Fill profile fAverageMultiplicity to hold average multiplicities and 
 // number of events for events with nRP>=0, nRP>=1, ... , and nRP>=8
 
 // Binning of fAverageMultiplicity is organized as follows:
 //  1st bin: all events (including the empty ones)
 //  2nd bin: event with # of RPs greater or equal to 1
 //  3rd bin: event with # of RPs greater or equal to 2
 //  4th bin: event with # of RPs greater or equal to 3
 //  5th bin: event with # of RPs greater or equal to 4
 //  6th bin: event with # of RPs greater or equal to 5
 //  7th bin: event with # of RPs greater or equal to 6
 //  8th bin: event with # of RPs greater or equal to 7
 //  9th bin: event with # of RPs greater or equal to 8
 
 if(nRP<0)
 {
  cout<<endl;
  cout<<" WARNING (QC): nRP<0 in in AFAWQC::FAM() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 
 for(Int_t i=0;i<9;i++)
 {
  if(nRP>=i){fAvMultiplicity->Fill(i+0.5,nRP,1);}
 }
 
} // end of AliFlowAnalysisWithQCumulants::FillAverageMultiplicities(nRP)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateCumulantsIntFlow()
{ 
 // a) Calculate Q-cumulants from the measured multiparticle correlations;
 // b) Propagate the statistical errors from measured multiparticle correlations to statistical errors of Q-cumulants;  
 // c) Remark: Q-cumulants calculated in this method are biased by non-uniform acceptance of detector !!!! 
 //            Method CalculateQcumulantsCorrectedForNUAIntFlow() is called afterwards to correct for this bias;
 // d) Store the results and statistical error of Q-cumulants in histogram fIntFlowQcumulants.
 //    Binning of fIntFlowQcumulants is organized as follows:
 //
 //            1st bin: QC{2}
 //            2nd bin: QC{4}
 //            3rd bin: QC{6}
 //            4th bin: QC{8}
 //
 
 // Correlations:
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>> 
 Double_t four = fIntFlowCorrelationsHist->GetBinContent(2); // <<4>>  
 Double_t six = fIntFlowCorrelationsHist->GetBinContent(3); // <<6>> 
 Double_t eight = fIntFlowCorrelationsHist->GetBinContent(4); // <<8>>  
 // Statistical errors of average 2-, 4-, 6- and 8-particle azimuthal correlations:
 Double_t twoError = fIntFlowCorrelationsHist->GetBinError(1); // statistical error of <2>  
 Double_t fourError = fIntFlowCorrelationsHist->GetBinError(2); // statistical error of <4>   
 Double_t sixError = fIntFlowCorrelationsHist->GetBinError(3); // statistical error of <6> 
 Double_t eightError = fIntFlowCorrelationsHist->GetBinError(4); // statistical error of <8> 
 // Covariances (multiplied by prefactor depending on weights - see comments in CalculateCovariancesIntFlow()):
 Double_t wCov24 = 0.; // Cov(<2>,<4>) * prefactor(w_<2>,w_<4>)
 Double_t wCov26 = 0.; // Cov(<2>,<6>) * prefactor(w_<2>,w_<6>)
 Double_t wCov28 = 0.; // Cov(<2>,<8>) * prefactor(w_<2>,w_<8>)
 Double_t wCov46 = 0.; // Cov(<4>,<6>) * prefactor(w_<4>,w_<6>)
 Double_t wCov48 = 0.; // Cov(<4>,<8>) * prefactor(w_<4>,w_<8>)
 Double_t wCov68 = 0.; // Cov(<6>,<8>) * prefactor(w_<6>,w_<8>)  
 if(!fForgetAboutCovariances)
 {
  wCov24 = fIntFlowCovariances->GetBinContent(1); // Cov(<2>,<4>) * prefactor(w_<2>,w_<4>)
  wCov26 = fIntFlowCovariances->GetBinContent(2); // Cov(<2>,<6>) * prefactor(w_<2>,w_<6>)
  wCov28 = fIntFlowCovariances->GetBinContent(3); // Cov(<2>,<8>) * prefactor(w_<2>,w_<8>)
  wCov46 = fIntFlowCovariances->GetBinContent(4); // Cov(<4>,<6>) * prefactor(w_<4>,w_<6>)
  wCov48 = fIntFlowCovariances->GetBinContent(5); // Cov(<4>,<8>) * prefactor(w_<4>,w_<8>)
  wCov68 = fIntFlowCovariances->GetBinContent(6); // Cov(<6>,<8>) * prefactor(w_<6>,w_<8>) 
 }
 // Q-cumulants: 
 Double_t qc2 = 0.; // QC{2}
 Double_t qc4 = 0.; // QC{4}
 Double_t qc6 = 0.; // QC{6}
 Double_t qc8 = 0.; // QC{8}
 if(TMath::Abs(two) > 0.){qc2 = two;} 
 if(TMath::Abs(four) > 0.){qc4 = four-2.*pow(two,2.);} 
 if(TMath::Abs(six) > 0.){qc6 = six-9.*two*four+12.*pow(two,3.);} 
 if(TMath::Abs(eight) > 0.){qc8 = eight-16.*two*six-18.*pow(four,2.)+144.*pow(two,2.)*four-144.*pow(two,4.);} 
 // Statistical errors of Q-cumulants:       
 Double_t qc2Error = 0.;
 Double_t qc4Error = 0.;
 Double_t qc6Error = 0.;
 Double_t qc8Error = 0.; 
 // Squared statistical errors of Q-cumulants:       
 //Double_t qc2ErrorSquared = 0.;
 Double_t qc4ErrorSquared = 0.;
 Double_t qc6ErrorSquared = 0.;
 Double_t qc8ErrorSquared = 0.;        
 // Statistical error of QC{2}:              
 qc2Error = twoError;                                                 
 // Statistical error of QC{4}:              
 qc4ErrorSquared = 16.*pow(two,2.)*pow(twoError,2)+pow(fourError,2.)
                 - 8.*two*wCov24;                     
 if(qc4ErrorSquared>0.)
 {
  qc4Error = pow(qc4ErrorSquared,0.5);
 } else 
   {
    cout<<" WARNING (QC): Statistical error of QC{4} is imaginary !!!!"<<endl;
   }                                           
 // Statistical error of QC{6}:              
 qc6ErrorSquared = 81.*pow(4.*pow(two,2.)-four,2.)*pow(twoError,2.)
                 + 81.*pow(two,2.)*pow(fourError,2.)
                 + pow(sixError,2.)
                 - 162.*two*(4.*pow(two,2.)-four)*wCov24
                 + 18.*(4.*pow(two,2.)-four)*wCov26
                 - 18.*two*wCov46;                     
 if(qc6ErrorSquared>0.)
 {
  qc6Error = pow(qc6ErrorSquared,0.5);
 } else 
   {
    cout<<" WARNING (QC): Statistical error of QC{6} is imaginary !!!!"<<endl;
   }                       
 // Statistical error of QC{8}:              
 qc8ErrorSquared = 256.*pow(36.*pow(two,3.)-18.*four*two+six,2.)*pow(twoError,2.)
                 + 1296.*pow(4.*pow(two,2.)-four,2.)*pow(fourError,2.)
                 + 256.*pow(two,2.)*pow(sixError,2.)
                 + pow(eightError,2.)
                 - 1152.*(36.*pow(two,3.)-18.*four*two+six)*(4.*pow(two,2.)-four)*wCov24
                 + 512.*two*(36.*pow(two,3.)-18.*four*two+six)*wCov26
                 - 32.*(36.*pow(two,3.)-18.*four*two+six)*wCov28
                 - 1152.*two*(4.*pow(two,2.)-four)*wCov46
                 + 72.*(4.*pow(two,2.)-four)*wCov48
                 - 32.*two*wCov68;      
 if(qc8ErrorSquared>0.)
 {
  qc8Error = pow(qc8ErrorSquared,0.5);
 } else 
   {
    cout<<"WARNING (QC): Statistical error of QC{8} is imaginary !!!!"<<endl;
   }
 // Store the results and statistical errors for Q-cumulants:
 if(TMath::Abs(qc2)>0.)
 {
  fIntFlowQcumulants->SetBinContent(1,qc2);
  fIntFlowQcumulants->SetBinError(1,qc2Error);
 }
 if(TMath::Abs(qc4)>0.)
 {
  fIntFlowQcumulants->SetBinContent(2,qc4);
  fIntFlowQcumulants->SetBinError(2,qc4Error);
 }
 if(TMath::Abs(qc6)>0.)
 {
  fIntFlowQcumulants->SetBinContent(3,qc6);
  fIntFlowQcumulants->SetBinError(3,qc6Error);
 }
 if(TMath::Abs(qc8)>0.)
 {
  fIntFlowQcumulants->SetBinContent(4,qc8); 
  fIntFlowQcumulants->SetBinError(4,qc8Error);
 } 
 
 // Versus multiplicity: 
 if(!fCalculateCumulantsVsM){return;}
 Int_t nBins = fIntFlowCorrelationsVsMPro[0]->GetNbinsX(); // to be improved (hardwired 0) 
 Double_t value[4] = {0.}; // QCs vs M
 Double_t error[4] = {0.}; // error of QCs vs M
 Double_t dSum1[4] = {0.}; // sum value_i/(error_i)^2
 Double_t dSum2[4] = {0.}; // sum 1/(error_i)^2
 for(Int_t b=1;b<=nBins;b++)
 {
  // Correlations:
  two = fIntFlowCorrelationsVsMHist[0]->GetBinContent(b); // <<2>> 
  four = fIntFlowCorrelationsVsMHist[1]->GetBinContent(b); // <<4>>  
  six = fIntFlowCorrelationsVsMHist[2]->GetBinContent(b); // <<6>> 
  eight = fIntFlowCorrelationsVsMHist[3]->GetBinContent(b); // <<8>>  
  // Statistical errors of average 2-, 4-, 6- and 8-particle azimuthal correlations:
  twoError = fIntFlowCorrelationsVsMHist[0]->GetBinError(b); // statistical error of <2>  
  fourError = fIntFlowCorrelationsVsMHist[1]->GetBinError(b); // statistical error of <4>   
  sixError = fIntFlowCorrelationsVsMHist[2]->GetBinError(b); // statistical error of <6> 
  eightError = fIntFlowCorrelationsVsMHist[3]->GetBinError(b); // statistical error of <8> 
  // Covariances (multiplied by prefactor depending on weights - see comments in CalculateCovariancesIntFlow()):
  if(!fForgetAboutCovariances)
  {
   wCov24 = fIntFlowCovariancesVsM[0]->GetBinContent(b); // Cov(<2>,<4>) * prefactor(w_<2>,w_<4>)
   wCov26 = fIntFlowCovariancesVsM[1]->GetBinContent(b); // Cov(<2>,<6>) * prefactor(w_<2>,w_<6>)
   wCov28 = fIntFlowCovariancesVsM[2]->GetBinContent(b); // Cov(<2>,<8>) * prefactor(w_<2>,w_<8>)
   wCov46 = fIntFlowCovariancesVsM[3]->GetBinContent(b); // Cov(<4>,<6>) * prefactor(w_<4>,w_<6>)
   wCov48 = fIntFlowCovariancesVsM[4]->GetBinContent(b); // Cov(<4>,<8>) * prefactor(w_<4>,w_<8>)
   wCov68 = fIntFlowCovariancesVsM[5]->GetBinContent(b); // Cov(<6>,<8>) * prefactor(w_<6>,w_<8>) 
  }
  // Q-cumulants: 
  qc2 = 0.; // QC{2}
  qc4 = 0.; // QC{4}
  qc6 = 0.; // QC{6}
  qc8 = 0.; // QC{8}
  if(TMath::Abs(two) > 0.){qc2 = two;} 
  if(TMath::Abs(four) > 0.){qc4 = four-2.*pow(two,2.);} 
  if(TMath::Abs(six) > 0.){qc6 = six-9.*two*four+12.*pow(two,3.);} 
  if(TMath::Abs(eight) > 0.){qc8 = eight-16.*two*six-18.*pow(four,2.)+144.*pow(two,2.)*four-144.*pow(two,4.);}  
  // Statistical errors of Q-cumulants:       
  qc2Error = 0.;
  qc4Error = 0.;
  qc6Error = 0.;
  qc8Error = 0.; 
  // Squared statistical errors of Q-cumulants:       
  //Double_t qc2ErrorSquared = 0.;
  qc4ErrorSquared = 0.;
  qc6ErrorSquared = 0.;
  qc8ErrorSquared = 0.;    
  // Statistical error of QC{2}:              
  qc2Error = twoError;                                             
  // Statistical error of QC{4}:              
  qc4ErrorSquared = 16.*pow(two,2.)*pow(twoError,2)+pow(fourError,2.)
                  - 8.*two*wCov24;                     
  if(qc4ErrorSquared>0.)
  {
   qc4Error = pow(qc4ErrorSquared,0.5);
  } else 
    {
     // cout<<"WARNING: Statistical error of QC{4} is imaginary in multiplicity bin "<<b<<" !!!!"<<endl;
    }                                       
  // Statistical error of QC{6}:              
  qc6ErrorSquared = 81.*pow(4.*pow(two,2.)-four,2.)*pow(twoError,2.)
                  + 81.*pow(two,2.)*pow(fourError,2.)
                  + pow(sixError,2.)
                  - 162.*two*(4.*pow(two,2.)-four)*wCov24
                  + 18.*(4.*pow(two,2.)-four)*wCov26
                  - 18.*two*wCov46;                     
  if(qc6ErrorSquared>0.)
  {
   qc6Error = pow(qc6ErrorSquared,0.5);
  } else 
    {
     // cout<<"WARNING: Statistical error of QC{6} is imaginary in multiplicity bin "<<b<<" !!!!"<<endl;
    }                            
  // Statistical error of QC{8}:              
  qc8ErrorSquared = 256.*pow(36.*pow(two,3.)-18.*four*two+six,2.)*pow(twoError,2.)
                  + 1296.*pow(4.*pow(two,2.)-four,2.)*pow(fourError,2.)
                  + 256.*pow(two,2.)*pow(sixError,2.)
                  + pow(eightError,2.)
                  - 1152.*(36.*pow(two,3.)-18.*four*two+six)*(4.*pow(two,2.)-four)*wCov24
                  + 512.*two*(36.*pow(two,3.)-18.*four*two+six)*wCov26
                  - 32.*(36.*pow(two,3.)-18.*four*two+six)*wCov28
                  - 1152.*two*(4.*pow(two,2.)-four)*wCov46
                  + 72.*(4.*pow(two,2.)-four)*wCov48
                  - 32.*two*wCov68;      
  if(qc8ErrorSquared>0.)
  {
   qc8Error = pow(qc8ErrorSquared,0.5);
  } else 
    {
     // cout<<"WARNING: Statistical error of QC{8} is imaginary in multiplicity bin "<<b<<" !!!!"<<endl;
    }
  // Store the results and statistical errors for Q-cumulants:
  if(TMath::Abs(qc2)>0.)
  {
   fIntFlowQcumulantsVsM[0]->SetBinContent(b,qc2);
   fIntFlowQcumulantsVsM[0]->SetBinError(b,qc2Error);  
  }
  if(TMath::Abs(qc4)>0.)
  {
   fIntFlowQcumulantsVsM[1]->SetBinContent(b,qc4);  
   fIntFlowQcumulantsVsM[1]->SetBinError(b,qc4Error);
  }
  if(TMath::Abs(qc6)>0.)
  {
   fIntFlowQcumulantsVsM[2]->SetBinContent(b,qc6); 
   fIntFlowQcumulantsVsM[2]->SetBinError(b,qc6Error);
  }
  if(TMath::Abs(qc8)>0.)
  {  
   fIntFlowQcumulantsVsM[3]->SetBinContent(b,qc8);
   fIntFlowQcumulantsVsM[3]->SetBinError(b,qc8Error);
  } 
  // Rebin in M:
  for(Int_t co=0;co<4;co++)
  {
   if(fIntFlowCorrelationsVsMPro[co]->GetBinEffectiveEntries(b)<2){continue;}
   value[co] = fIntFlowQcumulantsVsM[co]->GetBinContent(b);
   error[co] = fIntFlowQcumulantsVsM[co]->GetBinError(b);
   if(error[co]>0.)
   {
    dSum1[co]+=value[co]/(error[co]*error[co]);
    dSum2[co]+=1./(error[co]*error[co]);
   }
  } // end of for(Int_t co=0;co<4;co++) 
 } // end of for(Int_t b=1;b<=nBins;b++)
 // Store rebinned Q-cumulants:
 for(Int_t co=0;co<4;co++)
 {
  if(dSum2[co]>0.)
  {
   fIntFlowQcumulantsRebinnedInM->SetBinContent(co+1,dSum1[co]/dSum2[co]);
   fIntFlowQcumulantsRebinnedInM->SetBinError(co+1,pow(1./dSum2[co],0.5));
  }
 } // end of for(Int_t co=0;co<4;co++)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateCumulantsIntFlow()

//================================================================================================================================ 

void AliFlowAnalysisWithQCumulants::CalculateReferenceFlow()
{
 // a) Calculate the final results for reference flow estimates from Q-cumulants;
 // b) Propagate the statistical errors to reference flow estimates from statistical error of Q-cumulants; 
 // c) Store the results and statistical errors of reference flow estimates in histogram fIntFlow.
 //    Binning of fIntFlow is organized as follows:
 //
 //            1st bin: v{2,QC}
 //            2nd bin: v{4,QC}
 //            3rd bin: v{6,QC}
 //            4th bin: v{8,QC}
 //
 
 // Reference flow estimates:
 Double_t v2 = 0.; // v{2,QC}  
 Double_t v4 = 0.; // v{4,QC}  
 Double_t v6 = 0.; // v{6,QC}  
 Double_t v8 = 0.; // v{8,QC}
 // Reference flow's statistical errors:
 Double_t v2Error = 0.; // v{2,QC} stat. error 
 Double_t v4Error = 0.; // v{4,QC} stat. error
 Double_t v6Error = 0.; // v{6,QC} stat. error
 Double_t v8Error = 0.; // v{8,QC} stat. error
  
 // Q-cumulants:
 Double_t qc2 = fIntFlowQcumulants->GetBinContent(1); // QC{2}  
 Double_t qc4 = fIntFlowQcumulants->GetBinContent(2); // QC{4}  
 Double_t qc6 = fIntFlowQcumulants->GetBinContent(3); // QC{6}  
 Double_t qc8 = fIntFlowQcumulants->GetBinContent(4); // QC{8}
 // Q-cumulants's statistical errors: 
 Double_t qc2Error = fIntFlowQcumulants->GetBinError(1); // QC{2} stat. error  
 Double_t qc4Error = fIntFlowQcumulants->GetBinError(2); // QC{4} stat. error  
 Double_t qc6Error = fIntFlowQcumulants->GetBinError(3); // QC{6} stat. error  
 Double_t qc8Error = fIntFlowQcumulants->GetBinError(4); // QC{8} stat. error
 // Calculate reference flow estimates from Q-cumulants: 
 if(qc2>=0.){v2 = pow(qc2,0.5);} 
 if(qc4<=0.){v4 = pow(-1.*qc4,1./4.);} 
 if(qc6>=0.){v6 = pow((1./4.)*qc6,1./6.);}
 if(qc8<=0.){v8 = pow((-1./33.)*qc8,1./8.);}  
 // Calculate stat. error for reference flow estimates from stat. error of Q-cumulants:  
 if(qc2>0.){v2Error = (1./2.)*pow(qc2,-0.5)*qc2Error;} 
 if(qc4<0.){v4Error = (1./4.)*pow(-qc4,-3./4.)*qc4Error;} 
 if(qc6>0.){v6Error = (1./6.)*pow(2.,-1./3.)*pow(qc6,-5./6.)*qc6Error;}   
 if(qc8<0.){v8Error = (1./8.)*pow(33.,-1./8.)*pow(-qc8,-7./8.)*qc8Error;}   
 // Print warnings for the 'wrong sign' cumulants: 
 if(TMath::Abs(v2) < 1.e-44)
 {
  cout<<" WARNING: Wrong sign QC{2}, couldn't calculate v{2,QC} !!!!"<<endl;
 }
 if(TMath::Abs(v4) < 1.e-44)
 {
  cout<<" WARNING: Wrong sign QC{4}, couldn't calculate v{4,QC} !!!!"<<endl;
 } 
 if(TMath::Abs(v6) < 1.e-44)
 {
  cout<<" WARNING: Wrong sign QC{6}, couldn't calculate v{6,QC} !!!!"<<endl; 
 }
 if(TMath::Abs(v8) < 1.e-44)
 {
  cout<<" WARNING: Wrong sign QC{8}, couldn't calculate v{8,QC} !!!!"<<endl;
 }                       
 // Store the results and statistical errors of integrated flow estimates:
 fIntFlow->SetBinContent(1,v2);
 fIntFlow->SetBinError(1,v2Error);
 fIntFlow->SetBinContent(2,v4);
 fIntFlow->SetBinError(2,v4Error);
 fIntFlow->SetBinContent(3,v6);
 fIntFlow->SetBinError(3,v6Error);
 fIntFlow->SetBinContent(4,v8);
 fIntFlow->SetBinError(4,v8Error);  
  
 // Versus multiplicity: 
 if(!fCalculateCumulantsVsM){return;} 
 Int_t nBins = fIntFlowCorrelationsVsMPro[0]->GetNbinsX(); // to be improved (hardwired 0) 
 for(Int_t b=1;b<=nBins;b++)
 {
  // Q-cumulants:
  Double_t qc2VsM = fIntFlowQcumulantsVsM[0]->GetBinContent(b); // QC{2}  
  Double_t qc4VsM = fIntFlowQcumulantsVsM[1]->GetBinContent(b); // QC{4}  
  Double_t qc6VsM = fIntFlowQcumulantsVsM[2]->GetBinContent(b); // QC{6}  
  Double_t qc8VsM = fIntFlowQcumulantsVsM[3]->GetBinContent(b); // QC{8}
  // Q-cumulants's statistical errors: 
  Double_t qc2ErrorVsM = fIntFlowQcumulantsVsM[0]->GetBinError(b); // QC{2} stat. error  
  Double_t qc4ErrorVsM = fIntFlowQcumulantsVsM[1]->GetBinError(b); // QC{4} stat. error  
  Double_t qc6ErrorVsM = fIntFlowQcumulantsVsM[2]->GetBinError(b); // QC{6} stat. error  
  Double_t qc8ErrorVsM = fIntFlowQcumulantsVsM[3]->GetBinError(b); // QC{8} stat. error
  // Reference flow estimates:
  Double_t v2VsM = 0.; // v{2,QC}  
  Double_t v4VsM = 0.; // v{4,QC}  
  Double_t v6VsM = 0.; // v{6,QC}  
  Double_t v8VsM = 0.; // v{8,QC}
  // Reference flow estimates errors:
  Double_t v2ErrorVsM = 0.; // v{2,QC} stat. error  
  Double_t v4ErrorVsM = 0.; // v{4,QC} stat. error
  Double_t v6ErrorVsM = 0.; // v{6,QC} stat. error  
  Double_t v8ErrorVsM = 0.; // v{8,QC} stat. error
  // Calculate reference flow estimates from Q-cumulants: 
  if(qc2VsM>=0.){v2VsM = pow(qc2VsM,0.5);} 
  if(qc4VsM<=0.){v4VsM = pow(-1.*qc4VsM,1./4.);} 
  if(qc6VsM>=0.){v6VsM = pow((1./4.)*qc6VsM,1./6.);}
  if(qc8VsM<=0.){v8VsM = pow((-1./33.)*qc8VsM,1./8.);}  
  // Calculate stat. error for reference flow estimates from stat. error of Q-cumulants: 
  if(qc2VsM>0.){v2ErrorVsM = (1./2.)*pow(qc2VsM,-0.5)*qc2ErrorVsM;} 
  if(qc4VsM<0.){v4ErrorVsM = (1./4.)*pow(-qc4VsM,-3./4.)*qc4ErrorVsM;} 
  if(qc6VsM>0.){v6ErrorVsM = (1./6.)*pow(2.,-1./3.)*pow(qc6VsM,-5./6.)*qc6ErrorVsM;}   
  if(qc8VsM<0.){v8ErrorVsM = (1./8.)*pow(33.,-1./8.)*pow(-qc8VsM,-7./8.)*qc8ErrorVsM;}                       
  // Store the results and statistical errors of integrated flow estimates:
  fIntFlowVsM[0]->SetBinContent(b,v2VsM);
  fIntFlowVsM[0]->SetBinError(b,v2ErrorVsM);
  fIntFlowVsM[1]->SetBinContent(b,v4VsM);
  fIntFlowVsM[1]->SetBinError(b,v4ErrorVsM);
  fIntFlowVsM[2]->SetBinContent(b,v6VsM);
  fIntFlowVsM[2]->SetBinError(b,v6ErrorVsM);
  fIntFlowVsM[3]->SetBinContent(b,v8VsM);
  fIntFlowVsM[3]->SetBinError(b,v8ErrorVsM);
 } // end of for(Int_t b=1;b<=nBins;b++)
 
 // 'Rebinned in M' calculation: // to be improved - this can be implemented better:   
 // Reference flow estimates:
 Double_t v2RebinnedInM = 0.; // v{2,QC}  
 Double_t v4RebinnedInM = 0.; // v{4,QC}  
 Double_t v6RebinnedInM = 0.; // v{6,QC}  
 Double_t v8RebinnedInM = 0.; // v{8,QC}
 // Reference flow's statistical errors:
 Double_t v2ErrorRebinnedInM = 0.; // v{2,QC} stat. error 
 Double_t v4ErrorRebinnedInM = 0.; // v{4,QC} stat. error
 Double_t v6ErrorRebinnedInM = 0.; // v{6,QC} stat. error
 Double_t v8ErrorRebinnedInM = 0.; // v{8,QC} stat. error
 // Q-cumulants:
 Double_t qc2RebinnedInM = fIntFlowQcumulantsRebinnedInM->GetBinContent(1); // QC{2}  
 Double_t qc4RebinnedInM = fIntFlowQcumulantsRebinnedInM->GetBinContent(2); // QC{4}  
 Double_t qc6RebinnedInM = fIntFlowQcumulantsRebinnedInM->GetBinContent(3); // QC{6}  
 Double_t qc8RebinnedInM = fIntFlowQcumulantsRebinnedInM->GetBinContent(4); // QC{8}
 // Q-cumulants's statistical errors: 
 Double_t qc2ErrorRebinnedInM = fIntFlowQcumulantsRebinnedInM->GetBinError(1); // QC{2} stat. error  
 Double_t qc4ErrorRebinnedInM = fIntFlowQcumulantsRebinnedInM->GetBinError(2); // QC{4} stat. error  
 Double_t qc6ErrorRebinnedInM = fIntFlowQcumulantsRebinnedInM->GetBinError(3); // QC{6} stat. error  
 Double_t qc8ErrorRebinnedInM = fIntFlowQcumulantsRebinnedInM->GetBinError(4); // QC{8} stat. error
 // Calculate reference flow estimates from Q-cumulants: 
 if(qc2RebinnedInM>=0.){v2RebinnedInM = pow(qc2RebinnedInM,0.5);} 
 if(qc4RebinnedInM<=0.){v4RebinnedInM = pow(-1.*qc4RebinnedInM,1./4.);} 
 if(qc6RebinnedInM>=0.){v6RebinnedInM = pow((1./4.)*qc6RebinnedInM,1./6.);}
 if(qc8RebinnedInM<=0.){v8RebinnedInM = pow((-1./33.)*qc8RebinnedInM,1./8.);}  
 // Calculate stat. error for reference flow estimates from stat. error of Q-cumulants: 
 if(qc2RebinnedInM>0.){v2ErrorRebinnedInM = (1./2.)*pow(qc2RebinnedInM,-0.5)*qc2ErrorRebinnedInM;} 
 if(qc4RebinnedInM<0.){v4ErrorRebinnedInM = (1./4.)*pow(-qc4RebinnedInM,-3./4.)*qc4ErrorRebinnedInM;} 
 if(qc6RebinnedInM>0.){v6ErrorRebinnedInM = (1./6.)*pow(2.,-1./3.)*pow(qc6RebinnedInM,-5./6.)*qc6ErrorRebinnedInM;}   
 if(qc8RebinnedInM<0.){v8ErrorRebinnedInM = (1./8.)*pow(33.,-1./8.)*pow(-qc8RebinnedInM,-7./8.)*qc8ErrorRebinnedInM;}   
 // Print warnings for the 'wrong sign' cumulants: 
 if(TMath::Abs(v2RebinnedInM) < 1.e-44)
 {
  cout<<" WARNING: Wrong sign QC{2} rebinned in M, couldn't calculate v{2,QC} !!!!"<<endl;
 }
 if(TMath::Abs(v4RebinnedInM) < 1.e-44)
 {
  cout<<" WARNING: Wrong sign QC{4} rebinned in M, couldn't calculate v{4,QC} !!!!"<<endl;
 }
 if(TMath::Abs(v6RebinnedInM) < 1.e-44)
 {
  cout<<" WARNING: Wrong sign QC{6} rebinned in M, couldn't calculate v{6,QC} !!!!"<<endl;
 }
 if(TMath::Abs(v8RebinnedInM) < 1.e-44)
 {
  cout<<" WARNING: Wrong sign QC{8} rebinned in M, couldn't calculate v{8,QC} !!!!"<<endl;
 }                       
 // Store the results and statistical errors of integrated flow estimates:
 fIntFlowRebinnedInM->SetBinContent(1,v2RebinnedInM);
 fIntFlowRebinnedInM->SetBinError(1,v2ErrorRebinnedInM);
 fIntFlowRebinnedInM->SetBinContent(2,v4RebinnedInM);
 fIntFlowRebinnedInM->SetBinError(2,v4ErrorRebinnedInM);
 fIntFlowRebinnedInM->SetBinContent(3,v6RebinnedInM);
 fIntFlowRebinnedInM->SetBinError(3,v6ErrorRebinnedInM);
 fIntFlowRebinnedInM->SetBinContent(4,v8RebinnedInM);
 fIntFlowRebinnedInM->SetBinError(4,v8ErrorRebinnedInM);    
  
} // end of AliFlowAnalysisWithQCumulants::CalculateReferenceFlow()

//================================================================================================================================ 

void AliFlowAnalysisWithQCumulants::FillCommonHistResultsIntFlow()
{
 // Fill in AliFlowCommonHistResults histograms relevant for reference flow.
 
 // There are two possibilities here:
 // a) Store minimum bias reference flow - use SetMinimumBiasReferenceFlow(kTRUE). This result is 
 //    biased by the interplay between nonflow correlations and multiplicity fluctuations and is 
 //    also stored in local histogram fIntFlow; 
 // b) Store reference flow obtained from flow analysis performed at fixed multiplicity and 
 //    rebinned only at the end of the day - use SetMinimumBiasReferenceFlow(kFALSE). This result
 //    is also stored in local histogram fIntFlowRebinnedInM.
 
 // Reference flow estimates:
 Double_t v[4] = {0.};
 // Statistical errors of reference flow estimates:
 Double_t vError[4] = {0.};
  
 for(Int_t b=0;b<4;b++)
 {
  if(fMinimumBiasReferenceFlow)
  { 
   v[b] = fIntFlow->GetBinContent(b+1);
   vError[b] = fIntFlow->GetBinError(b+1);
  } else
    {
     v[b] = fIntFlowRebinnedInM->GetBinContent(b+1);
     vError[b] = fIntFlowRebinnedInM->GetBinError(b+1);
    }
 } // end of for(Int_t b=0;b<4;b++)
  
 // Fill AliFlowCommonHistResults histogram:
 fCommonHistsResults2nd->FillIntegratedFlow(v[0],vError[0]); // to be improved (hardwired 2nd in the name)  
 fCommonHistsResults4th->FillIntegratedFlow(v[1],vError[1]); // to be improved (hardwired 4th in the name)
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)) // to be improved (calculate also 6th and 8th order)
 {
  fCommonHistsResults6th->FillIntegratedFlow(v[2],vError[2]); // to be improved (hardwired 6th in the name)
  fCommonHistsResults8th->FillIntegratedFlow(v[3],vError[3]); // to be improved (hardwired 8th in the name) 
 }
 
} // end of AliFlowAnalysisWithQCumulants::FillCommonHistResultsIntFlow()

//================================================================================================================================ 

void AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrelationsUsingParticleWeights()
{
 // Calculate all correlations needed for integrated flow using particle weights.
  
 // Remark 1: When particle weights are used the binning of fIntFlowCorrelationAllPro is organized as follows:
 //
 //  1st bin: <2>_{1n|1n} = two1n1nW1W1 = <w1 w2 cos(n*(phi1-phi2))>
 //  2nd bin: <2>_{2n|2n} = two2n2nW2W2 = <w1^2 w2^2 cos(2n*(phi1-phi2))>
 //  3rd bin: <2>_{3n|3n} = two3n3nW3W3 = <w1^3 w2^3 cos(3n*(phi1-phi2))> 
 //  4th bin: <2>_{4n|4n} = two4n4nW4W4 = <w1^4 w2^4 cos(4n*(phi1-phi2))>
 //  5th bin:           ----  EMPTY ----
 //  6th bin: <3>_{2n|1n,1n} = three2n1n1nW2W1W1 = <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))>
 //  7th bin: <3>_{3n|2n,1n} = ...
 //  8th bin: <3>_{4n|2n,2n} = ...
 //  9th bin: <3>_{4n|3n,1n} = ...
 // 10th bin:           ----  EMPTY ----
 // 11th bin: <4>_{1n,1n|1n,1n} = four1n1n1n1nW1W1W1W1 = <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
 // 12th bin: <4>_{2n,1n|2n,1n} = ...
 // 13th bin: <4>_{2n,2n|2n,2n} = ...
 // 14th bin: <4>_{3n|1n,1n,1n} = ... 
 // 15th bin: <4>_{3n,1n|3n,1n} = ...
 // 16th bin: <4>_{3n,1n|2n,2n} = ...
 // 17th bin: <4>_{4n|2n,1n,1n} = ... 
 // 18th bin:           ----  EMPTY ----
 // 19th bin: <5>_{2n|1n,1n,1n,1n} = ...
 // 20th bin: <5>_{2n,2n|2n,1n,1n} = ...
 // 21st bin: <5>_{3n,1n|2n,1n,1n} = ...
 // 22nd bin: <5>_{4n|1n,1n,1n,1n} = ...
 // 23rd bin:           ----  EMPTY ----
 // 24th bin: <6>_{1n,1n,1n|1n,1n,1n} = ...
 // 25th bin: <6>_{2n,1n,1n|2n,1n,1n} = ...
 // 26th bin: <6>_{2n,2n|1n,1n,1n,1n} = ...
 // 27th bin: <6>_{3n,1n|1n,1n,1n,1n} = ...
 // 28th bin:           ----  EMPTY ----
 // 29th bin: <7>_{2n,1n,1n|1n,1n,1n,1n} = ...
 // 30th bin:           ----  EMPTY ----
 // 31st bin: <8>_{1n,1n,1n,1n|1n,1n,1n,1n} = ...
 
 // Remark 2: When particle weights are used there are some extra correlations. They are stored in 
 // fIntFlowExtraCorrelationsPro binning of which is organized as follows:
 
 // 1st bin: two1n1nW3W1 = <w1^3 w2 cos(n*(phi1-phi2))>
 // 2nd bin: two1n1nW1W1W2 = <w1 w2 w3^2 cos(n*(phi1-phi2))>  
 
 // multiplicity (number of particles used to determine the reaction plane)
 Double_t dMult = (*fSpk)(0,0);
 
 // real and imaginary parts of weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n1k = (*fReQ)(0,1);
 Double_t dReQ2n2k = (*fReQ)(1,2);
 Double_t dReQ3n3k = (*fReQ)(2,3);
 Double_t dReQ4n4k = (*fReQ)(3,4);
 Double_t dReQ1n3k = (*fReQ)(0,3);
 Double_t dImQ1n1k = (*fImQ)(0,1);
 Double_t dImQ2n2k = (*fImQ)(1,2);
 Double_t dImQ3n3k = (*fImQ)(2,3);
 Double_t dImQ4n4k = (*fImQ)(3,4);
 Double_t dImQ1n3k = (*fImQ)(0,3);

 // dMs are variables introduced in order to simplify some Eqs. bellow:
 //..............................................................................................
 Double_t dM11 = (*fSpk)(1,1)-(*fSpk)(0,2); // dM11 = sum_{i,j=1,i!=j}^M w_i w_j
 Double_t dM22 = (*fSpk)(1,2)-(*fSpk)(0,4); // dM22 = sum_{i,j=1,i!=j}^M w_i^2 w_j^2
 Double_t dM33 = (*fSpk)(1,3)-(*fSpk)(0,6); // dM33 = sum_{i,j=1,i!=j}^M w_i^3 w_j^3
 Double_t dM44 = (*fSpk)(1,4)-(*fSpk)(0,8); // dM44 = sum_{i,j=1,i!=j}^M w_i^4 w_j^4
 Double_t dM31 = (*fSpk)(0,3)*(*fSpk)(0,1)-(*fSpk)(0,4); // dM31 = sum_{i,j=1,i!=j}^M w_i^3 w_j
 Double_t dM211 = (*fSpk)(0,2)*(*fSpk)(1,1)-2.*(*fSpk)(0,3)*(*fSpk)(0,1)
                - (*fSpk)(1,2)+2.*(*fSpk)(0,4); // dM211 = sum_{i,j,k=1,i!=j!=k}^M w_i^2 w_j w_k
 Double_t dM1111 = (*fSpk)(3,1)-6.*(*fSpk)(0,2)*(*fSpk)(1,1)  
                 + 8.*(*fSpk)(0,3)*(*fSpk)(0,1)
                 + 3.*(*fSpk)(1,2)-6.*(*fSpk)(0,4); // dM1111 = sum_{i,j,k,l=1,i!=j!=k!=l}^M w_i w_j w_k w_l
 //..............................................................................................

 // 2-particle correlations:
 Double_t two1n1nW1W1 = 0.; // <w1 w2 cos(n*(phi1-phi2))>
 Double_t two2n2nW2W2 = 0.; // <w1^2 w2^2 cos(2n*(phi1-phi2))>
 Double_t two3n3nW3W3 = 0.; // <w1^3 w2^3 cos(3n*(phi1-phi2))>
 Double_t two4n4nW4W4 = 0.; // <w1^4 w2^4 cos(4n*(phi1-phi2))>
 if(dMult>1) 
 { 
  if(dM11)
  {
   two1n1nW1W1 = (pow(dReQ1n1k,2)+pow(dImQ1n1k,2)-(*fSpk)(0,2))/dM11;    
   // average correlation <w1 w2 cos(n*(phi1-phi2))> for single event: 
   fIntFlowCorrelationsEBE->SetBinContent(1,two1n1nW1W1);
   fIntFlowEventWeightsForCorrelationsEBE->SetBinContent(1,dM11);
   // average correlation <w1 w2 cos(n*(phi1-phi2))> for all events:
   fIntFlowCorrelationsPro->Fill(0.5,two1n1nW1W1,dM11);  
   // average squared correlation <w1 w2 cos(n*(phi1-phi2))> for all events:
   fIntFlowSquaredCorrelationsPro->Fill(0.5,two1n1nW1W1*two1n1nW1W1,dM11); 
   fIntFlowCorrelationsAllPro->Fill(0.5,two1n1nW1W1,dM11);   
  }
  if(dM22)
  {
   two2n2nW2W2 = (pow(dReQ2n2k,2)+pow(dImQ2n2k,2)-(*fSpk)(0,4))/dM22; 
   // ...
   // average correlation <w1^2 w2^2 cos(2n*(phi1-phi2))> for all events:
   fIntFlowCorrelationsAllPro->Fill(1.5,two2n2nW2W2,dM22);   
  }
  if(dM33)
  {
   two3n3nW3W3 = (pow(dReQ3n3k,2)+pow(dImQ3n3k,2)-(*fSpk)(0,6))/dM33;
   // ...
   // average correlation <w1^3 w2^3 cos(3n*(phi1-phi2))> for all events:
   fIntFlowCorrelationsAllPro->Fill(2.5,two3n3nW3W3,dM33);   
  }
  if(dM44)
  {
   two4n4nW4W4 = (pow(dReQ4n4k,2)+pow(dImQ4n4k,2)-(*fSpk)(0,8))/dM44; 
   // ...
   // average correlation <w1^4 w2^4 cos(4n*(phi1-phi2))> for all events:
   fIntFlowCorrelationsAllPro->Fill(3.5,two4n4nW4W4,dM44);      
  }
 } // end of if(dMult>1) 

 // extra 2-particle correlations:
 Double_t two1n1nW3W1 = 0.; // <w1^3 w2 cos(n*(phi1-phi2))>
 Double_t two1n1nW1W1W2 = 0.; // <w1 w2 w3^2 cos(n*(phi1-phi2))> 
 if(dMult>1) 
 {    
  if(dM31)
  {
   two1n1nW3W1 = (dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k-(*fSpk)(0,4))/dM31; 
   fIntFlowExtraCorrelationsPro->Fill(0.5,two1n1nW3W1,dM31);  
  } 
  if(dM211)
  {
   two1n1nW1W1W2 = ((*fSpk)(0,2)*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2)-(*fSpk)(0,2))
                 - 2.*(dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k
                 - (*fSpk)(0,4)))/dM211;
   fIntFlowExtraCorrelationsPro->Fill(1.5,two1n1nW1W1W2,dM211);  
  }  
 } // end of if(dMult>1)
 //..............................................................................................
 
 //..............................................................................................
 // 3-particle correlations:
 Double_t three2n1n1nW2W1W1 = 0.; // <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))>
 
 if(dMult>2) 
 { 
  if(dM211)
  {                                                       
   three2n1n1nW2W1W1 = (pow(dReQ1n1k,2.)*dReQ2n2k+2.*dReQ1n1k*dImQ1n1k*dImQ2n2k-pow(dImQ1n1k,2.)*dReQ2n2k
                     - 2.*(dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k)
                     - pow(dReQ2n2k,2)-pow(dImQ2n2k,2)
                     + 2.*(*fSpk)(0,4))/dM211;                                                                               
   fIntFlowCorrelationsAllPro->Fill(5.5,three2n1n1nW2W1W1,dM211);
  } 
 } // end of if(dMult>2) 
 //..............................................................................................
 
 //..............................................................................................
 // 4-particle correlations:
 Double_t four1n1n1n1nW1W1W1W1 = 0.; // <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
 if(dMult>3) 
 { 
  if(dM1111)
  {      
   four1n1n1n1nW1W1W1W1 = (pow(pow(dReQ1n1k,2.)+pow(dImQ1n1k,2.),2)
                        - 2.*(pow(dReQ1n1k,2.)*dReQ2n2k+2.*dReQ1n1k*dImQ1n1k*dImQ2n2k-pow(dImQ1n1k,2.)*dReQ2n2k)
                        + 8.*(dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k)
                        + (pow(dReQ2n2k,2)+pow(dImQ2n2k,2))
                        - 4.*(*fSpk)(0,2)*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2))
                        - 6.*(*fSpk)(0,4)+2.*(*fSpk)(1,2))/dM1111;  
                          
   // average correlation <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> for single event: 
   fIntFlowCorrelationsEBE->SetBinContent(2,four1n1n1n1nW1W1W1W1);
   fIntFlowEventWeightsForCorrelationsEBE->SetBinContent(2,dM1111);
   // average correlation <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> for all events:
   fIntFlowCorrelationsPro->Fill(1.5,four1n1n1n1nW1W1W1W1,dM1111);   
   // average squared correlation <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> for all events:
   fIntFlowSquaredCorrelationsPro->Fill(1.5,four1n1n1n1nW1W1W1W1*four1n1n1n1nW1W1W1W1,dM1111);      
   fIntFlowCorrelationsAllPro->Fill(10.5,four1n1n1n1nW1W1W1W1,dM1111);   
  } 
 } // end of if(dMult>3) 
 //..............................................................................................
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrelationsUsingParticleWeights()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::InitializeArraysForIntFlow()
{
 // Initialize all arrays used to calculate integrated flow.
 
 for(Int_t sc=0;sc<2;sc++) // sin or cos terms
 {
  fIntFlowCorrectionTermsForNUAEBE[sc] = NULL;
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[sc] = NULL;
  fIntFlowCorrectionTermsForNUAPro[sc] = NULL;
  fIntFlowCorrectionTermsForNUAHist[sc] = NULL;
  for(Int_t ci=0;ci<4;ci++) // correction term index (to be improved - hardwired 4)
  {
   fIntFlowCorrectionTermsForNUAVsMPro[sc][ci] = NULL;
  }
  for(Int_t power=0;power<2;power++) // linear or quadratic 
  {
   fIntFlowSumOfEventWeightsNUA[sc][power] = NULL;
  }
 }
 for(Int_t power=0;power<2;power++) // linear or quadratic 
 {
  fIntFlowSumOfEventWeights[power] = NULL;    
 }
 for(Int_t i=0;i<4;i++) // print on the screen the final results (0=RF, 1=RP, 2=POI, 3=RF (rebbined in M))
 {
  fPrintFinalResults[i] = kTRUE;
 }
 for(Int_t ci=0;ci<4;ci++) // correlation index or cumulant order
 {
  fIntFlowCorrelationsVsMPro[ci] = NULL;
  fIntFlowSquaredCorrelationsVsMPro[ci] = NULL;
  fIntFlowCorrelationsVsMHist[ci] = NULL;
  fIntFlowQcumulantsVsM[ci] = NULL;
  fIntFlowVsM[ci] = NULL;
  fIntFlowDetectorBiasVsM[ci] = NULL;
  for(Int_t lc=0;lc<2;lc++)
  {
   fIntFlowSumOfEventWeightsVsM[ci][lc] = NULL;
  }
 } 
 for(Int_t pi=0;pi<6;pi++) // product or covariance index
 {
  fIntFlowProductOfCorrelationsVsMPro[pi] = NULL;
  fIntFlowCovariancesVsM[pi] = NULL;
  fIntFlowSumOfProductOfEventWeightsVsM[pi] = NULL;
 } 
 for(Int_t ci=0;ci<64;ci++) // correlation index for all correlations vs M profiles (to be improved - hardwired 64)
 {
  fIntFlowCorrelationsAllVsMPro[ci] = NULL;
 } 

} // end of void AliFlowAnalysisWithQCumulants::InitializeArraysForIntFlow()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::InitializeArraysForDiffFlow()
{
 // Initialize all arrays needed to calculate differential flow.
 //  a) Initialize lists holding profiles;
 //  b) Initialize lists holding histograms;
 //  c) Initialize event-by-event quantities;
 //  d) Initialize profiles;
 //  e) Initialize histograms holding final results.
 
 // a) Initialize lists holding profiles;
 for(Int_t t=0;t<2;t++) // type (RP, POI)
 {
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  {
   fDiffFlowCorrelationsProList[t][pe] = NULL;
   fDiffFlowProductOfCorrelationsProList[t][pe] = NULL;
   fDiffFlowCorrectionsProList[t][pe] = NULL;
  }
  // 2D:
  f2DDiffFlowCorrelationsProList[t] = NULL;
 }  
 
 // b) Initialize lists holding histograms;
 for(Int_t t=0;t<2;t++) // type (RP, POI)
 {
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  {
   fDiffFlowCorrelationsHistList[t][pe] = NULL;
   for(Int_t power=0;power<2;power++)
   {
    fDiffFlowSumOfEventWeightsHistList[t][pe][power] = NULL;
   } // end of for(Int_t power=0;power<2;power++)  
   fDiffFlowSumOfProductOfEventWeightsHistList[t][pe] = NULL;
   fDiffFlowCorrectionsHistList[t][pe] = NULL;
   fDiffFlowCovariancesHistList[t][pe] = NULL;
   fDiffFlowCumulantsHistList[t][pe] = NULL;
   fDiffFlowDetectorBiasHistList[t][pe] = NULL;
   fDiffFlowHistList[t][pe] = NULL;
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
 } // enf of for(Int_t t=0;t<2;t++) // type (RP, POI) 
 
 // c) Initialize event-by-event quantities:
 // 1D:
 for(Int_t t=0;t<3;t++) // type (RP, POI, POI&&RP)
 {
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  { 
   for(Int_t m=0;m<4;m++) // multiple of harmonic
   {
    for(Int_t k=0;k<9;k++) // power of weight
    {
     fReRPQ1dEBE[t][pe][m][k] = NULL;
     fImRPQ1dEBE[t][pe][m][k] = NULL;
     fs1dEBE[t][pe][k] = NULL; // to be improved (this doesn't need to be within loop over m)
    }   
   }
  }
 }
 // 1D:
 for(Int_t t=0;t<2;t++) // type (RP or POI)
 {
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  { 
   for(Int_t sc=0;sc<2;sc++) // sin or cos terms
   {
    for(Int_t cti=0;cti<9;cti++) // correction term index
    {
     fDiffFlowCorrectionTermsForNUAEBE[t][pe][sc][cti] = NULL;
    }   
   }
  }
 }
 // 2D:  
 for(Int_t t=0;t<3;t++) // type (RP, POI, POI&&RP)
 {
  for(Int_t m=0;m<4;m++) // multiple of harmonic
  {
   for(Int_t k=0;k<9;k++) // power of weight
   {
    fReRPQ2dEBE[t][m][k] = NULL;
    fImRPQ2dEBE[t][m][k] = NULL;
    fs2dEBE[t][k] = NULL; // to be improved (this doesn't need to be within loop over m)
   }   
  }
 }
 
 // d) Initialize profiles:
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  {
   for(Int_t ci=0;ci<4;ci++) // correlation index
   {
    fDiffFlowCorrelationsPro[t][pe][ci] = NULL;
    fDiffFlowSquaredCorrelationsPro[t][pe][ci] = NULL;
   } // end of for(Int_t ci=0;ci<4;ci++)   
   for(Int_t mci1=0;mci1<8;mci1++) // mixed correlation index
   {
    for(Int_t mci2=0;mci2<8;mci2++) // mixed correlation index
    {
     fDiffFlowProductOfCorrelationsPro[t][pe][mci1][mci2] = NULL;
    } // end of for(Int_t mci2=0;mci2<8;mci2++) // mixed correlation index
   } // end of for(Int_t mci1=0;mci1<8;mci1++) // mixed correlation index  
   // correction terms for nua:
   for(Int_t sc=0;sc<2;sc++) // sin or cos terms
   {
    for(Int_t cti=0;cti<9;cti++) // correction term index
    {
     fDiffFlowCorrectionTermsForNUAPro[t][pe][sc][cti] = NULL;
    }   
   }
   // other differential correlators:
   for(Int_t sc=0;sc<2;sc++) // sin or cos terms
   {
    for(Int_t ci=0;ci<1;ci++) // correction term index
    {
     fOtherDiffCorrelators[t][pe][sc][ci] = NULL;
    }   
   }
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
  for(Int_t ci=0;ci<4;ci++) // correlation index
  {
   f2DDiffFlowCorrelationsPro[t][ci] = NULL;
  }
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
  
 // e) Initialize histograms holding final results.
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  {
   for(Int_t ci=0;ci<4;ci++) // correlation index
   {
    fDiffFlowCorrelationsHist[t][pe][ci] = NULL;
    fDiffFlowCumulants[t][pe][ci] = NULL;
    fDiffFlowDetectorBias[t][pe][ci] = NULL;
    fDiffFlow[t][pe][ci] = NULL;
   } // end of for(Int_t ci=0;ci<4;ci++)    
   for(Int_t covarianceIndex=0;covarianceIndex<5;covarianceIndex++) 
   {
    fDiffFlowCovariances[t][pe][covarianceIndex] = NULL;     
   } // end of for(Int_t covarianceIndex=0;covarianceIndex<5;covarianceIndex++) 
   // correction terms for nua:
   for(Int_t sc=0;sc<2;sc++) // sin or cos terms
   {
    for(Int_t cti=0;cti<9;cti++) // correction term index
    {
     fDiffFlowCorrectionTermsForNUAHist[t][pe][sc][cti] = NULL;
    }   
   }
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
  for(Int_t ci=0;ci<4;ci++) // correlation index
  {
   f2DDiffFlowCumulants[t][ci] = NULL;
   f2DDiffFlow[t][ci] = NULL;
  }
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
 
 // sum of event weights for reduced correlations:
 for(Int_t t=0;t<2;t++) // type = RP or POI
 {
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  {
   for(Int_t p=0;p<2;p++) // power of weight is 1 or 2
   {
    for(Int_t ew=0;ew<4;ew++) // event weight index for reduced correlations
    {
     fDiffFlowSumOfEventWeights[t][pe][p][ew] = NULL;
    } 
   }   
  }
 }
 // product of event weights for both types of correlations:
 for(Int_t t=0;t<2;t++) // type = RP or POI
 {
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  {
   for(Int_t mci1=0;mci1<8;mci1++) // mixed correlation index
   {
    for(Int_t mci2=0;mci2<8;mci2++) // mixed correlation index
    {
     fDiffFlowSumOfProductOfEventWeights[t][pe][mci1][mci2] = NULL;
    } 
   }   
  }
 }
    
} // end of AliFlowAnalysisWithQCumulants::InitializeArraysForDiffFlow()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCumulants(TString type, TString ptOrEta)
{
 // Calculate differential flow cumulants from measured multiparticle correlations.
 
 // REMARK: Cumulants calculated in this method are NOT corrected for non-uniform acceptance. 
 // This correction, if enabled via setter SetApplyCorrectionForNUA(Bool_t), is applied 
 // in the method CalculateDiffFlowCumulantsCorrectedForNUA(TString type, TString ptOrEta)
 
 Int_t t = 0;
 Int_t pe = 0;

 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   } 
     
 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   } 
     
 // Common:
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 
 // Correlation <<2>>: 
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1);
 Double_t twoError = fIntFlowCorrelationsHist->GetBinError(1);
 
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // Reduced correlations:   
  Double_t twoPrime = fDiffFlowCorrelationsHist[t][pe][0]->GetBinContent(b); // <<2'>>
  Double_t twoPrimeError = fDiffFlowCorrelationsHist[t][pe][0]->GetBinError(b); // stat. error of <<2'>>
  Double_t fourPrime = fDiffFlowCorrelationsHist[t][pe][1]->GetBinContent(b); // <<4'>>
  Double_t fourPrimeError = fDiffFlowCorrelationsHist[t][pe][1]->GetBinError(b); // stat. error of <<4'>>
  // Covariances:
  Double_t wCovTwoTwoReduced = fDiffFlowCovariances[t][pe][0]->GetBinContent(b); // Cov(<2>,<2'>) * prefactor(<2>,<2'>)
  Double_t wCovTwoFourReduced = fDiffFlowCovariances[t][pe][1]->GetBinContent(b); // Cov(<2>,<4'>) * prefactor(<2>,<4'>)
  Double_t wCovTwoReducedFourReduced = fDiffFlowCovariances[t][pe][4]->GetBinContent(b); // Cov(<2'>,<4'>) * prefactor(<2'>,<4'>)
  // QC{2'}:
  Double_t qc2Prime = twoPrime; // QC{2'}
  Double_t qc2PrimeError = twoPrimeError; // stat. error of QC{2'}
  fDiffFlowCumulants[t][pe][0]->SetBinContent(b,qc2Prime); 
  fDiffFlowCumulants[t][pe][0]->SetBinError(b,qc2PrimeError); 
  // QC{4'}:
  Double_t qc4Prime = fourPrime - 2.*twoPrime*two; // QC{4'} = <<4'>> - 2*<<2'>><<2>>
  Double_t qc4PrimeError = 0.; // stat. error of QC{4'}
  Double_t qc4PrimeErrorSquared = 4.*pow(twoPrime,2.)*pow(twoError,2.)
                                + 4.*pow(two,2.)*pow(twoPrimeError,2.)
                                + pow(fourPrimeError,2.)
                                + 8.*two*twoPrime*wCovTwoTwoReduced
                                - 4.*twoPrime*wCovTwoFourReduced
                                - 4.*two*wCovTwoReducedFourReduced;
  if(qc4PrimeErrorSquared>0.)
  {
   qc4PrimeError = pow(qc4PrimeErrorSquared,0.5);
  } 
  fDiffFlowCumulants[t][pe][1]->SetBinContent(b,qc4Prime); 
  fDiffFlowCumulants[t][pe][1]->SetBinError(b,qc4PrimeError); 
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
    
} // end of void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCumulants(TString type, Bool_t useParticleWeights, TString eventWeights); 

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::Calculate2DDiffFlowCumulants(TString type)
{
 // Calculate 2D differential cumulants.
 
 // Remark: correction for detector effects and error propagation not implemented yet for 2D differential cumulants.
 
 Int_t t = 0;

 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }
       
 // Reference correlation <<2>>: 
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1);
 
 // Looping over all (pt,eta) bins and calculating differential flow cumulants: 
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   // Reduced correlations:   
   Double_t twoPrime = f2DDiffFlowCorrelationsPro[t][0]->GetBinContent(f2DDiffFlowCorrelationsPro[t][0]->GetBin(p,e)); // <<2'>>(pt,eta)
   Double_t fourPrime = f2DDiffFlowCorrelationsPro[t][1]->GetBinContent(f2DDiffFlowCorrelationsPro[t][1]->GetBin(p,e)); // <<4'>>(pt,eta)
   // Cumulants:
   Double_t qc2Prime = twoPrime; // QC{2'} = <<2'>>   
   f2DDiffFlowCumulants[t][0]->SetBinContent(f2DDiffFlowCumulants[t][0]->GetBin(p,e),qc2Prime); 
   Double_t qc4Prime = fourPrime - 2.*twoPrime*two; // QC{4'} = <<4'>> - 2*<<2'>><<2>>
   f2DDiffFlowCumulants[t][1]->SetBinContent(f2DDiffFlowCumulants[t][1]->GetBin(p,e),qc4Prime); 
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
} // end of void AliFlowAnalysisWithQCumulants::Calculate2DDiffFlowCumulants(TString type)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateFinalResultsForRPandPOIIntegratedFlow(TString type)
{
 // Calculate final results for integrated flow of RPs and POIs. 
  
 // to be improved - check if the integrated flow calculation here is actually correct 
  
 Int_t t = 0; // RP = 0, POI = 1

 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }
     
 // pt yield:    
 TH1F *yield2ndPt = NULL;
 TH1F *yield4thPt = NULL;
 TH1F *yield6thPt = NULL;
 TH1F *yield8thPt = NULL;
 
 if(type == "POI")
 {
  if(fFillMultipleControlHistograms)
  {
   yield2ndPt = (TH1F*)(fCommonHists2nd->GetHistPtPOI())->Clone();
   yield4thPt = (TH1F*)(fCommonHists4th->GetHistPtPOI())->Clone();
   yield6thPt = (TH1F*)(fCommonHists6th->GetHistPtPOI())->Clone();
   yield8thPt = (TH1F*)(fCommonHists8th->GetHistPtPOI())->Clone();  
  } else
    {
     yield2ndPt = (TH1F*)(fCommonHists->GetHistPtPOI())->Clone();
     yield4thPt = (TH1F*)(fCommonHists->GetHistPtPOI())->Clone();
     yield6thPt = (TH1F*)(fCommonHists->GetHistPtPOI())->Clone();
     yield8thPt = (TH1F*)(fCommonHists->GetHistPtPOI())->Clone();     
    }
 } 
 else if(type == "RP")
 {
  if(fFillMultipleControlHistograms)
  {
   yield2ndPt = (TH1F*)(fCommonHists2nd->GetHistPtRP())->Clone();
   yield4thPt = (TH1F*)(fCommonHists4th->GetHistPtRP())->Clone();
   yield6thPt = (TH1F*)(fCommonHists6th->GetHistPtRP())->Clone();
   yield8thPt = (TH1F*)(fCommonHists8th->GetHistPtRP())->Clone();  
  } else
    {
     yield2ndPt = (TH1F*)(fCommonHists->GetHistPtRP())->Clone();
     yield4thPt = (TH1F*)(fCommonHists->GetHistPtRP())->Clone();
     yield6thPt = (TH1F*)(fCommonHists->GetHistPtRP())->Clone();
     yield8thPt = (TH1F*)(fCommonHists->GetHistPtRP())->Clone();    
    } 
 } 
 
 if(!yield2ndPt){return;}
 if(!yield4thPt){return;}
 if(!yield6thPt){return;}
 if(!yield8thPt){return;} 

 Int_t nBinsPt = yield2ndPt->GetNbinsX();
 
 TH1D *flow2ndPt = NULL;
 TH1D *flow4thPt = NULL;
 TH1D *flow6thPt = NULL;
 TH1D *flow8thPt = NULL;
 
 // to be improved (hardwired pt index)
 flow2ndPt = (TH1D*)fDiffFlow[t][0][0]->Clone();
 flow4thPt = (TH1D*)fDiffFlow[t][0][1]->Clone();
 flow6thPt = (TH1D*)fDiffFlow[t][0][2]->Clone();
 flow8thPt = (TH1D*)fDiffFlow[t][0][3]->Clone(); 

 if(!flow2ndPt){return;}
 if(!flow4thPt){return;}
 if(!flow6thPt){return;}
 if(!flow8thPt){return;} 
   
 Double_t dvn2nd = 0., dvn4th = 0., dvn6th = 0., dvn8th = 0.; // differential flow
 Double_t dErrvn2nd = 0., dErrvn4th = 0., dErrvn6th = 0., dErrvn8th = 0.; // error on differential flow
 
 Double_t dVn2nd = 0., dVn4th = 0., dVn6th = 0., dVn8th = 0.; // integrated flow 
 Double_t dErrVn2nd = 0., dErrVn4th = 0., dErrVn6th = 0., dErrVn8th = 0.; // error on integrated flow

 Double_t dYield2nd = 0., dYield4th = 0., dYield6th = 0., dYield8th = 0.; // pt yield 
 Double_t dSum2nd = 0., dSum4th = 0., dSum6th = 0., dSum8th = 0.; // needed for normalizing integrated flow
 
 // looping over pt bins:
 for(Int_t p=1;p<nBinsPt+1;p++)
 {
  dvn2nd = flow2ndPt->GetBinContent(p);
  dvn4th = flow4thPt->GetBinContent(p);
  dvn6th = flow6thPt->GetBinContent(p);
  dvn8th = flow8thPt->GetBinContent(p);
  
  dErrvn2nd = flow2ndPt->GetBinError(p);
  dErrvn4th = flow4thPt->GetBinError(p);
  dErrvn6th = flow6thPt->GetBinError(p);
  dErrvn8th = flow8thPt->GetBinError(p);

  dYield2nd = yield2ndPt->GetBinContent(p);  
  dYield4th = yield4thPt->GetBinContent(p);
  dYield6th = yield6thPt->GetBinContent(p);
  dYield8th = yield8thPt->GetBinContent(p);
  
  dVn2nd += dvn2nd*dYield2nd;
  dVn4th += dvn4th*dYield4th;
  dVn6th += dvn6th*dYield6th;
  dVn8th += dvn8th*dYield8th;
  
  dSum2nd += dYield2nd;
  dSum4th += dYield4th;
  dSum6th += dYield6th;
  dSum8th += dYield8th;
  
  dErrVn2nd += dYield2nd*dYield2nd*dErrvn2nd*dErrvn2nd; // ro be improved (check this relation)
  dErrVn4th += dYield4th*dYield4th*dErrvn4th*dErrvn4th;
  dErrVn6th += dYield6th*dYield6th*dErrvn6th*dErrvn6th;
  dErrVn8th += dYield8th*dYield8th*dErrvn8th*dErrvn8th;
    
 } // end of for(Int_t p=1;p<nBinsPt+1;p++)

 // normalizing the results for integrated flow:
 if(dSum2nd) 
 {
  dVn2nd /= dSum2nd;
  dErrVn2nd /= (dSum2nd*dSum2nd);
  dErrVn2nd = TMath::Sqrt(dErrVn2nd);
 } 
 if(dSum4th) 
 {
  dVn4th /= dSum4th;
  dErrVn4th /= (dSum4th*dSum4th);
  dErrVn4th = TMath::Sqrt(dErrVn4th);
 } 
 //if(dSum6th) dVn6th/=dSum6th;
 //if(dSum8th) dVn8th/=dSum8th;
  
 // storing the results for integrated flow in common histos: (to be improved: new method for this?)
 if(type == "POI")
 {
  fCommonHistsResults2nd->FillIntegratedFlowPOI(dVn2nd,dErrVn2nd); 
  fCommonHistsResults4th->FillIntegratedFlowPOI(dVn4th,dErrVn4th); 
  fCommonHistsResults6th->FillIntegratedFlowPOI(dVn6th,0.); // to be improved (errors)
  fCommonHistsResults8th->FillIntegratedFlowPOI(dVn8th,0.); // to be improved (errors)
 }
 else if (type == "RP")
 {
  fCommonHistsResults2nd->FillIntegratedFlowRP(dVn2nd,dErrVn2nd); 
  fCommonHistsResults4th->FillIntegratedFlowRP(dVn4th,dErrVn4th);
  fCommonHistsResults6th->FillIntegratedFlowRP(dVn6th,0.); // to be improved (errors)
  fCommonHistsResults8th->FillIntegratedFlowRP(dVn8th,0.); // to be improved (errors)
 }
 
 delete flow2ndPt;
 delete flow4thPt;
 //delete flow6thPt;
 //delete flow8thPt;
 
 delete yield2ndPt;
 delete yield4thPt;
 delete yield6thPt;
 delete yield8thPt;
           
} // end of AliFlowAnalysisWithQCumulants::CalculateFinalResultsForRPandPOIIntegratedFlow(TString type)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::InitializeArraysForDistributions()
{
 // Initialize all arrays used for distributions.
 
 // a) Initialize arrays of histograms used to hold distributions of correlations; 
 // b) Initialize array to hold min and max values of correlations.
 
 // a) Initialize arrays of histograms used to hold distributions of correlations:
 for(Int_t di=0;di<4;di++) // distribution index
 {
  fDistributions[di] = NULL;
 }
 
 // b) Initialize default min and max values of correlations:
 //    (Remark: The default values bellow were chosen for v2=5% and M=500)
 fMinValueOfCorrelation[0] = -0.01; // <2>_min 
 fMaxValueOfCorrelation[0] = 0.04; // <2>_max 
 fMinValueOfCorrelation[1] = -0.00002; // <4>_min 
 fMaxValueOfCorrelation[1] = 0.00015; // <4>_max  
 fMinValueOfCorrelation[2] = -0.0000003; // <6>_min 
 fMaxValueOfCorrelation[2] = 0.0000006; // <6>_max  
 fMinValueOfCorrelation[3] = -0.000000006; // <8>_min 
 fMaxValueOfCorrelation[3] = 0.000000003; // <8>_max 
 
} // end of void AliFlowAnalysisWithQCumulants::InitializeArraysForDistributions()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::InitializeArraysForVarious()
{
 // Initialize all arrays used for various unclassified objects.
 
 for(Int_t p=0;p<4;p++) // [v_min,v_max,refMult_min,refMult_max]
 {
  fPhiDistributionForOneEventSettings[p] = 0.;
 } 
   
} //  end of void AliFlowAnalysisWithQCumulants::InitializeArraysForVarious()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::BookEverythingForDistributions()
{
 // a) Book profile to hold all flags for distributions of correlations;
 // b) Book all histograms to hold distributions of correlations.
 
 TString correlationIndex[4] = {"<2>","<4>","<6>","<8>"}; // to be improved (should I promote this to data members?)
  
 // a) Book profile to hold all flags for distributions of correlations:
 TString distributionsFlagsName = "fDistributionsFlags";
 distributionsFlagsName += fAnalysisLabel->Data();
 fDistributionsFlags = new TProfile(distributionsFlagsName.Data(),"Flags for Distributions of Correlations",9,0,9);
 fDistributionsFlags->SetTickLength(-0.01,"Y");
 fDistributionsFlags->SetMarkerStyle(25);
 fDistributionsFlags->SetLabelSize(0.05);
 fDistributionsFlags->SetLabelOffset(0.02,"Y");
 fDistributionsFlags->GetXaxis()->SetBinLabel(1,"Store or not?");
 fDistributionsFlags->GetXaxis()->SetBinLabel(2,"<2>_{min}");
 fDistributionsFlags->GetXaxis()->SetBinLabel(3,"<2>_{max}");
 fDistributionsFlags->GetXaxis()->SetBinLabel(4,"<4>_{min}");
 fDistributionsFlags->GetXaxis()->SetBinLabel(5,"<4>_{max}");
 fDistributionsFlags->GetXaxis()->SetBinLabel(6,"<6>_{min}");
 fDistributionsFlags->GetXaxis()->SetBinLabel(7,"<6>_{max}");
 fDistributionsFlags->GetXaxis()->SetBinLabel(8,"<8>_{min}");
 fDistributionsFlags->GetXaxis()->SetBinLabel(9,"<8>_{max}");
 fDistributionsList->Add(fDistributionsFlags);
 
 // b) Book all histograms to hold distributions of correlations.
 if(fStoreDistributions)
 { 
  TString distributionsName = "fDistributions";
  distributionsName += fAnalysisLabel->Data();
  for(Int_t di=0;di<4;di++) // distribution index
  {
   fDistributions[di] = new TH1D(Form("Distribution of %s",correlationIndex[di].Data()),Form("Distribution of %s",correlationIndex[di].Data()),10000,fMinValueOfCorrelation[di],fMaxValueOfCorrelation[di]); 
   fDistributions[di]->SetXTitle(correlationIndex[di].Data());
   fDistributionsList->Add(fDistributions[di]);
  } // end of for(Int_t di=0;di<4;di++) // distribution index
 } // end of if(fStoreDistributions)
 
} // end of void AliFlowAnalysisWithQCumulants::BookEverythingForDistributions()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::BookEverythingForVarious()
{
 // Book all objects for various unclassified quantities.
 
 if(!fStorePhiDistributionForOneEvent){return;}
 
 // a) Book histogram holding phi distribution for single event to illustrate flow.
 
 // a) Book histogram holding phi distribution for single event to illustrate flow:
 fPhiDistributionForOneEvent = new TH1D("fPhiDistributionForOneEvent","",360,0.,TMath::TwoPi());
 fPhiDistributionForOneEvent->GetXaxis()->SetTitle("#phi");
 fVariousList->Add(fPhiDistributionForOneEvent);
 
} // end of void AliFlowAnalysisWithQCumulants::BookEverythingForVarious()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::StoreFlagsForDistributions()
{
 // Store all flags for distributiuons of correlations in profile fDistributionsFlags.
 
 if(!fDistributionsFlags)
 {
  cout<<"WARNING: fDistributionsFlags is NULL in AFAWQC::SDF() !!!!"<<endl;
  exit(0);
 } 

 fDistributionsFlags->Fill(0.5,(Int_t)fStoreDistributions); // histos with distributions of correlations stored or not in the output file
 // store min and max values of correlations:
 for(Int_t di=0;di<4;di++) // distribution index
 {
  fDistributionsFlags->Fill(1.5+2.*(Double_t)di,fMinValueOfCorrelation[di]);
  fDistributionsFlags->Fill(2.5+2.*(Double_t)di,fMaxValueOfCorrelation[di]);
 }
     
} // end of void AliFlowAnalysisWithQCumulants::StoreFlagsForDistributions()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::StoreDistributionsOfCorrelations()
{
 // Store distributions of correlations.
 
 if(!(fIntFlowCorrelationsEBE && fIntFlowEventWeightsForCorrelationsEBE))
 {
  cout<<"WARNING: fIntFlowCorrelationsEBE && fIntFlowEventWeightsForCorrelationsEBE"<<endl; 
  cout<<"         is NULL in AFAWQC::SDOC() !!!!"<<endl;
  exit(0);
 }

 for(Int_t di=0;di<4;di++) // distribution index
 {
  if(!fDistributions[di])
  { 
   cout<<"WARNING: fDistributions[di] is NULL in AFAWQC::SDOC() !!!!"<<endl;
   cout<<"di = "<<di<<endl;
   exit(0);
  } else 
    {
     fDistributions[di]->Fill(fIntFlowCorrelationsEBE->GetBinContent(di+1),fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(di+1)); 
    } 
 } // end of for(Int_t di=0;di<4;di++) // distribution index

} // end of void AliFlowAnalysisWithQCumulants::StoreDistributionsOfCorrelations()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.
 //  a) Book and nest lists for integrated flow;
 //  b) Book and nest lists for differential flow;
 //  c) Book and nest list for particle weights;
 //  d) Book and nest list for distributions;
 //  e) Book and nest list for various unclassified objects; 
 //  f) Book and nest list for nested loops.
 
 // a) Book and nest all lists for integrated flow:
 //  Base list for integrated flow:
 fIntFlowList = new TList();
 fIntFlowList->SetName("Integrated Flow");
 fIntFlowList->SetOwner(kTRUE);
 fHistList->Add(fIntFlowList);
 //  List holding profiles: 
 fIntFlowProfiles = new TList();
 fIntFlowProfiles->SetName("Profiles");
 fIntFlowProfiles->SetOwner(kTRUE);
 fIntFlowList->Add(fIntFlowProfiles);
 //  List holding all profiles with results for correlations vs M:
 if(fCalculateAllCorrelationsVsM)
 {
  fIntFlowAllCorrelationsVsM = new TList();
  fIntFlowAllCorrelationsVsM->SetName("Correlations vs M");
  fIntFlowAllCorrelationsVsM->SetOwner(kTRUE); 
  fIntFlowProfiles->Add(fIntFlowAllCorrelationsVsM);
 } // end of if(fCalculateAllCorrelationsVsM)
 //  List holding histograms with results:
 fIntFlowResults = new TList();
 fIntFlowResults->SetName("Results");
 fIntFlowResults->SetOwner(kTRUE);
 fIntFlowList->Add(fIntFlowResults);
 
 // b) Book and nest lists for differential flow:
 this->BookAndNestListsForDifferentialFlow();
 
 // c) Book and nest list for particle weights:
 fWeightsList->SetName("Weights");
 fWeightsList->SetOwner(kTRUE);   
 fHistList->Add(fWeightsList); 

 // d) Book and nest list for distributions:
 fDistributionsList = new TList();
 fDistributionsList->SetName("Distributions");
 fDistributionsList->SetOwner(kTRUE);
 fHistList->Add(fDistributionsList);
 
 // e) Book and nest list for various unclassified objects:
 if(fStorePhiDistributionForOneEvent)
 {
  fVariousList = new TList();
  fVariousList->SetName("Various");
  fVariousList->SetOwner(kTRUE);
  fHistList->Add(fVariousList);
 }
  
 // f) Book and nest list for other differential correlators:
 fOtherDiffCorrelatorsList = new TList();
 fOtherDiffCorrelatorsList->SetName("Other differential correlators");
 fOtherDiffCorrelatorsList->SetOwner(kTRUE);
 if(fCalculateDiffFlow){fHistList->Add(fOtherDiffCorrelatorsList);} // TBI: Use another flag here instead of fCalculateDiffFlow 
  
 // g) Book and nest list for nested loops:
 fNestedLoopsList = new TList();
 fNestedLoopsList->SetName("Nested Loops");
 fNestedLoopsList->SetOwner(kTRUE);
 fHistList->Add(fNestedLoopsList);
 
} // end of void AliFlowAnalysisWithQCumulants::BookAndNestAllLists()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::BookAndNestListsForDifferentialFlow()
{
 // Book and nest lists for differential flow.

 // Base list for differential flow objects:
 fDiffFlowList = new TList();
 fDiffFlowList->SetName("Differential Flow");
 fDiffFlowList->SetOwner(kTRUE); 
 fHistList->Add(fDiffFlowList);
 
 // Local flags: 
 TString typeFlag[2] = {"RP","POI"};  
 TString ptEtaFlag[2] = {"p_{T}","#eta"}; 
 TString powerFlag[2] = {"linear","quadratic"};   

 // 2D:
 if(fCalculate2DDiffFlow)
 {
  fDiffFlow2D = new TList(); 
  fDiffFlow2D->SetName("2D");
  fDiffFlow2D->SetOwner(kTRUE);
  fDiffFlowList->Add(fDiffFlow2D); 
  for(Int_t t=0;t<2;t++)
  {
   f2DDiffFlowCorrelationsProList[t] = new TList();
   f2DDiffFlowCorrelationsProList[t]->SetOwner(kTRUE);
   f2DDiffFlowCorrelationsProList[t]->SetName(Form("Profiles with 2D correlations (%s)",typeFlag[t].Data()));
   fDiffFlow2D->Add(f2DDiffFlowCorrelationsProList[t]);
  } // end of for(Int_t t=0;t<2;t++)
 } // end of if(fCalculate2DDiffFlow)

 // What follows bellow in this method is relevant only for 1D differential flow:
 if(!fCalculateDiffFlow){return;}
 
 // List holding profiles: 
 fDiffFlowProfiles = new TList(); 
 fDiffFlowProfiles->SetName("Profiles");
 fDiffFlowProfiles->SetOwner(kTRUE);
 fDiffFlowList->Add(fDiffFlowProfiles);
 // List holding histograms with results: 
 fDiffFlowResults = new TList();
 fDiffFlowResults->SetName("Results");
 fDiffFlowResults->SetOwner(kTRUE);
 fDiffFlowList->Add(fDiffFlowResults);
 // Flags used for naming nested lists in list fDiffFlowProfiles and fDiffFlowResults:  
 TList list;
 list.SetOwner(kTRUE);
 // Nested lists in fDiffFlowProfiles (~/Differential Flow/Profiles):
 for(Int_t t=0;t<2;t++) // type: RP or POI
 {
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   // list holding profiles with correlations:
   fDiffFlowCorrelationsProList[t][pe] = (TList*)list.Clone();
   fDiffFlowCorrelationsProList[t][pe]->SetName(Form("Profiles with correlations (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowProfiles->Add(fDiffFlowCorrelationsProList[t][pe]);
   // list holding profiles with products of correlations:
   fDiffFlowProductOfCorrelationsProList[t][pe] = (TList*)list.Clone();
   fDiffFlowProductOfCorrelationsProList[t][pe]->SetName(Form("Profiles with products of correlations (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowProfiles->Add(fDiffFlowProductOfCorrelationsProList[t][pe]);
   // list holding profiles with corrections:
   fDiffFlowCorrectionsProList[t][pe] = (TList*)list.Clone();
   fDiffFlowCorrectionsProList[t][pe]->SetName(Form("Profiles with correction terms for NUA (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowProfiles->Add(fDiffFlowCorrectionsProList[t][pe]);   
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta 
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI   
 // nested lists in fDiffFlowResults (~/Differential Flow/Results):
 for(Int_t t=0;t<2;t++) // type: RP or POI
 {
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   // list holding histograms with correlations:
   fDiffFlowCorrelationsHistList[t][pe] = (TList*)list.Clone();
   fDiffFlowCorrelationsHistList[t][pe]->SetName(Form("Correlations (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowResults->Add(fDiffFlowCorrelationsHistList[t][pe]);
   // list holding histograms with corrections:
   fDiffFlowCorrectionsHistList[t][pe] = (TList*)list.Clone();
   fDiffFlowCorrectionsHistList[t][pe]->SetName(Form("Histograms with correction terms for NUA (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowResults->Add(fDiffFlowCorrectionsHistList[t][pe]);   
   for(Int_t power=0;power<2;power++)
   {
    // list holding histograms with sums of event weights:
    fDiffFlowSumOfEventWeightsHistList[t][pe][power] = (TList*)list.Clone();
    fDiffFlowSumOfEventWeightsHistList[t][pe][power]->SetName(Form("Sum of %s event weights (%s, %s)",powerFlag[power].Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data()));
    fDiffFlowResults->Add(fDiffFlowSumOfEventWeightsHistList[t][pe][power]);    
   } // end of for(Int_t power=0;power<2;power++)
   // list holding histograms with sums of products of event weights:
   fDiffFlowSumOfProductOfEventWeightsHistList[t][pe] = (TList*)list.Clone();
   fDiffFlowSumOfProductOfEventWeightsHistList[t][pe]->SetName(Form("Sum of products of event weights (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowResults->Add(fDiffFlowSumOfProductOfEventWeightsHistList[t][pe]);
   // list holding histograms with covariances of correlations:
   fDiffFlowCovariancesHistList[t][pe] = (TList*)list.Clone();
   fDiffFlowCovariancesHistList[t][pe]->SetName(Form("Covariances of correlations (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowResults->Add(fDiffFlowCovariancesHistList[t][pe]);
   // list holding histograms with differential Q-cumulants:
   fDiffFlowCumulantsHistList[t][pe] = (TList*)list.Clone();
   fDiffFlowCumulantsHistList[t][pe]->SetName(Form("Differential Q-cumulants (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowResults->Add(fDiffFlowCumulantsHistList[t][pe]);   
   // list holding histograms which quantify detector bias to differential Q-cumulants:
   fDiffFlowDetectorBiasHistList[t][pe] = (TList*)list.Clone();
   fDiffFlowDetectorBiasHistList[t][pe]->SetName(Form("Detector bias (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowResults->Add(fDiffFlowDetectorBiasHistList[t][pe]);   
   // list holding histograms with differential flow estimates from Q-cumulants:
   fDiffFlowHistList[t][pe] = (TList*)list.Clone();
   fDiffFlowHistList[t][pe]->SetName(Form("Differential flow (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowResults->Add(fDiffFlowHistList[t][pe]);      
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
  
} // end of void AliFlowAnalysisWithQCumulants::BookAndNestListsForDifferentialFlow()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::FillCommonHistResultsDiffFlow(TString type)
{
 // Fill common result histograms for differential flow.
 
 Int_t t = 0; 

 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   } 
  
 // to be improved - check all pointers used in this method
     
 if(!(fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th))
 {
  cout<<"WARNING: fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th"<<endl; 
  cout<<"         is NULL in AFAWQC::FCHRIF() !!!!"<<endl;
  exit(0);
 }
 
 // pt:
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  Double_t v2 = fDiffFlow[t][0][0]->GetBinContent(p);
  Double_t v4 = fDiffFlow[t][0][1]->GetBinContent(p);
  Double_t v6 = fDiffFlow[t][0][2]->GetBinContent(p);
  Double_t v8 = fDiffFlow[t][0][3]->GetBinContent(p);
  
  Double_t v2Error = fDiffFlow[t][0][0]->GetBinError(p);
  Double_t v4Error = fDiffFlow[t][0][1]->GetBinError(p);
  //Double_t v6Error = fFinalFlow1D[t][pW][nua][0][2]->GetBinError(p);
  //Double_t v8Error = fFinalFlow1D[t][pW][nua][0][3]->GetBinError(p);
 
  if(type == "RP")
  {
   fCommonHistsResults2nd->FillDifferentialFlowPtRP(p,v2,v2Error);
   fCommonHistsResults4th->FillDifferentialFlowPtRP(p,v4,v4Error);
   fCommonHistsResults6th->FillDifferentialFlowPtRP(p,v6,0.);
   fCommonHistsResults8th->FillDifferentialFlowPtRP(p,v8,0.);
  } else if(type == "POI")
    {
     fCommonHistsResults2nd->FillDifferentialFlowPtPOI(p,v2,v2Error);
     fCommonHistsResults4th->FillDifferentialFlowPtPOI(p,v4,v4Error);
     fCommonHistsResults6th->FillDifferentialFlowPtPOI(p,v6,0.);
     fCommonHistsResults8th->FillDifferentialFlowPtPOI(p,v8,0.);
    }
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)   
 
 // eta:
 if(!fCalculateDiffFlowVsEta){return;}
 for(Int_t e=1;e<=fnBinsEta;e++)
 {
  Double_t v2 = fDiffFlow[t][1][0]->GetBinContent(e);
  Double_t v4 = fDiffFlow[t][1][1]->GetBinContent(e);
  Double_t v6 = fDiffFlow[t][1][2]->GetBinContent(e);
  Double_t v8 = fDiffFlow[t][1][3]->GetBinContent(e);
  
  Double_t v2Error = fDiffFlow[t][1][0]->GetBinError(e);
  Double_t v4Error = fDiffFlow[t][1][1]->GetBinError(e);
  //Double_t v6Error = fDiffFlow[t][1][2]->GetBinError(e);
  //Double_t v8Error = fDiffFlow[t][1][3]->GetBinError(e);
 
  if(type == "RP")
  {
   fCommonHistsResults2nd->FillDifferentialFlowEtaRP(e,v2,v2Error);
   fCommonHistsResults4th->FillDifferentialFlowEtaRP(e,v4,v4Error);
   fCommonHistsResults6th->FillDifferentialFlowEtaRP(e,v6,0.);
   fCommonHistsResults8th->FillDifferentialFlowEtaRP(e,v8,0.);
  } else if(type == "POI")
    {
     fCommonHistsResults2nd->FillDifferentialFlowEtaPOI(e,v2,v2Error);
     fCommonHistsResults4th->FillDifferentialFlowEtaPOI(e,v4,v4Error);
     fCommonHistsResults6th->FillDifferentialFlowEtaPOI(e,v6,0.);
     fCommonHistsResults8th->FillDifferentialFlowEtaPOI(e,v8,0.);
    }
 } // end of for(Int_t e=1;e<=fnBinsEta;e++)    
 
} // end of void AliFlowAnalysisWithQCumulants::FillCommonHistResultsDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights, Bool_t correctedForNUA)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CommonConstants(TString method)
{
 // Access and store common constants.
 
 // a) If this method was called in Init() access common constants from AliFlowCommonConstants;
 // b) If this method was called in Init() book and fill TProfile to hold constants accessed in a);
 // c) If this method was called in Finish() access common constants from TProfile booked and filled in b).

 if(method == "Init")
 {
  // a) If this method was called in Init() access common constants from AliFlowCommonConstants:
  fnBinsPhi = AliFlowCommonConstants::GetMaster()->GetNbinsPhi();
  fPhiMin = AliFlowCommonConstants::GetMaster()->GetPhiMin();	     
  fPhiMax = AliFlowCommonConstants::GetMaster()->GetPhiMax();
  if(fnBinsPhi){fPhiBinWidth = (fPhiMax-fPhiMin)/fnBinsPhi;}  
  fnBinsPt = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
  fPtMin = AliFlowCommonConstants::GetMaster()->GetPtMin();	     
  fPtMax = AliFlowCommonConstants::GetMaster()->GetPtMax();
  if(fnBinsPt){fPtBinWidth = (fPtMax-fPtMin)/fnBinsPt;}  
  fnBinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
  fEtaMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();	     
  fEtaMax = AliFlowCommonConstants::GetMaster()->GetEtaMax();
  if(fnBinsEta){fEtaBinWidth = (fEtaMax-fEtaMin)/fnBinsEta;}  
  
  // b) If this method was called in Init() book and fill TProfile to hold constants accessed in a):
  TString fCommonConstantsName = "fCommonConstants";
  fCommonConstantsName += fAnalysisLabel->Data();
  fCommonConstants = new TProfile(fCommonConstantsName.Data(),"Common constants",9,0.,9.);
  fCommonConstants->SetLabelSize(0.05);
  fCommonConstants->GetXaxis()->SetBinLabel(1,"nBins (#phi)");
  fCommonConstants->Fill(0.5,fnBinsPhi);
  fCommonConstants->GetXaxis()->SetBinLabel(2,"#phi_{min}");
  fCommonConstants->Fill(1.5,fPhiMin);
  fCommonConstants->GetXaxis()->SetBinLabel(3,"#phi_{max}");
  fCommonConstants->Fill(2.5,fPhiMax);
  fCommonConstants->GetXaxis()->SetBinLabel(4,"nBins (p_{t})");
  fCommonConstants->Fill(3.5,fnBinsPt);
  fCommonConstants->GetXaxis()->SetBinLabel(5,"(p_{t})_{min}");
  fCommonConstants->Fill(4.5,fPtMin);
  fCommonConstants->GetXaxis()->SetBinLabel(6,"(p_{t})_{max}");
  fCommonConstants->Fill(5.5,fPtMax);
  fCommonConstants->GetXaxis()->SetBinLabel(7,"nBins (#eta)");
  fCommonConstants->Fill(6.5,fnBinsEta);
  fCommonConstants->GetXaxis()->SetBinLabel(8,"#eta_{min}");
  fCommonConstants->Fill(7.5,fEtaMin);
  fCommonConstants->GetXaxis()->SetBinLabel(9,"#eta_{max}");
  fCommonConstants->Fill(8.5,fEtaMax);
  fHistList->Add(fCommonConstants); 
 } // end of if(method == "Init")
 else if(method == "Finish")
 {
  // c) If this method was called in Finish() access common constants from TProfile booked and filled in b):
  if(!fCommonConstants)
  {
   printf("\n WARNING (QC): fCommonConstants is NULL in AFAWQC::AC(\"%s\") !!!!\n\n",method.Data());
   exit(0);
  } 
  fnBinsPhi = (Int_t)fCommonConstants->GetBinContent(1);
  fPhiMin = fCommonConstants->GetBinContent(2);	     
  fPhiMax = fCommonConstants->GetBinContent(3);
  if(fnBinsPhi){fPhiBinWidth = (fPhiMax-fPhiMin)/fnBinsPhi;}  
  fnBinsPt = (Int_t)fCommonConstants->GetBinContent(4);
  fPtMin = fCommonConstants->GetBinContent(5);	     
  fPtMax = fCommonConstants->GetBinContent(6);
  if(fnBinsPt){fPtBinWidth = (fPtMax-fPtMin)/fnBinsPt;}  
  fnBinsEta = (Int_t)fCommonConstants->GetBinContent(7);
  fEtaMin = fCommonConstants->GetBinContent(8);	     
  fEtaMax = fCommonConstants->GetBinContent(9);
  if(fnBinsEta){fEtaBinWidth = (fEtaMax-fEtaMin)/fnBinsEta;}  
 } // end of else if(method == "Finish")

} // end of void AliFlowAnalysisWithQCumulants::CommonConstants(TString method)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CrossCheckSettings()
{
 // a) Cross check if the choice for multiplicity weights make sense.
 
 // a) Cross check if the choice for multiplicity weights make sense:
 if(strcmp(fMultiplicityWeight->Data(),"combinations") && 
    strcmp(fMultiplicityWeight->Data(),"unit") &&
    strcmp(fMultiplicityWeight->Data(),"multiplicity"))
 {
  cout<<"WARNING (QC): Multiplicity weight can be either \"combinations\", \"unit\""<<endl;
  cout<<"              or \"multiplicity\". Certainly not \""<<fMultiplicityWeight->Data()<<"\"."<<endl;
  exit(0);
 }   
 
} // end of void AliFlowAnalysisWithQCumulants::CrossCheckSettings()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateIntFlowSumOfEventWeights()
{
 // Calculate sum of linear and quadratic event weights for correlations.
 
 // multiplicity:
 Double_t dMult = (*fSpk)(0,0);
                        
 for(Int_t p=0;p<2;p++) // power-1
 {
  for(Int_t ci=0;ci<4;ci++) // correlation index
  { 
   fIntFlowSumOfEventWeights[p]->Fill(ci+0.5,pow(fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci+1),p+1)); 
   if(fCalculateCumulantsVsM)
   {
    fIntFlowSumOfEventWeightsVsM[ci][p]->Fill(dMult+0.5,pow(fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci+1),p+1)); // to be improved: dMult => sum of weights?
   }
  }
 }
  
} // end of void AliFlowAnalysisWithQCumulants::CalculateIntFlowSumOfEventWeights()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateIntFlowSumOfEventWeightsNUA()
{
 // Calculate sum of linear and quadratic event weights for NUA terms.
                       
 for(Int_t sc=0;sc<2;sc++) // sin or cos terms
 {
  for(Int_t p=0;p<2;p++) // power-1
  {
   for(Int_t ci=0;ci<4;ci++) // nua term index
   { 
    fIntFlowSumOfEventWeightsNUA[sc][p]->Fill(ci+0.5,pow(fIntFlowEventWeightForCorrectionTermsForNUAEBE[sc]->GetBinContent(ci+1),p+1)); 
   }
  }
 }
  
} // end of void AliFlowAnalysisWithQCumulants::CalculateIntFlowSumOfEventWeightsNUA()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateIntFlowSumOfProductOfEventWeights()
{
 // Calculate sum of product of event weights for correlations.
  
 // multiplicity:
 Double_t dMult = (*fSpk)(0,0);
  
 Int_t counter = 0;
 
 for(Int_t ci1=1;ci1<4;ci1++)
 {
  for(Int_t ci2=ci1+1;ci2<=4;ci2++)
  {
   fIntFlowSumOfProductOfEventWeights->Fill(0.5+counter,
                                            fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci1)*
                                            fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci2));
   if(fCalculateCumulantsVsM)
   {                                                                                    
    fIntFlowSumOfProductOfEventWeightsVsM[counter]->Fill(dMult+0.5, // to be improved: dMult => sum of weights?
                                                         fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci1)*
                                                         fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci2));
   } // end of if(fCalculateCumulantsVsM)
   counter++;                                         
  }
 }

} // end of void AliFlowAnalysisWithQCumulants::CalculateIntFlowSumOfProductOfEventWeights()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateIntFlowSumOfProductOfEventWeightsNUA()
{
 // Calculate sum of product of event weights for NUA terms.
  
 // w_{<2>} * w_{<cos(#phi)>}:
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(0.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)*
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1));
 // w_{<2>} * w_{<sin(#phi)>}:
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(1.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)*
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1));
 // w_{<cos(#phi)> * w_{<sin(#phi)>}:
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(2.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1)*
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1));
 // w_{<2>} * w{<cos(phi1+phi2)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(3.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)*
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)); 
 // w_{<2>} * w{<sin(phi1+phi2)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(4.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)*
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2));
 // w_{<2>} * w{<cos(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(5.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)*
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3));
 // w_{<2>} * w{<sin(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(6.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1)*
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3));  
 // w_{<4>} * w{<cos(phi1)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(7.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)*
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1));
 // w_{<4>} * w{<sin(phi1)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(8.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)*
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1));
 // w_{<4>} * w{<cos(phi1+phi2)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(9.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)*
                                                 fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)); 
 // w_{<4>} * w{<sin(phi1+phi2)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(10.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2));
 // w_{<4>} * w{<cos(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(11.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3));
 // w_{<4>} * w{<sin(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(12.5,fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3));
 // w_{<cos(phi1)>} * w{<cos(phi1+phi2)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(13.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2));
 // w_{<cos(phi1)>} * w{<sin(phi1+phi2)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(14.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2)); 
 // w_{<cos(phi1)>} * w{<cos(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(15.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3));
 // w_{<cos(phi1)>} * w{<sin(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(16.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(1)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3));
 // w_{<sin(phi1)>} * w{<cos(phi1+phi2)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(17.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2));
 // w_{<sin(phi1)>} * w{<sin(phi1+phi2)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(18.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2));
 // w_{<sin(phi1)>} * w{<cos(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(19.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3));
 // w_{<sin(phi1)>} * w{<sin(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(20.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(1)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3)); 
 // w_{<cos(phi1+phi2)>} * w{<sin(phi1+phi2))>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(21.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2)); 
 // w_{<cos(phi1+phi2)>} * w{<cos(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(22.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3)); 
 // w_{<cos(phi1+phi2)>} * w{<sin(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(23.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(2)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3)); 
 // w_{<sin(phi1+phi2)>} * w{<cos(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(24.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3)); 
 // w_{<sin(phi1+phi2)>} * w{<sin(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(25.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(2)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3)); 
 // w_{<cos(phi1-phi2-phi3)>} * w{<sin(phi1-phi2-phi3)>}
 fIntFlowSumOfProductOfEventWeightsNUA->Fill(26.5,fIntFlowEventWeightForCorrectionTermsForNUAEBE[1]->GetBinContent(3)*
                                                  fIntFlowEventWeightForCorrectionTermsForNUAEBE[0]->GetBinContent(3));

} // end of void AliFlowAnalysisWithQCumulants::CalculateIntFlowIntFlowSumOfProductOfEventWeightsNUA()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrelations(TString type, TString ptOrEta)
{
 // Calculate reduced correlations for RPs or POIs for all pt and eta bins.

 // Multiplicity:
 Double_t dMult = (*fSpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 //Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 //Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);

 // reduced correlations are stored in fDiffFlowCorrelationsPro[0=RP,1=POI][0=pt,1=eta][correlation index]. Correlation index runs as follows:
 // 
 // 0: <<2'>>
 // 1: <<4'>>
 // 2: <<6'>>
 // 3: <<8'>>
 
 Int_t t = 0; // type flag 
 Int_t pe = 0; // ptEta flag
 
 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }

 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   }
    
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};

 // looping over all bins and calculating reduced correlations: 
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // real and imaginary parts of p_{m*n,0} (non-weighted Q-vector evaluated for POIs in particular pt or eta bin): 
  Double_t p1n0kRe = 0.;
  Double_t p1n0kIm = 0.;

  // number of POIs in particular pt or eta bin:
  Double_t mp = 0.;

  // real and imaginary parts of q_{m*n,0} (non-weighted Q-vector evaluated for particles which are both RPs and POIs in particular pt or eta bin):
  Double_t q1n0kRe = 0.;
  Double_t q1n0kIm = 0.;
  Double_t q2n0kRe = 0.;
  Double_t q2n0kIm = 0.;

  // number of particles which are both RPs and POIs in particular pt or eta bin:
  Double_t mq = 0.;
   
  if(type == "POI")
  {
   // q_{m*n,0}:
   q1n0kRe = fReRPQ1dEBE[2][pe][0][0]->GetBinContent(fReRPQ1dEBE[2][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(fReRPQ1dEBE[2][pe][0][0]->GetBin(b));
   q1n0kIm = fImRPQ1dEBE[2][pe][0][0]->GetBinContent(fImRPQ1dEBE[2][pe][0][0]->GetBin(b))
           * fImRPQ1dEBE[2][pe][0][0]->GetBinEntries(fImRPQ1dEBE[2][pe][0][0]->GetBin(b));
   q2n0kRe = fReRPQ1dEBE[2][pe][1][0]->GetBinContent(fReRPQ1dEBE[2][pe][1][0]->GetBin(b))
           * fReRPQ1dEBE[2][pe][1][0]->GetBinEntries(fReRPQ1dEBE[2][pe][1][0]->GetBin(b));
   q2n0kIm = fImRPQ1dEBE[2][pe][1][0]->GetBinContent(fImRPQ1dEBE[2][pe][1][0]->GetBin(b))
           * fImRPQ1dEBE[2][pe][1][0]->GetBinEntries(fImRPQ1dEBE[2][pe][1][0]->GetBin(b));         
                 
   mq = fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(fReRPQ1dEBE[2][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
  } 
  else if(type == "RP")
  {
   // q_{m*n,0}:
   q1n0kRe = fReRPQ1dEBE[0][pe][0][0]->GetBinContent(fReRPQ1dEBE[0][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b));
   q1n0kIm = fImRPQ1dEBE[0][pe][0][0]->GetBinContent(fImRPQ1dEBE[0][pe][0][0]->GetBin(b))
           * fImRPQ1dEBE[0][pe][0][0]->GetBinEntries(fImRPQ1dEBE[0][pe][0][0]->GetBin(b));
   q2n0kRe = fReRPQ1dEBE[0][pe][1][0]->GetBinContent(fReRPQ1dEBE[0][pe][1][0]->GetBin(b))
           * fReRPQ1dEBE[0][pe][1][0]->GetBinEntries(fReRPQ1dEBE[0][pe][1][0]->GetBin(b));
   q2n0kIm = fImRPQ1dEBE[0][pe][1][0]->GetBinContent(fImRPQ1dEBE[0][pe][1][0]->GetBin(b))
           * fImRPQ1dEBE[0][pe][1][0]->GetBinEntries(fImRPQ1dEBE[0][pe][1][0]->GetBin(b));         
                 
   mq = fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)  
  }
      
   if(type == "POI")
   {
    // p_{m*n,0}:
    p1n0kRe = fReRPQ1dEBE[1][pe][0][0]->GetBinContent(fReRPQ1dEBE[1][pe][0][0]->GetBin(b))
            * fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b));
    p1n0kIm = fImRPQ1dEBE[1][pe][0][0]->GetBinContent(fImRPQ1dEBE[1][pe][0][0]->GetBin(b))  
            * fImRPQ1dEBE[1][pe][0][0]->GetBinEntries(fImRPQ1dEBE[1][pe][0][0]->GetBin(b));
            
    mp = fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
    
    t = 1; // typeFlag = RP or POI
   }
   else if(type == "RP")
   {
    // p_{m*n,0} = q_{m*n,0}:
    p1n0kRe = q1n0kRe; 
    p1n0kIm = q1n0kIm; 
            
    mp = mq; 
    
    t = 0; // typeFlag = RP or POI
   }
      
   // 2'-particle correlation for particular pt or eta bin:
   Double_t two1n1nPtEta = 0.;
   Double_t mWeight2pPrime = 0.; // multiplicity weight for <2'>
   if(mp*dMult-mq)
   {
    two1n1nPtEta = (p1n0kRe*dReQ1n+p1n0kIm*dImQ1n-mq)
                 / (mp*dMult-mq);
    // determine multiplicity weight:
    if(!strcmp(fMultiplicityWeight->Data(),"combinations"))
    {
     mWeight2pPrime = mp*dMult-mq;
    } else if(!strcmp(fMultiplicityWeight->Data(),"unit"))
      {
       mWeight2pPrime = 1.;    
      } 
    if(type == "POI") // to be improved (I do not this if)
    { 
     // fill profile to get <<2'>> for POIs
     fDiffFlowCorrelationsPro[1][pe][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],two1n1nPtEta,mWeight2pPrime);
     // fill profile to get <<2'>^2> for POIs
     fDiffFlowSquaredCorrelationsPro[1][pe][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],two1n1nPtEta*two1n1nPtEta,mWeight2pPrime);   
     // histogram to store <2'> for POIs e-b-e (needed in some other methods):
     fDiffFlowCorrelationsEBE[1][pe][0]->SetBinContent(b,two1n1nPtEta);      
     fDiffFlowEventWeightsForCorrelationsEBE[1][pe][0]->SetBinContent(b,mWeight2pPrime);      
    }
    else if(type == "RP") // to be improved (I do not this if)
    {
     // profile to get <<2'>> for RPs:
     fDiffFlowCorrelationsPro[0][pe][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],two1n1nPtEta,mWeight2pPrime);     
     // profile to get <<2'>^2> for RPs:
     fDiffFlowSquaredCorrelationsPro[0][pe][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],two1n1nPtEta*two1n1nPtEta,mWeight2pPrime);          
     // histogram to store <2'> for RPs e-b-e (needed in some other methods):
     fDiffFlowCorrelationsEBE[0][pe][0]->SetBinContent(b,two1n1nPtEta); 
     fDiffFlowEventWeightsForCorrelationsEBE[0][pe][0]->SetBinContent(b,mWeight2pPrime); 
    }
   } // end of if(mp*dMult-mq)
  
   // 4'-particle correlation:
   Double_t four1n1n1n1nPtEta = 0.;
   Double_t mWeight4pPrime = 0.; // multiplicity weight for <4'>
   if((mp-mq)*dMult*(dMult-1.)*(dMult-2.)
       + mq*(dMult-1.)*(dMult-2.)*(dMult-3.)) // to be improved (introduce a new variable for this expression)
   {
    four1n1n1n1nPtEta = ((pow(dReQ1n,2.)+pow(dImQ1n,2.))*(p1n0kRe*dReQ1n+p1n0kIm*dImQ1n)
                      - q2n0kRe*(pow(dReQ1n,2.)-pow(dImQ1n,2.))
                      - 2.*q2n0kIm*dReQ1n*dImQ1n
                      - p1n0kRe*(dReQ1n*dReQ2n+dImQ1n*dImQ2n)
                      + p1n0kIm*(dImQ1n*dReQ2n-dReQ1n*dImQ2n)
                      - 2.*dMult*(p1n0kRe*dReQ1n+p1n0kIm*dImQ1n)
                      - 2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*mq                      
                      + 6.*(q1n0kRe*dReQ1n+q1n0kIm*dImQ1n)                                            
                      + 1.*(q2n0kRe*dReQ2n+q2n0kIm*dImQ2n)                      
                      + 2.*(p1n0kRe*dReQ1n+p1n0kIm*dImQ1n)                       
                      + 2.*mq*dMult                      
                      - 6.*mq)        
                      / ((mp-mq)*dMult*(dMult-1.)*(dMult-2.)
                          + mq*(dMult-1.)*(dMult-2.)*(dMult-3.)); 
    // determine multiplicity weight:
    if(!strcmp(fMultiplicityWeight->Data(),"combinations"))
    {
     mWeight4pPrime = (mp-mq)*dMult*(dMult-1.)*(dMult-2.) + mq*(dMult-1.)*(dMult-2.)*(dMult-3.);
    } else if(!strcmp(fMultiplicityWeight->Data(),"unit"))
      {
       mWeight4pPrime = 1.;    
      }     
    if(type == "POI")
    {
     // profile to get <<4'>> for POIs:
     fDiffFlowCorrelationsPro[1][pe][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],four1n1n1n1nPtEta,mWeight4pPrime);      
     // profile to get <<4'>^2> for POIs:
     fDiffFlowSquaredCorrelationsPro[1][pe][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],four1n1n1n1nPtEta*four1n1n1n1nPtEta,mWeight4pPrime); 
     // histogram to store <4'> for POIs e-b-e (needed in some other methods):
     fDiffFlowCorrelationsEBE[1][pe][1]->SetBinContent(b,four1n1n1n1nPtEta);                               
     fDiffFlowEventWeightsForCorrelationsEBE[1][pe][1]->SetBinContent(b,mWeight4pPrime);                               
    }
    else if(type == "RP")
    {
     // profile to get <<4'>> for RPs:
     fDiffFlowCorrelationsPro[0][pe][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],four1n1n1n1nPtEta,mWeight4pPrime);    
     // profile to get <<4'>^2> for RPs:
     fDiffFlowSquaredCorrelationsPro[0][pe][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],four1n1n1n1nPtEta*four1n1n1n1nPtEta,mWeight4pPrime);    
     // histogram to store <4'> for RPs e-b-e (needed in some other methods):
     fDiffFlowCorrelationsEBE[0][pe][1]->SetBinContent(b,four1n1n1n1nPtEta);                   
     fDiffFlowEventWeightsForCorrelationsEBE[0][pe][1]->SetBinContent(b,mWeight4pPrime);                   
    }
   } // end of if((mp-mq)*dMult*(dMult-1.)*(dMult-2.)
     //            +mq*(dMult-1.)*(dMult-2.)*(dMult-3.))
   
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 
   
} // end of void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrelations(TString type, TString ptOrEta);

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateOtherDiffCorrelators(TString type, TString ptOrEta)
{
 // Calculate other differential correlators for RPs or POIs for all pt and eta bins.
 
 // Multiplicity:
 Double_t dMult = (*fSpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);

 // Other correlations are stored in fOtherDiffCorrelators[2][2][2][1], [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correlator index]
 // Correlation index runs as follows:
 // 
 //  0: <exp[in(psi1-3phi2+2phi3)]>
 
 Int_t t = 0; // type flag 
 Int_t pe = 0; // ptEta flag
 
 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }

 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   }
    
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};

 // looping over all bins and calculating reduced correlations: 
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // real and imaginary parts of p_{m*n,0} (non-weighted Q-vector evaluated for POIs in particular pt or eta bin): 
  Double_t p1n0kRe = 0.;
  Double_t p1n0kIm = 0.;

  // number of POIs in particular pt or eta bin:
  Double_t mp = 0.;

  // real and imaginary parts of q_{m*n,0} (non-weighted Q-vector evaluated for particles which are both RPs and POIs in particular pt or eta bin):
  Double_t q1n0kRe = 0.;
  Double_t q1n0kIm = 0.;
  Double_t q2n0kRe = 0.;
  Double_t q2n0kIm = 0.;
  Double_t q3n0kRe = 0.;
  Double_t q3n0kIm = 0.;

  // number of particles which are both RPs and POIs in particular pt or eta bin:
  Double_t mq = 0.;
   
  if(type == "POI")
  {
   // q_{m*n,0}:
   q1n0kRe = fReRPQ1dEBE[2][pe][0][0]->GetBinContent(fReRPQ1dEBE[2][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(fReRPQ1dEBE[2][pe][0][0]->GetBin(b));
   q1n0kIm = fImRPQ1dEBE[2][pe][0][0]->GetBinContent(fImRPQ1dEBE[2][pe][0][0]->GetBin(b))
           * fImRPQ1dEBE[2][pe][0][0]->GetBinEntries(fImRPQ1dEBE[2][pe][0][0]->GetBin(b));
   q2n0kRe = fReRPQ1dEBE[2][pe][1][0]->GetBinContent(fReRPQ1dEBE[2][pe][1][0]->GetBin(b))
           * fReRPQ1dEBE[2][pe][1][0]->GetBinEntries(fReRPQ1dEBE[2][pe][1][0]->GetBin(b));
   q2n0kIm = fImRPQ1dEBE[2][pe][1][0]->GetBinContent(fImRPQ1dEBE[2][pe][1][0]->GetBin(b))
           * fImRPQ1dEBE[2][pe][1][0]->GetBinEntries(fImRPQ1dEBE[2][pe][1][0]->GetBin(b));                         
   q3n0kRe = fReRPQ1dEBE[2][pe][2][0]->GetBinContent(fReRPQ1dEBE[2][pe][2][0]->GetBin(b))
           * fReRPQ1dEBE[2][pe][2][0]->GetBinEntries(fReRPQ1dEBE[2][pe][2][0]->GetBin(b));
   q3n0kIm = fImRPQ1dEBE[2][pe][2][0]->GetBinContent(fImRPQ1dEBE[2][pe][2][0]->GetBin(b))
           * fImRPQ1dEBE[2][pe][2][0]->GetBinEntries(fImRPQ1dEBE[2][pe][2][0]->GetBin(b));         

   mq = fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(fReRPQ1dEBE[2][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
  } 
  else if(type == "RP")
  {
   // q_{m*n,0}:
   q1n0kRe = fReRPQ1dEBE[0][pe][0][0]->GetBinContent(fReRPQ1dEBE[0][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b));
   q1n0kIm = fImRPQ1dEBE[0][pe][0][0]->GetBinContent(fImRPQ1dEBE[0][pe][0][0]->GetBin(b))
           * fImRPQ1dEBE[0][pe][0][0]->GetBinEntries(fImRPQ1dEBE[0][pe][0][0]->GetBin(b));
   q2n0kRe = fReRPQ1dEBE[0][pe][1][0]->GetBinContent(fReRPQ1dEBE[0][pe][1][0]->GetBin(b))
           * fReRPQ1dEBE[0][pe][1][0]->GetBinEntries(fReRPQ1dEBE[0][pe][1][0]->GetBin(b));
   q2n0kIm = fImRPQ1dEBE[0][pe][1][0]->GetBinContent(fImRPQ1dEBE[0][pe][1][0]->GetBin(b))
           * fImRPQ1dEBE[0][pe][1][0]->GetBinEntries(fImRPQ1dEBE[0][pe][1][0]->GetBin(b));         
   q3n0kRe = fReRPQ1dEBE[0][pe][2][0]->GetBinContent(fReRPQ1dEBE[0][pe][2][0]->GetBin(b))
           * fReRPQ1dEBE[0][pe][2][0]->GetBinEntries(fReRPQ1dEBE[0][pe][2][0]->GetBin(b));
   q3n0kIm = fImRPQ1dEBE[0][pe][2][0]->GetBinContent(fImRPQ1dEBE[0][pe][2][0]->GetBin(b))
           * fImRPQ1dEBE[0][pe][2][0]->GetBinEntries(fImRPQ1dEBE[0][pe][2][0]->GetBin(b));         
                 
   mq = fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)  
  }
      
   if(type == "POI")
   {
    // p_{m*n,0}:
    p1n0kRe = fReRPQ1dEBE[1][pe][0][0]->GetBinContent(fReRPQ1dEBE[1][pe][0][0]->GetBin(b))
            * fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b));
    p1n0kIm = fImRPQ1dEBE[1][pe][0][0]->GetBinContent(fImRPQ1dEBE[1][pe][0][0]->GetBin(b))  
            * fImRPQ1dEBE[1][pe][0][0]->GetBinEntries(fImRPQ1dEBE[1][pe][0][0]->GetBin(b));
            
    mp = fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
    
    t = 1; // typeFlag = RP or POI
   }
   else if(type == "RP")
   {
    // p_{m*n,0} = q_{m*n,0}:
    p1n0kRe = q1n0kRe; 
    p1n0kIm = q1n0kIm; 
            
    mp = mq; 
    
    t = 0; // typeFlag = RP or POI
   }
      
   // 3'-particle correlators:
   //  Taeney-Yan correlator:
   Double_t dTaeneyYan = 0.;
   Double_t mWeightTaeneyYan = 0.; // multiplicity weight for Taeney-Yan correlator
   if((mp*dMult-2.*mq)*(dMult-1.) > 0.) // to be improved - is this condition fully justified?
   {
    dTaeneyYan = (dReQ3n*(p1n0kRe*dReQ2n-p1n0kIm*dImQ2n)+dImQ3n*(p1n0kIm*dReQ2n+p1n0kRe*dImQ2n)
               - p1n0kRe*dReQ1n - p1n0kIm*dImQ1n
               - q2n0kRe*dReQ2n - q2n0kIm*dImQ2n              
               - q3n0kRe*dReQ3n - q3n0kIm*dImQ3n
               + 2.*mq)
               / ((mp*dMult-2.*mq)*(dMult-1.));
    // determine multiplicity weight:
    if(!strcmp(fMultiplicityWeight->Data(),"combinations"))
    {
     mWeightTaeneyYan = (mp*dMult-2.*mq)*(dMult-1.);
    } else if(!strcmp(fMultiplicityWeight->Data(),"unit"))
      {
       mWeightTaeneyYan = 1.;    
      } 
    // Fill profiles:
    fOtherDiffCorrelators[t][pe][1][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dTaeneyYan,mWeightTaeneyYan);
   } // end of if((mp*dMult-2.*mq)*(dMult-1.) > 0.)  
   
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 
} // end of void AliFlowAnalysisWithQCumulants::CalculateOtherDiffCorrelators(TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::Calculate2DDiffFlowCorrelations(TString type)
{
 // Calculate all reduced correlations needed for 2D differential flow for each (pt,eta) bin. 
 
 // Multiplicity:
 Double_t dMult = (*fSpk)(0,0);
 // Real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 //Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 //Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);

 // 2D reduced correlations are stored in TProfile2D f2DDiffFlowCorrelationsPro[0=RP,1=POI][correlation index]. 
 // Correlation index runs as follows:
 //  0: <<2'>> 
 //  1: <<4'>>
 //  2: <<6'>>
 //  3: <<8'>>
 
 Int_t t = 0; // type flag  
 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }

 // Looping over all (pt,eta) bins and calculating correlations needed for differential flow: 
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   // Real and imaginary parts of p_{m*n,0} (non-weighted Q-vector evaluated for POIs in particular (pt,eta) bin): 
   Double_t p1n0kRe = 0.;
   Double_t p1n0kIm = 0.;
   // Number of POIs in particular pt or eta bin:
   Double_t mp = 0.;
   // Real and imaginary parts of q_{m*n,0} (non-weighted Q-vector evaluated for 'RP && POI particles' in particular pt or eta bin):
   Double_t q1n0kRe = 0.;
   Double_t q1n0kIm = 0.;
   Double_t q2n0kRe = 0.;
   Double_t q2n0kIm = 0.; 
   // Number of 'RP && POI particles' in particular pt or eta bin:
   Double_t mq = 0.;
   if(type == "POI")
   {
    // q_{m*n,0}:
    q1n0kRe = fReRPQ2dEBE[2][0][0]->GetBinContent(fReRPQ2dEBE[2][0][0]->GetBin(p,e))
            * fReRPQ2dEBE[2][0][0]->GetBinEntries(fReRPQ2dEBE[2][0][0]->GetBin(p,e));
    q1n0kIm = fImRPQ2dEBE[2][0][0]->GetBinContent(fImRPQ2dEBE[2][0][0]->GetBin(p,e))
            * fImRPQ2dEBE[2][0][0]->GetBinEntries(fImRPQ2dEBE[2][0][0]->GetBin(p,e));
    q2n0kRe = fReRPQ2dEBE[2][1][0]->GetBinContent(fReRPQ2dEBE[2][1][0]->GetBin(p,e))
            * fReRPQ2dEBE[2][1][0]->GetBinEntries(fReRPQ2dEBE[2][1][0]->GetBin(p,e));
    q2n0kIm = fImRPQ2dEBE[2][1][0]->GetBinContent(fImRPQ2dEBE[2][1][0]->GetBin(p,e))
            * fImRPQ2dEBE[2][1][0]->GetBinEntries(fImRPQ2dEBE[2][1][0]->GetBin(p,e));         
    // m_{q}:             
    mq = fReRPQ2dEBE[2][0][0]->GetBinEntries(fReRPQ2dEBE[2][0][0]->GetBin(p,e)); // to be improved (cross-checked by accessing other profiles here)
   } // end of if(type == "POI")
   else if(type == "RP")
   {
    // q_{m*n,0}:
    q1n0kRe = fReRPQ2dEBE[0][0][0]->GetBinContent(fReRPQ2dEBE[0][0][0]->GetBin(p,e))
            * fReRPQ2dEBE[0][0][0]->GetBinEntries(fReRPQ2dEBE[0][0][0]->GetBin(p,e));
    q1n0kIm = fImRPQ2dEBE[0][0][0]->GetBinContent(fImRPQ2dEBE[0][0][0]->GetBin(p,e))
            * fImRPQ2dEBE[0][0][0]->GetBinEntries(fImRPQ2dEBE[0][0][0]->GetBin(p,e));
    q2n0kRe = fReRPQ2dEBE[0][1][0]->GetBinContent(fReRPQ2dEBE[0][1][0]->GetBin(p,e))
            * fReRPQ2dEBE[0][1][0]->GetBinEntries(fReRPQ2dEBE[0][1][0]->GetBin(p,e));
    q2n0kIm = fImRPQ2dEBE[0][1][0]->GetBinContent(fImRPQ2dEBE[0][1][0]->GetBin(p,e))
            * fImRPQ2dEBE[0][1][0]->GetBinEntries(fImRPQ2dEBE[0][1][0]->GetBin(p,e));         
    // m_{q}:             
    mq = fReRPQ2dEBE[0][0][0]->GetBinEntries(fReRPQ2dEBE[0][0][0]->GetBin(p,e)); // to be improved (cross-checked by accessing other profiles here)  
   } // end of else if(type == "RP")
   if(type == "POI")
   {
    // p_{m*n,0}:
    p1n0kRe = fReRPQ2dEBE[1][0][0]->GetBinContent(fReRPQ2dEBE[1][0][0]->GetBin(p,e))
            * fReRPQ2dEBE[1][0][0]->GetBinEntries(fReRPQ2dEBE[1][0][0]->GetBin(p,e));
    p1n0kIm = fImRPQ2dEBE[1][0][0]->GetBinContent(fImRPQ2dEBE[1][0][0]->GetBin(p,e))  
            * fImRPQ2dEBE[1][0][0]->GetBinEntries(fImRPQ2dEBE[1][0][0]->GetBin(p,e));
    // m_{p}        
    mp = fReRPQ2dEBE[1][0][0]->GetBinEntries(fReRPQ2dEBE[1][0][0]->GetBin(p,e)); // to be improved (cross-checked by accessing other profiles here)
    
    t = 1; // typeFlag = RP or POI
   } // end of if(type == "POI")
   else if(type == "RP")
   {
    // p_{m*n,0} = q_{m*n,0}:
    p1n0kRe = q1n0kRe; 
    p1n0kIm = q1n0kIm; 
    // m_{p} = m_{q}:        
    mp = mq; 

    t = 0; // typeFlag = RP or POI
   } // end of if(type == "RP")

   // 2'-particle correlation for particular (pt,eta) bin:
   Double_t two1n1nPtEta = 0.;
   Double_t mWeight2pPrime = 0.; // multiplicity weight for <2'>
   if(mp*dMult-mq)
   {
    two1n1nPtEta = (p1n0kRe*dReQ1n+p1n0kIm*dImQ1n-mq)
                 / (mp*dMult-mq);
    // Determine multiplicity weight:
    if(!strcmp(fMultiplicityWeight->Data(),"combinations"))
    {
     mWeight2pPrime = mp*dMult-mq;
    } else if(!strcmp(fMultiplicityWeight->Data(),"unit"))
      {
       mWeight2pPrime = 1.;    
      } 
    // Fill 2D profile holding <<2'>>:     
    f2DDiffFlowCorrelationsPro[t][0]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,two1n1nPtEta,mWeight2pPrime);
   } // end of if(mp*dMult-mq)
   
   // 4'-particle correlation:
   Double_t four1n1n1n1nPtEta = 0.;
   Double_t mWeight4pPrime = 0.; // multiplicity weight for <4'>
   if((mp-mq)*dMult*(dMult-1.)*(dMult-2.)
       + mq*(dMult-1.)*(dMult-2.)*(dMult-3.)) // to be improved (introduce a new variable for this expression)
   {
    four1n1n1n1nPtEta = ((pow(dReQ1n,2.)+pow(dImQ1n,2.))*(p1n0kRe*dReQ1n+p1n0kIm*dImQ1n)
                      - q2n0kRe*(pow(dReQ1n,2.)-pow(dImQ1n,2.))
                      - 2.*q2n0kIm*dReQ1n*dImQ1n
                      - p1n0kRe*(dReQ1n*dReQ2n+dImQ1n*dImQ2n)
                      + p1n0kIm*(dImQ1n*dReQ2n-dReQ1n*dImQ2n)
                      - 2.*dMult*(p1n0kRe*dReQ1n+p1n0kIm*dImQ1n)
                      - 2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*mq                      
                      + 6.*(q1n0kRe*dReQ1n+q1n0kIm*dImQ1n)                                            
                      + 1.*(q2n0kRe*dReQ2n+q2n0kIm*dImQ2n)                      
                      + 2.*(p1n0kRe*dReQ1n+p1n0kIm*dImQ1n)                       
                      + 2.*mq*dMult                      
                      - 6.*mq)        
                      / ((mp-mq)*dMult*(dMult-1.)*(dMult-2.)
                          + mq*(dMult-1.)*(dMult-2.)*(dMult-3.)); 
    // Determine multiplicity weight:
    if(!strcmp(fMultiplicityWeight->Data(),"combinations"))
    {
     mWeight4pPrime = (mp-mq)*dMult*(dMult-1.)*(dMult-2.) + mq*(dMult-1.)*(dMult-2.)*(dMult-3.);
    } else if(!strcmp(fMultiplicityWeight->Data(),"unit"))
      {
       mWeight4pPrime = 1.;    
      }     
    // Fill 2D profile holding <<4'>>:
    f2DDiffFlowCorrelationsPro[t][1]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nPtEta,mWeight4pPrime);      
   } // end of if((mp-mq)*dMult*(dMult-1.)*(dMult-2.)
     //            +mq*(dMult-1.)*(dMult-2.)*(dMult-3.))
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)   
      
} // end of AliFlowAnalysisWithQCumulants::Calculate2DDiffFlowCorrelations(TString type)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateDiffFlowSumOfEventWeights(TString type, TString ptOrEta)
{
 // Calculate sums of various event weights for reduced correlations. 
 // (These quantitites are needed in expressions for unbiased estimators relevant for the statistical errors.)

 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;

 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   } 
     
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
   
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
 
 // binning:
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};
 
 for(Int_t rpq=0;rpq<3;rpq++)
 {
  for(Int_t m=0;m<4;m++)
  {
   for(Int_t k=0;k<9;k++)
   {
    if(!fReRPQ1dEBE[rpq][pe][m][k])
    {
     cout<<"WARNING: fReRPQ1dEBE[rpq][pe][m][k] is NULL in AFAWQC::CSAPOEWFDF() !!!!"<<endl;
     cout<<"pe  = "<<pe<<endl;
     cout<<"rpq = "<<rpq<<endl;
     cout<<"m   = "<<m<<endl;
     cout<<"k   = "<<k<<endl;
     exit(0); 
    }
   }
  }
 }  

 // multiplicities:
 Double_t dMult = (*fSpk)(0,0); // total event multiplicity
 //Double_t mr = 0.; // number of RPs in particular pt or eta bin
 Double_t mp = 0.; // number of POIs in particular pt or eta bin 
 Double_t mq = 0.; // number of particles which are both RPs and POIs in particular pt or eta bin
 
 // event weights for reduced correlations:
 Double_t dw2 = 0.; // event weight for <2'>
 Double_t dw4 = 0.; // event weight for <4'>
 //Double_t dw6 = 0.; // event weight for <6'>
 //Double_t dw8 = 0.; // event weight for <8'>

 // looping over bins:
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  if(type == "RP")
  {
   mq = fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(b);
   mp = mq; // trick to use the very same Eqs. bellow both for RP's and POI's diff. flow
  } else if(type == "POI")
    {
     mp = fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(b);
     mq = fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(b);    
    }
  
  // event weight for <2'>:
  dw2 = mp*dMult-mq;  
  fDiffFlowSumOfEventWeights[t][pe][0][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw2);
  fDiffFlowSumOfEventWeights[t][pe][1][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],pow(dw2,2.));
  
  // event weight for <4'>:
  dw4 = (mp-mq)*dMult*(dMult-1.)*(dMult-2.)
     + mq*(dMult-1.)*(dMult-2.)*(dMult-3.);  
  fDiffFlowSumOfEventWeights[t][pe][0][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw4);
  fDiffFlowSumOfEventWeights[t][pe][1][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],pow(dw4,2.));
  
  // event weight for <6'>:
  //dw6 = ...;  
  //fDiffFlowSumOfEventWeights[t][pe][0][2]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw6);
  //fDiffFlowSumOfEventWeights[t][pe][t][1][2]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],pow(dw6,2.));
  
  // event weight for <8'>:
  //dw8 = ...;  
  //fDiffFlowSumOfEventWeights[t][pe][0][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw8);
  //fDiffFlowSumOfEventWeights[t][pe][1][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],pow(dw8,2.));   
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++) 
 
} // end of void AliFlowAnalysisWithQCumulants::CalculateDiffFlowSumOfEventWeights()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateDiffFlowSumOfProductOfEventWeights(TString type, TString ptOrEta)
{
 // Calculate sum of products of various event weights for both types of correlations (the ones for int. and diff. flow). 
 // (These quantitites are needed in expressions for unbiased estimators relevant for the statistical errors.)
 //
 // Important: To fill fDiffFlowSumOfProductOfEventWeights[][][][] use bellow table (i,j) with following constraints: 
 // 1.) i<j  
 // 2.) do not store terms which DO NOT include reduced correlations;
 // Table:
 // [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>] x [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>]
  
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;

 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   } 
     
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
     
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
  
 // binning:
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};
 
 // protection:
 for(Int_t rpq=0;rpq<3;rpq++)
 {
  for(Int_t m=0;m<4;m++)
  {
   for(Int_t k=0;k<9;k++)
   {
    if(!fReRPQ1dEBE[rpq][pe][m][k])
    {
     cout<<"WARNING: fReRPQ1dEBE[rpq][pe][m][k] is NULL in AFAWQC::CSAPOEWFDF() !!!!"<<endl;
     cout<<"pe  = "<<pe<<endl;
     cout<<"rpq = "<<rpq<<endl;
     cout<<"m   = "<<m<<endl;
     cout<<"k   = "<<k<<endl;
     exit(0); 
    }
   }
  }
 }  
 
 // multiplicities:
 Double_t dMult = (*fSpk)(0,0); // total event multiplicity
 //Double_t mr = 0.; // number of RPs in particular pt or eta bin
 Double_t mp = 0.; // number of POIs in particular pt or eta bin 
 Double_t mq = 0.; // number of particles which are both RPs and POIs in particular pt or eta bin
 
 // event weights for correlations:
 Double_t dW2 = dMult*(dMult-1); // event weight for <2> 
 Double_t dW4 = dMult*(dMult-1)*(dMult-2)*(dMult-3); // event weight for <4> 
 Double_t dW6 = dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5); // event weight for <6> 
 Double_t dW8 = dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5)*(dMult-6)*(dMult-7); // event weight for <8> 

 // event weights for reduced correlations:
 Double_t dw2 = 0.; // event weight for <2'>
 Double_t dw4 = 0.; // event weight for <4'>
 //Double_t dw6 = 0.; // event weight for <6'>
 //Double_t dw8 = 0.; // event weight for <8'>
 
 // looping over bins:
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  if(type == "RP")
  {
   mq = fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(b);
   mp = mq; // trick to use the very same Eqs. bellow both for RP's and POI's diff. flow
  } else if(type == "POI")
    {
     mp = fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(b);
     mq = fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(b);    
    }
  
  // event weight for <2'>:
  dw2 = mp*dMult-mq;  
  fDiffFlowSumOfProductOfEventWeights[t][pe][0][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dW2*dw2); // storing product of even weights for <2> and <2'>
  fDiffFlowSumOfProductOfEventWeights[t][pe][1][2]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw2*dW4); // storing product of even weights for <4> and <2'>
  fDiffFlowSumOfProductOfEventWeights[t][pe][1][4]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw2*dW6); // storing product of even weights for <6> and <2'>
  fDiffFlowSumOfProductOfEventWeights[t][pe][1][6]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw2*dW8); // storing product of even weights for <8> and <2'>
  
  // event weight for <4'>:
  dw4 = (mp-mq)*dMult*(dMult-1.)*(dMult-2.)
     + mq*(dMult-1.)*(dMult-2.)*(dMult-3.);  
  fDiffFlowSumOfProductOfEventWeights[t][pe][0][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dW2*dw4); // storing product of even weights for <2> and <4'>
  fDiffFlowSumOfProductOfEventWeights[t][pe][1][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw2*dw4); // storing product of even weights for <2'> and <4'>
  fDiffFlowSumOfProductOfEventWeights[t][pe][2][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dW4*dw4); // storing product of even weights for <4> and <4'>
  fDiffFlowSumOfProductOfEventWeights[t][pe][3][4]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw4*dW6); // storing product of even weights for <6> and <4'> 
  fDiffFlowSumOfProductOfEventWeights[t][pe][3][6]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw4*dW8); // storing product of even weights for <8> and <4'>

  // event weight for <6'>:
  //dw6 = ...;  
  //fDiffFlowSumOfProductOfEventWeights[t][pe][0][5]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dW2*dw6); // storing product of even weights for <2> and <6'>
  //fDiffFlowSumOfProductOfEventWeights[t][pe][1][5]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw2*dw6); // storing product of even weights for <2'> and <6'>
  //fDiffFlowSumOfProductOfEventWeights[t][pe][2][5]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dW4*dw6); // storing product of even weights for <4> and <6'>
  //fDiffFlowSumOfProductOfEventWeights[t][pe][3][5]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw4*dw6); // storing product of even weights for <4'> and <6'> 
  //fDiffFlowSumOfProductOfEventWeights[t][pe][4][5]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dW6*dw6); // storing product of even weights for <6> and <6'>
  //fDiffFlowSumOfProductOfEventWeights[t][pe][5][6]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw6*dW8); // storing product of even weights for <6'> and <8>
  //fDiffFlowSumOfProductOfEventWeights[t][pe][5][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw6*dw8); // storing product of even weights for <6'> and <8'>

  // event weight for <8'>:
  //dw8 = ...;  
  //fDiffFlowSumOfProductOfEventWeights[t][pe][0][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dW2*dw8); // storing product of even weights for <2> and <8'>
  //fDiffFlowSumOfProductOfEventWeights[t][pe][1][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw2*dw8); // storing product of even weights for <2'> and <8'>
  //fDiffFlowSumOfProductOfEventWeights[t][pe][2][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dW4*dw8); // storing product of even weights for <4> and <8'>
  //fDiffFlowSumOfProductOfEventWeights[t][pe][3][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw4*dw8); // storing product of even weights for <4'> and <8'> 
  //fDiffFlowSumOfProductOfEventWeights[t][pe][4][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dW6*dw8); // storing product of even weights for <6> and <8'>
  //fDiffFlowSumOfProductOfEventWeights[t][pe][5][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dw6*dw8); // storing product of even weights for <6'> and <8'>
  //fDiffFlowSumOfProductOfEventWeights[t][pe][6][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],dW8*dw8); // storing product of even weights for <8> and <8'>
  
  // Table:
  // [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>] x [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>]
   
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 


} // end of void AliFlowAnalysisWithQCumulants::CalculateDiffFlowSumOfProductOfEventWeights(TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::FinalizeReducedCorrelations(TString type, TString ptOrEta)
{
 // Transfer profiles into histograms and calculate statistical errors correctly.

 Int_t t = 0; // RP or POI
 Int_t pe = 0; // pt or eta

 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   } 
     
 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   } 
               
 for(Int_t rci=0;rci<4;rci++) // to be improved - moved into the method CheckPointersUsedInFinish()
 {
  if(!fDiffFlowCorrelationsPro[t][pe][rci])
  {
   cout<<"WARNING: fDiffFlowCorrelationsPro[t][pe][rci] is NULL in AFAWQC::FRC() !!!!"<<endl;
   cout<<"t   = "<<t<<endl; 
   cout<<"pe  = "<<pe<<endl; 
   cout<<"rci = "<<rci<<endl;
   exit(0); 
  }
  if(!fDiffFlowSquaredCorrelationsPro[t][pe][rci])
  {
   cout<<"WARNING: fDiffFlowSquaredCorrelationsPro[t][pe][rci] is NULL in AFAWQC::FRC() !!!!"<<endl;
   cout<<"t   = "<<t<<endl; 
   cout<<"pe  = "<<pe<<endl; 
   cout<<"rci = "<<rci<<endl;
   exit(0); 
  }
  for(Int_t power=0;power<2;power++)
  {
   if(!fDiffFlowSumOfEventWeights[t][pe][power][rci])
   {
    cout<<"WARNING: fDiffFlowSumOfEventWeights[t][pe][power][rci] is NULL in AFAWQC::FRC() !!!!"<<endl;
    cout<<"t     = "<<t<<endl; 
    cout<<"pe    = "<<pe<<endl;
    cout<<"power = "<<power<<endl; 
    cout<<"rci   = "<<rci<<endl;
    exit(0); 
   }   
  } // end of for(Int_t power=0;power<2;power++)
 } // end of for(Int_t rci=0;rci<4;rci++)
    
 // common:
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta}; 
 // transfer 1D profile into 1D histogram:
 Double_t correlation = 0.;
 Double_t squaredCorrelation = 0.;
 Double_t spread = 0.;
 Double_t sumOfWeights = 0.; // sum of weights for particular reduced correlations for particular pt or eta bin
 Double_t sumOfSquaredWeights = 0.; // sum of squared weights for particular reduced correlations for particular pt or eta bin
 Double_t error = 0.; // error = termA * spread * termB
                      // termA = (sqrt(sumOfSquaredWeights)/sumOfWeights) 
                      // termB = 1/pow(1-termA^2,0.5)
 Double_t termA = 0.;                      
 Double_t termB = 0.;                      
 for(Int_t rci=0;rci<4;rci++) // index of reduced correlation
 {
  for(Int_t b=1;b<=nBinsPtEta[pe];b++) // number of pt or eta bins
  {
   if(fDiffFlowCorrelationsPro[t][pe][rci]->GetBinEffectiveEntries(b) < 2 || 
      fDiffFlowSquaredCorrelationsPro[t][pe][rci]->GetBinEffectiveEntries(b) < 2)
   {
    fDiffFlowCorrelationsPro[t][pe][rci]->SetBinError(b,0.);
    fDiffFlowSquaredCorrelationsPro[t][pe][rci]->SetBinError(b,0.);
    continue; // to be improved - should I ignore results in pt bins with one entry for reduced correlations or not?
   }  
   correlation = fDiffFlowCorrelationsPro[t][pe][rci]->GetBinContent(b); 
   squaredCorrelation = fDiffFlowSquaredCorrelationsPro[t][pe][rci]->GetBinContent(b); 
   if(squaredCorrelation-correlation*correlation >= 0.)
   {
    spread = pow(squaredCorrelation-correlation*correlation,0.5);
   } else
     {
      cout<<endl;
      cout<<Form(" WARNING: Imaginary 'spread' for rci = %d, pe = %d, bin = %d !!!!",rci,pe,b)<<endl;
      cout<<endl;
     }
   sumOfWeights = fDiffFlowSumOfEventWeights[t][pe][0][rci]->GetBinContent(b);
   sumOfSquaredWeights = fDiffFlowSumOfEventWeights[t][pe][1][rci]->GetBinContent(b);
   if(TMath::Abs(sumOfWeights)>0.){termA = (pow(sumOfSquaredWeights,0.5)/sumOfWeights);}
   if(1.-pow(termA,2.)>0.){termB = 1./pow(1.-pow(termA,2.),0.5);} 
   error = termA*spread*termB; // final error (unbiased estimator for standard deviation)
   fDiffFlowCorrelationsHist[t][pe][rci]->SetBinContent(b,correlation); 
   fDiffFlowCorrelationsHist[t][pe][rci]->SetBinError(b,error); 
  } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 } // end of for(Int_t rci=0;rci<4;rci++)
 
} // end of void AliFlowAnalysisWithQCumulants::FinalizeReducedCorrelations(TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateDiffFlowProductOfCorrelations(TString type, TString ptOrEta)
{
 // store products: <2><2'>, <2><4'>, <2><6'>, <2><8'>, <2'><4>, 
 //                 <2'><4'>, <2'><6>, <2'><6'>, <2'><8>, <2'><8'>,
 //                 <4><4'>, <4><6'>, <4><8'>, <4'><6>, <4'><6'>, 
 //                 <4'><8>, <4'><8'>, <6><6'>, <6><8'>, <6'><8>, 
 //                 <6'><8'>, <8><8'>.
  
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;

 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   } 
     
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
  
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
     
 // common:
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};
   
 // protections // to be improved (add protection for all pointers in this method)
 if(!fIntFlowCorrelationsEBE)
 {
  cout<<"WARNING: fIntFlowCorrelationsEBE is NULL in AFAWQC::CDFPOC() !!!!"<<endl;
  exit(0);
 } 
 
 /*    
 Double_t dMult = (*fSpk)(0,0); // multiplicity (number of particles used to determine the reaction plane)
 //Double_t mr = 0.; // number of RPs in particular pt or eta bin
 Double_t mp = 0.; // number of POIs in particular pt or eta bin 
 Double_t mq = 0.; // number of particles which are both RPs and POIs in particular pt or eta bin
 */

 // e-b-e correlations:
 Double_t twoEBE = fIntFlowCorrelationsEBE->GetBinContent(1); // <2>
 Double_t fourEBE = fIntFlowCorrelationsEBE->GetBinContent(2); // <4>
 Double_t sixEBE = fIntFlowCorrelationsEBE->GetBinContent(3); // <6>
 Double_t eightEBE = fIntFlowCorrelationsEBE->GetBinContent(4); // <8>
 
 // event weights for correlations:
 Double_t dW2 = fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(1); // event weight for <2> 
 Double_t dW4 = fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(2); // event weight for <4> 
 Double_t dW6 = fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(3); // event weight for <6> 
 Double_t dW8 = fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(4); // event weight for <8> 
  
 // e-b-e reduced correlations:
 Double_t twoReducedEBE = 0.; // <2'>
 Double_t fourReducedEBE = 0.; // <4'>
 Double_t sixReducedEBE = 0.; // <6'>
 Double_t eightReducedEBE = 0.; // <8'> 
 
 // event weights for reduced correlations:
 Double_t dw2 = 0.; // event weight for <2'>
 Double_t dw4 = 0.; // event weight for <4'>
 //Double_t dw6 = 0.; // event weight for <6'>
 //Double_t dw8 = 0.; // event weight for <8'>

 // looping over bins:
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // e-b-e reduced correlations:
  twoReducedEBE = fDiffFlowCorrelationsEBE[t][pe][0]->GetBinContent(b);
  fourReducedEBE = fDiffFlowCorrelationsEBE[t][pe][1]->GetBinContent(b);
  sixReducedEBE = fDiffFlowCorrelationsEBE[t][pe][2]->GetBinContent(b);
  eightReducedEBE = fDiffFlowCorrelationsEBE[t][pe][3]->GetBinContent(b);
  
  /*
  // to be improved (I should not do this here again)
  if(type == "RP")
  {
   mq = fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(b);
   mp = mq; // trick to use the very same Eqs. bellow both for RP's and POI's diff. flow
  } else if(type == "POI")
    {
     mp = fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(b);
     mq = fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(b);    
    }
  
  // event weights for reduced correlations:
  dw2 = mp*dMult-mq; // weight for <2'> 
  dw4 = (mp-mq)*dMult*(dMult-1.)*(dMult-2.)
     + mq*(dMult-1.)*(dMult-2.)*(dMult-3.); // weight for <4'>
  //dw6 = ...     
  //dw8 = ...     
  
  */
  
  dw2 = fDiffFlowEventWeightsForCorrelationsEBE[t][pe][0]->GetBinContent(b);
  dw4 = fDiffFlowEventWeightsForCorrelationsEBE[t][pe][1]->GetBinContent(b);
 
  // storing all products:
  fDiffFlowProductOfCorrelationsPro[t][pe][0][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],twoEBE*twoReducedEBE,dW2*dw2); // storing <2><2'>
  fDiffFlowProductOfCorrelationsPro[t][pe][1][2]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],fourEBE*twoReducedEBE,dW4*dw2); // storing <4><2'>
  fDiffFlowProductOfCorrelationsPro[t][pe][1][4]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sixEBE*twoReducedEBE,dW6*dw2); // storing <6><2'>
  fDiffFlowProductOfCorrelationsPro[t][pe][1][6]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],eightEBE*twoReducedEBE,dW8*dw2); // storing <8><2'>
  
  // event weight for <4'>:
  fDiffFlowProductOfCorrelationsPro[t][pe][0][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],twoEBE*fourReducedEBE,dW2*dw4); // storing <2><4'>
  fDiffFlowProductOfCorrelationsPro[t][pe][1][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],twoReducedEBE*fourReducedEBE,dw2*dw4); // storing <2'><4'>
  fDiffFlowProductOfCorrelationsPro[t][pe][2][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],fourEBE*fourReducedEBE,dW4*dw4); // storing <4><4'>
  fDiffFlowProductOfCorrelationsPro[t][pe][3][4]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sixEBE*fourReducedEBE,dW6*dw4); // storing <6><4'> 
  fDiffFlowProductOfCorrelationsPro[t][pe][3][6]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],eightEBE*fourReducedEBE,dW8*dw4); // storing <8><4'>

  // event weight for <6'>:
  //dw6 = ...;  
  //fDiffFlowProductOfCorrelationsPro[t][pe][0][5]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],twoEBE*sixReducedEBE,dW2*dw6); // storing <2><6'>
  //fDiffFlowProductOfCorrelationsPro[t][pe][1][5]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],twoReducedEBE*sixReducedEBE,dw2*dw6); // storing <2'><6'>
  //fDiffFlowProductOfCorrelationsPro[t][pe][2][5]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],fourEBE*sixReducedEBE,dW4*dw6); // storing <4><6'>
  //fDiffFlowProductOfCorrelationsPro[t][pe][3][5]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],fourReducedEBE*sixReducedEBE,dw4*dw6); // storing <4'><6'> 
  //fDiffFlowProductOfCorrelationsPro[t][pe][4][5]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sixEBE*sixReducedEBE,dW6*dw6); // storing <6><6'>
  //fDiffFlowProductOfCorrelationsPro[t][pe][5][6]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sixReducedEBE*eightEBE,dw6*dW8); // storing <6'><8>
  //fDiffFlowProductOfCorrelationsPro[t][pe][5][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sixReducedEBE*eightReducedEBE,dw6*dw8); // storing <6'><8'>

  // event weight for <8'>:
  //dw8 = ...;  
  //fDiffFlowProductOfCorrelationsPro[t][pe][0][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],twoEBE*eightReducedEBE,dW2*dw8); // storing <2><8'>
  //fDiffFlowProductOfCorrelationsPro[t][pe][1][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],twoReducedEBE*eightReducedEBE,dw2*dw8); // storing <2'><8'>
  //fDiffFlowProductOfCorrelationsPro[t][pe][2][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],fourEBE*eightReducedEBE,dW4*dw8); // storing <4><8'>
  //fDiffFlowProductOfCorrelationsPro[t][pe][3][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],fourReducedEBE*eightReducedEBE,dw4*dw8); // storing <4'><8'> 
  //fDiffFlowProductOfCorrelationsPro[t][pe][4][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sixEBE*eightReducedEBE,dW6*dw8); // storing <6><8'>
  //fDiffFlowProductOfCorrelationsPro[t][pe][5][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sixReducedEBE*eightReducedEBE,dw6*dw8); // storing <6'><8'>
  //fDiffFlowProductOfCorrelationsPro[t][pe][6][7]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],eightEBE*eightReducedEBE,dW8*dw8); // storing <8><8'> 
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++       
     
} // end of void AliFlowAnalysisWithQCumulants::CalculateDiffFlowProductOfCorrelations(TString type, TString ptOrEta)

//================================================================================================================================
    
void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCovariances(TString type, TString ptOrEta) // to be improved (reimplemented)
{
 // a) Calculate unbiased estimators Cov(<2>,<2'>), Cov(<2>,<4'>), Cov(<4>,<2'>), Cov(<4>,<4'>) and Cov(<2'>,<4'>)
 //    for covariances V(<2>,<2'>), V(<2>,<4'>), V(<4>,<2'>), V(<4>,<4'>) and V(<2'>,<4'>).  
 // b) Store in histogram fDiffFlowCovariances[t][pe][index] for instance the following: 
 //
 //             Cov(<2>,<2'>) * (sum_{i=1}^{N} w_{<2>}_i w_{<2'>}_i )/[(sum_{i=1}^{N} w_{<2>}_i) * (sum_{j=1}^{N} w_{<2'>}_j)]
 // 
 //     where N is the number of events, w_{<2>} is event weight for <2> and w_{<2'>} is event weight for <2'>.
 // c) Binning of fDiffFlowCovariances[t][pe][index] is organized as follows:
 // 
 //     1st bin: Cov(<2>,<2'>) * (sum_{i=1}^{N} w_{<2>}_i w_{<2'>}_i )/[(sum_{i=1}^{N} w_{<2>}_i) * (sum_{j=1}^{N} w_{<2'>}_j)] 
 //     2nd bin: Cov(<2>,<4'>) * (sum_{i=1}^{N} w_{<2>}_i w_{<4'>}_i )/[(sum_{i=1}^{N} w_{<2>}_i) * (sum_{j=1}^{N} w_{<4'>}_j)] 
 //     3rd bin: Cov(<4>,<2'>) * (sum_{i=1}^{N} w_{<4>}_i w_{<2'>}_i )/[(sum_{i=1}^{N} w_{<4>}_i) * (sum_{j=1}^{N} w_{<2'>}_j)] 
 //     4th bin: Cov(<4>,<4'>) * (sum_{i=1}^{N} w_{<4>}_i w_{<4'>}_i )/[(sum_{i=1}^{N} w_{<4>}_i) * (sum_{j=1}^{N} w_{<4'>}_j)] 
 //     5th bin: Cov(<2'>,<4'>) * (sum_{i=1}^{N} w_{<2'>}_i w_{<4'>}_i )/[(sum_{i=1}^{N} w_{<2'>}_i) * (sum_{j=1}^{N} w_{<4'>}_j)] 
 //     ...
  
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;

 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   } 
     
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
  
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
     
 // common:
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 //Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 //Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};
 
 // average correlations:
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>>
 Double_t four = fIntFlowCorrelationsHist->GetBinContent(2); // <<4>>
 //Double_t six = fIntFlowCorrelationsHist->GetBinContent(3); // <<6>>
 //Double_t eight = fIntFlowCorrelationsHist->GetBinContent(4); // <<8>>
 
 // sum of weights for correlation:
 Double_t sumOfWeightsForTwo = fIntFlowSumOfEventWeights[0]->GetBinContent(1); // sum_{i=1}^{N} w_{<2>}
 Double_t sumOfWeightsForFour = fIntFlowSumOfEventWeights[0]->GetBinContent(2); // sum_{i=1}^{N} w_{<4>}
 //Double_t sumOfWeightsForSix = fIntFlowSumOfEventWeights[0]->GetBinContent(3); // sum_{i=1}^{N} w_{<6>}
 //Double_t sumOfWeightsForEight = fIntFlowSumOfEventWeights[0]->GetBinContent(4); // sum_{i=1}^{N} w_{<8>}
 
 // average reduced correlations:
 Double_t twoReduced = 0.; // <<2'>> 
 Double_t fourReduced = 0.; // <<4'>>
 //Double_t sixReduced = 0.; // <<6'>>
 //Double_t eightReduced = 0.; // <<8'>>

 // sum of weights for reduced correlation:
 Double_t sumOfWeightsForTwoReduced = 0.; // sum_{i=1}^{N} w_{<2'>}
 Double_t sumOfWeightsForFourReduced = 0.; // sum_{i=1}^{N} w_{<4'>}
 //Double_t sumOfWeightsForSixReduced = 0.; // sum_{i=1}^{N} w_{<6'>}
 //Double_t sumOfWeightsForEightReduced = 0.; // sum_{i=1}^{N} w_{<8'>}
  
 // product of weights for reduced correlation:
 Double_t productOfWeightsForTwoTwoReduced = 0.; // sum_{i=1}^{N} w_{<2>}w_{<2'>}
 Double_t productOfWeightsForTwoFourReduced = 0.; // sum_{i=1}^{N} w_{<2>}w_{<4'>}
 Double_t productOfWeightsForFourTwoReduced = 0.; // sum_{i=1}^{N} w_{<4>}w_{<2'>}
 Double_t productOfWeightsForFourFourReduced = 0.; // sum_{i=1}^{N} w_{<4>}w_{<4'>}
 Double_t productOfWeightsForTwoReducedFourReduced = 0.; // sum_{i=1}^{N} w_{<2'>}w_{<4'>}
 // ...
 
 // products for differential flow:
 Double_t twoTwoReduced = 0; // <<2><2'>> 
 Double_t twoFourReduced = 0; // <<2><4'>> 
 Double_t fourTwoReduced = 0; // <<4><2'>> 
 Double_t fourFourReduced = 0; // <<4><4'>> 
 Double_t twoReducedFourReduced = 0; // <<2'><4'>> 

 // denominators in the expressions for the unbiased estimators for covariances:
 // denominator = 1 - term1/(term2*term3)
 // prefactor = term1/(term2*term3)
 Double_t denominator = 0.; 
 Double_t prefactor = 0.;
 Double_t term1 = 0.; 
 Double_t term2 = 0.; 
 Double_t term3 = 0.; 
 
 // unbiased estimators for covariances for differential flow:
 Double_t covTwoTwoReduced = 0.; // Cov(<2>,<2'>)
 Double_t wCovTwoTwoReduced = 0.; // Cov(<2>,<2'>) * prefactor(w_{<2>},w_{<2'>})
 Double_t covTwoFourReduced = 0.; // Cov(<2>,<4'>)
 Double_t wCovTwoFourReduced = 0.; // Cov(<2>,<4'>) * prefactor(w_{<2>},w_{<4'>})
 Double_t covFourTwoReduced = 0.; // Cov(<4>,<2'>)
 Double_t wCovFourTwoReduced = 0.; // Cov(<4>,<2'>) * prefactor(w_{<4>},w_{<2'>})
 Double_t covFourFourReduced = 0.; // Cov(<4>,<4'>)
 Double_t wCovFourFourReduced = 0.; // Cov(<4>,<4'>) * prefactor(w_{<4>},w_{<4'>})
 Double_t covTwoReducedFourReduced = 0.; // Cov(<2'>,<4'>)
 Double_t wCovTwoReducedFourReduced = 0.; // Cov(<2'>,<4'>) * prefactor(w_{<2'>},w_{<4'>})
 
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // average reduced corelations:
  twoReduced = fDiffFlowCorrelationsHist[t][pe][0]->GetBinContent(b);
  fourReduced = fDiffFlowCorrelationsHist[t][pe][1]->GetBinContent(b);
  // average products:
  twoTwoReduced = fDiffFlowProductOfCorrelationsPro[t][pe][0][1]->GetBinContent(b);
  twoFourReduced = fDiffFlowProductOfCorrelationsPro[t][pe][0][3]->GetBinContent(b);
  fourTwoReduced = fDiffFlowProductOfCorrelationsPro[t][pe][1][2]->GetBinContent(b);
  fourFourReduced = fDiffFlowProductOfCorrelationsPro[t][pe][2][3]->GetBinContent(b);
  twoReducedFourReduced = fDiffFlowProductOfCorrelationsPro[t][pe][1][3]->GetBinContent(b);  
  // sum of weights for reduced correlations:
  sumOfWeightsForTwoReduced = fDiffFlowSumOfEventWeights[t][pe][0][0]->GetBinContent(b);
  sumOfWeightsForFourReduced = fDiffFlowSumOfEventWeights[t][pe][0][1]->GetBinContent(b);
  // products of weights for correlations:
  productOfWeightsForTwoTwoReduced = fDiffFlowSumOfProductOfEventWeights[t][pe][0][1]->GetBinContent(b); 
  productOfWeightsForTwoFourReduced = fDiffFlowSumOfProductOfEventWeights[t][pe][0][3]->GetBinContent(b);
  productOfWeightsForFourTwoReduced = fDiffFlowSumOfProductOfEventWeights[t][pe][1][2]->GetBinContent(b);
  productOfWeightsForFourFourReduced = fDiffFlowSumOfProductOfEventWeights[t][pe][2][3]->GetBinContent(b);
  productOfWeightsForTwoReducedFourReduced = fDiffFlowSumOfProductOfEventWeights[t][pe][1][3]->GetBinContent(b);
  // denominator for the unbiased estimator for covariances: 1 - term1/(term2*term3) 
  // prefactor (multiplies Cov's) = term1/(term2*term3)       
  // <2>,<2'>:
  term1 = productOfWeightsForTwoTwoReduced;      
  term2 = sumOfWeightsForTwo;
  term3 = sumOfWeightsForTwoReduced;        
  if(term2*term3>0.)
  {
   denominator = 1.-term1/(term2*term3);
   prefactor = term1/(term2*term3);
   if(TMath::Abs(denominator)>1.e-6)
   {
    covTwoTwoReduced = (twoTwoReduced-two*twoReduced)/denominator;            
    wCovTwoTwoReduced = covTwoTwoReduced*prefactor; 
    fDiffFlowCovariances[t][pe][0]->SetBinContent(b,wCovTwoTwoReduced);
   }
  }
  // <2>,<4'>:
  term1 = productOfWeightsForTwoFourReduced;      
  term2 = sumOfWeightsForTwo;
  term3 = sumOfWeightsForFourReduced;        
  if(term2*term3>0.)
  {
   denominator = 1.-term1/(term2*term3);
   prefactor = term1/(term2*term3);
   if(TMath::Abs(denominator)>1.e-6)
   {
    covTwoFourReduced = (twoFourReduced-two*fourReduced)/denominator;            
    wCovTwoFourReduced = covTwoFourReduced*prefactor; 
    fDiffFlowCovariances[t][pe][1]->SetBinContent(b,wCovTwoFourReduced);
   }
  }
  // <4>,<2'>:
  term1 = productOfWeightsForFourTwoReduced;      
  term2 = sumOfWeightsForFour;
  term3 = sumOfWeightsForTwoReduced;        
  if(term2*term3>0.)
  {
   denominator = 1.-term1/(term2*term3);
   prefactor = term1/(term2*term3);
   if(TMath::Abs(denominator)>1.e-6)
   {
    covFourTwoReduced = (fourTwoReduced-four*twoReduced)/denominator;            
    wCovFourTwoReduced = covFourTwoReduced*prefactor; 
    fDiffFlowCovariances[t][pe][2]->SetBinContent(b,wCovFourTwoReduced);
   }
  }
  // <4>,<4'>:
  term1 = productOfWeightsForFourFourReduced;      
  term2 = sumOfWeightsForFour;
  term3 = sumOfWeightsForFourReduced;        
  if(term2*term3>0.)
  {
   denominator = 1.-term1/(term2*term3);
   prefactor = term1/(term2*term3);
   if(TMath::Abs(denominator)>1.e-6)
   {
    covFourFourReduced = (fourFourReduced-four*fourReduced)/denominator;            
    wCovFourFourReduced = covFourFourReduced*prefactor; 
    fDiffFlowCovariances[t][pe][3]->SetBinContent(b,wCovFourFourReduced);
   }
  }
  // <2'>,<4'>:
  term1 = productOfWeightsForTwoReducedFourReduced;      
  term2 = sumOfWeightsForTwoReduced;
  term3 = sumOfWeightsForFourReduced;        
  if(term2*term3>0.)
  {
   denominator = 1.-term1/(term2*term3);
   prefactor = term1/(term2*term3);
   if(TMath::Abs(denominator)>1.e-6)
   {
    covTwoReducedFourReduced = (twoReducedFourReduced-twoReduced*fourReduced)/denominator;            
    wCovTwoReducedFourReduced = covTwoReducedFourReduced*prefactor; 
    fDiffFlowCovariances[t][pe][4]->SetBinContent(b,wCovTwoReducedFourReduced);
   }
  }   
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)
  
} // end of void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCovariances(TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateDiffFlow(TString type, TString ptOrEta)
{
 // Calculate final results for differential flow.
 
 // REMARK: Differential flow calculated in this method is NOT corrected for non-uniform acceptance. 
 // This correction, if enabled via setter SetApplyCorrectionForNUA(Bool_t), is applied in the method 
 // CalculateDiffFlowCorrectedForNUA(TString type, TString ptOrEta)
  
 Int_t t = 0; // RP or POI
 Int_t pe = 0; // pt or eta

 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   } 
     
 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   } 
       
 // Common:
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 // Correlations:
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>>
 Double_t four = fIntFlowCorrelationsHist->GetBinContent(2); // <<4>> 
 // Statistical errors of correlations:
 Double_t twoError = fIntFlowCorrelationsHist->GetBinError(1);
 Double_t fourError = fIntFlowCorrelationsHist->GetBinError(2);   
 // Reduced correlations:
 Double_t twoReduced = 0.; // <<2'>>
 Double_t fourReduced = 0.; // <<4'>>
 // Statistical errors of reduced correlations:
 Double_t twoReducedError = 0.; 
 Double_t fourReducedError = 0.; 
 // Covariances:
 Double_t wCovTwoFour = 0.; // Cov(<2>,<4>) * prefactor(<2>,<4>)
 if(!fForgetAboutCovariances)
 {
  wCovTwoFour = fIntFlowCovariances->GetBinContent(1); // Cov(<2>,<4>) * prefactor(<2>,<4>)
 }
 Double_t wCovTwoTwoReduced = 0.; // Cov(<2>,<2'>) * prefactor(<2>,<2'>)
 Double_t wCovTwoFourReduced = 0.; // Cov(<2>,<4'>) * prefactor(<2>,<4'>)
 Double_t wCovFourTwoReduced = 0.; // Cov(<4>,<2'>) * prefactor(<4>,<2'>)
 Double_t wCovFourFourReduced = 0.; // Cov(<4>,<4'>) * prefactor(<4>,<4'>)
 Double_t wCovTwoReducedFourReduced = 0.; // Cov(<2'>,<4'>) * prefactor(<2'>,<4'>)
 // Differential flow:
 Double_t v2Prime = 0.; // v'{2}                   
 Double_t v4Prime = 0.; // v'{4}
 // Statistical error of differential flow:
 Double_t v2PrimeError = 0.;                    
 Double_t v4PrimeError = 0.; 
 // Squared statistical error of differential flow:
 Double_t v2PrimeErrorSquared = 0.;                    
 Double_t v4PrimeErrorSquared = 0.; 
 // Loop over pt or eta bins:
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // Reduced correlations and statistical errors:
  twoReduced = fDiffFlowCorrelationsHist[t][pe][0]->GetBinContent(b);
  twoReducedError = fDiffFlowCorrelationsHist[t][pe][0]->GetBinError(b);
  fourReduced = fDiffFlowCorrelationsHist[t][pe][1]->GetBinContent(b);
  fourReducedError = fDiffFlowCorrelationsHist[t][pe][1]->GetBinError(b);
  // Covariances:
  if(!fForgetAboutCovariances)
  {
   wCovTwoTwoReduced = fDiffFlowCovariances[t][pe][0]->GetBinContent(b);
   wCovTwoFourReduced = fDiffFlowCovariances[t][pe][1]->GetBinContent(b);
   wCovFourTwoReduced = fDiffFlowCovariances[t][pe][2]->GetBinContent(b);
   wCovFourFourReduced = fDiffFlowCovariances[t][pe][3]->GetBinContent(b);
   wCovTwoReducedFourReduced = fDiffFlowCovariances[t][pe][4]->GetBinContent(b);
  }
  // Differential flow:
  // v'{2}:
  if(two>0.) 
  {
   v2Prime = twoReduced/pow(two,0.5);
   v2PrimeErrorSquared = (1./4.)*pow(two,-3.)*(pow(twoReduced,2.)*pow(twoError,2.)
                       + 4.*pow(two,2.)*pow(twoReducedError,2.)
                       - 4.*two*twoReduced*wCovTwoTwoReduced);
   if(v2PrimeErrorSquared>0.){v2PrimeError = pow(v2PrimeErrorSquared,0.5);}
   if(TMath::Abs(v2Prime)>0.)
   {
    fDiffFlow[t][pe][0]->SetBinContent(b,v2Prime); 
    fDiffFlow[t][pe][0]->SetBinError(b,v2PrimeError);    
   }  
  } // end of if(two>0.) 
  // differential flow:
  // v'{4}
  if(2.*pow(two,2.)-four > 0.) 
  {
   v4Prime = (2.*two*twoReduced-fourReduced)/pow(2.*pow(two,2.)-four,3./4.);
   v4PrimeErrorSquared = pow(2.*pow(two,2.)-four,-7./2.)
                       * (pow(2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced,2.)*pow(twoError,2.)
                       + (9./16.)*pow(2.*two*twoReduced-fourReduced,2.)*pow(fourError,2.)
                       + 4.*pow(two,2.)*pow(2.*pow(two,2.)-four,2.)*pow(twoReducedError,2.)
                       + pow(2.*pow(two,2.)-four,2.)*pow(fourReducedError,2.)                          
                       - (3./2.)*(2.*two*twoReduced-fourReduced)
                       * (2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced)*wCovTwoFour
                       - 4.*two*(2.*pow(two,2.)-four)
                       * (2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced)*wCovTwoTwoReduced
                       + 2.*(2.*pow(two,2.)-four)
                       * (2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced)*wCovTwoFourReduced
                       + 3.*two*(2.*pow(two,2.)-four)*(2.*two*twoReduced-fourReduced)*wCovFourTwoReduced
                       - (3./2.)*(2.*pow(two,2.)-four)*(2.*two*twoReduced-fourReduced)*wCovFourFourReduced 
                       - 4.*two*pow(2.*pow(two,2.)-four,2.)*wCovTwoReducedFourReduced);  
   if(v4PrimeErrorSquared>0.){v4PrimeError = pow(v4PrimeErrorSquared,0.5);}        
   if(TMath::Abs(v4Prime)>0.)
   {
    fDiffFlow[t][pe][1]->SetBinContent(b,v4Prime);
    fDiffFlow[t][pe][1]->SetBinError(b,v4PrimeError);     
   }
  } // end of if(2.*pow(two,2.)-four > 0.)  
 } // end of for(Int_t b=1;b<=fnBinsPtEta[pe];b++)

} // end of AliFlowAnalysisWithQCumulants::CalculateDiffFlow(TString type, Bool_t useParticleWeights)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::Calculate2DDiffFlow(TString type)
{
 // Calculate final results for 2D diferential flow.

 // to be improved - check pointers used in this method

 Int_t t = 0; // RP or POI

 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   } 
 
 // Differential flow:
 Double_t v2Prime = 0.; // v'{2}                   
 Double_t v4Prime = 0.; // v'{4}
 // Differential cumulants:
 Double_t qc2Prime = 0.; // QC{2'}                   
 Double_t qc4Prime = 0.; // QC{4'}
 // Looping over all (pt,eta) bins and calculating differential flow: 
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   // QC{2'}:
   qc2Prime = f2DDiffFlowCumulants[t][0]->GetBinContent(f2DDiffFlowCumulants[t][0]->GetBin(p,e));
   if(qc2Prime>=0.)
   {
    v2Prime = pow(qc2Prime,0.5);
    f2DDiffFlow[t][0]->SetBinContent(f2DDiffFlow[t][0]->GetBin(p,e),v2Prime); 
   } 
   // QC{4'}:
   qc4Prime = f2DDiffFlowCumulants[t][1]->GetBinContent(f2DDiffFlowCumulants[t][1]->GetBin(p,e));
   if(qc4Prime<=0.)
   {
    v4Prime = pow(-1.*qc4Prime,1./4.);
    f2DDiffFlow[t][1]->SetBinContent(f2DDiffFlow[t][1]->GetBin(p,e),v4Prime); 
   }   
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
} // end of void AliFlowAnalysisWithQCumulants::Calculate2DDiffFlow(TString type)  

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::StoreIntFlowFlags()
{
 // a) Store all flags for integrated flow in profile fIntFlowFlags.
 
 if(!fIntFlowFlags)
 {
  cout<<"WARNING: fIntFlowFlags is NULL in AFAWQC::SFFIF() !!!!"<<endl;
  exit(0);
 } 

 // particle weights used or not:
 fIntFlowFlags->Fill(0.5,(Int_t)fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights);
 // which event weights were used:
 if(strcmp(fMultiplicityWeight->Data(),"combinations"))
 {
  fIntFlowFlags->Fill(1.5,0); // 0 = "combinations" (default)
 } else if(strcmp(fMultiplicityWeight->Data(),"unit"))
   {
    fIntFlowFlags->Fill(1.5,1); // 1 = "unit"   
   } else if(strcmp(fMultiplicityWeight->Data(),"multiplicity"))
     {
      fIntFlowFlags->Fill(1.5,2); // 2 = "multiplicity"        
     } 
 fIntFlowFlags->Fill(2.5,(Int_t)fApplyCorrectionForNUA);
 fIntFlowFlags->Fill(3.5,(Int_t)fPrintFinalResults[0]);
 fIntFlowFlags->Fill(4.5,(Int_t)fPrintFinalResults[1]);
 fIntFlowFlags->Fill(5.5,(Int_t)fPrintFinalResults[2]);
 fIntFlowFlags->Fill(6.5,(Int_t)fPrintFinalResults[3]);
 fIntFlowFlags->Fill(7.5,(Int_t)fApplyCorrectionForNUAVsM);
 fIntFlowFlags->Fill(8.5,(Int_t)fPropagateErrorAlsoFromNIT);
 fIntFlowFlags->Fill(9.5,(Int_t)fCalculateCumulantsVsM);
 fIntFlowFlags->Fill(10.5,(Int_t)fMinimumBiasReferenceFlow);
 fIntFlowFlags->Fill(11.5,(Int_t)fForgetAboutCovariances);
 fIntFlowFlags->Fill(12.5,(Int_t)fStorePhiDistributionForOneEvent); 
 fIntFlowFlags->Fill(13.5,(Int_t)fFillMultipleControlHistograms);  
 fIntFlowFlags->Fill(14.5,(Int_t)fCalculateAllCorrelationsVsM);  
} // end of void AliFlowAnalysisWithQCumulants::StoreIntFlowFlags()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::StoreDiffFlowFlags()
{
 // Store all flags for differential flow in the profile fDiffFlowFlags.
  
 if(!fDiffFlowFlags)
 {
  printf("\n WARNING (QC): fDiffFlowFlags is NULL in AFAWQC::SDFF() !!!!\n\n");
  exit(0);
 } 
 
 fDiffFlowFlags->Fill(0.5,fCalculateDiffFlow); // calculate differential flow
 fDiffFlowFlags->Fill(1.5,fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights); // particle weights used or not?
 //fDiffFlowFlags->Fill(2.5,""); // which event weight was used? ("combinations", "unit" or "multiplicity") to be improved - finalized
 fDiffFlowFlags->Fill(3.5,fApplyCorrectionForNUA); // corrected for non-uniform acceptance or not
 fDiffFlowFlags->Fill(4.5,fCalculate2DDiffFlow); // calculate also 2D differential flow vs (pt,eta) 
 fDiffFlowFlags->Fill(5.5,fCalculateDiffFlowVsEta); // if you set kFALSE only differential flow vs pt is calculated
     
} // end of void AliFlowAnalysisWithQCumulants::StoreDiffFlowFlags()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::GetPointersForCommonHistograms() 
{
 // Access all pointers to common control and common result histograms and profiles.
 
 TString sCommonConstantsName = "fCommonConstants";
 sCommonConstantsName += fAnalysisLabel->Data();
 fCommonConstants = dynamic_cast<TProfile*>(fHistList->FindObject(sCommonConstantsName.Data()));
 if(!fCommonConstants)
 {
  printf("\n WARNING (QC): fCommonConstants is NULL in AFAWQC::GPFCH() !!!!\n\n");
  exit(0);
 }
 
 // to be improved - lines bellow can be implemented better.
 
 TString commonHistsName = "AliFlowCommonHistQC";
 commonHistsName += fAnalysisLabel->Data();
 AliFlowCommonHist *commonHist = dynamic_cast<AliFlowCommonHist*>(fHistList->FindObject(commonHistsName.Data()));
 if(commonHist) 
 {
  this->SetCommonHists(commonHist); 
  if(fCommonHists->GetHarmonic())
  {
   fHarmonic = (Int_t)(fCommonHists->GetHarmonic())->GetBinContent(1);
  } 
 } // end of if(commonHist) 
 TString commonHists2ndOrderName = "AliFlowCommonHist2ndOrderQC";
 commonHists2ndOrderName += fAnalysisLabel->Data();
 AliFlowCommonHist *commonHist2nd = dynamic_cast<AliFlowCommonHist*>(fHistList->FindObject(commonHists2ndOrderName.Data()));
 if(commonHist2nd) this->SetCommonHists2nd(commonHist2nd);   
 TString commonHists4thOrderName = "AliFlowCommonHist4thOrderQC";
 commonHists4thOrderName += fAnalysisLabel->Data();
 AliFlowCommonHist *commonHist4th = dynamic_cast<AliFlowCommonHist*>(fHistList->FindObject(commonHists4thOrderName.Data()));
 if(commonHist4th) this->SetCommonHists4th(commonHist4th);  
 TString commonHists6thOrderName = "AliFlowCommonHist6thOrderQC";
 commonHists6thOrderName += fAnalysisLabel->Data();
 AliFlowCommonHist *commonHist6th = dynamic_cast<AliFlowCommonHist*>(fHistList->FindObject(commonHists6thOrderName.Data()));
 if(commonHist6th) this->SetCommonHists6th(commonHist6th);  
 TString commonHists8thOrderName = "AliFlowCommonHist8thOrderQC";
 commonHists8thOrderName += fAnalysisLabel->Data();
 AliFlowCommonHist *commonHist8th = dynamic_cast<AliFlowCommonHist*>(fHistList->FindObject(commonHists8thOrderName.Data()));
 if(commonHist8th) this->SetCommonHists8th(commonHist8th); 
  
 TString commonHistResults2ndOrderName = "AliFlowCommonHistResults2ndOrderQC"; 
 commonHistResults2ndOrderName += fAnalysisLabel->Data(); 
 AliFlowCommonHistResults *commonHistRes2nd = dynamic_cast<AliFlowCommonHistResults*>
                                              (fHistList->FindObject(commonHistResults2ndOrderName.Data()));
 if(commonHistRes2nd) this->SetCommonHistsResults2nd(commonHistRes2nd);   
 TString commonHistResults4thOrderName = "AliFlowCommonHistResults4thOrderQC";
 commonHistResults4thOrderName += fAnalysisLabel->Data();
 AliFlowCommonHistResults *commonHistRes4th = dynamic_cast<AliFlowCommonHistResults*>
                                              (fHistList->FindObject(commonHistResults4thOrderName.Data()));
 if(commonHistRes4th) this->SetCommonHistsResults4th(commonHistRes4th);  
 TString commonHistResults6thOrderName = "AliFlowCommonHistResults6thOrderQC";
 commonHistResults6thOrderName += fAnalysisLabel->Data();
 AliFlowCommonHistResults *commonHistRes6th = dynamic_cast<AliFlowCommonHistResults*>
                                              (fHistList->FindObject(commonHistResults6thOrderName.Data()));
 if(commonHistRes6th) this->SetCommonHistsResults6th(commonHistRes6th);  
 TString commonHistResults8thOrderName = "AliFlowCommonHistResults8thOrderQC";
 commonHistResults8thOrderName += fAnalysisLabel->Data();
 AliFlowCommonHistResults *commonHistRes8th = dynamic_cast<AliFlowCommonHistResults*>
                                              (fHistList->FindObject(commonHistResults8thOrderName.Data()));  
 if(commonHistRes8th) this->SetCommonHistsResults8th(commonHistRes8th);
       
} // end of void AliFlowAnalysisWithQCumulants::GetPointersForCommonHistograms() 

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::GetPointersForParticleWeightsHistograms() 
{
 // Get pointers for histograms with particle weights.

 TList *weightsList = dynamic_cast<TList*>(fHistList->FindObject("Weights"));
 if(!weightsList){printf("\n WARNING (QC): weightsList is NULL in AFAWQC::GPFPWH() !!!!\n");exit(0);}
 this->SetWeightsList(weightsList);
 TString fUseParticleWeightsName = "fUseParticleWeightsQC"; // to be improved (hirdwired label QC)
 fUseParticleWeightsName += fAnalysisLabel->Data();
 TProfile *useParticleWeights = dynamic_cast<TProfile*>(weightsList->FindObject(fUseParticleWeightsName.Data()));
 if(useParticleWeights)
 {
  this->SetUseParticleWeights(useParticleWeights);  
  fUsePhiWeights = (Int_t)fUseParticleWeights->GetBinContent(1); 
  fUsePtWeights = (Int_t)fUseParticleWeights->GetBinContent(2); 
  fUseEtaWeights = (Int_t)fUseParticleWeights->GetBinContent(3);  
  fUseTrackWeights = (Int_t)fUseParticleWeights->GetBinContent(4);  
 }
} // end of void AliFlowAnalysisWithQCumulants::GetPointersForParticleWeightsHistograms(); 

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::GetPointersForIntFlowHistograms() 
{
 // Get pointers for histograms and profiles relevant for integrated flow:
 //  a) Get pointer to base list for integrated flow holding profile fIntFlowFlags and lists fIntFlowProfiles and fIntFlowResults.
 //  b) Get pointer to profile fIntFlowFlags holding all flags for integrated flow.
 //  c) Get pointer to list fIntFlowProfiles and pointers to all objects that she holds. 
 //  d) Get pointer to list fIntFlowResults and pointers to all objects that she holds. 
  
 TString sinCosFlag[2] = {"sin","cos"}; // to be improved (should I promote this to data member?)
 TString powerFlag[2] = {"linear","quadratic"}; // to be improved (should I promote this to data member?)
 TString correlationFlag[4] = {"#LT#LT2#GT#GT","#LT#LT4#GT#GT","#LT#LT6#GT#GT","#LT#LT8#GT#GT"}; // to be improved (should I promote this to data member?)
 TString squaredCorrelationFlag[4] = {"#LT#LT2#GT^{2}#GT","#LT#LT4#GT^{2}#GT","#LT#LT6#GT^{2}#GT","#LT#LT8#GT^{2}#GT"}; // to be improved (should I promote this to data member?)
 
 // a) Get pointer to base list for integrated flow holding profile fIntFlowFlags and lists fIntFlowProfiles and fIntFlowResults:
 TList *intFlowList = NULL;
 intFlowList = dynamic_cast<TList*>(fHistList->FindObject("Integrated Flow"));
 if(!intFlowList) 
 {
  cout<<"WARNING: intFlowList is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
  exit(0); 
 }  
  
 // b) Get pointer to profile fIntFlowFlags holding all flags for integrated flow:
 TString intFlowFlagsName = "fIntFlowFlags";
 intFlowFlagsName += fAnalysisLabel->Data();
 TProfile *intFlowFlags = dynamic_cast<TProfile*>(intFlowList->FindObject(intFlowFlagsName.Data()));
 if(intFlowFlags)
 {
  this->SetIntFlowFlags(intFlowFlags);  
  fApplyCorrectionForNUA = (Bool_t)intFlowFlags->GetBinContent(3); 
  fApplyCorrectionForNUAVsM = (Bool_t)intFlowFlags->GetBinContent(8); 
  fCalculateCumulantsVsM = (Bool_t)intFlowFlags->GetBinContent(10);  
 } else 
   {
    cout<<"WARNING: intFlowFlags is NULL in FAWQC::GPFIFH() !!!!"<<endl;
   }
  
  // c) Get pointer to list fIntFlowProfiles and pointers to all objects that she holds:
  TList *intFlowProfiles = NULL;
  intFlowProfiles = dynamic_cast<TList*>(intFlowList->FindObject("Profiles"));
  if(intFlowProfiles)  
  {
   // average multiplicities:
   TString avMultiplicityName = "fAvMultiplicity";
   avMultiplicityName += fAnalysisLabel->Data();
   TProfile *avMultiplicity = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(avMultiplicityName.Data()));
   if(avMultiplicity) 
   {
    this->SetAvMultiplicity(avMultiplicity);
   } else 
     {
      cout<<"WARNING: avMultiplicity is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     }
   // average correlations <<2>>, <<4>>, <<6>> and <<8>> (with wrong errors!):
   TString intFlowCorrelationsProName = "fIntFlowCorrelationsPro";
   intFlowCorrelationsProName += fAnalysisLabel->Data();
   TProfile *intFlowCorrelationsPro = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(intFlowCorrelationsProName.Data()));
   if(intFlowCorrelationsPro) 
   {
    this->SetIntFlowCorrelationsPro(intFlowCorrelationsPro);
   } else 
     {
      cout<<"WARNING: intFlowCorrelationsPro is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     }      
   // average squared correlations <<2>^2>, <<4>^2>, <<6>^2> and <<8^2>>:
   TString intFlowSquaredCorrelationsProName = "fIntFlowSquaredCorrelationsPro";
   intFlowSquaredCorrelationsProName += fAnalysisLabel->Data();
   TProfile *intFlowSquaredCorrelationsPro = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(intFlowSquaredCorrelationsProName.Data()));
   if(intFlowSquaredCorrelationsPro) 
   {
    this->SetIntFlowSquaredCorrelationsPro(intFlowSquaredCorrelationsPro);
   } else 
     {
      cout<<"WARNING: intFlowSquaredCorrelationsPro is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     }             
   if(fCalculateCumulantsVsM)
   {
    // Average correlations <<2>>, <<4>>, <<6>> and <<8>> versus multiplicity for all events (error is wrong here):   
    TString intFlowCorrelationsVsMProName = "fIntFlowCorrelationsVsMPro";
    intFlowCorrelationsVsMProName += fAnalysisLabel->Data();
    for(Int_t ci=0;ci<4;ci++) // correlation index
    {
     TProfile *intFlowCorrelationsVsMPro = dynamic_cast<TProfile*>
                                        (intFlowProfiles->FindObject(Form("%s, %s",intFlowCorrelationsVsMProName.Data(),correlationFlag[ci].Data())));
     if(intFlowCorrelationsVsMPro)
     {
      this->SetIntFlowCorrelationsVsMPro(intFlowCorrelationsVsMPro,ci);
     } else
       {
        cout<<"WARNING: "<<Form("intFlowCorrelationsVsMPro[%d]",ci)<<" is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
       }   
    } // end of for(Int_t ci=0;ci<4;ci++) // correlation index 
    // Average squared correlations <<2>^2>, <<4>^2>, <<6>^2> and <<8>^2> versus multiplicity for all events:   
    TString intFlowSquaredCorrelationsVsMProName = "fIntFlowSquaredCorrelationsVsMPro";
    intFlowSquaredCorrelationsVsMProName += fAnalysisLabel->Data();
    for(Int_t ci=0;ci<4;ci++) // correlation index
    {
     TProfile *intFlowSquaredCorrelationsVsMPro = dynamic_cast<TProfile*>
                      (intFlowProfiles->FindObject(Form("%s, %s",intFlowSquaredCorrelationsVsMProName.Data(),squaredCorrelationFlag[ci].Data())));
     if(intFlowSquaredCorrelationsVsMPro)
     {
      this->SetIntFlowSquaredCorrelationsVsMPro(intFlowSquaredCorrelationsVsMPro,ci);
     } else
       {
        cout<<"WARNING: "<<Form("intFlowSquaredCorrelationsVsMPro[%d]",ci)<<" is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
       }   
    } // end of for(Int_t ci=0;ci<4;ci++) // correlation index 
   } // end of if(fCalculateCumulantsVsM)
   // average all correlations for integrated flow (with wrong errors!):
   TString intFlowCorrelationsAllProName = "fIntFlowCorrelationsAllPro";
   intFlowCorrelationsAllProName += fAnalysisLabel->Data();
   TProfile *intFlowCorrelationsAllPro = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(intFlowCorrelationsAllProName.Data()));
   if(intFlowCorrelationsAllPro) 
   {
    this->SetIntFlowCorrelationsAllPro(intFlowCorrelationsAllPro);
   } else 
     {
      cout<<"WARNING: intFlowCorrelationsAllPro is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     }     
   // average extra correlations for integrated flow (which appear only when particle weights are used):
   // (to be improved: Weak point in implementation, I am assuming here that method GetPointersForParticleWeightsHistograms() was called)
   if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)
   {
    TString intFlowExtraCorrelationsProName = "fIntFlowExtraCorrelationsPro";
    intFlowExtraCorrelationsProName += fAnalysisLabel->Data();
    TProfile *intFlowExtraCorrelationsPro = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(intFlowExtraCorrelationsProName.Data()));
    if(intFlowExtraCorrelationsPro) 
    {
     this->SetIntFlowExtraCorrelationsPro(intFlowExtraCorrelationsPro);
    } else 
      {
       cout<<"WARNING: intFlowExtraCorrelationsPro is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
      }
   } // end of if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)        
   // average products of correlations <2>, <4>, <6> and <8>:  
   TString intFlowProductOfCorrelationsProName = "fIntFlowProductOfCorrelationsPro";
   intFlowProductOfCorrelationsProName += fAnalysisLabel->Data();
   TProfile *intFlowProductOfCorrelationsPro = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(intFlowProductOfCorrelationsProName.Data()));
   if(intFlowProductOfCorrelationsPro) 
   {
    this->SetIntFlowProductOfCorrelationsPro(intFlowProductOfCorrelationsPro);
   } else 
     {
      cout<<"WARNING: intFlowProductOfCorrelationsPro is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     }               
   // average product of correlations <2>, <4>, <6> and <8> versus multiplicity  
   // [0=<<2><4>>,1=<<2><6>>,2=<<2><8>>,3=<<4><6>>,4=<<4><8>>,5=<<6><8>>]  
   if(fCalculateCumulantsVsM)
   {
    TString intFlowProductOfCorrelationsVsMProName = "fIntFlowProductOfCorrelationsVsMPro";
    intFlowProductOfCorrelationsVsMProName += fAnalysisLabel->Data();
    TString productFlag[6] = {"#LT#LT2#GT#LT4#GT#GT","#LT#LT2#GT#LT6#GT#GT","#LT#LT2#GT#LT8#GT#GT",
                              "#LT#LT4#GT#LT6#GT#GT","#LT#LT4#GT#LT8#GT#GT","#LT#LT6#GT#LT8#GT#GT"};
    for(Int_t pi=0;pi<6;pi++)
    { 
     TProfile *intFlowProductOfCorrelationsVsMPro = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(Form("%s, %s",intFlowProductOfCorrelationsVsMProName.Data(),productFlag[pi].Data())));
     if(intFlowProductOfCorrelationsVsMPro)
     {
      this->SetIntFlowProductOfCorrelationsVsMPro(intFlowProductOfCorrelationsVsMPro,pi);
     } else
       {
        cout<<"WARNING: "<<Form("intFlowProductOfCorrelationsVsMPro[%d]",pi)<<" is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
       }
    } // end of for(Int_t pi=0;pi<6;pi++)
   } // end of if(fCalculateCumulantsVsM)
   // average correction terms for non-uniform acceptance (with wrong errors!):
   for(Int_t sc=0;sc<2;sc++)
   {
    TString intFlowCorrectionTermsForNUAProName = "fIntFlowCorrectionTermsForNUAPro";
    intFlowCorrectionTermsForNUAProName += fAnalysisLabel->Data();
    TProfile *intFlowCorrectionTermsForNUAPro = dynamic_cast<TProfile*>(intFlowProfiles->FindObject((Form("%s: %s terms",intFlowCorrectionTermsForNUAProName.Data(),sinCosFlag[sc].Data()))));
    if(intFlowCorrectionTermsForNUAPro) 
    {
     this->SetIntFlowCorrectionTermsForNUAPro(intFlowCorrectionTermsForNUAPro,sc);
    } else 
      {
       cout<<"WARNING: intFlowCorrectionTermsForNUAPro is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
       cout<<"sc = "<<sc<<endl;
      } 
    // versus multiplicity:
    if(fCalculateCumulantsVsM)
    {
     TString correctionTermFlag[4] = {"(n(phi1))","(n(phi1+phi2))","(n(phi1-phi2-phi3))","(n(2phi1-phi2))"}; // to be improved - hardwired 4
     TString intFlowCorrectionTermsForNUAVsMProName = "fIntFlowCorrectionTermsForNUAVsMPro";
     intFlowCorrectionTermsForNUAVsMProName += fAnalysisLabel->Data();
     for(Int_t ci=0;ci<4;ci++) // correction term index (to be improved - hardwired 4)
     {
      TProfile *intFlowCorrectionTermsForNUAVsMPro = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(Form("%s: #LT#LT%s%s#GT#GT",intFlowCorrectionTermsForNUAVsMProName.Data(),sinCosFlag[sc].Data(),correctionTermFlag[ci].Data())));
      if(intFlowCorrectionTermsForNUAVsMPro) 
      {
       this->SetIntFlowCorrectionTermsForNUAVsMPro(intFlowCorrectionTermsForNUAVsMPro,sc,ci);
      } else 
        {
         cout<<"WARNING: intFlowCorrectionTermsForNUAVsMPro is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
         cout<<"sc = "<<sc<<endl;
         cout<<"ci = "<<ci<<endl;
        }       
     } // end of for(Int_t ci=0;ci<4;ci++) // correction term index (to be improved - hardwired 4)
    } // end of if(fCalculateCumulantsVsM)
   } // end of for(Int_t sc=0;sc<2;sc++)           
   // average products of correction terms for NUA:  
   TString intFlowProductOfCorrectionTermsForNUAProName = "fIntFlowProductOfCorrectionTermsForNUAPro";
   intFlowProductOfCorrectionTermsForNUAProName += fAnalysisLabel->Data();
   TProfile *intFlowProductOfCorrectionTermsForNUAPro = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(intFlowProductOfCorrectionTermsForNUAProName.Data()));
   if(intFlowProductOfCorrectionTermsForNUAPro) 
   {
    this->SetIntFlowProductOfCorrectionTermsForNUAPro(intFlowProductOfCorrectionTermsForNUAPro);
   } else 
     {
      cout<<"WARNING: intFlowProductOfCorrectionTermsForNUAPro is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     }     
  } else // to if(intFlowProfiles)  
    {
     cout<<"WARNING: intFlowProfiles is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
    }
   
  //  d) Get pointer to list fIntFlowResults and pointers to all objects that she holds. 
  TList *intFlowResults = NULL;
  intFlowResults = dynamic_cast<TList*>(intFlowList->FindObject("Results"));
  if(intFlowResults)
  {
   // average correlations <<2>>, <<4>>, <<6>> and <<8>> (with correct errors!):
   TString intFlowCorrelationsHistName = "fIntFlowCorrelationsHist";
   intFlowCorrelationsHistName += fAnalysisLabel->Data();
   TH1D *intFlowCorrelationsHist = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowCorrelationsHistName.Data()));
   if(intFlowCorrelationsHist) 
   {
    this->SetIntFlowCorrelationsHist(intFlowCorrelationsHist);
   } else 
     {
      cout<<"WARNING: intFlowCorrelationsHist is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     } 
   // average correlations <<2>>, <<4>>, <<6>> and <<8>> (with correct errors!) vs M:    
   if(fCalculateCumulantsVsM)
   {
    TString intFlowCorrelationsVsMHistName = "fIntFlowCorrelationsVsMHist";
    intFlowCorrelationsVsMHistName += fAnalysisLabel->Data();
    for(Int_t ci=0;ci<4;ci++) // correlation index
    {
     TH1D *intFlowCorrelationsVsMHist = dynamic_cast<TH1D*>
                                        (intFlowResults->FindObject(Form("%s, %s",intFlowCorrelationsVsMHistName.Data(),correlationFlag[ci].Data())));
     if(intFlowCorrelationsVsMHist)
     {
      this->SetIntFlowCorrelationsVsMHist(intFlowCorrelationsVsMHist,ci);
     } else
       {
        cout<<"WARNING: "<<Form("intFlowCorrelationsVsMHist[%d]",ci)<<" is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
       }   
    } // end of for(Int_t ci=0;ci<4;ci++) // correlation index   
   } // end of if(fCalculateCumulantsVsM)
   // average all correlations for integrated flow (with correct errors!):
   TString intFlowCorrelationsAllHistName = "fIntFlowCorrelationsAllHist";
   intFlowCorrelationsAllHistName += fAnalysisLabel->Data();
   TH1D *intFlowCorrelationsAllHist = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowCorrelationsAllHistName.Data()));
   if(intFlowCorrelationsAllHist) 
   {
    this->SetIntFlowCorrelationsAllHist(intFlowCorrelationsAllHist);
   } else 
     {
      cout<<"WARNING: intFlowCorrelationsAllHist is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     }  
   // average correction terms for non-uniform acceptance (with correct errors!):
   TString intFlowCorrectionTermsForNUAHistName = "fIntFlowCorrectionTermsForNUAHist";
   intFlowCorrectionTermsForNUAHistName += fAnalysisLabel->Data();
   for(Int_t sc=0;sc<2;sc++)
   {
    TH1D *intFlowCorrectionTermsForNUAHist = dynamic_cast<TH1D*>(intFlowResults->FindObject((Form("%s: %s terms",intFlowCorrectionTermsForNUAHistName.Data(),sinCosFlag[sc].Data()))));
    if(intFlowCorrectionTermsForNUAHist) 
    {
     this->SetIntFlowCorrectionTermsForNUAHist(intFlowCorrectionTermsForNUAHist,sc);
    } else 
      {
       cout<<"WARNING: intFlowCorrectionTermsForNUAHist is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
       cout<<"sc = "<<sc<<endl;
      } 
   } // end of for(Int_t sc=0;sc<2;sc++)           
   // covariances (multiplied with weight dependent prefactor):
   TString intFlowCovariancesName = "fIntFlowCovariances";
   intFlowCovariancesName += fAnalysisLabel->Data();
   TH1D *intFlowCovariances = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowCovariancesName.Data()));
   if(intFlowCovariances) 
   {
    this->SetIntFlowCovariances(intFlowCovariances); 
   } else 
     {
      cout<<"WARNING: intFlowCovariances is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     } 
   // sum of linear and quadratic event weights for <2>, <4>, <6> and <8>:
   TString intFlowSumOfEventWeightsName = "fIntFlowSumOfEventWeights";
   intFlowSumOfEventWeightsName += fAnalysisLabel->Data();
   for(Int_t power=0;power<2;power++)
   {
    TH1D *intFlowSumOfEventWeights = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("%s: %s",intFlowSumOfEventWeightsName.Data(),powerFlag[power].Data())));
    if(intFlowSumOfEventWeights) 
    {
     this->SetIntFlowSumOfEventWeights(intFlowSumOfEventWeights,power);
    } else 
      {
       cout<<"WARNING: intFlowSumOfEventWeights is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
       cout<<"power = "<<power<<endl;
      }                                   
   } // end of for(Int_t power=0;power<2;power++)                                                                  
   // sum of products of event weights for correlations <2>, <4>, <6> and <8>:  
   TString intFlowSumOfProductOfEventWeightsName = "fIntFlowSumOfProductOfEventWeights";
   intFlowSumOfProductOfEventWeightsName += fAnalysisLabel->Data();
   TH1D *intFlowSumOfProductOfEventWeights = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowSumOfProductOfEventWeightsName.Data()));
   if(intFlowSumOfProductOfEventWeights) 
   {
    this->SetIntFlowSumOfProductOfEventWeights(intFlowSumOfProductOfEventWeights);
   } else 
     {
      cout<<"WARNING: intFlowSumOfProductOfEventWeights is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     } 
   // final result for covariances of correlations (multiplied with weight dependent prefactor) versus M
   // [0=Cov(2,4),1=Cov(2,6),2=Cov(2,8),3=Cov(4,6),4=Cov(4,8),5=Cov(6,8)]:
   if(fCalculateCumulantsVsM)
   {
    TString intFlowCovariancesVsMName = "fIntFlowCovariancesVsM";
    intFlowCovariancesVsMName += fAnalysisLabel->Data();
    TString covarianceFlag[6] = {"Cov(<2>,<4>)","Cov(<2>,<6>)","Cov(<2>,<8>)","Cov(<4>,<6>)","Cov(<4>,<8>)","Cov(<6>,<8>)"};
    for(Int_t ci=0;ci<6;ci++)
    { 
     TH1D *intFlowCovariancesVsM = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("%s, %s",intFlowCovariancesVsMName.Data(),covarianceFlag[ci].Data())));
     if(intFlowCovariancesVsM)
     {
      this->SetIntFlowCovariancesVsM(intFlowCovariancesVsM,ci);
     } else
       {
        cout<<"WARNING: "<<Form("intFlowCovariancesVsM[%d]",ci)<<" is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
       }    
    } // end of for(Int_t ci=0;ci<6;ci++)
   } // end of if(fCalculateCumulantsVsM)
   // sum of linear and quadratic event weights for <2>, <4>, <6> and <8> versus multiplicity
   // [0=sum{w_{<2>}},1=sum{w_{<4>}},2=sum{w_{<6>}},3=sum{w_{<8>}}][0=linear 1,1=quadratic]:
   if(fCalculateCumulantsVsM)
   {
    TString intFlowSumOfEventWeightsVsMName = "fIntFlowSumOfEventWeightsVsM";
    intFlowSumOfEventWeightsVsMName += fAnalysisLabel->Data();
    TString sumFlag[2][4] = {{"#sum_{i=1}^{N} w_{<2>}","#sum_{i=1}^{N} w_{<4>}","#sum_{i=1}^{N} w_{<6>}","#sum_{i=1}^{N} w_{<8>}"},
                             {"#sum_{i=1}^{N} w_{<2>}^{2}","#sum_{i=1}^{N} w_{<4>}^{2}","#sum_{i=1}^{N} w_{<6>}^{2}","#sum_{i=1}^{N} w_{<8>}^{2}"}};
    for(Int_t si=0;si<4;si++)
    {
     for(Int_t power=0;power<2;power++)
     {
      TH1D *intFlowSumOfEventWeightsVsM = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("%s, %s",intFlowSumOfEventWeightsVsMName.Data(),sumFlag[power][si].Data())));
      if(intFlowSumOfEventWeightsVsM)
      {
       this->SetIntFlowSumOfEventWeightsVsM(intFlowSumOfEventWeightsVsM,si,power);
      } else
        {
         cout<<"WARNING: "<<Form("intFlowSumOfEventWeightsVsM[%d][%d]",si,power)<<" is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
        }    
     } // end of for(Int_t power=0;power<2;power++)
    } // end of for(Int_t si=0;si<4;si++)   
   } // end of if(fCalculateCumulantsVsM)
   // sum of products of event weights for correlations <2>, <4>, <6> and <8> vs M
   // [0=sum{w_{<2>}w_{<4>}},1=sum{w_{<2>}w_{<6>}},2=sum{w_{<2>}w_{<8>}},
   //  3=sum{w_{<4>}w_{<6>}},4=sum{w_{<4>}w_{<8>}},5=sum{w_{<6>}w_{<8>}}]:  
   if(fCalculateCumulantsVsM)
   {
    TString intFlowSumOfProductOfEventWeightsVsMName = "fIntFlowSumOfProductOfEventWeightsVsM";
    intFlowSumOfProductOfEventWeightsVsMName += fAnalysisLabel->Data();
    TString sopowFlag[6] = {"#sum_{i=1}^{N} w_{<2>} w_{<4>}","#sum_{i=1}^{N} w_{<2>} w_{<6>}","#sum_{i=1}^{N} w_{<2>} w_{<8>}",
                            "#sum_{i=1}^{N} w_{<4>} w_{<6>}","#sum_{i=1}^{N} w_{<4>} w_{<8>}","#sum_{i=1}^{N} w_{<6>} w_{<8>}"}; 
    for(Int_t pi=0;pi<6;pi++)
    {
     TH1D *intFlowSumOfProductOfEventWeightsVsM = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("%s, %s",intFlowSumOfProductOfEventWeightsVsMName.Data(),sopowFlag[pi].Data())));
     if(intFlowSumOfProductOfEventWeightsVsM)
     {
      this->SetIntFlowSumOfProductOfEventWeightsVsM(intFlowSumOfProductOfEventWeightsVsM,pi);
     } else
       {
        cout<<"WARNING: "<<Form("intFlowSumOfProductOfEventWeightsVsM[%d]",pi)<<" is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
       }
    } // end of for(Int_t pi=0;pi<6;pi++)        
   } // end of if(fCalculateCumulantsVsM)
   // covariances for NUA (multiplied with weight dependent prefactor):
   TString intFlowCovariancesNUAName = "fIntFlowCovariancesNUA";
   intFlowCovariancesNUAName += fAnalysisLabel->Data();
   TH1D *intFlowCovariancesNUA = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowCovariancesNUAName.Data()));
   if(intFlowCovariancesNUA) 
   {
    this->SetIntFlowCovariancesNUA(intFlowCovariancesNUA); 
   } else 
     {
      cout<<"WARNING: intFlowCovariancesNUA is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     } 
   // sum of linear and quadratic event weights NUA terms:
   TString intFlowSumOfEventWeightsNUAName = "fIntFlowSumOfEventWeightsNUA";
   intFlowSumOfEventWeightsNUAName += fAnalysisLabel->Data();
   for(Int_t sc=0;sc<2;sc++)
   {
    for(Int_t power=0;power<2;power++)
    {
     TH1D *intFlowSumOfEventWeightsNUA = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("%s: %s, %s",intFlowSumOfEventWeightsNUAName.Data(),powerFlag[power].Data(),sinCosFlag[sc].Data())));
     if(intFlowSumOfEventWeightsNUA) 
     {
      this->SetIntFlowSumOfEventWeightsNUA(intFlowSumOfEventWeightsNUA,sc,power);
     } else 
       {
        cout<<"WARNING: intFlowSumOfEventWeightsNUA is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
        cout<<"sc    = "<<sc<<endl;
        cout<<"power = "<<power<<endl;
       }                                   
    } // end of for(Int_t power=0;power<2;power++)                                                                  
   } // end of for(Int_t sc=0;sc<2;sc++)     
   // sum of products of event weights for NUA terms:  
   TString intFlowSumOfProductOfEventWeightsNUAName = "fIntFlowSumOfProductOfEventWeightsNUA";
   intFlowSumOfProductOfEventWeightsNUAName += fAnalysisLabel->Data();
   TH1D *intFlowSumOfProductOfEventWeightsNUA = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowSumOfProductOfEventWeightsNUAName.Data()));
   if(intFlowSumOfProductOfEventWeightsNUA) 
   {
    this->SetIntFlowSumOfProductOfEventWeightsNUA(intFlowSumOfProductOfEventWeightsNUA);
   } else 
     {
      cout<<"WARNING: intFlowSumOfProductOfEventWeightsNUA is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     } 
   // Final results for reference Q-cumulants:
   TString intFlowQcumulantsName = "fIntFlowQcumulants";
   intFlowQcumulantsName += fAnalysisLabel->Data();
   TH1D *intFlowQcumulants = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowQcumulantsName.Data()));
   if(intFlowQcumulants) 
   {
    this->SetIntFlowQcumulants(intFlowQcumulants);
   } else 
     {
      cout<<"WARNING: intFlowQcumulants is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     }  
   // Final results for reference Q-cumulants rebinned in M:
   if(fCalculateCumulantsVsM)
   {
    TString intFlowQcumulantsRebinnedInMName = "fIntFlowQcumulantsRebinnedInM";
    intFlowQcumulantsRebinnedInMName += fAnalysisLabel->Data();
    TH1D *intFlowQcumulantsRebinnedInM = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowQcumulantsRebinnedInMName.Data()));
    if(intFlowQcumulantsRebinnedInM) 
    {
     this->SetIntFlowQcumulantsRebinnedInM(intFlowQcumulantsRebinnedInM);
    } else 
      {
       cout<<"WARNING: intFlowQcumulantsRebinnedInM is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
      }  
   } // end of if(fCalculateCumulantsVsM)
   // Ratio between error squared: with/without non-isotropic terms:
   TString intFlowQcumulantsErrorSquaredRatioName = "fIntFlowQcumulantsErrorSquaredRatio";
   intFlowQcumulantsErrorSquaredRatioName += fAnalysisLabel->Data();
   TH1D *intFlowQcumulantsErrorSquaredRatio = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowQcumulantsErrorSquaredRatioName.Data()));
   if(intFlowQcumulantsErrorSquaredRatio) 
   {
    this->SetIntFlowQcumulantsErrorSquaredRatio(intFlowQcumulantsErrorSquaredRatio);
   } else 
     {
      cout<<" WARNING: intntFlowQcumulantsErrorSquaredRatio is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
     }  
   // final results for integrated Q-cumulants versus multiplicity:
   TString cumulantFlag[4] = {"QC{2}","QC{4}","QC{6}","QC{8}"};
   if(fCalculateCumulantsVsM)
   {
    TString intFlowQcumulantsVsMName = "fIntFlowQcumulantsVsM";
    intFlowQcumulantsVsMName += fAnalysisLabel->Data();
    for(Int_t co=0;co<4;co++) // cumulant order
    {
     TH1D *intFlowQcumulantsVsM = dynamic_cast<TH1D*>
                                  (intFlowResults->FindObject(Form("%s, %s",intFlowQcumulantsVsMName.Data(),cumulantFlag[co].Data())));
     if(intFlowQcumulantsVsM)
     {
      this->SetIntFlowQcumulantsVsM(intFlowQcumulantsVsM,co);
     } else
       {
        cout<<"WARNING: "<<Form("intFlowQcumulantsVsM[%d]",co)<<" is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
       }
    } // end of for(Int_t co=0;co<4;co++) // cumulant order
   } // end of if(fCalculateCumulantsVsM)
   // Final reference flow estimates from Q-cumulants:
   TString intFlowName = "fIntFlow";
   intFlowName += fAnalysisLabel->Data();
   TH1D *intFlow = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowName.Data()));
   if(intFlow) 
   {
    this->SetIntFlow(intFlow);
   } else 
     {
      cout<<"WARNING: intFlow is NULL in AFAWQC::GPFIFH() !!!!"<<endl; 
     } 
   // Final reference flow estimates from Q-cumulants vs M rebinned in M:
   if(fCalculateCumulantsVsM)
   {
    TString intFlowRebinnedInMName = "fIntFlowRebinnedInM";
    intFlowRebinnedInMName += fAnalysisLabel->Data();
    TH1D *intFlowRebinnedInM = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowRebinnedInMName.Data()));
    if(intFlowRebinnedInM) 
    {
     this->SetIntFlowRebinnedInM(intFlowRebinnedInM);
    } else 
      {
       cout<<"WARNING: intFlowRebinnedInM is NULL in AFAWQC::GPFIFH() !!!!"<<endl; 
      } 
   } // end of if(fCalculateCumulantsVsM)
   // integrated flow from Q-cumulants versus multiplicity:
   if(fCalculateCumulantsVsM)
   {
    TString intFlowVsMName = "fIntFlowVsM";
    intFlowVsMName += fAnalysisLabel->Data();
    TString flowFlag[4] = {Form("v_{%d}{2,QC}",fHarmonic),Form("v_{%d}{4,QC}",fHarmonic),Form("v_{%d}{6,QC}",fHarmonic),Form("v_{%d}{8,QC}",fHarmonic)};
    for(Int_t co=0;co<4;co++) // cumulant order
    {
     TH1D *intFlowVsM = dynamic_cast<TH1D*>
                        (intFlowResults->FindObject(Form("%s, %s",intFlowVsMName.Data(),flowFlag[co].Data())));            
     if(intFlowVsM)
     {
      this->SetIntFlowVsM(intFlowVsM,co);
     } else
       {
        cout<<"WARNING: "<<Form("intFlowVsM[%d]",co)<<" is NULL in AFAWQC::GPFIFH() !!!!"<<endl;      
       }
    } // end of for(Int_t co=0;co<4;co++) // cumulant order
   } // end of if(fCalculateCumulantsVsM)
   // quantifying detector effects effects to correlations:
   TString intFlowDetectorBiasName = "fIntFlowDetectorBias";
   intFlowDetectorBiasName += fAnalysisLabel->Data();
   TH1D *intFlowDetectorBias = dynamic_cast<TH1D*>(intFlowResults->FindObject(intFlowDetectorBiasName.Data()));
   if(intFlowDetectorBias) 
   {
    this->SetIntFlowDetectorBias(intFlowDetectorBias);
   } else 
     {
      cout<<"WARNING: intFlowDetectorBias is NULL in AFAWQC::GPFIFH() !!!!"<<endl; 
     } 
   // quantifying detector effects effects to correlations vs multiplicity:
   if(fCalculateCumulantsVsM)
   {
    TString intFlowDetectorBiasVsMName = "fIntFlowDetectorBiasVsM";
    intFlowDetectorBiasVsMName += fAnalysisLabel->Data();
    for(Int_t ci=0;ci<4;ci++) // correlation index
    {
     TH1D *intFlowDetectorBiasVsM = dynamic_cast<TH1D*>
                                    (intFlowResults->FindObject(Form("%s for %s",intFlowDetectorBiasVsMName.Data(),cumulantFlag[ci].Data())));
     if(intFlowDetectorBiasVsM)
     {
      this->SetIntFlowDetectorBiasVsM(intFlowDetectorBiasVsM,ci);
     } else
       {
        cout<<"WARNING: "<<Form("intFlowDetectorBiasVsM[%d]",ci)<<" is NULL in AFAWQC::GPFIFH() !!!!"<<endl;      
       }
    } // end of for(Int_t ci=0;ci<4;ci++) // correlation index   
   } // end of if(fCalculateCumulantsVsM)
  } else // to if(intFlowResults)
    {
     cout<<"WARNING: intFlowResults is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
    }
    
} // end of void AliFlowAnalysisWithQCumulants::GetPointersForIntFlowHistograms()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::GetPointersFor2DDiffFlowHistograms()
{
 // Get pointers for 2D differential flow histograms.
 //  a) Check pointers used in this method;
 //  b) Get pointers to 2D differential flow lists;
 //  c) Get pointers to 2D differential flow profiles;
 //  d) Get pointers to 2D differential flow histograms. 

 // a) Check pointers used in this method:
 if(!fDiffFlowList)
 { 
  printf("\n WARNING (QC): fDiffFlowList is NULL in AFAWQC::GPF2DDFH() !!!!\n");
  printf("               Call method GetPointersForDiffFlowHistograms() first.\n\n");
  exit(0);
 }
 if(!fDiffFlowFlags)
 { 
  printf("\n WARNING (QC): fDiffFlowFlags is NULL in AFAWQC::GPF2DDFH() !!!!\n\n");
  printf("               Call method GetPointersForDiffFlowHistograms() first.\n\n");
  exit(0);
 }
 
 // b) Get pointers to 2D differential flow lists:
 this->SetCalculate2DDiffFlow((Bool_t)fDiffFlowFlags->GetBinContent(5)); // to be improved - hardwired 5
 if(!fCalculate2DDiffFlow){return;}
 TString typeFlag[2] = {"RP","POI"}; 
 TString reducedCorrelationIndex[4] = {"<2'>","<4'>","<6'>","<8'>"};
 TString differentialCumulantIndex[4] = {"QC{2'}","QC{4'}","QC{6'}","QC{8'}"};  
 TString differentialFlowIndex[4] = {"v'{2}","v'{4}","v'{6}","v'{8}"};  
 // Base list: 
 TString diffFlow2DListName = "2D"; 
 diffFlow2DListName += fAnalysisLabel->Data();
 fDiffFlow2D = dynamic_cast<TList*>(fDiffFlowList->FindObject(diffFlow2DListName.Data()));
 if(!fDiffFlow2D)
 { 
  printf("\n WARNING (QC): fDiffFlow2D is NULL in AFAWQC::GPFDFH() !!!!\n\n");
  exit(0);
 }
 // Lists holding profiles with 2D correlations: 
 TString s2DDiffFlowCorrelationsProListName = "Profiles with 2D correlations"; 
 s2DDiffFlowCorrelationsProListName += fAnalysisLabel->Data(); // to be improved
 for(Int_t t=0;t<2;t++)
 {
  f2DDiffFlowCorrelationsProList[t] = dynamic_cast<TList*>(fDiffFlow2D->FindObject(Form("Profiles with 2D correlations (%s)",typeFlag[t].Data())));
  if(!f2DDiffFlowCorrelationsProList[t])
  { 
   printf("\n WARNING (QC): f2DDiffFlowCorrelationsProList[%i] is NULL in AFAWQC::GPF2DFH() !!!!\n\n",t);
   exit(0);
  }
 } // end of for(Int_t t=0;t<2;t++) 
 
 // c) Get pointers to 2D differential flow profiles:
 TString s2DDiffFlowCorrelationsProName = "f2DDiffFlowCorrelationsPro";
 s2DDiffFlowCorrelationsProName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t rci=0;rci<4;rci++) // reduced correlation index
  {
   f2DDiffFlowCorrelationsPro[t][rci] = dynamic_cast<TProfile2D*>(f2DDiffFlowCorrelationsProList[t]->FindObject(Form("%s, %s, %s",s2DDiffFlowCorrelationsProName.Data(),typeFlag[t].Data(),reducedCorrelationIndex[rci].Data())));
   if(!f2DDiffFlowCorrelationsPro[t][rci])
   {
    printf("\n WARNING (QC): f2DDiffFlowCorrelationsPro[%i][%i] is NULL in AFAWQC::GPF2DFH() !!!!\n\n",t,rci);
    exit(0);   
   } else
     {
      this->Set2DDiffFlowCorrelationsPro(f2DDiffFlowCorrelationsPro[t][rci],t,rci);
     } 
  } // end of for(Int_t rci=0;rci<4;rci++) // reduced correlation index
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI 

 // d) Get pointers to 2D differential flow histograms: 
 TString s2DDiffFlowCumulantsName = "f2DDiffFlowCumulants";
 s2DDiffFlowCumulantsName += fAnalysisLabel->Data();
 TString s2DDiffFlowName = "f2DDiffFlow";
 s2DDiffFlowName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t rci=0;rci<4;rci++) // reduced correlation index
  {
   // 2D differential cumulants:
   f2DDiffFlowCumulants[t][rci] = dynamic_cast<TH2D*>(f2DDiffFlowCorrelationsProList[t]->FindObject(Form("%s, %s, %s",s2DDiffFlowCumulantsName.Data(),typeFlag[t].Data(),differentialCumulantIndex[rci].Data())));
   if(!f2DDiffFlowCumulants[t][rci])
   {
    printf("\n WARNING (QC): f2DDiffFlowCumulants[%i][%i] is NULL in AFAWQC::GPF2DFH() !!!!\n\n",t,rci);
    exit(0);   
   } else
     {
      this->Set2DDiffFlowCumulants(f2DDiffFlowCumulants[t][rci],t,rci);
     } 
   // 2D differential flow:
   f2DDiffFlow[t][rci] = dynamic_cast<TH2D*>(f2DDiffFlowCorrelationsProList[t]->FindObject(Form("%s, %s, %s",s2DDiffFlowName.Data(),typeFlag[t].Data(),differentialFlowIndex[rci].Data())));
   if(!f2DDiffFlow[t][rci])
   {
    printf("\n WARNING (QC): f2DDiffFlow[%i][%i] is NULL in AFAWQC::GPF2DFH() !!!!\n\n",t,rci);
    exit(0);   
   } else
     {
      this->Set2DDiffFlow(f2DDiffFlow[t][rci],t,rci);
     } 
  } // end of for(Int_t rci=0;rci<4;rci++) // reduced correlation index
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI 
  
} // end of void AliFlowAnalysisWithQCumulants::GetPointersFor2DDiffFlowHistograms()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::GetPointersForOtherDiffCorrelators()
{
 // Get pointers for other differential correlators.
 //  a) Get pointer to list with other differential correlators;
 //  b) Declare local flags;
 //  c) Get pointers to other differential profiles.

 if(!fCalculateDiffFlow){return;} // TBI: This must eventually be moved somewhere else 

 // a) Get pointer to list with other differential correlators:
 fOtherDiffCorrelatorsList = dynamic_cast<TList*>(fHistList->FindObject("Other differential correlators"));  
 if(!fOtherDiffCorrelatorsList)
 { 
  printf("\n WARNING (QC): fOtherDiffCorrelatorsList is NULL in AFAWQC::GPFDFH() !!!!\n\n");
  exit(0);
 }
 
 // b) Declare local flags: // (to be improved - promoted to data members)
 TString typeFlag[2] = {"RP","POI"}; 
 TString ptEtaFlag[2] = {"p_{T}","#eta"};
 TString sinCosFlag[2] = {"sin","cos"}; 
  
 // c) Get pointers to other differential profiles:
 TString otherDiffCorrelatorsName = "fOtherDiffCorrelators";
 otherDiffCorrelatorsName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // typeFlag (0 = RP, 1 = POI)
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t sc=0;sc<2;sc++) // sin or cos
   {
    for(Int_t ci=0;ci<1;ci++) // correlator index
    {
     fOtherDiffCorrelators[t][pe][sc][ci] = dynamic_cast<TProfile*>(fOtherDiffCorrelatorsList->FindObject(Form("%s, %s, %s, %s, ci = %d",otherDiffCorrelatorsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),ci+1))); 
     if(!fOtherDiffCorrelators[t][pe][sc][ci])
     {
      printf("\n WARNING (QC): fOtherDiffCorrelators[%i][%i][%i][%i] is NULL in AFAWQC::GPFODC() !!!!\n\n",t,pe,sc,ci);
      exit(0);       
     } else
       {
        this->SetOtherDiffCorrelators(fOtherDiffCorrelators[t][pe][sc][ci],t,pe,sc,ci);     
       } 
    } // end of for(Int_t ci=0;ci<1;ci++) // correlator index
   } // end of for(Int_t sc=0;sc<2;sc++) // sin or cos
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
 } // end of for(Int_t t=0;t<2;t++) // typeFlag (0 = RP, 1 = POI)
  
} // end of void AliFlowAnalysisWithQCumulants::GetPointersForOtherDiffCorrelators()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::GetPointersForDiffFlowHistograms()
{
 // Get pointer to all objects relevant for differential flow.
 //  a) Get pointer to base list for differential flow fDiffFlowList;
 //  b) Get pointer to profile fDiffFlowFlags holding all flags for differential flow. Access and set some flags;
 //  c) Get pointers to nested lists fDiffFlowListProfiles and fDiffFlowListResults;
 //  d) Define flags locally (to be improved: should I promote these flags to data members?);
 //  e) Get pointers to all nested lists in fDiffFlowListProfiles and to profiles which they hold;
 //  f) Get pointers to all nested lists in fDiffFlowListResults and to histograms which they hold.
 
 // a) Get pointer to base list for differential flow fDiffFlowList:
 fDiffFlowList = dynamic_cast<TList*>(fHistList->FindObject("Differential Flow"));  
 if(!fDiffFlowList)
 { 
  printf("\n WARNING (QC): fDiffFlowList is NULL in AFAWQC::GPFDFH() !!!!\n\n");
  exit(0);
 }
 
 // b) Get pointer to profile fDiffFlowFlags holding all flags for differential flow. Access and set some flags:
 TString diffFlowFlagsName = "fDiffFlowFlags";
 diffFlowFlagsName += fAnalysisLabel->Data();
 fDiffFlowFlags = dynamic_cast<TProfile*>(fDiffFlowList->FindObject(diffFlowFlagsName.Data()));
 if(fDiffFlowFlags)
 {
  this->SetCalculateDiffFlow((Bool_t)fDiffFlowFlags->GetBinContent(1)); // to be improved - hardwired 1
  this->SetCalculateDiffFlowVsEta((Bool_t)fDiffFlowFlags->GetBinContent(6)); // to be improved - hardwired 6
 } else
   {
    printf("\n WARNING (QC): fDiffFlowFlags is NULL in AFAWQC::GPFDFH() !!!!\n\n");
    printf("\n             Flags in method Finish() are wrong.\n\n");
    exit(0);
   } 
   
 if(!fCalculateDiffFlow){return;} // IMPORTANT: do not move this anywhere above in this method (to be improved)   
  
 // c) Get pointers to nested lists fDiffFlowListProfiles and fDiffFlowListResults:
 //  List holding nested lists holding profiles:
 TList *diffFlowListProfiles = NULL;
 diffFlowListProfiles = dynamic_cast<TList*>(fDiffFlowList->FindObject("Profiles"));
 if(!diffFlowListProfiles)
 { 
  printf("\n WARNING (QC): diffFlowListProfiles is NULL in AFAWQC::GPFDFH() !!!!\n\n");
  exit(0);
 }
 //  List holding nested lists holding histograms with final results:
 TList *diffFlowListResults = NULL;
 diffFlowListResults = dynamic_cast<TList*>(fDiffFlowList->FindObject("Results"));
 if(!diffFlowListResults)
 { 
  printf("\n WARNING (QC): diffFlowListResults is NULL in AFAWQC::GPFDFH() !!!!\n\n");
  exit(0);
 }
 
 // d) Define flags locally (to be improved: should I promote these flags to data members?):
 TString typeFlag[2] = {"RP","POI"}; 
 TString ptEtaFlag[2] = {"p_{T}","#eta"};
 TString powerFlag[2] = {"linear","quadratic"};
 TString sinCosFlag[2] = {"sin","cos"};
 TString differentialCumulantIndex[4] = {"QC{2'}","QC{4'}","QC{6'}","QC{8'}"};  
 TString differentialFlowIndex[4] = {"v'{2}","v'{4}","v'{6}","v'{8}"};  
 TString reducedCorrelationIndex[4] = {"<2'>","<4'>","<6'>","<8'>"};
 TString reducedSquaredCorrelationIndex[4] = {"<2'>^{2}","<4'>^{2}","<6'>^{2}","<8'>^{2}"}; 
 TString mixedCorrelationIndex[8] = {"<2>","<2'>","<4>","<4'>","<6>","<6'>","<8>","<8'>"};
 TString covarianceName[5] = {"Cov(<2>,<2'>)","Cov(<2>,<4'>)","Cov(<4>,<2'>)","Cov(<4>,<4'>)","Cov(<2'>,<4'>)"}; 
  
 // e) Get pointers to all nested lists in fDiffFlowListProfiles and to profiles which they hold:
 // correlations:
 TList *diffFlowCorrelationsProList[2][2] = {{NULL}};
 TString diffFlowCorrelationsProName = "fDiffFlowCorrelationsPro";
 diffFlowCorrelationsProName += fAnalysisLabel->Data();
 TProfile *diffFlowCorrelationsPro[2][2][4] = {{{NULL}}}; 
 // squared correlations:  
 TString diffFlowSquaredCorrelationsProName = "fDiffFlowSquaredCorrelationsPro";
 diffFlowSquaredCorrelationsProName += fAnalysisLabel->Data(); 
 TProfile *diffFlowSquaredCorrelationsPro[2][2][4] = {{{NULL}}};  
 // products of correlations:
 TList *diffFlowProductOfCorrelationsProList[2][2] = {{NULL}};
 TString diffFlowProductOfCorrelationsProName = "fDiffFlowProductOfCorrelationsPro";
 diffFlowProductOfCorrelationsProName += fAnalysisLabel->Data();  
 TProfile *diffFlowProductOfCorrelationsPro[2][2][8][8] = {{{{NULL}}}};   
 // corrections:
 TList *diffFlowCorrectionsProList[2][2] = {{NULL}};
 TString diffFlowCorrectionTermsForNUAProName = "fDiffFlowCorrectionTermsForNUAPro";
 diffFlowCorrectionTermsForNUAProName += fAnalysisLabel->Data();  
 TProfile *diffFlowCorrectionTermsForNUAPro[2][2][2][10] = {{{{NULL}}}};   
 for(Int_t t=0;t<2;t++)
 {
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++)
  {
   diffFlowCorrelationsProList[t][pe] = dynamic_cast<TList*>(diffFlowListProfiles->FindObject(Form("Profiles with correlations (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data())));
   if(!diffFlowCorrelationsProList[t][pe])
   { 
    cout<<"WARNING: diffFlowCorrelationsProList[t][pe] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
    cout<<"t = "<<t<<endl;
    cout<<"pe = "<<pe<<endl;
    exit(0);
   }
   for(Int_t ci=0;ci<4;ci++) // correlation index
   {
    // reduced correlations:
    diffFlowCorrelationsPro[t][pe][ci] = dynamic_cast<TProfile*>(diffFlowCorrelationsProList[t][pe]->FindObject(Form("%s, %s, %s, %s",diffFlowCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[ci].Data())));
    if(diffFlowCorrelationsPro[t][pe][ci])
    {
     this->SetDiffFlowCorrelationsPro(diffFlowCorrelationsPro[t][pe][ci],t,pe,ci);
    } else
      {
       cout<<"WARNING: diffFlowCorrelationsPro[t][pe][ci] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
       cout<<"t  = "<<t<<endl;
       cout<<"pe = "<<pe<<endl;   
       cout<<"ci = "<<ci<<endl;
      }     
    // reduced squared correlations:
    diffFlowSquaredCorrelationsPro[t][pe][ci] = dynamic_cast<TProfile*>(diffFlowCorrelationsProList[t][pe]->FindObject(Form("%s, %s, %s, %s",diffFlowSquaredCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedSquaredCorrelationIndex[ci].Data())));
    if(diffFlowSquaredCorrelationsPro[t][pe][ci])
    {
     this->SetDiffFlowSquaredCorrelationsPro(diffFlowSquaredCorrelationsPro[t][pe][ci],t,pe,ci);
    } else
      {
       cout<<"WARNING: diffFlowSquaredCorrelationsPro[t][pe][ci] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
       cout<<"t  = "<<t<<endl;
       cout<<"pe = "<<pe<<endl;   
       cout<<"ci = "<<ci<<endl;
      }       
   } // end of for(Int_t ci=0;ci<4;ci++) // correlation index  
   // products of correlations:    
   diffFlowProductOfCorrelationsProList[t][pe] = dynamic_cast<TList*>(diffFlowListProfiles->FindObject(Form("Profiles with products of correlations (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()))); 
   if(!diffFlowProductOfCorrelationsProList[t][pe])
   { 
    cout<<"WARNING: ddiffFlowProductOfCorrelationsProList[t][pe] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
    cout<<"t = "<<t<<endl;
    cout<<"pe = "<<pe<<endl;
    exit(0);
   }
   for(Int_t mci1=0;mci1<8;mci1++) // mixed correlation index
   {
    for(Int_t mci2=mci1+1;mci2<8;mci2++) // mixed correlation index
    {
     diffFlowProductOfCorrelationsPro[t][pe][mci1][mci2] = dynamic_cast<TProfile*>(diffFlowProductOfCorrelationsProList[t][pe]->FindObject(Form("%s, %s, %s, %s, %s",diffFlowProductOfCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),mixedCorrelationIndex[mci1].Data(),mixedCorrelationIndex[mci2].Data())));
     if(diffFlowProductOfCorrelationsPro[t][pe][mci1][mci2])
     {
      this->SetDiffFlowProductOfCorrelationsPro(diffFlowProductOfCorrelationsPro[t][pe][mci1][mci2],t,pe,mci1,mci2);
     } else
       {
        cout<<"WARNING: diffFlowProductOfCorrelationsPro[t][pe][ci] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
        cout<<"t    = "<<t<<endl;
        cout<<"pe   = "<<pe<<endl;   
        cout<<"mci1 = "<<mci1<<endl;
        cout<<"mci2 = "<<mci2<<endl;
       }
     if(mci1%2 == 0) mci2++; // products which DO NOT include reduced correlations are not stored here
    } // end of for(Int_t mci2=mci1+1;mci2<8;mci2++) // mixed correlation index
   } // end of for(Int_t mci1=0;mci1<8;mci1++) // mixed correlation index    
   // corrections:
   diffFlowCorrectionsProList[t][pe] = dynamic_cast<TList*>(diffFlowListProfiles->FindObject(Form("Profiles with correction terms for NUA (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data())));
   if(!diffFlowCorrectionsProList[t][pe])
   { 
    cout<<"WARNING: diffFlowCorrectionsProList[t][pe] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
    cout<<"t = "<<t<<endl;
    cout<<"pe = "<<pe<<endl;
    exit(0);
   }
   // correction terms for NUA:
   for(Int_t sc=0;sc<2;sc++) // sin or cos
   {
    for(Int_t cti=0;cti<9;cti++) // correction term index
    {
     diffFlowCorrectionTermsForNUAPro[t][pe][sc][cti] = dynamic_cast<TProfile*>(diffFlowCorrectionsProList[t][pe]->FindObject(Form("%s, %s, %s, %s, cti = %d",diffFlowCorrectionTermsForNUAProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),cti+1)));
     if(diffFlowCorrectionTermsForNUAPro[t][pe][sc][cti])
     {
      this->SetDiffFlowCorrectionTermsForNUAPro(diffFlowCorrectionTermsForNUAPro[t][pe][sc][cti],t,pe,sc,cti);
     } else
       {
        cout<<"WARNING: diffFlowCorrectionTermsForNUAPro[t][pe][sc][cti] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
        cout<<"t   = "<<t<<endl;
        cout<<"pe  = "<<pe<<endl;   
        cout<<"sc  = "<<sc<<endl;
        cout<<"cti = "<<cti<<endl;
       }    
    } // end of for(Int_t cti=0;cti<9;cti++) // correction term index
   } // end of for(Int_t sc=0;sc<2;sc++) // sin or cos
   // ...
  } // end of for(Int_t pe=0;pe<2;pe++)
 } // end of for(Int_t t=0;t<2;t++)
  
 // f) Get pointers to all nested lists in fDiffFlowListResults and to histograms which they hold:
 // reduced correlations:
 TList *diffFlowCorrelationsHistList[2][2] = {{NULL}};
 TString diffFlowCorrelationsHistName = "fDiffFlowCorrelationsHist";
 diffFlowCorrelationsHistName += fAnalysisLabel->Data();  
 TH1D *diffFlowCorrelationsHist[2][2][4] = {{{NULL}}};
 // corrections for NUA:
 TList *diffFlowCorrectionsHistList[2][2] = {{NULL}};
 TString diffFlowCorrectionTermsForNUAHistName = "fDiffFlowCorrectionTermsForNUAHist";
 diffFlowCorrectionTermsForNUAHistName += fAnalysisLabel->Data();  
 TH1D *diffFlowCorrectionTermsForNUAHist[2][2][2][10] = {{{{NULL}}}};
 // differential Q-cumulants:
 TList *diffFlowCumulantsHistList[2][2] = {{NULL}};
 TString diffFlowCumulantsName = "fDiffFlowCumulants";
 diffFlowCumulantsName += fAnalysisLabel->Data();  
 TH1D *diffFlowCumulants[2][2][4] = {{{NULL}}};
 // detector bias to differential Q-cumulants:
 TList *diffFlowDetectorBiasHistList[2][2] = {{NULL}};
 TString diffFlowDetectorBiasName = "fDiffFlowDetectorBias";
 diffFlowDetectorBiasName += fAnalysisLabel->Data();  
 TH1D *diffFlowDetectorBias[2][2][4] = {{{NULL}}}; 
 // differential flow estimates from Q-cumulants:
 TList *diffFlowHistList[2][2] = {{NULL}};
 TString diffFlowName = "fDiffFlow";
 diffFlowName += fAnalysisLabel->Data();  
 TH1D *diffFlow[2][2][4] = {{{NULL}}};
 // differential covariances:
 TList *diffFlowCovariancesHistList[2][2] = {{NULL}};
 TString diffFlowCovariancesName = "fDiffFlowCovariances";
 diffFlowCovariancesName += fAnalysisLabel->Data();  
 TH1D *diffFlowCovariances[2][2][5] = {{{NULL}}};
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   // reduced correlations:
   diffFlowCorrelationsHistList[t][pe] = dynamic_cast<TList*>(diffFlowListResults->FindObject(Form("Correlations (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data())));
   if(!diffFlowCorrelationsHistList[t][pe])
   { 
    cout<<"WARNING: diffFlowCorrelationsHistList[t][pe] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
    cout<<"t = "<<t<<endl;
    cout<<"pe = "<<pe<<endl;
    exit(0);
   }
   for(Int_t index=0;index<4;index++) 
   {
    diffFlowCorrelationsHist[t][pe][index] = dynamic_cast<TH1D*>(diffFlowCorrelationsHistList[t][pe]->FindObject(Form("%s, %s, %s, %s",diffFlowCorrelationsHistName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[index].Data())));
    if(diffFlowCorrelationsHist[t][pe][index])
    {
     this->SetDiffFlowCorrelationsHist(diffFlowCorrelationsHist[t][pe][index],t,pe,index);
    } else 
      {
       cout<<"WARNING: diffFlowCorrelationsHist[t][pe][index] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
       cout<<"t     = "<<t<<endl;
       cout<<"pe    = "<<pe<<endl;
       cout<<"index = "<<index<<endl;
       exit(0);       
      } 
   } // end of for(Int_t index=0;index<4;index++)
   // corrections:
   diffFlowCorrectionsHistList[t][pe] = dynamic_cast<TList*>(diffFlowListResults->FindObject(Form("Histograms with correction terms for NUA (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data())));
   if(!diffFlowCorrectionsHistList[t][pe])
   { 
    cout<<"WARNING: diffFlowCorrectionsHistList[t][pe] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
    cout<<"t = "<<t<<endl;
    cout<<"pe = "<<pe<<endl;
    exit(0);
   }
   // correction terms for NUA:
   for(Int_t sc=0;sc<2;sc++) // sin or cos
   {
    for(Int_t cti=0;cti<9;cti++) // correction term index
    {
     diffFlowCorrectionTermsForNUAHist[t][pe][sc][cti] = dynamic_cast<TH1D*>(diffFlowCorrectionsHistList[t][pe]->FindObject(Form("%s, %s, %s, %s, cti = %d",diffFlowCorrectionTermsForNUAHistName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),cti+1)));
     if(diffFlowCorrectionTermsForNUAHist[t][pe][sc][cti])
     {
      this->SetDiffFlowCorrectionTermsForNUAHist(diffFlowCorrectionTermsForNUAHist[t][pe][sc][cti],t,pe,sc,cti);
     } else
       {
        cout<<"WARNING: diffFlowCorrectionTermsForNUAHist[t][pe][sc][cti] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
        cout<<"t   = "<<t<<endl;
        cout<<"pe  = "<<pe<<endl;   
        cout<<"sc  = "<<sc<<endl;
        cout<<"cti = "<<cti<<endl;
       }    
    } // end of for(Int_t cti=0;cti<9;cti++) // correction term index
   } // end of for(Int_t sc=0;sc<2;sc++) // sin or cos
   // ...
   // differential Q-cumulants:
   diffFlowCumulantsHistList[t][pe] = dynamic_cast<TList*>(diffFlowListResults->FindObject(Form("Differential Q-cumulants (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data())));
   if(!diffFlowCumulantsHistList[t][pe])
   { 
    cout<<"WARNING: diffFlowCumulantsHistList[t][pe] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
    cout<<"t  = "<<t<<endl;
    cout<<"pe = "<<pe<<endl;
    exit(0);
   }
   for(Int_t index=0;index<4;index++) 
   {
    diffFlowCumulants[t][pe][index] = dynamic_cast<TH1D*>(diffFlowCumulantsHistList[t][pe]->FindObject(Form("%s, %s, %s, %s",diffFlowCumulantsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),differentialCumulantIndex[index].Data())));
    if(diffFlowCumulants[t][pe][index])
    {
     this->SetDiffFlowCumulants(diffFlowCumulants[t][pe][index],t,pe,index);
    } else 
      {
       cout<<"WARNING: diffFlowCumulants[t][pe][index] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
       cout<<"t     = "<<t<<endl;
       cout<<"pe    = "<<pe<<endl;
       cout<<"index = "<<index<<endl;
       exit(0);       
      } 
   } // end of for(Int_t index=0;index<4;index++)
   // Detector bias to differential Q-cumulants:
   diffFlowDetectorBiasHistList[t][pe] = dynamic_cast<TList*>(diffFlowListResults->FindObject(Form("Detector bias (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data())));
   if(!diffFlowDetectorBiasHistList[t][pe])
   { 
    cout<<"WARNING: diffFlowDetectorBiasHistList[t][pe] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
    cout<<"t  = "<<t<<endl;
    cout<<"pe = "<<pe<<endl;
    exit(0);
   }
   for(Int_t index=0;index<4;index++) 
   {
    diffFlowDetectorBias[t][pe][index] = dynamic_cast<TH1D*>(diffFlowDetectorBiasHistList[t][pe]->FindObject(Form("%s, %s, %s, %s",diffFlowDetectorBiasName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),differentialCumulantIndex[index].Data())));
    if(diffFlowDetectorBias[t][pe][index])
    {
     this->SetDiffFlowDetectorBias(diffFlowDetectorBias[t][pe][index],t,pe,index);
    } else 
      {
       cout<<"WARNING: diffFlowDetectorBias[t][pe][index] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
       cout<<"t     = "<<t<<endl;
       cout<<"pe    = "<<pe<<endl;
       cout<<"index = "<<index<<endl;
       exit(0);       
      } 
   } // end of for(Int_t index=0;index<4;index++)
   // differential flow estimates from Q-cumulants:
   diffFlowHistList[t][pe] = dynamic_cast<TList*>(diffFlowListResults->FindObject(Form("Differential flow (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data())));
   if(!diffFlowHistList[t][pe])
   { 
    cout<<"WARNING: diffFlowHistList[t][pe] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
    cout<<"t  = "<<t<<endl;
    cout<<"pe = "<<pe<<endl;
    exit(0);
   }
   for(Int_t index=0;index<4;index++) 
   {
    diffFlow[t][pe][index] = dynamic_cast<TH1D*>(diffFlowHistList[t][pe]->FindObject(Form("%s, %s, %s, %s",diffFlowName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),differentialFlowIndex[index].Data())));
    if(diffFlow[t][pe][index])
    {
     this->SetDiffFlow(diffFlow[t][pe][index],t,pe,index);
    } else 
      {
       cout<<"WARNING: diffFlow[t][pe][index] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
       cout<<"t     = "<<t<<endl;
       cout<<"pe    = "<<pe<<endl;
       cout<<"index = "<<index<<endl;
       exit(0);       
      } 
   } // end of for(Int_t index=0;index<4;index++)
   // differential covariances:
   diffFlowCovariancesHistList[t][pe] = dynamic_cast<TList*>(diffFlowListResults->FindObject(Form("Covariances of correlations (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data())));
   if(!diffFlowCovariancesHistList[t][pe])
   { 
    cout<<"WARNING: diffFlowCovariancesHistList[t][pe] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
    cout<<"t  = "<<t<<endl;
    cout<<"pe = "<<pe<<endl;
    exit(0);
   }
   for(Int_t covIndex=0;covIndex<5;covIndex++) 
   {
    diffFlowCovariances[t][pe][covIndex] = dynamic_cast<TH1D*>(diffFlowCovariancesHistList[t][pe]->FindObject(Form("%s, %s, %s, %s",diffFlowCovariancesName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),covarianceName[covIndex].Data())));
    if(diffFlowCovariances[t][pe][covIndex])
    {
     this->SetDiffFlowCovariances(diffFlowCovariances[t][pe][covIndex],t,pe,covIndex);
    } else 
      {
       cout<<"WARNING: diffFlowCovariances[t][pe][covIndex] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
       cout<<"t        = "<<t<<endl;
       cout<<"pe       = "<<pe<<endl;
       cout<<"covIndex = "<<covIndex<<endl;
       exit(0);       
      } 
   } // end of for(Int_t covIndex=0;covIndex<5;covIndex++) // covariance index    
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI 
 // sum of event weights for reduced correlations:
 TList *diffFlowSumOfEventWeightsHistList[2][2][2] = {{{NULL}}};
 TString diffFlowSumOfEventWeightsName = "fDiffFlowSumOfEventWeights";
 diffFlowSumOfEventWeightsName += fAnalysisLabel->Data();  
 TH1D *diffFlowSumOfEventWeights[2][2][2][4] = {{{{NULL}}}};
 for(Int_t t=0;t<2;t++) // type is RP or POI
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  { 
   for(Int_t p=0;p<2;p++) // power of event weights is either 1 or 2
   {
    diffFlowSumOfEventWeightsHistList[t][pe][p] = dynamic_cast<TList*>(diffFlowListResults->FindObject(Form("Sum of %s event weights (%s, %s)",powerFlag[p].Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data())));
    if(!diffFlowSumOfEventWeightsHistList[t][pe][p])
    { 
     cout<<"WARNING: diffFlowSumOfEventWeightsHistList[t][pe][p] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
     cout<<"t     = "<<t<<endl;
     cout<<"pe    = "<<pe<<endl;
     cout<<"power = "<<p<<endl;
     exit(0);
    }
    for(Int_t ew=0;ew<4;ew++) // index of reduced correlation
    {
     diffFlowSumOfEventWeights[t][pe][p][ew] = dynamic_cast<TH1D*>(diffFlowSumOfEventWeightsHistList[t][pe][p]->FindObject(Form("%s, %s, %s, %s, %s",diffFlowSumOfEventWeightsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),powerFlag[p].Data(),reducedCorrelationIndex[ew].Data())));    
     if(diffFlowSumOfEventWeights[t][pe][p][ew])
     {
      this->SetDiffFlowSumOfEventWeights(diffFlowSumOfEventWeights[t][pe][p][ew],t,pe,p,ew);
     } else 
       {
        cout<<"WARNING: diffFlowSumOfEventWeights[t][pe][p][ew] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
        cout<<"t     = "<<t<<endl;
        cout<<"pe    = "<<pe<<endl;
        cout<<"power = "<<p<<endl;
        cout<<"ew    = "<<ew<<endl;
        exit(0);       
       } 
    }
   } // end of for(Int_t p=0;p<2;p++) // power of event weights is either 1 or 2
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
 } // end of for(Int_t t=0;t<2;t++) // type is RP or POI
 //  
 TList *diffFlowSumOfProductOfEventWeightsHistList[2][2] = {{NULL}};
 TString diffFlowSumOfProductOfEventWeightsName = "fDiffFlowSumOfProductOfEventWeights";
 diffFlowSumOfProductOfEventWeightsName += fAnalysisLabel->Data();  
 TH1D *diffFlowSumOfProductOfEventWeights[2][2][8][8] = {{{{NULL}}}};
 for(Int_t t=0;t<2;t++) // type is RP or POI
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  { 
   diffFlowSumOfProductOfEventWeightsHistList[t][pe] = dynamic_cast<TList*>(diffFlowListResults->FindObject(Form("Sum of products of event weights (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data())));
   if(!diffFlowSumOfProductOfEventWeightsHistList[t][pe])
   { 
    cout<<"WARNING: diffFlowSumOfProductOfEventWeightsHistList[t][pe] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
    cout<<"t     = "<<t<<endl;
    cout<<"pe    = "<<pe<<endl;
    exit(0);
   }
   for(Int_t mci1=0;mci1<8;mci1++) // mixed correlation index
   {
    for(Int_t mci2=mci1+1;mci2<8;mci2++) // mixed correlation index
    {
     diffFlowSumOfProductOfEventWeights[t][pe][mci1][mci2] = dynamic_cast<TH1D*>(diffFlowSumOfProductOfEventWeightsHistList[t][pe]->FindObject(Form("%s, %s, %s, %s, %s",diffFlowSumOfProductOfEventWeightsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),mixedCorrelationIndex[mci1].Data(),mixedCorrelationIndex[mci2].Data())));    
      if(diffFlowSumOfProductOfEventWeights[t][pe][mci1][mci2])
      {
       this->SetDiffFlowSumOfProductOfEventWeights(diffFlowSumOfProductOfEventWeights[t][pe][mci1][mci2],t,pe,mci1,mci2);
      } else 
        {
         cout<<"WARNING: diffFlowSumOfProductOfEventWeights[t][pe][mci1][mci2] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
         cout<<"t    = "<<t<<endl;
         cout<<"pe   = "<<pe<<endl;
         cout<<"mci1 = "<<mci1<<endl;
         cout<<"mci2 = "<<mci2<<endl;
         exit(0);       
        } 
     if(mci1%2 == 0) mci2++; // products which DO NOT include reduced correlations are not stored here
    } // end of for(Int_t mci2=mci1+1;mci2<8;mci2++) // mixed correlation index
   } // end of for(Int_t mci1=0;mci1<8;mci1++) // mixed correlation index
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
 } // end of for(Int_t t=0;t<2;t++) // type is RP or POI

} // end void AliFlowAnalysisWithQCumulants::GetPointersForDiffFlowHistograms()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::BookEverythingFor2DDifferentialFlow()
{
 // Book all objects needed for 2D differential flow.
 //  a) Define flags locally (to be improved: should I promote flags to data members?);
 //  b) Book e-b-e quantities;
 //  c) Book 2D profiles;
 //  d) Book 2D histograms.
 
 if(!fCalculate2DDiffFlow){return;}

 // a) Define flags locally (to be improved: should I promote flags to data members?):
 TString typeFlag[2] = {"RP","POI"}; 
 TString reducedCorrelationIndex[4] = {"<2'>","<4'>","<6'>","<8'>"};
 TString differentialCumulantIndex[4] = {"QC{2'}","QC{4'}","QC{6'}","QC{8'}"};  
 TString differentialFlowIndex[4] = {"v'{2}","v'{4}","v'{6}","v'{8}"};  
  
 // b) Book e-b-e quantities: 
 TProfile2D styleRe("typeMultiplePowerRe","typeMultiplePowerRe",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 TProfile2D styleIm("typeMultiplePowerIm","typeMultiplePowerIm",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 for(Int_t t=0;t<3;t++) // typeFlag (0 = RP, 1 = POI, 2 = RP&&POI )
 { 
  for(Int_t m=0;m<4;m++)
  {
   for(Int_t k=0;k<9;k++)
   {
    fReRPQ2dEBE[t][m][k] = (TProfile2D*)styleRe.Clone(Form("typeFlag%dmultiple%dpower%dRe",t,m,k)); 
    fImRPQ2dEBE[t][m][k] = (TProfile2D*)styleIm.Clone(Form("typeFlag%dmultiple%dpower%dIm",t,m,k));
   }
  } 
 } 
 TProfile2D styleS("typePower","typePower",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 for(Int_t t=0;t<3;t++) // typeFlag (0 = RP, 1 = POI, 2 = RP&&POI )
 { 
  for(Int_t k=0;k<9;k++)
  {
   fs2dEBE[t][k] = (TProfile2D*)styleS.Clone(Form("typeFlag%dpower%d",t,k));
  }
 }

 // c) Book 2D profiles:
 TString s2DDiffFlowCorrelationsProName = "f2DDiffFlowCorrelationsPro";
 s2DDiffFlowCorrelationsProName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t rci=0;rci<4;rci++) // reduced correlation index
  {
   f2DDiffFlowCorrelationsPro[t][rci] = new TProfile2D(Form("%s, %s, %s",s2DDiffFlowCorrelationsProName.Data(),typeFlag[t].Data(),reducedCorrelationIndex[rci].Data()),Form("%s, %s, %s",s2DDiffFlowCorrelationsProName.Data(),typeFlag[t].Data(),reducedCorrelationIndex[rci].Data()),fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax,"");
   f2DDiffFlowCorrelationsPro[t][rci]->Sumw2();
   f2DDiffFlowCorrelationsPro[t][rci]->SetXTitle("p_{t}");
   f2DDiffFlowCorrelationsPro[t][rci]->SetYTitle("#eta");
   f2DDiffFlowCorrelationsProList[t]->Add(f2DDiffFlowCorrelationsPro[t][rci]); 
  } // end of for(Int_t rci=0;rci<4;rci++) // correlation index
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POIs

 // d) Book 2D histograms:
 TString s2DDiffFlowCumulantsName = "f2DDiffFlowCumulants";
 s2DDiffFlowCumulantsName += fAnalysisLabel->Data();
 TString s2DDiffFlowName = "f2DDiffFlow";
 s2DDiffFlowName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t rci=0;rci<4;rci++) // reduced correlation index
  {
   // 2D diferential cumulants:
   f2DDiffFlowCumulants[t][rci] = new TH2D(Form("%s, %s, %s",s2DDiffFlowCumulantsName.Data(),typeFlag[t].Data(),differentialCumulantIndex[rci].Data()),Form("%s, %s, %s",s2DDiffFlowCumulantsName.Data(),typeFlag[t].Data(),differentialCumulantIndex[rci].Data()),fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
   f2DDiffFlowCumulants[t][rci]->SetXTitle("p_{t}");
   f2DDiffFlowCumulants[t][rci]->SetYTitle("#eta");
   f2DDiffFlowCorrelationsProList[t]->Add(f2DDiffFlowCumulants[t][rci]); //  to be improved - moved to another list 
   // 2D differential flow:
   f2DDiffFlow[t][rci] = new TH2D(Form("%s, %s, %s",s2DDiffFlowName.Data(),typeFlag[t].Data(),differentialFlowIndex[rci].Data()),Form("%s, %s, %s",s2DDiffFlowName.Data(),typeFlag[t].Data(),differentialFlowIndex[rci].Data()),fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
   f2DDiffFlow[t][rci]->SetXTitle("p_{t}");
   f2DDiffFlow[t][rci]->SetYTitle("#eta");
   f2DDiffFlowCorrelationsProList[t]->Add(f2DDiffFlow[t][rci]); //  to be improved - moved to another list 
  } // end of for(Int_t rci=0;rci<4;rci++) // correlation index
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POIs

} // void AliFlowAnalysisWithQCumulants::BookEverythingFor2DDifferentialFlow()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::BookEverythingForDifferentialFlow()
{
 // Book all histograms and profiles needed for differential flow.
 //  a) Book profile to hold all flags for differential flow;
 //  b) Define flags locally (to be improved: should I promote flags to data members?);
 //  c) Book e-b-e quantities;
 //  d) Book profiles;
 //  e) Book histograms holding final results. 
 
 // a) Book profile to hold all flags for differential flow:
 TString diffFlowFlagsName = "fDiffFlowFlags";
 diffFlowFlagsName += fAnalysisLabel->Data();
 fDiffFlowFlags = new TProfile(diffFlowFlagsName.Data(),"Flags for differential flow",6,0,6);
 fDiffFlowFlags->SetTickLength(-0.01,"Y");
 fDiffFlowFlags->SetMarkerStyle(25);
 fDiffFlowFlags->SetLabelSize(0.04,"X");
 fDiffFlowFlags->SetLabelOffset(0.02,"Y");
 fDiffFlowFlags->GetXaxis()->SetBinLabel(1,"Calculate diff. flow"); 
 fDiffFlowFlags->GetXaxis()->SetBinLabel(2,"Particle weights");
 fDiffFlowFlags->GetXaxis()->SetBinLabel(3,"Event weights");
 fDiffFlowFlags->GetXaxis()->SetBinLabel(4,"Correct for NUA");
 fDiffFlowFlags->GetXaxis()->SetBinLabel(5,"Calculate 2D diff. flow");
 fDiffFlowFlags->GetXaxis()->SetBinLabel(6,"Calculate diff. flow vs eta");
 fDiffFlowList->Add(fDiffFlowFlags);

 if(!fCalculateDiffFlow){return;}
  
 // b) Define flags locally (to be improved: should I promote flags to data members?): 
 TString typeFlag[2] = {"RP","POI"}; 
 TString ptEtaFlag[2] = {"p_{T}","#eta"};
 TString powerFlag[2] = {"linear","quadratic"};
 TString sinCosFlag[2] = {"sin","cos"};
 TString differentialCumulantIndex[4] = {"QC{2'}","QC{4'}","QC{6'}","QC{8'}"};  
 TString differentialFlowIndex[4] = {"v'{2}","v'{4}","v'{6}","v'{8}"};  
 TString reducedCorrelationIndex[4] = {"<2'>","<4'>","<6'>","<8'>"};
 TString reducedSquaredCorrelationIndex[4] = {"<2'>^{2}","<4'>^{2}","<6'>^{2}","<8'>^{2}"};
 TString mixedCorrelationIndex[8] = {"<2>","<2'>","<4>","<4'>","<6>","<6'>","<8>","<8'>"};
 TString covarianceName[5] = {"Cov(<2>,<2'>)","Cov(<2>,<4'>)","Cov(<4>,<2'>)","Cov(<4>,<4'>)","Cov(<2'>,<4'>)"}; 
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 Double_t maxPtEta[2] = {fPtMax,fEtaMax};
   
 // c) Book e-b-e quantities:
 // Event-by-event r_{m*n,k}(pt,eta), p_{m*n,k}(pt,eta) and q_{m*n,k}(pt,eta)
 // Explanantion of notation:
 //  1.) n is harmonic, m is multiple of harmonic;
 //  2.) k is power of particle weight;
 //  3.) r_{m*n,k}(pt,eta) = Q-vector evaluated in harmonic m*n for RPs in particular (pt,eta) bin (i-th RP is weighted with w_i^k);   
 //  4.) p_{m*n,k}(pt,eta) = Q-vector evaluated in harmonic m*n for POIs in particular (pt,eta) bin 
 //                          (if i-th POI is also RP, than it is weighted with w_i^k);   
 //  5.) q_{m*n,k}(pt,eta) = Q-vector evaluated in harmonic m*n for particles which are both RPs and POIs in particular (pt,eta) bin 
 //                          (i-th RP&&POI is weighted with w_i^k)            
  
 // 1D:
 for(Int_t t=0;t<3;t++) // typeFlag (0 = RP, 1 = POI, 2 = RP && POI )
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t m=0;m<4;m++) // multiple of harmonic
   {
    for(Int_t k=0;k<9;k++) // power of particle weight
    {
     fReRPQ1dEBE[t][pe][m][k] = new TProfile(Form("TypeFlag%dpteta%dmultiple%dpower%dRe",t,pe,m,k),
                                             Form("TypeFlag%dpteta%dmultiple%dpower%dRe",t,pe,m,k),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]); 
     fImRPQ1dEBE[t][pe][m][k] = new TProfile(Form("TypeFlag%dpteta%dmultiple%dpower%dIm",t,pe,m,k),
                                             Form("TypeFlag%dpteta%dmultiple%dpower%dIm",t,pe,m,k),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]); 
    }
   }
  }
 } 
 // to be improved (add explanation of fs1dEBE[t][pe][k]):   
 for(Int_t t=0;t<3;t++) // typeFlag (0 = RP, 1 = POI, 2 = RP&&POI )
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t k=0;k<9;k++) // power of particle weight
   {
    fs1dEBE[t][pe][k] = new TProfile(Form("TypeFlag%dpteta%dmultiple%d",t,pe,k),
                                     Form("TypeFlag%dpteta%dmultiple%d",t,pe,k),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]); 
   }
  }
 }
 // correction terms for nua:
 for(Int_t t=0;t<2;t++) // typeFlag (0 = RP, 1 = POI)
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t sc=0;sc<2;sc++) // sin or cos
   {
    for(Int_t cti=0;cti<9;cti++) // correction term index
    {
     fDiffFlowCorrectionTermsForNUAEBE[t][pe][sc][cti] = new TH1D(Form("typeFlag%d pteta%d sincos%d cti%d",t,pe,sc,cti),
                                             Form("typeFlag%d pteta%d sincos%d cti%d",t,pe,sc,cti),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]); 
    }
   }
  }
 } 
 // reduced correlations e-b-e:
 TString diffFlowCorrelationsEBEName = "fDiffFlowCorrelationsEBE";
 diffFlowCorrelationsEBEName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t rci=0;rci<4;rci++) // reduced correlation index
   {
    fDiffFlowCorrelationsEBE[t][pe][rci] = new TH1D(Form("%s, %s, %s, %s",diffFlowCorrelationsEBEName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[rci].Data()),Form("%s, %s, %s, %s",diffFlowCorrelationsEBEName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[rci].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
   } // end of for(Int_t ci=0;ci<4;ci++) // correlation index
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta 
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
 // event weights for reduced correlations e-b-e:
 TString diffFlowEventWeightsForCorrelationsEBEName = "fDiffFlowEventWeightsForCorrelationsEBE";
 diffFlowEventWeightsForCorrelationsEBEName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t rci=0;rci<4;rci++) // event weight for reduced correlation index
   {
    fDiffFlowEventWeightsForCorrelationsEBE[t][pe][rci] = new TH1D(Form("%s, %s, %s, eW for %s",diffFlowEventWeightsForCorrelationsEBEName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[rci].Data()),Form("%s, %s, %s, eW for %s",diffFlowEventWeightsForCorrelationsEBEName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[rci].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
   } // end of for(Int_t ci=0;ci<4;ci++) // correlation index
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta 
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
      
 // d) Book profiles;
 // reduced correlations:
 TString diffFlowCorrelationsProName = "fDiffFlowCorrelationsPro";
 diffFlowCorrelationsProName += fAnalysisLabel->Data();
 // reduced squared correlations:
 TString diffFlowSquaredCorrelationsProName = "fDiffFlowSquaredCorrelationsPro";
 diffFlowSquaredCorrelationsProName += fAnalysisLabel->Data();
 // corrections terms:
 TString diffFlowCorrectionTermsForNUAProName = "fDiffFlowCorrectionTermsForNUAPro";
 diffFlowCorrectionTermsForNUAProName += fAnalysisLabel->Data();
 // reduced correlations:
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t rci=0;rci<4;rci++) // reduced correlation index
   {
    fDiffFlowCorrelationsPro[t][pe][rci] = new TProfile(Form("%s, %s, %s, %s",diffFlowCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[rci].Data()),Form("%s, %s, %s, %s",diffFlowCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[rci].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe],"s");
    fDiffFlowCorrelationsPro[t][pe][rci]->Sumw2();
    fDiffFlowCorrelationsPro[t][pe][rci]->SetXTitle(ptEtaFlag[pe].Data());
    fDiffFlowCorrelationsProList[t][pe]->Add(fDiffFlowCorrelationsPro[t][pe][rci]); // to be improved (add dedicated list to hold reduced correlations)
   } // end of for(Int_t rci=0;rci<4;rci++) // correlation index
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta 
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
 // reduced squared correlations:
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t rci=0;rci<4;rci++) // reduced correlation index
   {
    fDiffFlowSquaredCorrelationsPro[t][pe][rci] = new TProfile(Form("%s, %s, %s, %s",diffFlowSquaredCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedSquaredCorrelationIndex[rci].Data()),Form("%s, %s, %s, %s",diffFlowSquaredCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedSquaredCorrelationIndex[rci].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe],"s");
    fDiffFlowSquaredCorrelationsPro[t][pe][rci]->Sumw2();
    fDiffFlowSquaredCorrelationsPro[t][pe][rci]->SetXTitle(ptEtaFlag[pe].Data());
    fDiffFlowCorrelationsProList[t][pe]->Add(fDiffFlowSquaredCorrelationsPro[t][pe][rci]); // to be improved (add dedicated list to hold reduced correlations)
   } // end of for(Int_t rci=0;rci<4;rci++) // correlation index
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta 
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
 // correction terms for nua:
 for(Int_t t=0;t<2;t++) // typeFlag (0 = RP, 1 = POI)
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t sc=0;sc<2;sc++) // sin or cos
   {
    for(Int_t cti=0;cti<9;cti++) // correction term index
    {
     fDiffFlowCorrectionTermsForNUAPro[t][pe][sc][cti] = new TProfile(Form("%s, %s, %s, %s, cti = %d",diffFlowCorrectionTermsForNUAProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),cti+1),Form("%s, %s, %s, %s, cti = %d",diffFlowCorrectionTermsForNUAProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),cti+1),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]); 
     fDiffFlowCorrectionsProList[t][pe]->Add(fDiffFlowCorrectionTermsForNUAPro[t][pe][sc][cti]);
    }
   }
  }
 } 
 // Other differential correlators:
 TString otherDiffCorrelatorsName = "fOtherDiffCorrelators";
 otherDiffCorrelatorsName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // typeFlag (0 = RP, 1 = POI)
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t sc=0;sc<2;sc++) // sin or cos
   {
    for(Int_t ci=0;ci<1;ci++) // correlator index
    {
     fOtherDiffCorrelators[t][pe][sc][ci] = new TProfile(Form("%s, %s, %s, %s, ci = %d",otherDiffCorrelatorsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),ci+1),Form("%s, %s, %s, %s, ci = %d",otherDiffCorrelatorsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),ci+1),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]); 
     fOtherDiffCorrelators[t][pe][sc][ci]->Sumw2();
     fOtherDiffCorrelatorsList->Add(fOtherDiffCorrelators[t][pe][sc][ci]);
    }
   }
  }
 }  
 // e) Book histograms holding final results. 
 // reduced correlations:
 TString diffFlowCorrelationsHistName = "fDiffFlowCorrelationsHist";
 diffFlowCorrelationsHistName += fAnalysisLabel->Data();
 // corrections terms:
 TString diffFlowCorrectionTermsForNUAHistName = "fDiffFlowCorrectionTermsForNUAHist";
 diffFlowCorrectionTermsForNUAHistName += fAnalysisLabel->Data();
 // differential covariances:
 TString diffFlowCovariancesName = "fDiffFlowCovariances";
 diffFlowCovariancesName += fAnalysisLabel->Data();
 // differential Q-cumulants:
 TString diffFlowCumulantsName = "fDiffFlowCumulants";
 diffFlowCumulantsName += fAnalysisLabel->Data();
 // Detector bias to differential Q-cumulants:
 TString diffFlowDetectorBiasName = "fDiffFlowDetectorBias";
 diffFlowDetectorBiasName += fAnalysisLabel->Data();
 // differential flow:
 TString diffFlowName = "fDiffFlow";
 diffFlowName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t index=0;index<4;index++) 
   {
    // reduced correlations:
    fDiffFlowCorrelationsHist[t][pe][index] = new TH1D(Form("%s, %s, %s, %s",diffFlowCorrelationsHistName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[index].Data()),Form("%s, %s, %s, %s",diffFlowCorrelationsHistName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[index].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
    fDiffFlowCorrelationsHist[t][pe][index]->SetXTitle(ptEtaFlag[pe].Data());
    fDiffFlowCorrelationsHistList[t][pe]->Add(fDiffFlowCorrelationsHist[t][pe][index]); 
    // differential Q-cumulants:
    fDiffFlowCumulants[t][pe][index] = new TH1D(Form("%s, %s, %s, %s",diffFlowCumulantsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),differentialCumulantIndex[index].Data()),Form("%s, %s, %s, %s",diffFlowCumulantsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),differentialCumulantIndex[index].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
    fDiffFlowCumulants[t][pe][index]->SetXTitle(ptEtaFlag[pe].Data());
    fDiffFlowCumulantsHistList[t][pe]->Add(fDiffFlowCumulants[t][pe][index]); 
    // Detector bias to differential Q-cumulants:
    fDiffFlowDetectorBias[t][pe][index] = new TH1D(Form("%s, %s, %s, %s",diffFlowDetectorBiasName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),differentialCumulantIndex[index].Data()),Form("%s, %s, %s, %s",diffFlowDetectorBiasName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),differentialCumulantIndex[index].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
    fDiffFlowDetectorBias[t][pe][index]->SetXTitle(ptEtaFlag[pe].Data());
    fDiffFlowDetectorBias[t][pe][index]->SetTitle(Form("#frac{corrected}{measured} %s",differentialCumulantIndex[index].Data()));
    fDiffFlowDetectorBiasHistList[t][pe]->Add(fDiffFlowDetectorBias[t][pe][index]); 
    // differential flow estimates from Q-cumulants:
    fDiffFlow[t][pe][index] = new TH1D(Form("%s, %s, %s, %s",diffFlowName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),differentialFlowIndex[index].Data()),Form("%s, %s, %s, %s",diffFlowName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),differentialFlowIndex[index].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
    fDiffFlow[t][pe][index]->SetXTitle(ptEtaFlag[pe].Data());
    fDiffFlowHistList[t][pe]->Add(fDiffFlow[t][pe][index]); 
   } // end of for(Int_t index=0;index<4;index++) 
   for(Int_t covIndex=0;covIndex<5;covIndex++) // covariance index 
   {
    // differential covariances:
    fDiffFlowCovariances[t][pe][covIndex] = new TH1D(Form("%s, %s, %s, %s",diffFlowCovariancesName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),covarianceName[covIndex].Data()),Form("%s, %s, %s, %s",diffFlowCovariancesName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),covarianceName[covIndex].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
    fDiffFlowCovariances[t][pe][covIndex]->SetXTitle(ptEtaFlag[pe].Data());
    fDiffFlowCovariancesHistList[t][pe]->Add(fDiffFlowCovariances[t][pe][covIndex]); 
   } // end of for(Int_t covIndex=0;covIndex<5;covIndex++) // covariance index
   // products of both types of correlations: 
   TString diffFlowProductOfCorrelationsProName = "fDiffFlowProductOfCorrelationsPro";
   diffFlowProductOfCorrelationsProName += fAnalysisLabel->Data();  
   for(Int_t mci1=0;mci1<8;mci1++) // mixed correlation index
   {
    for(Int_t mci2=mci1+1;mci2<8;mci2++) // mixed correlation index
    {
     fDiffFlowProductOfCorrelationsPro[t][pe][mci1][mci2] = new TProfile(Form("%s, %s, %s, %s, %s",diffFlowProductOfCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),mixedCorrelationIndex[mci1].Data(),mixedCorrelationIndex[mci2].Data()),Form("%s, %s, %s, %s #times %s",diffFlowProductOfCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),mixedCorrelationIndex[mci1].Data(),mixedCorrelationIndex[mci2].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]); 
     fDiffFlowProductOfCorrelationsPro[t][pe][mci1][mci2]->SetXTitle(ptEtaFlag[pe].Data());
     fDiffFlowProductOfCorrelationsProList[t][pe]->Add(fDiffFlowProductOfCorrelationsPro[t][pe][mci1][mci2]); 
     if(mci1%2 == 0) mci2++; // products which DO NOT include reduced correlations are not stored here
    } // end of for(Int_t mci2=mci1+1;mci2<8;mci2++) // mixed correlation index
   } // end of for(Int_t mci1=0;mci1<8;mci1++) // mixed correlation index    
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta 
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
 // sums of event weights for reduced correlations: 
 TString diffFlowSumOfEventWeightsName = "fDiffFlowSumOfEventWeights";
 diffFlowSumOfEventWeightsName += fAnalysisLabel->Data();  
 for(Int_t t=0;t<2;t++) // type is RP or POI
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  { 
   for(Int_t p=0;p<2;p++) // power of weights is either 1 or 2
   {
    for(Int_t ew=0;ew<4;ew++) // index of reduced correlation
    {
     fDiffFlowSumOfEventWeights[t][pe][p][ew] = new TH1D(Form("%s, %s, %s, %s, %s",diffFlowSumOfEventWeightsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),powerFlag[p].Data(),reducedCorrelationIndex[ew].Data()),Form("%s, %s, %s, power = %s, %s",diffFlowSumOfEventWeightsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),powerFlag[p].Data(),reducedCorrelationIndex[ew].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]); 
     fDiffFlowSumOfEventWeights[t][pe][p][ew]->SetXTitle(ptEtaFlag[pe].Data());
     fDiffFlowSumOfEventWeightsHistList[t][pe][p]->Add(fDiffFlowSumOfEventWeights[t][pe][p][ew]); // to be improved (add dedicated list to hold all this)
    }
   }
  }
 } 
 // sum of products of event weights for both types of correlations: 
 TString diffFlowSumOfProductOfEventWeightsName = "fDiffFlowSumOfProductOfEventWeights";
 diffFlowSumOfProductOfEventWeightsName += fAnalysisLabel->Data();  
 for(Int_t t=0;t<2;t++) // type is RP or POI
 {
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  { 
   for(Int_t mci1=0;mci1<8;mci1++) // mixed correlation index
   {
    for(Int_t mci2=mci1+1;mci2<8;mci2++) // mixed correlation index
    {
     fDiffFlowSumOfProductOfEventWeights[t][pe][mci1][mci2] = new TH1D(Form("%s, %s, %s, %s, %s",diffFlowSumOfProductOfEventWeightsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),mixedCorrelationIndex[mci1].Data(),mixedCorrelationIndex[mci2].Data()),Form("%s, %s, %s, %s #times %s",diffFlowSumOfProductOfEventWeightsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),mixedCorrelationIndex[mci1].Data(),mixedCorrelationIndex[mci2].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]); 
     fDiffFlowSumOfProductOfEventWeights[t][pe][mci1][mci2]->SetXTitle(ptEtaFlag[pe].Data());
     fDiffFlowSumOfProductOfEventWeightsHistList[t][pe]->Add(fDiffFlowSumOfProductOfEventWeights[t][pe][mci1][mci2]); 
     if(mci1%2 == 0) mci2++; // products which DO NOT include reduced correlations are not stored here
    }
   }
  }
 } 
 // correction terms for nua:
 for(Int_t t=0;t<2;t++) // typeFlag (0 = RP, 1 = POI)
 { 
  for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
  {
   for(Int_t sc=0;sc<2;sc++) // sin or cos
   {
    for(Int_t cti=0;cti<9;cti++) // correction term index
    {
     fDiffFlowCorrectionTermsForNUAHist[t][pe][sc][cti] = new TH1D(Form("%s, %s, %s, %s, cti = %d",diffFlowCorrectionTermsForNUAHistName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),cti+1),Form("%s, %s, %s, %s, cti = %d",diffFlowCorrectionTermsForNUAHistName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),cti+1),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]); 
     fDiffFlowCorrectionsHistList[t][pe]->Add(fDiffFlowCorrectionTermsForNUAHist[t][pe][sc][cti]);
    }
   }
  }
 } 
          
} // end of AliFlowAnalysisWithQCumulants::BookEverythingForDifferentialFlow()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateQcumulantsCorrectedForNUAIntFlow()
{
 // Calculate generalized Q-cumulants (cumulants corrected for non-unifom acceptance).
 
 // Isotropic cumulants:
 Double_t QC2 = fIntFlowQcumulants->GetBinContent(1);
 Double_t QC2Error = fIntFlowQcumulants->GetBinError(1);
 Double_t QC4 = fIntFlowQcumulants->GetBinContent(2);
 Double_t QC4Error = fIntFlowQcumulants->GetBinError(2);
 //Double_t QC6 = fIntFlowQcumulants->GetBinContent(3);
 //Double_t QC6Error = fIntFlowQcumulants->GetBinError(3);
 //Double_t QC8 = fIntFlowQcumulants->GetBinContent(4);
 //Double_t QC8Error = fIntFlowQcumulants->GetBinError(4);
 
 // Measured 2-, 4-, 6- and 8-particle correlations:
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>>
 Double_t twoError = fIntFlowCorrelationsHist->GetBinError(1); // statistical error of <<2>>
 Double_t four = fIntFlowCorrelationsHist->GetBinContent(2); // <<4>>
 Double_t fourError = fIntFlowCorrelationsHist->GetBinError(2); // statistical error of <<4>>
 //Double_t six = fIntFlowCorrelationsHist->GetBinContent(3); // <<6>>
 //Double_t sixError = fIntFlowCorrelationsHist->GetBinError(3); // statistical error of <<6>>
 //Double_t eight = fIntFlowCorrelationsHist->GetBinContent(4); // <<8>>
 //Double_t eightError = fIntFlowCorrelationsHist->GetBinError(4); // statistical error of <<8>>
  
 // Non-isotropic terms:
 Double_t c1 = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(1); // <<cos(n*phi1)>>
 Double_t c1Error = fIntFlowCorrectionTermsForNUAHist[1]->GetBinError(1); // statistical error of <<cos(n*phi1)>>
 Double_t c2 = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(2); // <<cos(n*(phi1+phi2))>>
 Double_t c2Error = fIntFlowCorrectionTermsForNUAHist[1]->GetBinError(2); // statistical error of <<cos(n*(phi1+phi2))>>
 Double_t c3 = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(3); // <<cos(n*(phi1-phi2-phi3))>>
 Double_t c3Error = fIntFlowCorrectionTermsForNUAHist[1]->GetBinError(3); // statistical error of <<cos(n*(phi1-phi2-phi3))>>
 Double_t s1 = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(1); // <<sin(n*phi1)>>
 Double_t s1Error = fIntFlowCorrectionTermsForNUAHist[0]->GetBinError(1); // statistical error of <<sin(n*phi1)>>
 Double_t s2 = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(2); // <<sin(n*(phi1+phi2))>>
 Double_t s2Error = fIntFlowCorrectionTermsForNUAHist[0]->GetBinError(2); // statistical error of <<sin(n*(phi1+phi2))>>
 Double_t s3 = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(3); // <<sin(n*(phi1-phi2-phi3))>>
 Double_t s3Error = fIntFlowCorrectionTermsForNUAHist[0]->GetBinError(3); // statistical error of <<sin(n*(phi1-phi2-phi3))>>
 
 // Shortcuts:
 Double_t a1 = 2.*pow(c1,2.)+2.*pow(s1,2.)-two;
 Double_t a2 = 6.*pow(c1,3.)-2.*c1*c2+c3+6.*c1*pow(s1,2.)-2.*s1*s2-4.*c1*two;
 Double_t a3 = 2.*pow(s1,2.)-2.*pow(c1,2.)+c2;
 Double_t a4 = 6.*pow(s1,3.)+6.*pow(c1,2.)*s1+2.*c2*s1-2.*c1*s2-s3-4.*s1*two;
 Double_t a5 = 4.*c1*s1-s2;
 
 // Covariances (including weight dependent prefactor):
 Double_t wCov1 = 0.; // w*Cov(<2>,<cos(phi)) 
 Double_t wCov2 = 0.; // w*Cov(<2>,<sin(phi))
 Double_t wCov3 = 0.; // w*Cov(<cos(phi),<sin(phi))
 Double_t wCov4 = 0.; // w*Cov(<2>,<4>) 
 Double_t wCov5 = 0.; // w*Cov(<2>,<cos(#phi_{1}+#phi_{2})>)
 Double_t wCov6 = 0.; // w*Cov(<2>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov7 = 0.; // w*Cov(<2>,<sin(#phi_{1}+#phi_{2})>)
 Double_t wCov8 = 0.; // w*Cov(<2>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov9 = 0.; // w*Cov(<4>,<cos(#phi)>
 Double_t wCov10 = 0.; // w*Cov(<4>,<cos(#phi_{1}+#phi_{2})>)
 Double_t wCov11 = 0.; // w*Cov(<4>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov12 = 0.; // w*Cov(<4>,<sin(#phi)>
 Double_t wCov13 = 0.; // w*Cov(<4>,<sin(#phi_{1}+#phi_{2})>)
 Double_t wCov14 = 0.; // w*Cov(<4>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov15 = 0.; // w*Cov(<cos(#phi)>,<cos(#phi_{1}+#phi_{2})>)
 Double_t wCov16 = 0.; // w*Cov(<cos(#phi)>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov17 = 0.; // w*Cov(<cos(#phi)>,<sin(#phi_{1}+#phi_{2})>)
 Double_t wCov18 = 0.; // w*Cov(<cos(#phi)>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov19 = 0.; // w*Cov(<cos(#phi_{1}+#phi_{2})>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov20 = 0.; // w*Cov(<sin(#phi)>,<cos(#phi_{1}+#phi_{2})>)
 Double_t wCov21 = 0.; // w*Cov(<cos(#phi_{1}+#phi_{2})>,<sin(#phi_{1}+#phi_{2})>)
 Double_t wCov22 = 0.; // w*Cov(<cos(#phi_{1}+#phi_{2})>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov23 = 0.; // w*Cov(<sin(#phi)>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov24 = 0.; // w*Cov(<sin(#phi_{1}+#phi_{2})>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov25 = 0.; // w*Cov(<cos(#phi_{1}-#phi_{2}-#phi_{3}>,<sin(#phi_{1}-#phi_{2}-#phi_{3}>)
 Double_t wCov26 = 0.; // w*Cov(<sin(#phi)>,<sin(#phi_{1}+#phi_{2})>)
 Double_t wCov27 = 0.; // w*Cov(<sin(#phi)>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
 Double_t wCov28 = 0.; // w*Cov(<sin(#phi_{1}+#phi_{2})>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
 if(!fForgetAboutCovariances)
 {
  wCov1 = fIntFlowCovariancesNUA->GetBinContent(1); // w*Cov(<2>,<cos(phi)) 
  wCov2 = fIntFlowCovariancesNUA->GetBinContent(2); // w*Cov(<2>,<sin(phi))
  wCov3 = fIntFlowCovariancesNUA->GetBinContent(3); // w*Cov(<cos(phi),<sin(phi))
  wCov4 = fIntFlowCovariances->GetBinContent(1); // w*Cov(<2>,<4>) 
  wCov5 = fIntFlowCovariancesNUA->GetBinContent(4); // w*Cov(<2>,<cos(#phi_{1}+#phi_{2})>)
  wCov6 = fIntFlowCovariancesNUA->GetBinContent(6); // w*Cov(<2>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov7 = fIntFlowCovariancesNUA->GetBinContent(5); // w*Cov(<2>,<sin(#phi_{1}+#phi_{2})>)
  wCov8 = fIntFlowCovariancesNUA->GetBinContent(7); // w*Cov(<2>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov9 = fIntFlowCovariancesNUA->GetBinContent(8); // w*Cov(<4>,<cos(#phi)>
  wCov10 = fIntFlowCovariancesNUA->GetBinContent(10); // w*Cov(<4>,<cos(#phi_{1}+#phi_{2})>)
  wCov11 = fIntFlowCovariancesNUA->GetBinContent(12); // w*Cov(<4>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov12 = fIntFlowCovariancesNUA->GetBinContent(9); // w*Cov(<4>,<sin(#phi)>
  wCov13 = fIntFlowCovariancesNUA->GetBinContent(11); // w*Cov(<4>,<sin(#phi_{1}+#phi_{2})>)
  wCov14 = fIntFlowCovariancesNUA->GetBinContent(13); // w*Cov(<4>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov15 = fIntFlowCovariancesNUA->GetBinContent(14); // w*Cov(<cos(#phi)>,<cos(#phi_{1}+#phi_{2})>)
  wCov16 = fIntFlowCovariancesNUA->GetBinContent(16); // w*Cov(<cos(#phi)>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov17 = fIntFlowCovariancesNUA->GetBinContent(15); // w*Cov(<cos(#phi)>,<sin(#phi_{1}+#phi_{2})>)
  wCov18 = fIntFlowCovariancesNUA->GetBinContent(17); // w*Cov(<cos(#phi)>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov19 = fIntFlowCovariancesNUA->GetBinContent(23); // w*Cov(<cos(#phi_{1}+#phi_{2})>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov20 = fIntFlowCovariancesNUA->GetBinContent(18); // w*Cov(<sin(#phi)>,<cos(#phi_{1}+#phi_{2})>)
  wCov21 = fIntFlowCovariancesNUA->GetBinContent(22); // w*Cov(<cos(#phi_{1}+#phi_{2})>,<sin(#phi_{1}+#phi_{2})>)
  wCov22 = fIntFlowCovariancesNUA->GetBinContent(24); // w*Cov(<cos(#phi_{1}+#phi_{2})>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov23 = fIntFlowCovariancesNUA->GetBinContent(20); // w*Cov(<sin(#phi)>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov24 = fIntFlowCovariancesNUA->GetBinContent(25); // w*Cov(<sin(#phi_{1}+#phi_{2})>,<cos(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov25 = fIntFlowCovariancesNUA->GetBinContent(27); // w*Cov(<cos(#phi_{1}-#phi_{2}-#phi_{3}>,<sin(#phi_{1}-#phi_{2}-#phi_{3}>)
  wCov26 = fIntFlowCovariancesNUA->GetBinContent(19); // w*Cov(<sin(#phi)>,<sin(#phi_{1}+#phi_{2})>)
  wCov27 = fIntFlowCovariancesNUA->GetBinContent(21); // w*Cov(<sin(#phi)>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
  wCov28 = fIntFlowCovariancesNUA->GetBinContent(26); // w*Cov(<sin(#phi_{1}+#phi_{2})>,<sin(#phi_{1}-#phi_{2}-#phi_{3})>)
 } // end of if(!fForgetAboutCovariances)
 
 // Calculating generalized QC{2}:
 //  Generalized QC{2}:
 Double_t gQC2 = two - pow(c1,2.) - pow(s1,2.);
 if(fApplyCorrectionForNUA){fIntFlowQcumulants->SetBinContent(1,gQC2);} 
 //  Statistical error of generalized QC{2}:
 Double_t gQC2ErrorSquared = pow(twoError,2.)+4.*pow(c1,2.)*pow(c1Error,2.)
                           + 4.*pow(s1,2.)*pow(s1Error,2.)
                           - 4*c1*wCov1-4*s1*wCov2 
                           + 8.*c1*s1*wCov3;
 //  Store ratio of error squared - with/without NUA terms:
 Double_t ratioErrorSquaredQC2 = 0.;
 if(fIntFlowQcumulants->GetBinError(1)>0.)
 { 
  ratioErrorSquaredQC2 = (gQC2ErrorSquared/pow(fIntFlowQcumulants->GetBinError(1),2.));
  fIntFlowQcumulantsErrorSquaredRatio->SetBinContent(1,ratioErrorSquaredQC2);
 }
 //  If enabled, store error by including non-isotropic terms:                         
 if(fApplyCorrectionForNUA && fPropagateErrorAlsoFromNIT)
 {
  if(gQC2ErrorSquared>=0.)
  {
   fIntFlowQcumulants->SetBinError(1,pow(gQC2ErrorSquared,0.5));
  } else
    {
     fIntFlowQcumulants->SetBinError(1,0.);
     cout<<endl;
     cout<<" WARNING (QC): Statistical error of generalized QC{2} is imaginary !!!!"<<endl;
     cout<<endl;
    }   
 } // end of if(fApplyCorrectionForNUA && fPropagateErrorAlsoFromNIT)
 // Quantify detector bias to QC{2}:
 if(TMath::Abs(QC2)>0.)
 {
  fIntFlowDetectorBias->SetBinContent(1,gQC2/QC2); 
  if(QC2Error>0.)
  {
   Double_t errorSquared = gQC2ErrorSquared/pow(QC2,2.)+pow(gQC2,2.)*pow(QC2Error,2.)/pow(QC2,4.);
   if(errorSquared>0.)
   {
    fIntFlowDetectorBias->SetBinError(1,pow(errorSquared,0.5));  
   }
  }
 } // end of if(TMath::Abs(QC2)>0.)

 // Calculating generalized QC{4}:
 //  Generalized QC{4}:
 Double_t gQC4 = four-2.*pow(two,2.)
               - 4.*c1*c3+4.*s1*s3-pow(c2,2.)-pow(s2,2.)
               + 4.*c2*(pow(c1,2.)-pow(s1,2.))+8.*s2*s1*c1
               + 8.*two*(pow(c1,2.)+pow(s1,2.))-6.*pow((pow(c1,2.)+pow(s1,2.)),2.);
 if(fApplyCorrectionForNUA){fIntFlowQcumulants->SetBinContent(2,gQC4);}   
 //  Statistical error of generalized QC{4}:
 Double_t gQC4ErrorSquared = 16.*pow(a1,2.)*pow(twoError,2.)+pow(fourError,2.)+16.*pow(a2,2.)*pow(c1Error,2.)
                           + 4.*pow(a3,2.)*pow(c2Error,2.)+16.*pow(c1,2.)*pow(c3Error,2.)
                           + 16.*pow(a4,2.)*pow(s1Error,2.)+4.*pow(a5,2.)*pow(s2Error,2.)
                           + 16.*pow(s1,2.)*pow(s3Error,2.)+8.*a1*wCov4-32.*a1*a2*wCov1
                           - 16.*a3*a1*wCov5-32.*c1*a1*wCov6-32.*a1*a4*wCov2+16.*a5*a1*wCov7
                           + 32.*s1*a1*wCov8-8.*a2*wCov9-4.*a3*wCov10-8.*c1*wCov11-8.*a4*wCov12
                           + 4.*a5*wCov13+8.*s1*wCov14+16.*a3*a2*wCov15+32.*c1*a2*wCov16+32.*a2*a4*wCov3
                           - 16.*a5*a2*wCov17-32.*s1*a2*wCov18+16.*c1*a3*wCov19+16.*a3*a4*wCov20
                           - 8.*a3*a5*wCov21-16.*s1*a3*wCov22+32.*c1*a4*wCov23-16.*c1*a5*wCov24
                           - 32.*c1*s1*wCov25-16.*a5*a4*wCov26-32.*s1*a4*wCov27+16.*s1*a5*wCov28;
 //  Store ratio of error squared - with/without NUA terms:
 Double_t ratioErrorSquaredQC4 = 0.;
 if(fIntFlowQcumulants->GetBinError(2)>0.)
 { 
  ratioErrorSquaredQC4 = (gQC4ErrorSquared/pow(fIntFlowQcumulants->GetBinError(2),2.));
  fIntFlowQcumulantsErrorSquaredRatio->SetBinContent(2,ratioErrorSquaredQC4);
 }                          
 if(fApplyCorrectionForNUA && fPropagateErrorAlsoFromNIT)
 {
  if(gQC4ErrorSquared>=0.)
  {
   fIntFlowQcumulants->SetBinError(2,pow(gQC4ErrorSquared,0.5));
  } else
    {
     fIntFlowQcumulants->SetBinError(2,0.);
     cout<<endl;
     cout<<" WARNING (QC): Statistical error of generalized QC{4} is imaginary !!!!"<<endl;
     cout<<endl;
    }   
 } // end of if(fApplyCorrectionForNUA && fPropagateErrorAlsoFromNIT)
 // Quantify detector bias to QC{4}:
 if(TMath::Abs(QC4)>0.)
 {
  fIntFlowDetectorBias->SetBinContent(2,gQC4/QC4); 
  if(QC4Error>0.)
  {
   Double_t errorSquared = gQC4ErrorSquared/pow(QC4,2.)+pow(gQC4,2.)*pow(QC4Error,2.)/pow(QC4,4.);
   if(errorSquared>0.)
   {
    fIntFlowDetectorBias->SetBinError(2,pow(errorSquared,0.5));  
   }
  }
 } // end of if(TMath::Abs(QC4)>0.)


 // .... to be improved (continued for 6th and 8th order) ....            
 
     
 // versus multiplicity:
 if(fCalculateCumulantsVsM) // to be improved - propagate error for nua terms vs M
 { 
  Int_t nBins = fIntFlowCorrelationsVsMPro[0]->GetNbinsX(); // to be improved (hardwired 0) 
  Double_t value[4] = {0.}; // QCs vs M
  Double_t error[4] = {0.}; // error of QCs vs M
  Double_t dSum1[4] = {0.}; // sum value_i/(error_i)^2
  Double_t dSum2[4] = {0.}; // sum 1/(error_i)^2
  for(Int_t b=1;b<=nBins;b++)
  {
   // Measured correlations:
   two = fIntFlowCorrelationsVsMHist[0]->GetBinContent(b); // <<2>> vs M
   four = fIntFlowCorrelationsVsMHist[1]->GetBinContent(b); // <<4>> vs M
   // Isotropic cumulants:
   QC2 = two;
   QC4 = four-2.*pow(two,2.);
   // Non-isotropic terms:
   c1 = fIntFlowCorrectionTermsForNUAVsMPro[1][0]->GetBinContent(b); // <<cos(n*phi1)>>
   c2 = fIntFlowCorrectionTermsForNUAVsMPro[1][1]->GetBinContent(b); // <<cos(n*(phi1+phi2))>>
   c3 = fIntFlowCorrectionTermsForNUAVsMPro[1][2]->GetBinContent(b); // <<cos(n*(phi1-phi2-phi3))>>
   s1 = fIntFlowCorrectionTermsForNUAVsMPro[0][0]->GetBinContent(b); // <<sin(n*phi1)>>
   s2 = fIntFlowCorrectionTermsForNUAVsMPro[0][1]->GetBinContent(b); // <<sin(n*(phi1+phi2))>>
   s3 = fIntFlowCorrectionTermsForNUAVsMPro[0][2]->GetBinContent(b); // <<sin(n*(phi1-phi2-phi3))>>
   // Generalized QC{2} vs M:
   gQC2 = two - pow(c1,2.) - pow(s1,2.); 
   if(fApplyCorrectionForNUAVsM){fIntFlowQcumulantsVsM[0]->SetBinContent(b,gQC2);}   
   // Generalized QC{4} vs M:  
   gQC4 = four-2.*pow(two,2.)
                 - 4.*c1*c3+4.*s1*s3-pow(c2,2.)-pow(s2,2.)
                 + 4.*c2*(pow(c1,2.)-pow(s1,2.))+8.*s2*s1*c1
                 + 8.*two*(pow(c1,2.)+pow(s1,2.))-6.*pow((pow(c1,2.)+pow(s1,2.)),2.);
   if(fApplyCorrectionForNUAVsM){fIntFlowQcumulantsVsM[1]->SetBinContent(b,gQC4);}   
   // Detector bias vs M:
   if(TMath::Abs(QC2)>0.)
   {
    fIntFlowDetectorBiasVsM[0]->SetBinContent(b,gQC2/QC2); 
   } // end of if(TMath::Abs(QC2)>0.)
   if(TMath::Abs(QC4)>0.)
   {
    fIntFlowDetectorBiasVsM[1]->SetBinContent(b,gQC4/QC4); 
   } // end of if(TMath::Abs(QC4)>0.)  
   // Rebin in M:
   for(Int_t co=0;co<4;co++)
   {
    value[co] = fIntFlowQcumulantsVsM[co]->GetBinContent(b);
    error[co] = fIntFlowQcumulantsVsM[co]->GetBinError(b);
    if(error[co]>0.)
    {
     dSum1[co]+=value[co]/(error[co]*error[co]);
     dSum2[co]+=1./(error[co]*error[co]);
    }
   } // end of for(Int_t co=0;co<4;co++) 
  } // end of for(Int_t b=1;b<=nBins;b++)
  // Store rebinned Q-cumulants:
  if(fApplyCorrectionForNUAVsM)
  {
   for(Int_t co=0;co<4;co++)
   {
    if(dSum2[co]>0.)
    {
     fIntFlowQcumulantsRebinnedInM->SetBinContent(co+1,dSum1[co]/dSum2[co]);
     fIntFlowQcumulantsRebinnedInM->SetBinError(co+1,pow(1./dSum2[co],0.5));
    }
   } // end of for(Int_t co=0;co<4;co++)
  } // end of if(fApplyCorrectionForNUAVsM)
 } // end of if(fCalculateCumulantsVsM) 
     
} // end of void AliFlowAnalysisWithQCumulants::CalculateQcumulantsCorrectedForNUAIntFlow()
 
//================================================================================================================================

void AliFlowAnalysisWithQCumulants::FinalizeCorrectionTermsForNUAIntFlow() 
{
 // From profile fIntFlowCorrectionTermsForNUAPro[sc] access measured correction terms for NUA
 // and their spread, correctly calculate the statistical errors and store the final 
 // results and statistical errors for correction terms for NUA in histogram fIntFlowCorrectionTermsForNUAHist[sc].
 //
 // Remark: Statistical error of correction temrs is calculated as:
 //
 //          statistical error = termA * spread * termB:
 //          termA = sqrt{sum_{i=1}^{N} w^2}/(sum_{i=1}^{N} w)
 //          termB = 1/sqrt(1-termA^2)   
 
 TString sinCosFlag[2] = {"sin","cos"}; // to be improved - promore this to data member?
 TString nonisotropicTermFlag[4] = {"(n(phi1))","(n(phi1+phi2))","(n(phi1-phi2-phi3))","(n(2phi1-phi2))"}; // to be improved - hardwired 4
    
 for(Int_t sc=0;sc<2;sc++) // sin or cos correction terms 
 {
  for(Int_t ci=1;ci<=4;ci++) // correction term index (to be improved - hardwired 4)
  {
   Double_t correction = fIntFlowCorrectionTermsForNUAPro[sc]->GetBinContent(ci);
   Double_t spread = fIntFlowCorrectionTermsForNUAPro[sc]->GetBinError(ci);
   Double_t sumOfLinearEventWeights = fIntFlowSumOfEventWeightsNUA[sc][0]->GetBinContent(ci);
   Double_t sumOfQuadraticEventWeights = fIntFlowSumOfEventWeightsNUA[sc][1]->GetBinContent(ci);
   Double_t termA = 0.;
   Double_t termB = 0.;
   if(TMath::Abs(sumOfLinearEventWeights)>1.e-44)
   {
    termA = pow(sumOfQuadraticEventWeights,0.5)/sumOfLinearEventWeights;
   } else
     {
      cout<<" WARNING (QC): sumOfLinearEventWeights == 0 in AFAWQC::FCTFNIF() !!!!"<<endl;
      cout<<Form("               (for <<%s[%s]>> non-isotropic term)",sinCosFlag[sc].Data(),nonisotropicTermFlag[ci-1].Data())<<endl;
     }
   if(1.-pow(termA,2.) > 0.)
   {
    termB = 1./pow(1-pow(termA,2.),0.5);
   } else
     {
      cout<<" WARNING (QC): 1.-pow(termA,2.) <= 0 in AFAWQC::FCTFNIF() !!!!"<<endl;   
      cout<<Form("               (for <<%s[%s]>> non-isotropic term)",sinCosFlag[sc].Data(),nonisotropicTermFlag[ci-1].Data())<<endl;
     }     
   Double_t statisticalError = termA * spread * termB;
   fIntFlowCorrectionTermsForNUAHist[sc]->SetBinContent(ci,correction);
   fIntFlowCorrectionTermsForNUAHist[sc]->SetBinError(ci,statisticalError);
  } // end of for(Int_t ci=1;ci<=4;ci++) // correction term index
 } // end of for(Int sc=0;sc<2;sc++) // sin or cos correction terms 
                                                                                                                                                                                               
} // end of void AliFlowAnalysisWithQCumulants::FinalizeCorrectionTermsForNUAIntFlow()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::GetPointersForNestedLoopsHistograms()
{
 // Get pointers to all objects relevant for calculations with nested loops.
   
 TList *nestedLoopsList = dynamic_cast<TList*>(fHistList->FindObject("Nested Loops"));
 if(nestedLoopsList) 
 {
  this->SetNestedLoopsList(nestedLoopsList);
 } else
   {
    cout<<"WARNING: nestedLoopsList is NULL in AFAWQC::GPFNLH() !!!!"<<endl;
    exit(0);
   }
    
  TString sinCosFlag[2] = {"sin","cos"}; // to be improved (should I promote this to data members?)
  TString typeFlag[2] = {"RP","POI"}; // to be improved (should I promote this to data members?)
  TString ptEtaFlag[2] = {"p_{T}","#eta"}; // to be improved (should I promote this to data members?)
  TString reducedCorrelationIndex[4] = {"<2'>","<4'>","<6'>","<8'>"}; // to be improved (should I promote this to data members?)
   
  TString evaluateNestedLoopsName = "fEvaluateNestedLoops";
  evaluateNestedLoopsName += fAnalysisLabel->Data();  
  TProfile *evaluateNestedLoops = dynamic_cast<TProfile*>(nestedLoopsList->FindObject(evaluateNestedLoopsName.Data()));
  Bool_t bEvaluateIntFlowNestedLoops = kFALSE;
  Bool_t bEvaluateDiffFlowNestedLoops = kFALSE;
  if(evaluateNestedLoops)
  {
   this->SetEvaluateNestedLoops(evaluateNestedLoops);
   bEvaluateIntFlowNestedLoops = (Int_t)evaluateNestedLoops->GetBinContent(1);
   bEvaluateDiffFlowNestedLoops = (Int_t)evaluateNestedLoops->GetBinContent(2);
  }
  // nested loops relevant for integrated flow:  
  if(bEvaluateIntFlowNestedLoops)
  {
   // correlations:
   TString intFlowDirectCorrelationsName = "fIntFlowDirectCorrelations";
   intFlowDirectCorrelationsName += fAnalysisLabel->Data();
   TProfile *intFlowDirectCorrelations = dynamic_cast<TProfile*>(nestedLoopsList->FindObject(intFlowDirectCorrelationsName.Data()));
   if(intFlowDirectCorrelations) 
   { 
    this->SetIntFlowDirectCorrelations(intFlowDirectCorrelations);
   } else
     {
      cout<<"WARNING: intFlowDirectCorrelations is NULL in AFAWQC::GPFNLH() !!!!"<<endl;
      exit(0);
     }
   if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)  
   {
    TString intFlowExtraDirectCorrelationsName = "fIntFlowExtraDirectCorrelations";
    intFlowExtraDirectCorrelationsName += fAnalysisLabel->Data();
    TProfile *intFlowExtraDirectCorrelations = dynamic_cast<TProfile*>(nestedLoopsList->FindObject(intFlowExtraDirectCorrelationsName.Data()));
    if(intFlowExtraDirectCorrelations) 
    { 
     this->SetIntFlowExtraDirectCorrelations(intFlowExtraDirectCorrelations);
    } else
      {
       cout<<"WARNING: intFlowExtraDirectCorrelations is NULL in AFAWQC::GPFNLH() !!!!"<<endl;
       exit(0);
      }       
   } // end of if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)  
   // correction terms for non-uniform acceptance:
   TString intFlowDirectCorrectionTermsForNUAName = "fIntFlowDirectCorrectionTermsForNUA";
   intFlowDirectCorrectionTermsForNUAName += fAnalysisLabel->Data();
   TProfile *intFlowDirectCorrectionTermsForNUA[2] = {NULL};
   for(Int_t sc=0;sc<2;sc++) // sin or cos terms
   {
    intFlowDirectCorrectionTermsForNUA[sc] = dynamic_cast<TProfile*>(nestedLoopsList->FindObject(Form("%s: %s terms",intFlowDirectCorrectionTermsForNUAName.Data(),sinCosFlag[sc].Data())));
    if(intFlowDirectCorrectionTermsForNUA[sc]) 
    { 
     this->SetIntFlowDirectCorrectionTermsForNUA(intFlowDirectCorrectionTermsForNUA[sc],sc);
    } else
      {
       cout<<"WARNING: intFlowDirectCorrectionTermsForNUA[sc] is NULL in AFAWQC::GPFNLH() !!!!"<<endl;
       cout<<"sc = "<<sc<<endl;
       exit(0);
      }
   } // end of for(Int_t sc=0;sc<2;sc++) 
  } // end of if(bEvaluateIntFlowNestedLoops)
    
  // nested loops relevant for differential flow:  
  if(bEvaluateDiffFlowNestedLoops)
  {
   // correlations:
   TString diffFlowDirectCorrelationsName = "fDiffFlowDirectCorrelations";
   diffFlowDirectCorrelationsName += fAnalysisLabel->Data();
   TProfile *diffFlowDirectCorrelations[2][2][4] = {{{NULL}}};
   for(Int_t t=0;t<2;t++)
   {
    for(Int_t pe=0;pe<2;pe++)
    {
     for(Int_t ci=0;ci<4;ci++) // correlation index
     {
      diffFlowDirectCorrelations[t][pe][ci] = dynamic_cast<TProfile*>(nestedLoopsList->FindObject(Form("%s, %s, %s, %s",diffFlowDirectCorrelationsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[ci].Data())));
      if(diffFlowDirectCorrelations[t][pe][ci])
      {
       this->SetDiffFlowDirectCorrelations(diffFlowDirectCorrelations[t][pe][ci],t,pe,ci);
      } else
        {
         cout<<"WARNING: diffFlowDirectCorrelations[t][pe][ci] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
         cout<<"t  = "<<t<<endl;
         cout<<"pe = "<<pe<<endl;   
         cout<<"ci = "<<ci<<endl;
        }     
     } // end of for(Int_t ci=0;ci<4;ci++) // correlation index  
    } // end of for(Int_t pe=0;pe<2;pe++)
   } // end of for(Int_t t=0;t<2;t++)   
   // correction terms for non-uniform acceptance:
   TString diffFlowDirectCorrectionTermsForNUAName = "fDiffFlowDirectCorrectionTermsForNUA";
   diffFlowDirectCorrectionTermsForNUAName += fAnalysisLabel->Data();  
   TProfile *diffFlowDirectCorrectionTermsForNUA[2][2][2][10] = {{{{NULL}}}};   
   for(Int_t t=0;t<2;t++)
   {
    for(Int_t pe=0;pe<2;pe++)
    {
     // correction terms for NUA:
     for(Int_t sc=0;sc<2;sc++) // sin or cos
     {
      for(Int_t cti=0;cti<9;cti++) // correction term index
      {
       diffFlowDirectCorrectionTermsForNUA[t][pe][sc][cti] = dynamic_cast<TProfile*>(nestedLoopsList->FindObject(Form("%s, %s, %s, %s, cti = %d",diffFlowDirectCorrectionTermsForNUAName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),cti+1)));
       if(diffFlowDirectCorrectionTermsForNUA[t][pe][sc][cti])
       {
        this->SetDiffFlowDirectCorrectionTermsForNUA(diffFlowDirectCorrectionTermsForNUA[t][pe][sc][cti],t,pe,sc,cti);
       } else
         {
          cout<<"WARNING: diffFlowDirectCorrectionTermsForNUA[t][pe][sc][cti] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
          cout<<"t   = "<<t<<endl;
          cout<<"pe  = "<<pe<<endl;   
          cout<<"sc  = "<<sc<<endl;
          cout<<"cti = "<<cti<<endl;
         }    
      } // end of for(Int_t cti=0;cti<9;cti++) // correction term index
     } // end of for(Int_t sc=0;sc<2;sc++) // sin or cos
    } // end of for(Int_t pe=0;pe<2;pe++)
   } // end of for(Int_t t=0;t<2;t++)
   // other differential correlators:
   TString otherDirectDiffCorrelatorsName = "fOtherDirectDiffCorrelators";
   otherDirectDiffCorrelatorsName += fAnalysisLabel->Data();  
   TProfile *otherDirectDiffCorrelators[2][2][2][1] = {{{{NULL}}}};   
   for(Int_t t=0;t<2;t++)
   {
    for(Int_t pe=0;pe<2;pe++)
    {
     // correction terms for NUA:
     for(Int_t sc=0;sc<2;sc++) // sin or cos
     {
      for(Int_t ci=0;ci<1;ci++) // correlator index
      {
       otherDirectDiffCorrelators[t][pe][sc][ci] = dynamic_cast<TProfile*>(nestedLoopsList->FindObject(Form("%s, %s, %s, %s, ci = %d",otherDirectDiffCorrelatorsName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),sinCosFlag[sc].Data(),ci+1)));
       if(otherDirectDiffCorrelators[t][pe][sc][ci])
       {
        this->SetOtherDirectDiffCorrelators(otherDirectDiffCorrelators[t][pe][sc][ci],t,pe,sc,ci);
       } else
         {
          cout<<"WARNING: otherDirectDiffCorrelators[t][pe][sc][ci] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
          cout<<"t   = "<<t<<endl;
          cout<<"pe  = "<<pe<<endl;   
          cout<<"sc  = "<<sc<<endl;
          cout<<"ci = "<<ci<<endl;
         }    
      } // end of for(Int_t ci=0;ci<9;ci++) // correction term index
     } // end of for(Int_t sc=0;sc<2;sc++) // sin or cos
    } // end of for(Int_t pe=0;pe<2;pe++)
   } // end of for(Int_t t=0;t<2;t++)
   // number of RPs and POIs in selected pt and eta bins for cross-checkings:
   TString noOfParticlesInBinName = "fNoOfParticlesInBin";
   TH1D *noOfParticlesInBin = NULL;
   noOfParticlesInBin = dynamic_cast<TH1D*>(nestedLoopsList->FindObject(noOfParticlesInBinName.Data()));
   if(noOfParticlesInBin)
   {
    this->SetNoOfParticlesInBin(noOfParticlesInBin);
   } else
     {
      cout<<endl;
      cout<<" WARNING (QC): noOfParticlesInBin is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
      cout<<endl;
     }
  } // end of if(bEvaluateDiffFlowNestedLoops)

} // end of void AliFlowAnalysisWithQCumulants::GetPointersForNestedLoopsHistograms()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::StoreHarmonic()
{
 // Store flow harmonic in common control histograms.

 (fCommonHists->GetHarmonic())->Fill(0.5,fHarmonic);
 if(fFillMultipleControlHistograms)
 {
  (fCommonHists2nd->GetHarmonic())->Fill(0.5,fHarmonic);
  (fCommonHists4th->GetHarmonic())->Fill(0.5,fHarmonic);
  (fCommonHists6th->GetHarmonic())->Fill(0.5,fHarmonic);
  (fCommonHists8th->GetHarmonic())->Fill(0.5,fHarmonic);
 }
 
} // end of void AliFlowAnalysisWithQCumulants::StoreHarmonic()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrelationsUsingParticleWeights(TString type, TString ptOrEta) // type = RP or POI 
{
 // Calculate all correlations needed for differential flow using particle weights.
 
 Int_t t = 0; // type flag 
 Int_t pe = 0; // ptEta flag
 
 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }

 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   }
    
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};

 // real and imaginary parts of weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n1k = (*fReQ)(0,1);
 Double_t dReQ2n2k = (*fReQ)(1,2);
 Double_t dReQ1n3k = (*fReQ)(0,3);
 //Double_t dReQ4n4k = (*fReQ)(3,4);
 Double_t dImQ1n1k = (*fImQ)(0,1);
 Double_t dImQ2n2k = (*fImQ)(1,2);
 Double_t dImQ1n3k = (*fImQ)(0,3);
 //Double_t dImQ4n4k = (*fImQ)(3,4);
 
 // S^M_{p,k} (see .h file for the definition of fSpk):
 Double_t dSM1p1k = (*fSpk)(0,1);
 Double_t dSM1p2k = (*fSpk)(0,2);
 Double_t dSM1p3k = (*fSpk)(0,3);
 Double_t dSM2p1k = (*fSpk)(1,1);
 Double_t dSM3p1k = (*fSpk)(2,1);
 
 // looping over all bins and calculating reduced correlations: 
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // real and imaginary parts of p_{m*n,0} (non-weighted Q-vector evaluated for POIs in particular (pt,eta) bin):  
  Double_t p1n0kRe = 0.;
  Double_t p1n0kIm = 0.;

  // number of POIs in particular (pt,eta) bin):
  Double_t mp = 0.;

  // real and imaginary parts of q_{m*n,k}: 
  // (weighted Q-vector evaluated for particles which are both RPs and POIs in particular (pt,eta) bin)
  Double_t q1n2kRe = 0.;
  Double_t q1n2kIm = 0.;
  Double_t q2n1kRe = 0.;
  Double_t q2n1kIm = 0.;

  // s_{1,1}, s_{1,2} and s_{1,3} // to be improved (add explanation)  
  Double_t s1p1k = 0.; 
  Double_t s1p2k = 0.; 
  Double_t s1p3k = 0.; 
   
  // M0111 from Eq. (118) in QC2c (to be improved (notation))
  Double_t dM0111 = 0.;
 
  if(type == "POI")
  {
   p1n0kRe = fReRPQ1dEBE[1][pe][0][0]->GetBinContent(fReRPQ1dEBE[1][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b));
   p1n0kIm = fImRPQ1dEBE[1][pe][0][0]->GetBinContent(fImRPQ1dEBE[1][pe][0][0]->GetBin(b))  
           * fImRPQ1dEBE[1][pe][0][0]->GetBinEntries(fImRPQ1dEBE[1][pe][0][0]->GetBin(b));
            
   mp = fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
    
   t = 1; // typeFlag = RP or POI
    
   // q_{m*n,k}: (Remark: m=1 is 0, k=0 iz zero (to be improved!)) 
   q1n2kRe = fReRPQ1dEBE[2][pe][0][2]->GetBinContent(fReRPQ1dEBE[2][pe][0][2]->GetBin(b))
           * fReRPQ1dEBE[2][pe][0][2]->GetBinEntries(fReRPQ1dEBE[2][pe][0][2]->GetBin(b));
   q1n2kIm = fImRPQ1dEBE[2][pe][0][2]->GetBinContent(fImRPQ1dEBE[2][pe][0][2]->GetBin(b))
           * fImRPQ1dEBE[2][pe][0][2]->GetBinEntries(fImRPQ1dEBE[2][pe][0][2]->GetBin(b));
   q2n1kRe = fReRPQ1dEBE[2][pe][1][1]->GetBinContent(fReRPQ1dEBE[2][pe][1][1]->GetBin(b))
           * fReRPQ1dEBE[2][pe][1][1]->GetBinEntries(fReRPQ1dEBE[2][pe][1][1]->GetBin(b));
   q2n1kIm = fImRPQ1dEBE[2][pe][1][1]->GetBinContent(fImRPQ1dEBE[2][pe][1][1]->GetBin(b))
           * fImRPQ1dEBE[2][pe][1][1]->GetBinEntries(fImRPQ1dEBE[2][pe][1][1]->GetBin(b));
       
   // s_{1,1}, s_{1,2} and s_{1,3} // to be improved (add explanation)  
   s1p1k = pow(fs1dEBE[2][pe][1]->GetBinContent(b)*fs1dEBE[2][pe][1]->GetBinEntries(b),1.); 
   s1p2k = pow(fs1dEBE[2][pe][2]->GetBinContent(b)*fs1dEBE[2][pe][2]->GetBinEntries(b),1.); 
   s1p3k = pow(fs1dEBE[2][pe][3]->GetBinContent(b)*fs1dEBE[2][pe][3]->GetBinEntries(b),1.); 
     
   // M0111 from Eq. (118) in QC2c (to be improved (notation)):
   dM0111 = mp*(dSM3p1k-3.*dSM1p1k*dSM1p2k+2.*dSM1p3k)
          - 3.*(s1p1k*(dSM2p1k-dSM1p2k)
          + 2.*(s1p3k-s1p2k*dSM1p1k));
  }
   else if(type == "RP")
   {
    // q_{m*n,k}: (Remark: m=1 is 0, k=0 iz zero (to be improved!)) 
    q1n2kRe = fReRPQ1dEBE[0][pe][0][2]->GetBinContent(fReRPQ1dEBE[0][pe][0][2]->GetBin(b))
            * fReRPQ1dEBE[0][pe][0][2]->GetBinEntries(fReRPQ1dEBE[0][pe][0][2]->GetBin(b));
    q1n2kIm = fImRPQ1dEBE[0][pe][0][2]->GetBinContent(fImRPQ1dEBE[0][pe][0][2]->GetBin(b))
            * fImRPQ1dEBE[0][pe][0][2]->GetBinEntries(fImRPQ1dEBE[0][pe][0][2]->GetBin(b));
    q2n1kRe = fReRPQ1dEBE[0][pe][1][1]->GetBinContent(fReRPQ1dEBE[0][pe][1][1]->GetBin(b))
            * fReRPQ1dEBE[0][pe][1][1]->GetBinEntries(fReRPQ1dEBE[0][pe][1][1]->GetBin(b));
    q2n1kIm = fImRPQ1dEBE[0][pe][1][1]->GetBinContent(fImRPQ1dEBE[0][pe][1][1]->GetBin(b))
            * fImRPQ1dEBE[0][pe][1][1]->GetBinEntries(fImRPQ1dEBE[0][pe][1][1]->GetBin(b));

    // s_{1,1}, s_{1,2} and s_{1,3} // to be improved (add explanation)  
    s1p1k = pow(fs1dEBE[0][pe][1]->GetBinContent(b)*fs1dEBE[0][pe][1]->GetBinEntries(b),1.); 
    s1p2k = pow(fs1dEBE[0][pe][2]->GetBinContent(b)*fs1dEBE[0][pe][2]->GetBinEntries(b),1.); 
    s1p3k = pow(fs1dEBE[0][pe][3]->GetBinContent(b)*fs1dEBE[0][pe][3]->GetBinEntries(b),1.); 
    
    // to be improved (cross-checked):
    p1n0kRe = fReRPQ1dEBE[0][pe][0][0]->GetBinContent(fReRPQ1dEBE[0][pe][0][0]->GetBin(b))
            * fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b));
    p1n0kIm = fImRPQ1dEBE[0][pe][0][0]->GetBinContent(fImRPQ1dEBE[0][pe][0][0]->GetBin(b))  
            * fImRPQ1dEBE[0][pe][0][0]->GetBinEntries(fImRPQ1dEBE[0][pe][0][0]->GetBin(b));
            
    mp = fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
     
    t = 0; // typeFlag = RP or POI
    
    // M0111 from Eq. (118) in QC2c (to be improved (notation)):
    dM0111 = mp*(dSM3p1k-3.*dSM1p1k*dSM1p2k+2.*dSM1p3k)
           - 3.*(s1p1k*(dSM2p1k-dSM1p2k)
           + 2.*(s1p3k-s1p2k*dSM1p1k));
    //...............................................................................................   
   }
   
   // 2'-particle correlation:
   Double_t two1n1nW0W1 = 0.;
   if(mp*dSM1p1k-s1p1k)
   {
    two1n1nW0W1 = (p1n0kRe*dReQ1n1k+p1n0kIm*dImQ1n1k-s1p1k)
                / (mp*dSM1p1k-s1p1k);
   
    // fill profile to get <<2'>>     
    fDiffFlowCorrelationsPro[t][pe][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],two1n1nW0W1,mp*dSM1p1k-s1p1k);    
    // fill profile to get <<2'>^2>     
    fDiffFlowSquaredCorrelationsPro[t][pe][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],two1n1nW0W1*two1n1nW0W1,mp*dSM1p1k-s1p1k);        
    // histogram to store <2'> e-b-e (needed in some other methods):
    fDiffFlowCorrelationsEBE[t][pe][0]->SetBinContent(b,two1n1nW0W1);      
    fDiffFlowEventWeightsForCorrelationsEBE[t][pe][0]->SetBinContent(b,mp*dSM1p1k-s1p1k);      
   } // end of if(mp*dSM1p1k-s1p1k)
   
   // 4'-particle correlation:
   Double_t four1n1n1n1nW0W1W1W1 = 0.;
   if(dM0111)
   {
    four1n1n1n1nW0W1W1W1 = ((pow(dReQ1n1k,2.)+pow(dImQ1n1k,2.))*(p1n0kRe*dReQ1n1k+p1n0kIm*dImQ1n1k)
                         - q2n1kRe*(pow(dReQ1n1k,2.)-pow(dImQ1n1k,2.))
                         - 2.*q2n1kIm*dReQ1n1k*dImQ1n1k
                         - p1n0kRe*(dReQ1n1k*dReQ2n2k+dImQ1n1k*dImQ2n2k)
                         + p1n0kIm*(dImQ1n1k*dReQ2n2k-dReQ1n1k*dImQ2n2k)
                         - 2.*dSM1p2k*(p1n0kRe*dReQ1n1k+p1n0kIm*dImQ1n1k)
                         - 2.*(pow(dReQ1n1k,2.)+pow(dImQ1n1k,2.))*s1p1k                                            
                         + 6.*(q1n2kRe*dReQ1n1k+q1n2kIm*dImQ1n1k)                                           
                         + 1.*(q2n1kRe*dReQ2n2k+q2n1kIm*dImQ2n2k)                         
                         + 2.*(p1n0kRe*dReQ1n3k+p1n0kIm*dImQ1n3k)                      
                         + 2.*s1p1k*dSM1p2k                                      
                         - 6.*s1p3k)        
                         / dM0111; // to be improved (notation of dM0111)
   
    // fill profile to get <<4'>>     
    fDiffFlowCorrelationsPro[t][pe][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],four1n1n1n1nW0W1W1W1,dM0111);    
    // fill profile to get <<4'>^2>     
    fDiffFlowSquaredCorrelationsPro[t][pe][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],four1n1n1n1nW0W1W1W1*four1n1n1n1nW0W1W1W1,dM0111);        
    // histogram to store <4'> e-b-e (needed in some other methods):
    fDiffFlowCorrelationsEBE[t][pe][1]->SetBinContent(b,four1n1n1n1nW0W1W1W1);      
    fDiffFlowEventWeightsForCorrelationsEBE[t][pe][1]->SetBinContent(b,dM0111);      
   } // end of if(dM0111)
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)

} // end of void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrelationsUsingParticleWeights(TString type, TString ptOrEta); // type = RP or POI 

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::FillCommonControlHistograms(AliFlowEventSimple *anEvent)
{
 // Fill common control histograms.
 
 Int_t nRP = anEvent->GetEventNSelTracksRP(); // number of RPs (i.e. number of particles used to determine the reaction plane)
 fCommonHists->FillControlHistograms(anEvent); 
 if(fFillMultipleControlHistograms)
 {
  if(nRP>1)
  {
   fCommonHists2nd->FillControlHistograms(anEvent);                                        
   if(nRP>3)
   {
    fCommonHists4th->FillControlHistograms(anEvent);                                        
    if(nRP>5)
    {
     fCommonHists6th->FillControlHistograms(anEvent);                                        
     if(nRP>7)
     {
      fCommonHists8th->FillControlHistograms(anEvent);                                        
     } // end of if(nRP>7)  
    } // end of if(nRP>5) 
   } // end of if(nRP>3)                                                                                                                      
  } // end of if(nRP>1) 
 } // end of if(fFillMultipleControlHistograms)
 
} // end of void AliFlowAnalysisWithQCumulants::FillCommonControlHistograms(AliFlowEventSimple *anEvent)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::ResetEventByEventQuantities()
{
 // Reset all event by event quantities.
 
 // Reference flow:
 fReQ->Zero();
 fImQ->Zero();
 fSpk->Zero();
 fIntFlowCorrelationsEBE->Reset();
 fIntFlowEventWeightsForCorrelationsEBE->Reset();
 fIntFlowCorrelationsAllEBE->Reset();
 
 for(Int_t sc=0;sc<2;sc++)
 {
  fIntFlowCorrectionTermsForNUAEBE[sc]->Reset();
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[sc]->Reset(); 
 }
    
 // Differential flow:
 if(fCalculateDiffFlow)
 {
  for(Int_t t=0;t<3;t++) // type (RP, POI, POI&&RP)
  {
   for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // 1D in pt or eta
   {
    for(Int_t m=0;m<4;m++) // multiple of harmonic
    {
     for(Int_t k=0;k<9;k++) // power of weight
     {
      if(fReRPQ1dEBE[t][pe][m][k]) fReRPQ1dEBE[t][pe][m][k]->Reset();
      if(fImRPQ1dEBE[t][pe][m][k]) fImRPQ1dEBE[t][pe][m][k]->Reset();
     }   
    } 
   }
  } 
  for(Int_t t=0;t<3;t++) // type (0 = RP, 1 = POI, 2 = RP&&POI )
  { 
   for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // 1D in pt or eta
   {
    for(Int_t k=0;k<9;k++)
    {
     if(fs1dEBE[t][pe][k]) fs1dEBE[t][pe][k]->Reset();
    }
   }
  }
  // e-b-e reduced correlations:
  for(Int_t t=0;t<2;t++) // type (0 = RP, 1 = POI)
  {  
   for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
   {
    for(Int_t rci=0;rci<4;rci++) // reduced correlation index
    {
     if(fDiffFlowCorrelationsEBE[t][pe][rci]) fDiffFlowCorrelationsEBE[t][pe][rci]->Reset();
     if(fDiffFlowEventWeightsForCorrelationsEBE[t][pe][rci]) fDiffFlowEventWeightsForCorrelationsEBE[t][pe][rci]->Reset();
    }
   }
  }  
  // correction terms for NUA:
  for(Int_t t=0;t<2;t++) // type (0 = RP, 1 = POI)
  {  
   for(Int_t pe=0;pe<1+(Int_t)fCalculateDiffFlowVsEta;pe++) // pt or eta
   {
    for(Int_t sc=0;sc<2;sc++) // sin or cos
    {
     for(Int_t cti=0;cti<9;cti++) // correction term index
     {
     fDiffFlowCorrectionTermsForNUAEBE[t][pe][sc][cti]->Reset();  
     }
    }
   }      
  }
 } // end of if(fCalculateDiffFlow)   

 // 2D (pt,eta)
 if(fCalculate2DDiffFlow)
 {
  for(Int_t t=0;t<3;t++) // type (RP, POI, POI&&RP)
  {
   for(Int_t m=0;m<4;m++) // multiple of harmonic
   {
    for(Int_t k=0;k<9;k++) // power of weight
    {
     if(fReRPQ2dEBE[t][m][k]){fReRPQ2dEBE[t][m][k]->Reset();}
     if(fImRPQ2dEBE[t][m][k]){fImRPQ2dEBE[t][m][k]->Reset();}
    }   
   }
  }
  for(Int_t t=0;t<3;t++) // type (0 = RP, 1 = POI, 2 = RP&&POI )
  { 
   for(Int_t k=0;k<9;k++)
   {
    if(fs2dEBE[t][k]){fs2dEBE[t][k]->Reset();}
   }
  }  
 } // end of if(fCalculate2DDiffFlow) 

} // end of void AliFlowAnalysisWithQCumulants::ResetEventByEventQuantities();

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectionsForNUASinTerms(TString type, TString ptOrEta)
{
 // Calculate correction terms for non-uniform acceptance for differential flow (sin terms).
 
 // Results are stored in fDiffFlowCorrectionTermsForNUAPro[t][pe][0][cti], where cti runs as follows:
 //  0: <<sin n(psi1)>>
 //  1: <<sin n(psi1+phi2)>>
 //  2: <<sin n(psi1+phi2-phi3)>>
 //  3: <<sin n(psi1-phi2-phi3)>>:
 //  4:
 //  5:
 //  6:
 
 // multiplicity:
 Double_t dMult = (*fSpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 //Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 //Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);

 Int_t t = 0; // type flag 
 Int_t pe = 0; // ptEta flag
 
 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }

 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   }
    
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};

 // looping over all bins and calculating correction terms: 
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // real and imaginary parts of p_{m*n,0} (non-weighted Q-vector evaluated for POIs in particular pt or eta bin): 
  Double_t p1n0kRe = 0.;
  Double_t p1n0kIm = 0.;

  // number of POIs in particular pt or eta bin:
  Double_t mp = 0.;

  // real and imaginary parts of q_{m*n,0} (non-weighted Q-vector evaluated for particles which are both RPs and POIs in particular pt or eta bin):
  Double_t q1n0kRe = 0.;
  Double_t q1n0kIm = 0.;
  Double_t q2n0kRe = 0.;
  Double_t q2n0kIm = 0.;

  // number of particles which are both RPs and POIs in particular pt or eta bin:
  Double_t mq = 0.;
   
  if(type == "POI")
  {
   // q_{m*n,0}:
   q1n0kRe = fReRPQ1dEBE[2][pe][0][0]->GetBinContent(fReRPQ1dEBE[2][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(fReRPQ1dEBE[2][pe][0][0]->GetBin(b));
   q1n0kIm = fImRPQ1dEBE[2][pe][0][0]->GetBinContent(fImRPQ1dEBE[2][pe][0][0]->GetBin(b))
           * fImRPQ1dEBE[2][pe][0][0]->GetBinEntries(fImRPQ1dEBE[2][pe][0][0]->GetBin(b));
   q2n0kRe = fReRPQ1dEBE[2][pe][1][0]->GetBinContent(fReRPQ1dEBE[2][pe][1][0]->GetBin(b))
           * fReRPQ1dEBE[2][pe][1][0]->GetBinEntries(fReRPQ1dEBE[2][pe][1][0]->GetBin(b));
   q2n0kIm = fImRPQ1dEBE[2][pe][1][0]->GetBinContent(fImRPQ1dEBE[2][pe][1][0]->GetBin(b))
           * fImRPQ1dEBE[2][pe][1][0]->GetBinEntries(fImRPQ1dEBE[2][pe][1][0]->GetBin(b));         
                 
   mq = fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(fReRPQ1dEBE[2][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
  } 
  else if(type == "RP")
  {
   // q_{m*n,0}:
   q1n0kRe = fReRPQ1dEBE[0][pe][0][0]->GetBinContent(fReRPQ1dEBE[0][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b));
   q1n0kIm = fImRPQ1dEBE[0][pe][0][0]->GetBinContent(fImRPQ1dEBE[0][pe][0][0]->GetBin(b))
           * fImRPQ1dEBE[0][pe][0][0]->GetBinEntries(fImRPQ1dEBE[0][pe][0][0]->GetBin(b));
   q2n0kRe = fReRPQ1dEBE[0][pe][1][0]->GetBinContent(fReRPQ1dEBE[0][pe][1][0]->GetBin(b))
           * fReRPQ1dEBE[0][pe][1][0]->GetBinEntries(fReRPQ1dEBE[0][pe][1][0]->GetBin(b));
   q2n0kIm = fImRPQ1dEBE[0][pe][1][0]->GetBinContent(fImRPQ1dEBE[0][pe][1][0]->GetBin(b))
           * fImRPQ1dEBE[0][pe][1][0]->GetBinEntries(fImRPQ1dEBE[0][pe][1][0]->GetBin(b));         
                 
   mq = fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)  
  }    
  if(type == "POI")
  {
   // p_{m*n,0}:
   p1n0kRe = fReRPQ1dEBE[1][pe][0][0]->GetBinContent(fReRPQ1dEBE[1][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b));
   p1n0kIm = fImRPQ1dEBE[1][pe][0][0]->GetBinContent(fImRPQ1dEBE[1][pe][0][0]->GetBin(b))  
           * fImRPQ1dEBE[1][pe][0][0]->GetBinEntries(fImRPQ1dEBE[1][pe][0][0]->GetBin(b));
            
   mp = fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
    
   t = 1; // typeFlag = RP or POI
  }
  else if(type == "RP")
  {
   // p_{m*n,0} = q_{m*n,0}:
   p1n0kRe = q1n0kRe; 
   p1n0kIm = q1n0kIm; 
           
   mp = mq; 
   
   t = 0; // typeFlag = RP or POI
  }

  // <<sin n(psi1)>>:
  Double_t sinP1nPsi = 0.;
  if(mp)
  {
   sinP1nPsi = p1n0kIm/mp;
   // fill profile for <<sin n(psi1)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][0][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sinP1nPsi,mp);
   // histogram to store <sin n(psi1)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][0][0]->SetBinContent(b,sinP1nPsi);
  } // end of if(mp)   
  
  // <<sin n(psi1+phi2)>>:
  Double_t sinP1nPsiP1nPhi = 0.;
  if(mp*dMult-mq)
  {
   sinP1nPsiP1nPhi = (p1n0kRe*dImQ1n+p1n0kIm*dReQ1n-q2n0kIm)/(mp*dMult-mq);
   // fill profile for <<sin n(psi1+phi2)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][0][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sinP1nPsiP1nPhi,mp*dMult-mq);
   // histogram to store <sin n(psi1+phi2)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][0][1]->SetBinContent(b,sinP1nPsiP1nPhi);
  } // end of if(mp*dMult-mq)   
  
  // <<sin n(psi1+phi2-phi3)>>:
  Double_t sinP1nPsi1P1nPhi2MPhi3 = 0.;
  if(mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.))
  {
   sinP1nPsi1P1nPhi2MPhi3 = (p1n0kIm*(pow(dImQ1n,2.)+pow(dReQ1n,2.)-dMult)
                          - 1.*(q2n0kIm*dReQ1n-q2n0kRe*dImQ1n)  
                          - mq*dImQ1n+2.*q1n0kIm)
                          / (mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.));
   // fill profile for <<sin n(psi1+phi2)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][0][2]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sinP1nPsi1P1nPhi2MPhi3,mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.));
   // histogram to store <sin n(psi1+phi2)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][0][2]->SetBinContent(b,sinP1nPsi1P1nPhi2MPhi3);
  } // end of if(mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.))   
  
  // <<sin n(psi1-phi2-phi3)>>:
  Double_t sinP1nPsi1M1nPhi2MPhi3 = 0.;
  if(mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.))
  {
   sinP1nPsi1M1nPhi2MPhi3 = (p1n0kIm*(pow(dReQ1n,2.)-pow(dImQ1n,2.))-2.*p1n0kRe*dReQ1n*dImQ1n
                          - 1.*(p1n0kIm*dReQ2n-p1n0kRe*dImQ2n)
                          + 2.*mq*dImQ1n-2.*q1n0kIm)
                          / (mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.));
   // fill profile for <<sin n(psi1+phi2)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][0][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sinP1nPsi1M1nPhi2MPhi3,mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.));
   // histogram to store <sin n(psi1+phi2)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][0][3]->SetBinContent(b,sinP1nPsi1M1nPhi2MPhi3);
  } // end of if(mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.))   
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectionsForNUASinTerms(TString type, TString ptOrEta)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectionsForNUACosTerms(TString type, TString ptOrEta)
{
 // Calculate correction terms for non-uniform acceptance for differential flow (cos terms).
 
 // Results are stored in fDiffFlowCorrectionTermsForNUAPro[t][pe][1][cti], where cti runs as follows:
 //  0: <<cos n(psi)>>
 //  1: <<cos n(psi1+phi2)>>
 //  2: <<cos n(psi1+phi2-phi3)>>
 //  3: <<cos n(psi1-phi2-phi3)>>
 //  4:
 //  5:
 //  6:
 
 // multiplicity:
 Double_t dMult = (*fSpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 //Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 //Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);

 Int_t t = 0; // type flag 
 Int_t pe = 0; // ptEta flag
 
 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }

 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   }
    
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};

 // looping over all bins and calculating correction terms: 
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // real and imaginary parts of p_{m*n,0} (non-weighted Q-vector evaluated for POIs in particular pt or eta bin): 
  Double_t p1n0kRe = 0.;
  Double_t p1n0kIm = 0.;

  // number of POIs in particular pt or eta bin:
  Double_t mp = 0.;

  // real and imaginary parts of q_{m*n,0} (non-weighted Q-vector evaluated for particles which are both RPs and POIs in particular pt or eta bin):
  Double_t q1n0kRe = 0.;
  Double_t q1n0kIm = 0.;
  Double_t q2n0kRe = 0.;
  Double_t q2n0kIm = 0.;

  // number of particles which are both RPs and POIs in particular pt or eta bin:
  Double_t mq = 0.;
   
  if(type == "POI")
  {
   // q_{m*n,0}:
   q1n0kRe = fReRPQ1dEBE[2][pe][0][0]->GetBinContent(fReRPQ1dEBE[2][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(fReRPQ1dEBE[2][pe][0][0]->GetBin(b));
   q1n0kIm = fImRPQ1dEBE[2][pe][0][0]->GetBinContent(fImRPQ1dEBE[2][pe][0][0]->GetBin(b))
           * fImRPQ1dEBE[2][pe][0][0]->GetBinEntries(fImRPQ1dEBE[2][pe][0][0]->GetBin(b));
   q2n0kRe = fReRPQ1dEBE[2][pe][1][0]->GetBinContent(fReRPQ1dEBE[2][pe][1][0]->GetBin(b))
           * fReRPQ1dEBE[2][pe][1][0]->GetBinEntries(fReRPQ1dEBE[2][pe][1][0]->GetBin(b));
   q2n0kIm = fImRPQ1dEBE[2][pe][1][0]->GetBinContent(fImRPQ1dEBE[2][pe][1][0]->GetBin(b))
           * fImRPQ1dEBE[2][pe][1][0]->GetBinEntries(fImRPQ1dEBE[2][pe][1][0]->GetBin(b));         
                 
   mq = fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(fReRPQ1dEBE[2][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
  } 
  else if(type == "RP")
  {
   // q_{m*n,0}:
   q1n0kRe = fReRPQ1dEBE[0][pe][0][0]->GetBinContent(fReRPQ1dEBE[0][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b));
   q1n0kIm = fImRPQ1dEBE[0][pe][0][0]->GetBinContent(fImRPQ1dEBE[0][pe][0][0]->GetBin(b))
           * fImRPQ1dEBE[0][pe][0][0]->GetBinEntries(fImRPQ1dEBE[0][pe][0][0]->GetBin(b));
   q2n0kRe = fReRPQ1dEBE[0][pe][1][0]->GetBinContent(fReRPQ1dEBE[0][pe][1][0]->GetBin(b))
           * fReRPQ1dEBE[0][pe][1][0]->GetBinEntries(fReRPQ1dEBE[0][pe][1][0]->GetBin(b));
   q2n0kIm = fImRPQ1dEBE[0][pe][1][0]->GetBinContent(fImRPQ1dEBE[0][pe][1][0]->GetBin(b))
           * fImRPQ1dEBE[0][pe][1][0]->GetBinEntries(fImRPQ1dEBE[0][pe][1][0]->GetBin(b));         
                 
   mq = fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)  
  }    
  if(type == "POI")
  {
   // p_{m*n,0}:
   p1n0kRe = fReRPQ1dEBE[1][pe][0][0]->GetBinContent(fReRPQ1dEBE[1][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b));
   p1n0kIm = fImRPQ1dEBE[1][pe][0][0]->GetBinContent(fImRPQ1dEBE[1][pe][0][0]->GetBin(b))  
           * fImRPQ1dEBE[1][pe][0][0]->GetBinEntries(fImRPQ1dEBE[1][pe][0][0]->GetBin(b));
            
   mp = fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
    
   t = 1; // typeFlag = RP or POI
  }
  else if(type == "RP")
  {
   // p_{m*n,0} = q_{m*n,0}:
   p1n0kRe = q1n0kRe; 
   p1n0kIm = q1n0kIm; 
           
   mp = mq; 
   
   t = 0; // typeFlag = RP or POI
  }

  // <<cos n(psi1)>>:
  Double_t cosP1nPsi = 0.;
  if(mp)
  {
   cosP1nPsi = p1n0kRe/mp;
   
   // fill profile for <<cos n(psi1)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][1][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],cosP1nPsi,mp);
   // histogram to store <cos n(psi1)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][1][0]->SetBinContent(b,cosP1nPsi);
  } // end of if(mp)   
  
  // <<cos n(psi1+phi2)>>:
  Double_t cosP1nPsiP1nPhi = 0.;
  if(mp*dMult-mq)
  {
   cosP1nPsiP1nPhi = (p1n0kRe*dReQ1n-p1n0kIm*dImQ1n-q2n0kRe)/(mp*dMult-mq);
   // fill profile for <<sin n(psi1+phi2)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][1][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],cosP1nPsiP1nPhi,mp*dMult-mq);
   // histogram to store <sin n(psi1+phi2)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][1][1]->SetBinContent(b,cosP1nPsiP1nPhi);
  } // end of if(mp*dMult-mq)   
  
  // <<cos n(psi1+phi2-phi3)>>:
  Double_t cosP1nPsi1P1nPhi2MPhi3 = 0.;
  if(mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.))
  {
   cosP1nPsi1P1nPhi2MPhi3 = (p1n0kRe*(pow(dImQ1n,2.)+pow(dReQ1n,2.)-dMult)
                          - 1.*(q2n0kRe*dReQ1n+q2n0kIm*dImQ1n)  
                          - mq*dReQ1n+2.*q1n0kRe)
                          / (mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.));
   // fill profile for <<sin n(psi1+phi2)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][1][2]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],cosP1nPsi1P1nPhi2MPhi3,mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.));
   // histogram to store <sin n(psi1+phi2)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][1][2]->SetBinContent(b,cosP1nPsi1P1nPhi2MPhi3);
  } // end of if(mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.))   
  
  // <<cos n(psi1-phi2-phi3)>>:
  Double_t cosP1nPsi1M1nPhi2MPhi3 = 0.;
  if(mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.))
  {
   cosP1nPsi1M1nPhi2MPhi3 = (p1n0kRe*(pow(dReQ1n,2.)-pow(dImQ1n,2.))+2.*p1n0kIm*dReQ1n*dImQ1n
                          - 1.*(p1n0kRe*dReQ2n+p1n0kIm*dImQ2n)  
                          - 2.*mq*dReQ1n+2.*q1n0kRe)
                          / (mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.));
   // fill profile for <<sin n(psi1+phi2)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][1][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],cosP1nPsi1M1nPhi2MPhi3,mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.));
   // histogram to store <sin n(psi1+phi2)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][1][3]->SetBinContent(b,cosP1nPsi1M1nPhi2MPhi3);
  } // end of if(mq*(dMult-1.)*(dMult-2.)+(mp-mq)*dMult*(dMult-1.))   
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectionsForNUACosTerms(TString type, TString ptOrEta)

//==================================================================================================================================

void AliFlowAnalysisWithQCumulants::FinalizeCorrectionTermsForNUADiffFlow(TString type, TString ptOrEta)
{
 // Transfer profiles into histogams and correctly propagate the error.
 
 Int_t t = 0; // type flag 
 Int_t pe = 0; // ptEta flag
 
 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }

 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   }
    
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 //Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 //Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};

 for(Int_t sc=0;sc<2;sc++) // sin or cos
 {
  for(Int_t cti=0;cti<9;cti++) // correction term index
  {
   for(Int_t b=1;b<=nBinsPtEta[pe];b++)
   {
    Double_t correctionTerm = fDiffFlowCorrectionTermsForNUAPro[t][pe][sc][cti]->GetBinContent(b);
    fDiffFlowCorrectionTermsForNUAHist[t][pe][sc][cti]->SetBinContent(b,correctionTerm);
    // to be improved (propagate error correctly)
    // ...
   } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)
  } // correction term index
 } // end of for(Int_t sc=0;sc<2;sc++) // sin or cos

}// end of void AliFlowAnalysisWithQCumulants::FinalizeCorrectionTermsForNUADiffFlow(TString type, TString ptOrEta)

//==================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCumulantsCorrectedForNUA(TString type, TString ptOrEta)
{ 
 // Calculate generalized differential flow cumulants (corrected for non-uniform acceptance).
 
 // to be improved - propagate error also from non-isotropic terms
  
 Int_t t = 0; // RP = 0, POI = 1
 Int_t pe = 0; // pt = 0, eta = 1

 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   } 
     
 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   } 
       
 // Common:
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 // 2-particle correlation:
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>>
 // sinus terms coming from reference flow: 
 Double_t sinP1nPhi = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(1); // <<sin(n*phi1)>>
 Double_t sinP1nPhi1P1nPhi2 = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(2); // <<sin(n*(phi1+phi2))>>
 Double_t sinP1nPhi1M1nPhi2M1nPhi3 = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(3); // <<sin(n*(phi1-phi2-phi3))>>
 // cosinus terms coming from reference flow: 
 Double_t cosP1nPhi = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(1); // <<cos(n*phi1)>>
 Double_t cosP1nPhi1P1nPhi2 = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(2); // <<cos(n*(phi1+phi2))>>
 Double_t cosP1nPhi1M1nPhi2M1nPhi3 = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(3); // <<cos(n*(phi1-phi2-phi3))>>

 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  Double_t twoPrime = fDiffFlowCorrelationsHist[t][pe][0]->GetBinContent(b); // <<2'>>
  Double_t fourPrime = fDiffFlowCorrelationsHist[t][pe][1]->GetBinContent(b); // <<4'>>
  Double_t sinP1nPsi = fDiffFlowCorrectionTermsForNUAHist[t][pe][0][0]->GetBinContent(b); // <<sin n(Psi)>> 
  Double_t cosP1nPsi = fDiffFlowCorrectionTermsForNUAHist[t][pe][1][0]->GetBinContent(b); // <<cos n(Psi)>> 
  Double_t sinP1nPsi1P1nPhi2 = fDiffFlowCorrectionTermsForNUAHist[t][pe][0][1]->GetBinContent(b); // <<sin n(psi1+phi2)>> 
  Double_t cosP1nPsi1P1nPhi2 = fDiffFlowCorrectionTermsForNUAHist[t][pe][1][1]->GetBinContent(b); // <<cos n(psi1+phi2)>> 
  Double_t sinP1nPsi1P1nPhi2M1nPhi3 = fDiffFlowCorrectionTermsForNUAHist[t][pe][0][2]->GetBinContent(b); // <<sin n(psi1+phi2-phi3)>> 
  Double_t cosP1nPsi1P1nPhi2M1nPhi3 = fDiffFlowCorrectionTermsForNUAHist[t][pe][1][2]->GetBinContent(b); // <<cos n(psi1+phi2-phi3)>> 
  Double_t sinP1nPsi1M1nPhi2M1nPhi3 = fDiffFlowCorrectionTermsForNUAHist[t][pe][0][3]->GetBinContent(b); // <<sin n(psi1-phi2-phi3)>> 
  Double_t cosP1nPsi1M1nPhi2M1nPhi3 = fDiffFlowCorrectionTermsForNUAHist[t][pe][1][3]->GetBinContent(b); // <<cos n(psi1-phi2-phi3)>> 
  // Generalized QC{2'}:
  Double_t qc2Prime = twoPrime - sinP1nPsi*sinP1nPhi - cosP1nPsi*cosP1nPhi;
  if(fApplyCorrectionForNUA)
  {
   fDiffFlowCumulants[t][pe][0]->SetBinContent(b,qc2Prime);
  }
  if(TMath::Abs(twoPrime)>0.)
  {
   fDiffFlowDetectorBias[t][pe][0]->SetBinContent(b,qc2Prime/twoPrime); // detector bias = generalized/isotropic cumulant.   
  }
  // Generalized QC{4'}:
  Double_t qc4Prime = fourPrime-2.*twoPrime*two
                    - cosP1nPsi*cosP1nPhi1M1nPhi2M1nPhi3
                    + sinP1nPsi*sinP1nPhi1M1nPhi2M1nPhi3
                    - cosP1nPhi*cosP1nPsi1M1nPhi2M1nPhi3
                    + sinP1nPhi*sinP1nPsi1M1nPhi2M1nPhi3
                    - 2.*cosP1nPhi*cosP1nPsi1P1nPhi2M1nPhi3
                    - 2.*sinP1nPhi*sinP1nPsi1P1nPhi2M1nPhi3
                    - cosP1nPsi1P1nPhi2*cosP1nPhi1P1nPhi2
                    - sinP1nPsi1P1nPhi2*sinP1nPhi1P1nPhi2
                    + 2.*cosP1nPhi1P1nPhi2*(cosP1nPsi*cosP1nPhi-sinP1nPsi*sinP1nPhi)
                    + 2.*sinP1nPhi1P1nPhi2*(cosP1nPsi*sinP1nPhi+sinP1nPsi*cosP1nPhi)
                    + 4.*two*(cosP1nPsi*cosP1nPhi+sinP1nPsi*sinP1nPhi)
                    + 2.*cosP1nPsi1P1nPhi2*(pow(cosP1nPhi,2.)-pow(sinP1nPhi,2.))
                    + 4.*sinP1nPsi1P1nPhi2*cosP1nPhi*sinP1nPhi
                    + 4.*twoPrime*(pow(cosP1nPhi,2.)+pow(sinP1nPhi,2.))
                    - 6.*(pow(cosP1nPhi,2.)-pow(sinP1nPhi,2.)) 
                    * (cosP1nPsi*cosP1nPhi-sinP1nPsi*sinP1nPhi)
                    - 12.*cosP1nPhi*sinP1nPhi
                    * (sinP1nPsi*cosP1nPhi+cosP1nPsi*sinP1nPhi);
  if(fApplyCorrectionForNUA)
  {
   fDiffFlowCumulants[t][pe][1]->SetBinContent(b,qc4Prime);   
  }
  if(TMath::Abs(fourPrime-2.*twoPrime*two)>0.)
  {
   fDiffFlowDetectorBias[t][pe][1]->SetBinContent(b,qc4Prime/(fourPrime-2.*twoPrime*two)); // detector bias = generalized/isotropic cumulant.   
  }
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateDiffFlowCumulantsCorrectedForNUA(TString type, TString ptOrEta)

//==================================================================================================================================    

void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectedForNUA(TString type, TString ptOrEta)
{
 // Calculate differential flow corrected for non-uniform acceptance.
 
 // to be improved: eventually I will have to access here masured correlations and NUA terms
 //                 instead of cumulants in order to propagate statistical error correctly also 
 //                 to NUA terms (propagating errors directly from cumulants is WRONG for 
 //                 differential flow becuase that doesn't account at all cross-covariance terms) 
 
 // REMARK: When NUA correction is apllied error for differential flow DOES NOT get corrected,
 //         i.e. only value is being corrected, error is still the one relevant for isotropic
 //         case. This eventually will be resolved. 
  
 
 Int_t t = 0; // RP or POI
 Int_t pe = 0; // pt or eta

 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }     
 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   } 
  
 // Common:
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 // Reference Q-cumulants
 Double_t qc2 = fIntFlowQcumulants->GetBinContent(1); // QC{2} 
 Double_t qc4 = fIntFlowQcumulants->GetBinContent(2); // QC{4}
 // Loop over pt or eta bins:
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // Differential Q-cumulants:
  Double_t qc2Prime = fDiffFlowCumulants[t][pe][0]->GetBinContent(b); // QC{2'}
  Double_t qc4Prime = fDiffFlowCumulants[t][pe][1]->GetBinContent(b); // QC{4'}
  // v'{2}:
  if(qc2>0.)
  { 
   Double_t v2Prime = qc2Prime/pow(qc2,0.5);
   if(TMath::Abs(v2Prime)>0.){fDiffFlow[t][pe][0]->SetBinContent(b,v2Prime);} 
  }  
  // v'{4}:
  if(qc4<0.)
  { 
   Double_t v4Prime = -qc4Prime/pow(-qc4,3./4.);
   if(TMath::Abs(v4Prime)>0.){fDiffFlow[t][pe][1]->SetBinContent(b,v4Prime);} 
  }  
 } // end of for(Int_t b=1;b<=fnBinsPtEta[pe];b++)
  
} // end of void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectedForNUA(TString type, TString ptOrEta); 

//==================================================================================================================================

void AliFlowAnalysisWithQCumulants::EvaluateIntFlowCorrelationsWithNestedLoops(AliFlowEventSimple * const anEvent)
{
 // Evaluate with nested loops multiparticle correlations for integrated flow (without using the particle weights). 

 // Remark: Results are stored in profile fIntFlowDirectCorrelations whose binning is organized as follows:
 // 
 //  1st bin: <2>_{1n|1n} = two1n1n = cos(n*(phi1-phi2))>
 //  2nd bin: <2>_{2n|2n} = two2n2n = cos(2n*(phi1-phi2))>
 //  3rd bin: <2>_{3n|3n} = two3n3n = cos(3n*(phi1-phi2))> 
 //  4th bin: <2>_{4n|4n} = two4n4n = cos(4n*(phi1-phi2))>
 //  5th bin:           ----  EMPTY ----
 //  6th bin: <3>_{2n|1n,1n} = three2n1n1n = <cos(n*(2.*phi1-phi2-phi3))>
 //  7th bin: <3>_{3n|2n,1n} = three3n2n1n = <cos(n*(3.*phi1-2.*phi2-phi3))>
 //  8th bin: <3>_{4n|2n,2n} = three4n2n2n = <cos(n*(4.*phi1-2.*phi2-2.*phi3))>
 //  9th bin: <3>_{4n|3n,1n} = three4n3n1n = <cos(n*(4.*phi1-3.*phi2-phi3))>
 // 10th bin:           ----  EMPTY ----
 // 11th bin: <4>_{1n,1n|1n,1n} = four1n1n1n1n = <cos(n*(phi1+phi2-phi3-phi4))>
 // 12th bin: <4>_{2n,1n|2n,1n} = four2n1n2n1n = <cos(2.*n*(phi1+phi2-phi3-phi4))>
 // 13th bin: <4>_{2n,2n|2n,2n} = four2n2n2n2n = <cos(n*(2.*phi1+phi2-2.*phi3-phi4))>
 // 14th bin: <4>_{3n|1n,1n,1n} = four3n1n1n1n = <cos(n*(3.*phi1-phi2-phi3-phi4))> 
 // 15th bin: <4>_{3n,1n|3n,1n} = four3n1n3n1n = <cos(n*(4.*phi1-2.*phi2-phi3-phi4))>
 // 16th bin: <4>_{3n,1n|2n,2n} = four3n1n2n2n = <cos(n*(3.*phi1+phi2-2.*phi3-2.*phi4))>
 // 17th bin: <4>_{4n|2n,1n,1n} = four4n2n1n1n = <cos(n*(3.*phi1+phi2-3.*phi3-phi4))> 
 // 18th bin:           ----  EMPTY ----
 // 19th bin: <5>_{2n|1n,1n,1n,1n} = five2n1n1n1n1n = <cos(n*(2.*phi1+phi2-phi3-phi4-phi5))>
 // 20th bin: <5>_{2n,2n|2n,1n,1n} = five2n2n2n1n1n = <cos(n*(2.*phi1+2.*phi2-2.*phi3-phi4-phi5))>
 // 21st bin: <5>_{3n,1n|2n,1n,1n} = five3n1n2n1n1n = <cos(n*(3.*phi1+phi2-2.*phi3-phi4-phi5))>
 // 22nd bin: <5>_{4n|1n,1n,1n,1n} = five4n1n1n1n1n = <cos(n*(4.*phi1-phi2-phi3-phi4-phi5))>
 // 23rd bin:           ----  EMPTY ----
 // 24th bin: <6>_{1n,1n,1n|1n,1n,1n} = six1n1n1n1n1n1n = <cos(n*(phi1+phi2+phi3-phi4-phi5-phi6))>
 // 25th bin: <6>_{2n,1n,1n|2n,1n,1n} = six2n1n1n2n1n1n = <cos(n*(2.*phi1+2.*phi2-phi3-phi4-phi5-phi6))>
 // 26th bin: <6>_{2n,2n|1n,1n,1n,1n} = six2n2n1n1n1n1n = <cos(n*(3.*phi1+phi2-phi3-phi4-phi5-phi6))>
 // 27th bin: <6>_{3n,1n|1n,1n,1n,1n} = six3n1n1n1n1n1n = <cos(n*(2.*phi1+phi2+phi3-2.*phi4-phi5-phi6))>
 // 28th bin:           ----  EMPTY ----
 // 29th bin: <7>_{2n,1n,1n|1n,1n,1n,1n} = seven2n1n1n1n1n1n1n =  <cos(n*(2.*phi1+phi2+phi3-phi4-phi5-phi6-phi7))>
 // 30th bin:           ----  EMPTY ----
 // 31st bin: <8>_{1n,1n,1n,1n|1n,1n,1n,1n} = eight1n1n1n1n1n1n1n1n = <cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8))>
 // 32nd bin:           ----  EMPTY ----
  //  Extra correlations for 3p TY study: 
 // 33rd bin: <4>_{4n,2n|3n,3n}= four4n2n3n3n = <cos(n*(4.*phi1+2.*phi2-3.*phi3-3.*phi4))>
 // 34th bin: <5>_{3n,3n|2n,2n,2n} = five3n3n2n2n2n = <cos(n(3*phi1+3*phi2-2*phi3-2*phi4-2*phi5))> 
  //  Extra correlations for 6p TY study: 
 // 35th bin: <2>_{5n|5n} = two5n5n = <cos(5n*(phi1-phi2)> T
 // 36th bin: <2>_{6n|6n} = two6n6n = <cos(6n*(phi1-phi2)> T
 // 37th bin: <3>_{5n|3n,2n} = three5n3n2n = <cos(n*(5*phi1-3*phi2-2*phi3)> 
 // 38th bin: <3>_{5n|4n,1n} = three5n4n1n = <cos(n*(5*phi1-4*phi2-1*phi3)> 
 // 39th bin: <3>_{6n|3n,3n} = three6n3n3n = <cos(n*(6*phi1-3*phi2-3*phi3)> T 
 // 40th bin: <3>_{6n|4n,2n} = three6n4n2n = <cos(n*(6*phi1-4*phi2-2*phi3)> T
 // 41st bin: <3>_{6n|5n,1n} = three6n5n1n = <cos(n*(6*phi1-5*phi2-1*phi3)>
 // 42nd bin: <4>_{6n|3n,2n,1n} = four6n3n2n1n = <cos(n*(6*phi1-3*phi2-2*phi3-1*phi4)>
 // 43rd bin: <4>_{3n,2n|3n,2n} = four3n2n3n2n = <cos(n*(3*phi1+2*phi2-3*phi3-2*phi4)>
 // 44th bin: <4>_{4n,1n|3n,2n} = four4n1n3n2n = <cos(n*(4*phi1+1*phi2-3*phi3-2*phi4)>
 // 45th bin: <4>_{3n,3n|3n,3n} = four3n3n3n3n = <cos(3.*n*(phi1+phi2-phi3-phi4))> T
 // 46th bin: <4>_{4n,2n|3n,3n} = four4n2n3n3n = <cos(n*(4*phi1+2*phi2-3*phi3-3*phi4)>
 // 47th bin: <4>_{5n,1n|3n,3n} = four5n1n3n3n = <cos(n*(5*phi1+1*phi2-3*phi3-3*phi4)>
 // 48th bin: <4>_{4n,2n|4n,2n} = four4n2n4n2n = <cos(n*(4*phi1+2*phi2-4*phi3-2*phi4)> T
 // 49th bin: <4>_{5n,1n|4n,2n} = four5n1n4n2n = <cos(n*(5*phi1+1*phi2-4*phi3-2*phi4)>
 // 50th bin: <4>_{5n|3n,1n,1n} = four5n3n1n1n = <cos(n*(5*phi1-3*phi2-1*phi3-1*phi4)>
 // 51st bin: <4>_{5n|2n,2n,1n} = four5n2n2n1n = <cos(n*(5*phi1-2*phi2-2*phi3-1*phi4)>
 // 52nd bin: <4>_{5n,1n|5n,1n} = four5n1n5n1n = <cos(n*(5*phi1+1*phi2-5*phi3-1*phi4)>
 // 53rd bin: <5>_{3n,3n|3n,2n,1n} = four3n3n3n2n1n = <cos(n*(3*phi1+3*phi2-3*phi3-2*phi4-1*phi5)>
 // 54th bin: <5>_{4n,2n|3n,2n,1n} = four4n2n3n2n1n = <cos(n*(4*phi1+2*phi2-3*phi3-2*phi4-1*phi5)>
 // 55th bin: <5>_{3n,2n|3n,1n,1n} = four3n2n3n1n1n = <cos(n*(3*phi1+2*phi2-3*phi3-1*phi4-1*phi5)>
 // 56th bin: <5>_{3n,2n|2n,2n,1n} = four3n2n2n2n1n = <cos(n*(3*phi1+2*phi2-2*phi3-2*phi4-1*phi5)>
 // 57th bin: <5>_{5n,1n|3n,2n,1n} = four5n1n3n2n1n = <cos(n*(5*phi1+1*phi2-3*phi3-2*phi4-1*phi5)>
 // 58th bin: <6>_{3n,2n,1n|3n,2n,1n} = six3n2n1n3n2n1n = <cos(n*(3*phi1+2*phi2+1*phi3-3*phi4-2*phi5-1*phi6)>
  
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL; 
 Double_t phi1=0., phi2=0., phi3=0., phi4=0., phi5=0., phi6=0., phi7=0., phi8=0.; 
 Int_t n = fHarmonic; 
 Int_t eventNo = (Int_t)fAvMultiplicity->GetBinEntries(1); // to be improved (is this casting safe in general?)
 Double_t dMult = (*fSpk)(0,0);
 cout<<endl;
 cout<<"Multiparticle correlations: Event number: "<<eventNo<<", multiplicity is "<<dMult<<endl;
 if(dMult<2)
 {
  cout<<"... skipping this event (multiplicity too low) ..."<<endl;
 } else if (dMult>fMaxAllowedMultiplicity)
   {
    cout<<"... skipping this event (multiplicity too high) ..."<<endl;
   } else 
     { 
      cout<<"... evaluating nested loops (without using particle weights)..."<<endl;
     } 
 
 // 2-particle correlations:       
 if(nPrim>=2 && nPrim<=fMaxAllowedMultiplicity)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi(); 
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    if(nPrim==2) cout<<i1<<" "<<i2<<"\r"<<flush;
    // fill the profile with 2-p correlations: 
    fIntFlowDirectCorrelations->Fill(0.5,cos(n*(phi1-phi2)),1.);     // <cos(n*(phi1-phi2))>
    fIntFlowDirectCorrelations->Fill(1.5,cos(2.*n*(phi1-phi2)),1.);  // <cos(2n*(phi1-phi2))>
    fIntFlowDirectCorrelations->Fill(2.5,cos(3.*n*(phi1-phi2)),1.);  // <cos(3n*(phi1-phi2))>
    fIntFlowDirectCorrelations->Fill(3.5,cos(4.*n*(phi1-phi2)),1.);  // <cos(4n*(phi1-phi2))>   
    fIntFlowDirectCorrelations->Fill(34.5,cos(5.*n*(phi1-phi2)),1.); // <cos(5n*(phi1-phi2))>
    fIntFlowDirectCorrelations->Fill(35.5,cos(6.*n*(phi1-phi2)),1.); // <cos(6n*(phi1-phi2))>   
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=2)
 
 // 3-particle correlations:         
 if(nPrim>=3 && nPrim<=fMaxAllowedMultiplicity)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     if(nPrim==3) cout<<i1<<" "<<i2<<" "<<i3<<"\r"<<flush;
     // fill the profile with 3-p correlations:   
     fIntFlowDirectCorrelations->Fill(5.,cos(2.*n*phi1-n*(phi2+phi3)),1.);         //<3>_{2n|nn,n}
     fIntFlowDirectCorrelations->Fill(6.,cos(3.*n*phi1-2.*n*phi2-n*phi3),1.);      //<3>_{3n|2n,n}
     fIntFlowDirectCorrelations->Fill(7.,cos(4.*n*phi1-2.*n*phi2-2.*n*phi3),1.);   //<3>_{4n|2n,2n}
     fIntFlowDirectCorrelations->Fill(8.,cos(4.*n*phi1-3.*n*phi2-n*phi3),1.);      //<3>_{4n|3n,n}
     fIntFlowDirectCorrelations->Fill(36.5,cos(5.*n*phi1-3.*n*phi2-2.*n*phi3),1.); //<3>_{5n|3n,2n}
     fIntFlowDirectCorrelations->Fill(37.5,cos(5.*n*phi1-4.*n*phi2-1.*n*phi3),1.); //<3>_{5n|4n,1n}
     fIntFlowDirectCorrelations->Fill(38.5,cos(6.*n*phi1-3.*n*phi2-3.*n*phi3),1.); //<3>_{6n|3n,3n}
     fIntFlowDirectCorrelations->Fill(39.5,cos(6.*n*phi1-4.*n*phi2-2.*n*phi3),1.); //<3>_{6n|4n,2n}     
     fIntFlowDirectCorrelations->Fill(40.5,cos(6.*n*phi1-5.*n*phi2-1.*n*phi3),1.); //<3>_{6n|5n,1n}
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=3)

 // 4-particle correlations:
 if(nPrim>=4 && nPrim<=fMaxAllowedMultiplicity)
 {       
  for(Int_t i1=0;i1<nPrim;i1++)
  { 
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3)continue;
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())) continue;
      phi4=aftsTrack->Phi();
      if(nPrim==4) cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<"\r"<<flush;
      // fill the profile with 4-p correlations:   
      fIntFlowDirectCorrelations->Fill(10.,cos(n*phi1+n*phi2-n*phi3-n*phi4),1.);            // <4>_{n,n|n,n} 
      fIntFlowDirectCorrelations->Fill(11.,cos(2.*n*phi1+n*phi2-2.*n*phi3-n*phi4),1.);      // <4>_{2n,n|2n,n}
      fIntFlowDirectCorrelations->Fill(12.,cos(2.*n*phi1+2*n*phi2-2.*n*phi3-2.*n*phi4),1.); // <4>_{2n,2n|2n,2n}
      fIntFlowDirectCorrelations->Fill(13.,cos(3.*n*phi1-n*phi2-n*phi3-n*phi4),1.);         // <4>_{3n|n,n,n}
      fIntFlowDirectCorrelations->Fill(14.,cos(3.*n*phi1+n*phi2-3.*n*phi3-n*phi4),1.);      // <4>_{3n,n|3n,n}   
      fIntFlowDirectCorrelations->Fill(15.,cos(3.*n*phi1+n*phi2-2.*n*phi3-2.*n*phi4),1.);   // <4>_{3n,n|2n,2n}
      fIntFlowDirectCorrelations->Fill(16.,cos(4.*n*phi1-2.*n*phi2-n*phi3-n*phi4),1.);      // <4>_{4n|2n,n,n}     
      fIntFlowDirectCorrelations->Fill(32.,cos(n*(4.*phi1+2.*phi2-3.*phi3-3.*phi4)),1.);    // <4>_{4n,2n|3n,3n}        
      fIntFlowDirectCorrelations->Fill(41.5,cos(n*(6.*phi1-3.*phi2-2.*phi3-1.*phi4)),1.);    // <4>_{6n|3n,2n,1n}  
      fIntFlowDirectCorrelations->Fill(42.5,cos(n*(3.*phi1+2.*phi2-3.*phi3-2.*phi4)),1.);    // <4>_{3n,2n|3n,2n}  
      fIntFlowDirectCorrelations->Fill(43.5,cos(n*(4.*phi1+1.*phi2-3.*phi3-2.*phi4)),1.);    // <4>_{4n,1n|3n,2n}  
      fIntFlowDirectCorrelations->Fill(44.5,cos(n*(3.*phi1+3.*phi2-3.*phi3-3.*phi4)),1.);    // <4>_{3n,3n|3n,3n}  
      fIntFlowDirectCorrelations->Fill(45.5,cos(n*(4.*phi1+2.*phi2-3.*phi3-3.*phi4)),1.);    // <4>_{4n,2n|3n,3n}  
      fIntFlowDirectCorrelations->Fill(46.5,cos(n*(5.*phi1+1.*phi2-3.*phi3-3.*phi4)),1.);    // <4>_{5n,1n|3n,3n}  
      fIntFlowDirectCorrelations->Fill(47.5,cos(n*(4.*phi1+2.*phi2-4.*phi3-2.*phi4)),1.);    // <4>_{4n,2n|4n,2n}     
      fIntFlowDirectCorrelations->Fill(48.5,cos(n*(5.*phi1+1.*phi2-4.*phi3-2.*phi4)),1.);    // <4>_{5n,1n|4n,2n}  
      fIntFlowDirectCorrelations->Fill(49.5,cos(n*(5.*phi1-3.*phi2-1.*phi3-1.*phi4)),1.);    // <4>_{5n|3n,1n,1n}  
      fIntFlowDirectCorrelations->Fill(50.5,cos(n*(5.*phi1-2.*phi2-2.*phi3-1.*phi4)),1.);    // <4>_{5n|2n,2n,1n}  
      fIntFlowDirectCorrelations->Fill(51.5,cos(n*(5.*phi1+1.*phi2-5.*phi3-1.*phi4)),1.);    // <4>_{5n,1n|5n,1n}        
      fIntFlowDirectCorrelations->Fill(58.5,cos(n*(6.*phi1-4.*phi2-1.*phi3-1.*phi4)),1.);    // <4>_{6n|4n,1n,1n}  
      fIntFlowDirectCorrelations->Fill(59.5,cos(n*(6.*phi1-2.*phi2-2.*phi3-2.*phi4)),1.);    // <4>_{6n|2n,2n,2n}  
     } // end of for(Int_t i4=0;i4<nPrim;i4++) 
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=)

 // 5-particle correlations:      
 if(nPrim>=5 && nPrim<=fMaxAllowedMultiplicity)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;  
   phi1=aftsTrack->Phi();
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3)continue;
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())) continue;
      phi4=aftsTrack->Phi();
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())) continue;
       phi5=aftsTrack->Phi();
       if(nPrim==5) cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<i5<<"\r"<<flush;
       // fill the profile with 5-p correlations:   
       fIntFlowDirectCorrelations->Fill(18.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),1.);              // <5>_{2n,n|n,n,n}
       fIntFlowDirectCorrelations->Fill(19.,cos(2.*n*phi1+2.*n*phi2-2.*n*phi3-n*phi4-n*phi5),1.);        // <5>_{2n,2n|2n,n,n}
       fIntFlowDirectCorrelations->Fill(20.,cos(3.*n*phi1+n*phi2-2.*n*phi3-n*phi4-n*phi5),1.);           // <5>_{3n,n|2n,n,n}
       fIntFlowDirectCorrelations->Fill(21.,cos(4.*n*phi1-n*phi2-n*phi3-n*phi4-n*phi5),1.);              // <5>_{4n|n,n,n,n}
       fIntFlowDirectCorrelations->Fill(33.,cos(3.*n*phi1+3.*n*phi2-2.*n*phi3-2.*n*phi4-2.*n*phi5),1.);  // <5>_{3n,3n|2n,2n,2n}       
       fIntFlowDirectCorrelations->Fill(52.5,cos(3.*n*phi1+3.*n*phi2-3.*n*phi3-2.*n*phi4-1.*n*phi5),1.); // <5>_{3n,3n|3n,2n,1n}      
       fIntFlowDirectCorrelations->Fill(53.5,cos(4.*n*phi1+2.*n*phi2-3.*n*phi3-2.*n*phi4-1.*n*phi5),1.); // <5>_{4n,2n|3n,2n,1n}       
       fIntFlowDirectCorrelations->Fill(54.5,cos(3.*n*phi1+2.*n*phi2-3.*n*phi3-1.*n*phi4-1.*n*phi5),1.); // <5>_{3n,2n|3n,1n,1n}
       fIntFlowDirectCorrelations->Fill(55.5,cos(3.*n*phi1+2.*n*phi2-2.*n*phi3-2.*n*phi4-1.*n*phi5),1.); // <5>_{3n,2n|2n,2n,1n}
       fIntFlowDirectCorrelations->Fill(56.5,cos(5.*n*phi1+1.*n*phi2-3.*n*phi3-2.*n*phi4-1.*n*phi5),1.); // <5>_{5n,1n|3n,2n,1n}              
       fIntFlowDirectCorrelations->Fill(60.5,cos(6.*n*phi1-2.*n*phi2-2.*n*phi3-1.*n*phi4-1.*n*phi5),1.); // <5>_{6n|2n,2n,1n,1n}
       fIntFlowDirectCorrelations->Fill(61.5,cos(4.*n*phi1+1.*n*phi2+1.*n*phi3-3.*n*phi4-3.*n*phi5),1.); // <5>_{4n,1n,1n|3n,3n}              
      } // end of for(Int_t i5=0;i5<nPrim;i5++)
     } // end of for(Int_t i4=0;i4<nPrim;i4++)  
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=5)
  
 // 6-particle correlations:
 if(nPrim>=6 && nPrim<=fMaxAllowedMultiplicity)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3)continue;
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())) continue;
      phi4=aftsTrack->Phi();
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())) continue;
       phi5=aftsTrack->Phi();
       for(Int_t i6=0;i6<nPrim;i6++)
       {
        if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
        aftsTrack=anEvent->GetTrack(i6);
        if(!(aftsTrack->InRPSelection())) continue;
        phi6=aftsTrack->Phi(); 
        if(nPrim==6) cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<i5<<" "<<i6<<"\r"<<flush;
        // fill the profile with 6-p correlations:   
        fIntFlowDirectCorrelations->Fill(23.,cos(n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6),1.);                    // <6>_{1n,1n,1n|1n,1n,1n}
        fIntFlowDirectCorrelations->Fill(24.,cos(2.*n*phi1+n*phi2+n*phi3-2.*n*phi4-n*phi5-n*phi6),1.);              // <6>_{2n,1n,1n|2n,1n,1n}
        fIntFlowDirectCorrelations->Fill(25.,cos(2.*n*phi1+2.*n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),1.);              // <6>_{2n,2n|1n,1n,1n,1n}
        fIntFlowDirectCorrelations->Fill(26.,cos(3.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),1.);                 // <6>_{3n,1n|1n,1n,1n,1n}  
        fIntFlowDirectCorrelations->Fill(57.5,cos(3.*n*phi1+2.*n*phi2+1.*n*phi3-3.*n*phi4-2.*n*phi5-1.*n*phi6),1.); // <6>_{3n,2n,1n|3n,2n,1n}  
        fIntFlowDirectCorrelations->Fill(62.5,cos(3.*n*phi1+3.*n*phi2-2.*n*phi3-2.*n*phi4-1.*n*phi5-1.*n*phi6),1.); // <6>_{3n,3n|2n,2n,1n,1n}  
       } // end of for(Int_t i6=0;i6<nPrim;i6++)
      } // end of for(Int_t i5=0;i5<nPrim;i5++)
     } // end of for(Int_t i4=0;i4<nPrim;i4++)
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=6)
  
 // 7-particle correlations:
 if(nPrim>=7 && nPrim<=fMaxAllowedMultiplicity)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  { 
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3)continue;
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())) continue;
      phi4=aftsTrack->Phi();
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())) continue;
       phi5=aftsTrack->Phi();
       for(Int_t i6=0;i6<nPrim;i6++)
       {
        if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
        aftsTrack=anEvent->GetTrack(i6);
        if(!(aftsTrack->InRPSelection())) continue;
        phi6=aftsTrack->Phi(); 
        for(Int_t i7=0;i7<nPrim;i7++)
        {
         if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
         aftsTrack=anEvent->GetTrack(i7);
         if(!(aftsTrack->InRPSelection())) continue;
         phi7=aftsTrack->Phi(); 
         if(nPrim==7) cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<i5<<" "<<i6<<" "<<i7<<"\r"<<flush;
         // fill the profile with 7-p correlation:   
         fIntFlowDirectCorrelations->Fill(28.,cos(2.*n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6-n*phi7),1.); // <7>_{2n,n,n|n,n,n,n}
        } // end of for(Int_t i7=0;i7<nPrim;i7++)
       } // end of for(Int_t i6=0;i6<nPrim;i6++) 
      } // end of for(Int_t i5=0;i5<nPrim;i5++)
     } // end of for(Int_t i4=0;i4<nPrim;i4++)  
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=7)
 
 // 8-particle correlations:
 if(nPrim>=8 && nPrim<=fMaxAllowedMultiplicity)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3)continue;
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())) continue;
      phi4=aftsTrack->Phi();
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())) continue;
       phi5=aftsTrack->Phi();
       for(Int_t i6=0;i6<nPrim;i6++)
       {
        if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
        aftsTrack=anEvent->GetTrack(i6);
        if(!(aftsTrack->InRPSelection())) continue;
        phi6=aftsTrack->Phi();
        for(Int_t i7=0;i7<nPrim;i7++)
        {
         if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
         aftsTrack=anEvent->GetTrack(i7);
         if(!(aftsTrack->InRPSelection())) continue;
         phi7=aftsTrack->Phi();
         for(Int_t i8=0;i8<nPrim;i8++)
         {
          if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7)continue;
          aftsTrack=anEvent->GetTrack(i8);
          if(!(aftsTrack->InRPSelection())) continue;
          phi8=aftsTrack->Phi();
          cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<i5<<" "<<i6<<" "<<i7<<" "<<i8<<"\r"<<flush;
          // fill the profile with 8-p correlation:   
          fIntFlowDirectCorrelations->Fill(30.,cos(n*phi1+n*phi2+n*phi3+n*phi4-n*phi5-n*phi6-n*phi7-n*phi8),1.); // <8>_{n,n,n,n|n,n,n,n}
         } // end of for(Int_t i8=0;i8<nPrim;i8++)
        } // end of for(Int_t i7=0;i7<nPrim;i7++) 
       } // end of for(Int_t i6=0;i6<nPrim;i6++) 
      } // end of for(Int_t i5=0;i5<nPrim;i5++)
     } // end of for(Int_t i4=0;i4<nPrim;i4++)  
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=8)
 
 cout<<endl;

} // end of AliFlowAnalysisWithQCumulants::EvaluateIntFlowCorrelationsWithNestedLoops(AliFlowEventSimple* anEvent)


//==================================================================================================================================


void AliFlowAnalysisWithQCumulants::CrossCheckIntFlowCorrelations()
{
 // Cross-check results for multiparticle correlations needed for int. flow: results from Q-vectors vs results from nested loops.

 cout<<endl;
 cout<<endl;
 cout<<"   *****************************************"<<endl;
 cout<<"   **** cross-checking the correlations ****"<<endl;
 cout<<"   ****       for integrated flow       ****"<<endl;
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
 {
  cout<<"   ****   (particle weights not used)   ****"<<endl;
 } else
   {
    cout<<"   ****     (particle weights used)     ****"<<endl;
   } 
 cout<<"   *****************************************"<<endl;
 cout<<endl;
 cout<<endl;

 Int_t ciMax = 64; // to be improved (removed eventually when I calculate 6th and 8th order with particle weights)
 
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights)
 {
  ciMax = 11;
 }

 for(Int_t ci=1;ci<=ciMax;ci++)
 {
  if(strcmp((fIntFlowCorrelationsAllPro->GetXaxis())->GetBinLabel(ci), "") == 0) continue; // to be improved (access finalized histogram here)
  cout<<(fIntFlowCorrelationsAllPro->GetXaxis())->GetBinLabel(ci)<<":"<<endl; // to be improved (access finalized histogram here)
  cout<<"from Q-vectors    = "<<fIntFlowCorrelationsAllPro->GetBinContent(ci)<<endl; // to be improved (access finalized histogram here)
  cout<<"from nested loops = "<<fIntFlowDirectCorrelations->GetBinContent(ci)<<endl;
  cout<<endl;
 }
  
} // end of void AliFlowAnalysisWithQCumulants::CrossCheckIntFlowCorrelations()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CrossCheckIntFlowCorrectionTermsForNUA()
{
 // Cross-check results for corrections terms for non-uniform acceptance needed for int. flow: results from Q-vectors vs results from nested loops.

 cout<<endl;
 cout<<endl;
 cout<<"   *********************************************"<<endl;
 cout<<"   **** cross-checking the correction terms ****"<<endl;
 cout<<"   **** for non-uniform acceptance relevant ****"<<endl;
 cout<<"   ****         for integrated flow         ****"<<endl;
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
 {
  cout<<"   ****     (particle weights not used)     ****"<<endl;
 } else
   {
    cout<<"   ****       (particle weights used)       ****"<<endl;
   } 
 cout<<"   *********************************************"<<endl;
 cout<<endl;
 cout<<endl;

 for(Int_t ci=1;ci<=4;ci++) // correction term index (to be improved - hardwired 4)
 {
  for(Int_t sc=0;sc<2;sc++) // sin or cos term
  {
   if(strcmp((fIntFlowCorrectionTermsForNUAPro[sc]->GetXaxis())->GetBinLabel(ci), "") == 0) continue; // to be improved (access finalized histogram here)
   cout<<(fIntFlowCorrectionTermsForNUAPro[sc]->GetXaxis())->GetBinLabel(ci)<<":"<<endl; // to be improved (access finalized histogram here)
   cout<<"from Q-vectors    = "<<fIntFlowCorrectionTermsForNUAPro[sc]->GetBinContent(ci)<<endl; // to be improved (access finalized histogram here)
   cout<<"from nested loops = "<<fIntFlowDirectCorrectionTermsForNUA[sc]->GetBinContent(ci)<<endl;
   cout<<endl;
  } // end of for(Int_t sc=0;sc<2;sc++) // sin or cos term
 } // end of for(Int_t ci=1;ci<=10;ci++) // correction term index
  
} // end of void AliFlowAnalysisWithQCumulants::CrossCheckIntFlowCorrectionTermsForNUA() 

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::EvaluateIntFlowCorrelationsWithNestedLoopsUsingParticleWeights(AliFlowEventSimple * const anEvent)
{
 // Evaluate with nested loops multiparticle correlations for integrated flow (using the particle weights). 

 // Results are stored in profile fIntFlowDirectCorrelations. 
 // Remark 1: When particle weights are used the binning of fIntFlowDirectCorrelations is organized as follows:
 //
 //  1st bin: <2>_{1n|1n} = two1n1nW1W1 = <w1 w2 cos(n*(phi1-phi2))>
 //  2nd bin: <2>_{2n|2n} = two2n2nW2W2 = <w1^2 w2^2 cos(2n*(phi1-phi2))>
 //  3rd bin: <2>_{3n|3n} = two3n3nW3W3 = <w1^3 w2^3 cos(3n*(phi1-phi2))> 
 //  4th bin: <2>_{4n|4n} = two4n4nW4W4 = <w1^4 w2^4 cos(4n*(phi1-phi2))>
 //  5th bin:           ----  EMPTY ----
 //  6th bin: <3>_{2n|1n,1n} = three2n1n1nW2W1W1 = <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))>
 //  7th bin: <3>_{3n|2n,1n} = ...
 //  8th bin: <3>_{4n|2n,2n} = ...
 //  9th bin: <3>_{4n|3n,1n} = ...
 // 10th bin:           ----  EMPTY ----
 // 11th bin: <4>_{1n,1n|1n,1n} = four1n1n1n1nW1W1W1W1 = <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
 // 12th bin: <4>_{2n,1n|2n,1n} = ...
 // 13th bin: <4>_{2n,2n|2n,2n} = ...
 // 14th bin: <4>_{3n|1n,1n,1n} = ... 
 // 15th bin: <4>_{3n,1n|3n,1n} = ...
 // 16th bin: <4>_{3n,1n|2n,2n} = ...
 // 17th bin: <4>_{4n|2n,1n,1n} = ... 
 // 18th bin:           ----  EMPTY ----
 // 19th bin: <5>_{2n|1n,1n,1n,1n} = ...
 // 20th bin: <5>_{2n,2n|2n,1n,1n} = ...
 // 21st bin: <5>_{3n,1n|2n,1n,1n} = ...
 // 22nd bin: <5>_{4n|1n,1n,1n,1n} = ...
 // 23rd bin:           ----  EMPTY ----
 // 24th bin: <6>_{1n,1n,1n|1n,1n,1n} = ...
 // 25th bin: <6>_{2n,1n,1n|2n,1n,1n} = ...
 // 26th bin: <6>_{2n,2n|1n,1n,1n,1n} = ...
 // 27th bin: <6>_{3n,1n|1n,1n,1n,1n} = ...
 // 28th bin:           ----  EMPTY ----
 // 29th bin: <7>_{2n,1n,1n|1n,1n,1n,1n} = ...
 // 30th bin:           ----  EMPTY ----
 // 31st bin: <8>_{1n,1n,1n,1n|1n,1n,1n,1n} = ...
 
 // Remark 2: When particle weights are used there are some extra correlations. They are stored in 
 // fIntFlowExtraDirectCorrelations binning of which is organized as follows:
 
 // 1st bin: two1n1nW3W1 = <w1^3 w2 cos(n*(phi1-phi2))>
 // 2nd bin: two1n1nW1W1W2 = <w1 w2 w3^2 cos(n*(phi1-phi2))>  
 // ...
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 //Double_t phi1=0., phi2=0., phi3=0., phi4=0., phi5=0., phi6=0., phi7=0., phi8=0.;
 //Double_t wPhi1=1., wPhi2=1., wPhi3=1., wPhi4=1., wPhi5=1., wPhi6=1., wPhi7=1., wPhi8=1.;
 Double_t phi1=0., phi2=0., phi3=0., phi4=0.;
 Double_t wPhi1=1., wPhi2=1., wPhi3=1., wPhi4=1.;
 Int_t n = fHarmonic; 
 Int_t eventNo = (Int_t)fAvMultiplicity->GetBinEntries(1); // to be improved (is this casting safe in general?)
 Double_t dMult = (*fSpk)(0,0);
 cout<<endl;
 cout<<"Multiparticle correlations: Event number: "<<eventNo<<", multiplicity is "<<dMult<<endl;
 if(dMult<2)
 {
  cout<<"... skipping this event (multiplicity too low) ..."<<endl;
 } else if (dMult>fMaxAllowedMultiplicity)
   {
    cout<<"... skipping this event (multiplicity too high) ..."<<endl;
   } else 
     { 
      cout<<"... evaluating nested loops (using particle weights) ..."<<endl;
     } 
      
 // 2-particle correlations:       
 if(nPrim>=2 && nPrim<=fMaxAllowedMultiplicity)
 {
  // 2 nested loops multiparticle correlations using particle weights:       
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*fnBinsPhi/TMath::TwoPi())));
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));   
    if(nPrim==2) cout<<i1<<" "<<i2<<"\r"<<flush;
    // 2-p correlations using particle weights:
    if(fUsePhiWeights) fIntFlowDirectCorrelations->Fill(0.5,cos(n*(phi1-phi2)),wPhi1*wPhi2);                  // <w1   w2   cos( n*(phi1-phi2))>
    if(fUsePhiWeights) fIntFlowDirectCorrelations->Fill(1.5,cos(2.*n*(phi1-phi2)),pow(wPhi1,2)*pow(wPhi2,2)); // <w1^2 w2^2 cos(2n*(phi1-phi2))>
    if(fUsePhiWeights) fIntFlowDirectCorrelations->Fill(2.5,cos(3.*n*(phi1-phi2)),pow(wPhi1,3)*pow(wPhi2,3)); // <w1^3 w2^3 cos(3n*(phi1-phi2))>
    if(fUsePhiWeights) fIntFlowDirectCorrelations->Fill(3.5,cos(4.*n*(phi1-phi2)),pow(wPhi1,4)*pow(wPhi2,4)); // <w1^4 w2^4 cos(4n*(phi1-phi2))> 
    // extra correlations: 
    // 2-p extra correlations (do not appear if particle weights are not used):
    if(fUsePhiWeights) fIntFlowExtraDirectCorrelations->Fill(0.5,cos(n*(phi1-phi2)),pow(wPhi1,3)*wPhi2); // <w1^3 w2 cos(n*(phi1-phi2))>
    // ...
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=2)

 if(nPrim>=3 && nPrim<=fMaxAllowedMultiplicity)
 { 
  // 3 nested loops multiparticle correlations using particle weights:       
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*fnBinsPhi/TMath::TwoPi())));
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     if(fUsePhiWeights && fPhiWeights) wPhi3 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*fnBinsPhi/TMath::TwoPi())));
     if(nPrim==3) cout<<i1<<" "<<i2<<" "<<i3<<"\r"<<flush;
     // 3-p correlations using particle weights:
     if(fUsePhiWeights) fIntFlowDirectCorrelations->Fill(5.5,cos(2.*n*phi1-n*(phi2+phi3)),pow(wPhi1,2)*wPhi2*wPhi3); // <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))>
     // ...
     // extra correlations: 
     // 2-p extra correlations (do not appear if particle weights are not used):
      if(fUsePhiWeights) fIntFlowExtraDirectCorrelations->Fill(1.5,cos(n*(phi1-phi2)),wPhi1*wPhi2*pow(wPhi3,2)); // <w1 w2 w3^2 cos(n*(phi1-phi2))>
     // ...
     // 3-p extra correlations (do not appear if particle weights are not used):
     // ...
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=3)
 
 if(nPrim>=4 && nPrim<=fMaxAllowedMultiplicity)
 {
  // 4 nested loops multiparticle correlations using particle weights:       
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*fnBinsPhi/TMath::TwoPi())));
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     if(fUsePhiWeights && fPhiWeights) wPhi3 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*fnBinsPhi/TMath::TwoPi())));
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3)continue;
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())) continue;
      phi4=aftsTrack->Phi();
      if(fUsePhiWeights && fPhiWeights) wPhi4 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*fnBinsPhi/TMath::TwoPi())));
      if(nPrim>=4) cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<"\r"<<flush; // to be improved (replace eventually this if statement with if(nPrim==4))
      // 4-p correlations using particle weights:
      if(fUsePhiWeights) fIntFlowDirectCorrelations->Fill(10.5,cos(n*phi1+n*phi2-n*phi3-n*phi4),wPhi1*wPhi2*wPhi3*wPhi4); 
      // extra correlations: 
      // 2-p extra correlations (do not appear if particle weights are not used):
      // ...
      // 3-p extra correlations (do not appear if particle weights are not used):
      // ...
      // 4-p extra correlations (do not appear if particle weights are not used):
      // ...
     } // end of for(Int_t i4=0;i4<nPrim;i4++) 
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=4)

 cout<<endl; 

} // end of void AliFlowAnalysisWithQCumulants::EvaluateIntFlowCorrelationsWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* anEvent)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CrossCheckIntFlowExtraCorrelations()
{
 // Cross-check results for extra multiparticle correlations needed for int. flow 
 // which appear only when particle weights are used: results from Q-vectors vs results from nested loops.

 cout<<endl;
 cout<<endl;
 cout<<"   ***********************************************"<<endl;
 cout<<"   **** cross-checking the extra correlations ****"<<endl;
 cout<<"   ****          for integrated flow          ****"<<endl;
 cout<<"   ***********************************************"<<endl;
 cout<<endl;
 cout<<endl;
 
 for(Int_t eci=1;eci<=2;eci++) // to be improved (increased eciMax eventually when I calculate 6th and 8th)
 {
  if(strcmp((fIntFlowExtraCorrelationsPro->GetXaxis())->GetBinLabel(eci), "") == 0) continue;
  cout<<(fIntFlowExtraCorrelationsPro->GetXaxis())->GetBinLabel(eci)<<":"<<endl;
  cout<<"from Q-vectors    = "<<fIntFlowExtraCorrelationsPro->GetBinContent(eci)<<endl;
  cout<<"from nested loops = "<<fIntFlowExtraDirectCorrelations->GetBinContent(eci)<<endl;
  cout<<endl;
 }

} // end of void AliFlowAnalysisWithQCumulants::CrossCheckIntFlowExtraCorrelations()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::EvaluateIntFlowCorrectionsForNUAWithNestedLoops(AliFlowEventSimple * const anEvent)
{
 // Evaluate with nested loops correction terms for non-uniform acceptance relevant for NONAME integrated flow (to be improved (name)).
 //
 // Remark: Both sin and cos correction terms are calculated in this method. Sin terms are stored in fIntFlowDirectCorrectionTermsForNUA[0],
 // and cos terms in fIntFlowDirectCorrectionTermsForNUA[1]. Binning of fIntFlowDirectCorrectionTermsForNUA[sc] is organized as follows 
 // (sc stands for either sin or cos):
 
 //  1st bin: <<sc(n*(phi1))>> 
 //  2nd bin: <<sc(n*(phi1+phi2))>> 
 //  3rd bin: <<sc(n*(phi1-phi2-phi3))>>
 //  4th bin: <<sc(n*(2phi1-phi2))>>
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 Double_t phi1=0., phi2=0., phi3=0.;
 Int_t n = fHarmonic; 
 Int_t eventNo = (Int_t)fAvMultiplicity->GetBinEntries(1); // to be improved (is this casting safe in general?)
 Double_t dMult = (*fSpk)(0,0);
 cout<<endl;
 cout<<"Correction terms for non-uniform acceptance: Event number: "<<eventNo<<", multiplicity is "<<dMult<<endl;
 if(dMult<1)
 {
  cout<<"... skipping this event (multiplicity too low) ..."<<endl;
 } else if (dMult>fMaxAllowedMultiplicity)
   {
    cout<<"... skipping this event (multiplicity too high) ..."<<endl;
   } else 
     { 
      cout<<"... evaluating nested loops (without using particle weights)..."<<endl;
     }
 
 if(nPrim>=1 && nPrim<=fMaxAllowedMultiplicity)
 {
  // 1-particle correction terms for non-uniform acceptance:       
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   if(nPrim==1) cout<<i1<<"\r"<<flush;
   // sin terms:
   fIntFlowDirectCorrectionTermsForNUA[0]->Fill(0.5,sin(n*phi1),1.); // <sin(n*phi1)>  
   // cos terms:
   fIntFlowDirectCorrectionTermsForNUA[1]->Fill(0.5,cos(n*phi1),1.); // <cos(n*phi1)>
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=1) 
  
 if(nPrim>=2 && nPrim<=fMaxAllowedMultiplicity)
 {
  // 2-particle correction terms for non-uniform acceptance:       
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();  
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    if(nPrim==2) cout<<i1<<" "<<i2<<"\r"<<flush;
    // sin terms:
    fIntFlowDirectCorrectionTermsForNUA[0]->Fill(1.5,sin(n*(phi1+phi2)),1.); // <<sin(n*(phi1+phi2))>>
    fIntFlowDirectCorrectionTermsForNUA[0]->Fill(3.5,sin(n*(2*phi1-phi2)),1.); // <<sin(n*(2*phi1-phi2))>>
    // cos terms:
    fIntFlowDirectCorrectionTermsForNUA[1]->Fill(1.5,cos(n*(phi1+phi2)),1.); // <<cos(n*(phi1+phi2))>>
    fIntFlowDirectCorrectionTermsForNUA[1]->Fill(3.5,cos(n*(2*phi1-phi2)),1.); // <<cos(n*(2*phi1-phi2))>>
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=2)

 if(nPrim>=3 && nPrim<=fMaxAllowedMultiplicity)
 {
  // 3-particle correction terms for non-uniform acceptance:       
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     if(nPrim>=3) cout<<i1<<" "<<i2<<" "<<i3<<"\r"<<flush; // to be improved (eventually I will change this if statement)
     // sin terms:
     fIntFlowDirectCorrectionTermsForNUA[0]->Fill(2.5,sin(n*(phi1-phi2-phi3)),1.); // <<sin(n*(phi1-phi2-phi3))>>
     // cos terms:
     fIntFlowDirectCorrectionTermsForNUA[1]->Fill(2.5,cos(n*(phi1-phi2-phi3)),1.); // <<cos(n*(phi1-phi2-phi3))>>
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=3)

 cout<<endl;
}

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::EvaluateDiffFlowCorrelationsWithNestedLoops(AliFlowEventSimple * const anEvent, TString type, TString ptOrEta)
{
 // Evaluate reduced correlations with nested loops without using the particle weights.
 
 // Remark 1: Reduced correlations are evaluated in pt bin number fCrossCheckInPtBinNo and eta bin number fCrossCheckInEtaBinNo both for RPs and POIs.
 // Remark 2: Results are stored in 1 bin profiles fDiffFlowDirectCorrelations[t][pe][ci], where indices runs as follows:
 //           [0=RP,1=POI][0=Pt,1=Eta][0=<2'>,1=<4'>,2=<6'>,3=<8'>] 
 // Remark 3: <2'> = <cos(n*(psi1-phi2))>
 //           <4'> = <cos(n*(psi1+phi2-phi3-phi4))>
 // ...
 
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;
 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   }      
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
      
 Double_t lowerPtEtaEdge[2] = {fPtMin+(fCrossCheckInPtBinNo-1)*fPtBinWidth,fEtaMin+(fCrossCheckInEtaBinNo-1)*fEtaBinWidth};
 Double_t upperPtEtaEdge[2] = {fPtMin+fCrossCheckInPtBinNo*fPtBinWidth,fEtaMin+fCrossCheckInEtaBinNo*fEtaBinWidth};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 
 Double_t psi1=0., phi2=0., phi3=0., phi4=0.;// phi5=0., phi6=0., phi7=0., phi8=0.;
 
 Int_t n = fHarmonic; 
  
 // 2'-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
       
  psi1=aftsTrack->Phi(); 
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP):
   if(!(aftsTrack->InRPSelection()))continue;
   phi2=aftsTrack->Phi();   
   // 2'-particle correlations: 
   fDiffFlowDirectCorrelations[t][pe][0]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(1.*n*(psi1-phi2)),1.); // <cos(n*(psi1-phi2))  
  }//end of for(Int_t i2=0;i2<nPrim;i2++)
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 /*
 
 // 3'-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(ptOrEta == "Pt")
  { 
   if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
  } else if (ptOrEta == "Eta")
    {
     if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
    }
  psi1=aftsTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP):
   if(!(aftsTrack->InRPSelection())) continue;
   phi2=aftsTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2)continue;
    aftsTrack=anEvent->GetTrack(i3);
    // RP condition (!(first) particle in the correlator must be RP):
    if(!(aftsTrack->InRPSelection())) continue;
    phi3=aftsTrack->Phi();
    // to be improved : where to store it? ->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*(2.*phi1-phi2-phi3)),1.); // <w1 w2 w3 cos(n(2psi1-phi2-phi3))> 
   }//end of for(Int_t i3=0;i3<nPrim;i3++)  
  }//end of for(Int_t i2=0;i2<nPrim;i2++)  
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 */
 
 // 4'-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
    
  psi1=aftsTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1) continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP): 
   if(!(aftsTrack->InRPSelection())) continue;
   phi2=aftsTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2) continue;
    aftsTrack=anEvent->GetTrack(i3);
    // RP condition (!(first) particle in the correlator must be RP):
    if(!(aftsTrack->InRPSelection())) continue;
    phi3=aftsTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3) continue;
     aftsTrack=anEvent->GetTrack(i4);
     // RP condition (!(first) particle in the correlator must be RP):
     if(!(aftsTrack->InRPSelection())) continue;  
     phi4=aftsTrack->Phi();
     // 4'-particle correlations:
     fDiffFlowDirectCorrelations[t][pe][1]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*(psi1+phi2-phi3-phi4)),1.); // <cos(n(psi1+phi2-phi3-phi4))>     
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
      
 // count # of RPs and POIs in selected pt and eta bins for cross-checkings:
 for(Int_t i=0;i<nPrim;i++)
 {
  aftsTrack=anEvent->GetTrack(i); 
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  if(t==1)t++; 
  fNoOfParticlesInBin->Fill(t+pe+0.5);  
 }

} // end of void AliFlowAnalysisWithQCumulants::EvaluateDiffFlowCorrelationsWithNestedLoops(AliFlowEventSimple* anEvent, TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::EvaluateOtherDiffCorrelatorsWithNestedLoops(AliFlowEventSimple * const anEvent, TString type, TString ptOrEta)
{
 // Evaluate other differential correlators with nested loops without using the particle weights.

 // Remark 1: Other differential correlators are evaluated in pt bin number fCrossCheckInPtBinNo 
 //           and eta bin number fCrossCheckInEtaBinNo both for RPs and POIs.
 // Remark 2: Results are stored in 1 bin profiles fOtherDirectDiffCorrelators[t][pe][sc][ci], where indices runs as follows:
 //           [0=RP,1=POI][0=Pt,1=Eta][0=sin terms,1=cos terms][ci = correlator index] 
 // Remark 3: Correlator index 'ci' runs as follows:
 //            0: <exp(n*(psi1-3phi2+2phi3))> (Teaney-Yan correlator)
  
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;
 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   }      
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
      
 Double_t lowerPtEtaEdge[2] = {fPtMin+(fCrossCheckInPtBinNo-1)*fPtBinWidth,fEtaMin+(fCrossCheckInEtaBinNo-1)*fEtaBinWidth};
 Double_t upperPtEtaEdge[2] = {fPtMin+fCrossCheckInPtBinNo*fPtBinWidth,fEtaMin+fCrossCheckInEtaBinNo*fEtaBinWidth};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 
 Double_t psi1=0., phi2=0., phi3=0.;
 
 Int_t n = fHarmonic; 

 // 3-p correlators:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  psi1=aftsTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1) continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP): 
   if(!(aftsTrack->InRPSelection())) continue;
   phi2=aftsTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2) continue;
    aftsTrack=anEvent->GetTrack(i3);
    // RP condition (!(first) particle in the correlator must be RP):
    if(!(aftsTrack->InRPSelection())) continue;
    phi3=aftsTrack->Phi();
    // Fill 3-p correlators:
    fOtherDirectDiffCorrelators[t][pe][1][0]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*(psi1-3.*phi2+2.*phi3)),1.); // <cos(n(psi1-3.*phi2+2.*phi3))>     
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)   
} // end of void AliFlowAnalysisWithQCumulants::EvaluateOtherDiffCorrelatorsWithNestedLoops(AliFlowEventSimple * const anEvent, TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CrossCheckDiffFlowCorrelations(TString type, TString ptOrEta)
{
 // Compare correlations needed for diff. flow calculated with nested loops and those calculated from Q-vectors
 
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;
 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   }      
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
      
 TString rpORpoiString[2] = {"RP ","POI"}; // to be improved (name in the same way as in the other methods, eventually promote to data member) 
 TString ptORetaString[2] = {"pt","eta"}; // to be improved (name in the same way as in the other methods, eventually promote to data member) 
 TString reducedCorrelations[4] = {"<<cos(n(psi1-phi2))>>","<<cos(n(psi1+phi2-phi3-phi4))>>","",""}; // to be improved (access this from pro or hist)
 Double_t lowerPtEtaEdge[2] = {fPtMin+(fCrossCheckInPtBinNo-1)*fPtBinWidth,fEtaMin+(fCrossCheckInEtaBinNo-1)*fEtaBinWidth};
 Double_t upperPtEtaEdge[2] = {fPtMin+fCrossCheckInPtBinNo*fPtBinWidth,fEtaMin+fCrossCheckInEtaBinNo*fEtaBinWidth};
 
 Int_t crossCheckInPtEtaBinNo[2] = {fCrossCheckInPtBinNo,fCrossCheckInEtaBinNo};
 

 cout<<endl;
 cout<<"   *****************************************"<<endl;
 cout<<"   **** cross-checking the correlations ****"<<endl;
 cout<<"   ****   for differential flow ("<<rpORpoiString[t]<<")   ****"<<endl;
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
 {
  cout<<"   ****   (particle weights not used)   ****"<<endl;
 } else
   {
    cout<<"   ****    (particle weights used)      ****"<<endl;
   } 
 cout<<"   *****************************************"<<endl; 
 cout<<endl;
 cout<<"           "<<ptORetaString[pe]<<" bin: "<<lowerPtEtaEdge[pe]<<" <= "<<ptORetaString[pe]<<" < "<<upperPtEtaEdge[pe]<<endl;
 cout<<endl;
 
 for(Int_t rci=0;rci<2;rci++) // to be improved (calculate 6th and 8th order)
 {
  cout<<"      "<<reducedCorrelations[rci].Data()<<":"<<endl;
  cout<<"      from Q-vectors    = "<<fDiffFlowCorrelationsPro[t][pe][rci]->GetBinContent(crossCheckInPtEtaBinNo[pe])<<endl;
  cout<<"      from nested loops = "<<fDiffFlowDirectCorrelations[t][pe][rci]->GetBinContent(1)<<endl;
  cout<<endl;  
 } // end of for(Int_t rci=0;rci<4;rci++)
        
} // end of void AliFlowAnalysisWithQCumulants::CrossCheckDiffFlowCorrelations(TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CrossCheckOtherDiffCorrelators(TString type, TString ptOrEta)
{
 // Compare correlations needed for diff. flow calculated with nested loops and those calculated from Q-vectors
 
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;
 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   }      
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
      
 TString rpORpoiString[2] = {"RP ","POI"}; // to be improved (name in the same way as in the other methods, eventually promote to data member) 
 TString ptORetaString[2] = {"pt","eta"}; // to be improved (name in the same way as in the other methods, eventually promote to data member) 
 TString otherCorrelators[1] = {"<<cos(n(psi1-3phi2+2phi3))>>"}; // to be improved (access this from pro or hist)
 Double_t lowerPtEtaEdge[2] = {fPtMin+(fCrossCheckInPtBinNo-1)*fPtBinWidth,fEtaMin+(fCrossCheckInEtaBinNo-1)*fEtaBinWidth};
 Double_t upperPtEtaEdge[2] = {fPtMin+fCrossCheckInPtBinNo*fPtBinWidth,fEtaMin+fCrossCheckInEtaBinNo*fEtaBinWidth};
 
 Int_t crossCheckInPtEtaBinNo[2] = {fCrossCheckInPtBinNo,fCrossCheckInEtaBinNo};

 cout<<endl;
 cout<<"   *****************************************"<<endl;
 cout<<"   ****   cross-checking the other      ****"<<endl;
 cout<<"   ****   diff. correlators ("<<rpORpoiString[t]<<")       ****"<<endl;
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
 {
  cout<<"   ****   (particle weights not used)   ****"<<endl;
 } else
   {
    cout<<"   ****    (particle weights used)      ****"<<endl;
   } 
 cout<<"   *****************************************"<<endl; 
 cout<<endl;
 cout<<"           "<<ptORetaString[pe]<<" bin: "<<lowerPtEtaEdge[pe]<<" <= "<<ptORetaString[pe]<<" < "<<upperPtEtaEdge[pe]<<endl;
 cout<<endl;
 
 for(Int_t ci=0;ci<1;ci++) 
 {
  cout<<"      "<<otherCorrelators[ci].Data()<<":"<<endl;
  cout<<"      from Q-vectors    = "<<fOtherDiffCorrelators[t][pe][1][ci]->GetBinContent(crossCheckInPtEtaBinNo[pe])<<endl;
  cout<<"      from nested loops = "<<fOtherDirectDiffCorrelators[t][pe][1][ci]->GetBinContent(1)<<endl;
  cout<<endl;  
 } // end of for(Int_t ci=0;ci<1;ci++)
        
} // end of void AliFlowAnalysisWithQCumulants::CrossCheckOtherDiffCorrelators(TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::PrintNumberOfParticlesInSelectedBin()
{
 // Print on the screen number of RPs and POIs in selected pt and eta bin for cross checkings.
 
 cout<<endl;
 cout<<"Number of RPs in selected pt bin   = "<<fNoOfParticlesInBin->GetBinContent(1)<<endl;
 cout<<"Number of RPs in selected eta bin  = "<<fNoOfParticlesInBin->GetBinContent(2)<<endl;
 cout<<"Number of POIs in selected pt bin  = "<<fNoOfParticlesInBin->GetBinContent(3)<<endl;
 cout<<"Number of POIs in selected eta bin = "<<fNoOfParticlesInBin->GetBinContent(4)<<endl;
 
} // end of void AliFlowAnalysisWithQCumulants::PrintNumberOfParticlesInSelectedBin()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::EvaluateDiffFlowCorrelationsWithNestedLoopsUsingParticleWeights(AliFlowEventSimple * const anEvent, TString type, TString ptOrEta)
{
 // Evaluate reduced correlations with nested loops without using the particle weights.
 
 // Remark 1: Reduced correlations are evaluated in pt bin number fCrossCheckInPtBinNo and eta bin number fCrossCheckInEtaBinNo both for RPs and POIs.
 // Remark 2: Results are stored in 1 bin profiles fDiffFlowDirectCorrelations[t][pe][ci], where indices runs as follows:
 //           [0=RP,1=POI][0=Pt,1=Eta][0=<2'>,1=<4'>,2=<6'>,3=<8'>] 
 // Remark 3: <2'> = <w2 cos(n*(psi1-phi2))>
 //           <4'> = <w2 w3 w4 cos(n*(psi1+phi2-phi3-phi4))>
 // ...
  
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;
 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   }      
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
      
 Double_t lowerPtEtaEdge[2] = {fPtMin+(fCrossCheckInPtBinNo-1)*fPtBinWidth,fEtaMin+(fCrossCheckInEtaBinNo-1)*fEtaBinWidth};
 Double_t upperPtEtaEdge[2] = {fPtMin+fCrossCheckInPtBinNo*fPtBinWidth,fEtaMin+fCrossCheckInEtaBinNo*fEtaBinWidth};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 
 Double_t psi1=0., phi2=0., phi3=0., phi4=0.;// phi5=0., phi6=0., phi7=0., phi8=0.;
 Double_t wPhi2=1., wPhi3=1., wPhi4=1.;// wPhi5=1., wPhi6=1., wPhi7=1., wPhi8=1.;
 
 Int_t n = fHarmonic; 
 
 // 2'-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  psi1=aftsTrack->Phi(); 
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1) continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP):
   if(!(aftsTrack->InRPSelection())) continue;
   phi2=aftsTrack->Phi();   
   if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));
   // 2'-particle correlations: 
   fDiffFlowDirectCorrelations[t][pe][0]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(1.*n*(psi1-phi2)),wPhi2); // <w2 cos(n*(psi1-phi2))  
  }//end of for(Int_t i2=0;i2<nPrim;i2++)
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 // 4'-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  psi1=aftsTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1) continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP): 
   if(!(aftsTrack->InRPSelection())) continue;
   phi2=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2) continue;
    aftsTrack=anEvent->GetTrack(i3);
    // RP condition (!(first) particle in the correlator must be RP):
    if(!(aftsTrack->InRPSelection())) continue;
    phi3=aftsTrack->Phi();
    if(fUsePhiWeights && fPhiWeights) wPhi3 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*fnBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3) continue;
     aftsTrack=anEvent->GetTrack(i4);
     // RP condition (!(first) particle in the correlator must be RP):
     if(!(aftsTrack->InRPSelection())) continue;  
     phi4=aftsTrack->Phi();
     if(fUsePhiWeights && fPhiWeights) wPhi4 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*fnBinsPhi/TMath::TwoPi())));
     // 4'-particle correlations <w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))>:
     fDiffFlowDirectCorrelations[t][pe][1]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*(psi1+phi2-phi3-phi4)),wPhi2*wPhi3*wPhi4); 
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)      
 
 // count # of RPs and POIs in selected pt and eta bins for cross-checkings: (to be improved - moved to dedicated method)
 for(Int_t i=0;i<nPrim;i++)
 {
  aftsTrack=anEvent->GetTrack(i); 
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  if(t==1)t++; 
  fNoOfParticlesInBin->Fill(t+pe+0.5);  
 }
 
} // end of void AliFlowAnalysisWithQCumulants::EvaluateDiffFlowCorrelationsWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* anEvent, TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(AliFlowEventSimple * const anEvent, TString type, TString ptOrEta)
{
 // Evaluate with nested loops correction terms for non-uniform acceptance (both sin and cos terms) relevant for differential flow.
 
 // Remark 1: Reduced correction terms for non-uniform acceptance are evaluated in pt bin number fCrossCheckInPtBinNo 
 //           and eta bin number fCrossCheckInEtaBinNo both for RPs and POIs.
 // Remark 2: Results are stored in 1 bin profiles fDiffFlowDirectCorrections[t][pe][sc][cti], where first three indices runs as: 
 //           [0=RP,1=POI][0=Pt,1=Eta][0=sin terms,1=cos terms], whilst the cti (correction term index) runs as follows: 
 //  cti: 
 //    0: <<sc n(psi1)>>
 //    1: <<sc n(psi1+phi2)>> 
 //    2: <<sc n(psi1+phi2-phi3)>>
 //    3: <<sc n(psi1-phi2-phi3)>>
 //    4:
 //    5:
 //    6:
  
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;
 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   }      
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
      
 Double_t lowerPtEtaEdge[2] = {fPtMin+(fCrossCheckInPtBinNo-1)*fPtBinWidth,fEtaMin+(fCrossCheckInEtaBinNo-1)*fEtaBinWidth};
 Double_t upperPtEtaEdge[2] = {fPtMin+fCrossCheckInPtBinNo*fPtBinWidth,fEtaMin+fCrossCheckInEtaBinNo*fEtaBinWidth};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 
 Double_t psi1=0., phi2=0., phi3=0.;// phi4=0.;// phi5=0., phi6=0., phi7=0., phi8=0.;
 
 Int_t n = fHarmonic; 
 
 // 1-particle correction terms:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  psi1=aftsTrack->Phi(); 
  // sin terms: 
  fDiffFlowDirectCorrectionTermsForNUA[t][pe][0][0]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,sin(n*psi1),1.); // <<sin(n*(psi1))>>  
  // cos terms: 
  fDiffFlowDirectCorrectionTermsForNUA[t][pe][1][0]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*psi1),1.); // <<cos(n*(psi1))>>  
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
   
 // 2-particle correction terms:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
   // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  psi1=aftsTrack->Phi(); 
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1) continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP):
   if(!(aftsTrack->InRPSelection())) continue;
   phi2=aftsTrack->Phi();   
   // sin terms: 
   fDiffFlowDirectCorrectionTermsForNUA[t][pe][0][1]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,sin(n*(psi1+phi2)),1.); // <<sin(n*(psi1+phi2))>>  
   // cos terms: 
   fDiffFlowDirectCorrectionTermsForNUA[t][pe][1][1]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*(psi1+phi2)),1.); // <<cos(n*(psi1+phi2))>>  
  }//end of for(Int_t i2=0;i2<nPrim;i2++)
 }//end of for(Int_t i1=0;i1<nPrim;i1++)   
 
 // 3-particle correction terms:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
   // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  psi1=aftsTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1) continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP):
   if(!(aftsTrack->InRPSelection())) continue;
   phi2=aftsTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2) continue;
    aftsTrack=anEvent->GetTrack(i3);
    // RP condition (!(first) particle in the correlator must be RP):
    if(!(aftsTrack->InRPSelection())) continue;
    phi3=aftsTrack->Phi();
    // sin terms: 
    fDiffFlowDirectCorrectionTermsForNUA[t][pe][0][2]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,sin(n*(psi1+phi2-phi3)),1.); // <<sin(n*(psi1+phi2-phi3))>>  
    fDiffFlowDirectCorrectionTermsForNUA[t][pe][0][3]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,sin(n*(psi1-phi2-phi3)),1.); // <<sin(n*(psi1-phi2-phi3))>>  
    // cos terms: 
    fDiffFlowDirectCorrectionTermsForNUA[t][pe][1][2]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*(psi1+phi2-phi3)),1.); // <<cos(n*(psi1+phi2-phi3))>>  
    fDiffFlowDirectCorrectionTermsForNUA[t][pe][1][3]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*(psi1-phi2-phi3)),1.); // <<cos(n*(psi1-phi2-phi3))>>  
   }//end of for(Int_t i3=0;i3<nPrim;i3++)  
  }//end of for(Int_t i2=0;i2<nPrim;i2++)  
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
   
} // end of void AliFlowAnalysisWithQCumulants::EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(AliFlowEventSimple* anEvent, TString type, TString ptOrEta)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CrossCheckDiffFlowCorrectionTermsForNUA(TString type, TString ptOrEta)
{
 // Compare corrections temrs for non-uniform acceptance needed for diff. flow calculated with nested loops and those calculated from Q-vectors
 
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;
 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   }      
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
      
 TString rpORpoiString[2] = {"RP ","POI"}; // to be improved (name in the same way as in the other methods, eventually promote to data member) 
 TString ptORetaString[2] = {"pt","eta"}; // to be improved (name in the same way as in the other methods, eventually promote to data member) 
 //TString sinCosFlag[2] = {"sin","cos"}; // to be improved (eventually promote to data member)
 TString reducedCorrectionSinTerms[4] = {"<<sin(n(psi1))>>","<<sin(n(psi1+phi2))>>","<<sin(n*(psi1+phi2-phi3))>>","<<sin(n*(psi1-phi2-phi3))>>"}; // to be improved (access this from pro or hist)
 TString reducedCorrectionCosTerms[4] = {"<<cos(n(psi1))>>","<<cos(n(psi1+phi2))>>","<<cos(n*(psi1+phi2-phi3))>>","<<cos(n*(psi1-phi2-phi3))>>"}; // to be improved (access this from pro or hist)
 Double_t lowerPtEtaEdge[2] = {fPtMin+(fCrossCheckInPtBinNo-1)*fPtBinWidth,fEtaMin+(fCrossCheckInEtaBinNo-1)*fEtaBinWidth};
 Double_t upperPtEtaEdge[2] = {fPtMin+fCrossCheckInPtBinNo*fPtBinWidth,fEtaMin+fCrossCheckInEtaBinNo*fEtaBinWidth};
 
 Int_t crossCheckInPtEtaBinNo[2] = {fCrossCheckInPtBinNo,fCrossCheckInEtaBinNo};
 
 cout<<endl;
 cout<<"   ******************************************"<<endl;
 cout<<"   ****  cross-checking the correction   ****"<<endl;
 cout<<"   **** terms for non-uniform acceptance ****"<<endl; 
 cout<<"   ****    for differential flow ("<<rpORpoiString[t]<<")   ****"<<endl;
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights))
 {
  cout<<"   ****    (particle weights not used)   ****"<<endl;
 } else
   {
    cout<<"   ****     (particle weights used)      ****"<<endl;
   } 
 cout<<"   ******************************************"<<endl; 
 cout<<endl;
 cout<<"           "<<ptORetaString[pe]<<" bin: "<<lowerPtEtaEdge[pe]<<" <= "<<ptORetaString[pe]<<" < "<<upperPtEtaEdge[pe]<<endl;
 cout<<endl;
 
 for(Int_t cti=0;cti<4;cti++) // correction term index
 {
  for(Int_t sc=0;sc<2;sc++) // sin or cos terms
  {
   if(sc==0) // to be improved (this can be implemented better)
   { 
    cout<<"      "<<reducedCorrectionSinTerms[cti].Data()<<":"<<endl;
   } else
     {
      cout<<"      "<<reducedCorrectionCosTerms[cti].Data()<<":"<<endl;     
     }
   cout<<"      from Q-vectors    = "<<fDiffFlowCorrectionTermsForNUAPro[t][pe][sc][cti]->GetBinContent(crossCheckInPtEtaBinNo[pe])<<endl;
   cout<<"      from nested loops = "<<fDiffFlowDirectCorrectionTermsForNUA[t][pe][sc][cti]->GetBinContent(1)<<endl;
   cout<<endl;  
  } 
 } // end of for(Int_t rci=0;rci<4;rci++)

} // end of void AliFlowAnalysisWithQCumulants::CrossCheckDiffFlowCorrectionTermsForNUA(TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUACosTermsUsingParticleWeights()
{
 // Calculate corrections using particle weights for non-uniform acceptance of the detector for no-name integrated flow (cos terms).
 
 //                                  **********************************************************************
 //                                  **** weighted corrections for non-uniform acceptance (cos terms): ****
 //                                  **********************************************************************
 
 // Remark 1: When particle weights are used the binning of fIntFlowCorrectionTermsForNUAPro[1] is organized as follows:
 //
 // 1st bin: <<w1 cos(n*(phi1))>> = cosP1nW1
 // 2nd bin: <<w1 w2 cos(n*(phi1+phi2))>> = cosP1nP1nW1W1
 // 3rd bin: <<w1 w2 w3 cos(n*(phi1-phi2-phi3))>> = cosP1nM1nM1nW1W1W1 
 // ...

 // multiplicity (number of particles used to determine the reaction plane)
 Double_t dMult = (*fSpk)(0,0);
 
 // real and imaginary parts of weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n1k = (*fReQ)(0,1);
 Double_t dReQ2n2k = (*fReQ)(1,2);
 //Double_t dReQ3n3k = (*fReQ)(2,3);
 //Double_t dReQ4n4k = (*fReQ)(3,4);
 Double_t dReQ1n3k = (*fReQ)(0,3);
 Double_t dImQ1n1k = (*fImQ)(0,1);
 Double_t dImQ2n2k = (*fImQ)(1,2);
 //Double_t dImQ3n3k = (*fImQ)(2,3);
 //Double_t dImQ4n4k = (*fImQ)(3,4);
 //Double_t dImQ1n3k = (*fImQ)(0,3);

 // dMs are variables introduced in order to simplify some Eqs. bellow:
 //..............................................................................................
 Double_t dM11 = (*fSpk)(1,1)-(*fSpk)(0,2); // dM11 = sum_{i,j=1,i!=j}^M w_i w_j
 Double_t dM111 = (*fSpk)(2,1)-3.*(*fSpk)(0,2)*(*fSpk)(0,1)
                + 2.*(*fSpk)(0,3); // dM111 = sum_{i,j,k=1,i!=j!=k}^M w_i w_j w_k
 //..............................................................................................
         // 1-particle:
 Double_t cosP1nW1 = 0.; // <<w1 cos(n*(phi1))>>
   
 if(dMult>0 && TMath::Abs((*fSpk)(0,1))>1.e-6)
 {
  cosP1nW1 = dReQ1n1k/(*fSpk)(0,1); 
  
  // average weighted 1-particle correction (cos terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[1]->SetBinContent(1,cosP1nW1);
  
  // final average weighted 1-particle correction (cos terms) for non-uniform acceptance for all events:
  fIntFlowCorrectionTermsForNUAPro[1]->Fill(0.5,cosP1nW1,(*fSpk)(0,1));  
 } 
 
 // 2-particle:
 Double_t cosP1nP1nW1W1 = 0.; // <<w1 w2 cos(n*(phi1+phi2))>>
 
 if(dMult>1 && TMath::Abs(dM11)>1.e-6)
 {
  cosP1nP1nW1W1 = (pow(dReQ1n1k,2)-pow(dImQ1n1k,2)-dReQ2n2k)/dM11; 
  
  // average weighted 2-particle correction (cos terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[1]->SetBinContent(2,cosP1nP1nW1W1);
  
  // final average weighted 2-particle correction (cos terms) for non-uniform acceptance for all events:
  fIntFlowCorrectionTermsForNUAPro[1]->Fill(1.5,cosP1nP1nW1W1,dM11);  
 } 
 
 // 3-particle:
 Double_t cosP1nM1nM1nW1W1W1 = 0.; // <<w1 w2 w3 cos(n*(phi1-phi2-phi3))>>
 
 if(dMult>2 && TMath::Abs(dM111)>1.e-6)
 {
  cosP1nM1nM1nW1W1W1 = (dReQ1n1k*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2))
                     - dReQ1n1k*dReQ2n2k-dImQ1n1k*dImQ2n2k
                     - 2.*((*fSpk)(0,2))*dReQ1n1k
                     + 2.*dReQ1n3k) 
                     / dM111; 
  
  // average non-weighted 3-particle correction (cos terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[1]->SetBinContent(3,cosP1nM1nM1nW1W1W1);
  
  // final average non-weighted 3-particle correction (cos terms) for non-uniform acceptance for all events:
  fIntFlowCorrectionTermsForNUAPro[1]->Fill(2.5,cosP1nM1nM1nW1W1W1,dM111);  
 } 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUACosTermsUsingParticleWeights()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUASinTermsUsingParticleWeights()
{
 // calculate corrections using particle weights for non-uniform acceptance of the detector for no-name integrated flow (sin terms)
 
 //                                  **********************************************************************
 //                                  **** weighted corrections for non-uniform acceptance (sin terms): ****
 //                                  **********************************************************************
 
 // Remark 1: When particle weights are used the binning of fIntFlowCorrectionTermsForNUAPro[0] is organized as follows:
 //
 // 1st bin: <<w1 sin(n*(phi1))>> = sinP1nW1
 // 2nd bin: <<w1 w2 sin(n*(phi1+phi2))>> = sinP1nP1nW1W1
 // 3rd bin: <<w1 w2 w3 sin(n*(phi1-phi2-phi3))>> = sinP1nM1nM1nW1W1W1 
 // ...

 // multiplicity (number of particles used to determine the reaction plane)
 Double_t dMult = (*fSpk)(0,0);
 
 // real and imaginary parts of weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n1k = (*fReQ)(0,1);
 Double_t dReQ2n2k = (*fReQ)(1,2);
 //Double_t dReQ3n3k = (*fReQ)(2,3);
 //Double_t dReQ4n4k = (*fReQ)(3,4);
 //Double_t dReQ1n3k = (*fReQ)(0,3);
 Double_t dImQ1n1k = (*fImQ)(0,1);
 Double_t dImQ2n2k = (*fImQ)(1,2);
 //Double_t dImQ3n3k = (*fImQ)(2,3);
 //Double_t dImQ4n4k = (*fImQ)(3,4);
 Double_t dImQ1n3k = (*fImQ)(0,3);

 // dMs are variables introduced in order to simplify some Eqs. bellow:
 //..............................................................................................
 Double_t dM11 = (*fSpk)(1,1)-(*fSpk)(0,2); // dM11 = sum_{i,j=1,i!=j}^M w_i w_j
 Double_t dM111 = (*fSpk)(2,1)-3.*(*fSpk)(0,2)*(*fSpk)(0,1)
                + 2.*(*fSpk)(0,3); // dM111 = sum_{i,j,k=1,i!=j!=k}^M w_i w_j w_k
 //..............................................................................................
 
 // 1-particle:
 Double_t sinP1nW1 = 0.; // <<w1 sin(n*(phi1))>>
 
 if(dMult>0 && TMath::Abs((*fSpk)(0,1))>1.e-6)
 {
  sinP1nW1 = dImQ1n1k/((*fSpk)(0,1)); 
     
  // average weighted 1-particle correction (sin terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[0]->SetBinContent(1,sinP1nW1);
  
  // final average weighted 1-particle correction (sin terms) for non-uniform acceptance for all events:   
  fIntFlowCorrectionTermsForNUAPro[0]->Fill(0.5,sinP1nW1,(*fSpk)(0,1));  
 } 
 
 // 2-particle:
 Double_t sinP1nP1nW1W1 = 0.; // <<w1 w2 sin(n*(phi1+phi2))>>
 
 if(dMult>1 && TMath::Abs(dM11)>1.e-6)
 {
  sinP1nP1nW1W1 = (2.*dReQ1n1k*dImQ1n1k-dImQ2n2k)/dM11; 
     
  // average weighted 2-particle correction (sin terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[0]->SetBinContent(2,sinP1nP1nW1W1);
  
  // final average weighted 1-particle correction (sin terms) for non-uniform acceptance for all events:      
  fIntFlowCorrectionTermsForNUAPro[0]->Fill(1.5,sinP1nP1nW1W1,dM11);  
 } 
 
 // 3-particle:
 Double_t sinP1nM1nM1nW1W1W1 = 0.; // <<w1 w2 w3 sin(n*(phi1-phi2-phi3))>>
 
 if(dMult>2 && TMath::Abs(dM111)>1.e-6)
 {
  sinP1nM1nM1nW1W1W1 = (-dImQ1n1k*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2))
                     + dReQ1n1k*dImQ2n2k-dImQ1n1k*dReQ2n2k
                     + 2.*((*fSpk)(0,2))*dImQ1n1k
                     - 2.*dImQ1n3k)
                     / dM111; 
  
  // average weighted 3-particle correction (sin terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[0]->SetBinContent(3,sinP1nM1nM1nW1W1W1);
  
  // final average weighted 3-particle correction (sin terms) for non-uniform acceptance for all events:  
  fIntFlowCorrectionTermsForNUAPro[0]->Fill(2.5,sinP1nM1nM1nW1W1W1,dM111);  
 } 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUASinTermsUsingParticleWeights()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::EvaluateIntFlowCorrectionsForNUAWithNestedLoopsUsingParticleWeights(AliFlowEventSimple * const anEvent)
{
 // Evaluate with nested loops correction terms for non-uniform acceptance for integrated flow (using the particle weights). 

 // Results are stored in profiles fIntFlowDirectCorrectionTermsForNUA[0] (sin terms) and
 // fIntFlowDirectCorrectionTermsForNUA[1] (cos terms). 
 
 // Remark 1: When particle weights are used the binning of fIntFlowDirectCorrectionTermsForNUA[sc] is 
 // organized as follows (sc stands for either sin or cos):
 //
 // 1st bin: <<w1 sc(n*(phi1))>> = scP1nW1
 // 2nd bin: <<w1 w2 sc(n*(phi1+phi2))>> = scP1nP1nW1W1
 // 3rd bin: <<w1 w2 w3 sc(n*(phi1-phi2-phi3))>> = scP1nM1nM1nW1W1W1 
 // ...
  
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 //Double_t phi1=0., phi2=0., phi3=0., phi4=0., phi5=0., phi6=0., phi7=0., phi8=0.;
 //Double_t wPhi1=1., wPhi2=1., wPhi3=1., wPhi4=1., wPhi5=1., wPhi6=1., wPhi7=1., wPhi8=1.;
 Double_t phi1=0., phi2=0., phi3=0.;
 Double_t wPhi1=1., wPhi2=1., wPhi3=1.;
 Int_t n = fHarmonic; 
 Int_t eventNo = (Int_t)fAvMultiplicity->GetBinEntries(1); // to be improved (is this casting safe in general?)
 Double_t dMult = (*fSpk)(0,0);
 cout<<endl;
 cout<<"Correction terms for non-uniform acceptance: Event number: "<<eventNo<<", multiplicity is "<<dMult<<endl;
 if(dMult<1)
 {
  cout<<"... skipping this event (multiplicity too low) ..."<<endl;
 } else if (dMult>fMaxAllowedMultiplicity)
   {
    cout<<"... skipping this event (multiplicity too high) ..."<<endl;
   } else 
     { 
      cout<<"... evaluating nested loops (using particle weights) ..."<<endl;
     } 
      
 // 1-particle correction terms using particle weights:       
 if(nPrim>=1 && nPrim<=fMaxAllowedMultiplicity)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*fnBinsPhi/TMath::TwoPi())));
   // 1-particle correction terms using particle weights:
   if(fUsePhiWeights) fIntFlowDirectCorrectionTermsForNUA[0]->Fill(0.5,sin(n*phi1),wPhi1); // <w1 sin(n*phi1)>
   if(fUsePhiWeights) fIntFlowDirectCorrectionTermsForNUA[1]->Fill(0.5,cos(n*phi1),wPhi1); // <w1 cos(n*phi1)>
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=1 && nPrim<=fMaxAllowedMultiplicity) 
 
 // 2-particle correction terms using particle weights:       
 if(nPrim>=2 && nPrim<=fMaxAllowedMultiplicity)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*fnBinsPhi/TMath::TwoPi())));
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));   
    if(nPrim==2) cout<<i1<<" "<<i2<<"\r"<<flush;
    // 2-p correction terms using particle weights:    
    if(fUsePhiWeights) fIntFlowDirectCorrectionTermsForNUA[0]->Fill(1.5,sin(n*(phi1+phi2)),wPhi1*wPhi2); // <w1 w2 sin(n*(phi1+phi2))>
    if(fUsePhiWeights) fIntFlowDirectCorrectionTermsForNUA[1]->Fill(1.5,cos(n*(phi1+phi2)),wPhi1*wPhi2); // <w1 w2 cos(n*(phi1+phi2))>
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=2)

 // 3-particle correction terms using particle weights:       
 if(nPrim>=3 && nPrim<=fMaxAllowedMultiplicity)
 { 
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*fnBinsPhi/TMath::TwoPi())));
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     if(fUsePhiWeights && fPhiWeights) wPhi3 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*fnBinsPhi/TMath::TwoPi())));
     if(nPrim==3) cout<<i1<<" "<<i2<<" "<<i3<<"\r"<<flush;
     // 3-p correction terms using particle weights:    
     if(fUsePhiWeights) fIntFlowDirectCorrectionTermsForNUA[0]->Fill(2.5,sin(n*(phi1-phi2-phi3)),wPhi1*wPhi2*wPhi3); // <w1 w2 w3 sin(n*(phi1-phi2-phi3))>
     if(fUsePhiWeights) fIntFlowDirectCorrectionTermsForNUA[1]->Fill(2.5,cos(n*(phi1-phi2-phi3)),wPhi1*wPhi2*wPhi3); // <w1 w2 w3 cos(n*(phi1-phi2-phi3))>
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=3)
 
 /*
 
 if(nPrim>=4 && nPrim<=fMaxAllowedMultiplicity)
 {
  // 4 nested loops multiparticle correlations using particle weights:       
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())) continue;
   phi1=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*fnBinsPhi/TMath::TwoPi())));
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1)continue;
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())) continue;
    phi2=aftsTrack->Phi();
    if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2)continue;
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())) continue;
     phi3=aftsTrack->Phi();
     if(fUsePhiWeights && fPhiWeights) wPhi3 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*fnBinsPhi/TMath::TwoPi())));
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3)continue;
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())) continue;
      phi4=aftsTrack->Phi();
      if(fUsePhiWeights && fPhiWeights) wPhi4 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*fnBinsPhi/TMath::TwoPi())));
      if(nPrim>=4) cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<"\r"<<flush; // to be improved (replace eventually this if statement with if(nPrim==4))
      // 4-p correlations using particle weights:
      if(fUsePhiWeights) fIntFlowDirectCorrelations->Fill(10.5,cos(n*phi1+n*phi2-n*phi3-n*phi4),wPhi1*wPhi2*wPhi3*wPhi4); 
      // extra correlations: 
      // 2-p extra correlations (do not appear if particle weights are not used):
      // ...
      // 3-p extra correlations (do not appear if particle weights are not used):
      // ...
      // 4-p extra correlations (do not appear if particle weights are not used):
      // ...
     } // end of for(Int_t i4=0;i4<nPrim;i4++) 
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=4)

 */

 cout<<endl; 

} // end of void AliFlowAnalysisWithQCumulants::EvaluateIntFlowCorrectionsForNUAWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* anEvent)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights(TString type, TString ptOrEta)
{
 // Calculate correction terms for non-uniform acceptance for differential flow (cos terms) using particle weights.
 
 // Results are stored in fDiffFlowCorrectionTermsForNUAPro[t][pe][1][cti], where cti runs as follows:
 //
 //  0: <<cos n(psi)>>
 //  1: <<w2 cos n(psi1+phi2)>>
 //  2: <<w2 w3 cos n(psi1+phi2-phi3)>>
 //  3: <<w2 w3 cos n(psi1-phi2-phi3)>>
 //  4:
 //  5:
 //  6:
 
 // real and imaginary parts of weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n1k = (*fReQ)(0,1);
 Double_t dReQ2n2k = (*fReQ)(1,2);
 //Double_t dReQ1n3k = (*fReQ)(0,3);
 //Double_t dReQ4n4k = (*fReQ)(3,4);
 Double_t dImQ1n1k = (*fImQ)(0,1);
 Double_t dImQ2n2k = (*fImQ)(1,2);
 //Double_t dImQ1n3k = (*fImQ)(0,3);
 //Double_t dImQ4n4k = (*fImQ)(3,4);
 
 // S^M_{p,k} (see .h file for the definition of fSpk):
 Double_t dSM1p1k = (*fSpk)(0,1);
 Double_t dSM1p2k = (*fSpk)(0,2);
 Double_t dSM2p1k = (*fSpk)(1,1);

 Int_t t = 0; // type flag 
 Int_t pe = 0; // ptEta flag
 
 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }

 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   }
    
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};
 
 // looping over all bins and calculating correction terms: 
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // real and imaginary parts of p_{m*n,0} (non-weighted Q-vector evaluated for POIs in particular pt or eta bin): 
  Double_t p1n0kRe = 0.;
  Double_t p1n0kIm = 0.;

  // number of POIs in particular pt or eta bin:
  Double_t mp = 0.;

  // real and imaginary parts of q_{m*n,0} (weighted Q-vector evaluated for particles which are both RPs and POIs in particular pt or eta bin):
  Double_t q1n2kRe = 0.;
  Double_t q1n2kIm = 0.;
  Double_t q2n1kRe = 0.;
  Double_t q2n1kIm = 0.;
    
  // s_{1,1}, s_{1,2} // to be improved (add explanation)  
  Double_t s1p1k = 0.; 
  Double_t s1p2k = 0.; 
  
  // number of particles which are both RPs and POIs in particular pt or eta bin:
  Double_t mq = 0.;
  
  // M0111 from Eq. (118) in QC2c (to be improved (notation))
  Double_t dM01 = 0.;
  Double_t dM011 = 0.;
  
  if(type == "POI")
  {           
   // q_{m*n,k}:
   q1n2kRe = fReRPQ1dEBE[2][pe][0][2]->GetBinContent(fReRPQ1dEBE[2][pe][0][2]->GetBin(b))
           * fReRPQ1dEBE[2][pe][0][2]->GetBinEntries(fReRPQ1dEBE[2][pe][0][2]->GetBin(b));
   q1n2kIm = fImRPQ1dEBE[2][pe][0][2]->GetBinContent(fImRPQ1dEBE[2][pe][0][2]->GetBin(b))
           * fImRPQ1dEBE[2][pe][0][2]->GetBinEntries(fImRPQ1dEBE[2][pe][0][2]->GetBin(b));         
   q2n1kRe = fReRPQ1dEBE[2][pe][1][1]->GetBinContent(fReRPQ1dEBE[2][pe][1][1]->GetBin(b))
           * fReRPQ1dEBE[2][pe][1][1]->GetBinEntries(fReRPQ1dEBE[2][pe][1][1]->GetBin(b));
   q2n1kIm = fImRPQ1dEBE[2][pe][1][1]->GetBinContent(fImRPQ1dEBE[2][pe][1][1]->GetBin(b))
           * fImRPQ1dEBE[2][pe][1][1]->GetBinEntries(fImRPQ1dEBE[2][pe][1][1]->GetBin(b));         
   mq = fReRPQ1dEBE[2][pe][1][1]->GetBinEntries(fReRPQ1dEBE[2][pe][1][1]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
   
   s1p1k = pow(fs1dEBE[2][pe][1]->GetBinContent(b)*fs1dEBE[2][pe][1]->GetBinEntries(b),1.); 
   s1p2k = pow(fs1dEBE[2][pe][2]->GetBinContent(b)*fs1dEBE[2][pe][2]->GetBinEntries(b),1.); 
  }else if(type == "RP")
   {
    // q_{m*n,k}: (Remark: m=1 is 0, k=0 iz zero (to be improved!)) 
    q1n2kRe = fReRPQ1dEBE[0][pe][0][2]->GetBinContent(fReRPQ1dEBE[0][pe][0][2]->GetBin(b))
            * fReRPQ1dEBE[0][pe][0][2]->GetBinEntries(fReRPQ1dEBE[0][pe][0][2]->GetBin(b));
    q1n2kIm = fImRPQ1dEBE[0][pe][0][2]->GetBinContent(fImRPQ1dEBE[0][pe][0][2]->GetBin(b))
            * fImRPQ1dEBE[0][pe][0][2]->GetBinEntries(fImRPQ1dEBE[0][pe][0][2]->GetBin(b));
    q2n1kRe = fReRPQ1dEBE[0][pe][1][1]->GetBinContent(fReRPQ1dEBE[0][pe][1][1]->GetBin(b))
            * fReRPQ1dEBE[0][pe][1][1]->GetBinEntries(fReRPQ1dEBE[0][pe][1][1]->GetBin(b));
    q2n1kIm = fImRPQ1dEBE[0][pe][1][1]->GetBinContent(fImRPQ1dEBE[0][pe][1][1]->GetBin(b))
            * fImRPQ1dEBE[0][pe][1][1]->GetBinEntries(fImRPQ1dEBE[0][pe][1][1]->GetBin(b));
    // s_{1,1}, s_{1,2} and s_{1,3} // to be improved (add explanation)  
    s1p1k = pow(fs1dEBE[0][pe][1]->GetBinContent(b)*fs1dEBE[0][pe][1]->GetBinEntries(b),1.); 
    s1p2k = pow(fs1dEBE[0][pe][2]->GetBinContent(b)*fs1dEBE[0][pe][2]->GetBinEntries(b),1.); 
    //s1p3k = pow(fs1dEBE[0][pe][3]->GetBinContent(b)*fs1dEBE[0][pe][3]->GetBinEntries(b),1.);  
    
    mq = fReRPQ1dEBE[0][pe][1][1]->GetBinEntries(fReRPQ1dEBE[0][pe][1][1]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here) 
  }    
  
  if(type == "POI")
  {
   // p_{m*n,k}:   
   p1n0kRe = fReRPQ1dEBE[1][pe][0][0]->GetBinContent(fReRPQ1dEBE[1][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b));
   p1n0kIm = fImRPQ1dEBE[1][pe][0][0]->GetBinContent(fImRPQ1dEBE[1][pe][0][0]->GetBin(b))  
           * fImRPQ1dEBE[1][pe][0][0]->GetBinEntries(fImRPQ1dEBE[1][pe][0][0]->GetBin(b));
   mp = fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here) 
   // M01 from Eq. (118) in QC2c (to be improved (notation)):
   dM01 = mp*dSM1p1k-s1p1k;
   dM011 = mp*(dSM2p1k-dSM1p2k)
         - 2.*(s1p1k*dSM1p1k-s1p2k);
       
   // typeFlag = RP (0) or POI (1):   
   t = 1; 
  } else if(type == "RP")
    {  
     // to be improved (cross-checked):
     p1n0kRe = fReRPQ1dEBE[0][pe][0][0]->GetBinContent(fReRPQ1dEBE[0][pe][0][0]->GetBin(b))
             * fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b));
     p1n0kIm = fImRPQ1dEBE[0][pe][0][0]->GetBinContent(fImRPQ1dEBE[0][pe][0][0]->GetBin(b))  
             * fImRPQ1dEBE[0][pe][0][0]->GetBinEntries(fImRPQ1dEBE[0][pe][0][0]->GetBin(b));
     mp = fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
     // M01 from Eq. (118) in QC2c (to be improved (notation)):
     dM01 = mp*dSM1p1k-s1p1k;
     dM011 = mp*(dSM2p1k-dSM1p2k)
           - 2.*(s1p1k*dSM1p1k-s1p2k); 
     // typeFlag = RP (0) or POI (1): 
     t = 0;
    }
  
  // <<cos n(psi1)>>:
  Double_t cosP1nPsi = 0.;
  if(mp)
  {
   cosP1nPsi = p1n0kRe/mp;
   
   // fill profile for <<cos n(psi1)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][1][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],cosP1nPsi,mp);
   // histogram to store <cos n(psi1)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][1][0]->SetBinContent(b,cosP1nPsi);
  } // end of if(mp)   
  
  // <<w2 cos n(psi1+phi2)>>:
  Double_t cosP1nPsiP1nPhiW2 = 0.;
  if(dM01)
  {
   cosP1nPsiP1nPhiW2 = (p1n0kRe*dReQ1n1k-p1n0kIm*dImQ1n1k-q2n1kRe)/(dM01);
   // fill profile for <<w2 cos n(psi1+phi2)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][1][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],cosP1nPsiP1nPhiW2,dM01);
   // histogram to store <w2 cos n(psi1+phi2)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][1][1]->SetBinContent(b,cosP1nPsiP1nPhiW2);
  } // end of if(dM01)   
  
  // <<w2 w3 cos n(psi1+phi2-phi3)>>:
  Double_t cosP1nPsi1P1nPhi2MPhi3W2W3 = 0.;
  if(dM011)
  {
   cosP1nPsi1P1nPhi2MPhi3W2W3 = (p1n0kRe*(pow(dImQ1n1k,2.)+pow(dReQ1n1k,2.))
                              - p1n0kRe*dSM1p2k
                              - q2n1kRe*dReQ1n1k-q2n1kIm*dImQ1n1k
                              - s1p1k*dReQ1n1k
                              + 2.*q1n2kRe)
                              / dM011;  
   // fill profile for <<w1 w2 w3 cos n(psi1+phi2)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][1][2]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],cosP1nPsi1P1nPhi2MPhi3W2W3,dM011);
   // histogram to store <w1 w2 w3 cos n(psi1+phi2)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][1][2]->SetBinContent(b,cosP1nPsi1P1nPhi2MPhi3W2W3);
  } // end of if(dM011)   
  
  // <<w2 w3 cos n(psi1-phi2-phi3)>>:
  Double_t cosP1nPsi1M1nPhi2MPhi3W2W3 = 0.;
  if(dM011)
  {
   cosP1nPsi1M1nPhi2MPhi3W2W3 = (p1n0kRe*(pow(dReQ1n1k,2.)-pow(dImQ1n1k,2.))+2.*p1n0kIm*dReQ1n1k*dImQ1n1k
                              - 1.*(p1n0kRe*dReQ2n2k+p1n0kIm*dImQ2n2k)  
                              - 2.*s1p1k*dReQ1n1k
                              + 2.*q1n2kRe)
                              / dM011;
   // fill profile for <<w1 w2 w3 cos n(psi1+phi2)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][1][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],cosP1nPsi1M1nPhi2MPhi3W2W3,dM011);
   // histogram to store <w1 w2 w3 cos n(psi1+phi2)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][1][3]->SetBinContent(b,cosP1nPsi1M1nPhi2MPhi3W2W3);
  } // end of if(dM011)   
 
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)
   
} // end of AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights(TString type, TString ptOrEta)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights(TString type, TString ptOrEta)
{
 // Calculate correction terms for non-uniform acceptance for differential flow (sin terms).
  
 // Results are stored in fDiffFlowCorrectionTermsForNUAPro[t][pe][0][cti], where cti runs as follows:
 //  0: <<sin n(psi1)>>
 //  1: <<w2 sin n(psi1+phi2)>>
 //  2: <<w2 w3 sin n(psi1+phi2-phi3)>>
 //  3: <<w2 w3 sin n(psi1-phi2-phi3)>>:
 //  4:
 //  5:
 //  6:
 
 // real and imaginary parts of weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n1k = (*fReQ)(0,1);
 Double_t dReQ2n2k = (*fReQ)(1,2);
 //Double_t dReQ1n3k = (*fReQ)(0,3);
 //Double_t dReQ4n4k = (*fReQ)(3,4);
 Double_t dImQ1n1k = (*fImQ)(0,1);
 Double_t dImQ2n2k = (*fImQ)(1,2);
 //Double_t dImQ1n3k = (*fImQ)(0,3);
 //Double_t dImQ4n4k = (*fImQ)(3,4);
 
 // S^M_{p,k} (see .h file for the definition of fSpk):
 Double_t dSM1p1k = (*fSpk)(0,1);
 Double_t dSM1p2k = (*fSpk)(0,2);
 Double_t dSM2p1k = (*fSpk)(1,1);

 Int_t t = 0; // type flag 
 Int_t pe = 0; // ptEta flag
 
 if(type == "RP")
 {
  t = 0;
 } else if(type == "POI")
   {
    t = 1;
   }

 if(ptOrEta == "Pt")
 {
  pe = 0;
 } else if(ptOrEta == "Eta")
   {
    pe = 1;
   }
    
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 //Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};

 // looping over all bins and calculating correction terms: 
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // real and imaginary parts of p_{m*n,0} (non-weighted Q-vector evaluated for POIs in particular pt or eta bin): 
  Double_t p1n0kRe = 0.;
  Double_t p1n0kIm = 0.;

  // number of POIs in particular pt or eta bin:
  Double_t mp = 0.;

  // real and imaginary parts of q_{m*n,0} (weighted Q-vector evaluated for particles which are both RPs and POIs in particular pt or eta bin):
  Double_t q1n2kRe = 0.;
  Double_t q1n2kIm = 0.;
  Double_t q2n1kRe = 0.;
  Double_t q2n1kIm = 0.;
    
  // s_{1,1}, s_{1,2} and s_{1,3} // to be improved (add explanation)  
  Double_t s1p1k = 0.; 
  Double_t s1p2k = 0.; 
  
  // number of particles which are both RPs and POIs in particular pt or eta bin:
  Double_t mq = 0.;
  
  // M0111 from Eq. (118) in QC2c (to be improved (notation))
  Double_t dM01 = 0.;
  Double_t dM011 = 0.;

  if(type == "POI")
  {    
   // q_{m*n,k}:
   q1n2kRe = fReRPQ1dEBE[2][pe][0][2]->GetBinContent(fReRPQ1dEBE[2][pe][0][2]->GetBin(b))
           * fReRPQ1dEBE[2][pe][0][2]->GetBinEntries(fReRPQ1dEBE[2][pe][0][2]->GetBin(b));
   q1n2kIm = fImRPQ1dEBE[2][pe][0][2]->GetBinContent(fImRPQ1dEBE[2][pe][0][2]->GetBin(b))
           * fImRPQ1dEBE[2][pe][0][2]->GetBinEntries(fImRPQ1dEBE[2][pe][0][2]->GetBin(b));         
   q2n1kRe = fReRPQ1dEBE[2][pe][1][1]->GetBinContent(fReRPQ1dEBE[2][pe][1][1]->GetBin(b))
           * fReRPQ1dEBE[2][pe][1][1]->GetBinEntries(fReRPQ1dEBE[2][pe][1][1]->GetBin(b));
   q2n1kIm = fImRPQ1dEBE[2][pe][1][1]->GetBinContent(fImRPQ1dEBE[2][pe][1][1]->GetBin(b))
           * fImRPQ1dEBE[2][pe][1][1]->GetBinEntries(fImRPQ1dEBE[2][pe][1][1]->GetBin(b));         
   mq = fReRPQ1dEBE[2][pe][0][0]->GetBinEntries(fReRPQ1dEBE[2][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)
   
   s1p1k = pow(fs1dEBE[2][pe][1]->GetBinContent(b)*fs1dEBE[2][pe][1]->GetBinEntries(b),1.); 
   s1p2k = pow(fs1dEBE[2][pe][2]->GetBinContent(b)*fs1dEBE[2][pe][2]->GetBinEntries(b),1.); 
  }else if(type == "RP")
   {
    // q_{m*n,k}: (Remark: m=1 is 0, k=0 iz zero (to be improved!)) 
    q1n2kRe = fReRPQ1dEBE[0][pe][0][2]->GetBinContent(fReRPQ1dEBE[0][pe][0][2]->GetBin(b))
            * fReRPQ1dEBE[0][pe][0][2]->GetBinEntries(fReRPQ1dEBE[0][pe][0][2]->GetBin(b));
    q1n2kIm = fImRPQ1dEBE[0][pe][0][2]->GetBinContent(fImRPQ1dEBE[0][pe][0][2]->GetBin(b))
            * fImRPQ1dEBE[0][pe][0][2]->GetBinEntries(fImRPQ1dEBE[0][pe][0][2]->GetBin(b));
    q2n1kRe = fReRPQ1dEBE[0][pe][1][1]->GetBinContent(fReRPQ1dEBE[0][pe][1][1]->GetBin(b))
            * fReRPQ1dEBE[0][pe][1][1]->GetBinEntries(fReRPQ1dEBE[0][pe][1][1]->GetBin(b));
    q2n1kIm = fImRPQ1dEBE[0][pe][1][1]->GetBinContent(fImRPQ1dEBE[0][pe][1][1]->GetBin(b))
            * fImRPQ1dEBE[0][pe][1][1]->GetBinEntries(fImRPQ1dEBE[0][pe][1][1]->GetBin(b));
    // s_{1,1}, s_{1,2} and s_{1,3} // to be improved (add explanation)  
    s1p1k = pow(fs1dEBE[0][pe][1]->GetBinContent(b)*fs1dEBE[0][pe][1]->GetBinEntries(b),1.); 
    s1p2k = pow(fs1dEBE[0][pe][2]->GetBinContent(b)*fs1dEBE[0][pe][2]->GetBinEntries(b),1.); 
    //s1p3k = pow(fs1dEBE[0][pe][3]->GetBinContent(b)*fs1dEBE[0][pe][3]->GetBinEntries(b),1.); 
  }    
  
  if(type == "POI")
  {
   // p_{m*n,k}:   
   p1n0kRe = fReRPQ1dEBE[1][pe][0][0]->GetBinContent(fReRPQ1dEBE[1][pe][0][0]->GetBin(b))
           * fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b));
   p1n0kIm = fImRPQ1dEBE[1][pe][0][0]->GetBinContent(fImRPQ1dEBE[1][pe][0][0]->GetBin(b))  
           * fImRPQ1dEBE[1][pe][0][0]->GetBinEntries(fImRPQ1dEBE[1][pe][0][0]->GetBin(b));
   mp = fReRPQ1dEBE[1][pe][0][0]->GetBinEntries(fReRPQ1dEBE[1][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here) 
   // M01 from Eq. (118) in QC2c (to be improved (notation)):
   dM01 = mp*dSM1p1k-s1p1k;
   dM011 = mp*(dSM2p1k-dSM1p2k)
         - 2.*(s1p1k*dSM1p1k-s1p2k);  
   // typeFlag = RP (0) or POI (1):   
   t = 1;           
  } else if(type == "RP")
    { 
     // to be improved (cross-checked):
     p1n0kRe = fReRPQ1dEBE[0][pe][0][0]->GetBinContent(fReRPQ1dEBE[0][pe][0][0]->GetBin(b))
             * fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b));
     p1n0kIm = fImRPQ1dEBE[0][pe][0][0]->GetBinContent(fImRPQ1dEBE[0][pe][0][0]->GetBin(b))  
             * fImRPQ1dEBE[0][pe][0][0]->GetBinEntries(fImRPQ1dEBE[0][pe][0][0]->GetBin(b));
     mp = fReRPQ1dEBE[0][pe][0][0]->GetBinEntries(fReRPQ1dEBE[0][pe][0][0]->GetBin(b)); // to be improved (cross-checked by accessing other profiles here)    
     // M01 from Eq. (118) in QC2c (to be improved (notation)):
     dM01 = mp*dSM1p1k-s1p1k;
     dM011 = mp*(dSM2p1k-dSM1p2k)
           - 2.*(s1p1k*dSM1p1k-s1p2k); 
     // typeFlag = RP (0) or POI (1): 
     t = 0;
    }
  
  // <<sin n(psi1)>>:
  Double_t sinP1nPsi = 0.;
  if(mp)
  {
   sinP1nPsi = p1n0kIm/mp;
   
   // fill profile for <<sin n(psi1)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][0][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sinP1nPsi,mp);
   // histogram to store <sin n(psi1)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][0][0]->SetBinContent(b,sinP1nPsi);
  } // end of if(mp)   
  
  // <<w2 sin n(psi1+phi2)>>:
  Double_t sinP1nPsiP1nPhiW2 = 0.;
  if(dM01)
  {
   sinP1nPsiP1nPhiW2 = (p1n0kRe*dImQ1n1k+p1n0kIm*dReQ1n1k-q2n1kIm)/(dM01);
   // fill profile for <<w2 sin n(psi1+phi2)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][0][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sinP1nPsiP1nPhiW2,dM01);
   // histogram to store <w2 sin n(psi1+phi2)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][0][1]->SetBinContent(b,sinP1nPsiP1nPhiW2);
  } // end of if(mp*dMult-mq)   
  
  // <<w2 w3 sin n(psi1+phi2-phi3)>>:
  Double_t sinP1nPsi1P1nPhi2MPhi3W2W3 = 0.;
  if(dM011)
  {
   sinP1nPsi1P1nPhi2MPhi3W2W3 = (p1n0kIm*(pow(dImQ1n1k,2.)+pow(dReQ1n1k,2.))
                              - p1n0kIm*dSM1p2k
                              + q2n1kRe*dImQ1n1k-q2n1kIm*dReQ1n1k
                              - s1p1k*dImQ1n1k
                              + 2.*q1n2kIm)
                              / dM011;  
   // fill profile for <<w2 w3 sin n(psi1+phi2-phi3)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][0][2]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sinP1nPsi1P1nPhi2MPhi3W2W3,dM011);
   // histogram to store <w2 w3 sin n(psi1+phi2-phi3)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][0][2]->SetBinContent(b,sinP1nPsi1P1nPhi2MPhi3W2W3);
  } // end of if(dM011)   
  
  // <<w2 w3 sin n(psi1-phi2-phi3)>>:
  Double_t sinP1nPsi1M1nPhi2MPhi3W2W3 = 0.;
  if(dM011)
  {
   sinP1nPsi1M1nPhi2MPhi3W2W3 = (p1n0kIm*(pow(dReQ1n1k,2.)-pow(dImQ1n1k,2.))-2.*p1n0kRe*dReQ1n1k*dImQ1n1k
                              + 1.*(p1n0kRe*dImQ2n2k-p1n0kIm*dReQ2n2k)  
                              + 2.*s1p1k*dImQ1n1k
                              - 2.*q1n2kIm)
                              / dM011;
   // fill profile for <<w2 w3 sin n(psi1-phi2-phi3)>>:
   fDiffFlowCorrectionTermsForNUAPro[t][pe][0][3]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],sinP1nPsi1M1nPhi2MPhi3W2W3,dM011);
   // histogram to store <w2 w3 sin n(psi1-phi2-phi3)> e-b-e (needed in some other methods):
   fDiffFlowCorrectionTermsForNUAEBE[t][pe][0][3]->SetBinContent(b,sinP1nPsi1M1nPhi2MPhi3W2W3);
  } // end of if(dM011)   
  
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)

} // end of AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights(TString type, TString ptOrEta)

//================================================================================================================================
   
void AliFlowAnalysisWithQCumulants::EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoopsUsingParticleWeights(AliFlowEventSimple * const anEvent, TString type, TString ptOrEta)
{
 // Evaluate with nested loops correction terms for non-uniform acceptance 
 // with using particle weights (both sin and cos terms) relevant for differential flow.
 
 // Remark 1: "w1" in expressions bellow is a particle weight used only for particles which were 
 //           flagged both as POI and RP.
 // Remark 2: Reduced correction terms for non-uniform acceptance are evaluated in pt bin number fCrossCheckInPtBinNo 
 //           and eta bin number fCrossCheckInEtaBinNo both for RPs and POIs.
 // Remark 3: Results are stored in 1 bin profiles fDiffFlowDirectCorrections[t][pe][sc][cti], where first three indices runs as: 
 //           [0=RP,1=POI][0=Pt,1=Eta][0=sin terms,1=cos terms], whilst the cti (correction term index) runs as follows: 
 //  cti: 
 //    0: <<sc n(psi1)>>
 //    1: <<w2 sc n(psi1+phi2)>> 
 //    2: <<w2 w3 sc n(psi1+phi2-phi3)>>
 //    3: <<w2 w3 sc n(psi1-phi2-phi3)>>
 //    4:
 //    5:
 //    6:
     
 Int_t typeFlag = 0;
 Int_t ptEtaFlag = 0;
 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   }      
 if(ptOrEta == "Pt")
 {
  ptEtaFlag = 0;
 } else if(ptOrEta == "Eta")
   {
    ptEtaFlag = 1;
   } 
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pe = ptEtaFlag;
      
 Double_t lowerPtEtaEdge[2] = {fPtMin+(fCrossCheckInPtBinNo-1)*fPtBinWidth,fEtaMin+(fCrossCheckInEtaBinNo-1)*fEtaBinWidth};
 Double_t upperPtEtaEdge[2] = {fPtMin+fCrossCheckInPtBinNo*fPtBinWidth,fEtaMin+fCrossCheckInEtaBinNo*fEtaBinWidth};
 Double_t binWidthPtEta[2] = {fPtBinWidth,fEtaBinWidth};
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 
 Double_t psi1=0., phi2=0., phi3=0.;// phi4=0.;// phi5=0., phi6=0., phi7=0., phi8=0.;
 Double_t wPhi2=1., wPhi3=1.;
 
 Int_t n = fHarmonic; 
 
 // 1'-particle correction terms:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  psi1=aftsTrack->Phi(); 
  // sin terms: 
  fDiffFlowDirectCorrectionTermsForNUA[t][pe][0][0]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,sin(n*psi1),1.); // <<sin(n*(psi1))>>  
  // cos terms: 
  fDiffFlowDirectCorrectionTermsForNUA[t][pe][1][0]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*psi1),1.); // <<cos(n*(psi1))>>  
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
   
 // 2'-particle correction terms:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  psi1=aftsTrack->Phi(); 
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1) continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP):
   if(!(aftsTrack->InRPSelection())) continue;
   phi2=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));   
   // sin terms: 
   fDiffFlowDirectCorrectionTermsForNUA[t][pe][0][1]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,sin(n*(psi1+phi2)),wPhi2); // <<w2 sin(n*(psi1+phi2))>>  
   // cos terms: 
   fDiffFlowDirectCorrectionTermsForNUA[t][pe][1][1]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*(psi1+phi2)),wPhi2); // <<w2 cos(n*(psi1+phi2))>>  
  }//end of for(Int_t i2=0;i2<nPrim;i2++)
 }//end of for(Int_t i1=0;i1<nPrim;i1++)   
 
 // 3'-particle correction terms:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): // to be improved (this can be implemented much better)
  if(typeFlag==1) // this is diff flow of POIs 
  {
   if(ptOrEta == "Pt")
   { 
    if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;
   } else if (ptOrEta == "Eta")
     {
      if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InPOISelection())))continue;    
     }
  } else // this is diff flow of RPs 
    {
     if(ptOrEta == "Pt")
     { 
      if(!((aftsTrack->Pt()>=lowerPtEtaEdge[pe] && aftsTrack->Pt()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;
     } else if (ptOrEta == "Eta")
       {
        if(!((aftsTrack->Eta()>=lowerPtEtaEdge[pe] && aftsTrack->Eta()<upperPtEtaEdge[pe]) && (aftsTrack->InRPSelection())))continue;    
       }
    }
  psi1=aftsTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1) continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP):
   if(!(aftsTrack->InRPSelection())) continue;
   phi2=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));   
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2) continue;
    aftsTrack=anEvent->GetTrack(i3);
    // RP condition (!(first) particle in the correlator must be RP):
    if(!(aftsTrack->InRPSelection())) continue;
    phi3=aftsTrack->Phi();
    if(fUsePhiWeights && fPhiWeights) wPhi3 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*fnBinsPhi/TMath::TwoPi())));   
    // sin terms: 
    fDiffFlowDirectCorrectionTermsForNUA[t][pe][0][2]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,sin(n*(psi1+phi2-phi3)),wPhi2*wPhi3); // <<wPhi2*wPhi3 sin(n*(psi1+phi2-phi3))>>  
    fDiffFlowDirectCorrectionTermsForNUA[t][pe][0][3]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,sin(n*(psi1-phi2-phi3)),wPhi2*wPhi3); // <<wPhi2*wPhi3 sin(n*(psi1-phi2-phi3))>>  
    // cos terms: 
    fDiffFlowDirectCorrectionTermsForNUA[t][pe][1][2]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*(psi1+phi2-phi3)),wPhi2*wPhi3); // <<wPhi2*wPhi3 cos(n*(psi1+phi2-phi3))>>  
    fDiffFlowDirectCorrectionTermsForNUA[t][pe][1][3]->Fill(lowerPtEtaEdge[pe]+binWidthPtEta[pe]/2.,cos(n*(psi1-phi2-phi3)),wPhi2*wPhi3); // <<wPhi2*wPhi3 cos(n*(psi1-phi2-phi3))>>  
   }//end of for(Int_t i3=0;i3<nPrim;i3++)  
  }//end of for(Int_t i2=0;i2<nPrim;i2++)  
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
               
} // end of void AliFlowAnalysisWithQCumulants::EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* anEvent, TString type, TString ptOrEta)

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CheckPointersUsedInFinish()
{
 // Check all pointers used in method Finish().
 
 if(!fAvMultiplicity)
 {
  cout<<endl;
  cout<<" WARNING (QC): fAvMultiplicity is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fIntFlowCorrelationsPro)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowCorrelationsPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!fIntFlowSquaredCorrelationsPro)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowSquaredCorrelationsPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!fIntFlowCorrelationsHist)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowCorrelationsHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if((fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights) && !fIntFlowExtraCorrelationsPro) 
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowExtraCorrelationsPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 for(Int_t power=0;power<2;power++)
 { 
  if(!fIntFlowSumOfEventWeights[power]) 
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowSumOfEventWeights[%d] is NULL in CheckPointersUsedInFinish() !!!!",power)<<endl;
   cout<<endl;
   exit(0);
  }
 } // end of for(Int_t power=0;power<2;power++)
 if(!fIntFlowProductOfCorrelationsPro)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowProductOfCorrelationsPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!fIntFlowSumOfProductOfEventWeights)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowSumOfProductOfEventWeights is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!fIntFlowCovariances)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowCovariances is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }  
 if(!fIntFlowQcumulants)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowQcumulants is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }  
 if(!fIntFlow)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlow is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!fCommonHists)
 {
  cout<<endl;
  cout<<" WARNING (QC): fCommonHists is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!(fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th))
 {
  cout<<endl;
  cout<<" WARNING (QC): fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th"<<endl; 
  cout<<"               && fCommonHistsResults8th is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 
 
 // NUA stuff:
 for(Int_t sc=0;sc<2;sc++) // sin/cos
 { 
  if(!fIntFlowCorrectionTermsForNUAPro[sc]) 
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowCorrectionTermsForNUAPro[%d] is NULL in CheckPointersUsedInFinish() !!!!",sc)<<endl;
   cout<<endl;
   exit(0);
  }
  if(!fIntFlowCorrectionTermsForNUAHist[sc]) 
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowCorrectionTermsForNUAHist[%d] is NULL in CheckPointersUsedInFinish() !!!!",sc)<<endl;
   cout<<endl;
   exit(0);
  }
  for(Int_t lq=0;lq<2;lq++) // linear/quadratic
  {
   if(!fIntFlowSumOfEventWeightsNUA[sc][lq]) 
   {
    cout<<endl;
    cout<<Form(" WARNING (QC): fIntFlowSumOfEventWeightsNUA[%d][%d] is NULL in CheckPointersUsedInFinish() !!!!",sc,lq)<<endl;
    cout<<endl;
    exit(0);
   }
  } // end of for(Int_t lq=0;lq<2;lq++) // linear/quadratic
 } // end of for(Int_t power=0;power<2;power++) 
 if(!fIntFlowProductOfCorrectionTermsForNUAPro)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowProductOfCorrectionTermsForNUAPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!fIntFlowSumOfProductOfEventWeightsNUA)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowSumOfProductOfEventWeightsNUA is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!fIntFlowCovariancesNUA)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowCovariancesNUA is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!fIntFlowQcumulantsErrorSquaredRatio)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowQcumulantsErrorSquaredRatio is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!fIntFlowDetectorBias)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowDetectorBias is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 
 // Versus multiplicity:
 if(!fCalculateCumulantsVsM){return;}
 for(Int_t co=0;co<=3;co++) // cumulant order
 {
  if(!fIntFlowQcumulantsVsM[co])
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowQcumulantsVsM[%d] is NULL in CheckPointersUsedInFinish() !!!!",co)<<endl;
   cout<<endl;
   exit(0); 
  }
  if(!fIntFlowVsM[co])
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowVsM[%d] is NULL in CheckPointersUsedInFinish() !!!!",co)<<endl;
   cout<<endl;
   exit(0); 
  }
  if(!fIntFlowDetectorBiasVsM[co])
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowDetectorBiasVsM[%d] is NULL in CheckPointersUsedInFinish() !!!!",co)<<endl;
   cout<<endl;
   exit(0); 
  }
 } // end of for(Int_t c0=0;c0<=3;c0++) // cumulant order
 for(Int_t ci=0;ci<=3;ci++) // correlation index
 {
  if(!fIntFlowCorrelationsVsMPro[ci])
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowCorrelationsVsMPro[%d] is NULL in CheckPointersUsedInFinish() !!!!",ci)<<endl;
   cout<<endl;
   exit(0); 
  }
  if(!fIntFlowSquaredCorrelationsVsMPro[ci])
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowSquaredCorrelationsVsMPro[%d] is NULL in CheckPointersUsedInFinish() !!!!",ci)<<endl;
   cout<<endl;
   exit(0); 
  }  
  if(!fIntFlowCorrelationsVsMHist[ci])
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowCorrelationsVsMHist[%d] is NULL in CheckPointersUsedInFinish() !!!!",ci)<<endl;
   cout<<endl;
   exit(0); 
  }
  for(Int_t power=0;power<2;power++) 
  {
   if(!fIntFlowSumOfEventWeightsVsM[ci][power])
   {
    cout<<endl;
    cout<<Form(" WARNING (QC): fIntFlowSumOfEventWeightsVsM[%d][%d] is NULL in CheckPointersUsedInFinish() !!!!",ci,power)<<endl;
    cout<<endl;
    exit(0);   
   }
  } // end of for(Int_t power=0;power<2;power++) 
 } // end of for(Int_t ci=0;ci<=3;ci++) // correlation index
 for(Int_t i=0;i<6;i++)
 {
  if(!fIntFlowProductOfCorrelationsVsMPro[i])
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowProductOfCorrelationsVsMPro[%d] is NULL in CheckPointersUsedInFinish() !!!!",i)<<endl;
   cout<<endl;
   exit(0); 
  }
  if(!fIntFlowSumOfProductOfEventWeightsVsM[i])
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowSumOfProductOfEventWeightsVsM[%d] is NULL in CheckPointersUsedInFinish() !!!!",i)<<endl;
   cout<<endl;
   exit(0); 
  }
  if(!fIntFlowCovariancesVsM[i])
  {
   cout<<endl;
   cout<<Form(" WARNING (QC): fIntFlowCovariancesVsM[%d] is NULL in CheckPointersUsedInFinish() !!!!",i)<<endl;
   cout<<endl;
   exit(0); 
  }
 } // end of for(Int_t i=0;i<6;i++) 
 if(!fIntFlowRebinnedInM)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowRebinnedInM is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!fIntFlowQcumulantsRebinnedInM)
 {
  cout<<endl;
  cout<<" WARNING (QC): fIntFlowQcumulantsRebinnedInM is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }  
 
} // end of void AliFlowAnalysisWithQCumulants::CheckPointersUsedInFinish()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::CheckPointersUsedInMake()
{
 // Check all pointers used in method Make(). // to be improved - check other pointers as well
 
 if(!fAvMultiplicity)
 {
  printf("\n WARNING (QC): fAvMultiplicity is NULL in CheckPointersUsedInMake() !!!!\n\n");
  exit(0);
 }
 if((fUsePhiWeights||fUsePtWeights||fUseEtaWeights||fUseTrackWeights) && !fIntFlowExtraCorrelationsPro) 
 {
  printf("\n WARNING (QC): fIntFlowExtraCorrelationsPro is NULL in CheckPointersUsedInMake() !!!!\n\n");
  exit(0); 
 } 
 // 2D:
 if(fCalculate2DDiffFlow)
 {
  for(Int_t t=0;t<2;t++) // type = RP or POI
  { 
   for(Int_t rci=0;rci<4;rci++) // reduced correlation index
   {
    if(!f2DDiffFlowCorrelationsPro[t][rci])
    {
     printf("\n WARNING (QC): f2DDiffFlowCorrelationsPro[%i][%i] is NULL in CheckPointersUsedInMake() !!!!\n\n",t,rci);
     exit(0);     
    } // end of if(!f2DDiffFlowCorrelationsPro[t][rci])  
   } // end of for(Int_t rci=0;rci<4;rci++) // reduced correlation index
  } // end of for(Int_t t=0;t<2;t++)
 } // end of if(fCalculate2DDiffFlow)  

} // end of void AliFlowAnalysisWithQCumulants::CheckPointersUsedInMake()
 

