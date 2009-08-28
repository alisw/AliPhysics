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
 * author:  Ante Bilandzic        * 
 *           (anteb@nikhef.nl)    *
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
#include "TProfile3D.h"
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
 fHarmonic(2),
 fAnalysisLabel(NULL),
 // 2.) weights:
 fWeightsList(NULL),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fUseParticleWeights(NULL),
 fPhiWeights(NULL),
 fPtWeights(NULL),
 fEtaWeights(NULL),
 // 3.) integrated flow:
 fIntFlowList(NULL), 
 fIntFlowProfiles(NULL),
 fIntFlowResults(NULL),
 fReQ(NULL),
 fImQ(NULL),
 fSMpk(NULL),
 fAvMultiplicity(NULL),
 // 4.) differential flow:
 fDiffFlowList(NULL),
 fDiffFlowProfiles(NULL),
 fDiffFlowResults(NULL),
 // 5.) distributions:
 fDistributionsList(NULL),
 // x.) debugging and cross-checking:
 fNestedLoopsList(NULL),
 fEvaluateNestedLoopsForIntFlow(kFALSE),
 fEvaluateNestedLoopsForDiffFlow(kFALSE),  
 fEvaluateNestedLoops(NULL),
 fDirectCorrelations(NULL),
 fDirectCorrectionsCos(NULL),
 fDirectCorrectionsSin(NULL),
 fDirectCorrelationsDiffFlow(NULL),
 fDirectCorrectionsDiffFlowCos(NULL),
 fDirectCorrectionsDiffFlowSin(NULL),
 fDirectCorrelationsW(NULL),
 fDirectCorrectionsCosW(NULL),
 fDirectCorrectionsSinW(NULL),
 fDirectCorrelationsDiffFlowW(NULL),
 fDirectCorrectionsDiffFlowCosW(NULL),
 fDirectCorrectionsDiffFlowSinW(NULL)
 {
  // constructor  
  
  // base list to hold all output objects:
  fHistList = new TList();
  fHistList->SetName("cobjQC");
  fHistList->SetOwner(kTRUE);
  
  // list to hold histograms with phi, pt and eta weights:      
  fWeightsList = new TList();
    
  // analysis label;
  fAnalysisLabel = new TString();
      
  // initialize all arrays:  
  this->InitializeArraysForIntFlow();
  this->InitializeArraysForDiffFlow();
  this->InitializeArraysForDistributions();
  
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
 // initialize all constants and book everything
 
 // access constants:
 this->AccessConstants();
 
 // booking:
 this->BookAndFillWeightsHistograms();
 this->BookAndNestAllLists();
 this->BookCommonHistograms();
 this->BookEverythingForIntegratedFlow(); 
 this->BookEverythingForDifferentialFlow(); 
 this->BookEverythingForDistributions();
 this->BookEverythingForNestedLoops();
 
 // set harmonic in common control histograms (to be improved (should I do this somewhere else?)):
 (fCommonHists->GetHarmonic())->Fill(0.5,fHarmonic);
 (fCommonHists2nd->GetHarmonic())->Fill(0.5,fHarmonic);
 (fCommonHists4th->GetHarmonic())->Fill(0.5,fHarmonic);
 (fCommonHists6th->GetHarmonic())->Fill(0.5,fHarmonic);
 (fCommonHists8th->GetHarmonic())->Fill(0.5,fHarmonic);
 
} // end of void AliFlowAnalysisWithQCumulants::Init()


//================================================================================================================


void AliFlowAnalysisWithQCumulants::Make(AliFlowEventSimple* anEvent)
{
 // running over data only in this method
 
 // a) fill the common control histograms and call method to fill fAvMultiplicity;
 // b) loop over data to calculate e-b-e quantities;
 // c) call the methods;
 // d) debugging and cross-checking (evaluate nested loops);
 // e) reset e-b-e quantities.
 
 Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
 Double_t dPt  = 0.; // transverse momentum
 Double_t dEta = 0.; // pseudorapidity

 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight
                                                                                                                                
 // ********************************************
 // **** FILL THE COMMON CONTROL HISTOGRAMS ****
 // ********************************************
                                         
 Int_t nRP = anEvent->GetEventNSelTracksRP(); // number of RPs (i.e. number of particles used to determine the reaction plane)

 fCommonHists->FillControlHistograms(anEvent); 
 
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
                             
 this->FillAverageMultiplicities(nRP);                                                                  
                                                                                                                                            
 // *******************************************************
 // **** LOOP OVER DATA AND CALCULATE E-B-E QUANTITIES ****
 // *******************************************************
 
 Int_t nPrim = anEvent->NumberOfTracks();  // nPrim = total number of primary tracks, i.e. nPrim = nRP + nPOI + rest, where:
                                           // nRP   = # of particles used to determine the reaction plane;
                                           // nPOI  = # of particles of interest for a detailed flow analysis;
                                           // rest  = # of particles which are not niether RPs nor POIs.  
 
 AliFlowTrackSimple *aftsTrack = NULL;
 
 for(Int_t i=0;i<nPrim;i++) 
 { 
  aftsTrack=anEvent->GetTrack(i);
  if(aftsTrack)
  {
   if(!(aftsTrack->InRPSelection() || aftsTrack->InPOISelection())) continue; // consider only tracks which are RPs or POIs
   
   Int_t n = fHarmonic; // shortcut for the harmonic
 
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
    
    // integrated flow: 
    // calculate Re[Q_{m*n,k}] and Im[Q_{m*n,k}], m = 1,2,3,4, for this event:
    for(Int_t m=0;m<4;m++)
    {
     for(Int_t k=0;k<9;k++)
     {
      (*fReQ)(m,k)+=pow(wPhi*wPt*wEta,k)*TMath::Cos((m+1)*n*dPhi); 
      (*fImQ)(m,k)+=pow(wPhi*wPt*wEta,k)*TMath::Sin((m+1)*n*dPhi); 
     } 
    }
    // calculate S^{M}_{p,k} for this event 
    // Remark: final calculation of S^{M}_{p,k} follows after the loop over data bellow:
    for(Int_t p=0;p<8;p++)
    {
     for(Int_t k=0;k<9;k++)
     {     
      (*fSMpk)(p,k)+=pow(wPhi*wPt*wEta,k);
     }
    } 

    // differential flow:
    // (r_{m*m,k}(pt,eta)): 
    for(Int_t m=0;m<4;m++)
    {
     for(Int_t k=0;k<9;k++)
     {
      fReEBE[0][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k)*TMath::Cos((m+1.)*n*dPhi),1.);
      fImEBE[0][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k)*TMath::Sin((m+1.)*n*dPhi),1.);
     }
    } 
    // s_{k}(pt,eta) for RPs // to be improved (clarified)
    // Remark: final calculation of s_{p,k}(pt,eta) follows after the loop over data bellow:
    for(Int_t k=0;k<9;k++)
    {
     fs[0][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k),1.);
    } 
     
    if(aftsTrack->InPOISelection())
    {
     // (q_{m*m,k}(pt,eta)): 
     for(Int_t m=0;m<4;m++)
     {
      for(Int_t k=0;k<9;k++)
      {
       fReEBE[2][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k)*TMath::Cos((m+1.)*n*dPhi),1.);
       fImEBE[2][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k)*TMath::Sin((m+1.)*n*dPhi),1.);
      }
     } 
     // s_{k}(pt,eta) for RP&&POIs // to be improved (clarified)
     // Remark: final calculation of s_{p,k}(pt,eta) follows after the loop over data bellow:
     for(Int_t k=0;k<9;k++)
     {
      fs[2][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k),1.);
     }
      
    } // end of if(aftsTrack->InPOISelection())
   } // end of if(pTrack->InRPSelection())

   if(aftsTrack->InPOISelection())
   {
    dPhi = aftsTrack->Phi();
    dPt  = aftsTrack->Pt();
    dEta = aftsTrack->Eta();
       
    // p_n(m*n,0):   
    for(Int_t m=0;m<4;m++)
    {
     fReEBE[1][m][0]->Fill(dPt,dEta,TMath::Cos((m+1.)*n*dPhi),1.);
     fImEBE[1][m][0]->Fill(dPt,dEta,TMath::Sin((m+1.)*n*dPhi),1.);
    } 
    
   } // end of if(pTrack->InPOISelection() )   
 
  } else // to if(aftsTrack)
    {
     cout<<endl;
     cout<<" WARNING: no particle! (i.e. aftsTrack is a NULL pointer in AFAWQC::Make().)"<<endl;
     cout<<endl;       
    }
 } // end of for(Int_t i=0;i<nPrim;i++) 

 // calculate the final expressions for S^{M}_{p,k}:
 for(Int_t p=0;p<8;p++)
 {
  for(Int_t k=0;k<9;k++)
  {
   (*fSMpk)(p,k)=pow((*fSMpk)(p,k),p+1);
  }  
 } 
 
 // *****************************
 // **** CALL THE METHODS *******
 // *****************************

 if(!fEvaluateNestedLoopsForIntFlow)
 {
  // without weights:
  if(nRP>1) this->CalculateCorrelationsForIntegratedFlow();
  if(nRP>0) this->CalculateCorrectionsForNonUniformAcceptanceForIntFlowCosTerms();
  if(nRP>0) this->CalculateCorrectionsForNonUniformAcceptanceForIntFlowSinTerms();
  if(nRP>3) this->CalculateQProductsForIntFlow();
  if(nRP>1) this->CalculateSumAndProductOfEventWeights();
  // with weights:
  if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
  {
   if(nRP>1) this->CalculateWeightedCorrelationsForIntegratedFlow();
   if(nRP>3) this->CalculateWeightedQProductsForIntFlow();
  } 
 }

 if(!fEvaluateNestedLoopsForDiffFlow)
 {
  // without weights:
  if(nRP>1) this->CalculateCorrelationsForDifferentialFlow("RP");
  if(nRP>1) this->CalculateCorrelationsForDifferentialFlow("POI");
  
  // with weights:
  if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
  {
   if(nRP>1) this->CalculateWeightedCorrelationsForDifferentialFlow("RP");
   if(nRP>1) this->CalculateWeightedCorrelationsForDifferentialFlow("POI");
  } 
 }
 
 // **************************************************************
 // **** DEBUGGING AND CROSS-CHECKING (EVALUATE NESTED LOOPS) ****
 // **************************************************************

 if(fEvaluateNestedLoopsForIntFlow)
 {
  if(nPrim>0 && nPrim<15) // only for these multiplicities it is feasible to evaluate 8 nested loops in short time 
  {
   // without weights:
   if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
   {
    this->CalculateCorrelationsForIntegratedFlow();
    this->CalculateCorrectionsForNonUniformAcceptanceForIntFlowCosTerms();
    this->CalculateCorrectionsForNonUniformAcceptanceForIntFlowSinTerms();
   }
   // with weights:
   if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
   {
    this->CalculateWeightedCorrelationsForIntegratedFlow();
   }
    
   this->EvaluateNestedLoopsForIntegratedFlow(anEvent);  
  }
 } 
 
 if(fEvaluateNestedLoopsForDiffFlow)
 {
  if(nPrim>0 && nPrim<15) // only for these multiplicities it is feasible to evaluate 8 nested loops in short time 
  {
   // without weights:
   if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
   {
    this->CalculateCorrelationsForDifferentialFlow("RP");
    this->CalculateCorrelationsForDifferentialFlow("POI");
   }
   // with weights:
   if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
   {
    this->CalculateWeightedCorrelationsForDifferentialFlow("RP");
    this->CalculateWeightedCorrelationsForDifferentialFlow("POI");
   }
    
   this->EvaluateNestedLoopsForDifferentialFlow(anEvent);  
  }
 } 
 
 // ********************************
 // **** RESET E-B-E QUANTITIES ****
 // ********************************
 
 // integrated flow:
 fReQ->Zero();
 fImQ->Zero();
 fSMpk->Zero();
 for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++)
 {
  fQCorrelationsEBE[pW]->Reset();
  for(Int_t sc=0;sc<2;sc++)
  {
    fQCorrectionsEBE[pW][sc]->Reset();
  } 
 }  
 
 // differential flow:
 for(Int_t t=0;t<3;t++) // type (RP, POI, POI&&RP)
 {
  for(Int_t m=0;m<4;m++) // multiple of harmonic
  {
   for(Int_t k=0;k<9;k++) // power of weight
   {
    if(fReEBE[t][m][k]) fReEBE[t][m][k]->Reset();
    if(fImEBE[t][m][k]) fImEBE[t][m][k]->Reset();
   }   
  }
 }
 
 for(Int_t t=0;t<3;t++) // type (0 = RP, 1 = POI, 2 = RP&&POI )
 { 
  for(Int_t k=0;k<9;k++)
  {
   if(fs[t][k]) fs[t][k]->Reset();
  }
 }
 
} // end of AliFlowAnalysisWithQCumulants::Make(AliFlowEventSimple* anEvent)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::Finish()
{
 // calculate the final results

 // a) acces the constants;
 // b) access the flags;
 // c) calculate the final results for integrated flow (without and with weights);
 // d) store in AliFlowCommonHistResults and print the final results for integrated flow;
 // e) calculate the final results for differential flow (without and with weights);
 // f) print the final results for integrated flow obtained from differential flow (to be improved (terminology));
 // g) COMPARE RESULTS FROM NESTED LOOPS vs RESULTS FROM Q-VECTORS FOR INTEGRATED FLOW
 
 // ******************************
 // **** ACCESS THE CONSTANTS ****
 // ******************************
 
 this->AccessConstants();          
 
 // **************************
 // **** ACCESS THE FLAGS ****
 // **************************    

 fUsePhiWeights = (Int_t)fUseParticleWeights->GetBinContent(1); 
 fUsePtWeights = (Int_t)fUseParticleWeights->GetBinContent(2); 
 fUseEtaWeights = (Int_t)fUseParticleWeights->GetBinContent(3); 
 fEvaluateNestedLoopsForIntFlow = (Int_t)fEvaluateNestedLoops->GetBinContent(1);
 fEvaluateNestedLoopsForDiffFlow = (Int_t)fEvaluateNestedLoops->GetBinContent(2); 
    
 // *********************************************************
 // **** CALCULATE THE FINAL RESULTS FOR INTEGRATED FLOW ****
 // *********************************************************    
 
 // without weights:
 this->FinalizeCorrelationsForIntFlow(kFALSE,"exact");
 this->CalculateFinalCorrectionsForNonUniformAcceptanceForCumulantsForIntFlow(kFALSE,"exact");
 this->CalculateCovariancesForIntFlow(kFALSE,"exact");
 this->CalculateCumulantsForIntFlow(kFALSE,"exact");
 this->ApplyCorrectionForNonUniformAcceptanceToCumulantsForIntFlow(kFALSE,"exact");
 this->CalculateIntFlow(kFALSE,"exact",kFALSE); // pW = 0, eW = 0, not corrected for non-uniform acceptance
 this->CalculateIntFlow(kFALSE,"exact",kTRUE); // pW = 0, eW = 0, corrected for non-uniform acceptance
 
 // with weights:
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
 {
  this->FinalizeCorrelationsForIntFlow(kTRUE,"exact"); 
  // this->CalculateFinalCorrectionsForNonUniformAcceptanceForCumulantsForIntFlow(kTRUE,"exact");
  this->CalculateCovariancesForIntFlow(kTRUE,"exact");
  this->CalculateCumulantsForIntFlow(kTRUE,"exact");
  // this->ApplyCorrectionForNonUniformAcceptanceToCumulantsForIntFlow(kTRUE,"exact");
  this->CalculateIntFlow(kTRUE,"exact",kFALSE); // weighted and not corrected for non-uniform acceptance
  // this->CalculateIntFlow(kTRUE,"exact",kTRUE); // weighted and corrected for non-uniform acceptance
 }
 
 // ***************************************************************
 // **** STORE AND PRINT THE FINAL RESULTS FOR INTEGRATED FLOW ****
 // ***************************************************************
 
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
 {        
  this->FillCommonHistResultsIntFlow(kTRUE,"exact",kFALSE); // weighted and not corrected for non-uniform acceptance       
  // this->FillCommonHistResultsIntFlow(kTRUE,kTRUE); // weighted and corrected for non-uniform acceptance (to be improved (enabled))    
  // this->PrintQuantifyingCorrectionsForNonUniformAcceptance(kTRUE,"exact"); // (to be improved (enabled))
 } else 
   {
    this->FillCommonHistResultsIntFlow(kFALSE,"exact",kTRUE); // non-weighted and corrected for non-uniform acceptance 
    this->PrintQuantifyingCorrectionsForNonUniformAcceptance(kFALSE,"exact"); 
   }
 
 this->PrintFinalResultsForIntegratedFlow("NONAME"); // to be improved (name)
 
 // ***********************************************************
 // **** CALCULATE THE FINAL RESULTS FOR DIFFERENTIAL FLOW ****
 // ***********************************************************    
 
 // without weights:
 this->FinalizeCorrelationsForDiffFlow("RP",kFALSE,"exact");
 this->FinalizeCorrelationsForDiffFlow("POI",kFALSE,"exact");
 this->CalculateCumulantsForDiffFlow("RP",kFALSE,"exact");
 this->CalculateCumulantsForDiffFlow("POI",kFALSE,"exact");
 this->CalculateDiffFlow("RP",kFALSE,"exact");
 this->CalculateDiffFlow("POI",kFALSE,"exact");
 this->CalculateFinalResultsForRPandPOIIntegratedFlow("RP",kFALSE,"exact");
 this->CalculateFinalResultsForRPandPOIIntegratedFlow("POI",kFALSE,"exact");
 
 // with weights:
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
 {
  this->FinalizeCorrelationsForDiffFlow("RP",kTRUE,"exact");
  this->FinalizeCorrelationsForDiffFlow("POI",kTRUE,"exact");
  this->CalculateCumulantsForDiffFlow("RP",kTRUE,"exact");
  this->CalculateCumulantsForDiffFlow("POI",kTRUE,"exact");
  this->CalculateDiffFlow("RP",kTRUE,"exact");
  this->CalculateDiffFlow("POI",kTRUE,"exact");
  this->CalculateFinalResultsForRPandPOIIntegratedFlow("RP",kTRUE,"exact");
  this->CalculateFinalResultsForRPandPOIIntegratedFlow("POI",kTRUE,"exact");
 }
 
   
  //this->CalculateFinalCorrectionsForNonUniformAcceptanceForDifferentialFlow(kFALSE,"POI"); // to be improved (to calculate also when weights are used) 
  //this->CalculateFinalCorrectionsForNonUniformAcceptanceForDifferentialFlow(kFALSE,"RP"); // to be improved (to calculate also when weights are used)

 
 // *****************************************************************
 // **** STORE AND PRINT THE FINAL RESULTS FOR DIFFERENTIAL FLOW ****
 // *****************************************************************
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
 {
  this->FillCommonHistResultsDiffFlow("RP",kTRUE,"exact",kFALSE);
  this->FillCommonHistResultsDiffFlow("POI",kTRUE,"exact",kFALSE);
 } else
   {
    this->FillCommonHistResultsDiffFlow("RP",kFALSE,"exact",kFALSE);
    this->FillCommonHistResultsDiffFlow("POI",kFALSE,"exact",kFALSE);
   }
 
 this->PrintFinalResultsForIntegratedFlow("RP"); 
 this->PrintFinalResultsForIntegratedFlow("POI"); 
  
 // *****************************************************************************************
 // **** COMPARE RESULTS FROM NESTED LOOPS vs RESULTS FROM Q-VECTORS FOR INTEGRATED FLOW ****
 // *****************************************************************************************    
 
 if(fEvaluateNestedLoopsForIntFlow) 
 {
  this->CompareResultsFromNestedLoopsAndFromQVectorsForIntFlow(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);
 } 
 
 if(fEvaluateNestedLoopsForDiffFlow) 
 {
  this->CompareResultsFromNestedLoopsAndFromQVectorsForDiffFlow(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);
 } 
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithQCumulants::Finish()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCorrectionsForNonUniformAcceptanceForIntFlowCosTerms()
{
 // calculate corrections for non-uniform acceptance of the detector for no-name integrated flow (cos terms)
 
 // multiplicity:
 Double_t dMult = (*fSMpk)(0,0);
 
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
 // Remark 2: binning of fQCorrectionsCos is organized as follows:
 // --------------------------------------------------------------------------------------------------------------------
 // 1st bin: <<cos(n*(phi1))>> = cosP1n
 // 2nd bin: <<cos(n*(phi1+phi2))>> = cosP1nP1n
 // 3rd bin: <<cos(n*(phi1-phi2-phi3))>> = cosP1nM1nM1n
 // ...
 // --------------------------------------------------------------------------------------------------------------------
  
 // 1-particle:
 Double_t cosP1n = 0.; // <<cos(n*(phi1))>>
   
 if(dMult>0)
 {
  cosP1n = dReQ1n/dMult; 
  
  // average non-weighted 1-particle correction (cos terms) for non-uniform acceptance for single event:
  fQCorrectionsEBE[0][1]->SetBinContent(1,cosP1n);
  
  // final average non-weighted 1-particle correction (cos terms) for non-uniform acceptance for all events:
  fQCorrections[0][0][1]->Fill(0.5,cosP1n,dMult);  
 } 
 
 // 2-particle:
 Double_t cosP1nP1n = 0.; // <<cos(n*(phi1+phi2))>>
 
 if(dMult>1)
 {
  cosP1nP1n = (pow(dReQ1n,2)-pow(dImQ1n,2)-dReQ2n)/(dMult*(dMult-1)); 
  
  // average non-weighted 2-particle correction (cos terms) for non-uniform acceptance for single event:
  fQCorrectionsEBE[0][1]->SetBinContent(2,cosP1nP1n);
  
  // final average non-weighted 2-particle correction (cos terms) for non-uniform acceptance for all events:
  fQCorrections[0][0][1]->Fill(1.5,cosP1nP1n,dMult*(dMult-1));  
 } 
 
 // 3-particle:
 Double_t cosP1nM1nM1n = 0.; // <<cos(n*(phi1-phi2-phi3))>>
 
 if(dMult>2)
 {
  cosP1nM1nM1n = (dReQ1n*(pow(dReQ1n,2)+pow(dImQ1n,2))-dReQ1n*dReQ2n-dImQ1n*dImQ2n-2.*(dMult-1)*dReQ1n)
               / (dMult*(dMult-1)*(dMult-2)); 
  
  // average non-weighted 3-particle correction (cos terms) for non-uniform acceptance for single event:
  fQCorrectionsEBE[0][1]->SetBinContent(3,cosP1nM1nM1n);
  
  // final average non-weighted 3-particle correction (cos terms) for non-uniform acceptance for all events:
  fQCorrections[0][0][1]->Fill(2.5,cosP1nM1nM1n,dMult*(dMult-1)*(dMult-2));  
 } 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateCorrectionsForNonUniformAcceptanceForIntFlowCosTerms()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCorrectionsForNonUniformAcceptanceForIntFlowSinTerms()
{
 // calculate corrections for non-uniform acceptance of the detector for no-name integrated flow (sin terms)
 
 // multiplicity:
 Double_t dMult = (*fSMpk)(0,0);
 
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
 // Remark 2: binning of fQCorrectionsSin is organized as follows:
 // --------------------------------------------------------------------------------------------------------------------
 // 1st bin: <<sin(n*(phi1))>> = sinP1n
 // 2nd bin: <<sin(n*(phi1+phi2))>> = sinP1nP1n
 // 3rd bin: <<sin(n*(phi1-phi2-phi3))>> = sinP1nM1nM1n
 // ...
 // --------------------------------------------------------------------------------------------------------------------
 
 // 1-particle:
 Double_t sinP1n = 0.; // <sin(n*(phi1))>
 
 if(dMult>0)
 {
  sinP1n = dImQ1n/dMult; 
     
  // average non-weighted 1-particle correction (sin terms) for non-uniform acceptance for single event:
  fQCorrectionsEBE[0][0]->SetBinContent(1,sinP1n);
  
  // final average non-weighted 1-particle correction (sin terms) for non-uniform acceptance for all events:   
  fQCorrections[0][0][0]->Fill(0.5,sinP1n,dMult);  
 } 
 
 // 2-particle:
 Double_t sinP1nP1n = 0.; // <<sin(n*(phi1+phi2))>>
 
 if(dMult>1)
 {
  sinP1nP1n = (2.*dReQ1n*dImQ1n-dImQ2n)/(dMult*(dMult-1)); 
     
  // average non-weighted 2-particle correction (sin terms) for non-uniform acceptance for single event:
  fQCorrectionsEBE[0][0]->SetBinContent(2,sinP1nP1n);
  
  // final average non-weighted 1-particle correction (sin terms) for non-uniform acceptance for all events:      
  fQCorrections[0][0][0]->Fill(1.5,sinP1nP1n,dMult*(dMult-1));  
 } 
 
 // 3-particle:
 Double_t sinP1nM1nM1n = 0.; // <<sin(n*(phi1-phi2-phi3))>>
 
 if(dMult>2)
 {
  sinP1nM1nM1n = (-dImQ1n*(pow(dReQ1n,2)+pow(dImQ1n,2))+dReQ1n*dImQ2n-dImQ1n*dReQ2n+2.*(dMult-1)*dImQ1n)
               / (dMult*(dMult-1)*(dMult-2)); 
  
  // average non-weighted 3-particle correction (sin terms) for non-uniform acceptance for single event:
  fQCorrectionsEBE[0][0]->SetBinContent(3,sinP1nM1nM1n);
  
  // final average non-weighted 3-particle correction (sin terms) for non-uniform acceptance for all events:  
  fQCorrections[0][0][0]->Fill(2.5,sinP1nM1nM1n,dMult*(dMult-1)*(dMult-2));  
 } 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateCorrectionsForNonUniformAcceptanceForIntFlowSinTerms()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCorrectionsForNonUniformAcceptanceForDifferentialFlowCosTerms(TString type)
{
 // calculate corrections for non-uniform acceptance of the detector for differential flow (cos terms)
 
 // multiplicity:
 //Double_t dMult = (*fSMpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 //Double_t dReQ1n = (*fReQ)(0,0);
 //Double_t dReQ2n = (*fReQ)(1,0);
 //Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 //Double_t dImQ1n = (*fImQ)(0,0);
 //Double_t dImQ2n = (*fImQ)(1,0);
 //Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);

 // looping over all (pt,eta) bins and calculating correlations needed for differential flow: 
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   // real and imaginary parts of q_n (non-weighted Q-vector evaluated only for POIs in harmonic n for each (pt,eta) bin): 
   //Double_t dReqnPtEta = 0.;
   //Double_t dImqnPtEta = 0.;

   // number of POIs in each (pt,eta) bin:
   Double_t dmPtEta = 0.;

   // real and imaginary parts of q''_{n}, q''_{2n}, ... 
   // (non-weighted Q-vectors evaluated only for particles which are both RPs and POIs in harmonic n, 2n, ... for each (pt,eta) bin): 
   //Double_t dReqPrimePrime1nPtEta = 0.;
   //Double_t dImqPrimePrime1nPtEta = 0.;
   //Double_t dReqPrimePrime2nPtEta = 0.;
   //Double_t dImqPrimePrime2nPtEta = 0.;

   // number of particles which are both RPs and POIs in each (pt,eta) bin:
   //Double_t dmPrimePrimePtEta = 0.;
   
   if(type == "POI")
   {
    // q''_{n}, q''_{2n}:
    //...............................................................................................
    //dReqPrimePrime1nPtEta = fReqPrimePrime1nPtEta->GetBinContent(fReqPrimePrime1nPtEta->GetBin(p,e));
    //dImqPrimePrime1nPtEta = fImqPrimePrime1nPtEta->GetBinContent(fImqPrimePrime1nPtEta->GetBin(p,e));
    //dReqPrimePrime2nPtEta = fReqPrimePrime2nPtEta->GetBinContent(fReqPrimePrime2nPtEta->GetBin(p,e));
    //dImqPrimePrime2nPtEta = fImqPrimePrime2nPtEta->GetBinContent(fImqPrimePrime2nPtEta->GetBin(p,e));
    //...............................................................................................
   
    // m'':
    //dmPrimePrimePtEta = fmPrimePrimePtEta->GetBinContent(fmPrimePrimePtEta->GetBin(p,e));
   
    // q'_{n}: 
    //dReqnPtEta = fReqnPtEta->GetBinContent(fReqnPtEta->GetBin(p,e));
    //dImqnPtEta = fImqnPtEta->GetBinContent(fImqnPtEta->GetBin(p,e));
    //dmPtEta    = fmPtEta->GetBinContent(fmPtEta->GetBin(p,e));
   }
   else if(type == "RP")
   {
    // q_RP{n}, q_RP{2n}:
    //...............................................................................................
    //dReqPrimePrime1nPtEta = fReqRP1nPtEta->GetBinContent(fReqRP1nPtEta->GetBin(p,e));
    //dImqPrimePrime1nPtEta = fImqRP1nPtEta->GetBinContent(fImqRP1nPtEta->GetBin(p,e));
    //dReqPrimePrime2nPtEta = fReqRP2nPtEta->GetBinContent(fReqRP2nPtEta->GetBin(p,e));
    //dImqPrimePrime2nPtEta = fImqRP2nPtEta->GetBinContent(fImqRP2nPtEta->GetBin(p,e));
    //...............................................................................................
   
    // m'':
    //dmPrimePrimePtEta = fmRPPtEta->GetBinContent(fmRPPtEta->GetBin(p,e));
   
    //dReqnPtEta = fReqRP1nPtEta->GetBinContent(fReqRP1nPtEta->GetBin(p,e)); // not a bug ;-)
    //dImqnPtEta = fImqRP1nPtEta->GetBinContent(fImqRP1nPtEta->GetBin(p,e)); // not a bug ;-)
    //dmPtEta    = fmRPPtEta->GetBinContent(fmRPPtEta->GetBin(p,e));         // not a bug ;-) 
   }
   
   // 1'-p correction:
   //Double_t oneCosP1nPsiPtEta = 0.;
   
   if(dmPtEta)
   {
    //oneCosP1nPsiPtEta = dReqnPtEta/dmPtEta;
   
    // fill the 2D profile to get the average 1'-p correction for each (pt, eta) bin:
    if(type == "POI")
    { 
     //fCorrectionsCosP1nPsiPtEtaPOI->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,
     //                                    oneCosP1nPsiPtEta,dmPtEta);
    }
    else if(type == "RP")
    {
     //fCorrectionsCosP1nPsiPtEtaRP->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,
     //                                    oneCosP1nPsiPtEta,dmPtEta);
    }
   } // end of if(dmPtEta*dMult-dmPrimePrimePtEta)
   
   /*
   
   // 4'-particle correlation:
   Double_t four1n1n1n1nPtEta = 0.;
   if((dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
       + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.)) // to be improved (introduce a new variable for this expression)
   {
    four1n1n1n1nPtEta = ((pow(dReQ1n,2.)+pow(dImQ1n,2.))*(dReqnPtEta*dReQ1n+dImqnPtEta*dImQ1n)
                      - dReqPrimePrime2nPtEta*(pow(dReQ1n,2.)-pow(dImQ1n,2.))
                      - 2.*dImqPrimePrime2nPtEta*dReQ1n*dImQ1n
                      - dReqnPtEta*(dReQ1n*dReQ2n+dImQ1n*dImQ2n)
                      + dImqnPtEta*(dImQ1n*dReQ2n-dReQ1n*dImQ2n)
                      - 2.*dMult*(dReqnPtEta*dReQ1n+dImqnPtEta*dImQ1n)
                      - 2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*dmPrimePrimePtEta                      
                      + 6.*(dReqPrimePrime1nPtEta*dReQ1n+dImqPrimePrime1nPtEta*dImQ1n)                                            
                      + 1.*(dReqPrimePrime2nPtEta*dReQ2n+dImqPrimePrime2nPtEta*dImQ2n)                      
                      + 2.*(dReqnPtEta*dReQ1n+dImqnPtEta*dImQ1n)                       
                      + 2.*dmPrimePrimePtEta*dMult                      
                      - 6.*dmPrimePrimePtEta)        
                      / ((dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
                          + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.)); 
    
    // fill the 2D profile to get the average correlation for each (pt, eta) bin:
    if(type == "POI")
    {
     f4pPtEtaPOI->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nPtEta,
                       (dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
                        + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.));
    }
    else if(type == "RP")
    {
     f4pPtEtaRP->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nPtEta,
                      (dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
                       + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.));   
    }
   } // end of if((dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
     //            +dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.))
   
  */
   
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateCorrectionsForNonUniformAcceptanceForDifferentialFlowCosTerms(TString type)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCorrectionsForNonUniformAcceptanceForDifferentialFlowSinTerms(TString type)
{
 // calculate corrections for non-uniform acceptance of the detector for differential flow (sin terms)
 
 // multiplicity:
 //Double_t dMult = (*fSMpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 //Double_t dReQ1n = (*fReQ)(0,0);
 //Double_t dReQ2n = (*fReQ)(1,0);
 //Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 //Double_t dImQ1n = (*fImQ)(0,0);
 //Double_t dImQ2n = (*fImQ)(1,0);
 //Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);

 // looping over all (pt,eta) bins and calculating correlations needed for differential flow: 
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   // real and imaginary parts of q_n (non-weighted Q-vector evaluated only for POIs in harmonic n for each (pt,eta) bin): 
   //Double_t dReqnPtEta = 0.;
   //Double_t dImqnPtEta = 0.;

   // number of POIs in each (pt,eta) bin:
   Double_t dmPtEta = 0.;

   // real and imaginary parts of q''_{n}, q''_{2n}, ... 
   // (non-weighted Q-vectors evaluated only for particles which are both RPs and POIs in harmonic n, 2n, ... for each (pt,eta) bin): 
   //Double_t dReqPrimePrime1nPtEta = 0.;
   //Double_t dImqPrimePrime1nPtEta = 0.;
   //Double_t dReqPrimePrime2nPtEta = 0.;
   //Double_t dImqPrimePrime2nPtEta = 0.;

   // number of particles which are both RPs and POIs in each (pt,eta) bin:
   //Double_t dmPrimePrimePtEta = 0.;
   
   if(type == "POI")
   {
    // q''_{n}, q''_{2n}:
    //...............................................................................................
    //dReqPrimePrime1nPtEta = fReqPrimePrime1nPtEta->GetBinContent(fReqPrimePrime1nPtEta->GetBin(p,e));
    //dImqPrimePrime1nPtEta = fImqPrimePrime1nPtEta->GetBinContent(fImqPrimePrime1nPtEta->GetBin(p,e));
    //dReqPrimePrime2nPtEta = fReqPrimePrime2nPtEta->GetBinContent(fReqPrimePrime2nPtEta->GetBin(p,e));
    //dImqPrimePrime2nPtEta = fImqPrimePrime2nPtEta->GetBinContent(fImqPrimePrime2nPtEta->GetBin(p,e));
    //...............................................................................................
   
    // m'':
    //dmPrimePrimePtEta = fmPrimePrimePtEta->GetBinContent(fmPrimePrimePtEta->GetBin(p,e));
   
    // q'_{n}: 
    //dReqnPtEta = fReqnPtEta->GetBinContent(fReqnPtEta->GetBin(p,e));
    //dImqnPtEta = fImqnPtEta->GetBinContent(fImqnPtEta->GetBin(p,e));
    //dmPtEta    = fmPtEta->GetBinContent(fmPtEta->GetBin(p,e));
   }
   else if(type == "RP")
   {
    // q_RP{n}, q_RP{2n}:
    //...............................................................................................
    //dReqPrimePrime1nPtEta = fReqRP1nPtEta->GetBinContent(fReqRP1nPtEta->GetBin(p,e));
    //dImqPrimePrime1nPtEta = fImqRP1nPtEta->GetBinContent(fImqRP1nPtEta->GetBin(p,e));
    //dReqPrimePrime2nPtEta = fReqRP2nPtEta->GetBinContent(fReqRP2nPtEta->GetBin(p,e));
    //dImqPrimePrime2nPtEta = fImqRP2nPtEta->GetBinContent(fImqRP2nPtEta->GetBin(p,e));
    //...............................................................................................
   
    // m'':
    //dmPrimePrimePtEta = fmRPPtEta->GetBinContent(fmRPPtEta->GetBin(p,e));
   
    //dReqnPtEta = fReqRP1nPtEta->GetBinContent(fReqRP1nPtEta->GetBin(p,e)); // not a bug ;-)
    //dImqnPtEta = fImqRP1nPtEta->GetBinContent(fImqRP1nPtEta->GetBin(p,e)); // not a bug ;-)
    //dmPtEta    = fmRPPtEta->GetBinContent(fmRPPtEta->GetBin(p,e));         // not a bug ;-) 
   }
   
   // 1'-p correction:
   //Double_t oneSinP1nPsiPtEta = 0.;
   
   if(dmPtEta)
   {
    //oneSinP1nPsiPtEta = dImqnPtEta/dmPtEta;
   
    // fill the 2D profile to get the average 1'-p correction for each (pt, eta) bin:
    if(type == "POI")
    { 
     //fCorrectionsSinP1nPsiPtEtaPOI->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,
     //                                    oneSinP1nPsiPtEta,dmPtEta);
    }
    else if(type == "RP")
    {
     //fCorrectionsSinP1nPsiPtEtaRP->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,
     //                                    oneSinP1nPsiPtEta,dmPtEta);
    }
   } // end of if(dmPtEta*dMult-dmPrimePrimePtEta)
   
   /*
   
   // 4'-particle correlation:
   Double_t four1n1n1n1nPtEta = 0.;
   if((dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
       + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.)) // to be improved (introduce a new variable for this expression)
   {
    four1n1n1n1nPtEta = ((pow(dReQ1n,2.)+pow(dImQ1n,2.))*(dReqnPtEta*dReQ1n+dImqnPtEta*dImQ1n)
                      - dReqPrimePrime2nPtEta*(pow(dReQ1n,2.)-pow(dImQ1n,2.))
                      - 2.*dImqPrimePrime2nPtEta*dReQ1n*dImQ1n
                      - dReqnPtEta*(dReQ1n*dReQ2n+dImQ1n*dImQ2n)
                      + dImqnPtEta*(dImQ1n*dReQ2n-dReQ1n*dImQ2n)
                      - 2.*dMult*(dReqnPtEta*dReQ1n+dImqnPtEta*dImQ1n)
                      - 2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*dmPrimePrimePtEta                      
                      + 6.*(dReqPrimePrime1nPtEta*dReQ1n+dImqPrimePrime1nPtEta*dImQ1n)                                            
                      + 1.*(dReqPrimePrime2nPtEta*dReQ2n+dImqPrimePrime2nPtEta*dImQ2n)                      
                      + 2.*(dReqnPtEta*dReQ1n+dImqnPtEta*dImQ1n)                       
                      + 2.*dmPrimePrimePtEta*dMult                      
                      - 6.*dmPrimePrimePtEta)        
                      / ((dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
                          + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.)); 
    
    // fill the 2D profile to get the average correlation for each (pt, eta) bin:
    if(type == "POI")
    {
     f4pPtEtaPOI->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nPtEta,
                       (dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
                        + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.));
    }
    else if(type == "RP")
    {
     f4pPtEtaRP->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nPtEta,
                      (dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
                       + dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.));   
    }
   } // end of if((dmPtEta-dmPrimePrimePtEta)*dMult*(dMult-1.)*(dMult-2.)
     //            +dmPrimePrimePtEta*(dMult-1.)*(dMult-2.)*(dMult-3.))
   
  */
   
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateCorrectionsForNonUniformAcceptanceForDifferentialFlowSinTerms(TString type)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::EvaluateNestedLoopsForDifferentialFlow(AliFlowEventSimple* anEvent)
{
 // evaluate the nested loops relevant for differential flow (needed for cross-checking the results)
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 
 Double_t psi1=0., phi2=0., phi3=0., phi4=0.;// phi5=0., phi6=0., phi7=0., phi8=0.;
 Double_t wPhi1=1., wPhi2=1., wPhi3=1., wPhi4=1.;// wPhi5=1., wPhi6=1., wPhi7=1., wPhi8=1.;
 
 Int_t n = fHarmonic; // to be improved
 
 //                                          ********************************************
 //                                          **** NESTED LOOPS FOR DIFFERENTIAL FLOW ****
 //                                          ******************************************** 
 
 // Remark 1: (pt,eta) bin in which the cross-checking will be performed is given by 1.1 < pt < 1.2 GeV and -0.55 < eta < -0.525 
 
 // Remark 2: multi-particle correlations needed for diff. flow calculated with nested loops without weights are stored in 1D profile  
 //           fDirectCorrelationsDiffFlow
 
 // Remark 3: multi-particle correlations needed for diff. flow calculated with nested loops with weights are stored in 1D profile  
 //           fDirectCorrelationsDiffFlowW;
 
 // Remark 4: binning of fDirectCorrelationsDiffFlow is organized as follows:
 //......................................................................................
 //       ---- bins 1-20: 2-particle correlations ----
 //  1st bin: <2'>_{1n|1n} = twoPrime1n1n = <cos(n*(psi1-phi2))>
 //       ---- bins 21-40: 3-particle correlations ----
 //       ---- bins 41-60: 4-particle correlations ----
 // 41st bin: <4'>_{1n,1n|1n,1n} = fourPrime1n1n1n1n  = <cos(n*(psi1+phi2-phi3-phi4))>
 //......................................................................................
 
 // Remark 5: binning of fDirectCorrelationsDiffFlow is organized as follows:
 //......................................................................................
 //       ---- bins 1-20: 2-particle correlations ----
 //  1st bin: twoPrime1n1nW0W1 = <w2 cos(n*(psi1-phi2))>
 //       ---- bins 21-40: 3-particle correlations ----
 //       ---- bins 41-60: 4-particle correlations ----
 // 41st bin: fourPrime1n1n1n1nW0W1W1W1 = <w2 w3 w4 cos(n*(psi1+phi2-phi3-phi4))>
 //......................................................................................
 
 // 2'-particle:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): 
  if(!((aftsTrack->Pt()>=1.1 && aftsTrack->Pt()<1.2) && (aftsTrack->Eta()>=-0.55 && aftsTrack->Eta()<-0.525) && (aftsTrack->InPOISelection())))continue;
  psi1=aftsTrack->Phi(); 
  if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(psi1*fnBinsPhi/TMath::TwoPi())));
  
  fDirectCorrectionsDiffFlowCos->Fill(0.,cos(1.*n*(psi1)),1.); // <<cos(n*(psi1))>>
  fDirectCorrectionsDiffFlowSin->Fill(0.,sin(1.*n*(psi1)),1.); // <<sin(n*(psi1))>>
  
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP):
   if(!(aftsTrack->InRPSelection()))continue;
   phi2=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));
    
   // non-weighted: 
   //.....................................................................................  
   fDirectCorrelationsDiffFlow->Fill(0.,cos(1.*n*(psi1-phi2)),1.); // <cos(n*(psi1-phi2))>
   //.....................................................................................  
   // weighted:
   //.....................................................................................   
   if(fUsePhiWeights) fDirectCorrelationsDiffFlowW->Fill(0.,cos(1.*n*(psi1-phi2)),wPhi2); // <w2 cos(n*(psi1-phi2))>
   //.....................................................................................  
   
   //fDirectCorrelations->Fill(103.,cos(1.*n*(phi1-phi2)),pow(wPhi1,2)*wPhi2);//<2'>_{n,n}
   //fDirectCorrelations->Fill(104.,cos(2.*n*(phi1-phi2)),wPhi1*pow(wPhi2,2));//<2'>_{n,n}
   //fDirectCorrelations->Fill(105.,cos(1.*n*(phi1-phi2)),pow(wPhi2,3));//<2'>_{n,n}  
   //fDirectCorrelations->Fill(41.,cos(2.*n*(phi1-phi2)),1);//<2'>_{2n,2n}
   //fDirectCorrelations->Fill(42.,cos(3.*n*(phi1-phi2)),1);//<2'>_{3n,3n}
   //fDirectCorrelations->Fill(43.,cos(4.*n*(phi1-phi2)),1);//<2'>_{4n,4n}   
    
  }//end of for(Int_t i2=0;i2<nPrim;i2++)
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 
 
 /*
 
 //<3'>_{2n|n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  if(!((aftsTrack->Pt()>=0.5&&aftsTrack->Pt()<0.6)&&(aftsTrack->InPOISelection())))continue;//POI condition
  psi1=aftsTrack->Phi();
  if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(psi1*fnBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   aftsTrack=anEvent->GetTrack(i2);
   if(!(aftsTrack->InRPSelection()))continue;//RP condition
   phi2=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2)continue;
    aftsTrack=anEvent->GetTrack(i3);
    if(!(aftsTrack->InRPSelection()))continue;//RP condition
    phi3=aftsTrack->Phi();
    if(fUsePhiWeights && fPhiWeights) wPhi3 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*fnBinsPhi/TMath::TwoPi())));
    //fill the fDirectCorrelations:     
    
    // 2-p
    //fDirectCorrelations->Fill(101.,cos(n*(phi2-phi3)),wPhi1*wPhi2*wPhi3); // <w1 w2 w3 cos(n(phi2-phi3))>
    //fDirectCorrelations->Fill(102.,cos(n*(phi1-phi3)),pow(wPhi2,2.)*wPhi3); // <w2^2 w3 cos(n(psi1-phi2))>
    
    // 3-p            
    //fDirectCorrelations->Fill(110.,cos(n*(2.*phi1-phi2-phi3)),wPhi1*wPhi2*wPhi3); // <w1 w2 w3 cos(n(2psi1-phi2-phi3))>
    //fDirectCorrelations->Fill(111.,cos(n*(phi1+phi2-2.*phi3)),wPhi2*pow(wPhi3,2.)); // <w2 w3^2 cos(n(psi1+phi2-2.*phi3))>
    
    
    //fDirectCorrelations->Fill(46.,cos(n*(phi1+phi2-2.*phi3)),1);//<3'>_{n,n|2n}    
   }//end of for(Int_t i3=0;i3<nPrim;i3++)  
  }//end of for(Int_t i2=0;i2<nPrim;i2++)  
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 */
 
 // 4'-particle:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // POI condition (first particle in the correlator must be POI): 
  if(!((aftsTrack->Pt()>=1.1 && aftsTrack->Pt()<1.2) && (aftsTrack->Eta()>=-0.55 && aftsTrack->Eta()<-0.525) && (aftsTrack->InPOISelection())))continue;
  psi1=aftsTrack->Phi();
  if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(psi1*fnBinsPhi/TMath::TwoPi())));
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   aftsTrack=anEvent->GetTrack(i2);
   // RP condition (!(first) particle in the correlator must be RP): 
   if(!(aftsTrack->InRPSelection()))continue;
   phi2=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    aftsTrack=anEvent->GetTrack(i3);
    // RP condition (!(first) particle in the correlator must be RP):
    if(!(aftsTrack->InRPSelection()))continue;
    phi3=aftsTrack->Phi();
    if(fUsePhiWeights && fPhiWeights) wPhi3 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi3*fnBinsPhi/TMath::TwoPi())));
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     aftsTrack=anEvent->GetTrack(i4);
     // RP condition (!(first) particle in the correlator must be RP):
     if(!(aftsTrack->InRPSelection()))continue;  
     phi4=aftsTrack->Phi();
     if(fUsePhiWeights && fPhiWeights) wPhi4 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi4*fnBinsPhi/TMath::TwoPi())));
     
     // non-weighted:
     //.........................................................................................................................
     fDirectCorrelationsDiffFlow->Fill(40.,cos(n*(psi1+phi2-phi3-phi4)),1.); // <cos(n(psi1+phi1-phi2-phi3))> 
     //.........................................................................................................................     
     // weighted:
     //...............................................................................................................................
     if(fUsePhiWeights) fDirectCorrelationsDiffFlowW->Fill(40.,cos(n*(psi1+phi2-phi3-phi4)),wPhi2*wPhi3*wPhi4); // <w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))> 
     //...............................................................................................................................     
          
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 /*                
 //<5'>_{2n,n|n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  if(!((aftsTrack->Pt()>=0.5&&aftsTrack->Pt()<0.6)&&(aftsTrack->InPOISelection())))continue;//POI condition
  phi1=aftsTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   aftsTrack=anEvent->GetTrack(i2);
   if(!(aftsTrack->InRPSelection()))continue;//RP condition   
   phi2=aftsTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    aftsTrack=anEvent->GetTrack(i3);
    if(!(aftsTrack->InRPSelection()))continue;//RP condition   
    phi3=aftsTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     aftsTrack=anEvent->GetTrack(i4);
     if(!(aftsTrack->InRPSelection()))continue;//RP condition  
     phi4=aftsTrack->Phi();//
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      aftsTrack=anEvent->GetTrack(i5);
      if(!(aftsTrack->InRPSelection()))continue;//RP condition  
      phi5=aftsTrack->Phi();    
      //fill the fDirectCorrelations:if(bNestedLoops)
      //fDirectCorrelations->Fill(55.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),1);//<5'>_{2n,n|n,n,n}
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 

  
 */
 /*
 
 
 
 //<6'>_{n,n,n|n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  if(!((aftsTrack->Pt()>=0.5&&aftsTrack->Pt()<0.6)&&(aftsTrack->InPOISelection())))continue;//POI condition
  phi1=aftsTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   aftsTrack=anEvent->GetTrack(i2);
   if(!(aftsTrack->InRPSelection()))continue;//RP condition   
   phi2=aftsTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    aftsTrack=anEvent->GetTrack(i3);
    if(!(aftsTrack->InRPSelection()))continue;//RP condition   
    phi3=aftsTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     aftsTrack=anEvent->GetTrack(i4);
     if(!(aftsTrack->InRPSelection()))continue;//RP condition  
     phi4=aftsTrack->Phi();
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      aftsTrack=anEvent->GetTrack(i5);
      if(!(aftsTrack->InRPSelection()))continue;//RP condition  
      phi5=aftsTrack->Phi();    
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       aftsTrack=anEvent->GetTrack(i6);
       if(!(aftsTrack->InRPSelection()))continue;//RP condition  
       phi6=aftsTrack->Phi();  
       //fill the fDirectCorrelations:
       //fDirectCorrelations->Fill(60.,cos(n*(phi1+phi2+phi3-phi4-phi5-phi6)),1);//<6'>_{n,n,n|n,n,n}
      }//end of for(Int_t i6=0;i6<nPrim;i6++)   
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)

 
 */
 /* 
   
     
 //<7'>_{2n,n,n|n,n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  if(!((aftsTrack->Pt()>=0.5&&aftsTrack->Pt()<0.6)&&(aftsTrack->InPOISelection())))continue;//POI condition
  phi1=aftsTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   aftsTrack=anEvent->GetTrack(i2);
   if(!(aftsTrack->InRPSelection()))continue;//RP condition   
   phi2=aftsTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    aftsTrack=anEvent->GetTrack(i3);
    if(!(aftsTrack->InRPSelection()))continue;//RP condition   
    phi3=aftsTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     aftsTrack=anEvent->GetTrack(i4);
     if(!(aftsTrack->InRPSelection()))continue;//RP condition  
     phi4=aftsTrack->Phi();
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      aftsTrack=anEvent->GetTrack(i5);
      if(!(aftsTrack->InRPSelection()))continue;//RP condition  
      phi5=aftsTrack->Phi();    
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       aftsTrack=anEvent->GetTrack(i6);
       if(!(aftsTrack->InRPSelection()))continue;//RP condition  
       phi6=aftsTrack->Phi();
       for(Int_t i7=0;i7<nPrim;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        aftsTrack=anEvent->GetTrack(i7);
        if(!(aftsTrack->InRPSelection()))continue;//RP condition  
        phi7=aftsTrack->Phi();   
        //fill the fDirectCorrelations:
        //fDirectCorrelations->Fill(65.,cos(2.*n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6-n*phi7),1);//<7'>_{2n,n,n|n,n,n,n}
       }//end of for(Int_t i7=0;i7<nPrim;i7++)  
      }//end of for(Int_t i6=0;i6<nPrim;i6++)   
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)

 
  
 */
 /*  
    
     
       
 //<8'>_{n,n,n,n|n,n,n,n}
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  if(!((aftsTrack->Pt()>=0.5&&aftsTrack->Pt()<0.6)&&(aftsTrack->InPOISelection())))continue;//POI condition
  phi1=aftsTrack->Phi();
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   aftsTrack=anEvent->GetTrack(i2);
   if(!(aftsTrack->InRPSelection()))continue;//RP condition   
   phi2=aftsTrack->Phi();
   for(Int_t i3=0;i3<nPrim;i3++)
   { 
    if(i3==i1||i3==i2)continue;
    aftsTrack=anEvent->GetTrack(i3);
    if(!(aftsTrack->InRPSelection()))continue;//RP condition   
    phi3=aftsTrack->Phi();
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3)continue;
     aftsTrack=anEvent->GetTrack(i4);
     if(!(aftsTrack->InRPSelection()))continue;//RP condition  
     phi4=aftsTrack->Phi();
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      aftsTrack=anEvent->GetTrack(i5);
      if(!(aftsTrack->InRPSelection()))continue;//RP condition  
      phi5=aftsTrack->Phi();    
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       aftsTrack=anEvent->GetTrack(i6);
       if(!(aftsTrack->InRPSelection()))continue;//RP condition  
       phi6=aftsTrack->Phi();
       for(Int_t i7=0;i7<nPrim;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        aftsTrack=anEvent->GetTrack(i7);
        if(!(aftsTrack->InRPSelection()))continue;//RP condition  
        phi7=aftsTrack->Phi();
        for(Int_t i8=0;i8<nPrim;i8++)
        {
         if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7)continue;
         aftsTrack=anEvent->GetTrack(i8);
         if(!(aftsTrack->InRPSelection()))continue;//RP condition  
         phi8=aftsTrack->Phi();           
         //fill the fDirectCorrelations:
         //fDirectCorrelations->Fill(70.,cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)),1);//<8'>_{n,n,n,n|n,n,n,n}
        }//end of for(Int_t i8=0;i8<nPrim;i8++) 
       }//end of for(Int_t i7=0;i7<nPrim;i7++)  
      }//end of for(Int_t i6=0;i6<nPrim;i6++)   
     }//end of for(Int_t i5=0;i5<nPrim;i5++)  
    }//end of for(Int_t i4=0;i4<nPrim;i4++)
   }//end of for(Int_t i3=0;i3<nPrim;i3++)
  }//end of for(Int_t i2=0;i2<nPrim;i2++) 
 }//end of for(Int_t i1=0;i1<nPrim;i1++)
 
 
 
 */ 
 
 
 
 
} // end of AliFlowAnalysisWithQCumulants::EvaluateNestedLoopsForDifferentialFlow(AliFlowEventSimple* anEvent)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::GetOutputHistograms(TList *outputListHistos)
{
 // get pointers to all output histograms (called before Finish())
 if(outputListHistos)
 {	
  // 1.) common control histograms and common histograms for final results:
  AliFlowCommonHist *commonHist = dynamic_cast<AliFlowCommonHist*>(outputListHistos->FindObject("AliFlowCommonHistQC"));
  if(commonHist) this->SetCommonHists(commonHist); 
  AliFlowCommonHist *commonHist2nd = dynamic_cast<AliFlowCommonHist*>(outputListHistos->FindObject("AliFlowCommonHist2ndOrderQC"));
  if(commonHist2nd) this->SetCommonHists2nd(commonHist2nd); 
  AliFlowCommonHist *commonHist4th = dynamic_cast<AliFlowCommonHist*>(outputListHistos->FindObject("AliFlowCommonHist4thOrderQC"));
  if(commonHist4th) this->SetCommonHists4th(commonHist4th);
  AliFlowCommonHist *commonHist6th = dynamic_cast<AliFlowCommonHist*>(outputListHistos->FindObject("AliFlowCommonHist6thOrderQC"));
  if(commonHist6th) this->SetCommonHists6th(commonHist6th);
  AliFlowCommonHist *commonHist8th = dynamic_cast<AliFlowCommonHist*>(outputListHistos->FindObject("AliFlowCommonHist8thOrderQC"));
  if(commonHist8th) this->SetCommonHists8th(commonHist8th);
  AliFlowCommonHistResults *commonHistRes2nd = dynamic_cast<AliFlowCommonHistResults*>
                                               (outputListHistos->FindObject("AliFlowCommonHistResults2ndOrderQC"));
  if(commonHistRes2nd) this->SetCommonHistsResults2nd(commonHistRes2nd); 
  AliFlowCommonHistResults *commonHistRes4th = dynamic_cast<AliFlowCommonHistResults*>
                                               (outputListHistos->FindObject("AliFlowCommonHistResults4thOrderQC"));
  if(commonHistRes4th) this->SetCommonHistsResults4th(commonHistRes4th);
  AliFlowCommonHistResults *commonHistRes6th = dynamic_cast<AliFlowCommonHistResults*>
                                               (outputListHistos->FindObject("AliFlowCommonHistResults6thOrderQC"));
  if(commonHistRes6th) this->SetCommonHistsResults6th(commonHistRes6th);
  AliFlowCommonHistResults *commonHistRes8th = dynamic_cast<AliFlowCommonHistResults*>
                                               (outputListHistos->FindObject("AliFlowCommonHistResults8thOrderQC"));  
  if(commonHistRes8th) this->SetCommonHistsResults8th(commonHistRes8th);
  
  // 2.) weights: 
  TList *weightsList = dynamic_cast<TList*>(outputListHistos->FindObject("Weights"));
  if(weightsList) this->SetWeightsList(weightsList);
  Bool_t bUsePhiWeights = kFALSE;
  Bool_t bUsePtWeights = kFALSE;
  Bool_t bUseEtaWeights = kFALSE;
  TString fUseParticleWeightsName = "fUseParticleWeightsQC";
  fUseParticleWeightsName += fAnalysisLabel->Data();
  TProfile *useParticleWeights = dynamic_cast<TProfile*>(weightsList->FindObject(fUseParticleWeightsName.Data()));
  if(useParticleWeights)
  {
   this->SetUseParticleWeights(useParticleWeights);  
   bUsePhiWeights = (Int_t)useParticleWeights->GetBinContent(1);
   bUsePtWeights = (Int_t)useParticleWeights->GetBinContent(2);
   bUseEtaWeights = (Int_t)useParticleWeights->GetBinContent(3);
  }
  
  // 3.) integrated flow:
  TList *intFlowList = NULL;
  TList *intFlowProfiles = NULL;
  TList *intFlowResults = NULL;
   
  intFlowList = dynamic_cast<TList*>(outputListHistos->FindObject("Integrated Flow"));
  if(intFlowList) 
  {
   intFlowProfiles = dynamic_cast<TList*>(intFlowList->FindObject("Profiles"));
   intFlowResults = dynamic_cast<TList*>(intFlowList->FindObject("Results"));
  } else
    {
     cout<<"WARNING: intFlowList is NULL in AFAWQC::GOH() !!!!"<<endl;
    }  
   
  // profiles:   
  if(intFlowProfiles)  
  {
   // average multiplicities:
   TProfile *avMultiplicity = dynamic_cast<TProfile*>(intFlowProfiles->FindObject("fAvMultiplicity"));
   if(avMultiplicity) 
   {
    this->SetAvMultiplicity(avMultiplicity);
   } else 
     {
      cout<<"WARNING: avMultiplicity is NULL in AFAWQC::GOH() !!!!"<<endl;
     }
   
   // flags: (to be improved (united with other flags in this method))
   TString pWeightsFlag[2] = {"pWeights not used","pWeights used"};
   TString eWeightsFlag[2] = {"exact eWeights","non-exact eWeights"};
   //TString nuaFlag[2] = {"not corrected","corrected"};
   TString sinCosFlag[2] = {"sin","cos"};
    
   for(Int_t pW=0;pW<1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights);pW++)
   {
    for(Int_t eW=0;eW<1;eW++)
    {
     // correlations (profiles):
     TProfile *qCorrelations = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(Form("fQCorrelations: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
     if(qCorrelations) 
     {
      this->SetQCorrelations(qCorrelations,pW,eW);
     } else 
       {
        cout<<"WARNING: qCorrelations is NULL in AFAWQC::GOH() !!!!"<<endl;
        cout<<"pW = "<<pW<<endl;
        cout<<"eW = "<<eW<<endl;
       } 
     // products (profiles):  
     TProfile *qProducts = dynamic_cast<TProfile*>(intFlowProfiles->FindObject(Form("fQProducts: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
     if(qProducts) 
     {
      this->SetQProducts(qProducts,pW,eW);
     } else 
       {
        cout<<"WARNING: qProducts is NULL in AFAWQC::GOH() !!!!"<<endl;
        cout<<"pW = "<<pW<<endl;
        cout<<"eW = "<<eW<<endl;
       }     
     // corrections (profiles):
     for(Int_t sc=0;sc<2;sc++)
     {
      TProfile *qCorrections = dynamic_cast<TProfile*>(intFlowProfiles->FindObject((Form("fQCorrections: %s, %s, %s terms",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),sinCosFlag[sc].Data()))));
      if(qCorrections) 
      {
       this->SetQCorrections(qCorrections,pW,eW,sc);
      } else 
        {
         cout<<"WARNING: qCorrections is NULL in AFAWQC::GOH() !!!!"<<endl;
         cout<<"pW = "<<pW<<endl;
         cout<<"eW = "<<eW<<endl;
         cout<<"sc = "<<sc<<endl;
        } 
     } // end of for(Int_t sc=0;sc<2;sc++)           
    } // end of for(Int_t eW=0;eW<1;eW++)
   } // end of for(Int_t pW=0;pW<1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights);pW++)
  } else // to if(intFlowProfiles)  
    {
     cout<<"WARNING: intFlowProfiles is NULL in AFAWQC::GOH() !!!!"<<endl;
    }
   
  // results:   
  if(intFlowResults)  
  {
   for(Int_t pW=0;pW<1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights);pW++)
   {
    for(Int_t eW=0;eW<2;eW++)
    {
     // flags: (to be improved (united with other flags in this method))
     TString pWeightsFlag[2] = {"pWeights not used","pWeights used"};
     TString eWeightsFlag[2] = {"exact eWeights","non-exact eWeights"};
     TString nuaFlag[2] = {"not corrected","corrected"};
     TString powerFlag[2] = {"linear","quadratic"};
     //TString sinCosFlag[2] = {"sin","cos"};     
     // correlations (results):
     TH1D *correlations = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("fCorrelations: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
     if(correlations) 
     {
      this->SetCorrelations(correlations,pW,eW);
     } else 
       {
        cout<<"WARNING: correlations is NULL in AFAWQC::GOH() !!!!"<<endl;
        cout<<"pW = "<<pW<<endl;
        cout<<"eW = "<<eW<<endl;
       }     
     // corrections (results):
     TH1D *corrections = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("fCorrections: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
     if(corrections) 
     {
      this->SetCorrections(corrections,pW,eW);
     } else 
       {
        cout<<"WARNING: corrections is NULL in AFAWQC::GOH() !!!!"<<endl;
        cout<<"pW = "<<pW<<endl;
        cout<<"eW = "<<eW<<endl;
       }
     // covariances (results):
     TH1D *covariances = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("fCovariances: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
     if(covariances) 
     {
      this->SetCovariances(covariances,pW,eW);
     } else 
       {
        cout<<"WARNING: covariances is NULL in AFAWQC::GOH() !!!!"<<endl;
        cout<<"pW = "<<pW<<endl;
        cout<<"eW = "<<eW<<endl;
       } 
     // sum of linear and quadratic event weights (results):
     for(Int_t power=0;power<2;power++)
     {
      TH1D *sumOfEventWeights = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("fSumOfEventWeights: %s, %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),powerFlag[power].Data())));
      if(sumOfEventWeights) 
      {
       this->SetSumOfEventWeights(sumOfEventWeights,pW,eW,power);
      } else 
        {
         cout<<"WARNING: sumOfEventWeights is NULL in AFAWQC::GOH() !!!!"<<endl;
         cout<<"pW    = "<<pW<<endl;
         cout<<"eW    = "<<eW<<endl;
         cout<<"power = "<<power<<endl;
        }                                   
     } // end of for(Int_t power=0;power<2;power++)                                                                  
     // products of event weights (results):
     TH1D *productOfEventWeights = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("fProductOfEventWeights: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
     if(productOfEventWeights) 
     {
      this->SetProductOfEventWeights(productOfEventWeights,pW,eW);
     } else 
       {
        cout<<"WARNING: productOfEventWeights is NULL in AFAWQC::GOH() !!!!"<<endl;
        cout<<"pW = "<<pW<<endl;
        cout<<"eW = "<<eW<<endl;
       } 
       
     for(Int_t nua=0;nua<2;nua++)
     {
      // integrated Q-cumulants:
      TH1D *cumulants = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("fCumulants: %s, %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())));
      if(cumulants) 
      {
       this->SetCumulants(cumulants,pW,eW,nua);
      } else 
        {
         cout<<"WARNING: cumulants is NULL in AFAWQC::GOH() !!!!"<<endl;
         cout<<"pW = "<<pW<<endl;
         cout<<"eW = "<<eW<<endl;
         cout<<"nua = "<<nua<<endl;
        }  
      // integrated flow estimates from Q-cumulants:
      TH1D *intFlow = dynamic_cast<TH1D*>(intFlowResults->FindObject(Form("fIntFlow: %s, %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())));
      if(intFlow) 
      {
       this->SetIntFlow(intFlow,pW,eW,nua);
      } else 
        {
         cout<<"WARNING: intFlow is NULL in AFAWQC::GOH() !!!!"<<endl;
         cout<<"pW = "<<pW<<endl;
         cout<<"eW = "<<eW<<endl;
         cout<<"nua = "<<nua<<endl;
        }   
     } // end of for(Int_t nua=0;nua<2;nua++)  
    } // end of for(Int_t eW=0;eW<1;eW++) 
   } // end of for(Int_t pW=0;pW<1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights);pW++)  
  } else // to if(intFlowResults)
    {
     cout<<"WARNING: intFlowResults is NULL in AFAWQC::GOH() !!!!"<<endl;
    }
          
  // differential flow:
  TString typeFlag[2] = {"RP","POI"}; 
  TString pWeightsFlag[2] = {"not used","used"};
  TString eWeightsFlag[2] = {"exact","non-exact"}; 
  TString sinCosFlag[2] = {"sin","cos"};
  TString nuaFlag[2] = {"not corrected","corrected"}; // nua = non-uniform acceptance
  TString ptEtaFlag[2] = {"p_{t}","#eta"};
  // base list fDiffFlowList "Differential Flow":
  TList *diffFlowList = NULL;
  diffFlowList = dynamic_cast<TList*>(outputListHistos->FindObject("Differential Flow"));  
  // list holding nested lists containing profiles:
  TList *diffFlowListProfiles = NULL;
  // list holding nested lists containing 2D and 1D histograms with final results:
  TList *diffFlowListResults = NULL;
  if(diffFlowList)
  {  
   diffFlowListProfiles = dynamic_cast<TList*>(diffFlowList->FindObject("Profiles"));
   diffFlowListResults = dynamic_cast<TList*>(diffFlowList->FindObject("Results"));
  } else
    {
     cout<<"WARNING: diffFlowList is NULL in AFAWQC::GOH() !!!!"<<endl;
    }      
  
  // nested list in the list of profiles fDiffFlowListProfiles "Profiles":
  TList *dfpType[2] = {NULL};
  TList *dfpParticleWeights[2][2] = {{NULL}};
  TList *dfpEventWeights[2][2][2] = {{{NULL}}};
  TList *diffFlowCorrelations[2][2][2] = {{{NULL}}};
  TList *diffFlowProductsOfCorrelations[2][2][2] = {{{NULL}}};
  TList *diffFlowCorrectionTerms[2][2][2][2] = {{{{NULL}}}};
  
  if(diffFlowListProfiles)
  {
   for(Int_t t=0;t<2;t++)
   {
    dfpType[t] = dynamic_cast<TList*>(diffFlowListProfiles->FindObject(typeFlag[t].Data()));
    if(dfpType[t])
    {
     for(Int_t pW=0;pW<(1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights));pW++)
     {
      dfpParticleWeights[t][pW] = dynamic_cast<TList*>(dfpType[t]->FindObject(Form("%s, pWeights %s",typeFlag[t].Data(),pWeightsFlag[pW].Data())));
      if(dfpParticleWeights[t][pW])
      {
       for(Int_t eW=0;eW<2;eW++)
       {
        dfpEventWeights[t][pW][eW] = dynamic_cast<TList*>(dfpParticleWeights[t][pW]->FindObject(Form("%s, pWeights %s, eWeights %s",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
        if(dfpEventWeights[t][pW][eW])
        {
         // correlations: 
         diffFlowCorrelations[t][pW][eW] = dynamic_cast<TList*>(dfpEventWeights[t][pW][eW]->FindObject(Form("Correlations (%s, pWeights %s, eWeights %s)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
         // products of correlations:
         diffFlowProductsOfCorrelations[t][pW][eW] = dynamic_cast<TList*>(dfpEventWeights[t][pW][eW]->FindObject(Form("Products of correlations (%s, pWeights %s, eWeights %s)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
         // correction terms:
         for(Int_t sc=0;sc<2;sc++)
         {
          diffFlowCorrectionTerms[t][pW][eW][sc] = dynamic_cast<TList*>(dfpEventWeights[t][pW][eW]->FindObject(Form("Corrections for NUA (%s, pWeights %s, eWeights %s, %s terms)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),sinCosFlag[sc].Data())));
          //this->SetDiffFlowCorrectionTerms(diffFlowCorrectionTerms[t][pW][sc],t,pW,sc);   
         }
        } else // to if(dfpEventWeights[t][pW][eW])
          {
           cout<<"WARNING: dfpEventWeights[t][pW][eW] is NULL in AFAWQC::GOH() !!!!"<<endl;
           cout<<"t = "<<t<<endl;
           cout<<"pW = "<<pW<<endl;
           cout<<"eW = "<<eW<<endl;
          }
       } // end of for(Int_t eW=0;eW<2;eW++)   
      } else // to if(dfpParticleWeights[t][pW])
        {
         cout<<"WARNING: dfpParticleWeights[t][pW] is NULL in AFAWQC::GOH() !!!!"<<endl;
         cout<<"t = "<<t<<endl;
         cout<<"pW = "<<pW<<endl;
        }
     } // end of for(Int_t pW=0;pW<(1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights));pW++) 
     } else // if(dfpType[t]) 
      {
       cout<<"WARNING: dfpType[t] is NULL in AFAWQC::GOH() !!!!"<<endl;
       cout<<"t = "<<t<<endl;
      }
   } // end of for(Int_t t=0;t<2;t++)
  } else // to if(diffFlowListProfiles)
    {
     cout<<"WARNING: diffFlowListProfiles is NULL in AFAWQC::GOH() !!!!"<<endl;
    }  
  
  TProfile2D *correlationsPro[2][2][2][4] = {{{{NULL}}}};
  TProfile2D *productsOfCorrelationsPro[2][2][2][5] = {{{{NULL}}}};
  TProfile2D *correctionTermsPro[2][2][2][2][2] = {{{{{NULL}}}}};
  
  TString correlationName[4] = {"<2'>","<4'>","<6'>","<8'>"};
  TString cumulantName[4] = {"QC{2'}","QC{4'}","QC{6'}","QC{8'}"};
  TString productOfCorrelationsName[5] = {"<2><2'>","<2><4'>","<4><2'>","<4><4'>","<2'><4'>"};
  TString correctionsSinTermsName[2] = {"sin(1)","sin(2)"};
  TString correctionsCosTermsName[2] = {"cos(1)","cos(2)"};
  TString correctionName[4] = {"corr. to QC{2'}","corr. to QC{4'}","corr. to QC{6'}","corr. to QC{8'}"};
  TString covarianceName[4] = {"Cov(2,2')","Cov(2,4')","Cov(2',4')","Cov(4,2')"}; // to be improved (reorganized)
  TString flowName[4] = {"v'{2}","v'{4}","v'{6}","v'{8}"}; 
  
  for(Int_t t=0;t<2;t++)
  {
   for(Int_t pW=0;pW<(1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights));pW++)
   {
    for(Int_t eW=0;eW<2;eW++)
    {
     // correlations:
     if(diffFlowCorrelations[t][pW][eW])
     {
      for(Int_t correlationIndex=0;correlationIndex<4;correlationIndex++)
      {
       correlationsPro[t][pW][eW][correlationIndex] = dynamic_cast<TProfile2D*>(diffFlowCorrelations[t][pW][eW]->FindObject(correlationName[correlationIndex].Data()));
       if(correlationsPro[t][pW][eW][correlationIndex])
       {
        this->SetCorrelationsPro(correlationsPro[t][pW][eW][correlationIndex],t,pW,eW,correlationIndex);
       } else // to if(correlationsPro[t][pW][ew][correlationIndex])
         {
          cout<<"WARNING: correlationsPro[t][pW][eW][correlationIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
          cout<<"t = "<<t<<endl;
          cout<<"pW = "<<pW<<endl;
          cout<<"eW = "<<eW<<endl;
          cout<<"ci = "<<correlationIndex<<endl;
         } 
      }
     } else // to if(diffFlowCorrelations[t][pW][eW])
       {
        cout<<"WARNING: diffFlowCorrelations[t][pW][eW] is NULL in AFAWQC::GOH() !!!!"<<endl;
        cout<<"t = "<<t<<endl;
        cout<<"pW = "<<pW<<endl;
        cout<<"eW = "<<eW<<endl;
       } 
     // products of correlations:
     if(diffFlowProductsOfCorrelations[t][pW][eW])
     {
      for(Int_t productOfCorrelationsIndex=0;productOfCorrelationsIndex<5;productOfCorrelationsIndex++)
      {
       productsOfCorrelationsPro[t][pW][eW][productOfCorrelationsIndex] = dynamic_cast<TProfile2D*>(diffFlowProductsOfCorrelations[t][pW][eW]->FindObject(productOfCorrelationsName[productOfCorrelationsIndex].Data()));
       if(productsOfCorrelationsPro[t][pW][eW][productOfCorrelationsIndex])
       {
        this->SetProductsOfCorrelationsPro(productsOfCorrelationsPro[t][pW][eW][productOfCorrelationsIndex],t,pW,eW,productOfCorrelationsIndex);
       } else // to if(productsOfCorrelationsPro[t][pW][eW][productOfCorrelationsIndex])
         {
          cout<<"WARNING: productsOfCorrelationsPro[t][pW][eW][productOfCorrelationsIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
          cout<<"t = "<<t<<endl;
          cout<<"pW = "<<pW<<endl;
          cout<<"eW = "<<eW<<endl;
          cout<<"ci = "<<productOfCorrelationsIndex<<endl;
         } 
      }
     } else // to if(diffFlowProductsOfCorrelations[t][pW][eW])
       {
        cout<<"WARNING: diffFlowProductsOfCorrelations[t][pW][eW] is NULL in AFAWQC::GOH() !!!!"<<endl;
        cout<<"t = "<<t<<endl;
        cout<<"pW = "<<pW<<endl;
        cout<<"eW = "<<eW<<endl;
       }
     // correction terms:
     for(Int_t sc=0;sc<2;sc++)
     {
      if(diffFlowCorrectionTerms[t][pW][eW][sc])
      {
       for(Int_t correctionIndex=0;correctionIndex<2;correctionIndex++)
       {
        if(sc==0)
        {
         correctionTermsPro[t][pW][eW][sc][correctionIndex] = dynamic_cast<TProfile2D*>(diffFlowCorrectionTerms[t][pW][eW][sc]->FindObject(correctionsSinTermsName[correctionIndex].Data())); 
         if(correctionTermsPro[t][pW][eW][sc][correctionIndex])
         {
          this->SetCorrectionTermsPro(correctionTermsPro[t][pW][eW][sc][correctionIndex],t,pW,eW,sc,correctionIndex);
         } else 
           {
            cout<<"WARNING: correctionTermsPro[t][pW][eW][sc][correctionIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
            cout<<"t = "<<t<<endl;
            cout<<"pW = "<<pW<<endl;
            cout<<"eW = "<<eW<<endl;
            cout<<"sc = "<<sc<<endl;
            cout<<"ci = "<<correctionIndex<<endl;
           }
        } 
        if(sc==1)
        {
         correctionTermsPro[t][pW][eW][sc][correctionIndex] = dynamic_cast<TProfile2D*>(diffFlowCorrectionTerms[t][pW][eW][sc]->FindObject(correctionsCosTermsName[correctionIndex].Data())); 
         if(correctionTermsPro[t][pW][eW][sc][correctionIndex])
         {
          this->SetCorrectionTermsPro(correctionTermsPro[t][pW][eW][sc][correctionIndex],t,pW,eW,sc,correctionIndex);
         } else 
           {
            cout<<"WARNING: correctionTermsPro[t][pW][eW][sc][correctionIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
            cout<<"t = "<<t<<endl;
            cout<<"pW = "<<pW<<endl;
            cout<<"eW = "<<eW<<endl;
            cout<<"sc = "<<sc<<endl;
            cout<<"ci = "<<correctionIndex<<endl;
           }         
        }  
       } // end of for(Int_t correctionIndex=0;correctionIndex<2;correctionIndex++)
      } else // to if(diffFlowCorrectionTerms[t][pW][eW][sc])
        {
         cout<<"WARNING: diffFlowCorrectionTerms[t][pW][eW][sc] is NULL in AFAWQC::GOH() !!!!"<<endl;
         cout<<"t = "<<t<<endl;
         cout<<"pW = "<<pW<<endl;
         cout<<"eW = "<<eW<<endl;
         cout<<"sc = "<<sc<<endl;
        }
     } // end of for(Int_t sc=0;sc<2;sc++)  
    } // end of for(Int_t eW=0;eW<2;eW++)
   }  // end of for(Int_t pW=0;pW<(1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights));pW++)
  } // end of for(Int_t t=0;t<2;t++)
  
  // nested list in the list of results fDiffFlowListResults "Results":
  TList *dfrType[2] = {NULL};
  TList *dfrParticleWeights[2][2] = {{NULL}};
  TList *dfrEventWeights[2][2][2] = {{{NULL}}};
  TList *dfrCorrections[2][2][2][2] = {{{{NULL}}}}; 
  TList *diffFlowFinalCorrelations[2][2][2] = {{{NULL}}}; 
  TList *diffFlowFinalCorrections[2][2][2] = {{{NULL}}}; 
  TList *diffFlowFinalCovariances[2][2][2] = {{{NULL}}};   
  TList *diffFlowFinalCumulants[2][2][2][2] = {{{{NULL}}}};  
  TList *diffFlowFinalFlow[2][2][2][2] = {{{{NULL}}}}; 
 
  if(diffFlowListResults)
  {
   for(Int_t t=0;t<2;t++)
   {
    dfrType[t] = dynamic_cast<TList*>(diffFlowListResults->FindObject(typeFlag[t].Data()));
    if(dfrType[t])
    {
     for(Int_t pW=0;pW<(1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights));pW++)
     {
      dfrParticleWeights[t][pW] = dynamic_cast<TList*>(dfrType[t]->FindObject(Form("%s, pWeights %s",typeFlag[t].Data(),pWeightsFlag[pW].Data())));
      if(dfrParticleWeights[t][pW])
      {
       for(Int_t eW=0;eW<2;eW++)
       {
        dfrEventWeights[t][pW][eW] = dynamic_cast<TList*>(dfrParticleWeights[t][pW]->FindObject(Form("%s, pWeights %s, eWeights %s",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
        if(dfrEventWeights[t][pW][eW])
        {
         diffFlowFinalCorrelations[t][pW][eW] = dynamic_cast<TList*>(dfrEventWeights[t][pW][eW]->FindObject(Form("Correlations (%s, pWeights %s, eWeights %s)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
         diffFlowFinalCorrections[t][pW][eW] = dynamic_cast<TList*>(dfrEventWeights[t][pW][eW]->FindObject(Form("Corrections (%s, pWeights %s, eWeights %s)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));
         diffFlowFinalCovariances[t][pW][eW] = dynamic_cast<TList*>(dfrEventWeights[t][pW][eW]->FindObject(Form("Covariances (%s, pWeights %s, eWeights %s)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())));      
         for(Int_t nua=0;nua<2;nua++)
         {
          dfrCorrections[t][pW][eW][nua] = dynamic_cast<TList*>(dfrEventWeights[t][pW][eW]->FindObject(Form("%s, pWeights %s, eWeights %s, %s for NUA",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())));
          if(dfrCorrections[t][pW][eW][nua])
          {
           diffFlowFinalCumulants[t][pW][eW][nua] = dynamic_cast<TList*>(dfrCorrections[t][pW][eW][nua]->FindObject(Form("Cumulants (%s, pWeights %s, eWeights %s, %s for NUA)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())));
           diffFlowFinalFlow[t][pW][eW][nua] = dynamic_cast<TList*>(dfrCorrections[t][pW][eW][nua]->FindObject(Form("Differential Flow (%s, pWeights %s, eWeights %s, %s for NUA)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())));
          } else // to if(dfrCorrections[t][pW][eW][nua])
            {
             cout<<"WARNING: dfrCorrections[t][pW][eW][nua] is NULL in AFAWQC::GOH() !!!!"<<endl;
             cout<<"t = "<<t<<endl;
             cout<<"pW = "<<pW<<endl;
             cout<<"eW = "<<eW<<endl;
             cout<<"nua = "<<nua<<endl;
            }
         } // end of for(Int_t nua=0;nua<2;nua++)
        } else // to if(dfrEventWeights[t][pW][eW])
          {
           cout<<"WARNING: dfrEventWeights[t][pW][eW] is NULL in AFAWQC::GOH() !!!!"<<endl;
           cout<<"t = "<<t<<endl;
           cout<<"pW = "<<pW<<endl;
           cout<<"eW = "<<eW<<endl;
          }
       } // end of for(Int_t eW=0;eW<2;eW++)
      } else // to if(dfrParticleWeights[t][pW])
        {
         cout<<"WARNING: dfrParticleWeights[t][pW] is NULL in AFAWQC::GOH() !!!!"<<endl;
         cout<<"t = "<<t<<endl;
         cout<<"pW = "<<pW<<endl;
        }
     } // end of for(Int_t pW=0;pW<(1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights));pW++)
    } else // to if(dfrType[t])
      {
       cout<<"WARNING: dfrType[t] is NULL in AFAWQC::GOH() !!!!"<<endl;
       cout<<"t = "<<t<<endl;
      }
   } // end of for(Int_t t=0;t<2;t++)
  } else // to if(diffFlowListResults)
    {
     cout<<"WARNING: diffFlowListResults is NULL in AFAWQC::GOH() !!!!"<<endl;
    }
 
 TH2D *finalCorrelations2D[2][2][2][4] = {{{{NULL}}}};
 TH1D *finalCorrelations1D[2][2][2][2][4] = {{{{{NULL}}}}}; 
 TH2D *finalCumulants2D[2][2][2][2][4] = {{{{{NULL}}}}};
 TH1D *finalCumulantsPt[2][2][2][2][4] = {{{{{NULL}}}}}; 
 TH1D *finalCumulantsEta[2][2][2][2][4] = {{{{{NULL}}}}}; 
 TH2D *finalCorrections2D[2][2][2][4] = {{{{NULL}}}};
 TH1D *finalCorrections1D[2][2][2][2][4] = {{{{{NULL}}}}}; 
 TH2D *finalCovariances2D[2][2][2][4] = {{{{NULL}}}};
 TH1D *finalCovariances1D[2][2][2][2][4] = {{{{{NULL}}}}}; 
 TH2D *finalFlow2D[2][2][2][2][4] = {{{{{NULL}}}}};
 TH1D *finalFlowPt[2][2][2][2][4] = {{{{{NULL}}}}}; 
 TH1D *finalFlowEta[2][2][2][2][4] = {{{{{NULL}}}}}; 
 TH2D *nonEmptyBins2D[2] = {NULL};
 TH1D *nonEmptyBins1D[2][2] = {{NULL}};
 
 for(Int_t t=0;t<2;t++)
 {
  // 2D:
  nonEmptyBins2D[t] = dynamic_cast<TH2D*>(dfrType[t]->FindObject(Form("%s, (p_{t},#eta)",typeFlag[t].Data())));
  if(nonEmptyBins2D[t])
  {
   this->SetNonEmptyBins2D(nonEmptyBins2D[t],t);
  } else
    { 
     cout<<"WARNING: nonEmptyBins2D[t] is NULL in AFAWQC::GOH() !!!!"<<endl;
     cout<<"t = "<<t<<endl;
    }
  // 1D:
  for(Int_t pe=0;pe<2;pe++)
  {
   nonEmptyBins1D[t][pe] = dynamic_cast<TH1D*>(dfrType[t]->FindObject(Form("%s, %s",typeFlag[t].Data(),ptEtaFlag[pe].Data())));
   if(nonEmptyBins1D[t][pe])
   {
    this->SetNonEmptyBins1D(nonEmptyBins1D[t][pe],t,pe);
   } else
     { 
      cout<<"WARNING: nonEmptyBins1D[t][pe] is NULL in AFAWQC::GOH() !!!!"<<endl;
      cout<<"t = "<<t<<endl;
      cout<<"pe = "<<pe<<endl;
     }
  }  
  
  for(Int_t pW=0;pW<(1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights));pW++)
  {
   for(Int_t eW=0;eW<2;eW++)
   {
    // correlations:
    if(diffFlowFinalCorrelations[t][pW][eW])
    {
     for(Int_t correlationIndex=0;correlationIndex<4;correlationIndex++)
     {
      // 2D:
      finalCorrelations2D[t][pW][eW][correlationIndex] = dynamic_cast<TH2D*>(diffFlowFinalCorrelations[t][pW][eW]->FindObject(Form("%s, (p_{t},#eta)",correlationName[correlationIndex].Data())));
      if(finalCorrelations2D[t][pW][eW][correlationIndex])
      {
       this->SetFinalCorrelations2D(finalCorrelations2D[t][pW][eW][correlationIndex],t,pW,eW,correlationIndex);
      } else
        {
         cout<<"WARNING: finalCorrelations2D[t][pW][eW][correlationIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
         cout<<"t = "<<t<<endl;
         cout<<"pW = "<<pW<<endl;
         cout<<"eW = "<<eW<<endl;
         cout<<"ci = "<<correlationIndex<<endl;
        }
       // 1D
       for(Int_t pe=0;pe<2;pe++)
       {
        if(pe==0)
        {
         finalCorrelations1D[t][pW][eW][pe][correlationIndex] = dynamic_cast<TH1D*>(diffFlowFinalCorrelations[t][pW][eW]->FindObject(Form("%s, (p_{t})",correlationName[correlationIndex].Data())));
        }
        if(pe==1)
        {
         finalCorrelations1D[t][pW][eW][pe][correlationIndex] = dynamic_cast<TH1D*>(diffFlowFinalCorrelations[t][pW][eW]->FindObject(Form("%s, (#eta)",correlationName[correlationIndex].Data())));
        }          
        if(finalCorrelations1D[t][pW][eW][pe][correlationIndex])
        {
         this->SetFinalCorrelations1D(finalCorrelations1D[t][pW][eW][pe][correlationIndex],t,pW,eW,pe,correlationIndex);
        } else
          {
           cout<<"WARNING: finalCorrelations1D[t][pW][eW][pe][correlationIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
           cout<<"t = "<<t<<endl;
           cout<<"pW = "<<pW<<endl;
           cout<<"eW = "<<eW<<endl;
           cout<<"pe = "<<pe<<endl;
           cout<<"ci = "<<correlationIndex<<endl;
          }   
       } // end of for(Int_t pe=0;pe<2;pe++)  
     } // end of for(Int_t correlationIndex=0;correlationIndex<4;correlationIndex++)
    } else // to if(diffFlowFinalCorrelations[t][pW][eW])
      { 
       cout<<"WARNING: diffFlowFinalCorrelations[t][pW] is NULL in AFAWQC::GOH() !!!!"<<endl;
       cout<<"t = "<<t<<endl;
       cout<<"pW = "<<pW<<endl;
       cout<<"eW = "<<eW<<endl;
      }
    // corrections:
    if(diffFlowFinalCorrections[t][pW][eW])
    {
     for(Int_t correctionIndex=0;correctionIndex<4;correctionIndex++)
     {
      // 2D:
      finalCorrections2D[t][pW][eW][correctionIndex] = dynamic_cast<TH2D*>(diffFlowFinalCorrections[t][pW][eW]->FindObject(Form("%s, (p_{t},#eta)",correctionName[correctionIndex].Data())));
      if(finalCorrections2D[t][pW][eW][correctionIndex])
      {
       this->SetFinalCorrections2D(finalCorrections2D[t][pW][eW][correctionIndex],t,pW,eW,correctionIndex);
      } else
        {
         cout<<"WARNING: finalCorrections2D[t][pW][eW][correctionIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
         cout<<"t = "<<t<<endl;
         cout<<"pW = "<<pW<<endl;
         cout<<"eW = "<<eW<<endl;
         cout<<"ci = "<<correctionIndex<<endl;
        }
      // 1D
      for(Int_t pe=0;pe<2;pe++)
      {
       if(pe==0)
       {
        finalCorrections1D[t][pW][eW][pe][correctionIndex] = dynamic_cast<TH1D*>(diffFlowFinalCorrections[t][pW][eW]->FindObject(Form("%s, (p_{t})",correctionName[correctionIndex].Data())));
       }
       if(pe==1)
       {
        finalCorrections1D[t][pW][eW][pe][correctionIndex] = dynamic_cast<TH1D*>(diffFlowFinalCorrections[t][pW][eW]->FindObject(Form("%s, (#eta)",correctionName[correctionIndex].Data())));
       }          
       if(finalCorrections1D[t][pW][eW][pe][correctionIndex])
       {
        this->SetFinalCorrections1D(finalCorrections1D[t][pW][eW][pe][correctionIndex],t,pW,eW,pe,correctionIndex);
       } else
         {
          cout<<"WARNING: finalCorrections1D[t][pW][eW][pe][correctionIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
          cout<<"t = "<<t<<endl;
          cout<<"pW = "<<pW<<endl;
          cout<<"eW = "<<eW<<endl;
          cout<<"pe = "<<pe<<endl;
          cout<<"ci = "<<correctionIndex<<endl;
         }   
      } // end of for(Int_t pe=0;pe<2;pe++)  
     } // end of for(Int_t correctionIndex=0;correctionIndex<4;correctionIndex++)
    } else // to if(diffFlowFinalCorrections[t][pW][eW])
      { 
       cout<<"WARNING: diffFlowFinalCorrections[t][pW][eW] is NULL in AFAWQC::GOH() !!!!"<<endl;
       cout<<"t = "<<t<<endl;
       cout<<"pW = "<<pW<<endl;
       cout<<"eW = "<<eW<<endl;
      }     
    // covariances:
    if(diffFlowFinalCovariances[t][pW][eW])
    {
     for(Int_t covarianceIndex=0;covarianceIndex<4;covarianceIndex++)
     {
      // 2D:
      finalCovariances2D[t][pW][eW][covarianceIndex] = dynamic_cast<TH2D*>(diffFlowFinalCovariances[t][pW][eW]->FindObject(Form("%s, (p_{t},#eta)",covarianceName[covarianceIndex].Data())));
      if(finalCovariances2D[t][pW][eW][covarianceIndex])
      {
       this->SetFinalCovariances2D(finalCovariances2D[t][pW][eW][covarianceIndex],t,pW,eW,covarianceIndex);
      } else
        {
         cout<<"WARNING: finalCovariances2D[t][pW][eW][nua][covarianceIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
         cout<<"t = "<<t<<endl;
         cout<<"pW = "<<pW<<endl;
         cout<<"eW = "<<eW<<endl;
         cout<<"ci = "<<covarianceIndex<<endl;
        }
      // 1D:
      for(Int_t pe=0;pe<2;pe++)
      {
       if(pe==0)
       {
        finalCovariances1D[t][pW][eW][pe][covarianceIndex] = dynamic_cast<TH1D*>(diffFlowFinalCovariances[t][pW][eW]->FindObject(Form("%s, (p_{t})",covarianceName[covarianceIndex].Data())));
       }
       if(pe==1)
       {
        finalCovariances1D[t][pW][eW][pe][covarianceIndex] = dynamic_cast<TH1D*>(diffFlowFinalCovariances[t][pW][eW]->FindObject(Form("%s, (#eta)",covarianceName[covarianceIndex].Data())));
       }          
       if(finalCovariances1D[t][pW][eW][pe][covarianceIndex])
       {
        this->SetFinalCovariances1D(finalCovariances1D[t][pW][eW][pe][covarianceIndex],t,pW,eW,pe,covarianceIndex);
       } else
         {
          cout<<"WARNING: finalCovariances1D[t][pW][eW][pe][covarianceIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
          cout<<"t = "<<t<<endl;
          cout<<"pW = "<<pW<<endl;
          cout<<"eW = "<<eW<<endl;
          cout<<"pe = "<<pe<<endl;
          cout<<"ci = "<<covarianceIndex<<endl;
         }   
      } // end of for(Int_t pe=0;pe<2;pe++)  
     } // end of for(Int_t covarianceIndex=0;covarianceIndex<4;covarianceIndex++)
    } else // to if(diffFlowFinalCovariances[t][pW][eW])
      { 
       cout<<"WARNING: diffFlowFinalCovariances[t][pW][eW] is NULL in AFAWQC::GOH() !!!!"<<endl;
       cout<<"t = "<<t<<endl;
       cout<<"pW = "<<pW<<endl;
       cout<<"eW = "<<eW<<endl;
      }        
       
    for(Int_t nua=0;nua<2;nua++)
    {  
     // cumulants:
     if(diffFlowFinalCumulants[t][pW][eW][nua])
     {
      for(Int_t cumulantIndex=0;cumulantIndex<4;cumulantIndex++)
      {
       // 2D:
       finalCumulants2D[t][pW][eW][nua][cumulantIndex] = dynamic_cast<TH2D*>(diffFlowFinalCumulants[t][pW][eW][nua]->FindObject(Form("%s, (p_{t},#eta)",cumulantName[cumulantIndex].Data())));
       if(finalCumulants2D[t][pW][eW][nua][cumulantIndex])
       {
        this->SetFinalCumulants2D(finalCumulants2D[t][pW][eW][nua][cumulantIndex],t,pW,eW,nua,cumulantIndex);
       } else
         {
          cout<<"WARNING: finalCumulants2D[t][pW][eW][nua][cumulantIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
          cout<<"t = "<<t<<endl;
          cout<<"pW = "<<pW<<endl;
          cout<<"eW = "<<eW<<endl;
          cout<<"nua = "<<nua<<endl;
          cout<<"ci = "<<cumulantIndex<<endl;
         }
       // pt:  
       finalCumulantsPt[t][pW][eW][nua][cumulantIndex] = dynamic_cast<TH1D*>(diffFlowFinalCumulants[t][pW][eW][nua]->FindObject(Form("%s, (p_{t})",cumulantName[cumulantIndex].Data())));
       if(finalCumulantsPt[t][pW][eW][nua][cumulantIndex])
       {
        this->SetFinalCumulantsPt(finalCumulantsPt[t][pW][eW][nua][cumulantIndex],t,pW,eW,nua,cumulantIndex);
       } else
         {
          cout<<"WARNING: finalCumulantsPt[t][pW][eW][nua][cumulantIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
          cout<<"t = "<<t<<endl;
          cout<<"pW = "<<pW<<endl;
          cout<<"eW = "<<eW<<endl;
          cout<<"nua = "<<nua<<endl;
          cout<<"ci = "<<cumulantIndex<<endl;
         } 
       // eta:
       finalCumulantsEta[t][pW][eW][nua][cumulantIndex] = dynamic_cast<TH1D*>(diffFlowFinalCumulants[t][pW][eW][nua]->FindObject(Form("%s, (#eta)",cumulantName[cumulantIndex].Data())));
       if(finalCumulantsEta[t][pW][eW][nua][cumulantIndex])
       {
        this->SetFinalCumulantsEta(finalCumulantsEta[t][pW][eW][nua][cumulantIndex],t,pW,eW,nua,cumulantIndex);
       } else
         {
          cout<<"WARNING: finalCumulantsEta[t][pW][eW][nua][cumulantIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
          cout<<"t = "<<t<<endl;
          cout<<"pW = "<<pW<<endl;
          cout<<"eW = "<<eW<<endl;
          cout<<"nua = "<<nua<<endl;
          cout<<"ci = "<<cumulantIndex<<endl;
         }   
      } // end of for(Int_t cumulantIndex=0;cumulantIndex<4;cumulantIndex++)
     } else // to if(diffFlowFinalCumulants[t][pW][eW][nua])
       { 
        cout<<"WARNING: diffFlowFinalCumulants[t][pW][eW][nua] is NULL in AFAWQC::GOH() !!!!"<<endl;
        cout<<"t = "<<t<<endl;
        cout<<"pW = "<<pW<<endl;
        cout<<"eW = "<<eW<<endl;
        cout<<"nua = "<<nua<<endl;
       }      
     // flow:
     if(diffFlowFinalFlow[t][pW][eW][nua])
     {
      for(Int_t flowIndex=0;flowIndex<4;flowIndex++)
      {
       // 2D:
       finalFlow2D[t][pW][eW][nua][flowIndex] = dynamic_cast<TH2D*>(diffFlowFinalFlow[t][pW][eW][nua]->FindObject(Form("%s, (p_{t},#eta)",flowName[flowIndex].Data())));
       if(finalFlow2D[t][pW][eW][nua][flowIndex])
       {
        this->SetFinalFlow2D(finalFlow2D[t][pW][eW][nua][flowIndex],t,pW,eW,nua,flowIndex);
       } else
         {
          cout<<"WARNING: finalFlow2D[t][pW][eW][nua][flowIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
          cout<<"t = "<<t<<endl;
          cout<<"pW = "<<pW<<endl;
          cout<<"eW = "<<eW<<endl;
          cout<<"nua = "<<nua<<endl;
          cout<<"ci = "<<flowIndex<<endl;
         }
       // pt:
       finalFlowPt[t][pW][eW][nua][flowIndex] = dynamic_cast<TH1D*>(diffFlowFinalFlow[t][pW][eW][nua]->FindObject(Form("%s, (p_{t})",flowName[flowIndex].Data())));
       if(finalFlowPt[t][pW][eW][nua][flowIndex])
       {
        this->SetFinalFlowPt(finalFlowPt[t][pW][eW][nua][flowIndex],t,pW,eW,nua,flowIndex);
       } else
         {
          cout<<"WARNING: finalFlow1D[t][pW][nua][flowIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
          cout<<"t = "<<t<<endl;
          cout<<"pW = "<<pW<<endl;
          cout<<"eW = "<<eW<<endl;
          cout<<"nua = "<<nua<<endl;
          cout<<"ci = "<<flowIndex<<endl;
         }   
       // eta: 
       finalFlowEta[t][pW][eW][nua][flowIndex] = dynamic_cast<TH1D*>(diffFlowFinalFlow[t][pW][eW][nua]->FindObject(Form("%s, (#eta)",flowName[flowIndex].Data())));          
       if(finalFlowEta[t][pW][eW][nua][flowIndex])
       {
        this->SetFinalFlowEta(finalFlowEta[t][pW][eW][nua][flowIndex],t,pW,eW,nua,flowIndex);
       } else
         {
          cout<<"WARNING: finalFlow1D[t][pW][nua][flowIndex] is NULL in AFAWQC::GOH() !!!!"<<endl;
          cout<<"t = "<<t<<endl;
          cout<<"pW = "<<pW<<endl;
          cout<<"eW = "<<eW<<endl;
          cout<<"nua = "<<nua<<endl;
          cout<<"ci = "<<flowIndex<<endl;
         }   
      } // end of for(Int_t flowIndex=0;flowIndex<4;flowIndex++)
     } else // to if(diffFlowFinalFlow[t][pW][eW][nua])
       { 
        cout<<"WARNING: diffFlowFinalFlow[t][pW][eW][nua] is NULL in AFAWQC::GOH() !!!!"<<endl;
        cout<<"t = "<<t<<endl;
        cout<<"pW = "<<pW<<endl;
        cout<<"eW = "<<eW<<endl;
        cout<<"nua = "<<nua<<endl;
       }               
    } // end of for(Int_t nua=0;nua<2;nua++)
   } // end of for(Int_t eW=0;eW<2;eW++) 
  } // end of for(Int_t pW=0;pW<(1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights));pW++)
 } // end of for(Int_t t=0;t<2;t++)
 
  // x.) nested loops:
  TList *nestedLoopsList = dynamic_cast<TList*>(outputListHistos->FindObject("Nested Loops"));
  if(nestedLoopsList) this->SetNestedLoopsList(nestedLoopsList);
  TProfile *evaluateNestedLoops = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fEvaluateNestedLoops"));
  Bool_t bEvaluateNestedLoopsForIntFlow = kFALSE;
  Bool_t bEvaluateNestedLoopsForDiffFlow = kFALSE;
  if(evaluateNestedLoops)
  {
   this->SetEvaluateNestedLoops(evaluateNestedLoops);
   bEvaluateNestedLoopsForIntFlow = (Int_t)evaluateNestedLoops->GetBinContent(1);
   bEvaluateNestedLoopsForDiffFlow = (Int_t)evaluateNestedLoops->GetBinContent(2);
  }
  
  if(bEvaluateNestedLoopsForIntFlow)
  {
   TProfile *directCorrelations = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrelation"));
   if(directCorrelations) this->SetDirectCorrelations(directCorrelations);
   TProfile *directCorrectionsCos = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrectionsCos"));
   if(directCorrectionsCos) this->SetDirectCorrectionsCos(directCorrectionsCos);
   TProfile *directCorrectionsSin = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrectionsSin"));
   if(directCorrectionsSin) this->SetDirectCorrectionsSin(directCorrectionsSin);
   if(bUsePhiWeights||bUsePtWeights||bUseEtaWeights)
   {
    TProfile *directCorrelationsW = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrelationW"));
    if(directCorrelationsW) this->SetDirectCorrelationsW(directCorrelationsW);
    TProfile *directCorrectionsCosW = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrectionsCosW"));
    if(directCorrectionsCosW) this->SetDirectCorrectionsCosW(directCorrectionsCosW);
    TProfile *directCorrectionsSinW = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrectionsSinW")); 
    if(directCorrectionsSinW) this->SetDirectCorrectionsSinW(directCorrectionsSinW);
   }
  }
  
  if(bEvaluateNestedLoopsForDiffFlow)
  {
   TProfile *directCorrelationsDiffFlow = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrelationsDiffFlow"));
   if(directCorrelationsDiffFlow) this->SetDirectCorrelationsDiffFlow(directCorrelationsDiffFlow);
   TProfile *directCorrectionsDiffFlowCos = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrectionsDiffFlowCos"));
   if(directCorrectionsDiffFlowCos) this->SetDirectCorrectionsDiffFlowCos(directCorrectionsDiffFlowCos);
   TProfile *directCorrectionsDiffFlowSin = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrectionsDiffFlowSin"));
   if(directCorrectionsDiffFlowSin) this->SetDirectCorrectionsDiffFlowSin(directCorrectionsDiffFlowSin);
   if(bUsePhiWeights||bUsePtWeights||bUseEtaWeights)
   {
    TProfile *directCorrelationsDiffFlowW = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrelationsDiffFlowW")); 
    if(directCorrelationsDiffFlowW) this->SetDirectCorrelationsDiffFlowW(directCorrelationsDiffFlowW);
    TProfile *directCorrectionsDiffFlowCosW = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrectionsDiffFlowCosW"));
    if(directCorrectionsDiffFlowCosW) this->SetDirectCorrectionsDiffFlowCosW(directCorrectionsDiffFlowCosW);
    TProfile *directCorrectionsDiffFlowSinW = dynamic_cast<TProfile*>(nestedLoopsList->FindObject("fDirectCorrectionsDiffFlowSinW"));
    if(directCorrectionsDiffFlowSinW) this->SetDirectCorrectionsDiffFlowSinW(directCorrectionsDiffFlowSinW);
   }
  }
  
 } 
}


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


void AliFlowAnalysisWithQCumulants::CalculateFinalCorrectionsForNonUniformAcceptanceForCumulantsForIntFlow(Bool_t useParticleWeights, TString eventWeights)
{
 // calculate final corrections for non-uniform acceptance for QC{2}, QC{4}, QC{6} and QC{8}
 
 // corrections for non-uniform acceptance (NUA) are stored in histogram fCorrectionsForNUA,
 // binning of fCorrectionsForNUA is organized as follows:
 //
 // 1st bin: correction to QC{2}
 // 2nd bin: correction to QC{4}
 // 3rd bin: correction to QC{6}
 // 4th bin: correction to QC{8}
 
 // shortcuts flags:
 Int_t pW = (Int_t)(useParticleWeights);
 
 Int_t eW = -1;
 
 if(eventWeights == "exact")
 {
  eW = 0;
 }

 for(Int_t sc=0;sc<2;sc++) // sin or cos terms flag
 {
  if(!(fQCorrelations[pW][eW] && fQCorrections[pW][eW][sc] && fCorrections[pW][eW]))
  {
   cout<<"WARNING: fQCorrelations[pW][eW] && fQCorrections[pW][eW][sc] && fCorrections[pW][eW] is NULL in AFAWQC::CFCFNUAFIF() !!!!"<<endl;
   cout<<"pW = "<<pW<<endl;
   cout<<"eW = "<<eW<<endl;
   cout<<"sc = "<<sc<<endl;
   exit(0);
  }
 }  

 // measured 2-, 4-, 6- and 8-particle azimuthal correlations (biased with non-uniform acceptance!):
 Double_t two = fQCorrelations[pW][eW]->GetBinContent(1); // <<2>>
 //Double_t four = fQCorrelations[pW][eW]->GetBinContent(11); // <<4>>
 //Double_t six = fQCorrelations[pW][eW]->GetBinContent(24); // <<6>>
 //Double_t eight = fQCorrelations[pW][eW]->GetBinContent(31); // <<8>>
 
 // correction terms to QC{2}:
 // <<cos(n*phi1)>>^2
 Double_t two1stTerm = pow(fQCorrections[pW][eW][1]->GetBinContent(1),2); 
 // <<sin(n*phi1)>>^2
 Double_t two2ndTerm = pow(fQCorrections[pW][eW][0]->GetBinContent(1),2); 
 // final corrections for non-uniform acceptance to QC{2}:
 Double_t correctionQC2 = -1.*two1stTerm-1.*two2ndTerm;
 fCorrections[pW][eW]->SetBinContent(1,correctionQC2); 
 
 // correction terms to QC{4}:
 // <<cos(n*phi1)>> <<cos(n*(phi1-phi2-phi3))>>
 Double_t four1stTerm = fQCorrections[pW][eW][1]->GetBinContent(1)*fQCorrections[pW][eW][1]->GetBinContent(3);  
 // <<sin(n*phi1)>> <<sin(n*(phi1-phi2-phi3))>>
 Double_t four2ndTerm = fQCorrections[pW][eW][0]->GetBinContent(1)*fQCorrections[pW][eW][0]->GetBinContent(3);  
 // <<cos(n*(phi1+phi2))>>^2
 Double_t four3rdTerm = pow(fQCorrections[pW][eW][1]->GetBinContent(2),2); 
 // <<sin(n*(phi1+phi2))>>^2
 Double_t four4thTerm = pow(fQCorrections[pW][eW][0]->GetBinContent(2),2); 
 // <<cos(n*(phi1+phi2))>> (<<cos(n*phi1)>>^2 - <<sin(n*phi1)>>^2)
 Double_t four5thTerm = fQCorrections[pW][eW][1]->GetBinContent(2)
                      * (pow(fQCorrections[pW][eW][1]->GetBinContent(1),2)-pow(fQCorrections[pW][eW][0]->GetBinContent(1),2));
 // <<sin(n*(phi1+phi2))>> <<cos(n*phi1)>> <<sin(n*phi1)>>
 Double_t four6thTerm = fQCorrections[pW][eW][0]->GetBinContent(2)
                      * fQCorrections[pW][eW][1]->GetBinContent(1)
                      * fQCorrections[pW][eW][0]->GetBinContent(1);         
 // <<cos(n*(phi1-phi2))>> (<<cos(n*phi1)>>^2 + <<sin(n*phi1)>>^2)
 Double_t four7thTerm = two*(pow(fQCorrections[pW][eW][1]->GetBinContent(1),2)+pow(fQCorrections[pW][eW][0]->GetBinContent(1),2));  
 // (<<cos(n*phi1)>>^2 + <<sin(n*phi1)>>^2)^2
 Double_t four8thTerm = pow(pow(fQCorrections[pW][eW][1]->GetBinContent(1),2)+pow(fQCorrections[pW][eW][0]->GetBinContent(1),2),2);      
 // final correction to QC{4}:
 Double_t correctionQC4 = -4.*four1stTerm+4.*four2ndTerm-four3rdTerm-four4thTerm
                        + 4.*four5thTerm+8.*four6thTerm+8.*four7thTerm-6.*four8thTerm;                            
 fCorrections[pW][eW]->SetBinContent(2,correctionQC4);   

 // ... to be improved (continued for 6th and 8th order)                                                    

} // end of AliFlowAnalysisWithQCumulants::CalculateFinalCorrectionsForNonUniformAcceptanceForCumulantsForIntFlow(Bool_t useParticleWeights, TString eventWeights)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateFinalCorrectionsForNonUniformAcceptanceForDifferentialFlow(Bool_t useParticleWeights,TString type)
{
 
 useParticleWeights=kFALSE;
 type="ac";
  
 // calculate final corrections due to non-uniform acceptance of the detector to reduced multi-particle correlations
 /*
 if(!(useParticleWeights))
 {
  if(type == "POI")
  { 
   // **** corrections for non-uniform acceptance for 2nd order QC' for POI's ****
   
   // 1st term: <<cos(n*psi)>><<cos(n*phi)>>:
   if(fCorrectionsCosP1nPsiPtEtaPOI && fQCorrectionsCos)
   {
    // pt,eta: 
    if(f2pFinalCorrectionsForNUAPtEtaPOI) f2pFinalCorrectionsForNUAPtEtaPOI->Reset(); // to be improved
    TH2D *correctionPtEta1stTerm = new TH2D(*(fCorrectionsCosP1nPsiPtEtaPOI->ProjectionXY("","e")));
    correctionPtEta1stTerm->Scale(fQCorrectionsCos->GetBinContent(1)); // to be improved: are errors propagated correctly here?   
    if(f2pFinalCorrectionsForNUAPtEtaPOI) f2pFinalCorrectionsForNUAPtEtaPOI->Add(correctionPtEta1stTerm); // to be improved (if condition goes somewhere else)
    delete correctionPtEta1stTerm;
    // pt:
    if(f2pFinalCorrectionsForNUAPtPOI) f2pFinalCorrectionsForNUAPtPOI->Reset(); // to be improved
    TH1D *correctionPt1stTerm = new TH1D(*((this->MakePtProjection(fCorrectionsCosP1nPsiPtEtaPOI))->ProjectionX("","e"))); // to be improved: are errors propagated correctly here? 
    correctionPt1stTerm->Scale(fQCorrectionsCos->GetBinContent(1)); // to be improved: are errors propagated correctly here? 
    if(f2pFinalCorrectionsForNUAPtPOI) f2pFinalCorrectionsForNUAPtPOI->Add(correctionPt1stTerm); // to be improved (if condition goes somewhere else)
    delete correctionPt1stTerm;
    // eta:
    if(f2pFinalCorrectionsForNUAEtaPOI) f2pFinalCorrectionsForNUAEtaPOI->Reset(); // to be improved    
    TH1D *correctionEta1stTerm = new TH1D(*((this->MakeEtaProjection(fCorrectionsCosP1nPsiPtEtaPOI))->ProjectionX("","e"))); // to be improved: are errors propagated correctly here? 
    correctionEta1stTerm->Scale(fQCorrectionsCos->GetBinContent(1)); // to be improved: are errors propagated correctly here? 
    if(f2pFinalCorrectionsForNUAEtaPOI) f2pFinalCorrectionsForNUAEtaPOI->Add(correctionEta1stTerm); // to be improved (if condition goes somewhere else)
    delete correctionEta1stTerm;    
   } else
     { 
      cout<<"WARNING: (fCorrectionsCosP1nPsiPtEtaPOI && fQCorrectionsCos && f2pFinalCorrectionsForNUAPtEtaPOI) is NULL in QC::CFCFNUAFDF() !!!!  "<<endl;
      cout<<"         Corrections for non-uniform acceptance for differential flow are not correct."<<endl;
     } 
     
   // 2nd term: <<sin(n*psi)>><<sin(n*phi)>>:  
   if(fCorrectionsSinP1nPsiPtEtaPOI && fQCorrectionsSin)
   {
    // pt,eta:
    TH2D *correctionPtEta2ndTerm = new TH2D(*(fCorrectionsSinP1nPsiPtEtaPOI->ProjectionXY("","e")));
    correctionPtEta2ndTerm->Scale(fQCorrectionsSin->GetBinContent(1)); // to be improved: are errors propagated correctly here?    
    if(f2pFinalCorrectionsForNUAPtEtaPOI) f2pFinalCorrectionsForNUAPtEtaPOI->Add(correctionPtEta2ndTerm); // to be improved (if condition goes somewhere else)
    delete correctionPtEta2ndTerm;
    // pt:
    TH1D *correctionPt2ndTerm = new TH1D(*((this->MakePtProjection(fCorrectionsSinP1nPsiPtEtaPOI))->ProjectionX("","e"))); // to be improved: are errors propagated correctly here? 
    correctionPt2ndTerm->Scale(fQCorrectionsSin->GetBinContent(1)); // to be improved: are errors propagated correctly here? 
    if(f2pFinalCorrectionsForNUAPtPOI) f2pFinalCorrectionsForNUAPtPOI->Add(correctionPt2ndTerm); // to be improved (if condition goes somewhere else)
    delete correctionPt2ndTerm;
    // eta:
    TH1D *correctionEta2ndTerm = new TH1D(*((this->MakeEtaProjection(fCorrectionsSinP1nPsiPtEtaPOI))->ProjectionX("","e"))); // to be improved: are errors propagated correctly here? 
    correctionEta2ndTerm->Scale(fQCorrectionsSin->GetBinContent(1)); // to be improved: are errors propagated correctly here? 
    if(f2pFinalCorrectionsForNUAEtaPOI) f2pFinalCorrectionsForNUAEtaPOI->Add(correctionEta2ndTerm); // to be improved (if condition goes somewhere else)
    delete correctionEta2ndTerm; 
   } else
     { 
      cout<<"WARNING: (fCorrectionsSinP1nPsiPtEtaPOI && fQCorrectionsSin) is NULL in QC::CFCFNUAFDF() !!!!  "<<endl;
      cout<<"         Corrections for non-uniform acceptance for differential flow are not correct."<<endl;
     } 
  } else if(type == "RP")
    {
     // **** corrections for non-uniform acceptance for 2nd order QC' for RP's ****
   
     // 1st term: <<cos(n*psi)>><<cos(n*phi)>>:
     if(fCorrectionsCosP1nPsiPtEtaRP && fQCorrectionsCos)
     {
      // pt,eta: 
      if(f2pFinalCorrectionsForNUAPtEtaRP) f2pFinalCorrectionsForNUAPtEtaRP->Reset(); // to be improved
      TH2D *correctionPtEta1stTerm = new TH2D(*(fCorrectionsCosP1nPsiPtEtaRP->ProjectionXY("","e")));
      correctionPtEta1stTerm->Scale(fQCorrectionsCos->GetBinContent(1)); // to be improved: are errors propagated correctly here?    
      if(f2pFinalCorrectionsForNUAPtEtaRP) f2pFinalCorrectionsForNUAPtEtaRP->Add(correctionPtEta1stTerm); // to be improved (if condition goes somewhere else)
      delete correctionPtEta1stTerm;
      // pt:
      if(f2pFinalCorrectionsForNUAPtRP) f2pFinalCorrectionsForNUAPtRP->Reset(); // to be improved
      TH1D *correctionPt1stTerm = new TH1D(*((this->MakePtProjection(fCorrectionsCosP1nPsiPtEtaRP))->ProjectionX("","e"))); // to be improved: are errors propagated correctly here? 
      correctionPt1stTerm->Scale(fQCorrectionsCos->GetBinContent(1)); // to be improved: are errors propagated correctly here? 
      if(f2pFinalCorrectionsForNUAPtRP) f2pFinalCorrectionsForNUAPtRP->Add(correctionPt1stTerm); // to be improved (if condition goes somewhere else)
      delete correctionPt1stTerm;
      // eta:
      if(f2pFinalCorrectionsForNUAEtaRP) f2pFinalCorrectionsForNUAEtaRP->Reset(); // to be improved
      TH1D *correctionEta1stTerm = new TH1D(*((this->MakeEtaProjection(fCorrectionsCosP1nPsiPtEtaRP))->ProjectionX("","e"))); // to be improved: are errors propagated correctly here? 
      correctionEta1stTerm->Scale(fQCorrectionsCos->GetBinContent(1)); // to be improved: are errors propagated correctly here? 
      if(f2pFinalCorrectionsForNUAEtaRP) f2pFinalCorrectionsForNUAEtaRP->Add(correctionEta1stTerm); // to be improved (if condition goes somewhere else)
      delete correctionEta1stTerm;    
     } else
       { 
        cout<<"WARNING: (fCorrectionsCosP1nPsiPtEtaRP && fQCorrectionsCos) is NULL in QC::CFCFNUAFDF() !!!!  "<<endl;
        cout<<"         Corrections for non-uniform acceptance for differential flow are not correct."<<endl;
       } 
     // 2nd term: <<sin(n*psi)>><<sin(n*phi)>>:  
     if(fCorrectionsSinP1nPsiPtEtaRP && fQCorrectionsSin)
     {
      // pt,eta: 
      TH2D *correctionPtEta2ndTerm = new TH2D(*(fCorrectionsSinP1nPsiPtEtaRP->ProjectionXY("","e")));
      correctionPtEta2ndTerm->Scale(fQCorrectionsSin->GetBinContent(1)); // to be improved: are errors propagated correctly here?    
      if(f2pFinalCorrectionsForNUAPtEtaRP) f2pFinalCorrectionsForNUAPtEtaRP->Add(correctionPtEta2ndTerm); // to be improved (if condition goes somewhere else)
      delete correctionPtEta2ndTerm;
      // pt:
      TH1D *correctionPt2ndTerm = new TH1D(*((this->MakePtProjection(fCorrectionsSinP1nPsiPtEtaRP))->ProjectionX("","e"))); // to be improved: are errors propagated correctly here? 
      correctionPt2ndTerm->Scale(fQCorrectionsSin->GetBinContent(1)); // to be improved: are errors propagated correctly here? 
      if(f2pFinalCorrectionsForNUAPtRP) f2pFinalCorrectionsForNUAPtRP->Add(correctionPt2ndTerm); // to be improved (if condition goes somewhere else)
      delete correctionPt2ndTerm;
      // eta:
      TH1D *correctionEta2ndTerm = new TH1D(*((this->MakeEtaProjection(fCorrectionsSinP1nPsiPtEtaRP))->ProjectionX("","e"))); // to be improved: are errors propagated correctly here? 
      correctionEta2ndTerm->Scale(fQCorrectionsSin->GetBinContent(1)); // to be improved: are errors propagated correctly here? 
      if(f2pFinalCorrectionsForNUAEtaRP) f2pFinalCorrectionsForNUAEtaRP->Add(correctionEta2ndTerm); // to be improved (if condition goes somewhere else)
      delete correctionEta2ndTerm; 
     } else
       { 
        cout<<"WARNING: (fCorrectionsSinP1nPsiPtEtaRP && fQCorrectionsSin) is NULL in QC::CFCFNUAFDF() !!!!  "<<endl;
        cout<<"         Corrections for non-uniform acceptance for differential flow are not correct."<<endl;
       }              
    } else // to else if(type == "RP")
      {
       cout<<"WARNING: Type must be either POI or RP in QC::CFCFNUAFDF() !!!!                           "<<endl;
       cout<<"         Corrections for non-uniform acceptance for differential flow were not calculated."<<endl;
      }  
 } else // to if(!(useParticleWeights))
   {
    // ...
   }
 */
} // end of AliFlowAnalysisWithQCumulants::CalculateFinalCorrectionsForNonUniformAcceptanceForDifferentialFlow(Bool_t useParticleWeights,TString type)


//==================================================================================================================================

/*
void AliFlowAnalysisWithQCumulants::CalculateFinalResultsForDifferentialFlow(
                                                        TH2D *flowPtEta, TH1D *flowPt, TH1D *flowEta, 
                                                        TProfile2D *profile2ndPtEta, TProfile2D *profile4thPtEta, 
                                                        TProfile2D *profile6thPtEta, TProfile2D *profile8thPtEta)
{
 // calculate and store the final results for integrated flow
 
 TString *namePtEta = new TString();
 TString *type = new TString();
 TString *order2nd = new TString();
 TString *order4th = new TString();
 TString *order6th = new TString();
 TString *order8th = new TString(); 
 TString *pW = new TString();

 if(profile2ndPtEta) *namePtEta = profile2ndPtEta->GetName();
 if(namePtEta->Contains("POI")) *type = "POI";
 if(namePtEta->Contains("RP")) *type  = "RP";
 if(namePtEta->Contains("W")) *pW      = "W";
 if(namePtEta->Contains("2")) *order2nd  = "2";
 if(profile4thPtEta) *namePtEta = profile4thPtEta->GetName();
 if(namePtEta->Contains("4")) *order4th = "4";

 if(profile6thPtEta) *namePtEta = profile6thPtEta->GetName();
 if(namePtEta->Contains("6")) *order6th = "6";
 
 if(profile8thPtEta) *namePtEta = profile8thPtEta->GetName();
 if(namePtEta->Contains("8")) *order8th = "8";

 TProfile *profile2ndPt = NULL;
 TProfile *profile4thPt = NULL;
 TProfile *profile6thPt = NULL;
 TProfile *profile8thPt = NULL;

 TProfile *profile2ndEta = NULL;
 TProfile *profile4thEta = NULL;
 TProfile *profile6thEta = NULL;
 TProfile *profile8thEta = NULL;
  
 if(*order2nd == "2")
 {
  profile2ndPt  = new TProfile(*(this->MakePtProjection(profile2ndPtEta))); 
  profile2ndEta = new TProfile(*(this->MakeEtaProjection(profile2ndPtEta))); 
  if(*order4th == "4")
  {
   profile4thPt  = new TProfile(*(this->MakePtProjection(profile4thPtEta))); 
   profile4thEta = new TProfile(*(this->MakeEtaProjection(profile4thPtEta))); 
   if(*order6th == "6")
   {
    profile6thPt  = new TProfile(*(this->MakePtProjection(profile6thPtEta))); 
    profile6thEta = new TProfile(*(this->MakeEtaProjection(profile6thPtEta))); 
    if(*order8th == "8")
    {
     profile8thPt  = new TProfile(*(this->MakePtProjection(profile8thPtEta))); 
     profile8thEta = new TProfile(*(this->MakeEtaProjection(profile8thPtEta))); 
    }     
   }    
  } 
 }
 
 Int_t nBinsPt  = profile2ndPt->GetNbinsX();
 Int_t nBinsEta = profile2ndEta->GetNbinsX();
 
 Double_t dV2 = 0.;
 Double_t dV4 = 0.;
 Double_t dV6 = 0.;
 Double_t dV8 = 0.; 
 
 if(!(*pW == "W"))
 {
  dV2 = fIntFlowResultsQC->GetBinContent(1); 
  dV4 = fIntFlowResultsQC->GetBinContent(2); 
  dV6 = fIntFlowResultsQC->GetBinContent(3); 
  dV8 = fIntFlowResultsQC->GetBinContent(4); 
 } 
 else if(*pW == "W")
 {
  dV2 = fIntFlowResultsQCW->GetBinContent(1);  
  dV4 = fIntFlowResultsQCW->GetBinContent(2); 
  dV6 = fIntFlowResultsQCW->GetBinContent(3); 
  dV8 = fIntFlowResultsQCW->GetBinContent(4); 
 }    
 
 // 3D (pt,eta): 
 Double_t twoPrimePtEta   = 0.; // <<2'>> (pt,eta) 
 Double_t fourPrimePtEta  = 0.; // <<4'>> (pt,eta)  
 //Double_t sixPrimePtEta   = 0.; // <<6'>> (pt,eta) 
 //Double_t eightPrimePtEta = 0.; // <<8'>> (pt,eta) 
 Double_t secondOrderDiffFlowCumulantPtEta = 0.; // d_n{2,Q} (pt,eta)
 Double_t fourthOrderDiffFlowCumulantPtEta = 0.; // d_n{4,Q} (pt,eta) 
 //Double_t sixthOrderDiffFlowCumulantPtEta = 0.; // d_n{6,Q} (pt,eta)
 //Double_t eightOrderDiffFlowCumulantPtEta = 0.; // d_n{8,Q} (pt,eta)2nd
 Double_t dv2PtEta = 0.; // v'_n{2} (pt,eta) 
 Double_t dv4PtEta = 0.; // v'_n{4} (pt,eta) 
 //Double_t dv6PtEta = 0.; // v'_n{6} (pt,eta) 
 //Double_t dv8PtEta = 0.; // v'_n{8} (pt,eta)  

 // 2D (pt):   
 Double_t twoPrimePt   = 0.; // <<2'>> (pt)  
 Double_t fourPrimePt  = 0.; // <<4'>> (pt) 
 //Double_t sixPrimePt   = 0.; // <<6'>> (pt) 
 //Double_t eightPrimePt = 0.; // <<8'>> (pt)          
 Double_t secondOrderDiffFlowCumulantPt = 0.; // d_n{2,Q} (pt) 
 Double_t fourthOrderDiffFlowCumulantPt = 0.; // d_n{4,Q} (pt)  
 //Double_t sixthOrderDiffFlowCumulantPt = 0.; // d_n{6,Q} (pt)
 //Double_t eightOrderDiffFlowCumulantPt = 0.; // d_n{8,Q} (pt)
 Double_t dv2Pt = 0.; // v'_n{2} (pt)
 Double_t dv4Pt = 0.; // v'_n{4} (pt)
 //Double_t dv6Pt = 0.; // v'_n{6} (pt) 
 //Double_t dv8Pt = 0.; // v'_n{8} (pt)  

 // 2D (eta):           
 Double_t twoPrimeEta   = 0.; // <<2'>> (eta)  
 Double_t fourPrimeEta  = 0.; // <<4>> (eta) 
 //Double_t sixPrimeEta   = 0.; // <<6>> (eta) 
 //Double_t eightPrimeEta = 0.; // <<8'>> (eta)  
 Double_t secondOrderDiffFlowCumulantEta = 0.; // d_n{2,Q} (eta)
 Double_t fourthOrderDiffFlowCumulantEta = 0.; // d_n{4,Q} (eta) 
 //Double_t sixthOrderDiffFlowCumulantEta = 0.; // d_n{6,Q} (eta) 
 //Double_t eightOrderDiffFlowCumulantEta = 0.; // d_n{8,Q} (eta) 
 Double_t dv2Eta = 0.; // v'_n{2} (eta)
 Double_t dv4Eta = 0.; // v'_n{4} (eta)
 //Double_t dv6Eta = 0.; // v'_n{6} (eta) 
 //Double_t dv8Eta = 0.; // v'_n{8} (eta)
 

 // looping over (pt,eta) bins to calculate v'(pt,eta) 
 for(Int_t p=1;p<nBinsPt+1;p++)
 {
  for(Int_t e=1;e<nBinsEta+1;e++)
  {
  
   // 2nd order: 
   twoPrimePtEta = profile2ndPtEta->GetBinContent(profile2ndPtEta->GetBin(p,e));
   secondOrderDiffFlowCumulantPtEta = twoPrimePtEta;
   

   //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   // to be improved (applying correction for NUA):
   if(namePtEta->Contains("POI"))
   {
    if(f2pFinalCorrectionsForNUAPtEtaPOI) secondOrderDiffFlowCumulantPtEta = twoPrimePtEta 
                                     - f2pFinalCorrectionsForNUAPtEtaPOI->GetBinContent(f2pFinalCorrectionsForNUAPtEtaPOI->GetBin(p,e)) ;
   } else if (namePtEta->Contains("RP"))
     {  
      if(f2pFinalCorrectionsForNUAPtEtaRP) secondOrderDiffFlowCumulantPtEta = twoPrimePtEta 
                                       - f2pFinalCorrectionsForNUAPtEtaRP->GetBinContent(f2pFinalCorrectionsForNUAPtEtaRP->GetBin(p,e));
     }
   //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   
   if(dV2)
   {
    dv2PtEta = secondOrderDiffFlowCumulantPtEta/dV2;
    if(*order2nd == "2") 
    {
     flowPtEta->SetBinContent(p,e,dv2PtEta);   
    } 
   }
   
   // 4th order: 
   if(*order4th == "4" || *order6th == "6" || *order8th == "8")
   {
    fourPrimePtEta = profile4thPtEta->GetBinContent(profile4thPtEta->GetBin(p,e));
    fourthOrderDiffFlowCumulantPtEta = fourPrimePtEta - 2.*twoPrimePtEta*pow(dV2,2.); // to be improved (correlations instead of pow(dV2,2.))
    if(dV4)
    {
     dv4PtEta = -fourthOrderDiffFlowCumulantPtEta/pow(dV4,3);
     if(*order4th == "4")
     {
      flowPtEta->SetBinContent(p,e,dv4PtEta);
     } 
    }
   }    
   
  } // end of for(Int_t e=1;e<nBinsEta+1;e++)
 } // end of for(Int_t p=1;p<nBinsPt+1;p++) 
   
   
 // looping over (pt) bins to calcualate v'(pt)
 for(Int_t p=1;p<nBinsPt+1;p++)
 {
 
  // 2nd order: 
  twoPrimePt = profile2ndPt->GetBinContent(p);
  secondOrderDiffFlowCumulantPt = twoPrimePt;
  
  
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  // to be improved (applying correction for NUA):
  if(namePtEta->Contains("POI"))
  {
   if(f2pFinalCorrectionsForNUAPtPOI) secondOrderDiffFlowCumulantPt = twoPrimePt
                                    - f2pFinalCorrectionsForNUAPtPOI->GetBinContent(p) ;
  } else if (namePtEta->Contains("RP"))
    {
     if(f2pFinalCorrectionsForNUAPtRP) secondOrderDiffFlowCumulantPt = twoPrimePt
                                      - f2pFinalCorrectionsForNUAPtRP->GetBinContent(p);
    }
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  
  if(dV2)
  {
   dv2Pt = secondOrderDiffFlowCumulantPt/dV2;
   if(*order2nd == "2") 
   {
    flowPt->SetBinContent(p,dv2Pt);
   }
   
   // common control histos: (to be improved fill only once. now they are filled first without weights and then with weights):
   if(namePtEta->Contains("POI") && *order2nd == "2")
   {
    fCommonHistsResults2nd->FillDifferentialFlowPtPOI(p,dv2Pt,0.); //to be improved (errors && bb or bb+1 ?)
   } 
   else if(namePtEta->Contains("RP") && *order2nd == "2")
   {
    fCommonHistsResults2nd->FillDifferentialFlowPtRP(p,dv2Pt,0.); //to be improved (errors && bb or bb+1 ?)
   }
   
  }
  
  // 4th order: 
  if(*order4th == "4" || *order6th == "6" || *order8th == "8")
  {
   fourPrimePt = profile4thPt->GetBinContent(profile4thPt->GetBin(p));
   fourthOrderDiffFlowCumulantPt = fourPrimePt - 2.*twoPrimePt*pow(dV2,2.); // to be improved (correlations instead of pow(dV2,2.))
   if(dV4)
   {
    dv4Pt = -fourthOrderDiffFlowCumulantPt/pow(dV4,3);
    if(*order4th == "4") 
    {
     flowPt->SetBinContent(p,dv4Pt);
    }
    
    // common control histos: (to be improved):
    if(namePtEta->Contains("POI") && *order4th == "4")
    {
     fCommonHistsResults4th->FillDifferentialFlowPtPOI(p,dv4Pt,0.); //to be improved (errors && bb or bb+1 ?)
    } 
    else if(namePtEta->Contains("RP") && *order4th == "4" )
    {
     fCommonHistsResults4th->FillDifferentialFlowPtRP(p,dv4Pt,0.); //to be improved (errors && bb or bb+1 ?)
    }
        
   }
  }    
  
 } // end of for(Int_t p=1;p<nBinsPt+1;p++)  
 
 
 // looping over (eta) bins to calcualate v'(eta)
 for(Int_t e=1;e<nBinsEta+1;e++)
 {
 
  // 2nd order: 
  twoPrimeEta = profile2ndEta->GetBinContent(e);
  secondOrderDiffFlowCumulantEta = twoPrimeEta;
  
  
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  // to be improved (applying correction for NUA):
  if(namePtEta->Contains("POI"))
  {
   if(f2pFinalCorrectionsForNUAEtaPOI) secondOrderDiffFlowCumulantEta = twoPrimeEta
                                    - f2pFinalCorrectionsForNUAEtaPOI->GetBinContent(e) ;
  } else if (namePtEta->Contains("RP"))
    {
     if(f2pFinalCorrectionsForNUAEtaRP) secondOrderDiffFlowCumulantEta = twoPrimeEta
                                      - f2pFinalCorrectionsForNUAEtaRP->GetBinContent(e);
    }
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  
  if(dV2)
  {
   dv2Eta = secondOrderDiffFlowCumulantEta/dV2;
   if(*order2nd == "2") 
   {
    flowEta->SetBinContent(e,dv2Eta);
   }
   
   // common control histos: (to be improved):
   if(namePtEta->Contains("POI") && *order2nd == "2")
   {
    fCommonHistsResults2nd->FillDifferentialFlowEtaPOI(e,dv2Eta,0.); //to be improved (errors && bb or bb+1 ?)
   } 
   else if(namePtEta->Contains("RP") && *order2nd == "2")
   {
    fCommonHistsResults2nd->FillDifferentialFlowEtaRP(e,dv2Eta,0.); //to be improved (errors && bb or bb+1 ?)
   }
     

  }
  
  // 4th order: 
  if(*order4th == "4" || *order6th == "6" || *order8th == "8")
  {
   fourPrimeEta = profile4thEta->GetBinContent(profile4thEta->GetBin(e));
   fourthOrderDiffFlowCumulantEta = fourPrimeEta - 2.*twoPrimeEta*pow(dV2,2.); // to be improved (correlations instead of pow(dV2,2.))
   if(dV4)
   {
    dv4Eta = -fourthOrderDiffFlowCumulantEta/pow(dV4,3);
    if(*order4th == "4")
    {
     flowEta->SetBinContent(e,dv4Eta);
    }
    
    // common control histos: (to be improved):
    if(namePtEta->Contains("POI") && *order4th == "4")
    {
     fCommonHistsResults4th->FillDifferentialFlowEtaPOI(e,dv4Eta,0.); //to be improved (errors && bb or bb+1 ?)
    } 
    else if(namePtEta->Contains("RP") && *order4th == "4")
    {
     fCommonHistsResults4th->FillDifferentialFlowEtaRP(e,dv4Eta,0.); //to be improved (errors && bb or bb+1 ?)
    }
   
   }
  }    
  
 } // end of for(Int_t e=1;e<nBinsEta+1;e++)    
    
 delete namePtEta;
 delete type;
 delete order2nd;
 delete order4th;
 delete order6th;
 delete order8th;
 delete pW;
 delete profile2ndPt;
 delete profile4thPt;
 delete profile6thPt;
 delete profile8thPt;
 delete profile2ndEta;
 delete profile4thEta;
 delete profile6thEta;
 delete profile8thEta;

} // end of AliFlowAnalysisWithQCumulants::CalculateFinalResultsForDifferentialFlow(Bool_t useParticleWeights, TString type)
*/

//================================================================================================================================


void AliFlowAnalysisWithQCumulants::PrintFinalResultsForIntegratedFlow(TString type)
{
 // printing on the screen the final results for integrated flow (NONAME, POI and RP) // to be improved (NONAME) 
 
 Int_t n = fHarmonic; 
 
 if(type == "NONAME" || type == "RP" || type == "POI")
 {
  if(!(fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th))
  {
   cout<<"WARNING: fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th"<<endl;
   cout<<"         is NULL in AFAWQC::PFRFIF() !!!!"<<endl;
  }
 } else
   {
    cout<<"WARNING: type in not from {NONAME, RP, POI} in AFAWQC::PFRFIF() !!!!"<<endl;
    exit(0);
   }
 
 Double_t dVn[4] = {0.}; // array to hold Vn{2}, Vn{4}, Vn{6} and Vn{8}   
 Double_t dVnErr[4] = {0.}; // array to hold errors of Vn{2}, Vn{4}, Vn{6} and Vn{8}   
 
 if(type == "NONAME")
 {
  dVn[0] = (fCommonHistsResults2nd->GetHistIntFlow())->GetBinContent(1); 
  dVnErr[0] = (fCommonHistsResults2nd->GetHistIntFlow())->GetBinError(1); 
  dVn[1] = (fCommonHistsResults4th->GetHistIntFlow())->GetBinContent(1); 
  dVnErr[1] = (fCommonHistsResults4th->GetHistIntFlow())->GetBinError(1); 
  dVn[2] = (fCommonHistsResults6th->GetHistIntFlow())->GetBinContent(1); 
  dVnErr[2] = (fCommonHistsResults6th->GetHistIntFlow())->GetBinError(1); 
  dVn[3] = (fCommonHistsResults8th->GetHistIntFlow())->GetBinContent(1); 
  dVnErr[3] = (fCommonHistsResults8th->GetHistIntFlow())->GetBinError(1); 
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
     }
 
 TString title = " flow estimates from Q-cumulants"; 
 TString subtitle = "    ("; 
 
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
 {
  subtitle.Append(type);
  subtitle.Append(", without weights)");
 } else  
   {
    subtitle.Append(type);
    subtitle.Append(", with weights)");
   }
  
 cout<<endl;
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<title.Data()<<endl; 
 cout<<subtitle.Data()<<endl; 
 cout<<endl;
  
 for(Int_t i=0;i<4;i++)
 {
  if(dVn[i]>=0.)
  {
   cout<<"  v_"<<n<<"{"<<2*(i+1)<<"} = "<<dVn[i]<<" +/- "<<dVnErr[i]<<endl;
  }
  else
  {
   cout<<"  v_"<<n<<"{"<<2*(i+1)<<"} = Im"<<endl;
  }  
 }

 cout<<endl;
 /*
 if(type == "NONAME")
 {
  cout<<"     nEvts = "<<nEvtsNoName<<", AvM = "<<dMultNoName<<endl; // to be improved
 }
 else if (type == "RP")
 {
  cout<<"     nEvts = "<<nEvtsRP<<", AvM = "<<dMultRP<<endl; // to be improved  
 } 
 else if (type == "POI")
 {
  cout<<"     nEvts = "<<nEvtsPOI<<", AvM = "<<dMultPOI<<endl; // to be improved  
 } 
 */
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<endl; 
  
}// end of AliFlowAnalysisWithQCumulants::PrintFinalResultsForIntegratedFlow(TString type="NONAME");


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CompareResultsFromNestedLoopsAndFromQVectorsForDiffFlow(Bool_t useParticleWeights)
{
 // compare correlations needed for diff. flow calculated with nested loops and those calculated from Q-vectors

 cout<<endl;
 cout<<"   *************************************"<<endl;
 cout<<"   **** cross-checking the formulas ****"<<endl;
 cout<<"   ****    for differential flow    ****"<<endl;
 cout<<"   ****                             ****"<<endl;
 cout<<"   ****        (pt,eta) bin:        ****"<<endl; 
 cout<<"   ****    1.1  < pt  <  1.2  GeV   ****"<<endl;  
 cout<<"   ****   -0.55 < eta < -0.525      ****"<<endl; 
 cout<<"   *************************************"<<endl;                             
 cout<<endl;
 
 if(!useParticleWeights)
 {                                        
  cout<<"<cos(n(psi1-phi2))> from Q-vectors    = "<<fCorrelationsPro[0][0][0][0]->GetBinContent(fCorrelationsPro[0][0][0][0]->GetBin(12,19))<<endl;
  cout<<"<cos(n(psi1-phi2))> from nested loops = "<<fDirectCorrelationsDiffFlow->GetBinContent(1)<<endl;
  cout<<endl;  
  cout<<"<cos(n(psi1+phi2-phi3-phi4))> from Q-vectors    = "<<fCorrelationsPro[0][0][0][1]->GetBinContent(fCorrelationsPro[0][0][0][1]->GetBin(12,19))<<endl;
  cout<<"<cos(n(psi1+phi2-phi3-phi4))> from nested loops = "<<fDirectCorrelationsDiffFlow->GetBinContent(41)<<endl;
  cout<<endl;   
  cout<<"****************************************************"<<endl;
  cout<<"****************************************************"<<endl;
  cout<<endl;
  /*
  cout<<"<cos(n(psi1))> from Q-vectors    = "<<fCorrectionsCosP1nPsiPtEtaPOI->GetBinContent(fCorrectionsCosP1nPsiPtEtaPOI->GetBin(12,19))<<endl;
  cout<<"<cos(n(psi1))> from nested loops = "<<fDirectCorrectionsDiffFlowCos->GetBinContent(1)<<endl;
  cout<<endl;  
  cout<<"<sin(n(psi1))> from Q-vectors    = "<<fCorrectionsSinP1nPsiPtEtaPOI->GetBinContent(fCorrectionsSinP1nPsiPtEtaPOI->GetBin(12,19))<<endl;
  cout<<"<sin(n(psi1))> from nested loops = "<<fDirectCorrectionsDiffFlowSin->GetBinContent(1)<<endl;
  cout<<endl;
  */
 }
 
 if(useParticleWeights)
 {
  cout<<"<w2 cos(n(psi1-phi2))> from Q-vectors (RP)   = "<<fCorrelationsPro[0][1][0][0]->GetBinContent(fCorrelationsPro[0][1][0][0]->GetBin(12,19))<<endl;
  cout<<"<w2 cos(n(psi1-phi2))> from Q-vectors (POI)  = "<<fCorrelationsPro[1][1][0][0]->GetBinContent(fCorrelationsPro[1][1][0][0]->GetBin(12,19))<<endl;
  cout<<"<w2 cos(n(psi1-phi2))> from nested loops     = "<<fDirectCorrelationsDiffFlowW->GetBinContent(1)<<endl;
  cout<<endl;  
  cout<<"<w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))> from Q-vectors (RP)  = "<<fCorrelationsPro[0][1][0][1]->GetBinContent(fCorrelationsPro[0][1][0][1]->GetBin(12,19))<<endl;
  cout<<"<w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))> from Q-vectors (POI) = "<<fCorrelationsPro[1][1][0][1]->GetBinContent(fCorrelationsPro[1][1][0][1]->GetBin(12,19))<<endl;
  cout<<"<w2 w3 w4 cos(n(psi1+phi2-phi3-phi4))> from nested loops    = "<<fDirectCorrelationsDiffFlowW->GetBinContent(41)<<endl;
  cout<<endl;  
 }
 
} // end of void AliFlowAnalysisWithQCumulants::CompareResultsFromNestedLoopsAndFromQVectorsForDiffFlow()


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


void AliFlowAnalysisWithQCumulants::BookCommonHistograms()
{
 // book common control histograms and common histograms for final results
 
 // common control histogram (ALL events)
 fCommonHists = new AliFlowCommonHist("AliFlowCommonHistQC");
 fHistList->Add(fCommonHists);  
 
 // common control histogram (for events with 2 and more particles)
 fCommonHists2nd = new AliFlowCommonHist("AliFlowCommonHist2ndOrderQC");
 fHistList->Add(fCommonHists2nd);  
 
 // common control histogram (for events with 4 and more particles)
 fCommonHists4th = new AliFlowCommonHist("AliFlowCommonHist4thOrderQC");
 fHistList->Add(fCommonHists4th);  
 
 // common control histogram (for events with 6 and more particles)
 fCommonHists6th = new AliFlowCommonHist("AliFlowCommonHist6thOrderQC");
 fHistList->Add(fCommonHists6th);  
 
 // common control histogram (for events with 8 and more particles)
 fCommonHists8th = new AliFlowCommonHist("AliFlowCommonHist8thOrderQC");
 fHistList->Add(fCommonHists8th);  
  
 // common histograms for final results (calculated for events with 2 and more particles)
 fCommonHistsResults2nd = new AliFlowCommonHistResults("AliFlowCommonHistResults2ndOrderQC");
 fHistList->Add(fCommonHistsResults2nd); 
 
 // common histograms for final results (calculated for events with 4 and more particles)
 fCommonHistsResults4th = new AliFlowCommonHistResults("AliFlowCommonHistResults4thOrderQC");
 fHistList->Add(fCommonHistsResults4th);
 
 // common histograms for final results (calculated for events with 6 and more particles)
 fCommonHistsResults6th = new AliFlowCommonHistResults("AliFlowCommonHistResults6thOrderQC");
 fHistList->Add(fCommonHistsResults6th); 
 
 // common histograms for final results (calculated for events with 8 and more particles)
 fCommonHistsResults8th = new AliFlowCommonHistResults("AliFlowCommonHistResults8thOrderQC");
 fHistList->Add(fCommonHistsResults8th); 
 
} // end of void AliFlowAnalysisWithQCumulants::BookCommonHistograms()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::BookAndFillWeightsHistograms()
{
 // book and fill histograms which hold phi, pt and eta weights

 if(!fWeightsList)
 {
  cout<<"WARNING: fWeightsList is NULL in AFAWQC::BAFWH() !!!!"<<endl;
  exit(0);  
 }
    
 TString fUseParticleWeightsName = "fUseParticleWeightsQC";
 fUseParticleWeightsName += fAnalysisLabel->Data();
 fUseParticleWeights = new TProfile(fUseParticleWeightsName.Data(),"0 = particle weight not used, 1 = particle weight used ",3,0,3);
 fUseParticleWeights->SetLabelSize(0.06);
 (fUseParticleWeights->GetXaxis())->SetBinLabel(1,"w_{#phi}");
 (fUseParticleWeights->GetXaxis())->SetBinLabel(2,"w_{p_{T}}");
 (fUseParticleWeights->GetXaxis())->SetBinLabel(3,"w_{#eta}");
 fUseParticleWeights->Fill(0.5,(Int_t)fUsePhiWeights);
 fUseParticleWeights->Fill(1.5,(Int_t)fUsePtWeights);
 fUseParticleWeights->Fill(2.5,(Int_t)fUseEtaWeights);
 fWeightsList->Add(fUseParticleWeights); 
  
 if(fUsePhiWeights)
 {
  if(fWeightsList->FindObject("phi_weights"))
  {
   fPhiWeights = dynamic_cast<TH1F*>(fWeightsList->FindObject("phi_weights"));
   if(fPhiWeights->GetBinWidth(1) != fPhiBinWidth)
   {
    cout<<"WARNING: fPhiWeights->GetBinWidth(1) != fPhiBinWidth in AFAWQC::BAFWH() !!!!        "<<endl;
    cout<<"         This indicates inconsistent binning in phi histograms throughout the code."<<endl;
    exit(0);
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
   if(fPtWeights->GetBinWidth(1) != fPtBinWidth)
   {
    cout<<"WARNING: fPtWeights->GetBinWidth(1) != fPtBinWidth in AFAWQC::BAFWH() !!!!         "<<endl;
    cout<<"         This indicates insconsistent binning in pt histograms throughout the code."<<endl;
    exit(0);
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
   if(fEtaWeights->GetBinWidth(1) != fEtaBinWidth)
   {
    cout<<"WARNING: fEtaWeights->GetBinWidth(1) != fEtaBinWidth in AFAWQC::BAFWH() !!!!        "<<endl;
    cout<<"         This indicates insconsistent binning in eta histograms throughout the code."<<endl;
    exit(0);
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
 // book all objects for integrated flow 
 
 // a) common;
 // b) profiles;
 // c) results;
 
 // ****************
 // **** COMMON ****
 // **************** 

 // Re[Q_{m*n,k}], Im[Q_{m*n,k}] and S_{p,k}^M: 
 fReQ  = new TMatrixD(4,9);
 fImQ  = new TMatrixD(4,9);
 fSMpk = new TMatrixD(8,9);
 
 // profile to hold average multiplicities and number of events for events with nRP>=0, nRP>=1, ... , and nRP>=8
 fAvMultiplicity = new TProfile("fAvMultiplicity","Average Multiplicities of RPs",9,0,9);
 fAvMultiplicity->SetTickLength(-0.01,"Y");
 fAvMultiplicity->SetMarkerStyle(25);
 fAvMultiplicity->SetLabelSize(0.05);
 fAvMultiplicity->SetLabelOffset(0.02,"Y");
 fAvMultiplicity->SetYTitle("Average Multiplicity");
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
  
 TString pWeightsFlag[2] = {"pWeights not used","pWeights used"};
 TString eWeightsFlag[2] = {"exact eWeights","non-exact eWeights"};
 TString nuaFlag[2] = {"not corrected","corrected"};
 TString sinCosFlag[2] = {"sin","cos"};
 TString powerFlag[2] = {"linear","quadratic"};
  
 for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++)
 {
  // ***************
  // **** E-B-E ****
  // ***************
  
  // average multiparticle correlations for single event calculated from Q-vectors
  // (Remark: binning is organized in the same way as in fQCorrelations[pW][eW] bellow):   
  fQCorrelationsEBE[pW] = new TH1D(Form("fQCorrelationsEBE: %s",pWeightsFlag[pW].Data()),Form("Average multi-particle correlations for single event calculated from Q-vectors (%s)",pWeightsFlag[pW].Data()),32,0,32);
  
  for(Int_t sc=0;sc<2;sc++)
  {
   // correction terms for non-uniform acceptance for single event calculated from Q-vectors
   // (Remark: binning is organized in the same way as in fQCorrections[pW][sc]):
   fQCorrectionsEBE[pW][sc] = new TH1D(Form("fQCorrectionsEBE: %s, %s terms",pWeightsFlag[pW].Data(),sinCosFlag[sc].Data()),Form("Correction terms for non-uniform acceptance for single event (%s, %s terms)",pWeightsFlag[pW].Data(),sinCosFlag[sc].Data()),10,0,10);  
  }
  
  for(Int_t eW=0;eW<2;eW++)
  {
   // ******************
   // **** PROFILES ****
   // ******************
  
   // final average multiparticle correlations for all events calculated from Q-vectors:
   fQCorrelations[pW][eW] = new TProfile(Form("fQCorrelations: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),Form("Final average multi-particle correlations calculated from Q-vectors (%s, %s)",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),32,0,32,"s");
   fQCorrelations[pW][eW]->SetTickLength(-0.01,"Y");
   fQCorrelations[pW][eW]->SetMarkerStyle(25);
   fQCorrelations[pW][eW]->SetLabelSize(0.03);
   fQCorrelations[pW][eW]->SetLabelOffset(0.01,"Y");
   // 2-p correlations:
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(1,"<<2>>_{n|n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(2,"<<2>>_{2n|2n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(3,"<<2>>_{3n|3n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(4,"<<2>>_{4n|4n}");
   // 3-p correlations:
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(6,"<<3>>_{2n|n,n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(7,"<<3>>_{3n|2n,n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(8,"<<3>>_{4n|2n,2n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(9,"<<3>>_{4n|3n,n}");
   // 4-p correlations:
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(11,"<<4>>_{n,n|n,n}"); 
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(12,"<<4>>_{2n,n|2n,n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(13,"<<4>>_{2n,2n|2n,2n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(14,"<<4>>_{3n|n,n,n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(15,"<<4>>_{3n,n|3n,n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(16,"<<4>>_{3n,n|2n,2n}"); 
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(17,"<<4>>_{4n|2n,n,n}");
   // 5-p correlations:
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(19,"<<5>>_{2n|n,n,n,n}"); 
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(20,"<<5>>_{2n,2n|2n,n,n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(21,"<<5>>_{3n,n|2n,n,n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(22,"<<5>>_{4n|n,n,n,n}");
   // 6-p correlations:
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(24,"<<6>>_{n,n,n|n,n,n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(25,"<<6>>_{2n,n,n|2n,n,n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(26,"<<6>>_{2n,2n|n,n,n,n}");
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(27,"<<6>>_{3n,n|n,n,n,n}");
   // 7-p correlations:
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(29,"<<7>>_{2n,n,n|n,n,n,n}");
   // 8-p correlations:
   (fQCorrelations[pW][eW]->GetXaxis())->SetBinLabel(31,"<<8>>_{n,n,n,n|n,n,n,n}");
   // add fQCorrelations[0] to the list fIntFlowList:
   fIntFlowProfiles->Add(fQCorrelations[pW][eW]);
  
   // averages <<2><4>>, <<2><6>>, <<4><6>>, etc,  needed to calculate covariances:
   fQProducts[pW][eW] = new TProfile(Form("fQProducts: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),Form("Averages of products (%s, %s)",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),6,0,6);
   fQProducts[pW][eW]->SetTickLength(-0.01,"Y");
   fQProducts[pW][eW]->SetMarkerStyle(25);
   fQProducts[pW][eW]->SetLabelSize(0.05);
   fQProducts[pW][eW]->SetLabelOffset(0.01,"Y");
   (fQProducts[pW][eW]->GetXaxis())->SetBinLabel(1,"<<2><4>>");
   (fQProducts[pW][eW]->GetXaxis())->SetBinLabel(2,"<<2><6>>");
   (fQProducts[pW][eW]->GetXaxis())->SetBinLabel(3,"<<2><8>>");
   (fQProducts[pW][eW]->GetXaxis())->SetBinLabel(4,"<<4><6>>");
   (fQProducts[pW][eW]->GetXaxis())->SetBinLabel(5,"<<4><8>>");
   (fQProducts[pW][eW]->GetXaxis())->SetBinLabel(6,"<<6><8>>");
   fIntFlowProfiles->Add(fQProducts[pW][eW]);
 
   for(Int_t sc=0;sc<2;sc++)
   {
    // final average correction terms for non-uniform acceptance calculated from Q-vectors: 
    fQCorrections[pW][eW][sc] = new TProfile(Form("fQCorrections: %s, %s, %s terms",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),sinCosFlag[sc].Data()),Form("Correction terms for non-uniform acceptance (%s, %s, %s terms)",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),sinCosFlag[sc].Data()),10,0,10,"s");
    fQCorrections[pW][eW][sc]->SetTickLength(-0.01,"Y");
    fQCorrections[pW][eW][sc]->SetMarkerStyle(25);
    fQCorrections[pW][eW][sc]->SetLabelSize(0.03);
    fQCorrections[pW][eW][sc]->SetLabelOffset(0.01,"Y");
    // ......................................................................... 
    // 1-p terms:
    (fQCorrections[pW][eW][sc]->GetXaxis())->SetBinLabel(1,Form("%s(n(#phi_{1}))>",sinCosFlag[sc].Data()));
    // 2-p terms:
    // 3-p terms:
    // ...
    // ......................................................................... 
    // add fQCorrectionsCos to the list fIntFlowList:
    fIntFlowProfiles->Add(fQCorrections[pW][eW][sc]);
   } // end of for(Int_t sc=0;sc<2;sc++) 
  
   // *****************
   // **** RESULTS ****
   // *****************
  
   // final results for average multi-particle correlations with correct statistical errors:
   fCorrelations[pW][eW] = new TH1D(Form("fCorrelations: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),Form("Average multi-particle correlations (%s, %s)",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),32,0,32);
   fCorrelations[pW][eW]->SetTickLength(-0.01,"Y");
   fCorrelations[pW][eW]->SetMarkerStyle(25);
   fCorrelations[pW][eW]->SetLabelSize(0.03);
   fCorrelations[pW][eW]->SetLabelOffset(0.01,"Y");
   // 2-p correlations:
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(1,"<<2>>_{n|n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(2,"<<2>>_{2n|2n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(3,"<<2>>_{3n|3n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(4,"<<2>>_{4n|4n}");
   // 3-p correlations:
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(6,"<<3>>_{2n|n,n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(7,"<<3>>_{3n|2n,n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(8,"<<3>>_{4n|2n,2n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(9,"<<3>>_{4n|3n,n}");
   // 4-p correlations:
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(11,"<<4>>_{n,n|n,n}"); 
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(12,"<<4>>_{2n,n|2n,n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(13,"<<4>>_{2n,2n|2n,2n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(14,"<<4>>_{3n|n,n,n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(15,"<<4>>_{3n,n|3n,n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(16,"<<4>>_{3n,n|2n,2n}"); 
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(17,"<<4>>_{4n|2n,n,n}");
   // 5-p correlations:
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(19,"<<5>>_{2n|n,n,n,n}"); 
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(20,"<<5>>_{2n,2n|2n,n,n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(21,"<<5>>_{3n,n|2n,n,n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(22,"<<5>>_{4n|n,n,n,n}");
   // 6-p correlations:
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(24,"<<6>>_{n,n,n|n,n,n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(25,"<<6>>_{2n,n,n|2n,n,n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(26,"<<6>>_{2n,2n|n,n,n,n}");
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(27,"<<6>>_{3n,n|n,n,n,n}");
   // 7-p correlations:
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(29,"<<7>>_{2n,n,n|n,n,n,n}");
   // 8-p correlations:
   (fCorrelations[pW][eW]->GetXaxis())->SetBinLabel(31,"<<8>>_{n,n,n,n|n,n,n,n}");
   // add fCorrelations to the list fIntFlowList:
   fIntFlowResults->Add(fCorrelations[pW][eW]);
 
   // final corrections for non-uniform acceptance for QC{2}, QC{4}, QC{6} and QC{8}:
   fCorrections[pW][eW] = new TH1D(Form("fCorrections: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),Form("Corrections for non-uniform acceptance for Q-cumulants (%s, %s)",eWeightsFlag[eW].Data()),4,0,4);
   fCorrections[pW][eW]->SetLabelSize(0.05);
   fCorrections[pW][eW]->SetMarkerStyle(25);
   (fCorrections[pW][eW]->GetXaxis())->SetBinLabel(1,"corr. for QC{2}");
   (fCorrections[pW][eW]->GetXaxis())->SetBinLabel(2,"corr. for QC{4}");
   (fCorrections[pW][eW]->GetXaxis())->SetBinLabel(3,"corr. for QC{6}");
   (fCorrections[pW][eW]->GetXaxis())->SetBinLabel(4,"corr. for QC{8}");
   // add fCorrections[pW] to list fIntFlowResults: 
   fIntFlowResults->Add(fCorrections[pW][eW]);
  
   // final results for covariances:
   fCovariances[pW][eW] = new TH1D(Form("fCovariances: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),Form("Covariances (%s, %s)",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),6,0,6);
   fCovariances[pW][eW]->SetLabelSize(0.05);
   fCovariances[pW][eW]->SetMarkerStyle(25);
   (fCovariances[pW][eW]->GetXaxis())->SetBinLabel(1,"Cov(2,4)");
   (fCovariances[pW][eW]->GetXaxis())->SetBinLabel(2,"Cov(2,6)");
   (fCovariances[pW][eW]->GetXaxis())->SetBinLabel(3,"Cov(2,8)");
   (fCovariances[pW][eW]->GetXaxis())->SetBinLabel(4,"Cov(4,6)");
   (fCovariances[pW][eW]->GetXaxis())->SetBinLabel(5,"Cov(4,8)");
   (fCovariances[pW][eW]->GetXaxis())->SetBinLabel(6,"Cov(6,8)");
   // add fCovariances[pW][eW] to list fIntFlowResults: 
   fIntFlowResults->Add(fCovariances[pW][eW]);
  
   // final results for sum of event weights:
   for(Int_t power=0;power<2;power++)
   {
    fSumOfEventWeights[pW][eW][power] = new TH1D(Form("fSumOfEventWeights: %s, %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),powerFlag[power].Data()),Form("SumOfEventWeights (%s, %s, %s)",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),powerFlag[power].Data()),4,0,4);
    fSumOfEventWeights[pW][eW][power]->SetLabelSize(0.05);
    fSumOfEventWeights[pW][eW][power]->SetMarkerStyle(25);
    (fSumOfEventWeights[pW][eW][power]->GetXaxis())->SetBinLabel(1,"#sum_{i=1}^{N} w_{<2>}");
    (fSumOfEventWeights[pW][eW][power]->GetXaxis())->SetBinLabel(2,"#sum_{i=1}^{N} w_{<4>}");
    (fSumOfEventWeights[pW][eW][power]->GetXaxis())->SetBinLabel(3,"#sum_{i=1}^{N} w_{<6>}");
    (fSumOfEventWeights[pW][eW][power]->GetXaxis())->SetBinLabel(4,"#sum_{i=1}^{N} w_{<8>}");
    // add fSumOfEventWeights[pW][eW] to list fIntFlowResults: 
    fIntFlowResults->Add(fSumOfEventWeights[pW][eW][power]);
   } 
   
   // final results for sum of product of event weights:
   fProductOfEventWeights[pW][eW] = new TH1D(Form("fProductOfEventWeights: %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),Form("ProductOfEventWeights (%s, %s)",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()),6,0,6);
   fProductOfEventWeights[pW][eW]->SetLabelSize(0.05);
   fProductOfEventWeights[pW][eW]->SetMarkerStyle(25);
   (fProductOfEventWeights[pW][eW]->GetXaxis())->SetBinLabel(1,"#sum_{i=1}^{N} w_{<2>} w_{<4>}");
   (fProductOfEventWeights[pW][eW]->GetXaxis())->SetBinLabel(2,"#sum_{i=1}^{N} w_{<2>} w_{<6>}");
   (fProductOfEventWeights[pW][eW]->GetXaxis())->SetBinLabel(3,"#sum_{i=1}^{N} w_{<2>} w_{<8>}");
   (fProductOfEventWeights[pW][eW]->GetXaxis())->SetBinLabel(4,"#sum_{i=1}^{N} w_{<4>} w_{<6>}");
   (fProductOfEventWeights[pW][eW]->GetXaxis())->SetBinLabel(5,"#sum_{i=1}^{N} w_{<4>} w_{<8>}");
   (fProductOfEventWeights[pW][eW]->GetXaxis())->SetBinLabel(6,"#sum_{i=1}^{N} w_{<6>} w_{<8>}");
   // add fProductOfEventWeights[pW][eW] to list fIntFlowResults: 
   fIntFlowResults->Add(fProductOfEventWeights[pW][eW]);

   for(Int_t nua=0;nua<2;nua++)
   {
    // integrated Q-cumulants:
    fCumulants[pW][eW][nua] = new TH1D(Form("fCumulants: %s, %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data()),Form("Q-cumulants (%s, %s, %s)",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data()),4,0,4);
    fCumulants[pW][eW][nua]->SetLabelSize(0.05);
    fCumulants[pW][eW][nua]->SetMarkerStyle(25);
    (fCumulants[pW][eW][nua]->GetXaxis())->SetBinLabel(1,"QC{2}");
    (fCumulants[pW][eW][nua]->GetXaxis())->SetBinLabel(2,"QC{4}");
    (fCumulants[pW][eW][nua]->GetXaxis())->SetBinLabel(3,"QC{6}");
    (fCumulants[pW][eW][nua]->GetXaxis())->SetBinLabel(4,"QC{8}");
    // add fCumulants[pW][nua] to list fIntFlowResults: 
    fIntFlowResults->Add(fCumulants[pW][eW][nua]);
   
    // integrated flow from Q-cumulants:
    fIntFlow[pW][eW][nua] = new TH1D(Form("fIntFlow: %s, %s, %s",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data()),Form("Integrated flow from Q-cumulants (%s, %s, %s)",pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data()),4,0,4);
    fIntFlow[pW][eW][nua]->SetLabelSize(0.05);
    fIntFlow[pW][eW][nua]->SetMarkerStyle(25);
    (fIntFlow[pW][eW][nua]->GetXaxis())->SetBinLabel(1,"v_{2}{2,QC}");
    (fIntFlow[pW][eW][nua]->GetXaxis())->SetBinLabel(2,"v_{2}{4,QC}");
    (fIntFlow[pW][eW][nua]->GetXaxis())->SetBinLabel(3,"v_{2}{6,QC}");
    (fIntFlow[pW][eW][nua]->GetXaxis())->SetBinLabel(4,"v_{2}{8,QC}");
    // add fIntFlow[pW][nua] to list fIntFlowResults: 
    fIntFlowResults->Add(fIntFlow[pW][eW][nua]);
   
   } // end of for(Int_t nua=0;nua<2;nua++)
  } // end of for(Int_t eW=0;eW<2;eW++)
 } // end of for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++)
  


 /*
  // final average weighted multi-particle correlations for all events calculated from Q-vectors
  fQCorrelations[1] = new TProfile("Weighted correlations","final average multi-particle correlations from weighted Q-vectors",200,0,200,"s");
  fQCorrelations[1]->SetTickLength(-0.01,"Y");
  fQCorrelations[1]->SetMarkerStyle(25);
  fQCorrelations[1]->SetLabelSize(0.03);
  fQCorrelations[1]->SetLabelOffset(0.01,"Y");
  // 2-particle correlations:
  (fQCorrelations[1]->GetXaxis())->SetBinLabel(1,"<w_{1}w_{2}cos(n(#phi_{1}-#phi_{2}))>");
  (fQCorrelations[1]->GetXaxis())->SetBinLabel(2,"<w_{1}^{2}w_{2}^{2}cos(2n(#phi_{1}-#phi_{2}))>");
  (fQCorrelations[1]->GetXaxis())->SetBinLabel(3,"<w_{1}^{3}w_{2}^{3}cos(3n(#phi_{1}-#phi_{2}))>");
  (fQCorrelations[1]->GetXaxis())->SetBinLabel(4,"<w_{1}^{4}w_{2}^{4}cos(4n(#phi_{1}-#phi_{2}))>");
  (fQCorrelations[1]->GetXaxis())->SetBinLabel(5,"<w_{1}^{3}w_{2}cos(n(#phi_{1}-#phi_{2}))>");
  (fQCorrelations[1]->GetXaxis())->SetBinLabel(6,"<w_{1}^{2}w_{2}w_{3}cos(n(#phi_{1}-#phi_{2}))>");
  // 3-particle correlations:
  (fQCorrelations[1]->GetXaxis())->SetBinLabel(21,"<w_{1}w_{2}w_{3}^{2}cos(n(2#phi_{1}-#phi_{2}-#phi_{3}))>");
  // 4-particle correlations:
  (fQCorrelations[1]->GetXaxis())->SetBinLabel(41,"<w_{1}w_{2}w_{3}w_{4}cos(n(#phi_{1}+#phi_{2}-#phi_{3}-#phi_{4}))>");
  // add fQCorrelations[1] to the list fIntFlowList:
  fIntFlowList->Add(fQCorrelations[1]); 
 */
  
 
} // end of AliFlowAnalysisWithQCumulants::BookEverythingForIntegratedFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::BookEverythingForNestedLoops()
{
 // book all profiles relevant for calculations with nested loops
 
 fEvaluateNestedLoops = new TProfile("fEvaluateNestedLoops","1 = evaluate, 0 = do not evaluate",2,0,2);
 fEvaluateNestedLoops->SetLabelSize(0.05);
 (fEvaluateNestedLoops->GetXaxis())->SetBinLabel(1,"Nested Loops (Int. Flow)");
 (fEvaluateNestedLoops->GetXaxis())->SetBinLabel(2,"Nested Loops (Diff. Flow)");
 fEvaluateNestedLoops->Fill(0.5,(Int_t)fEvaluateNestedLoopsForIntFlow);
 fEvaluateNestedLoops->Fill(1.5,(Int_t)fEvaluateNestedLoopsForDiffFlow);
 // add fEvaluateNestedLoops to the list fNestedLoopsList:
 fNestedLoopsList->Add(fEvaluateNestedLoops);
 
 if(fEvaluateNestedLoopsForIntFlow)
 {
  fDirectCorrelations = new TProfile("fDirectCorrelations","multi-particle correlations with nested loops",100,0,100,"s");
  fNestedLoopsList->Add(fDirectCorrelations);
  fDirectCorrectionsCos = new TProfile("fDirectCorrectionsCos"," corrections for non-uniform acceptance (cos terms)",100,0,100,"s");
  fNestedLoopsList->Add(fDirectCorrectionsCos);
  fDirectCorrectionsSin = new TProfile("fDirectCorrectionsSin"," corrections for non-uniform acceptance (sin terms)",100,0,100,"s");
  fNestedLoopsList->Add(fDirectCorrectionsSin); 
  if(fUsePhiWeights) // Remark: cross-checking performed only with phi-weights (this is sufficient)
  {
   fDirectCorrelationsW = new TProfile("fDirectCorrelationsW","multi-particle correlations with nested loops",200,0,200,"s"); 
   fNestedLoopsList->Add(fDirectCorrelationsW);
   fDirectCorrectionsCosW = new TProfile("fDirectCorrectionsCosW"," corrections for non-uniform acceptance (cos terms)",100,0,100,"s");
   fNestedLoopsList->Add(fDirectCorrectionsCosW);
   fDirectCorrectionsSinW = new TProfile("fDirectCorrectionsSinW"," corrections for non-uniform acceptance (sin terms)",100,0,100,"s");
   fNestedLoopsList->Add(fDirectCorrectionsSinW); 
  }
 }
 
 if(fEvaluateNestedLoopsForDiffFlow)
 {
  fDirectCorrelationsDiffFlow = new TProfile("fDirectCorrelationsDiffFlow","multi-particle correlations with nested loops",200,0,200,"s");
  fNestedLoopsList->Add(fDirectCorrelationsDiffFlow); 
  fDirectCorrectionsDiffFlowCos = new TProfile("fDirectCorrectionsDiffFlowCos",
                                               "corrections for non-uniform acceptance (cos terms) with nested loops",200,0,200,"s");
  fNestedLoopsList->Add(fDirectCorrectionsDiffFlowCos);   
  fDirectCorrectionsDiffFlowSin = new TProfile("fDirectCorrectionsDiffFlowSin",
                                               "corrections for non-uniform acceptance (sin terms) with nested loops",200,0,200,"s");
  fNestedLoopsList->Add(fDirectCorrectionsDiffFlowSin);           
  if(fUsePhiWeights) // Remark: cross-checking performed only with phi-weights (this is sufficient)
  {
   fDirectCorrelationsDiffFlowW = new TProfile("fDirectCorrelationsDiffFlowW","multi-particle correlations with nested loops",200,0,200,"s");
   fNestedLoopsList->Add(fDirectCorrelationsDiffFlowW);
   fDirectCorrectionsDiffFlowCosW = new TProfile("fDirectCorrectionsDiffFlowCosW",
                                               "corrections for non-uniform acceptance (cos terms) with nested loops",200,0,200,"s");
   fNestedLoopsList->Add(fDirectCorrectionsDiffFlowCosW);   
   fDirectCorrectionsDiffFlowSinW = new TProfile("fDirectCorrectionsDiffFlowSinW",
                                               "corrections for non-uniform acceptance (sin terms) with nested loops",200,0,200,"s");
   fNestedLoopsList->Add(fDirectCorrectionsDiffFlowSinW);   
  }
 }
 
} // end of AliFlowAnalysisWithQCumulants::BookEverythingForNestedLoops()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCorrelationsForIntegratedFlow()
{
 // calculate all correlations needed for integrated flow
 
 // multiplicity:
 Double_t dMult = (*fSMpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 Double_t dReQ3n = (*fReQ)(2,0);
 Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 Double_t dImQ3n = (*fImQ)(2,0);
 Double_t dImQ4n = (*fImQ)(3,0);
  
 // real and imaginary parts of some expressions involving various combinations of Q-vectors evaluated in harmonics n, 2n, 3n and 4n:
 // (these expression appear in the Eqs. for the multi-particle correlations bellow)
 
 // Re[Q_{2n} Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ1nstarQ1nstar = pow(dReQ1n,2.)*dReQ2n + 2.*dReQ1n*dImQ1n*dImQ2n - pow(dImQ1n,2.)*dReQ2n; 
 
 // Im[Q_{2n} Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ1nstarQ1nstar = pow(dReQ1n,2.)*dImQ2n-2.*dReQ1n*dImQ1n*dReQ2n-pow(dImQ1n,2.)*dImQ2n; 
 
 // Re[Q_{n} Q_{n} Q_{2n}^*] = Re[Q_{2n} Q_{n}^* Q_{n}^*]
 Double_t reQ1nQ1nQ2nstar = reQ2nQ1nstarQ1nstar; 
 
 // Re[Q_{3n} Q_{n} Q_{2n}^* Q_{2n}^*]
 Double_t reQ3nQ1nQ2nstarQ2nstar = (pow(dReQ2n,2.)-pow(dImQ2n,2.))*(dReQ3n*dReQ1n-dImQ3n*dImQ1n) 
                                 + 2.*dReQ2n*dImQ2n*(dReQ3n*dImQ1n+dImQ3n*dReQ1n);

 // Im[Q_{3n} Q_{n} Q_{2n}^* Q_{2n}^*]                                                                  
 //Double_t imQ3nQ1nQ2nstarQ2nstar = calculate and implement this (deleteMe)
  
 // Re[Q_{2n} Q_{2n} Q_{3n}^* Q_{1n}^*] = Re[Q_{3n} Q_{n} Q_{2n}^* Q_{2n}^*]
 Double_t reQ2nQ2nQ3nstarQ1nstar = reQ3nQ1nQ2nstarQ2nstar;
  
 // Re[Q_{4n} Q_{2n}^* Q_{2n}^*]
 Double_t reQ4nQ2nstarQ2nstar = pow(dReQ2n,2.)*dReQ4n+2.*dReQ2n*dImQ2n*dImQ4n-pow(dImQ2n,2.)*dReQ4n;

 // Im[Q_{4n} Q_{2n}^* Q_{2n}^*]
 //Double_t imQ4nQ2nstarQ2nstar = calculate and implement this (deleteMe)
 
 // Re[Q_{2n} Q_{2n} Q_{4n}^*] =  Re[Q_{4n} Q_{2n}^* Q_{2n}^*]
 Double_t reQ2nQ2nQ4nstar = reQ4nQ2nstarQ2nstar;
 
 // Re[Q_{4n} Q_{3n}^* Q_{n}^*]
 Double_t reQ4nQ3nstarQ1nstar = dReQ4n*(dReQ3n*dReQ1n-dImQ3n*dImQ1n)+dImQ4n*(dReQ3n*dImQ1n+dImQ3n*dReQ1n);
 
 // Re[Q_{3n} Q_{n} Q_{4n}^*] = Re[Q_{4n} Q_{3n}^* Q_{n}^*]
 Double_t reQ3nQ1nQ4nstar = reQ4nQ3nstarQ1nstar;
 
 // Im[Q_{4n} Q_{3n}^* Q_{n}^*]
 //Double_t imQ4nQ3nstarQ1nstar = calculate and implement this (deleteMe)

 // Re[Q_{3n} Q_{2n}^* Q_{n}^*]
 Double_t reQ3nQ2nstarQ1nstar = dReQ3n*dReQ2n*dReQ1n-dReQ3n*dImQ2n*dImQ1n+dImQ3n*dReQ2n*dImQ1n
                              + dImQ3n*dImQ2n*dReQ1n;
                              
 // Re[Q_{2n} Q_{n} Q_{3n}^*] = Re[Q_{3n} Q_{2n}^* Q_{n}^*]
 Double_t reQ2nQ1nQ3nstar = reQ3nQ2nstarQ1nstar;
 
 // Im[Q_{3n} Q_{2n}^* Q_{n}^*]
 //Double_t imQ3nQ2nstarQ1nstar; //calculate and implement this (deleteMe)
 
 // Re[Q_{3n} Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ3nQ1nstarQ1nstarQ1nstar = dReQ3n*pow(dReQ1n,3)-3.*dReQ1n*dReQ3n*pow(dImQ1n,2)
                                     + 3.*dImQ1n*dImQ3n*pow(dReQ1n,2)-dImQ3n*pow(dImQ1n,3);

 // Im[Q_{3n} Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ3nQ1nstarQ1nstarQ1nstar; //calculate and implement this (deleteMe)
 
 // |Q_{2n}|^2 |Q_{n}|^2
 Double_t dQ2nQ1nQ2nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.));
 
 // Re[Q_{4n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ4nQ2nstarQ1nstarQ1nstar = (dReQ4n*dReQ2n+dImQ4n*dImQ2n)*(pow(dReQ1n,2)-pow(dImQ1n,2))
                                     + 2.*dReQ1n*dImQ1n*(dImQ4n*dReQ2n-dReQ4n*dImQ2n); 
 
 // Im[Q_{4n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ4nQ2nstarQ1nstarQ1nstar; //calculate and implement this (deleteMe)
 
 // Re[Q_{2n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ1nQ1nstarQ1nstarQ1nstar = (dReQ2n*dReQ1n-dImQ2n*dImQ1n)*(pow(dReQ1n,3)-3.*dReQ1n*pow(dImQ1n,2))
                                        + (dReQ2n*dImQ1n+dReQ1n*dImQ2n)*(3.*dImQ1n*pow(dReQ1n,2)-pow(dImQ1n,3));

 // Im[Q_{2n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^*] 
 //Double_t imQ2nQ1nQ1nstarQ1nstarQ1nstar; //calculate and implement this (deleteMe)
 
 // Re[Q_{2n} Q_{2n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ2nQ2nstarQ1nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))
                                        * (dReQ2n*(pow(dReQ1n,2.)-pow(dImQ1n,2.)) + 2.*dImQ2n*dReQ1n*dImQ1n);

 // Im[Q_{2n} Q_{2n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ2nQ2nstarQ1nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))
 //                                       * (dImQ2n*(pow(dReQ1n,2.)-pow(dImQ1n,2.)) - 2.*dReQ2n*dReQ1n*dImQ1n);
 
 // Re[Q_{4n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ4nQ1nstarQ1nstarQ1nstarQ1nstar = pow(dReQ1n,4.)*dReQ4n-6.*pow(dReQ1n,2.)*dReQ4n*pow(dImQ1n,2.)
                                            + pow(dImQ1n,4.)*dReQ4n+4.*pow(dReQ1n,3.)*dImQ1n*dImQ4n
                                            - 4.*pow(dImQ1n,3.)*dReQ1n*dImQ4n;
                                            
 // Im[Q_{4n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ4nQ1nstarQ1nstarQ1nstarQ1nstar = pow(dReQ1n,4.)*dImQ4n-6.*pow(dReQ1n,2.)*dImQ4n*pow(dImQ1n,2.)
 //                                           + pow(dImQ1n,4.)*dImQ4n+4.*pow(dImQ1n,3.)*dReQ1n*dReQ4n
 //                                           - 4.*pow(dReQ1n,3.)*dImQ1n*dReQ4n;
 
 // Re[Q_{3n} Q_{n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ3nQ1nQ2nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)+pow(dImQ1n,2.))
                                        * (dReQ1n*dReQ2n*dReQ3n-dReQ3n*dImQ1n*dImQ2n+dReQ2n*dImQ1n*dImQ3n+dReQ1n*dImQ2n*dImQ3n);
 
 // Im[Q_{3n} Q_{n} Q_{2n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ3nQ1nQ2nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)+pow(dImQ1n,2.))
 //                                       * (-dReQ2n*dReQ3n*dImQ1n-dReQ1n*dReQ3n*dImQ2n+dReQ1n*dReQ2n*dImQ3n-dImQ1n*dImQ2n*dImQ3n);
 
 
 // Re[Q_{2n} Q_{2n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)*dReQ2n-2.*dReQ1n*dReQ2n*dImQ1n-dReQ2n*pow(dImQ1n,2.)
                                               + dImQ2n*pow(dReQ1n,2.)+2.*dReQ1n*dImQ1n*dImQ2n-pow(dImQ1n,2.)*dImQ2n)
                                               * (pow(dReQ1n,2.)*dReQ2n+2.*dReQ1n*dReQ2n*dImQ1n-dReQ2n*pow(dImQ1n,2.)
                                               - dImQ2n*pow(dReQ1n,2.)+2.*dReQ1n*dImQ1n*dImQ2n+pow(dImQ1n,2.)*dImQ2n);
 
 // Im[Q_{2n} Q_{2n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 //Double_t imQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar = 2.*(pow(dReQ1n,2.)*dReQ2n-dReQ2n*pow(dImQ1n,2.)
 //                                              + 2.*dReQ1n*dImQ1n*dImQ2n)*(pow(dReQ1n,2.)*dImQ2n
 //                                              - 2.*dReQ1n*dImQ1n*dReQ2n-pow(dImQ1n,2.)*dImQ2n);
 
 // Re[Q_{3n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)+pow(dImQ1n,2.))
                                               * (pow(dReQ1n,3.)*dReQ3n-3.*dReQ1n*dReQ3n*pow(dImQ1n,2.)
                                               + 3.*pow(dReQ1n,2.)*dImQ1n*dImQ3n-pow(dImQ1n,3.)*dImQ3n);
  
 // Im[Q_{3n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]                                                                                           
 //Double_t imQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = (pow(dReQ1n,2.)+pow(dImQ1n,2.))
 //                                              * (pow(dImQ1n,3.)*dReQ3n-3.*dImQ1n*dReQ3n*pow(dReQ1n,2.)
 //                                              - 3.*pow(dImQ1n,2.)*dReQ1n*dImQ3n+pow(dReQ1n,3.)*dImQ3n);
 
 // |Q_{2n}|^2 |Q_{n}|^4
 Double_t dQ2nQ1nQ1nQ2nstarQ1nstarQ1nstar = (pow(dReQ2n,2.)+pow(dImQ2n,2.))*pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.);
 
 // Re[Q_{2n} Q_{n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]
 Double_t reQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                                                  * (pow(dReQ1n,2.)*dReQ2n-dReQ2n*pow(dImQ1n,2.)
                                                  + 2.*dReQ1n*dImQ1n*dImQ2n);
                                                  
 // Im[Q_{2n} Q_{n} Q_{n} Q_{n}^* Q_{n}^* Q_{n}^* Q_{n}^*]                                                  
 //Double_t imQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar = pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
 //                                                 * (pow(dReQ1n,2.)*dImQ2n-dImQ2n*pow(dImQ1n,2.)
 //                                                 - 2.*dReQ1n*dReQ2n*dImQ1n);
 
  
 
       
 //                                        **************************************
 //                                        **** multi-particle correlations: ****
 //                                        **************************************
 //
 // Remark 1: multi-particle correlations calculated with non-weighted Q-vectors are stored in 1D profile fQCorrelations[0].
 // Remark 2: binning of fQCorrelations[0] is organized as follows:
 // --------------------------------------------------------------------------------------------------------------------
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
 // --------------------------------------------------------------------------------------------------------------------
    
 // 2-particle:
 Double_t two1n1n = 0.; // <cos(n*(phi1-phi2))>
 Double_t two2n2n = 0.; // <cos(2n*(phi1-phi2))>
 Double_t two3n3n = 0.; // <cos(3n*(phi1-phi2))>
 Double_t two4n4n = 0.; // <cos(4n*(phi1-phi2))>
 
 if(dMult>1)
 {
  two1n1n = (pow(dReQ1n,2.)+pow(dImQ1n,2.)-dMult)/(dMult*(dMult-1.)); 
  two2n2n = (pow(dReQ2n,2.)+pow(dImQ2n,2.)-dMult)/(dMult*(dMult-1.)); 
  two3n3n = (pow(dReQ3n,2.)+pow(dImQ3n,2.)-dMult)/(dMult*(dMult-1.)); 
  two4n4n = (pow(dReQ4n,2.)+pow(dImQ4n,2.)-dMult)/(dMult*(dMult-1.)); 
  
  // average non-weighted 2-particle correlations for single event: 
  fQCorrelationsEBE[0]->SetBinContent(1,two1n1n);
  fQCorrelationsEBE[0]->SetBinContent(2,two2n2n);
  fQCorrelationsEBE[0]->SetBinContent(3,two3n3n);
  fQCorrelationsEBE[0]->SetBinContent(4,two4n4n);
        
  // final average non-weighted 2-particle correlations for all events:      
  fQCorrelations[0][0]->Fill(0.5,two1n1n,dMult*(dMult-1.));  
  fQCorrelations[0][0]->Fill(1.5,two2n2n,dMult*(dMult-1.)); 
  fQCorrelations[0][0]->Fill(2.5,two3n3n,dMult*(dMult-1.)); 
  fQCorrelations[0][0]->Fill(3.5,two4n4n,dMult*(dMult-1.)); 
  
  // distribution of <cos(n*(phi1-phi2))>:
  //f2pDistribution->Fill(two1n1n,dMult*(dMult-1.)); 
 } // end of if(dMult>1)
 
 // 3-particle:
 Double_t three2n1n1n = 0.; // <cos(n*(2.*phi1-phi2-phi3))>
 Double_t three3n2n1n = 0.; // <cos(n*(3.*phi1-2.*phi2-phi3))>
 Double_t three4n2n2n = 0.; // <cos(n*(4.*phi1-2.*phi2-2.*phi3))>
 Double_t three4n3n1n = 0.; // <cos(n*(4.*phi1-3.*phi2-phi3))>
 
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
              
  // average non-weighted 3-particle correlations for single event: 
  fQCorrelationsEBE[0]->SetBinContent(6,three2n1n1n);
  fQCorrelationsEBE[0]->SetBinContent(7,three3n2n1n);
  fQCorrelationsEBE[0]->SetBinContent(8,three4n2n2n);
  fQCorrelationsEBE[0]->SetBinContent(9,three4n3n1n);
        
  // final average non-weighted 3-particle correlations for all events:                
  fQCorrelations[0][0]->Fill(5.5,three2n1n1n,dMult*(dMult-1.)*(dMult-2.)); 
  fQCorrelations[0][0]->Fill(6.5,three3n2n1n,dMult*(dMult-1.)*(dMult-2.));
  fQCorrelations[0][0]->Fill(7.5,three4n2n2n,dMult*(dMult-1.)*(dMult-2.)); 
  fQCorrelations[0][0]->Fill(8.5,three4n3n1n,dMult*(dMult-1.)*(dMult-2.));    
 } // end of if(dMult>2)
 
 // 4-particle:
 Double_t four1n1n1n1n = 0.; // <cos(n*(phi1+phi2-phi3-phi4))>
 Double_t four2n2n2n2n = 0.; // <cos(2.*n*(phi1+phi2-phi3-phi4))>
 Double_t four2n1n2n1n = 0.; // <cos(n*(2.*phi1+phi2-2.*phi3-phi4))> 
 Double_t four3n1n1n1n = 0.; // <cos(n*(3.*phi1-phi2-phi3-phi4))> 
 Double_t four4n2n1n1n = 0.; // <cos(n*(4.*phi1-2.*phi2-phi3-phi4))> 
 Double_t four3n1n2n2n = 0.; // <cos(n*(3.*phi1+phi2-2.*phi3-2.*phi4))> 
 Double_t four3n1n3n1n = 0.; // <cos(n*(3.*phi1+phi2-3.*phi3-phi4))>   
 
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
  four3n1n1n1n = (reQ3nQ1nstarQ1nstarQ1nstar-3.*reQ3nQ2nstarQ1nstar-3.*reQ2nQ1nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               + (2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+3.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
               + 6.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-6.*dMult)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  four4n2n1n1n = (reQ4nQ2nstarQ1nstarQ1nstar-2.*reQ4nQ3nstarQ1nstar-reQ4nQ2nstarQ2nstar-2.*reQ3nQ2nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - (reQ2nQ1nstarQ1nstar-2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))-2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
               - 3.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - 6./((dMult-1.)*(dMult-2.)*(dMult-3.));
  four3n1n2n2n = (reQ3nQ1nQ2nstarQ2nstar-reQ4nQ2nstarQ2nstar-reQ3nQ1nQ4nstar-2.*reQ3nQ2nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - (2.*reQ1nQ1nQ2nstar-(pow(dReQ4n,2.)+pow(dImQ4n,2.))-2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
               - 4.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))-4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               - 6./((dMult-1.)*(dMult-2.)*(dMult-3.)); 
  four3n1n3n1n = ((pow(dReQ3n,2.)+pow(dImQ3n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
               - 2.*reQ4nQ3nstarQ1nstar-2.*reQ3nQ2nstarQ1nstar)
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               + ((pow(dReQ4n,2.)+pow(dImQ4n,2.))-(dMult-4.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
               + (pow(dReQ2n,2.)+pow(dImQ2n,2.))-(dMult-4.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
               / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.))
               + (dMult-6.)/((dMult-1.)*(dMult-2.)*(dMult-3.));
               
  // average non-weighted 4-particle correlations for single event: 
  fQCorrelationsEBE[0]->SetBinContent(11,four1n1n1n1n);
  fQCorrelationsEBE[0]->SetBinContent(12,four2n1n2n1n);
  fQCorrelationsEBE[0]->SetBinContent(13,four2n2n2n2n);
  fQCorrelationsEBE[0]->SetBinContent(14,four3n1n1n1n);
  fQCorrelationsEBE[0]->SetBinContent(15,four3n1n3n1n);
  fQCorrelationsEBE[0]->SetBinContent(16,four3n1n2n2n);
  fQCorrelationsEBE[0]->SetBinContent(17,four4n2n1n1n);
        
  // final average non-weighted 4-particle correlations for all events:                
  fQCorrelations[0][0]->Fill(10.5,four1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations[0][0]->Fill(11.5,four2n1n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations[0][0]->Fill(12.5,four2n2n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations[0][0]->Fill(13.5,four3n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations[0][0]->Fill(14.5,four3n1n3n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fQCorrelations[0][0]->Fill(15.5,four3n1n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));  
  fQCorrelations[0][0]->Fill(16.5,four4n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)); 
  
  // distribution of <cos(n*(phi1+phi2-phi3-phi4))>
  //f4pDistribution->Fill(four1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  
 } // end of if(dMult>3)

 // 5-particle:
 Double_t five2n1n1n1n1n = 0.; // <cos(n*(2.*phi1+phi2-phi3-phi4-phi5))>
 Double_t five2n2n2n1n1n = 0.; // <cos(n*(2.*phi1+2.*phi2-2.*phi3-phi4-phi5))>
 Double_t five3n1n2n1n1n = 0.; // <cos(n*(3.*phi1+phi2-2.*phi3-phi4-phi5))>
 Double_t five4n1n1n1n1n = 0.; // <cos(n*(4.*phi1-phi2-phi3-phi4-phi5))>
 
 if(dMult>4)
 {
  five2n1n1n1n1n = (reQ2nQ1nQ1nstarQ1nstarQ1nstar-reQ3nQ1nstarQ1nstarQ1nstar+6.*reQ3nQ2nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (reQ2nQ1nQ3nstar+3.*(dMult-6.)*reQ2nQ1nstarQ1nstar+3.*reQ1nQ1nQ2nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 + 3.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))     
                 - 3.*(dMult-4.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.)))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - 3.*(pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                 - 2.*(2*dMult-5.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))+2.*dMult*(dMult-4.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
                 
  five2n2n2n1n1n = (reQ2nQ2nQ2nstarQ1nstarQ1nstar-reQ4nQ2nstarQ1nstarQ1nstar-2.*reQ2nQ2nQ3nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 + 2.*(reQ4nQ2nstarQ2nstar+4.*reQ3nQ2nstarQ1nstar+reQ3nQ1nQ4nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 + (reQ2nQ2nQ4nstar-2.*(dMult-5.)*reQ2nQ1nstarQ1nstar+2.*reQ1nQ1nQ2nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+4.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 + 1.*pow((pow(dReQ2n,2.)+pow(dImQ2n,2.)),2.)
                 - 2.*(3.*dMult-10.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.)))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (4.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 - 4.*(dMult-5.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.))+4.*dMult*(dMult-6.))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)); 

  five4n1n1n1n1n = (reQ4nQ1nstarQ1nstarQ1nstarQ1nstar-6.*reQ4nQ2nstarQ1nstarQ1nstar-4.*reQ3nQ1nstarQ1nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 + (8.*reQ4nQ3nstarQ1nstar+3.*reQ4nQ2nstarQ2nstar+12.*reQ3nQ2nstarQ1nstar+12.*reQ2nQ1nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (6.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))+8.*(pow(dReQ3n,2.)+pow(dImQ3n,2.))
                 + 12.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))+24.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-24.*dMult)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  
  five3n1n2n1n1n = (reQ3nQ1nQ2nstarQ1nstarQ1nstar-reQ4nQ2nstarQ1nstarQ1nstar-reQ3nQ1nstarQ1nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (reQ3nQ1nQ2nstarQ2nstar-3.*reQ4nQ3nstarQ1nstar-reQ4nQ2nstarQ2nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - ((2.*dMult-13.)*reQ3nQ2nstarQ1nstar-reQ3nQ1nQ4nstar-9.*reQ2nQ1nstarQ1nstar)
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - (2.*reQ1nQ1nQ2nstar+2.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
                 - 2.*(dMult-5.)*(pow(dReQ3n,2.)+pow(dImQ3n,2.))+2.*(pow(dReQ3n,2.)
                 + pow(dImQ3n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 + (2.*(dMult-6.)*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 - 2.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                 - pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)
                 + 2.*(3.*dMult-11.)*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
                 / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.))
                 - 4.*(dMult-6.)/((dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
                 
  // average non-weighted 5-particle correlations for single event: 
  fQCorrelationsEBE[0]->SetBinContent(19,five2n1n1n1n1n);
  fQCorrelationsEBE[0]->SetBinContent(20,five2n2n2n1n1n);
  fQCorrelationsEBE[0]->SetBinContent(21,five3n1n2n1n1n);
  fQCorrelationsEBE[0]->SetBinContent(22,five4n1n1n1n1n);
        
  // final average non-weighted 5-particle correlations for all events:                         
  fQCorrelations[0][0]->Fill(18.5,five2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)); 
  fQCorrelations[0][0]->Fill(19.5,five2n2n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fQCorrelations[0][0]->Fill(20.5,five3n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fQCorrelations[0][0]->Fill(21.5,five4n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
 } // end of if(dMult>4)
    
 // 6-particle:
 Double_t six1n1n1n1n1n1n = 0.; // <cos(n*(phi1+phi2+phi3-phi4-phi5-phi6))>
 Double_t six2n2n1n1n1n1n = 0.; // <cos(n*(2.*phi1+2.*phi2-phi3-phi4-phi5-phi6))>
 Double_t six3n1n1n1n1n1n = 0.; // <cos(n*(3.*phi1+phi2-phi3-phi4-phi5-phi6))>
 Double_t six2n1n1n2n1n1n = 0.; // <cos(n*(2.*phi1+phi2+phi3-2.*phi4-phi5-phi6))>
 
 if(dMult>5)
 {
  six1n1n1n1n1n1n = (pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),3.)+9.*dQ2nQ1nQ2nstarQ1nstar-6.*reQ2nQ1nQ1nstarQ1nstarQ1nstar)
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.))
                  + 4.*(reQ3nQ1nstarQ1nstarQ1nstar-3.*reQ3nQ2nstarQ1nstar)
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.))
                  + 2.*(9.*(dMult-4.)*reQ2nQ1nstarQ1nstar+2.*(pow(dReQ3n,2.)+pow(dImQ3n,2.)))
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.))
                  - 9.*(pow((pow(dReQ1n,2.)+pow(dImQ1n,2.)),2.)+(pow(dReQ2n,2.)+pow(dImQ2n,2.)))
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-5.))
                  + (18.*(pow(dReQ1n,2.)+pow(dImQ1n,2.)))
                  / (dMult*(dMult-1)*(dMult-3)*(dMult-4))
                  - 6./((dMult-1.)*(dMult-2.)*(dMult-3.));
                  
  six2n1n1n2n1n1n = (dQ2nQ1nQ1nQ2nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)
                  * (2.*five2n2n2n1n1n+4.*five2n1n1n1n1n+4.*five3n1n2n1n1n+4.*four2n1n2n1n+1.*four1n1n1n1n)
                  - dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(4.*four1n1n1n1n+4.*two1n1n
                  + 2.*three2n1n1n+2.*three2n1n1n+4.*four3n1n1n1n+8.*three2n1n1n+2.*four4n2n1n1n
                  + 4.*four2n1n2n1n+2.*two2n2n+8.*four2n1n2n1n+4.*four3n1n3n1n+8.*three3n2n1n
                  + 4.*four3n1n2n2n+4.*four1n1n1n1n+4.*four2n1n2n1n+1.*four2n2n2n2n)
                  - dMult*(dMult-1.)*(dMult-2.)*(2.*three2n1n1n+8.*two1n1n+4.*two1n1n+2.
                  + 4.*two1n1n+4.*three2n1n1n+2.*two2n2n+4.*three2n1n1n+8.*three3n2n1n
                  + 8.*two2n2n+4.*three4n3n1n+4.*two3n3n+4.*three3n2n1n+4.*two1n1n
                  + 8.*three2n1n1n+4.*two1n1n+4.*three3n2n1n+4.*three2n1n1n+2.*two2n2n
                  + 4.*three3n2n1n+2.*three4n2n2n)-dMult*(dMult-1.)
                  * (4.*two1n1n+4.+4.*two1n1n+2.*two2n2n+1.+4.*two1n1n+4.*two2n2n+4.*two3n3n
                  + 1.+2.*two2n2n+1.*two4n4n)-dMult)
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); // to be improved (direct formula needed)
 
  six2n2n1n1n1n1n = (reQ2nQ2nQ1nstarQ1nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)
                  * (five4n1n1n1n1n+8.*five2n1n1n1n1n+6.*five2n2n2n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)
                  * (4.*four3n1n1n1n+6.*four4n2n1n1n+12.*three2n1n1n+12.*four1n1n1n1n+24.*four2n1n2n1n
                  + 4.*four3n1n2n2n+3.*four2n2n2n2n)-dMult*(dMult-1.)*(dMult-2.)*(6.*three2n1n1n+12.*three3n2n1n
                  + 4.*three4n3n1n+3.*three4n2n2n+8.*three2n1n1n+24.*two1n1n+12.*two2n2n+12.*three2n1n1n+8.*three3n2n1n
                  + 1.*three4n2n2n)-dMult*(dMult-1.)*(4.*two1n1n+6.*two2n2n+4.*two3n3n+1.*two4n4n+2.*two2n2n+8.*two1n1n+6.)-dMult)
                  / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); // to be improved (direct formula needed)
   
  six3n1n1n1n1n1n = (reQ3nQ1nQ1nstarQ1nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)
                  * (five4n1n1n1n1n+4.*five2n1n1n1n1n+6.*five3n1n2n1n1n+4.*four3n1n1n1n)
                  - dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(4.*four3n1n1n1n+6.*four4n2n1n1n+6.*four1n1n1n1n
                  + 12.*three2n1n1n+12.*four2n1n2n1n+6.*four3n1n1n1n+12.*three3n2n1n+4.*four3n1n3n1n+3.*four3n1n2n2n)
                  - dMult*(dMult-1.)*(dMult-2.)*(6.*three2n1n1n+12.*three3n2n1n+4.*three4n3n1n+3.*three4n2n2n+4.*two1n1n
                  + 12.*two1n1n+6.*three2n1n1n+12.*three2n1n1n+4.*three3n2n1n+12.*two2n2n+4.*three3n2n1n+4.*two3n3n+1.*three4n3n1n
                  + 6.*three3n2n1n)-dMult*(dMult-1.)*(4.*two1n1n+6.*two2n2n+4.*two3n3n+1.*two4n4n+1.*two1n1n+4.+6.*two1n1n+4.*two2n2n
                  + 1.*two3n3n)-dMult)/(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); // to be improved (direct formula needed)
                                 
  // average non-weighted 6-particle correlations for single event: 
  fQCorrelationsEBE[0]->SetBinContent(24,six1n1n1n1n1n1n);
  fQCorrelationsEBE[0]->SetBinContent(25,six2n1n1n2n1n1n);
  fQCorrelationsEBE[0]->SetBinContent(26,six2n2n1n1n1n1n);
  fQCorrelationsEBE[0]->SetBinContent(27,six3n1n1n1n1n1n);
        
  // final average non-weighted 6-particle correlations for all events:         
  fQCorrelations[0][0]->Fill(23.5,six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  fQCorrelations[0][0]->Fill(24.5,six2n1n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  fQCorrelations[0][0]->Fill(25.5,six2n2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  fQCorrelations[0][0]->Fill(26.5,six3n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 

  // distribution of <cos(n*(phi1+phi2+phi3-phi4-phi5-phi6))>
  //f6pDistribution->Fill(six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
 } // end of if(dMult>5)
 
 // 7-particle:
 Double_t seven2n1n1n1n1n1n1n = 0.; // <cos(n*(2.*phi1+phi2+phi3-phi4-phi5-phi6-phi7))>
 
 if(dMult>6)
 {
  seven2n1n1n1n1n1n1n = (reQ2nQ1nQ1nQ1nstarQ1nstarQ1nstarQ1nstar-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)
                      * (2.*six3n1n1n1n1n1n+4.*six1n1n1n1n1n1n+1.*six2n2n1n1n1n1n+6.*six2n1n1n2n1n1n+8.*five2n1n1n1n1n)
                      - dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(1.*five4n1n1n1n1n +8.*five2n1n1n1n1n+8.*four3n1n1n1n
                      + 12.*five3n1n2n1n1n+4.*five2n1n1n1n1n+3.*five2n2n2n1n1n+6.*five2n2n2n1n1n+6.*four1n1n1n1n+24.*four1n1n1n1n
                      + 12.*five2n1n1n1n1n+12.*five2n1n1n1n1n+12.*three2n1n1n+24.*four2n1n2n1n+4.*five3n1n2n1n1n+4.*five2n1n1n1n1n)
                      - dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(4.*four3n1n1n1n+6.*four4n2n1n1n+12.*four1n1n1n1n+24.*three2n1n1n
                      + 24.*four2n1n2n1n+12.*four3n1n1n1n+24.*three3n2n1n+8.*four3n1n3n1n+6.*four3n1n2n2n+6.*three2n1n1n+12.*four1n1n1n1n
                      + 12.*four2n1n2n1n+6.*three2n1n1n+12.*four2n1n2n1n+4.*four3n1n2n2n+3.*four2n2n2n2n+4.*four1n1n1n1n+6.*three2n1n1n
                      + 24.*two1n1n+24.*four1n1n1n1n+4.*four3n1n1n1n+24.*two1n1n+24.*three2n1n1n+12.*two2n2n+24.*three2n1n1n+12.*four2n1n2n1n
                      + 8.*three3n2n1n+8.*four2n1n2n1n+1.*four4n2n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(6.*three2n1n1n+1.*three2n1n1n+8.*two1n1n
                      + 12.*three3n2n1n+24.*two1n1n+12.*three2n1n1n+4.*three2n1n1n+8.*two1n1n+4.*three4n3n1n+24.*three2n1n1n+8.*three3n2n1n
                      + 12.*two1n1n+12.*two1n1n+3.*three4n2n2n+24.*two2n2n+6.*two2n2n+12.+12.*three3n2n1n+8.*two3n3n+12.*three2n1n1n+24.*two1n1n
                      + 4.*three3n2n1n+8.*three3n2n1n+2.*three4n3n1n+12.*two1n1n+8.*three2n1n1n+4.*three2n1n1n+2.*three3n2n1n+6.*two2n2n+8.*two2n2n
                      + 1.*three4n2n2n+4.*three3n2n1n+6.*three2n1n1n)-dMult*(dMult-1.)*(4.*two1n1n+2.*two1n1n+6.*two2n2n+8.+1.*two2n2n+4.*two3n3n
                      + 12.*two1n1n+4.*two1n1n+1.*two4n4n+8.*two2n2n+6.+2.*two3n3n+4.*two1n1n+1.*two2n2n)-dMult)
                      / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)); // to be improved (direct formula needed)
        
  // average non-weighted 7-particle correlations for single event: 
  fQCorrelationsEBE[0]->SetBinContent(29,seven2n1n1n1n1n1n1n);
       
  // final average non-weighted 7-particle correlations for all events:                      
  fQCorrelations[0][0]->Fill(28.5,seven2n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.));
 } // end of if(dMult>6)
 
 // 8-particle:
 Double_t eight1n1n1n1n1n1n1n1n = 0.; // <cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8))>
 if(dMult>7)
 {
  eight1n1n1n1n1n1n1n1n = (pow(pow(dReQ1n,2.)+pow(dImQ1n,2.),4.)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)
                        * (12.*seven2n1n1n1n1n1n1n+16.*six1n1n1n1n1n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)
                        * (8.*six3n1n1n1n1n1n+48.*six1n1n1n1n1n1n+6.*six2n2n1n1n1n1n+96.*five2n1n1n1n1n+72.*four1n1n1n1n+36.*six2n1n1n2n1n1n)
                        - dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(2.*five4n1n1n1n1n+32.*five2n1n1n1n1n+36.*four1n1n1n1n
                        + 32.*four3n1n1n1n+48.*five2n1n1n1n1n+48.*five3n1n2n1n1n+144.*five2n1n1n1n1n+288.*four1n1n1n1n+36.*five2n2n2n1n1n
                        + 144.*three2n1n1n+96.*two1n1n+144.*four2n1n2n1n)-dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)
                        * (8.*four3n1n1n1n+48.*four1n1n1n1n+12.*four4n2n1n1n+96.*four2n1n2n1n+96.*three2n1n1n+72.*three2n1n1n+144.*two1n1n
                        + 16.*four3n1n3n1n+48.*four3n1n1n1n+144.*four1n1n1n1n+72.*four1n1n1n1n+96.*three3n2n1n+24.*four3n1n2n2n+144.*four2n1n2n1n
                        + 288.*two1n1n+288.*three2n1n1n+9.*four2n2n2n2n+72.*two2n2n+24.)-dMult*(dMult-1.)*(dMult-2.)*(12.*three2n1n1n+16.*two1n1n
                        + 24.*three3n2n1n+48.*three2n1n1n+96.*two1n1n+8.*three4n3n1n+32.*three3n2n1n+96.*three2n1n1n+144.*two1n1n+6.*three4n2n2n
                        + 96.*two2n2n+36.*two2n2n+72.+48.*three3n2n1n+16.*two3n3n+72.*three2n1n1n+144.*two1n1n)-dMult*(dMult-1.)*(8.*two1n1n
                        + 12.*two2n2n+16.+8.*two3n3n+48.*two1n1n+1.*two4n4n+16.*two2n2n+18.)-dMult)
                        / (dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.)); // to be improved (direct formula needed)
  
  // average non-weighted 8-particle correlations for single event: 
  fQCorrelationsEBE[0]->SetBinContent(31,eight1n1n1n1n1n1n1n1n);
       
  // final average non-weighted 8-particle correlations for all events:                       
  fQCorrelations[0][0]->Fill(30.5,eight1n1n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
 
  // distribution of <cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8))>
  //f8pDistribution->Fill(eight1n1n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
 } // end of if(dMult>7) 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateCorrelationsForIntegratedFlow()











































//================================================================================================================================


void AliFlowAnalysisWithQCumulants::EvaluateNestedLoopsForIntegratedFlow(AliFlowEventSimple* anEvent)
{
 // 1.) Evaluate with nested loops the relevant correlations for integrated flow without and with using the particle weights. 
 //     Results are stored in profiles fDirectCorrelations and fDirectCorrelationsW, respectively.
 
 // 2.) Evaluate with nested loops corrections for non-uniform acceptance relevant for integrated flow, 
 //     without and with using the particle weights.
 //     Without weights: cos terms are stored in profile fDirectCorrectionsCos, and sin terms in profile fDirectCorrectionsSin.
 //     With weights: cos terms are stored in profile fDirectCorrectionsCosW, and sin terms in profile fDirectCorrectionsSinW.
 
 // 3.) Binning of fDirectCorrelations is organized as follows:
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
 
 // 4.) Binning of fDirectCorrelationsW is organized as follows:
 // ..............................................................................................
 //       ---- bins 1-20: 2-particle correlations ----
 // 1st bin: two1n1nW1W1 = <w1 w2 cos(n*(phi1-phi2))>
 // 2nd bin: two2n2nW2W2 = <w1^2 w2^2 cos(2n*(phi1-phi2))>
 // 3rd bin: two3n3nW3W3 = <w1^3 w2^3 cos(3n*(phi1-phi2))>
 // 4th bin: two4n4nW4W4 = <w1^4 w2^4 cos(4n*(phi1-phi2))>
 // 5th bin: two1n1nW3W1 = <w1^3 w2 cos(n*(phi1-phi2))>
 // 6th bin: two1n1nW1W1W2 = <w1 w2 w3^2 cos(n*(phi1-phi2))>  
 //       ---- bins 21-40: 3-particle correlations ----
 // 21st bin: three2n1n1nW2W1W1 = <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))> 
 //       ---- bins 41-60: 4-particle correlations ----
 // 41st bin: four1n1n1n1nW1W1W1W1 = <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
 //       ---- bins 61-80: 5-particle correlations ---- 
 //       ---- bins 81-100: 6-particle correlations ----
 //       ---- bins 101-120: 7-particle correlations ----
 //       ---- bins 121-140: 8-particle correlations ----
 // ..............................................................................................
 
 // 5.) Binning of fDirectCorrectionsCos is organized as follows:
 // ..............................................................................................
 // 1st bin: <<cos(n*(phi1))>> = cosP1n
 // 2nd bin: <<cos(n*(phi1+phi2))>> = cosP1nP1n
 // 3rd bin: <<cos(n*(phi1-phi2-phi3))>> = cosP1nM1nM1n
 // ...
 // ..............................................................................................
              
 // 6.) Binning of fDirectCorrectionsSin is organized as follows:
 // ..............................................................................................
 // 1st bin: <<sin(n*(phi1))>> = sinP1n
 // 2nd bin: <<sin(n*(phi1+phi2))>> = sinP1nP1n
 // 3rd bin: <<sin(n*(phi1-phi2-phi3))>> = sinP1nM1nM1n
 // ...
 // ..............................................................................................
       
 // 7.) Binning of fDirectCorrectionsCosW is organized as follows:
 // ..............................................................................................     
 // ...
 // ..............................................................................................     
 
 // 8.) Binning of fDirectCorrectionsSinW is organized as follows:
 // ..............................................................................................     
 // ...
 // ..............................................................................................     
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 
 Double_t phi1=0., phi2=0., phi3=0., phi4=0., phi5=0., phi6=0., phi7=0., phi8=0.;
 Double_t wPhi1=1., wPhi2=1., wPhi3=1., wPhi4=1., wPhi5=1., wPhi6=1., wPhi7=1., wPhi8=1.;
 
 Int_t n = fHarmonic; 
 
 // 2-particle correlations and 1- and 2-particle correction terms:       
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  if(!(aftsTrack->InRPSelection())) continue;
  phi1=aftsTrack->Phi();
  if(fUsePhiWeights && fPhiWeights) wPhi1 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi1*fnBinsPhi/TMath::TwoPi())));
  
  // corrections for non-uniform acceptance:
  // non-weighted: 
  fDirectCorrectionsCos->Fill(0.5,cos(n*phi1),1.); // <cos(n*phi1)>
  fDirectCorrectionsSin->Fill(0.5,sin(n*phi1),1.); // <sin(n*phi1)>  
  // weighted:
  // fDirectCorrectionsCosW->Fill(0.5,???,1); // to be improved (continued)
  // fDirectCorrectionsSinW->Fill(0.5,???,1); // to be improved (continued)
  
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1)continue;
   aftsTrack=anEvent->GetTrack(i2);
   if(!(aftsTrack->InRPSelection())) continue;
   phi2=aftsTrack->Phi();
   if(fUsePhiWeights && fPhiWeights) wPhi2 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi2*fnBinsPhi/TMath::TwoPi())));
    
   // non-weighted correlations: 
   fDirectCorrelations->Fill(0.5,cos(n*(phi1-phi2)),1.);    // <cos(n*(phi1-phi2))>
   fDirectCorrelations->Fill(1.5,cos(2.*n*(phi1-phi2)),1.); // <cos(2n*(phi1-phi2))>
   fDirectCorrelations->Fill(2.5,cos(3.*n*(phi1-phi2)),1.); // <cos(3n*(phi1-phi2))>
   fDirectCorrelations->Fill(3.5,cos(4.*n*(phi1-phi2)),1.); // <cos(4n*(phi1-phi2))> 
   
   // weighted correlations:
   // ................................................................................................................
   if(fUsePhiWeights) fDirectCorrelationsW->Fill(0.5,cos(n*(phi1-phi2)),wPhi1*wPhi2);                  // <w1   w2   cos( n*(phi1-phi2))>
   if(fUsePhiWeights) fDirectCorrelationsW->Fill(1.5,cos(2.*n*(phi1-phi2)),pow(wPhi1,2)*pow(wPhi2,2)); // <w1^2 w2^2 cos(2n*(phi1-phi2))>
   if(fUsePhiWeights) fDirectCorrelationsW->Fill(2.5,cos(3.*n*(phi1-phi2)),pow(wPhi1,3)*pow(wPhi2,3)); // <w1^3 w2^3 cos(3n*(phi1-phi2))>
   if(fUsePhiWeights) fDirectCorrelationsW->Fill(3.5,cos(4.*n*(phi1-phi2)),pow(wPhi1,4)*pow(wPhi2,4)); // <w1^4 w2^4 cos(4n*(phi1-phi2))> 
   if(fUsePhiWeights) fDirectCorrelationsW->Fill(4.5,cos(n*(phi1-phi2)),pow(wPhi1,3)*wPhi2);           // <w1^3 w2 cos(n*(phi1-phi2))>
   // ...
   // ................................................................................................................
 
   // non-weighted corrections for non-uniform acceptance (cos terms)
   // ................................................................................................................
   fDirectCorrectionsCos->Fill(1.5,cos(n*(phi1+phi2)),1.); // <<cos(n*(phi1+phi2))>>
   // ...
   // ................................................................................................................
  
   // non-weighted corrections for non-uniform acceptance (sin terms)
   // ................................................................................................................
   fDirectCorrectionsSin->Fill(1.5,sin(n*(phi1+phi2)),1.); // <<sin(n*(phi1+phi2))>>
   // ...
   // ................................................................................................................
   
   // weighted corrections for non-uniform acceptance (cos terms)
   // ................................................................................................................
   // fDirectCorrectionsCosW->Fill(1.5,???,1.); // to be improved (continued)
   // ...
   // ................................................................................................................
  
   // non-weighted corrections for non-uniform acceptance (sin terms)
   // ................................................................................................................
   // fDirectCorrectionsSinW->Fill(1.5,???,1.); // to be improved (continued)
   // ...
   // ................................................................................................................

  } // end of for(Int_t i2=0;i2<nPrim;i2++)
 } // end of for(Int_t i1=0;i1<nPrim;i1++)
 
 // 3-particle correlations:         
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
    
    // non-weighted correlations:
    fDirectCorrelations->Fill(5.,cos(2.*n*phi1-n*(phi2+phi3)),1.);       //<3>_{2n|nn,n}
    fDirectCorrelations->Fill(6.,cos(3.*n*phi1-2.*n*phi2-n*phi3),1.);    //<3>_{3n|2n,n}
    fDirectCorrelations->Fill(7.,cos(4.*n*phi1-2.*n*phi2-2.*n*phi3),1.); //<3>_{4n|2n,2n}
    fDirectCorrelations->Fill(8.,cos(4.*n*phi1-3.*n*phi2-n*phi3),1.);    //<3>_{4n|3n,n}
    
    // weighted correlations:
    // ..............................................................................................................................
    // 2-p:
    if(fUsePhiWeights) fDirectCorrelationsW->Fill(5.,cos(n*(phi1-phi2)),wPhi1*wPhi2*pow(wPhi3,2)); // <w1 w2 w3^2 cos(n*(phi1-phi2))>
    // 3-p:
    if(fUsePhiWeights) fDirectCorrelationsW->Fill(20.,cos(2.*n*phi1-n*(phi2+phi3)),pow(wPhi1,2)*wPhi2*wPhi3); // <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))>
    // ...
    // ..............................................................................................................................
    
    // non-weighted corrections for non-uniform acceptance (cos terms)
    // ................................................................................................................
    fDirectCorrectionsCos->Fill(2.,cos(n*(phi1-phi2-phi3)),1.); // <<cos(n*(phi1-phi2-phi3))>>
    // ...
    // ................................................................................................................
  
    // non-weighted corrections for non-uniform acceptance (sin terms)
    // ................................................................................................................
    fDirectCorrectionsSin->Fill(2.,sin(n*(phi1-phi2-phi3)),1.); // <<sin(n*(phi1-phi2-phi3))>>
    // ...
    // ................................................................................................................
    
    // weighted corrections for non-uniform acceptance (cos terms)
    // ................................................................................................................
    // ...
    // ................................................................................................................
    
    // weighted corrections for non-uniform acceptance (sin terms)
    // ................................................................................................................
    // ...
    // ................................................................................................................

   } // end of for(Int_t i3=0;i3<nPrim;i3++)
  } // end of for(Int_t i2=0;i2<nPrim;i2++)
 } // end of for(Int_t i1=0;i1<nPrim;i1++)

 // 4-particle correlations:       
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
     
     // non-weighted:
     fDirectCorrelations->Fill(10.,cos(n*phi1+n*phi2-n*phi3-n*phi4),1.);            // <4>_{n,n|n,n} 
     fDirectCorrelations->Fill(11.,cos(2.*n*phi1+n*phi2-2.*n*phi3-n*phi4),1.);      // <4>_{2n,n|2n,n}
     fDirectCorrelations->Fill(12.,cos(2.*n*phi1+2*n*phi2-2.*n*phi3-2.*n*phi4),1.); // <4>_{2n,2n|2n,2n}
     fDirectCorrelations->Fill(13.,cos(3.*n*phi1-n*phi2-n*phi3-n*phi4),1.);         // <4>_{3n|n,n,n}
     fDirectCorrelations->Fill(14.,cos(3.*n*phi1+n*phi2-3.*n*phi3-n*phi4),1.);      // <4>_{3n,n|3n,n}   
     fDirectCorrelations->Fill(15.,cos(3.*n*phi1+n*phi2-2.*n*phi3-2.*n*phi4),1.);   // <4>_{3n,n|2n,2n}
     fDirectCorrelations->Fill(16.,cos(4.*n*phi1-2.*n*phi2-n*phi3-n*phi4),1.);      // <4>_{4n|2n,n,n}
     
     // weighted:
     //.......................................................................................
     // 4-p:
     if(fUsePhiWeights) fDirectCorrelationsW->Fill(40.,cos(n*phi1+n*phi2-n*phi3-n*phi4),wPhi1*wPhi2*wPhi3*wPhi4); 
     // ...             
     //.......................................................................................
     
    }  
   }
  }
 }

 // 5-particle correlations:      
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  //cout<<"i1 = "<<i1<<endl;
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
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      aftsTrack=anEvent->GetTrack(i5);
      if(!(aftsTrack->InRPSelection())) continue;
      phi5=aftsTrack->Phi();
      if(fUsePhiWeights && fPhiWeights) wPhi5 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*fnBinsPhi/TMath::TwoPi())));
      
      // non-weighted:
      //------------------------------------------------------------------------------------------------------
      fDirectCorrelations->Fill(18.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),1.);       //<5>_{2n,n|n,n,n}
      fDirectCorrelations->Fill(19.,cos(2.*n*phi1+2.*n*phi2-2.*n*phi3-n*phi4-n*phi5),1.); //<5>_{2n,2n|2n,n,n}
      fDirectCorrelations->Fill(20.,cos(3.*n*phi1+n*phi2-2.*n*phi3-n*phi4-n*phi5),1.);    //<5>_{3n,n|2n,n,n}
      fDirectCorrelations->Fill(21.,cos(4.*n*phi1-n*phi2-n*phi3-n*phi4-n*phi5),1.);       //<5>_{4n|n,n,n,n}
      //------------------------------------------------------------------------------------------------------
      
      // weighted:
      //..............................................................................................................
      // 5-p:
      if(fUsePhiWeights) fDirectCorrelationsW->Fill(60.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),pow(wPhi1,2)*wPhi2*wPhi3*wPhi4*wPhi5);     
      //..............................................................................................................
      
     }
    }  
   }
  }
 }
 
 // 6-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  //cout<<"i1 = "<<i1<<endl;
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
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      aftsTrack=anEvent->GetTrack(i5);
      if(!(aftsTrack->InRPSelection())) continue;
      phi5=aftsTrack->Phi();
      if(fUsePhiWeights && fPhiWeights) wPhi5 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*fnBinsPhi/TMath::TwoPi())));
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       aftsTrack=anEvent->GetTrack(i6);
       if(!(aftsTrack->InRPSelection())) continue;
       phi6=aftsTrack->Phi(); 
       if(fUsePhiWeights && fPhiWeights) wPhi6 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi6*fnBinsPhi/TMath::TwoPi())));
       
       // non-weighted:
       //-----------------------------------------------------------------------------------------------------------
       fDirectCorrelations->Fill(23.,cos(n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6),1.);       //<6>_{n,n,n|n,n,n}
       fDirectCorrelations->Fill(24.,cos(2.*n*phi1+n*phi2+n*phi3-2.*n*phi4-n*phi5-n*phi6),1.); //<6>_{2n,n,n|2n,n,n}
       fDirectCorrelations->Fill(25.,cos(2.*n*phi1+2.*n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),1.); //<6>_{2n,2n|n,n,n,n}
       fDirectCorrelations->Fill(26.,cos(3.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),1.);    //<6>_{3n,n|n,n,n,n}  
       //-----------------------------------------------------------------------------------------------------------

       // weighted:
       //.................................................................................................................
       // 6-p:
       if(fUsePhiWeights) fDirectCorrelationsW->Fill(80.,cos(n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6);
       //.................................................................................................................       
          
      } 
     }
    }  
   }
  }
 }
 
 // 7-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  //cout<<"i1 = "<<i1<<endl;
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
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      aftsTrack=anEvent->GetTrack(i5);
      if(!(aftsTrack->InRPSelection())) continue;
      phi5=aftsTrack->Phi();
      if(fUsePhiWeights && fPhiWeights) wPhi5 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*fnBinsPhi/TMath::TwoPi())));
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       aftsTrack=anEvent->GetTrack(i6);
       if(!(aftsTrack->InRPSelection())) continue;
       phi6=aftsTrack->Phi(); 
       if(fUsePhiWeights && fPhiWeights) wPhi6 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi6*fnBinsPhi/TMath::TwoPi())));
       for(Int_t i7=0;i7<nPrim;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        aftsTrack=anEvent->GetTrack(i7);
        if(!(aftsTrack->InRPSelection())) continue;
        phi7=aftsTrack->Phi(); 
        if(fUsePhiWeights && fPhiWeights) wPhi7 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi7*fnBinsPhi/TMath::TwoPi())));
        
        // non-weighted:
        //---------------------------------------------------------------------------------------------------------------
        fDirectCorrelations->Fill(28.,cos(2.*n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6-n*phi7),1.);//<7>_{2n,n,n|n,n,n,n}
        //---------------------------------------------------------------------------------------------------------------
        
        // weighted:
        //..........................................................................................................................................
        if(fUsePhiWeights) fDirectCorrelationsW->Fill(100.,cos(2.*n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6-n*phi7),
                                                           pow(wPhi1,2.)*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6*wPhi7);
        //..........................................................................................................................................
        
       } 
      } 
     }
    }  
   }
  }
 }
 
 cout<<endl;
 
 // 8-particle correlations:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  cout<<"i1 = "<<i1<<endl;
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
     for(Int_t i5=0;i5<nPrim;i5++)
     {
      if(i5==i1||i5==i2||i5==i3||i5==i4)continue;
      aftsTrack=anEvent->GetTrack(i5);
      if(!(aftsTrack->InRPSelection())) continue;
      phi5=aftsTrack->Phi();
      if(fUsePhiWeights && fPhiWeights) wPhi5 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi5*fnBinsPhi/TMath::TwoPi())));
      for(Int_t i6=0;i6<nPrim;i6++)
      {
       if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5)continue;
       aftsTrack=anEvent->GetTrack(i6);
       if(!(aftsTrack->InRPSelection())) continue;
       phi6=aftsTrack->Phi();
       if(fUsePhiWeights && fPhiWeights) wPhi6 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi6*fnBinsPhi/TMath::TwoPi()))); 
       for(Int_t i7=0;i7<nPrim;i7++)
       {
        if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6)continue;
        aftsTrack=anEvent->GetTrack(i7);
        if(!(aftsTrack->InRPSelection())) continue;
        phi7=aftsTrack->Phi();
        if(fUsePhiWeights && fPhiWeights) wPhi7 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi7*fnBinsPhi/TMath::TwoPi()))); 
        for(Int_t i8=0;i8<nPrim;i8++)
        {
         if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7)continue;
         aftsTrack=anEvent->GetTrack(i8);
         if(!(aftsTrack->InRPSelection())) continue;
         phi8=aftsTrack->Phi();
         if(fUsePhiWeights && fPhiWeights) wPhi8 = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(phi8*fnBinsPhi/TMath::TwoPi()))); 
          
         // non-weighted: 
         //--------------------------------------------------------------------------------------------------------------------
         fDirectCorrelations->Fill(30.,cos(n*phi1+n*phi2+n*phi3+n*phi4-n*phi5-n*phi6-n*phi7-n*phi8),1.);//<8>_{n,n,n,n|n,n,n,n}
         //--------------------------------------------------------------------------------------------------------------------
         
         // weighted: 
         //...........................................................................................................................................
         if(fUsePhiWeights) fDirectCorrelationsW->Fill(120.,cos(n*phi1+n*phi2+n*phi3+n*phi4-n*phi5-n*phi6-n*phi7-n*phi8),
                                                           wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6*wPhi7*wPhi8);
         //...........................................................................................................................................
     
        } 
       } 
      } 
     }
    }  
   }
  }
 }
 
} // end of AliFlowAnalysisWithQCumulants::EvaluateNestedLoopsForIntegratedFlow(AliFlowEventSimple* anEvent)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CompareResultsFromNestedLoopsAndFromQVectorsForIntFlow(Bool_t useParticleWeights)
{
 // compare results needed for int. flow calculated with nested loops and with those calculated from Q-vectors

 cout<<endl;
 cout<<endl;
 cout<<"   *************************************"<<endl;
 cout<<"   **** cross-checking the formulas ****"<<endl;
 cout<<"   ****     for integrated flow     ****"<<endl;
 cout<<"   *************************************"<<endl;
 cout<<endl;
 cout<<endl;
 
 if(!(useParticleWeights))
 {
  cout<<"   **** results for non-weighted correlations: ****"<<endl;
  cout<<endl;
  cout<<"<2>_{1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(1)<<endl;
  cout<<"<2>_{1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(1)<<endl;
  cout<<endl;
  cout<<"<2>_{2n,2n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(2)<<endl;
  cout<<"<2>_{2n,2n} from nested loops = "<<fDirectCorrelations->GetBinContent(2)<<endl;
  cout<<endl;
  cout<<"<2>_{3n,3n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(3)<<endl;
  cout<<"<2>_{3n,3n} from nested loops = "<<fDirectCorrelations->GetBinContent(3)<<endl;
  cout<<endl;
  cout<<"<2>_{4n,4n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(4)<<endl;
  cout<<"<2>_{4n,4n} from nested loops = "<<fDirectCorrelations->GetBinContent(4)<<endl;
  cout<<endl; 
  cout<<"<3>_{2n|1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(6)<<endl;
  cout<<"<3>_{2n|1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(6)<<endl;
  cout<<endl;
  cout<<"<3>_{3n|2n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(7)<<endl;
  cout<<"<3>_{3n|2n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(7)<<endl;
  cout<<endl;
  cout<<"<3>_{4n,2n,2n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(8)<<endl;
  cout<<"<3>_{4n,2n,2n} from nested loops = "<<fDirectCorrelations->GetBinContent(8)<<endl;
  cout<<endl;
  cout<<"<3>_{4n,3n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(9)<<endl;
  cout<<"<3>_{4n,3n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(9)<<endl;
  cout<<endl; 
  cout<<"<4>_{1n,1n|1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(11)<<endl;
  cout<<"<4>_{1n,1n|1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(11)<<endl;
  cout<<endl;
  cout<<"<4>_{2n,1n|2n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(12)<<endl;
  cout<<"<4>_{2n,1n|2n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(12)<<endl;
  cout<<endl;
  cout<<"<4>_{2n,2n|2n,2n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(13)<<endl;
  cout<<"<4>_{2n,2n|2n,2n} from nested loops = "<<fDirectCorrelations->GetBinContent(13)<<endl;
  cout<<endl;
  cout<<"<4>_{3n|1n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(14)<<endl;
  cout<<"<4>_{3n|1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(14)<<endl;
  cout<<endl;
  cout<<"<4>_{3n,1n|3n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(15)<<endl;
  cout<<"<4>_{3n,1n|3n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(15)<<endl;
  cout<<endl;
  cout<<"<4>_{3n,1n|2n,2n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(16)<<endl;
  cout<<"<4>_{3n,1n|2n,2n} from nested loops = "<<fDirectCorrelations->GetBinContent(16)<<endl;
  cout<<endl; 
  cout<<"<4>_{4n|2n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(17)<<endl;
  cout<<"<4>_{4n|2n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(17)<<endl;
  cout<<endl;
  cout<<"<5>_{2n,1n|1n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(19)<<endl;
  cout<<"<5>_{2n,1n|1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(19)<<endl;
  cout<<endl;
  cout<<"<5>_{2n,2n|2n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(20)<<endl;
  cout<<"<5>_{2n,2n|2n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(20)<<endl;
  cout<<endl;
  cout<<"<5>_{3n,1n|2n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(21)<<endl;
  cout<<"<5>_{3n,1n|2n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(21)<<endl;
  cout<<endl;
  cout<<"<5>_{4n|1n,1n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(22)<<endl;
  cout<<"<5>_{4n|1n,1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(22)<<endl;
  cout<<endl;
  cout<<"<6>_{1n,1n,1n|1n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(24)<<endl;
  cout<<"<6>_{1n,1n,1n|1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(24)<<endl;
  cout<<endl; 
  cout<<"<6>_{2n,1n,1n|2n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(25)<<endl;
  cout<<"<6>_{2n,1n,1n|2n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(25)<<endl;
  cout<<endl;
  cout<<"<6>_{2n,2n|1n,1n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(26)<<endl;
  cout<<"<6>_{2n,2n|1n,1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(26)<<endl;
  cout<<endl; 
  cout<<"<6>_{3n,1n|1n,1n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(27)<<endl;
  cout<<"<6>_{3n,1n|1n,1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(27)<<endl;
  cout<<endl; 
  cout<<"<7>_{2n,1n,1n|1n,1n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(29)<<endl;
  cout<<"<7>_{2n,1n,1n|1n,1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(29)<<endl;
  cout<<endl; 
  cout<<"<8>_{1n,1n,1n,1n|1n,1n,1n,1n} from Q-vectors    = "<<fQCorrelations[0][0]->GetBinContent(31)<<endl;
  cout<<"<8>_{1n,1n,1n,1n|1n,1n,1n,1n} from nested loops = "<<fDirectCorrelations->GetBinContent(31)<<endl;
  cout<<endl; 
  cout<<endl; 
  cout<<"   **** results for non-weighted correction terms: ****"<<endl;
  cout<<endl;
  //.........................................................................................
  cout<<"<cos(n*phi1)> from Q-vectors    = "<<fQCorrections[0][0][1]->GetBinContent(1)<<endl;
  cout<<"<cos(n*phi1)> from nested loops = "<<fDirectCorrectionsCos->GetBinContent(1)<<endl;
  cout<<endl;
  cout<<"<sin(n*phi1)> from Q-vectors    = "<<fQCorrections[0][0][0]->GetBinContent(1)<<endl;
  cout<<"<sin(n*phi1)> from nested loops = "<<fDirectCorrectionsSin->GetBinContent(1)<<endl;
  cout<<endl;  
  cout<<"<cos(n*(phi1+phi2))> from Q-vectors    = "<<fQCorrections[0][0][1]->GetBinContent(2)<<endl;
  cout<<"<cos(n*(phi1+phi2))> from nested loops = "<<fDirectCorrectionsCos->GetBinContent(2)<<endl;
  cout<<endl;
  cout<<"<sin(n*(phi1+phi2))> from Q-vectors    = "<<fQCorrections[0][0][0]->GetBinContent(2)<<endl;
  cout<<"<sin(n*(phi1+phi2))> from nested loops = "<<fDirectCorrectionsSin->GetBinContent(2)<<endl;
  cout<<endl; 
  cout<<"<cos(n*(phi1-phi2-phi3))> from Q-vectors    = "<<fQCorrections[0][0][1]->GetBinContent(3)<<endl;
  cout<<"<cos(n*(phi1-phi2-phi3))> from nested loops = "<<fDirectCorrectionsCos->GetBinContent(3)<<endl;
  cout<<endl;
  cout<<"<sin(n*(phi1-phi2-phi3))> from Q-vectors    = "<<fQCorrections[0][0][0]->GetBinContent(3)<<endl;
  cout<<"<sin(n*(phi1-phi2-phi3))> from nested loops = "<<fDirectCorrectionsSin->GetBinContent(3)<<endl;
  cout<<endl;  
  //.........................................................................................
 }
 
 if(useParticleWeights)
 {
  cout<<"   **** results for weighted correlations: ****"<<endl;
  cout<<endl;
  //.........................................................................................
  cout<<"<w1 w2 cos(n*(phi1-phi2))> from Q-vectors         = "<<fQCorrelations[1][0]->GetBinContent(1)<<endl;
  cout<<"<<w1 w2 cos(n*(phi1-phi2))> from nested loops     = "<<fDirectCorrelationsW->GetBinContent(1)<<endl;
  cout<<endl;
  cout<<"<w1^2 w2^2 cos(2n*(phi1-phi2))> from Q-vectors    = "<<fQCorrelations[1][0]->GetBinContent(2)<<endl;
  cout<<"<w1^2 w2^2 cos(2n*(phi1-phi2))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(2)<<endl;
  cout<<endl;
  cout<<"<w1^3 w2^3 cos(3n*(phi1-phi2))> from Q-vectors    = "<<fQCorrelations[1][0]->GetBinContent(3)<<endl;
  cout<<"<w1^3 w2^3 cos(3n*(phi1-phi2))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(3)<<endl;
  cout<<endl;
  cout<<"<w1^4 w2^4 cos(4n*(phi1-phi2))> from Q-vectors    = "<<fQCorrelations[1][0]->GetBinContent(4)<<endl;
  cout<<"<w1^4 w2^4 cos(4n*(phi1-phi2))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(4)<<endl;
  cout<<endl;  
  /*
  cout<<"<w1^3 w2 cos(n*(phi1-phi2))> from Q-vectors       = "<<fQCorrelationsW->GetBinContent(5)<<endl;
  cout<<"<w1^3 w2 cos(n*(phi1-phi2))> from nested loops    = "<<fDirectCorrelationsW->GetBinContent(5)<<endl;
  cout<<endl;
  cout<<"<w1 w2 w3^2 cos(n*(phi1-phi2))> from Q-vectors    = "<<fQCorrelationsW->GetBinContent(6)<<endl;
  cout<<"<w1 w2 w3^2 cos(n*(phi1-phi2))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(6)<<endl;
  cout<<endl;
  cout<<"<w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))> from Q-vectors    = "<<fQCorrelationsW->GetBinContent(21)<<endl;
  cout<<"<w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(21)<<endl;
  cout<<endl;
  */ 
  cout<<"<w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> from Q-vectors    = "<<fQCorrelations[1][0]->GetBinContent(11)<<endl;
  cout<<"<w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> from nested loops = "<<fDirectCorrelationsW->GetBinContent(41)<<endl;
  cout<<endl;
  //.........................................................................................
  //cout<<endl; 
  //cout<<endl; 
  //cout<<"   **** results for weighted correction terms: ****"<<endl;
  //cout<<endl;
  //.........................................................................................
 }
 
} // end of AliFlowAnalysisWithQCumulants::CompareResultsFromNestedLoopsAndFromQVectorsForIntFlow(Bool_t useParticleWeights)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateQProductsForIntFlow()
{
 // calculate averages like <<2><4>>, <<2><6>>, <<4><6>>, etc. which are needed to calculate covariances 

 // binning of fQProducts is organized as follows:
 // 
 // 1st bin: <2><4> 
 // 2nd bin: <2><6>
 // 3rd bin: <2><8>
 // 4th bin: <4><6>
 // 5th bin: <4><8>
 // 6th bin: <6><8>

 Double_t dMult = (*fSMpk)(0,0); // multiplicity (number of particles used to determine the reaction plane)

 Double_t twoEBE = 0.; // <2>
 Double_t fourEBE = 0.; // <4>
 Double_t sixEBE = 0.; // <6>
 Double_t eightEBE = 0.; // <8>
 
 twoEBE = fQCorrelationsEBE[0]->GetBinContent(1);
 fourEBE = fQCorrelationsEBE[0]->GetBinContent(11);
 if(dMult>5) 
 {
  sixEBE = fQCorrelationsEBE[0]->GetBinContent(24);
  if(dMult>7) 
  { 
   eightEBE = fQCorrelationsEBE[0]->GetBinContent(31);
  }
 }
 
 // <2><4>
 if(dMult>3)
 {
  fQProducts[0][0]->Fill(0.5,twoEBE*fourEBE,dMult*(dMult-1)*dMult*(dMult-1)*(dMult-2)*(dMult-3));
 }
 
 // <2><6>
 if(dMult>5)
 {
  fQProducts[0][0]->Fill(1.5,twoEBE*sixEBE,dMult*(dMult-1)*dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5));
 }
 
 // <2><8>
 if(dMult>7)
 {
  fQProducts[0][0]->Fill(2.5,twoEBE*eightEBE,dMult*(dMult-1)*dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5)*(dMult-6)*(dMult-7));
 }
 
 // <4><6>
 if(dMult>5)
 {
  fQProducts[0][0]->Fill(3.5,fourEBE*sixEBE,dMult*(dMult-1)*(dMult-2)*(dMult-3)*dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5));
 }
 
 // <4><8>
 if(dMult>7)
 {
  fQProducts[0][0]->Fill(4.5,fourEBE*eightEBE,dMult*(dMult-1)*(dMult-2)*(dMult-3)*
                                 dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5)*(dMult-6)*(dMult-7));
 }
 
 // <6><8>
 if(dMult>7)
 {
  fQProducts[0][0]->Fill(5.5,sixEBE*eightEBE,dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5)*
                                dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5)*(dMult-6)*(dMult-7));
 }
 
} // end of AliFlowAnalysisWithQCumulants::CalculateQProductsForIntFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCovariancesForIntFlow(Bool_t useParticleWeights, TString eventWeights)
{
 // calculate covariances Cov(2,4), Cov(2,6), etc (needed to propagate errors for integrated flow)
 
 // binning of fCovariances[pW] is organized as follows:
 // 
 // 1st bin: <<2><4>>-<<2>><<4>> 
 // 2nd bin: <<2><6>>-<<2>><<6>>
 // 3rd bin: <<2><8>>-<<2>><<8>>
 // 4th bin: <<4><6>>-<<4>><<6>>
 // 5th bin: <<4><8>>-<<4>><<8>>
 // 6th bin: <<6><8>>-<<6>><<8>>
 
 // shortcuts for flags:
 Int_t pW = (Int_t)(useParticleWeights);
 Int_t eW = -1;
 
 if(eventWeights == "exact")
 {
  eW = 0;
 }
  
 for(Int_t power=0;power<2;power++)
 { 
  if(!(fQCorrelations[pW][eW] && fQProducts[pW][eW] && fCovariances[pW][eW] && fSumOfEventWeights[pW][eW][power] && fProductOfEventWeights[pW][eW]))
  { 
   cout<<"WARNING: fQCorrelations[pW][eW] && fQProducts[pW][eW] && fCovariances[pW][eW] && fSumOfEventWeights[pW][eW][power] && fProductOfEventWeights[pW][eW] is NULL in AFAWQC::CCFIF() !!!!"<<endl;
   cout<<"pW    = "<<pW<<endl;
   cout<<"eW    = "<<eW<<endl;
   cout<<"power = "<<power<<endl;
   exit(0);
  }
 } 
 
 // average 2-, 4-, 6- and 8-particle azimuthal correlations for all events:
 Double_t two = fCorrelations[pW][eW]->GetBinContent(1); // <<2>>  
 Double_t four = fCorrelations[pW][eW]->GetBinContent(11); // <<4>>  
 Double_t six = fCorrelations[pW][eW]->GetBinContent(24); // <<6>>  
 Double_t eight = fCorrelations[pW][eW]->GetBinContent(31); // <<8>> 
 // average products of 2-, 4-, 6- and 8-particle azimuthal correlations:  
 Double_t product24 = fQProducts[pW][eW]->GetBinContent(1); // <<2><4>>  
 Double_t product26 = fQProducts[pW][eW]->GetBinContent(2); // <<2><6>>  
 Double_t product28 = fQProducts[pW][eW]->GetBinContent(3); // <<2><8>>  
 Double_t product46 = fQProducts[pW][eW]->GetBinContent(4); // <<4><6>> 
 Double_t product48 = fQProducts[pW][eW]->GetBinContent(5); // <<4><8>>  
 Double_t product68 = fQProducts[pW][eW]->GetBinContent(6); // <<6><8>>
 // denominator in the expression for the unbiased estimator for covariance:
 Double_t denom24 = 0.; 
 if(fSumOfEventWeights[pW][eW][0]->GetBinContent(1) && fSumOfEventWeights[pW][eW][0]->GetBinContent(2))
 {
  denom24 = 1-(fProductOfEventWeights[pW][eW]->GetBinContent(1))/
              (fSumOfEventWeights[pW][eW][0]->GetBinContent(1) * fSumOfEventWeights[pW][eW][0]->GetBinContent(2));
 }
 Double_t denom26 = 0.; 
 if(fSumOfEventWeights[pW][eW][0]->GetBinContent(1) && fSumOfEventWeights[pW][eW][0]->GetBinContent(3))
 {
  denom26 = 1-(fProductOfEventWeights[pW][eW]->GetBinContent(2))/
              (fSumOfEventWeights[pW][eW][0]->GetBinContent(1) * fSumOfEventWeights[pW][eW][0]->GetBinContent(3));
 }
 Double_t denom28 = 0.; 
 if(fSumOfEventWeights[pW][eW][0]->GetBinContent(1) && fSumOfEventWeights[pW][eW][0]->GetBinContent(4))
 {
  denom28 = 1-(fProductOfEventWeights[pW][eW]->GetBinContent(3))/
              (fSumOfEventWeights[pW][eW][0]->GetBinContent(1) * fSumOfEventWeights[pW][eW][0]->GetBinContent(4));
 }
 Double_t denom46 = 0.; 
 if(fSumOfEventWeights[pW][eW][0]->GetBinContent(2) && fSumOfEventWeights[pW][eW][0]->GetBinContent(3))
 {
  denom46 = 1-(fProductOfEventWeights[pW][eW]->GetBinContent(4))/
              (fSumOfEventWeights[pW][eW][0]->GetBinContent(2) * fSumOfEventWeights[pW][eW][0]->GetBinContent(3));
 }
 Double_t denom48 = 0.; 
 if(fSumOfEventWeights[pW][eW][0]->GetBinContent(2) && fSumOfEventWeights[pW][eW][0]->GetBinContent(4))
 {
  denom48 = 1-(fProductOfEventWeights[pW][eW]->GetBinContent(5))/
              (fSumOfEventWeights[pW][eW][0]->GetBinContent(2) * fSumOfEventWeights[pW][eW][0]->GetBinContent(4));
 }
 Double_t denom68 = 0.; 
 if(fSumOfEventWeights[pW][eW][0]->GetBinContent(3) && fSumOfEventWeights[pW][eW][0]->GetBinContent(4))
 {
  denom68 = 1-(fProductOfEventWeights[pW][eW]->GetBinContent(6))/
              (fSumOfEventWeights[pW][eW][0]->GetBinContent(3) * fSumOfEventWeights[pW][eW][0]->GetBinContent(4));
 }
 
 // final covariances: 
 Double_t cov24 = 0.;
 if(denom24)
 {
  cov24 = (product24-two*four)/denom24; // Cov(<2>,<4>) 
  fCovariances[pW][eW]->SetBinContent(1,cov24);
 } 
 Double_t cov26 = 0.;
 if(denom26)
 {
  cov26 = (product26-two*six)/denom26; // Cov(<2>,<6>) 
  fCovariances[pW][eW]->SetBinContent(2,cov26);
 } 
 Double_t cov28 = 0.;
 if(denom28)
 {
  cov28 = (product28-two*eight)/denom28; // Cov(<2>,<8>) 
  fCovariances[pW][eW]->SetBinContent(3,cov28);
 } 
 Double_t cov46 = 0.;
 if(denom46)
 {
  cov46 = (product46-four*six)/denom46; // Cov(<4>,<6>) 
  fCovariances[pW][eW]->SetBinContent(4,cov46);
 } 
 Double_t cov48 = 0.;
 if(denom48)
 {
  cov48 = (product48-four*eight)/denom48; // Cov(<4>,<8>) 
  fCovariances[pW][eW]->SetBinContent(5,cov48);
 } 
 Double_t cov68 = 0.;
 if(denom68)
 {
  cov68 = (product68-six*eight)/denom68; // Cov(<6>,<8>) 
  fCovariances[pW][eW]->SetBinContent(6,cov68);
 } 

} // end of AliFlowAnalysisWithQCumulants::CalculateCovariancesForIntFlow(Bool_t useParticleWeights, TString eventWeights)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::FinalizeCorrelationsForIntFlow(Bool_t useParticleWeights, TString eventWeights) // to be improved (there is better way to implement this method)
{
 // From fQCorrelations[pW][eW] access measured correlations and spread, calculate statistical errors and store the
 // final results and statistical errors for correlations in fCorrelations[pW][eW] (binning is the same as in fQCorrelations[pW][eW]).
 
 // shortcuts for flags:
 Int_t pW = (Int_t)(useParticleWeights);
 
 Int_t eW = -1;
 
 if(eventWeights == "exact")
 {
  eW = 0;
 }
 
 for(Int_t power=0;power<2;power++)
 {
  if(!(fQCorrelations[pW][eW] && fCorrelations[pW][eW] && fSumOfEventWeights[pW][eW][power])) 
  {
   cout<<"WARNING: fQCorrelations[pW][eW] && fCorrelations[pW][eW] && fSumOfEventWeights[pW][eW][power] is NULL in AFAWQC::FCFIF() !!!!"<<endl;
   cout<<"pW    = "<<pW<<endl;
   cout<<"eW    = "<<eW<<endl;
   cout<<"power = "<<power<<endl;
   exit(0);
  }
 }
 
 // <<2>>:
 Double_t correlation2p = fQCorrelations[pW][eW]->GetBinContent(1);  
 Double_t spread2p = fQCorrelations[pW][eW]->GetBinError(1); 
 Double_t sumOfEventWeightsLinear2p = fSumOfEventWeights[pW][eW][0]->GetBinContent(1);
 Double_t sumOfEventWeightsQuadratic2p = fSumOfEventWeights[pW][eW][1]->GetBinContent(1);
 // stat.error = termA * spread * termB:
 // termB = 1/sqrt(1-termA^2)
 Double_t termA2p = 0.;  
 Double_t termB2p = 0.;
 Double_t statError2p = 0.;
 if(sumOfEventWeightsLinear2p)
 {
  termA2p = pow(sumOfEventWeightsQuadratic2p,0.5)/sumOfEventWeightsLinear2p;
 } else
   {
    cout<<"WARNING: sumOfEventWeightsLinear2p == 0 in in AFAWQC::FCFIF() !!!!"<<endl;
   }
      
 if(1-pow(termA2p,2.) > 0)
 {
  termB2p = 1./pow(1-pow(termA2p,2.),0.5);
 } else
   {
    cout<<"WARNING: 1-pow(termA2p,2.) <= 0 in in AFAWQC::FCFIF() !!!!"<<endl;   
   }     
   
 statError2p = termA2p*spread2p*termB2p;
 fCorrelations[pW][eW]->SetBinContent(1,correlation2p);
 fCorrelations[pW][eW]->SetBinError(1,statError2p);
         
 // <<4>>:
 Double_t correlation4p = fQCorrelations[pW][eW]->GetBinContent(11);  
 Double_t spread4p = fQCorrelations[pW][eW]->GetBinError(11); 
 Double_t sumOfEventWeightsLinear4p = fSumOfEventWeights[pW][eW][0]->GetBinContent(2);
 Double_t sumOfEventWeightsQuadratic4p = fSumOfEventWeights[pW][eW][1]->GetBinContent(2);
 Double_t termA4p = 0.;  
 Double_t termB4p = 0.;
 Double_t statError4p = 0.;
 if(sumOfEventWeightsLinear4p)
 {
  termA4p = pow(sumOfEventWeightsQuadratic4p,0.5)/sumOfEventWeightsLinear4p;
 } else
   {
    cout<<"WARNING: sumOfEventWeightsLinear4p == 0 in in AFAWQC::FCFIF() !!!!"<<endl;
   }
 if(1-pow(termA4p,2.) > 0)
 {
  termB4p = 1./pow(1-pow(termA4p,2.),0.5);
 } else
   {
    cout<<"WARNING: 1-pow(termA4p,2.) <= 0 in in AFAWQC::FCFIF() !!!!"<<endl;   
   }  
   
 statError4p = termA4p*spread4p*termB4p;
 fCorrelations[pW][eW]->SetBinContent(11,correlation4p);
 fCorrelations[pW][eW]->SetBinError(11,statError4p);

 // <<6>>:
 Double_t correlation6p = fQCorrelations[pW][eW]->GetBinContent(24);  
 Double_t spread6p = fQCorrelations[pW][eW]->GetBinError(24); 
 Double_t sumOfEventWeightsLinear6p = fSumOfEventWeights[pW][eW][0]->GetBinContent(3);
 Double_t sumOfEventWeightsQuadratic6p = fSumOfEventWeights[pW][eW][1]->GetBinContent(3);
 Double_t termA6p = 0.;  
 Double_t termB6p = 0.;
 Double_t statError6p = 0.;
 if(sumOfEventWeightsLinear6p)
 {
  termA6p = pow(sumOfEventWeightsQuadratic6p,0.5)/sumOfEventWeightsLinear6p;
 } else
   {
    cout<<"WARNING: sumOfEventWeightsLinear6p == 0 in in AFAWQC::FCFIF() !!!!"<<endl;
   }
 if(1-pow(termA6p,2.) > 0)
 {
  termB6p = 1./pow(1-pow(termA6p,2.),0.5);
 } else
   {
    cout<<"WARNING: 1-pow(termA6p,2.) <= 0 in in AFAWQC::FCFIF() !!!!"<<endl;   
   }  
   
 statError6p = termA6p*spread6p*termB6p;
 fCorrelations[pW][eW]->SetBinContent(24,correlation6p);
 fCorrelations[pW][eW]->SetBinError(24,statError6p);
              
 // <<8>>             
 Double_t correlation8p = fQCorrelations[pW][eW]->GetBinContent(31);  
 // ...
 fCorrelations[pW][eW]->SetBinContent(31,correlation8p);
 // ... 
                          
} // end of AliFlowAnalysisWithQCumulants::FinalizeCorrelationsForIntFlow(Bool_t useParticleWeights, TString eventWeights)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::FillAverageMultiplicities(Int_t nRP)
{
 // Fill profile fAverageMultiplicity to hold average multiplicities and number of events for events with nRP>=0, nRP>=1, ... , and nRP>=8
 
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
 
 if(!fAvMultiplicity)
 {
  cout<<"WARNING: fAvMultiplicity is NULL in AFAWQC::FAM() !!!!"<<endl;
  exit(0);
 }
 
 if(nRP<0)
 {
  cout<<"WARNING: nRP<0 in in AFAWQC::FAM() !!!!"<<endl;
  exit(0);
 }
 
 for(Int_t i=0;i<9;i++)
 {
  if(nRP>=i) fAvMultiplicity->Fill(i+0.5,nRP,1);
 }
 
} // end of AliFlowAnalysisWithQCumulants::FillAverageMultiplicities(nRP)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCumulantsForIntFlow(Bool_t useParticleWeights, TString eventWeights)
{
 // calculate cumulants from measured correlations and store results in fCumulants[pW][eW][0]. 
 // (Remark: these cumulants are biased by non-uniform acceptance, corrected cumulants are stored in fCumulants[pW][eW][1].)
 
 // Binning of fCumulants[pW][nua] is organized as follows:
 //  1st bin: QC{2}
 //  2nd bin: QC{4}
 //  3rd bin: QC{6}
 //  4th bin: QC{8}
 
 // shortcuts for the flags:
 Int_t pW = (Int_t)(useParticleWeights); // (0=weights not used, 1=weights used)
 
 Int_t eW = -1;
 
 if(eventWeights == "exact")
 {
  eW = 0;
 }
 
 if(!(fCorrelations[pW][eW] && fCovariances[pW][eW] && fCumulants[pW][eW][0] && fAvMultiplicity))
 {
  cout<<"WARNING: fCorrelations[pW][eW] && fCovariances[pW][eW] && fCumulants[pW][eW][0] && fAvMultiplicity is NULL in AFAWQC::CCFIF() !!!!"<<endl;
  cout<<"pW = "<<pW<<endl;
  cout<<"eW = "<<eW<<endl;
  exit(0);
 }
 
 // correlations:
 Double_t two = fCorrelations[pW][eW]->GetBinContent(1); // <<2>> 
 Double_t four = fCorrelations[pW][eW]->GetBinContent(11); // <<4>>  
 Double_t six = fCorrelations[pW][eW]->GetBinContent(24); // <<6>> 
 Double_t eight = fCorrelations[pW][eW]->GetBinContent(31); // <<8>>  
 // stat. error of correlations:
 Double_t twoError = fCorrelations[pW][eW]->GetBinError(1); // stat. error of <<2>>   
 Double_t fourError = fCorrelations[pW][eW]->GetBinError(11); // stat. error of <<4>>  
 Double_t sixError = fCorrelations[pW][eW]->GetBinError(24); // stat. error of <<6>>  
 //Double_t eightError = fQCorrelations[pW]->GetBinError(31); // stat. error of <<8>>
 // spread of correlations:
 Double_t twoSpread = 0.; // spread of <<2>>
 Double_t fourSpread = 0.; // spread of <<6>>
 Double_t sixSpread = 0.; // spread of <<8>>
 //Double_t eightSpread = 0.;
 // stat. error = prefactor * spread:
 Double_t twoPrefactor = 0; 
 if(fSumOfEventWeights[pW][eW][0]->GetBinContent(1))
 {
  twoPrefactor = pow(fSumOfEventWeights[pW][eW][1]->GetBinContent(1),0.5)/fSumOfEventWeights[pW][eW][0]->GetBinContent(1);
  if(twoPrefactor) twoSpread = twoError/twoPrefactor;
 }
 Double_t fourPrefactor = 0; 
 if(fSumOfEventWeights[pW][eW][0]->GetBinContent(2))
 {
  fourPrefactor = pow(fSumOfEventWeights[pW][eW][1]->GetBinContent(2),0.5)/fSumOfEventWeights[pW][eW][0]->GetBinContent(2);
  if(fourPrefactor) fourSpread = fourError/fourPrefactor;
 }
 Double_t sixPrefactor = 0; 
 if(fSumOfEventWeights[pW][eW][0]->GetBinContent(3))
 {
  sixPrefactor = pow(fSumOfEventWeights[pW][eW][1]->GetBinContent(3),0.5)/fSumOfEventWeights[pW][eW][0]->GetBinContent(3);
  if(sixPrefactor) sixSpread = sixError/sixPrefactor;
 }
 // ... 8th
 // covariances:
 Double_t cov24 = fCovariances[pW][eW]->GetBinContent(1); // Cov(<2>,<4>)
 Double_t cov26 = fCovariances[pW][eW]->GetBinContent(2); // Cov(<2>,<6>) 
 //Double_t cov28 = fCovariances[pW]->GetBinContent(3); // Cov(<2>,<8>)
 Double_t cov46 = fCovariances[pW][eW]->GetBinContent(4); // Cov(<4>,<6>)
 //Double_t cov48 = fCovariances[pW]->GetBinContent(5); // Cov(<4>,<8>) 
 //Double_t cov68 = fCovariances[pW]->GetBinContent(6); // Cov(<6>,<8>)  
 // cumulants: 
 Double_t qc2 = 0.; // QC{2}
 Double_t qc4 = 0.; // QC{4}
 Double_t qc6 = 0.; // QC{6}
 Double_t qc8 = 0.; // QC{8}
 
 if(two) qc2 = two; 
 if(four) qc4 = four-2.*pow(two,2.); 
 if(six) qc6 = six-9.*two*four+12.*pow(two,3.); 
 if(eight) qc8 = eight-16.*two*six-18.*pow(four,2.)+144.*pow(two,2.)*four-144.*pow(two,4.); 
 
 // stat. error of cumulants:
 Double_t qc2Error = 0; // stat. error of QC{2}
 Double_t qc4Error = 0; // stat. error of QC{4}
 Double_t qc6Error = 0; // stat. error of QC{6}
 // Double_t qc8Error = 0; // stat. error of QC{8} // to be improved (calculated)
 // spread of cumulants:
 //Double_t qc2Spread = 0; // spread of QC{2}
 Double_t qc4Spread = 0; // spread of QC{4}
 Double_t qc6Spread = 0; // spread of QC{6}
 // Double_t qc8Spread = 0; // spread of QC{8} // to be improved (calculated)
 
 qc2Error = twoError; // final stat. error of QC{2}
 
 if(16.*pow(two,2.)*pow(twoSpread,2.)+pow(fourSpread,2.)-8.*two*cov24 >= 0.)
 {
  qc4Spread = pow(16.*pow(two,2.)*pow(twoSpread,2.)+pow(fourSpread,2.)-8.*two*cov24,0.5); 
 } else
   {
    cout<<"WARNING: spread of QC{4} is imaginary !!!!"<<endl;
   }
  
 qc4Error = fourPrefactor*qc4Spread; // final stat. error of QC{4}  
     
 if(81.*pow(4.*pow(two,2.)-four,2.)*pow(twoSpread,2.)   
    + 81.*pow(two,2.)*pow(fourSpread,2.)+pow(sixSpread,2.) 
    - 162.*(4.*pow(two,2.)-four)*two*cov24
    + 18.*(4.*pow(two,2.)-four)*cov26
    - 18.*two*cov46 >= 0.)
 {  
  qc6Spread = pow(81.*pow(4.*pow(two,2.)-four,2.)*pow(twoSpread,2.)   
                  + 81.*pow(two,2.)*pow(fourSpread,2.)+pow(sixSpread,2.) 
                  - 162.*(4.*pow(two,2.)-four)*two*cov24
                  + 18.*(4.*pow(two,2.)-four)*cov26
                  - 18.*two*cov46,0.5);
 } else 
   {
    cout<<"WARNING: stat. error of QC{6} is imaginary !!!!"<<endl;
   }          
   
 qc6Error = sixPrefactor*qc6Spread; // final stat. error of QC{6} 
   
 // store the results and statistical errors for cumulants:
 fCumulants[pW][eW][0]->SetBinContent(1,qc2);
 fCumulants[pW][eW][0]->SetBinError(1,qc2Error);
 fCumulants[pW][eW][0]->SetBinContent(2,qc4);
 fCumulants[pW][eW][0]->SetBinError(2,qc4Error);
 fCumulants[pW][eW][0]->SetBinContent(3,qc6);
 fCumulants[pW][eW][0]->SetBinError(3,qc6Error);
 fCumulants[pW][eW][0]->SetBinContent(4,qc8); 
 // fCumulants[pW]->SetBinError(4,qc8Error); // to be improved (calculated)
   
} // end of AliFlowAnalysisWithQCumulants::CalculateCumulantsForIntFlow(Bool_t useParticleWeights, TString eventWeights)


//================================================================================================================================ 


void AliFlowAnalysisWithQCumulants::CalculateIntFlow(Bool_t useParticleWeights, TString eventWeights, Bool_t correctedForNUA)
{
 // calculate final results for integrated flow 
 
 // Results for integrated flow are stored in fInfFlow[pW][nua]. Binning of fIntFlow[pW][nua] is organized as follows:
 //  1st bin: v{2}
 //  2nd bin: v{4}
 //  3rd bin: v{6}
 //  4th bin: v{8}
 
 // shortcuts for the flags:
 Int_t pW = (Int_t)(useParticleWeights); // 0=pWeights not useed, 1=pWeights used
 Int_t nua = (Int_t)(correctedForNUA); // 0=not corrected for NUA, 1=corrected for NUA
 Int_t eW = -1;
 
 if(eventWeights = "exact")
 {
  eW = 0;
 }
   
 if(!(fCumulants[pW][eW][nua] && fIntFlow[pW][eW][nua]))
 {
  cout<<"WARNING: fCumulants[pW][eW][nua] && fIntFlow[pW][eW][nua] is NULL in AFAWQC::CIF() !!!!"<<endl;
  cout<<"pW = "<<pW<<endl;
  cout<<"eW = "<<eW<<endl;
  cout<<"nua = "<<nua<<endl;
  exit(0);
 }
   
 // cumulants:
 Double_t qc2 = fCumulants[pW][eW][nua]->GetBinContent(1); // QC{2}  
 Double_t qc4 = fCumulants[pW][eW][nua]->GetBinContent(2); // QC{4}  
 Double_t qc6 = fCumulants[pW][eW][nua]->GetBinContent(3); // QC{6}  
 Double_t qc8 = fCumulants[pW][eW][nua]->GetBinContent(4); // QC{8}
   
 // cumulants' statistical errors:
 Double_t qc2Error = fCumulants[pW][eW][nua]->GetBinError(1); // error of QC{2}  
 Double_t qc4Error = fCumulants[pW][eW][nua]->GetBinError(2); // error of QC{4}  
 Double_t qc6Error = fCumulants[pW][eW][nua]->GetBinError(3); // error of QC{6}  
 //Double_t qc8Error = fCumulants[pW][nua]->GetBinError(4); // error of QC{8}
   
 // integrated flow estimates:
 Double_t v2 = 0.; // v{2}  
 Double_t v4 = 0.; // v{4}  
 Double_t v6 = 0.; // v{6}  
 Double_t v8 = 0.; // v{8}
 
 // calculate integrated flow estimates from cumulants: 
 if(qc2>=0.) v2 = pow(qc2,1./2.); 
 if(qc4<=0.) v4 = pow(-1.*qc4,1./4.); 
 if(qc6>=0.) v6 = pow((1./4.)*qc6,1./6.); 
 if(qc8<=0.) v8 = pow((-1./33.)*qc8,1./8.); 
   
 // statistical errors of integrated flow estimates:
 Double_t v2Error = 0.; // error of v{2}  
 Double_t v4Error = 0.; // error of v{4}  
 Double_t v6Error = 0.; // error of v{6}  
 Double_t v8Error = 0.; // error of v{8}
   
 // calculate statistical errors for integrated flow estimates:
 if(qc2>0.) v2Error = (1./(2.*pow(qc2,0.5)))*qc2Error; 
 if(qc4<0.) v4Error = (1./(4.*pow(-1.*qc4,3./4.)))*qc4Error; 
 if(qc6>0.) v6Error = (1./(6.*pow(2.,1./3.)*pow(qc6,5./6.)))*qc6Error;
 if(qc8<0.) v8Error = 0.; // to be improved (calculated)
  
 // store final results and statistical errors for integrated flow:
 fIntFlow[pW][eW][nua]->SetBinContent(1,v2);
 fIntFlow[pW][eW][nua]->SetBinError(1,v2Error);
 fIntFlow[pW][eW][nua]->SetBinContent(2,v4);
 fIntFlow[pW][eW][nua]->SetBinError(2,v4Error);
 fIntFlow[pW][eW][nua]->SetBinContent(3,v6);
 fIntFlow[pW][eW][nua]->SetBinError(3,v6Error);
 fIntFlow[pW][eW][nua]->SetBinContent(4,v8);
 fIntFlow[pW][eW][nua]->SetBinError(4,v8Error);
     
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlow(Bool_t useParticleWeights, TString eventWeights, Bool_t correctedForNUA)


//================================================================================================================================ 


void AliFlowAnalysisWithQCumulants::FillCommonHistResultsIntFlow(Bool_t useParticleWeights, TString eventWeights, Bool_t correctedForNUA)
{
 // fill in AliFlowCommonHistResults histograms relevant for 'NONAME' integrated flow (to be improved (name))
 
 // shortcuts for the flags:
 Int_t pW = (Int_t)(useParticleWeights); // 0=pWeights not useed, 1=pWeights used
 Int_t nua = (Int_t)(correctedForNUA); // 0=not corrected for NUA, 1=corrected for NUA
 Int_t eW = -1;
 
 if(eventWeights == "exact")
 {
  eW = 0;
 }
 
 if(!fIntFlow[pW][eW][nua])
 {
  cout<<"WARNING: fIntFlow[pW][eW][nua] is NULL in AFAWQC::FCHRIF() !!!!"<<endl;
  cout<<"pW = "<<pW<<endl;
  cout<<"eW = "<<eW<<endl;
  cout<<"nua = "<<nua<<endl;
  exit(0); 
 }  
    
 if(!(fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th))
 {
  cout<<"WARNING: fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th"<<endl; 
  cout<<"         is NULL in AFAWQC::FCHRIF() !!!!"<<endl;
  exit(0);
 }
 
 Double_t v2 = fIntFlow[pW][eW][nua]->GetBinContent(1);
 Double_t v4 = fIntFlow[pW][eW][nua]->GetBinContent(2);
 Double_t v6 = fIntFlow[pW][eW][nua]->GetBinContent(3);
 Double_t v8 = fIntFlow[pW][eW][nua]->GetBinContent(4);
  
 Double_t v2Error = fIntFlow[pW][eW][nua]->GetBinError(1);
 Double_t v4Error = fIntFlow[pW][eW][nua]->GetBinError(2);
 Double_t v6Error = fIntFlow[pW][eW][nua]->GetBinError(3);
 Double_t v8Error = fIntFlow[pW][eW][nua]->GetBinError(4);
 
 fCommonHistsResults2nd->FillIntegratedFlow(v2,v2Error); 
 fCommonHistsResults4th->FillIntegratedFlow(v4,v4Error); 
 fCommonHistsResults6th->FillIntegratedFlow(v6,v6Error); 
 fCommonHistsResults8th->FillIntegratedFlow(v8,v8Error);  

} // end of AliFlowAnalysisWithQCumulants::FillCommonHistResultsIntFlow(Bool_t useParticleWeights, TString eventWeights, Bool_t correctedForNUA)


//================================================================================================================================ 


void AliFlowAnalysisWithQCumulants::ApplyCorrectionForNonUniformAcceptanceToCumulantsForIntFlow(Bool_t useParticleWeights, TString eventWeights)
{
 // apply correction for non-uniform acceptance to cumulants for integrated flow 
 // (Remark: non-corrected cumulants are accessed from fCumulants[pW][0], corrected cumulants are stored in fCumulants[pW][1])
   
 // shortcuts for the flags:
 Int_t pW = (Int_t)(useParticleWeights); // 0=pWeights not used, 1=pWeights used
 Int_t eW = -1;
 
 if(eventWeights == "exact")
 {
  eW = 0;
 }
 
 if(!(fCumulants[pW][eW][0] && fCumulants[pW][eW][1] && fCorrections[pW][eW]))
 {
  cout<<"WARNING: fCumulants[pW][eW][0] && fCumulants[pW][eW][1] && fCorrections[pW][eW] is NULL in AFAWQC::ACFNUATCFIF() !!!!"<<endl;
  cout<<"pW = "<<pW<<endl;
  cout<<"eW = "<<eW<<endl;
  exit(0);
 } 
 
 // non-corrected cumulants:
 Double_t qc2 = fCumulants[pW][eW][0]->GetBinContent(1); 
 Double_t qc4 = fCumulants[pW][eW][0]->GetBinContent(2); 
 Double_t qc6 = fCumulants[pW][eW][0]->GetBinContent(3); 
 Double_t qc8 = fCumulants[pW][eW][0]->GetBinContent(4); 
 // statistical error of non-corrected cumulants:  
 Double_t qc2Error = fCumulants[pW][eW][0]->GetBinError(1); 
 Double_t qc4Error = fCumulants[pW][eW][0]->GetBinError(2); 
 Double_t qc6Error = fCumulants[pW][eW][0]->GetBinError(3); 
 Double_t qc8Error = fCumulants[pW][eW][0]->GetBinError(4); 
 // corrections for non-uniform acceptance:
 Double_t qc2Correction = fCorrections[pW][eW]->GetBinContent(1); 
 Double_t qc4Correction = fCorrections[pW][eW]->GetBinContent(2); 
 Double_t qc6Correction = fCorrections[pW][eW]->GetBinContent(3); 
 Double_t qc8Correction = fCorrections[pW][eW]->GetBinContent(4); 
 // corrected cumulants:
 Double_t qc2Corrected = qc2 + qc2Correction;
 Double_t qc4Corrected = qc4 + qc4Correction;
 Double_t qc6Corrected = qc6 + qc6Correction;
 Double_t qc8Corrected = qc8 + qc8Correction;
  
 // ... to be improved (I need here also to correct error of QCs for NUA. 
 // For simplicity sake I assume at the moment that this correction is negliglible, but it will be added eventually...)
 
 // store corrected results and statistical errors for cumulants:   
 fCumulants[pW][eW][1]->SetBinContent(1,qc2Corrected);
 fCumulants[pW][eW][1]->SetBinContent(2,qc4Corrected);
 fCumulants[pW][eW][1]->SetBinContent(3,qc6Corrected);
 fCumulants[pW][eW][1]->SetBinContent(4,qc8Corrected);
 fCumulants[pW][eW][1]->SetBinError(1,qc2Error); // to be improved (correct also qc2Error for NUA)
 fCumulants[pW][eW][1]->SetBinError(2,qc4Error); // to be improved (correct also qc4Error for NUA)
 fCumulants[pW][eW][1]->SetBinError(3,qc6Error); // to be improved (correct also qc6Error for NUA)
 fCumulants[pW][eW][1]->SetBinError(4,qc8Error); // to be improved (correct also qc8Error for NUA)  
 
} // end of AliFlowAnalysisWithQCumulants::ApplyCorrectionForNonUniformAcceptanceToCumulantsForIntFlow(Bool_t useParticleWeights, TString eventWeights)


//================================================================================================================================

  
void AliFlowAnalysisWithQCumulants::PrintQuantifyingCorrectionsForNonUniformAcceptance(Bool_t useParticleWeights, TString eventWeights)
{
 // print on the screen QC{n,biased}/QC{n,corrected}
 
 // shortcuts for the flags:
 Int_t pW = (Int_t)(useParticleWeights); // 0=pWeights not used, 1=pWeights used
 
 Int_t eW = -1;
 
 if(eventWeights == "exact")
 {
  eW = 0;
 } 
 
 if(!(fCumulants[pW][eW][0] && fCumulants[pW][eW][1]))
 {
  cout<<"WARNING: fCumulants[pW][eW][0] && fCumulants[pW][eW][1] is NULL in AFAWQC::PQCFNUA() !!!!"<<endl;
  cout<<"pW = "<<pW<<endl;
  cout<<"eW = "<<eW<<endl;
  exit(0);
 }
   
 cout<<endl;
 cout<<" Quantifying the bias to Q-cumulants from"<<endl;
 cout<<"  non-uniform acceptance of the detector:"<<endl;
 cout<<endl;
  
 if(fCumulants[pW][eW][1]->GetBinContent(1)) 
 { 
  cout<<"  QC{2,biased}/QC{2,corrected} = "<<(fCumulants[pW][eW][0]->GetBinContent(1))/(fCumulants[pW][eW][1]->GetBinContent(1))<<endl;   
 }
 if(fCumulants[pW][eW][1]->GetBinContent(2)) 
 { 
  cout<<"  QC{4,biased}/QC{4,corrected} = "<<fCumulants[pW][eW][0]->GetBinContent(2)/fCumulants[pW][eW][1]->GetBinContent(2)<<endl;   
 }
 
 cout<<endl;
   
} // end of AliFlowAnalysisWithQCumulants::PrintQuantifyingCorrectionsForNonUniformAcceptance(Bool_t useParticleWeights, TString eventWeights)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateWeightedCorrelationsForIntegratedFlow()
{
 // calculate all weighted correlations needed for integrated flow and store them in 1D profile fQCorrelations[1][eW] and fQExtraCorrelations[1][eW]
 
 for(Int_t eW=0;eW<2;eW++)
 {
  if(!(fQCorrelationsEBE[1] && fQCorrelations[1][eW]))
  {
   cout<<"WARNING: fQCorrelationsEBE[1] && fQCorrelations[1][eW] is NULL in AFAWQC::CWCFIF() !!!!"<<endl;
   cout<<"eW = "<<eW<<endl;
   exit(0);
  }
 }
 
 // Remark 1: binning of fQCorrelations[W] is organized as follows:
 //..............................................................................................
 //       ---- bins 1-20: 2-particle correlations ----
 // 1st bin: two1n1nW1W1 = <w1 w2 cos(n*(phi1-phi2))>
 // 2nd bin: two2n2nW2W2 = <w1^2 w2^2 cos(2n*(phi1-phi2))>
 // 3rd bin: two3n3nW3W3 = <w1^3 w2^3 cos(3n*(phi1-phi2))>
 // 4th bin: two4n4nW4W4 = <w1^4 w2^4 cos(4n*(phi1-phi2))>
 // 5th bin: two1n1nW3W1 = <w1^3 w2 cos(n*(phi1-phi2))>
 // 6th bin: two1n1nW1W1W2 = <w1 w2 w3^2 cos(n*(phi1-phi2))>  
 //       ---- bins 21-40: 3-particle correlations ----
 // 21st bin: three2n1n1nW2W1W1 = <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))> 
 //       ---- bins 41-60: 4-particle correlations ----
 // 41st bin: four1n1n1n1nW1W1W1W1 = <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
 //       ---- bins 61-80: 5-particle correlations ---- 
 //       ---- bins 81-100: 6-particle correlations ----
 //       ---- bins 101-120: 7-particle correlations ----
 //       ---- bins 121-140: 8-particle correlations ----
 //..............................................................................................
 
 // multiplicity (number of particles used to determine the reaction plane)
 Double_t dMult = (*fSMpk)(0,0);
 
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
 Double_t dM11 = (*fSMpk)(1,1)-(*fSMpk)(0,2); // dM11 = sum_{i,j=1,i!=j}^M w_i w_j
 Double_t dM22 = (*fSMpk)(1,2)-(*fSMpk)(0,4); // dM22 = sum_{i,j=1,i!=j}^M w_i^2 w_j^2
 Double_t dM33 = (*fSMpk)(1,3)-(*fSMpk)(0,6); // dM33 = sum_{i,j=1,i!=j}^M w_i^3 w_j^3
 Double_t dM44 = (*fSMpk)(1,4)-(*fSMpk)(0,8); // dM44 = sum_{i,j=1,i!=j}^M w_i^4 w_j^4
 Double_t dM31 = (*fSMpk)(0,3)*(*fSMpk)(0,1)-(*fSMpk)(0,4); // dM31 = sum_{i,j=1,i!=j}^M w_i^3 w_j
 Double_t dM211 = (*fSMpk)(0,2)*(*fSMpk)(1,1)-2.*(*fSMpk)(0,3)*(*fSMpk)(0,1)
                - (*fSMpk)(1,2)+2.*(*fSMpk)(0,4); // dM211 = sum_{i,j,k=1,i!=j!=k}^M w_i^2 w_j w_k
 Double_t dM1111 = (*fSMpk)(3,1)-6.*(*fSMpk)(0,2)*(*fSMpk)(1,1)  
                 + 8.*(*fSMpk)(0,3)*(*fSMpk)(0,1)
                 + 3.*(*fSMpk)(1,2)-6.*(*fSMpk)(0,4); // dM1111 = sum_{i,j,k,l=1,i!=j!=k!=l}^M w_i w_j w_k w_l
 //..............................................................................................

 //  ***********************************************
 //  **** weighted multi-particle correlations: ****
 //  ***********************************************
 //.............................................................................................. 
 // weighted 2-particle correlations:
 Double_t two1n1nW1W1 = 0.; // <w1 w2 cos(n*(phi1-phi2))>
 Double_t two2n2nW2W2 = 0.; // <w1^2 w2^2 cos(2n*(phi1-phi2))>
 Double_t two3n3nW3W3 = 0.; // <w1^3 w2^3 cos(3n*(phi1-phi2))>
 Double_t two4n4nW4W4 = 0.; // <w1^4 w2^4 cos(4n*(phi1-phi2))>
 Double_t two1n1nW3W1 = 0.; // <w1^3 w2 cos(n*(phi1-phi2))>
 Double_t two1n1nW1W1W2 = 0.; // <w1 w2 w3^2 cos(n*(phi1-phi2))> 
 
 if(dMult>1) 
 { 
  if(dM11)
  {
   two1n1nW1W1 = (pow(dReQ1n1k,2)+pow(dImQ1n1k,2)-(*fSMpk)(0,2))/dM11; 
   
   // average weighted correlation <w1 w2 cos(n*(phi1-phi2))> for single event: 
   fQCorrelationsEBE[1]->SetBinContent(1,two1n1nW1W1);

   // average weighted correlation <w1 w2 cos(n*(phi1-phi2))> for all events:
   //fQCorrelationsW->Fill(0.,two1n1nW1W1,dM11);
   fQCorrelations[1][0]->Fill(0.5,two1n1nW1W1,dM11);
  }
  if(dM22)
  {
   two2n2nW2W2 = (pow(dReQ2n2k,2)+pow(dImQ2n2k,2)-(*fSMpk)(0,4))/dM22; 
   //fQCorrelationsW->Fill(1.,two2n2nW2W2,dM22); 
   fQCorrelations[1][0]->Fill(1.5,two2n2nW2W2,dM22);
  }
  if(dM33)
  {
   two3n3nW3W3 = (pow(dReQ3n3k,2)+pow(dImQ3n3k,2)-(*fSMpk)(0,6))/dM33;
   //fQCorrelationsW->Fill(2.,two3n3nW3W3,dM33);
   fQCorrelations[1][0]->Fill(2.5,two3n3nW3W3,dM33);
  }
  if(dM44)
  {
   two4n4nW4W4 = (pow(dReQ4n4k,2)+pow(dImQ4n4k,2)-(*fSMpk)(0,8))/dM44; 
   //fQCorrelationsW->Fill(3.,two4n4nW4W4,dM44); 
   fQCorrelations[1][0]->Fill(3.5,two4n4nW4W4,dM44); 
  } 
  if(dM31)
  {
   two1n1nW3W1 = (dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k-(*fSMpk)(0,4))/dM31; 
   //fQCorrelationsW->Fill(4.,two1n1nW3W1,dM31);  
  } 
  if(dM211)
  {
   two1n1nW1W1W2 = ((*fSMpk)(0,2)*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2)-(*fSMpk)(0,2))
                 - 2.*(dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k
                 - (*fSMpk)(0,4)))/dM211;
   //fQCorrelationsW->Fill(5.,two1n1nW1W1W2,dM211);  
  }  
 } // end of if(dMult>1)
 //..............................................................................................
 
 //..............................................................................................
 // weighted 3-particle correlations:
 Double_t three2n1n1nW2W1W1 = 0.; // <w1^2 w2 w3 cos(n*(2phi1-phi2-phi3))>
 
 if(dMult>2) 
 { 
  if(dM211)
  {                                                       
   three2n1n1nW2W1W1 = (pow(dReQ1n1k,2.)*dReQ2n2k+2.*dReQ1n1k*dImQ1n1k*dImQ2n2k-pow(dImQ1n1k,2.)*dReQ2n2k
                     - 2.*(dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k)
                     - pow(dReQ2n2k,2)-pow(dImQ2n2k,2)
                     + 2.*(*fSMpk)(0,4))/dM211;                                                                               
   //fQCorrelationsW->Fill(20.,three2n1n1nW2W1W1,dM211);
  } 
 } // end of if(dMult>2) 
 //..............................................................................................
 
 //..............................................................................................
 // weighted 4-particle correlations:
 Double_t four1n1n1n1nW1W1W1W1 = 0.; // <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))>
 if(dMult>3) 
 { 
  if(dM1111)
  {      
   four1n1n1n1nW1W1W1W1 = (pow(pow(dReQ1n1k,2.)+pow(dImQ1n1k,2.),2)
                        - 2.*(pow(dReQ1n1k,2.)*dReQ2n2k+2.*dReQ1n1k*dImQ1n1k*dImQ2n2k-pow(dImQ1n1k,2.)*dReQ2n2k)
                        + 8.*(dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k)
                        + (pow(dReQ2n2k,2)+pow(dImQ2n2k,2))
                        - 4.*(*fSMpk)(0,2)*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2))
                        - 6.*(*fSMpk)(0,4)+2.*(*fSMpk)(1,2))/dM1111;  
                        
   // average weighted correlation <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> for single event: 
   fQCorrelationsEBE[1]->SetBinContent(11,four1n1n1n1nW1W1W1W1);

   // average weighted correlation <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> for all events:                        
   //fQCorrelationsW->Fill(40.,four1n1n1n1nW1W1W1W1,dM1111);
   
   fQCorrelations[1][0]->Fill(10.5,four1n1n1n1nW1W1W1W1,dM1111);
  } 
 } // end of if(dMult>3) 
 //..............................................................................................
 
} // end of AliFlowAnalysisWithQCumulants::CalculateWeightedCorrelationsForIntegratedFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateWeightedQProductsForIntFlow() // to be improved (completed)
{
 // calculate averages like <<2><4>>, <<2><6>>, <<4><6>>, etc. which are needed to calculate covariances 
 // Remark: here we take weighted correlations!
 
 // binning of fQProductsW is organized as follows:
 // 
 // 1st bin: <2><4> 
 // 2nd bin: <2><6>
 // 3rd bin: <2><8>
 // 4th bin: <4><6>
 // 5th bin: <4><8>
 // 6th bin: <6><8>
 
 Double_t dMult = (*fSMpk)(0,0); // multiplicity (number of particles used to determine the reaction plane)

 Double_t dM11 = (*fSMpk)(1,1)-(*fSMpk)(0,2); // dM11 = sum_{i,j=1,i!=j}^M w_i w_j
 Double_t dM1111 = (*fSMpk)(3,1)-6.*(*fSMpk)(0,2)*(*fSMpk)(1,1)  
                 + 8.*(*fSMpk)(0,3)*(*fSMpk)(0,1)
                 + 3.*(*fSMpk)(1,2)-6.*(*fSMpk)(0,4); // dM1111 = sum_{i,j,k,l=1,i!=j!=k!=l}^M w_i w_j w_k w_l

 Double_t twoEBEW = 0.; // <2>
 Double_t fourEBEW = 0.; // <4>
 
 twoEBEW = fQCorrelationsEBE[1]->GetBinContent(1);
 fourEBEW = fQCorrelationsEBE[1]->GetBinContent(11);
 
 // <2><4>
 if(dMult>3)
 {
  fQProducts[1][0]->Fill(0.5,twoEBEW*fourEBEW,dM11*dM1111);
 }
 
} // end of AliFlowAnalysisWithQCumulants::CalculateWeightedQProductsForIntFlow()  


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::InitializeArraysForIntFlow()
{
 // initialize all arrays needed to calculate the integrated flow
 
 for(Int_t pW=0;pW<2;pW++) // particle weights not used (0) or used (1)
 {
  fQCorrelationsEBE[pW] = NULL;
  for(Int_t sc=0;sc<2;sc++)
  {
   fQCorrectionsEBE[pW][sc] = NULL;
  }
  for(Int_t eW=0;eW<2;eW++)
  {
   // profiles:
   fQCorrelations[pW][eW] = NULL;
   fQProducts[pW][eW] = NULL;
   for(Int_t sc=0;sc<2;sc++) // sin or cos terms
   {
    fQCorrections[pW][eW][sc] = NULL;
   }
   // histograms with results:
   fCorrelations[pW][eW] = NULL;
   fCovariances[pW][eW] = NULL;
   fCorrections[pW][eW] = NULL;
   fProductOfEventWeights[pW][eW] = NULL;
   for(Int_t power=0;power<2;power++) 
   {
    fSumOfEventWeights[pW][eW][power] = NULL;    
   }
   for(Int_t nua=0;nua<2;nua++) // not corrected or corrected
   {
    fCumulants[pW][eW][nua] = NULL;
    fIntFlow[pW][eW][nua] = NULL;
   }
  } 
 }
 
} // end of void AliFlowAnalysisWithQCumulants::InitializeArraysForIntFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::InitializeArraysForDiffFlow()
{
 // initialize all arrays needed to calcualted differential flow
 
 // nested lists in fDiffFlowProfiles:
 for(Int_t t=0;t<2;t++)
 {
  fDFPType[t] = NULL;
  for(Int_t pW=0;pW<2;pW++) // particle weights not used (0) or used (1)
  {
   fDFPParticleWeights[t][pW] = NULL;
   for(Int_t eW=0;eW<2;eW++)
   {   
    fDFPEventWeights[t][pW][eW] = NULL;
    fDiffFlowCorrelations[t][pW][eW] = NULL;
    fDiffFlowProductsOfCorrelations[t][pW][eW] = NULL;
    for(Int_t sc=0;sc<2;sc++)
    {
     fDiffFlowCorrectionTerms[t][pW][eW][sc] = NULL;
    }
   } 
  }
 }  
 
 // profiles in nested lists in fDiffFlowProfiles:
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pW=0;pW<2;pW++) // particle weights not used (0) or used (1)
  {
   for(Int_t eW=0;eW<2;eW++)
   {
    // correlations:
    for(Int_t correlationIndex=0;correlationIndex<4;correlationIndex++)
    {
     fCorrelationsPro[t][pW][eW][correlationIndex] = NULL;
    }
    // products of correlations:
    for(Int_t productOfCorrelationsIndex=0;productOfCorrelationsIndex<6;productOfCorrelationsIndex++)
    {
     fProductsOfCorrelationsPro[t][pW][eW][productOfCorrelationsIndex] = NULL;
    }
    // correction terms:
    for(Int_t sc=0;sc<2;sc++)
    {
     for(Int_t correctionsIndex=0;correctionsIndex<2;correctionsIndex++)
     {
      fCorrectionTermsPro[t][pW][eW][sc][correctionsIndex] = NULL;
     } 
    } 
   }
  }
 }  
 
 // nested lists in fDiffFlowResults:
 for(Int_t t=0;t<2;t++)
 {
  fDFRType[t] = NULL;
  for(Int_t pW=0;pW<2;pW++) // particle weights not used (0) or used (1)
  {
   fDFRParticleWeights[t][pW] = NULL;
   for(Int_t eW=0;eW<2;eW++)
   {    
    fDFREventWeights[t][pW][eW] = NULL;
    fDiffFlowFinalCorrelations[t][pW][eW] = NULL;
    fDiffFlowFinalCorrections[t][pW][eW] = NULL;
    fDiffFlowFinalCovariances[t][pW][eW] = NULL;
    for(Int_t nua=0;nua<2;nua++)
    {
     fDFRCorrections[t][pW][eW][nua] = NULL;   
     fDiffFlowFinalCumulants[t][pW][eW][nua] = NULL;   
     fDiffFlowFinalFlow[t][pW][eW][nua] = NULL;
    }
   } 
  }
 }  
 
 // 2D and 1D histograms in nested lists in fDiffFlowResults:
 for(Int_t t=0;t<2;t++) 
 {
  fNonEmptyBins2D[t] = NULL;
  for(Int_t pe=0;pe<2;pe++)
  {
   fNonEmptyBins1D[t][pe] = NULL;
  }
  for(Int_t pW=0;pW<2;pW++) // particle weights not used (0) or used (1)
  {
   for(Int_t eW=0;eW<2;eW++)
   {
    // correlations:
    for(Int_t correlationIndex=0;correlationIndex<4;correlationIndex++) 
    {
     fFinalCorrelations2D[t][pW][eW][correlationIndex] = NULL;
     for(Int_t pe=0;pe<2;pe++)
     {
      fFinalCorrelations1D[t][pW][eW][pe][correlationIndex] = NULL;   
     }
    }
    // corrections:
    for(Int_t correctionIndex=0;correctionIndex<4;correctionIndex++) 
    {
     fFinalCorrections2D[t][pW][eW][correctionIndex] = NULL;
     for(Int_t pe=0;pe<2;pe++)
     {
      fFinalCorrections1D[t][pW][eW][pe][correctionIndex] = NULL;     
     }
    }
    // covariances:
    for(Int_t covarianceIndex=0;covarianceIndex<4;covarianceIndex++) 
    {
     fFinalCovariances2D[t][pW][eW][covarianceIndex] = NULL;
     for(Int_t pe=0;pe<2;pe++)
     {
      fFinalCovariances1D[t][pW][eW][pe][covarianceIndex] = NULL;     
     }
    }
    for(Int_t nua=0;nua<2;nua++) 
    {  
     // cumulants:
     for(Int_t cumulantIndex=0;cumulantIndex<4;cumulantIndex++) 
     {
      fFinalCumulants2D[t][pW][eW][nua][cumulantIndex] = NULL;
      fFinalCumulantsPt[t][pW][eW][nua][cumulantIndex] = NULL;     
      fFinalCumulantsEta[t][pW][eW][nua][cumulantIndex] = NULL;       
     }
     // final flow:
     for(Int_t finalFlowIndex=0;finalFlowIndex<4;finalFlowIndex++) 
     {
      fFinalFlow2D[t][pW][eW][nua][finalFlowIndex] = NULL;
      fFinalFlowPt[t][pW][eW][nua][finalFlowIndex] = NULL;   
      fFinalFlowEta[t][pW][eW][nua][finalFlowIndex] = NULL;   
     }
    } 
   }  
  }
 }
 
 for(Int_t t=0;t<3;t++) // type (RP, POI, POI&&RP)
 {
  for(Int_t m=0;m<4;m++) // multiple of harmonic
  {
   for(Int_t k=0;k<9;k++) // power of weight
   {
    fReEBE[t][m][k] = NULL;
    fImEBE[t][m][k] = NULL;
    fs[t][k] = NULL; // to be improved (this doesn't need to be within loop over m)
   }   
  }
 }
  
} // end of AliFlowAnalysisWithQCumulants::InitializeArraysForDiffFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::BookEverythingForDifferentialFlow()
{
 // organize and book everything needed for differential flow
 
 // book and nest all lists belonging to the list fDiffFlowProfiles "Profiles":
 TList list; // to be improved (do I need this here?)
 list.SetOwner(kTRUE); // to be improved (do I need this here?)
 TString typeFlag[2] = {"RP","POI"}; // to be improved (do I need this here?)
 TString pWeightsFlag[2] = {"not used","used"};  // to be improved (do I need this here?)
 TString eWeightsFlag[2] = {"exact","non-exact"};  // to be improved (do I need this here?)
 //TString sinCosFlag[2] = {"sin","cos"}; // to be improved (do I need this here?)
 TString nuaFlag[2] = {"not corrected","corrected"}; // nua = non-uniform acceptance // to be improved (do I need this here?)
 
 // book all 2D profiles and add to appropriate lists:
 TProfile2D styleProfile("","",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 styleProfile.SetXTitle("p_{t}");
 styleProfile.SetYTitle("#eta");
 TString correlationName[4] = {"<2'>","<4'>","<6'>","<8'>"};
 TString cumulantName[4] = {"QC{2'}","QC{4'}","QC{6'}","QC{8'}"};
 TString productOfCorrelationsName[5] = {"<2><2'>","<2><4'>","<4><2'>","<4><4'>","<2'><4'>"};
 TString correctionsSinTermsName[2] = {"sin(1)","sin(2)"};
 TString correctionsCosTermsName[2] = {"cos(1)","cos(2)"};
 TString correctionName[4] = {"corr. to QC{2'}","corr. to QC{4'}","corr. to QC{6'}","corr. to QC{8'}"};
 TString covarianceName[4] = {"Cov(2,2')","Cov(2,4')","Cov(2',4')","Cov(4,2')"}; // to be improved (reorganized)
 TString flowName[4] = {"v'{2}","v'{4}","v'{6}","v'{8}"}; 
 
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++) // weights: not used or used 
  {
   for(Int_t eW=0;eW<2;eW++)
   {
    // correlations:
    for(Int_t correlationIndex=0;correlationIndex<4;correlationIndex++)
    {
     fCorrelationsPro[t][pW][eW][correlationIndex] = (TProfile2D*)styleProfile.Clone(correlationName[correlationIndex].Data());
     fCorrelationsPro[t][pW][eW][correlationIndex]->SetTitle(correlationName[correlationIndex].Data()); 
     fDiffFlowCorrelations[t][pW][eW]->Add(fCorrelationsPro[t][pW][eW][correlationIndex]);
    }
    // products of correlations:
    for(Int_t productOfCorrelationsIndex=0;productOfCorrelationsIndex<5;productOfCorrelationsIndex++)
    {
     fProductsOfCorrelationsPro[t][pW][eW][productOfCorrelationsIndex] = (TProfile2D*)styleProfile.Clone(productOfCorrelationsName[productOfCorrelationsIndex].Data());
     fProductsOfCorrelationsPro[t][pW][eW][productOfCorrelationsIndex]->SetTitle(productOfCorrelationsName[productOfCorrelationsIndex].Data()); 
     fDiffFlowProductsOfCorrelations[t][pW][eW]->Add(fProductsOfCorrelationsPro[t][pW][eW][productOfCorrelationsIndex]);
    }
    // correction terms:
    for(Int_t sc=0;sc<2;sc++)
    {
     for(Int_t correctionsIndex=0;correctionsIndex<2;correctionsIndex++)
     {
      if(sc==0)
      {
       fCorrectionTermsPro[t][pW][eW][sc][correctionsIndex] = (TProfile2D*)styleProfile.Clone(correctionsSinTermsName[correctionsIndex].Data());
       fCorrectionTermsPro[t][pW][eW][sc][correctionsIndex]->SetTitle(correctionsSinTermsName[correctionsIndex].Data());
       fDiffFlowCorrectionTerms[t][pW][eW][sc]->Add(fCorrectionTermsPro[t][pW][eW][sc][correctionsIndex]);
      }
      if(sc==1)
      {
       fCorrectionTermsPro[t][pW][eW][sc][correctionsIndex] = (TProfile2D*)styleProfile.Clone(correctionsCosTermsName[correctionsIndex].Data());
       fCorrectionTermsPro[t][pW][eW][sc][correctionsIndex]->SetTitle(correctionsCosTermsName[correctionsIndex].Data());
       fDiffFlowCorrectionTerms[t][pW][eW][sc]->Add(fCorrectionTermsPro[t][pW][eW][sc][correctionsIndex]);
      }    
     }
    }
   }   
  }  
 }
   
 // book all 2D and 1D histograms for final results and store them in appropriate lists:
 TH2D styleHistPtEta("","",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 styleHistPtEta.SetXTitle("p_{t}");
 styleHistPtEta.SetYTitle("#eta");
 TH1D styleHistPt("","",fnBinsPt,fPtMin,fPtMax);
 styleHistPt.SetXTitle("p_{t}");
 TH1D styleHistEta("","",fnBinsEta,fEtaMin,fEtaMax);
 styleHistEta.SetXTitle("#eta");
 for(Int_t t=0;t<2;t++) // RP or POI
 {
  // 2D:
  fNonEmptyBins2D[t] = (TH2D*)styleHistPtEta.Clone(Form("%s, (p_{t},#eta)",typeFlag[t].Data()));
  fNonEmptyBins2D[t]->SetTitle(Form("Non-empty bins: %s (p_{t},#eta)",typeFlag[t].Data())); 
  fDFRType[t]->Add(fNonEmptyBins2D[t]); 
  // 1D:
  fNonEmptyBins1D[t][0] = (TH1D*)styleHistPt.Clone(Form("%s, p_{t}",typeFlag[t].Data()));
  fNonEmptyBins1D[t][0]->SetTitle(Form("Non-empty bins: %s (p_{t})",typeFlag[t].Data())); 
  fDFRType[t]->Add(fNonEmptyBins1D[t][0]); 
  fNonEmptyBins1D[t][1] = (TH1D*)styleHistEta.Clone(Form("%s, #eta",typeFlag[t].Data()));
  fNonEmptyBins1D[t][1]->SetTitle(Form("Non-empty bins: %s (#eta)",typeFlag[t].Data())); 
  fDFRType[t]->Add(fNonEmptyBins1D[t][1]); 
  
  for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++) // particle weights not used or used
  {
   for(Int_t eW=0;eW<2;eW++)
   {
    // correlations:
    for(Int_t correlationIndex=0;correlationIndex<4;correlationIndex++)
    {
     // 2D:
     fFinalCorrelations2D[t][pW][eW][correlationIndex] = (TH2D*)styleHistPtEta.Clone(Form("%s, (p_{t},#eta)",correlationName[correlationIndex].Data()));
     fFinalCorrelations2D[t][pW][eW][correlationIndex]->SetTitle(Form("%s (p_{t},#eta): %s, pWeights %s, eWeights %s",correlationName[correlationIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())); 
     fDiffFlowFinalCorrelations[t][pW][eW]->Add(fFinalCorrelations2D[t][pW][eW][correlationIndex]);   
     // 1D:
     for(Int_t pe=0;pe<2;pe++) // pt or eta:
     {
      if(pe==0) 
      { 
       fFinalCorrelations1D[t][pW][eW][pe][correlationIndex] = (TH1D*)styleHistPt.Clone(Form("%s, (p_{t})",correlationName[correlationIndex].Data()));
       fFinalCorrelations1D[t][pW][eW][pe][correlationIndex]->SetTitle(Form("%s (p_{t}): %s, pWeights %s, eWeights %s",correlationName[correlationIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())); 
       fDiffFlowFinalCorrelations[t][pW][eW]->Add(fFinalCorrelations1D[t][pW][eW][pe][correlationIndex]);   
      }
      if(pe==1)
      {
       fFinalCorrelations1D[t][pW][eW][pe][correlationIndex] = (TH1D*)styleHistEta.Clone(Form("%s, (#eta)",correlationName[correlationIndex].Data()));
       fFinalCorrelations1D[t][pW][eW][pe][correlationIndex]->SetTitle(Form("%s (#eta): %s, pWeights %s, eWeights %s",correlationName[correlationIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())); 
       fDiffFlowFinalCorrelations[t][pW][eW]->Add(fFinalCorrelations1D[t][pW][eW][pe][correlationIndex]); 
      } 
     } 
    }    
    // corrections:
    for(Int_t correctionIndex=0;correctionIndex<4;correctionIndex++)
    {
     // 2D:
     fFinalCorrections2D[t][pW][eW][correctionIndex] = (TH2D*)styleHistPtEta.Clone(Form("%s, (p_{t},#eta)",correctionName[correctionIndex].Data()));
     fFinalCorrections2D[t][pW][eW][correctionIndex]->SetTitle(Form("%s (p_{t},#eta): %s, pWeights %s, eWeights %s",correctionName[correctionIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())); 
     fDiffFlowFinalCorrections[t][pW][eW]->Add(fFinalCorrections2D[t][pW][eW][correctionIndex]);   
     // 1D:
     for(Int_t pe=0;pe<2;pe++) // pt or eta:
     {
      if(pe==0) 
      { 
       fFinalCorrections1D[t][pW][eW][pe][correctionIndex] = (TH1D*)styleHistPt.Clone(Form("%s, (p_{t})",correctionName[correctionIndex].Data()));
       fFinalCorrections1D[t][pW][eW][pe][correctionIndex]->SetTitle(Form("%s (p_{t}): %s, pWeights %s, eWeights %s",correctionName[correctionIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())); 
       fDiffFlowFinalCorrections[t][pW][eW]->Add(fFinalCorrections1D[t][pW][eW][pe][correctionIndex]);   
      }
      if(pe==1)
      {
       fFinalCorrections1D[t][pW][eW][pe][correctionIndex] = (TH1D*)styleHistEta.Clone(Form("%s, (#eta)",correctionName[correctionIndex].Data()));
       fFinalCorrections1D[t][pW][eW][pe][correctionIndex]->SetTitle(Form("%s (#eta): %s, pWeights %s, eWeights %s",correctionName[correctionIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())); 
       fDiffFlowFinalCorrections[t][pW][eW]->Add(fFinalCorrections1D[t][pW][eW][pe][correctionIndex]); 
      } 
     } 
    }
    // covariances:
    for(Int_t covarianceIndex=0;covarianceIndex<4;covarianceIndex++)
    {
     // 2D:
     fFinalCovariances2D[t][pW][eW][covarianceIndex] = (TH2D*)styleHistPtEta.Clone(Form("%s, (p_{t},#eta)",covarianceName[covarianceIndex].Data()));
     fFinalCovariances2D[t][pW][eW][covarianceIndex]->SetTitle(Form("%s (p_{t},#eta): %s, pWeights %s, eWeights %s",covarianceName[covarianceIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())); 
     fDiffFlowFinalCovariances[t][pW][eW]->Add(fFinalCovariances2D[t][pW][eW][covarianceIndex]);  
     // 1D:
     for(Int_t pe=0;pe<2;pe++) // pt or eta:
     {
      if(pe==0) 
      { 
       fFinalCovariances1D[t][pW][eW][pe][covarianceIndex] = (TH1D*)styleHistPt.Clone(Form("%s, (p_{t})",covarianceName[covarianceIndex].Data()));
       fFinalCovariances1D[t][pW][eW][pe][covarianceIndex]->SetTitle(Form("%s (p_{t}): %s, pWeights %s, eWeights %s",covarianceName[covarianceIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())); 
       fDiffFlowFinalCovariances[t][pW][eW]->Add(fFinalCovariances1D[t][pW][eW][pe][covarianceIndex]);   
      }
      if(pe==1)
      {
       fFinalCovariances1D[t][pW][eW][pe][covarianceIndex] = (TH1D*)styleHistEta.Clone(Form("%s, (#eta)",covarianceName[covarianceIndex].Data()));
       fFinalCovariances1D[t][pW][eW][pe][covarianceIndex]->SetTitle(Form("%s (#eta): %s, pWeights %s, eWeights %s",covarianceName[covarianceIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data())); 
       fDiffFlowFinalCovariances[t][pW][eW]->Add(fFinalCovariances1D[t][pW][eW][pe][covarianceIndex]); 
      } 
     } 
    }
    for(Int_t nua=0;nua<2;nua++) // corrected or not 
    {
     // cumulants:
     for(Int_t cumulantIndex=0;cumulantIndex<4;cumulantIndex++)
     {
      // 2D:
      fFinalCumulants2D[t][pW][eW][nua][cumulantIndex] = (TH2D*)styleHistPtEta.Clone(Form("%s, (p_{t},#eta)",cumulantName[cumulantIndex].Data()));
      fFinalCumulants2D[t][pW][eW][nua][cumulantIndex]->SetTitle(Form("%s (p_{t},#eta): %s, pWeights %s, eWeights %s, %s for NUA",cumulantName[cumulantIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())); 
      fDiffFlowFinalCumulants[t][pW][eW][nua]->Add(fFinalCumulants2D[t][pW][eW][nua][cumulantIndex]);   
      // pt:
      fFinalCumulantsPt[t][pW][eW][nua][cumulantIndex] = (TH1D*)styleHistPt.Clone(Form("%s, (p_{t})",cumulantName[cumulantIndex].Data()));
      fFinalCumulantsPt[t][pW][eW][nua][cumulantIndex]->SetTitle(Form("%s (p_{t}): %s, pWeights %s, eWeights %s, %s for NUA",cumulantName[cumulantIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())); 
      fDiffFlowFinalCumulants[t][pW][eW][nua]->Add(fFinalCumulantsPt[t][pW][eW][nua][cumulantIndex]);   
      // eta:     
      fFinalCumulantsEta[t][pW][eW][nua][cumulantIndex] = (TH1D*)styleHistEta.Clone(Form("%s, (#eta)",cumulantName[cumulantIndex].Data()));
      fFinalCumulantsEta[t][pW][eW][nua][cumulantIndex]->SetTitle(Form("%s (#eta): %s, pWeights %s, eWeights %s, %s for NUA",cumulantName[cumulantIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())); 
      fDiffFlowFinalCumulants[t][pW][eW][nua]->Add(fFinalCumulantsEta[t][pW][eW][nua][cumulantIndex]);  
     }       
     // flow:
     for(Int_t flowIndex=0;flowIndex<4;flowIndex++)
     {
      // 2D:
      fFinalFlow2D[t][pW][eW][nua][flowIndex] = (TH2D*)styleHistPtEta.Clone(Form("%s, (p_{t},#eta)",flowName[flowIndex].Data()));
      fFinalFlow2D[t][pW][eW][nua][flowIndex]->SetTitle(Form("%s (p_{t},#eta): %s, pWeights %s, eWeights %s, %s for NUA",flowName[flowIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())); 
      fDiffFlowFinalFlow[t][pW][eW][nua]->Add(fFinalFlow2D[t][pW][eW][nua][flowIndex]);   
      // pt: 
      fFinalFlowPt[t][pW][eW][nua][flowIndex] = (TH1D*)styleHistPt.Clone(Form("%s, (p_{t})",flowName[flowIndex].Data()));
      fFinalFlowPt[t][pW][eW][nua][flowIndex]->SetTitle(Form("%s (p_{t}): %s, pWeights %s, eWeights %s, %s for NUA",flowName[flowIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())); 
      fDiffFlowFinalFlow[t][pW][eW][nua]->Add(fFinalFlowPt[t][pW][eW][nua][flowIndex]);   
      // eta:
      fFinalFlowEta[t][pW][eW][nua][flowIndex] = (TH1D*)styleHistEta.Clone(Form("%s, (#eta)",flowName[flowIndex].Data()));
      fFinalFlowEta[t][pW][eW][nua][flowIndex]->SetTitle(Form("%s (#eta): %s, pWeights %s, eWeights %s, %s for NUA",flowName[flowIndex].Data(),typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data())); 
      fDiffFlowFinalFlow[t][pW][eW][nua]->Add(fFinalFlowEta[t][pW][eW][nua][flowIndex]);   
     }
    }   
   } 
  } 
 } 
   
 // Event-by-event r_{m*n,k}(pt,eta), p_{m*n,k}(pt,eta) and q_{m*n,k}(pt,eta)
 // Explanantion of notation:
 //  1.) n is harmonic, m is multiple of harmonic;
 //  2.) k is power of weight;
 //  3.) r_{m*n,k}(pt,eta) = Q-vector evaluated in harmonic m*n for RPs in particular (pt,eta) bin (i-th RP is weighted with w_i^k);   
 //  4.) p_{m*n,k}(pt,eta) = Q-vector evaluated in harmonic m*n for POIs in particular (pt,eta) bin 
 //                          (if i-th POI is also RP, than it is weighted with w_i^k);   
 //  5.) q_{m*n,k}(pt,eta) = Q-vector evaluated in harmonic m*n for particles which are both RPs and POIs in particular (pt,eta) bin 
 //                          (i-th RP&&POI is weighted with w_i^k)            
  
 TProfile2D styleRe("typeMultiplePowerRe","typeMultiplePowerRe",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 TProfile2D styleIm("typeMultiplePowerIm","typeMultiplePowerIm",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);
 
 for(Int_t t=0;t<3;t++) // typeFlag (0 = RP, 1 = POI, 2 = RP&&POI )
 { 
  for(Int_t m=0;m<4;m++)
  {
   for(Int_t k=0;k<9;k++)
   {
    fReEBE[t][m][k] = (TProfile2D*)styleRe.Clone(Form("typeFlag%dmultiple%dpower%dRe",t,m,k)); 
    fImEBE[t][m][k] = (TProfile2D*)styleIm.Clone(Form("typeFlag%dmultiple%dpower%dIm",t,m,k));
   }
  } 
 } 

 TProfile2D styleS("typePower","typePower",fnBinsPt,fPtMin,fPtMax,fnBinsEta,fEtaMin,fEtaMax);

 for(Int_t t=0;t<3;t++) // typeFlag (0 = RP, 1 = POI, 2 = RP&&POI )
 { 
  for(Int_t k=0;k<9;k++)
  {
   fs[t][k] = (TProfile2D*)styleS.Clone(Form("typeFlag%dpower%d",t,k));
  }
 }
   
} // end of AliFlowAnalysisWithQCumulants::BookEverythingForDifferentialFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCorrelationsForDifferentialFlow(TString type)
{
 // calculate all reduced correlations needed for differential flow for each (pt,eta) bin:
 
 Int_t typeFlag = -1; 
  
 // reduced correlations ares stored in fCorrelationsPro[t][pW][index] and are indexed as follows:
 // index:
 // 0: <2'>
 // 1: <4'>

 // multiplicity:
 Double_t dMult = (*fSMpk)(0,0);
 
 // real and imaginary parts of non-weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n = (*fReQ)(0,0);
 Double_t dReQ2n = (*fReQ)(1,0);
 //Double_t dReQ3n = (*fReQ)(2,0);
 //Double_t dReQ4n = (*fReQ)(3,0);
 Double_t dImQ1n = (*fImQ)(0,0);
 Double_t dImQ2n = (*fImQ)(1,0);
 //Double_t dImQ3n = (*fImQ)(2,0);
 //Double_t dImQ4n = (*fImQ)(3,0);

 // looping over all (pt,eta) bins and calculating correlations needed for differential flow: 
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   // real and imaginary parts of p_{m*n,0} (non-weighted Q-vector evaluated for POIs in particular (pt,eta) bin): 
   Double_t p1n0kRe = 0.;
   Double_t p1n0kIm = 0.;

   // number of POIs in particular (pt,eta) bin:
   Double_t mp = 0.;

   // real and imaginary parts of q_{m*n,0} (non-weighted Q-vector evaluated for particles which are both RPs and POIs in particular (pt,eta) bin):
   Double_t q1n0kRe = 0.;
   Double_t q1n0kIm = 0.;
   Double_t q2n0kRe = 0.;
   Double_t q2n0kIm = 0.;

   // number of particles which are both RPs and POIs in particular (pt,eta) bin:
   Double_t mq = 0.;
   
   // q_{m*n,0}:
   q1n0kRe = fReEBE[2][0][0]->GetBinContent(fReEBE[2][0][0]->GetBin(p,e))
           * fReEBE[2][0][0]->GetBinEntries(fReEBE[2][0][0]->GetBin(p,e));
   q1n0kIm = fImEBE[2][0][0]->GetBinContent(fImEBE[2][0][0]->GetBin(p,e))
           * fImEBE[2][0][0]->GetBinEntries(fImEBE[2][0][0]->GetBin(p,e));
   q2n0kRe = fReEBE[2][1][0]->GetBinContent(fReEBE[2][1][0]->GetBin(p,e))
           * fReEBE[2][1][0]->GetBinEntries(fReEBE[2][1][0]->GetBin(p,e));
   q2n0kIm = fImEBE[2][1][0]->GetBinContent(fImEBE[2][1][0]->GetBin(p,e))
           * fImEBE[2][1][0]->GetBinEntries(fImEBE[2][1][0]->GetBin(p,e));
           
   mq = fReEBE[2][0][0]->GetBinEntries(fReEBE[2][0][0]->GetBin(p,e)); // to be improved (cross-checked by accessing other profiles here)
   
   if(type == "POI")
   {
    // p_{m*n,0}:
    p1n0kRe = fReEBE[1][0][0]->GetBinContent(fReEBE[1][0][0]->GetBin(p,e))
            * fReEBE[1][0][0]->GetBinEntries(fReEBE[1][0][0]->GetBin(p,e));
    p1n0kIm = fImEBE[1][0][0]->GetBinContent(fImEBE[1][0][0]->GetBin(p,e))  
            * fImEBE[1][0][0]->GetBinEntries(fImEBE[1][0][0]->GetBin(p,e));
            
    mp = fReEBE[1][0][0]->GetBinEntries(fReEBE[1][0][0]->GetBin(p,e)); // to be improved (cross-checked by accessing other profiles here)
    
    typeFlag = 1;
   }
   else if(type == "RP")
   {
    // p_{m*n,0} = q_{m*n,0}:
    p1n0kRe = q1n0kRe; 
    p1n0kIm = q1n0kIm; 
    mp = mq; 
    
    typeFlag = 0;
   }
   
   // count events with non-empty (pt,eta) bin:
   if(mp>0)
   {
    fNonEmptyBins2D[typeFlag]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,1);
   }
   
   // 2'-particle correlation for particular (pt,eta) bin:
   Double_t two1n1nPtEta = 0.;
   if(mp*dMult-mq)
   {
    two1n1nPtEta = (p1n0kRe*dReQ1n+p1n0kIm*dImQ1n-mq)
                 / (mp*dMult-mq);
   
    // fill the 2D profile to get the average correlation for each (pt,eta) bin:
    if(type == "POI")
    { 
     //f2pPtEtaPOI->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,two1n1nPtEta,mp*dMult-mq);
     
     fCorrelationsPro[1][0][0][0]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,two1n1nPtEta,mp*dMult-mq);
    }
    else if(type == "RP")
    {
     //f2pPtEtaRP->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,two1n1nPtEta,mp*dMult-mq);   
     fCorrelationsPro[0][0][0][0]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,two1n1nPtEta,mp*dMult-mq);
    }
   } // end of if(mp*dMult-mq)
  
   // 4'-particle correlation:
   Double_t four1n1n1n1nPtEta = 0.;
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
    
    // fill the 2D profile to get the average correlation for each (pt, eta) bin:
    if(type == "POI")
    {
     //f4pPtEtaPOI->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nPtEta,
     //                  (mp-mq)*dMult*(dMult-1.)*(dMult-2.)
     //                   + mq*(dMult-1.)*(dMult-2.)*(dMult-3.));
     
     fCorrelationsPro[1][0][0][1]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nPtEta,
                                     (mp-mq)*dMult*(dMult-1.)*(dMult-2.)
                                     + mq*(dMult-1.)*(dMult-2.)*(dMult-3.));
    }
    else if(type == "RP")
    {
     //f4pPtEtaRP->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nPtEta,
     //                 (mp-mq)*dMult*(dMult-1.)*(dMult-2.)
     //                  + mq*(dMult-1.)*(dMult-2.)*(dMult-3.));   
                       
     fCorrelationsPro[0][0][0][1]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nPtEta,
                                       (mp-mq)*dMult*(dMult-1.)*(dMult-2.)
                                       + mq*(dMult-1.)*(dMult-2.)*(dMult-3.));                   
    }
   } // end of if((mp-mq)*dMult*(dMult-1.)*(dMult-2.)
     //            +mq*(dMult-1.)*(dMult-2.)*(dMult-3.))
   
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateCorrelationsForDifferentialFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateWeightedCorrelationsForDifferentialFlow(TString type)
{
 // calculate all weighted correlations needed for differential flow 
 
 // real and imaginary parts of weighted Q-vectors evaluated in harmonics n, 2n, 3n and 4n: 
 Double_t dReQ1n1k = (*fReQ)(0,1);
 Double_t dReQ2n2k = (*fReQ)(1,2);
 Double_t dReQ1n3k = (*fReQ)(0,3);
 //Double_t dReQ4n4k = (*fReQ)(3,4);
 Double_t dImQ1n1k = (*fImQ)(0,1);
 Double_t dImQ2n2k = (*fImQ)(1,2);
 Double_t dImQ1n3k = (*fImQ)(0,3);
 //Double_t dImQ4n4k = (*fImQ)(3,4);
 
 // S^M_{p,k} (see .h file for the definition of fSMpk):
 Double_t dSM1p1k = (*fSMpk)(0,1);
 Double_t dSM1p2k = (*fSMpk)(0,2);
 Double_t dSM1p3k = (*fSMpk)(0,3);
 Double_t dSM2p1k = (*fSMpk)(1,1);
 Double_t dSM3p1k = (*fSMpk)(2,1);
 
 // looping over all (pt,eta) bins and calculating weighted correlations needed for differential flow: 
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
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
    // p_{m*n,0}:
    p1n0kRe = fReEBE[1][0][0]->GetBinContent(fReEBE[1][0][0]->GetBin(p,e))
            * fReEBE[1][0][0]->GetBinEntries(fReEBE[1][0][0]->GetBin(p,e));
    p1n0kIm = fImEBE[1][0][0]->GetBinContent(fImEBE[1][0][0]->GetBin(p,e))
            * fImEBE[1][0][0]->GetBinEntries(fImEBE[1][0][0]->GetBin(p,e)); 
            
    mp = fReEBE[1][0][0]->GetBinEntries(fReEBE[1][0][0]->GetBin(p,e));
    
    // q_{m*n,k}: 
    q1n2kRe = fReEBE[2][0][2]->GetBinContent(fReEBE[2][0][2]->GetBin(p,e))
            * fReEBE[2][0][2]->GetBinEntries(fReEBE[2][0][2]->GetBin(p,e));
    q1n2kIm = fImEBE[2][0][2]->GetBinContent(fImEBE[2][0][2]->GetBin(p,e))
            * fImEBE[2][0][2]->GetBinEntries(fImEBE[2][0][2]->GetBin(p,e));
    q2n1kRe = fReEBE[2][1][1]->GetBinContent(fReEBE[2][1][1]->GetBin(p,e))
            * fReEBE[2][1][1]->GetBinEntries(fReEBE[2][1][1]->GetBin(p,e)); 
    q2n1kIm = fImEBE[2][1][1]->GetBinContent(fImEBE[2][1][1]->GetBin(p,e))
            * fImEBE[2][1][1]->GetBinEntries(fImEBE[2][1][1]->GetBin(p,e));
       
    // s_{1,1}, s_{1,2} and s_{1,3} // to be improved (add explanation)  
    s1p1k = pow(fs[2][1]->GetBinContent(fs[2][1]->GetBin(p,e)),1.); 
    s1p2k = pow(fs[2][2]->GetBinContent(fs[2][2]->GetBin(p,e)),1.); 
    s1p3k = pow(fs[2][3]->GetBinContent(fs[2][3]->GetBin(p,e)),1.); 
   
    // M0111 from Eq. (118) in QC2c (to be improved (notation)):
    dM0111 = mp*(dSM3p1k-3.*dSM1p1k*dSM1p2k+2.*dSM1p3k)
           - 3.*(s1p1k*(dSM2p1k-dSM1p2k)
           + 2.*(s1p3k-s1p2k*dSM1p1k));
   }
   else if(type == "RP")
   {
    p1n0kRe = fReEBE[0][0][0]->GetBinContent(fReEBE[0][0][0]->GetBin(p,e))
            * fReEBE[0][0][0]->GetBinEntries(fReEBE[0][0][0]->GetBin(p,e));
    p1n0kIm = fImEBE[0][0][0]->GetBinContent(fImEBE[0][0][0]->GetBin(p,e))
            * fImEBE[0][0][0]->GetBinEntries(fImEBE[0][0][0]->GetBin(p,e));
            
    mp = fReEBE[0][0][0]->GetBinEntries(fReEBE[0][0][0]->GetBin(p,e));
    
    // q_{m*n,k}: 
    q1n2kRe = fReEBE[0][0][2]->GetBinContent(fReEBE[0][0][2]->GetBin(p,e))
            * fReEBE[0][0][2]->GetBinEntries(fReEBE[0][0][2]->GetBin(p,e));
    q1n2kIm = fImEBE[0][0][2]->GetBinContent(fImEBE[0][0][2]->GetBin(p,e))
            * fImEBE[0][0][2]->GetBinEntries(fImEBE[0][0][2]->GetBin(p,e));
    q2n1kRe = fReEBE[0][1][1]->GetBinContent(fReEBE[0][1][1]->GetBin(p,e))
            * fReEBE[0][1][1]->GetBinEntries(fReEBE[0][1][1]->GetBin(p,e));
    q2n1kIm = fImEBE[0][1][1]->GetBinContent(fImEBE[0][1][1]->GetBin(p,e))
            * fImEBE[0][1][1]->GetBinEntries(fImEBE[0][1][1]->GetBin(p,e));
   
    // s_{1,1}, s_{1,2} and s_{1,3} // to be improved (add explanation)  
    s1p1k = pow(fs[0][1]->GetBinContent(fs[0][1]->GetBin(p,e)),1.); 
    s1p2k = pow(fs[0][2]->GetBinContent(fs[0][2]->GetBin(p,e)),1.); 
    s1p3k = pow(fs[0][3]->GetBinContent(fs[0][3]->GetBin(p,e)),1.); 
   
    // M0111 from Eq. (118) in QC2c (to be improved (notation)):
    dM0111 = mp*(dSM3p1k-3.*dSM1p1k*dSM1p2k+2.*dSM1p3k)
           - 3.*(s1p1k*(dSM2p1k-dSM1p2k)
           + 2.*(s1p3k-s1p2k*dSM1p1k));
    //...............................................................................................   
   }
   
   // 2'-particle correlation:
   Double_t two1n1nW0W1PtEta = 0.;
   if(mp*dSM1p1k-s1p1k)
   {
    two1n1nW0W1PtEta = (p1n0kRe*dReQ1n1k+p1n0kIm*dImQ1n1k-s1p1k)
                 / (mp*dSM1p1k-s1p1k);
   
    // fill the 2D profile to get the average correlation for each (pt, eta) bin:
    if(type == "POI")
    {
     //f2pPtEtaPOIW->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,two1n1nW0W1PtEta,
     //                   mp*dSM1p1k-s1p1k);
     fCorrelationsPro[1][1][0][0]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,two1n1nW0W1PtEta,mp*dSM1p1k-s1p1k);
    }
    else if(type == "RP")
    {
     //f2pPtEtaRPW->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,two1n1nW0W1PtEta,
     //                  mp*dSM1p1k-s1p1k); 
     fCorrelationsPro[0][1][0][0]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,two1n1nW0W1PtEta,mp*dSM1p1k-s1p1k);  
    }
   } // end of if(mp*dMult-dmPrimePrimePtEta)
   
   // 4'-particle correlation:
   Double_t four1n1n1n1nW0W1W1W1PtEta = 0.;
   if(dM0111)
   {
    four1n1n1n1nW0W1W1W1PtEta = ((pow(dReQ1n1k,2.)+pow(dImQ1n1k,2.))*(p1n0kRe*dReQ1n1k+p1n0kIm*dImQ1n1k)
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
                      / dM0111; // to be imropoved (notation of dM0111)
   
    // fill the 2D profile to get the average correlation for each (pt, eta) bin:
    if(type == "POI")
    {
     //f4pPtEtaPOIW->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nW0W1W1W1PtEta,dM0111);
     fCorrelationsPro[1][1][0][1]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nW0W1W1W1PtEta,dM0111);
    }
    else if(type == "RP")
    {
     //f4pPtEtaRPW->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nW0W1W1W1PtEta,dM0111); 
     fCorrelationsPro[0][1][0][1]->Fill(fPtMin+(p-1)*fPtBinWidth,fEtaMin+(e-1)*fEtaBinWidth,four1n1n1n1nW0W1W1W1PtEta,dM0111); 
    }
   } // end of if(dM0111)
  
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
  
} // end of AliFlowAnalysisWithQCumulants::CalculateWeightedCorrelationsForDifferentialFlow(TString type)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::FinalizeCorrelationsForDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights)
{
 // 1.) Access average for 2D correlations from profiles and store them in 2D final results histograms;
 // 2.) Access spread for 2D correlations from profiles, calculate error and store it in 2D final results histograms;
 // 3.) Make projections along pt and eta axis and store results and errors in 1D final results histograms.
 
 Int_t typeFlag = -1;
 Int_t pWeightsFlag = -1;
 Int_t eWeightsFlag = -1;

 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   } else 
     {
      cout<<"WARNING: type must be either RP or POI in AFAWQC::FCFDF() !!!!"<<endl;
      exit(0);
     }
     
 if(!useParticleWeights)
 {
  pWeightsFlag = 0;
 } else 
   {
    pWeightsFlag = 1;   
   }   
   
 if(eventWeights == "exact")
 {
  eWeightsFlag = 0;
 }          
  
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pW = pWeightsFlag;
 Int_t eW = eWeightsFlag;
 
 // from 2D histogram fNonEmptyBins2D make two 1D histograms fNonEmptyBins1D in pt and eta (to be improved (i.e. moved somewhere else))  
 // pt:
 for(Int_t p=1;p<fnBinsPt;p++)
 {
  Double_t contentPt = 0.;
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   contentPt += (fNonEmptyBins2D[t]->GetBinContent(fNonEmptyBins2D[t]->GetBin(p,e)));          
  }
  fNonEmptyBins1D[t][0]->SetBinContent(p,contentPt);
 }
 // eta:
 for(Int_t e=1;e<fnBinsEta;e++)
 {
  Double_t contentEta = 0.;
  for(Int_t p=1;p<=fnBinsPt;p++)
  {
   contentEta += (fNonEmptyBins2D[t]->GetBinContent(fNonEmptyBins2D[t]->GetBin(p,e)));          
  }
  fNonEmptyBins1D[t][1]->SetBinContent(e,contentEta);
 }
 
 // from 2D profile in (pt,eta) make two 1D profiles in (pt) and (eta):
 TProfile *profile[2][4]; // [0=pt,1=eta][correlation index] // to be improved (do not hardwire the correlation index)
 
 for(Int_t pe=0;pe<2;pe++) // pt or eta
 {
  for(Int_t ci=0;ci<4;ci++) // correlation index
  {
   if(pe==0) profile[pe][ci] = this->MakePtProjection(fCorrelationsPro[t][pW][eW][ci]);
   if(pe==1) profile[pe][ci] = this->MakeEtaProjection(fCorrelationsPro[t][pW][eW][ci]);
  }
 }
  
 // transfer 2D profile into 2D histogram:
 // to be improved (see in documentation if there is a method to transfer values from 2D profile into 2D histogram)    
 for(Int_t ci=0;ci<4;ci++)
 {
  for(Int_t p=1;p<=fnBinsPt;p++)
  {
   for(Int_t e=1;e<=fnBinsEta;e++)
   {
    Double_t correlation = fCorrelationsPro[t][pW][eW][ci]->GetBinContent(fCorrelationsPro[t][pW][eW][ci]->GetBin(p,e)); 
    Double_t spread = fCorrelationsPro[t][pW][eW][ci]->GetBinError(fCorrelationsPro[t][pW][eW][ci]->GetBin(p,e));
    Double_t nEvts = fNonEmptyBins2D[t]->GetBinContent(fNonEmptyBins2D[t]->GetBin(p,e));
    Double_t error = 0.;
    fFinalCorrelations2D[t][pW][eW][ci]->SetBinContent(fFinalCorrelations2D[t][pW][eW][ci]->GetBin(p,e),correlation);          
    if(nEvts>0)
    {
     error = spread/pow(nEvts,0.5);
     fFinalCorrelations2D[t][pW][eW][ci]->SetBinError(fFinalCorrelations2D[t][pW][eW][ci]->GetBin(p,e),error);
    }
   } // end of for(Int_t e=1;e<=fnBinsEta;e++)
  } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 } // end of for(Int_t ci=0;ci<4;ci++)
 
 // transfer 1D profile into 1D histogram (pt):
 // to be improved (see in documentation if there is a method to transfer values from 1D profile into 1D histogram)    
 for(Int_t ci=0;ci<4;ci++)
 {
  for(Int_t p=1;p<=fnBinsPt;p++)
  {
   if(profile[0][ci])
   {
    Double_t correlation = profile[0][ci]->GetBinContent(p); 
    Double_t spread = profile[0][ci]->GetBinError(p);
    Double_t nEvts = fNonEmptyBins1D[t][0]->GetBinContent(p);
    Double_t error = 0.;
    fFinalCorrelations1D[t][pW][eW][0][ci]->SetBinContent(p,correlation); 
    if(nEvts>0)
    {
     error = spread/pow(nEvts,0.5);
     fFinalCorrelations1D[t][pW][eW][0][ci]->SetBinError(p,error);
    }  
   }   
  } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 } // end of for(Int_t ci=0;ci<4;ci++)
 
 // transfer 1D profile into 1D histogram (eta):
 // to be improved (see in documentation if there is a method to transfer values from 1D profile into 1D histogram)    
 for(Int_t ci=0;ci<4;ci++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   if(profile[1][ci])
   {
    Double_t correlation = profile[1][ci]->GetBinContent(e); 
    fFinalCorrelations1D[t][pW][eW][1][ci]->SetBinContent(e,correlation);      
   }    
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t ci=0;ci<4;ci++)
  
} // end of void AliFlowAnalysisWithQCumulants::FinalizeCorrelationsForDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateCumulantsForDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights)
{
 // calcualate cumulants for differential flow from measured correlations
 // Remark: cumulants calculated here are NOT corrected for non-uniform acceptance. This correction is applied in the method ...
  
 Int_t typeFlag = -1;
 Int_t pWeightsFlag = -1;
 Int_t eWeightsFlag = -1;

 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   } else 
     {
      cout<<"WARNING: type must be either RP or POI in AFAWQC::CCFDF() !!!!"<<endl;
      exit(0);
     }
     
 if(!useParticleWeights)
 {
  pWeightsFlag = 0;
 } else 
   {
    pWeightsFlag = 1;   
   }       
  
 if(eventWeights == "exact")
 {
  eWeightsFlag = 0;
 } 
  
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pW = pWeightsFlag; 
 Int_t eW = eWeightsFlag; 
 
 // correlation <<2>>: 
 Double_t two = fQCorrelations[pW][eW]->GetBinContent(1);
  
 // 2D (pt,eta):
 // to be improved (see documentation if I can do all this without looping)
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++) 
  {  
   // reduced correlations:   
   Double_t twoPrime = fFinalCorrelations2D[t][pW][eW][0]->GetBinContent(fFinalCorrelations2D[t][pW][eW][0]->GetBin(p,e)); // <<2'>>(pt,eta)
   Double_t fourPrime = fFinalCorrelations2D[t][pW][eW][1]->GetBinContent(fFinalCorrelations2D[t][pW][eW][1]->GetBin(p,e)); // <<4'>>(pt,eta)
   for(Int_t nua=0;nua<2;nua++)
   {
    // QC{2'}:
    Double_t qc2Prime = twoPrime; // QC{2'} = <<2'>>
    fFinalCumulants2D[t][pW][eW][nua][0]->SetBinContent(fFinalCumulants2D[t][pW][eW][nua][0]->GetBin(p,e),qc2Prime);    
    // QC{4'}:
    Double_t qc4Prime = fourPrime - 2.*twoPrime*two; // QC{4'} = <<4'>> - 2*<<2'>><<2>>
    fFinalCumulants2D[t][pW][eW][nua][1]->SetBinContent(fFinalCumulants2D[t][pW][eW][nua][1]->GetBin(p,e),qc4Prime);   
   } // end of for(Int_t nua=0;nua<2;nua++)   
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
 // 1D (pt):
 // to be improved (see documentation if I can do all this without looping)
 // to be improved (treat pt and eta in one go)
 // to be improved (combine loops over nua for 2D and 1D)
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  // reduced correlations:   
  Double_t twoPrime = fFinalCorrelations1D[t][pW][eW][0][0]->GetBinContent(p); // <<2'>>(pt)
  Double_t fourPrime = fFinalCorrelations1D[t][pW][eW][0][1]->GetBinContent(p); // <<4'>>(pt)
  // spread of reduced correlations:
  Double_t twoPrimeError = fFinalCorrelations1D[t][pW][eW][0][0]->GetBinError(p); // sigma_2'/sqrt{N_2'}(pt)
  for(Int_t nua=0;nua<2;nua++)
  {
   // QC{2'}:
   Double_t qc2Prime = twoPrime; // QC{2'}
   Double_t qc2PrimeError = twoPrimeError; // sigma_{d_n{2}}/sqrt{N_2'} // to be improved
   fFinalCumulantsPt[t][pW][eW][nua][0]->SetBinContent(p,qc2Prime); 
   fFinalCumulantsPt[t][pW][eW][nua][0]->SetBinError(p,qc2PrimeError);   
   // QC{4'}:
   Double_t qc4Prime = fourPrime - 2.*twoPrime*two; // QC{4'} = <<4'>> - 2*<<2'>><<2>>
   fFinalCumulantsPt[t][pW][eW][nua][1]->SetBinContent(p,qc4Prime); 
  } // end of for(Int_t nua=0;nua<2;nua++) 
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
 
 // 1D (eta):
 // to be improved (see documentation if I can do all this without looping)
 // to be improved (treat pt and eta in one go)
 // to be improved (combine loops over nua for 2D and 1D)
 for(Int_t e=1;e<=fnBinsEta;e++)
 {
  // reduced correlations:   
  Double_t twoPrime = fFinalCorrelations1D[t][pW][eW][1][0]->GetBinContent(e); // <<2'>>(eta)
  Double_t fourPrime = fFinalCorrelations1D[t][pW][eW][1][1]->GetBinContent(e); // <<4'>>(eta)
  for(Int_t nua=0;nua<2;nua++)
  {
   // QC{2'}:
   Double_t qc2Prime = twoPrime; // QC{2'}
   fFinalCumulantsEta[t][pW][eW][nua][0]->SetBinContent(e,qc2Prime);    
   // QC{4'}:
   Double_t qc4Prime = fourPrime - 2.*twoPrime*two; // QC{4'} = <<4'>> - 2*<<2'>><<2>>
   fFinalCumulantsEta[t][pW][eW][nua][1]->SetBinContent(e,qc4Prime);
  } // end of for(Int_t nua=0;nua<2;nua++)   
 } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 
  
} // end of void AliFlowAnalysisWithQCumulants::CalculateCumulantsForDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights); 


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights)
{
 // calculate differential flow from differential cumulants and previously obtained integrated flow:
 
 Int_t typeFlag = -1;
 Int_t pWeightsFlag = -1;
 Int_t eWeightsFlag = -1;

 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   } else 
     {
      cout<<"WARNING: type must be either RP or POI in AFAWQC::CDF() !!!!"<<endl;
      exit(0);
     }
     
 if(!useParticleWeights)
 {
  pWeightsFlag = 0;
 } else 
   {
    pWeightsFlag = 1;   
   }       
  
 if(eventWeights == "exact")
 {
  eWeightsFlag = 0;
 }  
  
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pW = pWeightsFlag; 
 Int_t eW = eWeightsFlag; 
 
 // integrated flow:
 Double_t v2 = 0.;
 Double_t v4 = 0.;  
 
 if(pW == 0)
 {
  v2 = fIntFlow[pW][eW][1]->GetBinContent(1);   
  v4 = fIntFlow[pW][eW][1]->GetBinContent(2);    
 }
     
 if(pW == 1)
 {
  v2 = fIntFlow[pW][eW][0]->GetBinContent(1);   
  v4 = fIntFlow[pW][eW][0]->GetBinContent(2);    
 }
 
 // 2D:
 for(Int_t nua=0;nua<2;nua++)
 {
  for(Int_t p=1;p<=fnBinsPt;p++)
  {
   for(Int_t e=1;e<=fnBinsEta;e++) 
   { 
    // differential cumulants:
    Double_t qc2Prime = fFinalCumulants2D[t][pW][eW][nua][0]->GetBinContent(fFinalCumulants2D[t][pW][eW][nua][0]->GetBin(p,e)); // QC{2'}                    
    Double_t qc4Prime = fFinalCumulants2D[t][pW][eW][nua][1]->GetBinContent(fFinalCumulants2D[t][pW][eW][nua][1]->GetBin(p,e)); // QC{4'}
    // differential flow:
    Double_t v2Prime = 0.;                    
    Double_t v4Prime = 0.; 
    if(v2) 
    {
     v2Prime = qc2Prime/v2;
     fFinalFlow2D[t][pW][eW][nua][0]->SetBinContent(fFinalFlow2D[t][pW][eW][nua][0]->GetBin(p,e),v2Prime);  
    }                   
    if(v4)
    {
     v4Prime = -qc4Prime/pow(v4,3.); 
     fFinalFlow2D[t][pW][eW][nua][1]->SetBinContent(fFinalFlow2D[t][pW][eW][nua][1]->GetBin(p,e),v4Prime);  
    }                    
   } // end of for(Int_t e=1;e<=fnBinsEta;e++)
  } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 } // end of for(Int_t nua=0;nua<2;nua++)
 
 // 1D (pt): // to be improved (combined with eta, combined nua loop with 2D)
 for(Int_t nua=0;nua<2;nua++)
 {
  for(Int_t p=1;p<=fnBinsPt;p++)
  {
   // differential cumulants:
   Double_t qc2Prime = fFinalCumulantsPt[t][pW][eW][nua][0]->GetBinContent(p); // QC{2'}                    
   Double_t qc4Prime = fFinalCumulantsPt[t][pW][eW][nua][1]->GetBinContent(p); // QC{4'}
   // differential flow:
   Double_t v2Prime = 0.;                    
   Double_t v4Prime = 0.; 
   if(v2) 
   {
    v2Prime = qc2Prime/v2;
    fFinalFlowPt[t][pW][eW][nua][0]->SetBinContent(p,v2Prime);  
   }                   
   if(v4)
   {
    v4Prime = -qc4Prime/pow(v4,3.); 
    fFinalFlowPt[t][pW][eW][nua][1]->SetBinContent(p,v4Prime);  
   }                    
  } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 } // end of for(Int_t nua=0;nua<2;nua++)
 
 // 1D (eta): // to be improved (combined with pt, combined nua loop with 2D)
 for(Int_t nua=0;nua<2;nua++)
 {
  for(Int_t e=1;e<=fnBinsEta;e++)
  {
   // differential cumulants:
   Double_t qc2Prime = fFinalCumulantsEta[t][pW][eW][nua][0]->GetBinContent(e); // QC{2'}                    
   Double_t qc4Prime = fFinalCumulantsEta[t][pW][eW][nua][1]->GetBinContent(e); // QC{4'}
   // differential flow:
   Double_t v2Prime = 0.;                    
   Double_t v4Prime = 0.; 
   if(v2) 
   {
    v2Prime = qc2Prime/v2;
    fFinalFlowEta[t][pW][eW][nua][0]->SetBinContent(e,v2Prime);  
   }                   
   if(v4)
   {
    v4Prime = -qc4Prime/pow(v4,3.); 
    fFinalFlowEta[t][pW][eW][nua][1]->SetBinContent(e,v4Prime);  
   }                    
  } // end of for(Int_t e=1;e<=fnBinsEta;e++)
 } // end of for(Int_t nua=0;nua<2;nua++)

} // end of AliFlowAnalysisWithQCumulants::CalculateDiffFlow(TString type, Bool_t useParticleWeights)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateFinalResultsForRPandPOIIntegratedFlow(TString type, Bool_t useParticleWeights, TString eventWeights)
{
 // calculate final results for integrated flow of RPs and POIs 
 
 Int_t typeFlag = -1;
 Int_t pWeightsFlag = -1;
 Int_t eWeightsFlag = -1;

 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   } else 
     {
      cout<<"WARNING: type must be either RP or POI in AFAWQC::CDF() !!!!"<<endl;
      exit(0);
     }
     
 if(!useParticleWeights)
 {
  pWeightsFlag = 0;
 } else 
   {
    pWeightsFlag = 1;   
   }       
  
 if(eventWeights == "exact")
 {
  eWeightsFlag = 0;
 } 
  
 // shortcuts:
 Int_t t = typeFlag;
 Int_t pW = pWeightsFlag; 
 Int_t eW = eWeightsFlag; 
  
 // pt yield:    
 TH1F *yield2ndPt = NULL;
 TH1F *yield4thPt = NULL;
 TH1F *yield6thPt = NULL;
 TH1F *yield8thPt = NULL;
 
 if(type == "POI")
 {
  yield2ndPt = (TH1F*)(fCommonHists2nd->GetHistPtPOI())->Clone();
  yield4thPt = (TH1F*)(fCommonHists4th->GetHistPtPOI())->Clone();
  yield6thPt = (TH1F*)(fCommonHists6th->GetHistPtPOI())->Clone();
  yield8thPt = (TH1F*)(fCommonHists8th->GetHistPtPOI())->Clone();  
 } 
 else if(type == "RP")
 {
  yield2ndPt = (TH1F*)(fCommonHists2nd->GetHistPtRP())->Clone();
  yield4thPt = (TH1F*)(fCommonHists4th->GetHistPtRP())->Clone();
  yield6thPt = (TH1F*)(fCommonHists6th->GetHistPtRP())->Clone();
  yield8thPt = (TH1F*)(fCommonHists8th->GetHistPtRP())->Clone();  
 } 
 
 Int_t nBinsPt = yield2ndPt->GetNbinsX();
 
 TH1D *flow2ndPt = NULL;
 TH1D *flow4thPt = NULL;
 TH1D *flow6thPt = NULL;
 TH1D *flow8thPt = NULL;
 
 // to be improve (nua = 0 is hardwired)
 flow2ndPt = (TH1D*)fFinalFlowPt[t][pW][eW][0][0]->Clone();
 flow4thPt = (TH1D*)fFinalFlowPt[t][pW][eW][0][1]->Clone();
 flow6thPt = (TH1D*)fFinalFlowPt[t][pW][eW][0][2]->Clone();
 flow8thPt = (TH1D*)fFinalFlowPt[t][pW][eW][0][3]->Clone(); 
   
 Double_t dvn2nd = 0., dvn4th = 0., dvn6th = 0., dvn8th = 0.; // differential flow
 Double_t dVn2nd = 0., dVn4th = 0., dVn6th = 0., dVn8th = 0.; // integrated flow 
 //Double_t dSd2nd = 0., dSd4th = 0., dSd6th = 0., dSd8th = 0.; // error on integrated flow (to be improved - calculation needed)

 Double_t dYield2nd = 0., dYield4th = 0., dYield6th = 0., dYield8th = 0.; // pt yield 
 Double_t dSum2nd = 0., dSum4th = 0., dSum6th = 0., dSum8th = 0.; // needed for normalizing integrated flow
 
 // looping over pt bins:
 for(Int_t p=1;p<nBinsPt+1;p++)
 {
  dvn2nd = flow2ndPt->GetBinContent(p);
  dvn4th = flow4thPt->GetBinContent(p);
  dvn6th = flow6thPt->GetBinContent(p);
  dvn8th = flow8thPt->GetBinContent(p);

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
  
  // ... to be improved - errors needed to be calculated   
  
 } // end of for(Int_t p=1;p<nBinsPt+1;p++)

 // normalizing the results for integrated flow:
 
 
 
 
 
 if(dSum2nd) dVn2nd/=dSum2nd;
 if(dSum4th) dVn4th/=dSum4th;
 //if(dSum6th) dVn6th/=dSum6th;
 //if(dSum8th) dVn8th/=dSum8th;
 
 /*
 // storing the results for integrated flow:
 if(!(useParticleWeights))
 {
  if(type == "POI")
  {
   // 2nd:
   fIntFlowResultsPOIQC->SetBinContent(1,dVn2nd);
   fIntFlowResultsPOIQC->SetBinError(1,dSd2nd);
   // 4th:
   fIntFlowResultsPOIQC->SetBinContent(2,dVn4th);
   fIntFlowResultsPOIQC->SetBinError(2,dSd4th);
   // 6th:
   fIntFlowResultsPOIQC->SetBinContent(3,dVn6th);
   fIntFlowResultsPOIQC->SetBinError(3,dSd6th);
   // 8th:
   fIntFlowResultsPOIQC->SetBinContent(4,dVn8th);
   fIntFlowResultsPOIQC->SetBinError(4,dSd8th);
  }
  else if (type == "RP")
  {
   // 2nd:
   fIntFlowResultsRPQC->SetBinContent(1,dVn2nd);
   fIntFlowResultsRPQC->SetBinError(1,dSd2nd);
   // 4th:
   fIntFlowResultsRPQC->SetBinContent(2,dVn4th);
   fIntFlowResultsRPQC->SetBinError(2,dSd4th);
   // 6th:
   fIntFlowResultsRPQC->SetBinContent(3,dVn6th);
   fIntFlowResultsRPQC->SetBinError(3,dSd6th);
   // 8th:
   fIntFlowResultsRPQC->SetBinContent(4,dVn8th);
   fIntFlowResultsRPQC->SetBinError(4,dSd8th);
  }
 } 
 else if (useParticleWeights)
 {
  if(type == "POI")
  {
   // 2nd:
   fIntFlowResultsPOIQCW->SetBinContent(1,dVn2nd);
   fIntFlowResultsPOIQCW->SetBinError(1,dSd2nd);
   // 4th:
   fIntFlowResultsPOIQCW->SetBinContent(2,dVn4th);
   fIntFlowResultsPOIQCW->SetBinError(2,dSd4th);
   // 6th:
   fIntFlowResultsPOIQCW->SetBinContent(3,dVn6th);
   fIntFlowResultsPOIQCW->SetBinError(3,dSd6th);
   // 8th:
   fIntFlowResultsPOIQCW->SetBinContent(4,dVn8th);
   fIntFlowResultsPOIQCW->SetBinError(4,dSd8th);
  }
  else if (type == "RP")
  {
   // 2nd:
   fIntFlowResultsRPQCW->SetBinContent(1,dVn2nd);
   fIntFlowResultsRPQCW->SetBinError(1,dSd2nd);
   // 4th:
   fIntFlowResultsRPQCW->SetBinContent(2,dVn4th);
   fIntFlowResultsRPQCW->SetBinError(2,dSd4th);
   // 6th:
   fIntFlowResultsRPQCW->SetBinContent(3,dVn6th);
   fIntFlowResultsRPQCW->SetBinError(3,dSd6th);
   // 8th:
   fIntFlowResultsRPQCW->SetBinContent(4,dVn8th);
   fIntFlowResultsRPQCW->SetBinError(4,dSd8th);
  }
 }
 */
 
 // storing the results for integrated flow in common histos:
 // to be improved - now they are being filled twice ...
 if(type == "POI")
 {
  fCommonHistsResults2nd->FillIntegratedFlowPOI(dVn2nd,0.); // to be improved (errors)
  fCommonHistsResults4th->FillIntegratedFlowPOI(dVn4th,0.); // to be improved (errors)
  fCommonHistsResults6th->FillIntegratedFlowPOI(dVn6th,0.); // to be improved (errors)
  fCommonHistsResults8th->FillIntegratedFlowPOI(dVn8th,0.); // to be improved (errors)
 }
 else if (type == "RP")
 {
  fCommonHistsResults2nd->FillIntegratedFlowRP(dVn2nd,0.); // to be improved (errors)
  fCommonHistsResults4th->FillIntegratedFlowRP(dVn4th,0.); // to be improved (errors)
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
  
} // end of AliFlowAnalysisWithQCumulants::CalculateFinalResultsForRPandPOIIntegratedFlow(TString type, Bool_t useParticleWeights, TString eventWeights)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::InitializeArraysForDistributions()
{
 // initialize arrays used for distributions:
 
 for(Int_t pW=0;pW<2;pW++) // particle weights not used (0) or used (1)
 {
  for(Int_t eW=0;eW<2;eW++)
  {
   for(Int_t di=0;di<4;di++) // distribution index
   {
    fDistributions[pW][eW][di] = NULL;
   }
  } 
 }
 
} // end of void AliFlowAnalysisWithQCumulants::InitializeArraysForDistributions()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::BookEverythingForDistributions()
{
 // book all histograms for distributions
 
 /*
 //weighted <2>_{n|n} distribution
 f2pDistribution = new TH1D("f2pDistribution","<2>_{n|n} distribution",100000,-0.02,0.1);
 f2pDistribution->SetXTitle("<2>_{n|n}");
 f2pDistribution->SetYTitle("Counts");
 fHistList->Add(f2pDistribution);

 //weighted <4>_{n,n|n,n} distribution
 f4pDistribution = new TH1D("f4pDistribution","<4>_{n,n|n,n} distribution",100000,-0.00025,0.002);
 f4pDistribution->SetXTitle("<4>_{n,n|n,n}");
 f4pDistribution->SetYTitle("Counts");
 fHistList->Add(f4pDistribution); 
 
 //weighted <6>_{n,n,n|n,n,n} distribution
 f6pDistribution = new TH1D("f6pDistribution","<6>_{n,n,n|n,n,n} distribution",100000,-0.000005,0.000025);
 f6pDistribution->SetXTitle("<6>_{n,n,n|n,n,n}");
 f6pDistribution->SetYTitle("Counts");
 fHistList->Add(f6pDistribution);
 
 //weighted <8>_{n,n,n,n|n,n,n,n} distribution
 f8pDistribution = new TH1D("f8pDistribution","<8>_{n,n,n,n|n,n,n,n} distribution",100000,-0.000000001,0.00000001);
 f8pDistribution->SetXTitle("<8>_{n,n,n,n|n,n,n,n}");
 f8pDistribution->SetYTitle("Counts");
 fHistList->Add(f8pDistribution);
 */
 
} // end of void AliFlowAnalysisWithQCumulants::BookEverythingForDistributions()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::BookAndNestAllLists()
{
 // book and nest all lists in fHistList 
 
 // booking and nesting:  
 fWeightsList->SetName("Weights");
 fWeightsList->SetOwner(kTRUE);   
 fHistList->Add(fWeightsList); 
 
 fIntFlowList = new TList();
 fIntFlowList->SetName("Integrated Flow");
 fIntFlowList->SetOwner(kTRUE);
 fHistList->Add(fIntFlowList);
 
 fIntFlowResults = new TList();
 fIntFlowResults->SetName("Results");
 fIntFlowResults->SetOwner(kTRUE);
 fIntFlowList->Add(fIntFlowResults);
 
 fIntFlowProfiles = new TList();
 fIntFlowProfiles->SetName("Profiles");
 fIntFlowProfiles->SetOwner(kTRUE);
 fIntFlowList->Add(fIntFlowProfiles);
 
 fDiffFlowList = new TList();
 fDiffFlowList->SetName("Differential Flow");
 fDiffFlowList->SetOwner(kTRUE); 
 fHistList->Add(fDiffFlowList);
   
 fDiffFlowProfiles = new TList(); 
 fDiffFlowProfiles->SetName("Profiles");
 fDiffFlowProfiles->SetOwner(kTRUE);
 fDiffFlowList->Add(fDiffFlowProfiles);
 
 fDiffFlowResults = new TList();
 fDiffFlowResults->SetName("Results");
 fDiffFlowResults->SetOwner(kTRUE);
 fDiffFlowList->Add(fDiffFlowResults);
 
 fDistributionsList = new TList();
 fDistributionsList->SetName("Distributions");
 fDistributionsList->SetOwner(kTRUE);
 fHistList->Add(fDistributionsList);
 
 fNestedLoopsList = new TList();
 fNestedLoopsList->SetName("Nested Loops");
 fNestedLoopsList->SetOwner(kTRUE);
 fHistList->Add(fNestedLoopsList);
 
 // flags for naming nested lists:  
 TList list;
 list.SetOwner(kTRUE);
 TString typeFlag[2] = {"RP","POI"}; 
 TString pWeightsFlag[2] = {"not used","used"}; 
 TString eWeightsFlag[2] = {"exact","non-exact"};
 TString nuaFlag[2] = {"not corrected","corrected"}; // nua = non-uniform acceptance
 TString sinCosFlag[2] = {"sin","cos"};

 // nested lists in fDiffFlowProfiles (Differential Flow/Profiles): 
 for(Int_t t=0;t<2;t++) // type: RP or POI
 {
  fDFPType[t] = (TList*)list.Clone();
  fDFPType[t]->SetName(Form("%s",typeFlag[t].Data()));
  fDiffFlowProfiles->Add(fDFPType[t]);  
  for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++) // weights: not used or used  
  {
   fDFPParticleWeights[t][pW] = (TList*)list.Clone();
   fDFPParticleWeights[t][pW]->SetName(Form("%s, pWeights %s",typeFlag[t].Data(),pWeightsFlag[pW].Data()));
   fDFPType[t]->Add(fDFPParticleWeights[t][pW]);
   for(Int_t eW=0;eW<2;eW++)
   { 
    fDFPEventWeights[t][pW][eW] = (TList*)list.Clone();
    fDFPEventWeights[t][pW][eW]->SetName(Form("%s, pWeights %s, eWeights %s",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()));
    fDFPParticleWeights[t][pW]->Add(fDFPEventWeights[t][pW][eW]);
    // correlations: 
    fDiffFlowCorrelations[t][pW][eW] = (TList*)list.Clone();
    fDiffFlowCorrelations[t][pW][eW]->SetName(Form("Correlations (%s, pWeights %s, eWeights %s)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()));
    fDFPEventWeights[t][pW][eW]->Add(fDiffFlowCorrelations[t][pW][eW]);
    // products of correlations:
    fDiffFlowProductsOfCorrelations[t][pW][eW] = (TList*)list.Clone();
    fDiffFlowProductsOfCorrelations[t][pW][eW]->SetName(Form("Products of correlations (%s, pWeights %s, eWeights %s)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()));
    fDFPEventWeights[t][pW][eW]->Add(fDiffFlowProductsOfCorrelations[t][pW][eW]);
    // correction terms:   
    for(Int_t sc=0;sc<2;sc++) // corrections for NUA: sin or cos terms 
    {
     fDiffFlowCorrectionTerms[t][pW][eW][sc] = (TList*)list.Clone();
     fDiffFlowCorrectionTerms[t][pW][eW][sc]->SetName(Form("Corrections for NUA (%s, pWeights %s, eWeights %s, %s terms)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),sinCosFlag[sc].Data()));
     fDFPEventWeights[t][pW][eW]->Add(fDiffFlowCorrectionTerms[t][pW][eW][sc]);
    }
   }
  } 
 } 
 
 // nested lists in fDiffFlowResults (Differential Flow/Results):
 for(Int_t t=0;t<2;t++) // RP or POI
 {
  fDFRType[t] = (TList*)list.Clone();
  fDFRType[t]->SetName(Form("%s",typeFlag[t].Data()));
  fDiffFlowResults->Add(fDFRType[t]);  
  for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++) // particle weights not used or used
  {
   fDFRParticleWeights[t][pW] = (TList*)list.Clone();
   fDFRParticleWeights[t][pW]->SetName(Form("%s, pWeights %s",typeFlag[t].Data(),pWeightsFlag[pW].Data()));
   fDFRType[t]->Add(fDFRParticleWeights[t][pW]);
   for(Int_t eW=0;eW<2;eW++) // event weights: exact ot non-exact // to be improved (terminology)
   {
    fDFREventWeights[t][pW][eW] = (TList*)list.Clone();
    fDFREventWeights[t][pW][eW]->SetName(Form("%s, pWeights %s, eWeights %s",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()));
    fDFRParticleWeights[t][pW]->Add(fDFREventWeights[t][pW][eW]);
    // final correlations:    
    fDiffFlowFinalCorrelations[t][pW][eW] = (TList*)list.Clone();
    fDiffFlowFinalCorrelations[t][pW][eW]->SetName(Form("Correlations (%s, pWeights %s, eWeights %s)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()));
    fDFREventWeights[t][pW][eW]->Add(fDiffFlowFinalCorrelations[t][pW][eW]); 
    // final covariances:    
    fDiffFlowFinalCovariances[t][pW][eW] = (TList*)list.Clone();
    fDiffFlowFinalCovariances[t][pW][eW]->SetName(Form("Covariances (%s, pWeights %s, eWeights %s)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()));
    fDFREventWeights[t][pW][eW]->Add(fDiffFlowFinalCovariances[t][pW][eW]);
    // final corrections:    
    fDiffFlowFinalCorrections[t][pW][eW] = (TList*)list.Clone();
    fDiffFlowFinalCorrections[t][pW][eW]->SetName(Form("Corrections (%s, pWeights %s, eWeights %s)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data()));
    fDFREventWeights[t][pW][eW]->Add(fDiffFlowFinalCorrections[t][pW][eW]);
    for(Int_t nua=0;nua<2;nua++)
    {
     fDFRCorrections[t][pW][eW][nua] = (TList*)list.Clone();
     fDFRCorrections[t][pW][eW][nua]->SetName(Form("%s, pWeights %s, eWeights %s, %s for NUA",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data()));
     fDFREventWeights[t][pW][eW]->Add(fDFRCorrections[t][pW][eW][nua]); 
     // final cumulants:    
     fDiffFlowFinalCumulants[t][pW][eW][nua] = (TList*)list.Clone();
     fDiffFlowFinalCumulants[t][pW][eW][nua]->SetName(Form("Cumulants (%s, pWeights %s, eWeights %s, %s for NUA)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data()));
     fDFRCorrections[t][pW][eW][nua]->Add(fDiffFlowFinalCumulants[t][pW][eW][nua]); 
     // final differential flow:
     fDiffFlowFinalFlow[t][pW][eW][nua] = (TList*)list.Clone();
     fDiffFlowFinalFlow[t][pW][eW][nua]->SetName(Form("Differential Flow (%s, pWeights %s, eWeights %s, %s for NUA)",typeFlag[t].Data(),pWeightsFlag[pW].Data(),eWeightsFlag[eW].Data(),nuaFlag[nua].Data()));
     fDFRCorrections[t][pW][eW][nua]->Add(fDiffFlowFinalFlow[t][pW][eW][nua]);          
    }
   } 
  }
 } 
 
} // end of void AliFlowAnalysisWithQCumulants::BookAndNestAllLists()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::FillCommonHistResultsDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights, Bool_t correctedForNUA)
{
 // fill common result histograms for differential flow
 
 // shortcuts for the flags:
 Int_t t = 0; // 0=RP, 1=POI
 Int_t pW = (Int_t)(useParticleWeights); // 0=weights not useed, 1=weights used
 Int_t eW = -1; 
 Int_t nua = (Int_t)(correctedForNUA); // 0=not corrected for NUA, 1=corrected for NUA
 
 if(type == "POI")
 {
  t = 1;
 }
 
 if(eventWeights == "exact")
 {
  eW = 0;
 }
  
 for(Int_t o=0;o<4;o++) // order
 {
  if(!fFinalFlowPt[t][pW][eW][nua][o])
  {
    cout<<"WARNING: fFinalFlowPt[t][pW][eW][nua][o] is NULL in AFAWQC::FCHRIF() !!!!"<<endl;
    cout<<"t = "<<t<<endl;
    cout<<"pW = "<<pW<<endl;
    cout<<"eW = "<<eW<<endl;
    cout<<"nua = "<<nua<<endl;
    cout<<"o = "<<o<<endl;
    exit(0); 
  }
  if(!fFinalFlowEta[t][pW][eW][nua][o])
  {
    cout<<"WARNING: fFinalFlowEta[t][pW][eW][nua][o] is NULL in AFAWQC::FCHRIF() !!!!"<<endl;
    cout<<"t = "<<t<<endl;
    cout<<"pW = "<<pW<<endl;
    cout<<"eW = "<<eW<<endl;
    cout<<"nua = "<<nua<<endl;
    cout<<"o = "<<o<<endl;
    exit(0); 
  }  
 } 
   
 if(!(fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th))
 {
  cout<<"WARNING: fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th"<<endl; 
  cout<<"         is NULL in AFAWQC::FCHRIF() !!!!"<<endl;
  exit(0);
 }
 
 // pt:
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  Double_t v2 = fFinalFlowPt[t][pW][eW][nua][0]->GetBinContent(p);
  Double_t v4 = fFinalFlowPt[t][pW][eW][nua][1]->GetBinContent(p);
  Double_t v6 = fFinalFlowPt[t][pW][eW][nua][2]->GetBinContent(p);
  Double_t v8 = fFinalFlowPt[t][pW][eW][nua][3]->GetBinContent(p);
  
  //Double_t v2Error = fFinalFlow1D[t][pW][nua][0][0]->GetBinError(p);
  //Double_t v4Error = fFinalFlow1D[t][pW][nua][0][1]->GetBinError(p);
  //Double_t v6Error = fFinalFlow1D[t][pW][nua][0][2]->GetBinError(p);
  //Double_t v8Error = fFinalFlow1D[t][pW][nua][0][3]->GetBinError(p);
 
  if(type == "RP")
  {
   fCommonHistsResults2nd->FillDifferentialFlowPtRP(p,v2,0.);
   fCommonHistsResults4th->FillDifferentialFlowPtRP(p,v4,0.);
   fCommonHistsResults6th->FillDifferentialFlowPtRP(p,v6,0.);
   fCommonHistsResults8th->FillDifferentialFlowPtRP(p,v8,0.);
  } else if(type == "POI")
    {
     fCommonHistsResults2nd->FillDifferentialFlowPtPOI(p,v2,0.);
     fCommonHistsResults4th->FillDifferentialFlowPtPOI(p,v4,0.);
     fCommonHistsResults6th->FillDifferentialFlowPtPOI(p,v6,0.);
     fCommonHistsResults8th->FillDifferentialFlowPtPOI(p,v8,0.);
    }
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)   
 
 // eta:
 for(Int_t e=1;e<=fnBinsEta;e++)
 {
  Double_t v2 = fFinalFlowEta[t][pW][eW][nua][0]->GetBinContent(e);
  Double_t v4 = fFinalFlowEta[t][pW][eW][nua][1]->GetBinContent(e);
  Double_t v6 = fFinalFlowEta[t][pW][eW][nua][2]->GetBinContent(e);
  Double_t v8 = fFinalFlowEta[t][pW][eW][nua][3]->GetBinContent(e);
  
  //Double_t v2Error = fFinalFlow1D[t][pW][nua][1][0]->GetBinError(e);
  //Double_t v4Error = fFinalFlow1D[t][pW][nua][1][1]->GetBinError(e);
  //Double_t v6Error = fFinalFlow1D[t][pW][nua][1][2]->GetBinError(e);
  //Double_t v8Error = fFinalFlow1D[t][pW][nua][1][3]->GetBinError(e);
 
  if(type == "RP")
  {
   fCommonHistsResults2nd->FillDifferentialFlowEtaRP(e,v2,0.);
   fCommonHistsResults4th->FillDifferentialFlowEtaRP(e,v4,0.);
   fCommonHistsResults6th->FillDifferentialFlowEtaRP(e,v6,0.);
   fCommonHistsResults8th->FillDifferentialFlowEtaRP(e,v8,0.);
  } else if(type == "POI")
    {
     fCommonHistsResults2nd->FillDifferentialFlowEtaPOI(e,v2,0.);
     fCommonHistsResults4th->FillDifferentialFlowEtaPOI(e,v4,0.);
     fCommonHistsResults6th->FillDifferentialFlowEtaPOI(e,v6,0.);
     fCommonHistsResults8th->FillDifferentialFlowEtaPOI(e,v8,0.);
    }
 } // end of for(Int_t e=1;e<=fnBinsEta;e++)    
 
} // end of void AliFlowAnalysisWithQCumulants::FillCommonHistResultsDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights, Bool_t correctedForNUA)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::AccessConstants()
{
 // access needed common constants from AliFlowCommonConstants
 
 fnBinsPhi = AliFlowCommonConstants::GetNbinsPhi();
 fPhiMin = AliFlowCommonConstants::GetPhiMin();	     
 fPhiMax = AliFlowCommonConstants::GetPhiMax();
 if(fnBinsPhi) fPhiBinWidth = (fPhiMax-fPhiMin)/fnBinsPhi;  
 fnBinsPt = AliFlowCommonConstants::GetNbinsPt();
 fPtMin = AliFlowCommonConstants::GetPtMin();	     
 fPtMax = AliFlowCommonConstants::GetPtMax();
 if(fnBinsPt) fPtBinWidth = (fPtMax-fPtMin)/fnBinsPt;  
 fnBinsEta = AliFlowCommonConstants::GetNbinsEta();
 fEtaMin = AliFlowCommonConstants::GetEtaMin();	     
 fEtaMax = AliFlowCommonConstants::GetEtaMax();
 if(fnBinsEta) fEtaBinWidth = (fEtaMax-fEtaMin)/fnBinsEta;  
 
} // end of void AliFlowAnalysisWithQCumulants::AccessConstants()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateSumAndProductOfEventWeights()
{
 // 1.) calculate sum of linear and quadratic event weights;
 // 2.) calculate products of event weights
 
 Double_t dMult = (*fSMpk)(0,0); // multiplicity (number of particles used to determine the reaction plane)

 Double_t eventWeight[4] = {0}; 
 eventWeight[0] = dMult*(dMult-1); // event weight for <2> 
 eventWeight[1] = dMult*(dMult-1)*(dMult-2)*(dMult-3); // event weight for <4> 
 eventWeight[2] = dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5); // event weight for <6> 
 eventWeight[3] = dMult*(dMult-1)*(dMult-2)*(dMult-3)*(dMult-4)*(dMult-5)*(dMult-6)*(dMult-7); // event weight for <8> 
 
 for(Int_t p=0;p<2;p++)
 {
  for(Int_t c=0;c<4;c++)
  { 
   fSumOfEventWeights[0][0][p]->Fill(c+0.5,pow(eventWeight[c],p+1)); 
  }
 }
  
 // to be improved (hardwired pW and eW):
 fProductOfEventWeights[0][0]->Fill(0.5,eventWeight[0]*eventWeight[1]); 
 fProductOfEventWeights[0][0]->Fill(1.5,eventWeight[0]*eventWeight[2]); 
 fProductOfEventWeights[0][0]->Fill(2.5,eventWeight[0]*eventWeight[3]); 
 fProductOfEventWeights[0][0]->Fill(3.5,eventWeight[1]*eventWeight[2]); 
 fProductOfEventWeights[0][0]->Fill(4.5,eventWeight[1]*eventWeight[3]); 
 fProductOfEventWeights[0][0]->Fill(5.5,eventWeight[2]*eventWeight[3]); 
       
} // end of void AliFlowAnalysisWithQCumulants::CalculateSumAndProductOfEventWeights()


//================================================================================================================================

