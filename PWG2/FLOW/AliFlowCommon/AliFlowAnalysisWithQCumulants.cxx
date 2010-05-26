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
 // 2a.) particle weights:
 fWeightsList(NULL),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
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
 fIntFlowFlags(NULL),
 fApplyCorrectionForNUA(kTRUE), 
 fReQ(NULL),
 fImQ(NULL),
 fSMpk(NULL),
 fIntFlowCorrelationsEBE(NULL),
 fIntFlowEventWeightsForCorrelationsEBE(NULL),
 fIntFlowCorrelationsAllEBE(NULL),
 fAvMultiplicity(NULL),
 fIntFlowCorrelationsPro(NULL),
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
 fIntFlow(NULL),
 // 4.) differential flow:
 fDiffFlowList(NULL),
 fDiffFlowProfiles(NULL),
 fDiffFlowResults(NULL),
 fDiffFlowFlags(NULL),
 fCalculate2DFlow(kFALSE),
 // 5.) distributions:
 fDistributionsList(NULL),
 fDistributionsFlags(NULL),
 fStoreDistributions(kFALSE),
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
 // b) Access all common constants:
 this->AccessConstants();
 // c) Book all objects:
 this->BookAndFillWeightsHistograms();
 this->BookAndNestAllLists();
 this->BookCommonHistograms();
 this->BookEverythingForIntegratedFlow(); 
 this->BookEverythingForDifferentialFlow(); 
 this->BookEverythingForDistributions();
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
 
 //  a) Fill the common control histograms and call the method to fill fAvMultiplicity;
 //  b) Loop over data and calculate e-b-e quantities;
 //  c) Call all the methods;
 //  d) Debugging and cross-checking (evaluate nested loops);
 //  e) Reset all event by event quantities. 
 
 Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
 Double_t dPt  = 0.; // transverse momentum
 Double_t dEta = 0.; // pseudorapidity

 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight
 
 Int_t nRP = anEvent->GetEventNSelTracksRP(); // number of RPs (i.e. number of particles used to determine the reaction plane)
 
 // a) Fill the common control histograms and call the method to fill fAvMultiplicity:
 this->FillCommonControlHistograms(anEvent);                                                               
 this->FillAverageMultiplicities(nRP);                                                                  
                                                                                                                                                                                                                                                                                        
 // b) Loop over data and calculate e-b-e quantities:
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
    // 1D (pt):
    // (r_{m*m,k}(pt)): 
    for(Int_t m=0;m<4;m++)
    {
     for(Int_t k=0;k<9;k++)
     {
      fReRPQ1dEBE[0][0][m][k]->Fill(dPt,pow(wPhi*wPt*wEta,k)*TMath::Cos((m+1.)*n*dPhi),1.);
      fImRPQ1dEBE[0][0][m][k]->Fill(dPt,pow(wPhi*wPt*wEta,k)*TMath::Sin((m+1.)*n*dPhi),1.);
     }
    }
           
    // s_{k}(pt) for RPs // to be improved (clarified)
    // Remark: final calculation of s_{p,k}(pt) follows after the loop over data bellow:
    for(Int_t k=0;k<9;k++)
    {
     fs1dEBE[0][0][k]->Fill(dPt,pow(wPhi*wPt*wEta,k),1.);
    }
    // 1D (eta):
    // (r_{m*m,k}(eta)): 
    for(Int_t m=0;m<4;m++)
    {
     for(Int_t k=0;k<9;k++)
     {
      fReRPQ1dEBE[0][1][m][k]->Fill(dEta,pow(wPhi*wPt*wEta,k)*TMath::Cos((m+1.)*n*dPhi),1.);
      fImRPQ1dEBE[0][1][m][k]->Fill(dEta,pow(wPhi*wPt*wEta,k)*TMath::Sin((m+1.)*n*dPhi),1.);
     }
    }   
    // s_{k}(eta) for RPs // to be improved (clarified)
    // Remark: final calculation of s_{p,k}(eta) follows after the loop over data bellow:
    for(Int_t k=0;k<9;k++)
    {
     fs1dEBE[0][1][k]->Fill(dEta,pow(wPhi*wPt*wEta,k),1.);
    }
    
    
    
    /*
    // 2D (pt,eta):
    if(fCalculate2DFlow)
    {
     // (r_{m*m,k}(pt,eta)): 
     for(Int_t m=0;m<4;m++)
     {
      for(Int_t k=0;k<9;k++)
      {
       fReRPQ2dEBE[0][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k)*TMath::Cos((m+1.)*n*dPhi),1.);
       fImRPQ2dEBE[0][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k)*TMath::Sin((m+1.)*n*dPhi),1.);
      }
     }    
     // s_{k}(pt,eta) for RPs // to be improved (clarified)
     // Remark: final calculation of s_{p,k}(pt,eta) follows after the loop over data bellow:
     for(Int_t k=0;k<9;k++)
     {
      fs2dEBE[0][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k),1.);
     }
    } // end of if(fCalculate2DFlow)  
    */ 
    
      
     
    if(aftsTrack->InPOISelection())
    {
     // 1D (pt): 
     // (q_{m*m,k}(pt)): 
     for(Int_t m=0;m<4;m++)
     {
      for(Int_t k=0;k<9;k++)
      {
       fReRPQ1dEBE[2][0][m][k]->Fill(dPt,pow(wPhi*wPt*wEta,k)*TMath::Cos((m+1.)*n*dPhi),1.);
       fImRPQ1dEBE[2][0][m][k]->Fill(dPt,pow(wPhi*wPt*wEta,k)*TMath::Sin((m+1.)*n*dPhi),1.);
      }
     } 
     // s_{k}(pt) for RP&&POIs // to be improved (clarified)
     // Remark: final calculation of s_{p,k}(pt,eta) follows after the loop over data bellow:
     for(Int_t k=0;k<9;k++)
     {
      fs1dEBE[2][0][k]->Fill(dPt,pow(wPhi*wPt*wEta,k),1.);
     }
     // 1D (eta): 
     // (q_{m*m,k}(eta)): 
     for(Int_t m=0;m<4;m++)
     {
      for(Int_t k=0;k<9;k++)
      {
       fReRPQ1dEBE[2][1][m][k]->Fill(dEta,pow(wPhi*wPt*wEta,k)*TMath::Cos((m+1.)*n*dPhi),1.);
       fImRPQ1dEBE[2][1][m][k]->Fill(dEta,pow(wPhi*wPt*wEta,k)*TMath::Sin((m+1.)*n*dPhi),1.);
      }
     } 
     // s_{k}(eta) for RP&&POIs // to be improved (clarified)
     // Remark: final calculation of s_{p,k}(pt,eta) follows after the loop over data bellow:
     for(Int_t k=0;k<9;k++)
     {
      fs1dEBE[2][1][k]->Fill(dEta,pow(wPhi*wPt*wEta,k),1.);
     }
     
     /*
     // 2D (pt,eta) 
     if(fCalculate2DFlow)
     {
      // (q_{m*m,k}(pt,eta)): 
      for(Int_t m=0;m<4;m++)
      {
       for(Int_t k=0;k<9;k++)
       {
        fReRPQ2dEBE[2][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k)*TMath::Cos((m+1.)*n*dPhi),1.);
        fImRPQ2dEBE[2][m][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k)*TMath::Sin((m+1.)*n*dPhi),1.);
       }
      } 
      // s_{k}(pt,eta) for RP&&POIs // to be improved (clarified)
      // Remark: final calculation of s_{p,k}(pt,eta) follows after the loop over data bellow:
      for(Int_t k=0;k<9;k++)
      {
       fs2dEBE[2][k]->Fill(dPt,dEta,pow(wPhi*wPt*wEta,k),1.);
      }
     } // end of if(fCalculate2DFlow) 
     */
      
    } // end of if(aftsTrack->InPOISelection())
    

     
   } // end of if(pTrack->InRPSelection())

  
  
   if(aftsTrack->InPOISelection())
   {
    dPhi = aftsTrack->Phi();
    dPt  = aftsTrack->Pt();
    dEta = aftsTrack->Eta();
    
    // 1D (pt)
    // p_n(m*n,0):   
    for(Int_t m=0;m<4;m++)
    {
     fReRPQ1dEBE[1][0][m][0]->Fill(dPt,TMath::Cos((m+1.)*n*dPhi),1.);
     fImRPQ1dEBE[1][0][m][0]->Fill(dPt,TMath::Sin((m+1.)*n*dPhi),1.);
    }
    // 1D (eta)
    // p_n(m*n,0):   
    for(Int_t m=0;m<4;m++)
    {
     fReRPQ1dEBE[1][1][m][0]->Fill(dEta,TMath::Cos((m+1.)*n*dPhi),1.);
     fImRPQ1dEBE[1][1][m][0]->Fill(dEta,TMath::Sin((m+1.)*n*dPhi),1.);
    }
    
    
    /*
    // 2D (pt,eta):
    if(fCalculate2DFlow)
    {      
     // p_n(m*n,0):   
     for(Int_t m=0;m<4;m++)
     {
      fReRPQ2dEBE[1][m][0]->Fill(dPt,dEta,TMath::Cos((m+1.)*n*dPhi),1.);
      fImRPQ2dEBE[1][m][0]->Fill(dPt,dEta,TMath::Sin((m+1.)*n*dPhi),1.);
     }
    } // end of if(fCalculate2DFlow)  
    */
    
    
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
 // integrated flow:
 if(!fEvaluateIntFlowNestedLoops)
 {
  if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
  {
   if(nRP>1) this->CalculateIntFlowCorrelations(); // without using particle weights
  } else // to if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
    {
     if(nRP>1) this->CalculateIntFlowCorrelationsUsingParticleWeights(); // with using particle weights   
    } 
       
  if(nRP>3) this->CalculateIntFlowProductOfCorrelations();
  if(nRP>1) this->CalculateIntFlowSumOfEventWeights();
  if(nRP>1) this->CalculateIntFlowSumOfProductOfEventWeights();
  if(fApplyCorrectionForNUA)
  {
   if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
   {
    if(nRP>0) this->CalculateIntFlowCorrectionsForNUASinTerms();
    if(nRP>0) this->CalculateIntFlowCorrectionsForNUACosTerms();
   } else // to if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
     {
      if(nRP>0) this->CalculateIntFlowCorrectionsForNUASinTermsUsingParticleWeights();
      if(nRP>0) this->CalculateIntFlowCorrectionsForNUACosTermsUsingParticleWeights();     
     }  
     
   if(nRP>0) this->CalculateIntFlowProductOfCorrectionTermsForNUA();     
   if(nRP>0) this->CalculateIntFlowSumOfEventWeightsNUA();     
   if(nRP>0) this->CalculateIntFlowSumOfProductOfEventWeightsNUA();     
  } // end of if(fApplyCorrectionForNUA)
 } // end of if(!fEvaluateIntFlowNestedLoops)

 // differential flow:
 if(!fEvaluateDiffFlowNestedLoops)
 {
  if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
  {
   // without using particle weights:
   this->CalculateDiffFlowCorrelations("RP","Pt"); 
   this->CalculateDiffFlowCorrelations("RP","Eta");
   this->CalculateDiffFlowCorrelations("POI","Pt");
   this->CalculateDiffFlowCorrelations("POI","Eta");
   if(fApplyCorrectionForNUA)
   {
    this->CalculateDiffFlowCorrectionsForNUASinTerms("RP","Pt");
    this->CalculateDiffFlowCorrectionsForNUASinTerms("RP","Eta");
    this->CalculateDiffFlowCorrectionsForNUASinTerms("POI","Pt");
    this->CalculateDiffFlowCorrectionsForNUASinTerms("POI","Eta");
    this->CalculateDiffFlowCorrectionsForNUACosTerms("RP","Pt");
    this->CalculateDiffFlowCorrectionsForNUACosTerms("RP","Eta");
    this->CalculateDiffFlowCorrectionsForNUACosTerms("POI","Pt");
    this->CalculateDiffFlowCorrectionsForNUACosTerms("POI","Eta");   
   } // end of if(fApplyCorrectionForNUA)  
  } else // to if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
    {
     // with using particle weights:   
     this->CalculateDiffFlowCorrelationsUsingParticleWeights("RP","Pt"); 
     this->CalculateDiffFlowCorrelationsUsingParticleWeights("RP","Eta"); 
     this->CalculateDiffFlowCorrelationsUsingParticleWeights("POI","Pt"); 
     this->CalculateDiffFlowCorrelationsUsingParticleWeights("POI","Eta"); 
     if(fApplyCorrectionForNUA)
     {
      this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("RP","Pt");
      this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("RP","Eta");
      this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("POI","Pt");
      this->CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights("POI","Eta");
      this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("RP","Pt");
      this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("RP","Eta");
      this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("POI","Pt");
      this->CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights("POI","Eta");   
     } // end of if(fApplyCorrectionForNUA)  
    } 
    
  // whether or not using particle weights the following is calculated in the same way:  
  this->CalculateDiffFlowProductOfCorrelations("RP","Pt");
  this->CalculateDiffFlowProductOfCorrelations("RP","Eta");
  this->CalculateDiffFlowProductOfCorrelations("POI","Pt");
  this->CalculateDiffFlowProductOfCorrelations("POI","Eta");
  this->CalculateDiffFlowSumOfEventWeights("RP","Pt");
  this->CalculateDiffFlowSumOfEventWeights("RP","Eta");
  this->CalculateDiffFlowSumOfEventWeights("POI","Pt");
  this->CalculateDiffFlowSumOfEventWeights("POI","Eta");
  this->CalculateDiffFlowSumOfProductOfEventWeights("RP","Pt");
  this->CalculateDiffFlowSumOfProductOfEventWeights("RP","Eta");
  this->CalculateDiffFlowSumOfProductOfEventWeights("POI","Pt");
  this->CalculateDiffFlowSumOfProductOfEventWeights("POI","Eta");   
 } // end of if(!fEvaluateDiffFlowNestedLoops)


   
  // with weights:
  // ... 
  
  /*
  // 2D differential flow
  if(fCalculate2DFlow)
  {
   // without weights:
   if(nRP>1) this->CalculateCorrelationsForDifferentialFlow2D("RP");
   if(nRP>1) this->CalculateCorrelationsForDifferentialFlow2D("POI");
  
   // with weights:
   if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
   {
    if(nRP>1) this->CalculateWeightedCorrelationsForDifferentialFlow2D("RP");
    if(nRP>1) this->CalculateWeightedCorrelationsForDifferentialFlow2D("POI");
   } 
  } // end of if(fCalculate2DFlow)
  */
  
 // distributions of correlations:
 if(fStoreDistributions)
 {
  this->StoreDistributionsOfCorrelations();
 }
  
 // d) Debugging and cross-checking (evaluate nested loops):
 //  d1) cross-checking results for integrated flow:
 if(fEvaluateIntFlowNestedLoops)
 {
  if(nPrim>0 && nPrim<=fMaxAllowedMultiplicity) // by default fMaxAllowedMultiplicity = 10 
  {
   // without using particle weights:
   if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
   {
    // correlations:
    this->CalculateIntFlowCorrelations(); // from Q-vectors
    this->EvaluateIntFlowCorrelationsWithNestedLoops(anEvent); // from nested loops (to be improved: do I have to pass here anEvent or not?)
    // correction for non-uniform acceptance:
    this->CalculateIntFlowCorrectionsForNUASinTerms(); // from Q-vectors (sin terms)
    this->CalculateIntFlowCorrectionsForNUACosTerms(); // from Q-vectors (cos terms)
    this->EvaluateIntFlowCorrectionsForNUAWithNestedLoops(anEvent); // from nested loops (both sin and cos terms)
   }
   // using particle weights:
   if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
   {
    // correlations:
    this->CalculateIntFlowCorrelationsUsingParticleWeights(); // from Q-vectors
    this->EvaluateIntFlowCorrelationsWithNestedLoopsUsingParticleWeights(anEvent); // from nested loops (to be improved: do I have to pass here anEvent or not?)
    // correction for non-uniform acceptance:
    this->CalculateIntFlowCorrectionsForNUASinTermsUsingParticleWeights(); // from Q-vectors (sin terms)
    this->CalculateIntFlowCorrectionsForNUACosTermsUsingParticleWeights(); // from Q-vectors (cos terms)
    this->EvaluateIntFlowCorrectionsForNUAWithNestedLoopsUsingParticleWeights(anEvent); // from nested loops (both sin and cos terms)   
   }
  } else if (nPrim>fMaxAllowedMultiplicity) // to if(nPrim>0 && nPrim<=fMaxAllowedMultiplicity)
    {
     cout<<endl;
     cout<<"Skipping the event because multiplicity is "<<nPrim<<". Too high to evaluate nested loops!"<<endl;
    } else
      {
       cout<<endl;
       cout<<"Skipping the event because multiplicity is "<<nPrim<<"."<<endl;      
      } 
 } // end of if(fEvaluateIntFlowNestedLoops) 
 
 //  d2) cross-checking results for differential flow:
 if(fEvaluateDiffFlowNestedLoops)
 {
  if(nPrim>0 && nPrim<=fMaxAllowedMultiplicity) // by default fMaxAllowedMultiplicity = 10
  {
   // without using particle weights:
   if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
   {
    // reduced correlations:
    // Q-vectors:
    this->CalculateDiffFlowCorrelations("RP","Pt");
    this->CalculateDiffFlowCorrelations("RP","Eta");
    this->CalculateDiffFlowCorrelations("POI","Pt");
    this->CalculateDiffFlowCorrelations("POI","Eta");
    // nested loops:
    this->EvaluateDiffFlowCorrelationsWithNestedLoops(anEvent,"RP","Pt"); 
    this->EvaluateDiffFlowCorrelationsWithNestedLoops(anEvent,"RP","Eta"); 
    this->EvaluateDiffFlowCorrelationsWithNestedLoops(anEvent,"POI","Pt"); 
    this->EvaluateDiffFlowCorrelationsWithNestedLoops(anEvent,"POI","Eta"); 
    // reduced corrections for non-uniform acceptance:
    // Q-vectors:
    this->CalculateDiffFlowCorrectionsForNUASinTerms("RP","Pt");
    this->CalculateDiffFlowCorrectionsForNUASinTerms("RP","Eta");
    this->CalculateDiffFlowCorrectionsForNUASinTerms("POI","Pt");
    this->CalculateDiffFlowCorrectionsForNUASinTerms("POI","Eta");
    this->CalculateDiffFlowCorrectionsForNUACosTerms("RP","Pt");
    this->CalculateDiffFlowCorrectionsForNUACosTerms("RP","Eta");
    this->CalculateDiffFlowCorrectionsForNUACosTerms("POI","Pt");
    this->CalculateDiffFlowCorrectionsForNUACosTerms("POI","Eta");
    // nested loops:
    this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(anEvent,"RP","Pt");
    this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(anEvent,"RP","Eta");
    this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(anEvent,"POI","Pt"); 
    this->EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(anEvent,"POI","Eta"); 
   } // end of if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
   // using particle weights:
   if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
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
   } // end of if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
  } // end of if(nPrim>0 && nPrim<=fMaxAllowedMultiplicity) // by default fMaxAllowedMultiplicity = 10
 } // end of if(fEvaluateDiffFlowNestedLoops) 
 
 // e) Reset all event by event quantities: 
 this->ResetEventByEventQuantities();
 
} // end of AliFlowAnalysisWithQCumulants::Make(AliFlowEventSimple* anEvent)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::Finish()
{
 // Calculate the final results.
 //  a) acces the constants;
 //  b) access the flags;
 //  c) calculate the final results for integrated flow (without and with weights);
 //  d) store in AliFlowCommonHistResults and print the final results for integrated flow;
 //  e) calculate the final results for differential flow (without and with weights);
 //  f) print the final results for integrated flow obtained from differential flow (to be improved (terminology));
 //  g) cross-check the results: results from Q-vectors vs results from nested loops
 
 // ******************************
 // **** ACCESS THE CONSTANTS ****
 // ******************************
 
 this->AccessConstants();          
 
 if(fCommonHists && fCommonHists->GetHarmonic())
 {
  fHarmonic = (Int_t)(fCommonHists->GetHarmonic())->GetBinContent(1); // to be improved (moved somewhere else)
 } 

 // **************************
 // **** ACCESS THE FLAGS **** // to be improved (moved somewhere else)
 // **************************    
 fUsePhiWeights = (Int_t)fUseParticleWeights->GetBinContent(1); 
 fUsePtWeights = (Int_t)fUseParticleWeights->GetBinContent(2); 
 fUseEtaWeights = (Int_t)fUseParticleWeights->GetBinContent(3);  
 fApplyCorrectionForNUA = (Int_t)fIntFlowFlags->GetBinContent(3); 
 fPrintFinalResults[0] = (Int_t)fIntFlowFlags->GetBinContent(4);
 fPrintFinalResults[1] = (Int_t)fIntFlowFlags->GetBinContent(5);
 fPrintFinalResults[2] = (Int_t)fIntFlowFlags->GetBinContent(6);
 fEvaluateIntFlowNestedLoops = (Int_t)fEvaluateNestedLoops->GetBinContent(1);
 fEvaluateDiffFlowNestedLoops = (Int_t)fEvaluateNestedLoops->GetBinContent(2); 
 fCrossCheckInPtBinNo = (Int_t)fEvaluateNestedLoops->GetBinContent(3);
 fCrossCheckInEtaBinNo = (Int_t)fEvaluateNestedLoops->GetBinContent(4); 
    
 // *********************************************************
 // **** CALCULATE THE FINAL RESULTS FOR INTEGRATED FLOW ****
 // *********************************************************     
 
 this->FinalizeCorrelationsIntFlow();
 this->CalculateCovariancesIntFlow();
 this->CalculateCumulantsIntFlow();
 this->CalculateIntFlow(); 

 if(fApplyCorrectionForNUA) // to be improved (reorganized, etc)
 {
  this->FinalizeCorrectionTermsForNUAIntFlow();
  this->CalculateCovariancesNUAIntFlow();
  this->CalculateQcumulantsCorrectedForNUAIntFlow();   
  this->CalculateIntFlowCorrectedForNUA(); 
 }
  
 // ***************************************************************
 // **** STORE AND PRINT THE FINAL RESULTS FOR INTEGRATED FLOW ****
 // ***************************************************************
 
 this->FillCommonHistResultsIntFlow();  
  
 if(fPrintFinalResults[0])
 {
  this->PrintFinalResultsForIntegratedFlow("NONAME"); // to be improved (name)
 }
 
 // ***********************************************************
 // **** CALCULATE THE FINAL RESULTS FOR DIFFERENTIAL FLOW ****
 // ***********************************************************    
 
 this->FinalizeReducedCorrelations("RP","Pt"); 
 this->FinalizeReducedCorrelations("RP","Eta"); 
 this->FinalizeReducedCorrelations("POI","Pt"); 
 this->FinalizeReducedCorrelations("POI","Eta");
 this->CalculateDiffFlowCovariances("RP","Pt");
 this->CalculateDiffFlowCovariances("RP","Eta");
 this->CalculateDiffFlowCovariances("POI","Pt");
 this->CalculateDiffFlowCovariances("POI","Eta");
 this->CalculateDiffFlowCumulants("RP","Pt");
 this->CalculateDiffFlowCumulants("RP","Eta");
 this->CalculateDiffFlowCumulants("POI","Pt");
 this->CalculateDiffFlowCumulants("POI","Eta");
 this->CalculateDiffFlow("RP","Pt");
 this->CalculateDiffFlow("RP","Eta");
 this->CalculateDiffFlow("POI","Pt");
 this->CalculateDiffFlow("POI","Eta");
 
 if(fApplyCorrectionForNUA) // to be improved (reorganized, etc)
 {
  this->FinalizeCorrectionTermsForNUADiffFlow("RP","Pt");
  this->FinalizeCorrectionTermsForNUADiffFlow("RP","Eta");
  this->FinalizeCorrectionTermsForNUADiffFlow("POI","Pt");
  this->FinalizeCorrectionTermsForNUADiffFlow("POI","Eta");      
  this->CalculateDiffFlowCumulantsCorrectedForNUA("RP","Pt");   
  this->CalculateDiffFlowCumulantsCorrectedForNUA("RP","Eta");   
  this->CalculateDiffFlowCumulantsCorrectedForNUA("POI","Pt");   
  this->CalculateDiffFlowCumulantsCorrectedForNUA("POI","Eta");  
  this->CalculateDiffFlowCorrectedForNUA("RP","Pt"); 
  this->CalculateDiffFlowCorrectedForNUA("RP","Eta"); 
  this->CalculateDiffFlowCorrectedForNUA("POI","Pt"); 
  this->CalculateDiffFlowCorrectedForNUA("POI","Eta"); 
 }
 
 this->CalculateFinalResultsForRPandPOIIntegratedFlow("RP");
 this->CalculateFinalResultsForRPandPOIIntegratedFlow("POI");

 // *****************************************************************
 // **** STORE AND PRINT THE FINAL RESULTS FOR DIFFERENTIAL FLOW ****
 // *****************************************************************
 this->FillCommonHistResultsDiffFlow("RP");
 this->FillCommonHistResultsDiffFlow("POI");

 if(fPrintFinalResults[1])
 {
  this->PrintFinalResultsForIntegratedFlow("RP"); 
 } 
 if(fPrintFinalResults[2])
 {
  this->PrintFinalResultsForIntegratedFlow("POI"); 
 } 
 // g) cross-check the results: results from Q-vectors vs results from nested loops
 
 //  g1) integrated flow:
 if(fEvaluateIntFlowNestedLoops)
 {
  this->CrossCheckIntFlowCorrelations();
  this->CrossCheckIntFlowCorrectionTermsForNUA(); 
  if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights) this->CrossCheckIntFlowExtraCorrelations();     
 } // end of if(fEvaluateIntFlowNestedLoops)  
 
 //  g2) differential flow: 
 if(fEvaluateDiffFlowNestedLoops) 
 {
  // correlations:
  this->PrintNumberOfParticlesInSelectedBin();
  this->CrossCheckDiffFlowCorrelations("RP","Pt");  
  this->CrossCheckDiffFlowCorrelations("RP","Eta"); 
  this->CrossCheckDiffFlowCorrelations("POI","Pt");  
  this->CrossCheckDiffFlowCorrelations("POI","Eta");
  // correction terms for non-uniform acceptance:
  this->CrossCheckDiffFlowCorrectionTermsForNUA("RP","Pt");      
  this->CrossCheckDiffFlowCorrectionTermsForNUA("RP","Eta");       
  this->CrossCheckDiffFlowCorrectionTermsForNUA("POI","Pt");      
  this->CrossCheckDiffFlowCorrectionTermsForNUA("POI","Eta");       
 } // end of if(fEvaluateDiffFlowNestedLoops)
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithQCumulants::Finish()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUACosTerms()
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
 } 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUACosTerms()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUASinTerms()
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
 } 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectionsForNUASinTerms()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::GetOutputHistograms(TList *outputListHistos)
{
 // a) Get pointers for common control and common result histograms and profiles.
 // b) Get pointers for histograms with particle weights.
 // c) Get pointers for histograms and profiles relevant for integrated flow.
 // d) Get pointers for histograms and profiles relevant for differental flow.
 // e) Get pointers for histograms and profiles holding results obtained with nested loops.
 
 if(outputListHistos)
 {	
  this->SetHistList(outputListHistos);
  if(!fHistList)
  {
   cout<<endl;
   cout<<" WARNING (QC): fHistList is NULL in AFAWQC::GOH() !!!!"<<endl;
   cout<<endl;
   exit(0);
  }
  this->GetPointersForCommonHistograms(); 
  this->GetPointersForParticleWeightsHistograms(); 
  this->GetPointersForIntFlowHistograms();
  this->GetPointersForDiffFlowHistograms(); 
  this->GetPointersForNestedLoopsHistograms(); 
 } else 
   {
    cout<<endl;
    cout<<" WARNING (QC): outputListHistos is NULL in AFAWQC::GOH() !!!!"<<endl;
    cout<<endl;
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
 // common control histogram (ALL events)
 TString commonHistsName = "AliFlowCommonHistQC";
 commonHistsName += fAnalysisLabel->Data();
 fCommonHists = new AliFlowCommonHist(commonHistsName.Data());
 fHistList->Add(fCommonHists);  
 // common control histogram (for events with 2 and more particles)
 TString commonHists2ndOrderName = "AliFlowCommonHist2ndOrderQC";
 commonHists2ndOrderName += fAnalysisLabel->Data();
 fCommonHists2nd = new AliFlowCommonHist(commonHists2ndOrderName.Data());
 fHistList->Add(fCommonHists2nd);  
 // common control histogram (for events with 4 and more particles)
 TString commonHists4thOrderName = "AliFlowCommonHist4thOrderQC";
 commonHists4thOrderName += fAnalysisLabel->Data();
 fCommonHists4th = new AliFlowCommonHist(commonHists4thOrderName.Data());
 fHistList->Add(fCommonHists4th);  
 // common control histogram (for events with 6 and more particles)
 TString commonHists6thOrderName = "AliFlowCommonHist6thOrderQC";
 commonHists6thOrderName += fAnalysisLabel->Data();
 fCommonHists6th = new AliFlowCommonHist(commonHists6thOrderName.Data());
 fHistList->Add(fCommonHists6th);  
 // common control histogram (for events with 8 and more particles)
 TString commonHists8thOrderName = "AliFlowCommonHist8thOrderQC";
 commonHists8thOrderName += fAnalysisLabel->Data();
 fCommonHists8th = new AliFlowCommonHist(commonHists8thOrderName.Data());
 fHistList->Add(fCommonHists8th);    
 // common histograms for final results (calculated for events with 2 and more particles)
 TString commonHistResults2ndOrderName = "AliFlowCommonHistResults2ndOrderQC";
 commonHistResults2ndOrderName += fAnalysisLabel->Data();
 fCommonHistsResults2nd = new AliFlowCommonHistResults(commonHistResults2ndOrderName.Data());
 fHistList->Add(fCommonHistsResults2nd);  
 // common histograms for final results (calculated for events with 4 and more particles)
 TString commonHistResults4thOrderName = "AliFlowCommonHistResults4thOrderQC";
 commonHistResults4thOrderName += fAnalysisLabel->Data();
 fCommonHistsResults4th = new AliFlowCommonHistResults(commonHistResults4thOrderName.Data());
 fHistList->Add(fCommonHistsResults4th); 
 // common histograms for final results (calculated for events with 6 and more particles)
 TString commonHistResults6thOrderName = "AliFlowCommonHistResults6thOrderQC";
 commonHistResults6thOrderName += fAnalysisLabel->Data();
 fCommonHistsResults6th = new AliFlowCommonHistResults(commonHistResults6thOrderName.Data());
 fHistList->Add(fCommonHistsResults6th);  
 // common histograms for final results (calculated for events with 8 and more particles)
 TString commonHistResults8thOrderName = "AliFlowCommonHistResults8thOrderQC";
 commonHistResults8thOrderName += fAnalysisLabel->Data();
 fCommonHistsResults8th = new AliFlowCommonHistResults(commonHistResults8thOrderName.Data());
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
   if(TMath::Abs(fPhiWeights->GetBinWidth(1)-fPhiBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (QC): Inconsistent binning in histograms for phi-weights throughout the code."<<endl;
    cout<<endl;
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
   if(TMath::Abs(fPtWeights->GetBinWidth(1)-fPtBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (QC): Inconsistent binning in histograms for pt-weights throughout the code."<<endl;
    cout<<endl;
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
   if(TMath::Abs(fEtaWeights->GetBinWidth(1)-fEtaBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (QC): Inconsistent binning in histograms for eta-weights throughout the code."<<endl;
    cout<<endl;
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
 // Book all objects for integrated flow:
 //  a) Book profile to hold all flags for integrated flow.
 //  b) Book event-by-event quantities.
 //  c) Book profiles. // to be improved (comment)
 //  d) Book histograms holding the final results.
 
 TString sinCosFlag[2] = {"sin","cos"}; // to be improved (should I promote this to data members?)
 TString powerFlag[2] = {"linear","quadratic"}; // to be improved (should I promote this to data members?)
 
 // a) Book profile to hold all flags for integrated flow:
 TString intFlowFlagsName = "fIntFlowFlags";
 intFlowFlagsName += fAnalysisLabel->Data();
 fIntFlowFlags = new TProfile(intFlowFlagsName.Data(),"Flags for Integrated Flow",6,0,6);
 fIntFlowFlags->SetTickLength(-0.01,"Y");
 fIntFlowFlags->SetMarkerStyle(25);
 fIntFlowFlags->SetLabelSize(0.05);
 fIntFlowFlags->SetLabelOffset(0.02,"Y");
 fIntFlowFlags->GetXaxis()->SetBinLabel(1,"Particle Weights");
 fIntFlowFlags->GetXaxis()->SetBinLabel(2,"Event Weights");
 fIntFlowFlags->GetXaxis()->SetBinLabel(3,"Corrected for NUA?");
 fIntFlowFlags->GetXaxis()->SetBinLabel(4,"Print NONAME results");
 fIntFlowFlags->GetXaxis()->SetBinLabel(5,"Print RP results");
 fIntFlowFlags->GetXaxis()->SetBinLabel(6,"Print POI results");
 fIntFlowList->Add(fIntFlowFlags);

 // b) Book event-by-event quantities:
 // Re[Q_{m*n,k}], Im[Q_{m*n,k}] and S_{p,k}^M: 
 fReQ  = new TMatrixD(4,9);
 fImQ  = new TMatrixD(4,9);
 fSMpk = new TMatrixD(8,9);
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
 fIntFlowCorrelationsAllEBE = new TH1D(intFlowCorrelationsAllEBEName.Data(),intFlowCorrelationsAllEBEName.Data(),32,0,32);
 // average correction terms for non-uniform acceptance for single event 
 // (binning is the same as in fIntFlowCorrectionTermsForNUAPro[2] and fIntFlowCorrectionTermsForNUAHist[2]):
 TString fIntFlowCorrectionTermsForNUAEBEName = "fIntFlowCorrectionTermsForNUAEBE";
 fIntFlowCorrectionTermsForNUAEBEName += fAnalysisLabel->Data();
 for(Int_t sc=0;sc<2;sc++) // sin or cos terms
 {
  fIntFlowCorrectionTermsForNUAEBE[sc] = new TH1D(Form("%s: %s terms",fIntFlowCorrectionTermsForNUAEBEName.Data(),sinCosFlag[sc].Data()),Form("Correction terms for non-uniform acceptance (%s terms)",sinCosFlag[sc].Data()),10,0,10);  
 }
 // event weights for terms for non-uniform acceptance: 
 TString fIntFlowEventWeightForCorrectionTermsForNUAEBEName = "fIntFlowEventWeightForCorrectionTermsForNUAEBE";
 fIntFlowEventWeightForCorrectionTermsForNUAEBEName += fAnalysisLabel->Data();
 for(Int_t sc=0;sc<2;sc++) // sin or cos terms
 {
  fIntFlowEventWeightForCorrectionTermsForNUAEBE[sc] = new TH1D(Form("%s: %s terms",fIntFlowEventWeightForCorrectionTermsForNUAEBEName.Data(),sinCosFlag[sc].Data()),Form("Event weights for terms for non-uniform acceptance (%s terms)",sinCosFlag[sc].Data()),10,0,10);  
 }
 // c) Book profiles: // to be improved (comment)
 // profile to hold average multiplicities and number of events for events with nRP>=0, nRP>=1, ... , and nRP>=8:
 TString avMultiplicityName = "fAvMultiplicity";
 avMultiplicityName += fAnalysisLabel->Data();
 fAvMultiplicity = new TProfile(avMultiplicityName.Data(),"Average Multiplicities of RPs",9,0,9);
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
 // average correlations <<2>>, <<4>>, <<6>> and <<8>> for all events (with wrong errors!):
 TString intFlowCorrelationsProName = "fIntFlowCorrelationsPro";
 intFlowCorrelationsProName += fAnalysisLabel->Data();
 fIntFlowCorrelationsPro = new TProfile(intFlowCorrelationsProName.Data(),"Average correlations for all events",4,0,4,"s");
 fIntFlowCorrelationsPro->SetTickLength(-0.01,"Y");
 fIntFlowCorrelationsPro->SetMarkerStyle(25);
 fIntFlowCorrelationsPro->SetLabelSize(0.06);
 fIntFlowCorrelationsPro->SetLabelOffset(0.01,"Y");
 (fIntFlowCorrelationsPro->GetXaxis())->SetBinLabel(1,"<<2>>");
 (fIntFlowCorrelationsPro->GetXaxis())->SetBinLabel(2,"<<4>>");
 (fIntFlowCorrelationsPro->GetXaxis())->SetBinLabel(3,"<<6>>");
 (fIntFlowCorrelationsPro->GetXaxis())->SetBinLabel(4,"<<8>>");
 fIntFlowProfiles->Add(fIntFlowCorrelationsPro);
 // average correlations <<2>>, <<4>>, <<6>> and <<8>> versus multiplicity for all events (error is biased estimator):
 TString correlationFlag[4] = {"<<2>>","<<4>>","<<6>>","<<8>>"};
 Int_t fMinMultiplicity = 0; // to be improved (setter to this?)
 Int_t fMaxMultiplicity = 100; // to be improved (setter to this?)
 for(Int_t ci=0;ci<4;ci++) // correlation index
 {
  TString intFlowCorrelationsVsMProName = "fIntFlowCorrelationsVsMPro";
  intFlowCorrelationsVsMProName += fAnalysisLabel->Data();
  fIntFlowCorrelationsVsMPro[ci] = new TProfile(Form("%s, %s",intFlowCorrelationsVsMProName.Data(),correlationFlag[ci].Data()),
                                                Form("%s vs multiplicity",correlationFlag[ci].Data()),
                                                fMaxMultiplicity-fMinMultiplicity,fMinMultiplicity,fMaxMultiplicity,"s");                                            
  fIntFlowCorrelationsVsMPro[ci]->GetYaxis()->SetTitle(correlationFlag[ci].Data());
  fIntFlowCorrelationsVsMPro[ci]->GetXaxis()->SetTitle("M");
  fIntFlowProfiles->Add(fIntFlowCorrelationsVsMPro[ci]);
 } // end of for(Int_t ci=0;ci<4;ci++) // correlation index  
 // averaged all correlations for all events (with wrong errors!):
 TString intFlowCorrelationsAllProName = "fIntFlowCorrelationsAllPro";
 intFlowCorrelationsAllProName += fAnalysisLabel->Data();
 fIntFlowCorrelationsAllPro = new TProfile(intFlowCorrelationsAllProName.Data(),"Average correlations for all events",32,0,32,"s");
 fIntFlowCorrelationsAllPro->SetTickLength(-0.01,"Y");
 fIntFlowCorrelationsAllPro->SetMarkerStyle(25);
 fIntFlowCorrelationsAllPro->SetLabelSize(0.03);
 fIntFlowCorrelationsAllPro->SetLabelOffset(0.01,"Y");
 // 2-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(1,"<<2>>_{n|n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(2,"<<2>>_{2n|2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(3,"<<2>>_{3n|3n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(4,"<<2>>_{4n|4n}");
 // 3-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(6,"<<3>>_{2n|n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(7,"<<3>>_{3n|2n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(8,"<<3>>_{4n|2n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(9,"<<3>>_{4n|3n,n}");
 // 4-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(11,"<<4>>_{n,n|n,n}"); 
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(12,"<<4>>_{2n,n|2n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(13,"<<4>>_{2n,2n|2n,2n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(14,"<<4>>_{3n|n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(15,"<<4>>_{3n,n|3n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(16,"<<4>>_{3n,n|2n,2n}"); 
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(17,"<<4>>_{4n|2n,n,n}");
 // 5-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(19,"<<5>>_{2n|n,n,n,n}"); 
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(20,"<<5>>_{2n,2n|2n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(21,"<<5>>_{3n,n|2n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(22,"<<5>>_{4n|n,n,n,n}");
 // 6-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(24,"<<6>>_{n,n,n|n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(25,"<<6>>_{2n,n,n|2n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(26,"<<6>>_{2n,2n|n,n,n,n}");
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(27,"<<6>>_{3n,n|n,n,n,n}");
 // 7-p correlations:  
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(29,"<<7>>_{2n,n,n|n,n,n,n}");
 // 8-p correlations:
 (fIntFlowCorrelationsAllPro->GetXaxis())->SetBinLabel(31,"<<8>>_{n,n,n,n|n,n,n,n}");
 fIntFlowProfiles->Add(fIntFlowCorrelationsAllPro);
 // when particle weights are used some extra correlations appear:
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights) 
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
 } // end of if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
 // average product of correlations <2>, <4>, <6> and <8>:  
 TString intFlowProductOfCorrelationsProName = "fIntFlowProductOfCorrelationsPro";
 intFlowProductOfCorrelationsProName += fAnalysisLabel->Data();
 fIntFlowProductOfCorrelationsPro = new TProfile(intFlowProductOfCorrelationsProName.Data(),"Average products of correlations",6,0,6);
 fIntFlowProductOfCorrelationsPro->SetTickLength(-0.01,"Y");
 fIntFlowProductOfCorrelationsPro->SetMarkerStyle(25); 
 fIntFlowProductOfCorrelationsPro->SetLabelSize(0.05);
 fIntFlowProductOfCorrelationsPro->SetLabelOffset(0.01,"Y");
 (fIntFlowProductOfCorrelationsPro->GetXaxis())->SetBinLabel(1,"<<2><4>>");
 (fIntFlowProductOfCorrelationsPro->GetXaxis())->SetBinLabel(2,"<<2><6>>");
 (fIntFlowProductOfCorrelationsPro->GetXaxis())->SetBinLabel(3,"<<2><8>>");
 (fIntFlowProductOfCorrelationsPro->GetXaxis())->SetBinLabel(4,"<<4><6>>");
 (fIntFlowProductOfCorrelationsPro->GetXaxis())->SetBinLabel(5,"<<4><8>>");
 (fIntFlowProductOfCorrelationsPro->GetXaxis())->SetBinLabel(6,"<<6><8>>");
 fIntFlowProfiles->Add(fIntFlowProductOfCorrelationsPro);
 // average product of correlations <2>, <4>, <6> and <8> versus multiplicity
 // [0=<<2><4>>,1=<<2><6>>,2=<<2><8>>,3=<<4><6>>,4=<<4><8>>,5=<<6><8>>]  
 TString intFlowProductOfCorrelationsVsMProName = "fIntFlowProductOfCorrelationsVsMPro";
 intFlowProductOfCorrelationsVsMProName += fAnalysisLabel->Data();
 TString productFlag[6] = {"<<2><4>>","<<2><6>>","<<2><8>>","<<4><6>>","<<4><8>>","<<6><8>>"};
 for(Int_t pi=0;pi<6;pi++)
 { 
  fIntFlowProductOfCorrelationsVsMPro[pi] = new TProfile(Form("%s, %s",intFlowProductOfCorrelationsVsMProName.Data(),productFlag[pi].Data()),
                                                     Form("%s versus multiplicity",productFlag[pi].Data()),
                                                     fMaxMultiplicity-fMinMultiplicity,fMinMultiplicity,fMaxMultiplicity);             
  fIntFlowProductOfCorrelationsVsMPro[pi]->GetXaxis()->SetTitle("M");
  fIntFlowProfiles->Add(fIntFlowProductOfCorrelationsVsMPro[pi]);
 } // end of for(Int_t pi=0;pi<6;pi++)
 // average product of correction terms for NUA:  
 TString intFlowProductOfCorrectionTermsForNUAProName = "fIntFlowProductOfCorrectionTermsForNUAPro";
 intFlowProductOfCorrectionTermsForNUAProName += fAnalysisLabel->Data();
 fIntFlowProductOfCorrectionTermsForNUAPro = new TProfile(intFlowProductOfCorrectionTermsForNUAProName.Data(),"Average products of correction terms for NUA",27,0,27);
 fIntFlowProductOfCorrectionTermsForNUAPro->SetTickLength(-0.01,"Y");
 fIntFlowProductOfCorrectionTermsForNUAPro->SetMarkerStyle(25); 
 fIntFlowProductOfCorrectionTermsForNUAPro->SetLabelSize(0.05);
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
  fIntFlowCorrectionTermsForNUAPro[sc] = new TProfile(Form("%s: %s terms",intFlowCorrectionTermsForNUAProName.Data(),sinCosFlag[sc].Data()),Form("Correction terms for non-uniform acceptance (%s terms)",sinCosFlag[sc].Data()),10,0,10,"s");
  fIntFlowCorrectionTermsForNUAPro[sc]->SetTickLength(-0.01,"Y");
  fIntFlowCorrectionTermsForNUAPro[sc]->SetMarkerStyle(25);
  fIntFlowCorrectionTermsForNUAPro[sc]->SetLabelSize(0.03);
  fIntFlowCorrectionTermsForNUAPro[sc]->SetLabelOffset(0.01,"Y");
  (fIntFlowCorrectionTermsForNUAPro[sc]->GetXaxis())->SetBinLabel(1,Form("<<%s(n(phi1))>>",sinCosFlag[sc].Data()));
  (fIntFlowCorrectionTermsForNUAPro[sc]->GetXaxis())->SetBinLabel(2,Form("<<%s(n(phi1+phi2))>>",sinCosFlag[sc].Data()));  
  (fIntFlowCorrectionTermsForNUAPro[sc]->GetXaxis())->SetBinLabel(3,Form("<<%s(n(phi1-phi2-phi3))>>",sinCosFlag[sc].Data()));  
  (fIntFlowCorrectionTermsForNUAPro[sc]->GetXaxis())->SetBinLabel(4,Form("<<%s(n(2phi1-phi2))>>",sinCosFlag[sc].Data()));  
  fIntFlowProfiles->Add(fIntFlowCorrectionTermsForNUAPro[sc]);
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
 (fIntFlowCorrelationsHist->GetXaxis())->SetBinLabel(1,"<<2>>");
 (fIntFlowCorrelationsHist->GetXaxis())->SetBinLabel(2,"<<4>>");
 (fIntFlowCorrelationsHist->GetXaxis())->SetBinLabel(3,"<<6>>");
 (fIntFlowCorrelationsHist->GetXaxis())->SetBinLabel(4,"<<8>>");
 fIntFlowResults->Add(fIntFlowCorrelationsHist);
 // average correlations <<2>>, <<4>>, <<6>> and <<8>> for all events (with correct errors!) vs M:
 for(Int_t ci=0;ci<4;ci++) // correlation index
 {
  TString intFlowCorrelationsVsMHistName = "fIntFlowCorrelationsVsMHist";
  intFlowCorrelationsVsMHistName += fAnalysisLabel->Data();
  fIntFlowCorrelationsVsMHist[ci] = new TH1D(Form("%s, %s",intFlowCorrelationsVsMHistName.Data(),correlationFlag[ci].Data()),
                                             Form("%s vs multiplicity",correlationFlag[ci].Data()),
                                             fMaxMultiplicity-fMinMultiplicity,fMinMultiplicity,fMaxMultiplicity);                                            
  fIntFlowCorrelationsVsMHist[ci]->GetYaxis()->SetTitle(correlationFlag[ci].Data());
  fIntFlowCorrelationsVsMHist[ci]->GetXaxis()->SetTitle("M");
  fIntFlowResults->Add(fIntFlowCorrelationsVsMHist[ci]);
 } // end of for(Int_t ci=0;ci<4;ci++) // correlation index   
 // average all correlations for all events (with correct errors!):
 TString intFlowCorrelationsAllHistName = "fIntFlowCorrelationsAllHist";
 intFlowCorrelationsAllHistName += fAnalysisLabel->Data();
 fIntFlowCorrelationsAllHist = new TH1D(intFlowCorrelationsAllHistName.Data(),"Average correlations for all events",32,0,32);
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
  fIntFlowCorrectionTermsForNUAHist[sc] = new TH1D(Form("%s: %s terms",intFlowCorrectionTermsForNUAHistName.Data(),sinCosFlag[sc].Data()),Form("Correction terms for non-uniform acceptance (%s terms)",sinCosFlag[sc].Data()),10,0,10);
  fIntFlowCorrectionTermsForNUAHist[sc]->SetTickLength(-0.01,"Y");
  fIntFlowCorrectionTermsForNUAHist[sc]->SetMarkerStyle(25);
  fIntFlowCorrectionTermsForNUAHist[sc]->SetLabelSize(0.03);
  fIntFlowCorrectionTermsForNUAHist[sc]->SetLabelOffset(0.01,"Y");
  // ......................................................................... 
  // 1-p terms:
  (fIntFlowCorrectionTermsForNUAHist[sc]->GetXaxis())->SetBinLabel(1,Form("%s(n(#phi_{1}))>",sinCosFlag[sc].Data()));
  // 2-p terms:
  // 3-p terms:
  // ...
  // ......................................................................... 
  fIntFlowResults->Add(fIntFlowCorrectionTermsForNUAHist[sc]);
 } // end of for(Int_t sc=0;sc<2;sc++) 
 // covariances (multiplied with weight dependent prefactor):
 TString intFlowCovariancesName = "fIntFlowCovariances";
 intFlowCovariancesName += fAnalysisLabel->Data();
 fIntFlowCovariances = new TH1D(intFlowCovariancesName.Data(),"Covariances (multiplied with weight dependent prefactor)",6,0,6);
 fIntFlowCovariances->SetLabelSize(0.04);
 fIntFlowCovariances->SetMarkerStyle(25);
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(1,"Cov(<2>,<4>)");
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(2,"Cov(<2>,<6>)");
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(3,"Cov(<2>,<8>)");
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(4,"Cov(<4>,<6>)");
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(5,"Cov(<4>,<8>)");
 (fIntFlowCovariances->GetXaxis())->SetBinLabel(6,"Cov(<6>,<8>)");  
 fIntFlowResults->Add(fIntFlowCovariances);
 // sum of linear and quadratic event weights for <2>, <4>, <6> and <8>:
 TString intFlowSumOfEventWeightsName = "fIntFlowSumOfEventWeights";
 intFlowSumOfEventWeightsName += fAnalysisLabel->Data();
 for(Int_t power=0;power<2;power++)
 {
  fIntFlowSumOfEventWeights[power] = new TH1D(Form("%s: %s",intFlowSumOfEventWeightsName.Data(),powerFlag[power].Data()),Form("Sum of %s event weights for correlations",powerFlag[power].Data()),4,0,4);
  fIntFlowSumOfEventWeights[power]->SetLabelSize(0.05);
  fIntFlowSumOfEventWeights[power]->SetMarkerStyle(25);
  if(power == 0)
  {
   (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(1,"#sum_{i=1}^{N} w_{<2>}");
   (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(2,"#sum_{i=1}^{N} w_{<4>}");
   (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(3,"#sum_{i=1}^{N} w_{<6>}");
   (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(4,"#sum_{i=1}^{N} w_{<8>}");
  } else if (power == 1) 
    {
     (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(1,"#sum_{i=1}^{N} w_{<2>}^{2}");
     (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(2,"#sum_{i=1}^{N} w_{<4>}^{2}");
     (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(3,"#sum_{i=1}^{N} w_{<6>}^{2}");
     (fIntFlowSumOfEventWeights[power]->GetXaxis())->SetBinLabel(4,"#sum_{i=1}^{N} w_{<8>}^{2}");
    }
  fIntFlowResults->Add(fIntFlowSumOfEventWeights[power]);
 } 
 // sum of products of event weights for correlations <2>, <4>, <6> and <8>:  
 TString intFlowSumOfProductOfEventWeightsName = "fIntFlowSumOfProductOfEventWeights";
 intFlowSumOfProductOfEventWeightsName += fAnalysisLabel->Data();
 fIntFlowSumOfProductOfEventWeights = new TH1D(intFlowSumOfProductOfEventWeightsName.Data(),"Sum of product of event weights for correlations",6,0,6);
 fIntFlowSumOfProductOfEventWeights->SetLabelSize(0.05);
 fIntFlowSumOfProductOfEventWeights->SetMarkerStyle(25);
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(1,"#sum_{i=1}^{N} w_{<2>} w_{<4>}");
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(2,"#sum_{i=1}^{N} w_{<2>} w_{<6>}");
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(3,"#sum_{i=1}^{N} w_{<2>} w_{<8>}");
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(4,"#sum_{i=1}^{N} w_{<4>} w_{<6>}");
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(5,"#sum_{i=1}^{N} w_{<4>} w_{<8>}");
 (fIntFlowSumOfProductOfEventWeights->GetXaxis())->SetBinLabel(6,"#sum_{i=1}^{N} w_{<6>} w_{<8>}");
 fIntFlowResults->Add(fIntFlowSumOfProductOfEventWeights);
 // final result for covariances of correlations (multiplied with weight dependent prefactor) versus M
 // [0=Cov(2,4),1=Cov(2,6),2=Cov(2,8),3=Cov(4,6),4=Cov(4,8),5=Cov(6,8)]:
 TString intFlowCovariancesVsMName = "fIntFlowCovariancesVsM";
 intFlowCovariancesVsMName += fAnalysisLabel->Data();
 TString covarianceFlag[6] = {"Cov(<2>,<4>)","Cov(<2>,<6>)","Cov(<2>,<8>)","Cov(<4>,<6>)","Cov(<4>,<8>)","Cov(<6>,<8>)"};
 for(Int_t ci=0;ci<6;ci++)
 {
  fIntFlowCovariancesVsM[ci] = new TH1D(Form("%s, %s",intFlowCovariancesVsMName.Data(),covarianceFlag[ci].Data()),
                                        Form("%s vs multiplicity",covarianceFlag[ci].Data()),
                                        fMaxMultiplicity-fMinMultiplicity,fMinMultiplicity,fMaxMultiplicity);
  fIntFlowCovariancesVsM[ci]->GetYaxis()->SetTitle(covarianceFlag[ci].Data());
  fIntFlowCovariancesVsM[ci]->GetXaxis()->SetTitle("M");
  fIntFlowResults->Add(fIntFlowCovariancesVsM[ci]);
 }
 // sum of linear and quadratic event weights for <2>, <4>, <6> and <8> versus multiplicity
 // [0=sum{w_{<2>}},1=sum{w_{<4>}},2=sum{w_{<6>}},3=sum{w_{<8>}}][0=linear 1,1=quadratic]:
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
                          fMaxMultiplicity-fMinMultiplicity,fMinMultiplicity,fMaxMultiplicity);    
   fIntFlowSumOfEventWeightsVsM[si][power]->GetYaxis()->SetTitle(sumFlag[power][si].Data());  
   fIntFlowSumOfEventWeightsVsM[si][power]->GetXaxis()->SetTitle("M");  
   fIntFlowResults->Add(fIntFlowSumOfEventWeightsVsM[si][power]);
  } // end of for(Int_t power=0;power<2;power++)
 } // end of for(Int_t si=0;si<4;si++)   
 // sum of products of event weights for correlations <2>, <4>, <6> and <8> vs M
 // [0=sum{w_{<2>}w_{<4>}},1=sum{w_{<2>}w_{<6>}},2=sum{w_{<2>}w_{<8>}},
 //  3=sum{w_{<4>}w_{<6>}},4=sum{w_{<4>}w_{<8>}},5=sum{w_{<6>}w_{<8>}}]:  
 TString intFlowSumOfProductOfEventWeightsVsMName = "fIntFlowSumOfProductOfEventWeightsVsM";
 intFlowSumOfProductOfEventWeightsVsMName += fAnalysisLabel->Data();
 TString sopowFlag[6] = {"#sum_{i=1}^{N} w_{<2>} w_{<4>}","#sum_{i=1}^{N} w_{<2>} w_{<6>}","#sum_{i=1}^{N} w_{<2>} w_{<8>}",
                         "#sum_{i=1}^{N} w_{<4>} w_{<6>}","#sum_{i=1}^{N} w_{<4>} w_{<8>}","#sum_{i=1}^{N} w_{<6>} w_{<8>}"}; 
 for(Int_t pi=0;pi<6;pi++)
 {
  fIntFlowSumOfProductOfEventWeightsVsM[pi] = new TH1D(Form("%s, %s",intFlowSumOfProductOfEventWeightsVsMName.Data(),sopowFlag[pi].Data()),
                  Form("%s versus multiplicity",sopowFlag[pi].Data()),
                  fMaxMultiplicity-fMinMultiplicity,fMinMultiplicity,fMaxMultiplicity); 
  fIntFlowSumOfProductOfEventWeightsVsM[pi]->GetXaxis()->SetTitle("M");
  fIntFlowSumOfProductOfEventWeightsVsM[pi]->GetYaxis()->SetTitle(sopowFlag[pi].Data()); 
  fIntFlowResults->Add(fIntFlowSumOfProductOfEventWeightsVsM[pi]);
 } // end of for(Int_t pi=0;pi<6;pi++) 
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
   fIntFlowSumOfEventWeightsNUA[sc][power] = new TH1D(Form("%s: %s, %s",intFlowSumOfEventWeightsNUAName.Data(),powerFlag[power].Data(),sinCosFlag[sc].Data()),Form("Sum of %s event weights for NUA %s terms",powerFlag[power].Data(),sinCosFlag[sc].Data()),3,0,3);
   fIntFlowSumOfEventWeightsNUA[sc][power]->SetLabelSize(0.05);
   fIntFlowSumOfEventWeightsNUA[sc][power]->SetMarkerStyle(25);
   if(power == 0)
   {
    (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(1,Form("#sum_{i=1}^{N} w_{<%s(#phi)>}",sinCosFlag[sc].Data()));
    (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(2,Form("#sum_{i=1}^{N} w_{<%s(#phi_{1}+#phi_{2})>}",sinCosFlag[sc].Data()));
    (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(3,Form("#sum_{i=1}^{N} w_{<%s(#phi_{1}-#phi_{2}-#phi_{3})>}",sinCosFlag[sc].Data()));
   } else if(power == 1) 
     {
      (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(1,Form("#sum_{i=1}^{N} w_{<%s(#phi)>}^{2}",sinCosFlag[sc].Data()));
      (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(2,Form("#sum_{i=1}^{N} w_{<%s(#phi_{1}+#phi_{2})>}^{2}",sinCosFlag[sc].Data()));
      (fIntFlowSumOfEventWeightsNUA[sc][power]->GetXaxis())->SetBinLabel(3,Form("#sum_{i=1}^{N} w_{<%s(#phi_{1}-#phi_{2}-#phi_{3})>}^{2}",sinCosFlag[sc].Data()));
     }
   fIntFlowResults->Add(fIntFlowSumOfEventWeightsNUA[sc][power]);
  }
 }  
 // sum of products of event weights for NUA terms:  
 TString intFlowSumOfProductOfEventWeightsNUAName = "fIntFlowSumOfProductOfEventWeightsNUA";
 intFlowSumOfProductOfEventWeightsNUAName += fAnalysisLabel->Data();
 fIntFlowSumOfProductOfEventWeightsNUA = new TH1D(intFlowSumOfProductOfEventWeightsNUAName.Data(),"Sum of product of event weights for NUA terms",27,0,27);
 fIntFlowSumOfProductOfEventWeightsNUA->SetLabelSize(0.05);
 fIntFlowSumOfProductOfEventWeightsNUA->SetMarkerStyle(25);
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(1,"#sum_{i=1}^{N} w_{<2>} w_{<cos(#phi)>}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(2,"#sum_{i=1}^{N} w_{<2>} w_{<sin(#phi)>}");
 (fIntFlowSumOfProductOfEventWeightsNUA->GetXaxis())->SetBinLabel(3,"#sum_{i=1}^{N} w_{<cos(#phi)>} w_{<sin(#phi)>}");
 // ....
 // to be improved - add labels for remaining bins
 // ....
 fIntFlowResults->Add(fIntFlowSumOfProductOfEventWeightsNUA);
 // final results for integrated Q-cumulants:
 TString intFlowQcumulantsName = "fIntFlowQcumulants";
 intFlowQcumulantsName += fAnalysisLabel->Data();
 fIntFlowQcumulants = new TH1D(intFlowQcumulantsName.Data(),"Integrated Q-cumulants",4,0,4);
 fIntFlowQcumulants->SetLabelSize(0.05);
 fIntFlowQcumulants->SetMarkerStyle(25);
 (fIntFlowQcumulants->GetXaxis())->SetBinLabel(1,"QC{2}");
 (fIntFlowQcumulants->GetXaxis())->SetBinLabel(2,"QC{4}");
 (fIntFlowQcumulants->GetXaxis())->SetBinLabel(3,"QC{6}");
 (fIntFlowQcumulants->GetXaxis())->SetBinLabel(4,"QC{8}");
 fIntFlowResults->Add(fIntFlowQcumulants);
 // final results for integrated Q-cumulants versus multiplicity:
 TString intFlowQcumulantsVsMName = "fIntFlowQcumulantsVsM";
 intFlowQcumulantsVsMName += fAnalysisLabel->Data();
 TString cumulantFlag[4] = {"QC{2}","QC{4}","QC{6}","QC{8}"};
 for(Int_t co=0;co<4;co++) // cumulant order
 {
  fIntFlowQcumulantsVsM[co] = new TH1D(Form("%s, %s",intFlowQcumulantsVsMName.Data(),cumulantFlag[co].Data()),
                                       Form("%s vs multipicity",cumulantFlag[co].Data()),
                                       fMaxMultiplicity-fMinMultiplicity,fMinMultiplicity,fMaxMultiplicity);
  fIntFlowQcumulantsVsM[co]->GetXaxis()->SetTitle("M");                                     
  fIntFlowQcumulantsVsM[co]->GetYaxis()->SetTitle(cumulantFlag[co].Data());  
  fIntFlowResults->Add(fIntFlowQcumulantsVsM[co]);                                    
 } // end of for(Int_t co=0;co<4;co++) // cumulant order
 // final integrated flow estimates from Q-cumulants:
 TString intFlowName = "fIntFlow";
 intFlowName += fAnalysisLabel->Data();  
 // integrated flow from Q-cumulants:
 fIntFlow = new TH1D(intFlowName.Data(),"Integrated flow estimates from Q-cumulants",4,0,4);
 fIntFlow->SetLabelSize(0.05);
 fIntFlow->SetMarkerStyle(25);
 (fIntFlow->GetXaxis())->SetBinLabel(1,"v_{2}{2,QC}"); // to be improved (harwired harmonic)
 (fIntFlow->GetXaxis())->SetBinLabel(2,"v_{2}{4,QC}"); // to be improved (harwired harmonic)
 (fIntFlow->GetXaxis())->SetBinLabel(3,"v_{2}{6,QC}"); // to be improved (harwired harmonic)
 (fIntFlow->GetXaxis())->SetBinLabel(4,"v_{2}{8,QC}"); // to be improved (harwired harmonic)
 fIntFlowResults->Add(fIntFlow); 
 // integrated flow from Q-cumulants: versus multiplicity:
 TString intFlowVsMName = "fIntFlowVsM";
 intFlowVsMName += fAnalysisLabel->Data();
 TString flowFlag[4] = {"v_{2}{2,QC}","v_{2}{4,QC}","v_{2}{6,QC}","v_{2}{8,QC}"}; // to be improved (harwired harmonic)
 for(Int_t co=0;co<4;co++) // cumulant order
 {
  fIntFlowVsM[co] = new TH1D(Form("%s, %s",intFlowVsMName.Data(),flowFlag[co].Data()),
                                       Form("%s vs multipicity",flowFlag[co].Data()),
                                       fMaxMultiplicity-fMinMultiplicity,fMinMultiplicity,fMaxMultiplicity);
  fIntFlowVsM[co]->GetXaxis()->SetTitle("M");                                     
  fIntFlowVsM[co]->GetYaxis()->SetTitle(flowFlag[co].Data());  
  fIntFlowResults->Add(fIntFlowVsM[co]);                                    
 } // end of for(Int_t co=0;co<4;co++) // cumulant order

 /* // to be improved (removed):
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
  fIntFlowDirectCorrelations = new TProfile(intFlowDirectCorrelationsName.Data(),"Multiparticle correlations calculated with nested loops (for int. flow)",32,0,32,"s");
  fNestedLoopsList->Add(fIntFlowDirectCorrelations);
  if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
  {
   TString intFlowExtraDirectCorrelationsName = "fIntFlowExtraDirectCorrelations";
   intFlowExtraDirectCorrelationsName += fAnalysisLabel->Data();
   fIntFlowExtraDirectCorrelations = new TProfile(intFlowExtraDirectCorrelationsName.Data(),"Extra multiparticle correlations calculated with nested loops (for int. flow)",100,0,100,"s");
   fNestedLoopsList->Add(fIntFlowExtraDirectCorrelations);  
  } // end of if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
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
   for(Int_t pe=0;pe<2;pe++) // pt or eta
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
   for(Int_t pe=0;pe<2;pe++) // pt or eta
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
 // Remark 1: multi-particle correlations calculated with non-weighted Q-vectors are stored in 1D profile fQCorrelations[0]. // to be improved (wrong profiles)
 // Remark 2: binning of fQCorrelations[0] is organized as follows: // to be improved (wrong profiles)
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
  
  // average 2-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(1,two1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(2,two2n2n);
  fIntFlowCorrelationsAllEBE->SetBinContent(3,two3n3n);
  fIntFlowCorrelationsAllEBE->SetBinContent(4,two4n4n);
          
  // average 2-particle correlations for all events:      
  fIntFlowCorrelationsAllPro->Fill(0.5,two1n1n,dMult*(dMult-1.));  
  fIntFlowCorrelationsAllPro->Fill(1.5,two2n2n,dMult*(dMult-1.)); 
  fIntFlowCorrelationsAllPro->Fill(2.5,two3n3n,dMult*(dMult-1.)); 
  fIntFlowCorrelationsAllPro->Fill(3.5,two4n4n,dMult*(dMult-1.)); 
  
  // store separetately <2> (to be improved: do I really need this?)
  fIntFlowCorrelationsEBE->SetBinContent(1,two1n1n); // <2>
  
  // to be improved (this can be implemented better):
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
  fIntFlowCorrelationsVsMPro[0]->Fill(dMult+0.5,two1n1n,mWeight2p); 
    
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
              
  // average 3-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(6,three2n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(7,three3n2n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(8,three4n2n2n);
  fIntFlowCorrelationsAllEBE->SetBinContent(9,three4n3n1n);
        
  // average 3-particle correlations for all events:                
  fIntFlowCorrelationsAllPro->Fill(5.5,three2n1n1n,dMult*(dMult-1.)*(dMult-2.)); 
  fIntFlowCorrelationsAllPro->Fill(6.5,three3n2n1n,dMult*(dMult-1.)*(dMult-2.));
  fIntFlowCorrelationsAllPro->Fill(7.5,three4n2n2n,dMult*(dMult-1.)*(dMult-2.)); 
  fIntFlowCorrelationsAllPro->Fill(8.5,three4n3n1n,dMult*(dMult-1.)*(dMult-2.));    
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
               
  // average 4-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(11,four1n1n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(12,four2n1n2n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(13,four2n2n2n2n);
  fIntFlowCorrelationsAllEBE->SetBinContent(14,four3n1n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(15,four3n1n3n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(16,four3n1n2n2n);
  fIntFlowCorrelationsAllEBE->SetBinContent(17,four4n2n1n1n);
        
  // average 4-particle correlations for all events:                
  fIntFlowCorrelationsAllPro->Fill(10.5,four1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(11.5,four2n1n2n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(12.5,four2n2n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(13.5,four3n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(14.5,four3n1n3n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));
  fIntFlowCorrelationsAllPro->Fill(15.5,four3n1n2n2n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.));  
  fIntFlowCorrelationsAllPro->Fill(16.5,four4n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)); 
  
  // store separetately <4> (to be improved: do I really need this?)
  fIntFlowCorrelationsEBE->SetBinContent(2,four1n1n1n1n); // <4>
  
  // to be improved (this can be implemented better):
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
  fIntFlowCorrelationsVsMPro[1]->Fill(dMult+0.5,four1n1n1n1n,mWeight4p); 
  
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
                 
  // average 5-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(19,five2n1n1n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(20,five2n2n2n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(21,five3n1n2n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(22,five4n1n1n1n1n);
        
  // average 5-particle correlations for all events:                         
  fIntFlowCorrelationsAllPro->Fill(18.5,five2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)); 
  fIntFlowCorrelationsAllPro->Fill(19.5,five2n2n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(20.5,five3n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  fIntFlowCorrelationsAllPro->Fill(21.5,five4n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
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
                                 
  // average 6-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(24,six1n1n1n1n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(25,six2n1n1n2n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(26,six2n2n1n1n1n1n);
  fIntFlowCorrelationsAllEBE->SetBinContent(27,six3n1n1n1n1n1n);
        
  // average 6-particle correlations for all events:         
  fIntFlowCorrelationsAllPro->Fill(23.5,six1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  fIntFlowCorrelationsAllPro->Fill(24.5,six2n1n1n2n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 
  fIntFlowCorrelationsAllPro->Fill(25.5,six2n2n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.));
  fIntFlowCorrelationsAllPro->Fill(26.5,six3n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)); 

  // store separetately <6> (to be improved: do I really need this?)
  fIntFlowCorrelationsEBE->SetBinContent(3,six1n1n1n1n1n1n); // <6>
  
  // to be improved (this can be implemented better):
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
  fIntFlowCorrelationsVsMPro[2]->Fill(dMult+0.5,six1n1n1n1n1n1n,mWeight6p);   
 
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
        
  // average 7-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(29,seven2n1n1n1n1n1n1n);
       
  // average 7-particle correlations for all events:                      
  fIntFlowCorrelationsAllPro->Fill(28.5,seven2n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.));
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
  
  // average 8-particle correlations for single event: 
  fIntFlowCorrelationsAllEBE->SetBinContent(31,eight1n1n1n1n1n1n1n1n);
       
  // average 8-particle correlations for all events:                       
  fIntFlowCorrelationsAllPro->Fill(30.5,eight1n1n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
 
  // store separetately <8> (to be improved: do I really need this?)
  fIntFlowCorrelationsEBE->SetBinContent(4,eight1n1n1n1n1n1n1n1n); // <8>
  
  // to be improved (this can be implemented better):
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
  fIntFlowCorrelationsVsMPro[3]->Fill(dMult+0.5,eight1n1n1n1n1n1n1n1n,mWeight8p);    
  
  // distribution of <cos(n*(phi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8))>
  //f8pDistribution->Fill(eight1n1n1n1n1n1n1n1n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.)*(dMult-5.)*(dMult-6.)*(dMult-7.));
 } // end of if(dMult>7) 
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrelations()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateIntFlowProductOfCorrelations()
{
 // Calculate averages of products of correlations for integrated flow.
 
 // multiplicity:
 Double_t dMult = (*fSMpk)(0,0);
 
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
   fIntFlowProductOfCorrelationsVsMPro[counter]->Fill(dMult+0.5,
                                                      fIntFlowCorrelationsEBE->GetBinContent(ci1)*
                                                      fIntFlowCorrelationsEBE->GetBinContent(ci2),
                                                      fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci1)*
                                                      fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci2));
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
    
 for(Int_t power=0;power<2;power++)
 { 
  if(!(fIntFlowCorrelationsPro && fIntFlowProductOfCorrelationsPro 
       && fIntFlowSumOfEventWeights[power] && fIntFlowSumOfProductOfEventWeights
       && fIntFlowCovariances)) 
  {
   cout<<"WARNING: fIntFlowCorrelationsPro && fIntFlowProductOfCorrelationsPro "<<endl;
   cout<<"         && fIntFlowSumOfEventWeights[power] && fIntFlowSumOfProductOfEventWeights"<<endl;
   cout<<"         && fIntFlowCovariances is NULL in AFAWQC::FCIF() !!!!"<<endl;
   cout<<"power = "<<power<<endl;
   exit(0);
  }
 }
   
 // average 2-, 4-, 6- and 8-particle correlations for all events:
 Double_t correlation[4] = {0.};
 for(Int_t ci=0;ci<4;ci++)
 {
  correlation[ci] = fIntFlowCorrelationsPro->GetBinContent(ci+1);
 } 
 // average products of 2-, 4-, 6- and 8-particle correlations: 
 Double_t productOfCorrelations[4][4] = {{0.}};
 Int_t productOfCorrelationsLabel = 1;
 // denominators in the expressions for the unbiased estimator for covariance:
 Double_t denominator[4][4] = {{0.}};
 Int_t sumOfProductOfEventWeightsLabel1 = 1;
 // weight dependent prefactor which multiply unbiased estimators for covariances:
 Double_t wPrefactor[4][4] = {{0.}}; 
 Int_t sumOfProductOfEventWeightsLabel2 = 1;
 for(Int_t c1=0;c1<4;c1++)
 {
  for(Int_t c2=c1+1;c2<4;c2++)
  {
   productOfCorrelations[c1][c2] = fIntFlowProductOfCorrelationsPro->GetBinContent(productOfCorrelationsLabel);
   if(fIntFlowSumOfEventWeights[0]->GetBinContent(c1+1) && fIntFlowSumOfEventWeights[0]->GetBinContent(c2+1))
   {
    denominator[c1][c2] = 1.-(fIntFlowSumOfProductOfEventWeights->GetBinContent(sumOfProductOfEventWeightsLabel1))/
                             (fIntFlowSumOfEventWeights[0]->GetBinContent(c1+1) 
                              * fIntFlowSumOfEventWeights[0]->GetBinContent(c2+1));
                              
    wPrefactor[c1][c2] =  fIntFlowSumOfProductOfEventWeights->GetBinContent(sumOfProductOfEventWeightsLabel2)/ 
                          (fIntFlowSumOfEventWeights[0]->GetBinContent(c1+1)
                            * fIntFlowSumOfEventWeights[0]->GetBinContent(c2+1));
                          
                              
   }
   productOfCorrelationsLabel++;
   sumOfProductOfEventWeightsLabel1++;
   sumOfProductOfEventWeightsLabel2++;  
  }
 }
 
 // covariance label:
 Int_t covarianceLabel = 1;
 for(Int_t c1=0;c1<4;c1++)
 {
  for(Int_t c2=c1+1;c2<4;c2++)
  {
   if(denominator[c1][c2])
   {
    // covariances:
    Double_t cov = (productOfCorrelations[c1][c2]-correlation[c1]*correlation[c2])/denominator[c1][c2]; 
    // covarianced multiplied with weight dependent prefactor:
    Double_t wCov = cov * wPrefactor[c1][c2];
    fIntFlowCovariances->SetBinContent(covarianceLabel,wCov);
   }
   covarianceLabel++;
  }
 }
 
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
 Double_t denominator1 = 1.-sumOfWW1/(sumOfW1st1*sumOfW2nd1);
 // covariance:
 Double_t covariance1 = numerator1/denominator1;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor1 = sumOfWW1/(sumOfW1st1*sumOfW2nd1);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(1,wPrefactor1*covariance1);
  
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
 Double_t denominator2 = 1.-sumOfWW2/(sumOfW1st2*sumOfW2nd2);
 // covariance:
 Double_t covariance2 = numerator2/denominator2;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor2 = sumOfWW2/(sumOfW1st2*sumOfW2nd2);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(2,wPrefactor2*covariance2);
 
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
 Double_t denominator3 = 1.-sumOfWW3/(sumOfW1st3*sumOfW2nd3);
 // covariance:
 Double_t covariance3 = numerator3/denominator3;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor3 = sumOfWW3/(sumOfW1st3*sumOfW2nd3);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(3,wPrefactor3*covariance3);
 
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
 Double_t denominator4 = 1.-sumOfWW4/(sumOfW1st4*sumOfW2nd4);
 // covariance:
 Double_t covariance4 = numerator4/denominator4;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor4 = sumOfWW4/(sumOfW1st4*sumOfW2nd4);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(4,wPrefactor4*covariance4);

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
 Double_t denominator5 = 1.-sumOfWW5/(sumOfW1st5*sumOfW2nd5);
 // covariance:
 Double_t covariance5 = numerator5/denominator5;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor5 = sumOfWW5/(sumOfW1st5*sumOfW2nd5);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(5,wPrefactor5*covariance5);

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
 Double_t denominator6 = 1.-sumOfWW6/(sumOfW1st6*sumOfW2nd6);
 // covariance:
 Double_t covariance6 = numerator6/denominator6;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor6 = sumOfWW6/(sumOfW1st6*sumOfW2nd6);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(6,wPrefactor6*covariance6);

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
 Double_t denominator7 = 1.-sumOfWW7/(sumOfW1st7*sumOfW2nd7);
 // covariance:
 Double_t covariance7 = numerator7/denominator7;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor7 = sumOfWW7/(sumOfW1st7*sumOfW2nd7);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(7,wPrefactor7*covariance7);

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
 Double_t denominator8 = 1.-sumOfWW8/(sumOfW1st8*sumOfW2nd8);
 // covariance:
 Double_t covariance8 = numerator8/denominator8;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor8 = sumOfWW8/(sumOfW1st8*sumOfW2nd8);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(8,wPrefactor8*covariance8);

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
 Double_t denominator9 = 1.-sumOfWW9/(sumOfW1st9*sumOfW2nd9);
 // covariance:
 Double_t covariance9 = numerator9/denominator9;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor9 = sumOfWW9/(sumOfW1st9*sumOfW2nd9);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(9,wPrefactor9*covariance9);

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
 Double_t denominator10 = 1.-sumOfWW10/(sumOfW1st10*sumOfW2nd10);
 // covariance:
 Double_t covariance10 = numerator10/denominator10;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor10 = sumOfWW10/(sumOfW1st10*sumOfW2nd10);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(10,wPrefactor10*covariance10);

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
 Double_t denominator11 = 1.-sumOfWW11/(sumOfW1st11*sumOfW2nd11);
 // covariance:
 Double_t covariance11 = numerator11/denominator11;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor11 = sumOfWW11/(sumOfW1st11*sumOfW2nd11);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(11,wPrefactor11*covariance11);

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
 Double_t denominator12 = 1.-sumOfWW12/(sumOfW1st12*sumOfW2nd12);
 // covariance:
 Double_t covariance12 = numerator12/denominator12;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor12 = sumOfWW12/(sumOfW1st12*sumOfW2nd12);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(12,wPrefactor12*covariance12);

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
 Double_t denominator13 = 1.-sumOfWW13/(sumOfW1st13*sumOfW2nd13);
 // covariance:
 Double_t covariance13 = numerator13/denominator13;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor13 = sumOfWW13/(sumOfW1st13*sumOfW2nd13);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(13,wPrefactor13*covariance13);

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
 Double_t denominator14 = 1.-sumOfWW14/(sumOfW1st14*sumOfW2nd14);
 // covariance:
 Double_t covariance14 = numerator14/denominator14;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor14 = sumOfWW14/(sumOfW1st14*sumOfW2nd14);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(14,wPrefactor14*covariance14);

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
 Double_t denominator15 = 1.-sumOfWW15/(sumOfW1st15*sumOfW2nd15);
 // covariance:
 Double_t covariance15 = numerator15/denominator15;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor15 = sumOfWW15/(sumOfW1st15*sumOfW2nd15);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(15,wPrefactor15*covariance15);

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
 Double_t denominator16 = 1.-sumOfWW16/(sumOfW1st16*sumOfW2nd16);
 // covariance:
 Double_t covariance16 = numerator16/denominator16;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor16 = sumOfWW16/(sumOfW1st16*sumOfW2nd16);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(16,wPrefactor16*covariance16);

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
 Double_t denominator17 = 1.-sumOfWW17/(sumOfW1st17*sumOfW2nd17);
 // covariance:
 Double_t covariance17 = numerator17/denominator17;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor17 = sumOfWW17/(sumOfW1st17*sumOfW2nd17);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(17,wPrefactor17*covariance17);

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
 Double_t denominator18 = 1.-sumOfWW18/(sumOfW1st18*sumOfW2nd18);
 // covariance:
 Double_t covariance18 = numerator18/denominator18;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor18 = sumOfWW18/(sumOfW1st18*sumOfW2nd18);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(18,wPrefactor18*covariance18);

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
 Double_t denominator19 = 1.-sumOfWW19/(sumOfW1st19*sumOfW2nd19);
 // covariance:
 Double_t covariance19 = numerator19/denominator19;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor19 = sumOfWW19/(sumOfW1st19*sumOfW2nd19);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(19,wPrefactor19*covariance19);

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
 Double_t denominator20 = 1.-sumOfWW20/(sumOfW1st20*sumOfW2nd20);
 // covariance:
 Double_t covariance20 = numerator20/denominator20;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor20 = sumOfWW20/(sumOfW1st20*sumOfW2nd20);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(20,wPrefactor20*covariance20);

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
 Double_t denominator21 = 1.-sumOfWW21/(sumOfW1st21*sumOfW2nd21);
 // covariance:
 Double_t covariance21 = numerator21/denominator21;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor21 = sumOfWW21/(sumOfW1st21*sumOfW2nd21);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(21,wPrefactor21*covariance21);

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
 Double_t denominator22 = 1.-sumOfWW22/(sumOfW1st22*sumOfW2nd22);
 // covariance:
 Double_t covariance22 = numerator22/denominator22;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor22 = sumOfWW22/(sumOfW1st22*sumOfW2nd22);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(22,wPrefactor22*covariance22);

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
 Double_t denominator23 = 1.-sumOfWW23/(sumOfW1st23*sumOfW2nd23);
 // covariance:
 Double_t covariance23 = numerator23/denominator23;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor23 = sumOfWW23/(sumOfW1st23*sumOfW2nd23);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(23,wPrefactor23*covariance23);

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
 Double_t denominator24 = 1.-sumOfWW24/(sumOfW1st24*sumOfW2nd24);
 // covariance:
 Double_t covariance24 = numerator24/denominator24;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor24 = sumOfWW24/(sumOfW1st24*sumOfW2nd24);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(24,wPrefactor24*covariance24);

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
 Double_t denominator25 = 1.-sumOfWW25/(sumOfW1st25*sumOfW2nd25);
 // covariance:
 Double_t covariance25 = numerator25/denominator25;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor25 = sumOfWW25/(sumOfW1st25*sumOfW2nd25);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(25,wPrefactor25*covariance25);

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
 Double_t denominator26 = 1.-sumOfWW26/(sumOfW1st26*sumOfW2nd26);
 // covariance:
 Double_t covariance26 = numerator26/denominator26;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor26 = sumOfWW26/(sumOfW1st26*sumOfW2nd26);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(26,wPrefactor26*covariance26);

 // Cov(<cos(phi1-phi2-phi3)>,<sin(phi1-phi2-phi3)>):
 Double_t product27 = fIntFlowProductOfCorrectionTermsForNUAPro->GetBinContent(27); // <<cos(phi1-phi2-phi3)><sin(phi1-phi2-phi3)>> 
 Double_t term1st27 = fIntFlowCorrectionTermsForNUAPro[1]->GetBinContent(2); // <<cos(phi1-phi2-phi3)>>
 Double_t term2nd27 = fIntFlowCorrectionTermsForNUAPro[0]->GetBinContent(3); // <<sin(phi1-phi2-phi3)>>
 Double_t sumOfW1st27 = fIntFlowSumOfEventWeightsNUA[1][0]->GetBinContent(2); // W_{<cos(phi1-phi2-phi3)>}
 Double_t sumOfW2nd27 = fIntFlowSumOfEventWeightsNUA[0][0]->GetBinContent(3); // W_{<sin(phi1-phi2-phi3)>}
 Double_t sumOfWW27 = fIntFlowSumOfProductOfEventWeightsNUA->GetBinContent(27); // W_{<cos(phi1-phi2-phi3)>} * W_{<sin(phi1-phi2-phi3)>}
 // numerator in the expression for the the unbiased estimator for covariance:
 Double_t numerator27 = product27 - term1st27*term2nd27; 
 // denominator in the expression for the the unbiased estimator for covariance:
 Double_t denominator27 = 1.-sumOfWW27/(sumOfW1st27*sumOfW2nd27);
 // covariance:
 Double_t covariance27 = numerator27/denominator27;
 // weight dependent prefactor for covariance:
 Double_t wPrefactor27 = sumOfWW27/(sumOfW1st27*sumOfW2nd27);
 // finally, store "weighted" covariance:
 fIntFlowCovariancesNUA->SetBinContent(27,wPrefactor27*covariance27);

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
 
 for(Int_t power=0;power<2;power++)
 { 
  if(!(fIntFlowCorrelationsHist && fIntFlowCorrelationsPro && fIntFlowSumOfEventWeights[power])) 
  {
   cout<<"WARNING: fIntFlowCorrelationsHist && fIntFlowCorrelationsPro && fIntFlowSumOfEventWeights[power] is NULL in AFAWQC::FCIF() !!!!"<<endl;
   cout<<"power = "<<power<<endl;
   exit(0);
  }
 }
  
 for(Int_t ci=1;ci<=4;ci++) // correlation index
 {
  Double_t correlation = fIntFlowCorrelationsPro->GetBinContent(ci);
  Double_t spread = fIntFlowCorrelationsPro->GetBinError(ci);
  Double_t sumOfLinearEventWeights = fIntFlowSumOfEventWeights[0]->GetBinContent(ci);
  Double_t sumOfQuadraticEventWeights = fIntFlowSumOfEventWeights[1]->GetBinContent(ci);
  Double_t termA = 0.;
  Double_t termB = 0.;
  if(sumOfLinearEventWeights)
  {
   termA = pow(sumOfQuadraticEventWeights,0.5)/sumOfLinearEventWeights;
  } else
    {
     cout<<"WARNING: sumOfLinearEventWeights == 0 in AFAWQC::FCIF() !!!!"<<endl;
     cout<<"         (for "<<2*ci<<"-particle correlation)"<<endl;
    }
  if(1.-pow(termA,2.) > 0.)
  {
   termB = 1./pow(1-pow(termA,2.),0.5);
  } else
    {
     cout<<"WARNING: 1.-pow(termA,2.) <= 0 in AFAWQC::FCIF() !!!!"<<endl;   
     cout<<"         (for "<<2*ci<<"-particle correlation)"<<endl;
    }     
  Double_t statisticalError = termA * spread * termB;
  fIntFlowCorrelationsHist->SetBinContent(ci,correlation);
  fIntFlowCorrelationsHist->SetBinError(ci,statisticalError);
 } // end of for(Int_t ci=1;ci<=4;ci++) // correlation index     
 
 // versus multiplicity: 
 for(Int_t ci=0;ci<=3;ci++) // correlation index
 {
  Int_t nBins = fIntFlowCorrelationsVsMPro[ci]->GetNbinsX(); 
  for(Int_t b=1;b<=nBins;b++) // looping over multiplicity bins
  {
   Double_t correlationVsM = fIntFlowCorrelationsVsMPro[ci]->GetBinContent(b);
   Double_t spreadVsM = fIntFlowCorrelationsVsMPro[ci]->GetBinError(b);
   Double_t sumOfLinearEventWeightsVsM = fIntFlowSumOfEventWeightsVsM[ci][0]->GetBinContent(b);
   Double_t sumOfQuadraticEventWeightsVsM = fIntFlowSumOfEventWeightsVsM[ci][1]->GetBinContent(b);
   Double_t termAVsM = 0.;
   Double_t termBVsM = 0.;
   if(sumOfLinearEventWeightsVsM)
   {
    termAVsM = pow(sumOfQuadraticEventWeightsVsM,0.5)/sumOfLinearEventWeightsVsM;
   } else
     {
      cout<<"WARNING: sumOfLinearEventWeightsVsM == 0 in AFAWQC::FCIF() !!!!"<<endl;
      cout<<"         (for "<<2*ci<<"-particle correlation versus multiplicity)"<<endl;
     }
   if(1.-pow(termAVsM,2.) > 0.)
   {
    termBVsM = 1./pow(1-pow(termAVsM,2.),0.5);
   } else
     {
      cout<<"WARNING: 1.-pow(termAVsM,2.) <= 0 in AFAWQC::FCIF() !!!!"<<endl;   
      cout<<"         (for "<<2*ci<<"-particle correlation versus multiplicity)"<<endl;
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


void AliFlowAnalysisWithQCumulants::CalculateCumulantsIntFlow()
{
 // a) Calculate Q-cumulants from the measured multiparticle correlations.
 // b) Propagate the statistical errors of measured multiparticle correlations to statistical errors of Q-cumulants.  
 // c) REMARK: Q-cumulants calculated in this method are biased by non-uniform acceptance of detector !!!! 
 //            Method ApplyCorrectionForNonUniformAcceptance* (to be improved: finalize the name here)
 //            is called afterwards to correct for this bias.   
 // d) Store the results and statistical error of Q-cumulants in histogram fCumulants.
 //    Binning of fCumulants is organized as follows:
 //
 //     1st bin: QC{2}
 //     2nd bin: QC{4}
 //     3rd bin: QC{6}
 //     4th bin: QC{8}
 
 if(!(fIntFlowCorrelationsHist && fIntFlowCovariances && fIntFlowQcumulants))
 {
  cout<<"WARNING: fIntFlowCorrelationsHist && fIntFlowCovariances && fIntFlowQcumulants is NULL in AFAWQC::CCIF() !!!!"<<endl;
  exit(0);
 }
 
 // correlations:
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>> 
 Double_t four = fIntFlowCorrelationsHist->GetBinContent(2); // <<4>>  
 Double_t six = fIntFlowCorrelationsHist->GetBinContent(3); // <<6>> 
 Double_t eight = fIntFlowCorrelationsHist->GetBinContent(4); // <<8>>  
 
 // statistical errors of average 2-, 4-, 6- and 8-particle azimuthal correlations:
 Double_t twoError = fIntFlowCorrelationsHist->GetBinError(1); // statistical error of <2>  
 Double_t fourError = fIntFlowCorrelationsHist->GetBinError(2); // statistical error of <4>   
 Double_t sixError = fIntFlowCorrelationsHist->GetBinError(3); // statistical error of <6> 
 Double_t eightError = fIntFlowCorrelationsHist->GetBinError(4); // statistical error of <8> 
 
 // covariances (multiplied by prefactor depending on weights - see comments in CalculateCovariancesIntFlow()):
 Double_t wCov24 = fIntFlowCovariances->GetBinContent(1); // Cov(<2>,<4>) * prefactor(w_<2>,w_<4>)
 Double_t wCov26 = fIntFlowCovariances->GetBinContent(2); // Cov(<2>,<6>) * prefactor(w_<2>,w_<6>)
 Double_t wCov28 = fIntFlowCovariances->GetBinContent(3); // Cov(<2>,<8>) * prefactor(w_<2>,w_<8>)
 Double_t wCov46 = fIntFlowCovariances->GetBinContent(4); // Cov(<4>,<6>) * prefactor(w_<4>,w_<6>)
 Double_t wCov48 = fIntFlowCovariances->GetBinContent(5); // Cov(<4>,<8>) * prefactor(w_<4>,w_<8>)
 Double_t wCov68 = fIntFlowCovariances->GetBinContent(6); // Cov(<6>,<8>) * prefactor(w_<6>,w_<8>)
 
 // Q-cumulants: 
 Double_t qc2 = 0.; // QC{2}
 Double_t qc4 = 0.; // QC{4}
 Double_t qc6 = 0.; // QC{6}
 Double_t qc8 = 0.; // QC{8}
 if(two) qc2 = two; 
 if(four) qc4 = four-2.*pow(two,2.); 
 if(six) qc6 = six-9.*two*four+12.*pow(two,3.); 
 if(eight) qc8 = eight-16.*two*six-18.*pow(four,2.)+144.*pow(two,2.)*four-144.*pow(two,4.); 
 
 // statistical errors of Q-cumulants:       
 Double_t qc2Error = 0.;
 Double_t qc4Error = 0.;
 Double_t qc6Error = 0.;
 Double_t qc8Error = 0.;
 
 // squared statistical errors of Q-cumulants:       
 //Double_t qc2ErrorSquared = 0.;
 Double_t qc4ErrorSquared = 0.;
 Double_t qc6ErrorSquared = 0.;
 Double_t qc8ErrorSquared = 0.;
        
 // statistical error of QC{2}:              
 qc2Error = twoError;                     
                             
 // statistical error of QC{4}:              
 qc4ErrorSquared = 16.*pow(two,2.)*pow(twoError,2)+pow(fourError,2.)
                 - 8.*two*wCov24;                     
 if(qc4ErrorSquared>0.)
 {
  qc4Error = pow(qc4ErrorSquared,0.5);
 } else 
   {
    cout<<"WARNING: Statistical error of QC{4} is imaginary !!!!"<<endl;
   }
                                           
 // statistical error of QC{6}:              
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
    cout<<"WARNING: Statistical error of QC{6} is imaginary !!!!"<<endl;
   }
                            
 // statistical error of QC{8}:              
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
    cout<<"WARNING: Statistical error of QC{8} is imaginary !!!!"<<endl;
   }

 // store the results and statistical errors for Q-cumulants:
 fIntFlowQcumulants->SetBinContent(1,qc2);
 fIntFlowQcumulants->SetBinError(1,qc2Error);
 fIntFlowQcumulants->SetBinContent(2,qc4);
 fIntFlowQcumulants->SetBinError(2,qc4Error);
 fIntFlowQcumulants->SetBinContent(3,qc6);
 fIntFlowQcumulants->SetBinError(3,qc6Error);
 fIntFlowQcumulants->SetBinContent(4,qc8); 
 fIntFlowQcumulants->SetBinError(4,qc8Error); 
  
} // end of AliFlowAnalysisWithQCumulants::CalculateCumulantsIntFlow()


//================================================================================================================================ 


void AliFlowAnalysisWithQCumulants::CalculateIntFlow()
{
 // a) Calculate the final results for reference flow estimates from Q-cumulants.
 // b) Propagate the statistical errors of measured multiparticle correlations to statistical errors of reference flow estimates.  
 // c) Store the results and statistical errors of reference flow estimates in histogram fIntFlow.
 //    Binning of fIntFlow is organized as follows:
 //
 //     1st bin: v{2,QC}
 //     2nd bin: v{4,QC}
 //     3rd bin: v{6,QC}
 //     4th bin: v{8,QC}
 
 if(!(fIntFlowCorrelationsHist && fIntFlowCovariances && fIntFlowQcumulants && fIntFlow))
 {
  cout<<"WARNING: fIntFlowCorrelationsHist && fIntFlowCovariances && fIntFlowQcumulants && fIntFlow is NULL in AFAWQC::CCIF() !!!!"<<endl;
  exit(0);
 }
   
 // Q-cumulants:
 Double_t qc2 = fIntFlowQcumulants->GetBinContent(1); // QC{2}  
 Double_t qc4 = fIntFlowQcumulants->GetBinContent(2); // QC{4}  
 Double_t qc6 = fIntFlowQcumulants->GetBinContent(3); // QC{6}  
 Double_t qc8 = fIntFlowQcumulants->GetBinContent(4); // QC{8}
  
 // correlations:
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>> 
 Double_t four = fIntFlowCorrelationsHist->GetBinContent(2); // <<4>>  
 Double_t six = fIntFlowCorrelationsHist->GetBinContent(3); // <<6>> 
 Double_t eight = fIntFlowCorrelationsHist->GetBinContent(4); // <<8>>  
 
 // statistical errors of average 2-, 4-, 6- and 8-particle azimuthal correlations:
 Double_t twoError = fIntFlowCorrelationsHist->GetBinError(1); // statistical error of <2>  
 Double_t fourError = fIntFlowCorrelationsHist->GetBinError(2); // statistical error of <4>   
 Double_t sixError = fIntFlowCorrelationsHist->GetBinError(3); // statistical error of <6> 
 Double_t eightError = fIntFlowCorrelationsHist->GetBinError(4); // statistical error of <8> 
 
 // covariances (multiplied by prefactor depending on weights - see comments in CalculateCovariancesIntFlow()):
 Double_t wCov24 = fIntFlowCovariances->GetBinContent(1); // Cov(<2>,<4>) * prefactor(w_<2>,w_<4>)
 Double_t wCov26 = fIntFlowCovariances->GetBinContent(2); // Cov(<2>,<6>) * prefactor(w_<2>,w_<6>)
 Double_t wCov28 = fIntFlowCovariances->GetBinContent(3); // Cov(<2>,<8>) * prefactor(w_<2>,w_<8>)
 Double_t wCov46 = fIntFlowCovariances->GetBinContent(4); // Cov(<4>,<6>) * prefactor(w_<4>,w_<6>)
 Double_t wCov48 = fIntFlowCovariances->GetBinContent(5); // Cov(<4>,<8>) * prefactor(w_<4>,w_<8>)
 Double_t wCov68 = fIntFlowCovariances->GetBinContent(6); // Cov(<6>,<8>) * prefactor(w_<6>,w_<8>)
  
 // integrated flow estimates:
 Double_t v2 = 0.; // v{2,QC}  
 Double_t v4 = 0.; // v{4,QC}  
 Double_t v6 = 0.; // v{6,QC}  
 Double_t v8 = 0.; // v{8,QC}
 
 // calculate integrated flow estimates from Q-cumulants: 
 if(qc2>=0.) v2 = pow(qc2,1./2.); 
 if(qc4<=0.) v4 = pow(-1.*qc4,1./4.); 
 if(qc6>=0.) v6 = pow((1./4.)*qc6,1./6.); 
 if(qc8<=0.) v8 = pow((-1./33.)*qc8,1./8.); 
   
 // statistical errors of integrated flow estimates:
 Double_t v2Error = 0.; // statistical error of v{2,QC}  
 Double_t v4Error = 0.; // statistical error of v{4,QC}  
 Double_t v6Error = 0.; // statistical error of v{6,QC}  
 Double_t v8Error = 0.; // statistical error of v{8,QC}
   
 // squares of statistical errors of integrated flow estimates:
 Double_t v2ErrorSquared = 0.; // squared statistical error of v{2,QC} 
 Double_t v4ErrorSquared = 0.; // squared statistical error of v{4,QC}   
 Double_t v6ErrorSquared = 0.; // squared statistical error of v{6,QC}   
 Double_t v8ErrorSquared = 0.; // squared statistical error of v{8,QC} 
 
 // calculate squared statistical errors of integrated flow estimates:
 if(two > 0.) 
 { 
  v2ErrorSquared = (1./(4.*two))*pow(twoError,2.);
 } 
 if(2.*pow(two,2.)-four > 0.)
 {
  v4ErrorSquared = (1./pow(2.*pow(two,2.)-four,3./2.))*
                   (pow(two,2.)*pow(twoError,2.)+(1./16.)*pow(fourError,2.)-(1./2.)*two*wCov24);
 }
 if(six-9.*four*two+12.*pow(two,3.) > 0.) 
 {
  v6ErrorSquared = ((1./2.)*(1./pow(2.,2./3.))*(1./pow(six-9.*four*two+12.*pow(two,3.),5./3.)))*
                   ((9./2.)*pow(4.*pow(two,2.)-four,2.)*pow(twoError,2.) 
                    + (9./2.)*pow(two,2.)*pow(fourError,2.)+(1./18.)*pow(sixError,2.)
                    - 9.*two*(4.*pow(two,2.)-four)*wCov24+(4.*pow(two,2.)-four)*wCov26-two*wCov46); 
 }
 if(-1.*eight+16.*six*two+18.*pow(four,2.)-144.*four*pow(two,2.)+144.*pow(two,4.) > 0.) 
 {
  v8ErrorSquared = (4./pow(33,1./4.))*(1./pow(-1.*eight+16.*six*two+18.*pow(four,2.)-144.*four*pow(two,2.)+144.*pow(two,4.),7./4.))*
                   (pow(36.*pow(two,3.)-18.*four*two+six,2.)*pow(twoError,2.)
                    + (81./16.)*pow(4.*pow(two,2.)-four,2.)*pow(fourError,2.)
                    + pow(two,2.)*pow(sixError,2.)
                    + (1./256.)*pow(eightError,2.)
                    - (9./2.)*(36.*pow(two,3.)-18.*four*two+six)*(4.*pow(two,2.)-four)*wCov24
                    + 2.*two*(36.*pow(two,3.)-18.*four*two+six)*wCov26
                    - (1./8.)*(36.*pow(two,3.)-18.*four*two+six)*wCov28                    
                    - (9./2.)*two*(4.*pow(two,2.)-four)*wCov46                   
                    + (9./32.)*(4.*pow(two,2.)-four)*wCov48                    
                    - (1./8.)*two*wCov68);
 } 

 // calculate statistical errors of integrated flow estimates: 
 if(v2ErrorSquared > 0.)
 {
  v2Error = pow(v2ErrorSquared,0.5);
 } else
   {
    cout<<"WARNING: Statistical error of v{2,QC} is imaginary !!!!"<<endl;
   }    
 if(v4ErrorSquared > 0.)
 {
  v4Error = pow(v4ErrorSquared,0.5);
 } else
   {
    cout<<"WARNING: Statistical error of v{4,QC} is imaginary !!!!"<<endl;
   }     
 if(v6ErrorSquared > 0.)
 {
  v6Error = pow(v6ErrorSquared,0.5);
 } else
   {
    cout<<"WARNING: Statistical error of v{6,QC} is imaginary !!!!"<<endl;
   }     
 if(v8ErrorSquared > 0.)
 {
  v8Error = pow(v8ErrorSquared,0.5);
 } else
   {
    cout<<"WARNING: Statistical error of v{8,QC} is imaginary !!!!"<<endl;
   }    
                     
 // store the results and statistical errors of integrated flow estimates:
 fIntFlow->SetBinContent(1,v2);
 fIntFlow->SetBinError(1,v2Error);
 fIntFlow->SetBinContent(2,v4);
 fIntFlow->SetBinError(2,v4Error);
 fIntFlow->SetBinContent(3,v6);
 fIntFlow->SetBinError(3,v6Error);
 fIntFlow->SetBinContent(4,v8);
 fIntFlow->SetBinError(4,v8Error);
       
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlow()


//================================================================================================================================ 


void AliFlowAnalysisWithQCumulants::FillCommonHistResultsIntFlow()
{
 // Fill in AliFlowCommonHistResults histograms relevant for 'NONAME' integrated flow (to be improved (name))
 
 if(!fIntFlow)
 {
  cout<<"WARNING: fIntFlow is NULL in AFAWQC::FCHRIF() !!!!"<<endl;
  exit(0); 
 }  
    
 if(!(fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th))
 {
  cout<<"WARNING: fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th"<<endl; 
  cout<<"         is NULL in AFAWQC::FCHRIF() !!!!"<<endl;
  exit(0);
 }
 
 Double_t v2 = fIntFlow->GetBinContent(1);
 Double_t v4 = fIntFlow->GetBinContent(2);
 Double_t v6 = fIntFlow->GetBinContent(3);
 Double_t v8 = fIntFlow->GetBinContent(4);
  
 Double_t v2Error = fIntFlow->GetBinError(1);
 Double_t v4Error = fIntFlow->GetBinError(2);
 Double_t v6Error = fIntFlow->GetBinError(3);
 Double_t v8Error = fIntFlow->GetBinError(4);
 
 fCommonHistsResults2nd->FillIntegratedFlow(v2,v2Error); // to be improved (hardwired 2nd in the name)  
 fCommonHistsResults4th->FillIntegratedFlow(v4,v4Error); // to be improved (hardwired 4th in the name)
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)) // to be improved (calculate also 6th and 8th order)
 {
  fCommonHistsResults6th->FillIntegratedFlow(v6,v6Error); // to be improved (hardwired 6th in the name)
  fCommonHistsResults8th->FillIntegratedFlow(v8,v8Error); // to be improved (hardwired 8th in the name) 
 }
 
} // end of AliFlowAnalysisWithQCumulants::FillCommonHistResultsIntFlow()


//================================================================================================================================ 


/*
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
*/


//================================================================================================================================


/*  
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
*/


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

 // 2-particle correlations:
 Double_t two1n1nW1W1 = 0.; // <w1 w2 cos(n*(phi1-phi2))>
 Double_t two2n2nW2W2 = 0.; // <w1^2 w2^2 cos(2n*(phi1-phi2))>
 Double_t two3n3nW3W3 = 0.; // <w1^3 w2^3 cos(3n*(phi1-phi2))>
 Double_t two4n4nW4W4 = 0.; // <w1^4 w2^4 cos(4n*(phi1-phi2))>
 if(dMult>1) 
 { 
  if(dM11)
  {
   two1n1nW1W1 = (pow(dReQ1n1k,2)+pow(dImQ1n1k,2)-(*fSMpk)(0,2))/dM11;    
   // average correlation <w1 w2 cos(n*(phi1-phi2))> for single event: 
   fIntFlowCorrelationsEBE->SetBinContent(1,two1n1nW1W1);
   fIntFlowEventWeightsForCorrelationsEBE->SetBinContent(1,dM11);
   // average correlation <w1 w2 cos(n*(phi1-phi2))> for all events:
   fIntFlowCorrelationsPro->Fill(0.5,two1n1nW1W1,dM11);   
   fIntFlowCorrelationsAllPro->Fill(0.5,two1n1nW1W1,dM11);   
  }
  if(dM22)
  {
   two2n2nW2W2 = (pow(dReQ2n2k,2)+pow(dImQ2n2k,2)-(*fSMpk)(0,4))/dM22; 
   // ...
   // average correlation <w1^2 w2^2 cos(2n*(phi1-phi2))> for all events:
   fIntFlowCorrelationsAllPro->Fill(1.5,two2n2nW2W2,dM22);   
  }
  if(dM33)
  {
   two3n3nW3W3 = (pow(dReQ3n3k,2)+pow(dImQ3n3k,2)-(*fSMpk)(0,6))/dM33;
   // ...
   // average correlation <w1^3 w2^3 cos(3n*(phi1-phi2))> for all events:
   fIntFlowCorrelationsAllPro->Fill(2.5,two3n3nW3W3,dM33);   
  }
  if(dM44)
  {
   two4n4nW4W4 = (pow(dReQ4n4k,2)+pow(dImQ4n4k,2)-(*fSMpk)(0,8))/dM44; 
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
   two1n1nW3W1 = (dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k-(*fSMpk)(0,4))/dM31; 
   fIntFlowExtraCorrelationsPro->Fill(0.5,two1n1nW3W1,dM31);  
  } 
  if(dM211)
  {
   two1n1nW1W1W2 = ((*fSMpk)(0,2)*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2)-(*fSMpk)(0,2))
                 - 2.*(dReQ1n3k*dReQ1n1k+dImQ1n3k*dImQ1n1k
                 - (*fSMpk)(0,4)))/dM211;
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
                     + 2.*(*fSMpk)(0,4))/dM211;                                                                               
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
                        - 4.*(*fSMpk)(0,2)*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2))
                        - 6.*(*fSMpk)(0,4)+2.*(*fSMpk)(1,2))/dM1111;  
                          
   // average correlation <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> for single event: 
   fIntFlowCorrelationsEBE->SetBinContent(2,four1n1n1n1nW1W1W1W1);
   fIntFlowEventWeightsForCorrelationsEBE->SetBinContent(2,dM1111);
   // average correlation <w1 w2 w3 w4 cos(n*(phi1+phi2-phi3-phi4))> for all events:
   fIntFlowCorrelationsPro->Fill(1.5,four1n1n1n1nW1W1W1W1,dM1111);   
   fIntFlowCorrelationsAllPro->Fill(10.5,four1n1n1n1nW1W1W1W1,dM1111);   
  } 
 } // end of if(dMult>3) 
 //..............................................................................................
 
} // end of AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrelationsUsingParticleWeights()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateWeightedQProductsForIntFlow() // to be improved (completed)
{
 // calculate averages like <<2><4>>, <<2><6>>, <<4><6>>, etc. which are needed to calculate covariances 
 // Remark: here we take weighted correlations!
 
 /*
 
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
 
 */
 
} // end of AliFlowAnalysisWithQCumulants::CalculateWeightedQProductsForIntFlow()  


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
  for(Int_t power=0;power<2;power++) // linear or quadratic 
  {
   fIntFlowSumOfEventWeightsNUA[sc][power] = NULL;
  }
 }
 for(Int_t power=0;power<2;power++) // linear or quadratic 
 {
  fIntFlowSumOfEventWeights[power] = NULL;    
 }
 for(Int_t i=0;i<3;i++) // print on the screen the final results (0=NONAME, 1=RP, 2=POI)
 {
  fPrintFinalResults[i] = kTRUE;
 }
 for(Int_t ci=0;ci<4;ci++) // correlation index or cumulant order
 {
  fIntFlowCorrelationsVsMPro[ci] = NULL;
  fIntFlowCorrelationsVsMHist[ci] = NULL;
  fIntFlowQcumulantsVsM[ci] = NULL;
  fIntFlowVsM[ci] = NULL;
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
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
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

 
 
 
 /*
 
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
 
 
 */
 
  
  
  /*
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
  */
    
} // end of AliFlowAnalysisWithQCumulants::InitializeArraysForDiffFlow()


//================================================================================================================================
 /*


void AliFlowAnalysisWithQCumulants::CalculateCorrelationsForDifferentialFlow2D(TString type)
{
 // calculate all reduced correlations needed for differential flow for each (pt,eta) bin: 
 
 if(type == "RP") // to be improved (removed)
 {
  cout<<endl;
 }
 // ... 
 
 
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
   q1n0kRe = fReEBE2D[2][0][0]->GetBinContent(fReEBE2D[2][0][0]->GetBin(p,e))
           * fReEBE2D[2][0][0]->GetBinEntries(fReEBE2D[2][0][0]->GetBin(p,e));
   q1n0kIm = fImEBE2D[2][0][0]->GetBinContent(fImEBE2D[2][0][0]->GetBin(p,e))
           * fImEBE2D[2][0][0]->GetBinEntries(fImEBE2D[2][0][0]->GetBin(p,e));
   q2n0kRe = fReEBE2D[2][1][0]->GetBinContent(fReEBE2D[2][1][0]->GetBin(p,e))
           * fReEBE2D[2][1][0]->GetBinEntries(fReEBE2D[2][1][0]->GetBin(p,e));
   q2n0kIm = fImEBE2D[2][1][0]->GetBinContent(fImEBE2D[2][1][0]->GetBin(p,e))
           * fImEBE2D[2][1][0]->GetBinEntries(fImEBE2D[2][1][0]->GetBin(p,e));
           
   mq = fReEBE2D[2][0][0]->GetBinEntries(fReEBE2D[2][0][0]->GetBin(p,e)); // to be improved (cross-checked by accessing other profiles here)
   
   if(type == "POI")
   {
    // p_{m*n,0}:
    p1n0kRe = fReEBE2D[1][0][0]->GetBinContent(fReEBE2D[1][0][0]->GetBin(p,e))
            * fReEBE2D[1][0][0]->GetBinEntries(fReEBE2D[1][0][0]->GetBin(p,e));
    p1n0kIm = fImEBE2D[1][0][0]->GetBinContent(fImEBE2D[1][0][0]->GetBin(p,e))  
            * fImEBE2D[1][0][0]->GetBinEntries(fImEBE2D[1][0][0]->GetBin(p,e));
            
    mp = fReEBE2D[1][0][0]->GetBinEntries(fReEBE2D[1][0][0]->GetBin(p,e)); // to be improved (cross-checked by accessing other profiles here)
    
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

 
   
    
      
} // end of AliFlowAnalysisWithQCumulants::CalculateCorrelationsForDifferentialFlow2D()



 
 

//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateWeightedCorrelationsForDifferentialFlow2D(TString type)
{
 // calculate all weighted correlations needed for differential flow 
 
  if(type == "RP") // to be improved (removed)
 {
  cout<<endl;
 }
 // ... 
 
 
 
 
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
    p1n0kRe = fReEBE2D[1][0][0]->GetBinContent(fReEBE2D[1][0][0]->GetBin(p,e))
            * fReEBE2D[1][0][0]->GetBinEntries(fReEBE2D[1][0][0]->GetBin(p,e));
    p1n0kIm = fImEBE2D[1][0][0]->GetBinContent(fImEBE2D[1][0][0]->GetBin(p,e))
            * fImEBE2D[1][0][0]->GetBinEntries(fImEBE2D[1][0][0]->GetBin(p,e)); 
            
    mp = fReEBE2D[1][0][0]->GetBinEntries(fReEBE2D[1][0][0]->GetBin(p,e));
    
    // q_{m*n,k}: 
    q1n2kRe = fReEBE2D[2][0][2]->GetBinContent(fReEBE2D[2][0][2]->GetBin(p,e))
            * fReEBE2D[2][0][2]->GetBinEntries(fReEBE2D[2][0][2]->GetBin(p,e));
    q1n2kIm = fImEBE2D[2][0][2]->GetBinContent(fImEBE2D[2][0][2]->GetBin(p,e))
            * fImEBE2D[2][0][2]->GetBinEntries(fImEBE2D[2][0][2]->GetBin(p,e));
    q2n1kRe = fReEBE2D[2][1][1]->GetBinContent(fReEBE2D[2][1][1]->GetBin(p,e))
            * fReEBE2D[2][1][1]->GetBinEntries(fReEBE2D[2][1][1]->GetBin(p,e)); 
    q2n1kIm = fImEBE2D[2][1][1]->GetBinContent(fImEBE2D[2][1][1]->GetBin(p,e))
            * fImEBE2D[2][1][1]->GetBinEntries(fImEBE2D[2][1][1]->GetBin(p,e));
       
    // s_{1,1}, s_{1,2} and s_{1,3} // to be improved (add explanation)  
    s1p1k = pow(fs2D[2][1]->GetBinContent(fs2D[2][1]->GetBin(p,e)),1.); 
    s1p2k = pow(fs2D[2][2]->GetBinContent(fs2D[2][2]->GetBin(p,e)),1.); 
    s1p3k = pow(fs2D[2][3]->GetBinContent(fs2D[2][3]->GetBin(p,e)),1.); 
   
    // M0111 from Eq. (118) in QC2c (to be improved (notation)):
    dM0111 = mp*(dSM3p1k-3.*dSM1p1k*dSM1p2k+2.*dSM1p3k)
           - 3.*(s1p1k*(dSM2p1k-dSM1p2k)
           + 2.*(s1p3k-s1p2k*dSM1p1k));
   }
   else if(type == "RP")
   {
    p1n0kRe = fReEBE2D[0][0][0]->GetBinContent(fReEBE2D[0][0][0]->GetBin(p,e))
            * fReEBE2D[0][0][0]->GetBinEntries(fReEBE2D[0][0][0]->GetBin(p,e));
    p1n0kIm = fImEBE2D[0][0][0]->GetBinContent(fImEBE2D[0][0][0]->GetBin(p,e))
            * fImEBE2D[0][0][0]->GetBinEntries(fImEBE2D[0][0][0]->GetBin(p,e));
            
    mp = fReEBE2D[0][0][0]->GetBinEntries(fReEBE2D[0][0][0]->GetBin(p,e));
    
    // q_{m*n,k}: 
    q1n2kRe = fReEBE2D[0][0][2]->GetBinContent(fReEBE2D[0][0][2]->GetBin(p,e))
            * fReEBE2D[0][0][2]->GetBinEntries(fReEBE2D[0][0][2]->GetBin(p,e));
    q1n2kIm = fImEBE2D[0][0][2]->GetBinContent(fImEBE2D[0][0][2]->GetBin(p,e))
            * fImEBE2D[0][0][2]->GetBinEntries(fImEBE2D[0][0][2]->GetBin(p,e));
    q2n1kRe = fReEBE2D[0][1][1]->GetBinContent(fReEBE2D[0][1][1]->GetBin(p,e))
            * fReEBE2D[0][1][1]->GetBinEntries(fReEBE2D[0][1][1]->GetBin(p,e));
    q2n1kIm = fImEBE2D[0][1][1]->GetBinContent(fImEBE2D[0][1][1]->GetBin(p,e))
            * fImEBE2D[0][1][1]->GetBinEntries(fImEBE2D[0][1][1]->GetBin(p,e));
   
    // s_{1,1}, s_{1,2} and s_{1,3} // to be improved (add explanation)  
    s1p1k = pow(fs2D[0][1]->GetBinContent(fs2D[0][1]->GetBin(p,e)),1.); 
    s1p2k = pow(fs2D[0][2]->GetBinContent(fs2D[0][2]->GetBin(p,e)),1.); 
    s1p3k = pow(fs2D[0][3]->GetBinContent(fs2D[0][3]->GetBin(p,e)),1.); 
   
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
 
  
    
      
} // end of AliFlowAnalysisWithQCumulants::CalculateWeightedCorrelationsForDifferentialFlow2D(TString type)


//================================================================================================================================

 */  

/*
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
*/


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCumulants(TString type, TString ptOrEta)
{
 // calcualate cumulants for differential flow from measured correlations
 // Remark: cumulants calculated here are NOT corrected for non-uniform acceptance. This correction is applied in the method ...
 // to be improved (description) 
 
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;

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
 
 // correlation <<2>>: 
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1);
 
 // 1D:
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // reduced correlations:   
  Double_t twoPrime = fDiffFlowCorrelationsHist[t][pe][0]->GetBinContent(b); // <<2'>>(pt)
  Double_t fourPrime = fDiffFlowCorrelationsHist[t][pe][1]->GetBinContent(b); // <<4'>>(pt)
  // final statistical error of reduced correlations:
  //Double_t twoPrimeError = fFinalCorrelations1D[t][pW][eW][0][0]->GetBinError(p); 
  // QC{2'}:
  Double_t qc2Prime = twoPrime; // QC{2'}
  //Double_t qc2PrimeError = twoPrimeError; // final stat. error of QC{2'}
  fDiffFlowCumulants[t][pe][0]->SetBinContent(b,qc2Prime); 
  //fFinalCumulantsPt[t][pW][eW][nua][0]->SetBinError(p,qc2PrimeError);   
  // QC{4'}:
  Double_t qc4Prime = fourPrime - 2.*twoPrime*two; // QC{4'} = <<4'>> - 2*<<2'>><<2>>
  fDiffFlowCumulants[t][pe][1]->SetBinContent(b,qc4Prime); 
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
    
 /* 
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
 */
   
} // end of void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCumulants(TString type, Bool_t useParticleWeights, TString eventWeights); 


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateFinalResultsForRPandPOIIntegratedFlow(TString type)
{
 // calculate final results for integrated flow of RPs and POIs 
  
 Int_t typeFlag = -1;

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
     
 // shortcuts:
 Int_t t = typeFlag;
  
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
 
 // to be improved (hardwired pt index)
 flow2ndPt = (TH1D*)fDiffFlow[t][0][0]->Clone();
 flow4thPt = (TH1D*)fDiffFlow[t][0][1]->Clone();
 flow6thPt = (TH1D*)fDiffFlow[t][0][2]->Clone();
 flow8thPt = (TH1D*)fDiffFlow[t][0][3]->Clone(); 
   
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
 //  e) Book and nest list for nested loops;
 
 // a) Book and nest all lists for integrated flow:
 // base list for integrated flow:
 fIntFlowList = new TList();
 fIntFlowList->SetName("Integrated Flow");
 fIntFlowList->SetOwner(kTRUE);
 fHistList->Add(fIntFlowList);
 // list holding profiles: 
 fIntFlowProfiles = new TList();
 fIntFlowProfiles->SetName("Profiles");
 fIntFlowProfiles->SetOwner(kTRUE);
 fIntFlowList->Add(fIntFlowProfiles);
 // list holding histograms with results:
 fIntFlowResults = new TList();
 fIntFlowResults->SetName("Results");
 fIntFlowResults->SetOwner(kTRUE);
 fIntFlowList->Add(fIntFlowResults);
 
 // b) Book and nest lists for differential flow;
 fDiffFlowList = new TList();
 fDiffFlowList->SetName("Differential Flow");
 fDiffFlowList->SetOwner(kTRUE); 
 fHistList->Add(fDiffFlowList);
 // list holding profiles: 
 fDiffFlowProfiles = new TList(); 
 fDiffFlowProfiles->SetName("Profiles");
 fDiffFlowProfiles->SetOwner(kTRUE);
 fDiffFlowList->Add(fDiffFlowProfiles);
 // list holding histograms with results: 
 fDiffFlowResults = new TList();
 fDiffFlowResults->SetName("Results");
 fDiffFlowResults->SetOwner(kTRUE);
 fDiffFlowList->Add(fDiffFlowResults);
 // flags used for naming nested lists in list fDiffFlowProfiles and fDiffFlowResults:  
 TList list;
 list.SetOwner(kTRUE);
 TString typeFlag[2] = {"RP","POI"};  
 TString ptEtaFlag[2] = {"p_{T}","#eta"}; 
 TString powerFlag[2] = {"linear","quadratic"};   
 // nested lists in fDiffFlowProfiles (~/Differential Flow/Profiles):
 for(Int_t t=0;t<2;t++) // type: RP or POI
 {
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
   // list holding histograms with differential flow estimates from Q-cumulants:
   fDiffFlowHistList[t][pe] = (TList*)list.Clone();
   fDiffFlowHistList[t][pe]->SetName(Form("Differential flow (%s, %s)",typeFlag[t].Data(),ptEtaFlag[pe].Data()));
   fDiffFlowResults->Add(fDiffFlowHistList[t][pe]);      
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
  
 // c) Book and nest list for particle weights:
 fWeightsList->SetName("Weights");
 fWeightsList->SetOwner(kTRUE);   
 fHistList->Add(fWeightsList); 

 // d) Book and nest list for distributions:
 fDistributionsList = new TList();
 fDistributionsList->SetName("Distributions");
 fDistributionsList->SetOwner(kTRUE);
 fHistList->Add(fDistributionsList);
 
 // e) Book and nest list for nested loops:
 fNestedLoopsList = new TList();
 fNestedLoopsList->SetName("Nested Loops");
 fNestedLoopsList->SetOwner(kTRUE);
 fHistList->Add(fNestedLoopsList);
 
} // end of void AliFlowAnalysisWithQCumulants::BookAndNestAllLists()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::FillCommonHistResultsDiffFlow(TString type)
{
 // fill common result histograms for differential flow
 
 Int_t typeFlag = -1;
 //Int_t ptEtaFlag = -1;

 if(type == "RP")
 {
  typeFlag = 0;
 } else if(type == "POI")
   {
    typeFlag = 1;
   } 
  
 // shortcuts:
 Int_t t = typeFlag;
 //Int_t pe = ptEtaFlag;

 // to be improved (implement protection here)
     
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


void AliFlowAnalysisWithQCumulants::AccessConstants()
{
 // Access needed common constants from AliFlowCommonConstants
 
 fnBinsPhi = AliFlowCommonConstants::GetMaster()->GetNbinsPhi();
 fPhiMin = AliFlowCommonConstants::GetMaster()->GetPhiMin();	     
 fPhiMax = AliFlowCommonConstants::GetMaster()->GetPhiMax();
 if(fnBinsPhi) fPhiBinWidth = (fPhiMax-fPhiMin)/fnBinsPhi;  
 fnBinsPt = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
 fPtMin = AliFlowCommonConstants::GetMaster()->GetPtMin();	     
 fPtMax = AliFlowCommonConstants::GetMaster()->GetPtMax();
 if(fnBinsPt) fPtBinWidth = (fPtMax-fPtMin)/fnBinsPt;  
 fnBinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
 fEtaMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();	     
 fEtaMax = AliFlowCommonConstants::GetMaster()->GetEtaMax();
 if(fnBinsEta) fEtaBinWidth = (fEtaMax-fEtaMin)/fnBinsEta;  
 
} // end of void AliFlowAnalysisWithQCumulants::AccessConstants()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CrossCheckSettings()
{
 // a) Cross check if the choice for multiplicity weights make sense;
 
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
 Double_t dMult = (*fSMpk)(0,0);
                       
 for(Int_t p=0;p<2;p++) // power-1
 {
  for(Int_t ci=0;ci<4;ci++) // correlation index
  { 
   fIntFlowSumOfEventWeights[p]->Fill(ci+0.5,pow(fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci+1),p+1)); 
   fIntFlowSumOfEventWeightsVsM[ci][p]->Fill(dMult+0.5,pow(fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci+1),p+1));
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
   for(Int_t ci=0;ci<3;ci++) // nua term index
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
 Double_t dMult = (*fSMpk)(0,0);
  
 Int_t counter = 0;
 
 for(Int_t ci1=1;ci1<4;ci1++)
 {
  for(Int_t ci2=ci1+1;ci2<=4;ci2++)
  {
   fIntFlowSumOfProductOfEventWeights->Fill(0.5+counter,
                                            fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci1)*
                                            fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci2));                                           
   fIntFlowSumOfProductOfEventWeightsVsM[counter]->Fill(dMult+0.5,
                                                        fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci1)*
                                                        fIntFlowEventWeightsForCorrelationsEBE->GetBinContent(ci2));
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
 // calculate reduced correlations for RPs or POIs in pt or eta bins

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

 // reduced correlations are stored in fDiffFlowCorrelationsPro[0=RP,1=POI][0=pt,1=eta][correlation index]. Correlation index runs as follows:
 // 
 // 0: <<2'>>
 // 1: <<4'>>
 // 2: <<6'>>
 // 3: <<8'>>
 
 Int_t t = -1; // type flag 
 Int_t pe = -1; // ptEta flag
 
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
      
   // 2'-particle correlation for particular (pt,eta) bin:
   Double_t two1n1nPtEta = 0.;
   if(mp*dMult-mq)
   {
    two1n1nPtEta = (p1n0kRe*dReQ1n+p1n0kIm*dImQ1n-mq)
                 / (mp*dMult-mq);
   
    if(type == "POI") // to be improved (I do not this if)
    { 
     // fill profile to get <<2'>> for POIs
     fDiffFlowCorrelationsPro[1][pe][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],two1n1nPtEta,mp*dMult-mq);
     // histogram to store <2'> for POIs e-b-e (needed in some other methods):
     fDiffFlowCorrelationsEBE[1][pe][0]->SetBinContent(b,two1n1nPtEta);      
     fDiffFlowEventWeightsForCorrelationsEBE[1][pe][0]->SetBinContent(b,mp*dMult-mq);      
    }
    else if(type == "RP") // to be improved (I do not this if)
    {
     // profile to get <<2'>> for RPs:
     fDiffFlowCorrelationsPro[0][pe][0]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],two1n1nPtEta,mp*dMult-mq);
     // histogram to store <2'> for RPs e-b-e (needed in some other methods):
     fDiffFlowCorrelationsEBE[0][pe][0]->SetBinContent(b,two1n1nPtEta); 
     fDiffFlowEventWeightsForCorrelationsEBE[0][pe][0]->SetBinContent(b,mp*dMult-mq); 
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
    
    if(type == "POI")
    {
     // profile to get <<4'>> for POIs:
     fDiffFlowCorrelationsPro[1][pe][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],four1n1n1n1nPtEta,
                                     (mp-mq)*dMult*(dMult-1.)*(dMult-2.)
                                     + mq*(dMult-1.)*(dMult-2.)*(dMult-3.)); 
     // histogram to store <4'> for POIs e-b-e (needed in some other methods):
     fDiffFlowCorrelationsEBE[1][pe][1]->SetBinContent(b,four1n1n1n1nPtEta);                               
     fDiffFlowEventWeightsForCorrelationsEBE[1][pe][1]->SetBinContent(b,(mp-mq)*dMult*(dMult-1.)*(dMult-2.)
                                                                        + mq*(dMult-1.)*(dMult-2.)*(dMult-3.));                               
    }
    else if(type == "RP")
    {
     // profile to get <<4'>> for RPs:
     fDiffFlowCorrelationsPro[0][pe][1]->Fill(minPtEta[pe]+(b-1)*binWidthPtEta[pe],four1n1n1n1nPtEta,
                                       (mp-mq)*dMult*(dMult-1.)*(dMult-2.)
                                       + mq*(dMult-1.)*(dMult-2.)*(dMult-3.));                   
     // histogram to store <4'> for RPs e-b-e (needed in some other methods):
     fDiffFlowCorrelationsEBE[0][pe][1]->SetBinContent(b,four1n1n1n1nPtEta);                   
     fDiffFlowEventWeightsForCorrelationsEBE[0][pe][1]->SetBinContent(b,(mp-mq)*dMult*(dMult-1.)*(dMult-2.)
                                                                        + mq*(dMult-1.)*(dMult-2.)*(dMult-3.));                   
    }
   } // end of if((mp-mq)*dMult*(dMult-1.)*(dMult-2.)
     //            +mq*(dMult-1.)*(dMult-2.)*(dMult-3.))
   
 } // end of for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 
   
} // end of void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrelations(TString type, TString ptOrEta);


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateDiffFlowSumOfEventWeights(TString type, TString ptOrEta)
{
 // Calculate sums of various event weights for reduced correlations. 
 // (These quantitites are needed in expressions for unbiased estimators relevant for the statistical errors.)

 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;

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
 Double_t dMult = (*fSMpk)(0,0); // total event multiplicity
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
  
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;

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
 Double_t dMult = (*fSMpk)(0,0); // total event multiplicity
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

 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;

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
             
 for(Int_t rci=0;rci<4;rci++)
 {
  if(!fDiffFlowCorrelationsPro[t][pe][rci])
  {
   cout<<"WARNING: fDiffFlowCorrelationsPro[t][pe][rci] is NULL in AFAWQC::FRC() !!!!"<<endl;
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
   correlation = fDiffFlowCorrelationsPro[t][pe][rci]->GetBinContent(b); 
   spread = fDiffFlowCorrelationsPro[t][pe][rci]->GetBinError(b);
   sumOfWeights = fDiffFlowSumOfEventWeights[t][pe][0][rci]->GetBinContent(b);
   sumOfSquaredWeights = fDiffFlowSumOfEventWeights[t][pe][1][rci]->GetBinContent(b);
   if(sumOfWeights) termA = (pow(sumOfSquaredWeights,0.5)/sumOfWeights);
   if(1.-pow(termA,2.)>0.) termB = 1./pow(1.-pow(termA,2.),0.5); 
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
  
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;

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
 Double_t dMult = (*fSMpk)(0,0); // multiplicity (number of particles used to determine the reaction plane)
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
  
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;

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
   if(TMath::Abs(denominator)>1e-6)
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
   if(TMath::Abs(denominator)>1e-6)
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
   if(TMath::Abs(denominator)>1e-6)
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
   if(TMath::Abs(denominator)>1e-6)
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
   if(TMath::Abs(denominator)>1e-6)
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
 // calculate differential flow from differential cumulants and previously obtained integrated flow: (to be improved: description)
 
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;

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
   
 // correlations:
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>>
 Double_t four = fIntFlowCorrelationsHist->GetBinContent(2); // <<4>>
 
 // statistical errors of correlations:
 Double_t twoError = fIntFlowCorrelationsHist->GetBinError(1);
 Double_t fourError = fIntFlowCorrelationsHist->GetBinError(2);   
    
 // reduced correlations:
 Double_t twoReduced = 0.; // <<2'>>
 Double_t fourReduced = 0.; // <<4'>>
 
 // statistical errors of reduced correlations:
 Double_t twoReducedError = 0.; 
 Double_t fourReducedError = 0.; 

 // covariances:
 Double_t wCovTwoFour = fIntFlowCovariances->GetBinContent(1);// // Cov(<2>,<4>) * prefactor(<2>,<4>)
 Double_t wCovTwoTwoReduced = 0.; // Cov(<2>,<2'>) * prefactor(<2>,<2'>)
 Double_t wCovTwoFourReduced = 0.; // Cov(<2>,<4'>) * prefactor(<2>,<4'>)
 Double_t wCovFourTwoReduced = 0.; // Cov(<4>,<2'>) * prefactor(<4>,<2'>)
 Double_t wCovFourFourReduced = 0.; // Cov(<4>,<4'>) * prefactor(<4>,<4'>)
 Double_t wCovTwoReducedFourReduced = 0.; // Cov(<2'>,<4'>) * prefactor(<2'>,<4'>)
 
 // differential flow:
 Double_t v2Prime = 0.; // v'{2}                   
 Double_t v4Prime = 0.; // v'{4}
 
 // statistical error of differential flow:
 Double_t v2PrimeError = 0.;                    
 Double_t v4PrimeError = 0.; 
 
 // squared statistical error of differential flow:
 Double_t v2PrimeErrorSquared = 0.;                    
 Double_t v4PrimeErrorSquared = 0.; 
 
 // loop over pt or eta bins:
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // reduced correlations and statistical errors:
  twoReduced = fDiffFlowCorrelationsHist[t][pe][0]->GetBinContent(b);
  twoReducedError = fDiffFlowCorrelationsHist[t][pe][0]->GetBinError(b);
  fourReduced = fDiffFlowCorrelationsHist[t][pe][1]->GetBinContent(b);
  fourReducedError = fDiffFlowCorrelationsHist[t][pe][1]->GetBinError(b);
  // covariances:
  wCovTwoTwoReduced = fDiffFlowCovariances[t][pe][0]->GetBinContent(b);
  wCovTwoFourReduced = fDiffFlowCovariances[t][pe][1]->GetBinContent(b);
  wCovFourTwoReduced = fDiffFlowCovariances[t][pe][2]->GetBinContent(b);
  wCovFourFourReduced = fDiffFlowCovariances[t][pe][3]->GetBinContent(b);
  wCovTwoReducedFourReduced = fDiffFlowCovariances[t][pe][4]->GetBinContent(b);
  // differential flow:
  // v'{2}:
  if(two>0.) 
  {
   v2Prime = twoReduced/pow(two,0.5);
   v2PrimeErrorSquared = (1./4.)*pow(two,-3.)*
                         (pow(twoReduced,2.)*pow(twoError,2.)
                          + 4.*pow(two,2.)*pow(twoReducedError,2.)
                          - 4.*two*twoReduced*wCovTwoTwoReduced);
     
                                                            
   if(v2PrimeErrorSquared>0.) v2PrimeError = pow(v2PrimeErrorSquared,0.5);
   fDiffFlow[t][pe][0]->SetBinContent(b,v2Prime); 
   if(TMath::Abs(v2Prime)>1.e-44)fDiffFlow[t][pe][0]->SetBinError(b,v2PrimeError);     
  }
  // differential flow:
  // v'{4}
  if(2.*pow(two,2.)-four > 0.) 
  {
   v4Prime = (2.*two*twoReduced-fourReduced)/pow(2.*pow(two,2.)-four,3./4.);
   v4PrimeErrorSquared = pow(2.*pow(two,2.)-four,-7./2.)*
                         (pow(2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced,2.)*pow(twoError,2.)
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
   if(v4PrimeErrorSquared>0.) v4PrimeError = pow(v4PrimeErrorSquared,0.5);        
   fDiffFlow[t][pe][1]->SetBinContent(b,v4Prime);
   if(TMath::Abs(v4Prime)>1.e-44)fDiffFlow[t][pe][1]->SetBinError(b,v4PrimeError);     
  }
  
 } // end of for(Int_t b=1;b<=fnBinsPtEta[pe];b++)
 
   
 
 
 /*
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
 */

} // end of AliFlowAnalysisWithQCumulants::CalculateDiffFlow(TString type, Bool_t useParticleWeights)


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
 fIntFlowFlags->Fill(0.5,(Int_t)fUsePhiWeights||fUsePtWeights||fUseEtaWeights);
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
 // corrected for non-uniform acceptance or not:
 fIntFlowFlags->Fill(2.5,(Int_t)fApplyCorrectionForNUA);
 fIntFlowFlags->Fill(3.5,(Int_t)fPrintFinalResults[0]);
 fIntFlowFlags->Fill(4.5,(Int_t)fPrintFinalResults[1]);
 fIntFlowFlags->Fill(5.5,(Int_t)fPrintFinalResults[2]);
  
} // end of void AliFlowAnalysisWithQCumulants::StoreIntFlowFlags()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::StoreDiffFlowFlags()
{
 // Store all flags for differential flow in the profile fDiffFlowFlags.
  
 if(!fDiffFlowFlags)
 {
  cout<<"WARNING: fDiffFlowFlags is NULL in AFAWQC::SFFDF() !!!!"<<endl;
  exit(0);
 } 
 
 fDiffFlowFlags->Fill(0.5,fUsePhiWeights||fUsePtWeights||fUseEtaWeights); // particle weights used or not
 //fDiffFlowFlags->Fill(1.5,""); // which event weight was used? // to be improved
 fDiffFlowFlags->Fill(2.5,fApplyCorrectionForNUA); // corrected for non-uniform acceptance or not
 fDiffFlowFlags->Fill(3.5,fCalculate2DFlow); // calculate also 2D differential flow in (pt,eta) or not
    
} // end of void AliFlowAnalysisWithQCumulants::StoreDiffFlowFlags()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::GetPointersForCommonHistograms() 
{
 // Access all pointers to common control and common result histograms and profiles.
 
 TString commonHistsName = "AliFlowCommonHistQC";
 commonHistsName += fAnalysisLabel->Data();
 AliFlowCommonHist *commonHist = dynamic_cast<AliFlowCommonHist*>(fHistList->FindObject(commonHistsName.Data()));
 if(commonHist) this->SetCommonHists(commonHist); 
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
 AliFlowCommonHistResults *commonHistRes2nd = dynamic_cast<AliFlowCommonHistResults*>                                              (fHistList->FindObject(commonHistResults2ndOrderName.Data()));
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
 if(weightsList) this->SetWeightsList(weightsList);
 TString fUseParticleWeightsName = "fUseParticleWeightsQC"; // to be improved (hirdwired label QC)
 fUseParticleWeightsName += fAnalysisLabel->Data();
 TProfile *useParticleWeights = dynamic_cast<TProfile*>(weightsList->FindObject(fUseParticleWeightsName.Data()));
 if(useParticleWeights)
 {
  this->SetUseParticleWeights(useParticleWeights);  
  fUsePhiWeights = (Int_t)fUseParticleWeights->GetBinContent(1); 
  fUsePtWeights = (Int_t)fUseParticleWeights->GetBinContent(2); 
  fUseEtaWeights = (Int_t)fUseParticleWeights->GetBinContent(3);  
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
 TString correlationFlag[4] = {"<<2>>","<<4>>","<<6>>","<<8>>"}; // to be improved (should I promote this to data member?)
 
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
  Bool_t bApplyCorrectionForNUA = kFALSE;
  if(intFlowFlags)
  {
   this->SetIntFlowFlags(intFlowFlags);  
   bApplyCorrectionForNUA = (Int_t)intFlowFlags->GetBinContent(3); 
   this->SetApplyCorrectionForNUA(bApplyCorrectionForNUA);      
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
   // average correlations <<2>>, <<4>>, <<6>> and <<8>> versus multiplicity for all events (error is wrong here):
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
   if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
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
   } // end of if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)        
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
   TString intFlowProductOfCorrelationsVsMProName = "fIntFlowProductOfCorrelationsVsMPro";
   intFlowProductOfCorrelationsVsMProName += fAnalysisLabel->Data();
   TString productFlag[6] = {"<<2><4>>","<<2><6>>","<<2><8>>","<<4><6>>","<<4><8>>","<<6><8>>"};
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
   // sum of linear and quadratic event weights for <2>, <4>, <6> and <8> versus multiplicity
   // [0=sum{w_{<2>}},1=sum{w_{<4>}},2=sum{w_{<6>}},3=sum{w_{<8>}}][0=linear 1,1=quadratic]:
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
   // sum of products of event weights for correlations <2>, <4>, <6> and <8> vs M
   // [0=sum{w_{<2>}w_{<4>}},1=sum{w_{<2>}w_{<6>}},2=sum{w_{<2>}w_{<8>}},
   //  3=sum{w_{<4>}w_{<6>}},4=sum{w_{<4>}w_{<8>}},5=sum{w_{<6>}w_{<8>}}]:  
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
   // final results for integrated Q-cumulants:
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
   // final results for integrated Q-cumulants versus multiplicity:
   TString intFlowQcumulantsVsMName = "fIntFlowQcumulantsVsM";
   intFlowQcumulantsVsMName += fAnalysisLabel->Data();
   TString cumulantFlag[4] = {"QC{2}","QC{4}","QC{6}","QC{8}"};
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
   // final integrated flow estimates from Q-cumulants:
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
   // integrated flow from Q-cumulants versus multiplicity:
   TString intFlowVsMName = "fIntFlowVsM";
   intFlowVsMName += fAnalysisLabel->Data();
   TString flowFlag[4] = {"v_{2}{2,QC}","v_{2}{4,QC}","v_{2}{6,QC}","v_{2}{8,QC}"}; // to be improved (harwired harmonic)
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
  } else // to if(intFlowResults)
    {
     cout<<"WARNING: intFlowResults is NULL in AFAWQC::GPFIFH() !!!!"<<endl;
    }
    
} // end of void AliFlowAnalysisWithQCumulants::GetPointersForIntFlowHistograms()

//================================================================================================================================

void AliFlowAnalysisWithQCumulants::GetPointersForDiffFlowHistograms()
{
 // Get pointer to all objects relevant for differential flow.
 //  a) Define flags locally (to be improved: should I promote flags to data members?);
 //  b) Get pointer to base list for differential flow fDiffFlowList and nested lists fDiffFlowListProfiles and fDiffFlowListResults;
 //  c) Get pointer to profile fDiffFlowFlags holding all flags for differential flow;
 //  d) Get pointers to all nested lists in fDiffFlowListProfiles and to profiles which they hold;
 //  e) Get pointers to all nested lists in fDiffFlowListResults and to histograms which they hold.
 
 // a) Define flags locally (to be improved: should I promote flags to data members?): 
 TString typeFlag[2] = {"RP","POI"}; 
 TString ptEtaFlag[2] = {"p_{T}","#eta"};
 TString powerFlag[2] = {"linear","quadratic"};
 TString sinCosFlag[2] = {"sin","cos"};
 TString differentialCumulantIndex[4] = {"QC{2'}","QC{4'}","QC{6'}","QC{8'}"};  
 TString differentialFlowIndex[4] = {"v'{2}","v'{4}","v'{6}","v'{8}"};  
 TString reducedCorrelationIndex[4] = {"<2'>","<4'>","<6'>","<8'>"};
 TString mixedCorrelationIndex[8] = {"<2>","<2'>","<4>","<4'>","<6>","<6'>","<8>","<8'>"};
 TString covarianceName[5] = {"Cov(<2>,<2'>)","Cov(<2>,<4'>)","Cov(<4>,<2'>)","Cov(<4>,<4'>)","Cov(<2'>,<4'>)"}; 
  
 // b) Get pointer to base list for differential flow fDiffFlowList and nested lists fDiffFlowListProfiles and fDiffFlowListResults:
 TList *diffFlowList = NULL;
 diffFlowList = dynamic_cast<TList*>(fHistList->FindObject("Differential Flow"));  
 if(!diffFlowList)
 { 
  cout<<"WARNING: diffFlowList is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
  exit(0);
 }
 // list holding nested lists containing profiles:
 TList *diffFlowListProfiles = NULL;
 diffFlowListProfiles = dynamic_cast<TList*>(diffFlowList->FindObject("Profiles"));
 if(!diffFlowListProfiles)
 { 
  cout<<"WARNING: diffFlowListProfiles is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
  exit(0);
 }
 // list holding nested lists containing 2D and 1D histograms with final results:
 TList *diffFlowListResults = NULL;
 diffFlowListResults = dynamic_cast<TList*>(diffFlowList->FindObject("Results"));
 if(!diffFlowListResults)
 { 
  cout<<"WARNING: diffFlowListResults is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
  exit(0);
 }
 
 // c) Get pointer to profile holding all flags for differential flow;
 TString diffFlowFlagsName = "fDiffFlowFlags";
 diffFlowFlagsName += fAnalysisLabel->Data();
 TProfile *diffFlowFlags = dynamic_cast<TProfile*>(diffFlowList->FindObject(diffFlowFlagsName.Data()));
 Bool_t bCalculate2DFlow = kFALSE;
 if(diffFlowFlags)
 {
  this->SetDiffFlowFlags(diffFlowFlags);  
  bCalculate2DFlow = (Int_t)diffFlowFlags->GetBinContent(4);
  this->SetCalculate2DFlow(bCalculate2DFlow); // to be improved (shoul I call this setter somewhere else?)     
 }
  
 // d) Get pointers to all nested lists in fDiffFlowListProfiles and to profiles which they hold;
 // correlations:
 TList *diffFlowCorrelationsProList[2][2] = {{NULL}};
 TString diffFlowCorrelationsProName = "fDiffFlowCorrelationsPro";
 diffFlowCorrelationsProName += fAnalysisLabel->Data();
 TProfile *diffFlowCorrelationsPro[2][2][4] = {{{NULL}}};   
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
  for(Int_t pe=0;pe<2;pe++)
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
        cout<<"WARNING: diffFlowCorrelationsPro[t][pe][ci] is NULL in AFAWQC::GPFDFH() !!!!"<<endl;
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
  
 // e) Get pointers to all nested lists in fDiffFlowListResults and to histograms which they hold.
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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


void AliFlowAnalysisWithQCumulants::BookEverythingForDifferentialFlow()
{
 // Book all histograms and profiles needed for differential flow.
 //  a) Define flags locally (to be improved: should I promote flags to data members?);
 //  b) Book profile to hold all flags for differential flow;
 //  c) Book e-b-e quantities;
 //  d) Book profiles;
 //  e) Book histograms holding final results. 
 
 // a) Define flags locally (to be improved: should I promote flags to data members?): 
 TString typeFlag[2] = {"RP","POI"}; 
 TString ptEtaFlag[2] = {"p_{T}","#eta"};
 TString powerFlag[2] = {"linear","quadratic"};
 TString sinCosFlag[2] = {"sin","cos"};
 TString differentialCumulantIndex[4] = {"QC{2'}","QC{4'}","QC{6'}","QC{8'}"};  
 TString differentialFlowIndex[4] = {"v'{2}","v'{4}","v'{6}","v'{8}"};  
 TString reducedCorrelationIndex[4] = {"<2'>","<4'>","<6'>","<8'>"};
 TString mixedCorrelationIndex[8] = {"<2>","<2'>","<4>","<4'>","<6>","<6'>","<8>","<8'>"};
 TString covarianceName[5] = {"Cov(<2>,<2'>)","Cov(<2>,<4'>)","Cov(<4>,<2'>)","Cov(<4>,<4'>)","Cov(<2'>,<4'>)"}; 
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 Double_t maxPtEta[2] = {fPtMax,fEtaMax};
  
 // b) Book profile to hold all flags for differential flow:
 TString diffFlowFlagsName = "fDiffFlowFlags";
 diffFlowFlagsName += fAnalysisLabel->Data();
 fDiffFlowFlags = new TProfile(diffFlowFlagsName.Data(),"Flags for Differential Flow",4,0,4);
 fDiffFlowFlags->SetTickLength(-0.01,"Y");
 fDiffFlowFlags->SetMarkerStyle(25);
 fDiffFlowFlags->SetLabelSize(0.05);
 fDiffFlowFlags->SetLabelOffset(0.02,"Y");
 (fDiffFlowFlags->GetXaxis())->SetBinLabel(1,"Particle Weights");
 (fDiffFlowFlags->GetXaxis())->SetBinLabel(2,"Event Weights");
 (fDiffFlowFlags->GetXaxis())->SetBinLabel(3,"Corrected for NUA?");
 (fDiffFlowFlags->GetXaxis())->SetBinLabel(4,"Calculated 2D flow?");
 fDiffFlowList->Add(fDiffFlowFlags);

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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
 // 2D:
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
 // reduced correlations e-b-e:
 TString diffFlowCorrelationsEBEName = "fDiffFlowCorrelationsEBE";
 diffFlowCorrelationsEBEName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
 // corrections terms:
 TString diffFlowCorrectionTermsForNUAProName = "fDiffFlowCorrectionTermsForNUAPro";
 diffFlowCorrectionTermsForNUAProName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<2;pe++) // pt or eta
  {
   for(Int_t rci=0;rci<4;rci++) // reduced correlation index
   {
    // reduced correlations:
    fDiffFlowCorrelationsPro[t][pe][rci] = new TProfile(Form("%s, %s, %s, %s",diffFlowCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[rci].Data()),Form("%s, %s, %s, %s",diffFlowCorrelationsProName.Data(),typeFlag[t].Data(),ptEtaFlag[pe].Data(),reducedCorrelationIndex[rci].Data()),nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe],"s");
    fDiffFlowCorrelationsPro[t][pe][rci]->SetXTitle(ptEtaFlag[pe].Data());
    fDiffFlowCorrelationsProList[t][pe]->Add(fDiffFlowCorrelationsPro[t][pe][rci]); // to be improved (add dedicated list to hold reduced correlations)
   } // end of for(Int_t rci=0;rci<4;rci++) // correlation index
  } // end of for(Int_t pe=0;pe<2;pe++) // pt or eta 
 } // end of for(Int_t t=0;t<2;t++) // type: RP or POI
 // correction terms for nua:
 for(Int_t t=0;t<2;t++) // typeFlag (0 = RP, 1 = POI)
 { 
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
 // differential flow:
 TString diffFlowName = "fDiffFlow";
 diffFlowName += fAnalysisLabel->Data();
 for(Int_t t=0;t<2;t++) // type: RP or POI
 { 
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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

/*
void AliFlowAnalysisWithQCumulants::CalculateCorrectionsForNUAForIntQcumulants() // to be improved (do I really need this method?)
{
 // Calculate final corrections for non-uniform acceptance for Q-cumulants.
  
 // Corrections for non-uniform acceptance are stored in histogram fCorrectionsForNUA,
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


} // end of AliFlowAnalysisWithQCumulants::CalculateCorrectionsForNUAForIntQcumulants()
*/

//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateQcumulantsCorrectedForNUAIntFlow()
{
 // Calculate generalized Q-cumulants (cumulants corrected for non-unifom acceptance).
 
 // measured 2-, 4-, 6- and 8-particle correlations (biased by non-uniform acceptance!):
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>>
 Double_t four = fIntFlowCorrelationsHist->GetBinContent(2); // <<4>>
 //Double_t six = fIntFlowCorrelationsHist->GetBinContent(3); // <<6>>
 //Double_t eight = fIntFlowCorrelationsHist->GetBinContent(4); // <<8>>
 
 // statistical error of measured 2-, 4-, 6- and 8-particle correlations:
 //Double_t twoError = fIntFlowCorrelationsHist->GetBinError(1); // statistical error of <<2>>
 //Double_t fourError = fIntFlowCorrelationsHist->GetBinError(2); // statistical error of <<4>>
 //Double_t sixError = fIntFlowCorrelationsHist->GetBinError(3); // statistical error of <<6>>
 //Double_t eightError = fIntFlowCorrelationsHist->GetBinError(4); // statistical error of <<8>>

 // QC{2}:
 // <<cos(n*phi1)>>^2
 Double_t two1stTerm = pow(fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(1),2); 
 //Double_t two1stTermErrorSquared = pow(fIntFlowCorrectionTermsForNUAHist[1]->GetBinError(1),2); 
 // <<sin(n*phi1)>>^2
 Double_t two2ndTerm = pow(fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(1),2); 
 //Double_t two2ndTermErrorSquared = pow(fIntFlowCorrectionTermsForNUAHist[0]->GetBinError(1),2); 
 // generalized QC{2}:
 Double_t gQC2 = two - two1stTerm - two2ndTerm; // to be improved (terminology, notation)
 fIntFlowQcumulants->SetBinContent(1,gQC2); 
 //fIntFlowQcumulants->SetBinError(1,0.); // to be improved (propagate error) 
 
 // QC{4}:
 // <<cos(n*phi1)>> <<cos(n*(phi1-phi2-phi3))>>
 Double_t four1stTerm = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(1)
                      * fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(3);  
 // <<sin(n*phi1)>> <<sin(n*(phi1-phi2-phi3))>>
 Double_t four2ndTerm = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(1)
                      * fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(3);  
 // <<cos(n*(phi1+phi2))>>^2
 Double_t four3rdTerm = pow(fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(2),2); 
 // <<sin(n*(phi1+phi2))>>^2
 Double_t four4thTerm = pow(fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(2),2); 
 // <<cos(n*(phi1+phi2))>> (<<cos(n*phi1)>>^2 - <<sin(n*phi1)>>^2)
 Double_t four5thTerm = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(2)
                      * (pow(fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(1),2)
                      - pow(fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(1),2));
 // <<sin(n*(phi1+phi2))>> <<cos(n*phi1)>> <<sin(n*phi1)>>
 Double_t four6thTerm = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(2)
                      * fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(1)
                      * fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(1);         
 // <<cos(n*(phi1-phi2))>> (<<cos(n*phi1)>>^2 + <<sin(n*phi1)>>^2)
 Double_t four7thTerm = two*(pow(fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(1),2)
                      + pow(fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(1),2));  
 // (<<cos(n*phi1)>>^2 + <<sin(n*phi1)>>^2)^2
 Double_t four8thTerm = pow(pow(fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(1),2)
                      + pow(fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(1),2),2);      
 // generalized QC{4}:
 Double_t gQC4 = four-2.*pow(two,2.)-4.*four1stTerm+4.*four2ndTerm-four3rdTerm
               - four4thTerm+4.*four5thTerm+8.*four6thTerm+8.*four7thTerm-6.*four8thTerm;                            
 fIntFlowQcumulants->SetBinContent(2,gQC4);   
 //fIntFlowQcumulants->SetBinError(2,0.); // to be improved (propagate error) 

 // ... to be improved (continued for 6th and 8th order)                                                    
    
} // end of void AliFlowAnalysisWithQCumulants::CalculateQcumulantsCorrectedForNUAIntFlow()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectedForNUA()
{
 // Calculate integrated flow from generalized Q-cumulants (corrected for non-uniform acceptance).
 
 // to be improved: add protection for NULL pointers, propagate statistical errors from 
 // measured correlations and correction terms
 
 // generalized Q-cumulants:
 Double_t qc2 = fIntFlowQcumulants->GetBinContent(1); // QC{2}  
 Double_t qc4 = fIntFlowQcumulants->GetBinContent(2); // QC{4}  
 //Double_t qc6 = fIntFlowQcumulants->GetBinContent(3); // QC{6}  
 //Double_t qc8 = fIntFlowQcumulants->GetBinContent(4); // QC{8}
 
 // integrated flow estimates:
 Double_t v2 = 0.; // v{2,QC}  
 Double_t v4 = 0.; // v{4,QC}  
 //Double_t v6 = 0.; // v{6,QC}  
 //Double_t v8 = 0.; // v{8,QC}

 // calculate integrated flow estimates from generalized Q-cumulants: 
 if(qc2>=0.) v2 = pow(qc2,1./2.); 
 if(qc4<=0.) v4 = pow(-1.*qc4,1./4.); 
 //if(qc6>=0.) v6 = pow((1./4.)*qc6,1./6.); 
 //if(qc8<=0.) v8 = pow((-1./33.)*qc8,1./8.); 

 // store integrated flow estimates from generalized Q-cumulants:
 fIntFlow->SetBinContent(1,v2);
 fIntFlow->SetBinContent(2,v4);
 //fIntFlow->SetBinContent(3,v6);
 //fIntFlow->SetBinContent(4,v8);
 
 /*
 // propagate correctly error by including non-isotropic terms (to be improved - moved somewhere else):
 // correlations:
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>> 
 //Double_t four = fIntFlowCorrelationsHist->GetBinContent(2); // <<4>>  
 //Double_t six = fIntFlowCorrelationsHist->GetBinContent(3); // <<6>> 
 //Double_t eight = fIntFlowCorrelationsHist->GetBinContent(4); // <<8>>  
 // statistical errors of average 2-, 4-, 6- and 8-particle azimuthal correlations:
 Double_t twoError = fIntFlowCorrelationsHist->GetBinError(1); // statistical error of <2>  
 Double_t fourError = fIntFlowCorrelationsHist->GetBinError(2); // statistical error of <4>   
 //Double_t sixError = fIntFlowCorrelationsHist->GetBinError(3); // statistical error of <6> 
 //Double_t eightError = fIntFlowCorrelationsHist->GetBinError(4); // statistical error of <8> 
 // nua terms: 
 Double_t c1 = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(1); // <<cos(phi1)>>
 Double_t c2 = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(2); // <<cos(phi1+phi2)>>
 Double_t c3 = fIntFlowCorrectionTermsForNUAHist[1]->GetBinContent(3); // <<cos(phi1-phi2-phi3)>>
 Double_t s1 = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(1); // <<sin(phi1)>>
 Double_t s2 = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(2); // <<sin(phi1+phi2)>>
 Double_t s3 = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(3); // <<sin(phi1-phi2-phi3)>>
 // statistical errors of nua terms:
 Double_t c1Error = fIntFlowCorrectionTermsForNUAHist[1]->GetBinError(1); // statistical error of <cos(phi1)>
 Double_t c2Error = fIntFlowCorrectionTermsForNUAHist[1]->GetBinError(2); // statistical error of <cos(phi1+phi2)>
 Double_t c3Error = fIntFlowCorrectionTermsForNUAHist[1]->GetBinError(3); // statistical error of <cos(phi1-phi2-phi3)>
 Double_t s1Error = fIntFlowCorrectionTermsForNUAHist[0]->GetBinError(1); // statistical error of <sin(phi1)>
 Double_t s2Error = fIntFlowCorrectionTermsForNUAHist[0]->GetBinError(2); // statistical error of <sin(phi1+phi2)>
 Double_t s3Error = fIntFlowCorrectionTermsForNUAHist[0]->GetBinError(3); // statistical error of <sin(phi1-phi2-phi3)>
 
 // covariances for nua:
 Double_t wCov24 = fIntFlowCovariances->GetBinContent(1); // Cov(<2>,<4>) * prefactor(w_<2>,w_<4>)
 Double_t wCov2c1 = fIntFlowCovariancesNUA->GetBinContent(1);
 Double_t wCov2s1 = fIntFlowCovariancesNUA->GetBinContent(2);
 Double_t wCovc1s1 = fIntFlowCovariancesNUA->GetBinContent(3);
 Double_t wCov2c2 = fIntFlowCovariancesNUA->GetBinContent(4);
 Double_t wCov2s2 = fIntFlowCovariancesNUA->GetBinContent(5); 
 Double_t wCov2c3 = fIntFlowCovariancesNUA->GetBinContent(6); 
 Double_t wCov2s3 = fIntFlowCovariancesNUA->GetBinContent(7); 
 Double_t wCov4c1 = fIntFlowCovariancesNUA->GetBinContent(8);  
 Double_t wCov4s1 = fIntFlowCovariancesNUA->GetBinContent(9); 
 Double_t wCov4c2 = fIntFlowCovariancesNUA->GetBinContent(10); 
 Double_t wCov4s2 = fIntFlowCovariancesNUA->GetBinContent(11); 
 Double_t wCov4c3 = fIntFlowCovariancesNUA->GetBinContent(12);
 Double_t wCov4s3 = fIntFlowCovariancesNUA->GetBinContent(13);
 Double_t wCovc1c2 = fIntFlowCovariancesNUA->GetBinContent(14);
 Double_t wCovc1s2 = fIntFlowCovariancesNUA->GetBinContent(15);
 Double_t wCovc1c3 = fIntFlowCovariancesNUA->GetBinContent(16);
 Double_t wCovc1s3 = fIntFlowCovariancesNUA->GetBinContent(17);
 Double_t wCovs1c2 = fIntFlowCovariancesNUA->GetBinContent(18);
 Double_t wCovs1s2 = fIntFlowCovariancesNUA->GetBinContent(19); 
 Double_t wCovs1c3 = fIntFlowCovariancesNUA->GetBinContent(20);
 Double_t wCovs1s3 = fIntFlowCovariancesNUA->GetBinContent(21); 
 Double_t wCovc2s2 = fIntFlowCovariancesNUA->GetBinContent(22);  
 Double_t wCovc2c3 = fIntFlowCovariancesNUA->GetBinContent(23);  
 Double_t wCovc2s3 = fIntFlowCovariancesNUA->GetBinContent(24); 
 Double_t wCovs2c3 = fIntFlowCovariancesNUA->GetBinContent(25); 
 Double_t wCovs2s3 = fIntFlowCovariancesNUA->GetBinContent(26);
 Double_t wCovc3s3 = fIntFlowCovariancesNUA->GetBinContent(27); 
 */
 
 /*
 // 2nd order:
 Double_t err2ndSquared = (1./(4.*pow(v2,2.)))
                        * (pow(twoError,2.)+4.*pow(s1*s1Error,2.)+4.*pow(c1*c1Error,2.)
                        // to be improved (add eventually also covariance terms)
                        //- 4.*c1*wCov2c1-4.*s1*wCov2s1+8.*c1*s1*wCovc1s1 
                        );
 if(err2ndSquared>=0.) 
 {
  fIntFlow->SetBinError(1,pow(err2ndSquared,0.5)); // to be improved (enabled eventually)
 } else
   {
    cout<<"WARNING: Statistical error of v{2,QC} (with non-isotropic terms included) is imaginary !!!! "<<endl;
   }
 // 4th order:
 Double_t err4thSquared = (1./(16.*pow(v4,6.)))
                        * (pow(4.*pow(two,2.)-8.*(pow(c1,2.)+pow(s1,2.)),2.)*pow(twoError,2.)
                        + pow(fourError,2.)
                        + 16.*pow(6.*pow(c1,3.)-2.*c1*c2+c3+6.*c1*pow(s1,2.)-2.*s1*s2-4.*c1*two,2.)*pow(c1Error,2.)
                        + 16.*pow(6.*pow(c1,2.)*s1+2.*c2*s1+6.*pow(s1,3.)-2.*c1*s2-s3-4.*s1*two,2.)*pow(s1Error,2.)
                        + 4.*pow(c2-2.*(pow(c1,2.)-pow(s1,2.)),2.)*pow(c2Error,2.)
                        + 4.*pow(4*c1*s1-s2,2.)*pow(s2Error,2.)
                        + 16.*pow(c1,2.)*pow(c3Error,2.)
                        + 16.*pow(s1,2.)*pow(s3Error,2.)
                        // to be improved (add eventually also covariance terms)
                        // ...
                        );
 if(err4thSquared>=0.) 
 {
  fIntFlow->SetBinError(2,pow(err4thSquared,0.5)); // to be improved (enabled eventually)
 } else
   {
    cout<<"WARNING: Statistical error of v{4,QC} (with non-isotropic terms included) is imaginary !!!! "<<endl;
   }
 */
 
} // end of void AliFlowAnalysisWithQCumulants::CalculateIntFlowCorrectedForNUA()

   
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
 
 for(Int_t sc=0;sc<2;sc++) // sin or cos correction terms 
 {
  for(Int_t ci=1;ci<=3;ci++) // correction term index
  {
   Double_t correction = fIntFlowCorrectionTermsForNUAPro[sc]->GetBinContent(ci);
   Double_t spread = fIntFlowCorrectionTermsForNUAPro[sc]->GetBinError(ci);
   Double_t sumOfLinearEventWeights = fIntFlowSumOfEventWeightsNUA[sc][0]->GetBinContent(ci);
   Double_t sumOfQuadraticEventWeights = fIntFlowSumOfEventWeightsNUA[sc][1]->GetBinContent(ci);
   Double_t termA = 0.;
   Double_t termB = 0.;
   if(sumOfLinearEventWeights)
   {
    termA = pow(sumOfQuadraticEventWeights,0.5)/sumOfLinearEventWeights;
   } else
     {
      cout<<"WARNING: sumOfLinearEventWeights == 0 in AFAWQC::FCTFNIF() !!!!"<<endl;
      cout<<"         (for "<<ci<<"-th correction term)"<<endl;
     }
   if(1.-pow(termA,2.) > 0.)
   {
    termB = 1./pow(1-pow(termA,2.),0.5);
   } else
     {
      cout<<"WARNING: 1.-pow(termA,2.) <= 0 in AFAWQC::FCTFNIF() !!!!"<<endl;   
      cout<<"         (for "<<ci<<"-th correction term)"<<endl;
     }     
   Double_t statisticalError = termA * spread * termB;
   fIntFlowCorrectionTermsForNUAHist[sc]->SetBinContent(ci,correction);
   fIntFlowCorrectionTermsForNUAHist[sc]->SetBinError(ci,statisticalError);
  } // end of for(Int_t ci=1;ci<=10;ci++) // correction term index
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
   if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)  
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
   } // end of if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)  
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
 (fCommonHists2nd->GetHarmonic())->Fill(0.5,fHarmonic);
 (fCommonHists4th->GetHarmonic())->Fill(0.5,fHarmonic);
 (fCommonHists6th->GetHarmonic())->Fill(0.5,fHarmonic);
 (fCommonHists8th->GetHarmonic())->Fill(0.5,fHarmonic);

} // end of void AliFlowAnalysisWithQCumulants::StoreHarmonic()


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrelationsUsingParticleWeights(TString type, TString ptOrEta) // type = RP or POI 
{
 // Calculate all correlations needed for differential flow using particle weights.
 
 Int_t t = -1; // type flag 
 Int_t pe = -1; // ptEta flag
 
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
 
 // S^M_{p,k} (see .h file for the definition of fSMpk):
 Double_t dSM1p1k = (*fSMpk)(0,1);
 Double_t dSM1p2k = (*fSMpk)(0,2);
 Double_t dSM1p3k = (*fSMpk)(0,3);
 Double_t dSM2p1k = (*fSMpk)(1,1);
 Double_t dSM3p1k = (*fSMpk)(2,1);
 
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
 
} // end of void AliFlowAnalysisWithQCumulants::FillCommonControlHistograms(AliFlowEventSimple *anEvent)


//================================================================================================================================


void AliFlowAnalysisWithQCumulants::ResetEventByEventQuantities()
{
 // Reset all event by event quantities.
 
 // integrated flow:
 fReQ->Zero();
 fImQ->Zero();
 fSMpk->Zero();
 fIntFlowCorrelationsEBE->Reset();
 fIntFlowEventWeightsForCorrelationsEBE->Reset();
 fIntFlowCorrelationsAllEBE->Reset();
 
 if(fApplyCorrectionForNUA)  
 {
  for(Int_t sc=0;sc<2;sc++)
  {
   fIntFlowCorrectionTermsForNUAEBE[sc]->Reset();
   fIntFlowEventWeightForCorrectionTermsForNUAEBE[sc]->Reset();
  } 
 }
    
 // differential flow:
 // 1D:
 for(Int_t t=0;t<3;t++) // type (RP, POI, POI&&RP)
 {
  for(Int_t pe=0;pe<2;pe++) // 1D in pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // 1D in pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
  for(Int_t pe=0;pe<2;pe++) // pt or eta
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
    
 // 2D (pt,eta)
 if(fCalculate2DFlow)
 {
  for(Int_t t=0;t<3;t++) // type (RP, POI, POI&&RP)
  {
   for(Int_t m=0;m<4;m++) // multiple of harmonic
   {
    for(Int_t k=0;k<9;k++) // power of weight
    {
     if(fReRPQ2dEBE[t][m][k]) fReRPQ2dEBE[t][m][k]->Reset();
     if(fImRPQ2dEBE[t][m][k]) fImRPQ2dEBE[t][m][k]->Reset();
    }   
   }
  }
  for(Int_t t=0;t<3;t++) // type (0 = RP, 1 = POI, 2 = RP&&POI )
  { 
   for(Int_t k=0;k<9;k++)
   {
    if(fs2dEBE[t][k]) fs2dEBE[t][k]->Reset();
   }
  }  
 } // end of if(fCalculate2DFlow) 

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

 Int_t t = -1; // type flag 
 Int_t pe = -1; // ptEta flag
 
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

 Int_t t = -1; // type flag 
 Int_t pe = -1; // ptEta flag
 
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
 // Transfer prolfiles into histogams and correctly propagate the error (to be improved: description)
 
 // to be improved: debugged - I do not correctly transfer all profiles into histos (bug appears only after merging) 
  
 Int_t t = -1; // type flag 
 Int_t pe = -1; // ptEta flag
 
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
 // Calculate generalized differential flow Q-cumulants (corrected for non-uniform acceptance)
  
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;

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
 
 // 2-particle correlation:
 Double_t two = fIntFlowCorrelationsHist->GetBinContent(1); // <<2>>
 // sin term coming from integrated flow: 
 Double_t sinP1nPhi = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(1); // <<sin(n*phi1)>>
 Double_t sinP1nPhi1P1nPhi2 = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(2); // <<sin(n*(phi1+phi2))>>
 Double_t sinP1nPhi1M1nPhi2M1nPhi3 = fIntFlowCorrectionTermsForNUAHist[0]->GetBinContent(3); // <<sin(n*(phi1-phi2-phi3))>>
 // cos term coming from integrated flow: 
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
  // generalized QC{2'}:
  Double_t qc2Prime = twoPrime - sinP1nPsi*sinP1nPhi - cosP1nPsi*cosP1nPhi;
  fDiffFlowCumulants[t][pe][0]->SetBinContent(b,qc2Prime);
  // generalized QC{4'}:
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
  fDiffFlowCumulants[t][pe][1]->SetBinContent(b,qc4Prime);   
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)
 
} // end of AliFlowAnalysisWithQCumulants::CalculateDiffFlowCumulantsCorrectedForNUA(TString type, TString ptOrEta)


//==================================================================================================================================
    

void AliFlowAnalysisWithQCumulants::CalculateDiffFlowCorrectedForNUA(TString type, TString ptOrEta)
{
 // Calculate differential flow corrected for non-uniform acceptance.
 
 // to be improved (rewritten completely)
 
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;

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
   
 // to be improved: access here generalized QC{2} and QC{4} instead: 
 Double_t dV2 = fIntFlow->GetBinContent(1); 
 Double_t dV4 = fIntFlow->GetBinContent(2); 
 
 // loop over pt or eta bins:
 for(Int_t b=1;b<=nBinsPtEta[pe];b++)
 {
  // generalized QC{2'}:
  Double_t gQC2Prime = fDiffFlowCumulants[t][pe][0]->GetBinContent(b);
  // v'{2}:
  if(dV2>0)
  { 
   Double_t v2Prime = gQC2Prime/dV2;
   fDiffFlow[t][pe][0]->SetBinContent(b,v2Prime); 
  }  
  // generalized QC{4'}:
  Double_t gQC4Prime = fDiffFlowCumulants[t][pe][1]->GetBinContent(b);
  // v'{4}:
  if(dV4>0)
  { 
   Double_t v4Prime = -gQC4Prime/pow(dV4,3.);
   fDiffFlow[t][pe][1]->SetBinContent(b,v4Prime); 
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
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL; 
 Double_t phi1=0., phi2=0., phi3=0., phi4=0., phi5=0., phi6=0., phi7=0., phi8=0.; 
 Int_t n = fHarmonic; 
 Int_t eventNo = (Int_t)fAvMultiplicity->GetBinEntries(1); // to be improved (is this casting safe in general?)
 Double_t dMult = (*fSMpk)(0,0);
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
    fIntFlowDirectCorrelations->Fill(0.5,cos(n*(phi1-phi2)),1.);    // <cos(n*(phi1-phi2))>
    fIntFlowDirectCorrelations->Fill(1.5,cos(2.*n*(phi1-phi2)),1.); // <cos(2n*(phi1-phi2))>
    fIntFlowDirectCorrelations->Fill(2.5,cos(3.*n*(phi1-phi2)),1.); // <cos(3n*(phi1-phi2))>
    fIntFlowDirectCorrelations->Fill(3.5,cos(4.*n*(phi1-phi2)),1.); // <cos(4n*(phi1-phi2))>   
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
     fIntFlowDirectCorrelations->Fill(5.,cos(2.*n*phi1-n*(phi2+phi3)),1.);       //<3>_{2n|nn,n}
     fIntFlowDirectCorrelations->Fill(6.,cos(3.*n*phi1-2.*n*phi2-n*phi3),1.);    //<3>_{3n|2n,n}
     fIntFlowDirectCorrelations->Fill(7.,cos(4.*n*phi1-2.*n*phi2-2.*n*phi3),1.); //<3>_{4n|2n,2n}
     fIntFlowDirectCorrelations->Fill(8.,cos(4.*n*phi1-3.*n*phi2-n*phi3),1.);    //<3>_{4n|3n,n}
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
       fIntFlowDirectCorrelations->Fill(18.,cos(2.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5),1.);       //<5>_{2n,n|n,n,n}
       fIntFlowDirectCorrelations->Fill(19.,cos(2.*n*phi1+2.*n*phi2-2.*n*phi3-n*phi4-n*phi5),1.); //<5>_{2n,2n|2n,n,n}
       fIntFlowDirectCorrelations->Fill(20.,cos(3.*n*phi1+n*phi2-2.*n*phi3-n*phi4-n*phi5),1.);    //<5>_{3n,n|2n,n,n}
       fIntFlowDirectCorrelations->Fill(21.,cos(4.*n*phi1-n*phi2-n*phi3-n*phi4-n*phi5),1.);       //<5>_{4n|n,n,n,n}
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
        fIntFlowDirectCorrelations->Fill(23.,cos(n*phi1+n*phi2+n*phi3-n*phi4-n*phi5-n*phi6),1.);       //<6>_{n,n,n|n,n,n}
        fIntFlowDirectCorrelations->Fill(24.,cos(2.*n*phi1+n*phi2+n*phi3-2.*n*phi4-n*phi5-n*phi6),1.); //<6>_{2n,n,n|2n,n,n}
        fIntFlowDirectCorrelations->Fill(25.,cos(2.*n*phi1+2.*n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),1.); //<6>_{2n,2n|n,n,n,n}
        fIntFlowDirectCorrelations->Fill(26.,cos(3.*n*phi1+n*phi2-n*phi3-n*phi4-n*phi5-n*phi6),1.);    //<6>_{3n,n|n,n,n,n}  
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
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
 {
  cout<<"   ****   (particle weights not used)   ****"<<endl;
 } else
   {
    cout<<"   ****     (particle weights used)     ****"<<endl;
   } 
 cout<<"   *****************************************"<<endl;
 cout<<endl;
 cout<<endl;

 Int_t ciMax = 32; // to be improved (removed eventually when I calculate 6th and 8th order with particle weights)
 
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
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
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
 {
  cout<<"   ****     (particle weights not used)     ****"<<endl;
 } else
   {
    cout<<"   ****       (particle weights used)       ****"<<endl;
   } 
 cout<<"   *********************************************"<<endl;
 cout<<endl;
 cout<<endl;

 for(Int_t ci=1;ci<=10;ci++) // correction term index
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
 Double_t dMult = (*fSMpk)(0,0);
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
 Double_t dMult = (*fSMpk)(0,0);
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
 
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;
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


void AliFlowAnalysisWithQCumulants::CrossCheckDiffFlowCorrelations(TString type, TString ptOrEta)
{
 // Compare correlations needed for diff. flow calculated with nested loops and those calculated from Q-vectors
 
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;
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
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
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
  
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;
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
  
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;
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
 
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;
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
 if(!(fUsePhiWeights||fUsePtWeights||fUseEtaWeights))
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
 Double_t dMult = (*fSMpk)(0,0);
 
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
 Double_t dM11 = (*fSMpk)(1,1)-(*fSMpk)(0,2); // dM11 = sum_{i,j=1,i!=j}^M w_i w_j
 Double_t dM111 = (*fSMpk)(2,1)-3.*(*fSMpk)(0,2)*(*fSMpk)(0,1)
                + 2.*(*fSMpk)(0,3); // dM111 = sum_{i,j,k=1,i!=j!=k}^M w_i w_j w_k
 //..............................................................................................
         // 1-particle:
 Double_t cosP1nW1 = 0.; // <<w1 cos(n*(phi1))>>
   
 if(dMult>0 && TMath::Abs((*fSMpk)(0,1))>1e-6)
 {
  cosP1nW1 = dReQ1n1k/(*fSMpk)(0,1); 
  
  // average weighted 1-particle correction (cos terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[1]->SetBinContent(1,cosP1nW1);
  
  // final average weighted 1-particle correction (cos terms) for non-uniform acceptance for all events:
  fIntFlowCorrectionTermsForNUAPro[1]->Fill(0.5,cosP1nW1,(*fSMpk)(0,1));  
 } 
 
 // 2-particle:
 Double_t cosP1nP1nW1W1 = 0.; // <<w1 w2 cos(n*(phi1+phi2))>>
 
 if(dMult>1 && TMath::Abs(dM11)>1e-6)
 {
  cosP1nP1nW1W1 = (pow(dReQ1n1k,2)-pow(dImQ1n1k,2)-dReQ2n2k)/dM11; 
  
  // average weighted 2-particle correction (cos terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[1]->SetBinContent(2,cosP1nP1nW1W1);
  
  // final average weighted 2-particle correction (cos terms) for non-uniform acceptance for all events:
  fIntFlowCorrectionTermsForNUAPro[1]->Fill(1.5,cosP1nP1nW1W1,dM11);  
 } 
 
 // 3-particle:
 Double_t cosP1nM1nM1nW1W1W1 = 0.; // <<w1 w2 w3 cos(n*(phi1-phi2-phi3))>>
 
 if(dMult>2 && TMath::Abs(dM111)>1e-6)
 {
  cosP1nM1nM1nW1W1W1 = (dReQ1n1k*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2))
                     - dReQ1n1k*dReQ2n2k-dImQ1n1k*dImQ2n2k
                     - 2.*((*fSMpk)(0,2))*dReQ1n1k
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
 Double_t dMult = (*fSMpk)(0,0);
 
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
 Double_t dM11 = (*fSMpk)(1,1)-(*fSMpk)(0,2); // dM11 = sum_{i,j=1,i!=j}^M w_i w_j
 Double_t dM111 = (*fSMpk)(2,1)-3.*(*fSMpk)(0,2)*(*fSMpk)(0,1)
                + 2.*(*fSMpk)(0,3); // dM111 = sum_{i,j,k=1,i!=j!=k}^M w_i w_j w_k
 //..............................................................................................
 
 // 1-particle:
 Double_t sinP1nW1 = 0.; // <<w1 sin(n*(phi1))>>
 
 if(dMult>0 && TMath::Abs((*fSMpk)(0,1))>1e-6)
 {
  sinP1nW1 = dImQ1n1k/((*fSMpk)(0,1)); 
     
  // average weighted 1-particle correction (sin terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[0]->SetBinContent(1,sinP1nW1);
  
  // final average weighted 1-particle correction (sin terms) for non-uniform acceptance for all events:   
  fIntFlowCorrectionTermsForNUAPro[0]->Fill(0.5,sinP1nW1,(*fSMpk)(0,1));  
 } 
 
 // 2-particle:
 Double_t sinP1nP1nW1W1 = 0.; // <<w1 w2 sin(n*(phi1+phi2))>>
 
 if(dMult>1 && TMath::Abs(dM11)>1e-6)
 {
  sinP1nP1nW1W1 = (2.*dReQ1n1k*dImQ1n1k-dImQ2n2k)/dM11; 
     
  // average weighted 2-particle correction (sin terms) for non-uniform acceptance for single event:
  fIntFlowCorrectionTermsForNUAEBE[0]->SetBinContent(2,sinP1nP1nW1W1);
  
  // final average weighted 1-particle correction (sin terms) for non-uniform acceptance for all events:      
  fIntFlowCorrectionTermsForNUAPro[0]->Fill(1.5,sinP1nP1nW1W1,dM11);  
 } 
 
 // 3-particle:
 Double_t sinP1nM1nM1nW1W1W1 = 0.; // <<w1 w2 w3 sin(n*(phi1-phi2-phi3))>>
 
 if(dMult>2 && TMath::Abs(dM111)>1e-6)
 {
  sinP1nM1nM1nW1W1W1 = (-dImQ1n1k*(pow(dReQ1n1k,2)+pow(dImQ1n1k,2))
                     + dReQ1n1k*dImQ2n2k-dImQ1n1k*dReQ2n2k
                     + 2.*((*fSMpk)(0,2))*dImQ1n1k
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
 Double_t dMult = (*fSMpk)(0,0);
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
 
 // S^M_{p,k} (see .h file for the definition of fSMpk):
 Double_t dSM1p1k = (*fSMpk)(0,1);
 Double_t dSM1p2k = (*fSMpk)(0,2);
 Double_t dSM2p1k = (*fSMpk)(1,1);

 Int_t t = -1; // type flag 
 Int_t pe = -1; // ptEta flag
 
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
 
 // S^M_{p,k} (see .h file for the definition of fSMpk):
 Double_t dSM1p1k = (*fSMpk)(0,1);
 Double_t dSM1p2k = (*fSMpk)(0,2);
 Double_t dSM2p1k = (*fSMpk)(1,1);

 Int_t t = -1; // type flag 
 Int_t pe = -1; // ptEta flag
 
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
     
 Int_t typeFlag = -1;
 Int_t ptEtaFlag = -1;
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





