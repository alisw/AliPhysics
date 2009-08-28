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

/******************************** 
 * integrated flow estimate by  *
 *   fitting q-distribution     * 
 *                              *
 * author: Ante Bilandzic       * 
 *          (anteb@nikhef.nl)   *
 *                              *  
 *  based on the macro written  *
 *     by Sergei Voloshin       *
 *******************************/  

#define AliFlowAnalysisWithFittingQDistribution_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h" 
#include "TF1.h"
#include "TLegend.h"
#include "TParticle.h"
#include "TProfile.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithFittingQDistribution.h"

class TH1;
class TGraph;
class TPave;
class TLatex;
class TMarker;
class TObjArray;
class TList;
class TSystem;
class TROOT;
class AliFlowVector;
class TVector;

//================================================================================================================

ClassImp(AliFlowAnalysisWithFittingQDistribution)

AliFlowAnalysisWithFittingQDistribution::AliFlowAnalysisWithFittingQDistribution():  
 fHistList(NULL),
 fCommonHists(NULL),
 fCommonHistsResults(NULL),
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
 fWeightsList(NULL),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fUseParticleWeights(NULL),
 fPhiWeights(NULL),
 fPtWeights(NULL),
 fEtaWeights(NULL),
 // fitting parameters with default values harwired here (use dedicated macro fqd.C to change them):
 fFittingParameters(NULL), 
 fTreshold(5),
 fvStart(0.05),
 fvMin(0.0),
 fvMax(0.25),
 fSigma2Start(0.75),
 fSigma2Min(0.5), 
 fSigma2Max(2.5),
 fPlotResults(kFALSE),
 // rest:
 fLegend(NULL)
 {
  // constructor 
  
  // base list to hold all output objects:
  fHistList = new TList();
  fHistList->SetName("cobjFQD");
  fHistList->SetOwner(kTRUE);
  
  // analysis label;
  fAnalysisLabel = new TString();
 
  // list to hold histograms with phi, pt and eta weights:      
  fWeightsList = new TList();

  // initialize all arrays:  
  this->InitializeArrays();

 } // end of constructor
 
 
//================================================================================================================
 

AliFlowAnalysisWithFittingQDistribution::~AliFlowAnalysisWithFittingQDistribution()
{
 // desctructor
 delete fHistList; 
}


//================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::Init()
{
 // access constants and book everything
 
 // access constants:
 this->AccessConstants();
 
 // booking:
 this->BookCommonHistograms();
 this->BookAndFillWeightsHistograms();
 this->BookEverythingForDistributions();
 
 // store fitting parameters:  
 this->StoreFittingParameters();
 
 // nest lists:
 fWeightsList->SetName("Weights");
 fWeightsList->SetOwner(kTRUE);   
 fHistList->Add(fWeightsList);
 
 // set harmonic in common control histograms (to be improved (should I do this somewhere else?)):
 (fCommonHists->GetHarmonic())->Fill(0.5,fHarmonic); 
 
} // end of void AliFlowAnalysisWithFittingQDistribution::Init()


//================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::Make(AliFlowEventSimple* anEvent)
{
 // loop over data
 
 // a) fill the common control histograms;
 // b) loop over data and calculate non-weighted and weighted Q-vector and sum of particle weights;
 // c) fill histograms for q-distribution;
 // d) reset e-b-e quantities.
 
 // a) fill the common control histograms:
 fCommonHists->FillControlHistograms(anEvent); 
 
 Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
 Double_t dPt  = 0.; // transverse momentum
 Double_t dEta = 0.; // pseudorapidity

 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight
 
 Int_t nPrim = anEvent->NumberOfTracks();  // nPrim = total number of primary tracks, i.e. nPrim = nRP + nPOI + rest, where:
                                           // nRP   = # of particles used to determine the reaction plane;
                                           // nPOI  = # of particles of interest for a detailed flow analysis;
                                           // rest  = # of particles which are not niether RPs nor POIs.  
 
 Int_t n = fHarmonic; // shortcut for the harmonic
 
 Double_t dReQ[2] = {0.}; // real part of Q-vector [0=particle weights not used, 1=particle weights used]
 Double_t dImQ[2] = {0.}; // imaginary part of Q-vector [0=particle weights not used, 1=particle weights used]
 Double_t dSumOfParticleWeights[2] = {0.}; // [0=particle weights not used, 1=particle weights used] 
                                                                                                                               
 AliFlowTrackSimple *aftsTrack = NULL;                                          
 
 // b) loop over data and calculate non-weighted and weighted Q-vector and sum of particle weights:                                          
 for(Int_t i=0;i<nPrim;i++) 
 { 
  aftsTrack=anEvent->GetTrack(i);
  if(aftsTrack)
  {
   if(!(aftsTrack->InRPSelection())) continue; // consider only tracks which are RPs
     
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
    
   // calculate real and imaginary part of non-weighted and weighted Q-vector and sum of particle weights for this event:
   for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++) // particle weights not used (0) or used (1)
   {
    // Q-vector:
    dReQ[pW]+=pow(wPhi*wPt*wEta,pW)*TMath::Cos(n*dPhi); 
    dImQ[pW]+=pow(wPhi*wPt*wEta,pW)*TMath::Sin(n*dPhi);
    // sum of particle weights:
    dSumOfParticleWeights[pW] += pow(wPhi*wPt*wEta,pW); // if pW = 0, this sum gives # of RPs, i.e. multiplicity
   } 
   
  } // end of if(aftsTrack)
 } // end of for(Int_t i=0;i<nPrim;i++)                                      
                                           
 // c) fill histograms for q-distribution:
 // calculating first q = Q\sqrt{sum of particle weights} (Remark: if particle weights are unit than sum of particle weights = multiplicity)
 Double_t q=0;                                          
 for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++) // particle weights not used (0) or used (1)
 {
  if(dSumOfParticleWeights[pW])
  {
   q = pow(dReQ[pW]*dReQ[pW]+dImQ[pW]*dImQ[pW],0.5)/pow(dSumOfParticleWeights[pW],0.5);
   // fill histograms:
   fqDistribution[pW]->Fill(q,1.);
   fSumOfParticleWeights[pW]->Fill(dSumOfParticleWeights[pW],1.);
  }
 } 
 
 // d) reset e-b-e quantities:
 for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++) // particle weights not used (0) or used (1)
 {
  dReQ[pW] = 0.;
  dImQ[pW] = 0.;
  dSumOfParticleWeights[pW] = 0.;
 }

} // end of Make()


//================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::Finish(Bool_t doFit)
{
 // calculate the final results
 
 // a) acces the constants and all fitting paremeters;
 // b) access the flags for particle weights;
 // c) do final fit;
 // d) fill common hist results;
 // e) print on the screen the final results.
 
 // a) access the constants and all fitting paremeters:
 this->AccessConstants();
 this->AccessFittingParameters();
 
 // b) access the flags for particle weights: 
 fUsePhiWeights = (Bool_t)fUseParticleWeights->GetBinContent(1); 
 fUsePtWeights = (Bool_t)fUseParticleWeights->GetBinContent(2); 
 fUseEtaWeights = (Bool_t)fUseParticleWeights->GetBinContent(3);

 // to be improved (moved somewhere else):
 if(fPlotResults)
 {
  fLegend = new TLegend(0.6,0.55,0.85,0.7); 
 }

 // c) do final fit:             
 if(doFit) 
 {
  // particle weights not used:
  // 1) sigma^2 not fitted (fixed to 0.5):
  this->DoFit(kFALSE,kFALSE);
  // 2) sigma^2 fitted:
  this->DoFit(kFALSE,kTRUE);
  // particle weights used:
  if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)      
  {
   // 1) sigma^2 not fitted (fixed to 0.5):
   this->DoFit(kTRUE,kFALSE);  
   // 2) sigma^2 fitted:
   this->DoFit(kTRUE,kTRUE);  
  }
  
  // d) fill common hist results (by default fill results obtained with sigma^2 fitted):
  if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
  {
   this->FillCommonHistResultsIntFlow(kTRUE,kTRUE); 
  } else
    {
     this->FillCommonHistResultsIntFlow(kFALSE,kTRUE);    
    } 
  
  // e) print on the screen the final results:
  this->PrintFinalResultsForIntegratedFlow();  
  
 } // end of if(doFit)
   
} // end of void AliFlowAnalysisWithFittingQDistribution::Finish(Bool_t doFit)


//================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::GetOutputHistograms(TList *outputListHistos) 
{
 // get pointers to all output histograms (called before Finish()) 
 
 if(outputListHistos)   
 {   
  // 1.) common control histograms and common histograms for final results:
  TString commonHistName = "AliFlowCommonHistFQD";
  commonHistName += fAnalysisLabel->Data();
  AliFlowCommonHist *commonHist = dynamic_cast<AliFlowCommonHist*>(outputListHistos->FindObject(commonHistName.Data()));
  if(commonHist) this->SetCommonHists(commonHist); 
  
  TString commonHistResName = "AliFlowCommonHistResultsFQD";
  commonHistResName += fAnalysisLabel->Data();
  AliFlowCommonHistResults *commonHistRes = dynamic_cast<AliFlowCommonHistResults*>
                                            (outputListHistos->FindObject(commonHistResName.Data()));
  if(commonHistRes) this->SetCommonHistsResults(commonHistRes); 
  
  // 2.) weights: 
  TList *weightsList = dynamic_cast<TList*>(outputListHistos->FindObject("Weights"));
  if(weightsList) this->SetWeightsList(weightsList);
  Bool_t bUsePhiWeights = kFALSE;
  Bool_t bUsePtWeights = kFALSE;
  Bool_t bUseEtaWeights = kFALSE;
  TString fUseParticleWeightsName = "fUseParticleWeightsFQD";
  fUseParticleWeightsName += fAnalysisLabel->Data();
  TProfile *useParticleWeights = dynamic_cast<TProfile*>(weightsList->FindObject(fUseParticleWeightsName.Data()));
  if(useParticleWeights)
  {
   this->SetUseParticleWeights(useParticleWeights);  
   bUsePhiWeights = (Int_t)useParticleWeights->GetBinContent(1);
   bUsePtWeights = (Int_t)useParticleWeights->GetBinContent(2);
   bUseEtaWeights = (Int_t)useParticleWeights->GetBinContent(3);
  }
  
  // 3.) distributions and 4.) final results of fitting:
  TString pWeightsFlag[2] = {"pWeights not used","pWeights used"};
  TString sigmaFlag[2] = {"#sigma^{2} not fitted","#sigma^{2} fitted"};
  
  // q-distribution:
  TString qDistributionName = "fqDistribution";
  qDistributionName += fAnalysisLabel->Data();
  // sum of particle weights:
  TString sumOfParticleWeightsName = "fSumOfParticleWeightsName"; 
  sumOfParticleWeightsName += fAnalysisLabel->Data();
  // final results for integrated flow:
  TString intFlowName = "fIntFlowFQD";
  intFlowName += fAnalysisLabel->Data();
  // sigma^2:
  TString sigma2Name = "fSigma2";
  sigma2Name += fAnalysisLabel->Data();
  // chi^2:
  TString chi2Name = "fChi2";
  chi2Name += fAnalysisLabel->Data();
  
  TH1D *qDistribution[2] = {NULL};
  TH1D *sumOfParticleWeights[2] = {NULL};
  TH1D *intFlow[2][2] = {{NULL}};
  TH1D *sigma2[2][2] = {{NULL}};
  TH1D *chi2[2][2] = {{NULL}};
   
  for(Int_t pW=0;pW<1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights);pW++)
  {
   // q-distribution:
   qDistribution[pW] = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("%s, %s",qDistributionName.Data(),pWeightsFlag[pW].Data())));
   if(qDistribution[pW])
   {
    this->SetqDistribution(qDistribution[pW],pW);
   } else
     {
      cout<<"WARNING: qDistribution[pW] is NULL in AFAWFQD::GOH() !!!!"<<endl;
      cout<<"pW = "<<pW<<endl;
     }
   // sum of particle weights:
   sumOfParticleWeights[pW] = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("%s, %s",sumOfParticleWeightsName.Data(),pWeightsFlag[pW].Data())));
   if(sumOfParticleWeights[pW])
   {
    this->SetSumOfParticleWeights(sumOfParticleWeights[pW],pW);
   } else
     {
      cout<<"WARNING: sumOfParticleWeights[pW] is NULL in AFAWFQD::GOH() !!!!"<<endl;
      cout<<"pW = "<<pW<<endl;
     }
   // final results:
   for(Int_t f=0;f<2;f++)
   {
    // final results for integrated flow:
    intFlow[pW][f] = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("%s, %s, %s",intFlowName.Data(),pWeightsFlag[pW].Data(),sigmaFlag[f].Data())));
    if(intFlow[pW][f])
    {
     this->SetIntFlow(intFlow[pW][f],pW,f);
    } else 
      {
       cout<<"WARNING: intFlow[pW][f] is NULL in AFAWFQD::GOH() !!!!"<<endl;
       cout<<"pW = "<<pW<<endl;
       cout<<"f  = "<<f<<endl;
      }
    // sigma^2:
    sigma2[pW][f] = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("%s, %s, %s",sigma2Name.Data(),pWeightsFlag[pW].Data(),sigmaFlag[f].Data())));
    if(sigma2[pW][f])
    {
     this->SetSigma2(sigma2[pW][f],pW,f);
    } else 
      {
       cout<<"WARNING: sigma2[pW][f] is NULL in AFAWFQD::GOH() !!!!"<<endl;
       cout<<"pW = "<<pW<<endl;
       cout<<"f  = "<<f<<endl;
      } 
    // chi^2:
    chi2[pW][f] = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("%s, %s, %s",chi2Name.Data(),pWeightsFlag[pW].Data(),sigmaFlag[f].Data())));
    if(chi2[pW][f])
    {
     this->SetChi2(chi2[pW][f],pW,f);
    } else 
      {
       cout<<"WARNING: chi2[pW][f] is NULL in AFAWFQD::GOH() !!!!"<<endl;
       cout<<"pW = "<<pW<<endl;
       cout<<"f  = "<<f<<endl;
      }  
      
   } // end of for(Int_t f=0;f<2;f++)
  } // end of for(Int_t pW=0;pW<1+(Int_t)(bUsePhiWeights||bUsePtWeights||bUseEtaWeights);pW++)
  
  // 5.) fitting parameters:
  // q-distribution:
  TString fittingParametersName = "fFittingParameters";
  fittingParametersName += fAnalysisLabel->Data();
  TProfile *fittingParameters = NULL;
  fittingParameters = dynamic_cast<TProfile*>(outputListHistos->FindObject(fittingParametersName.Data()));
  if(fittingParameters)
  {
   this->SetFittingParameters(fittingParameters);
  } else
    {
     cout<<"WARNING:fittingParameters is NULL in AFAWFQD::GOH() !!!!"<<endl;
    }
  
 } else // to if(outputListHistos)
   {
    cout<<"WARNING: outputListHistos is NULL in AFAWFQD::GOH() !!!!"<<endl;
    exit(0);
   } 
    
   
} // end of void AliFlowAnalysisWithFittingQDistribution::GetOutputHistograms(TList *outputListHistos) 


//================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::WriteHistograms(TString* outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName->Data(),"RECREATE");
 //output->WriteObject(fHistList, "cobjFQD","SingleKey");
 fHistList->SetName("cobjFQD");
 fHistList->SetOwner(kTRUE);
 fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
 delete output;
}

//================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::WriteHistograms(TString outputFileName)
{
 //store the final results in output .root file
 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 //output->WriteObject(fHistList, "cobjFQD","SingleKey");
 fHistList->SetName("cobjFQD");
 fHistList->SetOwner(kTRUE);
 fHistList->Write(fHistList->GetName(), TObject::kSingleKey); 
 delete output;
}


//================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::InitializeArrays()
{
 // initialize all arrays
 
 for(Int_t pW=0;pW<2;pW++) // particle weights not used (0) or used (1)
 {
  fSumOfParticleWeights[pW] = NULL;
  fqDistribution[pW] = NULL; 
  for(Int_t f=0;f<2;f++) // sigma^2 not fitted (0) or fitted (1)
  {
   fIntFlow[pW][f] = NULL;
   fSigma2[pW][f] = NULL;
   fChi2[pW][f] = NULL;
  }
 } 

} // end of void AliFlowAnalysisWithFittingQDistribution::InitializeArrays()


//================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::BookCommonHistograms()
{
 // book common histograms
 
 // common control histogram: 
 TString commonHistName = "AliFlowCommonHistFQD";
 commonHistName += fAnalysisLabel->Data();
 fCommonHists = new AliFlowCommonHist(commonHistName.Data());
 fHistList->Add(fCommonHists);  

 // common histograms for final results:
 TString commonHistResName = "AliFlowCommonHistResultsFQD";
 commonHistResName += fAnalysisLabel->Data();
 fCommonHistsResults = new AliFlowCommonHistResults(commonHistResName.Data());
 fHistList->Add(fCommonHistsResults); 

} // end of void AliFlowAnalysisWithFittingQDistribution::BookCommonHistograms()


//================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::BookAndFillWeightsHistograms()
{
 // book and fill histograms which hold phi, pt and eta weights

 if(!fWeightsList)
 {
  cout<<"WARNING: fWeightsList is NULL in AFAWFQD::BAFWH() !!!!"<<endl;
  exit(0);  
 }
    
 TString fUseParticleWeightsName = "fUseParticleWeightsFQD";
 fUseParticleWeightsName += fAnalysisLabel->Data();
 fUseParticleWeights = new TProfile(fUseParticleWeightsName.Data(),"0 = particle weight not used, 1 = particle weight used ",3,0,3);
 fUseParticleWeights->SetLabelSize(0.08);
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
    cout<<"WARNING: fPhiWeights->GetBinWidth(1) != fPhiBinWidth in AFAWFQD::BAFWH() !!!!        "<<endl;
    cout<<"         This indicates inconsistent binning in phi histograms throughout the code."<<endl;
    exit(0);
   }
  } else 
    {
     cout<<"WARNING: fWeightsList->FindObject(\"phi_weights\") is NULL in AFAWFQD::BAFWH() !!!!"<<endl;
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
    cout<<"WARNING: fPtWeights->GetBinWidth(1) != fPtBinWidth in AFAWFQD::BAFWH() !!!!         "<<endl;
    cout<<"         This indicates insconsistent binning in pt histograms throughout the code."<<endl;
    exit(0);
   }
  } else 
    {
     cout<<"WARNING: fWeightsList->FindObject(\"pt_weights\") is NULL in AFAWFQD::BAFWH() !!!!"<<endl;
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
    cout<<"WARNING: fEtaWeights->GetBinWidth(1) != fEtaBinWidth in AFAWFQD::BAFWH() !!!!        "<<endl;
    cout<<"         This indicates insconsistent binning in eta histograms throughout the code."<<endl;
    exit(0);
   }
  } else 
    {
     cout<<"WARNING: fUseEtaWeights && fWeightsList->FindObject(\"eta_weights\") is NULL in AFAWFQD::BAFWH() !!!!"<<endl;
     exit(0);
    }
 } // end of if(fUseEtaWeights)
 
} // end of AliFlowAnalysisWithFittingQDistribution::BookAndFillWeightsHistograms()


//================================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::AccessConstants()
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
 
} // end of void AliFlowAnalysisWithFittingQDistribution::AccessConstants()


//================================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::BookEverythingForDistributions()
{
 // book histograms for distributions
 
 TString pWeightsFlag[2] = {"pWeights not used","pWeights used"};
 TString sigmaFlag[2] = {"#sigma^{2} not fitted","#sigma^{2} fitted"};
 // q-distribution:
 TString fqDistributionName = "fqDistribution";
 fqDistributionName += fAnalysisLabel->Data();
 // sum of particle weights:
 TString fSumOfParticleWeightsName = "fSumOfParticleWeightsName";
 fSumOfParticleWeightsName += fAnalysisLabel->Data();
 // final results for integrated flow:
 TString fIntFlowName = "fIntFlowFQD";
 fIntFlowName += fAnalysisLabel->Data();
 // sigma^2:
 TString fSigma2Name = "fSigma2";
 fSigma2Name += fAnalysisLabel->Data();
 // chi^2:
 TString fChi2Name = "fChi2";
 fChi2Name += fAnalysisLabel->Data();
 
 for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++) // particle weights not used (0) or used (1)
 {
  // q-distribution:
  fqDistribution[pW] = new TH1D(Form("%s, %s",fqDistributionName.Data(),pWeightsFlag[pW].Data()),"q-distribution",10000,0,1000);
  fqDistribution[pW]->SetXTitle("q_{n}=Q_{n}/#sqrt{M}");
  fqDistribution[pW]->SetYTitle("Counts");
  fHistList->Add(fqDistribution[pW]);
  // sum of particle weights: 
  fSumOfParticleWeights[pW] = new TH1D(Form("%s, %s",fSumOfParticleWeightsName.Data(),pWeightsFlag[pW].Data()),"Sum of particle weights",10000,0,10000);
  fSumOfParticleWeights[pW]->SetXTitle("#sum_{i=1}^{N} w_{i}");
  fSumOfParticleWeights[pW]->SetYTitle("Counts");
  fHistList->Add(fSumOfParticleWeights[pW]);
  
  for(Int_t f=0;f<2;f++) // sigma^2 not fitted (0) or fitted (1)
  {
   // final results for integrated flow:
   fIntFlow[pW][f] = new TH1D(Form("%s, %s, %s",fIntFlowName.Data(),pWeightsFlag[pW].Data(),sigmaFlag[f].Data()),"Integrated Flow",1,0,1);
   fIntFlow[pW][f]->SetLabelSize(0.08);
   (fIntFlow[pW][f]->GetXaxis())->SetBinLabel(1,"v_{n}");
   fHistList->Add(fIntFlow[pW][f]);
   // sigma^2:
   fSigma2[pW][f] = new TH1D(Form("%s, %s, %s",fSigma2Name.Data(),pWeightsFlag[pW].Data(),sigmaFlag[f].Data()),"#sigma^{2}",1,0,1);
   fSigma2[pW][f]->SetLabelSize(0.08);
   (fSigma2[pW][f]->GetXaxis())->SetBinLabel(1,"#sigma^{2}");
   fHistList->Add(fSigma2[pW][f]);
   // chi^2:
   fChi2[pW][f] = new TH1D(Form("%s, %s, %s",fChi2Name.Data(),pWeightsFlag[pW].Data(),sigmaFlag[f].Data()),"#chi^{2} (Minuit)",1,0,1);
   fChi2[pW][f]->SetLabelSize(0.08);
   (fChi2[pW][f]->GetXaxis())->SetLabelOffset(0.01);
   (fChi2[pW][f]->GetXaxis())->SetBinLabel(1,"#chi^{2}");
   fHistList->Add(fChi2[pW][f]);
  } // end of for(Int_t f=0;f<2;f++) // sigma^2 not fitted or fitted
  
 } // end of for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++) // particle weights not used (0) or used (1)
 
 // book profile fFittingParameters which will hold all fitting parameters:
 TString fFittingParametersName = "fFittingParameters";
 fFittingParametersName += fAnalysisLabel->Data(); 
 fFittingParameters = new TProfile(fFittingParametersName.Data(),"Parameters for fitting q-distribution",8,0,8);
 fFittingParameters->SetLabelSize(0.05);
 (fFittingParameters->GetXaxis())->SetBinLabel(1,"treshold");
 (fFittingParameters->GetXaxis())->SetBinLabel(2,"starting v_{n}");
 (fFittingParameters->GetXaxis())->SetBinLabel(3,"min. v_{n}");
 (fFittingParameters->GetXaxis())->SetBinLabel(4,"max. v_{n}");
 (fFittingParameters->GetXaxis())->SetBinLabel(5,"starting #sigma^{2}");
 (fFittingParameters->GetXaxis())->SetBinLabel(6,"min. #sigma^{2}");
 (fFittingParameters->GetXaxis())->SetBinLabel(7,"max. #sigma^{2}");
 (fFittingParameters->GetXaxis())->SetBinLabel(8,"plot or not?");
 fHistList->Add(fFittingParameters);
 
} // end of void AliFlowAnalysisWithFittingQDistribution::BookEverythingForDistributions()


//================================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::DoFit(Bool_t useParticleWeights, Bool_t sigma2Fitted)
{
 // do final fit for q-distribution
 
 // shortcuts for flags:
 Int_t pW = (Int_t)(useParticleWeights);
 Int_t s2F = (Int_t)(sigma2Fitted);
 
 if(!(fqDistribution[pW] && fSumOfParticleWeights[pW] && fIntFlow[pW][s2F] && fSigma2[pW][s2F] && fChi2[pW][s2F])) 
 { 
  cout<<"WARNING: fqDistribution[pW] && fSumOfParticleWeights[pW] && fIntFlow[pW][s2F] && fSigma2[pW][s2F] && fChi2[pW][s2F] is NULL in AFAWFQD::DoFit() !!!!"<<endl;
  cout<<"pW  = "<<pW<<endl;
  cout<<"s2F = "<<s2F<<endl;
  exit(0);
 }
 
 // average multiplicity and number of events:
 Double_t AvM = fSumOfParticleWeights[pW]->GetMean(1);
 //Int_t nEvts = (Int_t)fSumOfParticleWeights[pW]->GetEntries();
 
 // start fitting from the bin with at least fTreshold entries, 
 // finish fitting at the bin with at least fTreshold entries:
 Int_t binMin = fqDistribution[pW]->FindFirstBinAbove(fTreshold);  
 Int_t binMax = fqDistribution[pW]->FindLastBinAbove(fTreshold);
 Double_t binWidth = fqDistribution[pW]->GetBinWidth(4); // assuming that all bins have the same width 
 if(binWidth == 0) 
 {
  cout<<"WARNING: binWidth == 0 in AFAWFQD::DoFit()"<<endl;
  exit(0);
 }
 Double_t qmin = (binMin-1)*binWidth; 
 Double_t qmax = (binMax)*binWidth;
 Double_t ent = 0.; // number of entries between binMin and binMax:
 for(Int_t b=binMin;b<=binMax;b++)
 {
  ent += fqDistribution[pW]->GetBinContent(b);
 }
 Double_t norm = binWidth*ent; // norm (assuming that all bins have the same width)

 // fitting function:
 TF1 *fittingFun = new TF1("fittingFun","[2]*(x/[1])*exp(-(x*x+[0]*[0])/(2.*[1]))*TMath::BesselI0(x*[0]/[1])",qmin,qmax); 
 
 fittingFun->SetParNames("v*sqrt{sum of particle weights}","sigma^2","norm");
 fittingFun->SetParameters(fvStart*pow(AvM,0.5),fSigma2Start,norm);         
 fittingFun->SetParLimits(0,fvMin*pow(AvM,0.5),fvMax*pow(AvM,0.5)); 
 
 if(s2F == 0)
 {
  fittingFun->FixParameter(1,0.5);
 } else
   {
    fittingFun->SetParLimits(1,fSigma2Min,fSigma2Max);          
   }
 fittingFun->FixParameter(2,norm);  

 // fitting (do it only if # of entries >50): // to be improved (this is only a pragmatics fix to avoid TMinuit crash)
 if(ent > 50)
 {
  fqDistribution[pW]->Fit("fittingFun","NQ","",qmin,qmax);
 }
 // results:
 Double_t v = 0.; // integrated flow
 Double_t vError = 0.; // error of integrated flow 
 Double_t sigma2 = 0.; // sigma^2
 Double_t sigma2Error = 0.; // error of sigma^2
 Double_t chi2 = 0; // chi^2 from Minuit
 
 if(AvM)
 { 
  // integrated flow:
  v = fittingFun->GetParameter(0)/pow(AvM,0.5);
  vError = fittingFun->GetParError(0)/pow(AvM,0.5);
  fIntFlow[pW][s2F]->SetBinContent(1,v); // s2F is shortcut for "sigma^2 fitted"
  fIntFlow[pW][s2F]->SetBinError(1,vError); // s2F is shortcut for "sigma^2 fitted"
 }
 
 if(s2F == 0) // sigma^2 not fitted, but fixed to 0.5
 {
  // sigma^2:
  sigma2 = 0.5;
  fSigma2[pW][0]->SetBinContent(1,sigma2);  
  fSigma2[pW][0]->SetBinError(1,0.);
  // chi^2:
  chi2 = fittingFun->GetChisquare();
  fChi2[pW][0]->SetBinContent(1,chi2);  
  //fChi2[pW][0]->SetBinError(1,0.);  
 } else // sigma^2 fitted
   {
    // sigma^2:
    sigma2 = fittingFun->GetParameter(1);
    sigma2Error = fittingFun->GetParError(1);
    fSigma2[pW][1]->SetBinContent(1,sigma2);  
    fSigma2[pW][1]->SetBinError(1,sigma2Error);    
    // chi^2:
    chi2 = fittingFun->GetChisquare();
    fChi2[pW][1]->SetBinContent(1,chi2);  
    //fChi2[pW][1]->SetBinError(1,0.);  
   }
 
 if(fPlotResults && !(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)) // to be improved (plot also the plot when particle weights are used)
 {
  // set ranges: // to be improved (there is certainly a better way to implement this)
  Int_t firstNonEmptyBin = fqDistribution[pW]->FindFirstBinAbove(0);
  Double_t lowRange = fqDistribution[pW]->GetBinLowEdge(firstNonEmptyBin);
  Int_t lastNonEmptyBin = fqDistribution[pW]->FindLastBinAbove(0);
  Double_t upperRange = fqDistribution[pW]->GetBinLowEdge(lastNonEmptyBin+10);
  (fqDistribution[pW]->GetXaxis())->SetRangeUser(lowRange,upperRange);
  
  if(s2F == 0)
  { 
   // to be improved (there is certainly a better way to implement this)
   fqDistribution[pW]->SetFillColor(16);  
   fqDistribution[pW]->SetTitle("Fitted q-distribution");
   fqDistribution[pW]->Draw("");
   fLegend->AddEntry(fqDistribution[pW],"q-distribution","f");
   TF1 *fittingFunTemp = (TF1*)fittingFun->Clone("fittingFunTemp");
   fittingFunTemp->SetLineColor(4); // 4 = blue color
   fittingFunTemp->Draw("SAME"); 
   fLegend->AddEntry("fittingFunTemp","#sigma^{2} fixed","l");
   fLegend->Draw("SAME");      
  } else
    {    
     fittingFun->SetLineColor(2); // 2 = red color   
     fittingFun->Draw("SAME");
     fLegend->AddEntry("fittingFun","#sigma^{2} fitted","l");    
    } 
 } // end of if(fPlotResults)

} // end of void AliFlowAnalysisWithFittingQDistribution::DoFit(Bool_t useParticleWeights)


//================================================================================================================================ 


void AliFlowAnalysisWithFittingQDistribution::FillCommonHistResultsIntFlow(Bool_t useParticleWeights, Bool_t sigma2Fitted)
{
 // fill in AliFlowCommonHistResults histograms relevant for 'NONAME' integrated flow (to be improved (name))
 
 // shortcuts for the flags:
 Int_t pW = (Int_t)(useParticleWeights); // particle weights not used (0) or used (1)
 Int_t s2F = (Int_t)(sigma2Fitted); // 0 = sigma^2 not fitted (but fixed to 0.5), 1 = sigma^2 fitted
 
 if(!fIntFlow[pW][s2F])
 {
  cout<<"WARNING: fIntFlow[pW][s2F] is NULL in AFAWFQD::FCHRIF() !!!!"<<endl;
  cout<<"pW  = "<<pW<<endl;
  cout<<"s2F = "<<s2F<<endl;
  exit(0); 
 }  
 
 if(!fSumOfParticleWeights[pW])
 {
  cout<<"WARNING: fSumOfParticleWeights[pW] is NULL in AFAWFQD::FCHRIF() !!!!"<<endl;
  cout<<"pW = "<<pW<<endl;
  exit(0);
 }
 
 if(!(fCommonHistsResults))
 {
  cout<<"WARNING: fCommonHistsResults is NULL in AFAWFQD::FCHRIF() !!!!"<<endl; 
  exit(0);
 }
  
 // fill integrated flow:
 Double_t v = fIntFlow[pW][s2F]->GetBinContent(1); 
 Double_t vError = fIntFlow[pW][s2F]->GetBinError(1);
 
 fCommonHistsResults->FillIntegratedFlow(v,vError);   
 
 // fill chi (this chi stands for resolution, not to be confused with chi2 used before):
 Double_t AvM = fSumOfParticleWeights[pW]->GetMean(1);
 Double_t chi = AvM*pow(v,2); 
 if(chi>=0)
 {
  fCommonHistsResults->FillChi(pow(chi,0.5));   
  fCommonHistsResults->FillChiRP(pow(chi,0.5));   
 }
   
} // end of void AliFlowAnalysisWithFittingQDistribution::FillCommonHistResultsIntFlow(Bool_t useParticleWeights, Bool_t sigma2NotFixed) 


//================================================================================================================================ 


void AliFlowAnalysisWithFittingQDistribution::PrintFinalResultsForIntegratedFlow()
{
 // print the final results for integrated flow on the screen
 
 // shortcuts: pW  = particle weights 
 //            s2F = sigma^2 fitted 
 
 for(Int_t pW=0;pW<1+(Int_t)(fUsePhiWeights||fUsePtWeights||fUseEtaWeights);pW++)
 {
  if(!fSumOfParticleWeights[pW])
  {
   cout<<"WARNING: fSumOfParticleWeights[pW] is NULL in AFAWFQD::FCHRIF() !!!!"<<endl;
   cout<<"pW = "<<pW<<endl;
   exit(0);
  }
  for(Int_t s2F=0;s2F<2;s2F++)
  {
   if(!fIntFlow[pW][s2F])
   {
    cout<<"WARNING: fIntFlow[pW][s2F] is NULL in AFAWFQD::FCHRIF() !!!!"<<endl;
    cout<<"pW  = "<<pW<<endl;
    cout<<"s2F = "<<s2F<<endl;
    exit(0); 
   }
  }  
 }  
 
 if(!(fCommonHistsResults))
 {
  cout<<"WARNING: fCommonHistsResults is NULL in AFAWFQD::FCHRIF() !!!!"<<endl; 
  exit(0);
 }
 
 // shortcut for the harmonic:
 Int_t n = (Int_t)(fCommonHists->GetHarmonic())->GetBinContent(1); 

 // printing:
 cout<<" "<<endl;
 cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl;
 cout<<"      integrated flow by fitting "<<endl;
 cout<<"           q-distribution:      "<<endl;
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
 {
  cout<<"           (with weights)       "<<endl;
 } else
   {
    cout<<"          (without weights)       "<<endl;
   }   
 cout<<endl;

 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
 {
  cout<<"1.) sigma^2 not fitted: "<<endl;
  cout<<"  v_"<<n<<"{FQD} = "<<fIntFlow[1][0]->GetBinContent(1)<<" +/- "<<fIntFlow[1][0]->GetBinError(1)<<endl;
  cout<<"  sigma^2 = 0.5 +/- 0 "<<endl; 
  cout<<"  chi^2 = "<<fChi2[1][0]->GetBinContent(1)<<" (Minuit)"<<endl; 
  cout<<" "<<endl;   
  cout<<"2.) sigma^2 fitted: "<<endl;
  cout<<"  v_"<<n<<"{FQD} = "<<fIntFlow[1][1]->GetBinContent(1)<<" +/- "<<fIntFlow[1][1]->GetBinError(1)<<endl;
  cout<<"  sigma^2 = "<<fSigma2[1][1]->GetBinContent(1)<<" +/- "<<fSigma2[1][1]->GetBinError(1)<<endl; 
  cout<<"  chi^2 = "<<fChi2[1][1]->GetBinContent(1)<<" (Minuit)"<<endl; 
  cout<<" "<<endl; 
  cout<<"      nEvts = "<<fSumOfParticleWeights[1]->GetEntries()<<", AvM = "<<fSumOfParticleWeights[1]->GetMean()<<endl;
  cout<<" "<<endl;
 } else
   { 
    cout<<"1.) sigma^2 not fitted: "<<endl;
    cout<<endl;
    cout<<"  v_"<<n<<"{FQD} = "<<fIntFlow[0][0]->GetBinContent(1)<<" +/- "<<fIntFlow[0][0]->GetBinError(1)<<endl;
    cout<<"  sigma^2 = 0.5 +/- 0 "<<endl; 
    cout<<"  chi^2 = "<<fChi2[0][0]->GetBinContent(1)<<" (Minuit)"<<endl;  
    cout<<" "<<endl;   
    cout<<"2.) sigma^2 fitted: "<<endl;
    cout<<endl;
    cout<<"  v_"<<n<<"{FQD} = "<<fIntFlow[0][1]->GetBinContent(1)<<" +/- "<<fIntFlow[0][1]->GetBinError(1)<<endl;
    cout<<"  sigma^2 = "<<fSigma2[0][1]->GetBinContent(1)<<" +/- "<<fSigma2[0][1]->GetBinError(1)<<endl; 
    cout<<"  chi^2 = "<<fChi2[0][1]->GetBinContent(1)<<" (Minuit)"<<endl; 
    cout<<" "<<endl;  
    cout<<"      nEvts = "<<fSumOfParticleWeights[0]->GetEntries()<<", AvM = "<<fSumOfParticleWeights[0]->GetMean()<<endl;
    cout<<" "<<endl;
   }
    
 cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl; 
 cout<<endl;
  
} // end of void AliFlowAnalysisWithFittingQDistribution::PrintFinalResultsForIntegratedFlow()


//================================================================================================================================ 


void AliFlowAnalysisWithFittingQDistribution::StoreFittingParameters()
{
 // store fitting parameters in profile fFittingParameters
 
 // Binning of fFittingParameters is organized as follows:
 // 1st bin: fTreshold
 // 2nd bin: fvStart
 // 3rd bin: fvMin
 // 4th bin: fvMax
 // 5th bin: fSigma2Start
 // 6th bin: fSigma2Min
 // 7th bin: fSigma2Max
 // 8th bin: fPlotResults
 
 if(!fFittingParameters)
 {
  cout<<"WARNING: fFittingParameters is NULL in AFAWFQD::SFP() !!!!"<<endl;
  exit(0);
 }
 
 fFittingParameters->Reset();
 fFittingParameters->Fill(0.5,fTreshold);
 fFittingParameters->Fill(1.5,fvStart);
 fFittingParameters->Fill(2.5,fvMin);
 fFittingParameters->Fill(3.5,fvMax);
 fFittingParameters->Fill(4.5,fSigma2Start);
 fFittingParameters->Fill(5.5,fSigma2Min);
 fFittingParameters->Fill(6.5,fSigma2Max);
 fFittingParameters->Fill(7.5,(Int_t)fPlotResults);
 
} // end of void AliFlowAnalysisWithFittingQDistribution::StoreFittingParameters()


//================================================================================================================================ 


void AliFlowAnalysisWithFittingQDistribution::AccessFittingParameters()
{
 // access fitting parameters:
 
 if(!fFittingParameters)
 {
  cout<<"WARNING: fFittingParameters is NULL in AFAWFQD::AFP() !!!!"<<endl;
  exit(0);
 }
 
 fTreshold = fFittingParameters->GetBinContent(1);
 fvStart = fFittingParameters->GetBinContent(2);
 fvMin = fFittingParameters->GetBinContent(3);
 fvMax = fFittingParameters->GetBinContent(4);
 fSigma2Start = fFittingParameters->GetBinContent(5);
 fSigma2Min = fFittingParameters->GetBinContent(6);
 fSigma2Max = fFittingParameters->GetBinContent(7);
 fPlotResults = (Bool_t) fFittingParameters->GetBinContent(8);
 
} // end of void AliFlowAnalysisWithFittingQDistribution::AccessFittingParameters()
