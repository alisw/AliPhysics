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
 * estimating reference flow by *
 *   fitting q-distribution     * 
 *                              *
 * author: Ante Bilandzic       * 
 *       (abilandzic@gmail.com) *
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
 fSumOfParticleWeights(NULL),
 fqDistribution(NULL),
 fqMin(0.),
 fqMax(1000.),
 fqNbins(10000),
 fFittingParameters(NULL), 
 fTreshold(5.),
 fvStart(0.05),
 fvMin(0.0),
 fvMax(0.25),
 fSigma2Start(0.75),
 fSigma2Min(0.5), 
 fSigma2Max(2.5),
 fFinalResultIsFromSigma2Fitted(kTRUE),
 fPrintOnTheScreen(kTRUE)
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
 // destructor
 
 delete fHistList; 
 
} // end of destructor

//================================================================================================================


void AliFlowAnalysisWithFittingQDistribution::Init()
{
 // Access constants and book everything.
 
 //save old value and prevent histograms from being added to directory
 //to avoid name clashes in case multiple analaysis objects are used
 //in an analysis
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
 TH1::AddDirectory(kFALSE);
 
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
 
 //restore old status
 TH1::AddDirectory(oldHistAddStatus);
} // end of void AliFlowAnalysisWithFittingQDistribution::Init()

//================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::Make(AliFlowEventSimple* anEvent)
{
 // Loop over data only in this method.
 
 // a) Check all pointers used in this method;
 // b) Fill the common control histograms;
 // c) Loop over data and calculate Q-vector and sum of particle weights;
 // d) Fill the histogram for q-distribution and sum of particle weights.
 
 // a) Check all pointers used in this method:
 this->CheckPointersUsedInMake();
 
 // b) fill the common control histograms:
 fCommonHists->FillControlHistograms(anEvent);                                         
 
 // c) Loop over data and fill histogram for q-distribution:                                          
 Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
 Double_t dPt  = 0.; // transverse momentum
 Double_t dEta = 0.; // pseudorapidity
 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight 
 Double_t dReQ = 0.; // real part of Q-vector 
 Double_t dImQ = 0.; // imaginary part of Q-vector
 Int_t n = fHarmonic; // shortcut for the harmonic 
 Double_t dSumOfParticleWeights = 0.; // when particle weights are not used dSumOfParticleWeights is equal to multiplicity
 AliFlowTrackSimple *aftsTrack = NULL;  
 Int_t nPrim = anEvent->NumberOfTracks(); // nPrim = total number of primary tracks, i.e. nPrim = nRP + nPOI + rest, where:
                                          // nRP   = # of particles used to determine the reaction plane;
                                          // nPOI  = # of particles of interest for a detailed flow analysis;
                                          // rest  = # of particles which are not niether RPs nor POIs.   
 // Start loop over particles:
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
   if(fUsePtWeights && fPtWeights && fPtBinWidth) // determine pt weight for this particle:
   {
    wPt = fPtWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/fPtBinWidth))); 
   }            
   if(fUseEtaWeights && fEtaWeights && fEtaBinWidth) // determine eta weight for this particle: 
   {
    wEta = fEtaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/fEtaBinWidth))); 
   }   
   // Calculate real and imaginary part of Q-vector and sum of particle weights for this event:
   // Q-vector:
   dReQ += wPhi*wPt*wEta*TMath::Cos(n*dPhi); 
   dImQ += wPhi*wPt*wEta*TMath::Sin(n*dPhi);
   // sum of particle weights:
   dSumOfParticleWeights += wPhi*wPt*wEta; // when particle weights are not used this sum gives # of RPs, i.e. multiplicity   
  } // end of if(aftsTrack)
 } // end of for(Int_t i=0;i<nPrim;i++)                                      
                                           
 // d) Fill the histogram for q-distribution and sum of particle weights:
 Double_t q = 0.; // q = Q\sqrt{sum of particle weights}                                         
 if(dSumOfParticleWeights)
 {
  q = pow(dReQ*dReQ+dImQ*dImQ,0.5)/pow(dSumOfParticleWeights,0.5);
  fqDistribution->Fill(q,1.);
  fSumOfParticleWeights->Fill(dSumOfParticleWeights,1.);
 } else
   {
    cout<<endl;
    cout<<" WARNING (FQD::Make()): dSumOfParticleWeights == 0. !!!!"<<endl;
    cout<<endl;
   }

} // end of Make()

//================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::Finish(Bool_t doFit)
{
 // Calculate the final results.
 
 // a) Check all pointers used in this method;
 // b) Acces common constants;
 // c) Access fitting paremeters;
 // d) Do the final fit of q-distribution;
 // e) Fill common hist results;
 // f) Print on the screen the final results.
 
 // a) Check all pointers used in this method:
 this->CheckPointersUsedInFinish(); 
 
 // b) Access common constants:
 this->AccessConstants();
 
 // c) Access fitting paremeters:
 this->AccessFittingParameters();
 
 // d) Do the final fit of q-distribution:             
 if(doFit) 
 {
  this->DoFit(kFALSE); // sigma^2 not fitted (fixed to 0.5)
  this->DoFit(kTRUE); // sigma^2 fitted
 } 
   
 // e) Fill common hist results (by default fill with results obtained with sigma^2 fitted,
 //    alternatively use a setter SetFinalResultIsFromSigma2Fitted(kFALSE)):
 this->FillCommonHistResults(fFinalResultIsFromSigma2Fitted); 
  
 // f) Print on the screen the final results:
 if(fPrintOnTheScreen) this->PrintOnTheScreen();  
  
} // end of void AliFlowAnalysisWithFittingQDistribution::Finish(Bool_t doFit)
 
//================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::CheckPointersUsedInMake()
{
 // Check all pointers used in method Make().
 
 if(!fCommonHists)
 {
  cout<<endl;
  cout<<" WARNING (FQD::Make()): fCommonHists is NULL !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!fqDistribution)
 {
  cout<<endl;
  cout<<" WARNING (FQD::Make()): fqDistribution is NULL !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!fSumOfParticleWeights)
 {
  cout<<endl;
  cout<<" WARNING (FQD::Make()): fSumOfParticleWeights is NULL !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(fUsePhiWeights && !fPhiWeights)
 {
  cout<<endl;
  cout<<" WARNING (FQD::Make()): fPhiWeights is NULL !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(fUsePtWeights && !fPtWeights)
 {
  cout<<endl;
  cout<<" WARNING (FQD::Make()): fPtWeights is NULL !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(fUseEtaWeights && !fEtaWeights)
 {
  cout<<endl;
  cout<<" WARNING (FQD::Make()): fEtaWeights is NULL !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 

} // end of void AliFlowAnalysisWithFittingQDistribution::CheckPointersUsedInMake()

//================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::CheckPointersUsedInFinish()
{
 // Check all pointers used in method Finish().
 
 if(!fFittingParameters)
 {
  cout<<endl;
  cout<<" WARNING (FQD::Finish()): fFittingParameters is NULL !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!fqDistribution)
 {
  cout<<endl;
  cout<<" WARNING (FQD::Finish()): fqDistribution is NULL !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!fSumOfParticleWeights)
 {
  cout<<endl;
  cout<<" WARNING (FQD::Finish()): fSumOfParticleWeights is NULL !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 for(Int_t s2F=0;s2F<2;s2F++)
 {
  if(!fIntFlow[s2F])
  {
   cout<<endl;
   cout<<" WARNING (FQD::Finish()): "<<Form("fIntFlow[%d] is NULL !!!!",s2F)<<endl;
   cout<<endl;
   exit(0); 
  }
  if(!fSigma2[s2F])
  {
   cout<<endl;
   cout<<" WARNING (FQD::Finish()): "<<Form("fSigma2[%d] is NULL !!!!",s2F)<<endl;
   cout<<endl;
   exit(0); 
  }
  if(!fChi2[s2F])
  {
   cout<<endl;
   cout<<" WARNING (FQD::Finish()): "<<Form("fChi2[%d] is NULL !!!!",s2F)<<endl;
   cout<<endl;
   exit(0); 
  }
 } // end of for(Int_t s2F=0;s2F<2;s2F++)
 if(!fCommonHistsResults)
 {
  cout<<endl;
  cout<<"WARNING (FQD::Finish()): fCommonHistsResults is NULL !!!!"<<endl; 
  cout<<endl;
  exit(0);
 }
 if(!fCommonHists)
 {
  cout<<endl;
  cout<<"WARNING (FQD::Finish()): fCommonHists is NULL !!!!"<<endl; 
  cout<<endl;
  exit(0);
 }
  
} // end of void AliFlowAnalysisWithFittingQDistribution::CheckPointersUsedInFinish()
 
//================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::GetOutputHistograms(TList *outputListHistos) 
{
 // Get pointers to all output histograms (called before Finish()). 
 
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
  if(weightsList){this->SetWeightsList(weightsList);}
  Bool_t bUsePhiWeights = kFALSE;
  Bool_t bUsePtWeights = kFALSE;
  Bool_t bUseEtaWeights = kFALSE;
  TString fUseParticleWeightsName = "fUseParticleWeightsFQD";
  fUseParticleWeightsName += fAnalysisLabel->Data();
  TProfile *useParticleWeights = NULL; 
  if(weightsList)
  {
   useParticleWeights = dynamic_cast<TProfile*>(weightsList->FindObject(fUseParticleWeightsName.Data()));
  }
  if(useParticleWeights)
  {
   this->SetUseParticleWeights(useParticleWeights);  
   bUsePhiWeights = (Int_t)useParticleWeights->GetBinContent(1);
   this->SetUsePhiWeights(bUsePhiWeights);
   bUsePtWeights = (Int_t)useParticleWeights->GetBinContent(2);
   this->SetUsePtWeights(bUsePtWeights);
   bUseEtaWeights = (Int_t)useParticleWeights->GetBinContent(3);
   this->SetUseEtaWeights(bUseEtaWeights);
  }
  
  // 3.) distributions and 4.) final results of fitting:
  TString sigmaFlag[2] = {"#sigma^{2} not fitted","#sigma^{2} fitted"};  
  // q-distribution:
  TString qDistributionName = "fqDistribution";
  qDistributionName += fAnalysisLabel->Data();
  // sum of particle weights:
  TString sumOfParticleWeightsName = "fSumOfParticleWeights"; 
  sumOfParticleWeightsName += fAnalysisLabel->Data();
  // final results for reference flow:
  TString intFlowName = "fIntFlowFQD";
  intFlowName += fAnalysisLabel->Data();
  // sigma^2:
  TString sigma2Name = "fSigma2";
  sigma2Name += fAnalysisLabel->Data();
  // chi^2:
  TString chi2Name = "fChi2";
  chi2Name += fAnalysisLabel->Data();
  // fitting function:
  TString fittingFunctionName = "fFittingFunction";
  fittingFunctionName += fAnalysisLabel->Data();
  
  TH1D *qDistribution = NULL;
  TH1D *sumOfParticleWeights = NULL;
  TH1D *intFlow[2] = {NULL};
  TH1D *sigma2[2] = {NULL};
  TH1D *chi2[2] = {NULL};
  TF1 *fittingFunction[2] = {NULL};
   
  // q-distribution:
  qDistribution = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("%s",qDistributionName.Data())));
  if(qDistribution)
  {
   this->SetqDistribution(qDistribution);
  } else
    {
     cout<<"WARNING: qDistribution is NULL in AFAWFQD::GOH() !!!!"<<endl;
    }
  // sum of particle weights:
  sumOfParticleWeights = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("%s",sumOfParticleWeightsName.Data())));
  if(sumOfParticleWeights)
  {
   this->SetSumOfParticleWeights(sumOfParticleWeights);
  } else
    {
     cout<<"WARNING: sumOfParticleWeights is NULL in AFAWFQD::GOH() !!!!"<<endl;
    }
  // final results:
  for(Int_t f=0;f<2;f++)
  {
   // final results for reference flow:
   intFlow[f] = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("%s, %s",intFlowName.Data(),sigmaFlag[f].Data())));
   if(intFlow[f])
   {
    this->SetIntFlow(intFlow[f],f);
   } else 
     {
      cout<<"WARNING: intFlow[f] is NULL in AFAWFQD::GOH() !!!!"<<endl;
      cout<<"f  = "<<f<<endl;
     }
   // sigma^2:
   sigma2[f] = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("%s, %s",sigma2Name.Data(),sigmaFlag[f].Data())));
   if(sigma2[f])
   {
    this->SetSigma2(sigma2[f],f);
   } else 
     { 
      cout<<"WARNING: sigma2[f] is NULL in AFAWFQD::GOH() !!!!"<<endl;
      cout<<"f  = "<<f<<endl;
     } 
   // chi^2:
   chi2[f] = dynamic_cast<TH1D*>(outputListHistos->FindObject(Form("%s, %s",chi2Name.Data(),sigmaFlag[f].Data())));
   if(chi2[f])
   {
    this->SetChi2(chi2[f],f);
   } else 
     {
      cout<<"WARNING: chi2[f] is NULL in AFAWFQD::GOH() !!!!"<<endl;
      cout<<"f  = "<<f<<endl;
     }   
   // fitting functions:
   fittingFunction[f] = dynamic_cast<TF1*>(outputListHistos->FindObject(Form("%s, %s",fittingFunctionName.Data(),sigmaFlag[f].Data())));
   if(fittingFunction[f])
   {
    this->SetFittingFunction(fittingFunction[f],f);
   } else 
     {
      cout<<"WARNING: fittingFunction[f] is NULL in AFAWFQD::GOH() !!!!"<<endl;
      cout<<"f  = "<<f<<endl;
     }       
  } // end of for(Int_t f=0;f<2;f++)
 
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

void AliFlowAnalysisWithFittingQDistribution::WriteHistograms(TDirectoryFile *outputFileName)
{
 //store the final results in output .root file
 fHistList->SetName("cobjFQD");
 fHistList->SetOwner(kTRUE);
 outputFileName->Add(fHistList);
 outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);
}

//================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::InitializeArrays()
{
 // Initialize all arrays.
 
 for(Int_t s2F=0;s2F<2;s2F++) // sigma^2 not fitted (0) or fitted (1)
 {
  fIntFlow[s2F] = NULL;
  fSigma2[s2F] = NULL;
  fChi2[s2F] = NULL;
  fFittingFunction[s2F] = NULL;
 } 

} // end of void AliFlowAnalysisWithFittingQDistribution::InitializeArrays()

//================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::BookCommonHistograms()
{
 // Book common histograms.
 
 // common control histogram: 
 TString commonHistName = "AliFlowCommonHistFQD";
 commonHistName += fAnalysisLabel->Data();
 fCommonHists = new AliFlowCommonHist(commonHistName.Data());
 fHistList->Add(fCommonHists);  

 // common histograms for final results:
 TString commonHistResName = "AliFlowCommonHistResultsFQD";
 commonHistResName += fAnalysisLabel->Data();
 fCommonHistsResults = new AliFlowCommonHistResults(commonHistResName.Data(),"",fHarmonic);
 fHistList->Add(fCommonHistsResults); 

} // end of void AliFlowAnalysisWithFittingQDistribution::BookCommonHistograms()

//================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::BookAndFillWeightsHistograms()
{
 // Book and fill histograms which hold phi, pt and eta weights.

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
 // phi weights:
 if(fUsePhiWeights)
 {
  if(fWeightsList->FindObject("phi_weights"))
  {
   fPhiWeights = dynamic_cast<TH1F*>(fWeightsList->FindObject("phi_weights"));
   if(!fPhiWeights){printf("\n WARNING (FQD): !fPhiWeights !!!!\n");exit(0);}   
   if(TMath::Abs(fPhiWeights->GetBinWidth(1)-fPhiBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (FQD): Inconsistent binning in histograms for phi-weights throughout the code."<<endl;
    cout<<endl;
    //exit(0);
   }
  } else 
    {
     cout<<"WARNING: fWeightsList->FindObject(\"phi_weights\") is NULL in AFAWFQD::BAFWH() !!!!"<<endl;
     exit(0);
    }
 } // end of if(fUsePhiWeights)
 // pt weights:
 if(fUsePtWeights) 
 {
  if(fWeightsList->FindObject("pt_weights"))
  {
   fPtWeights = dynamic_cast<TH1D*>(fWeightsList->FindObject("pt_weights"));
   if(!fPtWeights){printf("\n WARNING (FQD): !fPtWeights !!!!\n");exit(0);}   
   if(TMath::Abs(fPtWeights->GetBinWidth(1)-fPtBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (FQD): Inconsistent binning in histograms for pt-weights throughout the code."<<endl;
    cout<<endl;
    //exit(0);
   }
  } else 
    {
     cout<<"WARNING: fWeightsList->FindObject(\"pt_weights\") is NULL in AFAWFQD::BAFWH() !!!!"<<endl;
     exit(0);
    }
 } // end of if(fUsePtWeights)    
 // eta weights:
 if(fUseEtaWeights) 
 {
  if(fWeightsList->FindObject("eta_weights"))
  {
   fEtaWeights = dynamic_cast<TH1D*>(fWeightsList->FindObject("eta_weights"));
   if(!fEtaWeights){printf("\n WARNING (FQD): !fEtaWeights !!!!\n");exit(0);}   
   if(TMath::Abs(fEtaWeights->GetBinWidth(1)-fEtaBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (FQD): Inconsistent binning in histograms for eta-weights throughout the code."<<endl;
    cout<<endl;
    //exit(0);
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
 
} // end of void AliFlowAnalysisWithFittingQDistribution::AccessConstants()

//================================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::BookEverythingForDistributions()
{
 // Book all objects relevant for fitting of q-distributions.
 
 // Flags:
 TString sigmaFlag[2] = {"#sigma^{2} not fitted","#sigma^{2} fitted"};
 // q-distribution:
 TString fqDistributionName = "fqDistribution";
 fqDistributionName += fAnalysisLabel->Data();
 fqDistribution = new TH1D(Form("%s",fqDistributionName.Data()),"q-distribution",fqNbins,fqMin,fqMax);  
 fqDistribution->SetXTitle("q_{n}=|Q_{n}|/#sqrt{M}");
 fqDistribution->SetYTitle("Counts");
 fHistList->Add(fqDistribution);
 // Sum of particle weights: 
 TString fSumOfParticleWeightsName = "fSumOfParticleWeights";
 fSumOfParticleWeightsName += fAnalysisLabel->Data();
 fSumOfParticleWeights = new TH1D(Form("%s",fSumOfParticleWeightsName.Data()),"Sum of particle weights",10000,0,100000); // (to be improved - harwired limits and number of bins)
 fSumOfParticleWeights->SetXTitle("#sum_{i=1}^{N} w_{i}");
 fSumOfParticleWeights->GetXaxis()->SetTitleSize(0.03);
 fSumOfParticleWeights->SetYTitle("Counts");
 fHistList->Add(fSumOfParticleWeights); 
 // Final results for integrated flow:
 TString fIntFlowName = "fIntFlowFQD";
 fIntFlowName += fAnalysisLabel->Data();
 // sigma^2:
 TString fSigma2Name = "fSigma2";
 fSigma2Name += fAnalysisLabel->Data();
 // chi^2:
 TString fChi2Name = "fChi2";
 fChi2Name += fAnalysisLabel->Data();
 // fitting function: 
 TString fittingFunctionName = "fFittingFunction";
 fittingFunctionName += fAnalysisLabel->Data();
 
 for(Int_t f=0;f<2;f++) // sigma^2 not fitted (0) or fitted (1)
 {
  // final results for integrated flow:
  fIntFlow[f] = new TH1D(Form("%s, %s",fIntFlowName.Data(),sigmaFlag[f].Data()),Form("Reference Flow (%s)",sigmaFlag[f].Data()),1,0,1);
  fIntFlow[f]->SetLabelSize(0.08);
  (fIntFlow[f]->GetXaxis())->SetBinLabel(1,"v_{n}");
  fHistList->Add(fIntFlow[f]);
  // sigma^2:
  fSigma2[f] = new TH1D(Form("%s, %s",fSigma2Name.Data(),sigmaFlag[f].Data()),Form("#sigma^{2} (%s)",sigmaFlag[f].Data()),1,0,1);
  fSigma2[f]->SetLabelSize(0.08);
  (fSigma2[f]->GetXaxis())->SetBinLabel(1,"#sigma^{2}");
  fHistList->Add(fSigma2[f]);
  // chi^2:
  fChi2[f] = new TH1D(Form("%s, %s",fChi2Name.Data(),sigmaFlag[f].Data()),Form("#chi^{2} (%s)",sigmaFlag[f].Data()),1,0,1);
  fChi2[f]->SetLabelSize(0.08);
  (fChi2[f]->GetXaxis())->SetLabelOffset(0.01);
  (fChi2[f]->GetXaxis())->SetBinLabel(1,"#chi^{2}");
  fHistList->Add(fChi2[f]);
  // fitting functions:
  fFittingFunction[f] = new TF1(Form("%s, %s",fittingFunctionName.Data(),sigmaFlag[f].Data()),"[2]*(x/[1])*exp(-(x*x+[0]*[0])/(2.*[1]))*TMath::BesselI0(x*[0]/[1])");
  fHistList->Add(fFittingFunction[f]);
 } // end of for(Int_t f=0;f<2;f++) // sigma^2 not fitted or fitted 
 // Book profile fFittingParameters which will hold all fitting parameters:
 TString fFittingParametersName = "fFittingParameters";
 fFittingParametersName += fAnalysisLabel->Data(); 
 fFittingParameters = new TProfile(fFittingParametersName.Data(),"Parameters for fitting q-distribution",9,0,9);
 fFittingParameters->SetLabelSize(0.05);
 (fFittingParameters->GetXaxis())->SetBinLabel(1,"treshold");
 (fFittingParameters->GetXaxis())->SetBinLabel(2,"starting v_{n}");
 (fFittingParameters->GetXaxis())->SetBinLabel(3,"min. v_{n}");
 (fFittingParameters->GetXaxis())->SetBinLabel(4,"max. v_{n}");
 (fFittingParameters->GetXaxis())->SetBinLabel(5,"starting #sigma^{2}");
 (fFittingParameters->GetXaxis())->SetBinLabel(6,"min. #sigma^{2}");
 (fFittingParameters->GetXaxis())->SetBinLabel(7,"max. #sigma^{2}");
 (fFittingParameters->GetXaxis())->SetBinLabel(8,"Final result is from #sigma^{2} fitted?"); 
 (fFittingParameters->GetXaxis())->SetBinLabel(9,"Print results on the screen?"); 
 fHistList->Add(fFittingParameters);
 
} // end of void AliFlowAnalysisWithFittingQDistribution::BookEverythingForDistributions()

//================================================================================================================================

void AliFlowAnalysisWithFittingQDistribution::DoFit(Bool_t sigma2Fitted)
{
 // Do the final fit of q-distribution.
 
 Int_t s2F = (Int_t)(sigma2Fitted); // shortcut
 Double_t AvM = fSumOfParticleWeights->GetMean(1); // average multiplicity
 //Int_t nEvts = (Int_t)fSumOfParticleWeights->GetEntries(); // number of events:
 
 // Start fitting from the bin with at least fTreshold entries, 
 // finish fitting at the bin with at least fTreshold entries:
 Int_t binMin = fqDistribution->FindFirstBinAbove(fTreshold);  
 Int_t binMax = fqDistribution->FindLastBinAbove(fTreshold);
 Double_t binWidth = fqDistribution->GetBinWidth(4); // assuming that all bins have the same width 
 if(TMath::Abs(binWidth) < 1.e-44) 
 {
  cout<<endl;
  cout<<"WARNING (FQD): binWidth == 0 in AFAWFQD::DoFit()"<<endl;
  cout<<endl;
  exit(0);
 }
 Double_t qmin = (binMin-1)*binWidth; 
 Double_t qmax = (binMax)*binWidth;
 Double_t ent = 0.; // number of entries between binMin and binMax:
 for(Int_t b=binMin;b<=binMax;b++)
 {
  ent += fqDistribution->GetBinContent(b);
 }
 Double_t norm = binWidth*ent; // norm (assuming that all bins have the same width)
 // Fitting function:
 fFittingFunction[s2F]->SetRange(qmin,qmax); 
 fFittingFunction[s2F]->SetParNames("v*sqrt{sum of particle weights}","sigma^2","norm");
 fFittingFunction[s2F]->SetParameters(fvStart*pow(AvM,0.5),fSigma2Start,norm);         
 fFittingFunction[s2F]->SetParLimits(0,fvMin*pow(AvM,0.5),fvMax*pow(AvM,0.5)); 
 if(s2F == 0)
 {
  fFittingFunction[s2F]->FixParameter(1,0.5);
 } else
   {
    fFittingFunction[s2F]->SetParLimits(1,fSigma2Min,fSigma2Max);          
   }
 fFittingFunction[s2F]->FixParameter(2,norm);  
 // Fitting is done here:
 // Remark: do fit only if # of entries > 50 - this is only a pragmatics fix to avoid TMinuit crash (to be improved)
 if(ent > 50)
 {
  fqDistribution->Fit(fFittingFunction[s2F]->GetName(),"NQ","",qmin,qmax);
 }
 // Final results:
 Double_t v = 0.; // reference flow
 Double_t vError = 0.; // error of reference flow 
 Double_t sigma2 = 0.; // sigma^2
 Double_t sigma2Error = 0.; // error of sigma^2
 Double_t chi2 = 0; // chi^2 from Minuit
 // Reference flow:
 if(AvM)
 { 
  v = fFittingFunction[s2F]->GetParameter(0)/pow(AvM,0.5);
  vError = fFittingFunction[s2F]->GetParError(0)/pow(AvM,0.5);
  fIntFlow[s2F]->SetBinContent(1,v); // s2F is shortcut for "sigma^2 fitted"
  fIntFlow[s2F]->SetBinError(1,vError); // s2F is shortcut for "sigma^2 fitted"
 } else
   {
    cout<<endl;
    cout<<"WARNING (FQD): AvM == 0 in AFAWFQD::DoFit()"<<endl;
    cout<<endl;
   }    
 // sigma^2:: 
 if(s2F == 0) // sigma^2 not fitted, but fixed to 0.5
 {
  sigma2 = 0.5;
  fSigma2[0]->SetBinContent(1,sigma2);  
  fSigma2[0]->SetBinError(1,0.);
 } else // sigma^2 fitted
   {
    sigma2 = fFittingFunction[s2F]->GetParameter(1);
    sigma2Error = fFittingFunction[s2F]->GetParError(1);
    fSigma2[1]->SetBinContent(1,sigma2);  
    fSigma2[1]->SetBinError(1,sigma2Error);    
   }
 // chi^2:
 chi2 = fFittingFunction[s2F]->GetChisquare();
 fChi2[s2F]->SetBinContent(1,chi2);     
 
} // end of void AliFlowAnalysisWithFittingQDistribution::DoFit()

//================================================================================================================================ 

void AliFlowAnalysisWithFittingQDistribution::FillCommonHistResults(Bool_t sigma2Fitted)
{
 // Fill common result histogram for reference flow and resolution.
  
 // Remark: by default the result obtained with sigma^2 fitted is being stored. 
 // In order to store the result obtained with sigma^2 fixed use a setter SetFinalResultIsFromSigma2Fitted(kFALSE).
  
 Int_t s2F = (Int_t)sigma2Fitted; // shortcut
  
 // Reference flow:
 Double_t v = fIntFlow[s2F]->GetBinContent(1); 
 Double_t vError = fIntFlow[s2F]->GetBinError(1);
 fCommonHistsResults->FillIntegratedFlow(v,vError);   
 // Resolution:
 Double_t AvM = fSumOfParticleWeights->GetMean(1);
 Double_t chi2 = AvM*pow(v,2.); // chi^2
 if(chi2>=0.)
 {
  fCommonHistsResults->FillChi(pow(chi2,0.5));   
 }
   
} // end of void AliFlowAnalysisWithFittingQDistribution::FillCommonHistResultsIntFlow(Bool_t useParticleWeights, Bool_t sigma2NotFixed) 

//================================================================================================================================ 

void AliFlowAnalysisWithFittingQDistribution::PrintOnTheScreen()
{
 // Print the final results on the screen.
 
 // shortcut for the harmonic:
 Int_t n = 0;
 if(fCommonHists->GetHarmonic())
 {
  n = (Int_t)(fCommonHists->GetHarmonic())->GetBinContent(1); 
 }
 
 // printing:
 cout<<endl;
 cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl;
 cout<<"      reference flow by fitting "<<endl;
 cout<<"           q-distribution:      "<<endl;
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
 {
  cout<<"           (with weights)       "<<endl;
 } else
   {
    cout<<"          (without weights)       "<<endl;
   }   
 cout<<endl;
 cout<<"1.) sigma^2 not fitted: "<<endl;
 cout<<"  v_"<<n<<"{FQD} = "<<fIntFlow[0]->GetBinContent(1)<<" +/- "<<fIntFlow[0]->GetBinError(1)<<endl;
 cout<<"  sigma^2 = 0.5 +/- 0 "<<endl; 
 cout<<"  chi^2 = "<<fChi2[0]->GetBinContent(1)<<" (Minuit)"<<endl; 
 cout<<endl;   
 if(TMath::Abs(fIntFlow[0]->GetBinContent(1)-fvMin)<1.e-10 || 
    TMath::Abs(fIntFlow[0]->GetBinContent(1)-fvMax)<1.e-10)
 {
  cout<<"  WARNING: value of v_"<<n<<"{FQD}"<<" is on the boundary"<<endl;
  cout<<"           of fitting interval. Redo the fit."<< endl;
  cout<<endl;
 }
 cout<<"2.) sigma^2 fitted: "<<endl;
 cout<<"  v_"<<n<<"{FQD} = "<<fIntFlow[1]->GetBinContent(1)<<" +/- "<<fIntFlow[1]->GetBinError(1)<<endl;
 cout<<"  sigma^2 = "<<fSigma2[1]->GetBinContent(1)<<" +/- "<<fSigma2[1]->GetBinError(1)<<endl; 
 cout<<"  chi^2 = "<<fChi2[1]->GetBinContent(1)<<" (Minuit)"<<endl; 
 cout<<endl; 
 if(TMath::Abs(fIntFlow[0]->GetBinContent(1)-fvMin)<1.e-10 || 
    TMath::Abs(fIntFlow[0]->GetBinContent(1)-fvMax)<1.e-10)
 {
  cout<<"  WARNING: value of v_"<<n<<"{FQD}"<<" is on the boundary"<<endl;
  cout<<"           of fitting interval. Redo the fit."<< endl;
  cout<<endl;
 }     
 if(TMath::Abs(fSigma2[1]->GetBinContent(1)-fSigma2Min)<1.e-10 || 
    TMath::Abs(fSigma2[1]->GetBinContent(1)-fSigma2Max)<1.e-10)
 {
  cout<<"  WARNING: value of sigma^2 is on the boundary"<<endl;
  cout<<"           of fitting interval. Redo the fit."<< endl;
  cout<<endl;
 }     
 cout<<"      nEvts = "<<fSumOfParticleWeights->GetEntries()<<", AvM = "<<fSumOfParticleWeights->GetMean()<<endl;
 cout<<endl;
 cout<<"***************************************"<<endl;
 cout<<"***************************************"<<endl; 
 cout<<endl;
  
} // end of void AliFlowAnalysisWithFittingQDistribution::PrintOnTheScreen()

//================================================================================================================================ 

void AliFlowAnalysisWithFittingQDistribution::StoreFittingParameters()
{
 // Store fitting parameters for the fit of q-distribution in profile fFittingParameters.
 
 if(!fFittingParameters)
 {
  cout<<endl;
  cout<<"WARNING (FQD): fFittingParameters is NULL in AFAWFQD::SFP() !!!!"<<endl;
  cout<<endl;
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
 fFittingParameters->Fill(7.5,fFinalResultIsFromSigma2Fitted); 
 fFittingParameters->Fill(8.5,fPrintOnTheScreen); 
 
} // end of void AliFlowAnalysisWithFittingQDistribution::StoreFittingParameters()

//================================================================================================================================ 

void AliFlowAnalysisWithFittingQDistribution::AccessFittingParameters()
{
 // Access fitting parameters for the fit of q-distribution.

 fTreshold = fFittingParameters->GetBinContent(1);
 fvStart = fFittingParameters->GetBinContent(2);
 fvMin = fFittingParameters->GetBinContent(3);
 fvMax = fFittingParameters->GetBinContent(4);
 fSigma2Start = fFittingParameters->GetBinContent(5);
 fSigma2Min = fFittingParameters->GetBinContent(6);
 fSigma2Max = fFittingParameters->GetBinContent(7);
 fFinalResultIsFromSigma2Fitted = (Bool_t)fFittingParameters->GetBinContent(8);
 fPrintOnTheScreen = (Bool_t)fFittingParameters->GetBinContent(9);
 
} // end of void AliFlowAnalysisWithFittingQDistribution::AccessFittingParameters()
