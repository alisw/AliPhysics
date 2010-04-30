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

/* $Id$ */

/*************************************************************** 
 * Only in this class nested loops are used for flow analysis. *
 * Nested loops are used to evaluate:                          *
 *                                                             *  
 *  a) Distribution of relative angle difference (phi1-phi2);  *
 *  b) Cross-check the results for mixed harmonics.            *
 *                                                             *
 *       Author: Ante Bilandzic (abilandzic@gmail.com)         *
 ***************************************************************/ 

#define AliFlowAnalysisWithNestedLoops_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TProfile.h"

#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithNestedLoops.h"

class TH1;
class TList;
ClassImp(AliFlowAnalysisWithNestedLoops)

//================================================================================================================

AliFlowAnalysisWithNestedLoops::AliFlowAnalysisWithNestedLoops(): 
fHistList(NULL),
fHistListName(NULL),
fAnalysisLabel(NULL),
fAnalysisSettings(NULL),
fPrintOnTheScreen(kTRUE),
fCommonHists(NULL),
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
fWeightsList(NULL),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fUseParticleWeights(NULL),
fPhiWeights(NULL),
fPtWeights(NULL),
fEtaWeights(NULL),
fListRAD(NULL),
fEvaluateNestedLoopsForRAD(kTRUE),
fRelativeAngleDistribution(NULL),
fListMH(NULL),
fEvaluateNestedLoopsForMH(kFALSE),
fCorrelatorIntegerMH(1),
fCrossCheckInPtSumBinNo(4),  
fCrossCheckInPtDiffBinNo(5) 
{
 // Constructor. 
 
 // Base list to hold all output objects:
 fHistList = new TList();
 fHistListName = new TString("cobjNL");
 fHistList->SetName(fHistListName->Data());
 fHistList->SetOwner(kTRUE);
 
 // List to hold histograms with phi, pt and eta weights:      
 fWeightsList = new TList();
 
 // List to hold objects relevant for relative angle distributions:      
 fListRAD = new TList();
 
 // List holding objects relevant for debugging and cross-checking of MH class: 
 fListMH = new TList();
 
 // Initialize all arrays: 
 this->InitializeArraysForMH();
 
} // AliFlowAnalysisWithNestedLoops::AliFlowAnalysisWithNestedLoops()
 
//================================================================================================================  

AliFlowAnalysisWithNestedLoops::~AliFlowAnalysisWithNestedLoops()
{
 // Destructor.
 
 delete fHistList;

} // end of AliFlowAnalysisWithNestedLoops::~AliFlowAnalysisWithNestedLoops()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::Init()
{
 // Initialize and book all objects. 
 
 // a) Cross check if the user settings make sense before starting; 
 // b) Access all common constants;
 // c) Book and nest all lists in the base list fHistList;
 // d) Book profile holding seetings for analysis with nested loops;
 // e) Book common control histograms;
 // f) Book all objects relevant for distributions;
 // g) Book and fill histograms to hold phi, pt and eta weights.

 //save old value and prevent histograms from being added to directory
 //to avoid name clashes in case multiple analaysis objects are used
 //in an analysis
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
 TH1::AddDirectory(kFALSE);
 
 TH1::SetDefaultSumw2();
  
 this->CrossCheckSettings();
 this->AccessConstants();
 this->BookAndNestAllLists();
 this->BookAndFillProfileHoldingSettings();
 this->BookCommonHistograms();
 this->BookEverythingForRAD();
 this->BookEverythingForMH();
 this->BookAndFillWeightsHistograms();

 //restore old status
 TH1::AddDirectory(oldHistAddStatus);
} // end of void AliFlowAnalysisWithNestedLoops::Init()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::Make(AliFlowEventSimple* anEvent)
{
 // Running over data only in this method.
 
 // a) Check all pointers used in this method;
 // b) Fill common control histograms;
 // c) Evaluate nested loops for relative angle distribution;
 // d) Evaluate nested loops for mixed harmonics.
 
 this->CheckPointersUsedInMake();
 fCommonHists->FillControlHistograms(anEvent);  
 if(fEvaluateNestedLoopsForRAD) this->EvaluateNestedLoopsForRAD(anEvent);
 if(fEvaluateNestedLoopsForMH) this->EvaluateNestedLoopsForMH(anEvent);
  
} // end of AliFlowAnalysisWithNestedLoops::Make(AliFlowEventSimple* anEvent)

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::Finish()
{
 // Calculate the final results.
 
 // a) Check all pointers used in this method;
 // b) Access settings for analysis with mixed harmonics;
 // c) Print on the screen.
 
 this->CheckPointersUsedInFinish();
 this->AccessSettings();
 if(fPrintOnTheScreen) this->PrintOnTheScreen();
                                                                                                                                                                                                                                                                                                               
} // end of AliFlowAnalysisWithNestedLoops::Finish()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::GetOutputHistograms(TList *outputListHistos)
{
 // Get pointers to all objects saved in the output file.
 
 // a) Get pointers for common control histograms. 
 if(outputListHistos)
 {	
  this->SetHistList(outputListHistos);
  if(!fHistList)
  {
   cout<<endl;
   cout<<" WARNING (NL): fHistList is NULL in AFAWNL::GOH() !!!!"<<endl;
   cout<<endl;
   exit(0);
  }
  this->GetPointersForBaseHistograms();
  this->GetPointersForCommonHistograms();
  Bool_t bEvaluateNestedLoopsForRAD = (Bool_t) fAnalysisSettings->GetBinContent(1); // to be improved (not needed here?)
  Bool_t bEvaluateNestedLoopsForMH = (Bool_t) fAnalysisSettings->GetBinContent(2); // to be improved (not needed here?)
  if(bEvaluateNestedLoopsForRAD) this->GetPointersForRAD();
  if(bEvaluateNestedLoopsForMH) this->GetPointersForMH();
 } else 
   {
    cout<<endl;
    cout<<" WARNING (NL): outputListHistos is NULL in AFAWNL::GOH() !!!!"<<endl;
    cout<<endl;
    exit(0);
   }
   
} // end of void AliFlowAnalysisWithNestedLoops::GetOutputHistograms(TList *outputListHistos)

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::GetPointersForBaseHistograms() 
{
 // Get pointers to base histograms.
 
 TString analysisSettingsName = "fAnalysisSettings";
 TProfile *analysisSettings = dynamic_cast<TProfile*>(fHistList->FindObject(analysisSettingsName.Data()));
 if(analysisSettings) 
 {
  this->SetAnalysisSettings(analysisSettings); 
 } else
   {
    cout<<endl;
    cout<<" WARNING (NL): analysisSettings is NULL in AFAWNL::GPFBH() !!!!"<<endl;
    cout<<endl;
    exit(0);  
   }
 
} // end of void AliFlowAnalysisWithNestedLoops::GetPointersForBaseHistograms()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::GetPointersForCommonHistograms() 
{
 // Get pointers to common control histograms.
 
 TString commonHistsName = "AliFlowCommonHistNL";
 AliFlowCommonHist *commonHist = dynamic_cast<AliFlowCommonHist*>(fHistList->FindObject(commonHistsName.Data()));
 if(commonHist) 
 {
  this->SetCommonHists(commonHist); 
 } else
   {
    cout<<endl;
    cout<<" WARNING (NL): commonHist is NULL in AFAWNL::GPFCH() !!!!"<<endl;
    cout<<endl;
    exit(0);  
   }
 
} // end of void AliFlowAnalysisWithNestedLoops::GetPointersForCommonHistograms()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::GetPointersForRAD() 
{
 // Get pointers to objects relevant for relative angle distributions.
 
 TList *listRAD = NULL;
 listRAD = dynamic_cast<TList*>(fHistList->FindObject("Relative Angle Distribution"));
 if(!listRAD) 
 {
  cout<<"WARNING: listRAD is NULL in AFAWNL::GPFRAD() !!!!"<<endl;
  exit(0); 
 }  

 TString relativeAngleDistributionName = "fRelativeAngleDistribution";
 TH1D *relativeAngleDistribution = dynamic_cast<TH1D*>(listRAD->FindObject(relativeAngleDistributionName.Data()));
 if(relativeAngleDistribution)
 {
  this->SetRelativeAngleDistribution(relativeAngleDistribution);  
 }
  
} // end of void AliFlowAnalysisWithNestedLoops::GetPointersForRAD()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::GetPointersForMH() 
{
 // Get pointers to objects evaluated with nested loops.
 
 TList *listMH = NULL;
 listMH = dynamic_cast<TList*>(fHistList->FindObject("Mixed Harmonics"));
 if(!listMH) 
 {
  cout<<"WARNING: listMH is NULL in AFAWNL::GPFMH() !!!!"<<endl;
  exit(0); 
 }  
  
 TString psdFlag[2] = {"PtSum","PtDiff"};
 for(Int_t sd=0;sd<2;sd++)
 {
  TProfile *p3pCorrelatorVsPtSumDiffDirectPro = dynamic_cast<TProfile*>(listMH->FindObject(Form("f3pCorrelatorDirectVs%s",psdFlag[sd].Data())));
  if(p3pCorrelatorVsPtSumDiffDirectPro)
  {
   this->Set3pCorrelatorVsPtSumDiffDirectPro(p3pCorrelatorVsPtSumDiffDirectPro,sd);  
  }  
 } // end of for(Int_t sd=0;sd<2;sd++)
 
} // end of void AliFlowAnalysisWithNestedLoops::GetPointersForMH() 

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::WriteHistograms(TString outputFileName)
{
 // Store the final results in output .root file.
 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);
 delete output;
}

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::WriteHistograms(TDirectoryFile *outputFileName)
{
 // Store the final results in output .root file.
 fHistList->SetName("cobjNL");
 fHistList->SetOwner(kTRUE);
 outputFileName->Add(fHistList);
 outputFileName->Write(outputFileName->GetName(),TObject::kSingleKey);
}

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::BookAndNestAllLists()
{
 // Book and nest all list in base list fHistList.

 // Weights:
 fWeightsList->SetName("Weights");
 fWeightsList->SetOwner(kTRUE);   
 fHistList->Add(fWeightsList); 
 // List for Relative Angle Distribution:
 fListRAD->SetName("Relative Angle Distribution");
 fListRAD->SetOwner(kTRUE);   
 if(fEvaluateNestedLoopsForRAD) fHistList->Add(fListRAD); 
 // List for Mixed Harmonics:
 fListMH->SetName("Mixed Harmonics");
 fListMH->SetOwner(kTRUE);   
 if(fEvaluateNestedLoopsForMH) fHistList->Add(fListMH); 

} // end of void AliFlowAnalysisWithNestedLoops::BookAndNestAllLists()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::BookAndFillProfileHoldingSettings()
{
 // Book profile to hold all analysis settings.

 TString analysisSettingsName = "fAnalysisSettings";
 fAnalysisSettings = new TProfile(analysisSettingsName.Data(),"Settings for analysis with nested loops",6,0,6);
 fAnalysisSettings->GetXaxis()->SetLabelSize(0.035);
 fAnalysisSettings->GetXaxis()->SetBinLabel(1,"Nested loops for RAD?");
 fAnalysisSettings->Fill(0.5,(Int_t)fEvaluateNestedLoopsForRAD);
 fAnalysisSettings->GetXaxis()->SetBinLabel(2,"Nested loops for MH?");
 fAnalysisSettings->Fill(1.5,(Int_t)fEvaluateNestedLoopsForMH);
 fAnalysisSettings->GetXaxis()->SetBinLabel(3,"Integer n in cos(n(2#phi_{1}-#psi_{2}-#psi_{3}))");
 fAnalysisSettings->Fill(2.5,(Int_t)fCorrelatorIntegerMH);
 fAnalysisSettings->GetXaxis()->SetBinLabel(4,"Print on the screen?");
 fAnalysisSettings->Fill(3.5,(Int_t)fPrintOnTheScreen); 
 fAnalysisSettings->GetXaxis()->SetBinLabel(5,"fCrossCheckInPtSumBinNo");
 fAnalysisSettings->Fill(4.5,fCrossCheckInPtSumBinNo);
 fAnalysisSettings->GetXaxis()->SetBinLabel(6,"fCrossCheckInPtDiffBinNo");
 fAnalysisSettings->Fill(5.5,fCrossCheckInPtDiffBinNo);
 fHistList->Add(fAnalysisSettings);
 
} // end of void AliFlowAnalysisWithNestedLoops::BookAndFillProfileHoldingSettings()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::BookCommonHistograms()
{
 // Book common control histograms and common histograms for final results.
 
 TString commonHistsName = "AliFlowCommonHistNL";
 fCommonHists = new AliFlowCommonHist(commonHistsName.Data());
 fHistList->Add(fCommonHists);  
 
} // end of void AliFlowAnalysisWithNestedLoops::BookCommonHistograms()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::BookEverythingForRAD()
{
 // Book all objects relevant calculation of relative angle distribution.
 
 TString relativeAngleDistributionName = "fRelativeAngleDistribution";
 fRelativeAngleDistribution = new TH1D(relativeAngleDistributionName.Data(),"Relative angle distribution",720,-TMath::TwoPi(),TMath::TwoPi());
 fRelativeAngleDistribution->GetYaxis()->SetTitle("#frac{dN}{#Delta #phi}"); 
 fRelativeAngleDistribution->GetXaxis()->SetTitle("#Delta #phi");
 fListRAD->Add(fRelativeAngleDistribution);

} // end fo void AliFlowAnalysisWithNestedLoops::BookEverythingForRAD()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::AccessConstants()
{
 // Access needed common constants from AliFlowCommonConstants.
 
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
 
} // end of void AliFlowAnalysisWithNestedLoops::AccessConstants()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::CrossCheckSettings()
{
 // Cross-check if the user settings make sense. 
 
} // end of void AliFlowAnalysisWithNestedLoops::CrossCheckSettings()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::BookAndFillWeightsHistograms()
{
 // Book and fill (by accessing file "weights.root") histograms which hold phi, pt and eta weights.

 if(!fWeightsList)
 {
  cout<<"WARNING: fWeightsList is NULL in AFAWNL::BAFWH() !!!!"<<endl;
  exit(0);  
 }
 // Profile to hold flags for weights:   
 TString fUseParticleWeightsName = "fUseParticleWeightsNL";
 fUseParticleWeights = new TProfile(fUseParticleWeightsName.Data(),"0 = particle weight not used, 1 = particle weight used ",3,0,3);
 fUseParticleWeights->SetLabelSize(0.06);
 (fUseParticleWeights->GetXaxis())->SetBinLabel(1,"w_{#phi}");
 (fUseParticleWeights->GetXaxis())->SetBinLabel(2,"w_{p_{T}}");
 (fUseParticleWeights->GetXaxis())->SetBinLabel(3,"w_{#eta}");
 fUseParticleWeights->Fill(0.5,(Int_t)fUsePhiWeights);
 fUseParticleWeights->Fill(1.5,(Int_t)fUsePtWeights);
 fUseParticleWeights->Fill(2.5,(Int_t)fUseEtaWeights);
 fWeightsList->Add(fUseParticleWeights); 
 // Phi-weights: 
 if(fUsePhiWeights)
 {
  if(fWeightsList->FindObject("phi_weights"))
  {
   fPhiWeights = dynamic_cast<TH1F*>(fWeightsList->FindObject("phi_weights"));
   if(TMath::Abs(fPhiWeights->GetBinWidth(1)-fPhiBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (NL): Inconsistent binning in histograms for phi-weights throughout the code."<<endl;
    cout<<endl;
    exit(0);
   }
  } else 
    {
     cout<<"WARNING (NL): fWeightsList->FindObject(\"phi_weights\") is NULL in AFAWNL::BAFWH() !!!!"<<endl;
     exit(0);
    }
 } // end of if(fUsePhiWeights)
 // Pt-weights:
 if(fUsePtWeights) 
 {
  if(fWeightsList->FindObject("pt_weights"))
  {
   fPtWeights = dynamic_cast<TH1D*>(fWeightsList->FindObject("pt_weights"));
   if(TMath::Abs(fPtWeights->GetBinWidth(1)-fPtBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (NL): Inconsistent binning in histograms for pt-weights throughout the code."<<endl;
    cout<<endl;
    exit(0);
   }
  } else 
    {
     cout<<"WARNING (NL): fWeightsList->FindObject(\"pt_weights\") is NULL in AFAWNL::BAFWH() !!!!"<<endl;
     exit(0);
    }
 } // end of if(fUsePtWeights)    
 // Eta-weights:
 if(fUseEtaWeights) 
 {
  if(fWeightsList->FindObject("eta_weights"))
  {
   fEtaWeights = dynamic_cast<TH1D*>(fWeightsList->FindObject("eta_weights"));
   if(TMath::Abs(fEtaWeights->GetBinWidth(1)-fEtaBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (NL): Inconsistent binning in histograms for eta-weights throughout the code."<<endl;
    cout<<endl;
    exit(0);
   }
  } else 
    {
     cout<<"WARNING: fUseEtaWeights && fWeightsList->FindObject(\"eta_weights\") is NULL in AFAWNL::BAFWH() !!!!"<<endl;
     exit(0);
    }
 } // end of if(fUseEtaWeights)
 
} // end of AliFlowAnalysisWithNestedLoops::BookAndFillWeightsHistograms()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::CheckPointersUsedInMake()
{
 // Check pointers used in method Make().
                        
 if(fEvaluateNestedLoopsForRAD) CheckPointersForRAD("Make");
 if(fEvaluateNestedLoopsForMH) CheckPointersForMH("Make"); 
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithNestedLoops::CheckPointersUsedInMake()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::CheckPointersUsedInFinish()
{
 // Check pointers used in method Finish().
 
 if(fEvaluateNestedLoopsForRAD) CheckPointersForRAD("Finish");
 if(fEvaluateNestedLoopsForMH) CheckPointersForMH("Finish"); 
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithNestedLoops::CheckPointersUsedInFinish()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::CheckPointersForRAD(TString where)
{
 // Check pointers relevant for calculation of relative angle distribution.
 
 if(!fRelativeAngleDistribution)
 {
  cout<<endl;
  cout<<" WARNING (NL): !fRelativeAngleDistribution is NULL in "<<where.Data()<<"() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 
 if(strcmp(where.Data(),"Make") == 0)
 {
  // Check pointers used only in method Make():
  // ...
 }
 else if(strcmp(where.Data(),"Finish") == 0)
 {
  // Check pointers used only in method Finish():
  // ...
 }

} // end of void AliFlowAnalysisWithNestedLoops::CheckPointersForRAD(TString where)

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::CheckPointersForMH(TString where)
{
 // Check pointers relevant for calculation of mixed harmonics.
 
 for(Int_t sd=0;sd<2;sd++)
 {
  if(!(f3pCorrelatorVsPtSumDiffDirectPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (NL): !"<<Form("f3pCorrelatorVsPtSumDiffDirectPro[%d]",sd)<<" is NULL in "<<where.Data()<<"() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
 }

 if(strcmp(where.Data(),"Make") == 0)
 {
  // Check pointers used only in method Make():
  // ...
 }
 else if(strcmp(where.Data(),"Finish") == 0)
 {
  // Check pointers used only in method Finish():
  // ...
 }

} // end of void AliFlowAnalysisWithNestedLoops::CheckPointersForMH(TString where)

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::AccessSettings()
{
 // Access the settings for analysis.
 
 fEvaluateNestedLoopsForRAD = (Bool_t)fAnalysisSettings->GetBinContent(1);
 fEvaluateNestedLoopsForMH = (Bool_t)fAnalysisSettings->GetBinContent(2);
 fCorrelatorIntegerMH = (Int_t)fAnalysisSettings->GetBinContent(3);
 fPrintOnTheScreen = (Bool_t)fAnalysisSettings->GetBinContent(4);
 fCrossCheckInPtSumBinNo = (Int_t)fAnalysisSettings->GetBinContent(5);
 fCrossCheckInPtDiffBinNo = (Int_t)fAnalysisSettings->GetBinContent(6);
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithNestedLoops::AccessSettings()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::InitializeArraysForMH()
{
 // Initialize arrays mixed harmonics calculations.
 
 for(Int_t sd=0;sd<2;sd++) // sum or difference
 {
  f3pCorrelatorVsPtSumDiffDirectPro[sd] = NULL;
 }
  
} // end of AliFlowAnalysisWithNestedLoops::InitializeArraysForMH()

//================================================================================================================  

void AliFlowAnalysisWithNestedLoops::BookEverythingForMH()
{
 // Book all objects relevant for mixed harmonics.
 
 if(fEvaluateNestedLoopsForMH)
 {
  TString psdFlag[2] = {"PtSum","PtDiff"};
  TString psdTitleFlag[2] = {"(p_{t,1}+p_{t,2})/2","#left|p_{t,1}-p_{t,2}#right|"};
  //TString s3pCorrelatorVsPtSumDiffDirectProName = "f3pCorrelatorVsPtSumDiffDirectPro";
  for(Int_t sd=0;sd<2;sd++)
  {
   // to be improved: hardwired ,fnBinsPt,0.,fPtMax):
   f3pCorrelatorVsPtSumDiffDirectPro[sd] = new TProfile(Form("f3pCorrelatorDirectVs%s",psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
   //f3pCorrelatorVsPtSumDiffDirectPro[sd]->SetLabelSize(0.05);
   //f3pCorrelatorVsPtSumDiffDirectPro[sd]->SetMarkerStyle(25);
   f3pCorrelatorVsPtSumDiffDirectPro[sd]->GetXaxis()->SetTitle(psdTitleFlag[sd].Data());
   fListMH->Add(f3pCorrelatorVsPtSumDiffDirectPro[sd]);
  }  
 } // end of if(fEvaluateNestedLoopsForMH)

} // end of AliFlowAnalysisWithNestedLoops::BookEverythingForMH()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::EvaluateNestedLoopsForRAD(AliFlowEventSimple *anEvent)
{
 // Evaluate nested loops needed for calculation of relative angle distribution.
 
 Double_t dPhi1=0., dPhi2=0.; // azimuthal angles in the laboratory frame
 AliFlowTrackSimple *aftsTrack = NULL; // simple track
  
 // Loop over data and store for each distinct pair phi1-phi2 in fRelativeAngleDistribution:
 Int_t nPrim = anEvent->NumberOfTracks();  // nPrim = total number of primary tracks, i.e. nPrim = nRP + nPOI + rest, where:
                                           // nRP   = # of particles used to determine the reaction plane ("Reference Particles");
                                           // nPOI  = # of particles of interest for a detailed flow analysis ("Particles of Interest");
                                           // rest  = # of particles which are not niether RPs nor POIs.  
 // Start nested loops over data:
 for(Int_t i=0;i<nPrim;i++) 
 { 
  aftsTrack=anEvent->GetTrack(i);
  if(aftsTrack)
  {
   if(!aftsTrack->InRPSelection()) continue; // consider only tracks which are RPs 
   dPhi1 = aftsTrack->Phi();
   for(Int_t j=0;j<nPrim;j++) 
   { 
    if(j==i) continue; // eliminating trivial contribution from autocorrelation
    aftsTrack=anEvent->GetTrack(j);
    if(aftsTrack)
    {
     if(!aftsTrack->InRPSelection()) continue; // consider only tracks which are RPs 
     dPhi2 = aftsTrack->Phi();
     // Fill the histogram:
     fRelativeAngleDistribution->Fill(dPhi1-dPhi2);
    }
   } // end of for(Int_t j=0;j<nPrim;j++)
  } else // to if(aftsTrack)
    {
     cout<<endl;
     cout<<" WARNING (NL): No particle! (i.e. aftsTrack is a NULL pointer in AFAWNL::Make().)"<<endl;
     cout<<endl;       
    }
 } // end of for(Int_t i=0;i<nPrim;i++) 
 
} // end of void AliFlowAnalysisWithNestedLoops::EvaluateNestedLoopsForRAD(AliFlowEventSimple *anEvent)

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::EvaluateNestedLoopsForMH(AliFlowEventSimple *anEvent)
{
 // Evaluate nested loops needed for calculation of mixed harmonics.
 // Remark: phi label azimuthal angle of RP particle and psi label azimuthal angle of POI particle.
 
 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL;
 Double_t phi1=0.,psi2=0.,psi3=0.; // angles in the correlator cos[n(2phi1-psi2-psi3)] 
 Double_t pt2=0.,pt3=0.; // transverse momenta of psi2 and psi3
 Int_t n = fCorrelatorIntegerMH; 
 // Evaluting differential correlator cos[n(2phi1-psi2-psi3)] with three nested loops:
 for(Int_t i1=0;i1<nPrim;i1++)
 {
  aftsTrack=anEvent->GetTrack(i1);
  // RP condition (first particle in the correlator must be RP): 
  if(!(aftsTrack->InRPSelection())) continue;
  phi1 = aftsTrack->Phi();  
  for(Int_t i2=0;i2<nPrim;i2++)
  {
   if(i2==i1) continue;
   aftsTrack = anEvent->GetTrack(i2);
   // POI condition (second particle in the correlator must be POI):
   if(!(aftsTrack->InPOISelection())) continue;
   psi2 = aftsTrack->Phi();
   pt2 = aftsTrack->Pt();
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2) continue;
    aftsTrack=anEvent->GetTrack(i3);
    // POI condition (third particle in the correlator must be POI):
    if(!(aftsTrack->InPOISelection())) continue;
    psi3 = aftsTrack->Phi();
    pt3 = aftsTrack->Pt();   
    // Evaluate and store differential correlator cos[n(2phi1-psi2-psi3)]:
    Double_t ptSum = (pt2+pt3)/2.;
    Double_t ptDiff = TMath::Abs(pt2-pt3);
    Double_t diff3pCorrelator = TMath::Cos(n*(2.*phi1-psi2-psi3));
    f3pCorrelatorVsPtSumDiffDirectPro[0]->Fill(ptSum,diff3pCorrelator,1.);
    f3pCorrelatorVsPtSumDiffDirectPro[1]->Fill(ptDiff,diff3pCorrelator,1.);
   } // end of for(Int_t i3=0;i3<nPrim;i3++)  
  } // end of for(Int_t i2=0;i2<nPrim;i2++)  
 } // end of for(Int_t i1=0;i1<nPrim;i1++)

} // end of void AliFlowAnalysisWithNestedLoops::EvaluateNestedLoopsForMH(AliFlowEventSimple *anEvent)

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::PrintOnTheScreen()
{
 // Print on the screen.
 
 cout<<endl;
 cout<<"****************************************************"<<endl;
 cout<<"****************************************************"<<endl;
 cout<<"                  Nested Loops                 "<<endl; 
 cout<<endl;
 
 if(fEvaluateNestedLoopsForRAD) 
 {
  cout<<"  Evaluated for relative angle distribution."<<endl;
 }
 
 if(fEvaluateNestedLoopsForMH) 
 {
  cout<<"  Evaluated for mixed harmonics."<<endl;
  if(fCorrelatorIntegerMH!=1)
  {
   cout<< "  cos["<<fCorrelatorIntegerMH<<"(2phi1-psi2-psi3)] = "<<endl;
  } else
    {
     cout<< "  cos(2phi1-psi2-psi3) = "<<endl;
    } 
  cout<< "  a) in pt sum bin "<<fCrossCheckInPtSumBinNo<<": "<<
  f3pCorrelatorVsPtSumDiffDirectPro[0]->GetBinContent(fCrossCheckInPtSumBinNo)<<
  " +/- "<<f3pCorrelatorVsPtSumDiffDirectPro[0]->GetBinError(fCrossCheckInPtSumBinNo)<<endl;
  cout<< "  b) in pt diff bin "<<fCrossCheckInPtDiffBinNo<<": "<<
  f3pCorrelatorVsPtSumDiffDirectPro[1]->GetBinContent(fCrossCheckInPtDiffBinNo)<<
  " +/- "<<f3pCorrelatorVsPtSumDiffDirectPro[1]->GetBinError(fCrossCheckInPtDiffBinNo)<<endl;
 }
 
 if(!fEvaluateNestedLoopsForRAD && !fEvaluateNestedLoopsForMH)
 {
  cout<<"  Not evaluated."<<endl;
 }
 cout<<endl;
 cout<<"****************************************************"<<endl;
 cout<<"****************************************************"<<endl;
 cout<<endl;
 
} // end of void AliFlowAnalysisWithNestedLoops::PrintOnTheScreen()


