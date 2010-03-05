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
 *  a) Distribution of relative angle difference (phi1-phi2).  *
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
fResultsList(NULL),
fRelativeAngleDistribution(NULL)
{
 // Constructor. 
  // Base list to hold all output objects:
 fHistList = new TList();
 fHistListName = new TString("cobjNL");
 fHistList->SetName(fHistListName->Data());
 fHistList->SetOwner(kTRUE);
 
 // List to hold histograms with phi, pt and eta weights:      
 fWeightsList = new TList();
 
 // List to hold objects with final results:      
 fResultsList = new TList();
 
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
 
 this->CrossCheckSettings();
 this->AccessConstants();
 this->BookAndNestAllLists();
 this->BookProfileHoldingSettings();
 this->BookCommonHistograms();
 this->BookEverythingForDistributions();
 this->BookAndFillWeightsHistograms();

} // end of void AliFlowAnalysisWithNestedLoops::Init()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::Make(AliFlowEventSimple* anEvent)
{
 // Running over data only in this method.
 
 // a) Check all pointers used in method Make();
 // b) Define local variables;
 // c) Fill common control histograms;
 // d) Loop over data and store for each distinct pair phi1-phi2 in fRelativeAngleDistribution.
 
 // a) Check all pointers used in method Make():
 this->CheckPointersUsedInMake();
 
 // b) Define local variables:
 Double_t dPhi1=0., dPhi2=0.; // azimuthal angles in the laboratory frame
 AliFlowTrackSimple *aftsTrack = NULL; // simple track
 
 // c) Fill common control histograms:
 fCommonHists->FillControlHistograms(anEvent);  
 
 // d) Loop over data and store for each distinct pair phi1-phi2 in fRelativeAngleDistribution:
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
  
} // end of AliFlowAnalysisWithNestedLoops::Make(AliFlowEventSimple* anEvent)
//================================================================================================================

void AliFlowAnalysisWithNestedLoops::Finish()
{
 // Calculate the final results.
 
 // a) Check all pointers used in this method;
 // b) Access settings for analysis with mixed harmonics;
 
 this->CheckPointersUsedInFinish();
 this->AccessSettings();
                                                                                                                                                                                                                                                                                                               
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
  this->GetPointersForCommonHistograms();
  this->GetPointersForResultsHistograms();
 } else 
   {
    cout<<endl;
    cout<<" WARNING (NL): outputListHistos is NULL in AFAWNL::GOH() !!!!"<<endl;
    cout<<endl;
    exit(0);
   }
   
} // end of void AliFlowAnalysisWithNestedLoops::GetOutputHistograms(TList *outputListHistos)

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

void AliFlowAnalysisWithNestedLoops::GetPointersForResultsHistograms() 
{
 // Get pointers to histograms holding final results.
 
 TList *resultsList = NULL;
 resultsList = dynamic_cast<TList*>(fHistList->FindObject("Results"));
 if(!resultsList) 
 {
  cout<<"WARNING: resultsList is NULL in AFAWNL::GPFRH() !!!!"<<endl;
  exit(0); 
 }  
 
 TString relativeAngleDistributionName = "fRelativeAngleDistribution";
 TH1D *relativeAngleDistribution = dynamic_cast<TH1D*>(resultsList->FindObject(relativeAngleDistributionName.Data()));
 if(relativeAngleDistribution)
 {
  this->SetRelativeAngleDistribution(relativeAngleDistribution);  
 }
  
} // end of void AliFlowAnalysisWithNestedLoops::GetPointersForResultsHistograms()

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
 // Results:
 fResultsList->SetName("Results");
 fResultsList->SetOwner(kTRUE);   
 fHistList->Add(fResultsList); 

} // end of void AliFlowAnalysisWithNestedLoops::BookAndNestAllLists()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::BookProfileHoldingSettings()
{
 // Book profile to hold all analysis settings.

 TString analysisSettingsName = "fAnalysisSettings";
 fAnalysisSettings = new TProfile(analysisSettingsName.Data(),"Settings for analysis with nested loops",1,0,1);
 //fAnalysisSettings->GetXaxis()->SetBinLabel(1," ... ");
 //fAnalysisSettings->Fill(0.5, ... );
 fHistList->Add(fAnalysisSettings);
 
} // end of void AliFlowAnalysisWithNestedLoops::BookProfileHoldingSettings()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::BookCommonHistograms()
{
 // Book common control histograms and common histograms for final results.
 
 TString commonHistsName = "AliFlowCommonHistNL";
 fCommonHists = new AliFlowCommonHist(commonHistsName.Data());
 fHistList->Add(fCommonHists);  
 
} // end of void AliFlowAnalysisWithNestedLoops::BookCommonHistograms()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::BookEverythingForDistributions()
{
 // Book all objects relevant for distributions.
 
 // Histogram to hold distribution of phi1-phi2:
 TString relativeAngleDistributionName = "fRelativeAngleDistribution";
 fRelativeAngleDistribution = new TH1D(relativeAngleDistributionName.Data(),"Relative angle distribution",720,-TMath::TwoPi(),TMath::TwoPi());
 fRelativeAngleDistribution->GetYaxis()->SetTitle("#frac{dN}{#Delta #phi}"); 
 fRelativeAngleDistribution->GetXaxis()->SetTitle("#Delta #phi");
 fResultsList->Add(fRelativeAngleDistribution);

} // end fo void AliFlowAnalysisWithNestedLoops::BookEverythingForDistributions()

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
                        
 if(!fRelativeAngleDistribution)
 {
  cout<<endl;
  cout<<" WARNING (NL): !fRelativeAngleDistribution is NULL !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithNestedLoops::CheckPointersUsedInMake()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::CheckPointersUsedInFinish()
{
 // Check pointers used in method Finish().
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithNestedLoops::CheckPointersUsedInFinish()

//================================================================================================================

void AliFlowAnalysisWithNestedLoops::AccessSettings()
{
 // Access the settings for analysis with mixed harmonics.
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithNestedLoops::AccessSettings()

//================================================================================================================
