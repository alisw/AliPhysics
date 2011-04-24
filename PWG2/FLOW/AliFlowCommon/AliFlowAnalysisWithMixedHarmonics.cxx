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

/********************************************************** 
 * In this class azimuthal correlators in mixed harmonics *
 * are implemented in terms of Q-vectors. This approach   *
 * doesn't require evaluation of nested loops. This class *
 * can be used to:                                        *
 *                                                        *  
 *  a) Extract subdominant harmonics (like v1 and v4);    *
 *  b) Study flow of two-particle resonances;             *
 *  c) Study strong parity violation.                     * 
 *                                                        * 
 * Author: Ante Bilandzic (abilandzic@gmail.com)          *
 *********************************************************/ 

#define AliFlowAnalysisWithMixedHarmonics_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithMixedHarmonics.h"

class TH1;
class TList;
ClassImp(AliFlowAnalysisWithMixedHarmonics)

//================================================================================================================
AliFlowAnalysisWithMixedHarmonics::AliFlowAnalysisWithMixedHarmonics(): 
fHistList(NULL),
fHistListName(NULL),
fHarmonic(1),
fAnalysisLabel(NULL),
fAnalysisSettings(NULL),
fNoOfMultipicityBins(100),
fMultipicityBinWidth(1),
fMinMultiplicity(3),
fOppositeChargesPOI(kFALSE),
fEvaluateDifferential3pCorrelator(kFALSE),
fCorrectForDetectorEffects(kFALSE),
fPrintOnTheScreen(kTRUE),
fCalculateVsM(kFALSE),
fShowBinLabelsVsM(kFALSE),
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
fReQnk(NULL),
fImQnk(NULL),
fSpk(NULL),
fProfileList(NULL),
f3pCorrelatorPro(NULL),
f5pCorrelatorPro(NULL),
fNonIsotropicTermsPro(NULL),
f3pCorrelatorVsMPro(NULL),
f3pPOICorrelatorVsM(NULL),
fNonIsotropicTermsVsMPro(NULL),
f2pCorrelatorCosPsiDiffPtDiff(NULL),
f2pCorrelatorCosPsiSumPtDiff(NULL),
f2pCorrelatorSinPsiDiffPtDiff(NULL),
f2pCorrelatorSinPsiSumPtDiff(NULL),
f2pCorrelatorCosPsiDiffPtSum(NULL),
f2pCorrelatorCosPsiSumPtSum(NULL),
f2pCorrelatorSinPsiDiffPtSum(NULL),
f2pCorrelatorSinPsiSumPtSum(NULL),
f2pCorrelatorCosPsiDiffEtaDiff(NULL),
f2pCorrelatorCosPsiSumEtaDiff(NULL),
f2pCorrelatorSinPsiDiffEtaDiff(NULL),
f2pCorrelatorSinPsiSumEtaDiff(NULL),
f2pCorrelatorCosPsiDiffEtaSum(NULL),
f2pCorrelatorCosPsiSumEtaSum(NULL),
f2pCorrelatorSinPsiDiffEtaSum(NULL),
f2pCorrelatorSinPsiSumEtaSum(NULL),
fResultsList(NULL),
f3pCorrelatorHist(NULL),
fDetectorBiasHist(NULL),
f3pCorrelatorVsMHist(NULL),
fDetectorBiasVsMHist(NULL)
{
 // Constructor. 
 
 // Base list to hold all output objects:
 fHistList = new TList();
 fHistListName = new TString("cobjMH");
 fHistList->SetName(fHistListName->Data());
 fHistList->SetOwner(kTRUE);
 
 // List to hold histograms with phi, pt and eta weights:      
 fWeightsList = new TList();
 
 // List to hold all all-event profiles:      
 fProfileList = new TList();
 
 // List to hold objects with final results:      
 fResultsList = new TList();

 // Initialize all arrays:  
 this->InitializeArrays();
 
} // AliFlowAnalysisWithMixedHarmonics::AliFlowAnalysisWithMixedHarmonics()
 
//================================================================================================================  

AliFlowAnalysisWithMixedHarmonics::~AliFlowAnalysisWithMixedHarmonics()
{
 // Destructor.
 
 delete fHistList;

} // end of AliFlowAnalysisWithMixedHarmonics::~AliFlowAnalysisWithMixedHarmonics()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::Init()
{
 // Initialize and book all objects. 
 
 // a) Cross check if the user settings make sense before starting; 
 // b) Access all common constants;
 // c) Book and nest all lists in the base list fHistList;
 // d) Book common control histograms;
 // e) Book all event-by-event quantities;
 // f) Book all all-event quantities;
 // g) Book and fill histograms to hold phi, pt and eta weights;
 // h) Store harmonic n used in cos[n*(phi1+phi2-2phi3)] and cos[n*(psi1+psi2-2phi3)].
  
 //save old value and prevent histograms from being added to directory
 //to avoid name clashes in case multiple analaysis objects are used
 //in an analysis
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
 TH1::AddDirectory(kFALSE);
 
 TH1::SetDefaultSumw2();
 
 this->CrossCheckSettings();
 this->AccessConstants();
 this->BookAndNestAllLists();
 this->BookProfileHoldingSettings();
 this->BookCommonHistograms();
 this->BookAllEventByEventQuantities();
 this->BookAllAllEventQuantities();
 this->BookAndFillWeightsHistograms();
 this->StoreHarmonic();

 TH1::AddDirectory(oldHistAddStatus);
 
} // end of void AliFlowAnalysisWithMixedHarmonics::Init()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::Make(AliFlowEventSimple* anEvent)
{
 // Running over data only in this method.
 
 // a) Check all pointers used in this method;
 // b) Define local variables;
 // c) Fill common control histograms;
 // d) Loop over data and calculate e-b-e quantities Q_{n,k} and S_{p,k};
 // e) Calculate 3-p azimuthal correlator cos[n(phi1+phi2-2*phi3)] and non-isotropic terms in terms of Q_{n,k} and S_{p,k};
 // f) Calculate differential 3-p azimuthal correlator cos[n(psi1+psi2-2*phi3)] in terms of Q_{2n} and p_{n}:
 // g) Reset all event-by-event quantities.
 
 // a) Check all pointers used in this method:
 this->CheckPointersUsedInMake();
 
 // b) Define local variables:
 Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
 Double_t dPt  = 0.; // transverse momentum
 Double_t dEta = 0.; // pseudorapidity
 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight
 AliFlowTrackSimple *aftsTrack = NULL; // simple track
 
 // c) Fill common control histograms:
 fCommonHists->FillControlHistograms(anEvent);  
 
 // d) Loop over data and calculate e-b-e quantities:
 Int_t nPrim = anEvent->NumberOfTracks();  // nPrim = total number of primary tracks, i.e. nPrim = nRP + nPOI
                                           // nRP   = # of particles used to determine the reaction plane ("Reference Particles");
                                           // nPOI  = # of particles of interest for a detailed flow analysis ("Particles of Interest");

 Int_t nRefMult = anEvent->GetReferenceMultiplicity();

 // Start loop over data:
 for(Int_t i=0;i<nPrim;i++) 
 { 
  aftsTrack=anEvent->GetTrack(i);
  if(aftsTrack)
  {
   if(!(aftsTrack->InRPSelection() || aftsTrack->InPOISelection())) continue; // consider only tracks which are either RPs or POIs
   Int_t n = fHarmonic; 
   if(aftsTrack->InRPSelection()) // checking RP condition:
   {    
    dPhi = aftsTrack->Phi();
    dPt  = aftsTrack->Pt();
    dEta = aftsTrack->Eta();
    if(fUsePhiWeights && fPhiWeights && fnBinsPhi) // determine phi-weight for this particle:
    {
     wPhi = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*fnBinsPhi/TMath::TwoPi())));
    }
    if(fUsePtWeights && fPtWeights && fnBinsPt) // determine pt-weight for this particle:
    {
     wPt = fPtWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/fPtBinWidth))); 
    }              
    if(fUseEtaWeights && fEtaWeights && fEtaBinWidth) // determine eta-weight for this particle: 
    {
     wEta = fEtaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/fEtaBinWidth))); 
    } 
    // Calculate Re[Q_{m,k}] and Im[Q_{m,k}], (m = 1,2,3,4,5,6 and k = 0,1,2,3) for this event:
    for(Int_t m=0;m<6;m++) 
    {
     for(Int_t k=0;k<4;k++) // to be improved (what is the maximum k that I need?)
     {
      (*fReQnk)(m,k)+=pow(wPhi*wPt*wEta,k)*TMath::Cos((m+1)*n*dPhi); 
      (*fImQnk)(m,k)+=pow(wPhi*wPt*wEta,k)*TMath::Sin((m+1)*n*dPhi); 
     } 
    }
    // Calculate partially S_{p,k} for this event (final calculation of S_{p,k} follows after the loop over data bellow):
    for(Int_t p=0;p<4;p++) // to be improved (what is maximum p that I need?)
    {
     for(Int_t k=0;k<4;k++) // to be improved (what is maximum k that I need?)
     {     
      (*fSpk)(p,k)+=pow(wPhi*wPt*wEta,k);
     }
    }    
   } // end of if(aftsTrack->InRPSelection())
   // POIs:
   if(fEvaluateDifferential3pCorrelator)
   {
    if(aftsTrack->InPOISelection()) // 1st POI
    {
     Double_t dPsi1 = aftsTrack->Phi();
     Double_t dPt1 = aftsTrack->Pt();
     Double_t dEta1 = aftsTrack->Eta();
     Int_t iCharge1 = aftsTrack->Charge();
     Bool_t b1stPOIisAlsoRP = kFALSE;
     if(aftsTrack->InRPSelection()){b1stPOIisAlsoRP = kTRUE;}
     for(Int_t j=0;j<nPrim;j++)
     {
      if(j==i){continue;}
      aftsTrack=anEvent->GetTrack(j);
      if(aftsTrack->InPOISelection()) // 2nd POI
      {
       Double_t dPsi2 = aftsTrack->Phi();
       Double_t dPt2 = aftsTrack->Pt(); 
       Double_t dEta2 = aftsTrack->Eta();
       Int_t iCharge2 = aftsTrack->Charge();
       if(fOppositeChargesPOI && iCharge1 == iCharge2){continue;}
       Bool_t b2ndPOIisAlsoRP = kFALSE;
       if(aftsTrack->InRPSelection()){b2ndPOIisAlsoRP = kTRUE;}

       // Fill:Pt
       fRePEBE[0]->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi1+dPsi2)),1.);
       fImPEBE[0]->Fill((dPt1+dPt2)/2.,TMath::Sin(n*(dPsi1+dPsi2)),1.);
       fRePEBE[1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi1+dPsi2)),1.);
       fImPEBE[1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Sin(n*(dPsi1+dPsi2)),1.);

       // Fill:Eta
       fReEtaEBE[0]->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi1+dPsi2)),1.);
       fImEtaEBE[0]->Fill((dEta1+dEta2)/2.,TMath::Sin(n*(dPsi1+dPsi2)),1.);
       fReEtaEBE[1]->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi1+dPsi2)),1.);
       fImEtaEBE[1]->Fill(TMath::Abs(dEta1-dEta2),TMath::Sin(n*(dPsi1+dPsi2)),1.);

       //=========================================================//
       //2particle correlator <cos(n*(psi1 - ps12))> vs |Pt1-Pt2|
       f2pCorrelatorCosPsiDiffPtDiff->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi1-dPsi2)));
       f2pCorrelatorCosPsiSumPtDiff->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi1+dPsi2)));
       f2pCorrelatorSinPsiDiffPtDiff->Fill(TMath::Abs(dPt1-dPt2),TMath::Sin(n*(dPsi1-dPsi2)));
       f2pCorrelatorSinPsiSumPtDiff->Fill(TMath::Abs(dPt1-dPt2),TMath::Sin(n*(dPsi1+dPsi2)));
       //_______________________________________________________//
       //2particle correlator <cos(n*(psi1 - ps12))> vs (Pt1+Pt2)/2
       f2pCorrelatorCosPsiDiffPtSum->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi1-dPsi2)));
       f2pCorrelatorCosPsiSumPtSum->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi1+dPsi2)));
       f2pCorrelatorSinPsiDiffPtSum->Fill((dPt1+dPt2)/2.,TMath::Sin(n*(dPsi1-dPsi2)));
       f2pCorrelatorSinPsiSumPtSum->Fill((dPt1+dPt2)/2.,TMath::Sin(n*(dPsi1+dPsi2)));
       //_______________________________________________________//
       //2particle correlator <cos(n*(psi1 - ps12))> vs |eta1-eta2|
       f2pCorrelatorCosPsiDiffEtaDiff->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi1-dPsi2)));
       f2pCorrelatorCosPsiSumEtaDiff->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi1+dPsi2)));
       f2pCorrelatorSinPsiDiffEtaDiff->Fill(TMath::Abs(dEta1-dEta2),TMath::Sin(n*(dPsi1-dPsi2)));
       f2pCorrelatorSinPsiSumEtaDiff->Fill(TMath::Abs(dEta1-dEta2),TMath::Sin(n*(dPsi1+dPsi2)));
       //_______________________________________________________//
       //2particle correlator <cos(n*(psi1 - ps12))> vs (Pt1+Pt2)/2
       f2pCorrelatorCosPsiDiffEtaSum->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi1-dPsi2)));
       f2pCorrelatorCosPsiSumEtaSum->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi1+dPsi2)));
       f2pCorrelatorSinPsiDiffEtaSum->Fill((dEta1+dEta2)/2.,TMath::Sin(n*(dPsi1-dPsi2)));
       f2pCorrelatorSinPsiSumEtaSum->Fill((dEta1+dEta2)/2.,TMath::Sin(n*(dPsi1+dPsi2)));
       //=========================================================//
       
       if(b1stPOIisAlsoRP)
       {
        fOverlapEBE[0][0]->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE[0][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE2[0][0]->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE2[0][1]->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi1-dPsi2)),1.);
       }
       if(b2ndPOIisAlsoRP)
       {
        fOverlapEBE[1][0]->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE[1][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE2[1][0]->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE2[1][1]->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi1-dPsi2)),1.);
       }
      } // end of if(aftsTrack->InPOISelection()) // 2nd POI
     } // end of for(Int_t j=i+1;j<nPrim;j++)
    } // end of if(aftsTrack->InPOISelection()) // 1st POI  
   } // end of if(fEvaluateDifferential3pCorrelator)
  } else // to if(aftsTrack)
    {
     cout<<endl;
     cout<<" WARNING (MH): No particle! (i.e. aftsTrack is a NULL pointer in Make().)"<<endl;
     cout<<endl;       
    }
 } // end of for(Int_t i=0;i<nPrim;i++) 

 // Calculate the final expressions for S_{p,k}:
 for(Int_t p=0;p<4;p++) // to be improved (what is maximum p that I need?)
 {
  for(Int_t k=0;k<4;k++) // to be improved (what is maximum k that I need?)
  {
   (*fSpk)(p,k)=pow((*fSpk)(p,k),p+1);
  }  
 } 
 
 // e) Calculate 3-p correlator cos[n(phi1+phi2-2*phi3)] in terms of Q_{n,k} and S_{p,k}:
 if(anEvent->GetEventNSelTracksRP() >= 3) 
 {
  this->Calculate3pCorrelator();
  this->CalculateNonIsotropicTerms();                          
  if(anEvent->GetEventNSelTracksRP() >= 5) 
  {
   this->Calculate5pCorrelator();
  } // end of if(anEvent->GetEventNSelTracksRP() >= 5) 
 } // end of if(anEvent->GetEventNSelTracksRP() >= 3)             
                 
 // f) Calculate differential 3-p azimuthal correlator cos[n(psi1+psi2-2*phi3)] in terms of Q_{2n} and p_{n}:
 if(fEvaluateDifferential3pCorrelator && anEvent->GetEventNSelTracksRP() >= 1)
 {
   Double_t gIntegrated3pCorrelator = 0.;
   this->CalculateDifferential3pCorrelator(gIntegrated3pCorrelator); // to be improved - add relevant if statements for the min # POIs as well

  //3particle correlator vs ref. mult
   if(fCalculateVsM)
     f3pPOICorrelatorVsM->Fill(nRefMult,gIntegrated3pCorrelator);
 }
 
 // g) Reset all event-by-event quantities: 
 this->ResetEventByEventQuantities();
   
} // end of AliFlowAnalysisWithMixedHarmonics::Make(AliFlowEventSimple* anEvent)

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::Finish()
{
 // Calculate the final results.
 
 // a) Check all pointers used in this method;
 // b) Access settings for analysis with mixed harmonics;
 // c) Correct for detector effects;
 // d) Print on the screen the final results.
 
 this->CheckPointersUsedInFinish();
 this->AccessSettings();
 this->CorrectForDetectorEffects();
 if(fCalculateVsM){this->CorrectForDetectorEffectsVsM();}
 if(fPrintOnTheScreen){this->PrintOnTheScreen();}
                                                                                                                                                                                                                                                                                                               
} // end of AliFlowAnalysisWithMixedHarmonics::Finish()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::GetOutputHistograms(TList *outputListHistos)
{
 // Get pointers to all objects saved in the output file.
 
 // a) Get pointers for common control histograms. 
 if(outputListHistos)
 {	
  this->SetHistList(outputListHistos);
  if(!fHistList)
  {
   cout<<endl;
   cout<<" WARNING (MH): fHistList is NULL in GetOutputHistograms() !!!!"<<endl;
   cout<<endl;
   exit(0);
  }
  this->GetPointersForBaseHistograms();
  this->GetPointersForCommonHistograms();
  this->GetPointersForAllEventProfiles();
  this->GetPointersForResultsHistograms();
 } else 
   {
    cout<<endl;
    cout<<" WARNING (MH): outputListHistos is NULL in GetOutputHistograms() !!!!"<<endl;
    cout<<endl;
    exit(0);
   }
   
} // end of void AliFlowAnalysisWithMixedHarmonics::GetOutputHistograms(TList *outputListHistos)

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::GetPointersForBaseHistograms() 
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
    cout<<" WARNING (MH): analysisSettings is NULL in GetPointersForBaseHistograms() !!!!"<<endl;
    cout<<endl;
    exit(0);  
   }
 
} // end of void AliFlowAnalysisWithMixedHarmonics::GetPointersForBaseHistograms()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::GetPointersForCommonHistograms() 
{
 // Get pointers to common control histograms.
 
 TString commonHistsName = "AliFlowCommonHistMH";
 AliFlowCommonHist *commonHist = dynamic_cast<AliFlowCommonHist*>(fHistList->FindObject(commonHistsName.Data()));
 if(commonHist) 
 {
  this->SetCommonHists(commonHist); 
 } else
   {
    cout<<endl;
    cout<<" WARNING (MH): commonHist is NULL in GetPointersForCommonHistograms() !!!!"<<endl;
    cout<<endl;
    exit(0);  
   }
 
} // end of void AliFlowAnalysisWithMixedHarmonics::GetPointersForCommonHistograms()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::GetPointersForAllEventProfiles() 
{
 // Get pointers to profiles holding final results.
 
 TList *profileList = NULL;
 profileList = dynamic_cast<TList*>(fHistList->FindObject("Profiles"));
 if(!profileList) 
 {
  cout<<endl;
  cout<<" WARNING (MH): profileList is NULL in GetPointersForAllEventProfiles() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }  
 
 TString s3pCorrelatorProName = "f3pCorrelatorPro";
 TProfile *p3pCorrelatorPro = dynamic_cast<TProfile*>(profileList->FindObject(s3pCorrelatorProName.Data()));
 if(p3pCorrelatorPro)
 {
  this->Set3pCorrelatorPro(p3pCorrelatorPro);  
 }
 TString s3pCorrelatorVsMProName = "f3pCorrelatorVsMPro";
 TProfile *p3pCorrelatorVsMPro = dynamic_cast<TProfile*>(profileList->FindObject(s3pCorrelatorVsMProName.Data()));
 if(p3pCorrelatorVsMPro)
 {
  this->Set3pCorrelatorVsMPro(p3pCorrelatorVsMPro);  
 }
 TString s3pPOICorrelatorVsMName = "f3pPOICorrelatorVsM";
 TProfile *p3pPOICorrelatorVsM = dynamic_cast<TProfile*>(profileList->FindObject(s3pPOICorrelatorVsMName.Data()));
 if(p3pPOICorrelatorVsM)
 {
  this->Set3pPOICorrelatorVsM(p3pPOICorrelatorVsM);  
 }
 TString nonIsotropicTermsProName = "fNonIsotropicTermsPro";
 TProfile *nonIsotropicTermsPro = dynamic_cast<TProfile*>(profileList->FindObject(nonIsotropicTermsProName.Data()));
 if(nonIsotropicTermsPro)
 {
  this->SetNonIsotropicTermsPro(nonIsotropicTermsPro);  
 }
 TString nonIsotropicTermsVsMProName = "fNonIsotropicTermsVsMPro";
 TProfile2D *nonIsotropicTermsVsMPro = dynamic_cast<TProfile2D*>(profileList->FindObject(nonIsotropicTermsVsMProName.Data()));
 if(nonIsotropicTermsVsMPro)
 {
  this->SetNonIsotropicTermsVsMPro(nonIsotropicTermsVsMPro);  
 } 
 TString psdFlag[2] = {"PtSum","PtDiff"}; 
 for(Int_t sd=0;sd<2;sd++)
 {
  TProfile *p3pCorrelatorVsPtSumDiffPro = dynamic_cast<TProfile*>(profileList->FindObject(Form("f3pCorrelatorVs%sPro",psdFlag[sd].Data())));
  if(p3pCorrelatorVsPtSumDiffPro)
  {
   this->Set3pCorrelatorVsPtSumDiffPro(p3pCorrelatorVsPtSumDiffPro,sd);  
  }
 }

 //2p correlator vs |Pt1-Pt2|
 TProfile *g2pCorrelatorCosPsiDiffPtDiff = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorCosPsiDiffPtDiff"));
 if(g2pCorrelatorCosPsiDiffPtDiff)
   this->Set2pCorrelatorCosPsiDiffPtDiff(g2pCorrelatorCosPsiDiffPtDiff);
 TProfile *g2pCorrelatorCosPsiSumPtDiff = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorCosPsiSumPtDiff"));
 if(g2pCorrelatorCosPsiSumPtDiff)
   this->Set2pCorrelatorCosPsiSumPtDiff(g2pCorrelatorCosPsiSumPtDiff);
 TProfile *g2pCorrelatorSinPsiDiffPtDiff = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorSinPsiDiffPtDiff"));
 if(g2pCorrelatorSinPsiDiffPtDiff)
   this->Set2pCorrelatorSinPsiDiffPtDiff(g2pCorrelatorSinPsiDiffPtDiff);
 TProfile *g2pCorrelatorSinPsiSumPtDiff = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorSinPsiSumPtDiff"));
 if(g2pCorrelatorSinPsiSumPtDiff)
   this->Set2pCorrelatorSinPsiSumPtDiff(g2pCorrelatorSinPsiSumPtDiff);

 //2p correlator vs (Pt1+Pt2)/2
 TProfile *g2pCorrelatorCosPsiDiffPtSum = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorCosPsiDiffPtSum"));
 if(g2pCorrelatorCosPsiDiffPtSum)
   this->Set2pCorrelatorCosPsiDiffPtSum(g2pCorrelatorCosPsiDiffPtSum);
 TProfile *g2pCorrelatorCosPsiSumPtSum = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorCosPsiSumPtSum"));
 if(g2pCorrelatorCosPsiSumPtSum)
   this->Set2pCorrelatorCosPsiSumPtSum(g2pCorrelatorCosPsiSumPtSum);
 TProfile *g2pCorrelatorSinPsiDiffPtSum = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorSinPsiDiffPtSum"));
 if(g2pCorrelatorSinPsiDiffPtSum)
   this->Set2pCorrelatorSinPsiDiffPtSum(g2pCorrelatorSinPsiDiffPtSum);
 TProfile *g2pCorrelatorSinPsiSumPtSum = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorSinPsiSumPtSum"));
 if(g2pCorrelatorSinPsiSumPtSum)
   this->Set2pCorrelatorSinPsiSumPtSum(g2pCorrelatorSinPsiSumPtSum);

 //2p correlator vs |eta1-eta2|
 TProfile *g2pCorrelatorCosPsiDiffEtaDiff = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorCosPsiDiffEtaDiff"));
 if(g2pCorrelatorCosPsiDiffEtaDiff)
   this->Set2pCorrelatorCosPsiDiffEtaDiff(g2pCorrelatorCosPsiDiffEtaDiff);
 TProfile *g2pCorrelatorCosPsiSumEtaDiff = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorCosPsiSumEtaDiff"));
 if(g2pCorrelatorCosPsiSumEtaDiff)
   this->Set2pCorrelatorCosPsiSumEtaDiff(g2pCorrelatorCosPsiSumEtaDiff);
 TProfile *g2pCorrelatorSinPsiDiffEtaDiff = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorSinPsiDiffEtaDiff"));
 if(g2pCorrelatorSinPsiDiffEtaDiff)
   this->Set2pCorrelatorSinPsiDiffEtaDiff(g2pCorrelatorSinPsiDiffEtaDiff);
 TProfile *g2pCorrelatorSinPsiSumEtaDiff = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorSinPsiSumEtaDiff"));
 if(g2pCorrelatorSinPsiSumEtaDiff)
   this->Set2pCorrelatorSinPsiSumEtaDiff(g2pCorrelatorSinPsiSumEtaDiff);

 //2p correlator vs (eta1+eta2)/2
 TProfile *g2pCorrelatorCosPsiDiffEtaSum = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorCosPsiDiffEtaSum"));
 if(g2pCorrelatorCosPsiDiffEtaSum)
   this->Set2pCorrelatorCosPsiDiffEtaSum(g2pCorrelatorCosPsiDiffEtaSum);
 TProfile *g2pCorrelatorCosPsiSumEtaSum = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorCosPsiSumEtaSum"));
 if(g2pCorrelatorCosPsiSumEtaSum)
   this->Set2pCorrelatorCosPsiSumEtaSum(g2pCorrelatorCosPsiSumEtaSum);
 TProfile *g2pCorrelatorSinPsiDiffEtaSum = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorSinPsiDiffEtaSum"));
 if(g2pCorrelatorSinPsiDiffEtaSum)
   this->Set2pCorrelatorSinPsiDiffEtaSum(g2pCorrelatorSinPsiDiffEtaSum);
 TProfile *g2pCorrelatorSinPsiSumEtaSum = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorSinPsiSumEtaSum"));
 if(g2pCorrelatorSinPsiSumEtaSum)
   this->Set2pCorrelatorSinPsiSumEtaSum(g2pCorrelatorSinPsiSumEtaSum);
   
 TString s5pCorrelatorProName = "f5pCorrelatorPro";
 TProfile *p5pCorrelatorPro = dynamic_cast<TProfile*>(profileList->FindObject(s5pCorrelatorProName.Data()));
 if(p5pCorrelatorPro)
 {
  this->Set5pCorrelatorPro(p5pCorrelatorPro);  
 }
  
} // end of void AliFlowAnalysisWithMixedHarmonics::GetPointersForAllEventProfiles()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::GetPointersForResultsHistograms() 
{
 // Get pointers to histograms holding final results.
 
 TList *resultsList = NULL;
 resultsList = dynamic_cast<TList*>(fHistList->FindObject("Results"));
 if(!resultsList) 
 {
  cout<<endl;
  cout<<" WARNING (MH): resultsList is NULL in GetPointersForResultsHistograms() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }  
 TString s3pCorrelatorHistName = "f3pCorrelatorHist";
 TH1D *h3pCorrelatorHist = dynamic_cast<TH1D*>(resultsList->FindObject(s3pCorrelatorHistName.Data()));
 if(h3pCorrelatorHist)
 {
  this->Set3pCorrelatorHist(h3pCorrelatorHist);  
 }
 TString s3pCorrelatorVsMHistName = "f3pCorrelatorVsMHist";
 TH1D *h3pCorrelatorVsMHist = dynamic_cast<TH1D*>(resultsList->FindObject(s3pCorrelatorVsMHistName.Data()));
 if(h3pCorrelatorVsMHist)
 {
  this->Set3pCorrelatorVsMHist(h3pCorrelatorVsMHist);  
 }
 TString detectorBiasHistName = "fDetectorBiasHist";
 TH1D *detectorBiasHist = dynamic_cast<TH1D*>(resultsList->FindObject(detectorBiasHistName.Data()));
 if(detectorBiasHist)
 {
  this->SetDetectorBiasHist(detectorBiasHist);  
 }
 TString detectorBiasVsMHistName = "fDetectorBiasVsMHist";
 TH1D *detectorBiasVsMHist = dynamic_cast<TH1D*>(resultsList->FindObject(detectorBiasVsMHistName.Data()));
 if(detectorBiasVsMHist)
 {
  this->SetDetectorBiasVsMHist(detectorBiasVsMHist);  
 }

} // end of void AliFlowAnalysisWithMixedHarmonics::GetPointersForResultsHistograms()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::WriteHistograms(TString outputFileName)
{
 // Store the final results in output .root file.
 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);
 delete output;
}

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::WriteHistograms(TDirectoryFile *outputFileName)
{
 // Store the final results in output .root file.
 fHistList->SetName("cobjMH");
 fHistList->SetOwner(kTRUE);
 outputFileName->Add(fHistList);
 outputFileName->Write(outputFileName->GetName(),TObject::kSingleKey);
}

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::StoreHarmonic()
{
 // Store harmonic n used in cos[n*(phi1+phi2-2phi3)] and cos[n*(psi1+psi2-2phi3)].

 (fCommonHists->GetHarmonic())->Fill(0.5,fHarmonic);

} // end of void AliFlowAnalysisWithMixedHarmonics::StoreHarmonic()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::InitializeArrays()
{
 // Initialize arrays.
 
 for(Int_t sd=0;sd<2;sd++)
 {
  fRePEBE[sd] = NULL;
  fImPEBE[sd] = NULL;
  fReEtaEBE[sd] = NULL;
  fImEtaEBE[sd] = NULL;
  f3pCorrelatorVsPtSumDiffPro[sd] = NULL;
  f3pCorrelatorVsEtaSumDiffPro[sd] = NULL;
 }
 for(Int_t fs=0;fs<2;fs++) // 1st/2nd POI which is also RP
 {
  for(Int_t sd=0;sd<2;sd++)
  {
   fOverlapEBE[fs][sd] = NULL;
   fOverlapEBE2[fs][sd] = NULL;
  }  
 } 
  
} // end of AliFlowAnalysisWithMixedHarmonics::InitializeArrays()

//================================================================================================================  

void AliFlowAnalysisWithMixedHarmonics::BookAndNestAllLists()
{
 // Book and nest all list in base list fHistList.

 // Weights:
 fWeightsList->SetName("Weights");
 fWeightsList->SetOwner(kTRUE);   
 fHistList->Add(fWeightsList); 
 // Profiles:
 fProfileList->SetName("Profiles");
 fProfileList->SetOwner(kTRUE);   
 fHistList->Add(fProfileList); 
 // Results:
 fResultsList->SetName("Results");
 fResultsList->SetOwner(kTRUE);   
 fHistList->Add(fResultsList); 

} // end of void AliFlowAnalysisWithMixedHarmonics::BookAndNestAllLists()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::BookProfileHoldingSettings()
{
 // Book profile to hold all analysis settings.

 TString analysisSettingsName = "fAnalysisSettings";
 fAnalysisSettings = new TProfile(analysisSettingsName.Data(),"Settings for analysis with mixed harmonics",10,0,10);
 fAnalysisSettings->SetStats(kFALSE);
 fAnalysisSettings->GetXaxis()->SetLabelSize(0.03);
 fAnalysisSettings->GetXaxis()->SetBinLabel(1,"Corr. for det. effects?");
 fAnalysisSettings->Fill(0.5,(Int_t)fCorrectForDetectorEffects); 
 fAnalysisSettings->GetXaxis()->SetBinLabel(2,"# of mult. bins");
 fAnalysisSettings->Fill(1.5,fNoOfMultipicityBins); 
 fAnalysisSettings->GetXaxis()->SetBinLabel(3,"Width of mult. bins");
 fAnalysisSettings->Fill(2.5,fMultipicityBinWidth);  
 fAnalysisSettings->GetXaxis()->SetBinLabel(4,"Minimal mult.");
 fAnalysisSettings->Fill(3.5,fMinMultiplicity);
 fAnalysisSettings->GetXaxis()->SetBinLabel(5,"Print on the screen?");
 fAnalysisSettings->Fill(4.5,(Int_t)fPrintOnTheScreen);
 fAnalysisSettings->GetXaxis()->SetBinLabel(6,"fHarmonic");
 fAnalysisSettings->Fill(5.5,(Int_t)fHarmonic);
 fAnalysisSettings->GetXaxis()->SetBinLabel(7,"fOppositeChargesPOI");
 fAnalysisSettings->Fill(6.5,(Int_t)fOppositeChargesPOI); 
 fAnalysisSettings->GetXaxis()->SetBinLabel(8,"fEvaluateDifferential3pCorrelator");
 fAnalysisSettings->Fill(7.5,(Int_t)fOppositeChargesPOI); 
 fAnalysisSettings->GetXaxis()->SetBinLabel(9,"fCalculateVsM");
 fAnalysisSettings->Fill(8.5,(Int_t)fCalculateVsM);  
 fAnalysisSettings->GetXaxis()->SetBinLabel(10,"fShowBinLabelsVsM");
 fAnalysisSettings->Fill(9.5,(Int_t)fShowBinLabelsVsM); 
 fHistList->Add(fAnalysisSettings);
 
} // end of void AliFlowAnalysisWithMixedHarmonics::BookProfileHoldingSettings()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::BookCommonHistograms()
{
 // Book common control histograms and common histograms for final results.
 
 TString commonHistsName = "AliFlowCommonHistMH";
 fCommonHists = new AliFlowCommonHist(commonHistsName.Data());
 fHistList->Add(fCommonHists);  
 
} // end of void AliFlowAnalysisWithMixedHarmonics::BookCommonHistograms()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::BookAllEventByEventQuantities()
{
 // Book all event-by-event quantitites.
 
  // Q_{n,k} and S{p,k}:
 fReQnk = new TMatrixD(6,9); // to be improved (check bound on k!)
 fImQnk = new TMatrixD(6,9); // to be improved (check bound on k!)
 fSpk = new TMatrixD(4,4); // to be improved (check bound on p and k!)
 
 // p_n vs [(p1+p2)/2,|p1-p2|]
 if(!fEvaluateDifferential3pCorrelator){return;} 
 TString psdFlag[2] = {"PtSum","PtDiff"};
 TString p2sdFlag[2] = {"PtSum","PtDiff"};
 TString fsFlag[2] = {"1st","2nd"};
 for(Int_t sd=0;sd<2;sd++)
 {
  fRePEBE[sd] = new TProfile(Form("fRePEBE%s",psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
  fImPEBE[sd] = new TProfile(Form("fImPEBE%s",psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
  fReEtaEBE[sd] = new TProfile(Form("fReEtaEBE%s",p2sdFlag[sd].Data()),"",fnBinsEta,fEtaMin,fEtaMax);
  fImEtaEBE[sd] = new TProfile(Form("fImEtaEBE%s",p2sdFlag[sd].Data()),"",fnBinsEta,fEtaMin,fEtaMax);
 }  
 for(Int_t fs=0;fs<2;fs++)
 {
  for(Int_t sd=0;sd<2;sd++)
  {
   fOverlapEBE[fs][sd] = new TProfile(Form("%s POI, %s",fsFlag[sd].Data(),psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
   fOverlapEBE2[fs][sd] = new TProfile(Form("%s POI 2, %s",fsFlag[sd].Data(),p2sdFlag[sd].Data()),"",fnBinsEta,fEtaMin,fEtaMax);
  }
 }  

} // end fo void AliFlowAnalysisWithMixedHarmonics::BookAllEventByEventQuantities()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::BookAllAllEventQuantities()
{
 // Book all all-event quantitites.
 
 // a) Book histos and profiles without any binning in multiplicity, pt or eta;
 // b) Book quantites with multiplicity binning;
 // c) Book quantites with binning in (p1+p2)/2 and |p1-p2|. 
  
 this->BookDefault();
 if(fCalculateVsM){this->BookVsM();}
 if(fEvaluateDifferential3pCorrelator){this->BookDifferential();}  
   
} // end of void AliFlowAnalysisWithMixedHarmonics::BookAllAllEventQuantities()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::BookDefault()
{
 // Book histos and profiles without any binning in multiplicity, pt or eta. 
 
 // a) 3-p correlator <<cos[n*(phi1+phi2-2phi3)]>> for all events (not corrected for detector effects);
 // b) Non-isotropic terms in the decomposition of <<cos[n(phi1+phi2-2phi3)]>>;
 // c) 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> corrected for detector effects;
 // d) Histogram which quantifies bias coming from detector inefficiencies to 3-p correlator <<cos[n(phi1+phi2-2phi3)]>>;
 // e) 5-p correlator <<cos[n*(2phi1+2phi2+2phi3-3phi4-3phi5)]>> for all events (not corrected for detector effects - not supported yet).

 // a) 3-p correlator <<cos[n*(phi1+phi2-2phi3)]>> for all events (not corrected for detector effects);
 TString s3pCorrelatorProName = "f3pCorrelatorPro";
 f3pCorrelatorPro = new TProfile(s3pCorrelatorProName.Data(),"",1,0,1);
 f3pCorrelatorPro->SetStats(kFALSE);
 f3pCorrelatorPro->GetXaxis()->SetLabelOffset(0.01);
 f3pCorrelatorPro->GetXaxis()->SetLabelSize(0.05);
 if(fHarmonic == 1)
 {
  f3pCorrelatorPro->GetXaxis()->SetBinLabel(1,"#LT#LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT#GT");
 } else
   {
    f3pCorrelatorPro->GetXaxis()->SetBinLabel(1,Form("#LT#LTcos[%i(#phi_{1}+#phi_{2}-2#phi_{3})]#GT#GT",fHarmonic));
   }
 fProfileList->Add(f3pCorrelatorPro);
 
 // b) Non-isotropic terms in the decomposition of <<cos[n(phi1+phi2-2phi3)]>>:
 TString nonIsotropicTermsProName = "fNonIsotropicTermsPro";
 fNonIsotropicTermsPro = new TProfile(nonIsotropicTermsProName.Data(),"",8,0,8);
 fNonIsotropicTermsPro->SetStats(kFALSE);
 if(fHarmonic == 1)
 {
  fNonIsotropicTermsPro->SetTitle("Non-isotropic terms in decomposition of #LT#LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT#GT");
  fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(1,"cos(#phi_{1})");
  fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(2,"sin(#phi_{1})");
  fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(3,"cos(2#phi_{1})");
  fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(4,"sin(2#phi_{1})");
  fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(5,"cos(#phi_{1}+#phi_{2})");
  fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(6,"sin(#phi_{1}+#phi_{2})");
  fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(7,"cos(2#phi_{1}-#phi_{2})");
  fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(8,"sin(2#phi_{1}-#phi_{2})");
  // fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(9,"cos(#phi_{1}-#phi_{2}-#phi_{3})"); // not needed
  // fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(10,"sin(#phi_{1}-#phi_{2}-#phi_{3})"); // not needed  
 } else
   {
    fNonIsotropicTermsPro->SetTitle(Form("Non-isotropic terms in decomposition of #LT#LTcos[%i(#phi_{1}+#phi_{2}-2#phi_{3})]#GT#GT",fHarmonic));
    fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(1,Form("cos(%d#phi_{1})",fHarmonic));
    fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(2,Form("sin(%d#phi_{1})",fHarmonic));
    fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(3,Form("cos(%d#phi_{1})",2*fHarmonic));
    fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(4,Form("sin(%d#phi_{1})",2*fHarmonic));
    fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(5,Form("cos[%d(#phi_{1}+#phi_{2})]",fHarmonic));
    fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(6,Form("sin[%d(#phi_{1}+#phi_{2})]",fHarmonic));
    fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(7,Form("cos[%d(2#phi_{1}-#phi_{2})]",fHarmonic));
    fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(8,Form("sin[%d(2#phi_{1}-#phi_{2})]",fHarmonic));
    // fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(9,Form("cos(%d(#phi_{1}-#phi_{2}-#phi_{3}))",fHarmonic)); // not needed
    // fNonIsotropicTermsPro->GetXaxis()->SetBinLabel(10,Form("sin(%d(#phi_{1}-#phi_{2}-#phi_{3}))",fHarmonic)); // not needed
   } 
 fProfileList->Add(fNonIsotropicTermsPro);
 
 // c) 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> corrected for detector effects:
 TString s3pCorrelatorHistName = "f3pCorrelatorHist";
 f3pCorrelatorHist = new TH1D(s3pCorrelatorHistName.Data(),"",1,0,1);
 f3pCorrelatorHist->SetStats(kFALSE);
 f3pCorrelatorHist->GetXaxis()->SetLabelOffset(0.01);
 f3pCorrelatorHist->GetXaxis()->SetLabelSize(0.05);
 if(fHarmonic == 1)
 {
  f3pCorrelatorHist->GetXaxis()->SetBinLabel(1,"#LT#LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT#GT");
 } else
   {
    f3pCorrelatorHist->GetXaxis()->SetBinLabel(1,Form("#LT#LTcos[%i(#phi_{1}+#phi_{2}-2#phi_{3})]#GT#GT",fHarmonic)); 
   }
 fResultsList->Add(f3pCorrelatorHist);
 
 // d) Histogram which quantifies bias coming from detector inefficiencies to 3-p correlator <<cos[n(phi1+phi2-2phi3)]>>:
 TString detectorBiasHistName = "fDetectorBiasHist";
 fDetectorBiasHist = new TH1D(detectorBiasHistName.Data(),"Bias coming from detector inefficiences",1,0,1);
 fDetectorBiasHist->SetStats(kFALSE);
 if(fHarmonic == 1)
 {
  fDetectorBiasHist->GetXaxis()->SetBinLabel(1,"#frac{corrected}{measured} #LT#LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT#GT");
 } else
   {
    fDetectorBiasHist->GetXaxis()->SetBinLabel(1,Form("#frac{corrected}{measured} #LT#LTcos[%i(#phi_{1}+#phi_{2}-2#phi_{3})]#GT#GT",fHarmonic));  
   }  
 fResultsList->Add(fDetectorBiasHist);
 
 // e) 5-p correlator <<cos[n*(2phi1+2phi2+2phi3-3phi4-3phi5)]>> for all events (not corrected for detector effects - not supported yet):
 TString s5pCorrelatorProName = "f5pCorrelatorPro";
 f5pCorrelatorPro = new TProfile(s5pCorrelatorProName.Data(),"",1,0,1);
 f5pCorrelatorPro->SetStats(kFALSE);
 f5pCorrelatorPro->GetXaxis()->SetLabelOffset(0.01);
 f5pCorrelatorPro->GetXaxis()->SetLabelSize(0.05);
 if(fHarmonic == 1)
 {
  f5pCorrelatorPro->GetXaxis()->SetBinLabel(1,"#LT#LTcos(2#phi_{1}+2#phi_{2}+2#phi_{3}-3#phi_{4}-3#phi_{5})#GT#GT");
 } else
   {
    f5pCorrelatorPro->GetXaxis()->SetBinLabel(1,Form("#LT#LTcos[%i(2#phi_{1}+2#phi_{2}+2#phi_{3}-3#phi_{4}-3#phi_{5})]#GT#GT",fHarmonic));
   }
 fProfileList->Add(f5pCorrelatorPro);
 
} // end of void AliFlowAnalysisWithMixedHarmonics::BookDefault()     
      
//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::BookVsM()
{
 // Book histos and profiles holding results vs multiplicity. 

 // a) 3-p correlator <<cos[n*(phi1+phi2-2phi3)]>> for all events (not corrected for detector effects) vs M;
 // b) Non-isotropic terms in the decomposition of <<cos[n(phi1+phi2-2phi3)]>> vs M;
 // c) 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> corrected for detector effects vs M;
 // d) Histogram which quantifies bias coming from detector inefficiencies to 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> vs M.

 // a) 3-p correlator <<cos[n*(phi1+phi2-2phi3)]>> for all events (not corrected for detector effects) vs M:
 TString s3pCorrelatorVsMProName = "f3pCorrelatorVsMPro";
 f3pCorrelatorVsMPro = new TProfile(s3pCorrelatorVsMProName.Data(),"",fNoOfMultipicityBins+2,0,fNoOfMultipicityBins+2);
 f3pCorrelatorVsMPro->SetStats(kFALSE); 
 if(fHarmonic == 1)
 {
  f3pCorrelatorVsMPro->SetTitle("#LT#LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT#GT #font[72]{vs} M");
 } else
   {
    f3pCorrelatorVsMPro->SetTitle(Form("#LT#LTcos[%d(#phi_{1}+#phi_{2}-2#phi_{3})]#GT#GT #font[72]{vs} M",fHarmonic)); 
   }
 if(fShowBinLabelsVsM)
 {
  f3pCorrelatorVsMPro->GetXaxis()->SetBinLabel(1,Form("M < %d",(Int_t)fMinMultiplicity));
  for(Int_t b=2;b<=fNoOfMultipicityBins+1;b++)
  {
   f3pCorrelatorVsMPro->GetXaxis()->SetBinLabel(b,Form("%d #leq M < %d",(Int_t)(fMinMultiplicity+(b-2)*fMultipicityBinWidth),(Int_t)(fMinMultiplicity+(b-1)*fMultipicityBinWidth)));
  }
  f3pCorrelatorVsMPro->GetXaxis()->SetBinLabel(fNoOfMultipicityBins+2,Form(" M #geq %d",(Int_t)(fMinMultiplicity+fNoOfMultipicityBins*fMultipicityBinWidth)));
 } else
   {
    f3pCorrelatorVsMPro->GetXaxis()->SetTitle("M");
   }
 fProfileList->Add(f3pCorrelatorVsMPro); 

 TString s3pPOICorrelatorVsMName = "f3pPOICorrelatorVsM";
 f3pPOICorrelatorVsM = new TProfile(s3pPOICorrelatorVsMName.Data(),"",fNoOfMultipicityBins+2,0,fNoOfMultipicityBins+2);
 f3pPOICorrelatorVsM->SetStats(kFALSE); 
 if(fHarmonic == 1)
 {
  f3pPOICorrelatorVsM->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{3})#GT#GT #font[72]{vs} M");
 } else
   {
    f3pPOICorrelatorVsM->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2}-2#phi_{3})]#GT#GT #font[72]{vs} M",fHarmonic)); 
   }
 if(fShowBinLabelsVsM)
 {
  f3pPOICorrelatorVsM->GetXaxis()->SetBinLabel(1,Form("M < %d",(Int_t)fMinMultiplicity));
  for(Int_t b=2;b<=fNoOfMultipicityBins+1;b++)
  {
   f3pPOICorrelatorVsM->GetXaxis()->SetBinLabel(b,Form("%d #leq M < %d",(Int_t)(fMinMultiplicity+(b-2)*fMultipicityBinWidth),(Int_t)(fMinMultiplicity+(b-1)*fMultipicityBinWidth)));
  }
  f3pPOICorrelatorVsM->GetXaxis()->SetBinLabel(fNoOfMultipicityBins+2,Form(" M #geq %d",(Int_t)(fMinMultiplicity+fNoOfMultipicityBins*fMultipicityBinWidth)));
 } else
   {
    f3pPOICorrelatorVsM->GetXaxis()->SetTitle("M");
   }
 fProfileList->Add(f3pPOICorrelatorVsM); 
 
 // b) Non-isotropic terms in the decomposition of <<cos[n(phi1+phi2-2phi3)]>> vs M:
 TString s3pCorrelatorVsMHistName = "f3pCorrelatorVsMHist";
 f3pCorrelatorVsMHist = new TH1D(s3pCorrelatorVsMHistName.Data(),"",fNoOfMultipicityBins+2,0,fNoOfMultipicityBins+2);
 f3pCorrelatorVsMHist->SetStats(kFALSE); 
 if(fHarmonic == 1)
 {
  f3pCorrelatorVsMHist->SetTitle("cos(#phi_{1}+#phi_{2}-2#phi_{3}) #font[72]{vs} M");
 } else
   {
    f3pCorrelatorVsMHist->SetTitle(Form("cos[%d(#phi_{1}+#phi_{2}-2#phi_{3})] #font[72]{vs} M",fHarmonic)); 
   }
 if(fShowBinLabelsVsM)
 {   
  f3pCorrelatorVsMHist->GetXaxis()->SetBinLabel(1,Form("M < %d",(Int_t)fMinMultiplicity));
  for(Int_t b=2;b<=fNoOfMultipicityBins+1;b++)
  {
   f3pCorrelatorVsMHist->GetXaxis()->SetBinLabel(b,Form("%d #leq M < %d",(Int_t)(fMinMultiplicity+(b-2)*fMultipicityBinWidth),(Int_t)(fMinMultiplicity+(b-1)*fMultipicityBinWidth)));
  }
  f3pCorrelatorVsMHist->GetXaxis()->SetBinLabel(fNoOfMultipicityBins+2,Form(" M #geq %d",(Int_t)(fMinMultiplicity+fNoOfMultipicityBins*fMultipicityBinWidth)));
 } else
   {
    f3pCorrelatorVsMHist->GetXaxis()->SetTitle("M");   
   }
 fResultsList->Add(f3pCorrelatorVsMHist);
 
 // c) 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> corrected for detector effects vs M:
 TString nonIsotropicTermsVsMProName = "fNonIsotropicTermsVsMPro";
 fNonIsotropicTermsVsMPro = new TProfile2D(nonIsotropicTermsVsMProName.Data(),"",8,0,8,fNoOfMultipicityBins+2,0,fNoOfMultipicityBins+2);
 fNonIsotropicTermsVsMPro->SetStats(kFALSE);
 if(fHarmonic == 1)
 {
  fNonIsotropicTermsVsMPro->SetTitle("Non-isotropic terms in decomposition of #LT#LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT#GT #font[72]{vs} M");
  fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(1,"cos(#phi_{1})");
  fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(2,"sin(#phi_{1})");
  fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(3,"cos(2#phi_{1})");
  fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(4,"sin(2#phi_{1})");
  fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(5,"cos(#phi_{1}+#phi_{2})");
  fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(6,"sin(#phi_{1}+#phi_{2})");
  fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(7,"cos(2#phi_{1}-#phi_{2})");
  fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(8,"sin(2#phi_{1}-#phi_{2})");
  // fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(9,"cos(#phi_{1}-#phi_{2}-#phi_{3})"); // not needed
  // fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(10,"sin(#phi_{1}-#phi_{2}-#phi_{3})"); // not needed  
 } else
   {
    fNonIsotropicTermsVsMPro->SetTitle(Form("Non-isotropic terms in decomposition of #LT#LTcos[%d(#phi_{1}+#phi_{2}-2#phi_{3})]#GT#GT",fHarmonic));
    fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(1,Form("cos(%d#phi_{1})",fHarmonic));
    fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(2,Form("sin(%d#phi_{1})",fHarmonic));
    fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(3,Form("cos(%d#phi_{1})",2*fHarmonic));
    fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(4,Form("sin(%d#phi_{1})",2*fHarmonic));
    fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(5,Form("cos[%d(#phi_{1}+#phi_{2})]",fHarmonic));
    fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(6,Form("sin[%d(#phi_{1}+#phi_{2})]",fHarmonic));
    fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(7,Form("cos[%d(2#phi_{1}-#phi_{2})]",fHarmonic));
    fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(8,Form("sin[%d(2#phi_{1}-#phi_{2})]",fHarmonic));
    // fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(9,Form("cos(%d(#phi_{1}-#phi_{2}-#phi_{3}))",fHarmonic)); // not needed
    // fNonIsotropicTermsVsMPro->GetXaxis()->SetBinLabel(10,Form("sin(%d(#phi_{1}-#phi_{2}-#phi_{3}))",fHarmonic)); // not needed 
   } 
 if(fShowBinLabelsVsM)
 {     
  fNonIsotropicTermsVsMPro->GetYaxis()->SetBinLabel(1,Form("M < %d",(Int_t)fMinMultiplicity));
  for(Int_t b=2;b<=fNoOfMultipicityBins+1;b++)
  {
   fNonIsotropicTermsVsMPro->GetYaxis()->SetBinLabel(b,Form("%d #leq M < %d",(Int_t)(fMinMultiplicity+(b-2)*fMultipicityBinWidth),(Int_t)(fMinMultiplicity+(b-1)*fMultipicityBinWidth)));
  }
  fNonIsotropicTermsVsMPro->GetYaxis()->SetBinLabel(fNoOfMultipicityBins+2,Form(" M #geq %d",(Int_t)(fMinMultiplicity+fNoOfMultipicityBins*fMultipicityBinWidth)));
 } else
   {
    fNonIsotropicTermsVsMPro->GetYaxis()->SetTitle("M");
   }
 fProfileList->Add(fNonIsotropicTermsVsMPro); 
 
 // d) Histogram which quantifies bias coming from detector inefficiencies to 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> vs M:
 TString detectorBiasVsMHistName = "fDetectorBiasVsMHist";
 fDetectorBiasVsMHist = new TH1D(detectorBiasVsMHistName.Data(),"",fNoOfMultipicityBins+2,0,fNoOfMultipicityBins+2);
 fDetectorBiasVsMHist->SetStats(kFALSE);
 if(fHarmonic == 1)
 {
  fDetectorBiasVsMHist->SetTitle("#frac{corrected}{measured} #LT#LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT#GT #font[72]{vs} M"); 
 } else
   {
    fDetectorBiasVsMHist->SetTitle(Form("#frac{corrected}{measured} #LT#LTcos[%d(#phi_{1}+#phi_{2}-2#phi_{3})]#GT#GT #font[72]{vs} M",fHarmonic)); 
   }
 if(fShowBinLabelsVsM)
 {   
  fDetectorBiasVsMHist->GetXaxis()->SetBinLabel(1,Form("M < %d",(Int_t)fMinMultiplicity));
  for(Int_t b=2;b<=fNoOfMultipicityBins+1;b++)
  {
   fDetectorBiasVsMHist->GetXaxis()->SetBinLabel(b,Form("%d #leq M < %d",(Int_t)(fMinMultiplicity+(b-2)*fMultipicityBinWidth),(Int_t)(fMinMultiplicity+(b-1)*fMultipicityBinWidth)));
  }
  fDetectorBiasVsMHist->GetXaxis()->SetBinLabel(fNoOfMultipicityBins+2,Form(" M #geq %d",(Int_t)(fMinMultiplicity+fNoOfMultipicityBins*fMultipicityBinWidth)));
 } else
   {
    fDetectorBiasVsMHist->GetXaxis()->SetTitle("M");
   }
 fResultsList->Add(fDetectorBiasVsMHist);

} // end of void AliFlowAnalysisWithMixedHarmonics::BookVsM()
      
//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::BookDifferential()
{
 // Book histos and profiles holding results vs (p1+p2)/2 and |p1-p2|. 
 
 TString psdFlag[2] = {"PtSum","PtDiff"};
 TString psdTitleFlag[2] = {"(p_{T,1}+ p_{T,2})/2","#left|p_{T,1}- p_{T,2}#right|"};
 TString psdFlag2[2] = {"EtaSum","EtaDiff"};
 TString psdTitleFlag2[2] = {"(#eta_{1}+ #eta_{2})/2","#left|#eta_{1}- #eta_{2}#right|"};
 //TString s3pCorrelatorVsPtSumDiffProName = "f3pCorrelatorVsPtSumDiffPro";
 for(Int_t sd=0;sd<2;sd++)
 {
  f3pCorrelatorVsPtSumDiffPro[sd] = new TProfile(Form("f3pCorrelatorVs%sPro",psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
  f3pCorrelatorVsPtSumDiffPro[sd]->SetStats(kFALSE);
  f3pCorrelatorVsEtaSumDiffPro[sd] = new TProfile(Form("f3pCorrelatorVs%sPro",psdFlag2[sd].Data()),"",fnBinsEta,fEtaMin,fEtaMax);
  f3pCorrelatorVsEtaSumDiffPro[sd]->SetStats(kFALSE);
  //f3pCorrelatorVsPtSumDiffPro[sd]->SetLabelSize(0.05);
  //f3pCorrelatorVsPtSumDiffPro[sd]->SetMarkerStyle(25);
  if(fHarmonic == 1)
  {
   f3pCorrelatorVsPtSumDiffPro[sd]->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{3})#GT#GT #font[72]{vs} %s",psdTitleFlag[sd].Data())); 
   f3pCorrelatorVsEtaSumDiffPro[sd]->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{3})#GT#GT #font[72]{vs} %s",psdTitleFlag2[sd].Data())); 
  } else
    {
     f3pCorrelatorVsPtSumDiffPro[sd]->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2}-2#phi_{3})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag[sd].Data())); 
     f3pCorrelatorVsEtaSumDiffPro[sd]->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2}-2#phi_{3})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag2[sd].Data())); 
    }   
  f3pCorrelatorVsPtSumDiffPro[sd]->GetXaxis()->SetTitle(psdTitleFlag[sd].Data());
  fProfileList->Add(f3pCorrelatorVsPtSumDiffPro[sd]);
  f3pCorrelatorVsEtaSumDiffPro[sd]->GetXaxis()->SetTitle(psdTitleFlag2[sd].Data());
  fProfileList->Add(f3pCorrelatorVsEtaSumDiffPro[sd]);
 }

 //2p correlator vs |Pt1-Pt2|
 f2pCorrelatorCosPsiDiffPtDiff = new TProfile("f2pCorrelatorCosPsiDiffPtDiff",";|p_{T,1}-p_{T,2}| (GeV/c);#LT #LT cos(#psi_{1} - #psi_{2}) #GT #GT",fnBinsPt,0.,fPtMax);
 f2pCorrelatorCosPsiDiffPtDiff->SetStats(kFALSE);
 f2pCorrelatorCosPsiSumPtDiff = new TProfile("f2pCorrelatorCosPsiSumPtDiff",";|p_{T,1}-p_{T,2}| (GeV/c);#LT #LT cos(#psi_{1} + #psi_{2}) #GT #GT",fnBinsPt,0.,fPtMax);
 f2pCorrelatorCosPsiSumPtDiff->SetStats(kFALSE);
 f2pCorrelatorSinPsiDiffPtDiff = new TProfile("f2pCorrelatorSinPsiDiffPtDiff",";|p_{T,1}-p_{T,2}| (GeV/c);#LT #LT sin(#psi_{1} - #psi_{2}) #GT #GT",fnBinsPt,0.,fPtMax);
 f2pCorrelatorSinPsiDiffPtDiff->SetStats(kFALSE);
 f2pCorrelatorSinPsiSumPtDiff = new TProfile("f2pCorrelatorSinPsiSumPtDiff",";|p_{T,1}-p_{T,2}| (GeV/c);#LT #LT sin(#psi_{1} + #psi_{2}) #GT #GT",fnBinsPt,0.,fPtMax);
 f2pCorrelatorSinPsiSumPtDiff->SetStats(kFALSE);
 if(fHarmonic == 1) {
   f2pCorrelatorCosPsiDiffPtDiff->SetTitle(Form("#LT#LTcos(#psi_{1}-#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag[1].Data())); 
   f2pCorrelatorCosPsiSumPtDiff->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag[1].Data())); 
   f2pCorrelatorSinPsiDiffPtDiff->SetTitle(Form("#LT#LTsin(#psi_{1}-#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag[1].Data())); 
   f2pCorrelatorSinPsiSumPtDiff->SetTitle(Form("#LT#LTsin(#psi_{1}+#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag[1].Data())); 
 }
 else {
   f2pCorrelatorCosPsiDiffPtDiff->SetTitle(Form("#LT#LTcos[%d(#psi_{1}-#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag[1].Data()));
   f2pCorrelatorCosPsiSumPtDiff->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag[1].Data()));
   f2pCorrelatorSinPsiDiffPtDiff->SetTitle(Form("#LT#LTsin[%d(#psi_{1}-#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag[1].Data()));
   f2pCorrelatorSinPsiSumPtDiff->SetTitle(Form("#LT#LTsin[%d(#psi_{1}+#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag[1].Data()));
 }
 fProfileList->Add(f2pCorrelatorCosPsiDiffPtDiff);
 fProfileList->Add(f2pCorrelatorCosPsiSumPtDiff);
 fProfileList->Add(f2pCorrelatorSinPsiDiffPtDiff);
 fProfileList->Add(f2pCorrelatorSinPsiSumPtDiff);

 //2p correlator vs (Pt1+Pt2)/2
 f2pCorrelatorCosPsiDiffPtSum = new TProfile("f2pCorrelatorCosPsiDiffPtSum",";(p_{T,1}+p_{T,2})/2 (GeV/c);#LT #LT cos(#psi_{1} - #psi_{2}) #GT #GT",fnBinsPt,0.,fPtMax);
 f2pCorrelatorCosPsiDiffPtSum->SetStats(kFALSE);
 f2pCorrelatorCosPsiSumPtSum = new TProfile("f2pCorrelatorCosPsiSumPtSum",";(p_{T,1}+p_{T,2})/2 (GeV/c);#LT #LT cos(#psi_{1} + #psi_{2}) #GT #GT",fnBinsPt,0.,fPtMax);
 f2pCorrelatorCosPsiSumPtSum->SetStats(kFALSE);
 f2pCorrelatorSinPsiDiffPtSum = new TProfile("f2pCorrelatorSinPsiDiffPtSum",";(p_{T,1}+p_{T,2})/2 (GeV/c);#LT #LT sin(#psi_{1} - #psi_{2}) #GT #GT",fnBinsPt,0.,fPtMax);
 f2pCorrelatorSinPsiDiffPtSum->SetStats(kFALSE);
 f2pCorrelatorSinPsiSumPtSum = new TProfile("f2pCorrelatorSinPsiSumPtSum",";(p_{T,1}+p_{T,2})/2 (GeV/c);#LT #LT sin(#psi_{1} + #psi_{2}) #GT #GT",fnBinsPt,0.,fPtMax);
 f2pCorrelatorSinPsiSumPtSum->SetStats(kFALSE);
 if(fHarmonic == 1) {
   f2pCorrelatorCosPsiDiffPtSum->SetTitle(Form("#LT#LTcos(#psi_{1}-#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag[0].Data())); 
   f2pCorrelatorCosPsiSumPtSum->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag[0].Data())); 
   f2pCorrelatorSinPsiDiffPtSum->SetTitle(Form("#LT#LTsin(#psi_{1}-#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag[0].Data())); 
   f2pCorrelatorSinPsiSumPtSum->SetTitle(Form("#LT#LTsin(#psi_{1}+#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag[0].Data())); 
 }
 else {
   f2pCorrelatorCosPsiDiffPtSum->SetTitle(Form("#LT#LTcos[%d(#psi_{1}-#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag[0].Data()));
   f2pCorrelatorCosPsiSumPtSum->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag[0].Data()));
   f2pCorrelatorSinPsiDiffPtSum->SetTitle(Form("#LT#LTsin[%d(#psi_{1}-#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag[0].Data()));
   f2pCorrelatorSinPsiSumPtSum->SetTitle(Form("#LT#LTsin[%d(#psi_{1}+#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag[0].Data()));
 }
 fProfileList->Add(f2pCorrelatorCosPsiDiffPtSum);
 fProfileList->Add(f2pCorrelatorCosPsiSumPtSum);
 fProfileList->Add(f2pCorrelatorSinPsiDiffPtSum);
 fProfileList->Add(f2pCorrelatorSinPsiSumPtSum);

 //2p correlator vs |eta1-eta2|
 f2pCorrelatorCosPsiDiffEtaDiff = new TProfile("f2pCorrelatorCosPsiDiffEtaDiff",";|#eta_{1}-#eta_{2}|;#LT #LT cos(#psi_{1} - #psi_{2}) #GT #GT",fnBinsEta,fEtaMin,fEtaMax);
 f2pCorrelatorCosPsiDiffEtaDiff->SetStats(kFALSE);
 f2pCorrelatorCosPsiSumEtaDiff = new TProfile("f2pCorrelatorCosPsiSumEtaDiff",";|#eta_{1}-#eta_{2}|;#LT #LT cos(#psi_{1} + #psi_{2}) #GT #GT",fnBinsEta,fEtaMin,fEtaMax);
 f2pCorrelatorCosPsiSumEtaDiff->SetStats(kFALSE);
 f2pCorrelatorSinPsiDiffEtaDiff = new TProfile("f2pCorrelatorSinPsiDiffEtaDiff",";|#eta_{1}-#eta_{2}|;#LT #LT sin(#psi_{1} - #psi_{2}) #GT #GT",fnBinsEta,fEtaMin,fEtaMax);
 f2pCorrelatorSinPsiDiffEtaDiff->SetStats(kFALSE);
 f2pCorrelatorSinPsiSumEtaDiff = new TProfile("f2pCorrelatorSinPsiSumEtaDiff",";|#eta_{1}-#eta_{2}|;#LT #LT sin(#psi_{1} + #psi_{2}) #GT #GT",fnBinsEta,fEtaMin,fEtaMax);
 f2pCorrelatorSinPsiSumEtaDiff->SetStats(kFALSE);
 if(fHarmonic == 1) {
   f2pCorrelatorCosPsiDiffEtaDiff->SetTitle(Form("#LT#LTcos(#psi_{1}-#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag2[1].Data())); 
   f2pCorrelatorCosPsiSumEtaDiff->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag2[1].Data())); 
   f2pCorrelatorSinPsiDiffEtaDiff->SetTitle(Form("#LT#LTsin(#psi_{1}-#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag2[1].Data())); 
   f2pCorrelatorSinPsiSumEtaDiff->SetTitle(Form("#LT#LTsin(#psi_{1}+#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag2[1].Data())); 
 }
 else {
   f2pCorrelatorCosPsiDiffEtaDiff->SetTitle(Form("#LT#LTcos[%d(#psi_{1}-#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag2[1].Data()));
   f2pCorrelatorCosPsiSumEtaDiff->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag2[1].Data()));
   f2pCorrelatorSinPsiDiffEtaDiff->SetTitle(Form("#LT#LTsin[%d(#psi_{1}-#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag2[1].Data()));
   f2pCorrelatorSinPsiSumEtaDiff->SetTitle(Form("#LT#LTsin[%d(#psi_{1}+#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag2[1].Data()));
 }
 fProfileList->Add(f2pCorrelatorCosPsiDiffEtaDiff);
 fProfileList->Add(f2pCorrelatorCosPsiSumEtaDiff);
 fProfileList->Add(f2pCorrelatorSinPsiDiffEtaDiff);
 fProfileList->Add(f2pCorrelatorSinPsiSumEtaDiff);

 //2p correlator vs (eta1+eta2)/2
 f2pCorrelatorCosPsiDiffEtaSum = new TProfile("f2pCorrelatorCosPsiDiffEtaSum",";(#eta_{1}+#eta_{,2})/2;#LT #LT cos(#psi_{1} - #psi_{2}) #GT #GT",fnBinsEta,fEtaMin,fEtaMax);
 f2pCorrelatorCosPsiDiffEtaSum->SetStats(kFALSE);
 f2pCorrelatorCosPsiSumEtaSum = new TProfile("f2pCorrelatorCosPsiSumEtaSum",";(#eta_{1}+#eta_{,2})/2;#LT #LT cos(#psi_{1} + #psi_{2}) #GT #GT",fnBinsEta,fEtaMin,fEtaMax);
 f2pCorrelatorCosPsiSumEtaSum->SetStats(kFALSE);
 f2pCorrelatorSinPsiDiffEtaSum = new TProfile("f2pCorrelatorSinPsiDiffEtaSum",";(#eta_{1}+#eta_{,2})/2;#LT #LT sin(#psi_{1} - #psi_{2}) #GT #GT",fnBinsEta,fEtaMin,fEtaMax);
 f2pCorrelatorSinPsiDiffEtaSum->SetStats(kFALSE);
 f2pCorrelatorSinPsiSumEtaSum = new TProfile("f2pCorrelatorSinPsiSumEtaSum",";(#eta_{1}+#eta_{,2})/2;#LT #LT sin(#psi_{1} + #psi_{2}) #GT #GT",fnBinsEta,fEtaMin,fEtaMax);
 f2pCorrelatorSinPsiSumEtaSum->SetStats(kFALSE);
 if(fHarmonic == 1) {
   f2pCorrelatorCosPsiDiffEtaSum->SetTitle(Form("#LT#LTcos(#psi_{1}-#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag2[0].Data())); 
   f2pCorrelatorCosPsiSumEtaSum->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag2[0].Data())); 
   f2pCorrelatorSinPsiDiffEtaSum->SetTitle(Form("#LT#LTsin(#psi_{1}-#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag2[0].Data())); 
   f2pCorrelatorSinPsiSumEtaSum->SetTitle(Form("#LT#LTsin(#psi_{1}+#psi_{2})#GT#GT #font[72]{vs} %s",psdTitleFlag2[0].Data())); 
 }
 else {
   f2pCorrelatorCosPsiDiffEtaSum->SetTitle(Form("#LT#LTcos[%d(#psi_{1}-#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag2[0].Data()));
   f2pCorrelatorCosPsiSumEtaSum->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag2[0].Data()));
   f2pCorrelatorSinPsiDiffEtaSum->SetTitle(Form("#LT#LTsin[%d(#psi_{1}-#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag2[0].Data()));
   f2pCorrelatorSinPsiSumEtaSum->SetTitle(Form("#LT#LTsin[%d(#psi_{1}+#psi_{2})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag2[0].Data()));
 }
 fProfileList->Add(f2pCorrelatorCosPsiDiffEtaSum);
 fProfileList->Add(f2pCorrelatorCosPsiSumEtaSum);
 fProfileList->Add(f2pCorrelatorSinPsiDiffEtaSum);
 fProfileList->Add(f2pCorrelatorSinPsiSumEtaSum);

} // end of void AliFlowAnalysisWithMixedHarmonics::BookDifferential()     
   
//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::AccessConstants()
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
 
} // end of void AliFlowAnalysisWithMixedHarmonics::AccessConstants()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::CrossCheckSettings()
{
 // Cross-check if the user settings make sense. 
 
 // ...
  
} // end of void AliFlowAnalysisWithMixedHarmonics::CrossCheckSettings()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::BookAndFillWeightsHistograms()
{
 // Book and fill (by accessing file "weights.root") histograms which hold phi, pt and eta weights.

 if(!fWeightsList)
 {
  cout<<endl;
  cout<<" WARNING (MH): fWeightsList is NULL in BookAndFillWeightsHistograms() !!!!"<<endl;
  cout<<endl;
  exit(0);  
 }
 // Profile to hold flags for weights:   
 TString fUseParticleWeightsName = "fUseParticleWeightsMH";
 fUseParticleWeights = new TProfile(fUseParticleWeightsName.Data(),"0 = particle weight not used, 1 = particle weight used ",3,0,3);
 fUseParticleWeights->SetStats(kFALSE);
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
   if (!fPhiWeights)
   {
     printf("WARNING: no phi weights. bye!\n");
     exit(0);
   }
   if(TMath::Abs(fPhiWeights->GetBinWidth(1)-fPhiBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<" WARNING (MH): Inconsistent binning in histograms for phi-weights throughout the code."<<endl;
    cout<<endl;
    exit(0);
   }
  } else 
    {
     cout<<endl;
     cout<<" WARNING (MH): fWeightsList->FindObject(\"phi_weights\") is NULL in BookAndFillWeightsHistograms() !!!!"<<endl;
     cout<<endl;
     exit(0);
    }
 } // end of if(fUsePhiWeights)
 // Pt-weights:
 if(fUsePtWeights) 
 {
  if(fWeightsList->FindObject("pt_weights"))
  {
   fPtWeights = dynamic_cast<TH1D*>(fWeightsList->FindObject("pt_weights"));
   if (!fPtWeights)
   {
     printf("WARNING: no pt weights. bye!\n");
     exit(0);
   }
   if(TMath::Abs(fPtWeights->GetBinWidth(1)-fPtBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<" WARNING (MH): Inconsistent binning in histograms for pt-weights throughout the code."<<endl;
    cout<<endl;
    exit(0);
   }
  } else 
    {
     cout<<endl;
     cout<<" WARNING (MH): fWeightsList->FindObject(\"pt_weights\") is NULL in BookAndFillWeightsHistograms() !!!!"<<endl;
     cout<<endl;
     exit(0);
    }
 } // end of if(fUsePtWeights)    
 // Eta-weights:
 if(fUseEtaWeights) 
 {
  if(fWeightsList->FindObject("eta_weights"))
  {
   fEtaWeights = dynamic_cast<TH1D*>(fWeightsList->FindObject("eta_weights"));
   if (!fEtaWeights)
   {
     printf("WARNING: no pt weights. bye!\n");
     exit(0);
   }
   if(TMath::Abs(fEtaWeights->GetBinWidth(1)-fEtaBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<" WARNING (MH): Inconsistent binning in histograms for eta-weights throughout the code."<<endl;
    cout<<endl;
    exit(0);
   }
  } else 
    {
     cout<<endl;
     cout<<" WARNING (MH): fUseEtaWeights && fWeightsList->FindObject(\"eta_weights\") is NULL in BookAndFillWeightsHistograms() !!!!"<<endl;
     cout<<endl;
     exit(0);
    }
 } // end of if(fUseEtaWeights)
 
} // end of AliFlowAnalysisWithMixedHarmonics::BookAndFillWeightsHistograms()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::CheckPointersUsedInMake()
{
 // Check pointers used in method Make().
                        
 if(!fReQnk || !fImQnk || !fSpk )
 {                        
  cout<<endl;
  cout<<" WARNING (MH): fReQnk || fImQnk || fSpk is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f3pCorrelatorPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f5pCorrelatorPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f5pCorrelatorPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!fNonIsotropicTermsPro)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): fNonIsotropicTermsPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 
 if(!f3pCorrelatorVsMPro && fCalculateVsM)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorVsMPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f3pPOICorrelatorVsM && fCalculateVsM)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): f3pPOICorrelatorVsM is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fNonIsotropicTermsVsMPro && fCalculateVsM)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): fNonIsotropicTermsVsMPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 
 // Differential correlators:
 if(!fEvaluateDifferential3pCorrelator){return;}
 for(Int_t sd=0;sd<2;sd++)
 {
  if(!(f3pCorrelatorVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorVsEtaSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorVsEtaSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
 }
 for(Int_t sd=0;sd<2;sd++)
 {
  if(!fRePEBE[sd]||!fImPEBE[sd])
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("!fRePEBE[%d]||!fImPEBE[%d]",sd,sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  }
  if(!fReEtaEBE[sd]||!fImEtaEBE[sd])
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("!fReEtaEBE[%d]||!fImEtaEBE[%d]",sd,sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  }
for(Int_t fs=0;fs<2;fs++)
  {
   if(!fOverlapEBE[fs][sd]||!fOverlapEBE[fs][sd])
   {
    cout<<endl;
    cout<<" WARNING (MH): "<<Form("!fOverlapEBE[%d][%d]||!fOverlapEBE[%d][%d]",fs,sd,fs,sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
    cout<<endl;
    exit(0);   
   }
   if(!fOverlapEBE2[fs][sd]||!fOverlapEBE2[fs][sd])
   {
    cout<<endl;
    cout<<" WARNING (MH): "<<Form("!fOverlapEBE2[%d][%d]||!fOverlapEBE2[%d][%d]",fs,sd,fs,sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
    cout<<endl;
    exit(0);   
   }
  } // end of for(Int_t fs=0;fs<2;fs++)
 } // end of for(Int_t sd=0;sd<2;sd++)
 
 //2p correlator vs |Pt1-Pt2|
 if(!f2pCorrelatorCosPsiDiffPtDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffPtDiff is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCosPsiSumPtDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiSumPtDiff is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiDiffPtDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiDiffPtDiff is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiSumPtDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiSumPtDiff is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

 //2p correlator vs (Pt1+Pt2)/2
 if(!f2pCorrelatorCosPsiDiffPtSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffPtSum is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCosPsiSumPtSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiSumPtSum is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiDiffPtSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiDiffPtSum is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiSumPtSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiSumPtSum is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

 //2p correlator vs |eta1-eta2|
 if(!f2pCorrelatorCosPsiDiffEtaDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffEtaDiff is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCosPsiSumEtaDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiSumEtaDiff is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiDiffEtaDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiDiffEtaDiff is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiSumEtaDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiSumEtaDiff is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

 //2p correlator vs (eta1+eta2)/2
 if(!f2pCorrelatorCosPsiDiffEtaSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffEtaSum is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCosPsiSumEtaSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiSumEtaSum is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiDiffEtaSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiDiffEtaSum is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiSumEtaSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiSumEtaSum is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

} // end of AliFlowAnalysisWithMixedHarmonics::CheckPointersUsedInMake()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::CheckPointersUsedInFinish()
{
 // Check pointers used in method Finish().
 
 if(!fAnalysisSettings)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): fAnalysisSettings is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f3pCorrelatorPro)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fNonIsotropicTermsPro)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): fNonIsotropicTermsPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 
 if(!f3pPOICorrelatorVsM && fCalculateVsM)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): f3pPOICorrelatorVsM is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f3pCorrelatorVsMPro && fCalculateVsM)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorVsMPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f3pCorrelatorVsMHist && fCalculateVsM)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorVsMHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fNonIsotropicTermsVsMPro && fCalculateVsM)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): fNonIsotropicTermsVsMPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

 if(!f3pCorrelatorHist)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }   
 if(!fDetectorBiasHist)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): fDetectorBiasHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }   
 /* to be improved - enabled eventually
 if(!fDetectorBiasVsMHist && fCalculateVsM)
 {                        
  cout<<endl;
  cout<<" WARNING (MH): !fDetectorBiasVsMHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 
 */  
 
 // Differential correlators:
 if(!fEvaluateDifferential3pCorrelator){return;} 
 for(Int_t sd=0;sd<2;sd++)
 {
  if(!(f3pCorrelatorVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorVsEtaSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorVsEtaSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
 }

 //2p correlator vs |Pt1-Pt2|
 if(!f2pCorrelatorCosPsiDiffPtDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffPtDiff is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCosPsiSumPtDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiSumPtDiff is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiDiffPtDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiDiffPtDiff is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiSumPtDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiSumPtDiff is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

 //2p correlator vs (Pt1+Pt2)/2
 if(!f2pCorrelatorCosPsiDiffPtSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffPtSum is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCosPsiSumPtSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiSumPtSum is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiDiffPtSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiDiffPtSum is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiSumPtSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiSumPtSum is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

 //2p correlator vs |eta1-eta2|
 if(!f2pCorrelatorCosPsiDiffEtaDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffEtaDiff is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCosPsiSumEtaDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiSumEtaDiff is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiDiffEtaDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiDiffEtaDiff is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiSumEtaDiff) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiSumEtaDiff is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

 //2p correlator vs (eta1+eta2)/2
 if(!f2pCorrelatorCosPsiDiffEtaSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffEtaSum is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCosPsiSumEtaSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiSumEtaSum is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiDiffEtaSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiDiffEtaSum is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorSinPsiSumEtaSum) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorSinPsiSumEtaSum is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

} // end of AliFlowAnalysisWithMixedHarmonics::CheckPointersUsedInFinish()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::PrintOnTheScreen()
{
 // Print the final results on the screen.
 
 // <<cos[n*(phi1+phi2-2phi3)]>>:
 Double_t d3pCorrelator = 0.;
 Double_t d3pCorrelatorError = 0.;
 if(!fCorrectForDetectorEffects)
 {
  d3pCorrelator = f3pCorrelatorPro->GetBinContent(1);
  d3pCorrelatorError = f3pCorrelatorPro->GetBinError(1);
 } else
   {
    d3pCorrelator = f3pCorrelatorHist->GetBinContent(1);
    d3pCorrelatorError = f3pCorrelatorHist->GetBinError(1); 
   }
 
 // <<cos[n*(psi1+psi2-2phi3)]>>:
 Double_t d3pCorrelatorPoi = 0.;
 Double_t d3pCorrelatorPoiError = 0.;
 if(fEvaluateDifferential3pCorrelator)
 {
  GetCorrelatorAndError(f3pCorrelatorVsPtSumDiffPro[0],
		       d3pCorrelatorPoi,
		       d3pCorrelatorPoiError);
 }		       
 cout<<endl;
 cout<<"*******************************************************"<<endl;
 cout<<"*******************************************************"<<endl;
 cout<<"                    Mixed Harmonics                      "<<endl; 
 cout<<endl;
 if(fHarmonic!=1)
 {
  cout<<"  cos["<<fHarmonic<<"(phi1+phi2-2phi3)] = "<<d3pCorrelator<<" +/- "<<d3pCorrelatorError<<endl;
  cout<<"  cos["<<fHarmonic<<"(psi1+psi2-2phi3)] = "<<d3pCorrelatorPoi<<" +/- "<<d3pCorrelatorPoiError<<endl;
 } else
   {
    cout<<"  cos(phi1+phi2-2phi3) = "<<d3pCorrelator<<" +/- "<<d3pCorrelatorError<<endl;
    cout<<"  cos(psi1+psi2-2phi3) = "<<d3pCorrelatorPoi<<" +/- "<<d3pCorrelatorPoiError<<endl;
   }
 if(!fCorrectForDetectorEffects)
 {
  cout<<"  Detector Bias = "<<fDetectorBiasHist->GetBinContent(1)<<" (not corrected for)"<<endl;
 } else
   {
    cout<<"  Detector Bias = "<<fDetectorBiasHist->GetBinContent(1)<<" (corrected for)"<<endl; 
   }
 cout<<endl;
 cout<<"             nEvts = "<<(Int_t)fCommonHists->GetHistMultRP()->GetEntries()<<", <M> = "<<(Double_t)fCommonHists->GetHistMultRP()->GetMean()<<endl; 
 cout<<"*******************************************************"<<endl;
 cout<<"*******************************************************"<<endl;

} // end of void AliFlowAnalysisWithMixedHarmonics::PrintOnTheScreen()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::AccessSettings()
{
 // Access the settings for analysis with mixed harmonics.
 
 fCorrectForDetectorEffects = (Bool_t)fAnalysisSettings->GetBinContent(1);
 fNoOfMultipicityBins = (Int_t)fAnalysisSettings->GetBinContent(2);
 fMultipicityBinWidth = (Double_t)fAnalysisSettings->GetBinContent(3);
 fMinMultiplicity = (Double_t)fAnalysisSettings->GetBinContent(4);
 fPrintOnTheScreen = (Bool_t)fAnalysisSettings->GetBinContent(5);
 fHarmonic = (Int_t)fAnalysisSettings->GetBinContent(6);
 fOppositeChargesPOI = (Bool_t)fAnalysisSettings->GetBinContent(7);      
 fEvaluateDifferential3pCorrelator = (Bool_t)fAnalysisSettings->GetBinContent(8);    
 fCalculateVsM = (Bool_t)fAnalysisSettings->GetBinContent(9);  
 fShowBinLabelsVsM = (Bool_t)fAnalysisSettings->GetBinContent(10); 
                                                                                                                                                                                                                                                                                                                                                                                      
} // end of AliFlowAnalysisWithMixedHarmonics::AccessSettings()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::CorrectForDetectorEffects()
{
 // Correct measured 3-p correlator cos[n(phi1+phi2-2phi3)] for detector effects.
  
 Double_t measured3pCorrelator = f3pCorrelatorPro->GetBinContent(1); // biased by detector effects
 Double_t corrected3pCorrelator = 0.; // corrected for detector effects
 Double_t nonIsotropicTerms[10] = {0.}; // there are 10 distinct non-isotropic terms
 for(Int_t nit=0;nit<10;nit++)
 {
  nonIsotropicTerms[nit] = fNonIsotropicTermsPro->GetBinContent(nit+1);
 }                    
 // Calculate corrected 3-p correlator:                     
 corrected3pCorrelator = measured3pCorrelator
                       - nonIsotropicTerms[2]*nonIsotropicTerms[4]                                                                                
                       - nonIsotropicTerms[3]*nonIsotropicTerms[5]                                                              
                       - 2.*nonIsotropicTerms[0]*nonIsotropicTerms[6]                                       
                       - 2.*nonIsotropicTerms[1]*nonIsotropicTerms[7]                                       
                       + 2.*nonIsotropicTerms[2]*(pow(nonIsotropicTerms[0],2.)-pow(nonIsotropicTerms[1],2.))                                       
                       + 4.*nonIsotropicTerms[3]*nonIsotropicTerms[0]*nonIsotropicTerms[1]; 
 // Store corrected correlator:
 if(fCorrectForDetectorEffects)
 {
  f3pCorrelatorHist->SetBinContent(1,corrected3pCorrelator);
  f3pCorrelatorHist->SetBinError(1,f3pCorrelatorPro->GetBinError(1)); // to be improved (propagate error for non-isotropic terms)
 }
 // Quantify bias from detector inefficiences to 3-p correlator. Remark: Bias is quantified as a 
 // ratio between corrected and measured 3-p correlator:
 //              bias = corrected/measured
 // This bias is stored in histogram fDetectorBias.
 Double_t bias = 0.;
 if(TMath::Abs(measured3pCorrelator)>1.e-44)
 {
  bias = corrected3pCorrelator/measured3pCorrelator;
  fDetectorBiasHist->SetBinContent(1,bias);                                                          
 }   
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithMixedHarmonics::CorrectForDetectorEffects()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::CorrectForDetectorEffectsVsM()
{
 // Correct measured 3-p correlator cos[n(phi1+phi2-2phi3)] vs M for detector effects.
  
 for(Int_t b=1;b<=fNoOfMultipicityBins+2;b++) 
 {  
  Double_t measured3pCorrelator = f3pCorrelatorVsMPro->GetBinContent(b); // biased by detector effects
  Double_t corrected3pCorrelator = 0.; // corrected for detector effects
  Double_t nonIsotropicTerms[10] = {0.}; // there are 10 distinct non-isotropic terms
  for(Int_t nit=0;nit<10;nit++)
  {
   nonIsotropicTerms[nit] = fNonIsotropicTermsVsMPro->GetBinContent(fNonIsotropicTermsVsMPro->GetBin(nit+1,b));
  }                    
  // Calculate corrected 3-p correlator:                     
  corrected3pCorrelator = measured3pCorrelator
                        - nonIsotropicTerms[2]*nonIsotropicTerms[4]                                                                                
                        - nonIsotropicTerms[3]*nonIsotropicTerms[5]                                                              
                        - 2.*nonIsotropicTerms[0]*nonIsotropicTerms[6]                                       
                        - 2.*nonIsotropicTerms[1]*nonIsotropicTerms[7]                                       
                        + 2.*nonIsotropicTerms[2]*(pow(nonIsotropicTerms[0],2.)-pow(nonIsotropicTerms[1],2.))                                       
                        + 4.*nonIsotropicTerms[3]*nonIsotropicTerms[0]*nonIsotropicTerms[1]; 
  // Store corrected correlator:
  if(fCorrectForDetectorEffects)
  {
   f3pCorrelatorVsMHist->SetBinContent(b,corrected3pCorrelator);
   f3pCorrelatorVsMHist->SetBinError(b,f3pCorrelatorVsMPro->GetBinError(b)); // to be improved (propagate error for non-isotropic terms)
  }
  // Quantify bias from detector inefficiences to 3-p correlator. Remark: Bias is quantified as a 
  // ratio between corrected and measured 3-p correlator:
  //              bias = corrected/measured
  // This bias is stored in histogram fDetectorBias.
  Double_t bias = 0.;
  if(measured3pCorrelator)
  {
   bias = corrected3pCorrelator/measured3pCorrelator;
   fDetectorBiasVsMHist->SetBinContent(b,bias);                                                          
  }   
 } // end of for(Int_t b=1;b<=fNoOfMultipicityBins;b++) 
                                                                                                                                                                                                                                                                                                                                   
} // end of AliFlowAnalysisWithMixedHarmonics::CorrectForDetectorEffectsVsM()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::ResetEventByEventQuantities()
{
 // Reset all event-by-event quantities.
 
 fReQnk->Zero();
 fImQnk->Zero();
 fSpk->Zero();
  
 if(!fEvaluateDifferential3pCorrelator){return;}
 for(Int_t sd=0;sd<2;sd++)
 {
  fRePEBE[sd]->Reset();
  fImPEBE[sd]->Reset();
  fReEtaEBE[sd]->Reset();
  fImEtaEBE[sd]->Reset();
 }
 for(Int_t fs=0;fs<2;fs++)
 {
  for(Int_t sd=0;sd<2;sd++)
  {
   fOverlapEBE[fs][sd]->Reset();
   fOverlapEBE2[fs][sd]->Reset();
  } 
 }
 
} // end of void AliFlowAnalysisWithMixedHarmonics::ResetEventByEventQuantities()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::Calculate3pCorrelator()
{
 // Calculate 3-p azimuthal correlator cos[n(phi1+phi2-2phi3)] in terms of Q_{n,k} and S_{p,k}.
 
 // a) Calculate 3-p correlator without using particle weights;
 // b) Calculate 3-p correlator with using particle weights.

 // a) Calculate 3-p correlator without using particle weights: 
 if(!(fUsePhiWeights || fUsePtWeights || fUseEtaWeights))
 {
  // Multiplicity (number of RPs):
  Double_t dMult = (*fSpk)(0,0);
  // Real and imaginary parts of non-weighted Q-vectors (Q_{n,0}) evaluated in harmonics n and 2n: 
  Double_t dReQ1n = (*fReQnk)(0,0);
  Double_t dReQ2n = (*fReQnk)(1,0);
  Double_t dImQ1n = (*fImQnk)(0,0);
  Double_t dImQ2n = (*fImQnk)(1,0);
  // 3-particle azimuthal correlator <cos(n*(phi1+phi2-2phi3))>:
  Double_t three1n1n2n = (pow(dReQ1n,2.)*dReQ2n + 2.*dReQ1n*dImQ1n*dImQ2n - pow(dImQ1n,2.)*dReQ2n
                       - 2.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))
                       - (pow(dReQ2n,2.)+pow(dImQ2n,2.))+2.*dMult)
                       / (dMult*(dMult-1.)*(dMult-2.));                 
  
  // Fill all-events profile:                     
  f3pCorrelatorPro->Fill(0.5,three1n1n2n,dMult*(dMult-1.)*(dMult-2.));
  
  // 3-particle azimuthal correlator <cos(n*(phi1+phi2-2phi3))> vs multiplicity:
  if(fCalculateVsM)
  {
   if(dMult<fMinMultiplicity) 
   {
    f3pCorrelatorVsMPro->Fill(0.5,three1n1n2n,dMult*(dMult-1.)*(dMult-2.));
   } else if(dMult>=fMinMultiplicity+fNoOfMultipicityBins*fMultipicityBinWidth)
     {
      f3pCorrelatorVsMPro->Fill(0.5+fNoOfMultipicityBins+1,three1n1n2n,dMult*(dMult-1.)*(dMult-2.));  
     } else
       {
        f3pCorrelatorVsMPro->Fill(1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),three1n1n2n,dMult*(dMult-1.)*(dMult-2.));
       }
  } // end of if(fCalculateVsM)
      
 } // end of if(!(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)) 

 // b) Calculate 3-p correlator with using particle weights: 
 if(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)
 {
  // ...
 } // end of if(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)
 
} // end of void AliFlowAnalysisWithMixedHarmonics::Calculate3pCorrelator() 

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::Calculate5pCorrelator()
{
 // Calculate 5-p azimuthal correlator cos[n(2phi1+2phi2+2phi3-3phi4-3phi5)] in terms of Q_{n,k} and S_{p,k}.
 
 // a) Calculate 5-p correlator without using particle weights;
 // b) Calculate 5-p correlator with using particle weights.

 // a) Calculate 5-p correlator without using particle weights: 
 if(!(fUsePhiWeights || fUsePtWeights || fUseEtaWeights))
 {
  // Multiplicity (number of RPs):
  Double_t dMult = (*fSpk)(0,0);
  // Real and imaginary parts of non-weighted Q-vectors (Q_{n,0}) evaluated in harmonics n,2n,...,6n: 
  Double_t dReQ1n = (*fReQnk)(0,0);
  Double_t dReQ2n = (*fReQnk)(1,0);
  Double_t dReQ3n = (*fReQnk)(2,0);
  Double_t dReQ4n = (*fReQnk)(3,0);
  //Double_t dReQ5n = (*fReQnk)(4,0); // not needed
  Double_t dReQ6n = (*fReQnk)(5,0);
  Double_t dImQ1n = (*fImQnk)(0,0);
  Double_t dImQ2n = (*fImQnk)(1,0);
  Double_t dImQ3n = (*fImQnk)(2,0);
  Double_t dImQ4n = (*fImQnk)(3,0);
  //Double_t dImQ5n = (*fImQnk)(4,0); // not needed
  Double_t dImQ6n = (*fImQnk)(5,0);
  
  // 5-particle azimuthal correlator:
  Double_t five2n2n2n3n3n = 0.; // <cos[n(2phi1+2phi2+2phi3-3phi4-3phi5)]>
  Double_t reQ2nQ2nQ2nQ3nstarQ3nstar = pow(dReQ2n,3.)*pow(dReQ3n,2.) 
                                     - 3.*dReQ2n*pow(dReQ3n,2.)*pow(dImQ2n,2.)
                                     + 6.*pow(dReQ2n,2.)*dReQ3n*dImQ2n*dImQ3n 
                                     - 2.*dReQ3n*pow(dImQ2n,3.)*dImQ3n-pow(dReQ2n,3.)*pow(dImQ3n,2.) 
                                     + 3.*dReQ2n*pow(dImQ2n,2.)*pow(dImQ3n,2.);
  Double_t reQ2nQ2nQ2nQ6nstar = dReQ6n*pow(dReQ2n,3)-3.*dReQ2n*dReQ6n*pow(dImQ2n,2)
                               + 3.*dImQ2n*dImQ6n*pow(dReQ2n,2)-dImQ6n*pow(dImQ2n,3); 
  Double_t reQ4nQ2nQ3nstarQ3nstar = (dReQ4n*dReQ2n-dImQ4n*dImQ2n)*(dReQ3n*dReQ3n-dImQ3n*dImQ3n)
                                  + 2.*(dReQ4n*dImQ2n+dImQ4n*dReQ2n)*dReQ3n*dImQ3n;
  Double_t reQ2nQ2nQ1nstarQ3nstar = (pow(dReQ2n,2.)-pow(dImQ2n,2.))*(dReQ3n*dReQ1n-dImQ3n*dImQ1n) 
                                  + 2.*dReQ2n*dImQ2n*(dReQ3n*dImQ1n+dImQ3n*dReQ1n);                            
  Double_t reQ6nQ3nstarQ3nstar = pow(dReQ3n,2.)*dReQ6n + 2.*dReQ3n*dImQ3n*dImQ6n 
                               - pow(dImQ3n,2.)*dReQ6n; 
  Double_t reQ4nQ2nQ6nstar = dReQ6n*dReQ4n*dReQ2n-dReQ6n*dImQ4n*dImQ2n+dImQ6n*dReQ4n*dImQ2n
                           + dImQ6n*dImQ4n*dReQ2n;
  Double_t reQ4nQ1nstarQ3nstar = dReQ4n*(dReQ3n*dReQ1n-dImQ3n*dImQ1n)+dImQ4n*(dReQ3n*dImQ1n+dImQ3n*dReQ1n); 
  Double_t reQ2nQ2nQ4nstar = pow(dReQ2n,2.)*dReQ4n+2.*dReQ2n*dImQ2n*dImQ4n-pow(dImQ2n,2.)*dReQ4n;                       
  Double_t reQ2nQ1nQ3nstar = dReQ3n*dReQ2n*dReQ1n-dReQ3n*dImQ2n*dImQ1n+dImQ3n*dReQ2n*dImQ1n
                           + dImQ3n*dImQ2n*dReQ1n;
  Double_t reQ2nQ1nstarQ1nstar = pow(dReQ1n,2.)*dReQ2n + 2.*dReQ1n*dImQ1n*dImQ2n - pow(dImQ1n,2.)*dReQ2n; 
  // Analytic expression for 5-particle azimuthal correlator:                                   
  five2n2n2n3n3n = (reQ2nQ2nQ2nQ3nstarQ3nstar-reQ2nQ2nQ2nQ6nstar-3.*reQ4nQ2nQ3nstarQ3nstar 
                 - 6.*reQ2nQ2nQ1nstarQ3nstar+2.*reQ6nQ3nstarQ3nstar+3.*reQ4nQ2nQ6nstar
                 + 6.*reQ4nQ1nstarQ3nstar+6.*reQ2nQ2nQ4nstar
                 + 12.*reQ2nQ1nQ3nstar+6.*reQ2nQ1nstarQ1nstar
                 - 2.*((pow(dReQ6n,2.)+pow(dImQ6n,2.)) 
                 + 3.*(pow(dReQ4n,2.)+pow(dImQ4n,2.))
                 + 6.*(pow(dReQ3n,2.)+pow(dImQ3n,2.)) 
                 + 9.*(pow(dReQ2n,2.)+pow(dImQ2n,2.))
                 + 6.*(pow(dReQ1n,2.)+pow(dImQ1n,2.))-12.*dMult))
                 /(dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
  // Fill all-events profile:                       
  f5pCorrelatorPro->Fill(0.5,five2n2n2n3n3n,dMult*(dMult-1.)*(dMult-2.)*(dMult-3.)*(dMult-4.));
 } // end of if(!(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)) 

 // b) Calculate 5-p correlator with using particle weights: 
 if(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)
 {
  // ...
 } // end of if(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)
 
} // end of void AliFlowAnalysisWithMixedHarmonics::Calculate5pCorrelator() 

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::CalculateNonIsotropicTerms()
{
 // Calculate non-isotropic terms which appear in the decomposition of 3-p correlator <cos[n(phi1+phi2-2phi3)]>.
 
 // a) Calculate without using particle weights;
 // b) Calculate using particle weights.
 
 // For detector with uniform acceptance all these terms vanish. These non-isotropic terms are stored in fNonIsotropicTermsPro.
 // Binning of fNonIsotropicTermsPro is organized as follows:
 //  1st bin: <<cos(n*phi1)>>
 //  2nd bin: <<sin(n*phi1)>>
 //  3rd bin: <<cos(2n*phi1)>>
 //  4th bin: <<sin(2n*phi1)>>
 //  5th bin: <<cos(n*(phi1+phi2)>>
 //  6th bin: <<sin(n*(phi1+phi2)>>
 //  7th bin: <<cos(n*(2phi1-phi2)>>
 //  8th bin: <<sin(n*(2phi1-phi2)>>
 //  9th bin: <<cos(n*(phi1-phi2-phi3)>> // not needed
 // 10th bin: <<sin(n*(phi1-phi2-phi3)>> // not needed 
 
 // a) Calculate without using particle weights:
 if(!(fUsePhiWeights || fUsePtWeights || fUseEtaWeights))
 {
  // Multiplicity (number of RPs):
  Double_t dMult = (*fSpk)(0,0);
  // Real and imaginary parts of non-weighted Q-vectors (Q_{n,0}) evaluated in harmonics n and 2n: 
  Double_t dReQ1n = (*fReQnk)(0,0);
  Double_t dReQ2n = (*fReQnk)(1,0);
  Double_t dImQ1n = (*fImQnk)(0,0);
  Double_t dImQ2n = (*fImQnk)(1,0);
  // 1-particle terms:
  Double_t cosP1n = 0.; // <cos(n*(phi1))>  
  Double_t sinP1n = 0.; // <sin(n*(phi1))>
  Double_t cosP2n = 0.; // <cos(2n*(phi1))>  
  Double_t sinP2n = 0.; // <sin(2n*(phi1))>
  if(dMult>0)
  { 
   cosP1n = dReQ1n/dMult; 
   sinP1n = dImQ1n/dMult;
   cosP2n = dReQ2n/dMult; 
   sinP2n = dImQ2n/dMult;   
   // All-event avarages:
   fNonIsotropicTermsPro->Fill(0.5,cosP1n,dMult); // <<cos(n*(phi1))>> 
   fNonIsotropicTermsPro->Fill(1.5,sinP1n,dMult); // <<sin(n*(phi1))>>   
   fNonIsotropicTermsPro->Fill(2.5,cosP2n,dMult); // <<cos(2n*(phi1))>> 
   fNonIsotropicTermsPro->Fill(3.5,sinP2n,dMult); // <<sin(2n*(phi1))>>   
   // All-event avarages vs M:
   if(fCalculateVsM)
   {
    if(dMult<fMinMultiplicity) 
    {
     fNonIsotropicTermsVsMPro->Fill(0.5,0.5,cosP1n,dMult); // <<cos(n*(phi1))>> 
     fNonIsotropicTermsVsMPro->Fill(1.5,0.5,sinP1n,dMult); // <<sin(n*(phi1))>>   
     fNonIsotropicTermsVsMPro->Fill(2.5,0.5,cosP2n,dMult); // <<cos(2n*(phi1))>> 
     fNonIsotropicTermsVsMPro->Fill(3.5,0.5,sinP2n,dMult); // <<sin(2n*(phi1))>>   
    } else if(dMult>=fMinMultiplicity+fNoOfMultipicityBins*fMultipicityBinWidth)
      {
       fNonIsotropicTermsVsMPro->Fill(0.5,0.5+fNoOfMultipicityBins+1,cosP1n,dMult); // <<cos(n*(phi1))>> 
       fNonIsotropicTermsVsMPro->Fill(1.5,0.5+fNoOfMultipicityBins+1,sinP1n,dMult); // <<sin(n*(phi1))>>   
       fNonIsotropicTermsVsMPro->Fill(2.5,0.5+fNoOfMultipicityBins+1,cosP2n,dMult); // <<cos(2n*(phi1))>> 
       fNonIsotropicTermsVsMPro->Fill(3.5,0.5+fNoOfMultipicityBins+1,sinP2n,dMult); // <<sin(2n*(phi1))>>       
      } else
        {
         fNonIsotropicTermsVsMPro->Fill(0.5,1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),cosP1n,dMult); // <<cos(n*(phi1))>> 
         fNonIsotropicTermsVsMPro->Fill(1.5,1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),sinP1n,dMult); // <<sin(n*(phi1))>>   
         fNonIsotropicTermsVsMPro->Fill(2.5,1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),cosP2n,dMult); // <<cos(2n*(phi1))>> 
         fNonIsotropicTermsVsMPro->Fill(3.5,1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),sinP2n,dMult); // <<sin(2n*(phi1))>>    
        }
   } // end of if(fCalculateVsM)         
  } // end of if(dMult>0) 
  // 2-particle terms:
  Double_t cosP1nP1n = 0.; // <cos(n*(phi1+phi2))>
  Double_t sinP1nP1n = 0.; // <sin(n*(phi1+phi2))>
  Double_t cosP2nM1n = 0.; // <cos(n*(2phi1-phi2))>
  Double_t sinP2nM1n = 0.; // <sin(n*(2phi1-phi2))>
  if(dMult>1)
  {
   cosP1nP1n = (pow(dReQ1n,2)-pow(dImQ1n,2)-dReQ2n)/(dMult*(dMult-1)); 
   sinP1nP1n = (2.*dReQ1n*dImQ1n-dImQ2n)/(dMult*(dMult-1)); 
   cosP2nM1n = (dReQ2n*dReQ1n+dImQ2n*dImQ1n-dReQ1n)/(dMult*(dMult-1)); 
   sinP2nM1n = (dImQ2n*dReQ1n-dReQ2n*dImQ1n-dImQ1n)/(dMult*(dMult-1)); 
   // All-event avarages:
   fNonIsotropicTermsPro->Fill(4.5,cosP1nP1n,dMult*(dMult-1.)); // <<cos(n*(phi1+phi2))>> 
   fNonIsotropicTermsPro->Fill(5.5,sinP1nP1n,dMult*(dMult-1.)); // <<sin(n*(phi1+phi2))>>   
   fNonIsotropicTermsPro->Fill(6.5,cosP2nM1n,dMult*(dMult-1.)); // <<cos(n*(2phi1-phi2))>> 
   fNonIsotropicTermsPro->Fill(7.5,sinP2nM1n,dMult*(dMult-1.)); // <<sin(n*(2phi1-phi2))>>   
   // All-event avarages vs M:  
   if(fCalculateVsM)
   {
    if(dMult<fMinMultiplicity) 
    {
     fNonIsotropicTermsVsMPro->Fill(4.5,0.5,cosP1nP1n,dMult*(dMult-1.)); // <<cos(n*(phi1+phi2))>> 
     fNonIsotropicTermsVsMPro->Fill(5.5,0.5,sinP1nP1n,dMult*(dMult-1.)); // <<sin(n*(phi1+phi2))>>   
     fNonIsotropicTermsVsMPro->Fill(6.5,0.5,cosP2nM1n,dMult*(dMult-1.)); // <<cos(n*(2phi1-phi2))>> 
     fNonIsotropicTermsVsMPro->Fill(7.5,0.5,sinP2nM1n,dMult*(dMult-1.)); // <<sin(n*(2phi1-phi2))>>   
    } else if(dMult>=fMinMultiplicity+fNoOfMultipicityBins*fMultipicityBinWidth)
      {
       fNonIsotropicTermsVsMPro->Fill(4.5,0.5+fNoOfMultipicityBins+1,cosP1nP1n,dMult*(dMult-1.)); // <<cos(n*(phi1+phi2))>> 
       fNonIsotropicTermsVsMPro->Fill(5.5,0.5+fNoOfMultipicityBins+1,sinP1nP1n,dMult*(dMult-1.)); // <<sin(n*(phi1+phi2))>>   
       fNonIsotropicTermsVsMPro->Fill(6.5,0.5+fNoOfMultipicityBins+1,cosP2nM1n,dMult*(dMult-1.)); // <<cos(n*(2phi1-phi2))>> 
       fNonIsotropicTermsVsMPro->Fill(7.5,0.5+fNoOfMultipicityBins+1,sinP2nM1n,dMult*(dMult-1.)); // <<sin(n*(2phi1-phi2))>>    
      } else
        {
         fNonIsotropicTermsVsMPro->Fill(4.5,1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),cosP1nP1n,dMult*(dMult-1.)); // <<cos(n*(phi1+phi2))>> 
         fNonIsotropicTermsVsMPro->Fill(5.5,1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),sinP1nP1n,dMult*(dMult-1.)); // <<sin(n*(phi1+phi2))>>   
         fNonIsotropicTermsVsMPro->Fill(6.5,1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),cosP2nM1n,dMult*(dMult-1.)); // <<cos(n*(2phi1-phi2))>> 
         fNonIsotropicTermsVsMPro->Fill(7.5,1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),sinP2nM1n,dMult*(dMult-1.)); // <<sin(n*(2phi1-phi2))>>    
        }
   } // end of if(fCalculateVsM)       
  } // end of if(dMult>1) 
  // 3-particle: correct and ready but not needed, hence commented out.
  /*
  Double_t cosP1nM1nM1n = 0.; // <cos(n*(phi1-phi2-phi3))>
  Double_t sinP1nM1nM1n = 0.; // <sin(n*(phi1-phi2-phi3))>
  if(dMult>2)
  {
   cosP1nM1nM1n = (dReQ1n*(pow(dReQ1n,2)+pow(dImQ1n,2))-dReQ1n*dReQ2n-dImQ1n*dImQ2n-2.*(dMult-1)*dReQ1n)
                / (dMult*(dMult-1)*(dMult-2)); 
   sinP1nM1nM1n = (-dImQ1n*(pow(dReQ1n,2)+pow(dImQ1n,2))+dReQ1n*dImQ2n-dImQ1n*dReQ2n+2.*(dMult-1)*dImQ1n)
                / (dMult*(dMult-1)*(dMult-2));              
   // All-events avarages:
   fNonIsotropicTermsPro->Fill(8.5,cosP1nM1nM1n,dMult*(dMult-1.)*(dMult-2.)); // <<cos(n*(phi1-phi2-phi3))>> 
   fNonIsotropicTermsPro->Fill(9.5,sinP1nM1nM1n,dMult*(dMult-1.)*(dMult-2.)); // <<sin(n*(phi1-phi2-phi3))>>    
   // All-events avarages vs M:  
   if(fCalculateVsM)
   {
    if(dMult<fMinMultiplicity) 
    {
     fNonIsotropicTermsVsMPro->Fill(8.5,0.5,cosP1nM1nM1n,dMult*(dMult-1.)*(dMult-2.)); // <<cos(n*(phi1-phi2-phi3))>> 
     fNonIsotropicTermsVsMPro->Fill(9.5,0.5,sinP1nM1nM1n,dMult*(dMult-1.)*(dMult-2.)); // <<sin(n*(phi1-phi2-phi3))>>    
    } else if(dMult>=fMinMultiplicity+fNoOfMultipicityBins*fMultipicityBinWidth)
      {
       fNonIsotropicTermsVsMPro->Fill(8.5,0.5+fNoOfMultipicityBins+1,cosP1nM1nM1n,dMult*(dMult-1.)*(dMult-2.)); // <<cos(n*(phi1-phi2-phi3))>> 
       fNonIsotropicTermsVsMPro->Fill(9.5,0.5+fNoOfMultipicityBins+1,sinP1nM1nM1n,dMult*(dMult-1.)*(dMult-2.)); // <<sin(n*(phi1-phi2-phi3))>>    
      } else
        {
         fNonIsotropicTermsVsMPro->Fill(8.5,1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),
                                        cosP1nM1nM1n,dMult*(dMult-1.)*(dMult-2.)); // <<cos(n*(phi1-phi2-phi3))>> 
         fNonIsotropicTermsVsMPro->Fill(9.5,1.5+(Int_t)((dMult-fMinMultiplicity)/fMultipicityBinWidth),
                                        sinP1nM1nM1n,dMult*(dMult-1.)*(dMult-2.)); // <<sin(n*(phi1-phi2-phi3))>>           
        }  
   } // end of if(fCalculateVsM)
  } // end of if(dMult>2)
  */
 } // end of if(!(fUsePhiWeights || fUsePtWeights || fUseEtaWeights))
 
 // b) Calculate using particle weights:
 if(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)
 {
  // ...
 } // end of if(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)
 
} // end of void AliFlowAnalysisWithMixedHarmonics::CalculateNonIsotropicTerms()

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::CalculateDifferential3pCorrelator(Double_t &gIntegratedValue)
{
 // Calculate differential 3-p azimuthal correlator cos[n(psi1+psi2-2phi3)] in terms of Q_{2n}, p_{n}, q1_{n} and q2_{n}.
 
 // a) Calculate differential 3-p correlator without using particle weights;
 // b) Calculate differential 3-p correlator with using particle weights.

 // a) Calculate differential 3-p correlator without using particle weights: 
 if(!(fUsePhiWeights || fUsePtWeights || fUseEtaWeights))
 {
   Int_t iBinCounter = 0;
   Double_t gSumBinContentTimesWeight = 0., gSumWeight = 0.;
   Double_t gSumBinContentTimesWeightSquared = 0.;

  // Multiplicity (number of RPs):
  Double_t dMult = (*fSpk)(0,0);
  // Real and imaginary parts of non-weighted Q-vectors (Q_{n,0}) evaluated in harmonic 2n: 
  Double_t dReQ2n = (*fReQnk)(1,0);
  Double_t dImQ2n = (*fImQnk)(1,0);
  for(Int_t sd=0;sd<2;sd++) 
  {
    // [(p1+p2)/2,|p1-p2|]
   // looping over all bins and calculating reduced correlations: 
   for(Int_t b=1;b<=fnBinsPt;b++)
   {
    // real and imaginary parts of p_{n}: 
    Double_t p1nRe = fRePEBE[sd]->GetBinContent(b)*fRePEBE[sd]->GetBinEntries(b);
    Double_t p1nIm = fImPEBE[sd]->GetBinContent(b)*fImPEBE[sd]->GetBinEntries(b);
    // overlap 1: to be improved (terminology)
    Double_t overlap1 = fOverlapEBE[0][sd]->GetBinContent(b)*fOverlapEBE[0][sd]->GetBinEntries(b);
    // overlap 2: to be improved (terminology)
    Double_t overlap2 = fOverlapEBE[1][sd]->GetBinContent(b)*fOverlapEBE[1][sd]->GetBinEntries(b);    
    // number of pairs of POIs in particular (p1+p2)/2 or |p1-p2| bin:
    Double_t mp = fRePEBE[sd]->GetBinEntries(b);
    // number of pairs of POI1/RP and POI2 in particular (p1+p2)/2 or |p1-p2| bin:
    Double_t mOverlap1 = fOverlapEBE[0][sd]->GetBinEntries(b);
    // number of pairs of POI2/RP and POI1 in particular (p1+p2)/2 or |p1-p2| bin:
    Double_t mOverlap2 = fOverlapEBE[1][sd]->GetBinEntries(b);
    // e-b-e weight for cos[n(psi1+psi2-2phi3)]:
    Double_t weight = mp*dMult-mOverlap1-mOverlap2;  
    
    Double_t cosP2nphi1M1npsi2M1npsi2 = 0; // cos[n(psi1+psi2-2phi3)]
    if(weight>0.)
    {
     cosP2nphi1M1npsi2M1npsi2 = (p1nRe*dReQ2n+p1nIm*dImQ2n-overlap1-overlap2)/(weight);
    }
    f3pCorrelatorVsPtSumDiffPro[sd]->Fill(fPtMin+(b-1)*fPtBinWidth,cosP2nphi1M1npsi2M1npsi2,weight);
    if(sd == 0) {
      iBinCounter += 1;

      gSumBinContentTimesWeight += f3pCorrelatorVsPtSumDiffPro[sd]->GetBinContent(b)*f3pCorrelatorVsPtSumDiffPro[sd]->GetBinEntries(b);
      gSumWeight += f3pCorrelatorVsPtSumDiffPro[sd]->GetBinEntries(b);
      gSumBinContentTimesWeightSquared += TMath::Power(gSumBinContentTimesWeight,2);
    }
   } // end of for(Int_t b=1;b<=fnBinsPt;b++)

    // [(eta1+eta2)/2,|eta1-eta2|]
   // looping over all bins and calculating reduced correlations: 
   for(Int_t k=1;k<=fnBinsEta;k++)
   {
    // real and imaginary parts of p_{n}: 
    Double_t p1nRe = fReEtaEBE[sd]->GetBinContent(k)*fReEtaEBE[sd]->GetBinEntries(k);
    Double_t p1nIm = fImEtaEBE[sd]->GetBinContent(k)*fImEtaEBE[sd]->GetBinEntries(k);
    // overlap 1: to be improved (terminology)
    Double_t overlap1 = fOverlapEBE2[0][sd]->GetBinContent(k)*fOverlapEBE2[0][sd]->GetBinEntries(k);
    // overlap 2: to be improved (terminology)
    Double_t overlap2 = fOverlapEBE2[1][sd]->GetBinContent(k)*fOverlapEBE2[1][sd]->GetBinEntries(k);    
    // number of pairs of POIs in particular (eta1+eta2)/2 or |eta1-eta2| bin:
    Double_t mp = fReEtaEBE[sd]->GetBinEntries(k);
    // number of pairs of POI1/RP and POI2 in particular (eta1+eta2)/2 or |eta1-eta2| bin:
    Double_t mOverlap1 = fOverlapEBE2[0][sd]->GetBinEntries(k);
    // number of pairs of POI2/RP and POI1 in particular (eta1+eta2)/2 or |eta1-eta2| bin:
    Double_t mOverlap2 = fOverlapEBE2[1][sd]->GetBinEntries(k);
    // e-b-e weight for cos[n(psi1+psi2-2phi3)]:
    Double_t weight = mp*dMult-mOverlap1-mOverlap2;  
    
    Double_t cosP2nphi1M1npsi2M1npsi2 = 0; // cos[n(psi1+psi2-2phi3)]
    if(weight>0.)
    {
     cosP2nphi1M1npsi2M1npsi2 = (p1nRe*dReQ2n+p1nIm*dImQ2n-overlap1-overlap2)/(weight);
    }
    f3pCorrelatorVsEtaSumDiffPro[sd]->Fill(fEtaMin+(k-1)*fEtaBinWidth,cosP2nphi1M1npsi2M1npsi2,weight);
   } // end of for(Int_t k=1;k<=fnBinsEta;k++)
  } // end of for(Int_t sd=0;sd<2;sd++)      

  gIntegratedValue = -1000.;
  if((gSumWeight)&&(iBinCounter)) 
    gIntegratedValue = gSumBinContentTimesWeight/gSumWeight;
 } // end of if(!(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)) 

 // b) Calculate differential 3-p correlator by using particle weights: 
 if(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)
 {
  // ...
 } // end of if(fUsePhiWeights || fUsePtWeights || fUseEtaWeights)
 
} // end of void AliFlowAnalysisWithMixedHarmonics::CalculateDifferential3pCorrelator() 

//================================================================================================================

void AliFlowAnalysisWithMixedHarmonics::GetCorrelatorAndError(TProfile *g3pCorrelatorVsPt, Double_t &g3pCorrelatorValue, Double_t &g3pCorrelatorError) {
  //Retrieves the 3p correlator <<cos[n(psi1+psi2-2phi3)]>> 
  //and its error
  Double_t gSumXi = 0.;
  Double_t gSumYi = 0.;
  Double_t gSumXiYi = 0.;
  Double_t gSumXiYi2 = 0.;
  Double_t gSumXi2Yi2 = 0.;
  Double_t gSumDeltaXi2 = 0.;
  Double_t gSumYi2DeltaXi2 = 0.;

  for(Int_t iBin = 1; iBin <= g3pCorrelatorVsPt->GetNbinsX(); iBin++) {
    gSumXi += g3pCorrelatorVsPt->GetBinEntries(iBin);
    gSumYi += g3pCorrelatorVsPt->GetBinContent(iBin);
    gSumXiYi += g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinContent(iBin);
    gSumXiYi2 += g3pCorrelatorVsPt->GetBinEntries(iBin)*TMath::Power(g3pCorrelatorVsPt->GetBinContent(iBin),2);
    gSumXi2Yi2 += TMath::Power(g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinContent(iBin),2);
    gSumDeltaXi2 += TMath::Power(g3pCorrelatorVsPt->GetBinError(iBin),2);
    gSumYi2DeltaXi2 += TMath::Power(g3pCorrelatorVsPt->GetBinContent(iBin),2) + TMath::Power(g3pCorrelatorVsPt->GetBinError(iBin),2);
  }

  g3pCorrelatorValue = -1000.;
  g3pCorrelatorError = 1000.;
  
  if(gSumXi != 0.)
    g3pCorrelatorValue = gSumXiYi/gSumXi;
  if((gSumXi != 0.)&&(gSumXiYi != 0.))
    g3pCorrelatorError = TMath::Abs((gSumXiYi/gSumXi))*TMath::Sqrt(TMath::Power((TMath::Sqrt(gSumYi2DeltaXi2)/gSumXiYi),2) + TMath::Power((gSumDeltaXi2/gSumXi),2));
}
//================================================================================================================

