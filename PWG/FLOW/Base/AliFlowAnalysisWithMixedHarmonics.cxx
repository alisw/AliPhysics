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

using std::endl;
using std::cout;
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
fCommonConstants(NULL),
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
fPsi2V0C(0),
fPsi2V0A(0),
fPsiZNCA(0),
fPsiZNCC(0),
fPsiZNCCA(0),
fProfileList(NULL),
f2pCorrelatorCos2PsiDiff2PsiV0Pro(NULL),
f2pCorrelatorCos2PsiDiff2PsiZDCPro(NULL),
f3pCorrelatorPro(NULL),
f5pCorrelatorPro(NULL),
fNonIsotropicTermsPro(NULL),
f3pCorrelatorVsMPro(NULL),
f3pPOICorrelatorVsM(NULL),
fNonIsotropicTermsVsMPro(NULL),
fNonIsotropicTermsList(NULL),
fNonIsotropicPOITermsPro(NULL),
f2pCorrelatorCosPsiDiffPro(NULL),
f2pCorrelatorCos2PsiDiffPro(NULL),
f2pCorrelatorCosPsiDiffHist(NULL),
f2pCorrelatorCos2PsiDiffHist(NULL),
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
fDetectorBiasVsMHist(NULL),
fQAZDCCAvsV0C(NULL),
fQAZDCCAvsV0A(NULL),
fQAV0CEvPldistribution(NULL),
fQAV0AEvPldistribution(NULL),
fQAZDCCAEvPldistribution(NULL)
{
 // Constructor. 
 
 // Base list to hold all output objects:
 fHistList = new TList();
 fHistListName = new TString("cobjMH");
 fHistList->SetName(fHistListName->Data());
 fHistList->SetOwner(kTRUE);
 
 // List to hold histograms with phi, pt and eta weights:      
 fWeightsList = new TList();
 
 // Initialize event-by-event V0 and ZDC variables
 for(Int_t c=0;c<2;c++) {
  fVZFlowVect2ndHar[c] = AliFlowVector();
  fZDCFlowVect[c] = AliFlowVector();
 }
 
 // List to hold all all-event profiles:      
 fProfileList = new TList();
 
 // List to hold profiles with all non-isotropic terms for diff. correlators:      
 fNonIsotropicTermsList = new TList();

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
 this->AccessConstants("Init");
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
 
 // get V0 q-vector
 AliFlowVector vVZQarray[2];
 anEvent->GetV02Qsub(vVZQarray,2);
 fVZFlowVect2ndHar[0] = vVZQarray[0];
 fVZFlowVect2ndHar[1] = vVZQarray[1];
 // V0-C eta < 0
 Double_t VZCIm = fVZFlowVect2ndHar[0].X();
 Double_t VZCRe = fVZFlowVect2ndHar[0].Y();
 Double_t VZCM  = fVZFlowVect2ndHar[0].GetMult();
 // V0-A eta > 0
 Double_t VZAIm = fVZFlowVect2ndHar[1].X();
 Double_t VZARe = fVZFlowVect2ndHar[1].Y();
 Double_t VZAM  = fVZFlowVect2ndHar[1].GetMult();
 // Get event plane for V0
 fPsi2V0C = (1./2)*TMath::ATan2(VZCRe,VZCIm);
 if(fPsi2V0C < 0) fPsi2V0C += TMath::TwoPi()/2;
 fPsi2V0A = (1./2)*TMath::ATan2(VZARe,VZAIm);
 if(fPsi2V0A < 0) fPsi2V0A += TMath::TwoPi()/2;
 
 fNonIsotropicTermsV0EvPlPro->Fill(0.5, TMath::Cos(fPsi2V0C));
 fNonIsotropicTermsV0EvPlPro->Fill(1.5, TMath::Sin(fPsi2V0C));
 fNonIsotropicTermsV0EvPlPro->Fill(2.5, TMath::Cos(2*fPsi2V0C));
 fNonIsotropicTermsV0EvPlPro->Fill(3.5, TMath::Sin(2*fPsi2V0C));
 fNonIsotropicTermsV0EvPlPro->Fill(4.5, TMath::Cos(fPsi2V0A));
 fNonIsotropicTermsV0EvPlPro->Fill(5.5, TMath::Sin(fPsi2V0A));
 fNonIsotropicTermsV0EvPlPro->Fill(6.5, TMath::Cos(2*fPsi2V0A));
 fNonIsotropicTermsV0EvPlPro->Fill(7.5, TMath::Sin(2*fPsi2V0A));

 // get ZDC q-vector
 AliFlowVector vZDCQarray[2];
 anEvent->GetZDC2Qsub(vZDCQarray);
 AliFlowVector fZDCFlowVect[2];
 fZDCFlowVect[0] = vZDCQarray[0];
 fZDCFlowVect[1] = vZDCQarray[1];
 
 // ZDC-C (eta < -8.8)
 Double_t ZCIm = fZDCFlowVect[0].X();
 Double_t ZCRe = fZDCFlowVect[0].Y();
 Double_t ZCM  = fZDCFlowVect[0].GetMult();
 // ZDC-A (eta > 8.8)
 Double_t ZAIm = fZDCFlowVect[1].X();
 Double_t ZARe = fZDCFlowVect[1].Y();
 Double_t ZAM  = fZDCFlowVect[1].GetMult();

 fPsiZNCA = TMath::ATan2(ZARe,ZAIm); // Psi_{1,A} spectator plane -pi to pi
 if (fPsiZNCA < 0) { // Psi_{1,A} should be differ to Psi_{1,C} by pi. 
  fPsiZNCA = fPsiZNCA + TMath::Pi();
 } else if (fPsiZNCA >= 0) {
  fPsiZNCA = fPsiZNCA - TMath::Pi();
 }
 fPsiZNCC = TMath::ATan2(ZCRe,ZCIm); // Psi_{1,C} spectator plane 
 fPsiZNCCA = TMath::ATan2((ZCRe-ZARe),(ZCIm-ZAIm));  

 fNonIsotropicTermsZDCSpPlPPro->Fill(0.5, TMath::Cos(fPsiZNCC));
 fNonIsotropicTermsZDCSpPlPPro->Fill(1.5, TMath::Sin(fPsiZNCC));
 fNonIsotropicTermsZDCSpPlPPro->Fill(2.5, TMath::Cos(2*fPsiZNCC));
 fNonIsotropicTermsZDCSpPlPPro->Fill(3.5, TMath::Sin(2*fPsiZNCC));
 fNonIsotropicTermsZDCSpPlPPro->Fill(4.5, TMath::Cos(fPsiZNCA));
 fNonIsotropicTermsZDCSpPlPPro->Fill(5.5, TMath::Sin(fPsiZNCA));
 fNonIsotropicTermsZDCSpPlPPro->Fill(6.5, TMath::Cos(2*fPsiZNCA));
 fNonIsotropicTermsZDCSpPlPPro->Fill(7.5, TMath::Sin(2*fPsiZNCA));
 fNonIsotropicTermsZDCSpPlPPro->Fill(8.5, TMath::Cos(fPsiZNCCA));
 fNonIsotropicTermsZDCSpPlPPro->Fill(9.5, TMath::Sin(fPsiZNCCA));
 fNonIsotropicTermsZDCSpPlPPro->Fill(10.5, TMath::Cos(2*fPsiZNCCA));
 fNonIsotropicTermsZDCSpPlPPro->Fill(11.5, TMath::Sin(2*fPsiZNCCA));
 
 fQAZDCCAvsV0C->Fill(fPsiZNCCA, fPsi2V0C);
 fQAZDCCAvsV0A->Fill(fPsiZNCCA, fPsi2V0A);
 fQAV0CEvPldistribution->Fill(fPsi2V0C);
 fQAV0AEvPldistribution->Fill(fPsi2V0A);
 fQAZDCCAEvPldistribution->Fill(fPsiZNCCA);
 
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
    // Calculate <<cos(2a-2Psi_V0)>> for RP 
    f2pCorrelatorCos2PsiDiff2PsiV0Pro->Fill(0.5, TMath::Cos(2*(dPhi-fPsi2V0C)));
    f2pCorrelatorCos2PsiDiff2PsiV0Pro->Fill(1.5, TMath::Cos(2*(dPhi-fPsi2V0A)));
    // Calculate <<cos(2a-2Psi_ZDC)>> for RP
    f2pCorrelatorCos2PsiDiff2PsiZDCPro->Fill(0.5, TMath::Cos(2*(dPhi-fPsiZNCC)));
    f2pCorrelatorCos2PsiDiff2PsiZDCPro->Fill(1.5, TMath::Cos(2*(dPhi-fPsiZNCA)));
    
    f2pCorrelatorCos2PsiDiff2PsiZDCPro->Fill(2.5, TMath::Cos(2*(dPhi-fPsiZNCCA)));
    
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
     
     fNonIsotropicTermsV0Pro->Fill(0.5, TMath::Cos(dPsi1-2*fPsi2V0C));
     fNonIsotropicTermsV0Pro->Fill(1.5, TMath::Sin(dPsi1-2*fPsi2V0C));
     fNonIsotropicTermsV0Pro->Fill(2.5, TMath::Cos(dPsi1-2*fPsi2V0A));
     fNonIsotropicTermsV0Pro->Fill(3.5, TMath::Sin(dPsi1-2*fPsi2V0A));

     fNonIsotropicTermsZDCPro->Fill(0.5, TMath::Cos(dPsi1-2*fPsiZNCC));
     fNonIsotropicTermsZDCPro->Fill(1.5, TMath::Sin(dPsi1-2*fPsiZNCC));
     fNonIsotropicTermsZDCPro->Fill(2.5, TMath::Cos(dPsi1-2*fPsiZNCA));
     fNonIsotropicTermsZDCPro->Fill(3.5, TMath::Sin(dPsi1-2*fPsiZNCA));
     fNonIsotropicTermsZDCPro->Fill(4.5, TMath::Cos(dPsi1-2*fPsiZNCCA));
     fNonIsotropicTermsZDCPro->Fill(5.5, TMath::Sin(dPsi1-2*fPsiZNCCA));
     
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
       //2particle correlator fixed to 2nd order <cos(2*(psi1 - psi2))> profile
       f2pCorrelatorCosPsiDiffPro->Fill(0.5, TMath::Cos(dPsi1-dPsi2));
       f2pCorrelatorCos2PsiDiffPro->Fill(0.5, TMath::Cos(2*(dPsi1-dPsi2)));
       
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
       
       // non-isotropic terms for 2particle correlator <cos(psi1 - psi2)> and <cos(2*(psi1 - psi2))>
       fNonIsotropicPOITermsPro->Fill(0.5, TMath::Cos(dPsi1),1.);
       fNonIsotropicPOITermsPro->Fill(1.5, TMath::Sin(dPsi1),1.);
       fNonIsotropicPOITermsPro->Fill(2.5, TMath::Cos(dPsi2),1.);
       fNonIsotropicPOITermsPro->Fill(3.5, TMath::Sin(dPsi2),1.);
       fNonIsotropicPOITermsPro->Fill(4.5, TMath::Cos(2*(dPsi1)),1.);
       fNonIsotropicPOITermsPro->Fill(5.5, TMath::Sin(2*(dPsi1)),1.);
       fNonIsotropicPOITermsPro->Fill(6.5, TMath::Cos(2*(dPsi2)),1.);
       fNonIsotropicPOITermsPro->Fill(7.5, TMath::Sin(2*(dPsi2)),1.);
       
       // non-isotropic terms, 1st POI:
       fReNITEBE[0][0][0]->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi1)),1.);
       fReNITEBE[0][0][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi1)),1.);
       fReNITEBE[0][0][2]->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi1)),1.);
       fReNITEBE[0][0][3]->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi1)),1.);
       fImNITEBE[0][0][0]->Fill((dPt1+dPt2)/2.,TMath::Sin(n*(dPsi1)),1.);
       fImNITEBE[0][0][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Sin(n*(dPsi1)),1.);
       fImNITEBE[0][0][2]->Fill((dEta1+dEta2)/2.,TMath::Sin(n*(dPsi1)),1.);
       fImNITEBE[0][0][3]->Fill(TMath::Abs(dEta1-dEta2),TMath::Sin(n*(dPsi1)),1.);
       // non-isotropic terms, 2nd POI:
       fReNITEBE[1][0][0]->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi2)),1.);
       fReNITEBE[1][0][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi2)),1.);
       fReNITEBE[1][0][2]->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi2)),1.);
       fReNITEBE[1][0][3]->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi2)),1.);
       fImNITEBE[1][0][0]->Fill((dPt1+dPt2)/2.,TMath::Sin(n*(dPsi2)),1.);
       fImNITEBE[1][0][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Sin(n*(dPsi2)),1.);
       fImNITEBE[1][0][2]->Fill((dEta1+dEta2)/2.,TMath::Sin(n*(dPsi2)),1.);
       fImNITEBE[1][0][3]->Fill(TMath::Abs(dEta1-dEta2),TMath::Sin(n*(dPsi2)),1.);
	   
	   
       if(b1stPOIisAlsoRP)
       {
        fOverlapEBE[0][0]->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi1-dPsi2)),1.); // overlap in <<cos(psi1+psi2-2phi3)>> = <<cos[(psi1-psi2)-(2psi1-2phi3)]>>, if psi1 is both POI1 and RP and phi3 is RP, 2psi1-2phi3=0
        fOverlapEBE[0][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE2[0][0]->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE2[0][1]->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi1-dPsi2)),1.);
        // non-isotropic terms, 1st POI:
        fReNITEBE[0][1][0]->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi1)),1.);
        fReNITEBE[0][1][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi1)),1.);
        fReNITEBE[0][1][2]->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi1)),1.);
        fReNITEBE[0][1][3]->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi1)),1.);
        fImNITEBE[0][1][0]->Fill((dPt1+dPt2)/2.,TMath::Sin(n*(dPsi1)),1.);
        fImNITEBE[0][1][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Sin(n*(dPsi1)),1.);
        fImNITEBE[0][1][2]->Fill((dEta1+dEta2)/2.,TMath::Sin(n*(dPsi1)),1.);
        fImNITEBE[0][1][3]->Fill(TMath::Abs(dEta1-dEta2),TMath::Sin(n*(dPsi1)),1.);       
       }
       if(b2ndPOIisAlsoRP)
       {
        fOverlapEBE[1][0]->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE[1][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE2[1][0]->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi1-dPsi2)),1.);
        fOverlapEBE2[1][1]->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi1-dPsi2)),1.);
        // non-isotropic terms, 2nd POI:
        fReNITEBE[1][1][0]->Fill((dPt1+dPt2)/2.,TMath::Cos(n*(dPsi2)),1.);
        fReNITEBE[1][1][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Cos(n*(dPsi2)),1.);
        fReNITEBE[1][1][2]->Fill((dEta1+dEta2)/2.,TMath::Cos(n*(dPsi2)),1.);
        fReNITEBE[1][1][3]->Fill(TMath::Abs(dEta1-dEta2),TMath::Cos(n*(dPsi2)),1.);
        fImNITEBE[1][1][0]->Fill((dPt1+dPt2)/2.,TMath::Sin(n*(dPsi2)),1.);
        fImNITEBE[1][1][1]->Fill(TMath::Abs(dPt1-dPt2),TMath::Sin(n*(dPsi2)),1.);
        fImNITEBE[1][1][2]->Fill((dEta1+dEta2)/2.,TMath::Sin(n*(dPsi2)),1.);
        fImNITEBE[1][1][3]->Fill(TMath::Abs(dEta1-dEta2),TMath::Sin(n*(dPsi2)),1.);       
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
   Double_t gIntegrated3pCorrelator = 0., gIntegrated3pCorrelatorV0C = 0, gIntegrated3pCorrelatorV0A = 0, gIntegrated3pCorrelatorZDCC = 0, gIntegrated3pCorrelatorZDCA = 0, gIntegrated3pCorrelatorZDCCA = 0;
   this->CalculateDifferential3pCorrelator(gIntegrated3pCorrelator, gIntegrated3pCorrelatorV0C, gIntegrated3pCorrelatorV0A, gIntegrated3pCorrelatorZDCC, gIntegrated3pCorrelatorZDCA, gIntegrated3pCorrelatorZDCCA); // to be improved - add relevant if statements for the min # POIs as well]
   f3pCorrelatorPOIIntegratedPro->Fill(0.5, gIntegrated3pCorrelator);
   f3pCorrelatorPOIm2V0CIntegratedPro->Fill(0.5, gIntegrated3pCorrelatorV0C); // integrated differential 3-p correlator <<cos[psi1+psi2-2phi_V0C)]>>
   f3pCorrelatorPOIm2V0AIntegratedPro->Fill(0.5, gIntegrated3pCorrelatorV0A); // integrated differential 3-p correlator <<cos[psi1+psi2-2phi_V0A)]>>
   f3pCorrelatorPOIm2ZDCCIntegratedPro->Fill(0.5, gIntegrated3pCorrelatorZDCC); // integrated differential 3-p correlator <<cos[psi1+psi2-2phi_ZDCC)]>>
   f3pCorrelatorPOIm2ZDCAIntegratedPro->Fill(0.5, gIntegrated3pCorrelatorZDCA); // integrated differential 3-p correlator <<cos[psi1+psi2-2phi_ZDCA)]>>
   f3pCorrelatorPOIm2ZDCCAIntegratedPro->Fill(0.5, gIntegrated3pCorrelatorZDCCA); // integrated differential 3-p correlator <<cos[psi1+psi2-2phi_ZDCCA)]>>
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
 // b) Access common constants;
 // c) Access settings for analysis with mixed harmonics;
 // d) Correct for detector effects;
 // e) Print on the screen the final results.
 
 printf("AliFlowAnalysisWithMixedHarmonics::Finish()\n");
 
 this->CheckPointersUsedInFinish();
 this->AccessConstants("Finish");          
 this->AccessSettings();
 fCorrectForDetectorEffects = kTRUE;
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
   
 TString sCommonConstantsName = "fCommonConstants";
 fCommonConstants = dynamic_cast<TProfile*>(fHistList->FindObject(sCommonConstantsName.Data()));
 if(!fCommonConstants)
 {
  printf("\n WARNING (MH): fCommonConstants is NULL in GetPointersForBaseHistograms() !!!!\n\n");
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
 TList *nonIsotropicTermsList = NULL;
 nonIsotropicTermsList = dynamic_cast<TList*>(profileList->FindObject("Nonisotropic Terms"));
 if(!nonIsotropicTermsList && fEvaluateDifferential3pCorrelator) 
 {
  cout<<endl;
  cout<<" WARNING (MH): nonIsotropicTerms is NULL in GetPointersForAllEventProfiles() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }  
 
 TString s2pCorrelatorCos2PsiDiff2PsiV0ProName = "f2pCorrelatorCos2PsiDiff2PsiV0Pro";
 TProfile *p2pCorrelatorCos2PsiDiff2PsiV0Pro = dynamic_cast<TProfile*>(profileList->FindObject(s2pCorrelatorCos2PsiDiff2PsiV0ProName.Data()));
 if(p2pCorrelatorCos2PsiDiff2PsiV0Pro)
 {
  this->Set2pCorrelatorCos2PsiDiff2PsiV0Pro(p2pCorrelatorCos2PsiDiff2PsiV0Pro);  
 }
 TString s2pCorrelatorCos2PsiDiff2PsiZDCProName = "f2pCorrelatorCos2PsiDiff2PsiZDCPro";
 TProfile *p2pCorrelatorCos2PsiDiff2PsiZDCPro = dynamic_cast<TProfile*>(profileList->FindObject(s2pCorrelatorCos2PsiDiff2PsiZDCProName.Data()));
 if(p2pCorrelatorCos2PsiDiff2PsiZDCPro)
 {
  this->Set2pCorrelatorCos2PsiDiff2PsiZDCPro(p2pCorrelatorCos2PsiDiff2PsiZDCPro);  
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
 TString psdFlag2[2] = {"EtaSum","EtaDiff"};
 TString nonIsotropicTerm[10] = {"#LT#LTcos(#psi_{POI_1})#GT#GT","#LT#LTsin(#psi_{POI_1})#GT#GT",
                                 "#LT#LTcos(#psi_{POI_2})#GT#GT","#LT#LTsin(#psi_{POI_2})#GT#GT",
                                 "#LT#LTcos(#psi_{POI_1}-2#phi_{RP})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#phi_{RP})#GT#GT",
                                 "#LT#LTcos(#psi_{POI_2}-2#phi_{RP})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#phi_{RP})#GT#GT",
                                 "#LT#LTcos(#psi_{POI_1}+#psi_{POI_2})#GT#GT","#LT#LTsin(#psi_{POI_1}+#psi_{POI_2})#GT#GT"}; 

 TString nonIsotropicV0Term[8] = {"#LT#LTcos(#psi_{POI_1}-2#Psi_{V0C})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#Psi_{V0C})#GT#GT",
	                              "#LT#LTcos(#psi_{POI_1}-2#Psi_{V0A})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#Psi_{V0A})#GT#GT",
                                  "#LT#LTcos(#psi_{POI_2}-2#Psi_{V0C})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#Psi_{V0C})#GT#GT",
                                  "#LT#LTcos(#psi_{POI_2}-2#Psi_{V0A})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#Psi_{V0A})#GT#GT"};
                                 
 TString nonIsotropicZDCTerm[12] = {"#LT#LTcos(#psi_{POI_1}-2#Psi_{ZDCC})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#Psi_{ZDCC})#GT#GT",
	                              "#LT#LTcos(#psi_{POI_1}-2#Psi_{ZDCA})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#Psi_{ZDCA})#GT#GT",
	                              "#LT#LTcos(#psi_{POI_1}-2#Psi_{ZDCCA})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#Psi_{ZDCCA})#GT#GT",
                                  "#LT#LTcos(#psi_{POI_2}-2#Psi_{ZDCC})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#Psi_{ZDCC})#GT#GT",
                                  "#LT#LTcos(#psi_{POI_2}-2#Psi_{ZDCA})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#Psi_{ZDCA})#GT#GT",
                                  "#LT#LTcos(#psi_{POI_2}-2#Psi_{ZDCCA})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#Psi_{ZDCCA})#GT#GT"};
                                  
 for(Int_t sd=0;sd<2;sd++)
 {
  TProfile *p3pCorrelatorVsPtSumDiffPro = dynamic_cast<TProfile*>(profileList->FindObject(Form("f3pCorrelatorVs%sPro",psdFlag[sd].Data())));
  if(p3pCorrelatorVsPtSumDiffPro)
  {
   this->Set3pCorrelatorVsPtSumDiffPro(p3pCorrelatorVsPtSumDiffPro,sd);  
  }
  TProfile *p3pCorrelatorV0CVsPtSumDiffPro = dynamic_cast<TProfile*>(profileList->FindObject(Form("f3pCorrelatorV0CVs%sPro",psdFlag[sd].Data())));
  if(p3pCorrelatorV0CVsPtSumDiffPro)
  {
   this->Set3pCorrelatorV0CVsPtSumDiffPro(p3pCorrelatorV0CVsPtSumDiffPro,sd);  
  }
  TProfile *p3pCorrelatorV0AVsPtSumDiffPro = dynamic_cast<TProfile*>(profileList->FindObject(Form("f3pCorrelatorV0AVs%sPro",psdFlag[sd].Data())));
  if(p3pCorrelatorV0AVsPtSumDiffPro)
  {
   this->Set3pCorrelatorV0AVsPtSumDiffPro(p3pCorrelatorV0AVsPtSumDiffPro,sd);  
  }
  TProfile *p3pCorrelatorZDCCVsPtSumDiffPro = dynamic_cast<TProfile*>(profileList->FindObject(Form("f3pCorrelatorZDCCVs%sPro",psdFlag[sd].Data())));
  if(p3pCorrelatorZDCCVsPtSumDiffPro)
  {
   this->Set3pCorrelatorZDCCVsPtSumDiffPro(p3pCorrelatorZDCCVsPtSumDiffPro,sd);  
  }
  TProfile *p3pCorrelatorZDCAVsPtSumDiffPro = dynamic_cast<TProfile*>(profileList->FindObject(Form("f3pCorrelatorZDCAVs%sPro",psdFlag[sd].Data())));
  if(p3pCorrelatorZDCAVsPtSumDiffPro)
  {
   this->Set3pCorrelatorZDCAVsPtSumDiffPro(p3pCorrelatorZDCAVsPtSumDiffPro,sd);  
  }
  TProfile *p3pCorrelatorZDCCAVsPtSumDiffPro = dynamic_cast<TProfile*>(profileList->FindObject(Form("f3pCorrelatorZDCCAVs%sPro",psdFlag[sd].Data())));
  if(p3pCorrelatorZDCCAVsPtSumDiffPro)
  {
   this->Set3pCorrelatorZDCCAVsPtSumDiffPro(p3pCorrelatorZDCCAVsPtSumDiffPro,sd);  
  }
  TProfile *p3pCorrelatorVsEtaSumDiffPro = dynamic_cast<TProfile*>(profileList->FindObject(Form("f3pCorrelatorVs%sPro",psdFlag2[sd].Data())));
  if(p3pCorrelatorVsEtaSumDiffPro)
  {
   this->Set3pCorrelatorVsEtaSumDiffPro(p3pCorrelatorVsEtaSumDiffPro,sd);  
  }
  if(nonIsotropicTermsList) 
  {
   for(Int_t t=0;t<10;t++)
   {   
    // Pt:
    TProfile *pNonIsotropicTermsVsPtSumDiffPro = dynamic_cast<TProfile*>
                      (nonIsotropicTermsList->FindObject(Form("fNonIsotropicTermsVs%sPro %s",psdFlag[sd].Data(),nonIsotropicTerm[t].Data())));
    if(pNonIsotropicTermsVsPtSumDiffPro)
    {
     this->SetNonIsotropicTermsVsPtSumDiffPro(pNonIsotropicTermsVsPtSumDiffPro,sd,t);  
    } 
    // Eta:
    TProfile *pNonIsotropicTermsVsEtaSumDiffPro = dynamic_cast<TProfile*>
                      (nonIsotropicTermsList->FindObject(Form("fNonIsotropicTermsVs%sPro %s",psdFlag2[sd].Data(),nonIsotropicTerm[t].Data())));
    if(pNonIsotropicTermsVsEtaSumDiffPro)
    {
     this->SetNonIsotropicTermsVsEtaSumDiffPro(pNonIsotropicTermsVsEtaSumDiffPro,sd,t);  
    }                        
   } // end of for(Int_t t=0;t<10;t++)
  } // end of if(nonIsotropicTermsList)
 } // end of for(Int_t sd=0;sd<2;sd++)
 if(nonIsotropicTermsList) {
   TProfile *pNonIsotropicTermsV0Pro = dynamic_cast<TProfile *>(nonIsotropicTermsList->FindObject("fNonIsotropicTermsV0Pro"));
   if(pNonIsotropicTermsV0Pro)
     this->SetNonIsotropicTermsV0Pro(pNonIsotropicTermsV0Pro);
     
   TProfile *pNonIsotropicTermsZDCPro = dynamic_cast<TProfile *>(nonIsotropicTermsList->FindObject("fNonIsotropicTermsZDCPro"));
   if(pNonIsotropicTermsZDCPro)
     this->SetNonIsotropicTermsZDCPro(pNonIsotropicTermsZDCPro);
        
   TProfile *pNonIsotropicTermsV0EvPlPro = dynamic_cast<TProfile *>(nonIsotropicTermsList->FindObject("fNonIsotropicTermsV0EvPlPro"));
   if(pNonIsotropicTermsV0EvPlPro)
     this->SetNonIsotropicTermsV0EvPlPro(pNonIsotropicTermsV0EvPlPro);
  
   TProfile *pNonIsotropicTermsZDCSpPlPPro = dynamic_cast<TProfile *>(nonIsotropicTermsList->FindObject("fNonIsotropicTermsZDCSpPlPPro"));
   if(pNonIsotropicTermsZDCSpPlPPro)
     this->SetNonIsotropicTermsZDCSpPlPPro(pNonIsotropicTermsZDCSpPlPPro);
 }
 // non-isotropic terms for <<cos(2(phi1-phi2))>>
 TProfile *pNonIsotropicPOITermsPro = dynamic_cast<TProfile *>(profileList->FindObject("fNonIsotropicPOITermsPro"));
 if(pNonIsotropicPOITermsPro)
   this->SetNonIsotropicPOITermsPro(pNonIsotropicPOITermsPro);
  
 TProfile *p3pCorrelatorPOIIntegratedPro = dynamic_cast<TProfile*>(profileList->FindObject("f3pCorrelatorPOIIntegratedPro"));
 if(p3pCorrelatorPOIIntegratedPro)
 {
  this->Set3pCorrelatorPOIIntegratedPro(p3pCorrelatorPOIIntegratedPro);  
 }
 
 TProfile *p3pCorrelatorPOIm2V0CIntegratedPro = dynamic_cast<TProfile*>(profileList->FindObject("f3pCorrelatorPOIm2V0CIntegratedPro"));
 if(p3pCorrelatorPOIm2V0CIntegratedPro)
 {
  this->Set3pCorrelatorPOIm2V0CIntegratedPro(p3pCorrelatorPOIm2V0CIntegratedPro);  
 }
 TProfile *p3pCorrelatorPOIm2V0AIntegratedPro = dynamic_cast<TProfile*>(profileList->FindObject("f3pCorrelatorPOIm2V0AIntegratedPro"));
 if(p3pCorrelatorPOIm2V0AIntegratedPro)
 {
  this->Set3pCorrelatorPOIm2V0AIntegratedPro(p3pCorrelatorPOIm2V0AIntegratedPro);  
 }
 TProfile *p3pCorrelatorPOIm2ZDCCIntegratedPro = dynamic_cast<TProfile*>(profileList->FindObject("f3pCorrelatorPOIm2ZDCCIntegratedPro"));
 if(p3pCorrelatorPOIm2ZDCCIntegratedPro)
 {
  this->Set3pCorrelatorPOIm2ZDCCIntegratedPro(p3pCorrelatorPOIm2ZDCCIntegratedPro);  
 }
 TProfile *p3pCorrelatorPOIm2ZDCAIntegratedPro = dynamic_cast<TProfile*>(profileList->FindObject("f3pCorrelatorPOIm2ZDCAIntegratedPro"));
 if(p3pCorrelatorPOIm2ZDCAIntegratedPro)
 {
  this->Set3pCorrelatorPOIm2ZDCAIntegratedPro(p3pCorrelatorPOIm2ZDCAIntegratedPro);  
 }
 TProfile *p3pCorrelatorPOIm2ZDCCAIntegratedPro = dynamic_cast<TProfile*>(profileList->FindObject("f3pCorrelatorPOIm2ZDCCAIntegratedPro"));
 if(p3pCorrelatorPOIm2ZDCCAIntegratedPro)
 {
  this->Set3pCorrelatorPOIm2ZDCCAIntegratedPro(p3pCorrelatorPOIm2ZDCCAIntegratedPro);  
 }
 
   
 //2p correlator <<cos(psi1-psi2)>> and <<cos(2*(psi1 - psi2))>> profile
 TProfile *g2pCorrelatorCosPsiDiffPro = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorCosPsiDiffPro"));
 if(g2pCorrelatorCosPsiDiffPro)
   this->Set2pCorrelatorCosPsiDiffPro(g2pCorrelatorCosPsiDiffPro);
 TProfile *g2pCorrelatorCos2PsiDiffPro = dynamic_cast<TProfile *>(profileList->FindObject("f2pCorrelatorCos2PsiDiffPro"));
 if(g2pCorrelatorCos2PsiDiffPro)
   this->Set2pCorrelatorCos2PsiDiffPro(g2pCorrelatorCos2PsiDiffPro);
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
 
 TString psdFlag[2] = {"PtSum","PtDiff"}; 
 TString psdFlag2[2] = {"EtaSum","EtaDiff"};
 for(Int_t sd=0;sd<2;sd++)
 {
  TH1D *h3pCorrelatorVsPtSumDiffHist = dynamic_cast<TH1D*>(resultsList->FindObject(Form("f3pCorrelatorVs%sHist",psdFlag[sd].Data())));
  if(h3pCorrelatorVsPtSumDiffHist)
  {
   this->Set3pCorrelatorVsPtSumDiffHist(h3pCorrelatorVsPtSumDiffHist,sd);  
  }
  
  
  TH1D *h3pCorrelatorVsEtaSumDiffHist = dynamic_cast<TH1D*>(resultsList->FindObject(Form("f3pCorrelatorVs%sHist",psdFlag2[sd].Data())));
  if(h3pCorrelatorVsEtaSumDiffHist)
  {
   this->Set3pCorrelatorVsEtaSumDiffHist(h3pCorrelatorVsEtaSumDiffHist,sd);  
  } 
 } // end of for(Int_t sd=0;sd<2;sd++)
 
 TH1D *h3pCorrelatorPOIIntegratedHist = dynamic_cast<TH1D*>(resultsList->FindObject("f3pCorrelatorPOIIntegratedHist"));
 if(h3pCorrelatorPOIIntegratedHist)
 {
  this->Set3pCorrelatorPOIIntegratedHist(h3pCorrelatorPOIIntegratedHist);  
 }
 
 TH1D *h3pCorrelatorPOIm2V0CIntegratedHist = dynamic_cast<TH1D*>(resultsList->FindObject("f3pCorrelatorPOIm2V0CIntegratedHist"));
 if(h3pCorrelatorPOIm2V0CIntegratedHist)
 {
  this->Set3pCorrelatorPOIm2V0CIntegratedHist(h3pCorrelatorPOIm2V0CIntegratedHist);  
 }
 
 TH1D *h3pCorrelatorPOIm2V0AIntegratedHist = dynamic_cast<TH1D*>(resultsList->FindObject("f3pCorrelatorPOIm2V0AIntegratedHist"));
 if(h3pCorrelatorPOIm2V0AIntegratedHist)
 {
  this->Set3pCorrelatorPOIm2V0AIntegratedHist(h3pCorrelatorPOIm2V0AIntegratedHist);  
 }
 
 TH1D *h3pCorrelatorPOIm2ZDCCIntegratedHist = dynamic_cast<TH1D*>(resultsList->FindObject("f3pCorrelatorPOIm2ZDCCIntegratedHist"));
 if(h3pCorrelatorPOIm2ZDCCIntegratedHist)
 {
  this->Set3pCorrelatorPOIm2ZDCCIntegratedHist(h3pCorrelatorPOIm2ZDCCIntegratedHist);  
 }
 
 TH1D *h3pCorrelatorPOIm2ZDCAIntegratedHist = dynamic_cast<TH1D*>(resultsList->FindObject("f3pCorrelatorPOIm2ZDCAIntegratedHist"));
 if(h3pCorrelatorPOIm2ZDCAIntegratedHist)
 {
  this->Set3pCorrelatorPOIm2ZDCAIntegratedHist(h3pCorrelatorPOIm2ZDCAIntegratedHist);  
 }
 
 TH1D *h3pCorrelatorPOIm2ZDCCAIntegratedHist = dynamic_cast<TH1D*>(resultsList->FindObject("f3pCorrelatorPOIm2ZDCCAIntegratedHist"));
 if(h3pCorrelatorPOIm2ZDCCAIntegratedHist)
 {
  this->Set3pCorrelatorPOIm2ZDCCAIntegratedHist(h3pCorrelatorPOIm2ZDCCAIntegratedHist);  
 }
 
 TH1D *h2pCorrelatorCosPsiDiffHist = dynamic_cast<TH1D*>(resultsList->FindObject("f2pCorrelatorCosPsiDiffHist"));
 if(h2pCorrelatorCosPsiDiffHist)
 {
  this->Set2pCorrelatorCosPsiDiffHist(h2pCorrelatorCosPsiDiffHist);  
 }
 
 TH1D *h2pCorrelatorCos2PsiDiffHist = dynamic_cast<TH1D*>(resultsList->FindObject("f2pCorrelatorCos2PsiDiffHist"));
 if(h2pCorrelatorCos2PsiDiffHist)
 {
  this->Set2pCorrelatorCos2PsiDiffHist(h2pCorrelatorCos2PsiDiffHist);  
 }
 
 //QA for V0 and ZDCCA event plane
 TH2D *hQAZDCCAvsV0C = dynamic_cast<TH2D*>(resultsList->FindObject("fQAZDCCAvsV0C"));
 if(hQAZDCCAvsV0C)
 {
  this->SetQAZDCCAvsV0C(hQAZDCCAvsV0C);  
 }
 
 TH2D *hQAZDCCAvsV0A = dynamic_cast<TH2D*>(resultsList->FindObject("fQAZDCCAvsV0A"));
 if(hQAZDCCAvsV0A)
 {
  this->SetQAZDCCAvsV0A(hQAZDCCAvsV0A);  
 }
 
 TH1D *hQAV0CEvPldistribution = dynamic_cast<TH1D*>(resultsList->FindObject("fQAV0CEvPldistribution"));
 if(hQAV0CEvPldistribution)
 {
  this->SetQAV0CEvPldistribution(hQAV0CEvPldistribution);  
 }
 
 TH1D *hQAV0AEvPldistribution = dynamic_cast<TH1D*>(resultsList->FindObject("fQAV0AEvPldistribution"));
 if(hQAV0AEvPldistribution)
 {
  this->SetQAV0AEvPldistribution(hQAV0AEvPldistribution);  
 }
 
 TH1D *hQAZDCCAEvPldistribution = dynamic_cast<TH1D*>(resultsList->FindObject("fQAZDCCAEvPldistribution"));
 if(hQAZDCCAEvPldistribution)
 {
  this->SetQAZDCCAEvPldistribution(hQAZDCCAEvPldistribution);  
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
 fNonIsotropicTermsV0EvPlPro = NULL;
 fNonIsotropicTermsZDCSpPlPPro = NULL;
 
 for(Int_t sd=0;sd<2;sd++)
 {
  fRePEBE[sd] = NULL;
  fImPEBE[sd] = NULL;
  fReEtaEBE[sd] = NULL;
  fImEtaEBE[sd] = NULL;
  f3pCorrelatorVsPtSumDiffPro[sd] = NULL;
  f3pCorrelatorV0CVsPtSumDiffPro[sd] = NULL;
  f3pCorrelatorV0AVsPtSumDiffPro[sd] = NULL;
  f3pCorrelatorZDCCVsPtSumDiffPro[sd] = NULL;
  f3pCorrelatorZDCAVsPtSumDiffPro[sd] = NULL;
  f3pCorrelatorZDCCAVsPtSumDiffPro[sd] = NULL;
  f3pCorrelatorVsEtaSumDiffPro[sd] = NULL;
  f3pCorrelatorVsPtSumDiffHist[sd] = NULL;
  f3pCorrelatorVsEtaSumDiffHist[sd] = NULL;

  for(Int_t t=0;t<10;t++) // non-isotropic terms for diff. correlators
  {
   fNonIsotropicTermsVsPtSumDiffPro[sd][t] = NULL;
   fNonIsotropicTermsVsEtaSumDiffPro[sd][t] = NULL;
  } 
 } // end of for(Int_t sd=0;sd<2;sd++)
 
 fNonIsotropicTermsV0Pro = NULL;
 fNonIsotropicTermsZDCPro = NULL;
 
 f3pCorrelatorPOIIntegratedPro = NULL;
 f3pCorrelatorPOIm2V0CIntegratedPro = NULL;
 f3pCorrelatorPOIm2V0AIntegratedPro = NULL;
 f3pCorrelatorPOIm2ZDCCIntegratedPro = NULL;
 f3pCorrelatorPOIm2ZDCAIntegratedPro = NULL;
 f3pCorrelatorPOIm2ZDCCAIntegratedPro = NULL;
 f3pCorrelatorPOIIntegratedHist = NULL;
 f3pCorrelatorPOIm2V0CIntegratedHist = NULL;
 f3pCorrelatorPOIm2V0AIntegratedHist = NULL;
 f3pCorrelatorPOIm2ZDCCIntegratedHist = NULL;
 f3pCorrelatorPOIm2ZDCAIntegratedHist = NULL;
 f3pCorrelatorPOIm2ZDCCAIntegratedHist = NULL;
 for(Int_t fs=0;fs<2;fs++) // 1st/2nd POI which is also RP
 {
  for(Int_t sd=0;sd<2;sd++)
  {
   fOverlapEBE[fs][sd] = NULL;
   fOverlapEBE2[fs][sd] = NULL;
  }  
 } // end of for(Int_t fs=0;fs<2;fs++) // 1st/2nd POI which is also RP
 for(Int_t p12=0;p12<2;p12++) // 1st/2nd POI 
 {
  for(Int_t ao=0;ao<2;ao++) // all/overlap
  {
   for(Int_t pe=0;pe<4;pe++) // [(p1+p2)/2,|p1-p2|,(eta1+eta2)/2,|eta1-eta2|]
   {
    fReNITEBE[p12][ao][pe] = NULL;
    fImNITEBE[p12][ao][pe] = NULL;
   } // end of for(Int_t pe=0;pe<4;pe++) // [(p1+p2)/2,|p1-p2|,(eta1+eta2)/2,|eta1-eta2|]
  } // end of for(Int_t ao=0;ao<2;ao++) // all/overlap
 } // end of for(Int_t p12=0;p12<2;p12++) // 1st/2nd POI 
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
 // Profiles with non-isotropic terms for diff. correlators:
 fNonIsotropicTermsList->SetName("Nonisotropic Terms");
 fNonIsotropicTermsList->SetOwner(kTRUE);   
 if(fEvaluateDifferential3pCorrelator){fProfileList->Add(fNonIsotropicTermsList);} 

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
 fAnalysisSettings->Fill(7.5,(Int_t)fEvaluateDifferential3pCorrelator); 
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
 Int_t nBinsPtEta[4] = {fnBinsPt,fnBinsPt,fnBinsEta,fnBinsEta};
 Double_t dPtEtaMin[4] = {0.,0.,fEtaMin,fEtaMin}; 
 Double_t dPtEtaMax[4] = {fPtMax,fPtMax,fEtaMax,fEtaMax}; 
 for(Int_t p12=0;p12<2;p12++) // 1st/2nd POI 
 {
  for(Int_t ao=0;ao<2;ao++) // all/overlap
  {
   for(Int_t pe=0;pe<4;pe++) // [(p1+p2)/2,|p1-p2|,(eta1+eta2)/2,|eta1-eta2|]
   {
    fReNITEBE[p12][ao][pe] = new TProfile(Form("fReNITEBE%d%d%d",p12,ao,pe),"",nBinsPtEta[pe],dPtEtaMin[pe],dPtEtaMax[pe]);
    fImNITEBE[p12][ao][pe] = new TProfile(Form("fImNITEBE%d%d%d",p12,ao,pe),"",nBinsPtEta[pe],dPtEtaMin[pe],dPtEtaMax[pe]);
   } // end of for(Int_t pe=0;pe<4;pe++) // [(p1+p2)/2,|p1-p2|,(eta1+eta2)/2,|eta1-eta2|]
  } // end of for(Int_t ao=0;ao<2;ao++) // all/overlap
 } // end of for(Int_t p12=0;p12<2;p12++) // 1st/2nd POI 

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
 
 // a) 2-p (RP) correlator <<cos(2*(phi1-Psi_{V0})>> and <<cos(2*(phi1-Psi_{ZDC})>>
 // b) 3-p correlator <<cos[n*(phi1+phi2-2phi3)]>> for all events (not corrected for detector effects);
 // c) Non-isotropic terms in the decomposition of <<cos[n(phi1+phi2-2phi3)]>>;
 // d) 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> corrected for detector effects;
 // d.2) 3-p correlator where phi1 and phi2 are POI1 and POI2 and phi3 is RP <<cos[n(phi1+phi2-2phi3)]>> corrected for detector effects;
 // e) Histogram which quantifies bias coming from detector inefficiencies to 3-p correlator <<cos[n(phi1+phi2-2phi3)]>>;
 // f) 5-p correlator <<cos[n*(2phi1+2phi2+2phi3-3phi4-3phi5)]>> for all events (not corrected for detector effects - not supported yet);
 // g) 2-p (POI) correlator <<cos(phi1-phi2)>> and <<cos[2(phi1-phi2)]>> and non-isotropic terms for POI (not corrected for detector effects);
 // h) 2-p (POI) correlator <<cos(phi1-phi2)>> and <<cos[2(phi1-phi2)]>> (corrected for detector effects).
 // i) QA for V0 and ZDC event plane
 
 // a) 2-p (RP) correlator <<cos(2*(phi1-phi_{V0})>> and <<cos(2*(phi1-Psi_{ZDC})>>
 TString s2pCorrelatorCos2PsiDiff2PsiV0ProName = "f2pCorrelatorCos2PsiDiff2PsiV0Pro";
 f2pCorrelatorCos2PsiDiff2PsiV0Pro = new TProfile(s2pCorrelatorCos2PsiDiff2PsiV0ProName.Data(),"",2,0,2);
 f2pCorrelatorCos2PsiDiff2PsiV0Pro->SetStats(kFALSE);
 f2pCorrelatorCos2PsiDiff2PsiV0Pro->GetXaxis()->SetLabelOffset(0.01);
 f2pCorrelatorCos2PsiDiff2PsiV0Pro->GetXaxis()->SetLabelSize(0.05);
 f2pCorrelatorCos2PsiDiff2PsiV0Pro->GetXaxis()->SetBinLabel(1,"#LT#LTcos(2(#phi_{1}-#Psi_{V0C}))#GT#GT");
 f2pCorrelatorCos2PsiDiff2PsiV0Pro->GetXaxis()->SetBinLabel(2,"#LT#LTcos(2(#phi_{1}-#Psi_{V0A}))#GT#GT");
 fProfileList->Add(f2pCorrelatorCos2PsiDiff2PsiV0Pro);
 
 TString s2pCorrelatorCos2PsiDiff2PsiZDCProName = "f2pCorrelatorCos2PsiDiff2PsiZDCPro";
 f2pCorrelatorCos2PsiDiff2PsiZDCPro = new TProfile(s2pCorrelatorCos2PsiDiff2PsiZDCProName.Data(),"",3,0,3);
 f2pCorrelatorCos2PsiDiff2PsiZDCPro->SetStats(kFALSE);
 f2pCorrelatorCos2PsiDiff2PsiZDCPro->GetXaxis()->SetLabelOffset(0.01);
 f2pCorrelatorCos2PsiDiff2PsiZDCPro->GetXaxis()->SetLabelSize(0.05);
 f2pCorrelatorCos2PsiDiff2PsiZDCPro->GetXaxis()->SetBinLabel(1,"#LT#LTcos(2(#phi_{1}-#Psi_{ZDCC}))#GT#GT");
 f2pCorrelatorCos2PsiDiff2PsiZDCPro->GetXaxis()->SetBinLabel(2,"#LT#LTcos(2(#phi_{1}-#Psi_{ZDCA}))#GT#GT");
 f2pCorrelatorCos2PsiDiff2PsiZDCPro->GetXaxis()->SetBinLabel(3,"#LT#LTcos(2(#phi_{1}-#Psi_{ZDCCA}))#GT#GT");
 fProfileList->Add(f2pCorrelatorCos2PsiDiff2PsiZDCPro);
    
 // b) 3-p correlator <<cos[n*(phi1+phi2-2phi3)]>> for all events (not corrected for detector effects);
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
 
 // c) Non-isotropic terms in the decomposition of <<cos[n(phi1+phi2-2phi3)]>>:
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
 
 // d) 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> corrected for detector effects:
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
 
 // d.2) 3-p correlator where phi1 and phi2 are POI1 and POI2 and phi3 is RP <<cos[n(phi1+phi2-2phi3)]>> corrected for detector effects:
 f3pCorrelatorPOIIntegratedHist = new TH1D("f3pCorrelatorPOIIntegratedHist","",1,0.,1);
 f3pCorrelatorPOIIntegratedHist->SetStats(kFALSE);
 if(fHarmonic == 1) {
   f3pCorrelatorPOIIntegratedHist->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{3})#GT#GT POI");
 } else {
   f3pCorrelatorPOIIntegratedHist->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2}-2#phi_{3})]#GT#GT POI",fHarmonic));
 }
 fResultsList->Add(f3pCorrelatorPOIIntegratedHist);
 
 f3pCorrelatorPOIm2V0CIntegratedHist = new TH1D("f3pCorrelatorPOIm2V0CIntegratedHist","",1,0.,1);
 f3pCorrelatorPOIm2V0CIntegratedHist->SetStats(kFALSE);
 f3pCorrelatorPOIm2V0CIntegratedHist->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#Psi_{V0C})#GT#GT POI");
 fResultsList->Add(f3pCorrelatorPOIm2V0CIntegratedHist);
 
 f3pCorrelatorPOIm2V0AIntegratedHist = new TH1D("f3pCorrelatorPOIm2V0AIntegratedHist","",1,0.,1);
 f3pCorrelatorPOIm2V0AIntegratedHist->SetStats(kFALSE);
 f3pCorrelatorPOIm2V0AIntegratedHist->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#Psi_{V0A})#GT#GT POI");
 fResultsList->Add(f3pCorrelatorPOIm2V0AIntegratedHist);
 
 f3pCorrelatorPOIm2ZDCCIntegratedHist = new TH1D("f3pCorrelatorPOIm2ZDCCIntegratedHist","",1,0.,1);
 f3pCorrelatorPOIm2ZDCCIntegratedHist->SetStats(kFALSE);
 f3pCorrelatorPOIm2ZDCCIntegratedHist->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#Psi_{ZDCC})#GT#GT POI");
 fResultsList->Add(f3pCorrelatorPOIm2ZDCCIntegratedHist);
 
 f3pCorrelatorPOIm2ZDCAIntegratedHist = new TH1D("f3pCorrelatorPOIm2ZDCAIntegratedHist","",1,0.,1);
 f3pCorrelatorPOIm2ZDCAIntegratedHist->SetStats(kFALSE);
 f3pCorrelatorPOIm2ZDCAIntegratedHist->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#Psi_{ZDCA})#GT#GT POI");
 fResultsList->Add(f3pCorrelatorPOIm2ZDCAIntegratedHist);
 
 f3pCorrelatorPOIm2ZDCCAIntegratedHist = new TH1D("f3pCorrelatorPOIm2ZDCCAIntegratedHist","",2,0.,2);
 f3pCorrelatorPOIm2ZDCCAIntegratedHist->SetStats(kFALSE);
 f3pCorrelatorPOIm2ZDCCAIntegratedHist->GetXaxis()->SetBinLabel(1,"Corrected");
 f3pCorrelatorPOIm2ZDCCAIntegratedHist->GetXaxis()->SetBinLabel(2,"over Corrected");
 f3pCorrelatorPOIm2ZDCCAIntegratedHist->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#Psi_{ZDCCA})#GT#GT POI");
 fResultsList->Add(f3pCorrelatorPOIm2ZDCCAIntegratedHist);

 // e) Histogram which quantifies bias coming from detector inefficiencies to 3-p correlator <<cos[n(phi1+phi2-2phi3)]>>:
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
 
 // f) 5-p correlator <<cos[n*(2phi1+2phi2+2phi3-3phi4-3phi5)]>> for all events (not corrected for detector effects - not supported yet):
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
 
 // g) 2-p (POI) correlator <<cos(phi1-phi2)>> and <<cos[2(phi1-phi2)]>> and non-isotropic terms for POI (not corrected for detector effects):
 // non-isotropic terms for <<cos(phi1-phi2)>> and <<cos(2(phi1-phi2))>>
 fNonIsotropicPOITermsPro = new TProfile("fNonIsotropicPOITermsPro","non-isotropic terms for #LT #LT cos(2(#psi_{1} - #psi_{2})) #GT #GT", 8, 0, 8);
 fNonIsotropicPOITermsPro->SetStats(kFALSE);
 fNonIsotropicPOITermsPro->GetXaxis()->SetBinLabel(1,"cos(#phi_{POI_1})");
 fNonIsotropicPOITermsPro->GetXaxis()->SetBinLabel(2,"sin(#phi_{POI_1})");
 fNonIsotropicPOITermsPro->GetXaxis()->SetBinLabel(3,"cos(#phi_{POI_2})");
 fNonIsotropicPOITermsPro->GetXaxis()->SetBinLabel(4,"sin(#phi_{POI_2})");	
 fNonIsotropicPOITermsPro->GetXaxis()->SetBinLabel(5,"cos(2#phi_{POI_1})");
 fNonIsotropicPOITermsPro->GetXaxis()->SetBinLabel(6,"sin(2#phi_{POI_1})");
 fNonIsotropicPOITermsPro->GetXaxis()->SetBinLabel(7,"cos(2#phi_{POI_2})");
 fNonIsotropicPOITermsPro->GetXaxis()->SetBinLabel(8,"sin(2#phi_{POI_2})");	
 fProfileList->Add(fNonIsotropicPOITermsPro);
 
 // 2p correlator <<cos(psi1-psi2)>> and <<cos(2*(psi1 - psi2))>> profile
 f2pCorrelatorCosPsiDiffPro = new TProfile("f2pCorrelatorCosPsiDiffPro","; ;#LT #LT cos(#psi_{1} - #psi_{2}) #GT #GT",1,0.,1);
 f2pCorrelatorCosPsiDiffPro->SetStats(kFALSE);
 f2pCorrelatorCosPsiDiffPro->SetTitle("#LT#LTcos(#psi_{1}-#psi_{2})#GT#GT #font[72]");
 f2pCorrelatorCos2PsiDiffPro = new TProfile("f2pCorrelatorCos2PsiDiffPro",";;#LT #LT cos(2(#psi_{1} - #psi_{2})) #GT #GT",1,0.,1);
 f2pCorrelatorCos2PsiDiffPro->SetStats(kFALSE);
 f2pCorrelatorCos2PsiDiffPro->SetTitle("#LT#LTcos(2(#psi_{1}-#psi_{2}))#GT#GT #font[72]");
 fProfileList->Add(f2pCorrelatorCosPsiDiffPro);
 fProfileList->Add(f2pCorrelatorCos2PsiDiffPro);
 
 // h) 2-p (POI) correlator <<cos(phi1-phi2)>> and <<cos[2(phi1-phi2)]>> (corrected for detector effects):
 f2pCorrelatorCosPsiDiffHist = new TH1D("f2pCorrelatorCosPsiDiffHist","; ;#LT #LT cos(#psi_{1} - #psi_{2}) #GT #GT",1,0.,1);
 f2pCorrelatorCosPsiDiffHist->SetStats(kFALSE);
 f2pCorrelatorCosPsiDiffHist->SetTitle("#LT#LTcos(#psi_{1}-#psi_{2})#GT#GT #font[72]");
 f2pCorrelatorCos2PsiDiffHist = new TH1D("f2pCorrelatorCos2PsiDiffHist",";;#LT #LT cos(2(#psi_{1} - #psi_{2})) #GT #GT",1,0.,1);
 f2pCorrelatorCos2PsiDiffHist->SetStats(kFALSE);
 f2pCorrelatorCos2PsiDiffHist->SetTitle("#LT#LTcos(2(#psi_{1}-#psi_{2}))#GT#GT #font[72]");
 fResultsList->Add(f2pCorrelatorCosPsiDiffHist);
 fResultsList->Add(f2pCorrelatorCos2PsiDiffHist);
 
 // i) QA for V0 and ZDC event plane
 fQAZDCCAvsV0C = new TH2D("fQAZDCCAvsV0C","fQAZDCCAvsV0C",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
 fQAZDCCAvsV0C->Sumw2();
 fQAZDCCAvsV0C->SetXTitle("#Psi_{ZDCCA}");
 fQAZDCCAvsV0C->SetYTitle("#Psi_{V0C}");
 fResultsList->Add(fQAZDCCAvsV0C);
  
 fQAZDCCAvsV0A = new TH2D("fQAZDCCAvsV0A","fQAZDCCAvsV0A",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
 fQAZDCCAvsV0A->Sumw2();
 fQAZDCCAvsV0A->SetXTitle("#Psi_{ZDCCA}");
 fQAZDCCAvsV0A->SetYTitle("#Psi_{V0A}");
 fResultsList->Add(fQAZDCCAvsV0A);
 
 fQAV0CEvPldistribution = new TH1D("fQAV0CEvPldistribution","fQAV0CEvPldistribution",100,0,TMath::Pi());
 fQAV0CEvPldistribution->Sumw2();
 fQAV0CEvPldistribution->SetXTitle("#Psi_{V0C}");
 fResultsList->Add(fQAV0CEvPldistribution);
 
 fQAV0AEvPldistribution = new TH1D("fQAV0AEvPldistribution","fQAV0AEvPldistribution",100,0,TMath::Pi());
 fQAV0AEvPldistribution->Sumw2();
 fQAV0AEvPldistribution->SetXTitle("#Psi_{V0A}");
 fResultsList->Add(fQAV0AEvPldistribution);
 
 fQAZDCCAEvPldistribution = new TH1D("fQAZDCCAEvPldistribution","fQAZDCCAEvPldistribution",100,-TMath::Pi(),TMath::Pi());
 fQAZDCCAEvPldistribution->Sumw2();
 fQAZDCCAEvPldistribution->SetXTitle("#Psi_{ZDCCA}");
 fResultsList->Add(fQAZDCCAEvPldistribution);


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
  f3pCorrelatorV0CVsPtSumDiffPro[sd] = new TProfile(Form("f3pCorrelatorV0CVs%sPro",psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
  f3pCorrelatorV0CVsPtSumDiffPro[sd]->SetStats(kFALSE);
  f3pCorrelatorV0AVsPtSumDiffPro[sd] = new TProfile(Form("f3pCorrelatorV0AVs%sPro",psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
  f3pCorrelatorV0AVsPtSumDiffPro[sd]->SetStats(kFALSE);
  f3pCorrelatorZDCCVsPtSumDiffPro[sd] = new TProfile(Form("f3pCorrelatorZDCCVs%sPro",psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
  f3pCorrelatorZDCCVsPtSumDiffPro[sd]->SetStats(kFALSE);
  f3pCorrelatorZDCAVsPtSumDiffPro[sd] = new TProfile(Form("f3pCorrelatorZDCAVs%sPro",psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
  f3pCorrelatorZDCAVsPtSumDiffPro[sd]->SetStats(kFALSE);
  f3pCorrelatorZDCCAVsPtSumDiffPro[sd] = new TProfile(Form("f3pCorrelatorZDCCAVs%sPro",psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
  f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->SetStats(kFALSE);
  f3pCorrelatorVsEtaSumDiffPro[sd] = new TProfile(Form("f3pCorrelatorVs%sPro",psdFlag2[sd].Data()),"",fnBinsEta,fEtaMin,fEtaMax);
  f3pCorrelatorVsEtaSumDiffPro[sd]->SetStats(kFALSE);
  //f3pCorrelatorVsPtSumDiffPro[sd]->SetLabelSize(0.05);
  //f3pCorrelatorVsPtSumDiffPro[sd]->SetMarkerStyle(25);
  f3pCorrelatorV0CVsPtSumDiffPro[sd]->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2}-2#Psi_{V0C})#GT#GT #font[72]{vs} %s",psdTitleFlag[sd].Data())); 
  f3pCorrelatorV0AVsPtSumDiffPro[sd]->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2}-2#Psi_{V0A})#GT#GT #font[72]{vs} %s",psdTitleFlag[sd].Data())); 
  f3pCorrelatorZDCCVsPtSumDiffPro[sd]->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2}-2#Psi_{ZDCC})#GT#GT #font[72]{vs} %s",psdTitleFlag[sd].Data())); 
  f3pCorrelatorZDCAVsPtSumDiffPro[sd]->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2}-2#Psi_{ZDCA})#GT#GT #font[72]{vs} %s",psdTitleFlag[sd].Data())); 
  f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2}-2#Psi_{ZDCCA})#GT#GT #font[72]{vs} %s",psdTitleFlag[sd].Data())); 
   
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
  f3pCorrelatorV0CVsPtSumDiffPro[sd]->GetXaxis()->SetTitle(psdTitleFlag[sd].Data());
  fProfileList->Add(f3pCorrelatorV0CVsPtSumDiffPro[sd]);
  f3pCorrelatorV0AVsPtSumDiffPro[sd]->GetXaxis()->SetTitle(psdTitleFlag[sd].Data());
  fProfileList->Add(f3pCorrelatorV0AVsPtSumDiffPro[sd]);
  f3pCorrelatorZDCCVsPtSumDiffPro[sd]->GetXaxis()->SetTitle(psdTitleFlag[sd].Data());
  fProfileList->Add(f3pCorrelatorZDCCVsPtSumDiffPro[sd]);
  f3pCorrelatorZDCAVsPtSumDiffPro[sd]->GetXaxis()->SetTitle(psdTitleFlag[sd].Data());
  fProfileList->Add(f3pCorrelatorZDCAVsPtSumDiffPro[sd]);
  f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->GetXaxis()->SetTitle(psdTitleFlag[sd].Data());
  fProfileList->Add(f3pCorrelatorZDCCAVsPtSumDiffPro[sd]);
  f3pCorrelatorVsEtaSumDiffPro[sd]->GetXaxis()->SetTitle(psdTitleFlag2[sd].Data());
  fProfileList->Add(f3pCorrelatorVsEtaSumDiffPro[sd]);
 }
 
 f3pCorrelatorPOIIntegratedPro = new TProfile("f3pCorrelatorPOIIntegratedPro","",1,0.,1);
 f3pCorrelatorPOIIntegratedPro->SetStats(kFALSE);
 if(fHarmonic == 1) {
   f3pCorrelatorPOIIntegratedPro->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{3})#GT#GT POI");
 } else {
   f3pCorrelatorPOIIntegratedPro->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2}-2#phi_{3})]#GT#GT POI",fHarmonic));
 }
 fProfileList->Add(f3pCorrelatorPOIIntegratedPro);
 
 f3pCorrelatorPOIm2V0CIntegratedPro = new TProfile("f3pCorrelatorPOIm2V0CIntegratedPro","",1,0.,1);
 f3pCorrelatorPOIm2V0CIntegratedPro->SetStats(kFALSE);
 f3pCorrelatorPOIm2V0CIntegratedPro->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{V0C})#GT#GT POI");
 fProfileList->Add(f3pCorrelatorPOIm2V0CIntegratedPro);
 f3pCorrelatorPOIm2V0AIntegratedPro = new TProfile("f3pCorrelatorPOIm2V0AIntegratedPro","",1,0.,1);
 f3pCorrelatorPOIm2V0AIntegratedPro->SetStats(kFALSE);
 f3pCorrelatorPOIm2V0AIntegratedPro->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{V0A})#GT#GT POI");
 fProfileList->Add(f3pCorrelatorPOIm2V0AIntegratedPro);
 f3pCorrelatorPOIm2ZDCCIntegratedPro = new TProfile("f3pCorrelatorPOIm2ZDCCIntegratedPro","",1,0.,1);
 f3pCorrelatorPOIm2ZDCCIntegratedPro->SetStats(kFALSE);
 f3pCorrelatorPOIm2ZDCCIntegratedPro->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{ZDCC})#GT#GT POI");
 fProfileList->Add(f3pCorrelatorPOIm2ZDCCIntegratedPro);
 f3pCorrelatorPOIm2ZDCAIntegratedPro = new TProfile("f3pCorrelatorPOIm2ZDCAIntegratedPro","",1,0.,1);
 f3pCorrelatorPOIm2ZDCAIntegratedPro->SetStats(kFALSE);
 f3pCorrelatorPOIm2ZDCAIntegratedPro->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{ZDCA})#GT#GT POI");
 fProfileList->Add(f3pCorrelatorPOIm2ZDCAIntegratedPro);
 f3pCorrelatorPOIm2ZDCCAIntegratedPro = new TProfile("f3pCorrelatorPOIm2ZDCCAIntegratedPro","",1,0.,1);
 f3pCorrelatorPOIm2ZDCCAIntegratedPro->SetStats(kFALSE);
 f3pCorrelatorPOIm2ZDCCAIntegratedPro->SetTitle("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{ZDCCA})#GT#GT POI");
 fProfileList->Add(f3pCorrelatorPOIm2ZDCCAIntegratedPro);
 
 // Corrected for detector effects:
 for(Int_t sd=0;sd<2;sd++)
 {
  f3pCorrelatorVsPtSumDiffHist[sd] = new TH1D(Form("f3pCorrelatorVs%sHist",psdFlag[sd].Data()),"",fnBinsPt,0.,fPtMax);
  f3pCorrelatorVsPtSumDiffHist[sd]->SetStats(kFALSE);
  f3pCorrelatorVsEtaSumDiffHist[sd] = new TH1D(Form("f3pCorrelatorVs%sHist",psdFlag2[sd].Data()),"",fnBinsEta,fEtaMin,fEtaMax);
  f3pCorrelatorVsEtaSumDiffHist[sd]->SetStats(kFALSE);
  //f3pCorrelatorVsPtSumDiffHist[sd]->SetLabelSize(0.05);
  //f3pCorrelatorVsPtSumDiffHist[sd]->SetMarkerStyle(25);
  if(fHarmonic == 1)
  {
   f3pCorrelatorVsPtSumDiffHist[sd]->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{3})#GT#GT #font[72]{vs} %s",psdTitleFlag[sd].Data())); 
   f3pCorrelatorVsEtaSumDiffHist[sd]->SetTitle(Form("#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{3})#GT#GT #font[72]{vs} %s",psdTitleFlag2[sd].Data())); 
  } else
    {
     f3pCorrelatorVsPtSumDiffHist[sd]->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2}-2#phi_{3})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag[sd].Data())); 
     f3pCorrelatorVsEtaSumDiffHist[sd]->SetTitle(Form("#LT#LTcos[%d(#psi_{1}+#psi_{2}-2#phi_{3})]#GT#GT #font[72]{vs} %s",fHarmonic,psdTitleFlag2[sd].Data())); 
    }   
  
  f3pCorrelatorVsPtSumDiffHist[sd]->GetXaxis()->SetTitle(psdTitleFlag[sd].Data());
  fResultsList->Add(f3pCorrelatorVsPtSumDiffHist[sd]);
  f3pCorrelatorVsEtaSumDiffHist[sd]->GetXaxis()->SetTitle(psdTitleFlag2[sd].Data());
  fResultsList->Add(f3pCorrelatorVsEtaSumDiffHist[sd]);
 }
 
 //TString psdFlag[2] = {"PtSum","PtDiff"};
 //TString psdTitleFlag[2] = {"(p_{T,1}+ p_{T,2})/2","#left|p_{T,1}- p_{T,2}#right|"};
 //TString psdFlag2[2] = {"EtaSum","EtaDiff"};
 //TString psdTitleFlag2[2] = {"(#eta_{1}+ #eta_{2})/2","#left|#eta_{1}- #eta_{2}#right|"};
 TString nonIsotropicTerm[10] = {"#LT#LTcos(#psi_{POI_1})#GT#GT","#LT#LTsin(#psi_{POI_1})#GT#GT",
                                 "#LT#LTcos(#psi_{POI_2})#GT#GT","#LT#LTsin(#psi_{POI_2})#GT#GT",
                                 "#LT#LTcos(#psi_{POI_1}-2#phi_{RP})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#phi_{RP})#GT#GT",
                                 "#LT#LTcos(#psi_{POI_2}-2#phi_{RP})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#phi_{RP})#GT#GT",
                                 "#LT#LTcos(#psi_{POI_1}+#psi_{POI_2})#GT#GT","#LT#LTsin(#psi_{POI_1}+#psi_{POI_2})#GT#GT"};
 TString nonIsotropicV0Term[8] = {"#LT#LTcos(#psi_{POI_1}-2#Psi_{V0C})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#Psi_{V0C})#GT#GT",
	                              "#LT#LTcos(#psi_{POI_1}-2#Psi_{V0A})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#Psi_{V0A})#GT#GT",
                                  "#LT#LTcos(#psi_{POI_2}-2#Psi_{V0C})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#Psi_{V0C})#GT#GT",
                                  "#LT#LTcos(#psi_{POI_2}-2#Psi_{V0A})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#Psi_{V0A})#GT#GT"};
                                 
 TString nonIsotropicZDCTerm[12] = {"#LT#LTcos(#psi_{POI_1}-2#Psi_{ZDCC})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#Psi_{ZDCC})#GT#GT",
	                              "#LT#LTcos(#psi_{POI_1}-2#Psi_{ZDCA})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#Psi_{ZDCA})#GT#GT",
	                              "#LT#LTcos(#psi_{POI_1}-2#Psi_{ZDCCA})#GT#GT","#LT#LTsin(#psi_{POI_1}-2#Psi_{ZDCCA})#GT#GT",
                                  "#LT#LTcos(#psi_{POI_2}-2#Psi_{ZDCC})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#Psi_{ZDCC})#GT#GT",
                                  "#LT#LTcos(#psi_{POI_2}-2#Psi_{ZDCA})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#Psi_{ZDCA})#GT#GT",
                                  "#LT#LTcos(#psi_{POI_2}-2#Psi_{ZDCCA})#GT#GT","#LT#LTsin(#psi_{POI_2}-2#Psi_{ZDCCA})#GT#GT"};
                                 
 fNonIsotropicTermsV0EvPlPro = new TProfile("fNonIsotropicTermsV0EvPlPro","fNonIsotropicTermsV0EvPlPro",8,0.,8);
 fNonIsotropicTermsV0EvPlPro->GetXaxis()->SetBinLabel(1,"#LT#LTcos(#Psi_{V0C})#GT#GT");
 fNonIsotropicTermsV0EvPlPro->GetXaxis()->SetBinLabel(2,"#LT#LTsin(#Psi_{V0C})#GT#GT");
 fNonIsotropicTermsV0EvPlPro->GetXaxis()->SetBinLabel(3,"#LT#LTcos(#2Psi_{V0C})#GT#GT");
 fNonIsotropicTermsV0EvPlPro->GetXaxis()->SetBinLabel(4,"#LT#LTsin(#2Psi_{V0C})#GT#GT");
 fNonIsotropicTermsV0EvPlPro->GetXaxis()->SetBinLabel(5,"#LT#LTcos(#Psi_{V0A})#GT#GT");
 fNonIsotropicTermsV0EvPlPro->GetXaxis()->SetBinLabel(6,"#LT#LTsin(#Psi_{V0A})#GT#GT");
 fNonIsotropicTermsV0EvPlPro->GetXaxis()->SetBinLabel(7,"#LT#LTcos(#2Psi_{V0A})#GT#GT");
 fNonIsotropicTermsV0EvPlPro->GetXaxis()->SetBinLabel(8,"#LT#LTsin(#2Psi_{V0A})#GT#GT");
 fNonIsotropicTermsV0EvPlPro->SetStats(kFALSE);
 fNonIsotropicTermsList->Add(fNonIsotropicTermsV0EvPlPro);
 
 fNonIsotropicTermsZDCSpPlPPro = new TProfile("fNonIsotropicTermsZDCSpPlPPro","fNonIsotropicTermsZDCSpPlPPro",12,0.,12);
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(1,"#LT#LTcos(#Psi_{ZDCC})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(2,"#LT#LTsin(#Psi_{ZDCC})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(3,"#LT#LTcos(#2Psi_{ZDCC})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(4,"#LT#LTsin(#2Psi_{ZDCC})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(5,"#LT#LTcos(#Psi_{ZDCA})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(6,"#LT#LTsin(#Psi_{ZDCA})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(7,"#LT#LTcos(#2Psi_{ZDCA})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(8,"#LT#LTsin(#2Psi_{ZDCA})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(9,"#LT#LTcos(#Psi_{ZDCCA})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(10,"#LT#LTsin(#Psi_{ZDCCA})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(11,"#LT#LTcos(#2Psi_{ZDCCA})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->GetXaxis()->SetBinLabel(12,"#LT#LTsin(#2Psi_{ZDCCA})#GT#GT");
 fNonIsotropicTermsZDCSpPlPPro->SetStats(kFALSE);
 fNonIsotropicTermsList->Add(fNonIsotropicTermsZDCSpPlPPro);
 
 fNonIsotropicTermsV0Pro = new TProfile("fNonIsotropicTermsV0Pro","fNonIsotropicTermsV0Pro",4,0.,4);
 fNonIsotropicTermsV0Pro->GetXaxis()->SetBinLabel(1,"#LT#LTcos(#phi_{POI}-2#Psi_{V0C})#GT#GT");
 fNonIsotropicTermsV0Pro->GetXaxis()->SetBinLabel(2,"#LT#LTsin(#phi_{POI}-2#Psi_{V0C})#GT#GT");
 fNonIsotropicTermsV0Pro->GetXaxis()->SetBinLabel(3,"#LT#LTcos(#phi_{POI}-2#Psi_{V0A})#GT#GT");
 fNonIsotropicTermsV0Pro->GetXaxis()->SetBinLabel(4,"#LT#LTsin(#phi_{POI}-2#Psi_{V0A})#GT#GT");
 fNonIsotropicTermsV0Pro->SetStats(kFALSE);
 fNonIsotropicTermsList->Add(fNonIsotropicTermsV0Pro);
 
 fNonIsotropicTermsZDCPro = new TProfile("fNonIsotropicTermsZDCPro","fNonIsotropicTermsZDCPro",6,0.,6);
 fNonIsotropicTermsZDCPro->GetXaxis()->SetBinLabel(1,"#LT#LTcos(#phi_{POI}-2#Psi_{ZDCC})#GT#GT");
 fNonIsotropicTermsZDCPro->GetXaxis()->SetBinLabel(2,"#LT#LTsin(#phi_{POI}-2#Psi_{ZDCC})#GT#GT");
 fNonIsotropicTermsZDCPro->GetXaxis()->SetBinLabel(3,"#LT#LTcos(#phi_{POI}-2#Psi_{ZDCA})#GT#GT");
 fNonIsotropicTermsZDCPro->GetXaxis()->SetBinLabel(4,"#LT#LTsin(#phi_{POI}-2#Psi_{ZDCA})#GT#GT");
 fNonIsotropicTermsZDCPro->GetXaxis()->SetBinLabel(5,"#LT#LTcos(#phi_{POI}-2#Psi_{ZDCCA})#GT#GT");
 fNonIsotropicTermsZDCPro->GetXaxis()->SetBinLabel(6,"#LT#LTsin(#phi_{POI}-2#Psi_{ZDCCA})#GT#GT");
 fNonIsotropicTermsZDCPro->SetStats(kFALSE);
 fNonIsotropicTermsList->Add(fNonIsotropicTermsZDCPro);
  
 for(Int_t sd=0;sd<2;sd++)
 {
  for(Int_t t=0;t<10;t++)
  { 
   // Pt:
   fNonIsotropicTermsVsPtSumDiffPro[sd][t] = new TProfile(Form("fNonIsotropicTermsVs%sPro %s",psdFlag[sd].Data(),nonIsotropicTerm[t].Data()),"",fnBinsPt,0.,fPtMax);
   fNonIsotropicTermsVsPtSumDiffPro[sd][t]->SetTitle(Form("%s vs %s",nonIsotropicTerm[t].Data(),psdTitleFlag[sd].Data()));
   fNonIsotropicTermsVsPtSumDiffPro[sd][t]->SetStats(kFALSE);
   fNonIsotropicTermsVsPtSumDiffPro[sd][t]->GetXaxis()->SetTitle(psdTitleFlag[sd].Data());
   fNonIsotropicTermsList->Add(fNonIsotropicTermsVsPtSumDiffPro[sd][t]);
   // Eta:
   fNonIsotropicTermsVsEtaSumDiffPro[sd][t] = new TProfile(Form("fNonIsotropicTermsVs%sPro %s",psdFlag2[sd].Data(),nonIsotropicTerm[t].Data()),"",fnBinsEta,fEtaMin,fEtaMax);
   fNonIsotropicTermsVsEtaSumDiffPro[sd][t]->SetTitle(Form("%s vs %s",nonIsotropicTerm[t].Data(),psdTitleFlag2[sd].Data()));
   fNonIsotropicTermsVsEtaSumDiffPro[sd][t]->SetStats(kFALSE);
   fNonIsotropicTermsVsEtaSumDiffPro[sd][t]->GetXaxis()->SetTitle(psdTitleFlag2[sd].Data());
   fNonIsotropicTermsList->Add(fNonIsotropicTermsVsEtaSumDiffPro[sd][t]);
  } // end of for(Int_t t=0;t<10;t++)
 } // end of for(Int_t sd=0;sd<2;sd++)

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

void AliFlowAnalysisWithMixedHarmonics::AccessConstants(TString method)
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
   printf("\n WARNING (MH): fCommonConstants is NULL in AFAWMH::AC(\"%s\") !!!!\n\n",method.Data());
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

} // end of void AliFlowAnalysisWithMixedHarmonics::AccessConstants(TString method)

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
 if(!f2pCorrelatorCos2PsiDiff2PsiV0Pro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCos2PsiDiff2PsiV0Pro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f2pCorrelatorCos2PsiDiff2PsiZDCPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCos2PsiDiff2PsiZDCPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
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
 // non-isotropic terms for V0 and ZDC
 if(!(fNonIsotropicTermsV0EvPlPro))
 {
  cout<<endl;
  cout<<" WARNING (MH): "<<"fNonIsotropicTermsV0EvPlPro"<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);   
 } 
 if(!(fNonIsotropicTermsZDCSpPlPPro))
 {
  cout<<endl;
  cout<<" WARNING (MH): "<<"fNonIsotropicTermsZDCSpPlPPro"<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);   
 } 
 for(Int_t sd=0;sd<2;sd++)
 {
  if(!(f3pCorrelatorVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorV0CVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorV0CVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorV0AVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorV0AVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorZDCCVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorZDCCVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorZDCAVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorZDCAVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorZDCCAVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorZDCCAVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
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
  for(Int_t t=0;t<10;t++)
  { 
   if(!(fNonIsotropicTermsVsPtSumDiffPro[sd][t]))
   {
    cout<<endl;
    cout<<" WARNING (MH): "<<Form("fNonIsotropicTermsVsPtSumDiffPro[%d][%d]",sd,t)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
    cout<<endl;
    exit(0);   
   } 
   if(!(fNonIsotropicTermsVsEtaSumDiffPro[sd][t]))
   {
    cout<<endl;
    cout<<" WARNING (MH): "<<Form("fNonIsotropicTermsVsEtaSumDiffPro[%d][%d]",sd,t)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
    cout<<endl;
    exit(0);   
   } 
  } // end of for(Int_t t=0;t<10;t++) 
 } // end of for(Int_t sd=0;sd<2;sd++)
 
 if(!fNonIsotropicTermsV0Pro)
 {
  cout<<endl;
  cout<<" WARNING (MH): fNonIsotropicTermsV0Pro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 
 if(!fNonIsotropicTermsZDCPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): fNonIsotropicTermsZDCPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
  
 if(!f3pCorrelatorPOIIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIIntegratedPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIm2V0CIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2V0CIntegratedPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIm2V0AIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2V0AIntegratedPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIm2ZDCCIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2ZDCCIntegratedPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIm2ZDCAIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2ZDCAIntegratedPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIm2ZDCCAIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2ZDCCAIntegratedPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0); 
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
 for(Int_t p12=0;p12<2;p12++) // 1st/2nd POI 
 {
  for(Int_t ao=0;ao<2;ao++) // all/overlap
  {
   for(Int_t pe=0;pe<4;pe++) // [(p1+p2)/2,|p1-p2|,(eta1+eta2)/2,|eta1-eta2|]
   {
    if(!fReNITEBE[p12][ao][pe]||!fImNITEBE[p12][ao][pe])
    {
     cout<<endl;
     cout<<" WARNING (MH): "<<Form("!fReNITEBE[%d][%d][%d]||!fImNITEBE[%d][%d][%d]",p12,ao,pe,p12,ao,pe)<<" is NULL in CheckPointersUsedInMake() !!!!"<<endl;
     cout<<endl;
     exit(0);   
    }
   } // end of for(Int_t pe=0;pe<4;pe++) // [(p1+p2)/2,|p1-p2|,(eta1+eta2)/2,|eta1-eta2|]
  } // end of for(Int_t ao=0;ao<2;ao++) // all/overlap
 } // end of for(Int_t p12=0;p12<2;p12++) // 1st/2nd POI 
  
 //2p correlator <<cos(psi1-psi2)>> and <<cos(2*(psi1 - psi2))>> profile
 if(!fNonIsotropicPOITermsPro) {
  cout<<endl;
  cout<<" WARNING (MH): fNonIsotropicPOITermsPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 
 //2p correlator <<cos(psi1-psi2)>> and <<cos(2*(psi1 - psi2))>> profile
 if(!f2pCorrelatorCosPsiDiffPro) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCos2PsiDiffPro) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCos2PsiDiffPro is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 
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
 
 //QA for V0 and ZDC event plane
 if(!fQAZDCCAvsV0C) {
  cout<<endl;
  cout<<" WARNING (MH): fQAZDCCAvsV0C is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fQAZDCCAvsV0A) {
  cout<<endl;
  cout<<" WARNING (MH): fQAZDCCAvsV0A is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fQAV0CEvPldistribution) {
  cout<<endl;
  cout<<" WARNING (MH): fQAV0CEvPldistribution is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fQAV0AEvPldistribution) {
  cout<<endl;
  cout<<" WARNING (MH): fQAV0AEvPldistribution is NULL in CheckPointersUsedInMake() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fQAZDCCAEvPldistribution) {
  cout<<endl;
  cout<<" WARNING (MH): fQAZDCCAEvPldistribution is NULL in CheckPointersUsedInMake() !!!!"<<endl;
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
 if(!f2pCorrelatorCos2PsiDiff2PsiV0Pro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCos2PsiDiff2PsiV0Pro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f2pCorrelatorCos2PsiDiff2PsiZDCPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCos2PsiDiff2PsiZDCPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
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
 
 cout<<"======> in Finish() fEvaluateDifferential3pCorrelator = "<<fEvaluateDifferential3pCorrelator<<endl;
 fEvaluateDifferential3pCorrelator = kTRUE;
 
 if(!fEvaluateDifferential3pCorrelator){return;} 
 // non-isotropic terms for V0 and ZDC
 if(!(fNonIsotropicTermsV0EvPlPro))
 {
  cout<<endl;
  cout<<" WARNING (MH): "<<"fNonIsotropicTermsV0EvPlPro"<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);   
 } 
 if(!(fNonIsotropicTermsZDCSpPlPPro))
 {
  cout<<endl;
  cout<<" WARNING (MH): "<<"fNonIsotropicTermsZDCSpPlPPro"<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);   
 } 
 
 for(Int_t sd=0;sd<2;sd++)
 {
  if(!(f3pCorrelatorVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorV0CVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorV0CVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorV0AVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorV0AVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorZDCCVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorZDCCVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorZDCAVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorZDCAVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorZDCCAVsPtSumDiffPro[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorZDCCAVsPtSumDiffPro[%d]",sd)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
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
  if(!(f3pCorrelatorVsPtSumDiffHist[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorVsPtSumDiffHist[%d]",sd)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  } 
  if(!(f3pCorrelatorVsEtaSumDiffHist[sd]))
  {
   cout<<endl;
   cout<<" WARNING (MH): "<<Form("f3pCorrelatorVsEtaSumDiffHist[%d]",sd)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
   cout<<endl;
   exit(0);   
  }
  for(Int_t t=0;t<10;t++)
  { 
   if(!(fNonIsotropicTermsVsPtSumDiffPro[sd][t]))
   {
    cout<<endl;
    cout<<" WARNING (MH): "<<Form("fNonIsotropicTermsVsPtSumDiffPro[%d][%d]",sd,t)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
    cout<<endl;
    exit(0);   
   } 
   if(!(fNonIsotropicTermsVsEtaSumDiffPro[sd][t]))
   {
    cout<<endl;
    cout<<" WARNING (MH): "<<Form("fNonIsotropicTermsVsEtaSumDiffPro[%d][%d]",sd,t)<<" is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
    cout<<endl;
    exit(0);   
   } 
  } // end of for(Int_t t=0;t<10;t++) 
 } // end of for(Int_t sd=0;sd<2;sd++)
 
 if(!fNonIsotropicTermsV0Pro)
 {
  cout<<endl;
  cout<<" WARNING (MH): fNonIsotropicTermsV0Pro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 
 if(!fNonIsotropicTermsZDCPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): fNonIsotropicTermsZDCPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 
 if(!f3pCorrelatorPOIIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIIntegratedPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIm2V0CIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2V0CIntegratedPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIm2V0AIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2V0AIntegratedPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIm2ZDCCIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2ZDCCIntegratedPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIm2ZDCAIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2ZDCAIntegratedPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIm2ZDCCAIntegratedPro)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2ZDCCAIntegratedPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 }
 if(!f3pCorrelatorPOIIntegratedHist)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIIntegratedHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 

 if(!f3pCorrelatorPOIm2V0CIntegratedHist)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2V0CIntegratedHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!f3pCorrelatorPOIm2V0AIntegratedHist)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2V0AIntegratedHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!f3pCorrelatorPOIm2ZDCCIntegratedHist)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2ZDCCIntegratedHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!f3pCorrelatorPOIm2ZDCAIntegratedHist)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2ZDCAIntegratedHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 if(!f3pCorrelatorPOIm2ZDCCAIntegratedHist)
 {
  cout<<endl;
  cout<<" WARNING (MH): f3pCorrelatorPOIm2ZDCCAIntegratedHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 //2p correlator <<cos(psi1-psi2)>> and <<cos(2*(psi1 - psi2))>> profile
 if(!fNonIsotropicPOITermsPro) {
  cout<<endl;
  cout<<" WARNING (MH): fNonIsotropicPOITermsPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 
 //2p correlator <<cos(psi1-psi2)>> and <<cos(2*(psi1 - psi2))>> profile
 if(!f2pCorrelatorCosPsiDiffPro) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCos2PsiDiffPro) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCos2PsiDiffPro is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 
 //2p correlator <<cos(psi1-psi2)>> and <<cos(2*(psi1 - psi2))>> hist
 if(!f2pCorrelatorCosPsiDiffHist) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCosPsiDiffHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!f2pCorrelatorCos2PsiDiffHist) {
  cout<<endl;
  cout<<" WARNING (MH): f2pCorrelatorCos2PsiDiffHist is NULL in CheckPointersUsedInFinish() !!!!"<<endl;
  cout<<endl;
  exit(0);
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
 // a.) Correct integrated 3-p correlator cos[n(phi1+phi2-2phi3)] for detector effects;
 // b.) Correct differential 3-p correlator cos[n(psi1+psi2-2phi3)] for detector effects.
 // c.) Correct 2-p correlator cos[2(phi1-phi2)] for detector effects. This is needed for DeltaGamma/v2
 
 // a.) Correct integrated 3-p correlator cos[n(phi1+phi2-2phi3)] for detector effects:  
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
 
 if(!fEvaluateDifferential3pCorrelator){return;}
 // b.) Correct differential 3-p correlator cos[n(psi1+psi2-2phi3)] for detector effects:
 Int_t iBinCounter = 0;
   Double_t gSumBinContentTimesWeight = 0., gSumWeight = 0.;
   Double_t gSumBinContentTimesWeightSquared = 0.;
   
   Double_t gIntegratedValue = -1000.;
   Double_t gIntegratedValueV0C = -1000.;
   Double_t gIntegratedValueV0A = -1000.;
   Double_t gIntegratedValueZDCC = -1000.;
   Double_t gIntegratedValueZDCA = -1000.;
   Double_t gIntegratedValueZDCCA = -1000.;
   Double_t gIntegratedValueZDCCAoverCorrect = -1000.;
   
 for(Int_t sd=0;sd<2;sd++) 
 {
  // [(p1+p2)/2,|p1-p2|]
  // looping over all bins and calculating reduced correlations: 
  for(Int_t b=1;b<=fnBinsPt;b++)
  {
   Double_t measured = f3pCorrelatorVsPtSumDiffPro[sd]->GetBinContent(b); // cos[n(psi1+psi2-2phi3)]
   Double_t measuredErr = f3pCorrelatorVsPtSumDiffPro[sd]->GetBinError(b); 
  
   Double_t measuredZDCCA = f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->GetBinContent(b); // cos[n(psi1+psi2-2phi_ZDCCA)]
   Double_t measuredZDCCAErr = f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->GetBinError(b);
	
   Double_t corrected = 0.; // 3-p correlator corrected for detector effects
   Double_t correctedErr = measuredErr; // to be improved - propagate error also for non-isotropic terms
   
   
   // non-isotropic terms:
   Double_t cosPsiPOI1 = fNonIsotropicTermsVsPtSumDiffPro[sd][0]->GetBinContent(b); // <<cos(#psi_{POI_1})>>
   Double_t sinPsiPOI1 = fNonIsotropicTermsVsPtSumDiffPro[sd][1]->GetBinContent(b); // <<sin(#psi_{POI_1})>>
   Double_t cosPsiPOI2 = fNonIsotropicTermsVsPtSumDiffPro[sd][2]->GetBinContent(b); // <<cos(#psi_{POI_2})>>
   Double_t sinPsiPOI2 = fNonIsotropicTermsVsPtSumDiffPro[sd][3]->GetBinContent(b); // <<sin(#psi_{POI_2})>>
   Double_t cosPsiPOI1m2PhiRP = fNonIsotropicTermsVsPtSumDiffPro[sd][4]->GetBinContent(b); // <<cos(#psi_{POI_1}-2*phi_{RP})>>
   Double_t sinPsiPOI1m2PhiRP = fNonIsotropicTermsVsPtSumDiffPro[sd][5]->GetBinContent(b); // <<sin(#psi_{POI_1}-2*phi_{RP})>>
   Double_t cosPsiPOI2m2PhiRP = fNonIsotropicTermsVsPtSumDiffPro[sd][6]->GetBinContent(b); // <<cos(#psi_{POI_2}-2*phi_{RP})>>
   Double_t sinPsiPOI2m2PhiRP = fNonIsotropicTermsVsPtSumDiffPro[sd][7]->GetBinContent(b); // <<sin(#psi_{POI_2}-2*phi_{RP})>>
   Double_t cosPsiPOI1pPsiPOI2 = fNonIsotropicTermsVsPtSumDiffPro[sd][8]->GetBinContent(b); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
   Double_t sinPsiPOI1pPsiPOI2 = fNonIsotropicTermsVsPtSumDiffPro[sd][9]->GetBinContent(b); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>
   Double_t cos2PhiRP = fNonIsotropicTermsPro->GetBinContent(3); // <<cos(2n*(phi_{RP}))>>
   Double_t sin2PhiRP = fNonIsotropicTermsPro->GetBinContent(4); // <<sin(2n*(phi_{RP}))>>
   
   
	   
   // apply correction:
   corrected = measured 
             - (cosPsiPOI1*cosPsiPOI2m2PhiRP-sinPsiPOI1*sinPsiPOI2m2PhiRP 
             + cosPsiPOI2*cosPsiPOI1m2PhiRP-sinPsiPOI2*sinPsiPOI1m2PhiRP    
             + cos2PhiRP*cosPsiPOI1pPsiPOI2+sin2PhiRP*sinPsiPOI1pPsiPOI2)
             + 2.*cos2PhiRP*(cosPsiPOI1*cosPsiPOI2-sinPsiPOI1*sinPsiPOI2)
             + 2.*sin2PhiRP*(cosPsiPOI1*sinPsiPOI2+sinPsiPOI1*cosPsiPOI2); 
             

                
   if(fCorrectForDetectorEffects)
   {
    f3pCorrelatorVsPtSumDiffHist[sd]->SetBinContent(b,corrected);
    f3pCorrelatorVsPtSumDiffHist[sd]->SetBinError(b,correctedErr);
   }
   
   if(sd == 0) {
      iBinCounter += 1;

      gSumBinContentTimesWeight += f3pCorrelatorVsPtSumDiffHist[sd]->GetBinContent(b)*f3pCorrelatorVsPtSumDiffPro[sd]->GetBinEntries(b);
      gSumWeight += f3pCorrelatorVsPtSumDiffPro[sd]->GetBinEntries(b);
      gSumBinContentTimesWeightSquared += TMath::Power(gSumBinContentTimesWeight,2);
   }
  } // end of for(Int_t b=1;b<=fnBinsPt;b++)

  
  // [(eta1+eta2)/2,|eta1-eta2|]
  // looping over all bins and calculating reduced correlations: 
  for(Int_t b=1;b<=fnBinsEta;b++)
  {
   Double_t measured = f3pCorrelatorVsEtaSumDiffPro[sd]->GetBinContent(b);
   Double_t measuredErr = f3pCorrelatorVsEtaSumDiffPro[sd]->GetBinError(b);   
   Double_t corrected = 0.; // 3-p correlator corrected for detector effects
   Double_t correctedErr = measuredErr; // to be improved - propagate error also for non-isotropic terms
   // non-isotropic terms:
   Double_t cosPsiPOI1 = fNonIsotropicTermsVsEtaSumDiffPro[sd][0]->GetBinContent(b); // <<cos(#psi_{POI_1})>>
   Double_t sinPsiPOI1 = fNonIsotropicTermsVsEtaSumDiffPro[sd][1]->GetBinContent(b); // <<sin(#psi_{POI_1})>>
   Double_t cosPsiPOI2 = fNonIsotropicTermsVsEtaSumDiffPro[sd][2]->GetBinContent(b); // <<cos(#psi_{POI_2})>>
   Double_t sinPsiPOI2 = fNonIsotropicTermsVsEtaSumDiffPro[sd][3]->GetBinContent(b); // <<sin(#psi_{POI_2})>>
   Double_t cosPsiPOI1m2PhiRP = fNonIsotropicTermsVsEtaSumDiffPro[sd][4]->GetBinContent(b); // <<cos(#psi_{POI_1}-2*phi_{RP})>>
   Double_t sinPsiPOI1m2PhiRP = fNonIsotropicTermsVsEtaSumDiffPro[sd][5]->GetBinContent(b); // <<sin(#psi_{POI_1}-2*phi_{RP})>>
   Double_t cosPsiPOI2m2PhiRP = fNonIsotropicTermsVsEtaSumDiffPro[sd][6]->GetBinContent(b); // <<cos(#psi_{POI_2}-2*phi_{RP})>>
   Double_t sinPsiPOI2m2PhiRP = fNonIsotropicTermsVsEtaSumDiffPro[sd][7]->GetBinContent(b); // <<sin(#psi_{POI_2}-2*phi_{RP})>>
   Double_t cosPsiPOI1pPsiPOI2 = fNonIsotropicTermsVsEtaSumDiffPro[sd][8]->GetBinContent(b); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
   Double_t sinPsiPOI1pPsiPOI2 = fNonIsotropicTermsVsEtaSumDiffPro[sd][9]->GetBinContent(b); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>
   Double_t cos2PhiRP = fNonIsotropicTermsPro->GetBinContent(3); // <<cos(2n*(phi_{RP}))>>
   Double_t sin2PhiRP = fNonIsotropicTermsPro->GetBinContent(4); // <<sin(2n*(phi_{RP}))>>
   // apply correction:
   corrected = measured 
             - (cosPsiPOI1*cosPsiPOI2m2PhiRP-sinPsiPOI1*sinPsiPOI2m2PhiRP 
             + cosPsiPOI2*cosPsiPOI1m2PhiRP-sinPsiPOI2*sinPsiPOI1m2PhiRP    
             + cos2PhiRP*cosPsiPOI1pPsiPOI2+sin2PhiRP*sinPsiPOI1pPsiPOI2)
             + 2.*cos2PhiRP*(cosPsiPOI1*cosPsiPOI2-sinPsiPOI1*sinPsiPOI2)
             + 2.*sin2PhiRP*(cosPsiPOI1*sinPsiPOI2+sinPsiPOI1*cosPsiPOI2);
   if(fCorrectForDetectorEffects)
   {
    f3pCorrelatorVsEtaSumDiffHist[sd]->SetBinContent(b,corrected);
    f3pCorrelatorVsEtaSumDiffHist[sd]->SetBinError(b,correctedErr);
   }
  } // end of for(Int_t b=1;b<=fnBinsEta;b++)
 } // end of for(Int_t sd=0;sd<2;sd++) 
 
 
 if (kTRUE) {
  Double_t gSumBinContentTimesWeightV0C = 0., gSumWeightV0C = 0.;
   Double_t gSumBinContentTimesWeightV0CErr = 0.;
   
   Double_t gSumBinContentTimesWeightV0A = 0., gSumWeightV0A = 0.;
   Double_t gSumBinContentTimesWeightV0AErr = 0.;
   
   Double_t gSumBinContentTimesWeightZDCC = 0., gSumWeightZDCC = 0.;
   Double_t gSumBinContentTimesWeightZDCCErr = 0.;
   
   Double_t gSumBinContentTimesWeightZDCA = 0., gSumWeightZDCA = 0.;
   Double_t gSumBinContentTimesWeightZDCAErr = 0.;
   
   Double_t gSumBinContentTimesWeightZDCCA = 0., gSumWeightZDCCA = 0.;
   Double_t gSumBinContentTimesWeightZDCCAErr = 0.;
  
  Double_t cosPsiPOI1 = 0; // <<cos(#psi_{POI_1})>>
  Double_t sinPsiPOI1 = 0; // <<sin(#psi_{POI_1})>>
  Double_t cosPsiPOI2 = 0; // <<cos(#psi_{POI_2})>>
  Double_t sinPsiPOI2 = 0; // <<sin(#psi_{POI_2})>>
  Double_t cosPsiPOI1m2PhiRP = 0; // <<cos(#psi_{POI_1}-2*phi_{RP})>>
  Double_t sinPsiPOI1m2PhiRP = 0; // <<sin(#psi_{POI_1}-2*phi_{RP})>>
  Double_t cosPsiPOI2m2PhiRP = 0; // <<cos(#psi_{POI_2}-2*phi_{RP})>>
  Double_t sinPsiPOI2m2PhiRP = 0; // <<sin(#psi_{POI_2}-2*phi_{RP})>>
  Double_t cosPsiPOI1pPsiPOI2 = 0; // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
  Double_t sinPsiPOI1pPsiPOI2 = 0; // <<sin(#psi_{POI_1}+#psi_{POI_2})>>
  Double_t cos2PhiRP = 0; // <<cos(2n*(phi_{RP}))>>
  Double_t sin2PhiRP = 0; // <<sin(2n*(phi_{RP}))>>
  
  Double_t cosPsiPOI1SumWeight = 0; // <<cos(#psi_{POI_1})>>
  Double_t sinPsiPOI1SumWeight = 0; // <<sin(#psi_{POI_1})>>
  Double_t cosPsiPOI2SumWeight = 0; // <<cos(#psi_{POI_2})>>
  Double_t sinPsiPOI2SumWeight = 0; // <<sin(#psi_{POI_2})>>
  Double_t cosPsiPOI1m2PhiRPSumWeight = 0; // <<cos(#psi_{POI_1}-2*phi_{RP})>>
  Double_t sinPsiPOI1m2PhiRPSumWeight = 0; // <<sin(#psi_{POI_1}-2*phi_{RP})>>
  Double_t cosPsiPOI2m2PhiRPSumWeight = 0; // <<cos(#psi_{POI_2}-2*phi_{RP})>>
  Double_t sinPsiPOI2m2PhiRPSumWeight = 0; // <<sin(#psi_{POI_2}-2*phi_{RP})>>
  Double_t cosPsiPOI1pPsiPOI2SumWeight = 0; // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
  Double_t sinPsiPOI1pPsiPOI2SumWeight = 0; // <<sin(#psi_{POI_1}+#psi_{POI_2})>>
   
   
  Double_t cosPsiPOIm2PhiV0C = 0, sinPsiPOIm2PhiV0C = 0, cosPsiPOIm2PhiV0A = 0, sinPsiPOIm2PhiV0A = 0;
  
  Double_t cosPsiPOIm2PhiZDCC = 0, sinPsiPOIm2PhiZDCC = 0, cosPsiPOIm2PhiZDCA = 0, sinPsiPOIm2PhiZDCA = 0, cosPsiPOIm2PhiZDCCA = 0, sinPsiPOIm2PhiZDCCA = 0;
  
  
  // V0 and ZDC non isotropic terms do not fit for the same binning structure
  for(Int_t b=1;b<=fnBinsPt;b++)
  {
	  Int_t sd = 0;
	  gSumBinContentTimesWeightV0C += f3pCorrelatorV0CVsPtSumDiffPro[sd]->GetBinContent(b)*f3pCorrelatorV0CVsPtSumDiffPro[sd]->GetBinEntries(b); // cos[n(psi1+psi2-2phi_V0C)]
	  gSumWeightV0C += f3pCorrelatorV0CVsPtSumDiffPro[sd]->GetBinEntries(b);
	  gSumBinContentTimesWeightV0CErr += f3pCorrelatorV0CVsPtSumDiffPro[sd]->GetBinError(b)*f3pCorrelatorV0CVsPtSumDiffPro[sd]->GetBinEntries(b);
	  
      gSumBinContentTimesWeightV0A += f3pCorrelatorV0AVsPtSumDiffPro[sd]->GetBinContent(b)*f3pCorrelatorV0AVsPtSumDiffPro[sd]->GetBinEntries(b); // cos[n(psi1+psi2-2phi_V0A)]
      gSumWeightV0A += f3pCorrelatorV0AVsPtSumDiffPro[sd]->GetBinEntries(b);
	  gSumBinContentTimesWeightV0AErr += f3pCorrelatorV0AVsPtSumDiffPro[sd]->GetBinError(b)*f3pCorrelatorV0AVsPtSumDiffPro[sd]->GetBinEntries(b);
	  
	  gSumBinContentTimesWeightZDCC += f3pCorrelatorZDCCVsPtSumDiffPro[sd]->GetBinContent(b)*f3pCorrelatorZDCCVsPtSumDiffPro[sd]->GetBinEntries(b); // cos[n(psi1+psi2-2phi_ZDCC)]
      gSumWeightZDCC += f3pCorrelatorZDCCVsPtSumDiffPro[sd]->GetBinEntries(b);
	  gSumBinContentTimesWeightZDCCErr += f3pCorrelatorZDCCVsPtSumDiffPro[sd]->GetBinError(b)*f3pCorrelatorZDCCVsPtSumDiffPro[sd]->GetBinEntries(b);
	  
	  gSumBinContentTimesWeightZDCA += f3pCorrelatorZDCAVsPtSumDiffPro[sd]->GetBinContent(b)*f3pCorrelatorZDCAVsPtSumDiffPro[sd]->GetBinEntries(b); // cos[n(psi1+psi2-2phi_ZDCA)]
      gSumWeightZDCA += f3pCorrelatorZDCAVsPtSumDiffPro[sd]->GetBinEntries(b);
	  gSumBinContentTimesWeightZDCAErr += f3pCorrelatorZDCAVsPtSumDiffPro[sd]->GetBinError(b)*f3pCorrelatorZDCAVsPtSumDiffPro[sd]->GetBinEntries(b);
	  
	  gSumBinContentTimesWeightZDCCA += f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->GetBinContent(b)*f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->GetBinEntries(b); // cos[n(psi1+psi2-2phi_ZDCCA)]
      gSumWeightZDCCA += f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->GetBinEntries(b);
	  gSumBinContentTimesWeightZDCCAErr += f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->GetBinError(b)*f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->GetBinEntries(b);
	  
      

      cosPsiPOI1 += fNonIsotropicTermsVsPtSumDiffPro[sd][0]->GetBinContent(b)*fNonIsotropicTermsVsPtSumDiffPro[sd][0]->GetBinEntries(b); // <<cos(#psi_{POI_1})>>
      sinPsiPOI1 += fNonIsotropicTermsVsPtSumDiffPro[sd][1]->GetBinContent(b)*fNonIsotropicTermsVsPtSumDiffPro[sd][1]->GetBinEntries(b); // <<sin(#psi_{POI_1})>>
      cosPsiPOI2 += fNonIsotropicTermsVsPtSumDiffPro[sd][2]->GetBinContent(b)*fNonIsotropicTermsVsPtSumDiffPro[sd][2]->GetBinEntries(b); // <<cos(#psi_{POI_2})>>
      sinPsiPOI2 += fNonIsotropicTermsVsPtSumDiffPro[sd][3]->GetBinContent(b)*fNonIsotropicTermsVsPtSumDiffPro[sd][3]->GetBinEntries(b); // <<sin(#psi_{POI_2})>>
      cosPsiPOI1m2PhiRP += fNonIsotropicTermsVsPtSumDiffPro[sd][4]->GetBinContent(b)*fNonIsotropicTermsVsPtSumDiffPro[sd][4]->GetBinEntries(b); // <<cos(#psi_{POI_1}-2*phi_{RP})>>
      sinPsiPOI1m2PhiRP += fNonIsotropicTermsVsPtSumDiffPro[sd][5]->GetBinContent(b)*fNonIsotropicTermsVsPtSumDiffPro[sd][5]->GetBinEntries(b); // <<sin(#psi_{POI_1}-2*phi_{RP})>>
      cosPsiPOI2m2PhiRP += fNonIsotropicTermsVsPtSumDiffPro[sd][6]->GetBinContent(b)*fNonIsotropicTermsVsPtSumDiffPro[sd][6]->GetBinEntries(b); // <<cos(#psi_{POI_2}-2*phi_{RP})>>
      sinPsiPOI2m2PhiRP += fNonIsotropicTermsVsPtSumDiffPro[sd][7]->GetBinContent(b)*fNonIsotropicTermsVsPtSumDiffPro[sd][7]->GetBinEntries(b); // <<sin(#psi_{POI_2}-2*phi_{RP})>>
      cosPsiPOI1pPsiPOI2 += fNonIsotropicTermsVsPtSumDiffPro[sd][8]->GetBinContent(b)*fNonIsotropicTermsVsPtSumDiffPro[sd][8]->GetBinEntries(b); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
      sinPsiPOI1pPsiPOI2 += fNonIsotropicTermsVsPtSumDiffPro[sd][9]->GetBinContent(b)*fNonIsotropicTermsVsPtSumDiffPro[sd][9]->GetBinEntries(b); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>
      cos2PhiRP = fNonIsotropicTermsPro->GetBinContent(3); // <<cos(2n*(phi_{RP}))>>
      sin2PhiRP = fNonIsotropicTermsPro->GetBinContent(4); // <<sin(2n*(phi_{RP}))>>
      
      cosPsiPOI1SumWeight += fNonIsotropicTermsVsPtSumDiffPro[sd][0]->GetBinEntries(b); // <<cos(#psi_{POI_1})>>
      sinPsiPOI1SumWeight += fNonIsotropicTermsVsPtSumDiffPro[sd][1]->GetBinEntries(b); // <<sin(#psi_{POI_1})>>
      cosPsiPOI2SumWeight += fNonIsotropicTermsVsPtSumDiffPro[sd][2]->GetBinEntries(b); // <<cos(#psi_{POI_2})>>
      sinPsiPOI2SumWeight += fNonIsotropicTermsVsPtSumDiffPro[sd][3]->GetBinEntries(b); // <<sin(#psi_{POI_2})>>
      cosPsiPOI1m2PhiRPSumWeight += fNonIsotropicTermsVsPtSumDiffPro[sd][4]->GetBinEntries(b); // <<cos(#psi_{POI_1}-2*phi_{RP})>>
      sinPsiPOI1m2PhiRPSumWeight += fNonIsotropicTermsVsPtSumDiffPro[sd][5]->GetBinEntries(b); // <<sin(#psi_{POI_1}-2*phi_{RP})>>
      cosPsiPOI2m2PhiRPSumWeight += fNonIsotropicTermsVsPtSumDiffPro[sd][6]->GetBinEntries(b); // <<cos(#psi_{POI_2}-2*phi_{RP})>>
      sinPsiPOI2m2PhiRPSumWeight += fNonIsotropicTermsVsPtSumDiffPro[sd][7]->GetBinEntries(b); // <<sin(#psi_{POI_2}-2*phi_{RP})>>
      cosPsiPOI1pPsiPOI2SumWeight += fNonIsotropicTermsVsPtSumDiffPro[sd][8]->GetBinEntries(b); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
      sinPsiPOI1pPsiPOI2SumWeight += fNonIsotropicTermsVsPtSumDiffPro[sd][9]->GetBinEntries(b); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>
      
      
	  
  }
  
  Double_t measuredV0C = gSumBinContentTimesWeightV0C/gSumWeightV0C;
  //Double_t measuredV0CErr = 0;
  Double_t measuredV0A = gSumBinContentTimesWeightV0A/gSumWeightV0A;
  //Double_t measuredV0AErr = 0; 
  Double_t measuredZDCC = gSumBinContentTimesWeightZDCC/gSumWeightZDCC;
  //Double_t measuredZDCCErr = 0; 
  Double_t measuredZDCA = gSumBinContentTimesWeightZDCA/gSumWeightZDCA;
  //Double_t measuredZDCAErr = 0; 
  Double_t measuredZDCCA = gSumBinContentTimesWeightZDCCA/gSumWeightZDCCA;
  //Double_t measuredZDCCAErr = 0;
  
  Double_t correctedV0C = 0.; // 3-p correlator corrected for detector effects
  //Double_t correctedV0CErr = gSumBinContentTimesWeightV0CErr/gSumWeightV0C; // to be improved - propagate error also for non-isotropic terms
  Double_t correctedV0A = 0.; // 3-p correlator corrected for detector effects
  //Double_t correctedV0AErr = gSumBinContentTimesWeightV0AErr/gSumWeightV0A; // to be improved - propagate error also for non-isotropic terms
  Double_t correctedZDCC = 0.; // 3-p correlator corrected for detector effects
  //Double_t correctedZDCCErr = gSumBinContentTimesWeightZDCCErr/gSumWeightZDCC; // to be improved - propagate error also for non-isotropic terms
  Double_t correctedZDCA = 0.; // 3-p correlator corrected for detector effects
  //Double_t correctedZDCAErr = gSumBinContentTimesWeightZDCAErr/gSumWeightZDCA; // to be improved - propagate error also for non-isotropic terms
  Double_t correctedZDCCA = 0.; // 3-p correlator corrected for detector effects
  //Double_t correctedZDCCAErr = gSumBinContentTimesWeightZDCCAErr/gSumWeightZDCCA; // to be improved - propagate error also for non-isotropic terms
  Double_t correctedZDCCAoverCorrect = 0;
  
  cosPsiPOI1 = cosPsiPOI1/cosPsiPOI1SumWeight; // <<cos(#psi_{POI_1})>>
  sinPsiPOI1 = sinPsiPOI1/sinPsiPOI1SumWeight; // <<sin(#psi_{POI_1})>>
  cosPsiPOI2 = cosPsiPOI2/cosPsiPOI2SumWeight; // <<cos(#psi_{POI_2})>>
  sinPsiPOI2 = sinPsiPOI2/sinPsiPOI2SumWeight; // <<sin(#psi_{POI_2})>>
  cosPsiPOI1m2PhiRP = cosPsiPOI1m2PhiRP/cosPsiPOI1m2PhiRPSumWeight; // <<cos(#psi_{POI_1}-2*phi_{RP})>>
  sinPsiPOI1m2PhiRP = sinPsiPOI1m2PhiRP/sinPsiPOI1m2PhiRPSumWeight; // <<sin(#psi_{POI_1}-2*phi_{RP})>>
  cosPsiPOI2m2PhiRP = cosPsiPOI2m2PhiRP/cosPsiPOI2m2PhiRPSumWeight; // <<cos(#psi_{POI_2}-2*phi_{RP})>>
  sinPsiPOI2m2PhiRP = sinPsiPOI2m2PhiRP/sinPsiPOI2m2PhiRPSumWeight; // <<sin(#psi_{POI_2}-2*phi_{RP})>>
  cosPsiPOI1pPsiPOI2 = cosPsiPOI1pPsiPOI2/cosPsiPOI1pPsiPOI2SumWeight; // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
  sinPsiPOI1pPsiPOI2 = sinPsiPOI1pPsiPOI2/sinPsiPOI1pPsiPOI2SumWeight; // <<sin(#psi_{POI_1}+#psi_{POI_2})>>
  
  cosPsiPOIm2PhiV0C = fNonIsotropicTermsV0Pro->GetBinContent(1); //<<cos(psi_{POI}-2phi_V0C)>>, POI1 is same as POI2
  sinPsiPOIm2PhiV0C = fNonIsotropicTermsV0Pro->GetBinContent(2); //<<sin(psi_{POI}-2phi_V0C)>>, POI1 is same as POI2
  cosPsiPOIm2PhiV0A = fNonIsotropicTermsV0Pro->GetBinContent(3); //<<cos(psi_{POI}-2phi_V0A)>>, POI1 is same as POI2
  sinPsiPOIm2PhiV0A = fNonIsotropicTermsV0Pro->GetBinContent(4); //<<sin(psi_{POI}-2phi_V0A)>>, POI1 is same as POI2
      
  cosPsiPOIm2PhiZDCC = fNonIsotropicTermsZDCPro->GetBinContent(1); //<<cos(psi_{POI}-2phi_ZDCC)>>, POI1 is same as POI2
  sinPsiPOIm2PhiZDCC = fNonIsotropicTermsZDCPro->GetBinContent(2); //<<sin(psi_{POI}-2phi_ZDCC)>>, POI1 is same as POI2
  cosPsiPOIm2PhiZDCA = fNonIsotropicTermsZDCPro->GetBinContent(3); //<<cos(psi_{POI}-2phi_ZDCA)>>, POI1 is same as POI2
  sinPsiPOIm2PhiZDCA = fNonIsotropicTermsZDCPro->GetBinContent(4); //<<sin(psi_{POI}-2phi_ZDCA)>>, POI1 is same as POI2
  cosPsiPOIm2PhiZDCCA = fNonIsotropicTermsZDCPro->GetBinContent(5); //<<cos(psi_{POI}-2phi_ZDCCA)>>, POI1 is same as POI2
  sinPsiPOIm2PhiZDCCA = fNonIsotropicTermsZDCPro->GetBinContent(6); //<<sin(psi_{POI}-2phi_ZDCCA)>>, POI1 is same as POI2
  
  Double_t cos2PhiZDCCA = fNonIsotropicTermsZDCSpPlPPro->GetBinContent(11); // <<cos(2#psi_{ZDCCA})>>
  Double_t sin2PhiZDCCA = fNonIsotropicTermsZDCSpPlPPro->GetBinContent(12); // <<sin(2#psi_{ZDCCA})>>
  
  correctedV0C = measuredV0C 
                - (cosPsiPOI1*cosPsiPOIm2PhiV0C-sinPsiPOI1*sinPsiPOIm2PhiV0C 
                + cosPsiPOI2*cosPsiPOIm2PhiV0C-sinPsiPOI2*sinPsiPOIm2PhiV0C); 
                
  correctedV0A = measuredV0A 
                - (cosPsiPOI1*cosPsiPOIm2PhiV0A-sinPsiPOI1*sinPsiPOIm2PhiV0A 
                + cosPsiPOI2*cosPsiPOIm2PhiV0A-sinPsiPOI2*sinPsiPOIm2PhiV0A); 
                
  correctedZDCC = measuredZDCC 
                - (cosPsiPOI1*cosPsiPOIm2PhiZDCC-sinPsiPOI1*sinPsiPOIm2PhiZDCC 
                + cosPsiPOI2*cosPsiPOIm2PhiZDCC-sinPsiPOI2*sinPsiPOIm2PhiZDCC); 
  correctedZDCA = measuredZDCA 
                - (cosPsiPOI1*cosPsiPOIm2PhiZDCA-sinPsiPOI1*sinPsiPOIm2PhiZDCA 
                + cosPsiPOI2*cosPsiPOIm2PhiZDCA-sinPsiPOI2*sinPsiPOIm2PhiZDCA); 
  correctedZDCCA = measuredZDCCA 
                - (cosPsiPOI1*cosPsiPOIm2PhiZDCCA-sinPsiPOI1*sinPsiPOIm2PhiZDCCA 
                + cosPsiPOI2*cosPsiPOIm2PhiZDCCA-sinPsiPOI2*sinPsiPOIm2PhiZDCCA); 
                 
  correctedZDCCAoverCorrect = measuredZDCCA 
                - (cosPsiPOI1*cosPsiPOIm2PhiZDCCA-sinPsiPOI1*sinPsiPOIm2PhiZDCCA 
                + cosPsiPOI2*cosPsiPOIm2PhiZDCCA-sinPsiPOI2*sinPsiPOIm2PhiZDCCA
                + cos2PhiZDCCA*cosPsiPOI1pPsiPOI2+sin2PhiZDCCA*sinPsiPOI1pPsiPOI2)
                + 2.*cos2PhiZDCCA*(cosPsiPOI1*cosPsiPOI2-sinPsiPOI1*sinPsiPOI2)
                + 2.*sin2PhiZDCCA*(cosPsiPOI1*sinPsiPOI2+sinPsiPOI1*cosPsiPOI2); 
                
  gIntegratedValueV0C = correctedV0C;

  gIntegratedValueV0A = correctedV0A;

  gIntegratedValueZDCC = correctedZDCC;

  gIntegratedValueZDCA = correctedZDCA;
  
  gIntegratedValueZDCCA = correctedZDCCAoverCorrect;
  
  gIntegratedValueZDCCAoverCorrect = correctedZDCCAoverCorrect;
    
  }
  
  if((gSumWeight)&&(iBinCounter)) 
    gIntegratedValue = gSumBinContentTimesWeight/gSumWeight;
  

  f3pCorrelatorPOIIntegratedHist->SetBinContent(1, gIntegratedValue);
  f3pCorrelatorPOIIntegratedHist->SetBinError(1, f3pCorrelatorPOIIntegratedPro->GetBinError(1)); // to be improved later 

  f3pCorrelatorPOIm2V0CIntegratedHist->SetBinContent(1, gIntegratedValueV0C);
  f3pCorrelatorPOIm2V0CIntegratedHist->SetBinError(1, f3pCorrelatorPOIm2V0CIntegratedPro->GetBinError(1)); // to be improved later 
  f3pCorrelatorPOIm2V0AIntegratedHist->SetBinContent(1, gIntegratedValueV0A);
  f3pCorrelatorPOIm2V0AIntegratedHist->SetBinError(1, f3pCorrelatorPOIm2V0AIntegratedPro->GetBinError(1)); // to be improved later 
  f3pCorrelatorPOIm2ZDCCIntegratedHist->SetBinContent(1, gIntegratedValueZDCC);
  f3pCorrelatorPOIm2ZDCCIntegratedHist->SetBinError(1, f3pCorrelatorPOIm2ZDCCIntegratedPro->GetBinError(1)); // to be improved later 
  f3pCorrelatorPOIm2ZDCAIntegratedHist->SetBinContent(1, gIntegratedValueZDCA);
  f3pCorrelatorPOIm2ZDCAIntegratedHist->SetBinError(1, f3pCorrelatorPOIm2ZDCAIntegratedPro->GetBinError(1)); // to be improved later 
  f3pCorrelatorPOIm2ZDCCAIntegratedHist->SetBinContent(1, gIntegratedValueZDCCA);
  f3pCorrelatorPOIm2ZDCCAIntegratedHist->SetBinError(1, f3pCorrelatorPOIm2ZDCCAIntegratedPro->GetBinError(1)); // to be improved later 
  f3pCorrelatorPOIm2ZDCCAIntegratedHist->SetBinContent(2, gIntegratedValueZDCCAoverCorrect);
  f3pCorrelatorPOIm2ZDCCAIntegratedHist->SetBinError(2, f3pCorrelatorPOIm2ZDCCAIntegratedPro->GetBinError(1)); // to be improved later 
  
  // c.) Correct 2-p correlator cos[2(phi1-phi2)] for detector effects. This is needed for DeltaGamma/v2
  Double_t measured2pPsiCorrelator = f2pCorrelatorCosPsiDiffPro->GetBinContent(1); // cos[2(phi1-phi2)] biased by detector effects
  Double_t measured2p2PsiCorrelator = f2pCorrelatorCos2PsiDiffPro->GetBinContent(1); // cos[2(phi1-phi2)] biased by detector effects
  Double_t corrected2pPsiCorrelator = 0.;
  Double_t corrected2p2PsiCorrelator = 0.;
  
  Double_t nonIsotropicTermCosPsiPOI1 = fNonIsotropicPOITermsPro->GetBinContent(1); // cos(#phi_{POI_1})
  Double_t nonIsotropicTermSinPsiPOI1 = fNonIsotropicPOITermsPro->GetBinContent(2); // sin(#phi_{POI_1})
  Double_t nonIsotropicTermCosPsiPOI2 = fNonIsotropicPOITermsPro->GetBinContent(3); // cos(#phi_{POI_2})
  Double_t nonIsotropicTermSinPsiPOI2 = fNonIsotropicPOITermsPro->GetBinContent(4); // sin(#phi_{POI_2})
  Double_t nonIsotropicTermCos2PsiPOI1 = fNonIsotropicPOITermsPro->GetBinContent(5); // cos(2#phi_{POI_1})
  Double_t nonIsotropicTermSin2PsiPOI1 = fNonIsotropicPOITermsPro->GetBinContent(6); // sin(2#phi_{POI_1})
  Double_t nonIsotropicTermCos2PsiPOI2 = fNonIsotropicPOITermsPro->GetBinContent(7); // cos(2#phi_{POI_2})
  Double_t nonIsotropicTermSin2PsiPOI2 = fNonIsotropicPOITermsPro->GetBinContent(8); // sin(2#phi_{POI_2})

 
  // Calculate corrected 2-p correlator:
  corrected2pPsiCorrelator = measured2pPsiCorrelator - (nonIsotropicTermCosPsiPOI1*nonIsotropicTermCosPsiPOI2+nonIsotropicTermSinPsiPOI1*nonIsotropicTermSinPsiPOI2); // <<cos(phi1-phi2)>> - (<<cos(phi1)>><<cos(phi2)>>+<<sin(phi1)>><<sin(phi2)>>)
  f2pCorrelatorCosPsiDiffHist->SetBinContent(1, corrected2pPsiCorrelator);
  f2pCorrelatorCosPsiDiffHist->SetBinError(1, f2pCorrelatorCosPsiDiffPro->GetBinError(1));
  
  corrected2p2PsiCorrelator = measured2p2PsiCorrelator - (nonIsotropicTermCos2PsiPOI1*nonIsotropicTermCos2PsiPOI2+nonIsotropicTermSin2PsiPOI1*nonIsotropicTermSin2PsiPOI2); // <<cos(2(phi1-phi2))>> - (<<cos(2phi1)>><<cos(2phi2)>>+<<sin(2phi1)>><<sin(2phi2)>>)
  f2pCorrelatorCos2PsiDiffHist->SetBinContent(1, corrected2p2PsiCorrelator);
  f2pCorrelatorCos2PsiDiffHist->SetBinError(1, f2pCorrelatorCos2PsiDiffPro->GetBinError(1));
  
  
  
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
 
 fPsi2V0C = 0;
 fPsi2V0A = 0;
 fPsiZNCA = 0;
 fPsiZNCC = 0;
 fPsiZNCCA = 0;

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
 for(Int_t p12=0;p12<2;p12++) // 1st/2nd POI 
 {
  for(Int_t ao=0;ao<2;ao++) // all/overlap 
  {
   for(Int_t pe=0;pe<4;pe++) // [(p1+p2)/2,|p1-p2|,(eta1+eta2)/2,|eta1-eta2|]
   {
    fReNITEBE[p12][ao][pe]->Reset();
    fImNITEBE[p12][ao][pe]->Reset();   
   } // end of for(Int_t pe=0;pe<4;pe++) // [(p1+p2)/2,|p1-p2|,(eta1+eta2)/2,|eta1-eta2|]
  } 
 } // end of for(Int_t p12=0;p12<2;p12++) // 1st/2nd POI 
 
 
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
 // b) Calculate using particle weights;
 
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

void AliFlowAnalysisWithMixedHarmonics::CalculateDifferential3pCorrelator(Double_t &gIntegratedValue, Double_t &gIntegratedValueV0C, Double_t &gIntegratedValueV0A, Double_t &gIntegratedValueZDCC, Double_t &gIntegratedValueZDCA, Double_t &gIntegratedValueZDCCA)
{
 // Calculate differential 3-p azimuthal correlator cos[n(psi1+psi2-2phi3)] in terms of Q_{2n}, p_{n}, q1_{n} and q2_{n}.
 
 // a) Calculate differential 3-p correlator without using particle weights;
 // b) Calculate differential 3-p correlator with using particle weights;
 // c) Calculate non-isotropic terms for 3-p correlator.

 // a) Calculate differential 3-p correlator without using particle weights: 
 if(!(fUsePhiWeights || fUsePtWeights || fUseEtaWeights))
 {
   Int_t iBinCounter = 0;
   Double_t gSumBinContentTimesWeight = 0., gSumWeight = 0.;
   Double_t gSumBinContentTimesWeightSquared = 0.;

   Double_t gSumBinContentTimesWeightV0C = 0., gSumWeightV0C = 0.;
   Double_t gSumBinContentTimesWeightSquaredV0C = 0.;
   
   Double_t gSumBinContentTimesWeightV0A = 0., gSumWeightV0A = 0.;
   Double_t gSumBinContentTimesWeightSquaredV0A = 0.;
   
   Double_t gSumBinContentTimesWeightZDCC = 0., gSumWeightZDCC = 0.;
   Double_t gSumBinContentTimesWeightSquaredZDCC = 0.;
   
   Double_t gSumBinContentTimesWeightZDCA = 0., gSumWeightZDCA = 0.;
   Double_t gSumBinContentTimesWeightSquaredZDCA = 0.;
   
   Double_t gSumBinContentTimesWeightZDCCA = 0., gSumWeightZDCCA = 0.;
   Double_t gSumBinContentTimesWeightSquaredZDCCA = 0.;
   
  // Multiplicity (number of RPs):
  Double_t dMult = (*fSpk)(0,0);
  // Real and imaginary parts of non-weighted Q-vectors (Q_{n,0}) evaluated in harmonic 2n: 
  Double_t dReQ2n = (*fReQnk)(1,0);
  Double_t dImQ2n = (*fImQnk)(1,0);
  
  Double_t dReQ2nV0C = cos(2*fPsi2V0C);
  Double_t dImQ2nV0C = sin(2*fPsi2V0C);
  Double_t dReQ2nV0A = cos(2*fPsi2V0A);
  Double_t dImQ2nV0A = sin(2*fPsi2V0A);
  Double_t dReQ2nZDCC = cos(2*fPsiZNCC);
  Double_t dImQ2nZDCC = sin(2*fPsiZNCC);
  Double_t dReQ2nZDCA = cos(2*fPsiZNCA);
  Double_t dImQ2nZDCA = sin(2*fPsiZNCA);
  Double_t dReQ2nZDCCA = cos(2*fPsiZNCCA);
  Double_t dImQ2nZDCCA = sin(2*fPsiZNCCA);
  
  Double_t Intgp1nRe = 0, Intgp1nIm = 0, Intgmp = 0, Intgoverlap1 = 0, Intgoverlap2 = 0, Intgweight = 0;
  
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
    Double_t weight = mp*dMult-mOverlap1-mOverlap2;  // if POI1, POI2 and RP selection is same except for charge, mp=mOverlap1=mOverlap2
    Double_t weightForV0andZDC = mp; // no overlap between psi1, psi2 and V0,ZDC. dMult for V0 and ZDC is just 1
    
    Double_t cosP2nphi1M1npsi2M1npsi2 = 0; // cos[n(psi1+psi2-2phi3)]
    Double_t cosP2nphi1M1npsi2M1npsiV0C = 0, cosP2nphi1M1npsi2M1npsiV0A = 0; // cos[n(psi1+psi2-2phi_V0)]
    Double_t cosP2nphi1M1npsi2M1npsiZDCC = 0, cosP2nphi1M1npsi2M1npsiZDCA = 0, cosP2nphi1M1npsi2M1npsiZDCCA = 0; // cos[n(psi1+psi2-2phi_ZDC)]
    
    if(weight>0.)
    {
     cosP2nphi1M1npsi2M1npsi2 = (p1nRe*dReQ2n+p1nIm*dImQ2n-overlap1-overlap2)/(weight);
     f3pCorrelatorVsPtSumDiffPro[sd]->Fill(fPtMin+(b-1)*fPtBinWidth,cosP2nphi1M1npsi2M1npsi2,weight);
    }
    
    if(weightForV0andZDC>0.)
    {
     // Calculate differential 3-p correlator for V0
     cosP2nphi1M1npsi2M1npsiV0C = (p1nRe*dReQ2nV0C+p1nIm*dImQ2nV0C)/(weightForV0andZDC);
     cosP2nphi1M1npsi2M1npsiV0A = (p1nRe*dReQ2nV0A+p1nIm*dImQ2nV0A)/(weightForV0andZDC);
     // Calculate differential 3-p correlator for ZDC
     cosP2nphi1M1npsi2M1npsiZDCC = (p1nRe*dReQ2nZDCC+p1nIm*dImQ2nZDCC)/(weightForV0andZDC);
     cosP2nphi1M1npsi2M1npsiZDCA = (p1nRe*dReQ2nZDCA+p1nIm*dImQ2nZDCA)/(weightForV0andZDC);
     cosP2nphi1M1npsi2M1npsiZDCCA = (p1nRe*dReQ2nZDCCA+p1nIm*dImQ2nZDCCA)/(weightForV0andZDC);
     f3pCorrelatorV0CVsPtSumDiffPro[sd]->Fill(fPtMin+(b-1)*fPtBinWidth,cosP2nphi1M1npsi2M1npsiV0C,weightForV0andZDC);
     f3pCorrelatorV0AVsPtSumDiffPro[sd]->Fill(fPtMin+(b-1)*fPtBinWidth,cosP2nphi1M1npsi2M1npsiV0A,weightForV0andZDC);
     f3pCorrelatorZDCCVsPtSumDiffPro[sd]->Fill(fPtMin+(b-1)*fPtBinWidth,cosP2nphi1M1npsi2M1npsiZDCC,weightForV0andZDC);
     f3pCorrelatorZDCAVsPtSumDiffPro[sd]->Fill(fPtMin+(b-1)*fPtBinWidth,cosP2nphi1M1npsi2M1npsiZDCA,weightForV0andZDC);
     f3pCorrelatorZDCCAVsPtSumDiffPro[sd]->Fill(fPtMin+(b-1)*fPtBinWidth,cosP2nphi1M1npsi2M1npsiZDCCA,weightForV0andZDC);
	}
  
    if(sd == 0) {
      iBinCounter += 1;

	  Intgp1nRe += fRePEBE[sd]->GetBinContent(b)*fRePEBE[sd]->GetBinEntries(b);
	  Intgp1nIm += fImPEBE[sd]->GetBinContent(b)*fImPEBE[sd]->GetBinEntries(b);
	  Intgmp += fRePEBE[sd]->GetBinEntries(b);
	  Intgoverlap1 += fOverlapEBE[0][sd]->GetBinContent(b)*fOverlapEBE[0][sd]->GetBinEntries(b);
	  Intgoverlap2 += fOverlapEBE[1][sd]->GetBinContent(b)*fOverlapEBE[1][sd]->GetBinEntries(b);   
	  Intgweight += fRePEBE[sd]->GetBinEntries(b)*dMult-fOverlapEBE[0][sd]->GetBinEntries(b)-fOverlapEBE[1][sd]->GetBinEntries(b);  
	   
	  
      //gSumBinContentTimesWeight += f3pCorrelatorVsPtSumDiffPro[sd]->GetBinContent(b)*f3pCorrelatorVsPtSumDiffPro[sd]->GetBinEntries(b);
      gSumWeight += f3pCorrelatorVsPtSumDiffPro[sd]->GetBinEntries(b);
      //gSumBinContentTimesWeightSquared += TMath::Power(gSumBinContentTimesWeight,2);
    }
       
       
    // non-isotropic terms, 1st POI:
    Double_t p1nRePOI1 = fReNITEBE[0][0][sd]->GetBinContent(b)*fReNITEBE[0][0][sd]->GetBinEntries(b); //Cos(dPsi1) POI1
    Double_t p1nImPOI1 = fImNITEBE[0][0][sd]->GetBinContent(b)*fImNITEBE[0][0][sd]->GetBinEntries(b); //Sin(dPsi1) POI1
    Double_t mpPOI1 = fReNITEBE[0][0][sd]->GetBinEntries(b);
    Double_t q1nRePOI1 = fReNITEBE[0][1][sd]->GetBinContent(b)*fReNITEBE[0][1][sd]->GetBinEntries(b); //Cos(dPsi1) POI1&RP
    Double_t q1nImPOI1 = fImNITEBE[0][1][sd]->GetBinContent(b)*fImNITEBE[0][1][sd]->GetBinEntries(b); //Sin(dPsi1) POI1&RP
    Double_t mqPOI1 = fReNITEBE[0][1][sd]->GetBinEntries(b);   
    // non-isotropic terms, 2nd POI:
    Double_t p1nRePOI2 = fReNITEBE[1][0][sd]->GetBinContent(b)*fReNITEBE[1][0][sd]->GetBinEntries(b); //Cos(dPsi1) POI2
    Double_t p1nImPOI2 = fImNITEBE[1][0][sd]->GetBinContent(b)*fImNITEBE[1][0][sd]->GetBinEntries(b); //Sin(dPsi1) POI2
    Double_t mpPOI2 = fReNITEBE[1][0][sd]->GetBinEntries(b);
    Double_t q1nRePOI2 = fReNITEBE[1][1][sd]->GetBinContent(b)*fReNITEBE[1][1][sd]->GetBinEntries(b); //Cos(dPsi1) POI2&RP
    Double_t q1nImPOI2 = fImNITEBE[1][1][sd]->GetBinContent(b)*fImNITEBE[1][1][sd]->GetBinEntries(b); //Sin(dPsi1) POI2&RP
    Double_t mqPOI2 = fReNITEBE[1][1][sd]->GetBinEntries(b);
    // Fill all-event profiles:   
    if(weight>0. && mpPOI1>0.)
    { 
     fNonIsotropicTermsVsPtSumDiffPro[sd][0]->Fill(fPtMin+(b-1)*fPtBinWidth,p1nRePOI1/mpPOI1,mpPOI1); // <<cos(#psi_{POI_1})>>
     fNonIsotropicTermsVsPtSumDiffPro[sd][1]->Fill(fPtMin+(b-1)*fPtBinWidth,p1nImPOI1/mpPOI1,mpPOI1); // <<sin(#psi_{POI_1})>>
    }
    if(weight>0. && mpPOI2>0.)    
    {
     fNonIsotropicTermsVsPtSumDiffPro[sd][2]->Fill(fPtMin+(b-1)*fPtBinWidth,p1nRePOI2/mpPOI2,mpPOI2); // <<cos(#psi_{POI_2})>>
     fNonIsotropicTermsVsPtSumDiffPro[sd][3]->Fill(fPtMin+(b-1)*fPtBinWidth,p1nImPOI2/mpPOI2,mpPOI2); // <<sin(#psi_{POI_2})>>
    }
    if(weight>0. && mpPOI1*dMult-mqPOI1>0.)
    { 
     fNonIsotropicTermsVsPtSumDiffPro[sd][4]->Fill(fPtMin+(b-1)*fPtBinWidth,
      (p1nRePOI1*dReQ2n+p1nImPOI1*dImQ2n-q1nRePOI1)/(mpPOI1*dMult-mqPOI1),mpPOI1*dMult-mqPOI1); // <<cos(#psi_{POI_1}-2*phi)>>
     fNonIsotropicTermsVsPtSumDiffPro[sd][5]->Fill(fPtMin+(b-1)*fPtBinWidth,
      (p1nImPOI1*dReQ2n-p1nRePOI1*dImQ2n+q1nImPOI1)/(mpPOI1*dMult-mqPOI1),mpPOI1*dMult-mqPOI1); // <<sin(#psi_{POI_1}-2*phi)>>
    }
    if(weight>0. && mpPOI2*dMult-mqPOI2>0.)
    { 
     fNonIsotropicTermsVsPtSumDiffPro[sd][6]->Fill(fPtMin+(b-1)*fPtBinWidth,
      (p1nRePOI2*dReQ2n+p1nImPOI2*dImQ2n-q1nRePOI2)/(mpPOI2*dMult-mqPOI2),mpPOI2*dMult-mqPOI2); // <<cos(#psi_{POI_2}-2*phi)>>
     fNonIsotropicTermsVsPtSumDiffPro[sd][7]->Fill(fPtMin+(b-1)*fPtBinWidth,
      (p1nImPOI2*dReQ2n-p1nRePOI2*dImQ2n+q1nImPOI2)/(mpPOI2*dMult-mqPOI2),mpPOI2*dMult-mqPOI2); // <<sin(#psi_{POI_2}-2*phi)>>
    }
    if(weight>0. && mp>0.)
    {
     fNonIsotropicTermsVsPtSumDiffPro[sd][8]->Fill(fPtMin+(b-1)*fPtBinWidth,p1nRe/mp,mp); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
     fNonIsotropicTermsVsPtSumDiffPro[sd][9]->Fill(fPtMin+(b-1)*fPtBinWidth,p1nIm/mp,mp); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>   
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
   
    // non-isotropic terms, 1st POI:
    Double_t p1nRePOI1 = fReNITEBE[0][0][sd+2]->GetBinContent(k)*fReNITEBE[0][0][sd+2]->GetBinEntries(k);
    Double_t p1nImPOI1 = fImNITEBE[0][0][sd+2]->GetBinContent(k)*fImNITEBE[0][0][sd+2]->GetBinEntries(k);
    Double_t mpPOI1 = fReNITEBE[0][0][sd+2]->GetBinEntries(k);
    Double_t q1nRePOI1 = fReNITEBE[0][1][sd+2]->GetBinContent(k)*fReNITEBE[0][1][sd+2]->GetBinEntries(k);
    Double_t q1nImPOI1 = fImNITEBE[0][1][sd+2]->GetBinContent(k)*fImNITEBE[0][1][sd+2]->GetBinEntries(k);
    Double_t mqPOI1 = fReNITEBE[0][1][sd+2]->GetBinEntries(k);   
    // non-isotropic terms, 2nd POI:
    Double_t p1nRePOI2 = fReNITEBE[1][0][sd+2]->GetBinContent(k)*fReNITEBE[1][0][sd+2]->GetBinEntries(k);
    Double_t p1nImPOI2 = fImNITEBE[1][0][sd+2]->GetBinContent(k)*fImNITEBE[1][0][sd+2]->GetBinEntries(k);
    Double_t mpPOI2 = fReNITEBE[1][0][sd+2]->GetBinEntries(k);
    Double_t q1nRePOI2 = fReNITEBE[1][1][sd+2]->GetBinContent(k)*fReNITEBE[1][1][sd+2]->GetBinEntries(k);
    Double_t q1nImPOI2 = fImNITEBE[1][1][sd+2]->GetBinContent(k)*fImNITEBE[1][1][sd+2]->GetBinEntries(k);
    Double_t mqPOI2 = fReNITEBE[1][1][sd+2]->GetBinEntries(k);
    // Fill all-event profiles:    
    if(weight>0. && mpPOI1>0.)
    {
     fNonIsotropicTermsVsEtaSumDiffPro[sd][0]->Fill(fEtaMin+(k-1)*fEtaBinWidth,p1nRePOI1/mpPOI1,mpPOI1); // <<cos(#psi_{POI_1})>>
     fNonIsotropicTermsVsEtaSumDiffPro[sd][1]->Fill(fEtaMin+(k-1)*fEtaBinWidth,p1nImPOI1/mpPOI1,mpPOI1); // <<sin(#psi_{POI_1})>>
    }
    if(weight>0. && mpPOI2>0.)
    {     
     fNonIsotropicTermsVsEtaSumDiffPro[sd][2]->Fill(fEtaMin+(k-1)*fEtaBinWidth,p1nRePOI2/mpPOI2,mpPOI2); // <<cos(#psi_{POI_2})>>
     fNonIsotropicTermsVsEtaSumDiffPro[sd][3]->Fill(fEtaMin+(k-1)*fEtaBinWidth,p1nImPOI2/mpPOI2,mpPOI2); // <<sin(#psi_{POI_2})>>
    }
    if(weight>0. && mpPOI1*dMult-mqPOI1>0.)
    {
     fNonIsotropicTermsVsEtaSumDiffPro[sd][4]->Fill(fEtaMin+(k-1)*fEtaBinWidth,
      (p1nRePOI1*dReQ2n+p1nImPOI1*dImQ2n-q1nRePOI1)/(mpPOI1*dMult-mqPOI1),mpPOI1*dMult-mqPOI1); // <<cos(#psi_{POI_1}-2*phi)>>
     fNonIsotropicTermsVsEtaSumDiffPro[sd][5]->Fill(fEtaMin+(k-1)*fEtaBinWidth,
      (p1nImPOI1*dReQ2n-p1nRePOI1*dImQ2n+q1nImPOI1)/(mpPOI1*dMult-mqPOI1),mpPOI1*dMult-mqPOI1); // <<sin(#psi_{POI_1}-2*phi)>>
    }
    if(weight>0. && mpPOI2*dMult-mqPOI2>0.)
    {
     fNonIsotropicTermsVsEtaSumDiffPro[sd][6]->Fill(fEtaMin+(k-1)*fEtaBinWidth,
      (p1nRePOI2*dReQ2n+p1nImPOI2*dImQ2n-q1nRePOI2)/(mpPOI2*dMult-mqPOI2),mpPOI2*dMult-mqPOI2); // <<cos(#psi_{POI_2}-2*phi)>>
     fNonIsotropicTermsVsEtaSumDiffPro[sd][7]->Fill(fEtaMin+(k-1)*fEtaBinWidth,
      (p1nImPOI2*dReQ2n-p1nRePOI2*dImQ2n+q1nImPOI2)/(mpPOI2*dMult-mqPOI2),mpPOI2*dMult-mqPOI2); // <<sin(#psi_{POI_2}-2*phi)>>
    }
    if(weight>0. && mp>0.)
    {
     fNonIsotropicTermsVsEtaSumDiffPro[sd][8]->Fill(fEtaMin+(k-1)*fEtaBinWidth,p1nRe/mp,mp); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
     fNonIsotropicTermsVsEtaSumDiffPro[sd][9]->Fill(fEtaMin+(k-1)*fEtaBinWidth,p1nIm/mp,mp); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>    
    }
   } // end of for(Int_t k=1;k<=fnBinsEta;k++)  
  } // end of for(Int_t sd=0;sd<2;sd++)      

  gIntegratedValue = -1000.;
  gIntegratedValueV0C = -1000.;
  gIntegratedValueV0A = -1000.;
  gIntegratedValueZDCC = -1000.;
  gIntegratedValueZDCA = -1000.;
  gIntegratedValueZDCCA = -1000.;
  
  //if((gSumWeight)&&(iBinCounter)) 
    //gIntegratedValue = gSumBinContentTimesWeight/gSumWeight;
    
  if((gSumWeight)&&(iBinCounter)) 
    gIntegratedValue = (Intgp1nRe*dReQ2n+Intgp1nIm*dImQ2n-Intgoverlap1-Intgoverlap2)/(Intgweight);
  if((gSumWeight)&&(iBinCounter)) 
    gIntegratedValueV0C = (Intgp1nRe*dReQ2nV0C+Intgp1nIm*dImQ2nV0C)/(Intgmp);
  if((gSumWeight)&&(iBinCounter)) 
    gIntegratedValueV0A = (Intgp1nRe*dReQ2nV0A+Intgp1nIm*dImQ2nV0A)/(Intgmp);
  if((gSumWeight)&&(iBinCounter)) 
    gIntegratedValueZDCC = (Intgp1nRe*dReQ2nZDCC+Intgp1nIm*dImQ2nZDCC)/(Intgmp);
  if((gSumWeight)&&(iBinCounter)) 
    gIntegratedValueZDCA = (Intgp1nRe*dReQ2nZDCA+Intgp1nIm*dImQ2nZDCA)/(Intgmp);
  if((gSumWeight)&&(iBinCounter)) 
    gIntegratedValueZDCCA = (Intgp1nRe*dReQ2nZDCCA+Intgp1nIm*dImQ2nZDCCA)/(Intgmp);


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
