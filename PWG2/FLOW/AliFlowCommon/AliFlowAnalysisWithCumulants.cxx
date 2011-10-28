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

/************************************************* 
 * Flow analysis with cumulants. In this class   *
 * cumulants are calculated by making use of the *
 * formalism of generating functions proposed by *
 * Ollitrault et al.                             *
 *                                               * 
 *      Author: Ante Bilandzic                   * 
 *              (abilandzic@gmail.com)           *
 *************************************************/ 

#define AliFlowAnalysisWithCumulants_cxx

#include "Riostream.h"
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TProfile.h"
#include "TProfile2D.h" 
#include "TProfile3D.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowVector.h"

//================================================================================================================

ClassImp(AliFlowAnalysisWithCumulants)

AliFlowAnalysisWithCumulants::AliFlowAnalysisWithCumulants(): 
 fHistList(NULL),
 fHistListName(NULL),
 fAnalysisSettings(NULL),
 fCommonHists(NULL),
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
 fMultiple(1),
 fR0(2.2),
 fWeightsList(NULL),
 fUsePhiWeights(kFALSE),
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fPhiWeights(NULL),
 fPtWeights(NULL),
 fEtaWeights(NULL),
 fMultiplicityWeight(NULL),
 fReferenceFlowList(NULL),
 fReferenceFlowProfiles(NULL),
 fReferenceFlowResults(NULL),
 fReferenceFlowFlags(NULL),
 fCalculateVsMultiplicity(kFALSE),
 fnBinsMult(10000),  
 fMinMult(0.),   
 fMaxMult(10000.),
 fGEBE(NULL), 
 fReferenceFlowGenFun(NULL),
 fQvectorComponents(NULL),
 fAverageOfSquaredWeight(NULL),
 fReferenceFlowGenFunVsM(NULL),
 fQvectorComponentsVsM(NULL),
 fAverageOfSquaredWeightVsM(NULL),
 fAvMVsM(NULL),
 fAvM(0.),
 fnEvts(0), 
 fReferenceFlowCumulants(NULL), 
 fReferenceFlow(NULL),
 fChi(NULL),
 fDiffFlowList(NULL),
 fDiffFlowProfiles(NULL),
 fDiffFlowResults(NULL),
 fDiffFlowFlags(NULL),
 fTuningList(NULL), 
 fTuningProfiles(NULL),
 fTuningResults(NULL),
 fTuningFlags(NULL),
 fTuneParameters(kFALSE),
 fTuningAvM(NULL)   
 {
  // Constructor. 
  
  // Base list to hold all output objects:
  fHistList = new TList();
  fHistListName = new TString("cobjGFC");
  fHistList->SetName(fHistListName->Data());
  fHistList->SetOwner(kTRUE);
 
  // Multiplicity weight:
  fMultiplicityWeight = new TString("unit");
   
  // Initialize all arrays:
  this->InitializeArrays();
 
} // end of AliFlowAnalysisWithCumulants::AliFlowAnalysisWithCumulants()

//================================================================================================================

AliFlowAnalysisWithCumulants::~AliFlowAnalysisWithCumulants()
{
 // Desctructor.
 
 delete fHistList; 

} // end of AliFlowAnalysisWithCumulants::~AliFlowAnalysisWithCumulants()

//================================================================================================================

void AliFlowAnalysisWithCumulants::Init()
{
 // Initialize and book all objects. 
 
 // a) Cross check if the user settings make sense before starting; 
 // b) Access all common constants;
 // c) Book and fill weights histograms;
 // d) Book and nest all lists in the base list fHistList;
 // e) Book and fill profile holding analysis settings;
 // f) Book common control and results histograms;
 // g) Store flags for reference flow;
 // h) Store flags for differential flow;
 // i) Book all objects needed for tuning;
 // j) Book all objects needed for calculation versus multiplicity.
 
 //save old value and prevent histograms from being added to directory
 //to avoid name clashes in case multiple analaysis objects are used
 //in an analysis
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
 TH1::AddDirectory(kFALSE);
 
 this->CrossCheckSettings();
 this->AccessConstants();
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights){this->BookAndFillWeightsHistograms();}
 this->BookAndNestAllLists();
 this->BookProfileHoldingSettings();
 this->BookCommonHistograms();
 this->BookEverythingForReferenceFlow();
 this->BookEverythingForDiffFlow();
 this->StoreReferenceFlowFlags();
 this->StoreDiffFlowFlags();
 if(fTuneParameters){this->BookEverythingForTuning();}
 if(fCalculateVsMultiplicity){this->BookEverythingForCalculationVsMultiplicity();}
 
 (fCommonHists->GetHarmonic())->Fill(0.5,fHarmonic); // to be improved (moved somewhere else?)
 
 TH1::AddDirectory(oldHistAddStatus);

} // end of void AliFlowAnalysisWithCumulants::Init()

//================================================================================================================

void AliFlowAnalysisWithCumulants::Make(AliFlowEventSimple* anEvent)
{
 // Running over data only in this method.
 
 // a) Check all pointers used in this method;
 // b) If tuning enabled, fill generating functions for different values of tuning parameters;
 // c) For default values of tuning parameters (r0 = 2.2 and cutoff at 10th order):
 //  c1) Fill common control histograms;
 //  c2) Fill generating function for reference flow;
 //  c3) Fill profile holding average values of various Q-vector components;   
 //  c4) Fill generating function for differential flow.
 
 this->CheckPointersUsedInMake();
 if(fTuneParameters) {this->FillGeneratingFunctionsForDifferentTuningParameters(anEvent);} 
 Int_t nRP = anEvent->GetEventNSelTracksRP(); // number of RPs (i.e. number of particles used to determine the reaction plane)
 if(nRP<10) {return;} // generating function formalism make sense only for nRPs >= 10 for default settings 
 fCommonHists->FillControlHistograms(anEvent);                                                               
 this->FillGeneratingFunctionForReferenceFlow(anEvent);
 this->FillQvectorComponents(anEvent);
 this->FillGeneratingFunctionForDiffFlow(anEvent);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
} // end of void AliFlowAnalysisWithCumulants::Make()

//================================================================================================================

void AliFlowAnalysisWithCumulants::Finish()
{
 // Calculate the final results.
 
 // a) Check all pointers used in this method;
 // b) Access all common constants;
 // c) Access settings for analysis with Generating Function Cumulants;
 // d) From relevant common control histogram get average multiplicity of RPs and number of events;
 // e) Calculate cumulants for reference flow;
 // f) Calculate from isotropic cumulants reference flow;
 // g) Calculate error for reference flow estimates;
 // h) Store the final results for reference flow in common hist results;
 // i) Print on the screen the final results for reference flow;
 // j) Calculate cumulants for differential flow;
 // k) Calculate differential flow for RPs/POIs vs pt/eta from cumulants;
 // l) Calculate integrated flow of RPs and POIs;
 // m) Print on the screen the final results for integrated flow of RPs and POIs;
 // n) If tuning enabled, calculate results for different tuning parameters.
 
 this->CheckPointersUsedInFinish();
 this->AccessConstants(); 
 this->AccessSettings();
 this->GetAvMultAndNoOfEvts();
 this->CalculateCumulantsForReferenceFlow(); 
 this->CalculateReferenceFlow();
 this->CalculateReferenceFlowError();
 this->FillCommonHistResultsForReferenceFlow();
 if(fPrintFinalResults[0]){this->PrintFinalResults("RF");}
 this->CalculateCumulantsForDiffFlow("RP","pt");
 this->CalculateCumulantsForDiffFlow("RP","eta");
 this->CalculateCumulantsForDiffFlow("POI","pt");
 this->CalculateCumulantsForDiffFlow("POI","eta");
 this->CalculateDifferentialFlow("RP","pt");
 this->CalculateDifferentialFlow("RP","eta");
 this->CalculateDifferentialFlow("POI","pt");
 this->CalculateDifferentialFlow("POI","eta");
 this->CalculateDifferentialFlowErrors("RP","pt");
 this->CalculateDifferentialFlowErrors("RP","eta");
 this->CalculateDifferentialFlowErrors("POI","pt");
 this->CalculateDifferentialFlowErrors("POI","eta");
 this->FillCommonHistResultsForDifferentialFlow("RP");
 this->FillCommonHistResultsForDifferentialFlow("POI");
 this->CalculateIntegratedFlow("RP");
 this->CalculateIntegratedFlow("POI");
 if(fPrintFinalResults[1]){this->PrintFinalResults("RP");}
 if(fPrintFinalResults[2]){this->PrintFinalResults("POI");}
 if(fTuneParameters){this->FinalizeTuning();}
     
} // end of void AliFlowAnalysisWithCumulants::Finish()

//================================================================================================================

void AliFlowAnalysisWithCumulants::FinalizeTuning()
{
 // Finalize results with tuned inerpolating parameters.

 for(Int_t r=0;r<10;r++)
 {
  if(TMath::Abs(fTuningR0[r])<1.e-10) continue; // protection against division by r0 bellow
  for(Int_t pq=0;pq<5;pq++)
  {
   Int_t pMax = fTuningGenFun[r][pq]->GetXaxis()->GetNbins();
   Int_t qMax = fTuningGenFun[r][pq]->GetYaxis()->GetNbins(); 
   fAvM = fTuningAvM->GetBinContent(pq+1);
   // <G[p][q]>
   TMatrixD dAvG(pMax,qMax); 
   dAvG.Zero();
   Bool_t someAvGEntryIsNegative = kFALSE;
   for(Int_t p=0;p<pMax;p++)
   {
    for(Int_t q=0;q<qMax;q++)
    {
     dAvG(p,q) = fTuningGenFun[r][pq]->GetBinContent(fTuningGenFun[r][pq]->GetBin(p+1,q+1));
     if(dAvG(p,q)<0.)
     {
      someAvGEntryIsNegative = kTRUE;
      cout<<endl; 
      cout<<" WARNING: "<<Form("<G[%d][%d]> is negative !!!! GFC results are meaningless for r0 = %f, pq = %i.",p,q,fTuningR0[r],pq)<<endl; 
      cout<<endl; 
     }
    }  
   } 
   // C[p][q] (generating function for the cumulants)    
   TMatrixD dC(pMax,qMax);
   dC.Zero();
   if(fAvM>0. && !someAvGEntryIsNegative)
   {
    for(Int_t p=0;p<pMax;p++)
    {
     for(Int_t q=0;q<qMax;q++)
     {
      dC(p,q) = fAvM*(pow(dAvG(p,q),(1./fAvM))-1.); 
     }
    }
   }
   // Averaging the generating function for cumulants over azimuth
   // in order to eliminate detector effects.
   // <C[p][q]> (Remark: here <> stands for average over azimuth):
   TVectorD dAvC(pMax);
   dAvC.Zero();
   for(Int_t p=0;p<pMax;p++)
   {
    Double_t temp = 0.; 
    for(Int_t q=0;q<qMax;q++)
    {
     temp += 1.*dC(p,q);
    } 
    dAvC[p] = temp/qMax;
   }  
   // Finally, the isotropic cumulants for reference flow and reference flow itself:
   TVectorD cumulant(pMax);
   cumulant.Zero(); 
   TVectorD flow(pMax);
   flow.Zero(); 
   if(pMax==2)
   {
    cumulant[0]=(1./(fTuningR0[r]*fTuningR0[r]))*(2.*dAvC[0]-(1./2.)*dAvC[1]);
    cumulant[1]=(2./pow(fTuningR0[r],4.))*((-2.)*dAvC[0]+1.*dAvC[1]);
    if(cumulant[0]>=0.) {flow[0] = pow(cumulant[0],1./2.);}
    if(cumulant[1]<=0.) {flow[1] = pow(-1.*cumulant[1],1./4.);}
   }
   else if(pMax==3)
   {
    cumulant[0] = (1./(fTuningR0[r]*fTuningR0[r]))*(3.*dAvC[0]-(3./2.)*dAvC[1]+(1./3.)*dAvC[2]);
    cumulant[1] = (2./pow(fTuningR0[r],4.))*((-5.)*dAvC[0]+4.*dAvC[1]-1.*dAvC[2]);
    cumulant[2] = (6./pow(fTuningR0[r],6.))*(3.*dAvC[0]-3.*dAvC[1]+1.*dAvC[2]);
    if(cumulant[0]>=0.) {flow[0] = pow(cumulant[0],1./2.);}
    if(cumulant[1]<=0.) {flow[1] = pow(-1.*cumulant[1],1./4.);}
    if(cumulant[2]>=0.) {flow[2] = pow((1./4.)*cumulant[2],1./6.);}
   }
   else if(pMax==4)
   {
    cumulant[0] = (1./(fTuningR0[r]*fTuningR0[r]))*(4.*dAvC[0]-3.*dAvC[1]+(4./3.)*dAvC[2]-(1./4.)*dAvC[3]);
    cumulant[1] = (1./pow(fTuningR0[r],4.))*((-52./3.)*dAvC[0]+19.*dAvC[1]-(28./3.)*dAvC[2]+(11./6.)*dAvC[3]);
    cumulant[2] = (3./pow(fTuningR0[r],6.))*(18.*dAvC[0]-24.*dAvC[1]+14.*dAvC[2]-3.*dAvC[3]);
    cumulant[3] = (24./pow(fTuningR0[r],8.))*((-4.)*dAvC[0]+6.*dAvC[1]-4.*dAvC[2]+1.*dAvC[3]);
    if(cumulant[0]>=0.) {flow[0] = pow(cumulant[0],1./2.);}
    if(cumulant[1]<=0.) {flow[1] = pow(-1.*cumulant[1],1./4.);}
    if(cumulant[2]>=0.) {flow[2] = pow((1./4.)*cumulant[2],1./6.);}
    if(cumulant[3]<=0.) {flow[3] = pow((-1./33.)*cumulant[3],1./8.);}
   }
   else if(pMax==5)
   {
    cumulant[0] = (-1./(60*fTuningR0[r]*fTuningR0[r]))*((-300.)*dAvC[0]+300.*dAvC[1]-200.*dAvC[2]+75.*dAvC[3]-12.*dAvC[4]);
    cumulant[1] = (-1./(6.*pow(fTuningR0[r],4.)))*(154.*dAvC[0]-214.*dAvC[1]+156.*dAvC[2]-61.*dAvC[3]+10.*dAvC[4]);
    cumulant[2] = (3./(2.*pow(fTuningR0[r],6.)))*(71.*dAvC[0]-118.*dAvC[1]+98.*dAvC[2]-41.*dAvC[3]+7.*dAvC[4]);
    cumulant[3] = (-24./pow(fTuningR0[r],8.))*(14.*dAvC[0]-26.*dAvC[1]+24.*dAvC[2]-11.*dAvC[3]+2.*dAvC[4]);
    cumulant[4] = (120./pow(fTuningR0[r],10.))*(5.*dAvC[0]-10.*dAvC[1]+10.*dAvC[2]-5.*dAvC[3]+1.*dAvC[4]); 
    if(cumulant[0]>=0.) {flow[0] = pow(cumulant[0],1./2.);}
    if(cumulant[1]<=0.) {flow[1] = pow(-1.*cumulant[1],1./4.);}
    if(cumulant[2]>=0.) {flow[2] = pow((1./4.)*cumulant[2],1./6.);}
    if(cumulant[3]<=0.) {flow[3] = pow((-1./33.)*cumulant[3],1./8.);}
    if(cumulant[4]>=0.) {flow[4] = pow((1./456.)*cumulant[4],1./10.);}
   }   
   else if(pMax==8)
   {  
    cumulant[0] = (1./(fTuningR0[r]*fTuningR0[r]))*(8.*dAvC[0]-14.*dAvC[1]+(56./3.)*dAvC[2]-(35./2.)*dAvC[3] 
                + (56./5.)*dAvC[4]-(14./3.)*dAvC[5]+(8./7.)*dAvC[6]-(1./8.)*dAvC[7]);
    cumulant[1] = (1./pow(fTuningR0[r],4.))*((-1924./35.)*dAvC[0]+(621./5.)*dAvC[1]-(8012./45.)*dAvC[2] 
                + (691./4.)*dAvC[3]-(564./5.)*dAvC[4]+(2143./45.)*dAvC[5]-(412./35.)*dAvC[6]+(363./280.)*dAvC[7]);
    cumulant[2] = (1./pow(fTuningR0[r],6.))*(349.*dAvC[0]-(18353./20.)*dAvC[1]+(7173./5.)*dAvC[2]
                - 1457.*dAvC[3]+(4891./5.)*dAvC[4]-(1683./4.)*dAvC[5]+(527./5.)*dAvC[6]-(469./40.)*dAvC[7]);
    cumulant[3] = (1./pow(fTuningR0[r],8.))*((-10528./5.)*dAvC[0]+(30578./5.)*dAvC[1]-(51456./5.)*dAvC[2]
                + 10993.*dAvC[3]-(38176./5.)*dAvC[4]+(16818./5.)*dAvC[5]-(4288./5.)*dAvC[6]+(967./10.)*dAvC[7]);
    cumulant[4] = (1./pow(fTuningR0[r],10.))*(11500.*dAvC[0]-35800.*dAvC[1]+63900.*dAvC[2]-71600.*dAvC[3] 
	    + 51620.*dAvC[4]-23400.*dAvC[5]+6100.*dAvC[6]-700.*dAvC[7]);
    cumulant[5] = (1./pow(fTuningR0[r],12.))*(-52560.*dAvC[0]+172080.*dAvC[1]-321840.*dAvC[2]+376200.*dAvC[3]  
                - 281520.*dAvC[4]+131760.*dAvC[5]-35280.*dAvC[6]+4140.*dAvC[7]);
    cumulant[6] = (1./pow(fTuningR0[r],14.))*(176400.*dAvC[0]-599760.*dAvC[1]+1164240.*dAvC[2]-1411200.*dAvC[3] 
	    + 1093680.*dAvC[4]-529200.*dAvC[5]+146160.*dAvC[6]-17640.*dAvC[7]);
    cumulant[7] = (1./pow(fTuningR0[r],16.))*(-322560*dAvC[0]+1128960.*dAvC[1]-2257920.*dAvC[2]+2822400.*dAvC[3] 
         	    - 2257920.*dAvC[4]+1128960.*dAvC[5]-322560.*dAvC[6]+40320.*dAvC[7]);
    if(cumulant[0]>=0.) {flow[0] = pow(cumulant[0],1./2.);}
    if(cumulant[1]<=0.) {flow[1] = pow(-1.*cumulant[1],1./4.);}
    if(cumulant[2]>=0.) {flow[2] = pow((1./4.)*cumulant[2],1./6.);}
    if(cumulant[3]<=0.) {flow[3] = pow((-1./33.)*cumulant[3],1./8.);}
    if(cumulant[4]>=0.) {flow[4] = pow((1./456.)*cumulant[4],1./10.);}
    if(cumulant[5]<=0.) {flow[5] = pow((-1./9460.)*cumulant[5],1./12.);}
    if(cumulant[6]>=0.) {flow[6] = pow((1./274800.)*cumulant[6],1./14.);}
    if(cumulant[7]<=0.) {flow[7] = pow((-1./10643745.)*cumulant[7],1./16.);}
   }
   // Store cumulants and reference flow:
   for(Int_t co=0;co<pMax;co++) // cumulant order
   {
    fTuningCumulants[r][pq]->SetBinContent(co+1,cumulant[co]);
    fTuningFlow[r][pq]->SetBinContent(co+1,flow[co]);
   }
  } // end of for(Int_t pq=0;pq<5;pq++) 
 } // end of for(Int_t r=0;r<10;r++)

} // end of void AliFlowAnalysisWithCumulants::FinalizeTuning()

//================================================================================================================

void AliFlowAnalysisWithCumulants::FillGeneratingFunctionsForDifferentTuningParameters(AliFlowEventSimple *anEvent)
{
 // Fill generating function for reference flow evaluated for different tuning parameters.
 
 Int_t pMax[5] = {2,3,4,5,8};      
 Int_t qMax[5] = {5,7,9,11,17};  
 
 // Particle variables and weights:
 Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
 Double_t dPt  = 0.; // transverse momentum
 Double_t dEta = 0.; // pseudorapidity
 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight

 Int_t nPrim = anEvent->NumberOfTracks(); // nPrim = total number of primary tracks, i.e. nPrim = nRP + nPOI, where:
                                          // nRP   = # of particles used to determine the reaction plane;
                                          // nPOI  = # of particles of interest for a detailed flow analysis.
 
 Int_t nRP = anEvent->GetEventNSelTracksRP(); // nRP = # of particles used to determine the reaction plane;
 for(Int_t pq=0;pq<5;pq++)
 { 
  if(nRP<2.*pMax[pq]) continue; // results doesn't make sense if nRP is smaller than serie's cutoff
  fTuningAvM->Fill(pq+0.5,nRP,1.); // <M> for different classes of events }
 }
  
 Double_t tuningGenFunEBE[10][5][8][17] = {{{{0.}}}}; 
 for(Int_t r=0;r<10;r++)
 {
  for(Int_t pq=0;pq<5;pq++)
  {
   for(Int_t p=0;p<pMax[pq];p++)
   {
    for(Int_t q=0;q<qMax[pq];q++)
    {
     tuningGenFunEBE[r][pq][p][q] = 1.;
    }
   }
  } 
 }
 
 // Looping over tracks:
 for(Int_t i=0;i<nPrim;i++)
 {
  AliFlowTrackSimple *aftsTrack = anEvent->GetTrack(i);
  if(aftsTrack && aftsTrack->InRPSelection())
  {
   // Access particle variables and weights:
   dPhi = aftsTrack->Phi();
   dPt  = aftsTrack->Pt();
   dEta = aftsTrack->Eta();
   if(fUsePhiWeights && fnBinsPhi) // determine phi weight for this particle:
   {
    wPhi = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*fnBinsPhi/TMath::TwoPi())));
   }
   if(fUsePtWeights && fnBinsPt) // determine pt weight for this particle:
   {
    wPt = fPtWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/fPtBinWidth))); 
   }              
   if(fUseEtaWeights && fEtaBinWidth) // determine eta weight for this particle: 
   {
    wEta = fEtaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/fEtaBinWidth))); 
   }    
   // Fill the generating functions:
   for(Int_t r=0;r<10;r++) // 10 different values for interpolating parameter r0
   {
    if(TMath::Abs(fTuningR0[r])<1.e-10) continue;
    for(Int_t pq=0;pq<5;pq++) // 5 different values for set (pMax,qMax)
    {
     if(nRP<2.*pMax[pq]) continue; // results doesn't make sense if nRP is smaller than serie's cutoff
     for(Int_t p=0;p<pMax[pq];p++)
     {
      for(Int_t q=0;q<qMax[pq];q++)
      {
       tuningGenFunEBE[r][pq][p][q] *= (1.+wPhi*wPt*wEta*(2.*fTuningR0[r]*sqrt(p+1.)/nRP)*cos(fHarmonic*dPhi-2.*q*TMath::Pi()/qMax[pq]));
      } // end of for(Int_t q=0;q<qMax[pq];q++) 
     } // end of for(Int_t p=0;p<pMax[pq];p++)
    } // end for(Int_t pq=0;pq<5;pq++) // 5 different values for set (pMax,qMax)
   } // end of for(Int_t r=0;r<10;r++) // 10 different values for interpolating parameter r0  
  } // end of if(aftsTrack && aftsTrack->InRPSelection())
 } // end of for(Int_t i=0;i<nPrim;i++) 
  
 // Store G[p][q]:
 for(Int_t r=0;r<10;r++) 
 {
  for(Int_t pq=0;pq<5;pq++) 
  {
   if(nRP<2.*pMax[pq]) continue; // results doesn't make sense if nRP is smaller than serie's cutoff
   for(Int_t p=0;p<pMax[pq];p++)
   {
    for(Int_t q=0;q<qMax[pq];q++)
    {
     if(fTuningGenFun[r][pq]) {fTuningGenFun[r][pq]->Fill((Double_t)p,(Double_t)q,tuningGenFunEBE[r][pq][p][q],1.);} 
    }
   } 
  }
 } 
  
} // end of void AliFlowAnalysisWithCumulants::FillGeneratingFunctionsForDifferentTuningParameters(AliFlowEventSimple *anEvent)

//================================================================================================================

void AliFlowAnalysisWithCumulants::FillGeneratingFunctionForReferenceFlow(AliFlowEventSimple *anEvent)
{
 // Fill generating function for reference flow for current event.
 
 if(!anEvent)
 {
  printf(" WARNING (GFC): anEvent is NULL !!!!");
  return;
 }
 
 // Particle variables and weights:
 Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
 Double_t dPt  = 0.; // transverse momentum
 Double_t dEta = 0.; // pseudorapidity
 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight
  
 Int_t nPrim = anEvent->NumberOfTracks(); // nPrim = total number of primary tracks, i.e. nPrim = nRP + nPOI, where:
                                          // nRP   = # of particles used to determine the reaction plane;
                                          // nPOI  = # of particles of interest for a detailed flow analysis.
 
 Int_t nRP = anEvent->GetEventNSelTracksRP(); // nRP = # of particles used to determine the reaction plane;
 if(fCalculateVsMultiplicity){fAvMVsM->Fill(nRP+0.5,nRP,1.);}
 
 // Initializing the generating function G[p][q] for reference flow for current event: 
 Int_t pMax = fGEBE->GetNrows();
 Int_t qMax = fGEBE->GetNcols();
 for(Int_t p=0;p<pMax;p++)
 {
  for(Int_t q=0;q<qMax;q++)
  {
   (*fGEBE)(p,q) = 1.;
  }   
 }
    
 // Cross-checking the number of RPs in current event:
 Int_t crossCheckRP = 0; 
 
 // Looping over tracks:
 for(Int_t i=0;i<nPrim;i++)
 {
  AliFlowTrackSimple *aftsTrack = anEvent->GetTrack(i);
  if(aftsTrack && aftsTrack->InRPSelection())
  {
   crossCheckRP++;
   // Access particle variables and weights:
   dPhi = aftsTrack->Phi();
   dPt  = aftsTrack->Pt();
   dEta = aftsTrack->Eta();
   if(fUsePhiWeights && fnBinsPhi) // determine phi weight for this particle:
   {
    wPhi = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*fnBinsPhi/TMath::TwoPi())));
   }
   if(fUsePtWeights && fnBinsPt) // determine pt weight for this particle:
   {
    wPt = fPtWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/fPtBinWidth))); 
   }              
   if(fUseEtaWeights && fEtaBinWidth) // determine eta weight for this particle: 
   {
    wEta = fEtaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/fEtaBinWidth))); 
   }    
   // Fill the generating function:
   for(Int_t p=0;p<pMax;p++)
   {
    for(Int_t q=0;q<qMax;q++)
    {
     (*fGEBE)(p,q) *= (1.+wPhi*wPt*wEta*(2.*fR0*sqrt(p+1.)/nRP)*cos(fHarmonic*dPhi-2.*q*TMath::Pi()/qMax));
    }
   }
   // Fill the profile to calculate <<w^2>>: 
   fAverageOfSquaredWeight->Fill(0.5,pow(wPhi*wPt*wEta,2.),1.); 
  } // end of if(aftsTrack && aftsTrack->InRPSelection())
 } // end of for(Int_t i=0;i<nPrim;i++) 
 
 // Cross check # of RPs:
 if(anEvent && (crossCheckRP != anEvent->GetEventNSelTracksRP()))
 {
  cout<<endl; 
  cout<<"WARNING (GFC): crossCheckRP != nRP in GFC::Make(). Something is wrong with RP flagging !!!!"<<endl;
  cout<<endl; 
  exit(0);
 }
 
 // Storing the value of G[p][q] in 2D profile in order to get eventually the avarage <G[p][q]>:
 // Determine first the event weight for G[p][q]:
 // (to be improved - this can be implemented much better, this shall be executed only once out of Make(), eventWeight should be a data member)
 Double_t eventWeight = 0.;
 if(!strcmp(fMultiplicityWeight->Data(),"unit"))
 {
  eventWeight = 1.;
 } else if(!strcmp(fMultiplicityWeight->Data(),"multiplicity"))
   {
    eventWeight = anEvent->GetEventNSelTracksRP();           
   }
 // Store G[p][q] weighted appropriately:
 for(Int_t p=0;p<pMax;p++)
 {
  for(Int_t q=0;q<qMax;q++)
  {
   fReferenceFlowGenFun->Fill((Double_t)p,(Double_t)q,(*fGEBE)(p,q),eventWeight); 
   if(fCalculateVsMultiplicity){fReferenceFlowGenFunVsM->Fill(nRP+0.5,(Double_t)p,(Double_t)q,(*fGEBE)(p,q),eventWeight);}
  }
 } 
 
} // end of void AliFlowAnalysisWithCumulants::FillGeneratingFunctionForReferenceFlow(AliFlowEventSimple* anEvent)
 
//================================================================================================================

void AliFlowAnalysisWithCumulants::FillQvectorComponents(AliFlowEventSimple* anEvent)
{
 // Fill components of Q-vector for current event (needed for error calculation).
 
 // Remark: Components are stored in profile fQvectorComponents whose binning is organized as follows:
 //  1st bin: Q_x
 //  2nd bin: Q_y
 //  3rd bin: (Q_x)^2
 //  4th bin: (Q_y)^2
 
 AliFlowVector afv;
 afv.Set(0.,0.);
 afv.SetMult(0);
 
 Int_t n = 2; // to be removed
 
 if(anEvent)
 {
  afv = anEvent->GetQ(1*n,fWeightsList,fUsePhiWeights,fUsePtWeights,fUseEtaWeights); // get the Q-vector for this event
  fQvectorComponents->Fill(0.5,afv.X(),1.); // in the 1st bin fill Q_x
  fQvectorComponents->Fill(1.5,afv.Y(),1.); // in the 2nd bin fill Q_y
  fQvectorComponents->Fill(2.5,pow(afv.X(),2.),1.); // in the 3rd bin fill (Q_x)^2
  fQvectorComponents->Fill(3.5,pow(afv.Y(),2.),1.); // in the 4th bin fill (Q_y)^2
 }
 
} // end of void AliFlowAnalysisWithCumulants::FillQvectorComponents(AliFlowEventSimple* anEvent)

//================================================================================================================

void AliFlowAnalysisWithCumulants::FillGeneratingFunctionForDiffFlow(AliFlowEventSimple* anEvent)
{ 
 // Fill generating function for differential flow for the current event.
 
 // Remark 0: Generating function D[b][p][q] is a complex number => real and imaginary part are calculated separately
 //           (b denotes pt or eta bin);
 // Remark 1: Note that bellow G[p][q] is needed, the value of generating function for reference flow for the CURRENT event.
 //           This values is obtained in method FillGeneratingFunctionForReferenceFlow() as TMatrixD fGEBE;
 // Remark 2: Results for D[b][p][q] are stored in 3D profiles fDiffFlowGenFun[0=Re,1=Im][0=RP,1=POI][0=pt,1=eta] in order to 
 //           automatically get <Re(D[b][p][q])> and <Im(D[b][p][q])> at the end of the day.
 
 // Particle variables and weights:
 Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
 Double_t dPt  = 0.; // transverse momentum
 Double_t dEta = 0.; // pseudorapidity
 Double_t wPhi = 1.; // phi weight
 Double_t wPt  = 1.; // pt weight
 Double_t wEta = 1.; // eta weight
 
 // pMax and qMax:
 Int_t pMax = fGEBE->GetNrows();
 Int_t qMax = fGEBE->GetNcols(); 

 Int_t nPrim = anEvent->NumberOfTracks(); // nPrim = total number of primary tracks, i.e. nPrim = nRP + nPOI, where:
                                          // nRP   = # of particles used to determine the reaction plane;
                                          // nPOI  = # of particles of interest for a detailed flow analysis.
 
 Int_t nRP = anEvent->GetEventNSelTracksRP(); // nRP = # of particles used to determine the reaction plane
       
 // Start the second loop over event in order to evaluate the generating function D[b][p][q] for differential flow: 
 for(Int_t i=0;i<nPrim;i++)
 {
  AliFlowTrackSimple *aftsTrack = anEvent->GetTrack(i);
  if(aftsTrack)
  {
   if(!(aftsTrack->InRPSelection() || aftsTrack->InPOISelection())) continue;
   // Differential flow of POIs:
   if(aftsTrack->InPOISelection())
   {
    // Get azimuthal angle, momentum and pseudorapidity of a particle:
    dPhi = aftsTrack->Phi();
    dPt  = aftsTrack->Pt();
    dEta = aftsTrack->Eta();
    Double_t ptEta[2] = {dPt,dEta};    
   
    // Count number of POIs in pt/eta bin:
    for(Int_t pe=0;pe<2;pe++)
    { 
     fNoOfParticlesInBin[1][pe]->Fill(ptEta[pe],ptEta[pe],1.);
    }
  
    if(!(aftsTrack->InRPSelection())) // particle was flagged only as POI 
    {
     // Fill generating function:
     for(Int_t p=0;p<pMax;p++)
     {
      for(Int_t q=0;q<qMax;q++)
      {
       for(Int_t ri=0;ri<2;ri++)
       {
        for(Int_t pe=0;pe<2;pe++)
        {
         if(ri==0) // Real part (to be improved - this can be implemented better)
         {
          fDiffFlowGenFun[ri][1][pe]->Fill(ptEta[pe],(Double_t)p,(Double_t)q, // to be improved - hardwired weight 1. in the line bellow
                                       (*fGEBE)(p,q)*cos(fMultiple*fHarmonic*dPhi),1.);
         } 
         else if(ri==1) // Imaginary part (to be improved - this can be implemented better)
         {
          fDiffFlowGenFun[ri][1][pe]->Fill(ptEta[pe],(Double_t)p,(Double_t)q, // to be improved - hardwired weight 1. in the line bellow
                                       (*fGEBE)(p,q)*sin(fMultiple*fHarmonic*dPhi),1.);
         }
        } // end of for(Int_t pe=0;pe<2;pe++)
       } // end of for(Int_t ri=0;ri<2;ri++) 
      } // end of for(Int_t q=0;q<qMax;q++)
     } // end of for(Int_t p=0;p<pMax;p++)       
    } // end of if(!(aftsTrack->InRPSelection())) // particle was flagged only as POI 
    else if(aftsTrack->InRPSelection()) // particle was flagged both as RP and POI 
    {
     // If particle weights were used, get them:
     if(fUsePhiWeights && fnBinsPhi) // determine phi weight for this particle:
     {
      wPhi = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*fnBinsPhi/TMath::TwoPi())));
     }
     if(fUsePtWeights && fnBinsPt) // determine pt weight for this particle:
     {
      wPt = fPtWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/fPtBinWidth))); 
     }              
     if(fUseEtaWeights && fEtaBinWidth) // determine eta weight for this particle: 
     {
      wEta = fEtaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/fEtaBinWidth))); 
     }    
     // Fill generating function:
     for(Int_t p=0;p<pMax;p++)
     {
      for(Int_t q=0;q<qMax;q++)
      {
       for(Int_t ri=0;ri<2;ri++)
       {
        for(Int_t pe=0;pe<2;pe++)
        {
         if(ri==0) // Real part (to be improved - this can be implemented better)
         {
          fDiffFlowGenFun[ri][1][pe]->Fill(ptEta[pe],(Double_t)p,(Double_t)q, // to be improved - hardwired weight 1. in the line bellow
                                      (*fGEBE)(p,q)*cos(fMultiple*fHarmonic*dPhi)/(1.+wPhi*wPt*wEta*(2.*fR0*sqrt(p+1.)/nRP)*cos(fHarmonic*dPhi-2.*q*TMath::Pi()/qMax)),1.);
         } 
         else if(ri==1) // Imaginary part (to be improved - this can be implemented better)
         {
          fDiffFlowGenFun[ri][1][pe]->Fill(ptEta[pe],(Double_t)p,(Double_t)q, // to be improved - hardwired weight 1. in the line bellow
                                      (*fGEBE)(p,q)*sin(fMultiple*fHarmonic*dPhi)/(1.+wPhi*wPt*wEta*(2.*fR0*sqrt(p+1.)/nRP)*cos(fHarmonic*dPhi-2.*q*TMath::Pi()/qMax)),1.);
         }
        } // end of for(Int_t pe=0;pe<2;pe++)
       } // end of for(Int_t ri=0;ri<2;ri++) 
      } // end of for(Int_t q=0;q<qMax;q++)
     } // end of for(Int_t p=0;p<pMax;p++)
    } // end of else if (aftsTrack->InRPSelection()) // particle was flagged both as RP and POI 
   } // end of if(aftsTrack->InPOISelection())
   // Differential flow of RPs:
   if(aftsTrack->InRPSelection()) 
   {
    // Get azimuthal angle, momentum and pseudorapidity of a particle:
    dPhi = aftsTrack->Phi();
    dPt  = aftsTrack->Pt();
    dEta = aftsTrack->Eta();
    Double_t ptEta[2] = {dPt,dEta}; 
    
    // Count number of RPs in pt/eta bin:
    for(Int_t pe=0;pe<2;pe++)
    { 
     fNoOfParticlesInBin[0][pe]->Fill(ptEta[pe],ptEta[pe],1.);
    }
    
    // If particle weights were used, get them:
    if(fUsePhiWeights && fnBinsPhi) // determine phi weight for this particle:
    {
     wPhi = fPhiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*fnBinsPhi/TMath::TwoPi())));
    }
    if(fUsePtWeights && fnBinsPt) // determine pt weight for this particle: 
    {
     wPt = fPtWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-fPtMin)/fPtBinWidth))); 
    }              
    if(fUseEtaWeights && fEtaBinWidth) // determine eta weight for this particle: 
    {
     wEta = fEtaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-fEtaMin)/fEtaBinWidth))); 
    }    
    // Fill generating function:
    for(Int_t p=0;p<pMax;p++)
    {
     for(Int_t q=0;q<qMax;q++)
     {
      for(Int_t ri=0;ri<2;ri++)
      {
       for(Int_t pe=0;pe<2;pe++)
       {
        if(ri==0) // Real part (to be improved - this can be implemented better)
        {
         fDiffFlowGenFun[ri][0][pe]->Fill(ptEta[pe],(Double_t)p,(Double_t)q, // to be improved - hardwired weight 1. in the line bellow
                                     (*fGEBE)(p,q)*cos(fMultiple*fHarmonic*dPhi)/(1.+wPhi*wPt*wEta*(2.*fR0*sqrt(p+1.)/nRP)*cos(fHarmonic*dPhi-2.*q*TMath::Pi()/qMax)),1.);
        } 
        else if(ri==1) // Imaginary part (to be improved - this can be implemented better)
        {
         fDiffFlowGenFun[ri][0][pe]->Fill(ptEta[pe],(Double_t)p,(Double_t)q, // to be improved - hardwired weight 1. in the line bellow
                                     (*fGEBE)(p,q)*sin(fMultiple*fHarmonic*dPhi)/(1.+wPhi*wPt*wEta*(2.*fR0*sqrt(p+1.)/nRP)*cos(fHarmonic*dPhi-2.*q*TMath::Pi()/qMax)),1.);
        }
       } // end of for(Int_t pe=0;pe<2;pe++)
      } // end of for(Int_t ri=0;ri<2;ri++) 
     } // end of for(Int_t q=0;q<qMax;q++)
    } // end of for(Int_t p=0;p<pMax;p++)
   } // end of if(aftsTrack->InRPSelection()) 
  } // end of if(aftsTrack)  
 } // end of for(Int_t i=0;i<nPrim;i++)
 
} // end of void AliFlowAnalysisWithCumulants::FillGeneratingFunctionForDiffFlow(AliFlowEventSimple* anEvent)

//================================================================================================================

void AliFlowAnalysisWithCumulants::GetOutputHistograms(TList *outputListHistos) 
{
 // Get pointers to all objects saved in the output file.
 
 if(outputListHistos) 
 {
  this->SetHistList(outputListHistos);
  if(!fHistList)
  {
   cout<<endl;
   cout<<" WARNING (GFC): fHistList is NULL in AFAWGFC::GOH() !!!!"<<endl;
   cout<<endl;
   exit(0);
  }
  this->GetPointersForBaseHistograms();
  this->AccessSettings();
  this->GetPointersForCommonControlHistograms();
  this->GetPointersForCommonResultsHistograms();
  this->GetPointersForReferenceFlowObjects();
  this->GetPointersForDiffFlowObjects();
  if(fTuneParameters){this->GetPointersForTuningObjects();}
 } else 
   {
    cout<<endl;
    cout<<" WARNING (GFC): outputListHistos is NULL in AFAWGFC::GOH() !!!!"<<endl;
    cout<<endl;
    exit(0);
   }

} // end of void AliFlowAnalysisWithCumulants::GetOutputHistograms(TList *outputListHistos) 

//================================================================================================================

void AliFlowAnalysisWithCumulants::GetPointersForBaseHistograms() 
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
    cout<<" WARNING (GFC): analysisSettings is NULL in AFAWGFC::GPFBH() !!!!"<<endl;
    cout<<endl;
    exit(0);  
   }
   
} // end of void AliFlowAnalysisWithCumulants::GetPointersForBaseHistograms()

//================================================================================================================

void AliFlowAnalysisWithCumulants::GetPointersForCommonControlHistograms() 
{
 // Get pointers for common control histograms.
 
 TString commonHistsName = "AliFlowCommonHistGFC";
 AliFlowCommonHist *commonHist = dynamic_cast<AliFlowCommonHist*>(fHistList->FindObject(commonHistsName.Data()));
 if(commonHist) 
 {
  this->SetCommonHists(commonHist); 
 } else
   {
    cout<<endl;
    cout<<" WARNING (GFC): commonHist is NULL in AFAWGFC::GPFCH() !!!!"<<endl;
    cout<<endl;
    exit(0);  
   }
 
} // end of void AliFlowAnalysisWithCumulants::GetPointersForCommonControlHistograms()

//================================================================================================================

void AliFlowAnalysisWithCumulants::GetPointersForCommonResultsHistograms() 
{
 // Get pointers for common results histograms.

 TString commonHistResults2ndOrderName = "AliFlowCommonHistResults2ndOrderGFC"; 
 AliFlowCommonHistResults *commonHistRes2nd = dynamic_cast<AliFlowCommonHistResults*> 
                                              (fHistList->FindObject(commonHistResults2ndOrderName.Data()));
 if(commonHistRes2nd) 
 {
  this->SetCommonHistsResults2nd(commonHistRes2nd);   
 } else
   {
    cout<<endl;
    cout<<" WARNING (GFC): commonHistRes2nd is NULL in AFAWGFC::GPFCRH() !!!!"<<endl;
    cout<<endl;
    exit(0);  
   }
 TString commonHistResults4thOrderName = "AliFlowCommonHistResults4thOrderGFC";
 AliFlowCommonHistResults *commonHistRes4th = dynamic_cast<AliFlowCommonHistResults*>
                                              (fHistList->FindObject(commonHistResults4thOrderName.Data()));
 if(commonHistRes4th)
 { 
  this->SetCommonHistsResults4th(commonHistRes4th);  
 } else
   {
    cout<<endl;
    cout<<" WARNING (GFC): commonHistRes4th is NULL in AFAWGFC::GPFCRH() !!!!"<<endl;
    cout<<endl;
    exit(0);  
   }
 TString commonHistResults6thOrderName = "AliFlowCommonHistResults6thOrderGFC";
 AliFlowCommonHistResults *commonHistRes6th = dynamic_cast<AliFlowCommonHistResults*>
                                              (fHistList->FindObject(commonHistResults6thOrderName.Data()));
 if(commonHistRes6th)
 { 
  this->SetCommonHistsResults6th(commonHistRes6th);  
 } else
   {
    cout<<endl;
    cout<<" WARNING (GFC): commonHistRes6th is NULL in AFAWGFC::GPFCRH() !!!!"<<endl;
    cout<<endl;
    exit(0);  
   }
 TString commonHistResults8thOrderName = "AliFlowCommonHistResults8thOrderGFC";
 AliFlowCommonHistResults *commonHistRes8th = dynamic_cast<AliFlowCommonHistResults*>
                                              (fHistList->FindObject(commonHistResults8thOrderName.Data()));  
 if(commonHistRes8th)
 {
  this->SetCommonHistsResults8th(commonHistRes8th);
 } else
   {
    cout<<endl;
    cout<<" WARNING (GFC): commonHistRes8th is NULL in AFAWGFC::GPFCRH() !!!!"<<endl;
    cout<<endl;
    exit(0);  
   } 

} // end of void AliFlowAnalysisWithCumulants::GetPointersForCommonResultsHistograms() 

//================================================================================================================

void AliFlowAnalysisWithCumulants::GetPointersForTuningObjects()
{
 // Get pointers to all objects used for tuning.
 
 // a) Get pointers to all lists relevant for tuning;
 // b) Get pointer to profile holding flags for tuning;
 // c) Get pointers to all objects in the list fTuningProfiles;
 // d) Get pointers to all objects in the list fTuningResults.

 // a) Get pointers to all lists relevant for tuning:
 TList *tuningList = dynamic_cast<TList*>(fHistList->FindObject("Tuning"));
 if(!tuningList) 
 {
  cout<<endl; 
  cout<<"WARNING (GFC): uningList is NULL in AFAWGFC::GPFTO() !!!!"<<endl;
  cout<<endl; 
  exit(0); 
 }  
 TList *tuningProfiles = dynamic_cast<TList*>(tuningList->FindObject("Profiles"));
 if(!tuningProfiles)  
 {
  cout<<endl; 
  cout<<"WARNING (GFC): tuningProfiles is NULL in AFAWGFC::GPFTO() !!!!"<<endl;
  cout<<endl; 
  exit(0); 
 }
 TList *tuningResults = dynamic_cast<TList*>(tuningList->FindObject("Results"));
 if(!tuningResults)  
 {
  cout<<endl; 
  cout<<"WARNING (GFC): tuningResults is NULL in AFAWGFC::GPFTO() !!!!"<<endl;
  cout<<endl; 
  exit(0); 
 }

 // b) Get pointer to profile holding flags for tuning:
 TString tuningFlagsName = "fTuningFlags";
 TProfile *tuningFlags = dynamic_cast<TProfile*>(tuningList->FindObject(tuningFlagsName.Data())); 
 if(tuningFlags)
 {
  this->SetTuningFlags(tuningFlags);  
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): tuningFlags is NULL in AFAWGFC::GPFTO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   }
   
 // c) Get pointers to all objects in the list fTuningProfiles:
 // Generating function for different tuning parameters:
 TProfile2D *tuningGenFun[10][5] = {{NULL}};
 for(Int_t r=0;r<10;r++)
 { 
  for(Int_t pq=0;pq<5;pq++)
  {
   tuningGenFun[r][pq] = dynamic_cast<TProfile2D*>(tuningProfiles->FindObject(Form("fTuningGenFun (r_{0,%i}, pq set %i)",r,pq)));    
   if(tuningGenFun[r][pq])
   {
    this->SetTuningGenFun(tuningGenFun[r][pq],r,pq);  
   } else 
     {
      cout<<endl; 
      cout<<"WARNING (GFC): "<<Form("tuningGenFun[%i][%i]",r,pq)<<" is NULL in AFAWGFC::GPFTO() !!!!"<<endl;
      cout<<endl; 
      exit(0); 
     }
  } // end of for(Int_t pq=0;pq<5;pq++)
 } // end of for(Int_t r=0;r<10;r++)
 // Average multiplicities for events with nRPs >= cuttof 
 TProfile *tuningAvM = dynamic_cast<TProfile*>(tuningProfiles->FindObject("fTuningAvM"));   
 if(tuningAvM)
 {
  this->SetTuningAvM(tuningAvM);  
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): tuningAvM is NULL in AFAWGFC::GPFTO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   }
 
 // d) Get pointers to all objects in the list fTuningResults.
 // Cumulants for reference flow for 10 different r0s and 5 different sets of (pmax,qmax): 
 TH1D *tuningCumulants[10][5] = {{NULL}};
 for(Int_t r=0;r<10;r++)
 { 
  for(Int_t pq=0;pq<5;pq++)
  {
   tuningCumulants[r][pq] = dynamic_cast<TH1D*>(tuningResults->FindObject(Form("fTuningCumulants (r_{0,%i}, pq set %i)",r,pq)));       
   if(tuningCumulants[r][pq])
   {
    this->SetTuningCumulants(tuningCumulants[r][pq],r,pq);  
   } else 
     {
      cout<<endl; 
      cout<<"WARNING (GFC): "<<Form("tuningCumulants[%i][%i]",r,pq)<<" is NULL in AFAWGFC::GPFTO() !!!!"<<endl;
      cout<<endl; 
      exit(0); 
     }
  } // end of for(Int_t pq=0;pq<5;pq++)
 } // end of for(Int_t r=0;r<10;r++)
 // Reference flow for 10 different r0s and 5 different sets of (pmax,qmax):    
 TH1D *tuningFlow[10][5] = {{NULL}};
 for(Int_t r=0;r<10;r++)
 { 
  for(Int_t pq=0;pq<5;pq++)
  {
   tuningFlow[r][pq] = dynamic_cast<TH1D*>(tuningResults->FindObject(Form("fTuningFlow (r_{0,%i}, pq set %i)",r,pq))); 
   if(tuningFlow[r][pq])
   {
    this->SetTuningFlow(tuningFlow[r][pq],r,pq);  
   } else 
     {
      cout<<endl; 
      cout<<"WARNING (GFC): "<<Form("tuningFlow[%i][%i]",r,pq)<<" is NULL in AFAWGFC::GPFTO() !!!!"<<endl;
      cout<<endl; 
      exit(0); 
     }
  } // end of for(Int_t pq=0;pq<5;pq++)
 } // end of for(Int_t r=0;r<10;r++)

} // end of void AliFlowAnalysisWithCumulants::GetPointersForTuningObjects()

//================================================================================================================

void AliFlowAnalysisWithCumulants::GetPointersForReferenceFlowObjects()
{
 // Get pointers for all objects relevant for calculation of reference flow.
  
 // a) Get pointers to all lists relevant for reference flow;
 // b) Get pointer to profile holding flags;
 // c) Get pointers to all objects in the list fReferenceFlowProfiles;
 // d) Get pointers to all objects in the list fReferenceFlowResults;
 // e) Get pointers for all objects relevant for calculation of reference flow versus multiplicity.

 // a) Get pointers to all lists relevant for reference flow:
 TList *referenceFlowList = dynamic_cast<TList*>(fHistList->FindObject("Reference Flow"));
 if(!referenceFlowList) 
 {
  cout<<endl; 
  cout<<"WARNING (GFC): referenceFlowList is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
  cout<<endl; 
  exit(0); 
 }  
 TList *referenceFlowProfiles = dynamic_cast<TList*>(referenceFlowList->FindObject("Profiles"));
 if(!referenceFlowProfiles)  
 {
  cout<<endl; 
  cout<<"WARNING (GFC): referenceFlowProfiles is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
  cout<<endl; 
  exit(0); 
 }
 TList *referenceFlowResults = dynamic_cast<TList*>(referenceFlowList->FindObject("Results"));
 if(!referenceFlowResults)  
 {
  cout<<endl; 
  cout<<"WARNING (GFC): referenceFlowResults is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
  cout<<endl; 
  exit(0); 
 }
 
 // b) Get pointer to profile holding flags:
 TString referenceFlowFlagsName = "fReferenceFlowFlags";
 TProfile *referenceFlowFlags = dynamic_cast<TProfile*>(referenceFlowList->FindObject(referenceFlowFlagsName.Data()));
 if(referenceFlowFlags)
 {
  this->SetReferenceFlowFlags(referenceFlowFlags);  
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): referenceFlowFlags is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   }
 
 // c) Get pointers to all objects in the list fReferenceFlowProfiles:
 TString referenceFlowGenFunName = "fReferenceFlowGenFun";
 TProfile2D *referenceFlowGenFun = dynamic_cast<TProfile2D*>(referenceFlowProfiles->FindObject(referenceFlowGenFunName.Data())); 
 if(referenceFlowGenFun)
 {
  this->SetReferenceFlowGenFun(referenceFlowGenFun);  
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): referenceFlowGenFun is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   }
 // Averages of various Q-vector components:
 TString qvectorComponentsName = "fQvectorComponents";
 TProfile *qvectorComponents = dynamic_cast<TProfile*>(referenceFlowProfiles->FindObject(qvectorComponentsName.Data()));
 if(qvectorComponents)
 {
  this->SetQvectorComponents(qvectorComponents);  
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): qvectorComponents is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   }
 // <<w^2>>, where w = wPhi*wPt*wEta:
 TString averageOfSquaredWeightName = "fAverageOfSquaredWeight";
 TProfile *averageOfSquaredWeight = dynamic_cast<TProfile*>(referenceFlowProfiles->FindObject(averageOfSquaredWeightName.Data()));
 if(averageOfSquaredWeight)
 {
  this->SetAverageOfSquaredWeight(averageOfSquaredWeight);
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): averageOfSquaredWeight is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   } 
  
 // d) Get pointers to all objects in the list fReferenceFlowResults:
 // Final results for isotropic cumulants for reference flow:
 TString referenceFlowCumulantsName = "fReferenceFlowCumulants";
 TH1D *referenceFlowCumulants = dynamic_cast<TH1D*>(referenceFlowResults->FindObject(referenceFlowCumulantsName.Data()));
 if(referenceFlowCumulants) 
 {
  this->SetReferenceFlowCumulants(referenceFlowCumulants);
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): referenceFlowCumulants is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   }  
 // Final results for reference flow:
 TString referenceFlowName = "fReferenceFlow";
 TH1D *referenceFlow = dynamic_cast<TH1D*>(referenceFlowResults->FindObject(referenceFlowName.Data()));
 if(referenceFlow)
 {
  this->SetReferenceFlow(referenceFlow); 
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): referenceFlow is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   } 
 // Final results for resolution:
 TString chiName = "fChi";
 TH1D *chi = dynamic_cast<TH1D*>(referenceFlowResults->FindObject(chiName.Data()));
 if(chi)
 {
  this->SetChi(chi); 
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): chi is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   } 
 
 // e) Get pointers for all objects relevant for calculation of reference flow versus multiplicity:
 if(!fCalculateVsMultiplicity) {return;}
 // All-event average of the generating function used to calculate reference flow vs multiplicity:
 TString referenceFlowGenFunVsMName = "fReferenceFlowGenFunVsM";
 TProfile3D *referenceFlowGenFunVsM = dynamic_cast<TProfile3D*>(referenceFlowProfiles->FindObject(referenceFlowGenFunVsMName.Data())); 
 if(referenceFlowGenFunVsM)
 {
  this->SetReferenceFlowGenFunVsM(referenceFlowGenFunVsM);  
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): referenceFlowGenFunVsM is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   }
 // Averages of various Q-vector components versus multiplicity:
 TString qvectorComponentsVsMName = "fQvectorComponentsVsM";
 TProfile2D *qvectorComponentsVsM = dynamic_cast<TProfile2D*>(referenceFlowProfiles->FindObject(qvectorComponentsVsMName.Data()));
 if(qvectorComponentsVsM)
 {
  this->SetQvectorComponentsVsM(qvectorComponentsVsM);  
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): qvectorComponentsVsM is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   }
 // <<w^2>>, where w = wPhi*wPt*wEta versus multiplicity:
 TString averageOfSquaredWeightVsMName = "fAverageOfSquaredWeightVsM";
 TProfile2D *averageOfSquaredWeightVsM = dynamic_cast<TProfile2D*>(referenceFlowProfiles->FindObject(averageOfSquaredWeightVsMName.Data()));
 if(averageOfSquaredWeightVsM)
 {
  this->SetAverageOfSquaredWeightVsM(averageOfSquaredWeightVsM);
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): averageOfSquaredWeightVsM is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   } 
 // Final results for reference GF-cumulants versus multiplicity:
 TString cumulantFlag[4] = {"GFC{2}","GFC{4}","GFC{6}","GFC{8}"};
 TString referenceFlowCumulantsVsMName = "fReferenceFlowCumulantsVsM";
 TH1D *referenceFlowCumulantsVsM[4] = {NULL};
 for(Int_t co=0;co<4;co++) // cumulant order
 {
  referenceFlowCumulantsVsM[co] = dynamic_cast<TH1D*>(referenceFlowResults->FindObject(Form("%s, %s",referenceFlowCumulantsVsMName.Data(),cumulantFlag[co].Data())));
  if(referenceFlowCumulantsVsM[co])
  {
   this->SetReferenceFlowCumulantsVsM(referenceFlowCumulantsVsM[co],co);
  } else
    {
     cout<<endl; 
     cout<<"WARNING (GFC): "<<Form("referenceFlowCumulantsVsM[%i]",co)<<" is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
     cout<<endl; 
     exit(0); 
    }   
 } // end of for(Int_t co=0;co<4;co++) // cumulant order
 // <M> vs multiplicity bin: 
 TProfile *avMVsM = dynamic_cast<TProfile*>(referenceFlowProfiles->FindObject("fAvMVsM"));
 if(avMVsM)
 {
  this->SetAvMVsM(avMVsM);
 } else
   {
    cout<<endl; 
    cout<<"WARNING (GFC): avMVsM is NULL in AFAWGFC::GPFRFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   }   

}  // end of void AliFlowAnalysisWithCumulants::GetPointersForReferenceFlowObjects()

//================================================================================================================

void AliFlowAnalysisWithCumulants::GetPointersForDiffFlowObjects()
{
 // Get pointers to all objects relevant for differential flow.
 
 //  a) Define flags locally (to be improved: should I promote flags to data members?);
 //  b) Get pointer to base list for differential flow fDiffFlowList and nested lists fDiffFlowListProfiles and fDiffFlowListResults;
 //  c) Get pointer to profile fDiffFlowFlags holding all flags for differential flow;
 //  d) Get pointers to all profiles in the list fDiffFlowProfiles;
 //  e) Get pointers to all profiles in the list fDiffFlowResults.
  
 // a) Define flags locally (to be improved: should I promote flags to data members?): 
 TString reIm[2] = {"Re","Im"};
 TString rpPoi[2] = {"RP","POI"}; 
 TString ptEta[2] = {"p_{t}","#eta"};
 TString order[4] = {"2nd order","4th order","6th order","8th order"}; 
 //TString differentialFlowIndex[4] = {"v'{2}","v'{4}","v'{6}","v'{8}"};  
  
 // b) Get pointer to base list for differential flow fDiffFlowList and nested lists fDiffFlowListProfiles and fDiffFlowListResults:
 TList *diffFlowList = dynamic_cast<TList*>(fHistList->FindObject("Differential Flow")); // to be improved (hardwired name) 
 if(!diffFlowList)
 { 
  cout<<endl;
  cout<<"WARNING: diffFlowList is NULL in AFAWC::GPFDFO() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 TList *diffFlowProfiles = dynamic_cast<TList*>(diffFlowList->FindObject("Profiles")); // to be improved (hardwired name)
 if(!diffFlowProfiles)
 { 
  cout<<endl;
  cout<<"WARNING: diffFlowProfiles is NULL in AFAWC::GPFDFO() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 TList *diffFlowResults = dynamic_cast<TList*>(diffFlowList->FindObject("Results")); // to be improved (hardwired name)
 if(!diffFlowResults)
 { 
  cout<<endl;
  cout<<"WARNING: diffFlowResults is NULL in AFAWC::GPFDFO() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

 // c) Get pointer to profile holding flags:
 TString diffFlowFlagsName = "fDiffFlowFlags";
 TProfile *diffFlowFlags = dynamic_cast<TProfile*>(diffFlowList->FindObject(diffFlowFlagsName.Data()));
 if(diffFlowFlags)
 {
  this->SetDiffFlowFlags(diffFlowFlags);  
 } else 
   {
    cout<<endl; 
    cout<<"WARNING (GFC): diffFlowFlags is NULL in AFAWGFC::GPFDFO() !!!!"<<endl;
    cout<<endl; 
    exit(0); 
   }
    
 // d) Get pointers to all profiles in the list fDiffFlowListProfiles:
 // Generating functions for differential flow:
 TProfile3D *diffFlowGenFun[2][2][2] = {{{NULL}}};
 for(Int_t ri=0;ri<2;ri++)
 {
  for(Int_t rp=0;rp<2;rp++)
  {
   for(Int_t pe=0;pe<2;pe++)
   {
    diffFlowGenFun[ri][rp][pe] = dynamic_cast<TProfile3D*> // to be improved - harwired name fDiffFlowGenFun in the line bellow
                                 (diffFlowProfiles->FindObject(Form("fDiffFlowGenFun (%s, %s, %s)",reIm[ri].Data(),rpPoi[rp].Data(),ptEta[pe].Data())));
    if(diffFlowGenFun[ri][rp][pe])
    {
     this->SetDiffFlowGenFun(diffFlowGenFun[ri][rp][pe],ri,rp,pe);
    } else 
      {
       cout<<endl; 
       cout<<"WARNING (GFC): "<<Form("diffFlowGenFun[%d][%d][%d]",ri,rp,pe)<<" is NULL in AFAWGFC::GPFDFO() !!!!"<<endl;
       cout<<endl; 
       exit(0); 
      }   
   }
  }   
 } 
 // Number of particles in pt/eta bin for RPs/POIs:
 TProfile *noOfParticlesInBin[2][2] = {{NULL}};
 for(Int_t rp=0;rp<2;rp++)
 {
  for(Int_t pe=0;pe<2;pe++)
  {
   noOfParticlesInBin[rp][pe] = dynamic_cast<TProfile*> // to be improved - harwired name fNoOfParticlesInBin in the line bellow
                                 (diffFlowProfiles->FindObject(Form("fNoOfParticlesInBin (%s, %s)",rpPoi[rp].Data(),ptEta[pe].Data())));  
   if(noOfParticlesInBin[rp][pe])
   {
    this->SetNoOfParticlesInBin(noOfParticlesInBin[rp][pe],rp,pe);
   } else 
     {
      cout<<endl; 
      cout<<"WARNING (GFC): "<<Form("noOfParticlesInBin[%d][%d]",rp,pe)<<" is NULL in AFAWGFC::GPFDFO() !!!!"<<endl;
      cout<<endl; 
      exit(0); 
     }   
  }
 }   
 // Differential cumulants per pt/eta bin for RPs/POIs:
 TH1D *diffFlowCumulants[2][2][4] = {{{NULL}}};
 for(Int_t rp=0;rp<2;rp++)
 {
  for(Int_t pe=0;pe<2;pe++)
  {
   for(Int_t co=0;co<4;co++)
   {
    diffFlowCumulants[rp][pe][co] = dynamic_cast<TH1D*> // to be improved - hardwired name fDiffFlowCumulants in the line bellow
                                    (diffFlowResults->FindObject(Form("fDiffFlowCumulants (%s, %s, %s)",rpPoi[rp].Data(),ptEta[pe].Data(),order[co].Data())));  
    if(diffFlowCumulants[rp][pe][co])
    {
     this->SetDiffFlowCumulants(diffFlowCumulants[rp][pe][co],rp,pe,co);
    } else 
      {
       cout<<endl; 
       cout<<"WARNING (GFC): "<<Form("diffFlowCumulants[%d][%d][%d]",rp,pe,co)<<" is NULL in AFAWGFC::GPFDFO() !!!!"<<endl;
       cout<<endl; 
       exit(0); 
      }
   }      
  }
 }   
 // Differential flow per pt/eta bin for RPs/POIs:
 TH1D *diffFlow[2][2][4] = {{{NULL}}};
 for(Int_t rp=0;rp<2;rp++)
 {
  for(Int_t pe=0;pe<2;pe++)
  {
   for(Int_t co=0;co<4;co++)
   {
    diffFlow[rp][pe][co] = dynamic_cast<TH1D*> // to be improved - hardwired name fDiffFlow in the line bellow
                           (diffFlowResults->FindObject(Form("fDiffFlow (%s, %s, %s)",rpPoi[rp].Data(),ptEta[pe].Data(),order[co].Data())));  
    if(diffFlow[rp][pe][co])
    {
     this->SetDiffFlow(diffFlow[rp][pe][co],rp,pe,co);
    } else 
      {
       cout<<endl; 
       cout<<"WARNING (GFC): "<<Form("diffFlow[%d][%d][%d]",rp,pe,co)<<" is NULL in AFAWGFC::GPFDFO() !!!!"<<endl;
       cout<<endl; 
       exit(0); 
      }
   }      
  }
 }   

} // end of void AliFlowAnalysisWithCumulants::GetPointersForDiffFlowObjects()

//================================================================================================================

void AliFlowAnalysisWithCumulants::CalculateIntegratedFlow(TString rpPoi)
{
 // Calculate final results for integrated flow of RPs and POIs. 
 // (to be improved - this method can be implemented much better)
  
 Int_t rp = 0;

 if(rpPoi == "RP")
 {
  rp = 0;
 } else if(rpPoi == "POI")
   {
    rp = 1;
   } 
          
 // pt yield:    
 TH1F *yieldPt = NULL;
 
 if(rpPoi == "POI")
 {
  yieldPt = (TH1F*)(fCommonHists->GetHistPtPOI())->Clone();
 } else if(rpPoi == "RP")
   {
    yieldPt = (TH1F*)(fCommonHists->GetHistPtRP())->Clone();
   } 
    
 if(!yieldPt)
 {
  printf("\n WARNING (GFC): yieldPt is NULL in AFAWC::CIF() !!!!\n");
  return;
 } 
  
 TH1D *flow2ndPt = (TH1D*)fDiffFlow[rp][0][0]->Clone();
 TH1D *flow4thPt = (TH1D*)fDiffFlow[rp][0][1]->Clone();
 TH1D *flow6thPt = (TH1D*)fDiffFlow[rp][0][2]->Clone();
 TH1D *flow8thPt = (TH1D*)fDiffFlow[rp][0][3]->Clone();
 Double_t dvn2nd = 0., dvn4th = 0., dvn6th = 0., dvn8th = 0.; // differential flow
 Double_t dErrvn2nd = 0., dErrvn4th = 0., dErrvn6th = 0., dErrvn8th = 0.; // error on differential flow
 Double_t dVn2nd = 0., dVn4th = 0., dVn6th = 0., dVn8th = 0.; // integrated flow 
 Double_t dErrVn2nd = 0., dErrVn4th = 0., dErrVn6th = 0., dErrVn8th = 0.; // error on integrated flow
 Double_t dYield = 0.; // pt yield 
 Double_t dSum2nd = 0., dSum4th = 0., dSum6th = 0., dSum8th = 0.; // needed for normalizing integrated flow
 fnBinsPt = flow2ndPt->GetXaxis()->GetNbins(); 
 // looping over pt bins:
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  dvn2nd = flow2ndPt->GetBinContent(p);
  dvn4th = flow4thPt->GetBinContent(p);
  dvn6th = flow6thPt->GetBinContent(p);
  dvn8th = flow8thPt->GetBinContent(p);
  
  dErrvn2nd = flow2ndPt->GetBinError(p);
  dErrvn4th = flow4thPt->GetBinError(p);
  dErrvn6th = flow6thPt->GetBinError(p);
  dErrvn8th = flow8thPt->GetBinError(p);

  dYield = yieldPt->GetBinContent(p);  
  
  dVn2nd += dvn2nd*dYield;
  dVn4th += dvn4th*dYield;
  dVn6th += dvn6th*dYield;
  dVn8th += dvn8th*dYield;
  
  dSum2nd += dYield;
  dSum4th += dYield;
  dSum6th += dYield;
  dSum8th += dYield;
  
  dErrVn2nd += dYield*dYield*dErrvn2nd*dErrvn2nd; 
  dErrVn4th += dYield*dYield*dErrvn4th*dErrvn4th;
  dErrVn6th += dYield*dYield*dErrvn6th*dErrvn6th;
  dErrVn8th += dYield*dYield*dErrvn8th*dErrvn8th;
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)

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
 if(dSum6th) 
 {
  dVn6th /= dSum6th;
  dErrVn6th /= (dSum6th*dSum6th);
  dErrVn6th = TMath::Sqrt(dErrVn6th);
 } 
 if(dSum8th) 
 {
  dVn8th /= dSum8th;
  dErrVn8th /= (dSum8th*dSum8th);
  dErrVn8th = TMath::Sqrt(dErrVn8th);
 } 
  
 // storing the results for integrated flow in common hist results: 
 if(rpPoi == "POI")
 {
  fCommonHistsResults2nd->FillIntegratedFlowPOI(dVn2nd,dErrVn2nd); 
  fCommonHistsResults4th->FillIntegratedFlowPOI(dVn4th,dErrVn4th); 
  fCommonHistsResults6th->FillIntegratedFlowPOI(dVn6th,dErrVn6th); 
  fCommonHistsResults8th->FillIntegratedFlowPOI(dVn8th,dErrVn8th); 
 }
 else if(rpPoi == "RP")
 {
  fCommonHistsResults2nd->FillIntegratedFlowRP(dVn2nd,dErrVn2nd); 
  fCommonHistsResults4th->FillIntegratedFlowRP(dVn4th,dErrVn4th);
  fCommonHistsResults6th->FillIntegratedFlowRP(dVn6th,dErrVn6th); 
  fCommonHistsResults8th->FillIntegratedFlowRP(dVn8th,dErrVn8th); 
 }
 
 delete flow2ndPt;
 delete flow4thPt;
 delete flow6thPt;
 delete flow8thPt; 
 delete yieldPt;
 
} // end of void AliFlowAnalysisWithCumulants::CalculateIntegratedFlow(TString rpPoi)

//================================================================================================================

void AliFlowAnalysisWithCumulants::FillCommonHistResultsForDifferentialFlow(TString rpPoi)
{
 // Fill common result histograms for differential flow.
 // (to be improved - this method can be implemented much better)
 
 Int_t rp = 0;

 if(rpPoi == "RP")
 {
  rp = 0;
 } else if(rpPoi == "POI")
   {
    rp = 1;
   } 
   
 // pt:
 for(Int_t p=1;p<=fnBinsPt;p++)
 {
  // Result:
  Double_t v2 = fDiffFlow[rp][0][0]->GetBinContent(p);
  Double_t v4 = fDiffFlow[rp][0][1]->GetBinContent(p);
  Double_t v6 = fDiffFlow[rp][0][2]->GetBinContent(p);
  Double_t v8 = fDiffFlow[rp][0][3]->GetBinContent(p);
  // Error:
  Double_t v2Error = fDiffFlow[rp][0][0]->GetBinError(p);
  Double_t v4Error = fDiffFlow[rp][0][1]->GetBinError(p);
  Double_t v6Error = fDiffFlow[rp][0][2]->GetBinError(p);
  Double_t v8Error = fDiffFlow[rp][0][3]->GetBinError(p);
  // Fill common hist results:
  if(rpPoi == "RP")
  {
   fCommonHistsResults2nd->FillDifferentialFlowPtRP(p,v2,v2Error);
   fCommonHistsResults4th->FillDifferentialFlowPtRP(p,v4,v4Error);
   fCommonHistsResults6th->FillDifferentialFlowPtRP(p,v6,v6Error);
   fCommonHistsResults8th->FillDifferentialFlowPtRP(p,v8,v8Error);
  } else if(rpPoi == "POI")
    {
     fCommonHistsResults2nd->FillDifferentialFlowPtPOI(p,v2,v2Error);
     fCommonHistsResults4th->FillDifferentialFlowPtPOI(p,v4,v4Error);
     fCommonHistsResults6th->FillDifferentialFlowPtPOI(p,v6,v6Error);
     fCommonHistsResults8th->FillDifferentialFlowPtPOI(p,v8,v8Error);
    }
 } // end of for(Int_t p=1;p<=fnBinsPt;p++)   
 
 // eta:
 for(Int_t e=1;e<=fnBinsEta;e++)
 {
  // Results:
  Double_t v2 = fDiffFlow[rp][1][0]->GetBinContent(e);
  Double_t v4 = fDiffFlow[rp][1][1]->GetBinContent(e);
  Double_t v6 = fDiffFlow[rp][1][2]->GetBinContent(e);
  Double_t v8 = fDiffFlow[rp][1][3]->GetBinContent(e);
  // Errors:
  Double_t v2Error = fDiffFlow[rp][1][0]->GetBinError(e);
  Double_t v4Error = fDiffFlow[rp][1][1]->GetBinError(e);
  Double_t v6Error = fDiffFlow[rp][1][2]->GetBinError(e);
  Double_t v8Error = fDiffFlow[rp][1][3]->GetBinError(e); 
  // Fill common hist results:
  if(rpPoi == "RP")
  {
   fCommonHistsResults2nd->FillDifferentialFlowEtaRP(e,v2,v2Error);
   fCommonHistsResults4th->FillDifferentialFlowEtaRP(e,v4,v4Error);
   fCommonHistsResults6th->FillDifferentialFlowEtaRP(e,v6,v6Error);
   fCommonHistsResults8th->FillDifferentialFlowEtaRP(e,v8,v8Error);
  } else if(rpPoi == "POI")
    {
     fCommonHistsResults2nd->FillDifferentialFlowEtaPOI(e,v2,v2Error);
     fCommonHistsResults4th->FillDifferentialFlowEtaPOI(e,v4,v4Error);
     fCommonHistsResults6th->FillDifferentialFlowEtaPOI(e,v6,v6Error);
     fCommonHistsResults8th->FillDifferentialFlowEtaPOI(e,v8,v8Error);
    }
 } // end of for(Int_t e=1;e<=fnBinsEta;e++)    
 
} // end of void AliFlowAnalysisWithCumulants::FillCommonHistResultsForDifferentialFlow(TString rpPoi)

//================================================================================================================ 
   
void AliFlowAnalysisWithCumulants::CalculateDifferentialFlow(TString rpPoi, TString ptEta)
{
 // Calculate differential flow for RPs/POIs vs pt/eta from cumulants.
 
 Int_t rp = 0; // RP or POI
 Int_t pe = 0; // pt or eta

 if(rpPoi == "RP")
 {
  rp = 0;
 } else if(rpPoi == "POI")
   {
    rp = 1;
   } 
 if(ptEta == "pt")
 {
  pe = 0;
 } else if(ptEta == "eta")
   {
    pe = 1;
   } 
 
 // Reference cumulants:
 Double_t gfc2 = fReferenceFlowCumulants->GetBinContent(1); // reference 2nd order cumulant  
 Double_t gfc4 = fReferenceFlowCumulants->GetBinContent(2); // reference 4th order cumulant   
 Double_t gfc6 = fReferenceFlowCumulants->GetBinContent(3); // reference 6th order cumulant   
 Double_t gfc8 = fReferenceFlowCumulants->GetBinContent(4); // reference 8th order cumulant 
 
 Int_t nBins = fDiffFlowCumulants[rp][pe][0]->GetXaxis()->GetNbins();
 
 for(Int_t b=1;b<=nBins;b++)
 { 
  // Differential cumulants:
  Double_t gfd2 = fDiffFlowCumulants[rp][pe][0]->GetBinContent(b); // differential 2nd order cumulant
  Double_t gfd4 = fDiffFlowCumulants[rp][pe][1]->GetBinContent(b); // differential 4th order cumulant
  Double_t gfd6 = fDiffFlowCumulants[rp][pe][2]->GetBinContent(b); // differential 6th order cumulant
  Double_t gfd8 = fDiffFlowCumulants[rp][pe][3]->GetBinContent(b); // differential 8th order cumulant
  // Differential flow:
  Double_t v2 = 0.; // v'{2,GFC}
  Double_t v4 = 0.; // v'{4,GFC}
  Double_t v6 = 0.; // v'{6,GFC}
  Double_t v8 = 0.; // v'{8,GFC}  
  // 2nd order:
  if(gfc2>0.)
  {
   v2 = gfd2/pow(gfc2,0.5);
   fDiffFlow[rp][pe][0]->SetBinContent(b,v2);
  }
  // 4th order:
  if(gfc4<0.)
  {
   v4 = -gfd4/pow(-gfc4,.75);
   fDiffFlow[rp][pe][1]->SetBinContent(b,v4);
  }
  // 6th order:
  if(gfc6>0.)
  {
   v6 = gfd6/(4.*pow((1./4.)*gfc6,(5./6.)));
   fDiffFlow[rp][pe][2]->SetBinContent(b,v6);   
  }
  // 8th order:
  if(gfc8<0.)
  {
   v8 = -gfd8/(33.*pow(-(1./33.)*gfc8,(7./8.))); 
   fDiffFlow[rp][pe][3]->SetBinContent(b,v8);   
  }
 } // end of for(Int_t b=1;b<=nBins;b++)
 
} // end of void AliFlowAnalysisWithCumulants::CalculateDifferentialFlow(TString rpPoi,TString ptEta)

//================================================================================================================

void AliFlowAnalysisWithCumulants::CalculateDifferentialFlowErrors(TString rpPoi,TString ptEta)
{
 // Calculate errors of differential flow.

 Int_t rp = 0; // RP or POI
 Int_t pe = 0; // pt or eta

 if(rpPoi == "RP")
 {
  rp = 0;
 } else if(rpPoi == "POI")
   {
    rp = 1;
   } 
 if(ptEta == "pt")
 {
  pe = 0;
 } else if(ptEta == "eta")
   {
    pe = 1;
   } 
 
 // Resolution chi:
 Double_t chi2 = fChi->GetBinContent(1);
 Double_t chi4 = fChi->GetBinContent(2);
 //Double_t chi6 = fChi->GetBinContent(3);
 //Double_t chi8 = fChi->GetBinContent(4);
 
 Int_t nBins = fNoOfParticlesInBin[rp][pe]->GetXaxis()->GetNbins();
 for(Int_t b=1;b<=nBins;b++)
 {
  Int_t nParticles = (Int_t)fNoOfParticlesInBin[rp][pe]->GetBinEntries(b);
  // Error of 2nd order estimate:
  if(chi2>0. &&  nParticles>0.)
  {
   Double_t v2Error = pow((1./(2.*nParticles))*((1.+pow(chi2,2.))/pow(chi2,2.)),0.5); 
   fDiffFlow[rp][pe][0]->SetBinError(b,v2Error);
  }
  // Error of 4th order estimate:
  if(chi4>0. &&  nParticles>0.)
  {
   Double_t v4Error = pow((1./(2.*nParticles))*((2.+6.*pow(chi4,2.)+pow(chi4,4.)+pow(chi4,6.))/pow(chi4,6.)),0.5);
   fDiffFlow[rp][pe][1]->SetBinError(b,v4Error);
  }
  // Error of 6th order estimate:
  //if(chi6>0. &&  nParticles>0.)
  //{
  // Double_t v6Error = ... // to be improved - yet to be calculated
   fDiffFlow[rp][pe][2]->SetBinError(b,0.);
  //}
  // Error of 8th order estimate:
  //if(chi8>0. &&  nParticles>0.)
  //{
  // Double_t v8Error = ... // to be improved - yet to be calculated
   fDiffFlow[rp][pe][3]->SetBinError(b,0.);
  //}
 } // end of for(Int_t b=1;b<=nBins;b++)

} // end of void AliFlowAnalysisWithCumulants::CalculateDifferentialFlowErrors(TString rpPoi,TString ptEta)

//================================================================================================================

void AliFlowAnalysisWithCumulants::CalculateCumulantsForDiffFlow(TString rpPoi,TString ptEta)
{
 // Calculate cumulants for differential flow. 
 
 Int_t rp = 0; // RP or POI
 Int_t pe = 0; // pt or eta

 if(rpPoi == "RP")
 {
  rp = 0;
 } else if(rpPoi == "POI")
   {
    rp = 1;
   } 
 if(ptEta == "pt")
 {
  pe = 0;
 } else if(ptEta == "eta")
   {
    pe = 1;
   } 
         
 // [nBins][pMax][qMax]:
 Int_t nBins = fDiffFlowGenFun[0][rp][pe]->GetXaxis()->GetNbins();
 Int_t pMax  = fDiffFlowGenFun[0][rp][pe]->GetYaxis()->GetNbins();
 Int_t qMax  = fDiffFlowGenFun[0][rp][pe]->GetZaxis()->GetNbins();
 // <G[p][q]>
 TMatrixD dAvG(pMax,qMax); 
 dAvG.Zero();
 for(Int_t p=0;p<pMax;p++)
 {
  for(Int_t q=0;q<qMax;q++)
  {
   dAvG(p,q) = fReferenceFlowGenFun->GetBinContent(fReferenceFlowGenFun->GetBin(p+1,q+1));
  }  
 } 
 // Loop over pt/eta bins and calculate differential cumulants:
 for(Int_t b=0;b<nBins;b++)
 {
  Double_t gfc[5] = {0.}; // to be improved (hardwired 5)
  Double_t D[5] = {0.}; // D_{p} in Eq. (11) in Practical guide // to be improved (hardwired 5)
  // ptBinRPNoOfParticles[b]=fPtBinRPNoOfParticles->GetBinEntries(b+1); 
  for(Int_t p=0;p<pMax;p++)
  {
   Double_t tempSum = 0.;
   for(Int_t q=0;q<qMax;q++)
   {
    if(TMath::Abs(dAvG(p,q))>1.e-44)
    {   
     Double_t X = fDiffFlowGenFun[0][rp][pe]->GetBinContent(fDiffFlowGenFun[0][rp][pe]->GetBin(b+1,p+1,q+1))/dAvG(p,q); // see Ollitrault's Practical guide (Eq. 11)
     Double_t Y = fDiffFlowGenFun[1][rp][pe]->GetBinContent(fDiffFlowGenFun[0][rp][pe]->GetBin(b+1,p+1,q+1))/dAvG(p,q); // see Ollitrault's Practical guide (Eq. 11)
     tempSum += cos(fMultiple*2.*q*TMath::Pi()/qMax)*X 
              + sin(fMultiple*2.*q*TMath::Pi()/qMax)*Y;
    }
   }   
   D[p] = (pow(fR0*pow(p+1.0,0.5),fMultiple)/qMax)*tempSum;   
  }   
  gfc[0] = (1./(fR0*fR0))*(5.*D[0]-5.*D[1]+(10./3.)*D[2]-(5./4.)*D[3]+(1./5.)*D[4]); 
  gfc[1] = (1./pow(fR0,4.))*((-77./6.)*D[0]+(107./6.)*D[1]-(13./1.)*D[2]+(61./12.)*D[3]-(5./6.)*D[4]);
  gfc[2] = (1./pow(fR0,6.))*((71./2.)*D[0]-59.*D[1]+49.*D[2]-(41./2.)*D[3]+(7./2.)*D[4]);
  gfc[3] = (1./pow(fR0,8.))*(-84.*D[0]+156.*D[1]-144.*D[2]+66.*D[3]-12.*D[4]);
  // gfc[4] = (1./pow(fR0,10.))*(120.*D[0]-240.*D[1]+240.*D[2]-120.*D[3]+24.*D[4]); // 10th order cumulant (to be improved - where to store it?)
  // Store cumulants:
  for(Int_t co=0;co<4;co++)
  {
   fDiffFlowCumulants[rp][pe][co]->SetBinContent(b+1,gfc[co]);
  } 
 }   
 
} // end of void AliFlowAnalysisWithCumulants::CalculateCumulantsForDiffFlow(TString rpPoi, TString ptEta)

//================================================================================================================

void AliFlowAnalysisWithCumulants::PrintFinalResults(TString type)
{
 // Printing on the screen the final results for reference flow and for integrated flow of RPs and POIs.
  
 Int_t n = fHarmonic; 
  
 Double_t dVn[4] = {0.}; // array to hold Vn{2}, Vn{4}, Vn{6} and Vn{8}   
 Double_t dVnErr[4] = {0.}; // array to hold errors of Vn{2}, Vn{4}, Vn{6} and Vn{8}   
 
 if(type == "RF")
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
     } else
       {
        cout<<endl;
        cout<<" WARNING: Impossible type (can be RF, RP or POI) !!!!"<<endl;
        cout<<"          Results will not be printed on the screen."<<endl;
        cout<<endl;
        exit(0);
       }
 
 TString title = " flow estimates from GF-cumulants"; 
 TString subtitle = "      ("; 
 
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
  cout<<"  v_"<<n<<"{"<<2*(i+1)<<"} = "<<dVn[i]<<" +/- "<<dVnErr[i]<<endl;
 }
  
 cout<<endl;
 if(type == "RF")
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
  
} // end of AliFlowAnalysisWithCumulants::PrintFinalResults(TString type);

//================================================================================================================

void AliFlowAnalysisWithCumulants::FillCommonHistResultsForReferenceFlow()
{
 // Fill in AliFlowCommonHistResults dedicated histograms for reference flow. 
 
 // Results: 
 Double_t v2 = fReferenceFlow->GetBinContent(1);
 Double_t v4 = fReferenceFlow->GetBinContent(2);
 Double_t v6 = fReferenceFlow->GetBinContent(3);
 Double_t v8 = fReferenceFlow->GetBinContent(4);
 // Errors: 
 Double_t v2Error = fReferenceFlow->GetBinError(1);
 Double_t v4Error = fReferenceFlow->GetBinError(2);
 Double_t v6Error = fReferenceFlow->GetBinError(3);
 Double_t v8Error = fReferenceFlow->GetBinError(4);
 // Fill results end errors in common hist results:
 fCommonHistsResults2nd->FillIntegratedFlow(v2,v2Error);  
 fCommonHistsResults4th->FillIntegratedFlow(v4,v4Error); 
 fCommonHistsResults6th->FillIntegratedFlow(v6,v6Error); 
 fCommonHistsResults8th->FillIntegratedFlow(v8,v8Error); 
 // Chi:
 Double_t chi2 = fChi->GetBinContent(1);
 Double_t chi4 = fChi->GetBinContent(2);
 Double_t chi6 = fChi->GetBinContent(3);
 Double_t chi8 = fChi->GetBinContent(4);
 // Fill resolution chi in common hist results:
 fCommonHistsResults2nd->FillChi(chi2);
 fCommonHistsResults4th->FillChi(chi4);
 fCommonHistsResults6th->FillChi(chi6);
 fCommonHistsResults8th->FillChi(chi8); 
  
} // end of AliFlowAnalysisWithCumulants::FillCommonHistResultsForReferenceFlow()

//================================================================================================================

void AliFlowAnalysisWithCumulants::CalculateReferenceFlowError()
{
 // Calculate error of reference flow harmonics.
  
 // Generating Function Cumulants:
 Double_t gfc2 = fReferenceFlowCumulants->GetBinContent(1); // GFC{2}  
 Double_t gfc4 = fReferenceFlowCumulants->GetBinContent(2); // GFC{4}  
 Double_t gfc6 = fReferenceFlowCumulants->GetBinContent(3); // GFC{6}  
 Double_t gfc8 = fReferenceFlowCumulants->GetBinContent(4); // GFC{8}
 // Reference flow estimates:
 Double_t v2 = fReferenceFlow->GetBinContent(1); // v{2,GFC}  
 Double_t v4 = fReferenceFlow->GetBinContent(2); // v{4,GFC}  
 Double_t v6 = fReferenceFlow->GetBinContent(3); // v{6,GFC}  
 Double_t v8 = fReferenceFlow->GetBinContent(4); // v{8,GFC}
 // Statistical errors of reference flow estimates:
 Double_t v2Error = 0.; // statistical error of v{2,GFC}  
 Double_t v4Error = 0.; // statistical error of v{4,GFC}  
 Double_t v6Error = 0.; // statistical error of v{6,GFC}  
 Double_t v8Error = 0.; // statistical error of v{8,GFC}
 // Chi:
 Double_t chi2 = 0.;
 Double_t chi4 = 0.;
 Double_t chi6 = 0.;
 Double_t chi8 = 0.;
 // <Q-vector stuff>:
 Double_t dAvQx  = fQvectorComponents->GetBinContent(1); // <Q_x>
 Double_t dAvQy  = fQvectorComponents->GetBinContent(2); // <Q_y>
 Double_t dAvQ2x = fQvectorComponents->GetBinContent(3); // <(Q_x)^2>
 Double_t dAvQ2y = fQvectorComponents->GetBinContent(4); // <(Q_y)^2>
 // <w^2>:
 Double_t dAvw2 = 1.;
 if(fnEvts>0)
 { 
  dAvw2 = fAverageOfSquaredWeight->GetBinContent(1); 
  if(TMath::Abs(dAvw2)<1.e-44) 
  {
   cout<<endl;
   cout<<" WARNING (GFC): Average of squared weight is 0 in GFC. Most probably one of the histograms"<<endl; 
   cout<<"                in the file \"weights.root\" was empty. Nothing will be calculated !!!!"<<endl;
   cout<<endl;
  }
 } 
 // Calculating statistical error of v{2,GFC}:
 if(fnEvts>0. && fAvM>0. && dAvw2>0. && gfc2>=0.)
 { 
  if((dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(gfc2,(1./2.))*(fAvM/dAvw2),2.)>0.))       
  {
   chi2 = (fAvM*v2)/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(v2*fAvM/dAvw2,2.),0.5); 
  } 
  if(TMath::Abs(chi2)>1.e-44)
  {  
   v2Error = pow(((1./(2.*fAvM*fnEvts))*((1.+2.*pow(chi2,2))/(2.*pow(chi2,2)))),0.5);
  }
 } 
 // Calculating statistical error of v{4,GFC}:
 if(fnEvts>0 && fAvM>0 && dAvw2>0 && gfc4<=0.)
 {
  if((dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(-gfc4,(1./4.))*(fAvM/dAvw2),2.)>0.))
  {
   chi4 = (fAvM*v4)/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(v4*fAvM/dAvw2,2.),0.5);
  } 
  if(TMath::Abs(chi4)>1.e-44)
  {
   v4Error = (1./(pow(2.*fAvM*fnEvts,0.5)))*pow((1.+4.*pow(chi4,2)+1.*pow(chi4,4.)+2.*pow(chi4,6.))/(2.*pow(chi4,6.)),0.5);
  }
 } 
 // Calculating statistical error of v{6,GFC}:
 if(fnEvts>0 && fAvM>0 && dAvw2>0 && gfc6>=0.)
 {
  if((dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow((1./4.)*gfc6,(1./6.))*(fAvM/dAvw2),2.)>0.))
  {
   chi6 = (fAvM*v6)/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(v6*fAvM/dAvw2,2.),0.5);
  } 
  if(TMath::Abs(chi6)>1.e-44)
  {
   v6Error = (1./(pow(2.*fAvM*fnEvts,0.5)))*pow((3.+18.*pow(chi6,2)+9.*pow(chi6,4.)+28.*pow(chi6,6.)
              +12.*pow(chi6,8.)+24.*pow(chi6,10.))/(24.*pow(chi6,10.)),0.5);
  } 
 } 
 // Calculating statistical error of v{8,GFC}:
 if(fnEvts>0 && fAvM>0 && dAvw2>0 && gfc8<=0.)
 { 
  if((dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(pow(-(1./33.)*gfc8,(1./8.))*(fAvM/dAvw2),2.)>0.))
  { 
   chi8=(fAvM*v8)/pow(dAvQ2x+dAvQ2y-pow(dAvQx,2.)-pow(dAvQy,2.)-pow(v8*fAvM/dAvw2,2.),0.5);
  }
  if(TMath::Abs(chi8)>1.e-44)
  {
   v8Error = (1./(pow(2.*fAvM*fnEvts,0.5)))*pow((12.+96.*pow(chi8,2.)+72.*pow(chi8,4.)+304.*pow(chi8,6.)
              +257.*pow(chi8,8.)+804.*pow(chi8,10.)+363.*pow(chi8,12.)+726.*pow(chi8,14.))/(726.*pow(chi8,14.)),0.5);
  } 
 } 
  
 // Store errors for reference flow:
 fReferenceFlow->SetBinError(1,v2Error);
 fReferenceFlow->SetBinError(2,v4Error);
 fReferenceFlow->SetBinError(3,v6Error);
 fReferenceFlow->SetBinError(4,v8Error);
 // Store resolution chi:
 fChi->SetBinContent(1,chi2);
 fChi->SetBinContent(2,chi4);
 fChi->SetBinContent(3,chi6);
 fChi->SetBinContent(4,chi8);

} // end of void AliFlowAnalysisWithCumulants::CalculateReferenceFlowError()

//================================================================================================================

void AliFlowAnalysisWithCumulants::CalculateReferenceFlow()
{
 // Calculate from isotropic cumulants reference flow.

 // Generating Function Cumulants:
 Double_t gfc2 = fReferenceFlowCumulants->GetBinContent(1); // GFC{2}  
 Double_t gfc4 = fReferenceFlowCumulants->GetBinContent(2); // GFC{4}  
 Double_t gfc6 = fReferenceFlowCumulants->GetBinContent(3); // GFC{6}  
 Double_t gfc8 = fReferenceFlowCumulants->GetBinContent(4); // GFC{8}
 // Reference flow estimates:
 Double_t v2 = 0.; // v{2,GFC}  
 Double_t v4 = 0.; // v{4,GFC}  
 Double_t v6 = 0.; // v{6,GFC}  
 Double_t v8 = 0.; // v{8,GFC}
 // Calculate reference flow estimates from Q-cumulants: 
 if(gfc2>=0.) v2 = pow(gfc2,1./2.); 
 if(gfc4<=0.) v4 = pow(-1.*gfc4,1./4.); 
 if(gfc6>=0.) v6 = pow((1./4.)*gfc6,1./6.); 
 if(gfc8<=0.) v8 = pow((-1./33.)*gfc8,1./8.); 
 // Store results for reference flow:
 fReferenceFlow->SetBinContent(1,v2);
 fReferenceFlow->SetBinContent(2,v4);
 fReferenceFlow->SetBinContent(3,v6);
 fReferenceFlow->SetBinContent(4,v8);
 
} // end of void AliFlowAnalysisWithCumulants::CalculateReferenceFlow()

//================================================================================================================

void AliFlowAnalysisWithCumulants::CalculateCumulantsForReferenceFlow()
{
 // Calculate cumulants for reference flow.

 Int_t pMax = fReferenceFlowGenFun->GetXaxis()->GetNbins();
 Int_t qMax = fReferenceFlowGenFun->GetYaxis()->GetNbins();
 
 // <G[p][q]>
 TMatrixD dAvG(pMax,qMax); 
 dAvG.Zero();
 Bool_t someAvGEntryIsNegative = kFALSE;
 for(Int_t p=0;p<pMax;p++)
 {
  for(Int_t q=0;q<qMax;q++)
  {
   dAvG(p,q) = fReferenceFlowGenFun->GetBinContent(fReferenceFlowGenFun->GetBin(p+1,q+1));
   if(dAvG(p,q)<0.)
   {
    someAvGEntryIsNegative = kTRUE;
    cout<<endl; 
    cout<<" WARNING: "<<Form("<G[%d][%d]> is negative !!!! GFC results are meaningless.",p,q)<<endl; 
    cout<<endl; 
   }
  }  
 } 
      
 // C[p][q] (generating function for the cumulants)    
 TMatrixD dC(pMax,qMax);
 dC.Zero();
 if(fAvM>0. && !someAvGEntryIsNegative)
 {
  for(Int_t p=0;p<pMax;p++)
  {
   for(Int_t q=0;q<qMax;q++)
   {
    dC(p,q) = fAvM*(pow(dAvG(p,q),(1./fAvM))-1.); 
   }
  }
 }
    
 // Averaging the generating function for cumulants over azimuth
 // in order to eliminate detector effects.
 // <C[p][q]> (Remark: here <> stands for average over azimuth):
 TVectorD dAvC(pMax);
 dAvC.Zero();
 for(Int_t p=0;p<pMax;p++)
 {
  Double_t temp = 0.; 
  for(Int_t q=0;q<qMax;q++)
  {
   temp += 1.*dC(p,q);
  } 
  dAvC[p] = temp/qMax;
 }
 
 // Finally, the isotropic cumulants for reference flow:
 TVectorD cumulant(pMax);
 cumulant.Zero(); 
 cumulant[0] = (-1./(60*fR0*fR0))*((-300.)*dAvC[0]+300.*dAvC[1]-200.*dAvC[2]+75.*dAvC[3]-12.*dAvC[4]);
 cumulant[1] = (-1./(6.*pow(fR0,4.)))*(154.*dAvC[0]-214.*dAvC[1]+156.*dAvC[2]-61.*dAvC[3]+10.*dAvC[4]);
 cumulant[2] = (3./(2.*pow(fR0,6.)))*(71.*dAvC[0]-118.*dAvC[1]+98.*dAvC[2]-41.*dAvC[3]+7.*dAvC[4]);
 cumulant[3] = (-24./pow(fR0,8.))*(14.*dAvC[0]-26.*dAvC[1]+24.*dAvC[2]-11.*dAvC[3]+2.*dAvC[4]);
 cumulant[4] = (120./pow(fR0,10.))*(5.*dAvC[0]-10.*dAvC[1]+10.*dAvC[2]-5.*dAvC[3]+1.*dAvC[4]);
 
 // Store cumulants:
 // Remark: the highest order cumulant is on purpose in the overflow.
 for(Int_t co=0;co<pMax;co++) // cumulant order
 {
  fReferenceFlowCumulants->SetBinContent(co+1,cumulant[co]);
 }
 
 // Calculation versus multiplicity:
 if(!fCalculateVsMultiplicity){return;}
 for(Int_t b=0;b<fnBinsMult;b++)
 {
  fAvM = fAvMVsM->GetBinContent(b+1);
  // <G[p][q]>
  TMatrixD dAvGVsM(pMax,qMax); 
  dAvGVsM.Zero();
  Bool_t someAvGEntryIsNegativeVsM = kFALSE;
  for(Int_t p=0;p<pMax;p++)
  {
   for(Int_t q=0;q<qMax;q++)
   {
    dAvGVsM(p,q) = fReferenceFlowGenFunVsM->GetBinContent(fReferenceFlowGenFunVsM->GetBin(b+1,p+1,q+1));
    if(dAvGVsM(p,q)<0.)
    {
     someAvGEntryIsNegativeVsM = kTRUE;
     cout<<endl; 
     cout<<" WARNING: "<<Form("<G[%d][%d]> is negative !!!! GFC vs multiplicity results are meaningless.",p,q)<<endl; 
     cout<<endl; 
    }
   }  
  } 
      
  // C[p][q] (generating function for the cumulants)    
  TMatrixD dCVsM(pMax,qMax);
  dCVsM.Zero();
  if(fAvM>0. && !someAvGEntryIsNegativeVsM)
  {
   for(Int_t p=0;p<pMax;p++)
   {
    for(Int_t q=0;q<qMax;q++)
    {
     dCVsM(p,q) = fAvM*(pow(dAvGVsM(p,q),(1./fAvM))-1.); 
    }
   }
  }
    
  // Averaging the generating function for cumulants over azimuth
  // in order to eliminate detector effects.
  // <C[p][q]> (Remark: here <> stands for average over azimuth):
  TVectorD dAvCVsM(pMax);
  dAvCVsM.Zero();
  for(Int_t p=0;p<pMax;p++)
  {
   Double_t tempVsM = 0.; 
   for(Int_t q=0;q<qMax;q++)
   {
    tempVsM += 1.*dCVsM(p,q);
   } 
   dAvCVsM[p] = tempVsM/qMax;
  }
  
  // Finally, the isotropic cumulants for reference flow:
  TVectorD cumulantVsM(pMax);
  cumulantVsM.Zero(); 
  cumulantVsM[0] = (-1./(60*fR0*fR0))*((-300.)*dAvCVsM[0]+300.*dAvCVsM[1]-200.*dAvCVsM[2]+75.*dAvCVsM[3]-12.*dAvCVsM[4]);
  cumulantVsM[1] = (-1./(6.*pow(fR0,4.)))*(154.*dAvCVsM[0]-214.*dAvCVsM[1]+156.*dAvCVsM[2]-61.*dAvCVsM[3]+10.*dAvCVsM[4]);
  cumulantVsM[2] = (3./(2.*pow(fR0,6.)))*(71.*dAvCVsM[0]-118.*dAvCVsM[1]+98.*dAvCVsM[2]-41.*dAvCVsM[3]+7.*dAvCVsM[4]);
  cumulantVsM[3] = (-24./pow(fR0,8.))*(14.*dAvCVsM[0]-26.*dAvCVsM[1]+24.*dAvCVsM[2]-11.*dAvCVsM[3]+2.*dAvCVsM[4]);
  cumulantVsM[4] = (120./pow(fR0,10.))*(5.*dAvCVsM[0]-10.*dAvCVsM[1]+10.*dAvCVsM[2]-5.*dAvCVsM[3]+1.*dAvCVsM[4]);
  
  // Store cumulants:
  for(Int_t co=0;co<pMax-1;co++) // cumulant order
  {
   fReferenceFlowCumulantsVsM[co]->SetBinContent(b+1,cumulantVsM[co]);
  } 
 } // end of for(Int_t b=0;b<fnBinsMult;b++)

} // end of void AliFlowAnalysisWithCumulants::CalculateCumulantsForReferenceFlow()

//================================================================================================================

void AliFlowAnalysisWithCumulants::GetAvMultAndNoOfEvts()
{
 // From relevant common control histogram get average multiplicity of RPs and number of events.
 
 fAvM = (Double_t)fCommonHists->GetHistMultRP()->GetMean(); 
 fnEvts = (Int_t)fCommonHists->GetHistMultRP()->GetEntries(); 
 
} // end of void AliFlowAnalysisWithCumulants::GetAvMultAndNoOfEvts()

//================================================================================================================

void AliFlowAnalysisWithCumulants::InitializeArrays()
{
 // Initialize all arrays.
 
 for(Int_t ri=0;ri<2;ri++)
 {
  for(Int_t rp=0;rp<2;rp++)
  {
   for(Int_t pe=0;pe<2;pe++)
   {
    fDiffFlowGenFun[ri][rp][pe] = NULL;
   }
  }   
 } 
 for(Int_t rp=0;rp<2;rp++)
 {
  for(Int_t pe=0;pe<2;pe++)
  {
   fNoOfParticlesInBin[rp][pe] = NULL;
  }
 }   
 for(Int_t rp=0;rp<2;rp++)
 {
  for(Int_t pe=0;pe<2;pe++)
  {
   for(Int_t co=0;co<4;co++)
   {
    fDiffFlowCumulants[rp][pe][co] = NULL;
    fDiffFlow[rp][pe][co] = NULL;
   }
  }
 }   
 for(Int_t i=0;i<3;i++)
 {
  fPrintFinalResults[i] = kTRUE;
 }
 for(Int_t r=0;r<10;r++)
 {
  fTuningR0[r] = 0.;
  for(Int_t pq=0;pq<5;pq++)
  {
   fTuningGenFun[r][pq] = NULL;
   fTuningCumulants[r][pq] = NULL;
   fTuningFlow[r][pq] = NULL;
  }
 }
 for(Int_t co=0;co<4;co++)
 {
  fReferenceFlowCumulantsVsM[co] = NULL;
 }
 
} // end of void AliFlowAnalysisWithCumulants::InitializeArrays()

//================================================================================================================

void AliFlowAnalysisWithCumulants::CrossCheckSettings()
{
 // Cross-check the user settings before starting.
 
 // a) Cross check if the choice for multiplicity weight make sense.
 
 // a) Cross check if the choice for multiplicity weight make sense:
 if(strcmp(fMultiplicityWeight->Data(),"unit") &&
    strcmp(fMultiplicityWeight->Data(),"multiplicity"))
 {
  cout<<endl;
  cout<<"WARNING (GFC): Multiplicity weight can be either \"unit\" or \"multiplicity\"."<<endl;
  cout<<"               Certainly not \""<<fMultiplicityWeight->Data()<<"\"."<<endl;
  cout<<endl;
  exit(0);
 }   

 
 
} // end of void AliFlowAnalysisWithCumulants::CrossCheckSettings()

//================================================================================================================

void AliFlowAnalysisWithCumulants::AccessConstants()
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
  
} // end of void AliFlowAnalysisWithCumulants::AccessConstants()

//================================================================================================================

void AliFlowAnalysisWithCumulants::BookAndFillWeightsHistograms()
{
 // Book and fill histograms which hold phi, pt and eta weights.

 if(!fWeightsList)
 {
  cout<<"WARNING (GFC): fWeightsList is NULL in AFAWGFC::BAFWH() !!!!"<<endl;
  exit(0);  
 }
    
 if(fUsePhiWeights)
 {
  if(fWeightsList->FindObject("phi_weights"))
  {
   fPhiWeights = dynamic_cast<TH1F*>(fWeightsList->FindObject("phi_weights"));
   if(!fPhiWeights){printf("\n WARNING (GFC): !fPhiWeights !!!!\n");exit(0);}
   if(TMath::Abs(fPhiWeights->GetBinWidth(1)-fPhiBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (GFC): Inconsistent binning in histograms for phi-weights throughout the code."<<endl;
    cout<<endl;
    //exit(0);
   }
  } else 
    {
     cout<<endl;
     cout<<"WARNING (GFC): fWeightsList->FindObject(\"phi_weights\") is NULL in AFAWGFC::BAFWH() !!!!"<<endl;
     cout<<endl;
     exit(0);
    }
 } // end of if(fUsePhiWeights)
 
 if(fUsePtWeights) 
 {
  if(fWeightsList->FindObject("pt_weights"))
  {
   fPtWeights = dynamic_cast<TH1D*>(fWeightsList->FindObject("pt_weights"));
   if(!fPtWeights){printf("\n WARNING (GFC): !fPtWeights !!!!\n");exit(0);}
   if(TMath::Abs(fPtWeights->GetBinWidth(1)-fPtBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (GFC): Inconsistent binning in histograms for pt-weights throughout the code."<<endl;
    cout<<endl;
    //exit(0);
   }
  } else 
    {
     cout<<endl;
     cout<<"WARNING (GFC): fWeightsList->FindObject(\"pt_weights\") is NULL in AFAWGFC::BAFWH() !!!!"<<endl;
     cout<<endl;
     exit(0);
    }
 } // end of if(fUsePtWeights)    
 
 if(fUseEtaWeights) 
 {
  if(fWeightsList->FindObject("eta_weights"))
  {
   fEtaWeights = dynamic_cast<TH1D*>(fWeightsList->FindObject("eta_weights"));
   if(!fEtaWeights){printf("\n WARNING (GFC): !fEtaWeights !!!!\n");exit(0);}
   if(TMath::Abs(fEtaWeights->GetBinWidth(1)-fEtaBinWidth)>pow(10.,-6.))
   {
    cout<<endl;
    cout<<"WARNING (GFC): Inconsistent binning in histograms for eta-weights throughout the code."<<endl;
    cout<<endl;
    //exit(0);
   }
  } else 
    {
     cout<<endl;
     cout<<"WARNING (GFC): fUseEtaWeights && fWeightsList->FindObject(\"eta_weights\") is NULL in AFAWGFC::BAFWH() !!!!"<<endl;
     cout<<endl;
     exit(0);
    }
 } // end of if(fUseEtaWeights)
 
} // end of AliFlowAnalysisWithCumulants::BookAndFillWeightsHistograms()

//================================================================================================================

void AliFlowAnalysisWithCumulants::BookEverythingForCalculationVsMultiplicity()
{
 // Book all objects relevant for flow analysis versus multiplicity.
 
 // a) Define constants;
 // b) Book all profiles;
 // c) Book all results.
 
 // a) Define constants and local flags:
 Int_t pMax = 5;     
 Int_t qMax = 11;  
 TString cumulantFlag[4] = {"GFC{2}","GFC{4}","GFC{6}","GFC{8}"};

 // b) Book all profiles:
 // Average of the generating function for reference flow <G[p][q]> versus multiplicity:
 fReferenceFlowGenFunVsM = new TProfile3D("fReferenceFlowGenFunVsM","#LTG[p][q]#GT vs M",fnBinsMult,fMinMult,fMaxMult,pMax,0.,(Double_t)pMax,qMax,0.,(Double_t)qMax);
 fReferenceFlowGenFunVsM->SetXTitle("M");
 fReferenceFlowGenFunVsM->SetYTitle("p");
 fReferenceFlowGenFunVsM->SetZTitle("q");
 fReferenceFlowProfiles->Add(fReferenceFlowGenFunVsM);
 // Averages of Q-vector components versus multiplicity:
 fQvectorComponentsVsM = new TProfile2D("fQvectorComponentsVsM","Averages of Q-vector components",fnBinsMult,fMinMult,fMaxMult,4,0.,4.);
 //fQvectorComponentsVsM->SetLabelSize(0.06);
 fQvectorComponentsVsM->SetMarkerStyle(25);
 fQvectorComponentsVsM->SetXTitle("M");
 fQvectorComponentsVsM->GetYaxis()->SetBinLabel(1,"#LTQ_{x}#GT"); // Q_{x}
 fQvectorComponentsVsM->GetYaxis()->SetBinLabel(2,"#LTQ_{y}#GT"); // Q_{y}
 fQvectorComponentsVsM->GetYaxis()->SetBinLabel(3,"#LTQ_{x}^{2}#GT"); // Q_{x}^{2}
 fQvectorComponentsVsM->GetYaxis()->SetBinLabel(4,"#LTQ_{y}^{2}#GT"); // Q_{y}^{2}
 fReferenceFlowProfiles->Add(fQvectorComponentsVsM);
 // <<w^2>>, where w = wPhi*wPt*wEta versus multiplicity:
 fAverageOfSquaredWeightVsM = new TProfile2D("fAverageOfSquaredWeightVsM","#LT#LTw^{2}#GT#GT",fnBinsMult,fMinMult,fMaxMult,1,0,1);
 fAverageOfSquaredWeightVsM->SetLabelSize(0.06);
 fAverageOfSquaredWeightVsM->SetMarkerStyle(25);
 fAverageOfSquaredWeightVsM->SetLabelOffset(0.01);
 fAverageOfSquaredWeightVsM->GetXaxis()->SetBinLabel(1,"#LT#LTw^{2}#GT#GT");
 fReferenceFlowProfiles->Add(fAverageOfSquaredWeightVsM); 
 // <M> vs multiplicity bin: 
 fAvMVsM = new TProfile("fAvMVsM","#LTM#GT vs M",fnBinsMult,fMinMult,fMaxMult);
 //fAvMVsM->SetLabelSize(0.06);
 fAvMVsM->SetMarkerStyle(25);
 fAvMVsM->SetLabelOffset(0.01);
 fAvMVsM->SetXTitle("M");
 fAvMVsM->SetYTitle("#LTM#GT");
 fReferenceFlowProfiles->Add(fAvMVsM); 

 // c) Book all results:
 // Final results for reference GF-cumulants versus multiplicity:
 TString referenceFlowCumulantsVsMName = "fReferenceFlowCumulantsVsM";
 for(Int_t co=0;co<4;co++) // cumulant order
 {
  fReferenceFlowCumulantsVsM[co] = new TH1D(Form("%s, %s",referenceFlowCumulantsVsMName.Data(),cumulantFlag[co].Data()),
                                            Form("%s vs multipicity",cumulantFlag[co].Data()),
                                            fnBinsMult,fMinMult,fMaxMult);
  fReferenceFlowCumulantsVsM[co]->SetMarkerStyle(25);              
  fReferenceFlowCumulantsVsM[co]->GetXaxis()->SetTitle("M");                                     
  fReferenceFlowCumulantsVsM[co]->GetYaxis()->SetTitle(cumulantFlag[co].Data());  
  fReferenceFlowResults->Add(fReferenceFlowCumulantsVsM[co]);                                    
 } // end of for(Int_t co=0;co<4;co++) // cumulant order
 
} // end of void AliFlowAnalysisWithCumulants::BookEverythingForCalculationVsMultiplicity()

//================================================================================================================

void AliFlowAnalysisWithCumulants::BookEverythingForReferenceFlow()
{
 // Book all objects relevant for calculation of reference flow.
 
 // a) Define static constants for array's boundaries;
 // b) Book profile to hold all flags for reference flow;
 // c) Book all event-by-event quantities;
 // d) Book all profiles;
 // e) Book all histograms.
 
 // a) Define static constants for array's boundaries:
 static const Int_t pMax = 5;     
 static const Int_t qMax = 11;  
 
 // b) Book profile to hold all flags for reference flow:
 TString referenceFlowFlagsName = "fReferenceFlowFlags";
 fReferenceFlowFlags = new TProfile(referenceFlowFlagsName.Data(),"Flags for Reference Flow",2,0,2);
 fReferenceFlowFlags->SetTickLength(-0.01,"Y");
 fReferenceFlowFlags->SetMarkerStyle(25);
 fReferenceFlowFlags->SetLabelSize(0.05);
 fReferenceFlowFlags->SetLabelOffset(0.02,"Y");
 fReferenceFlowFlags->GetXaxis()->SetBinLabel(1,"Particle weights");
 fReferenceFlowFlags->GetXaxis()->SetBinLabel(2,"Event weights");
 fReferenceFlowList->Add(fReferenceFlowFlags);
 
 // c) Book all event-by-event quantities: 
 fGEBE = new TMatrixD(pMax,qMax);

 // d) Book all profiles:
 // Average of the generating function for reference flow <G[p][q]>:
 fReferenceFlowGenFun = new TProfile2D("fReferenceFlowGenFun","#LTG[p][q]#GT",pMax,0.,(Double_t)pMax,qMax,0.,(Double_t)qMax);
 fReferenceFlowGenFun->SetXTitle("p");
 fReferenceFlowGenFun->SetYTitle("q");
 fReferenceFlowProfiles->Add(fReferenceFlowGenFun);
 // Averages of Q-vector components:
 fQvectorComponents = new TProfile("fQvectorComponents","Averages of Q-vector components",4,0.,4.);
 fQvectorComponents->SetLabelSize(0.06);
 fQvectorComponents->SetMarkerStyle(25);
 fQvectorComponents->GetXaxis()->SetBinLabel(1,"#LTQ_{x}#GT"); // Q_{x}
 fQvectorComponents->GetXaxis()->SetBinLabel(2,"#LTQ_{y}#GT"); // Q_{y}
 fQvectorComponents->GetXaxis()->SetBinLabel(3,"#LTQ_{x}^{2}#GT"); // Q_{x}^{2}
 fQvectorComponents->GetXaxis()->SetBinLabel(4,"#LTQ_{y}^{2}#GT"); // Q_{y}^{2}
 fReferenceFlowProfiles->Add(fQvectorComponents);
 // <<w^2>>, where w = wPhi*wPt*wEta:
 fAverageOfSquaredWeight = new TProfile("fAverageOfSquaredWeight","#LT#LTw^{2}#GT#GT",1,0,1);
 fAverageOfSquaredWeight->SetLabelSize(0.06);
 fAverageOfSquaredWeight->SetMarkerStyle(25);
 fAverageOfSquaredWeight->SetLabelOffset(0.01);
 fAverageOfSquaredWeight->GetXaxis()->SetBinLabel(1,"#LT#LTw^{2}#GT#GT");
 fReferenceFlowProfiles->Add(fAverageOfSquaredWeight);
 
 // e) Book all histograms:
 // Final results for isotropic cumulants for reference flow:
 TString referenceFlowCumulantsName = "fReferenceFlowCumulants";
 fReferenceFlowCumulants = new TH1D(referenceFlowCumulantsName.Data(),"Isotropic Generating Function Cumulants for reference flow",4,0,4); // to be improved (hw 4)
 fReferenceFlowCumulants->SetLabelSize(0.05);
 fReferenceFlowCumulants->SetMarkerStyle(25);
 fReferenceFlowCumulants->GetXaxis()->SetBinLabel(1,"GFC{2}");
 fReferenceFlowCumulants->GetXaxis()->SetBinLabel(2,"GFC{4}");
 fReferenceFlowCumulants->GetXaxis()->SetBinLabel(3,"GFC{6}");
 fReferenceFlowCumulants->GetXaxis()->SetBinLabel(4,"GFC{8}");
 fReferenceFlowResults->Add(fReferenceFlowCumulants); 
 // Final results for reference flow:  
 fReferenceFlow = new TH1D("fReferenceFlow","Reference flow",4,0,4); // to be improved (hardwired 4)
 fReferenceFlow->SetLabelSize(0.05);
 fReferenceFlow->SetMarkerStyle(25);
 fReferenceFlow->GetXaxis()->SetBinLabel(1,"v_{n}{2,GFC}");
 fReferenceFlow->GetXaxis()->SetBinLabel(2,"v_{n}{4,GFC}");
 fReferenceFlow->GetXaxis()->SetBinLabel(3,"v_{n}{6,GFC}");
 fReferenceFlow->GetXaxis()->SetBinLabel(4,"v_{n}{8,GFC}");
 fReferenceFlowResults->Add(fReferenceFlow);
 // Final results for resolution:  
 fChi = new TH1D("fChi","Resolution",4,0,4); // to be improved (hardwired 4)
 fChi->SetLabelSize(0.06);
 fChi->SetMarkerStyle(25);
 fChi->GetXaxis()->SetBinLabel(1,"#chi_{2}");
 fChi->GetXaxis()->SetBinLabel(2,"#chi_{4}");
 fChi->GetXaxis()->SetBinLabel(3,"#chi_{6}");
 fChi->GetXaxis()->SetBinLabel(4,"#chi_{8}");
 fReferenceFlowResults->Add(fChi);

} // end of void AliFlowAnalysisWithCumulants::BookEverythingForReferenceFlow()

//================================================================================================================

void AliFlowAnalysisWithCumulants::BookEverythingForTuning()
{
 // Book all objects relevant for tuning.
 
 // a) Define pMax's and qMax's:
 // b) Book profile to hold all tuning parameters and flags;
 // c) Book all profiles;
 // d) Book all histograms.
 
 // a) Define pMax's and qMax's:
 Int_t pMax[5] = {2,3,4,5,8};      
 Int_t qMax[5] = {5,7,9,11,17};  
 
 // b) Book profile to hold all tuning parameters and flags:
 TString tuningFlagsName = "fTuningFlags";
 fTuningFlags = new TProfile(tuningFlagsName.Data(),"Tuning parameters",10,0,10);
 // fTuningFlags->SetTickLength(-0.01,"Y");
 fTuningFlags->SetMarkerStyle(25);
 fTuningFlags->SetLabelSize(0.05);
 fTuningFlags->SetLabelOffset(0.02,"X");
 for(Int_t r=1;r<=10;r++)
 {
  fTuningFlags->GetXaxis()->SetBinLabel(r,Form("r_{0,%d}",r-1));
  fTuningFlags->Fill(r-0.5,fTuningR0[r-1],1.);
 }
 fTuningList->Add(fTuningFlags);
  
 // c) Book all profiles:
 // Average of the generating function for reference flow <G[p][q]> for different tuning parameters:
 for(Int_t r=0;r<10;r++)
 { 
  for(Int_t pq=0;pq<5;pq++)
  {
   fTuningGenFun[r][pq] = new TProfile2D(Form("fTuningGenFun (r_{0,%i}, pq set %i)",r,pq),
                                         Form("#LTG[p][q]#GT for r_{0} = %f, p_{max} = %i, q_{max} = %i",fTuningR0[r],pMax[pq],qMax[pq]),
                                         pMax[pq],0.,(Double_t)pMax[pq],qMax[pq],0.,(Double_t)qMax[pq]);
   fTuningGenFun[r][pq]->SetXTitle("p");
   fTuningGenFun[r][pq]->SetYTitle("q");
   fTuningProfiles->Add(fTuningGenFun[r][pq]);
  }
 }
 // Average multiplicities for events with nRPs >= cuttof:
 fTuningAvM = new TProfile("fTuningAvM","Average multiplicity",5,0,5);
 fTuningAvM->SetMarkerStyle(25);
 for(Int_t b=1;b<=5;b++)
 {
  fTuningAvM->GetXaxis()->SetBinLabel(b,Form("nRP #geq %i",2*pMax[b-1]));
 }
 fTuningProfiles->Add(fTuningAvM); 

 // d) Book all histograms:
 // Final results for isotropic cumulants for reference flow for different tuning parameters:
 for(Int_t r=0;r<10;r++)
 { 
  for(Int_t pq=0;pq<5;pq++)
  {
   fTuningCumulants[r][pq] = new TH1D(Form("fTuningCumulants (r_{0,%i}, pq set %i)",r,pq),
                                      Form("GFC for r_{0} = %f, p_{max} = %i, q_{max} = %i",fTuningR0[r],pMax[pq],qMax[pq]),
                                      pMax[pq],0,pMax[pq]);
   // fTuningCumulants[r][pq]->SetLabelSize(0.05);
   fTuningCumulants[r][pq]->SetMarkerStyle(25);
   for(Int_t b=1;b<=pMax[pq];b++)
   {
    fTuningCumulants[r][pq]->GetXaxis()->SetBinLabel(b,Form("GFC{%i}",2*b));
   }
   fTuningResults->Add(fTuningCumulants[r][pq]); 
  }
 } 
 // Final results for reference flow for different tuning parameters: 
 for(Int_t r=0;r<10;r++)
 { 
  for(Int_t pq=0;pq<5;pq++)
  {
   fTuningFlow[r][pq] = new TH1D(Form("fTuningFlow (r_{0,%i}, pq set %i)",r,pq),
                                 Form("Reference flow for r_{0} = %f, p_{max} = %i, q_{max} = %i",fTuningR0[r],pMax[pq],qMax[pq]),
                                 pMax[pq],0,pMax[pq]);
   // fTuningFlow[r][pq]->SetLabelSize(0.06);
   fTuningFlow[r][pq]->SetMarkerStyle(25);
   for(Int_t b=1;b<=pMax[pq];b++)
   {
    fTuningFlow[r][pq]->GetXaxis()->SetBinLabel(b,Form("v{%i,GFC}",2*b));
   }
   fTuningResults->Add(fTuningFlow[r][pq]); 
  }
 }  
   
} // end of void AliFlowAnalysisWithCumulants::BookEverythingForTuning()

//================================================================================================================

void AliFlowAnalysisWithCumulants::BookEverythingForDiffFlow()
{
 // Book all objects relevant for calculation of differential flow.
 
 // a) Define static constants for array's boundaries;
 // b) Define local variables and local flags for booking;
 // c) Book profile to hold all flags for differential flow;
 // d) Book all event-by-event quantities;
 // e) Book all profiles;
 // f) Book all histograms.
 
 // a) Define static constants for array's boundaries:
 static const Int_t pMax = 5;     
 static const Int_t qMax = 11;  
 
 // b) Define local variables and local flags for booking:
 Int_t nBinsPtEta[2] = {fnBinsPt,fnBinsEta};
 Double_t minPtEta[2] = {fPtMin,fEtaMin};
 Double_t maxPtEta[2] = {fPtMax,fEtaMax};
 TString reIm[2] = {"Re","Im"};
 TString rpPoi[2] = {"RP","POI"};
 TString ptEta[2] = {"p_{t}","#eta"}; 
 TString order[4] = {"2nd order","4th order","6th order","8th order"}; 
 
 // c) Book profile to hold all flags for differential flow:
 TString diffFlowFlagsName = "fDiffFlowFlags";
 fDiffFlowFlags = new TProfile(diffFlowFlagsName.Data(),"Flags for Differential Flow",1,0,1);
 fDiffFlowFlags->SetTickLength(-0.01,"Y");
 fDiffFlowFlags->SetMarkerStyle(25);
 fDiffFlowFlags->SetLabelSize(0.05);
 fDiffFlowFlags->SetLabelOffset(0.02,"Y");
 fDiffFlowFlags->GetXaxis()->SetBinLabel(1,"...");
 fDiffFlowList->Add(fDiffFlowFlags);
 
 // d) Book all event-by-event quantities:
 // ... (to be improved - perhaps not needed)
 
 // e) Book all profiles:
 // Generating functions for differential flow:
 for(Int_t ri=0;ri<2;ri++)
 {
  for(Int_t rp=0;rp<2;rp++)
  {
   for(Int_t pe=0;pe<2;pe++)
   {
    fDiffFlowGenFun[ri][rp][pe] = new TProfile3D(Form("fDiffFlowGenFun (%s, %s, %s)",reIm[ri].Data(),rpPoi[rp].Data(),ptEta[pe].Data()),
                                                 Form("#LT%s[D[%s-bin][p][q]]#GT for %ss",reIm[ri].Data(),ptEta[pe].Data(),rpPoi[rp].Data()),
                                                 nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe],pMax,0.,(Double_t)pMax,qMax,0.,(Double_t)qMax);
    fDiffFlowGenFun[ri][rp][pe]->SetXTitle(ptEta[pe].Data());
    fDiffFlowGenFun[ri][rp][pe]->SetYTitle("p");
    fDiffFlowGenFun[ri][rp][pe]->SetZTitle("q");
    fDiffFlowGenFun[ri][rp][pe]->SetTitleOffset(1.44,"X");
    fDiffFlowGenFun[ri][rp][pe]->SetTitleOffset(1.44,"Y");
    fDiffFlowProfiles->Add(fDiffFlowGenFun[ri][rp][pe]);
    // to be improved - alternative // nBinsPtEta[pe],(Double_t)(fPtMin/fPtBinWidth),(Double_t)(fPtMax/fPtBinWidth),pMax,0.,(Double_t)pMax,qMax,0.,(Double_t)qMax);                                                 
   }
  }   
 } 
 // Number of particles in pt/eta bin for RPs/POIs:
 for(Int_t rp=0;rp<2;rp++)
 {
  for(Int_t pe=0;pe<2;pe++)
  {
   fNoOfParticlesInBin[rp][pe] = new TProfile(Form("fNoOfParticlesInBin (%s, %s)",rpPoi[rp].Data(),ptEta[pe].Data()),
                                              Form("Number of %ss per %s bin",rpPoi[rp].Data(),ptEta[pe].Data()),
                                              nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
   fNoOfParticlesInBin[rp][pe]->SetXTitle(ptEta[pe].Data());
   fDiffFlowProfiles->Add(fNoOfParticlesInBin[rp][pe]);
  }
 }   
 // Differential cumulants per pt/eta bin for RPs/POIs:
 for(Int_t rp=0;rp<2;rp++)
 {
  for(Int_t pe=0;pe<2;pe++)
  {
   for(Int_t co=0;co<4;co++)
   {
    fDiffFlowCumulants[rp][pe][co] = new TH1D(Form("fDiffFlowCumulants (%s, %s, %s)",rpPoi[rp].Data(),ptEta[pe].Data(),order[co].Data()),
                                              Form("Differential %s cumulant for %ss vs %s",order[co].Data(),rpPoi[rp].Data(),ptEta[pe].Data()),
                                              nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
    fDiffFlowCumulants[rp][pe][co]->SetXTitle(ptEta[pe].Data());
    fDiffFlowResults->Add(fDiffFlowCumulants[rp][pe][co]);
   }
  }
 }   
 // Differential flow per pt/eta bin for RPs/POIs:
 for(Int_t rp=0;rp<2;rp++)
 {
  for(Int_t pe=0;pe<2;pe++)
  {
   for(Int_t co=0;co<4;co++)
   {
    fDiffFlow[rp][pe][co] = new TH1D(Form("fDiffFlow (%s, %s, %s)",rpPoi[rp].Data(),ptEta[pe].Data(),order[co].Data()),
                                     Form("Differential flow from %s cumulant for %ss vs %s",order[co].Data(),rpPoi[rp].Data(),ptEta[pe].Data()),
                                     nBinsPtEta[pe],minPtEta[pe],maxPtEta[pe]);
    fDiffFlow[rp][pe][co]->SetXTitle(ptEta[pe].Data());
    fDiffFlowResults->Add(fDiffFlow[rp][pe][co]);
   }
  }
 }   

}// end of void AliFlowAnalysisWithCumulants::BookEverythingForDiffFlow()

//================================================================================================================

void AliFlowAnalysisWithCumulants::StoreReferenceFlowFlags()
{
 // Store all flags for reference flow in profile fReferenceFlowFlags.
 
 if(!fReferenceFlowFlags)
 {
  cout<<endl;
  cout<<"WARNING: !fReferenceFlowFlags is NULL in AFAWC::SRFF() !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 

 // Particle weights used or not:
 fReferenceFlowFlags->Fill(0.5,(Double_t)fUsePhiWeights||fUsePtWeights||fUseEtaWeights);
 // Which event weight was used to weight generating function event-by-event:
 if(strcmp(fMultiplicityWeight->Data(),"unit"))
 {
  fReferenceFlowFlags->Fill(1.5,0.); // 0 = "unit" (default)
 } else if(strcmp(fMultiplicityWeight->Data(),"multiplicity"))
   {
    fReferenceFlowFlags->Fill(1.5,1.); // 1 = "multiplicity"        
   } 
 fReferenceFlowFlags->Fill(2.5,fCalculateVsMultiplicity); // evaluate vs M?          

} // end of void AliFlowAnalysisWithCumulants::StoreReferenceFlowFlags()

//================================================================================================================

void AliFlowAnalysisWithCumulants::StoreDiffFlowFlags()
{
 // Store all flags for differential flow in profile fDiffFlowFlags.
 
 if(!fDiffFlowFlags)
 {
  cout<<endl;
  cout<<"WARNING: !fDiffFlowFlags is NULL in AFAWC::SRFF() !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 

 // fDiffFlags->Fill(0.5,(Double_t) ... );

} // end of void AliFlowAnalysisWithCumulants::StoreDiffFlowFlags()

//================================================================================================================

void AliFlowAnalysisWithCumulants::BookAndNestAllLists()
{
 // Book and nest all list in base list fHistList.
 
 // a) Book and nest lists for reference flow;
 // b) Book and nest lists for differential flow;
 // c) Book and nest lists for tuning;
 // d) If used, nest list for particle weights.   
  
 // a) Book and nest all lists for reference flow:
 fReferenceFlowList = new TList();
 fReferenceFlowList->SetName("Reference Flow");
 fReferenceFlowList->SetOwner(kTRUE);
 fHistList->Add(fReferenceFlowList);
 fReferenceFlowProfiles = new TList();
 fReferenceFlowProfiles->SetName("Profiles");
 fReferenceFlowProfiles->SetOwner(kTRUE);
 fReferenceFlowList->Add(fReferenceFlowProfiles);
 fReferenceFlowResults = new TList();
 fReferenceFlowResults->SetName("Results");
 fReferenceFlowResults->SetOwner(kTRUE);
 fReferenceFlowList->Add(fReferenceFlowResults);
 // b) Book and nest lists for differential flow:
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
 // c) Book and nest lists for tuning:
 if(fTuneParameters)
 {
  fTuningList = new TList();
  fTuningList->SetName("Tuning");
  fTuningList->SetOwner(kTRUE);
  fHistList->Add(fTuningList);
  fTuningProfiles = new TList();
  fTuningProfiles->SetName("Profiles");
  fTuningProfiles->SetOwner(kTRUE);
  fTuningList->Add(fTuningProfiles);
  fTuningResults = new TList();
  fTuningResults->SetName("Results");
  fTuningResults->SetOwner(kTRUE);
  fTuningList->Add(fTuningResults);
 } 
  
 // d) If used, nest list for particle weights.   
 if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights)
 {
  // Remark: pointer to this list is coming from the macro, no need to "new" it. 
  fWeightsList->SetName("Weights");
  fWeightsList->SetOwner(kTRUE);
  fHistList->Add(fWeightsList); 
 }

} // end of void AliFlowAnalysisWithCumulants::BookAndNestAllLists()

//================================================================================================================

void AliFlowAnalysisWithCumulants::BookProfileHoldingSettings()
{
 // Book profile to hold all analysis settings.

 TString analysisSettingsName = "fAnalysisSettings";
 fAnalysisSettings = new TProfile(analysisSettingsName.Data(),"Settings for analysis with Generating Function Cumulants",11,0.,11.);
 fAnalysisSettings->GetXaxis()->SetLabelSize(0.035);
 fAnalysisSettings->GetXaxis()->SetBinLabel(1,"Harmonic");
 fAnalysisSettings->Fill(0.5,fHarmonic);
 fAnalysisSettings->GetXaxis()->SetBinLabel(2,"Multiple");
 fAnalysisSettings->Fill(1.5,fMultiple); 
 fAnalysisSettings->GetXaxis()->SetBinLabel(3,"r_{0}");
 fAnalysisSettings->Fill(2.5,fR0);   
 fAnalysisSettings->GetXaxis()->SetBinLabel(4,"Use w_{#phi}?");
 fAnalysisSettings->Fill(3.5,fUsePhiWeights);
 fAnalysisSettings->GetXaxis()->SetBinLabel(5,"Use w_{p_{t}}?");
 fAnalysisSettings->Fill(4.5,fUsePtWeights);
 fAnalysisSettings->GetXaxis()->SetBinLabel(6,"Use w_{#eta}?");
 fAnalysisSettings->Fill(5.5,fUsePhiWeights);
 fAnalysisSettings->GetXaxis()->SetBinLabel(7,"Tune parameters?");
 fAnalysisSettings->Fill(6.5,fTuneParameters); 
 fAnalysisSettings->GetXaxis()->SetBinLabel(8,"Print RF results");
 fAnalysisSettings->Fill(7.5,fPrintFinalResults[0]); 
 fAnalysisSettings->GetXaxis()->SetBinLabel(9,"Print RP results");
 fAnalysisSettings->Fill(8.5,fPrintFinalResults[1]); 
 fAnalysisSettings->GetXaxis()->SetBinLabel(10,"Print POI results");
 fAnalysisSettings->Fill(9.5,fPrintFinalResults[2]);
 fAnalysisSettings->GetXaxis()->SetBinLabel(11,"Evaluate vs M?");
 fAnalysisSettings->Fill(10.5,fCalculateVsMultiplicity);
 fHistList->Add(fAnalysisSettings);

} // end of void AliFlowAnalysisWithCumulants::BookProfileHoldingSettings()

//================================================================================================================

void AliFlowAnalysisWithCumulants::BookCommonHistograms()
{
 // Book common control histograms and common histograms for final results.
 
 // Common control histogram:
 TString commonHistsName = "AliFlowCommonHistGFC";
 fCommonHists = new AliFlowCommonHist(commonHistsName.Data());
 fHistList->Add(fCommonHists);  
 // Common histograms for final results from 2nd order GFC:
 TString commonHistResults2ndOrderName = "AliFlowCommonHistResults2ndOrderGFC";
 fCommonHistsResults2nd = new AliFlowCommonHistResults(commonHistResults2ndOrderName.Data(),"",fHarmonic);
 fHistList->Add(fCommonHistsResults2nd);  
 // Common histograms for final results from 4th order GFC:
 TString commonHistResults4thOrderName = "AliFlowCommonHistResults4thOrderGFC";
 fCommonHistsResults4th = new AliFlowCommonHistResults(commonHistResults4thOrderName.Data(),"",fHarmonic);
 fHistList->Add(fCommonHistsResults4th); 
 // Common histograms for final results from 6th order GFC:
 TString commonHistResults6thOrderName = "AliFlowCommonHistResults6thOrderGFC";
 fCommonHistsResults6th = new AliFlowCommonHistResults(commonHistResults6thOrderName.Data(),"",fHarmonic);
 fHistList->Add(fCommonHistsResults6th);  
 // Common histograms for final results from 8th order GFC:
 TString commonHistResults8thOrderName = "AliFlowCommonHistResults8thOrderGFC";
 fCommonHistsResults8th = new AliFlowCommonHistResults(commonHistResults8thOrderName.Data(),"",fHarmonic);
 fHistList->Add(fCommonHistsResults8th);
 
} // end of void AliFlowAnalysisWithCumulants::BookCommonHistograms()

//================================================================================================================

void AliFlowAnalysisWithCumulants::CheckPointersUsedInMake()
{
 // Check pointers used in method Make().
 
 if(!fCommonHists)
 {                        
  cout<<endl;
  cout<<" WARNING (GFC): fCommonHists is NULL in CPUIM() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(fUsePhiWeights && !fPhiWeights)
 {                        
  cout<<endl;
  cout<<" WARNING (GFC): fPhiWeights is NULL in CPUIM() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(fUsePtWeights && !fPtWeights)
 {                        
  cout<<endl;
  cout<<" WARNING (GFC): fPtWeights is NULL in CPUIM() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(fUseEtaWeights && !fEtaWeights)
 {                        
  cout<<endl;
  cout<<" WARNING (GFC): fEtaWeights is NULL in CPUIM() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fAverageOfSquaredWeight)
 {                        
  cout<<endl;
  cout<<" WARNING (GFC): fAverageOfSquaredWeight is NULL in CPUIM() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fReferenceFlowGenFun)
 {                        
  cout<<endl;
  cout<<" WARNING (GFC): fReferenceFlowGenFun is NULL in CPUIM() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fQvectorComponents)
 {                        
  cout<<endl;
  cout<<" WARNING (GFC): fQvectorComponents is NULL in CPUIM() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!fGEBE)
 {
  cout<<endl; 
  cout<<"WARNING (GFC): fGEBE is NULL in CPUIM() !!!!"<<endl;
  cout<<endl; 
  exit(0);
 }
 // Checking pointers for vs multiplicity calculation: 
 if(fCalculateVsMultiplicity)
 {
  if(!fReferenceFlowGenFunVsM)
  {
   cout<<endl; 
   cout<<"WARNING (GFC): fReferenceFlowGenFunVsM is NULL in CPUIM() !!!!"<<endl;
   cout<<endl; 
   exit(0);
  }
  if(!fQvectorComponentsVsM)
  {
   cout<<endl; 
   cout<<"WARNING (GFC): fQvectorComponentsVsM is NULL in CPUIM() !!!!"<<endl;
   cout<<endl; 
   exit(0);
  }
  if(!fAverageOfSquaredWeightVsM)
  {
   cout<<endl; 
   cout<<"WARNING (GFC): fAverageOfSquaredWeightVsM is NULL in CPUIM() !!!!"<<endl;
   cout<<endl; 
   exit(0);
  }
  if(!fAvMVsM)
  {
   cout<<endl; 
   cout<<"WARNING (GFC): fAvMVsM is NULL in CPUIM() !!!!"<<endl;
   cout<<endl; 
   exit(0);
  }
 } // end of if(fCalculateVsMultiplicity) 
 
} // end of void AliFlowAnalysisWithCumulants::CheckPointersUsedInMake()

//================================================================================================================

void AliFlowAnalysisWithCumulants::CheckPointersUsedInFinish()
{
 // Check pointers used in method Finish().
 
 if(!fAnalysisSettings)
 {                        
  cout<<endl;
  cout<<" WARNING (GFC): fAnalysisSettings is NULL in CPUIF() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!(fCommonHists && fCommonHists->GetHistMultRP()))
 { 
  cout<<endl;
  cout<<" WARNING (GFC): (fCommonHists && fCommonHists->GetHistMultRP) is NULL in CPUIF() !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 
 if(!fReferenceFlowGenFun)
 { 
  cout<<endl;
  cout<<" WARNING (GFC): fReferenceFlowGenFun is NULL in CPUIF() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }  
 if(!fReferenceFlowCumulants)
 { 
  cout<<endl;
  cout<<" WARNING (GFC): fReferenceFlowCumulants is NULL in CPUIF() !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 
 if(!fQvectorComponents)
 { 
  cout<<endl;
  cout<<" WARNING (GFC): fQvectorComponents is NULL in CPUIF() !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 
 if(!fAverageOfSquaredWeight)
 { 
  cout<<endl;
  cout<<" WARNING (GFC): fAverageOfSquaredWeight is NULL in CPUIF() !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 
 if(!(fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && fCommonHistsResults8th))
 {
  cout<<endl;
  cout<<" WARNING (GFC): fCommonHistsResults2nd && fCommonHistsResults4th && fCommonHistsResults6th && "<<endl;
  cout<<"                fCommonHistsResults8th is NULL in CPUIF() !!!!"<<endl; 
  cout<<endl;
  exit(0);
 }
 if(!fReferenceFlow)
 { 
  cout<<endl;
  cout<<" WARNING (GFC): fReferenceFlow is NULL in CPUIF() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }  
 if(!fChi)
 { 
  cout<<endl;
  cout<<" WARNING (GFC): fChi is NULL in CPUIF() !!!!"<<endl;
  cout<<endl;
  exit(0);
 }  
 for(Int_t ri=0;ri<2;ri++)
 {
  for(Int_t rp=0;rp<2;rp++)
  {
   for(Int_t pe=0;pe<2;pe++)
   {
    if(!fDiffFlowGenFun[ri][rp][pe])
    {
     cout<<endl;
     cout<<" WARNING (GFC): "<<Form("fDiffFlowGenFun[%d][%d][%d]",ri,rp,pe)<<" is NULL in CPUIF() !!!!"<<endl;
     cout<<endl;
     exit(0);    
    }
   }
  }
 }
 for(Int_t rp=0;rp<2;rp++)
 {
  for(Int_t pe=0;pe<2;pe++)
  {
   for(Int_t co=0;co<4;co++)
   {
    if(!fDiffFlowCumulants[rp][pe][co])
    {
     cout<<endl;
     cout<<" WARNING (GFC): "<<Form("fDiffFlowCumulants[%d][%d][%d]",rp,pe,co)<<" is NULL in CPUIF() !!!!"<<endl;
     cout<<endl;
     exit(0);    
    }
    if(!fDiffFlow[rp][pe][co])
    {
     cout<<endl;
     cout<<" WARNING (GFC): "<<Form("fDiffFlow[%d][%d][%d]",rp,pe,co)<<" is NULL in CPUIF() !!!!"<<endl;
     cout<<endl;
     exit(0);    
    }
   }
  }
 }
 for(Int_t rp=0;rp<2;rp++)
 {
  for(Int_t pe=0;pe<2;pe++)
  {
   if(!fNoOfParticlesInBin[rp][pe])
   {
    cout<<endl;
    cout<<" WARNING (GFC): "<<Form("fNoOfParticlesInBin[%d][%d]",rp,pe)<<" is NULL in CPUIF() !!!!"<<endl;
    cout<<endl;
    exit(0);    
   } 
  }
 }  
 // Checking pointers for vs multiplicity calculation: 
 if(fCalculateVsMultiplicity)
 {
  if(!fReferenceFlowGenFunVsM)
  {
   cout<<endl; 
   cout<<"WARNING (GFC): fReferenceFlowGenFunVsM is NULL in CPUIF() !!!!"<<endl;
   cout<<endl; 
   exit(0);
  }
  if(!fQvectorComponentsVsM)
  {
   cout<<endl; 
   cout<<"WARNING (GFC): fQvectorComponentsVsM is NULL in CPUIF() !!!!"<<endl;
   cout<<endl; 
   exit(0);
  }
  if(!fAverageOfSquaredWeightVsM)
  {
   cout<<endl; 
   cout<<"WARNING (GFC): fAverageOfSquaredWeightVsM is NULL in CPUIF() !!!!"<<endl;
   cout<<endl; 
   exit(0);
  }
  if(!fAvMVsM)
  {
   cout<<endl; 
   cout<<"WARNING (GFC): fAvMVsM is NULL in CPUIF() !!!!"<<endl;
   cout<<endl; 
   exit(0);
  }
 } // end of if(fCalculateVsMultiplicity) 
 
} // end of void AliFlowAnalysisWithCumulants::CheckPointersUsedInFinish()

//================================================================================================================

void AliFlowAnalysisWithCumulants::AccessSettings()
{
 // Access the settings for analysis with Generating Function Cumulants.
 
 fHarmonic = (Int_t)fAnalysisSettings->GetBinContent(1); 
 fMultiple = (Int_t)fAnalysisSettings->GetBinContent(2); 
 fR0 = (Double_t)fAnalysisSettings->GetBinContent(3); 
 fUsePhiWeights = (Bool_t)fAnalysisSettings->GetBinContent(4);
 fUsePtWeights = (Bool_t)fAnalysisSettings->GetBinContent(5);
 fUseEtaWeights = (Bool_t)fAnalysisSettings->GetBinContent(6);
 fTuneParameters = (Bool_t)fAnalysisSettings->GetBinContent(7);
 fPrintFinalResults[0] = (Bool_t)fAnalysisSettings->GetBinContent(8);
 fPrintFinalResults[1] = (Bool_t)fAnalysisSettings->GetBinContent(9);
 fPrintFinalResults[2] = (Bool_t)fAnalysisSettings->GetBinContent(10);
 fCalculateVsMultiplicity = (Bool_t)fAnalysisSettings->GetBinContent(11);
 
} // end of AliFlowAnalysisWithCumulants::AccessSettings()

//================================================================================================================

void AliFlowAnalysisWithCumulants::WriteHistograms(TString* outputFileName)
{
 // Store the final results in output .root file.
 
 TFile *output = new TFile(outputFileName->Data(),"RECREATE");
 //output->WriteObject(fHistList, "cobjGFC","SingleKey");
 fHistList->SetName("cobjGFC");
 fHistList->SetOwner(kTRUE);
 fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
 delete output;

} // end of void AliFlowAnalysisWithCumulants::WriteHistograms(TString* outputFileName)

//================================================================================================================

void AliFlowAnalysisWithCumulants::WriteHistograms(TString outputFileName)
{
 // Store the final results in output .root file.
 
 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 //output->WriteObject(fHistList, "cobjGFC","SingleKey");
 fHistList->SetName("cobjGFC");
 fHistList->SetOwner(kTRUE);
 fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
 delete output;

} // end of void AliFlowAnalysisWithCumulants::WriteHistograms(TString outputFileName)

//================================================================================================================

void AliFlowAnalysisWithCumulants::WriteHistograms(TDirectoryFile *outputFileName)
{
 // Store the final results in output .root file.
 
 fHistList->SetName("cobjGFC");
 fHistList->SetOwner(kTRUE);
 outputFileName->Add(fHistList);
 outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);

} // end of void AliFlowAnalysisWithCumulants::WriteHistograms(TDirectoryFile *outputFileName)

//================================================================================================================




