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

 /************************************ 
 * flow analysis with multi-particle *
 *           correlations            * 
 *                                   * 
 * author: Ante Bilandzic            * 
 *        (abilandzic@gmail.com)     *
 ************************************/ 

#define AliFlowAnalysisWithMultiparticleCorrelations_cxx

#include "AliFlowAnalysisWithMultiparticleCorrelations.h"

using std::endl;
using std::cout;
using std::flush;

//================================================================================================================

ClassImp(AliFlowAnalysisWithMultiparticleCorrelations)

AliFlowAnalysisWithMultiparticleCorrelations::AliFlowAnalysisWithMultiparticleCorrelations(): 
 // 0.) Base:
 fHistList(NULL),
 fInternalFlagsPro(NULL),
 fUseInternalFlags(kFALSE),
 fMinNoRPs(-44),
 fMaxNoRPs(-44),
 fExactNoRPs(-44),
 // 1.) Control histograms:
 fControlHistogramsList(NULL),
 fControlHistogramsFlagsPro(NULL),
 fFillControlHistograms(kTRUE),
 // 2.) Q-vector:
 // ...
 // 3.) Correlations:
 fCorrelationsList(NULL),
 fCorrelationsFlagsPro(NULL),
 fCalculateCorrelations(kTRUE),
 fMaxHarmonic(6),
 fMaxCorrelator(9),
 // 4.) Cumulants:
 fCumulantsList(NULL),
 fCumulantsFlagsPro(NULL),
 f2pCumulantsPro(NULL),
 fCalculateCumulants(kTRUE),
 // 5.) Weights:
 fWeightsList(NULL),
 fWeightsFlagsPro(NULL),
 fUsePhiWeights(kFALSE),  
 fUsePtWeights(kFALSE),
 fUseEtaWeights(kFALSE),
 fPhiWeightsHist(NULL),
 fPtWeightsHist(NULL),
 fEtaWeightsHist(NULL),
 // 6.) Nested loops:
 fNestedLoopsList(NULL),
 fNestedLoopsFlagsPro(NULL),
 fCrossCheckWithNestedLoops(kFALSE),
 fNestedLoopsResultsCosPro(NULL),
 fNestedLoopsResultsSinPro(NULL),
 // 7.) 'Standard candles':
 fStandardCandlesList(NULL),
 fStandardCandlesFlagsPro(NULL),
 fCalculateStandardCandles(kFALSE),
 fStandardCandlesHist(NULL)
 {
  // Constructor.  
  
  // a) Book grandmother of all lists;
  // b) Initialize all arrays.

  // a) Book grandmother of all lists:
  fHistList = new TList();
  fHistList->SetName("cobjMPC");
  fHistList->SetOwner(kTRUE);

  // b) Initialize all arrays:
  this->InitializeArraysForControlHistograms();
  this->InitializeArraysForQvector();
  this->InitializeArraysForCorrelations();

 } // end of AliFlowAnalysisWithMultiparticleCorrelations::AliFlowAnalysisWithMultiparticleCorrelations()
 
//================================================================================================================  

AliFlowAnalysisWithMultiparticleCorrelations::~AliFlowAnalysisWithMultiparticleCorrelations()
{
 // Destructor.
 
 delete fHistList;

} // end of AliFlowAnalysisWithMultiparticleCorrelations::~AliFlowAnalysisWithMultiparticleCorrelations()

//================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::Init()
{
 // Well, this is method Init().

 // a) Trick to avoid name clashes, part 1; 
 // b) Cross-check the initial settings before starting this adventure;
 // c) Book all objects;
 // d) Set all flags;
 // *) Trick to avoid name clashes, part 2. 

 // a) Trick to avoid name clashes, part 1: 
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 TH1::AddDirectory(kFALSE);
 
 // b) Cross-check the initial settings before starting this adventure:
 this->CrossCheckSettings();

 // c) Book all objects:
 this->BookAndNestAllLists();
 this->BookEverythingForBase();
 this->BookEverythingForControlHistograms();
 this->BookEverythingForWeights();
 this->BookEverythingForCorrelations();
 this->BookEverythingForCumulants();
 this->BookEverythingForNestedLoops();
 this->BookEverythingForStandardCandles();

 // d) Set all flags:
 // ... 

 // *) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::Init()

//================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::Make(AliFlowEventSimple *anEvent)
{
 // Running over data only in this method.
 
 // a) Cross-check internal flags;
 // b) Fill control histograms;
 // c) Fill Q-vector components;
 // d) Calculate multi-particle correlations from Q-vector components; 
 // e) Calculate e-b-e cumulants; 
 // f) Reset Q-vector components;
 // g) Cross-check results with nested loops.

 // a) Cross-check internal flags:
 if(fUseInternalFlags){if(!this->CrossCheckInternalFlags(anEvent)){return;}}

 // b) Fill control histograms:
 if(fFillControlHistograms){this->FillControlHistograms(anEvent);}
 
 // c) Fill Q-vector components: // TBI add setter to disable it when needed
 this->FillQvector(anEvent);

 // d) Calculate multi-particle correlations from Q-vector components: 
 if(fCalculateCorrelations){this->CalculateCorrelations();}

 // e) Calculate e-b-e cumulants: 
 if(fCalculateCumulants){this->CalculateCumulants();}

 // f) Reset Q-vector components: // TBI add setter to disable it when needed
 this->ResetQvector();

 // g) Cross-check results with nested loops:
 if(fCrossCheckWithNestedLoops){this->CrossCheckWithNestedLoops(anEvent);}

} // end of AliFlowAnalysisWithMultiparticleCorrelations::Make(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::Finish()
{
 // Closing the curtains. 

 // a) Cross-check TBI;
 // b) Calculate 'standard candles'.

 // a) Cross-check TBI:
 // ...

 // b) Calculate 'standard candles':
 if(fCalculateStandardCandles){this->CalculateStandardCandles();}


 // ...
 printf("\n ... Closing the curtains ... \n");

} // end of AliFlowAnalysisWithMultiparticleCorrelations::Finish()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateStandardCandles()
{
 // Calculate and store 'standard candles'.

 // ...

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateStandardCandles()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForCorrelations()
{
 // Initialize all arrays for correlations.

 for(Int_t c=0;c<8;c++)
 {
  fCorrelationsPro[c] = NULL;
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateCorrelations()
{
 // Calculate multi-particle correlations from Q-vector components.

 // TBI mess

 // One-particle stuff;
 // Two-particle correlations;
 // 3-p; TBI
 // 4-p; TBI
 // ...

 // One-particle stuff:
 for(Int_t ci=1;ci<=fMaxHarmonic;ci++)
 {
  Double_t oneNum = One(ci).Re(); // numerator
  Double_t oneDen = One(0).Re(); // denominator
  if(TMath::Abs(oneDen)>0) 
  {
   Double_t oneWeight = oneDen; // TBI add choices for the weights
   fCorrelationsPro[0]->Fill(ci-0.5,oneNum/oneDen,oneWeight); 
  } // if(TMath::Abs(oneDen)>0) 
 } // for(Int_t ci=1;ci<=fMaxHarmonic;ci++)
 
 // Two-particle correlations:
 for(Int_t ci=1;ci<=fMaxHarmonic;ci++)
 {
  Double_t twoNum = Two(ci,-ci).Re(); // numerator
  Double_t twoDen = Two(0,0).Re(); // denominator
  if(TMath::Abs(twoDen)>0) 
  {
   Double_t twoWeight = twoDen; // TBI add choices for the weights
   fCorrelationsPro[1]->Fill(ci-0.5,twoNum/twoDen,twoWeight); 
  } // if(TMath::Abs(twoDen)>0)
 } // for(Int_t ci=1;ci<=fMaxHarmonic;ci++)

 // 3-p:
 // ...

 // 4-p:
 Double_t four3232Num = Four(3,2,-3,-2).Re(); // TBI hardwired stuff
 Double_t four4242Num = Four(4,2,-4,-2).Re(); // TBI hardwired stuff
 Double_t four4343Num = Four(4,3,-4,-3).Re(); // TBI hardwired stuff
 Double_t four5353Num = Four(5,3,-5,-3).Re(); // TBI hardwired stuff
 Double_t fourDen = Four(0,0,0,0).Re(); // TBI
 if(TMath::Abs(fourDen)>0) // protection against zero TBI
 {
  fCorrelationsPro[3]->Fill(0.5,four3232Num/fourDen,fourDen); // TBI 
  fCorrelationsPro[3]->Fill(1.5,four4242Num/fourDen,fourDen); // TBI 
  fCorrelationsPro[3]->Fill(2.5,four4343Num/fourDen,fourDen); // TBI 
  fCorrelationsPro[3]->Fill(3.5,four5353Num/fourDen,fourDen); // TBI 
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateCumulants()
{
 // Calculate e-b-e cumulants from Q-vector components.

 // 2-p.

 // 2-p:
 for(Int_t ci=1;ci<=6;ci++)
 {
  Double_t twoNum = Two(ci,-ci).Re(); // numerator TBI
  Double_t twoDen = Two(0,0).Re(); // denominator TBI

  TComplex oneFirstNum = One(ci); // numerator TBI
  Double_t oneFirstDen = One(0).Re(); // denominator TBI

  TComplex oneSecondNum = One(-ci); // numerator TBI
  Double_t oneSecondDen = One(0).Re(); // denominator TBI

  if(TMath::Abs(twoDen)>0 && TMath::Abs(oneFirstDen*oneSecondDen) > 0)  // protection against zero TBI
  {
   Double_t cumulant = twoNum/twoDen - (oneFirstNum*oneSecondNum).Re()/(oneFirstDen*oneSecondDen);
   f2pCumulantsPro->Fill(ci-0.5,cumulant,1.); // TBI add choices for the weights
  }
 } // for(Int_t ci=1;ci<=6;ci++)

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateCumulants()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckWithNestedLoops(AliFlowEventSimple *anEvent)
{
 // Cross-check results for multi-particle correlations with nested loops.

 // TBI add few comments here, there and over there
 // TBI this method is rather messy :'(

 Int_t h1 = fNestedLoopsFlagsPro->GetBinContent(2);
 Int_t h2 = fNestedLoopsFlagsPro->GetBinContent(3);
 Int_t h3 = fNestedLoopsFlagsPro->GetBinContent(4);
 Int_t h4 = fNestedLoopsFlagsPro->GetBinContent(5);
 Int_t h5 = fNestedLoopsFlagsPro->GetBinContent(6);
 Int_t h6 = fNestedLoopsFlagsPro->GetBinContent(7);
 Int_t h7 = fNestedLoopsFlagsPro->GetBinContent(8);
 Int_t h8 = fNestedLoopsFlagsPro->GetBinContent(9);

 this->ResetQvector();
 this->FillQvector(anEvent);

 fNestedLoopsResultsCosPro->Fill(1.5,One(h1).Re()/One(0).Re(),One(0).Re()); 
 fNestedLoopsResultsSinPro->Fill(1.5,One(h1).Im()/One(0).Re(),One(0).Re());  
 fNestedLoopsResultsCosPro->Fill(3.5,Two(h1,h2).Re()/Two(0,0).Re(),Two(0,0).Re()); 
 fNestedLoopsResultsSinPro->Fill(3.5,Two(h1,h2).Im()/Two(0,0).Re(),Two(0,0).Re()); 
 fNestedLoopsResultsCosPro->Fill(5.5,Three(h1,h2,h3).Re()/Three(0,0,0).Re(),Three(0,0,0).Re()); 
 fNestedLoopsResultsSinPro->Fill(5.5,Three(h1,h2,h3).Im()/Three(0,0,0).Re(),Three(0,0,0).Re()); 
 fNestedLoopsResultsCosPro->Fill(7.5,Four(h1,h2,h3,h4).Re()/Four(0,0,0,0).Re(),Four(0,0,0,0).Re()); 
 fNestedLoopsResultsSinPro->Fill(7.5,Four(h1,h2,h3,h4).Im()/Four(0,0,0,0).Re(),Four(0,0,0,0).Re()); 
 fNestedLoopsResultsCosPro->Fill(9.5,Five(h1,h2,h3,h4,h5).Re()/Five(0,0,0,0,0).Re(),Five(0,0,0,0,0).Re()); 
 fNestedLoopsResultsSinPro->Fill(9.5,Five(h1,h2,h3,h4,h5).Im()/Five(0,0,0,0,0).Re(),Five(0,0,0,0,0).Re()); 
 fNestedLoopsResultsCosPro->Fill(11.5,Six(h1,h2,h3,h4,h5,h6).Re()/Six(0,0,0,0,0,0).Re(),Six(0,0,0,0,0,0).Re()); 
 fNestedLoopsResultsSinPro->Fill(11.5,Six(h1,h2,h3,h4,h5,h6).Im()/Six(0,0,0,0,0,0).Re(),Six(0,0,0,0,0,0).Re()); 
 //fNestedLoopsResultsCosPro->Fill(13.5,Seven(h1,h2,h3,h4,h5,h6,h7).Re()/Seven(0,0,0,0,0,0,0).Re(),Seven(0,0,0,0,0,0,0).Re()); 
 //fNestedLoopsResultsSinPro->Fill(13.5,Seven(h1,h2,h3,h4,h5,h6,h7).Im()/Seven(0,0,0,0,0,0,0).Re(),Seven(0,0,0,0,0,0,0).Re()); 
 //fNestedLoopsResultsCosPro->Fill(15.5,Eight(h1,h2,h3,h4,h5,h6,h7,h8).Re()/Eight(0,0,0,0,0,0,0,0).Re(),Eight(0,0,0,0,0,0,0,0).Re()); 
 //fNestedLoopsResultsSinPro->Fill(15.5,Eight(h1,h2,h3,h4,h5,h6,h7,h8).Im()/Eight(0,0,0,0,0,0,0,0).Re(),Eight(0,0,0,0,0,0,0,0).Re()); 

 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL; 
 Double_t dPhi1=0.,dPhi2=0.,dPhi3=0.,dPhi4=0.,dPhi5=0.,dPhi6=0.,dPhi7=0.,dPhi8=0.; 
 Double_t wPhi1=1.,wPhi2=1.,wPhi3=1.,wPhi4=1.,wPhi5=1.,wPhi6=1.,wPhi7=1.,wPhi8=1.; 

 // 1-particle stuff: TBI       
 if(nPrim>=1)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack = anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   dPhi1 = aftsTrack->Phi(); 
   if(fUsePhiWeights){wPhi1 = PhiWeight(dPhi1);}
   // Fill:
   fNestedLoopsResultsCosPro->Fill(0.5,TMath::Cos(h1*dPhi1),wPhi1); 
   fNestedLoopsResultsSinPro->Fill(0.5,TMath::Sin(h1*dPhi1),wPhi1); 
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=1) 

 // 2-particle correlations:       
 if(nPrim>=2)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack = anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   dPhi1 = aftsTrack->Phi(); 
   if(fUsePhiWeights){wPhi1 = PhiWeight(dPhi1);}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack = anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    dPhi2 = aftsTrack->Phi();
    if(fUsePhiWeights){wPhi2 = PhiWeight(dPhi2);}
    // Fill:
    fNestedLoopsResultsCosPro->Fill(2.5,TMath::Cos(h1*dPhi1+h2*dPhi2),wPhi1*wPhi2); 
    fNestedLoopsResultsSinPro->Fill(2.5,TMath::Sin(h1*dPhi1+h2*dPhi2),wPhi1*wPhi2); 
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=2)

 // 3-particle correlations:         
 if(nPrim>=3)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   dPhi1=aftsTrack->Phi();
   if(fUsePhiWeights){wPhi1 = PhiWeight(dPhi1);}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    dPhi2=aftsTrack->Phi();
    if(fUsePhiWeights){wPhi2 = PhiWeight(dPhi2);}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     dPhi3=aftsTrack->Phi();
     if(fUsePhiWeights){wPhi3 = PhiWeight(dPhi3);}
     // Fill:
     fNestedLoopsResultsCosPro->Fill(4.5,TMath::Cos(h1*dPhi1+h2*dPhi2+h3*dPhi3),wPhi1*wPhi2*wPhi3);
     fNestedLoopsResultsSinPro->Fill(4.5,TMath::Sin(h1*dPhi1+h2*dPhi2+h3*dPhi3),wPhi1*wPhi2*wPhi3);
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=3)

 // 4-particle correlations:
 if(nPrim>=4)
 {       
  for(Int_t i1=0;i1<nPrim;i1++)
  { 
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   dPhi1=aftsTrack->Phi();
   if(fUsePhiWeights){wPhi1 = PhiWeight(dPhi1);}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    dPhi2=aftsTrack->Phi();
    if(fUsePhiWeights){wPhi2 = PhiWeight(dPhi2);}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     dPhi3=aftsTrack->Phi();
     if(fUsePhiWeights){wPhi3 = PhiWeight(dPhi3);}
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3){continue;}
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())){continue;}
      dPhi4=aftsTrack->Phi();
      if(fUsePhiWeights){wPhi4 = PhiWeight(dPhi4);}
      // Fill:
      fNestedLoopsResultsCosPro->Fill(6.5,TMath::Cos(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4),wPhi1*wPhi2*wPhi3*wPhi4);
      fNestedLoopsResultsSinPro->Fill(6.5,TMath::Sin(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4),wPhi1*wPhi2*wPhi3*wPhi4);
     } // end of for(Int_t i4=0;i4<nPrim;i4++) 
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=)

 // 5-particle correlations:      
 if(nPrim>=5)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}  
   dPhi1=aftsTrack->Phi();
   if(fUsePhiWeights){wPhi1 = PhiWeight(dPhi1);}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    dPhi2=aftsTrack->Phi();
    if(fUsePhiWeights){wPhi2 = PhiWeight(dPhi2);}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     dPhi3=aftsTrack->Phi();
     if(fUsePhiWeights){wPhi3 = PhiWeight(dPhi3);}
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3){continue;}
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())){continue;}
      dPhi4=aftsTrack->Phi();
      if(fUsePhiWeights){wPhi4 = PhiWeight(dPhi4);}
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())){continue;}
       dPhi5=aftsTrack->Phi();
       if(fUsePhiWeights){wPhi5 = PhiWeight(dPhi5);}
       // Fill:   
       fNestedLoopsResultsCosPro->Fill(8.5,TMath::Cos(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4+h5*dPhi5),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5);
       fNestedLoopsResultsSinPro->Fill(8.5,TMath::Sin(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4+h5*dPhi5),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5);
      } // end of for(Int_t i5=0;i5<nPrim;i5++)
     } // end of for(Int_t i4=0;i4<nPrim;i4++)  
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=5)
  
 // 6-particle correlations:
 if(nPrim>=6)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   dPhi1=aftsTrack->Phi();
   if(fUsePhiWeights){wPhi1 = PhiWeight(dPhi1);}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    dPhi2=aftsTrack->Phi();
    if(fUsePhiWeights){wPhi2 = PhiWeight(dPhi2);}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     dPhi3=aftsTrack->Phi();
     if(fUsePhiWeights){wPhi3 = PhiWeight(dPhi3);}
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3){continue;}
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())){continue;}
      dPhi4=aftsTrack->Phi();
      if(fUsePhiWeights){wPhi4 = PhiWeight(dPhi4);}
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())){continue;}
       dPhi5=aftsTrack->Phi();
       if(fUsePhiWeights){wPhi5=PhiWeight(dPhi5);}
       for(Int_t i6=0;i6<nPrim;i6++)
       {
        if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
        aftsTrack=anEvent->GetTrack(i6);
        if(!(aftsTrack->InRPSelection())){continue;}
        dPhi6=aftsTrack->Phi(); 
        if(fUsePhiWeights){wPhi6=PhiWeight(dPhi6);}
        // Fill:   
        fNestedLoopsResultsCosPro->Fill(10.5,TMath::Cos(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4+h5*dPhi5+h6*dPhi6),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6);
        fNestedLoopsResultsSinPro->Fill(10.5,TMath::Sin(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4+h5*dPhi5+h6*dPhi6),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6);
       } // end of for(Int_t i6=0;i6<nPrim;i6++)
      } // end of for(Int_t i5=0;i5<nPrim;i5++)
     } // end of for(Int_t i4=0;i4<nPrim;i4++)
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=6)
  
 // 7-particle correlations:
 if(nPrim>=7)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  { 
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   dPhi1=aftsTrack->Phi();
   if(fUsePhiWeights){wPhi1=PhiWeight(dPhi1);}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    dPhi2=aftsTrack->Phi();
    if(fUsePhiWeights){wPhi2=PhiWeight(dPhi2);}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     dPhi3=aftsTrack->Phi();
     if(fUsePhiWeights){wPhi3=PhiWeight(dPhi3);}
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3){continue;}
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())){continue;}
      dPhi4=aftsTrack->Phi();
      if(fUsePhiWeights){wPhi4=PhiWeight(dPhi4);}
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())){continue;}
       dPhi5=aftsTrack->Phi();
       if(fUsePhiWeights){wPhi5=PhiWeight(dPhi5);}
       for(Int_t i6=0;i6<nPrim;i6++)
       {
        if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
        aftsTrack=anEvent->GetTrack(i6);
        if(!(aftsTrack->InRPSelection())){continue;}
        dPhi6=aftsTrack->Phi(); 
        if(fUsePhiWeights){wPhi6=PhiWeight(dPhi6);}
        for(Int_t i7=0;i7<nPrim;i7++)
        {
         if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6){continue;}
         aftsTrack=anEvent->GetTrack(i7);
         if(!(aftsTrack->InRPSelection())){continue;}
         dPhi7=aftsTrack->Phi(); 
         if(fUsePhiWeights){wPhi7=PhiWeight(dPhi7);}
         // Fill:   
         fNestedLoopsResultsCosPro->Fill(12.5,TMath::Cos(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4+h5*dPhi5+h6*dPhi6+h7*dPhi7),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6*wPhi7);
         fNestedLoopsResultsSinPro->Fill(12.5,TMath::Sin(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4+h5*dPhi5+h6*dPhi6+h7*dPhi7),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6*wPhi7);
        } // end of for(Int_t i7=0;i7<nPrim;i7++)
       } // end of for(Int_t i6=0;i6<nPrim;i6++) 
      } // end of for(Int_t i5=0;i5<nPrim;i5++)
     } // end of for(Int_t i4=0;i4<nPrim;i4++)  
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=7)
 
 // 8-particle correlations:
 if(nPrim>=8)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   dPhi1=aftsTrack->Phi();
   if(fUsePhiWeights){wPhi1=PhiWeight(dPhi1);}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    dPhi2=aftsTrack->Phi();
    if(fUsePhiWeights){wPhi2=PhiWeight(dPhi2);}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     dPhi3=aftsTrack->Phi();
     if(fUsePhiWeights){wPhi3=PhiWeight(dPhi3);}
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3){continue;}
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())){continue;}
      dPhi4=aftsTrack->Phi();
      if(fUsePhiWeights){wPhi4=PhiWeight(dPhi4);}
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())){continue;}
       dPhi5=aftsTrack->Phi();
       if(fUsePhiWeights){wPhi5=PhiWeight(dPhi5);}
       for(Int_t i6=0;i6<nPrim;i6++)
       {
        if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
        aftsTrack=anEvent->GetTrack(i6);
        if(!(aftsTrack->InRPSelection())){continue;}
        dPhi6=aftsTrack->Phi();
        if(fUsePhiWeights){wPhi6=PhiWeight(dPhi6);}
        for(Int_t i7=0;i7<nPrim;i7++)
        {
         if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6){continue;}
         aftsTrack=anEvent->GetTrack(i7);
         if(!(aftsTrack->InRPSelection())){continue;}
         dPhi7=aftsTrack->Phi();
         if(fUsePhiWeights){wPhi7=PhiWeight(dPhi7);}
         for(Int_t i8=0;i8<nPrim;i8++)
         {
          if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7){continue;}
          aftsTrack=anEvent->GetTrack(i8);
          if(!(aftsTrack->InRPSelection())){continue;}
          dPhi8=aftsTrack->Phi();
          if(fUsePhiWeights){wPhi8=PhiWeight(dPhi8);}
          // Fill:   
          fNestedLoopsResultsCosPro->Fill(14.5,TMath::Cos(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4+h5*dPhi5+h6*dPhi6+h7*dPhi7+h8*dPhi8),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6*wPhi7*wPhi8);
          fNestedLoopsResultsSinPro->Fill(14.5,TMath::Sin(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4+h5*dPhi5+h6*dPhi6+h7*dPhi7+h8*dPhi8),wPhi1*wPhi2*wPhi3*wPhi4*wPhi5*wPhi6*wPhi7*wPhi8);
         } // end of for(Int_t i8=0;i8<nPrim;i8++)
        } // end of for(Int_t i7=0;i7<nPrim;i7++) 
       } // end of for(Int_t i6=0;i6<nPrim;i6++) 
      } // end of for(Int_t i5=0;i5<nPrim;i5++)
     } // end of for(Int_t i4=0;i4<nPrim;i4++)  
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=8)
 
 // *) Printout: TBI move somewhere else
 printf("\n cosine:");
 printf("\n  1-p => Q-vector:     %.12f",fNestedLoopsResultsCosPro->GetBinContent(2));
 printf("\n  1-p => Nested loops: %.12f",fNestedLoopsResultsCosPro->GetBinContent(1));
 printf("\n  2-p => Q-vector:     %.12f",fNestedLoopsResultsCosPro->GetBinContent(4));
 printf("\n  2-p => Nested loops: %.12f",fNestedLoopsResultsCosPro->GetBinContent(3));
 printf("\n  3-p => Q-vector:     %.12f",fNestedLoopsResultsCosPro->GetBinContent(6));
 printf("\n  3-p => Nested loops: %.12f",fNestedLoopsResultsCosPro->GetBinContent(5));
 printf("\n  4-p => Q-vector:     %.12f",fNestedLoopsResultsCosPro->GetBinContent(8));
 printf("\n  4-p => Nested loops: %.12f",fNestedLoopsResultsCosPro->GetBinContent(7));
 printf("\n  5-p => Q-vector:     %.12f",fNestedLoopsResultsCosPro->GetBinContent(10));
 printf("\n  5-p => Nested loops: %.12f",fNestedLoopsResultsCosPro->GetBinContent(9));
 printf("\n  6-p => Q-vector:     %.12f",fNestedLoopsResultsCosPro->GetBinContent(12));
 printf("\n  6-p => Nested loops: %.12f",fNestedLoopsResultsCosPro->GetBinContent(11));
 printf("\n  7-p => Q-vector:     %.12f",fNestedLoopsResultsCosPro->GetBinContent(14));
 printf("\n  7-p => Nested loops: %.12f",fNestedLoopsResultsCosPro->GetBinContent(13));
 printf("\n  8-p => Q-vector:     %.12f",fNestedLoopsResultsCosPro->GetBinContent(16));
 printf("\n  8-p => Nested loops: %.12f",fNestedLoopsResultsCosPro->GetBinContent(15));

 printf("\n\n sinus:");
 printf("\n  1-p => Q-vector:     %.12f",fNestedLoopsResultsSinPro->GetBinContent(2));
 printf("\n  1-p => Nested loops: %.12f",fNestedLoopsResultsSinPro->GetBinContent(1));
 printf("\n  2-p => Q-vector:     %.12f",fNestedLoopsResultsSinPro->GetBinContent(4));
 printf("\n  2-p => Nested loops: %.12f",fNestedLoopsResultsSinPro->GetBinContent(3));
 printf("\n  3-p => Q-vector:     %.12f",fNestedLoopsResultsSinPro->GetBinContent(6));
 printf("\n  3-p => Nested loops: %.12f",fNestedLoopsResultsSinPro->GetBinContent(5));
 printf("\n  4-p => Q-vector:     %.12f",fNestedLoopsResultsSinPro->GetBinContent(8));
 printf("\n  4-p => Nested loops: %.12f",fNestedLoopsResultsSinPro->GetBinContent(7));
 printf("\n  5-p => Q-vector:     %.12f",fNestedLoopsResultsSinPro->GetBinContent(10));
 printf("\n  5-p => Nested loops: %.12f",fNestedLoopsResultsSinPro->GetBinContent(9));
 printf("\n  6-p => Q-vector:     %.12f",fNestedLoopsResultsSinPro->GetBinContent(12));
 printf("\n  6-p => Nested loops: %.12f",fNestedLoopsResultsSinPro->GetBinContent(11));
 printf("\n  7-p => Q-vector:     %.12f",fNestedLoopsResultsSinPro->GetBinContent(14));
 printf("\n  7-p => Nested loops: %.12f",fNestedLoopsResultsSinPro->GetBinContent(13));
 printf("\n  8-p => Q-vector:     %.12f",fNestedLoopsResultsSinPro->GetBinContent(16));
 printf("\n  8-p => Nested loops: %.12f",fNestedLoopsResultsSinPro->GetBinContent(15));

 printf("\n\n"); 

} // void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckWithNestedLoops(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::FillQvector(AliFlowEventSimple *anEvent)
{
 // Fill Q-vector components.

 Int_t nTracks = anEvent->NumberOfTracks(); // TBI shall I promote this to data member?
 Double_t dPhi = 0., wPhi = 1.; // azimuthal angle and corresponding phi weight
 Double_t dPt = 0., wPt = 1.; // transverse momentum and corresponding pT weight
 Double_t dEta = 0., wEta = 1.; // pseudorapidity and corresponding eta weight
 Double_t wToPowerP = 1.; // weight raised to power p
 for(Int_t t=0;t<nTracks;t++) // loop over all tracks
 {
  AliFlowTrackSimple *pTrack = anEvent->GetTrack(t);
  if(!pTrack){printf("\n AAAARGH: pTrack is NULL in MPC::FillQvector(...) !!!!"); continue;}
  if(pTrack && pTrack->InRPSelection()) // fill Q-vector components only with reference particles
  {
   dPhi = pTrack->Phi(); // azimuthal angle
   if(fUsePhiWeights){wPhi = PhiWeight(dPhi);} // corresponding phi weight
   //if(dPhi < 0.){dPhi += TMath::TwoPi();} TBI
   //if(dPhi > TMath::TwoPi()){dPhi -= TMath::TwoPi();} TBI
   dPt = pTrack->Pt();
   if(fUsePtWeights){wPt = PtWeight(dPt);} // corresponding pT weight
   dEta = pTrack->Eta();
   if(fUseEtaWeights){wEta = EtaWeight(dEta);} // corresponding eta weight
   // Calculate Q-vector components:
   for(Int_t h=0;h<49;h++) // TBI hardwired 49 = maxHarmonic*maxCorrelator+1
   {
    for(Int_t p=0;p<9;p++) // TBI hardwired 9 = maxCorrelator+1
    {
     if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights){wToPowerP = pow(wPhi*wPt*wEta,p);} // TBI should I do something with the normalization of the product wPhi*wPt*wEta
     fQvector[h][p] += TComplex(wToPowerP*TMath::Cos(h*dPhi),wToPowerP*TMath::Sin(h*dPhi));
    } //  for(Int_t p=0;p<9;p++)
   } // for(Int_t h=0;h<49;h++)
  } // if(pTrack && pTrack->InRPSelection()) // fill Q-vector components only with reference particles
 } // for(Int_t t=0;t<nTracks;t++) // loop over all tracks

} // void AliFlowAnalysisWithMultiparticleCorrelations::FillQvector(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckSettings()
{
 // Cross-check all initial settings in this method. 
 
 // ...

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckSettings()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for weights;
 // c) Book and nest lists for correlations;
 // d) Book and nest lists for cumulants;
 // e) Book and nest lists for nested loops;
 // f) Book and nest lists for 'standard candles'.

 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("Control Histograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);

 // b) Book and nest lists for weights:
 fWeightsList = new TList();
 fWeightsList->SetName("Weights");
 fWeightsList->SetOwner(kTRUE);
 fHistList->Add(fWeightsList);

 // c) Book and nest lists for correlations:
 fCorrelationsList = new TList();
 fCorrelationsList->SetName("Correlations");
 fCorrelationsList->SetOwner(kTRUE);
 fHistList->Add(fCorrelationsList);

 // d) Book and nest lists for cumulants:
 fCumulantsList = new TList();
 fCumulantsList->SetName("Cumulants");
 fCumulantsList->SetOwner(kTRUE);
 fHistList->Add(fCumulantsList);

 // e) Book and nest lists for nested loops:
 fNestedLoopsList = new TList();
 fNestedLoopsList->SetName("Nested Loops");
 fNestedLoopsList->SetOwner(kTRUE);
 fHistList->Add(fNestedLoopsList);

 // f) Book and nest lists for 'standard candles':
 fStandardCandlesList = new TList();
 fStandardCandlesList->SetName("Standard Candles");
 fStandardCandlesList->SetOwner(kTRUE);
 fHistList->Add(fStandardCandlesList);

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::BookAndNestAllLists()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::WriteHistograms(TString outputFileName)
{
 // Store the final results in output file <outputFileName>.root.

 TFile *output = new TFile(outputFileName.Data(),"RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);

 delete output;

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::WriteHistograms(TString outputFileName)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::WriteHistograms(TDirectoryFile *outputFileName)
{
 // Store the final results in output file <outputFileName>.root.

 outputFileName->Add(fHistList);
 outputFileName->Write(outputFileName->GetName(),TObject::kSingleKey);

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::WriteHistograms(TDirectoryFile *outputFileName)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForControlHistograms()
{
 // Book all the stuff for control histograms.

 // a) Book the profile holding all the flags for control histograms;
 // b) Book all control histograms;
 //  b0) Book TH1D *fKinematicsHist[2][3]; 
 //  b1) Book TH1D *fMultDistributionsHist[3];
 //  b2) Book TH2D *fMultCorrelationsHist[3].

 // a) Book the profile holding all the flags for control histograms: TBI stil incomplete 
 fControlHistogramsFlagsPro = new TProfile("fControlHistogramsFlagsPro","Flags and settings for control histograms",1,0,1);
 fControlHistogramsFlagsPro->SetTickLength(-0.01,"Y");
 fControlHistogramsFlagsPro->SetMarkerStyle(25);
 fControlHistogramsFlagsPro->SetLabelSize(0.04);
 fControlHistogramsFlagsPro->SetLabelOffset(0.02,"Y");
 fControlHistogramsFlagsPro->SetStats(kFALSE);
 fControlHistogramsFlagsPro->SetFillColor(kGray);
 fControlHistogramsFlagsPro->SetLineColor(kBlack);
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(1,"fFillControlHistograms"); fControlHistogramsFlagsPro->Fill(0.5,fFillControlHistograms);
 fControlHistogramsList->Add(fControlHistogramsFlagsPro);

 if(!fFillControlHistograms){return;} // TBI is this safe? 

 // b) Book all control histograms: // TBI add setters for all these values
 //  b0) Book TH1D *fKinematicsHist[2][3]:
 Int_t nBins[2][3] = {{360,1000,1000},{360,1000,1000}}; // [RP,POI][phi,pt,eta]
 Double_t min[2][3] = {{0.,0.,-1.},{0.,0.,-1.}}; // [RP,POI][phi,pt,eta]
 Double_t max[2][3] = {{TMath::TwoPi(),10.,1.},{TMath::TwoPi(),10.,1.}}; // [RP,POI][phi,pt,eta]
 TString name[2][3] = {{"RP,phi","RP,pt","RP,eta"},{"POI,phi","POI,pt","POI,eta"}}; // [RP,POI][phi,pt,eta]
 TString title[2] = {"Reference particles (RP)","Particles of interest (POI)"}; // [RP,POI]
 Int_t lineColor[2] = {kBlue,kRed}; // [RP,POI]
 Int_t fillColor[2] = {kBlue-10,kRed-10}; // [RP,POI]
 TString xAxisTitle[3] = {"#phi","p_{T}","#eta"}; // [phi,pt,eta]
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   fKinematicsHist[rp][ppe] = new TH1D(name[rp][ppe].Data(),title[rp].Data(),nBins[rp][ppe],min[rp][ppe],max[rp][ppe]);
   fKinematicsHist[rp][ppe]->GetXaxis()->SetTitle(xAxisTitle[ppe].Data());
   fKinematicsHist[rp][ppe]->SetLineColor(lineColor[rp]);
   fKinematicsHist[rp][ppe]->SetFillColor(fillColor[rp]);
   fKinematicsHist[rp][ppe]->SetMinimum(0.); 
   fControlHistogramsList->Add(fKinematicsHist[rp][ppe]);
  }
 }

 //  b1) Book TH1D *fMultDistributionsHist[3]: // TBI add setters for all these values
 Int_t nBinsMult[3] = {3000,3000,3000}; // [RP,POI,reference multiplicity]
 Double_t minMult[3] = {0.,0.,0.}; // [RP,POI,reference multiplicity]
 Double_t maxMult[3] = {3000.,3000.,3000.}; // [RP,POI,reference multiplicity]
 TString nameMult[3] = {"Multiplicity (RP)","Multiplicity (POI)","Multiplicity (REF)"}; // [RP,POI,reference multiplicity]
 TString titleMult[3] = {"Reference particles (RP)","Particles of interest (POI)",""}; // [RP,POI,reference multiplicity]
 Int_t lineColorMult[3] = {kBlue,kRed,kGreen+2}; // [RP,POI,reference multiplicity]
 Int_t fillColorMult[3] = {kBlue-10,kRed-10,kGreen-10}; // [RP,POI,reference multiplicity]
 TString xAxisTitleMult[3] = {"Multiplicity (RP)","Multiplicity (POI)","Multiplicity (REF)"}; // [phi,pt,eta]
 for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 {
  fMultDistributionsHist[rprm] = new TH1D(nameMult[rprm].Data(),titleMult[rprm].Data(),nBinsMult[rprm],minMult[rprm],maxMult[rprm]);
  fMultDistributionsHist[rprm]->GetXaxis()->SetTitle(xAxisTitleMult[rprm].Data());
  fMultDistributionsHist[rprm]->SetLineColor(lineColorMult[rprm]);
  fMultDistributionsHist[rprm]->SetFillColor(fillColorMult[rprm]);
  fControlHistogramsList->Add(fMultDistributionsHist[rprm]);
 } // for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]

 //  b2) Book TH2D *fMultCorrelationsHist[3]: TBI too large objects to store in this way, perhaps,
 // ...
 fMultCorrelationsHist[0] = new TH2D("Multiplicity (RP vs. POI)","Multiplicity (RP vs. POI)",nBinsMult[0],minMult[0],maxMult[0],nBinsMult[1],minMult[1],maxMult[1]);
 fMultCorrelationsHist[0]->GetXaxis()->SetTitle(xAxisTitleMult[0].Data());
 fMultCorrelationsHist[0]->GetYaxis()->SetTitle(xAxisTitleMult[1].Data());
 fControlHistogramsList->Add(fMultCorrelationsHist[0]);

 // ...
 fMultCorrelationsHist[1] = new TH2D("Multiplicity (RP vs. REF)","Multiplicity (RP vs. REF)",nBinsMult[0],minMult[0],maxMult[0],nBinsMult[2],minMult[2],maxMult[2]);
 fMultCorrelationsHist[1]->GetXaxis()->SetTitle(xAxisTitleMult[0].Data());
 fMultCorrelationsHist[1]->GetYaxis()->SetTitle(xAxisTitleMult[2].Data());
 fControlHistogramsList->Add(fMultCorrelationsHist[1]);

 // ...
 fMultCorrelationsHist[2] = new TH2D("Multiplicity (POI vs. REF)","Multiplicity (POI vs. REF)",nBinsMult[1],minMult[1],maxMult[1],nBinsMult[2],minMult[2],maxMult[2]);
 fMultCorrelationsHist[2]->GetXaxis()->SetTitle(xAxisTitleMult[1].Data());
 fMultCorrelationsHist[2]->GetYaxis()->SetTitle(xAxisTitleMult[2].Data());
 fControlHistogramsList->Add(fMultCorrelationsHist[2]);

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForControlHistograms()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::FillControlHistograms(AliFlowEventSimple *anEvent)
{
 // Fill control histograms. 

 // a) Fill TH1D *fKinematicsHist[2][3];
 // b) Fill TH1D *fMultDistributionsHist[3]; 
 // c) Fill TH2D *fMultCorrelationsHist[3].  

 // a) Fill TH1D *fKinematicsHist[2][3]:
 Int_t nTracks = anEvent->NumberOfTracks(); // TBI shall I promote this to data member?
 for(Int_t t=0;t<nTracks;t++) // loop over all tracks
 {
  AliFlowTrackSimple *pTrack = anEvent->GetTrack(t);
  if(!pTrack){printf("\n AAAARGH: pTrack is NULL in MPC::FCH() !!!!");continue;}
  if(pTrack)
  {
   Double_t dPhi = pTrack->Phi(); 
   //if(dPhi < 0.){dPhi += TMath::TwoPi();} TBI
   //if(dPhi > TMath::TwoPi()){dPhi -= TMath::TwoPi();} TBI
   Double_t dPt = pTrack->Pt();
   Double_t dEta = pTrack->Eta();
   Double_t dPhiPtEta[3] = {dPhi,dPt,dEta};
   for(Int_t rp=0;rp<2;rp++) // [RP,POI]
   {
    for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
    {
     if((0==rp && pTrack->InRPSelection()) || (1==rp && pTrack->InPOISelection())) // TBI 
     { 
      fKinematicsHist[rp][ppe]->Fill(dPhiPtEta[ppe]);
     }
    } // for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
   } // for(Int_t rp=0;rp<2;rp++) // [RP,POI]
  } // if(pTrack)  
 } // for(Int_t t=0;t<nTracks;t++) // loop over all tracks

 // b) Fill TH1D *fMultDistributionsHist[3]: 
 Double_t dMultRP = anEvent->GetNumberOfRPs(); // TBI shall I promote these 3 variables into data members? 
 Double_t dMultPOI = anEvent->GetNumberOfPOIs();
 Double_t dMultREF = anEvent->GetReferenceMultiplicity();
 Double_t dMult[3] = {dMultRP,dMultPOI,dMultREF};
 for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 {
  fMultDistributionsHist[rprm]->Fill(dMult[rprm]);      
 } 

 // c) Fill TH2D *fMultCorrelationsHist[3]:  
 fMultCorrelationsHist[0]->Fill((Int_t)dMultRP,(Int_t)dMultPOI); // RP vs. POI
 fMultCorrelationsHist[1]->Fill((Int_t)dMultRP,(Int_t)dMultREF); // RP vs. refMult
 fMultCorrelationsHist[2]->Fill((Int_t)dMultPOI,(Int_t)dMultREF); // POI vs. refMult

} // void AliFlowAnalysisWithMultiparticleCorrelations::FillControlHistograms(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForControlHistograms()
{
 // Initialize all arrays for control histograms.

 // a) Initialize TH1D *fKinematicsHist[2][3];
 // b) Initialize TH1D *fMultDistributionsHist[3]; 
 // c) Initialize TH2D *fMultCorrelationsHist[3].  

 // a) Initialize TH1D *fKinematicsHist[2][3]:
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   fKinematicsHist[rp][ppe] = NULL;
  } 
 } 

 // b) Initialize TH1D *fMultDistributionsHist[3]:
 for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 {
  fMultDistributionsHist[rprm] = NULL;      
 } 

 // c) Initialize TH2D *fMultCorrelationsHist[3]: 
 for(Int_t r=0;r<3;r++) // [RP vs. POI, RP vs. refMult, POI vs. refMult]  
 {
  fMultCorrelationsHist[r] = NULL; 
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForControlHistograms()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForCorrelations()
{
 // Book all the stuff for correlations.

 // TBI this method can be implemented in a much more civilised way. 

 // a) Book the profile holding all the flags for correlations;
 // b) 1-p // TBI
 // c) 2-p // TBI
 // d) 3-p // TBI
 // e) 4-p // TBI
 // f) 5-p // TBI
 // g) 6-p // TBI
 // h) 7-p // TBI
 // i) 8-p // TBI

 // a) Book the profile holding all the flags for correlations:
 fCorrelationsFlagsPro = new TProfile("fCorrelationsFlagsPro","Flags for correlations",3,0,3);
 fCorrelationsFlagsPro->SetTickLength(-0.01,"Y");
 fCorrelationsFlagsPro->SetMarkerStyle(25);
 fCorrelationsFlagsPro->SetLabelSize(0.03);
 fCorrelationsFlagsPro->SetLabelOffset(0.02,"Y");
 fCorrelationsFlagsPro->SetStats(kFALSE);
 fCorrelationsFlagsPro->SetFillColor(kGray);
 fCorrelationsFlagsPro->SetLineColor(kBlack);
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateCorrelations"); fCorrelationsFlagsPro->Fill(0.5,fCalculateCorrelations); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(2,"fMaxHarmonic"); fCorrelationsFlagsPro->Fill(1.5,fMaxHarmonic); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(3,"fMaxCorrelator"); fCorrelationsFlagsPro->Fill(2.5,fMaxCorrelator); 
 fCorrelationsList->Add(fCorrelationsFlagsPro);

 if(!fCalculateCorrelations){return;} // TBI is this safe enough? 

 // b) 1-p // TBI
 fCorrelationsPro[0] = new TProfile("f1pCorrelationsPro","1-p correlations",6,0.,6.); // TBI better naming 
 fCorrelationsPro[0]->SetTickLength(-0.01,"Y");
 fCorrelationsPro[0]->SetMarkerStyle(25);
 fCorrelationsPro[0]->SetLabelSize(0.04);
 fCorrelationsPro[0]->SetLabelOffset(0.02,"Y");
 fCorrelationsPro[0]->SetStats(kFALSE);
 fCorrelationsPro[0]->Sumw2();
 for(Int_t ci=1;ci<=6;ci++) // correlation index
 {
  fCorrelationsPro[0]->GetXaxis()->SetBinLabel(ci,Form("#LT#LT1#GT#GT_{%d}",ci));
 } 
 fCorrelationsList->Add(fCorrelationsPro[0]);

 // c) 2-p // TBI
 fCorrelationsPro[1] = new TProfile("f2pCorrelationsPro","2-p correlations",6,0.,6.); // TBI better naming 
 fCorrelationsPro[1]->SetTickLength(-0.01,"Y");
 fCorrelationsPro[1]->SetMarkerStyle(25);
 fCorrelationsPro[1]->SetLabelSize(0.04);
 fCorrelationsPro[1]->SetLabelOffset(0.02,"Y");
 fCorrelationsPro[1]->SetStats(kFALSE);
 fCorrelationsPro[1]->Sumw2();
 for(Int_t ci=1;ci<=6;ci++) // correlation index
 {
  fCorrelationsPro[1]->GetXaxis()->SetBinLabel(ci,Form("#LT#LT2#GT#GT_{%d,%d}",ci,ci));
 } 
 fCorrelationsList->Add(fCorrelationsPro[1]);

 // d) 3-p // TBI
 fCorrelationsPro[2] = new TProfile("f3pCorrelationsPro","3-p correlations",1,0.,1.); // TBI better naming 
 fCorrelationsPro[2]->SetTickLength(-0.01,"Y");
 fCorrelationsPro[2]->SetMarkerStyle(25);
 fCorrelationsPro[2]->SetLabelSize(0.04);
 fCorrelationsPro[2]->SetLabelOffset(0.02,"Y");
 fCorrelationsPro[2]->SetStats(kFALSE);
 fCorrelationsPro[2]->Sumw2();
 //fCorrelationsPro[2]->GetXaxis()->SetBinLabel(1,"#LT#LT4#GT#GT_{3,2,-3,-2}"); // TBI automatize this
 fCorrelationsList->Add(fCorrelationsPro[2]);

 // e) 4-p // TBI
 fCorrelationsPro[3] = new TProfile("f4pCorrelationsPro","4-p correlations",4,0.,4.); // TBI better naming 
 fCorrelationsPro[3]->SetTickLength(-0.01,"Y");
 fCorrelationsPro[3]->SetMarkerStyle(25);
 fCorrelationsPro[3]->SetLabelSize(0.04);
 fCorrelationsPro[3]->SetLabelOffset(0.02,"Y");
 fCorrelationsPro[3]->SetStats(kFALSE);
 fCorrelationsPro[3]->Sumw2();
 fCorrelationsPro[3]->GetXaxis()->SetBinLabel(1,"#LT#LT4#GT#GT_{3,2,-3,-2}"); // TBI automatize this
 fCorrelationsPro[3]->GetXaxis()->SetBinLabel(2,"#LT#LT4#GT#GT_{4,2,-4,-2}"); // TBI automatize this
 fCorrelationsPro[3]->GetXaxis()->SetBinLabel(3,"#LT#LT4#GT#GT_{4,3,-4,-3}"); // TBI automatize this
 fCorrelationsPro[3]->GetXaxis()->SetBinLabel(4,"#LT#LT4#GT#GT_{5,3,-5,-3}"); // TBI automatize this
 fCorrelationsList->Add(fCorrelationsPro[3]);

 // f) 5-p // TBI
 fCorrelationsPro[4] = new TProfile("f5pCorrelationsPro","5-p correlations",1,0.,1.); // TBI better naming 
 fCorrelationsPro[4]->SetTickLength(-0.01,"Y");
 fCorrelationsPro[4]->SetMarkerStyle(25);
 fCorrelationsPro[4]->SetLabelSize(0.04);
 fCorrelationsPro[4]->SetLabelOffset(0.02,"Y");
 fCorrelationsPro[4]->SetStats(kFALSE);
 fCorrelationsPro[4]->Sumw2();
 //fCorrelationsPro[4]->GetXaxis()->SetBinLabel(1,"#LT#LT4#GT#GT_{3,2,-3,-2}"); // TBI automatize this
 fCorrelationsList->Add(fCorrelationsPro[4]);

 // g) 6-p // TBI
 fCorrelationsPro[5] = new TProfile("f6pCorrelationsPro","6-p correlations",1,0.,1.); // TBI better naming 
 fCorrelationsPro[5]->SetTickLength(-0.01,"Y");
 fCorrelationsPro[5]->SetMarkerStyle(25);
 fCorrelationsPro[5]->SetLabelSize(0.04);
 fCorrelationsPro[5]->SetLabelOffset(0.02,"Y");
 fCorrelationsPro[5]->SetStats(kFALSE);
 fCorrelationsPro[5]->Sumw2();
 //fCorrelationsPro[5]->GetXaxis()->SetBinLabel(1,"#LT#LT4#GT#GT_{3,2,-3,-2}"); // TBI automatize this
 fCorrelationsList->Add(fCorrelationsPro[5]);

 // h) 7-p // TBI
 fCorrelationsPro[6] = new TProfile("f7pCorrelationsPro","7-p correlations",1,0.,1.); // TBI better naming 
 fCorrelationsPro[6]->SetTickLength(-0.01,"Y");
 fCorrelationsPro[6]->SetMarkerStyle(25);
 fCorrelationsPro[6]->SetLabelSize(0.04);
 fCorrelationsPro[6]->SetLabelOffset(0.02,"Y");
 fCorrelationsPro[6]->SetStats(kFALSE);
 fCorrelationsPro[6]->Sumw2();
 //fCorrelationsPro[6]->GetXaxis()->SetBinLabel(1,"#LT#LT4#GT#GT_{3,2,-3,-2}"); // TBI automatize this
 fCorrelationsList->Add(fCorrelationsPro[6]);

 // i) 7-p // TBI
 fCorrelationsPro[7] = new TProfile("f8pCorrelationsPro","8-p correlations",1,0.,1.); // TBI better naming 
 fCorrelationsPro[7]->SetTickLength(-0.01,"Y");
 fCorrelationsPro[7]->SetMarkerStyle(25);
 fCorrelationsPro[7]->SetLabelSize(0.04);
 fCorrelationsPro[7]->SetLabelOffset(0.02,"Y");
 fCorrelationsPro[7]->SetStats(kFALSE);
 fCorrelationsPro[7]->Sumw2();
 //fCorrelationsPro[7]->GetXaxis()->SetBinLabel(1,"#LT#LT4#GT#GT_{3,2,-3,-2}"); // TBI automatize this
 fCorrelationsList->Add(fCorrelationsPro[7]);

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForCumulants()
{
 // Book all the stuff for cumualants.

 // a) Book the profile holding all the flags for cumulants;
 // b) 2p // TBI

 // a) Book the profile holding all the flags for cumulants:
 fCumulantsFlagsPro = new TProfile("fCumulantsFlagsPro","Flags for cumulants",1,0,1);
 fCumulantsFlagsPro->SetTickLength(-0.01,"Y");
 fCumulantsFlagsPro->SetMarkerStyle(25);
 fCumulantsFlagsPro->SetLabelSize(0.03);
 fCumulantsFlagsPro->SetLabelOffset(0.02,"Y");
 fCumulantsFlagsPro->SetStats(kFALSE);
 fCumulantsFlagsPro->SetFillColor(kGray);
 fCumulantsFlagsPro->SetLineColor(kBlack);
 fCumulantsFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateCumulants"); fCumulantsFlagsPro->Fill(0.5,fCalculateCumulants); 
 fCumulantsList->Add(fCumulantsFlagsPro);

 if(!fCalculateCumulants){return;} // TBI is this safe enough? 

 // b) 2p // TBI
 f2pCumulantsPro = new TProfile("f2pCumulantsPro","2-p cumulants",6,0.,6.);
 f2pCumulantsPro->SetTickLength(-0.01,"Y");
 f2pCumulantsPro->SetMarkerStyle(25);
 f2pCumulantsPro->SetLabelSize(0.04);
 f2pCumulantsPro->SetLabelOffset(0.02,"Y");
 f2pCumulantsPro->SetStats(kFALSE);
 f2pCumulantsPro->Sumw2();
 for(Int_t ci=1;ci<=6;ci++) // correlation index
 {
  f2pCumulantsPro->GetXaxis()->SetBinLabel(ci,Form("#LT#LT2#GT#GT_{C:%d,%d}",ci,ci)); // TBI
 } 
 fCumulantsList->Add(f2pCumulantsPro);
 
} // end of void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForNestedLoops()
{
 // Book all the stuff for nested loops.

 // a) Set default harmonic values; 
 // b) Book the profile holding all the flags for nested loops;
 // c) Book the profile holding all results for nested loops (cosine);
 // d) Book the profile holding all results for nested loops (sinus).

 // a) Set default harmonic values:
 Int_t h1=-1,h2=4,h3=-5,h4=6,h5=-6,h6=0,h7=-6,h8=2;
 // REMARK: This values can be overriden in a steering macro via 
 // mpc->GetNestedLoopsFlagsPro()->SetBinContent(<binNo>,<value>);

 // b) Book the profile holding all the flags for nested loops:
 fNestedLoopsFlagsPro = new TProfile("fNestedLoopsFlagsPro","Flags for nested loops",9,0,9);
 fNestedLoopsFlagsPro->SetTickLength(-0.01,"Y");
 fNestedLoopsFlagsPro->SetMarkerStyle(25);
 fNestedLoopsFlagsPro->SetLabelSize(0.03);
 fNestedLoopsFlagsPro->SetLabelOffset(0.02,"Y");
 fNestedLoopsFlagsPro->SetStats(kFALSE);
 fNestedLoopsFlagsPro->SetFillColor(kGray);
 fNestedLoopsFlagsPro->SetLineColor(kBlack);
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(1,"fCrossCheckWithNestedLoops"); fNestedLoopsFlagsPro->Fill(0.5,fCrossCheckWithNestedLoops);
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(2,"h_{1}"); fNestedLoopsFlagsPro->Fill(1.5,h1);
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(3,"h_{2}"); fNestedLoopsFlagsPro->Fill(2.5,h2);
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(4,"h_{3}"); fNestedLoopsFlagsPro->Fill(3.5,h3);
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(5,"h_{4}"); fNestedLoopsFlagsPro->Fill(4.5,h4);
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(6,"h_{5}"); fNestedLoopsFlagsPro->Fill(5.5,h5);
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(7,"h_{6}"); fNestedLoopsFlagsPro->Fill(6.5,h6);
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(8,"h_{7}"); fNestedLoopsFlagsPro->Fill(7.5,h7);
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(9,"h_{8}"); fNestedLoopsFlagsPro->Fill(8.5,h8);
 fNestedLoopsList->Add(fNestedLoopsFlagsPro);

 if(!fCrossCheckWithNestedLoops){return;} // TBI is this safe like this? 

 // c) Book the profile holding all results for nested loops (cosine):
 fNestedLoopsResultsCosPro = new TProfile("fNestedLoopsResultsCosPro","Nested loops results (cosine)",16,0,16);
 fNestedLoopsResultsCosPro->SetTickLength(-0.01,"Y");
 fNestedLoopsResultsCosPro->SetMarkerStyle(25);
 fNestedLoopsResultsCosPro->SetLabelSize(0.02);
 fNestedLoopsResultsCosPro->SetLabelOffset(0.02,"Y");
 fNestedLoopsResultsCosPro->SetStats(kFALSE);
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(1,Form("N: 1p(%d)",h1));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(2,Form("Q: 1p(%d)",h1));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(3,Form("N: 2p(%d,%d)",h1,h2));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(4,Form("Q: 2p(%d,%d)",h1,h2));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(5,Form("N: 3p(%d,%d,%d)",h1,h2,h3));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(6,Form("Q: 3p(%d,%d,%d)",h1,h2,h3));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(7,Form("N: 4p(%d,%d,%d,%d)",h1,h2,h3,h4));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(8,Form("Q: 4p(%d,%d,%d,%d)",h1,h2,h3,h4));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(9,Form("N: 5p(%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(10,Form("Q: 5p(%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(11,Form("N: 6p(%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(12,Form("Q: 6p(%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(13,Form("N: 7p(%d,%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6,h7));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(14,Form("Q: 7p(%d,%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6,h7));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(15,Form("N: 8p(%d,%d,%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6,h7,h8));
 fNestedLoopsResultsCosPro->GetXaxis()->SetBinLabel(16,Form("Q: 8p(%d,%d,%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6,h7,h8));
 fNestedLoopsList->Add(fNestedLoopsResultsCosPro);

 // d) Book the profile holding all results for nested loops (sinus):
 fNestedLoopsResultsSinPro = new TProfile("fNestedLoopsResultsSinPro","Nested loops results (sinus)",16,0,16);
 fNestedLoopsResultsSinPro->SetTickLength(-0.01,"Y");
 fNestedLoopsResultsSinPro->SetMarkerStyle(25);
 fNestedLoopsResultsSinPro->SetLabelSize(0.02);
 fNestedLoopsResultsSinPro->SetLabelOffset(0.02,"Y");
 fNestedLoopsResultsSinPro->SetStats(kFALSE);
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(1,Form("N: 1p(%d)",h1));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(2,Form("Q: 1p(%d)",h1));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(3,Form("N: 2p(%d,%d)",h1,h2));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(4,Form("Q: 2p(%d,%d)",h1,h2));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(5,Form("N: 3p(%d,%d,%d)",h1,h2,h3));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(6,Form("Q: 3p(%d,%d,%d)",h1,h2,h3));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(7,Form("N: 4p(%d,%d,%d,%d)",h1,h2,h3,h4));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(8,Form("Q: 4p(%d,%d,%d,%d)",h1,h2,h3,h4));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(9,Form("N: 5p(%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(10,Form("Q: 5p(%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(11,Form("N: 6p(%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(12,Form("Q: 6p(%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(13,Form("N: 7p(%d,%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6,h7));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(14,Form("Q: 7p(%d,%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6,h7));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(15,Form("N: 8p(%d,%d,%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6,h7,h8));
 fNestedLoopsResultsSinPro->GetXaxis()->SetBinLabel(16,Form("Q: 8p(%d,%d,%d,%d,%d,%d,%d,%d)",h1,h2,h3,h4,h5,h6,h7,h8));
 fNestedLoopsList->Add(fNestedLoopsResultsSinPro);

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForNestedLoops()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForStandardCandles()
{
 // Book all the stuff for 'standard candles'.

 // a) Book the profile holding all the flags for 'standard candles';
 // b) Book the histogram holding all results for 'standard candles'.

 // a) Book the profile holding all the flags for 'standard candles':
 fStandardCandlesFlagsPro = new TProfile("fStandardCandlesFlagsPro","Flags for standard candles",1,0,1);
 fStandardCandlesFlagsPro->SetTickLength(-0.01,"Y");
 fStandardCandlesFlagsPro->SetMarkerStyle(25);
 fStandardCandlesFlagsPro->SetLabelSize(0.03);
 fStandardCandlesFlagsPro->SetLabelOffset(0.02,"Y");
 fStandardCandlesFlagsPro->SetStats(kFALSE);
 fStandardCandlesFlagsPro->SetFillColor(kGray);
 fStandardCandlesFlagsPro->SetLineColor(kBlack);
 fStandardCandlesFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateStandardCandles"); fStandardCandlesFlagsPro->Fill(0.5,fCalculateStandardCandles);
 fStandardCandlesList->Add(fStandardCandlesFlagsPro);

 if(!fCalculateStandardCandles){return;} // TBI is this safe like this? 

 // b) Book the histogram holding all results for 'standard candles':
 fStandardCandlesHist = new TH1D("fStandardCandlesHist","Standard candles",1,0.,1.); 
 fStandardCandlesHist->SetStats(kFALSE);
 fStandardCandlesHist->GetXaxis()->SetBinLabel(1,"TBI");
 fStandardCandlesList->Add(fStandardCandlesHist);

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForStandardCandles()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetOutputHistograms(TList *histList)
{
 // Get pointers for everything and everywhere from the base list "fHistList". 

 // a) Get pointer for base list fHistList;
 // b) Get pointer for profile holding internal flags and, well, set again all flags;
 // c) Get pointers for control histograms;
 // d) Get pointers for correlations;
 // e) Get pointers for 'standard candles'.

 // a) Get pointer for base list fHistList and profile holding internal flags;
 fHistList = histList; 
 if(!fHistList){Fatal("AFAWMPC::GOH()","fHistList");}

 // b) Get pointer for profile holding internal flags and, well, set again all flags:
 fInternalFlagsPro = dynamic_cast<TProfile*>(fHistList->FindObject("fInternalFlagsPro"));
 if(!fInternalFlagsPro){Fatal("AFAWMPC::GOH()","fInternalFlagsPro");}
 fUseInternalFlags = fInternalFlagsPro->GetBinContent(1);
 fMinNoRPs = fInternalFlagsPro->GetBinContent(2);
 fMaxNoRPs = fInternalFlagsPro->GetBinContent(3);
 fExactNoRPs = fInternalFlagsPro->GetBinContent(4);

 // c) Get pointers for control histograms:
 this->GetPointersForControlHistograms(); 

 // d) Get pointers for correlations:
 this->GetPointersForCorrelations(); 

 // e) Get pointers for 'standard candles':
 this->GetPointersForStandardCandles(); 
  
} // void AliFlowAnalysisWithMultiparticleCorrelations::GetOutputHistograms(TList *histList)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForStandardCandles()
{
 // Get pointers for 'standard candles'.

 // a) Get pointer for fStandardCandlesList; TBI
 // b) Get pointer for fStandardCandlesFlagsPro; TBI
 // c) Set again all flags; TBI
 // d) Get pointer TH1D *fStandardCandlesHist; TBI

 // a) Get pointer for fStandardCandlesList: TBI
 fStandardCandlesList = dynamic_cast<TList*>(fHistList->FindObject("Standard Candles"));
 if(!fStandardCandlesList){Fatal("AFAWMPC::GPFSC()","fStandardCandlesList");}

 // b) Get pointer for fStandardCandlesFlagsPro: TBI
 fStandardCandlesFlagsPro = dynamic_cast<TProfile*>(fStandardCandlesList->FindObject("fStandardCandlesFlagsPro"));
 if(!fStandardCandlesFlagsPro){Fatal("AFAWMPC::GPFSC()","fStandardCandlesFlagsPro");}

 // c) Set again all flags: TBI
 fCalculateStandardCandles = fStandardCandlesFlagsPro->GetBinContent(1);

 if(!fCalculateStandardCandles){return;} // TBI is this safe enough

 // d) Get pointer TH1D *fStandardCandlesHist; TBI
 fStandardCandlesHist = dynamic_cast<TH1D*>(fStandardCandlesList->FindObject("fStandardCandlesHist"));
 
} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForStandardCandles()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForControlHistograms()
{
 // Get pointers for control histograms.

 // a) Get pointer for fControlHistogramsList; TBI
 // b) Get pointer for fControlHistogramsFlagsPro; TBI
 // c) Set again all flags; TBI
 // d) Get pointers to TH1D *fKinematicsHist[2][3]; TBI
 // e) Get pointers to TH1D *fMultDistributionsHist[3]; TBI
 // f) Get pointers to TH2D *fMultCorrelationsHist[3]. TBI

 // a) Get pointer for fControlHistogramsList: TBI
 fControlHistogramsList = dynamic_cast<TList*>(fHistList->FindObject("Control Histograms"));
 if(!fControlHistogramsList){Fatal("AFAWMPC::GPFCH()","fControlHistogramsList");}

 // b) Get pointer for fControlHistogramsFlagsPro: TBI
 fControlHistogramsFlagsPro = dynamic_cast<TProfile*>(fControlHistogramsList->FindObject("fControlHistogramsFlagsPro"));
 if(!fControlHistogramsFlagsPro){Fatal("AFAWMPC::GPFCH()","fControlHistogramsFlagsPro");}

 // c) Set again all flags: TBI
 fFillControlHistograms = fControlHistogramsFlagsPro->GetBinContent(1);

 if(!fFillControlHistograms){return;} // TBI is this safe enough

 // d) Get pointers to fKinematicsHist[2][3]: TBI
 TString name[2][3] = {{"RP,phi","RP,pt","RP,eta"},{"POI,phi","POI,pt","POI,eta"}}; // [RP,POI][phi,pt,eta]
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   fKinematicsHist[rp][ppe] = dynamic_cast<TH1D*>(fControlHistogramsList->FindObject(name[rp][ppe].Data()));
   if(!fKinematicsHist[rp][ppe]){Fatal("AFAWMPC::GPFCH()","%s",name[rp][ppe].Data());} // TBI 
  }
 }

 // e) Get pointers to TH1D *fMultDistributionsHist[3]:
 TString nameMult[3] = {"Multiplicity (RP)","Multiplicity (POI)","Multiplicity (REF)"}; // [RP,POI,reference multiplicity]
 for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 {
  fMultDistributionsHist[rprm] = dynamic_cast<TH1D*>(fControlHistogramsList->FindObject(nameMult[rprm].Data()));
  if(!fMultDistributionsHist[rprm]){Fatal("AFAWMPC::GPFCH()","%s",nameMult[rprm].Data());} // TBI 
 } // for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]

 // f) Get pointers to TH2D *fMultCorrelationsHist[3]: TBI automatize the things here...
 fMultCorrelationsHist[0] = dynamic_cast<TH2D*>(fControlHistogramsList->FindObject("Multiplicity (RP vs. POI)"));
 if(!fMultCorrelationsHist[0]){Fatal("AFAWMPC::GPFCH()","Multiplicity (RP vs. POI)");} // TBI 
 fMultCorrelationsHist[1] = dynamic_cast<TH2D*>(fControlHistogramsList->FindObject("Multiplicity (RP vs. REF)"));
 if(!fMultCorrelationsHist[1]){Fatal("AFAWMPC::GPFCH()","Multiplicity (RP vs. REF)");} // TBI 
 fMultCorrelationsHist[2] = dynamic_cast<TH2D*>(fControlHistogramsList->FindObject("Multiplicity (POI vs. REF)"));
 if(!fMultCorrelationsHist[2]){Fatal("AFAWMPC::GPFCH()","Multiplicity (POI vs. REF)");} // TBI 

} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForControlHistograms()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForCorrelations()
{
 // Get pointers for correlations.

 // a) Get pointer for fCorrelationsList; TBI
 // b) Get pointer for fCorrelationsFlagsPro; TBI
 // c) Set again all flags; TBI
 // d) Get pointers to TProfile *fCorrelationsPro[8].

 // a) Get pointer for fCorrelationsList: TBI
 fCorrelationsList = dynamic_cast<TList*>(fHistList->FindObject("Correlations"));
 if(!fCorrelationsList){Fatal("AFAWMPC::GPFC()","fCorrelationsList");}

 // b) Get pointer for fCorrelationsFlagsPro: TBI
 fCorrelationsFlagsPro = dynamic_cast<TProfile*>(fCorrelationsList->FindObject("fCorrelationsFlagsPro"));

 fCorrelationsFlagsPro = NULL;

 if(!fCorrelationsFlagsPro){Fatal("AFAWMPC::GPFC()","fCorrelationsFlagsPro");}

 // c) Set again all flags: 
 fCalculateCorrelations = fCorrelationsFlagsPro->GetBinContent(1);
 fMaxHarmonic = fCorrelationsFlagsPro->GetBinContent(2);
 fMaxCorrelator = fCorrelationsFlagsPro->GetBinContent(3);

 if(!fCalculateCorrelations){return;} // TBI is this safe enough

 // d) Get pointers to TProfile *fCorrelationsPro[8]:   
 for(Int_t c=0;c<8;c++)
 {
  fCorrelationsPro[c] = dynamic_cast<TProfile*>(fCorrelationsList->FindObject(Form("f%dpCorrelationsPro",c+1)));
  if(!fCorrelationsPro[c]){Fatal("AFAWMPC::GPFC()","f%dpCorrelationsPro",c+1);} 
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForQvector()
{
 // Initialize all arrays for Q-vector.

 for(Int_t h=0;h<49;h++) // harmonic TBI hardwired 49 = maxHarmonic*maxCorrelator+1
 {
  for(Int_t p=0;p<9;p++) // power TBI hardwired 9 = maxCorrelator+1
  {
   fQvector[h][p] = TComplex(0.,0.);
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForQvector()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::ResetQvector()
{
 // Reset all Q-vector components to zero before starting a new event. 

 for(Int_t h=0;h<49;h++) // harmonic TBI hardwired 49 = maxHarmonic*maxCorrelator+1
 {
  for(Int_t p=0;p<9;p++) // power TBI hardwired 9 = maxCorrelator+1
  {
   fQvector[h][p] = TComplex(0.,0.);
  } //  for(Int_t p=0;p<9;p++)
 } // for(Int_t h=0;h<49;h++)

} // void AliFlowAnalysisWithMultiparticleCorrelations::ResetQvector()

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Q(Int_t n, Int_t p)
{
 // Using the fact that Q{-n,p} = Q{n,p}^*. 
 
 if(n>=0){return fQvector[n][p];} 
 return TComplex::Conjugate(fQvector[-n][p]);
 
} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Q(Int_t n, Int_t p)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::One(Int_t n1)
{
 // Generic expression <exp[i(n1*phi1)]>. TBI comment

 TComplex one = Q(n1,1);

 return one;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::One(Int_t n1)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Two(Int_t n1, Int_t n2)
{
 // Generic two-particle correlation <exp[i(n1*phi1+n2*phi2)]>.

 TComplex two = Q(n1,1)*Q(n2,1)-Q(n1+n2,2);

 return two;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Two(Int_t n1, Int_t n2)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Three(Int_t n1, Int_t n2, Int_t n3)
{
 // Generic three-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3)]>.

 TComplex three = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
                - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3); 

 return three;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Three(Int_t n1, Int_t n2, Int_t n3)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)
{
 // Generic four-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4)]>.

 TComplex four = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
               - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
               + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
               + 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
               + 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);

 return four;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5)
{
 // Generic five-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5)]>.

 TComplex five = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
               - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
               + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
               + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
               + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
               - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
               + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
               - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
               + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
               + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
               - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
               + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
               - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
               - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
               + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
               + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
               + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
               + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
               - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
               + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
               + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
               + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
               + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
               - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3) 
               - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
               - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
 
 return five;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6)
{
 // Generic six-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6)]>.

 TComplex six = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
              - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
              - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
              - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
              + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
              - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
              - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
              - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
              - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
              - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
              - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
              - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
              + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
              - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
              - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
              - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
              - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
              + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
              + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
              - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
              - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
              - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
              + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
              - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
              + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
              + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
              - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
              - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
              - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
              - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
              + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
              - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
              - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
              - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
              - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
              + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
              - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
              - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
              - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
              + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
              - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
              + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
              - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
              - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
              - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
              + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
              + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
              - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
              - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
              - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
              - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
              + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
              + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
              - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
              - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
              - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
              - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
              + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
              - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
              - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
              - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
              + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
              + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
              + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
              + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
              + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
              + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
              - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
              + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
              - 120.*Q(n1+n2+n3+n4+n5+n6,6);

 return six;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7)
{
 // Generic seven-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7)]>.

 TComplex seven = Q(n1-n2+n2-n3+n4-n5+n6-n7,1); // TBI

 return seven;

} // end of TComplex AliFlowAnalysisWithMultiparticleCorrelations::Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)
{
 // Generic eight-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8)]>.

 TComplex eight = Q(n1-n2+n2-n3+n4-n5+n6-n7+n8,1); // TBI

 return eight;

} // end of TComplex AliFlowAnalysisWithMultiparticleCorrelations::Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForWeights()
{
 // Book all objects for calculations with weights. 

 // a) Book profile to hold all flags for weights;
 // b) Store histograms holding phi, pt and eta weights. 
    
 // a) Book profile to hold all flags for weights:
 fWeightsFlagsPro = new TProfile("fWeightsFlagsPro","0 = weight not used, 1 = weight used ",3,0,3);
 fWeightsFlagsPro->SetLabelSize(0.06);
 fWeightsFlagsPro->SetStats(kFALSE);
 fWeightsFlagsPro->SetFillColor(kGray);
 fWeightsFlagsPro->SetLineColor(kBlack);
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(1,"w_{#phi}"); fWeightsFlagsPro->Fill(0.5,fUsePhiWeights);
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(2,"w_{p_{T}}"); fWeightsFlagsPro->Fill(1.5,fUsePtWeights);
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(3,"w_{#eta}"); fWeightsFlagsPro->Fill(2.5,fUseEtaWeights);
 fWeightsList->Add(fWeightsFlagsPro); 
  
 // b) Store histograms holding phi, pt and eta weights:
 //    REMARK: It is assumed that these histos are accessed from external file "weights.root" 
 if(fPhiWeightsHist){fWeightsList->Add(fPhiWeightsHist);}
 if(fPtWeightsHist){fWeightsList->Add(fPtWeightsHist);}
 if(fEtaWeightsHist){fWeightsList->Add(fEtaWeightsHist);}

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForWeights()

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::PhiWeight(const Double_t &dPhi)
{
 // Determine phi weight for a given phi. 

 if(!fPhiWeightsHist){Fatal("AFAWMPC::PhiWeight()","fPhiWeightsHist");}

 Double_t wPhi = fPhiWeightsHist->GetBinContent(fPhiWeightsHist->FindBin(dPhi));

 return wPhi;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::PhiWeight(const Double_t &dPhi)

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::PtWeight(const Double_t &dPt)
{
 // Determine pt weight for a given pt. 

 if(!fPtWeightsHist){Fatal("AFAWMPC::PtWeight()","fPtWeightsHist");}

 Double_t wPt = fPtWeightsHist->GetBinContent(fPtWeightsHist->FindBin(dPt));

 return wPt;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::PtWeight(const Double_t &dPt)

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::EtaWeight(const Double_t &dEta)
{
 // Determine eta weight for a given eta. 

 if(!fEtaWeightsHist){Fatal("AFAWMPC::EtaWeight()","fEtaWeightsHist");}

 Double_t wEta = fEtaWeightsHist->GetBinContent(fEtaWeightsHist->FindBin(dEta));

 return wEta;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::EtaWeight(const Double_t &dEta)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForBase()
{
 // Book all base objects. 

 fInternalFlagsPro = new TProfile("fInternalFlagsPro","Internal flags and settings",4,0,4);
 fInternalFlagsPro->SetLabelSize(0.05);
 fInternalFlagsPro->SetStats(kFALSE);
 fInternalFlagsPro->SetFillColor(kGray);
 fInternalFlagsPro->SetLineColor(kBlack);
 fInternalFlagsPro->GetXaxis()->SetBinLabel(1,"fUseInternalFlags"); fInternalFlagsPro->Fill(0.5,fUseInternalFlags);  
 fInternalFlagsPro->GetXaxis()->SetBinLabel(2,"fMinNoRPs"); fInternalFlagsPro->Fill(1.5,fMinNoRPs);  
 fInternalFlagsPro->GetXaxis()->SetBinLabel(3,"fMaxNoRPs"); fInternalFlagsPro->Fill(2.5,fMaxNoRPs); 
 fInternalFlagsPro->GetXaxis()->SetBinLabel(4,"fExactNoRPs"); fInternalFlagsPro->Fill(3.5,fExactNoRPs);  
 fHistList->Add(fInternalFlagsPro); 

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForBase()

//=======================================================================================================================

Bool_t AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckInternalFlags(AliFlowEventSimple *anEvent)
{
 // Cross-check in this method wether "anEvent" passes internal flags. 

 // a) Cross-check min. and max. number of RPs. 
 // b) Cross-check... 

 Bool_t bPassesInternalFlags = kTRUE;

 // a) Cross-check min. and max. number of RPs: 
 fMinNoRPs <= anEvent->GetNumberOfRPs() && anEvent->GetNumberOfRPs() < fMaxNoRPs ? 1 : bPassesInternalFlags = kFALSE; // TBI can I leave 1 like this? 

 // ...

 return bPassesInternalFlags; 

} // Bool_t AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckInternalFlags(AliFlowEventSimple *anEvent)

//=======================================================================================================================


