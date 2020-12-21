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
using std::ofstream;
using std::ios;

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
 fPropagateError(kTRUE),
 fAnalysisTag(""),
 fDumpThePoints(kFALSE),
 fMaxNoEventsPerFile(100),
 fSelectRandomlyRPs(kFALSE),
 fnSelectedRandomlyRPs(-44),
 fRandomIndicesRPs(NULL),
 // 1.) Control histograms:
 fControlHistogramsList(NULL),
 fControlHistogramsFlagsPro(NULL),
 fFillControlHistograms(kFALSE),
 fFillKinematicsHist(kFALSE),
 fFillMultDistributionsHist(kFALSE),   
 fFillMultCorrelationsHist(kFALSE),
 fSkipSomeIntervals(kFALSE),
 fNumberOfSkippedRPParticles(0),
 // 2.) Q-vector:
 fQvectorList(NULL),       
 fQvectorFlagsPro(NULL),
 fCalculateQvector(kFALSE),
 fCalculateDiffQvectors(kFALSE),
 // 3.) Correlations:
 fCorrelationsList(NULL),
 fCorrelationsFlagsPro(NULL),
 fCalculateCorrelations(kFALSE),
 fMaxHarmonic(6), // TBI this shall not be hardwired in the ideal world...
 fMaxCorrelator(8),
 fCalculateIsotropic(kFALSE),
 fCalculateSame(kFALSE),
 fSkipZeroHarmonics(kFALSE),
 fCalculateSameIsotropic(kFALSE),
 fCalculateAll(kFALSE),
 fDontGoBeyond(0),
 fCalculateOnlyForHarmonicQC(kFALSE),
 fCalculateOnlyForSC(kFALSE),
 fCalculateOnlyCos(kFALSE),
 fCalculateOnlySin(kFALSE),
 // 4.) Event-by-event cumulants:
 fEbECumulantsList(NULL),
 fEbECumulantsFlagsPro(NULL),
 fCalculateEbECumulants(kFALSE),
 // 5.) Weights:
 fWeightsList(NULL),
 fWeightsFlagsPro(NULL),
 // 6.) Nested loops:
 fNestedLoopsList(NULL),
 fNestedLoopsFlagsPro(NULL),
 fCrossCheckWithNestedLoops(kFALSE),
 fCrossCheckDiffWithNestedLoops(kFALSE),
 fNestedLoopsResultsCosPro(NULL),
 fNestedLoopsResultsSinPro(NULL),
 fNestedLoopsDiffResultsPro(NULL),
 // 7.) 'Standard candles':
 fStandardCandlesList(NULL),
 fStandardCandlesFlagsPro(NULL),
 fCalculateStandardCandles(kFALSE),
 fPropagateErrorSC(kTRUE),
 fStandardCandlesHist(NULL),
 fProductsSCPro(NULL), 
 // 8.) Q-cumulants:
 fQcumulantsList(NULL),
 fQcumulantsFlagsPro(NULL),
 fCalculateQcumulants(kFALSE),
 fHarmonicQC(2),
 fPropagateErrorQC(kTRUE),
 fQcumulantsHist(NULL),
 fReferenceFlowHist(NULL),
 fProductsQCPro(NULL),
 // 9.) Differential correlations:
 fDiffCorrelationsList(NULL),
 fDiffCorrelationsFlagsPro(NULL),
 fCalculateDiffCorrelations(kFALSE),
 fCalculateDiffCos(kTRUE),
 fCalculateDiffSin(kFALSE),
 fCalculateDiffCorrelationsVsPt(kTRUE),
 fUseDefaultBinning(kTRUE),
 fnDiffBins(-44),
 fRangesDiffBins(NULL),
 fDiffBinNo(-1),
 // 10.) Symmetry plane correlations:
 fSymmetryPlanesList(NULL),
 fSymmetryPlanesFlagsPro(NULL),
 fCalculateSymmetryPlanes(kFALSE),
 // 11.) Eta gaps:
 fEtaGapsList(NULL),
 fEtaGapsFlagsPro(NULL),
 fCalculateEtaGaps(kFALSE),
 fLowestHarmonicEtaGaps(1),
 fHighestHarmonicEtaGaps(6)
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
  this->InitializeArraysForEbECumulants();
  this->InitializeArraysForWeights();
  this->InitializeArraysForQcumulants();
  this->InitializeArraysForDiffCorrelations();
  this->InitializeArraysForSymmetryPlanes();
  this->InitializeArraysForNestedLoops(); 
  this->InitializeArraysForEtaGaps(); 
  // TBI test for the workflow

  // c) Determine seed for gRandom:
  delete gRandom;
  gRandom = new TRandom3(0); // since 0 is in the argument, the seed is determined uniquely in space and time via TUUID
                             // TBI synchronize this eventually with seed set 'on-the-fly'   

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
 this->BookEverythingForQvector();
 this->BookEverythingForWeights();
 this->BookEverythingForCorrelations();
 this->BookEverythingForEbECumulants();
 this->BookEverythingForNestedLoops();
 this->BookEverythingForStandardCandles();
 this->BookEverythingForQcumulants();
 this->BookEverythingForDiffCorrelations(); 
 this->BookEverythingForSymmetryPlanes();
 this->BookEverythingForEtaGaps();

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
 // b) Cross-check all pointers used in this method;
 // c) Determine random indices;
 // d) Fill control histograms;
 // e) Fill Q-vector components;
 // f) Calculate multi-particle correlations from Q-vector components; 
 // g) Calculate e-b-e cumulants; 
 // h) Calculate symmetry plane correlations;
 // i) Calculate 2-p correlations with eta gaps;
 // j) Reset Q-vector components;
 // k) Cross-check results with nested loops;
 // l) Dump the points;
 // m) Reset array holding shuffled indices for RPs.

 // a) Cross-check internal flags:
 if(fUseInternalFlags){if(!this->CrossCheckInternalFlags(anEvent)){return;}}

 // b) Cross-check all pointers used in this method:
 this->CrossCheckPointersUsedInMake(); // TBI shall I call this method first  

 // c) Determine random indices:
 if(fnSelectedRandomlyRPs)
 {
  if(anEvent->GetNumberOfRPs() < fnSelectedRandomlyRPs){return;}
  this->DetermineRandomIndices(anEvent); 
  if(!fRandomIndicesRPs){return;} 
 } // TBI hw RPs in flag, make it more general

 // TBI temp gym: Remove this code eventually
 if(fSkipSomeIntervals)
 {
  fNumberOfSkippedRPParticles = 0;
  Int_t nTracks = anEvent->NumberOfTracks(); // TBI shall I promote this to data member?
  for(Int_t t=0;t<nTracks;t++) // loop over all tracks
  {
   AliFlowTrackSimple *pTrack = anEvent->GetTrack(t);
   if(!pTrack){printf("\n Error: pTrack is NULL in MPC::FCH() !!!!");continue;}
   if(!pTrack->InRPSelection()){continue;}
   if(pTrack)
   {
    Double_t dPhi = pTrack->Phi(); 
    //if(dPhi > TMath::TwoPi()){dPhi -= TMath::TwoPi();} TBI
    //if(dPhi < 0.){dPhi += TMath::TwoPi();} TBI
    Double_t dPt = pTrack->Pt();
    Double_t dEta = pTrack->Eta();
    Double_t dPhiPtEta[3] = {dPhi,dPt,dEta};

    // Skip some intervals: TBI promote eventually to AFTC class 
    Bool_t bPasses = kTRUE; 
    Bool_t bAlreadyCounted = kFALSE;
    for(Int_t ppe=0;ppe<3;ppe++)
    {
     if(!bPasses){break;} // found one kinematic window which shall be skipped
     for(Int_t b=0;b<10;b+=2)
     {
      if(-44==(Int_t)fSkip[ppe][b]){continue;}
      if(dPhiPtEta[ppe]>=fSkip[ppe][b] && dPhiPtEta[ppe]<fSkip[ppe][b+1])
      {
       bPasses = kFALSE; 
       if(!bAlreadyCounted)
       {
        fNumberOfSkippedRPParticles++;
        bAlreadyCounted = kTRUE;
       } // if(bAlreadyCounted)
       break;
      } // TBI this is a clear bug when this setter is used multiple times...
     } // for(Int_t b=0;b<10;b++)
    } // for(Int_t ppe=0;ppe<3;ppe++)
    if(!bPasses){continue;}

   } // if(pTrack)
  } // for(Int_t t=0;t<nTracks;t++) // loop over all tracks
 } // if(fSkipSomeIntervals)
 if(fSkipSomeIntervals) // TBI tmp gym
 {
  if(anEvent->GetNumberOfRPs() - fNumberOfSkippedRPParticles < 8){return;} // TBI tmp gym
 }



 // d) Fill control histograms:
 if(fFillControlHistograms){this->FillControlHistograms(anEvent);}
 
 // e) Fill Q-vector components:
 if(fCalculateQvector||fCalculateDiffQvectors){this->FillQvector(anEvent);}

 // f) Calculate multi-particle correlations from Q-vector components:
 if(fCalculateCorrelations){this->CalculateCorrelations(anEvent);}
 if(fCalculateDiffCorrelations){this->CalculateDiffCorrelations(anEvent);}

 // g) Calculate e-b-e cumulants: 
 if(fCalculateEbECumulants){this->CalculateEbECumulants(anEvent);}

 // h) Calculate symmetry plane correlations:
 if(fCalculateSymmetryPlanes){this->CalculateSymmetryPlanes(anEvent);}

 // i) Calculate 2-p correlations with eta gaps:
 if(fCalculateEtaGaps){this->CalculateEtaGaps(anEvent);}

 // j) Reset Q-vector components:
 if(fCalculateQvector||fCalculateDiffQvectors){this->ResetQvector();}

 // k) Cross-check results with nested loops:
 if(fCrossCheckWithNestedLoops){this->CrossCheckWithNestedLoops(anEvent);}
 if(fCrossCheckDiffWithNestedLoops){this->CrossCheckDiffWithNestedLoops(anEvent);}

 // l) Dump the points:
 if(fDumpThePoints){this->DumpThePoints(anEvent);}

 // m) Reset array holding shuffled indices for RPs:
 if(fSelectRandomlyRPs && fRandomIndicesRPs){delete fRandomIndicesRPs;}

} // end of AliFlowAnalysisWithMultiparticleCorrelations::Make(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::Finish()
{
 // Closing the curtains. 

 // a) Cross-check pointers used in this method;
 // b) Calculate 'standard candles';
 // c) Calculate Q-cumulants.
 
 // a) Cross-check pointers used in this method:
 this->CrossCheckPointersUsedInFinish();

 // b) Calculate 'standard candles':
 if(fCalculateStandardCandles){this->CalculateStandardCandles();}

 // c) Calculate Q-cumulants:
 if(fCalculateQcumulants){this->CalculateQcumulants();this->CalculateReferenceFlow();}

 // ...

 printf("\n  ... Closing the curtains ... \n\n");

} // end of AliFlowAnalysisWithMultiparticleCorrelations::Finish()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInFinish()
{
 // Cross-check all pointers used in method Finish().

 // a) Correlations;
 // b) 'Standard candles';
 // c) Q-cumulants.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInFinish()"; 

 // a) Correlations:
 if(fCalculateCorrelations)
 {
  for(Int_t cs=0;cs<2;cs++) 
  {
   if(fCalculateOnlyCos && 1==cs){continue;}
   else if(fCalculateOnlySin && 0==cs){continue;}
   for(Int_t c=0;c<fMaxCorrelator;c++) 
   {
    if(c==fDontGoBeyond){continue;}
    if(fCalculateOnlyForHarmonicQC && c%2==0){continue;}
    if(!fCorrelationsPro[cs][c]){Fatal(sMethodName.Data(),"fCorrelationsPro[%d][%d] is NULL, for one reason or another...",cs,c);}
   }
  }
  if(fCalculateQcumulants && fPropagateErrorQC && !fProductsQCPro)
   {Fatal(sMethodName.Data(),"fCalculateQcumulants && fPropagateErrorQC && !fProductsQCPro");}
 } // if(fCalculateCorrelations)

 // b) 'Standard candles':
 if(fCalculateStandardCandles)
 {
  if(!fStandardCandlesHist){Fatal(sMethodName.Data(),"fStandardCandlesHist is NULL, for one reason or another...");}
  if(fPropagateErrorSC)
  {
   if(!fProductsSCPro){Fatal(sMethodName.Data(),"fProductsSCPro is NULL, for one reason or another...");}
  }
 } // if(fCalculateStandardCandles)

 // c) Q-cumulants:
 if(fCalculateQcumulants)
 {
  if(!fQcumulantsHist){Fatal(sMethodName.Data(),"fQcumulantsHist is NULL, for one reason or another...");}
  if(!fReferenceFlowHist){Fatal(sMethodName.Data(),"fReferenceFlowHist is NULL, for one reason or another...");}
  if(fPropagateErrorQC)
  {
   if(!fProductsQCPro){Fatal(sMethodName.Data(),"fProductsQCPro is NULL, for one reason or another...");}
  }
 } // if(fCalculateQcumulants)

} // void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInFinish()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInMake()
{
 // Cross-check all pointers used in method Make().

 // a) Correlations;
 // b) Event-by-event cumulants;
 // c) ...

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInMake()"; 

 // a) Correlations:
 if(fCalculateCorrelations)
 {
  for(Int_t cs=0;cs<2;cs++) 
  {
   if(fCalculateOnlyCos && 1==cs){continue;}
   else if(fCalculateOnlySin && 0==cs){continue;}
   for(Int_t c=0;c<fMaxCorrelator;c++) 
   {
    if(c==fDontGoBeyond){continue;}
    if(fCalculateOnlyForHarmonicQC && c%2==0){continue;}
    if(!fCorrelationsPro[cs][c]){Fatal(sMethodName.Data(),"fCorrelationsPro[%d][%d] is NULL, for one reason or another...",cs,c);}
   }
  }
  if(fCalculateQcumulants && fPropagateErrorQC && !fProductsQCPro)
   {Fatal(sMethodName.Data(),"fCalculateQcumulants && fPropagateErrorQC && !fProductsQCPro");}
 } // if(fCalculateCorrelations)

 // b) Event-by-event cumulants:
 if(fCalculateEbECumulants)
 {
  for(Int_t cs=0;cs<2;cs++) 
  {
   for(Int_t c=0;c<fMaxCorrelator;c++) 
   {
    if(!fEbECumulantsPro[cs][c]){Fatal(sMethodName.Data(),"fEbECumulantsPro[%d][%d] is NULL, for one reason or another...",cs,c);}
   }
  }
 } // if(fCalculateEbECumulants)

} // void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInMake()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::DetermineRandomIndices(AliFlowEventSimple *anEvent)
{
 // Determine random indices.

 // Fisher-Yates algorithm:
 Int_t nPrim = anEvent->NumberOfTracks();
 if(nPrim > 0)
 {
  fRandomIndicesRPs = new TArrayI(nPrim);
 }
 else
 {
  return;
 }

 for(Int_t i=0;i<nPrim;i++)
 {
  fRandomIndicesRPs->AddAt(i,i);
 }
 for(Int_t i=nPrim-1;i>=1;i--)
 {
  Int_t j = gRandom->Integer(i+1);
  Int_t temp = fRandomIndicesRPs->GetAt(j);
  fRandomIndicesRPs->AddAt(fRandomIndicesRPs->GetAt(i),j);
  fRandomIndicesRPs->AddAt(temp,i);
 } // end of for(Int_t i=nPrim-1;i>=1;i--) 

} // void AliFlowAnalysisWithMultiparticleCorrelations::DetermineRandomIndices(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateStandardCandles()
{
 // Calculate 'standard candles'.

 // 'Standard candle' (SC) is defined in terms of average (all-event!) correlations as follows: 
 //   SC(-n1,-n2,n2,n1) = <<Cos(-n1,-n2,n2,n1)>> - <<Cos(-n1,n1)>>*<<Cos(-n2,n2)>>, n1 > n2.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CalculateStandardCandles()"; 

 Int_t nBins = fStandardCandlesHist->GetNbinsX(); 
 Int_t nBins2p = fCorrelationsPro[0][1]->GetNbinsX();
 Int_t nBins4p = fCorrelationsPro[0][3]->GetNbinsX();
 for(Int_t b=1;b<=nBins;b++)
 {
  // Standard candle:
  Double_t dSCn1n2n2n1 = 0.; // SC(-n1,-n2,n2,n1)
  Double_t dSCn1n2n2n1Err = 0.; // SC(-n1,-n2,n2,n1) stat. error
  fPropagateError = kTRUE;
 
  // Correlations:
  Double_t dCosn1n2n2n1 = 0.; // <<Cos(-n1,-n2,n2,n1)>>
  Double_t dCosn1n2n2n1Err = 0.; // <<Cos(-n1,-n2,n2,n1)>> stat. error
  Double_t dCosn1n1 = 0.; // <<Cos(-n1,n1)>>
  Double_t dCosn1n1Err = 0.; // <<Cos(-n1,n1)>> stat. error
  Double_t dCosn2n2 = 0.; // <<Cos(-n2,n2)>>
  Double_t dCosn2n2Err = 0.; // <<Cos(-n2,n2)>> stat. error
  
  // Labels:
  TString labeln1n2n2n1 = TString(fStandardCandlesHist->GetXaxis()->GetBinLabel(b)).ReplaceAll("SC","Cos");
  TString n1 = TString(fStandardCandlesHist->GetXaxis()->GetBinLabel(b))(4);  
  TString n2 = TString(fStandardCandlesHist->GetXaxis()->GetBinLabel(b))(7);
  if(n1.EqualTo("-") || n1.EqualTo(",")){Fatal(sMethodName.Data(),"n1.EqualTo...");}
  if(n2.EqualTo("-") || n2.EqualTo(",")){Fatal(sMethodName.Data(),"n2.EqualTo...");}
  TString labeln1n1 = Form("Cos(-%s,%s)",n1.Data(),n1.Data());
  TString labeln2n2 = Form("Cos(-%s,%s)",n2.Data(),n2.Data());
  //cout<<labeln1n2n2n1.Data()<<endl;
  //cout<<labeln1n1.Data()<<endl;
  //cout<<labeln2n2.Data()<<endl;
  //cout<<endl;  

  // Access <<Cos(-n1,-n2,n2,n1)>>:
  for(Int_t b4p=1;b4p<=nBins4p;b4p++)
  {
   if(labeln1n2n2n1.EqualTo(fCorrelationsPro[0][3]->GetXaxis()->GetBinLabel(b4p)))   
   {
    //cout<<labeln1n2n2n1.Data()<<endl;
    dCosn1n2n2n1 = fCorrelationsPro[0][3]->GetBinContent(b4p);
    dCosn1n2n2n1Err = fCorrelationsPro[0][3]->GetBinError(b4p);
    break; 
   }
  } // for(Int_t b4p=1;b4p<=nBins4p;b4p++)
  if(TMath::Abs(dCosn1n2n2n1) < 1.e-44)
  {
   cout<<Form("labeln1n2n2n1 = %s",labeln1n2n2n1.Data())<<endl;
   Warning(sMethodName.Data(),"TMath::Abs(dCosn1n2n2n1) < 1.e-44 !!!!");
  }

  // Access <<Cos(-n1,n1)>> and <<Cos(-n2,n2)>>:
  for(Int_t b2p=1;b2p<=nBins2p;b2p++)
  {
   if(labeln1n1.EqualTo(fCorrelationsPro[0][1]->GetXaxis()->GetBinLabel(b2p)))   
   {
    //cout<<labeln1n1.Data()<<endl;
    dCosn1n1 = fCorrelationsPro[0][1]->GetBinContent(b2p);
    dCosn1n1Err = fCorrelationsPro[0][1]->GetBinError(b2p);
   }  
   else if(labeln2n2.EqualTo(fCorrelationsPro[0][1]->GetXaxis()->GetBinLabel(b2p)))   
   {
    //cout<<labeln2n2.Data()<<endl;
    dCosn2n2 = fCorrelationsPro[0][1]->GetBinContent(b2p);
    dCosn2n2Err = fCorrelationsPro[0][1]->GetBinError(b2p);
   }
   if(TMath::Abs(dCosn1n1) > 0. && TMath::Abs(dCosn2n2) > 0.){break;} // found 'em both!
  } // for(Int_t b2p=1;b2p<=nBins2p;b2p++)
  if(TMath::Abs(dCosn1n1) < 1.e-44)
  {
   cout<<Form("labeln1n1 = %s",labeln1n1.Data())<<endl;
   Warning(sMethodName.Data(),"TMath::Abs(dCosn1n1) < 1.e-44 !!!!");
  }
  if(TMath::Abs(dCosn2n2) < 1.e-44)
  {
   cout<<Form("labeln2n2 = %s",labeln2n2.Data())<<endl;
   Warning(sMethodName.Data(),"TMath::Abs(dCosn2n2) < 1.e-44 !!!!");
  }

  // Calculate standard candles:
  dSCn1n2n2n1 = dCosn1n2n2n1-dCosn1n1*dCosn2n2;

  // Store the final results:
  fStandardCandlesHist->SetBinContent(b,dSCn1n2n2n1);

  // Error propagation:
  if(!fPropagateErrorSC)
  {
   fStandardCandlesHist->SetBinError(b,0.);
   continue;
  }

  // Access covariances (multiplied by weight dependent prefactor):
  Double_t wCovCosn1n2n2n1Cosn1n1 = Covariance(labeln1n2n2n1.Data(),labeln1n1.Data(),fProductsSCPro); // weighted Cov(<Cos(-n1,-n2,n2,n1)>,<Cos(-n1,n1)>)
  Double_t wCovCosn1n2n2n1Cosn2n2 = Covariance(labeln1n2n2n1.Data(),labeln2n2.Data(),fProductsSCPro); // weighted Cov(<Cos(-n1,-n2,n2,n1)>,<Cos(-n2,n2)>)
  Double_t wCovCosn1n1Cosn2n2 = Covariance(labeln1n1.Data(),labeln2n2.Data(),fProductsSCPro); // weighted Cov(<Cos(-n1,n1)>,<Cos(-n2,n2)>)

  // Explicit error propagation:
  Double_t dSCn1n2n2n1ErrSquared = pow(dCosn1n1,2.)*pow(dCosn2n2Err,2.) + pow(dCosn2n2,2.)*pow(dCosn1n1Err,2.) 
                                 + pow(dCosn1n2n2n1Err,2.) + 2.*dCosn1n1*dCosn2n2*wCovCosn1n1Cosn2n2 
                                 - 2.*dCosn1n1*wCovCosn1n2n2n1Cosn2n2 - 2.*dCosn2n2*wCovCosn1n2n2n1Cosn1n1;
  if(dSCn1n2n2n1ErrSquared > 0.)
  {
   dSCn1n2n2n1Err = pow(dSCn1n2n2n1ErrSquared,0.5);
  } else
    {
     Warning(sMethodName.Data(),"dSCn1n2n2n1ErrSquared > 0. is not satisfied for %s !!!!",labeln1n2n2n1.ReplaceAll("Cos","SC").Data());
     fPropagateError = kFALSE;
    }

  // Store the final stat. error:
  if(fPropagateError)
  {
   fStandardCandlesHist->SetBinError(b,dSCn1n2n2n1Err);
  }
 } // for(Int_t b=1;b<=nBins;b++)

 fPropagateError = kTRUE;

 return; 

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateStandardCandles()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForCorrelations()
{
 // Initialize all arrays for correlations.

 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  for(Int_t c=0;c<8;c++) // [1p,2p,...,8p]
  {
   fCorrelationsPro[cs][c] = NULL;
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForEbECumulants()
{
 // Initialize all arrays for event-by-event cumulants.

 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  for(Int_t c=0;c<8;c++) // [1p,2p,...,8p]
  {
   fEbECumulantsPro[cs][c] = NULL;
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForEbECumulants()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForDiffCorrelations()
{
 // Initialize all arrays for differential correlations.

 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  for(Int_t c=0;c<4;c++) // [1p,2p,3p,4p]
  {
   fDiffCorrelationsPro[cs][c] = NULL;
   fDiffHarmonics[cs][c] = 0;
  }
 }

 // Default values:
 // Cos, 2p:
 fDiffHarmonics[1][0] = -2;
 fDiffHarmonics[1][1] = 2;
 // Cos, 3p:
 fDiffHarmonics[2][0] = -3;
 fDiffHarmonics[2][1] = 1;
 fDiffHarmonics[2][2] = 2;
 // Cos, 4p:
 fDiffHarmonics[3][0] = -2;
 fDiffHarmonics[3][1] = -2;
 fDiffHarmonics[3][2] = 2;
 fDiffHarmonics[3][3] = 2;

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForDiffCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForSymmetryPlanes()
{
 // Initialize all arrays for symmetry plane correlations.

 for(Int_t gc=0;gc<1;gc++) // 'generic correlator': [[0]:(Psi2n,Psi1n),[1]:...] TBI upper boundary will change
 {
  for(Int_t n=0;n<2;n++) // 'harmonic n for generic correlator': [[0]:n=1,[1]:n=2,...] TBI upper boundary will change
  {
   fSymmetryPlanesPro[gc][n] = NULL;
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForSymmetryPlanes()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForNestedLoops()
{
 // Initialize all arrays for nested loops.  

 fCrossCheckDiffCSCOBN[0] = 0; // cos/sin
 fCrossCheckDiffCSCOBN[1] = 2; // correlator order
 fCrossCheckDiffCSCOBN[2] = 4; // bin number

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForNestedLoops()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForEtaGaps()
{
 // Initialize all arrays for eta gaps.  

 for(Int_t h=0;h<6;h++) // harmonic
 {
  fEtaGapsPro[h] = NULL; 
 } 

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForEtaGaps()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateCorrelations(AliFlowEventSimple *anEvent)
{
 // Calculate multi-particle correlations from Q-vector components.

 // a) Calculate all booked multi-particle correlations;
 // b) Calculate products needed for QC error propagation;
 // c) Calculate products needed for SC error propagation.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CalculateCorrelations(AliFlowEventSimple *anEvent)"; 
 if(!anEvent){Fatal(sMethodName.Data(),"'anEvent'!?!? You again!!!!");}

 // a) Calculate all booked multi-particle correlations:
 Double_t dMultRP = fSelectRandomlyRPs ? fnSelectedRandomlyRPs : anEvent->GetNumberOfRPs(); // TBI shall I promote this variable into data member? 
 if(fSkipSomeIntervals){ dMultRP = dMultRP - fNumberOfSkippedRPParticles; }
 
 for(Int_t cs=0;cs<2;cs++) // cos/sin 
 {
  if(fCalculateOnlyCos && 1==cs){continue;}
  else if(fCalculateOnlySin && 0==cs){continue;}
  for(Int_t co=0;co<8;co++) // correlator order (TBI hardwired 8) 
  {
   if(dMultRP < co+1){break;} // defines min. number of particles in an event for a certain correlator to make sense
   Int_t nBins = 0;
   if(fCorrelationsPro[cs][co]){nBins = fCorrelationsPro[cs][co]->GetNbinsX();}
   else{continue;}
   for(Int_t b=1;b<=nBins;b++)
   {
    TString sBinLabel = fCorrelationsPro[cs][co]->GetXaxis()->GetBinLabel(b);
    if(sBinLabel.EqualTo("")){break;} 
    Double_t num = CastStringToCorrelation(sBinLabel.Data(),kTRUE);
    Double_t den = CastStringToCorrelation(sBinLabel.Data(),kFALSE);
    Double_t weight = den; // TBI: add support for other options for the weight eventually
    if(den>0.) 
    {
     fCorrelationsPro[cs][co]->Fill(b-.5,num/den,weight);
    } else{Warning(sMethodName.Data(),"if(den>0.)");}
   } // for(Int_t b=1;b<=nBins;b++)
  } // for(Int_t co=0;co<8;co++) // correlator order (TBI hardwired 8) 
 } // for(Int_t cs=0;cs<=1;cs++) // cos/sin 

 // b) Calculate products needed for QC error propagation:
 if(fCalculateQcumulants && fPropagateErrorQC){this->CalculateProductsOfCorrelations(anEvent,fProductsQCPro);}

 // c) Calculate products needed for SC error propagation:
 if(fCalculateStandardCandles && fPropagateErrorSC){this->CalculateProductsOfCorrelations(anEvent,fProductsSCPro);}

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateCorrelations(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateDiffCorrelations(AliFlowEventSimple *anEvent)
{
 // Calculate differential multi-particle correlations from Q-, p- and q-vector components.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CalculateCorrelations(AliFlowEventSimple *anEvent)"; 
 if(!anEvent){Fatal(sMethodName.Data(),"'anEvent'!?!? You again!!!!");}

 Int_t nBins = 0; // TBI promote this to data member? 
 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  if(nBins != 0){break;}
  for(Int_t co=0;co<4;co++) // [1p,2p,3p,4p]
  {
   if(fDiffCorrelationsPro[cs][co] && 0==nBins)
   {
    nBins = fDiffCorrelationsPro[cs][co]->GetNbinsX(); 
   }
  } // for(Int_t co=0;co<4;co++) // [1p,2p,3p,4p]
 } // for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]

 // TBI: The lines below are genuine, most delicious, spaghetti ever... To be reimplemented (one day).
 if(fCalculateDiffCos)
 {
 for(Int_t b=1;b<=nBins;b++)
 {
  fDiffBinNo = b-1;
  // <2'>:  
  Double_t num2 = TwoDiff(fDiffHarmonics[1][0],fDiffHarmonics[1][1]).Re();
  Double_t den2 = TwoDiff(0,0).Re();
  Double_t w2 = den2; // TBI add support for other options for the weight
  if(den2>0.){fDiffCorrelationsPro[0][1]->Fill(fDiffCorrelationsPro[0][1]->GetBinCenter(b),num2/den2,w2);} 
  // <3'>:  
  Double_t num3 = ThreeDiff(fDiffHarmonics[2][0],fDiffHarmonics[2][1],fDiffHarmonics[2][2]).Re();
  Double_t den3 = ThreeDiff(0,0,0).Re();
  Double_t w3 = den3; // TBI add support for other options for the weight
  if(den3>0.){fDiffCorrelationsPro[0][2]->Fill(fDiffCorrelationsPro[0][2]->GetBinCenter(b),num3/den3,w3);} 
  // <4'>:  
  Double_t num4 = FourDiff(fDiffHarmonics[3][0],fDiffHarmonics[3][1],fDiffHarmonics[3][2],fDiffHarmonics[3][3]).Re();
  Double_t den4 = FourDiff(0,0,0,0).Re();
  Double_t w4 = den4; // TBI add support for other options for the weight
  if(den4>0.){fDiffCorrelationsPro[0][3]->Fill(fDiffCorrelationsPro[0][3]->GetBinCenter(b),num4/den4,w4);} 
 } // for(Int_t b=1;b<=nBins;b++)
 }
 // TBI: The lines below are genuine, most delicious, spaghetti ever... To be reimplemented (one day).
 if(fCalculateDiffSin)
 {
 for(Int_t b=1;b<=nBins;b++)
 {
  fDiffBinNo = b-1;
  // <2'>:  
  Double_t num2 = TwoDiff(fDiffHarmonics[1][0],fDiffHarmonics[1][1]).Im();
  Double_t den2 = TwoDiff(0,0).Re();
  Double_t w2 = den2; // TBI add support for other options for the weight
  if(den2>0.){fDiffCorrelationsPro[1][1]->Fill(fDiffCorrelationsPro[1][1]->GetBinCenter(b),num2/den2,w2);} 
  // <3'>:  
  Double_t num3 = ThreeDiff(fDiffHarmonics[2][0],fDiffHarmonics[2][1],fDiffHarmonics[2][2]).Im();
  Double_t den3 = ThreeDiff(0,0,0).Re();
  Double_t w3 = den3; // TBI add support for other options for the weight
  if(den3>0.){fDiffCorrelationsPro[1][2]->Fill(fDiffCorrelationsPro[1][2]->GetBinCenter(b),num3/den3,w3);} 
  // <4'>:  
  Double_t num4 = FourDiff(fDiffHarmonics[3][0],fDiffHarmonics[3][1],fDiffHarmonics[3][2],fDiffHarmonics[3][3]).Im();
  Double_t den4 = FourDiff(0,0,0,0).Re();
  Double_t w4 = den4; // TBI add support for other options for the weight
  if(den4>0.){fDiffCorrelationsPro[1][3]->Fill(fDiffCorrelationsPro[1][3]->GetBinCenter(b),num4/den4,w4);} 
 } // for(Int_t b=1;b<=nBins;b++)
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateDiffCorrelations(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateEtaGaps(AliFlowEventSimple *anEvent)
{
 // Calculate 2-p correlations with eta gaps.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CalculateEtaGaps(AliFlowEventSimple *anEvent)"; 
 if(!anEvent){Fatal(sMethodName.Data(),"'anEvent'? What's wrong with you today...");}

 TComplex Qa[6][11] = {{TComplex(0.,0.)}}; // -eta [harmonic][eta gap]
 Double_t Ma[6][11] = {{0.}}; // multiplicity for -eta TBI this shall not depend on harmonic, clearly
 TComplex Qb[6][11] = {{TComplex(0.,0.)}}; // +eta [harmonic][eta gap]
 Double_t Mb[6][11] = {{0.}}; // multiplicity for +eta TBI this shall not depend on harmonic, clearly
 Double_t dEtaGaps[11] = {1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0};

 Int_t nTracks = anEvent->NumberOfTracks(); // TBI shall I promote this to data member?
 Double_t dPhi=0.,dPt=0.,dEta=0.;
 Double_t wPhi=1.,wPt=1.,wEta=1.;
 Double_t wToPowerP=1.;
 for(Int_t t=0;t<nTracks;t++) // loop over all tracks
 {
  AliFlowTrackSimple *pTrack = anEvent->GetTrack(t);
  if(!pTrack){printf("\n pTrack is NULL in MPC::CalculateEtaGaps(AliFlowEventSimple *anEvent) !!!!"); continue;}
  if(!(pTrack->InRPSelection() || pTrack->InPOISelection())){printf("\n pTrack is neither RP nor POI !!!!"); continue;}

  if(pTrack->InRPSelection()) // fill Q-vector components only with reference particles
  {
   // Access kinematic variables for RP and corresponding weights:
   dPhi = pTrack->Phi(); // azimuthal angle
   if(fUseWeights[0][0]){wPhi = Weight(dPhi,"RP","phi");} // corresponding phi weight
   //if(dPhi < 0.){dPhi += TMath::TwoPi();} TBI
   //if(dPhi > TMath::TwoPi()){dPhi -= TMath::TwoPi();} TBI
   dPt = pTrack->Pt();
   if(fUseWeights[0][1]){wPt = Weight(dPt,"RP","pt");} // corresponding pT weight
   dEta = pTrack->Eta();
   if(fUseWeights[0][2]){wEta = Weight(dEta,"RP","eta");} // corresponding eta weight
   if(fUseWeights[0][0]||fUseWeights[0][1]||fUseWeights[0][2]){wToPowerP = wPhi*wPt*wEta;}
   // Calculate Qa and Qb vectors:
   if(dEta<0.) // Qa
   {
    for(Int_t h=fLowestHarmonicEtaGaps-1;h<=fHighestHarmonicEtaGaps-1;h++)
    {
     for(Int_t eg=0;eg<11;eg++) // eta gaps
     {  
      if(dEta<-1.*dEtaGaps[eg]/2.)  
      {
       Qa[h][eg] += TComplex(wToPowerP*TMath::Cos((h+1)*dPhi),wToPowerP*TMath::Sin((h+1)*dPhi));
       Ma[h][eg]+=wToPowerP;
      } 
     } // for(Int_t eg=0;eg<11;eg++) // eta gaps
    } // for(Int_t h=fLowestHarmonicEtaGaps-1;h<=fHighestHarmonicEtaGaps-1;h++)
   } // if(dEta<0.)
   else if(dEta>0.) // Qb
   {
    for(Int_t h=fLowestHarmonicEtaGaps-1;h<=fHighestHarmonicEtaGaps-1;h++)
    {
     for(Int_t eg=0;eg<11;eg++) // eta gaps
     {  
      if(dEta>dEtaGaps[eg]/2.)  
      {
       Qb[h][eg] += TComplex(wToPowerP*TMath::Cos((h+1)*dPhi),wToPowerP*TMath::Sin((h+1)*dPhi));
       Mb[h][eg]+=wToPowerP;
      } 
     } // for(Int_t eg=0;eg<11;eg++) // eta gaps
    } // for(Int_t h=fLowestHarmonicEtaGaps-1;h<=fHighestHarmonicEtaGaps-1;h++)
   }
  } // if(pTrack->InRPSelection())
 } // for(Int_t t=0;t<nTracks;t++) // loop over all tracks

 // Calculate 2-p correlations with eta gaps from Qa and Qb vectors:
 for(Int_t h=fLowestHarmonicEtaGaps-1;h<=fHighestHarmonicEtaGaps-1;h++)
 {
  for(Int_t eg=0;eg<11;eg++) // eta gaps
  {
   if(!(Qa[h][eg].Rho()>0. && Qb[h][eg].Rho()>0.)){continue;} 
   if(!(Ma[h][eg]>0. && Mb[h][eg]>0.)){continue;} 
   fEtaGapsPro[h]->Fill(eg+0.5,TComplex(Qa[h][eg]*TComplex::Conjugate(Qb[h][eg])).Re()/(Ma[h][eg]*Mb[h][eg]),Ma[h][eg]*Mb[h][eg]);
  } // for(Int_t eg=0;eg<11;eg++) // eta gaps
 } // for(Int_t h=fLowestHarmonicEtaGaps-1;h<=fHighestHarmonicEtaGaps-1;h++)

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateEtaGaps(AliFlowEventSimple *anEvent)

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::CastStringToCorrelation(const char *string, Bool_t numerator)
{
 // Cast string of the generic form Cos/Sin(-n_1,-n_2,...,n_{k-1},n_k) in the corresponding correlation value.
 // If you issue a call to this method with setting numerator = kFALSE, then you are getting back for free
 // the corresponding denumerator (a.k.a. weight 'number of combinations').

 // TBI:
 // a) add protection against cases a la:
 //     string = Cos(-3,-4,5,6,5,6,-3)
 //     method = Six(-3,-4,5,6,5,-3).Re()
 // b) cross-check with nested loops this method 

 Double_t dValue = 0.; // return value

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CastStringToCorrelation(const char *string, Bool_t numerator)"; 

 if(!(TString(string).BeginsWith("Cos") || TString(string).BeginsWith("Sin")))
 {
  cout<<Form("And the fatal string is... '%s'. Congratulations!!",string)<<endl; 
  Fatal(sMethodName.Data(),"!(TString(string).BeginsWith(...");
 }

 Bool_t bRealPart = kTRUE;
 if(TString(string).BeginsWith("Sin")){bRealPart = kFALSE;}

 Int_t n[8] = {0,0,0,0,0,0,0,0}; // harmonics, supporting up to 8p correlations
 UInt_t whichCorr = 0;   
 for(Int_t t=0;t<=TString(string).Length();t++)
 {
  if(TString(string[t]).EqualTo(",") || TString(string[t]).EqualTo(")")) // TBI this is just ugly
  {
   n[whichCorr] = string[t-1] - '0';
   if(TString(string[t-2]).EqualTo("-")){n[whichCorr] = -1*n[whichCorr];}
   if(!(TString(string[t-2]).EqualTo("-") 
      || TString(string[t-2]).EqualTo(",")
      || TString(string[t-2]).EqualTo("("))) // TBI relax this eventually to allow two-digits harmonics
   { 
    cout<<Form("And the fatal string is... '%s'. Congratulations!!",string)<<endl; 
    Fatal(sMethodName.Data(),"!(TString(string[t-2]).EqualTo(...");
   }
   whichCorr++;
   if(whichCorr>=9){Fatal(sMethodName.Data(),"whichCorr>=9");} // not supporting corr. beyond 8p 
  } // if(TString(string[t]).EqualTo(",") || TString(string[t]).EqualTo(")")) // TBI this is just ugly
 } // for(UInt_t t=0;t<=TString(string).Length();t++)

 switch(whichCorr)
 {
  case 1:
   if(!numerator){dValue = One(0).Re();}
   else if(bRealPart){dValue = One(n[0]).Re();} 
   else{dValue = One(n[0]).Im();}
  break;

  case 2: 
   if(!numerator){dValue = Two(0,0).Re();}
   else if(bRealPart){dValue = Two(n[0],n[1]).Re();}
   else{dValue = Two(n[0],n[1]).Im();}
  break;

  case 3: 
   if(!numerator){dValue = Three(0,0,0).Re();}
   else if(bRealPart){dValue = Three(n[0],n[1],n[2]).Re();}
   else{dValue = Three(n[0],n[1],n[2]).Im();}
  break;

  case 4: 
   if(!numerator){dValue = Four(0,0,0,0).Re();}
   else if(bRealPart){dValue = Four(n[0],n[1],n[2],n[3]).Re();}
   else{dValue = Four(n[0],n[1],n[2],n[3]).Im();}
  break;

  case 5: 
   if(!numerator){dValue = Five(0,0,0,0,0).Re();}
   else if(bRealPart){dValue = Five(n[0],n[1],n[2],n[3],n[4]).Re();}
   else{dValue = Five(n[0],n[1],n[2],n[3],n[4]).Im();}
  break;

  case 6: 
   if(!numerator){dValue = Six(0,0,0,0,0,0).Re();}
   else if(bRealPart){dValue = Six(n[0],n[1],n[2],n[3],n[4],n[5]).Re();}
   else{dValue = Six(n[0],n[1],n[2],n[3],n[4],n[5]).Im();}
  break;

  case 7: 
   if(!numerator){dValue = Seven(0,0,0,0,0,0,0).Re();}
   else if(bRealPart){dValue = Seven(n[0],n[1],n[2],n[3],n[4],n[5],n[6]).Re();}
   else{dValue = Seven(n[0],n[1],n[2],n[3],n[4],n[5],n[6]).Im();}
   break;

  case 8: 
   if(!numerator){dValue = Eight(0,0,0,0,0,0,0,0).Re();}
   else if(bRealPart){dValue = Eight(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]).Re();} 
   else{dValue = Eight(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]).Im();} 
  break;

  default:
   cout<<Form("And the fatal 'whichCorr' value is... %d. Congratulations!!",whichCorr)<<endl; 
   Fatal(sMethodName.Data(),"switch(whichCorr)"); 
 } // switch(whichCorr)
 
 return dValue;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::CastStringToCorrelation(const char *string, Bool_t numerator)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateProductsOfCorrelations(AliFlowEventSimple *anEvent, TProfile2D *profile2D)
{
 // Calculate products of multi-particle correlations (needed for error propagation).

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CalculateProductsOfCorrelations(AliFlowEventSimple *anEvent, TProfile2D *profile2D)"; 
 if(!anEvent){Fatal(sMethodName.Data(),"Sorry, 'anEvent' is on holidays.");} 
 if(!profile2D){Fatal(sMethodName.Data(),"Sorry, 'profile2D' is on holidays.");} 

 Int_t nBins = profile2D->GetXaxis()->GetNbins();
 for(Int_t bx=2;bx<=nBins;bx++)
 {
  for(Int_t by=1;by<bx;by++)
  {
   const char *binLabelX = profile2D->GetXaxis()->GetBinLabel(bx);
   const char *binLabelY = profile2D->GetYaxis()->GetBinLabel(by);
   Double_t numX = this->CastStringToCorrelation(binLabelX,kTRUE); // numerator
   Double_t denX = this->CastStringToCorrelation(binLabelX,kFALSE); // denominator
   Double_t wX = denX; // weight TBI add support for other options
   Double_t numY = this->CastStringToCorrelation(binLabelY,kTRUE); // numerator
   Double_t denY = this->CastStringToCorrelation(binLabelY,kFALSE); // denominator
   Double_t wY = denY; // weight TBI add support for other options
   if(TMath::Abs(denX) > 0. && TMath::Abs(denY) > 0.)
   {
    profile2D->Fill(bx-0.5,by-0.5,(numX/denX)*(numY/denY),wX*wY);
   } else
     {
      cout<<endl; 
      cout<<"Cannot calculate product for:"<<endl;    
      cout<<Form("binLabelX = %s",binLabelX)<<endl;
      cout<<Form("binLabelY = %s",binLabelY)<<endl;
      cout<<Form("anEvent->GetNumberOfRPs() = %d",anEvent->GetNumberOfRPs())<<endl; 
      cout<<Form("fNumberOfSkippedRPParticles = %d",fNumberOfSkippedRPParticles)<<endl; 
      Fatal(sMethodName.Data(),"if(TMath::Abs(denX) > 0. && TMath::Abs(denY) > 0.)");
     } // else
  } // for(Int_t by=1;by<bx;by++)
 } // for(Int_t bx=2;bx<=nBins;bx++)

} // void CalculateProductsOfCorrelations(AliFlowEventSimple *anEvent, TProfile2D *profile2D)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateEbECumulants(AliFlowEventSimple *anEvent)
{
 // Calculate e-b-e cumulants from Q-vector components.

 // TBI this mathod is far (very far, in fact) from being finalized :'(

 // a) Calculate and store e-b-e cumulants.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CalculateEbECumulants(AliFlowEventSimple *anEvent)"; 
 if(!anEvent){Fatal(sMethodName.Data(),"'anEvent'!?!? You again!!!!");}

 // a) Calculate and store e-b-e cumulants:
 Double_t dMultRP = anEvent->GetNumberOfRPs(); // TBI shall I promote this variable into data member? 

 if(fSkipSomeIntervals){ dMultRP = dMultRP - fNumberOfSkippedRPParticles; } // TBI tmp gym

 Int_t binNo[8]; for(Int_t c=0;c<8;c++){binNo[c]=1;} 
 // 1-p:
 for(Int_t n1=-fMaxHarmonic;n1<=fMaxHarmonic;n1++) 
 {
  if(fSkipZeroHarmonics && 0==n1){continue;}
  if(fCalculateAll)
  {
   TComplex oneN = One(n1); // numerator
   Double_t oneD = One(0).Re(); // denominator
   Double_t oneW = oneD; // weight TBI add other possibilities here for the weight
   if(oneD>0. && dMultRP>=1) 
   {
    fEbECumulantsPro[0][0]->Fill(binNo[0]-.5,oneN.Re()/oneD,oneW);
    fEbECumulantsPro[1][0]->Fill(binNo[0]++-.5,oneN.Im()/oneD,oneW);
   } else {Warning(sMethodName.Data(),"if(oneD>0. && dMultRP>=1) ");}
  } 
  if(1==fDontGoBeyond){continue;}
  // 2-p:
  for(Int_t n2=n1;n2<=fMaxHarmonic;n2++) 
  {
   if(fSkipZeroHarmonics && 0==n2){continue;}
   if(fCalculateAll 
      || (fCalculateIsotropic && 0==n1+n2) 
      || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2)) 
      || (fCalculateSameIsotropic && 0==n1+n2 &&  TMath::Abs(n1)==TMath::Abs(n2)))
   {
    Double_t cumulants2pCos = Two(n1,n2).Re()/Two(0,0).Re() 
                            - (One(n1).Re()/One(0).Re())*(One(n2).Re()/One(0).Re())
                            + (One(n1).Im()/One(0).Re())*(One(n2).Im()/One(0).Re());
                            
    Double_t cumulants2pSin = Two(n1,n2).Im()/Two(0,0).Re() 
                            - (One(n1).Re()/One(0).Re())*(One(n2).Im()/One(0).Re())
                            - (One(n2).Re()/One(0).Re())*(One(n1).Im()/One(0).Re());

    if(/*twoD>0. &&*/ dMultRP>=2) 
    {
     fEbECumulantsPro[0][1]->Fill(binNo[1]-.5,cumulants2pCos,1.);;
     fEbECumulantsPro[1][1]->Fill(binNo[1]++-.5,cumulants2pSin,1.);;
    } else {Warning(sMethodName.Data(),"/*twoD>0. &&*/ dMultRP>=2");} 
   } 
   if(2==fDontGoBeyond){continue;}
   
   /*

   // 3-p:
   for(Int_t n3=n2;n3<=fMaxHarmonic;n3++) 
   {
    if(fSkipZeroHarmonics && 0==n3){continue;}
    if(fCalculateAll 
       || (fCalculateIsotropic && 0==n1+n2+n3) 
       || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3))
       || (fCalculateSameIsotropic && 0==n1+n2+n3 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3)))
    {
     TComplex threeN = Three(n1,n2,n3); // numerator
     Double_t threeD = Three(0,0,0).Re(); // denominator
     Double_t threeW = threeD; // weight TBI add other possibilities here for the weight
     if(threeD>0. && dMultRP>=3) 
     {
      fEbECumulantsPro[0][2]->Fill(binNo[2]-.5,threeN.Re()/threeD,threeW);
      fEbECumulantsPro[1][2]->Fill(binNo[2]++-.5,threeN.Im()/threeD,threeW);
     } else {Warning(sMethodName.Data(),"threeD>0. && dMultRP>=3");} 
    }
    if(3==fDontGoBeyond){continue;}
    // 4-p:
    for(Int_t n4=n3;n4<=fMaxHarmonic;n4++) 
    {
     if(fSkipZeroHarmonics && 0==n4){continue;}
     if(fCalculateAll 
        || (fCalculateIsotropic && 0==n1+n2+n3+n4) 
        || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) 
            && TMath::Abs(n1)==TMath::Abs(n4))
        || (fCalculateSameIsotropic && 0==n1+n2+n3+n4 && TMath::Abs(n1)==TMath::Abs(n2) 
            && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4)))
     { 
      TComplex fourN = Four(n1,n2,n3,n4); // numerator
      Double_t fourD = Four(0,0,0,0).Re(); // denominator
      Double_t fourW = fourD; // weight TBI add other possibilities here for the weight
      if(fourD>0. && dMultRP>=4) 
      {
       fEbECumulantsPro[0][3]->Fill(binNo[3]-.5,fourN.Re()/fourD,fourW);
       fEbECumulantsPro[1][3]->Fill(binNo[3]++-.5,fourN.Im()/fourD,fourW);
      } else {Warning(sMethodName.Data(),"fourD>0. && dMultRP>=4");}
     }
     if(4==fDontGoBeyond){continue;}
     // 5-p:
     for(Int_t n5=n4;n5<=fMaxHarmonic;n5++) 
     {
      if(fSkipZeroHarmonics && 0==n5){continue;}
      if(fCalculateAll 
         || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5) 
         || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) 
             && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5))
         || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5 && TMath::Abs(n1)==TMath::Abs(n2) 
             && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5)))
      { 
       TComplex fiveN = Five(n1,n2,n3,n4,n5); // numerator
       Double_t fiveD = Five(0,0,0,0,0).Re(); // denominator
       Double_t fiveW = fiveD; // weight TBI add other possibilities here for the weight
       if(fiveD>0. && dMultRP>=5) 
       {
        fEbECumulantsPro[0][4]->Fill(binNo[4]-.5,fiveN.Re()/fiveD,fiveW);
        fEbECumulantsPro[1][4]->Fill(binNo[4]++-.5,fiveN.Im()/fiveD,fiveW);
       } else {Warning(sMethodName.Data(),"fiveD>0. && dMultRP>=5");}
      } 
      if(5==fDontGoBeyond){continue;}
      // 6-p:  
      for(Int_t n6=n5;n6<=fMaxHarmonic;n6++) 
      {
       if(fSkipZeroHarmonics && 0==n6){continue;}
       if(fCalculateAll 
          || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5+n6)  
          || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) 
              && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6))
          || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5+n6 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) 
              && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6)))
       { 
        TComplex sixN = Six(n1,n2,n3,n4,n5,n6); // numerator
        Double_t sixD = Six(0,0,0,0,0,0).Re(); // denominator
        Double_t sixW = sixD; // weight TBI add other possibilities here for the weight
        if(sixD>0. && dMultRP>=6) 
        {
         fEbECumulantsPro[0][5]->Fill(binNo[5]-.5,sixN.Re()/sixD,sixW);
         fEbECumulantsPro[1][5]->Fill(binNo[5]++-.5,sixN.Im()/sixD,sixW);
        } else {Warning(sMethodName.Data(),"sixD>0. && dMultRP>=6");}
       } 
       if(6==fDontGoBeyond){continue;}
       // 7-p:
       for(Int_t n7=n6;n7<=fMaxHarmonic;n7++) 
       {
        if(fSkipZeroHarmonics && 0==n7){continue;}
        if(fCalculateAll 
           || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5+n6+n7) 
           || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) 
               && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) && TMath::Abs(n1)==TMath::Abs(n7))
           || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5+n6+n7 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3)
               && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) 
               && TMath::Abs(n1)==TMath::Abs(n7)))
        { 
         TComplex sevenN = Seven(n1,n2,n3,n4,n5,n6,n7); // numerator
         Double_t sevenD = Seven(0,0,0,0,0,0,0).Re(); // denominator
         Double_t sevenW = sevenD; // weight TBI add other possibilities here for the weight
         if(sevenD>0. && dMultRP>=7) 
         {
          fEbECumulantsPro[0][6]->Fill(binNo[6]-.5,sevenN.Re()/sevenD,sevenW);
          fEbECumulantsPro[1][6]->Fill(binNo[6]++-.5,sevenN.Im()/sevenD,sevenW);
         } else {Warning(sMethodName.Data(),"sevenD>0. && dMultRP>=7");}
        } 
        if(7==fDontGoBeyond){continue;}
        // 8-p:
        for(Int_t n8=n7;n8<=fMaxHarmonic;n8++) 
        {
         if(fSkipZeroHarmonics && 0==n8){continue;}
         if(fCalculateAll 
            || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5+n6+n7+n8) 
            || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) 
                && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) && TMath::Abs(n1)==TMath::Abs(n7) 
                && TMath::Abs(n1)==TMath::Abs(n8))
            || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5+n6+n7+n8 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) 
                && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) 
                && TMath::Abs(n1)==TMath::Abs(n7) && TMath::Abs(n1)==TMath::Abs(n8)))
         { 
          TComplex eightN = Eight(n1,n2,n3,n4,n5,n6,n7,n8); // numerator
          Double_t eightD = Eight(0,0,0,0,0,0,0,0).Re(); // denominator
          Double_t eightW = eightD; // weight TBI add other possibilities here for the weight
          if(eightD>0. && dMultRP>=8) 
          {
           fEbECumulantsPro[0][7]->Fill(binNo[7]-.5,eightN.Re()/eightD,eightW);
           fEbECumulantsPro[1][7]->Fill(binNo[7]++-.5,eightN.Im()/eightD,eightW);
          }
         } 
        } // for(Int_t n8=n7;n8<=fMaxHarmonic;n8++)
       } // for(Int_t n7=n6;n7<=fMaxHarmonic;n7++) 
      } // for(Int_t n6=n5;n6<=fMaxHarmonic;n6++) 
     } // for(Int_t n5=n4;n5<=fMaxHarmonic;n5++) 
    } // for(Int_t n4=n3;n4<=fMaxHarmonic;n4++)   
   } // for(Int_t n3=n2;n3<=fMaxHarmonic;n3++) 
 
  */

  } // for(Int_t n2=n1;n2<=fMaxHarmonic;n2++)
 } // for(Int_t n1=-fMaxHarmonic;n1<=fMaxHarmonic;n1++) 

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateEbECumulants(AliFlowEventSimple *anEvent)

//=======================================================================================================================

Bool_t AliFlowAnalysisWithMultiparticleCorrelations::TrackIsInSpecifiedIntervals(AliFlowTrackSimple *pTrack)
{
 // TBI

 if(!pTrack){exit(0);} // TBI

 Double_t dPhi = pTrack->Phi();
 Double_t dPt = pTrack->Pt();
 Double_t dEta = pTrack->Eta();
 Double_t dPhiPtEta[3] = {dPhi,dPt,dEta};

 // Skip some intervals: TBI promote eventually to AFTC class 
 Bool_t bPasses = kTRUE; 
 for(Int_t ppe=0;ppe<3;ppe++)
 {
  if(!bPasses){break;} // found one kinematic window which shall be skipped
  for(Int_t b=0;b<10;b+=2)
  {
   if(-44==(Int_t)fSkip[ppe][b]){continue;}
   if(dPhiPtEta[ppe]>=fSkip[ppe][b] && dPhiPtEta[ppe]<fSkip[ppe][b+1])
   {
    bPasses = kFALSE; 
    break;
   } // TBI this is a clear bug when this setter is used multiple times...
  } // for(Int_t b=0;b<10;b++)
 } // for(Int_t ppe=0;ppe<3;ppe++)

 return bPasses;  

} // Bool_t AliFlowAnalysisWithMultiparticleCorrelations::TrackIsInSpecifiedIntervals(AliFlowTrackSimple *pTrack)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckWithNestedLoops(AliFlowEventSimple *anEvent)
{
 // Cross-check results for multi-particle correlations with nested loops.

 // TBI add few comments here, there and over there
 // TBI this method is rather messy :'(

 Int_t h1 = (Int_t)fNestedLoopsFlagsPro->GetBinContent(2);
 Int_t h2 = (Int_t)fNestedLoopsFlagsPro->GetBinContent(3);
 Int_t h3 = (Int_t)fNestedLoopsFlagsPro->GetBinContent(4);
 Int_t h4 = (Int_t)fNestedLoopsFlagsPro->GetBinContent(5);
 Int_t h5 = (Int_t)fNestedLoopsFlagsPro->GetBinContent(6);
 Int_t h6 = (Int_t)fNestedLoopsFlagsPro->GetBinContent(7);
 Int_t h7 = (Int_t)fNestedLoopsFlagsPro->GetBinContent(8);
 Int_t h8 = (Int_t)fNestedLoopsFlagsPro->GetBinContent(9);

 this->ResetQvector();
 this->FillQvector(anEvent);

 if(TMath::Abs(One(0).Re())>0.)
 {
  fNestedLoopsResultsCosPro->Fill(1.5,One(h1).Re()/One(0).Re(),One(0).Re()); 
  fNestedLoopsResultsSinPro->Fill(1.5,One(h1).Im()/One(0).Re(),One(0).Re());  
 } 
 if(TMath::Abs(Two(0,0).Re())>0.)
 {
  fNestedLoopsResultsCosPro->Fill(3.5,Two(h1,h2).Re()/Two(0,0).Re(),Two(0,0).Re()); 
  fNestedLoopsResultsSinPro->Fill(3.5,Two(h1,h2).Im()/Two(0,0).Re(),Two(0,0).Re()); 
 }
 if(TMath::Abs(Three(0,0,0).Re())>0.)
 {
  fNestedLoopsResultsCosPro->Fill(5.5,Three(h1,h2,h3).Re()/Three(0,0,0).Re(),Three(0,0,0).Re()); 
  fNestedLoopsResultsSinPro->Fill(5.5,Three(h1,h2,h3).Im()/Three(0,0,0).Re(),Three(0,0,0).Re()); 
 } 
 if(TMath::Abs(Four(0,0,0,0).Re())>0.)
 {
  fNestedLoopsResultsCosPro->Fill(7.5,Four(h1,h2,h3,h4).Re()/Four(0,0,0,0).Re(),Four(0,0,0,0).Re()); 
  fNestedLoopsResultsSinPro->Fill(7.5,Four(h1,h2,h3,h4).Im()/Four(0,0,0,0).Re(),Four(0,0,0,0).Re()); 
 } 
 if(TMath::Abs(Five(0,0,0,0,0).Re())>0.)
 {
  fNestedLoopsResultsCosPro->Fill(9.5,Five(h1,h2,h3,h4,h5).Re()/Five(0,0,0,0,0).Re(),Five(0,0,0,0,0).Re()); 
  fNestedLoopsResultsSinPro->Fill(9.5,Five(h1,h2,h3,h4,h5).Im()/Five(0,0,0,0,0).Re(),Five(0,0,0,0,0).Re()); 
 } 
 if(TMath::Abs(Six(0,0,0,0,0,0).Re())>0.)
 {
  fNestedLoopsResultsCosPro->Fill(11.5,Six(h1,h2,h3,h4,h5,h6).Re()/Six(0,0,0,0,0,0).Re(),Six(0,0,0,0,0,0).Re()); 
  fNestedLoopsResultsSinPro->Fill(11.5,Six(h1,h2,h3,h4,h5,h6).Im()/Six(0,0,0,0,0,0).Re(),Six(0,0,0,0,0,0).Re()); 
 }
 if(TMath::Abs(Seven(0,0,0,0,0,0,0).Re())>0.)
 {
  fNestedLoopsResultsCosPro->Fill(13.5,Seven(h1,h2,h3,h4,h5,h6,h7).Re()/Seven(0,0,0,0,0,0,0).Re(),Seven(0,0,0,0,0,0,0).Re()); 
  fNestedLoopsResultsSinPro->Fill(13.5,Seven(h1,h2,h3,h4,h5,h6,h7).Im()/Seven(0,0,0,0,0,0,0).Re(),Seven(0,0,0,0,0,0,0).Re()); 
 }
 if(TMath::Abs(Eight(0,0,0,0,0,0,0,0).Re())>0.)
 {
  fNestedLoopsResultsCosPro->Fill(15.5,Eight(h1,h2,h3,h4,h5,h6,h7,h8).Re()/Eight(0,0,0,0,0,0,0,0).Re(),Eight(0,0,0,0,0,0,0,0).Re()); 
  fNestedLoopsResultsSinPro->Fill(15.5,Eight(h1,h2,h3,h4,h5,h6,h7,h8).Im()/Eight(0,0,0,0,0,0,0,0).Re(),Eight(0,0,0,0,0,0,0,0).Re()); 
 }

 Int_t nPrim = anEvent->NumberOfTracks(); 
 Double_t dMultRP = anEvent->GetNumberOfRPs(); // TBI shall I promote this variable into data member? 
 AliFlowTrackSimple *aftsTrack = NULL; 
 Double_t dPhi1=0.,dPhi2=0.,dPhi3=0.,dPhi4=0.,dPhi5=0.,dPhi6=0.,dPhi7=0.,dPhi8=0.; 
 Double_t wPhi1=1.,wPhi2=1.,wPhi3=1.,wPhi4=1.,wPhi5=1.,wPhi6=1.,wPhi7=1.,wPhi8=1.; 

 // 1-particle stuff: TBI       
 if(dMultRP>=1)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack = anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi1 = aftsTrack->Phi(); 
   if(fUseWeights[0][0]){wPhi1 = Weight(dPhi1,"RP","phi");}
   // Fill: 
   fNestedLoopsResultsCosPro->Fill(0.5,TMath::Cos(h1*dPhi1),wPhi1); 
   fNestedLoopsResultsSinPro->Fill(0.5,TMath::Sin(h1*dPhi1),wPhi1); 
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=1) 

 // 2-particle correlations:       
 if(dMultRP>=2)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack = anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi1 = aftsTrack->Phi(); 
   if(fUseWeights[0][0]){wPhi1 = Weight(dPhi1,"RP","phi");}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack = anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
    dPhi2 = aftsTrack->Phi();
    if(fUseWeights[0][0]){wPhi2 = Weight(dPhi2,"RP","phi");}
    // Fill:
    fNestedLoopsResultsCosPro->Fill(2.5,TMath::Cos(h1*dPhi1+h2*dPhi2),wPhi1*wPhi2); 
    fNestedLoopsResultsSinPro->Fill(2.5,TMath::Sin(h1*dPhi1+h2*dPhi2),wPhi1*wPhi2); 
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=2)

 // 3-particle correlations:         
 if(dMultRP>=3)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi1=aftsTrack->Phi();
   if(fUseWeights[0][0]){wPhi1 = Weight(dPhi1,"RP","phi");}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;} 
    if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
    dPhi2=aftsTrack->Phi();
    if(fUseWeights[0][0]){wPhi2 = Weight(dPhi2,"RP","phi");}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
     dPhi3=aftsTrack->Phi();
     if(fUseWeights[0][0]){wPhi3 = Weight(dPhi3,"RP","phi");}
     // Fill:
     fNestedLoopsResultsCosPro->Fill(4.5,TMath::Cos(h1*dPhi1+h2*dPhi2+h3*dPhi3),wPhi1*wPhi2*wPhi3);
     fNestedLoopsResultsSinPro->Fill(4.5,TMath::Sin(h1*dPhi1+h2*dPhi2+h3*dPhi3),wPhi1*wPhi2*wPhi3);
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=3)

 // 4-particle correlations:
 if(dMultRP>=4)
 {       
  for(Int_t i1=0;i1<nPrim;i1++)
  { 
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi1=aftsTrack->Phi();
   if(fUseWeights[0][0]){wPhi1 = Weight(dPhi1,"RP","phi");}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
    dPhi2=aftsTrack->Phi();
    if(fUseWeights[0][0]){wPhi2 = Weight(dPhi2,"RP","phi");}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
     dPhi3=aftsTrack->Phi();
     if(fUseWeights[0][0]){wPhi3 = Weight(dPhi3,"RP","phi");}
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3){continue;}
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())){continue;}
      if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
      dPhi4=aftsTrack->Phi();
      if(fUseWeights[0][0]){wPhi4 = Weight(dPhi4,"RP","phi");}
      // Fill:
      fNestedLoopsResultsCosPro->Fill(6.5,TMath::Cos(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4),wPhi1*wPhi2*wPhi3*wPhi4);
      fNestedLoopsResultsSinPro->Fill(6.5,TMath::Sin(h1*dPhi1+h2*dPhi2+h3*dPhi3+h4*dPhi4),wPhi1*wPhi2*wPhi3*wPhi4);
     } // end of for(Int_t i4=0;i4<nPrim;i4++) 
    } // end of for(Int_t i3=0;i3<nPrim;i3++)
   } // end of for(Int_t i2=0;i2<nPrim;i2++)
  } // end of for(Int_t i1=0;i1<nPrim;i1++)
 } // end of if(nPrim>=)

 // 5-particle correlations:      
 if(dMultRP>=5)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}  
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi1=aftsTrack->Phi();
   if(fUseWeights[0][0]){wPhi1 = Weight(dPhi1,"RP","phi");}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
    dPhi2=aftsTrack->Phi();
    if(fUseWeights[0][0]){wPhi2 = Weight(dPhi2,"RP","phi");}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
     dPhi3=aftsTrack->Phi();
     if(fUseWeights[0][0]){wPhi3 = Weight(dPhi3,"RP","phi");}
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3){continue;}
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())){continue;}
      if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
      dPhi4=aftsTrack->Phi();
      if(fUseWeights[0][0]){wPhi4 = Weight(dPhi4,"RP","phi");}
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())){continue;}
       if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
       dPhi5=aftsTrack->Phi();
       if(fUseWeights[0][0]){wPhi5 = Weight(dPhi5,"RP","phi");}
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
 if(dMultRP>=6)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi1=aftsTrack->Phi();
   if(fUseWeights[0][0]){wPhi1 = Weight(dPhi1,"RP","phi");}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
    dPhi2=aftsTrack->Phi();
    if(fUseWeights[0][0]){wPhi2 = Weight(dPhi2,"RP","phi");}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
     dPhi3=aftsTrack->Phi();
     if(fUseWeights[0][0]){wPhi3 = Weight(dPhi3,"RP","phi");}
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3){continue;}
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())){continue;}
      if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
      dPhi4=aftsTrack->Phi();
      if(fUseWeights[0][0]){wPhi4 = Weight(dPhi4,"RP","phi");}
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())){continue;}
       if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
       dPhi5=aftsTrack->Phi();
       if(fUseWeights[0][0]){wPhi5=Weight(dPhi5,"RP","phi");}
       for(Int_t i6=0;i6<nPrim;i6++)
       {
        if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
        aftsTrack=anEvent->GetTrack(i6);
        if(!(aftsTrack->InRPSelection())){continue;}
        if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
        dPhi6=aftsTrack->Phi(); 
        if(fUseWeights[0][0]){wPhi6=Weight(dPhi6,"RP","phi");}
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
 if(dMultRP>=7)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  { 
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi1=aftsTrack->Phi();
   if(fUseWeights[0][0]){wPhi1=Weight(dPhi1,"RP","phi");}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
    dPhi2=aftsTrack->Phi();
    if(fUseWeights[0][0]){wPhi2=Weight(dPhi2,"RP","phi");}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
     dPhi3=aftsTrack->Phi();
     if(fUseWeights[0][0]){wPhi3=Weight(dPhi3,"RP","phi");}
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3){continue;}
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())){continue;}
      if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
      dPhi4=aftsTrack->Phi();
      if(fUseWeights[0][0]){wPhi4=Weight(dPhi4,"RP","phi");}
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())){continue;}
       if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
       dPhi5=aftsTrack->Phi();
       if(fUseWeights[0][0]){wPhi5=Weight(dPhi5,"RP","phi");}
       for(Int_t i6=0;i6<nPrim;i6++)
       {
        if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
        aftsTrack=anEvent->GetTrack(i6);
        if(!(aftsTrack->InRPSelection())){continue;}
        if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
        dPhi6=aftsTrack->Phi(); 
        if(fUseWeights[0][0]){wPhi6=Weight(dPhi6,"RP","phi");}
        for(Int_t i7=0;i7<nPrim;i7++)
        {
         if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6){continue;}
         aftsTrack=anEvent->GetTrack(i7);
         if(!(aftsTrack->InRPSelection())){continue;}
         if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
         dPhi7=aftsTrack->Phi(); 
         if(fUseWeights[0][0]){wPhi7=Weight(dPhi7,"RP","phi");}
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
 if(dMultRP>=8)
 {
  for(Int_t i1=0;i1<nPrim;i1++)
  {
   aftsTrack=anEvent->GetTrack(i1);
   if(!(aftsTrack->InRPSelection())){continue;}
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi1=aftsTrack->Phi();
   if(fUseWeights[0][0]){wPhi1=Weight(dPhi1,"RP","phi");}
   for(Int_t i2=0;i2<nPrim;i2++)
   {
    if(i2==i1){continue;}
    aftsTrack=anEvent->GetTrack(i2);
    if(!(aftsTrack->InRPSelection())){continue;}
    if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
    dPhi2=aftsTrack->Phi();
    if(fUseWeights[0][0]){wPhi2=Weight(dPhi2,"RP","phi");}
    for(Int_t i3=0;i3<nPrim;i3++)
    {
     if(i3==i1||i3==i2){continue;}
     aftsTrack=anEvent->GetTrack(i3);
     if(!(aftsTrack->InRPSelection())){continue;}
     if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
     dPhi3=aftsTrack->Phi();
     if(fUseWeights[0][0]){wPhi3=Weight(dPhi3,"RP","phi");}
     for(Int_t i4=0;i4<nPrim;i4++)
     {
      if(i4==i1||i4==i2||i4==i3){continue;}
      aftsTrack=anEvent->GetTrack(i4);
      if(!(aftsTrack->InRPSelection())){continue;}
      if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
      dPhi4=aftsTrack->Phi();
      if(fUseWeights[0][0]){wPhi4=Weight(dPhi4,"RP","phi");}
      for(Int_t i5=0;i5<nPrim;i5++)
      {
       if(i5==i1||i5==i2||i5==i3||i5==i4){continue;}
       aftsTrack=anEvent->GetTrack(i5);
       if(!(aftsTrack->InRPSelection())){continue;}
       if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
       dPhi5=aftsTrack->Phi();
       if(fUseWeights[0][0]){wPhi5=Weight(dPhi5,"RP","phi");}
       for(Int_t i6=0;i6<nPrim;i6++)
       {
        if(i6==i1||i6==i2||i6==i3||i6==i4||i6==i5){continue;}
        aftsTrack=anEvent->GetTrack(i6);
        if(!(aftsTrack->InRPSelection())){continue;}
        if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
        dPhi6=aftsTrack->Phi();
        if(fUseWeights[0][0]){wPhi6=Weight(dPhi6,"RP","phi");}
        for(Int_t i7=0;i7<nPrim;i7++)
        {
         if(i7==i1||i7==i2||i7==i3||i7==i4||i7==i5||i7==i6){continue;}
         aftsTrack=anEvent->GetTrack(i7);
         if(!(aftsTrack->InRPSelection())){continue;}
         if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
         dPhi7=aftsTrack->Phi();
         if(fUseWeights[0][0]){wPhi7=Weight(dPhi7,"RP","phi");}
         for(Int_t i8=0;i8<nPrim;i8++)
         {
          if(i8==i1||i8==i2||i8==i3||i8==i4||i8==i5||i8==i6||i8==i7){continue;}
          aftsTrack=anEvent->GetTrack(i8);
          if(!(aftsTrack->InRPSelection())){continue;}
          if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
          dPhi8=aftsTrack->Phi();
          if(fUseWeights[0][0]){wPhi8=Weight(dPhi8,"RP","phi");}
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

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckDiffWithNestedLoops(AliFlowEventSimple *anEvent)
{
 // Cross-check results for differential multi-particle correlations with nested loops.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckDiffWithNestedLoops(AliFlowEventSimple *anEvent)";

 Int_t nPrim = anEvent->NumberOfTracks(); 
 AliFlowTrackSimple *aftsTrack = NULL; 
 Double_t dPsi1=0.,dPhi2=0.,dPhi3=0.,dPhi4=0.; 
 Double_t wPsi1=1.,wPhi2=1.,wPhi3=1.,wPhi4=1.; 

 Int_t cs = fCrossCheckDiffCSCOBN[0]; // cos/sin

 // TBI reimplement lines below in a more civilised manner:
 Bool_t bCrossCheck2p = kFALSE;
 Bool_t bCrossCheck3p = kFALSE;
 Bool_t bCrossCheck4p = kFALSE;

 if(fCrossCheckDiffCSCOBN[1] == 2){bCrossCheck2p = kTRUE;}
 else if(fCrossCheckDiffCSCOBN[1] == 3){bCrossCheck3p = kTRUE;}
 else if(fCrossCheckDiffCSCOBN[1] == 4){bCrossCheck4p = kTRUE;}

 if(Int_t(bCrossCheck2p + bCrossCheck3p + bCrossCheck4p) > 1)
 {
  Fatal(sMethodName.Data(),"Int_t(bCrossCheck2p + bCrossCheck3p + bCrossCheck4p) > 1");
 }
 if(!(bCrossCheck2p || bCrossCheck3p || bCrossCheck4p))
 {
  Fatal(sMethodName.Data(),"!(bCrossCheck2p || bCrossCheck3p || bCrossCheck4p)");
 }
 Int_t nDiffBinNo = fCrossCheckDiffCSCOBN[2];
 Double_t dPt = 0., dEta = 0.;

 // <2'>:
 for(Int_t i1=0;i1<nPrim;i1++) // Loop over particles in a differential bin 
 {
  aftsTrack=anEvent->GetTrack(i1);
  if(!(aftsTrack->InPOISelection())){continue;}
  if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
  dPsi1=aftsTrack->Phi();
  if(fCalculateDiffCorrelationsVsPt)
  {
   dPt=aftsTrack->Pt();
   if(fDiffCorrelationsPro[0][1]->FindBin(dPt) != nDiffBinNo){continue;} // TBI spaghetti again 
  } else 
    {
     dEta=aftsTrack->Eta();
     if(fDiffCorrelationsPro[0][1]->FindBin(dEta) != nDiffBinNo){continue;} // TBI spaghetti again 
    }
  if(fUseWeights[1][0]){wPsi1=Weight(dPsi1,"POI","phi");}
  for(Int_t i2=0;i2<nPrim;i2++) // Loop over particles in an event
  {
   if(i2==i1){continue;} // get rid of autocorrelations
   aftsTrack=anEvent->GetTrack(i2);
   if(!(aftsTrack->InRPSelection())){continue;}
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi2=aftsTrack->Phi();
   if(fUseWeights[0][0]){wPhi2=Weight(dPhi2,"RP","phi");}
   // Fill profiles:
   if(bCrossCheck2p)
   {
    if(fCrossCheckDiffCSCOBN[0] == 0)
    {
     fNestedLoopsDiffResultsPro->Fill(0.5,TMath::Cos(fDiffHarmonics[1][0]*dPsi1+fDiffHarmonics[1][1]*dPhi2),wPsi1*wPhi2);
    } else {fNestedLoopsDiffResultsPro->Fill(0.5,TMath::Sin(fDiffHarmonics[1][0]*dPsi1+fDiffHarmonics[1][1]*dPhi2),wPsi1*wPhi2);}
   } // if(bCrossCheck2p) 
  } // for(Int_t i2=0;i2<nPrim;i2++)
 } // for(Int_t i1=0;i1<nPrim;i1++)

 // <3'>:
 for(Int_t i1=0;i1<nPrim;i1++) // Loop over particles in a differential bin
 {
  aftsTrack=anEvent->GetTrack(i1);
  if(!(aftsTrack->InPOISelection())){continue;}
  if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
  dPsi1=aftsTrack->Phi();
  if(fCalculateDiffCorrelationsVsPt)
  {
   dPt=aftsTrack->Pt();
   if(fDiffCorrelationsPro[0][1]->FindBin(dPt) != nDiffBinNo){continue;} // TBI spaghetti again 
  } else 
    {
     dEta=aftsTrack->Eta();
     if(fDiffCorrelationsPro[0][1]->FindBin(dEta) != nDiffBinNo){continue;} // TBI spaghetti again 
    }
  if(fUseWeights[1][0]){wPsi1=Weight(dPsi1,"POI","phi");}
  for(Int_t i2=0;i2<nPrim;i2++) // Loop over particles in an event
  {
   if(i2==i1){continue;} // get rid of autocorrelations
   aftsTrack=anEvent->GetTrack(i2);
   if(!(aftsTrack->InRPSelection())){continue;}
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi2=aftsTrack->Phi();
   if(fUseWeights[0][0]){wPhi2=Weight(dPhi2,"RP","phi");}
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2){continue;} // get rid of autocorrelations
    aftsTrack=anEvent->GetTrack(i3);
    if(!(aftsTrack->InRPSelection())){continue;}
    if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
    dPhi3=aftsTrack->Phi();
    if(fUseWeights[0][0]){wPhi3=Weight(dPhi3,"RP","phi");}
    // Fill the profiles:
    if(bCrossCheck3p)
    {
     if(fCrossCheckDiffCSCOBN[0] == 0)
     {
      fNestedLoopsDiffResultsPro->Fill(0.5,TMath::Cos(fDiffHarmonics[2][0]*dPsi1+fDiffHarmonics[2][1]*dPhi2+fDiffHarmonics[2][2]*dPhi3),wPsi1*wPhi2*wPhi3);  
     } else {fNestedLoopsDiffResultsPro->Fill(0.5,TMath::Sin(fDiffHarmonics[2][0]*dPsi1+fDiffHarmonics[2][1]*dPhi2+fDiffHarmonics[2][2]*dPhi3),wPsi1*wPhi2*wPhi3);}
    } // if(bCrossCheck3p)
   } // end of for(Int_t i3=0;i3<nPrim;i3++)  
  } // for(Int_t i2=0;i2<nPrim;i2++)
 } // for(Int_t i1=0;i1<nPrim;i1++)

 // <4'>:
 for(Int_t i1=0;i1<nPrim;i1++) // Loop over particles in a differential bin
 {
  aftsTrack=anEvent->GetTrack(i1);
  if(!(aftsTrack->InPOISelection())){continue;}
  if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
  dPsi1=aftsTrack->Phi();
  if(fCalculateDiffCorrelationsVsPt)
  {
   dPt=aftsTrack->Pt();
   if(fDiffCorrelationsPro[0][1]->FindBin(dPt) != nDiffBinNo){continue;} // TBI spaghetti again 
  } else 
    {
     dEta=aftsTrack->Eta();
     if(fDiffCorrelationsPro[0][1]->FindBin(dEta) != nDiffBinNo){continue;} // TBI spaghetti again 
    }
  if(fUseWeights[1][0]){wPsi1=Weight(dPsi1,"POI","phi");}
  for(Int_t i2=0;i2<nPrim;i2++) // Loop over particles in an event
  {
   if(i2==i1){continue;} // get rid of autocorrelations
   aftsTrack=anEvent->GetTrack(i2);
   if(!(aftsTrack->InRPSelection())){continue;}
   if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
   dPhi2=aftsTrack->Phi();
   if(fUseWeights[0][0]){wPhi2=Weight(dPhi2,"RP","phi");}
   for(Int_t i3=0;i3<nPrim;i3++)
   {
    if(i3==i1||i3==i2){continue;} // get rid of autocorrelations
    aftsTrack=anEvent->GetTrack(i3);
    if(!(aftsTrack->InRPSelection())){continue;}
    if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
    dPhi3=aftsTrack->Phi();
    if(fUseWeights[0][0]){wPhi3=Weight(dPhi3,"RP","phi");}
    for(Int_t i4=0;i4<nPrim;i4++)
    {
     if(i4==i1||i4==i2||i4==i3){continue;} // get rid of autocorrelations
     aftsTrack=anEvent->GetTrack(i4);
     if(!(aftsTrack->InRPSelection())){continue;}
     if(fSkipSomeIntervals && !TrackIsInSpecifiedIntervals(aftsTrack)){continue;} // TBI tmp gym
     dPhi4=aftsTrack->Phi();
     if(fUseWeights[0][0]){wPhi4=Weight(dPhi4,"RP","phi");}
     // Fill the profiles:
     if(bCrossCheck4p)
     {
      if(fCrossCheckDiffCSCOBN[0] == 0)
      {
       fNestedLoopsDiffResultsPro->Fill(0.5,TMath::Cos(fDiffHarmonics[3][0]*dPsi1+fDiffHarmonics[3][1]*dPhi2+fDiffHarmonics[3][2]*dPhi3+fDiffHarmonics[3][3]*dPhi4),wPsi1*wPhi2*wPhi3*wPhi4);
      } else {fNestedLoopsDiffResultsPro->Fill(0.5,TMath::Sin(fDiffHarmonics[3][0]*dPsi1+fDiffHarmonics[3][1]*dPhi2+fDiffHarmonics[3][2]*dPhi3+fDiffHarmonics[3][3]*dPhi4),wPsi1*wPhi2*wPhi3*wPhi4);} 
     } // if(bCrossCheck4p)
    } // end of for(Int_t i4=0;i4<nPrim;i4++) 
   } // end of for(Int_t i3=0;i3<nPrim;i3++)  
  } // for(Int_t i2=0;i2<nPrim;i2++)
 } // for(Int_t i1=0;i1<nPrim;i1++)

 // Printout:
 // 2-p:
 if(bCrossCheck2p)
 {
  printf("\n  2-p => Q-vector:     %.12f",fDiffCorrelationsPro[cs][1]->GetBinContent(nDiffBinNo));
  printf("\n  2-p => Nested loops: %.12f\n",fNestedLoopsDiffResultsPro->GetBinContent(1));
 }
 // 3-p:
 if(bCrossCheck3p)
 {
  printf("\n  3-p => Q-vector:     %.12f",fDiffCorrelationsPro[cs][2]->GetBinContent(nDiffBinNo));
  printf("\n  3-p => Nested loops: %.12f\n",fNestedLoopsDiffResultsPro->GetBinContent(1));
 } 
 // 4-p:
 if(bCrossCheck4p)
 {
  printf("\n  4-p => Q-vector:     %.12f",fDiffCorrelationsPro[cs][3]->GetBinContent(nDiffBinNo));
  printf("\n  4-p => Nested loops: %.12f\n",fNestedLoopsDiffResultsPro->GetBinContent(1));
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckDiffWithNestedLoops(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::FillQvector(AliFlowEventSimple *anEvent)
{
 // Fill Q-vector components.

 Int_t nTracks = anEvent->NumberOfTracks(); // TBI shall I promote this to data member?
 Double_t dPhi = 0., wPhi = 1.; // azimuthal angle and corresponding phi weight
 Double_t dPt = 0., wPt = 1.; // transverse momentum and corresponding pT weight
 Double_t dEta = 0., wEta = 1.; // pseudorapidity and corresponding eta weight
 Double_t wToPowerP = 1.; // weight raised to power p
 Int_t nCounterRPs = 0;
 for(Int_t t=0;t<nTracks;t++) // loop over all tracks
 {
  AliFlowTrackSimple *pTrack = NULL;
  if(!fSelectRandomlyRPs) // TBI hw RPs
  {
   pTrack = anEvent->GetTrack(t);
  }
  else
  {
   pTrack = anEvent->GetTrack((Int_t)fRandomIndicesRPs->GetAt(t));
  }

  if(!pTrack){printf("\n Error: pTrack is NULL in MPC::FillQvector(...) !!!!"); continue;}

  if(!TrackIsInSpecifiedIntervals(pTrack)){continue;} // TBI tmp gym

  if(!(pTrack->InRPSelection() || pTrack->InPOISelection())){printf("\n Error: pTrack is neither RP nor POI !!!!"); continue;}

  if(pTrack->InRPSelection()) // fill Q-vector components only with reference particles
  {
   nCounterRPs++;
   if(fSelectRandomlyRPs && nCounterRPs == fnSelectedRandomlyRPs){break;} // for(Int_t t=0;t<nTracks;t++) // loop over all tracks

   wPhi = 1.; wPt = 1.; wEta = 1.; wToPowerP = 1.; // TBI this shall go somewhere else, for performance sake

   // Access kinematic variables for RP and corresponding weights:
   dPhi = pTrack->Phi(); // azimuthal angle
   if(fUseWeights[0][0]){wPhi = Weight(dPhi,"RP","phi");} // corresponding phi weight
   //if(dPhi < 0.){dPhi += TMath::TwoPi();} TBI
   //if(dPhi > TMath::TwoPi()){dPhi -= TMath::TwoPi();} TBI
   dPt = pTrack->Pt();
   if(fUseWeights[0][1]){wPt = Weight(dPt,"RP","pt");} // corresponding pT weight
   dEta = pTrack->Eta();
   if(fUseWeights[0][2]){wEta = Weight(dEta,"RP","eta");} // corresponding eta weight

   // Calculate Q-vector components:
   for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)
   {
    for(Int_t wp=0;wp<fMaxCorrelator+1;wp++) // weight power
    {
     if(fUseWeights[0][0]||fUseWeights[0][1]||fUseWeights[0][2]){wToPowerP = pow(wPhi*wPt*wEta,wp);} 
     fQvector[h][wp] += TComplex(wToPowerP*TMath::Cos(h*dPhi),wToPowerP*TMath::Sin(h*dPhi));
    } // for(Int_t wp=0;wp<fMaxCorrelator+1;wp++)
   } // for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)
  } // if(pTrack->InRPSelection()) // fill Q-vector components only with reference particles

  // Differential Q-vectors (a.k.a. p-vector and q-vector):
  if(!fCalculateDiffQvectors){continue;}
  if(pTrack->InPOISelection()) 
  {
   wPhi = 1.; wPt = 1.; wEta = 1.; wToPowerP = 1.; // TBI this shall go somewhere else, for performance sake

   // Access kinematic variables for POI and corresponding weights:
   dPhi = pTrack->Phi(); // azimuthal angle
   if(fUseWeights[1][0]){wPhi = Weight(dPhi,"POI","phi");} // corresponding phi weight
   //if(dPhi < 0.){dPhi += TMath::TwoPi();} TBI
   //if(dPhi > TMath::TwoPi()){dPhi -= TMath::TwoPi();} TBI
   dPt = pTrack->Pt();
   if(fUseWeights[1][1]){wPt = Weight(dPt,"POI","pt");} // corresponding pT weight
   dEta = pTrack->Eta();
   if(fUseWeights[1][2]){wEta = Weight(dEta,"POI","eta");} // corresponding eta weight

   // Determine bin:
   Int_t binNo = -44;
   if(fCalculateDiffCorrelationsVsPt)
   { 
    binNo = fDiffCorrelationsPro[0][0]->FindBin(dPt); // TBI: hardwired [0][0]
   } else
     {
      binNo = fDiffCorrelationsPro[0][0]->FindBin(dEta); // TBI: hardwired [0][0]
     }
   // Calculate p-vector components:
   for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)
   {
    for(Int_t wp=0;wp<fMaxCorrelator+1;wp++) // weight power
    {
     if(fUseWeights[1][0]||fUseWeights[1][1]||fUseWeights[1][2]){wToPowerP = pow(wPhi*wPt*wEta,wp);} 
     fpvector[binNo-1][h][wp] += TComplex(wToPowerP*TMath::Cos(h*dPhi),wToPowerP*TMath::Sin(h*dPhi));

     if(pTrack->InRPSelection()) 
     {
      // Fill q-vector components:
      wPhi = 1.; wPt = 1.; wEta = 1.; wToPowerP = 1.; // TBI this shall go somewhere else, for performance sake

      if(fUseWeights[0][0]){wPhi = Weight(dPhi,"RP","phi");} // corresponding phi weight
      //if(dPhi < 0.){dPhi += TMath::TwoPi();} TBI
      //if(dPhi > TMath::TwoPi()){dPhi -= TMath::TwoPi();} TBI
      if(fUseWeights[0][1]){wPt = Weight(dPt,"RP","pt");} // corresponding pT weight
      if(fUseWeights[0][2]){wEta = Weight(dEta,"RP","eta");} // corresponding eta weight
      if(fUseWeights[1][0]){wPhi = Weight(dPhi,"POI","phi");} // corresponding phi weight
      //if(dPhi < 0.){dPhi += TMath::TwoPi();} TBI
      //if(dPhi > TMath::TwoPi()){dPhi -= TMath::TwoPi();} TBI
      if(fUseWeights[1][1]){wPt = Weight(dPt,"POI","pt");} // corresponding pT weight
      if(fUseWeights[1][2]){wEta = Weight(dEta,"POI","eta");} // corresponding eta weight
      if(fUseWeights[0][0]||fUseWeights[0][1]||fUseWeights[0][2]||fUseWeights[1][0]||fUseWeights[1][1]||fUseWeights[1][2]){wToPowerP = pow(wPhi*wPt*wEta,wp);} 
      fqvector[binNo-1][h][wp] += TComplex(wToPowerP*TMath::Cos(h*dPhi),wToPowerP*TMath::Sin(h*dPhi));
     } // if(pTrack->InRPSelection()) 

    } // for(Int_t wp=0;wp<fMaxCorrelator+1;wp++)
   } // for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)
  } // if(pTrack->InPOISelection()) 

 } // for(Int_t t=0;t<nTracks;t++) // loop over all tracks

} // void AliFlowAnalysisWithMultiparticleCorrelations::FillQvector(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckSettings()
{
 // Cross-check all initial settings in this method. 
 
 // a) Few cross-checks for control histograms;
 // b) Few cross-checks for flags for correlations;
 // c) 'Standard candles';
 // d) Q-cumulants;
 // e) Weights;
 // f) Differential correlations;
 // g) Nested loops;
 // h) Dump the points.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckSettings()";

 // a) Few cross-checks for control histograms: TBI the lines below are not really what they are supposed to be...
 /*
 if(fFillKinematicsHist && !fFillControlHistograms){Fatal(sMethodName.Data(),"fFillKinematicsHist && !fFillControlHistograms");}
 if(fFillMultDistributionsHist && !fFillControlHistograms){Fatal(sMethodName.Data(),"fFillMultDistributionsHist && !fFillControlHistograms");}
 if(fFillMultCorrelationsHist && !fFillControlHistograms){Fatal(sMethodName.Data(),"fFillMultCorrelationsHist && !fFillControlHistograms");}
 */ 

 // b) Few cross-checks for flags for correlations: // TBI the lines bellow can be civilized
 Int_t iSum = (Int_t)fCalculateIsotropic + (Int_t)fCalculateSame + (Int_t)fCalculateSameIsotropic;
 if(iSum>1){Fatal(sMethodName.Data(),"iSum is doing crazy things...");}
 if(fCalculateOnlyCos && fCalculateOnlySin){Fatal(sMethodName.Data(),"fCalculateOnlyCos && fCalculateOnlySin");}

 // c) 'Standard candles':
 if(fCalculateStandardCandles && !fCalculateCorrelations)
 {
  Fatal(sMethodName.Data(),"fCalculateStandardCandles && !fCalculateCorrelations");
 }
 if(fCalculateStandardCandles && fCalculateCorrelations && fCalculateSameIsotropic)
 {
  Fatal(sMethodName.Data(),"fCalculateStandardCandles && fCalculateCorrelations && fCalculateSameIsotropic");
 }
 if(fCalculateStandardCandles && fCalculateOnlyForHarmonicQC)
 {
  Fatal(sMethodName.Data(),"fCalculateStandardCandles && fCalculateOnlyForHarmonicQC");
 }
 if(fCalculateStandardCandles && fCalculateOnlyForSC && (4!=fDontGoBeyond))
 {
  Fatal(sMethodName.Data(),"fCalculateStandardCandles && fCalculateOnlyForSC && (4!=fDontGoBeyond)");
 }
 if(fCalculateStandardCandles && !fPropagateErrorSC)
 {
  Warning(sMethodName.Data(),"fCalculateStandardCandles && !fPropagateErrorSC");
 }
 if(fCalculateStandardCandles && fCalculateOnlySin)
 {
  Fatal(sMethodName.Data(),"fCalculateStandardCandles && fCalculateOnlySin");
 }
 if(fCalculateStandardCandles && fDontGoBeyond < 3)
 {
  Fatal(sMethodName.Data(),"fCalculateStandardCandles && fDontGoBeyond < 3");
 }

 // d) Q-cumulants:
 if(fCalculateQcumulants && !fCalculateCorrelations)
 {
  Fatal(sMethodName.Data(),"fCalculateQcumulants && !fCalculateCorrelations");
 }
 if(fCalculateQcumulants && !(fHarmonicQC > 0))
 {
  Fatal(sMethodName.Data(),"fCalculateQcumulants && !(fHarmonicQC > 0)");
 }
 if(fCalculateQcumulants && fCalculateOnlyForSC)
 {
  Fatal(sMethodName.Data(),"fCalculateQcumulants && fCalculateOnlyForSC");
 }
 if(fCalculateQcumulants && !fPropagateErrorQC)
 {
  Warning(sMethodName.Data(),"fCalculateQcumulants && !fPropagateErrorQC");
 }
 if(fCalculateQcumulants && fCalculateOnlySin)
 {
  Fatal(sMethodName.Data(),"fCalculateQcumulants && fCalculateOnlySin");
 }
 
 // e) Weights:
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   if(fUseWeights[rp][ppe] && !fWeightsHist[rp][ppe])
   {
    Fatal(sMethodName.Data(),"fUseWeights[rp][ppe] && !fWeightsHist[rp][ppe], rp = %d, ppe = %d",rp,ppe);
   }
  }
 }

 // f) Differential correlations:
 if(fCalculateDiffCorrelations && !fCalculateDiffQvectors)
 {
  Fatal(sMethodName.Data(),"fCalculateDiffCorrelations && !fCalculateDiffQvectors"); 
 }
 if(fCalculateDiffCorrelations && !fCalculateQvector)
 {
  Fatal(sMethodName.Data(),"fCalculateDiffCorrelations && !fCalculateQvector"); 
 }
 if(!fCalculateDiffCorrelations && fCalculateDiffQvectors)
 {
  Fatal(sMethodName.Data(),"!fCalculateDiffCorrelations && fCalculateDiffQvectors"); 
 }
 if(fCalculateDiffCorrelations && !fUseDefaultBinning && (fnDiffBins < 1 || !fRangesDiffBins))
 {
  Fatal(sMethodName.Data(),"fCalculateDiffCorrelations && !fUseDefaultBinning && (fnDiffBins < 1 || !fRangesDiffBins)"); 
 }
 if(fCalculateDiffCorrelations && !(fCalculateDiffCos || fCalculateDiffSin))
 {
  Fatal(sMethodName.Data(),"fCalculateDiffCorrelations && !(fCalculateDiffCos || fCalculateDiffSin)"); 
 }
 if(fCalculateDiffCorrelations && fDontFill[1])
 {
  Warning(sMethodName.Data(),"fCalculateDiffCorrelations && fDontFill[1]"); 
 }

 // g) Nested loops:
 if(fCrossCheckDiffWithNestedLoops && (1 == fCrossCheckDiffCSCOBN[0] && !fCalculateDiffSin))
 {
  Fatal(sMethodName.Data(),"fCrossCheckDiffWithNestedLoops && (1 == fCrossCheckDiffCSCOBN[0] && !CalculateDiffSin)"); 
 }
 if(fCrossCheckDiffWithNestedLoops && (0 == fCrossCheckDiffCSCOBN[0] && !fCalculateDiffCos))
 {
  Fatal(sMethodName.Data(),"fCrossCheckDiffWithNestedLoops && (0 == fCrossCheckDiffCSCOBN[0] && !CalculateDiffCos)"); 
 }

 // h) Dump the points:
 if(fDumpThePoints && !fFillMultDistributionsHist)
 {
  Fatal(sMethodName.Data(),"if(fDumpThePoints && !fFillMultDistributionsHist)"); 
 }
 if(fDumpThePoints && fMaxNoEventsPerFile <= 0)
 {
  Fatal(sMethodName.Data(),"if(fDumpThePoints && fMaxNoEventsPerFile <= 0)"); 
 }

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckSettings()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for Q-vectors;
 // c) Book and nest lists for correlations;
 // d) Book and nest lists for e-b-e cumulants;
 // e) Book and nest lists for weights;
 // f) Book and nest lists for nested loops;
 // g) Book and nest lists for 'standard candles';
 // h) Book and nest lists for Q-cumulants;
 // i) Book and nest lists for differential correlations;
 // j) Book and nest lists for symmetry plane correlations;
 // k) Book and nest lists for correlations with eta gaps.
 
 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("Control Histograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);

 // b) Book and nest lists for Q-vectors:
 fQvectorList = new TList();
 fQvectorList->SetName("Q-vectors");
 fQvectorList->SetOwner(kTRUE);
 fHistList->Add(fQvectorList);

 // c) Book and nest lists for correlations:
 fCorrelationsList = new TList();
 fCorrelationsList->SetName("Correlations");
 fCorrelationsList->SetOwner(kTRUE);
 fHistList->Add(fCorrelationsList);

 // d) Book and nest lists for e-b-e cumulants:
 fEbECumulantsList = new TList();
 fEbECumulantsList->SetName("E-b-e Cumulants");
 fEbECumulantsList->SetOwner(kTRUE);
 fHistList->Add(fEbECumulantsList);

 // e) Book and nest lists for weights:
 fWeightsList = new TList();
 fWeightsList->SetName("Weights");
 fWeightsList->SetOwner(kTRUE);
 fHistList->Add(fWeightsList);

 // f) Book and nest lists for nested loops:
 fNestedLoopsList = new TList();
 fNestedLoopsList->SetName("Nested Loops");
 fNestedLoopsList->SetOwner(kTRUE);
 fHistList->Add(fNestedLoopsList);

 // g) Book and nest lists for 'standard candles':
 fStandardCandlesList = new TList();
 fStandardCandlesList->SetName("Standard Candles");
 fStandardCandlesList->SetOwner(kTRUE);
 fHistList->Add(fStandardCandlesList);

 // h) Book and nest lists for Q-cumulants:
 fQcumulantsList = new TList();
 fQcumulantsList->SetName("Q-cumulants");
 fQcumulantsList->SetOwner(kTRUE);
 fHistList->Add(fQcumulantsList);

 // i) Book and nest lists for differential correlations:
 fDiffCorrelationsList = new TList();
 fDiffCorrelationsList->SetName("Differential Correlations");
 fDiffCorrelationsList->SetOwner(kTRUE);
 fHistList->Add(fDiffCorrelationsList);

 // j) Book and nest lists for symmetry plane correlations:
 fSymmetryPlanesList = new TList();
 fSymmetryPlanesList->SetName("Symmetry_Plane_Correlations");
 fSymmetryPlanesList->SetOwner(kTRUE);
 fHistList->Add(fSymmetryPlanesList);

 // k) Book and nest lists for correlations with eta gaps:
 fEtaGapsList = new TList();
 fEtaGapsList->SetName("Correlations_with_eta_gaps");
 fEtaGapsList->SetOwner(kTRUE);
 fHistList->Add(fEtaGapsList);

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
 fControlHistogramsFlagsPro = new TProfile("fControlHistogramsFlagsPro","Flags and settings for control histograms",7,0,7);
 fControlHistogramsFlagsPro->SetTickLength(-0.01,"Y");
 fControlHistogramsFlagsPro->SetMarkerStyle(25);
 fControlHistogramsFlagsPro->SetLabelSize(0.04);
 fControlHistogramsFlagsPro->SetLabelOffset(0.02,"Y");
 fControlHistogramsFlagsPro->SetStats(kFALSE);
 fControlHistogramsFlagsPro->SetFillColor(kGray);
 fControlHistogramsFlagsPro->SetLineColor(kBlack);
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(1,"fFillControlHistograms"); fControlHistogramsFlagsPro->Fill(0.5,fFillControlHistograms);
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(2,"fFillKinematicsHist"); fControlHistogramsFlagsPro->Fill(1.5,fFillKinematicsHist);
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(3,"fFillMultDistributionsHist"); fControlHistogramsFlagsPro->Fill(2.5,fFillMultDistributionsHist);
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(4,"fFillMultCorrelationsHist"); fControlHistogramsFlagsPro->Fill(3.5,fFillMultCorrelationsHist);
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(5,"fDontFill[0=RP]"); fControlHistogramsFlagsPro->Fill(4.5,fDontFill[0]);
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(6,"fDontFill[1=POI]"); fControlHistogramsFlagsPro->Fill(5.5,fDontFill[1]);
 fControlHistogramsFlagsPro->GetXaxis()->SetBinLabel(7,"fDontFill[2=REF]"); fControlHistogramsFlagsPro->Fill(6.5,fDontFill[2]);
 fControlHistogramsList->Add(fControlHistogramsFlagsPro);

 if(!fFillControlHistograms){return;} // TBI is this safe? Well, perhaps it is if I can't implement it better...

 // b) Book all control histograms: // TBI add setters for all these values
 //  b0) Book TH1D *fKinematicsHist[2][3]:
 TString name[2][3] = {{"RP,phi","RP,pt","RP,eta"},{"POI,phi","POI,pt","POI,eta"}}; // [RP,POI][phi,pt,eta]
 TString title[2] = {"Reference particles (RP)","Particles of interest (POI)"}; // [RP,POI]
 Int_t lineColor[2] = {kBlue,kRed}; // [RP,POI]
 Int_t fillColor[2] = {kBlue-10,kRed-10}; // [RP,POI]
 TString xAxisTitle[3] = {"#phi","p_{T}","#eta"}; // [phi,pt,eta]
 if(fFillKinematicsHist)
 {
  for(Int_t rp=0;rp<2;rp++) // [RP,POI]
  {
   if(fDontFill[rp]){continue;}
   for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
   {
    fKinematicsHist[rp][ppe] = new TH1D(name[rp][ppe].Data(),title[rp].Data(),fnBins[rp][ppe],fMin[rp][ppe],fMax[rp][ppe]);
    fKinematicsHist[rp][ppe]->GetXaxis()->SetTitle(xAxisTitle[ppe].Data());
    fKinematicsHist[rp][ppe]->SetLineColor(lineColor[rp]);
    fKinematicsHist[rp][ppe]->SetFillColor(fillColor[rp]);
    fKinematicsHist[rp][ppe]->SetMinimum(0.); 
    fControlHistogramsList->Add(fKinematicsHist[rp][ppe]);
   }
  }
 } // if(fFillKinematicsHist)

 //  b1) Book TH1D *fMultDistributionsHist[3]: // TBI add setters for all these values
 TString nameMult[3] = {"Multiplicity (RP)","Multiplicity (POI)","Multiplicity (REF)"}; // [RP,POI,reference multiplicity]
 TString titleMult[3] = {"Reference particles (RP)","Particles of interest (POI)",""}; // [RP,POI,reference multiplicity]
 Int_t lineColorMult[3] = {kBlue,kRed,kGreen+2}; // [RP,POI,reference multiplicity]
 Int_t fillColorMult[3] = {kBlue-10,kRed-10,kGreen-10}; // [RP,POI,reference multiplicity]
 TString xAxisTitleMult[3] = {"Multiplicity (RP)","Multiplicity (POI)","Multiplicity (REF)"}; // [phi,pt,eta]
 if(fFillMultDistributionsHist)
 {
  for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
  {
   if(fDontFill[rprm]){continue;}
   fMultDistributionsHist[rprm] = new TH1D(nameMult[rprm].Data(),titleMult[rprm].Data(),fnBinsMult[rprm],fMinMult[rprm],fMaxMult[rprm]);
   fMultDistributionsHist[rprm]->GetXaxis()->SetTitle(xAxisTitleMult[rprm].Data());
   fMultDistributionsHist[rprm]->SetLineColor(lineColorMult[rprm]);
   fMultDistributionsHist[rprm]->SetFillColor(fillColorMult[rprm]);
   fControlHistogramsList->Add(fMultDistributionsHist[rprm]);
  } // for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 } // if(fFillMultDistributionsHist)

 //  b2) Book TH2I *fMultCorrelationsHist[3]: 
 if(fFillMultCorrelationsHist)
 {
  if(!fDontFill[0] && !fDontFill[1])
  {
   fMultCorrelationsHist[0] = new TH2I("Multiplicity (RP vs. POI)","Multiplicity (RP vs. POI)",fnBinsMult[0],fMinMult[0],fMaxMult[0],fnBinsMult[1],fMinMult[1],fMaxMult[1]);
   fMultCorrelationsHist[0]->GetXaxis()->SetTitle(xAxisTitleMult[0].Data());
   fMultCorrelationsHist[0]->GetYaxis()->SetTitle(xAxisTitleMult[1].Data());
   fControlHistogramsList->Add(fMultCorrelationsHist[0]);
  }
  if(!fDontFill[0] && !fDontFill[2])
  {
   fMultCorrelationsHist[1] = new TH2I("Multiplicity (RP vs. REF)","Multiplicity (RP vs. REF)",fnBinsMult[0],fMinMult[0],fMaxMult[0],fnBinsMult[2],fMinMult[2],fMaxMult[2]);
   fMultCorrelationsHist[1]->GetXaxis()->SetTitle(xAxisTitleMult[0].Data());
   fMultCorrelationsHist[1]->GetYaxis()->SetTitle(xAxisTitleMult[2].Data());
   fControlHistogramsList->Add(fMultCorrelationsHist[1]);
  }
  if(!fDontFill[1] && !fDontFill[2])
  {
   fMultCorrelationsHist[2] = new TH2I("Multiplicity (POI vs. REF)","Multiplicity (POI vs. REF)",fnBinsMult[1],fMinMult[1],fMaxMult[1],fnBinsMult[2],fMinMult[2],fMaxMult[2]);
   fMultCorrelationsHist[2]->GetXaxis()->SetTitle(xAxisTitleMult[1].Data());
   fMultCorrelationsHist[2]->GetYaxis()->SetTitle(xAxisTitleMult[2].Data());
   fControlHistogramsList->Add(fMultCorrelationsHist[2]);
  }
 } // if(fFillMultCorrelationsHist){

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForControlHistograms()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::FillControlHistograms(AliFlowEventSimple *anEvent)
{
 // Fill control histograms. 
 // a) Fill TH1D *fKinematicsHist[2][3];
 // b) Fill TH1D *fMultDistributionsHist[3]; 
 // c) Fill TH2D *fMultCorrelationsHist[3].  

 // a) Fill TH1D *fKinematicsHist[2][3]:
 Int_t nCounterRPs = 0;
 if(fFillKinematicsHist)
 {
  Int_t nTracks = anEvent->NumberOfTracks(); // TBI shall I promote this to data member?
  for(Int_t t=0;t<nTracks;t++) // loop over all tracks
  {
   AliFlowTrackSimple *pTrack = NULL;
   if(!fSelectRandomlyRPs) // TBI hw RPs
   {
    pTrack = anEvent->GetTrack(t);
   }
   else
   {
    pTrack = anEvent->GetTrack((Int_t)fRandomIndicesRPs->GetAt(t));
   }
   if(!pTrack){printf("\n Error: pTrack is NULL in MPC::FCH() !!!!");continue;}
   if(pTrack)
   {

    if(!TrackIsInSpecifiedIntervals(pTrack)){continue;} // TBI tmp gym
    
    if(pTrack->InRPSelection())
    {
     nCounterRPs++;
     if(fSelectRandomlyRPs && nCounterRPs == fnSelectedRandomlyRPs){break;} // for(Int_t t=0;t<nTracks;t++) // loop over all tracks
    }
 
    Double_t dPhi = pTrack->Phi(); 
    //if(dPhi > TMath::TwoPi()){dPhi -= TMath::TwoPi();} TBI
    //if(dPhi < 0.){dPhi += TMath::TwoPi();} TBI
    Double_t dPt = pTrack->Pt();
    Double_t dEta = pTrack->Eta();
    Double_t dPhiPtEta[3] = {dPhi,dPt,dEta};

    for(Int_t rp=0;rp<2;rp++) // [RP,POI]
    {
     for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
     {
      if((0==rp && pTrack->InRPSelection()) || (1==rp && pTrack->InPOISelection())) // TBI 
      { 
       if(fKinematicsHist[rp][ppe]){fKinematicsHist[rp][ppe]->Fill(dPhiPtEta[ppe]);}
      }
     } // for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
    } // for(Int_t rp=0;rp<2;rp++) // [RP,POI]
   } // if(pTrack)  
  } // for(Int_t t=0;t<nTracks;t++) // loop over all tracks
 } // if(fFillKinematicsHist)

 // b) Fill TH1D *fMultDistributionsHist[3]:
 Double_t dMultRP = fSelectRandomlyRPs ? nCounterRPs : anEvent->GetNumberOfRPs();
 Double_t dMultPOI = anEvent->GetNumberOfPOIs(); // TBI reimplement when reshuffling is enabled, add support also for the POIs
 Double_t dMultREF = anEvent->GetReferenceMultiplicity();

 if(fSkipSomeIntervals) // TBI tmp gym
 {
  dMultRP = dMultRP - fNumberOfSkippedRPParticles; 
  dMultPOI = -44;
 }

 Double_t dMult[3] = {dMultRP,dMultPOI,dMultREF};
 for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 {
  if(fFillMultDistributionsHist && fMultDistributionsHist[rprm]){fMultDistributionsHist[rprm]->Fill(dMult[rprm]);}      
 } 

 // c) Fill TH2I *fMultCorrelationsHist[3]:  
 if(fFillMultCorrelationsHist)
 {
  if(fMultCorrelationsHist[0]){fMultCorrelationsHist[0]->Fill((Int_t)dMultRP,(Int_t)dMultPOI);} // RP vs. POI
  if(fMultCorrelationsHist[1]){fMultCorrelationsHist[1]->Fill((Int_t)dMultRP,(Int_t)dMultREF);} // RP vs. refMult
  if(fMultCorrelationsHist[2]){fMultCorrelationsHist[2]->Fill((Int_t)dMultPOI,(Int_t)dMultREF);} // POI vs. refMult
 } // if(fFillMultCorrelationsHist)

} // void AliFlowAnalysisWithMultiparticleCorrelations::FillControlHistograms(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForControlHistograms()
{
 // Initialize all arrays for control histograms.

 // a) Initialize TH1D *fKinematicsHist[2][3];
 // b) Initialize TH1D *fMultDistributionsHist[3]; 
 // c) Initialize TH2D *fMultCorrelationsHist[3];  
 // d) Initialize Bool_t fDontFill[3];   
 // e) Initialize default binning values for fKinematicsHist[2][3];
 // f) Initialize default binning values for fMultCorrelationsHist[3];
 // g) Initialize default rp, phi and eta intervals to be skipped.
 
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

 // c) Initialize TH2I *fMultCorrelationsHist[3]: 
 for(Int_t r=0;r<3;r++) // [RP vs. POI, RP vs. refMult, POI vs. refMult]  
 {
  fMultCorrelationsHist[r] = NULL; 
 }

 // d) Initialize Bool_t fDontFill[3]:
 for(Int_t rpr=0;rpr<3;rpr++) // [RP,POI,REF]
 {
  fDontFill[rpr] = kFALSE;
 }

 // e) Initialize default binning values for fKinematicsHist[2][3]:
 // nBins:
 fnBins[0][0] = 360;  // [RP][phi]
 fnBins[0][1] = 1000; // [RP][pt]
 fnBins[0][2] = 1000; // [RP][eta]
 fnBins[1][0] = 360;  // [POI][phi]
 fnBins[1][1] = 1000; // [POI][pt]
 fnBins[1][2] = 1000; // [POI][eta]
 // Min:
 fMin[0][0] = 0.;  // [RP][phi]
 fMin[0][1] = 0.;  // [RP][pt]
 fMin[0][2] = -1.; // [RP][eta]
 fMin[1][0] = 0.;  // [POI][phi]
 fMin[1][1] = 0.;  // [POI][pt]
 fMin[1][2] = -1.; // [POI][eta]
 // Max:
 fMax[0][0] = TMath::TwoPi(); // [RP][phi]
 fMax[0][1] = 10.;            // [RP][pt]
 fMax[0][2] = 1.;             // [RP][eta]
 fMax[1][0] = TMath::TwoPi(); // [POI][phi]
 fMax[1][1] = 10.;            // [POI][pt]
 fMax[1][2] = 1.;             // [POI][eta]

 // f) Initialize default binning values for fMultCorrelationsHist[3]:
 // nBins:
 fnBinsMult[0] = 3000; // [RP]
 fnBinsMult[1] = 3000; // [POI]
 fnBinsMult[2] = 3000; // [REF]
 // Min:
 fMinMult[0] = 0.; // [RP]
 fMinMult[1] = 0.; // [POI]
 fMinMult[2] = 0.; // [REF]
 // Max:
 fMaxMult[0] = 3000.; // [RP]
 fMaxMult[1] = 3000.; // [POI]
 fMaxMult[2] = 3000.; // [REF]

 // g) Initialize default rp, phi and eta intervals to be skipped:
 for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
 {
  for(Int_t i=0;i<10;i++) // interval boundaries, 10 boundaries at max
  { 
   fSkip[ppe][i] = -44.;
  } 
 } // for(Int_t ppe=0;ppe<3;ppe++)

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForControlHistograms()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForQvector()
{
 // Book all the stuff for Q-vector.

 // a) Book the profile holding all the flags for Q-vector;
 // ...

 // a) Book the profile holding all the flags for Q-vector:
 fQvectorFlagsPro = new TProfile("fQvectorFlagsPro","Flags for Q-vectors",2,0,2);
 fQvectorFlagsPro->SetTickLength(-0.01,"Y");
 fQvectorFlagsPro->SetMarkerStyle(25);
 fQvectorFlagsPro->SetLabelSize(0.03);
 fQvectorFlagsPro->SetLabelOffset(0.02,"Y");
 fQvectorFlagsPro->SetStats(kFALSE);
 fQvectorFlagsPro->SetFillColor(kGray);
 fQvectorFlagsPro->SetLineColor(kBlack);
 fQvectorFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateQvector"); fQvectorFlagsPro->Fill(0.5,fCalculateQvector); 
 fQvectorFlagsPro->GetXaxis()->SetBinLabel(2,"fCalculateDiffQvectors"); fQvectorFlagsPro->Fill(1.5,fCalculateDiffQvectors); 
 fQvectorList->Add(fQvectorFlagsPro);

 // ...

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForQvector()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForCorrelations()
{
 // Book all the stuff for correlations.

 // TBI this method can be implemented in a much more civilised way. 

 // a) Book the profile holding all the flags for correlations;
 // b) Book TProfile *fCorrelationsPro[2][8] ([0=cos,1=sin][1p,2p,...,8p]).

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForCorrelations()";

 // a) Book the profile holding all the flags for correlations:
 fCorrelationsFlagsPro = new TProfile("fCorrelationsFlagsPro","Flags for correlations",13,0,13);
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
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(4,"fCalculateIsotropic"); fCorrelationsFlagsPro->Fill(3.5,fCalculateIsotropic); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(5,"fCalculateSame"); fCorrelationsFlagsPro->Fill(4.5,fCalculateSame); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(6,"fSkipZeroHarmonics"); fCorrelationsFlagsPro->Fill(5.5,fSkipZeroHarmonics); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(7,"fCalculateSameIsotropic"); fCorrelationsFlagsPro->Fill(6.5,fCalculateSameIsotropic); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(8,"fCalculateAll"); fCorrelationsFlagsPro->Fill(7.5,fCalculateAll); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(9,"fDontGoBeyond"); fCorrelationsFlagsPro->Fill(8.5,fDontGoBeyond); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(10,"fCalculateOnlyForHarmonicQC"); fCorrelationsFlagsPro->Fill(9.5,fCalculateOnlyForHarmonicQC); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(11,"fCalculateOnlyForSC"); fCorrelationsFlagsPro->Fill(10.5,fCalculateOnlyForSC); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(12,"fCalculateOnlyCos"); fCorrelationsFlagsPro->Fill(11.5,fCalculateOnlyCos); 
 fCorrelationsFlagsPro->GetXaxis()->SetBinLabel(13,"fCalculateOnlySin"); fCorrelationsFlagsPro->Fill(12.5,fCalculateOnlySin);
 fCorrelationsList->Add(fCorrelationsFlagsPro);

 if(!fCalculateCorrelations){return;} // TBI is this safe enough? 

 // b) Book TProfile *fCorrelationsPro[2][8] ([0=cos,1=sin][1p,2p,...,8p]): // TBI hardwired 8, shall be fMaxCorrelator
 cout<<" => Booking TProfile *fCorrelationsPro[2][8]..."<<endl;
 TString sCosSin[2] = {"Cos","Sin"};
 Int_t markerColor[2] = {kBlue,kRed};
 Int_t markerStyle[2] = {kFullSquare,kFullSquare};
 Int_t nBins[8] = {1,1,1,1,1,1,1,1}; // TBI hardwired 8, shall be fMaxCorrelator
 Int_t nBinsTitle[8] = {1,1,1,1,1,1,1,1}; // TBI hardwired 8, shall be fMaxCorrelator
 Int_t nToBeFilled[8] = {0,0,0,0,0,0,0,0}; // TBI hardwired 8, shall be fMaxCorrelator
 for(Int_t c=0;c<fMaxCorrelator;c++) // [1p,2p,...,8p]
 {
  // Implementing \binom{n+k-1}{k}, which is the resulting number of sets obtained
  // after sampling n starting elements into k subsets, repetitions allowed.
  // In my case, n=2*fMaxHarmonic+1, k=c+1, hence:
  nBins[c] = (Int_t)(TMath::Factorial(2*fMaxHarmonic+1+c+1-1)
           / (TMath::Factorial(2*fMaxHarmonic+1-1)*TMath::Factorial(c+1)));
  nBinsTitle[c] = nBins[c];
  if(c>=fDontGoBeyond){nBins[c]=1;} // TBI is this really safe? 
 } // for(Int_t c=0;c<8;c++) // [1p,2p,...,8p]
 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  if(fCalculateOnlyCos && 1==cs){continue;}
  else if(fCalculateOnlySin && 0==cs){continue;}
  for(Int_t c=0;c<fMaxCorrelator;c++) // [1p,2p,...,8p]
  {
   if(c==fDontGoBeyond){continue;}
   if(fCalculateOnlyForHarmonicQC && c%2==0){continue;}
   fCorrelationsPro[cs][c] = new TProfile(Form("%dpCorrelations%s",c+1,sCosSin[cs].Data()),"",nBins[c],0.,1.*nBins[c]);
   fCorrelationsPro[cs][c]->Sumw2();
   fCorrelationsPro[cs][c]->SetStats(kFALSE);
   fCorrelationsPro[cs][c]->SetMarkerColor(markerColor[cs]);
   fCorrelationsPro[cs][c]->SetMarkerStyle(markerStyle[cs]);
   fCorrelationsList->Add(fCorrelationsPro[cs][c]);
  } // for(Int_t c=0;c<8;c++) // [1p,2p,...,8p]
 } // for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 // Set all bin labels: TBI this can be implemented better, most likely...
 Int_t binNo[2][8]; 
 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  if(fCalculateOnlyCos && 1==cs){continue;}
  else if(fCalculateOnlySin && 0==cs){continue;}
  for(Int_t c=0;c<fMaxCorrelator;c++)
  {
   binNo[cs][c] = 1;
  } 
 } // for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 
 for(Int_t n1=-fMaxHarmonic;n1<=fMaxHarmonic;n1++) 
 {
  cout<< Form("    Patience, this takes some time... n1 = %d/%d\r",n1+fMaxHarmonic,2*fMaxHarmonic)<<flush; // TBI
  if(fSkipZeroHarmonics && 0==n1){continue;}
  if(fCalculateOnlyForHarmonicQC && TMath::Abs(n1) != fHarmonicQC){continue;}
  if(fCalculateAll)
  {
   for(Int_t cs=0;cs<2;cs++) 
   {
    if(fCalculateOnlyCos && 1==cs){continue;}
    else if(fCalculateOnlySin && 0==cs){continue;}
    if(fCorrelationsPro[cs][0]){fCorrelationsPro[cs][0]->GetXaxis()->SetBinLabel(binNo[cs][0]++,Form("%s(%d)",sCosSin[cs].Data(),n1));}
   } // for(Int_t cs=0;cs<2;cs++) 
   nToBeFilled[0]++; 
  }
  if(1==fDontGoBeyond){continue;}
  for(Int_t n2=n1;n2<=fMaxHarmonic;n2++) 
  {
   if(fSkipZeroHarmonics && 0==n2){continue;}
   if(fCalculateOnlyForHarmonicQC && TMath::Abs(n2) != fHarmonicQC){continue;}
   if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2) || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2)) 
      || (fCalculateSameIsotropic && 0==n1+n2 && TMath::Abs(n1)==TMath::Abs(n2)) 
      || (fCalculateOnlyForHarmonicQC && 0==n1+n2)
      || (fCalculateOnlyForSC && 0==n1+n2))
   {  
    for(Int_t cs=0;cs<2;cs++) 
    {
     if(fCalculateOnlyCos && 1==cs){continue;}
     else if(fCalculateOnlySin && 0==cs){continue;}
     if(fCorrelationsPro[cs][1]){fCorrelationsPro[cs][1]->GetXaxis()->SetBinLabel(binNo[cs][1]++,Form("%s(%d,%d)",sCosSin[cs].Data(),n1,n2));}
    } // for(Int_t cs=0;cs<2;cs++) 
    nToBeFilled[1]++; 
   }
   if(2==fDontGoBeyond){continue;}
   for(Int_t n3=n2;n3<=fMaxHarmonic;n3++) 
   {
    if(fSkipZeroHarmonics && 0==n3){continue;}
    if(fCalculateOnlyForHarmonicQC && TMath::Abs(n3) != fHarmonicQC){continue;}
    if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3) || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3))
       || (fCalculateSameIsotropic && 0==n1+n2+n3 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3)) 
       || (fCalculateOnlyForHarmonicQC && 0==n1+n2+n3))
    {  
     for(Int_t cs=0;cs<2;cs++) 
     {
      if(fCalculateOnlyCos && 1==cs){continue;}
      else if(fCalculateOnlySin && 0==cs){continue;}
      if(fCorrelationsPro[cs][2]){fCorrelationsPro[cs][2]->GetXaxis()->SetBinLabel(binNo[cs][2]++,Form("%s(%d,%d,%d)",sCosSin[cs].Data(),n1,n2,n3));}
     } // for(Int_t cs=0;cs<2;cs++) 
     nToBeFilled[2]++; 
    }
    if(3==fDontGoBeyond){continue;}
    for(Int_t n4=n3;n4<=fMaxHarmonic;n4++) 
    {
     if(fSkipZeroHarmonics && 0==n4){continue;}
     if(fCalculateOnlyForHarmonicQC && TMath::Abs(n4) != fHarmonicQC){continue;}
     if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4) || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4))
       || (fCalculateSameIsotropic && 0==n1+n2+n3+n4 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4)) 
       || (fCalculateOnlyForHarmonicQC && 0==n1+n2+n3+n4)
       || (fCalculateOnlyForSC && (0==n1+n4 && 0==n2+n3 && n1 != n2 && n3 != n4)))
     {   
      for(Int_t cs=0;cs<2;cs++) 
      {
       if(fCalculateOnlyCos && 1==cs){continue;}
       else if(fCalculateOnlySin && 0==cs){continue;}
       if(fCorrelationsPro[cs][3]){fCorrelationsPro[cs][3]->GetXaxis()->SetBinLabel(binNo[cs][3]++,Form("%s(%d,%d,%d,%d)",sCosSin[cs].Data(),n1,n2,n3,n4));}
      } // for(Int_t cs=0;cs<2;cs++) 
      nToBeFilled[3]++; 
     } 
     if(4==fDontGoBeyond){continue;}
     for(Int_t n5=n4;n5<=fMaxHarmonic;n5++) 
     {
      if(fSkipZeroHarmonics && 0==n5){continue;}
      if(fCalculateOnlyForHarmonicQC && TMath::Abs(n5) != fHarmonicQC){continue;}
      if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5) 
         || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5))
         || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5)) 
         || (fCalculateOnlyForHarmonicQC && 0==n1+n2+n3+n4+n5))
      {   
       for(Int_t cs=0;cs<2;cs++) 
       {
        if(fCalculateOnlyCos && 1==cs){continue;}
        else if(fCalculateOnlySin && 0==cs){continue;}
        if(fCorrelationsPro[cs][4]){fCorrelationsPro[cs][4]->GetXaxis()->SetBinLabel(binNo[cs][4]++,Form("%s(%d,%d,%d,%d,%d)",sCosSin[cs].Data(),n1,n2,n3,n4,n5));}
       } // for(Int_t cs=0;cs<2;cs++) 
       nToBeFilled[4]++; 
      }
      if(5==fDontGoBeyond){continue;}
      for(Int_t n6=n5;n6<=fMaxHarmonic;n6++) 
      {
       if(fSkipZeroHarmonics && 0==n6){continue;}
       if(fCalculateOnlyForHarmonicQC && TMath::Abs(n6) != fHarmonicQC){continue;}
       if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5+n6)  
          || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) 
              && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6))
          || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5+n6 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) 
              && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6)) 
          || (fCalculateOnlyForHarmonicQC && 0==n1+n2+n3+n4+n5+n6))
       {   
        for(Int_t cs=0;cs<2;cs++) 
        {
         if(fCalculateOnlyCos && 1==cs){continue;}
         else if(fCalculateOnlySin && 0==cs){continue;}
         if(fCorrelationsPro[cs][5]){fCorrelationsPro[cs][5]->GetXaxis()->SetBinLabel(binNo[cs][5]++,Form("%s(%d,%d,%d,%d,%d,%d)",sCosSin[cs].Data(),n1,n2,n3,n4,n5,n6));}         
        } // for(Int_t cs=0;cs<2;cs++) 
        nToBeFilled[5]++; 
       }
       if(6==fDontGoBeyond){continue;}
       for(Int_t n7=n6;n7<=fMaxHarmonic;n7++) 
       {
        if(fSkipZeroHarmonics && 0==n7){continue;}
        if(fCalculateOnlyForHarmonicQC && TMath::Abs(n7) != fHarmonicQC){continue;}
        if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5+n6+n7) 
           || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) 
               && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) && TMath::Abs(n1)==TMath::Abs(n7))
           || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5+n6+n7 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) 
               && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) && TMath::Abs(n1)==TMath::Abs(n7)) 
           || (fCalculateOnlyForHarmonicQC && 0==n1+n2+n3+n4+n5+n6+n7))
        {   
         for(Int_t cs=0;cs<2;cs++) 
         {
          if(fCalculateOnlyCos && 1==cs){continue;}
          else if(fCalculateOnlySin && 0==cs){continue;}
          if(fCorrelationsPro[cs][6]){fCorrelationsPro[cs][6]->GetXaxis()->SetBinLabel(binNo[cs][6]++,Form("%s(%d,%d,%d,%d,%d,%d,%d)",sCosSin[cs].Data(),n1,n2,n3,n4,n5,n6,n7));}
         } // for(Int_t cs=0;cs<2;cs++) 
         nToBeFilled[6]++; 
        }
        if(7==fDontGoBeyond){continue;}
        for(Int_t n8=n7;n8<=fMaxHarmonic;n8++) 
        {
         if(fSkipZeroHarmonics && 0==n8){continue;}
         if(fCalculateOnlyForHarmonicQC && TMath::Abs(n8) != fHarmonicQC){continue;}
         if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5+n6+n7+n8) 
            || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) 
                && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) && TMath::Abs(n1)==TMath::Abs(n7) && TMath::Abs(n1)==TMath::Abs(n8))
            || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5+n6+n7+n8 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) 
                && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) && TMath::Abs(n1)==TMath::Abs(n7) 
                && TMath::Abs(n1)==TMath::Abs(n8))
            || (fCalculateOnlyForHarmonicQC && 0==n1+n2+n3+n4+n5+n6+n7+n8))
         {    
          for(Int_t cs=0;cs<2;cs++) 
          {
           if(fCalculateOnlyCos && 1==cs){continue;}
           else if(fCalculateOnlySin && 0==cs){continue;}
           if(fCorrelationsPro[cs][7]){fCorrelationsPro[cs][7]->GetXaxis()->SetBinLabel(binNo[cs][7]++,Form("%s(%d,%d,%d,%d,%d,%d,%d,%d)",sCosSin[cs].Data(),n1,n2,n3,n4,n5,n6,n7,n8));}
          } // for(Int_t cs=0;cs<2;cs++) 
          nToBeFilled[7]++; 
         }
        } // for(Int_t n8=n7;n8<=fMaxHarmonic;n8++)
       } // for(Int_t n7=n6;n7<=fMaxHarmonic;n7++) 
      } // for(Int_t n6=n5;n6<=fMaxHarmonic;n6++) 
     } // for(Int_t n5=n4;n5<=fMaxHarmonic;n5++) 
    } // for(Int_t n4=n3;n4<=fMaxHarmonic;n4++)   
   } // for(Int_t n3=n2;n3<=fMaxHarmonic;n3++) 
  } // for(Int_t n2=n1;n2<=fMaxHarmonic;n2++)
 } // for(Int_t n1=-fMaxHarmonic;n1<=fMaxHarmonic;n1++) 

 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  if(fCalculateOnlyCos && 1==cs){continue;}
  else if(fCalculateOnlySin && 0==cs){continue;}
  for(Int_t c=0;c<fMaxCorrelator;c++) // [1p,2p,...,8p]
  {
   if(!fCorrelationsPro[cs][c]){continue;}
   fCorrelationsPro[cs][c]->SetTitle(Form("%d-p correlations, %s terms, %d/%d in total",c+1,sCosSin[cs].Data(),nToBeFilled[c],nBinsTitle[c]));
   fCorrelationsPro[cs][c]->GetXaxis()->SetRangeUser(0.,fCorrelationsPro[cs][c]->GetBinLowEdge(nToBeFilled[c]+1));
  }
 } 
 cout<<"    Booked.                                           "<<endl; // TBI 

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForDiffCorrelations()
{
 // Book all the stuff for differential correlations.

 // a) Book the profile holding all the flags for differential correlations;
 // b) Book TProfile *fDiffCorrelationsPro[2][4] ([0=cos,1=sin][1p,2p,3p,4p]).

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForDiffCorrelations()";

 // a) Book the profile holding all the flags for differential correlations:
 fDiffCorrelationsFlagsPro = new TProfile("fDiffCorrelationsFlagsPro","Flags for differential correlations",5,0,5);
 fDiffCorrelationsFlagsPro->SetTickLength(-0.01,"Y");
 fDiffCorrelationsFlagsPro->SetMarkerStyle(25);
 fDiffCorrelationsFlagsPro->SetLabelSize(0.03);
 fDiffCorrelationsFlagsPro->SetLabelOffset(0.02,"Y");
 fDiffCorrelationsFlagsPro->SetStats(kFALSE);
 fDiffCorrelationsFlagsPro->SetFillColor(kGray);
 fDiffCorrelationsFlagsPro->SetLineColor(kBlack);
 fDiffCorrelationsFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateDiffCorrelations"); fDiffCorrelationsFlagsPro->Fill(0.5,fCalculateDiffCorrelations); 
 fDiffCorrelationsFlagsPro->GetXaxis()->SetBinLabel(2,"fCalculateDiffCos"); fDiffCorrelationsFlagsPro->Fill(1.5,fCalculateDiffCos); 
 fDiffCorrelationsFlagsPro->GetXaxis()->SetBinLabel(3,"fCalculateDiffSin"); fDiffCorrelationsFlagsPro->Fill(2.5,fCalculateDiffSin); 
 fDiffCorrelationsFlagsPro->GetXaxis()->SetBinLabel(4,"fCalculateDiffCorrelationsVsPt"); fDiffCorrelationsFlagsPro->Fill(3.5,fCalculateDiffCorrelationsVsPt); 
 fDiffCorrelationsFlagsPro->GetXaxis()->SetBinLabel(5,"fUseDefaultBinning"); fDiffCorrelationsFlagsPro->Fill(4.5,fUseDefaultBinning); 
 fDiffCorrelationsList->Add(fDiffCorrelationsFlagsPro);

 if(!fCalculateDiffCorrelations){return;}  

 // b) Book TProfile *fDiffCorrelationsPro[2][4] ([0=cos,1=sin][1p,2p,3p,4p]):
 Bool_t fDiffStore[2][4] = {{0,1,1,1},{0,0,0,0}}; // store or not TBI promote to data member, and implement setter perhaps  
 Int_t markerColor[2] = {kRed,kGreen};
 Int_t markerStyle[2] = {kFullSquare,kOpenSquare};
 TString sCosSin[2] = {"Cos","Sin"};
 TString sLabel[4] = {Form("%d",fDiffHarmonics[0][0]),
                      Form("%d,%d",fDiffHarmonics[1][0],fDiffHarmonics[1][1]),
                      Form("%d,%d,%d",fDiffHarmonics[2][0],fDiffHarmonics[2][1],fDiffHarmonics[2][2]),
                      Form("%d,%d,%d,%d",fDiffHarmonics[3][0],fDiffHarmonics[3][1],fDiffHarmonics[3][2],fDiffHarmonics[3][3])};

 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  if(!fCalculateDiffCos && 0==cs){continue;}
  if(!fCalculateDiffSin && 1==cs){continue;}

  for(Int_t c=0;c<4;c++) // [1p,2p,3p,4p]
  {
   if(fCalculateDiffCorrelationsVsPt)
   {
    if(fUseDefaultBinning)
    {
     // vs pt, default binning:  
     fDiffCorrelationsPro[cs][c] = new TProfile(Form("%s, %dp, %s",sCosSin[cs].Data(),c+1,"pt"),
                                                Form("%s(%s)",sCosSin[cs].Data(),sLabel[c].Data()),
                                                100,0.,10.);
    } else // if(fUseDefaultBinning)
      {
       // vs pt, non-default binning:
       fDiffCorrelationsPro[cs][c] = new TProfile(Form("%s, %dp, %s",sCosSin[cs].Data(),c+1,"pt"),
                                                  Form("%s(%s)",sCosSin[cs].Data(),sLabel[c].Data()),
                                                  fnDiffBins,fRangesDiffBins);
      }// else // if(fUseDefaultBinning) 
      fDiffCorrelationsPro[cs][c]->Sumw2();
      fDiffCorrelationsPro[cs][c]->SetStats(kFALSE);
      fDiffCorrelationsPro[cs][c]->SetMarkerColor(markerColor[cs]);
      fDiffCorrelationsPro[cs][c]->SetMarkerStyle(markerStyle[cs]);
      fDiffCorrelationsPro[cs][c]->GetXaxis()->SetTitle("p_{T}");
      if(fDiffStore[cs][c]){fDiffCorrelationsList->Add(fDiffCorrelationsPro[cs][c]);}
   } else // if(fCalculateDiffCorrelationsVsPt)
     {
      if(fUseDefaultBinning)
      {
       // vs eta, default binning:
       fDiffCorrelationsPro[cs][c] = new TProfile(Form("%s, %dp, %s",sCosSin[cs].Data(),c+1,"eta"),
                                                  Form("%s(%s)",sCosSin[cs].Data(),sLabel[c].Data()),
                                                  100,-1.,1.);
      } else // if(fUseDefaultBinning)
        {
         // vs eta, non-default binning:
         fDiffCorrelationsPro[cs][c] = new TProfile(Form("%s, %dp, %s",sCosSin[cs].Data(),c+1,"eta"),
                                                    Form("%s(%s)",sCosSin[cs].Data(),sLabel[c].Data()),
                                                    fnDiffBins,fRangesDiffBins);
        } // else // if(fUseDefaultBinning)
        fDiffCorrelationsPro[cs][c]->Sumw2();
        fDiffCorrelationsPro[cs][c]->SetStats(kFALSE);
        fDiffCorrelationsPro[cs][c]->SetMarkerColor(markerColor[cs]);
        fDiffCorrelationsPro[cs][c]->SetMarkerStyle(markerStyle[cs]);
        fDiffCorrelationsPro[cs][c]->GetXaxis()->SetTitle("#eta");
        if(fDiffStore[cs][c]){fDiffCorrelationsList->Add(fDiffCorrelationsPro[cs][c]);}
     } // else // if(fCalculateDiffCorrelationsVsPt)
  } // for(Int_t c=0;c<4;c++) // [1p,2p,3p,4p]
 } // for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForDiffCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForSymmetryPlanes()
{
 // Book all the stuff for symmetry plane correlations.

 // a) Book the profile holding all the flags for symmetry plane correlations;
 // b) Book TProfile *fSymmetryPlanesPro[?][?]. TBI check the exact dimensions in the header file.

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForSymmetryPlanes()";

 // a) Book the profile holding all the flags for symmetry plane correlations:
 fSymmetryPlanesFlagsPro = new TProfile("fSymmetryPlanesFlagsPro","Flags for Symmetry Plane Correlations (SPC)",1,0,1);
 fSymmetryPlanesFlagsPro->SetTickLength(-0.01,"Y");
 fSymmetryPlanesFlagsPro->SetMarkerStyle(25);
 fSymmetryPlanesFlagsPro->SetLabelSize(0.03);
 fSymmetryPlanesFlagsPro->SetLabelOffset(0.02,"Y");
 fSymmetryPlanesFlagsPro->SetStats(kFALSE);
 fSymmetryPlanesFlagsPro->SetFillColor(kGray);
 fSymmetryPlanesFlagsPro->SetLineColor(kBlack);
 fSymmetryPlanesFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateSymmetryPlanes"); fSymmetryPlanesFlagsPro->Fill(0.5,fCalculateSymmetryPlanes);
 //fSymmetryPlanesFlagsPro->GetXaxis()->SetBinLabel(2,"TBI"); fSymmetryPlanesFlagsPro->Fill(1.5,TBI);
 fSymmetryPlanesList->Add(fSymmetryPlanesFlagsPro);

 if(!fCalculateSymmetryPlanes){return;}

 // b) Book TProfile *fSymmetryPlanesPro[?][?]: TBI check the exact dimensions in the header file
 TString sTitle[1][2] = {"#LT#LTcos[4(#psi_{2}-#psi_{1})]#GT#GT","#LT#LTcos[4(#psi_{4}-#psi_{2})]#GT#GT"}; // TBI check h.w. harmonic 4, as well as h.w. indices on symmetry planes
 TString sLabels[4] = {"k = 0","k = 1","k = 2","k = 3"};
 for(Int_t gc=0;gc<1;gc++) // 'generic correlator': [[0]:(Psi2n,Psi1n),[1]:...] TBI upper boundary will change
 {
  for(Int_t n=0;n<2;n++) // 'harmonic n for generic correlator': [[0]:n=1,[1]:n=2,...] TBI upper boundary will change
  {
   fSymmetryPlanesPro[gc][n] = new TProfile(Form("%d,%d",gc,n),sTitle[gc][n].Data(),4,0.,4.); // TBI this is a landmine, introduce fnHighestOptimizerSPC eventually instead of '4'
   fSymmetryPlanesPro[gc][n]->Sumw2();
   fSymmetryPlanesPro[gc][n]->SetStats(kFALSE);
   fSymmetryPlanesPro[gc][n]->SetMarkerColor(kRed);
   fSymmetryPlanesPro[gc][n]->SetMarkerStyle(kFullSquare);
   fSymmetryPlanesPro[gc][n]->SetLineColor(kBlack);
   fSymmetryPlanesPro[gc][n]->SetLabelSize(0.0544);
   for(Int_t o=0;o<4;o++) // TBI h.w. '4'
   {
    fSymmetryPlanesPro[gc][n]->GetXaxis()->SetBinLabel(o+1,sLabels[o].Data());
   } // for(Int_t o=0;o<4;o++)
   fSymmetryPlanesList->Add(fSymmetryPlanesPro[gc][n]);
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForSymmetryPlanes()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForEtaGaps()
{
 // Book all the stuff for correlations with eta gaps.

 // a) Book the profile holding all the flags for correlations with eta gaps;  
 // b) Book TProfile *fEtaGapsPro[6][10]

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForEtaGaps()";

 // a) Book the profile holding all the flags for correlations with eta gaps:
 fEtaGapsFlagsPro = new TProfile("fEtaGapsFlagsPro","Flags for correlations with eta gaps",1,0,1);
 fEtaGapsFlagsPro->SetTickLength(-0.01,"Y");
 fEtaGapsFlagsPro->SetMarkerStyle(25);
 fEtaGapsFlagsPro->SetLabelSize(0.03);
 fEtaGapsFlagsPro->SetLabelOffset(0.02,"Y");
 fEtaGapsFlagsPro->SetStats(kFALSE);
 fEtaGapsFlagsPro->SetFillColor(kGray);
 fEtaGapsFlagsPro->SetLineColor(kBlack);
 fEtaGapsFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateEtaGaps"); fEtaGapsFlagsPro->Fill(0.5,fCalculateEtaGaps);
 fEtaGapsFlagsPro->GetXaxis()->SetBinLabel(2,"fLowestHarmonicEtaGaps"); fEtaGapsFlagsPro->Fill(1.5,fLowestHarmonicEtaGaps);
 fEtaGapsFlagsPro->GetXaxis()->SetBinLabel(3,"fHighestHarmonicEtaGaps"); fEtaGapsFlagsPro->Fill(2.5,fHighestHarmonicEtaGaps);
 fEtaGapsList->Add(fEtaGapsFlagsPro);

 if(!fCalculateEtaGaps){return;}

 // b) Book TProfile *fEtaGapsPro[6][10]:
 Int_t markerColor[6] = {kBlue,kRed,kGreen+2,kBlack,kBlack,kBlack};
 TString sEtaGaps[11] = {"|#Delta#eta| > 1.0","|#Delta#eta| > 0.9","|#Delta#eta| > 0.8","|#Delta#eta| > 0.7","|#Delta#eta| > 0.6","|#Delta#eta| > 0.5","|#Delta#eta| > 0.4","|#Delta#eta| > 0.3","|#Delta#eta| > 0.2","|#Delta#eta| > 0.1","|#Delta#eta| > 0.0"};
 for(Int_t h=fLowestHarmonicEtaGaps-1;h<=fHighestHarmonicEtaGaps-1;h++) // harmonics
 {
  fEtaGapsPro[h] = new TProfile(Form("EtaGaps_v%d",h+1),Form("2-p correlation for harmonic v_{%d} with eta gaps",h+1),11,0.,11.); 
  fEtaGapsPro[h]->Sumw2();
  fEtaGapsPro[h]->SetStats(kFALSE);
  fEtaGapsPro[h]->SetMarkerColor(markerColor[h]);
  fEtaGapsPro[h]->SetMarkerStyle(kFullSquare);
  fEtaGapsPro[h]->SetLineColor(markerColor[h]);
  fEtaGapsPro[h]->SetLabelSize(0.0544);
  fEtaGapsPro[h]->SetMinimum(0.0);
  for(Int_t eg=0;eg<11;eg++)
  {
   fEtaGapsPro[h]->GetXaxis()->SetBinLabel(eg+1,sEtaGaps[eg].Data());
  }
  fEtaGapsList->Add(fEtaGapsPro[h]);
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForEtaGaps()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForEbECumulants()
{
 // Book all the stuff for event-by-event cumulants.

 // a) Book the profile holding all the flags for e-b-e cumulants;
 // b) Book TProfile *fEbECumulantsPro[2][8] ([0=cos,1=sin][1p,2p,...,8p]).

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForEbECumulants()";

 // a) Book the profile holding all the flags for e-b-e cumulants:
 fEbECumulantsFlagsPro = new TProfile("fEbECumulantsFlagsPro","Flags for e-b-e cumulants",1,0,1);
 fEbECumulantsFlagsPro->SetTickLength(-0.01,"Y");
 fEbECumulantsFlagsPro->SetMarkerStyle(25);
 fEbECumulantsFlagsPro->SetLabelSize(0.03);
 fEbECumulantsFlagsPro->SetLabelOffset(0.02,"Y");
 fEbECumulantsFlagsPro->SetStats(kFALSE);
 fEbECumulantsFlagsPro->SetFillColor(kGray);
 fEbECumulantsFlagsPro->SetLineColor(kBlack);
 fEbECumulantsFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateEbECumulants"); fEbECumulantsFlagsPro->Fill(0.5,fCalculateEbECumulants); 
 fEbECumulantsList->Add(fEbECumulantsFlagsPro);

 if(!fCalculateEbECumulants){return;} // TBI is this safe enough? 

 // b) Book TProfile *fEbECumulantsPro[2][8] ([0=cos,1=sin][1p,2p,...,8p]):
 TString sCosSin[2] = {"Cos","Sin"};
 Int_t markerColor[2] = {kBlue,kRed};
 Int_t markerStyle[2] = {kFullSquare,kFullSquare};
 Int_t nBins[8] = {1,1,1,1,1,1,1,1}; // TBI hardwired 8, shall be fMaxCorrelator
 Int_t nBinsTitle[8] = {1,1,1,1,1,1,1,1}; // TBI hardwired 8, shall be fMaxCorrelator
 Int_t nToBeFilled[8] = {0,0,0,0,0,0,0,0}; // TBI hardwired 8, shall be fMaxCorrelator
 for(Int_t c=0;c<fMaxCorrelator;c++) // [1p,2p,...,8p]
 {
  // Implementing \binom{n+k-1}{k}, which is the resulting number of sets obtained
  // after sampling n starting elements into k subsets, repetitions allowed.
  // In my case, n=2*fMaxHarmonic+1, k=c+1, hence:
  nBins[c] = (Int_t)(TMath::Factorial(2*fMaxHarmonic+1+c+1-1)
           / (TMath::Factorial(2*fMaxHarmonic+1-1)*TMath::Factorial(c+1)));
  nBinsTitle[c] = nBins[c];
  if(c>=fDontGoBeyond){nBins[c]=1;} // TBI a bit of spaghetti here...
 } // for(Int_t c=0;c<8;c++) // [1p,2p,...,8p]
 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  for(Int_t c=0;c<fMaxCorrelator;c++) // [1p,2p,...,8p]
  {
   fEbECumulantsPro[cs][c] = new TProfile(Form("%dpEbECumulants%s",c+1,sCosSin[cs].Data()),"",nBins[c],0.,1.*nBins[c]);
   fEbECumulantsPro[cs][c]->Sumw2();
   fEbECumulantsPro[cs][c]->SetStats(kFALSE);
   fEbECumulantsPro[cs][c]->SetMarkerColor(markerColor[cs]);
   fEbECumulantsPro[cs][c]->SetMarkerStyle(markerStyle[cs]);
   fEbECumulantsList->Add(fEbECumulantsPro[cs][c]);
  } // for(Int_t c=0;c<8;c++) // [1p,2p,...,8p]
 } // for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 // Set all bin labels: TBI this can be implemented better, most likely...
 Int_t binNo[8]; for(Int_t c=0;c<fMaxCorrelator;c++){binNo[c]=1;} // TBI hardwired 8, shall be fMaxCorrelator
 for(Int_t n1=-fMaxHarmonic;n1<=fMaxHarmonic;n1++) 
 {
  if(fSkipZeroHarmonics && 0==n1){continue;}
  if(fCalculateAll)
  {
   fEbECumulantsPro[0][0]->GetXaxis()->SetBinLabel(binNo[0],Form("Cos(%d)",n1));
   fEbECumulantsPro[1][0]->GetXaxis()->SetBinLabel(binNo[0]++,Form("Sin(%d)",n1));
   nToBeFilled[0]++; 
  }
  if(1==fDontGoBeyond){continue;}
  for(Int_t n2=n1;n2<=fMaxHarmonic;n2++) 
  {
   if(fSkipZeroHarmonics && 0==n2){continue;}
   if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2) || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2)) 
      || (fCalculateSameIsotropic && 0==n1+n2 &&  TMath::Abs(n1)==TMath::Abs(n2)))
   {  
    fEbECumulantsPro[0][1]->GetXaxis()->SetBinLabel(binNo[1],Form("Cos(%d,%d)",n1,n2));
    fEbECumulantsPro[1][1]->GetXaxis()->SetBinLabel(binNo[1]++,Form("Sin(%d,%d)",n1,n2));
    nToBeFilled[1]++; 
   }
   if(2==fDontGoBeyond){continue;}
   for(Int_t n3=n2;n3<=fMaxHarmonic;n3++) 
   {
    if(fSkipZeroHarmonics && 0==n3){continue;}
    if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3) || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3))
       || (fCalculateSameIsotropic && 0==n1+n2+n3 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3)))
    {  
     fEbECumulantsPro[0][2]->GetXaxis()->SetBinLabel(binNo[2],Form("Cos(%d,%d,%d)",n1,n2,n3));
     fEbECumulantsPro[1][2]->GetXaxis()->SetBinLabel(binNo[2]++,Form("Sin(%d,%d,%d)",n1,n2,n3));
     nToBeFilled[2]++; 
    }
    if(3==fDontGoBeyond){continue;}
    for(Int_t n4=n3;n4<=fMaxHarmonic;n4++) 
    {
     if(fSkipZeroHarmonics && 0==n4){continue;}
     if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4) || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4))
       || (fCalculateSameIsotropic && 0==n1+n2+n3+n4 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4)))
     {   
      fEbECumulantsPro[0][3]->GetXaxis()->SetBinLabel(binNo[3],Form("Cos(%d,%d,%d,%d)",n1,n2,n3,n4));
      fEbECumulantsPro[1][3]->GetXaxis()->SetBinLabel(binNo[3]++,Form("Sin(%d,%d,%d,%d)",n1,n2,n3,n4)); 
      nToBeFilled[3]++; 
     } 
     if(4==fDontGoBeyond){continue;}
     for(Int_t n5=n4;n5<=fMaxHarmonic;n5++) 
     {
      if(fSkipZeroHarmonics && 0==n5){continue;}
      if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5) 
         || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5))
         || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5)))
      {   
       fEbECumulantsPro[0][4]->GetXaxis()->SetBinLabel(binNo[4],Form("Cos(%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5));
       fEbECumulantsPro[1][4]->GetXaxis()->SetBinLabel(binNo[4]++,Form("Sin(%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5));
       nToBeFilled[4]++; 
      }
      if(5==fDontGoBeyond){continue;}
      for(Int_t n6=n5;n6<=fMaxHarmonic;n6++) 
      {
       if(fSkipZeroHarmonics && 0==n6){continue;}
       if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5+n6)  
          || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) 
              && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6))
          || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5+n6 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) 
              && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6)))
       {   
        fEbECumulantsPro[0][5]->GetXaxis()->SetBinLabel(binNo[5],Form("Cos(%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6));
        fEbECumulantsPro[1][5]->GetXaxis()->SetBinLabel(binNo[5]++,Form("Sin(%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6));
        nToBeFilled[5]++; 
       }
       if(6==fDontGoBeyond){continue;}
       for(Int_t n7=n6;n7<=fMaxHarmonic;n7++) 
       {
        if(fSkipZeroHarmonics && 0==n7){continue;}
        if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5+n6+n7) 
           || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) 
               && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) && TMath::Abs(n1)==TMath::Abs(n7))
           || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5+n6+n7 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) 
               && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) && TMath::Abs(n1)==TMath::Abs(n7)))
        {   
         fEbECumulantsPro[0][6]->GetXaxis()->SetBinLabel(binNo[6],Form("Cos(%d,%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6,n7));
         fEbECumulantsPro[1][6]->GetXaxis()->SetBinLabel(binNo[6]++,Form("Sin(%d,%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6,n7));
         nToBeFilled[6]++; 
        }
        if(7==fDontGoBeyond){continue;}
        for(Int_t n8=n7;n8<=fMaxHarmonic;n8++) 
        {
         if(fSkipZeroHarmonics && 0==n8){continue;}
         if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4+n5+n6+n7+n8) 
            || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4) 
                && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) && TMath::Abs(n1)==TMath::Abs(n7) && TMath::Abs(n1)==TMath::Abs(n8))
            || (fCalculateSameIsotropic && 0==n1+n2+n3+n4+n5+n6+n7+n8 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) 
                && TMath::Abs(n1)==TMath::Abs(n4) && TMath::Abs(n1)==TMath::Abs(n5) && TMath::Abs(n1)==TMath::Abs(n6) && TMath::Abs(n1)==TMath::Abs(n7) 
                && TMath::Abs(n1)==TMath::Abs(n8)))
         {    
          fEbECumulantsPro[0][7]->GetXaxis()->SetBinLabel(binNo[7],Form("Cos(%d,%d,%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6,n7,n8));
          fEbECumulantsPro[1][7]->GetXaxis()->SetBinLabel(binNo[7]++,Form("Sin(%d,%d,%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6,n7,n8));
          nToBeFilled[7]++; 
         }
        } // for(Int_t n8=n7;n8<=fMaxHarmonic;n8++)
       } // for(Int_t n7=n6;n7<=fMaxHarmonic;n7++) 
      } // for(Int_t n6=n5;n6<=fMaxHarmonic;n6++) 
     } // for(Int_t n5=n4;n5<=fMaxHarmonic;n5++) 
    } // for(Int_t n4=n3;n4<=fMaxHarmonic;n4++)   
   } // for(Int_t n3=n2;n3<=fMaxHarmonic;n3++) 
  } // for(Int_t n2=n1;n2<=fMaxHarmonic;n2++)
 } // for(Int_t n1=-fMaxHarmonic;n1<=fMaxHarmonic;n1++) 

 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  for(Int_t c=0;c<fMaxCorrelator;c++) // [1p,2p,...,8p]
  {
   fEbECumulantsPro[cs][c]->SetTitle(Form("%d-p e-b-e cumulants, %s terms, %d/%d in total",c+1,sCosSin[cs].Data(),nToBeFilled[c],nBinsTitle[c]));
   fEbECumulantsPro[cs][c]->GetXaxis()->SetRangeUser(0.,fEbECumulantsPro[cs][c]->GetBinLowEdge(nToBeFilled[c]+1));
  }
 } 
 
} // end of void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForEbECumulants()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForNestedLoops()
{
 // Book all the stuff for nested loops.

 // TBI this method is just ugly, who implemented it like this... 

 // a) Set default harmonic values; 
 // b) Book the profile holding all the flags for nested loops;
 // c) Book the profile holding all results for nested loops (cosine);
 // d) Book the profile holding all results for nested loops (sinus);
 // e) Book the profile holding all results for differential nested loops.

 // a) Set default harmonic values: 
 //delete gRandom; // TBI this is not safe here, 
 //gRandom = new TRandom3(0);
 Int_t h1 = (Int_t)pow(-1.,1.*gRandom->Integer(fMaxHarmonic+1))*gRandom->Integer(fMaxHarmonic+1); // TBI reimplement all these lines eventually 
 Int_t h2 = (Int_t)pow(-1.,1.*gRandom->Integer(fMaxHarmonic+1))*gRandom->Integer(fMaxHarmonic+1);
 Int_t h3 = (Int_t)pow(-1.,1.*gRandom->Integer(fMaxHarmonic+1))*gRandom->Integer(fMaxHarmonic+1);
 Int_t h4 = (Int_t)pow(-1.,1.*gRandom->Integer(fMaxHarmonic+1))*gRandom->Integer(fMaxHarmonic+1);
 Int_t h5 = (Int_t)pow(-1.,1.*gRandom->Integer(fMaxHarmonic+1))*gRandom->Integer(fMaxHarmonic+1);
 Int_t h6 = (Int_t)pow(-1.,1.*gRandom->Integer(fMaxHarmonic+1))*gRandom->Integer(fMaxHarmonic+1);
 Int_t h7 = (Int_t)pow(-1.,1.*gRandom->Integer(fMaxHarmonic+1))*gRandom->Integer(fMaxHarmonic+1);
 Int_t h8 = (Int_t)pow(-1.,1.*gRandom->Integer(fMaxHarmonic+1))*gRandom->Integer(fMaxHarmonic+1);

 // REMARK: This values can be overriden in a steering macro via 
 // mpc->GetNestedLoopsFlagsPro()->SetBinContent(<binNo>,<value>);

 // b) Book the profile holding all the flags for nested loops:
 fNestedLoopsFlagsPro = new TProfile("fNestedLoopsFlagsPro","Flags for nested loops",10,0,10);
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
 fNestedLoopsFlagsPro->GetXaxis()->SetBinLabel(10,"fCrossCheckDiffWithNestedLoops"); fNestedLoopsFlagsPro->Fill(9.5,fCrossCheckDiffWithNestedLoops);
 fNestedLoopsList->Add(fNestedLoopsFlagsPro);

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
 if(fCrossCheckWithNestedLoops){fNestedLoopsList->Add(fNestedLoopsResultsCosPro);} else{delete fNestedLoopsResultsCosPro;}

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
 if(fCrossCheckWithNestedLoops){fNestedLoopsList->Add(fNestedLoopsResultsSinPro);} else{delete fNestedLoopsResultsSinPro;}

 // e) Book the profile holding all results for differential nested loops:
 fNestedLoopsDiffResultsPro = new TProfile("fNestedLoopsDiffResultsPro","Differential nested loops results",1,0.,1.);
 fNestedLoopsDiffResultsPro->SetStats(kFALSE);
 if(fCrossCheckDiffWithNestedLoops){fNestedLoopsList->Add(fNestedLoopsDiffResultsPro);} else{delete fNestedLoopsDiffResultsPro;}

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForNestedLoops()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForStandardCandles()
{
 // Book all the stuff for 'standard candles'.

 // a) Book the profile holding all the flags for 'standard candles';
 // b) Book the histogram holding all results for 'standard candles';
 // c) Book TProfile2D *fProductsSCPro. 

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForStandardCandles()";

 // a) Book the profile holding all the flags for 'standard candles':
 fStandardCandlesFlagsPro = new TProfile("fStandardCandlesFlagsPro","Flags for standard candles",2,0,2);
 fStandardCandlesFlagsPro->SetTickLength(-0.01,"Y");
 fStandardCandlesFlagsPro->SetMarkerStyle(25);
 fStandardCandlesFlagsPro->SetLabelSize(0.03);
 fStandardCandlesFlagsPro->SetLabelOffset(0.02,"Y");
 fStandardCandlesFlagsPro->SetStats(kFALSE);
 fStandardCandlesFlagsPro->SetFillColor(kGray);
 fStandardCandlesFlagsPro->SetLineColor(kBlack);
 fStandardCandlesFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateStandardCandles"); fStandardCandlesFlagsPro->Fill(0.5,fCalculateStandardCandles);
 fStandardCandlesFlagsPro->GetXaxis()->SetBinLabel(2,"fPropagateErrorSC"); fStandardCandlesFlagsPro->Fill(1.5,fPropagateErrorSC);
 fStandardCandlesList->Add(fStandardCandlesFlagsPro);

 if(!fCalculateStandardCandles){return;} // TBI is this safe like this? 

 // b) Book the histogram holding all results for 'standard candles':
 Int_t nBins = fMaxHarmonic*(fMaxHarmonic-1)/2;
 fStandardCandlesHist = new TH1D("fStandardCandlesHist","Standard candles (SC)",nBins,0.,1.*nBins); 
 fStandardCandlesHist->SetStats(kFALSE);
 fStandardCandlesHist->SetMarkerStyle(kFullSquare);
 fStandardCandlesHist->SetMarkerColor(kBlue);
 fStandardCandlesHist->SetLineColor(kBlue);
 Int_t binNo = 1;
 for(Int_t n1=-fMaxHarmonic;n1<=-2;n1++)
 {
  for(Int_t n2=n1+1;n2<=-1;n2++)
  {
   fStandardCandlesHist->GetXaxis()->SetBinLabel(binNo++,Form("SC(%d,%d,%d,%d)",n1,n2,-1*n2,-1*n1));
  }
 }
 if(binNo-1 != nBins){Fatal(sMethodName.Data(),"Well, binNo-1 != nBins ... :'( ");}
 fStandardCandlesList->Add(fStandardCandlesHist);

 if(!fPropagateErrorSC){return;} // TBI is this safe like this? 

 // c) Book TProfile2D *fProductsSCPro:
 Int_t nBins2D = fMaxHarmonic + fMaxHarmonic*(fMaxHarmonic-1)/2; // #2-p + #4-p distinct correlators in SC context 
 if(nBins2D<=0){Fatal(sMethodName.Data(),"nBins2D<=0");} // well, just in case...
 fProductsSCPro = new TProfile2D("fProductsSCPro","Products of correlations",nBins2D,0.,nBins2D,nBins2D,0.,nBins2D);
 fProductsSCPro->SetStats(kFALSE);
 fProductsSCPro->Sumw2();
 for(Int_t b=1;b<=fMaxHarmonic;b++) // 2-p correlators
 {
  fProductsSCPro->GetXaxis()->SetBinLabel(b,Form("Cos(%d,%d)",-(fMaxHarmonic+1-b),(fMaxHarmonic+1-b)));
  fProductsSCPro->GetYaxis()->SetBinLabel(b,Form("Cos(%d,%d)",-(fMaxHarmonic+1-b),(fMaxHarmonic+1-b)));
 } // for(Int_t b=1;b<=fMaxHarmonic;b++) // 2-p correlators
 for(Int_t b=fMaxHarmonic+1;b<=nBins2D;b++) // 4-p correlators
 {
  TString sBinLabel = fStandardCandlesHist->GetXaxis()->GetBinLabel(b-fMaxHarmonic);
  if(sBinLabel.EqualTo("")){Fatal(sMethodName.Data(),"sBinLabel.EqualTo...");}
  fProductsSCPro->GetXaxis()->SetBinLabel(b,sBinLabel.ReplaceAll("SC","Cos").Data());
  fProductsSCPro->GetYaxis()->SetBinLabel(b,sBinLabel.ReplaceAll("SC","Cos").Data());
 } // for(Int_t b=fMaxHarmonic+1;b<=nBins2D;b++) // 4-p correlators
 fStandardCandlesList->Add(fProductsSCPro);

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForStandardCandles()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetOutputHistograms(TList *histList)
{
 // Get pointers for everything and everywhere from the base list "fHistList". 

 // a) Get pointer for base list fHistList;
 // b) Get pointer for profile holding internal flags and, well, set again all flags;
 // c) Get pointers for control histograms;
 // d) Get pointers for Q-vector;
 // e) Get pointers for correlations;
 // f) Get pointers for 'standard candles';
 // g) Get pointers for Q-cumulants;
 // h) Get pointers for differential correlations;
 // i) Get pointers for symmetry planes.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::GetOutputHistograms(TList *histList)";

 // a) Get pointer for base list fHistList and profile holding internal flags;
 fHistList = histList; 
 if(!fHistList){Fatal(sMethodName.Data(),"fHistList is malicious today...");}

 // b) Get pointer for profile holding internal flags and, well, set again all flags:
 fInternalFlagsPro = dynamic_cast<TProfile*>(fHistList->FindObject("fInternalFlagsPro"));
 if(!fInternalFlagsPro){Fatal(sMethodName.Data(),"fInternalFlagsPro");}
 fUseInternalFlags = fInternalFlagsPro->GetBinContent(1);
 fMinNoRPs = (Int_t)fInternalFlagsPro->GetBinContent(2);
 fMaxNoRPs = (Int_t)fInternalFlagsPro->GetBinContent(3);
 fExactNoRPs = (Int_t)fInternalFlagsPro->GetBinContent(4);
 fPropagateError = (Bool_t)fInternalFlagsPro->GetBinContent(5);

 // c) Get pointers for control histograms:
 this->GetPointersForControlHistograms(); 

 // d) Get pointers for Q-vector:
 this->GetPointersForQvector(); 

 // e) Get pointers for correlations:
 this->GetPointersForCorrelations(); 

 // f) Get pointers for 'standard candles':
 this->GetPointersForStandardCandles(); 

 // g) Get pointers for Q-cumulants:
 this->GetPointersForQcumulants(); 

 // h) Get pointers for differential correlations:
 this->GetPointersForDiffCorrelations(); 

 // i) Get pointers for symmetry planes:
 this->GetPointersForSymmetryPlanes();
  
} // void AliFlowAnalysisWithMultiparticleCorrelations::GetOutputHistograms(TList *histList)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForQvector()
{
 // Get pointers for Q-vector objects.

 // a) Get pointer for fQvectorList; 
 // b) Get pointer for fQvectorFlagsPro;
 // c) Set again all flags.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForQvector()";

 // a) Get pointer for fQvectorList: 
 fQvectorList = dynamic_cast<TList*>(fHistList->FindObject("Q-vectors")); 
 if(!fQvectorList){Fatal(sMethodName.Data(),"fQvectorList");}

 // b) Get pointer for fQvectorFlagsPro: 
 fQvectorFlagsPro = dynamic_cast<TProfile*>(fQvectorList->FindObject("fQvectorFlagsPro"));
 if(!fQvectorFlagsPro){Fatal(sMethodName.Data(),"fQvectorFlagsPro");}

 // c) Set again all flags:
 fCalculateQvector = (Bool_t)fQvectorFlagsPro->GetBinContent(1);
 fCalculateDiffQvectors = (Bool_t)fQvectorFlagsPro->GetBinContent(2);

} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForQvector()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForStandardCandles()
{
 // Get pointers for 'standard candles'.

 // a) Get pointer for fStandardCandlesList; 
 // b) Get pointer for fStandardCandlesFlagsPro;
 // c) Set again all flags; 
 // d) Get pointer TH1D *fStandardCandlesHist;
 // e) Get pointer TProfile2D *fProductsSCPro.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForStandardCandles()";

 // a) Get pointer for fStandardCandlesList: 
 fStandardCandlesList = dynamic_cast<TList*>(fHistList->FindObject("Standard Candles"));
 if(!fStandardCandlesList){Fatal(sMethodName.Data(),"fStandardCandlesList");}

 // b) Get pointer for fStandardCandlesFlagsPro: 
 fStandardCandlesFlagsPro = dynamic_cast<TProfile*>(fStandardCandlesList->FindObject("fStandardCandlesFlagsPro"));
 if(!fStandardCandlesFlagsPro){Fatal(sMethodName.Data(),"fStandardCandlesFlagsPro");}

 // c) Set again all flags: 
 fCalculateStandardCandles = (Bool_t)fStandardCandlesFlagsPro->GetBinContent(1);
 fPropagateErrorSC = (Bool_t)fStandardCandlesFlagsPro->GetBinContent(2);

 if(!fCalculateStandardCandles){return;} // TBI is this safe enough

 // d) Get pointer TH1D *fStandardCandlesHist: 
 fStandardCandlesHist = dynamic_cast<TH1D*>(fStandardCandlesList->FindObject("fStandardCandlesHist"));

 if(!fPropagateErrorSC){return;} // TBI is this safe enough

 // e) Get pointer TProfile2D *fProductsSCPro:
 fProductsSCPro = dynamic_cast<TProfile2D*>(fStandardCandlesList->FindObject("fProductsSCPro"));
 
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

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForControlHistograms()";

 // a) Get pointer for fControlHistogramsList: TBI
 fControlHistogramsList = dynamic_cast<TList*>(fHistList->FindObject("Control Histograms"));
 if(!fControlHistogramsList){Fatal(sMethodName.Data(),"fControlHistogramsList");}

 // b) Get pointer for fControlHistogramsFlagsPro: TBI
 fControlHistogramsFlagsPro = dynamic_cast<TProfile*>(fControlHistogramsList->FindObject("fControlHistogramsFlagsPro"));
 if(!fControlHistogramsFlagsPro){Fatal(sMethodName.Data(),"fControlHistogramsFlagsPro");}

 // c) Set again all flags:
 fFillControlHistograms = fControlHistogramsFlagsPro->GetBinContent(1);
 fFillKinematicsHist = fControlHistogramsFlagsPro->GetBinContent(2);
 fFillMultDistributionsHist = fControlHistogramsFlagsPro->GetBinContent(3);
 fFillMultCorrelationsHist = fControlHistogramsFlagsPro->GetBinContent(4);
 fDontFill[0] = fControlHistogramsFlagsPro->GetBinContent(5);
 fDontFill[1] = fControlHistogramsFlagsPro->GetBinContent(6);
 fDontFill[2] = fControlHistogramsFlagsPro->GetBinContent(7);

 if(!fFillControlHistograms){return;} // TBI is this safe enough

 // d) Get pointers to fKinematicsHist[2][3]: TBI
 TString name[2][3] = {{"RP,phi","RP,pt","RP,eta"},{"POI,phi","POI,pt","POI,eta"}}; // [RP,POI][phi,pt,eta]
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  if(fDontFill[rp]){continue;}
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   fKinematicsHist[rp][ppe] = dynamic_cast<TH1D*>(fControlHistogramsList->FindObject(name[rp][ppe].Data()));
   if(!fKinematicsHist[rp][ppe] && fFillKinematicsHist){Fatal(sMethodName.Data(),"%s",name[rp][ppe].Data());} // TBI 
  }
 }

 // e) Get pointers to TH1D *fMultDistributionsHist[3]:
 TString nameMult[3] = {"Multiplicity (RP)","Multiplicity (POI)","Multiplicity (REF)"}; // [RP,POI,reference multiplicity]
 for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 {
  if(fDontFill[rprm]){continue;}
  fMultDistributionsHist[rprm] = dynamic_cast<TH1D*>(fControlHistogramsList->FindObject(nameMult[rprm].Data()));
  if(!fMultDistributionsHist[rprm] && fFillMultDistributionsHist){Fatal(sMethodName.Data(),"%s",nameMult[rprm].Data());} // TBI 
 } // for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]

 // f) Get pointers to TH2I *fMultCorrelationsHist[3]: TBI automatize the things here (at some point...)...
 if(!fDontFill[0] && !fDontFill[1])
 {
  fMultCorrelationsHist[0] = dynamic_cast<TH2I*>(fControlHistogramsList->FindObject("Multiplicity (RP vs. POI)"));
  if(!fMultCorrelationsHist[0] && fFillMultCorrelationsHist){Fatal(sMethodName.Data(),"Multiplicity (RP vs. POI)");} // TBI 
 }
 if(!fDontFill[0] && !fDontFill[2])
 {
  fMultCorrelationsHist[1] = dynamic_cast<TH2I*>(fControlHistogramsList->FindObject("Multiplicity (RP vs. REF)"));
  if(!fMultCorrelationsHist[1] && fFillMultCorrelationsHist){Fatal(sMethodName.Data(),"Multiplicity (RP vs. REF)");} // TBI 
 }
 if(!fDontFill[1] && !fDontFill[2])
 {
  fMultCorrelationsHist[2] = dynamic_cast<TH2I*>(fControlHistogramsList->FindObject("Multiplicity (POI vs. REF)"));
  if(!fMultCorrelationsHist[2] && fFillMultCorrelationsHist){Fatal(sMethodName.Data(),"Multiplicity (POI vs. REF)");} // TBI 
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForControlHistograms()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForCorrelations()
{
 // Get pointers for correlations.

 // a) Get pointer for fCorrelationsList; TBI
 // b) Get pointer for fCorrelationsFlagsPro; TBI
 // c) Set again all flags; TBI
 // d) Get pointers to TProfile *fCorrelationsPro[2][8].

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForCorrelations()";

 // a) Get pointer for fCorrelationsList: TBI
 fCorrelationsList = dynamic_cast<TList*>(fHistList->FindObject("Correlations"));
 if(!fCorrelationsList){Fatal(sMethodName.Data(),"fCorrelationsList");}

 // b) Get pointer for fCorrelationsFlagsPro: TBI
 fCorrelationsFlagsPro = dynamic_cast<TProfile*>(fCorrelationsList->FindObject("fCorrelationsFlagsPro"));

 if(!fCorrelationsFlagsPro){Fatal(sMethodName.Data(),"fCorrelationsFlagsPro");}

 // c) Set again all flags: 
 fCalculateCorrelations = fCorrelationsFlagsPro->GetBinContent(1);
 fMaxHarmonic = (Int_t)fCorrelationsFlagsPro->GetBinContent(2);
 fMaxCorrelator = (Int_t)fCorrelationsFlagsPro->GetBinContent(3);
 fCalculateIsotropic = (Bool_t)fCorrelationsFlagsPro->GetBinContent(4);
 fCalculateSame = (Bool_t)fCorrelationsFlagsPro->GetBinContent(5);
 fSkipZeroHarmonics = (Bool_t)fCorrelationsFlagsPro->GetBinContent(6);
 fCalculateSameIsotropic = (Bool_t)fCorrelationsFlagsPro->GetBinContent(7);
 fCalculateAll = (Bool_t)fCorrelationsFlagsPro->GetBinContent(8);
 fDontGoBeyond = (Int_t)fCorrelationsFlagsPro->GetBinContent(9);
 fCalculateOnlyForHarmonicQC = (Bool_t)fCorrelationsFlagsPro->GetBinContent(10);
 fCalculateOnlyForSC = (Bool_t)fCorrelationsFlagsPro->GetBinContent(11);
 fCalculateOnlyCos = (Bool_t)fCorrelationsFlagsPro->GetBinContent(12);
 fCalculateOnlySin = (Bool_t)fCorrelationsFlagsPro->GetBinContent(13);

 if(!fCalculateCorrelations){return;} // TBI is this safe enough, that is the question...

 // d) Get pointers to TProfile *fCorrelationsPro[2][8]:
 TString sCosSin[2] = {"Cos","Sin"}; 
 for(Int_t cs=0;cs<2;cs++)
 {
  if(fCalculateOnlyCos && 1==cs){continue;}
  else if(fCalculateOnlySin && 0==cs){continue;}
  for(Int_t c=0;c<fMaxCorrelator;c++)
  {
   if(c==fDontGoBeyond){continue;}
   if(fCalculateOnlyForHarmonicQC && c%2==0){continue;}
   fCorrelationsPro[cs][c] = dynamic_cast<TProfile*>(fCorrelationsList->FindObject(Form("%dpCorrelations%s",c+1,sCosSin[cs].Data())));
   if(!fCorrelationsPro[cs][c]){Fatal(sMethodName.Data(),"%dpCorrelations%s",c+1,sCosSin[cs].Data());} 
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForQcumulants()
{
 // Get pointers for Q-cumulants.

 // a) Get pointer for fQcumulantsList;
 // b) Get pointer for fQcumulantsFlagsPro;
 // c) Set again all flags; 
 // d) Get pointer for fQcumulantsHist;
 // e) Get pointer for fReferenceFlowHist;
 // f) Get pointer for fProductsQCPro.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForQcumulants()";

 // a) Get pointer for fQcumulantsList:
 fQcumulantsList = dynamic_cast<TList*>(fHistList->FindObject("Q-cumulants"));
 if(!fQcumulantsList){Fatal(sMethodName.Data(),"fQcumulantsList");}

 // b) Get pointer for fQcumulantsFlagsPro:
 fQcumulantsFlagsPro = dynamic_cast<TProfile*>(fQcumulantsList->FindObject("fQcumulantsFlagsPro"));
 if(!fQcumulantsFlagsPro){Fatal(sMethodName.Data(),"fQcumulantsFlagsPro");}

 // c) Set again all flags:
 fCalculateQcumulants = (Bool_t) fQcumulantsFlagsPro->GetBinContent(1);
 fHarmonicQC = (Int_t) fQcumulantsFlagsPro->GetBinContent(2);
 fPropagateErrorQC = (Bool_t) fQcumulantsFlagsPro->GetBinContent(3);

 if(!fCalculateQcumulants){return;} // TBI is this safe enough

 // d) Get pointer for fQcumulantsHist:
 fQcumulantsHist = dynamic_cast<TH1D*>(fQcumulantsList->FindObject("Q-cumulants"));
 if(!fQcumulantsHist){Fatal(sMethodName.Data(),"fQcumulantsHist");}

 // e) Get pointer for fReferenceFlowHist:
 fReferenceFlowHist = dynamic_cast<TH1D*>(fQcumulantsList->FindObject("Reference Flow"));
 if(!fReferenceFlowHist){Fatal(sMethodName.Data(),"fReferenceFlowHist");}

 if(!fPropagateErrorQC){return;} // TBI is this safe enough

 // f) Get pointer for fProductsQCPro:
 fProductsQCPro = dynamic_cast<TProfile2D*>(fQcumulantsList->FindObject("fProductsQCPro"));
 if(!fProductsQCPro){Fatal(sMethodName.Data(),"fProductsQCPro");}
 
} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForQcumulants()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForDiffCorrelations()
{
 // Get pointers for differential correlations.

 // a) Get pointer for fDiffCorrelationsList; 
 // b) Get pointer for fDiffCorrelationsFlagsPro; 
 // c) Set again all flags; 
 // d) Get pointers to TProfile *fDiffCorrelationsPro[2][4].

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForDiffCorrelations()";

 // a) Get pointer for fDiffCorrelationsList:
 fDiffCorrelationsList = dynamic_cast<TList*>(fHistList->FindObject("Differential Correlations"));
 if(!fDiffCorrelationsList){Fatal(sMethodName.Data(),"fDiffCorrelationsList");}

 // b) Get pointer for fDiffCorrelationsFlagsPro: 
 fDiffCorrelationsFlagsPro = dynamic_cast<TProfile*>(fDiffCorrelationsList->FindObject("fDiffCorrelationsFlagsPro"));
 if(!fDiffCorrelationsFlagsPro){Fatal(sMethodName.Data(),"fDiffCorrelationsFlagsPro");}

 // c) Set again all flags: 
 fCalculateDiffCorrelations = fDiffCorrelationsFlagsPro->GetBinContent(1);

 if(!fCalculateDiffCorrelations){return;} 

 // TBI get all pointers below for diff. stuff eventually, not needed for the time being. 

 // d) Get pointers to TProfile *fDiffCorrelationsPro[2][4]: // TBI
 /*
 TString sCosSin[2] = {"Cos","Sin"}; 
 for(Int_t cs=0;cs<2;cs++)
 {
  if(fCalculateOnlyCos && 1==cs){continue;}
  else if(fCalculateOnlySin && 0==cs){continue;}
  for(Int_t c=0;c<fMaxCorrelator;c++)
  {
   fCorrelationsPro[cs][c] = dynamic_cast<TProfile*>(fCorrelationsList->FindObject(Form("%dpCorrelations%s",c+1,sCosSin[cs].Data())));
   if(!fCorrelationsPro[cs][c]){Fatal(sMethodName.Data(),"%dpCorrelations%s",c+1,sCosSin[cs].Data());} 
  }
 }
 */

} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForDiffCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForSymmetryPlanes()
{
 // Get pointers for symmetry planes.

 // a) Get pointer for fSymmetryPlanesList;
 // b) Get pointer for fSymmetryPlanesFlagsPro;
 // c) Set again all flags;
 // d) Get pointers to TProfile *fSymmetryPlanesPro[?][?].

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForSymmetryPlanes()";

 // a) Get pointer for fSymmetryPlanesList:
 fSymmetryPlanesList = dynamic_cast<TList*>(fHistList->FindObject("Symmetry_Plane_Correlations"));
 if(!fSymmetryPlanesList){Fatal(sMethodName.Data(),"fSymmetryPlanesList");}

 // b) Get pointer for fSymmetryPlanesFlagsPro:
 fSymmetryPlanesFlagsPro = dynamic_cast<TProfile*>(fSymmetryPlanesList->FindObject("fSymmetryPlanesFlagsPro"));
 if(!fSymmetryPlanesFlagsPro){Fatal(sMethodName.Data(),"fSymmetryPlanesFlagsPro");}

 // c) Set again all flags:
 fCalculateSymmetryPlanes = fSymmetryPlanesFlagsPro->GetBinContent(1);

 if(!fCalculateSymmetryPlanes){return;}

 // d) Get pointers to TProfile *fSymmetryPlanesPro[?][?]: // TBI
 for(Int_t gc=0;gc<1;gc++) // 'generic correlator': [[0]:(Psi2n,Psi1n),[1]:...] TBI upper boundary will change
 {
  for(Int_t n=0;n<2;n++) // 'harmonic n for generic correlator': [[0]:n=1,[1]:n=2,...] TBI upper boundary will change
  {
   fSymmetryPlanesPro[gc][n] = dynamic_cast<TProfile*>(fSymmetryPlanesList->FindObject(Form("%d,%d",gc,n)));
   if(!fSymmetryPlanesPro[gc][n]){Fatal(sMethodName.Data(),"fSymmetryPlanesPro[%d][%d]",gc,n);}
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForSymmetryPlanes()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForQvector()
{
 // Initialize all arrays for Q-vector.

 for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++) 
 {
  for(Int_t wp=0;wp<fMaxCorrelator+1;wp++) // weight power
  {
   fQvector[h][wp] = TComplex(0.,0.);
   for(Int_t b=0;b<100;b++) // TBI hardwired 100 
   {  
    fpvector[b][h][wp] = TComplex(0.,0.); 
    fqvector[b][h][wp] = TComplex(0.,0.); 
   }
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForQvector()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::ResetQvector()
{
 // Reset all Q-vector components to zero before starting a new event. 

 for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++) 
 {
  for(Int_t wp=0;wp<fMaxCorrelator+1;wp++) // weight powe
  {
   fQvector[h][wp] = TComplex(0.,0.);
   if(!fCalculateDiffQvectors){continue;}
   for(Int_t b=0;b<100;b++) // TBI hardwired 100 
   {  
    fpvector[b][h][wp] = TComplex(0.,0.); 
    fqvector[b][h][wp] = TComplex(0.,0.); 
   }
  } 
 } 

} // void AliFlowAnalysisWithMultiparticleCorrelations::ResetQvector()

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Q(Int_t n, Int_t wp)
{
 // Using the fact that Q{-n,p} = Q{n,p}^*. 
 
 if(n>=0){return fQvector[n][wp];} 
 return TComplex::Conjugate(fQvector[-n][wp]);
 
} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Q(Int_t n, Int_t wp)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::p(Int_t n, Int_t wp)
{
 // Using the fact that p{-n,p} = p{n,p}^*.
 
 if(n>=0){return fpvector[fDiffBinNo][n][wp];} 
 return TComplex::Conjugate(fpvector[fDiffBinNo][-n][wp]);

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::p(Int_t n, Int_t p)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::q(Int_t n, Int_t wp)
{
 // Using the fact that q{-n,p} = q{n,p}^*.

 // When weights are used for RPs and not for POIs, and vice versa, some gymnastics is required here:
 // TBI rethink and reimplement the lines below:
 Int_t nUseWeightsForRP = (Int_t)(fUseWeights[0][0] || fUseWeights[0][1] || fUseWeights[0][2]); 
 Int_t nUseWeightsForPOI = (Int_t)(fUseWeights[1][0] || fUseWeights[1][1] || fUseWeights[1][2]); 
 if(nUseWeightsForPOI == 1 && nUseWeightsForRP == 0){wp=1;}
 else if(nUseWeightsForPOI == 0 && nUseWeightsForRP == 1){wp-=1;}

 if(n>=0){return fqvector[fDiffBinNo][n][wp];} 
 return TComplex::Conjugate(fqvector[fDiffBinNo][-n][wp]);

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::q(Int_t n, Int_t wp)

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

 Int_t harmonic[7] = {n1,n2,n3,n4,n5,n6,n7};

 TComplex seven = Recursion(7,harmonic); 

 return seven;

} // end of TComplex AliFlowAnalysisWithMultiparticleCorrelations::Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)
{
 // Generic eight-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8)]>.

 Int_t harmonic[8] = {n1,n2,n3,n4,n5,n6,n7,n8};

 TComplex eight = Recursion(8,harmonic); 

 return eight;

} // end of TComplex AliFlowAnalysisWithMultiparticleCorrelations::Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForWeights()
{
 // Book all objects for calculations with weights. 

 // a) Book profile to hold all flags for weights;
 // b) Store histograms holding phi, pt and eta weights.
    
 // a) Book profile to hold all flags for weights:
 fWeightsFlagsPro = new TProfile("fWeightsFlagsPro","0 = weight not used, 1 = weight used ",6,0,6);
 fWeightsFlagsPro->SetLabelSize(0.06);
 fWeightsFlagsPro->SetStats(kFALSE);
 fWeightsFlagsPro->SetFillColor(kGray);
 fWeightsFlagsPro->SetLineColor(kBlack);
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(1,"RP: w_{#phi}"); fWeightsFlagsPro->Fill(0.5,fUseWeights[0][0]);
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(2,"RP: w_{p_{T}}"); fWeightsFlagsPro->Fill(1.5,fUseWeights[0][1]);
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(3,"RP: w_{#eta}"); fWeightsFlagsPro->Fill(2.5,fUseWeights[0][2]);
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(4,"POI: w_{#phi}"); fWeightsFlagsPro->Fill(3.5,fUseWeights[1][0]);
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(5,"POI: w_{p_{T}}"); fWeightsFlagsPro->Fill(4.5,fUseWeights[1][1]);
 fWeightsFlagsPro->GetXaxis()->SetBinLabel(6,"POI: w_{#eta}"); fWeightsFlagsPro->Fill(5.5,fUseWeights[1][2]);
 fWeightsList->Add(fWeightsFlagsPro); 
  
 // b) Store histograms holding phi, pt and eta weights:
 //    REMARK: It is assumed that these histos are accessed from external file "weights.root" 
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   if(fWeightsHist[rp][ppe]){fWeightsList->Add(fWeightsHist[rp][ppe]);}
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForWeights()

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::Weight(const Double_t &value, const char *type, const char *variable) // value, [RP,POI], [phi,pt,eta]
{
 // Determine particle weight. 

 TString sMethodName = "Double_t AliFlowAnalysisWithMultiparticleCorrelations::Weight(const Double_t &value, const char *type, const char *variable)"; 

 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI"))){Fatal(sMethodName.Data(),"!(TString(type).EqualTo...");}
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))){Fatal(sMethodName.Data(),"!(TString(variable).EqualTo...");}

 Int_t rp = 0; // [RP,POI]
 if(TString(type).EqualTo("POI")){rp=1;} 

 Int_t ppe = 0; // [phi,pt,eta]
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 if(!fWeightsHist[rp][ppe]){Fatal(sMethodName.Data(),"!fWeightsHist[rp][ppe]");}

 Double_t weight = fWeightsHist[rp][ppe]->GetBinContent(fWeightsHist[rp][ppe]->FindBin(value));

 return weight;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::Weight(const Double_t &value, const char *type, const char *variable)

//=======================================================================================================================

/*
Double_t AliFlowAnalysisWithMultiparticleCorrelations::PhiWeight(const Double_t &dPhi, const char *type)
{
 // Determine phi weight for a given phi. 

 TString sMethodName = "Double_t AliFlowAnalysisWithMultiparticleCorrelations::PhiWeight(const Double_t &dPhi, const char *type)"; 

 // Basic protection:
 if(!(TString(type)::EqualTo("RP") || TString(type)::EqualTo("POI"))){Fatal(sMethodName.Data(),"!(TString(type)::EqualTo...");}

 Int_t rp = 0; // RP or POI
 if(TString(type)::EqualTo("POI")){rp=1;} 
 
 if(!fWeightsHist[rp][0]){Fatal("AliFlowAnalysisWithMultiparticleCorrelations::PhiWeight(const Double_t &dPhi)","fPhiWeightsHist");}
  








 Double_t wPhi = fPhiWeightsHist->GetBinContent(fPhiWeightsHist->FindBin(dPhi));
 

 wPhi = 0.;

 return wPhi;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::PhiWeight(const Double_t &dPhi)

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::PtWeight(const Double_t &dPt, const char *type)
{
 // Determine pt weight for a given pt. 

 if(!fPtWeightsHist){Fatal("AliFlowAnalysisWithMultiparticleCorrelations::PtWeight(const Double_t &dPt)","fPtWeightsHist");}

 Double_t wPt = fPtWeightsHist->GetBinContent(fPtWeightsHist->FindBin(dPt));

 return wPt;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::PtWeight(const Double_t &dPt)

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::EtaWeight(const Double_t &dEta, const char *type)
{
 // Determine eta weight for a given eta. 

 if(!fEtaWeightsHist){Fatal("AliFlowAnalysisWithMultiparticleCorrelations::EtaWeight(const Double_t &dEta)","fEtaWeightsHist");}

 Double_t wEta = fEtaWeightsHist->GetBinContent(fEtaWeightsHist->FindBin(dEta));

 return wEta;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::EtaWeight(const Double_t &dEta)

*/


//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForBase()
{
 // Book all base objects. 

 fInternalFlagsPro = new TProfile("fInternalFlagsPro","Internal flags and settings",10,0,10);
 fInternalFlagsPro->SetLabelSize(0.05);
 fInternalFlagsPro->SetStats(kFALSE);
 fInternalFlagsPro->SetFillColor(kGray);
 fInternalFlagsPro->SetLineColor(kBlack);
 fInternalFlagsPro->GetXaxis()->SetBinLabel(1,"fUseInternalFlags"); fInternalFlagsPro->Fill(0.5,fUseInternalFlags);  
 fInternalFlagsPro->GetXaxis()->SetBinLabel(2,"fMinNoRPs"); fInternalFlagsPro->Fill(1.5,fMinNoRPs);  
 fInternalFlagsPro->GetXaxis()->SetBinLabel(3,"fMaxNoRPs"); fInternalFlagsPro->Fill(2.5,fMaxNoRPs); 
 fInternalFlagsPro->GetXaxis()->SetBinLabel(4,"fExactNoRPs"); fInternalFlagsPro->Fill(3.5,fExactNoRPs);  
 fInternalFlagsPro->GetXaxis()->SetBinLabel(5,"fPropagateError"); fInternalFlagsPro->Fill(4.5,fPropagateError);  
 fInternalFlagsPro->GetXaxis()->SetBinLabel(6,Form("fAnalysisTag = %s",fAnalysisTag.Data())); 
 fInternalFlagsPro->GetXaxis()->SetBinLabel(7,"fDumpThePoints"); fInternalFlagsPro->Fill(6.5,fDumpThePoints);  
 fInternalFlagsPro->GetXaxis()->SetBinLabel(8,"fMaxNoEventsPerFile"); fInternalFlagsPro->Fill(7.5,fMaxNoEventsPerFile);  
 fInternalFlagsPro->GetXaxis()->SetBinLabel(9,"fSelectRandomlyRPs"); fInternalFlagsPro->Fill(8.5,fSelectRandomlyRPs);  
 fInternalFlagsPro->GetXaxis()->SetBinLabel(10,"fnSelectedRandomlyRPs"); fInternalFlagsPro->Fill(9.5,fnSelectedRandomlyRPs);  
 fHistList->Add(fInternalFlagsPro); 

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForBase()

//=======================================================================================================================

Bool_t AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckInternalFlags(AliFlowEventSimple *anEvent)
{
 // Cross-check in this method wether "anEvent" passes internal flags. 

 // a) Cross-check min. and max. number of RPs; 
 // b) Cross-check exact number of RPs. 

 Bool_t bPassesInternalFlags = kTRUE;

 // a) Cross-check min. and max. number of RPs: 
 if(-44 != fMinNoRPs)
 {
  fMinNoRPs <= anEvent->GetNumberOfRPs() ? 1 : bPassesInternalFlags = kFALSE; 
  if(!bPassesInternalFlags){return bPassesInternalFlags;}
 }
 if(-44 != fMaxNoRPs)
 {
  anEvent->GetNumberOfRPs() < fMaxNoRPs ? 1 : bPassesInternalFlags = kFALSE;  
  if(!bPassesInternalFlags){return bPassesInternalFlags;}
 }

 // b) Cross-check exact number of RPs:
 if(-44 != fExactNoRPs)
 {
  anEvent->GetNumberOfRPs() == fExactNoRPs ? 1 : bPassesInternalFlags = kFALSE;  
  if(!bPassesInternalFlags){return bPassesInternalFlags;}
 }

 return bPassesInternalFlags; 

} // Bool_t AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckInternalFlags(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::DumpPointsForDurham(TGraphErrors *ge)
{
 // Dump points from TGraphErrors object into Durham database format. 

 // Remark 1: format is <binCenter>  <value>  +-<stat.error>
 // Remark 2: the default precision is 6 significant digits

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::DumpPointsForDurham(TGraphErrors *ge)"; 

 if(!ge){Fatal(sMethodName.Data(),"ge is NULL, for one reason or another...");}

 ge->Draw("ap");

 Int_t nPoints = ge->GetN();
 Double_t x = 0.;
 //Double_t xErr = 0.;
 Double_t y = 0.;
 Double_t yErr = 0.;
 printf("\nbinCenter value     +-stat.error\n");
 for(Int_t p=0;p<nPoints;p++)
 { 
  ge->GetPoint(p,x,y);
  //xErr = ge->GetErrorX(p); 
  yErr = ge->GetErrorY(p); 
  printf("%f  %f  +-%f\n",x,y,yErr);
 } // end of for(Int_t p=0;p<nPoints;p++)
 cout<<endl;

} // void AliFlowAnalysisWithMultiparticleCorrelations::DumpPointsForDurham(TGraphErrors *ge)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::DumpPointsForDurham(TH1D *h)
{
 // Dump points from TH1D object into Durham database format.

 // Remark 1: format is <binCenter>  <value>  +-<stat.error>
 // Remark 2: the default precision is 6 significant digits

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::DumpPointsForDurham(TH1D *h)"; 

 if(!h){Fatal(sMethodName.Data(),"h is NULL, for one reason or another...");}

 h->Draw();

 Int_t nPoints = h->GetXaxis()->GetNbins();
 Double_t x = 0.;
 Double_t y = 0.;
 Double_t yErr = 0.;
 printf("\nbinCenter value     +-stat.error\n");
 for(Int_t p=1;p<=nPoints;p++)
 { 
  x = h->GetBinCenter(p);
  y = h->GetBinContent(p);
  yErr = h->GetBinError(p); 
  //printf("%f  %f  +-%f\n",x,y,yErr); 
  printf("%e  %e  +-%e\n",x,y,yErr); 
 } // end of for(Int_t p=0;p<nPoints;p++)
 cout<<endl;

} // void AliFlowAnalysisWithMultiparticleCorrelations::DumpPointsForDurham(TH1D *h)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::DumpPointsForDurham(TH1F *h)
{
 // Dump points from TH1F object into Durham database format.

 // Remark 1: format is <binCenter>  <value>  +-<stat.error>
 // Remark 2: the default precision is 6 significant digits

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::DumpPointsForDurham(TH1F *h)"; 

 if(!h){Fatal(sMethodName.Data(),"h is NULL, for one reason or another...");}

 h->Draw();

 Int_t nPoints = h->GetXaxis()->GetNbins();
 Double_t x = 0.;
 Double_t y = 0.;
 Double_t yErr = 0.;
 printf("\nbinCenter value     +-stat.error\n");
 for(Int_t p=1;p<=nPoints;p++)
 { 
  x = h->GetBinCenter(p);
  y = h->GetBinContent(p);
  yErr = h->GetBinError(p); 
  printf("%f  %f  +-%f\n",x,y,yErr);
 } // end of for(Int_t p=0;p<nPoints;p++)
 cout<<endl;

} // void AliFlowAnalysisWithMultiparticleCorrelations::DumpPointsForDurham(TH1F *h)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForQcumulants()
{
 // Initialize all arrays for Q-cumulants.

 // ... TBI 

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForQcumulants()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForQcumulants()
{
 // Book all the stuff for Q-cumulants.

 // a) Book the profile holding all the flags for Q-cumulants;
 // b) Book TH1D *fQcumulantsHist;
 // c) Book TH1D *fReferenceFlowHist;
 // d) Book TProfile2D *fProductsQCPro.

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForQcumulants()";

 // a) Book the profile holding all the flags for Q-cumulants:
 fQcumulantsFlagsPro = new TProfile("fQcumulantsFlagsPro","Flags for Q-cumulants",3,0,3);
 fQcumulantsFlagsPro->SetTickLength(-0.01,"Y");
 fQcumulantsFlagsPro->SetMarkerStyle(25);
 fQcumulantsFlagsPro->SetLabelSize(0.03);
 fQcumulantsFlagsPro->SetLabelOffset(0.02,"Y");
 fQcumulantsFlagsPro->SetStats(kFALSE);
 fQcumulantsFlagsPro->SetFillColor(kGray);
 fQcumulantsFlagsPro->SetLineColor(kBlack);
 fQcumulantsFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateQcumulants"); fQcumulantsFlagsPro->Fill(0.5,fCalculateQcumulants); 
 fQcumulantsFlagsPro->GetXaxis()->SetBinLabel(2,"fHarmonicQC"); fQcumulantsFlagsPro->Fill(1.5,fHarmonicQC); 
 fQcumulantsFlagsPro->GetXaxis()->SetBinLabel(3,"fPropagateErrorQC"); fQcumulantsFlagsPro->Fill(2.5,fPropagateErrorQC); 
 fQcumulantsList->Add(fQcumulantsFlagsPro);

 if(!fCalculateQcumulants){return;} // TBI is this safe enough? 

 // b) Book TH1D *fQcumulantsHist:
 fQcumulantsHist = new TH1D("Q-cumulants",Form("Q-cumulants (n = %d)",fHarmonicQC),4,0.,4.);
 fQcumulantsHist->SetStats(kFALSE);
 fQcumulantsHist->SetMarkerColor(kBlack);
 fQcumulantsHist->SetMarkerStyle(kFullSquare);
 fQcumulantsHist->GetXaxis()->SetLabelSize(0.045);
 fQcumulantsHist->GetXaxis()->SetLabelOffset(0.01);
 for(Int_t qc=1;qc<=4;qc++) // [QC{2},QC{4},QC{6},QC{8}]
 {
  fQcumulantsHist->GetXaxis()->SetBinLabel(qc,Form("QC{%d}",2*qc));
 } 
 fQcumulantsList->Add(fQcumulantsHist);

 // c) Book TH1D *fReferenceFlowHist:
 fReferenceFlowHist = new TH1D("Reference Flow","Reference flow from Q-cumulants",4,0.,4.);
 fReferenceFlowHist->SetStats(kFALSE);
 fReferenceFlowHist->SetMarkerColor(kBlack);
 fReferenceFlowHist->SetMarkerStyle(kFullSquare);
 fReferenceFlowHist->GetXaxis()->SetLabelSize(0.045);
 fReferenceFlowHist->GetXaxis()->SetLabelOffset(0.01);
 for(Int_t qc=1;qc<=4;qc++) // [vn{2},vn{4},vn{6},vn{8}]
 {
  fReferenceFlowHist->GetXaxis()->SetBinLabel(qc,Form("v_{%d}{%d}",fHarmonicQC,2*qc));
 } 
 fQcumulantsList->Add(fReferenceFlowHist);

 if(!fPropagateErrorQC){return;} // TBI is this safe enough? 

 // d) Book TProfile2D *fProductsQCPro:
 const Int_t nCorrelations = 4;
 Int_t n = fHarmonicQC; 
 TString sCorrelations[nCorrelations] = {Form("Cos(-%d,%d)",n,n),
                                         Form("Cos(-%d,-%d,%d,%d)",n,n,n,n),
                                         Form("Cos(-%d,-%d,-%d,%d,%d,%d)",n,n,n,n,n,n),
                                         Form("Cos(-%d,-%d,-%d,-%d,%d,%d,%d,%d)",n,n,n,n,n,n,n,n)}; 
 Int_t nBins2D = (Int_t)TMath::Floor(fDontGoBeyond/2.);
 if(fDontGoBeyond > 8){nBins2D = 4;}
 if(nBins2D < 1 || nBins2D > 4)
 {
  cout<<Form("nBins2D = %d",nBins2D)<<endl;
  cout<<Form("fDontGoBeyond = %d",fDontGoBeyond)<<endl;
  Fatal(sMethodName.Data(),"nBins2D < 1 || nBins2D > 4");
 }
 fProductsQCPro = new TProfile2D("fProductsQCPro","Products of correlations",nBins2D,0.,nBins2D,nBins2D,0.,nBins2D);
 fProductsQCPro->SetStats(kFALSE);
 fProductsQCPro->Sumw2();
 for(Int_t b=1;b<=nBins2D;b++)
 {
  fProductsQCPro->GetXaxis()->SetBinLabel(b,sCorrelations[b-1].Data());
  fProductsQCPro->GetYaxis()->SetBinLabel(b,sCorrelations[b-1].Data());
 } // for(Int_t b=1;b<=nBins2D;b++)
 fQcumulantsList->Add(fProductsQCPro);
 
} // end of void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForQcumulants()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateQcumulants()
{
 // Calculate Q-cumulants.

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::CalculateQcumulants()";

 fPropagateError = kTRUE;
 Int_t n = fHarmonicQC;
 fQcumulantsHist->SetTitle(Form("Q-cumulants (n = %d)",n));

 TString sCorrelations[4] = {Form("Cos(-%d,%d)",n,n),
                             Form("Cos(-%d,-%d,%d,%d)",n,n,n,n),
                             Form("Cos(-%d,-%d,-%d,%d,%d,%d)",n,n,n,n,n,n),
                             Form("Cos(-%d,-%d,-%d,-%d,%d,%d,%d,%d)",n,n,n,n,n,n,n,n)};

 Int_t nBins[4] = {fCorrelationsPro[0][1]->GetNbinsX(),fCorrelationsPro[0][3]->GetNbinsX(),
                   fCorrelationsPro[0][5]->GetNbinsX(),fCorrelationsPro[0][7]->GetNbinsX()};

 Double_t dCorrelation[4] = {0.};
 Double_t dCorrelationErr[4] = {0.};

 for(Int_t c=0;c<4;c++) // [<<2>>,<<4>>,<<6>>,<<8>>]
 {
  if(2*(c+1)>fDontGoBeyond){break;}
  for(Int_t b=1;b<=nBins[c];b++)
  {
   if(sCorrelations[c].EqualTo(fCorrelationsPro[0][2*c+1]->GetXaxis()->GetBinLabel(b)))   
   {
    dCorrelation[c] = fCorrelationsPro[0][2*c+1]->GetBinContent(b);
    dCorrelationErr[c] = fCorrelationsPro[0][2*c+1]->GetBinError(b);
    break; 
   }
  } // for(Int_t b=1;b<=nBins[c];b++)
 } // for(Int_t c=0;c<4;c++) // [<<2>>,<<4>>,<<6>>,<<8>>]

 // Correlations:
 Double_t two = dCorrelation[0]; // <<2>>
 Double_t four = dCorrelation[1]; // <<4>>  
 Double_t six = dCorrelation[2]; // <<6>> 
 Double_t eight = dCorrelation[3]; // <<8>>  

 // Q-cumulants: 
 Double_t qc2 = 0.; // QC{2}
 Double_t qc4 = 0.; // QC{4}
 Double_t qc6 = 0.; // QC{6}
 Double_t qc8 = 0.; // QC{8}
 if(TMath::Abs(two) > 0. && !(fDontGoBeyond < 2)){qc2 = two;} 
 if(TMath::Abs(four) > 0. && !(fDontGoBeyond < 4)){qc4 = four-2.*pow(two,2.);} 
 if(TMath::Abs(six) > 0. && !(fDontGoBeyond < 6)){qc6 = six-9.*two*four+12.*pow(two,3.);} 
 if(TMath::Abs(eight) > 0. && !(fDontGoBeyond < 8)){qc8 = eight-16.*two*six-18.*pow(four,2.)+144.*pow(two,2.)*four-144.*pow(two,4.);} 

 // Store the results for Q-cumulants:
 if(TMath::Abs(qc2)>0.)
 {
  fQcumulantsHist->SetBinContent(1,qc2);
 }
 if(TMath::Abs(qc4)>0.)
 {
  fQcumulantsHist->SetBinContent(2,qc4);
 }
 if(TMath::Abs(qc6)>0.)
 {
  fQcumulantsHist->SetBinContent(3,qc6);
 }
 if(TMath::Abs(qc8)>0.)
 {
  fQcumulantsHist->SetBinContent(4,qc8); 
 }

 if(!fPropagateErrorQC)
 {
  fQcumulantsHist->SetBinError(1,0.); 
  fQcumulantsHist->SetBinError(2,0.); 
  fQcumulantsHist->SetBinError(3,0.); 
  fQcumulantsHist->SetBinError(4,0.); 
  return;
 } // if(!fPropagateErrorQC)

 // Statistical errors of average 2-, 4-, 6- and 8-particle azimuthal correlations:
 Double_t twoError = dCorrelationErr[0]; // statistical error of <2>  
 Double_t fourError = dCorrelationErr[1]; // statistical error of <4>   
 Double_t sixError = dCorrelationErr[2]; // statistical error of <6> 
 Double_t eightError = dCorrelationErr[3]; // statistical error of <8> 

 // Covariances multiplied by a prefactor depending on weights: 
 Double_t wCov24 = 0.; // Cov(<2>,<4>) * prefactor(w_<2>,w_<4>) 
 Double_t wCov26 = 0.; // Cov(<2>,<6>) * prefactor(w_<2>,w_<6>)
 Double_t wCov28 = 0.; // Cov(<2>,<8>) * prefactor(w_<2>,w_<8>)
 Double_t wCov46 = 0.; // Cov(<4>,<6>) * prefactor(w_<4>,w_<6>)
 Double_t wCov48 = 0.; // Cov(<4>,<8>) * prefactor(w_<4>,w_<8>)
 Double_t wCov68 = 0.; // Cov(<6>,<8>) * prefactor(w_<6>,w_<8>)  
 if(!(fDontGoBeyond < 4)){wCov24 = Covariance(sCorrelations[0].Data(),sCorrelations[1].Data(),fProductsQCPro);} // Cov(<2>,<4>) * prefactor(w_<2>,w_<4>) 
 if(!(fDontGoBeyond < 6)){wCov26 = Covariance(sCorrelations[0].Data(),sCorrelations[2].Data(),fProductsQCPro);} // Cov(<2>,<6>) * prefactor(w_<2>,w_<6>)
 if(!(fDontGoBeyond < 8)){wCov28 = Covariance(sCorrelations[0].Data(),sCorrelations[3].Data(),fProductsQCPro);} // Cov(<2>,<8>) * prefactor(w_<2>,w_<8>)
 if(!(fDontGoBeyond < 6)){wCov46 = Covariance(sCorrelations[1].Data(),sCorrelations[2].Data(),fProductsQCPro);} // Cov(<4>,<6>) * prefactor(w_<4>,w_<6>)
 if(!(fDontGoBeyond < 8)){wCov48 = Covariance(sCorrelations[1].Data(),sCorrelations[3].Data(),fProductsQCPro);} // Cov(<4>,<8>) * prefactor(w_<4>,w_<8>)
 if(!(fDontGoBeyond < 8)){wCov68 = Covariance(sCorrelations[2].Data(),sCorrelations[3].Data(),fProductsQCPro);} // Cov(<6>,<8>) * prefactor(w_<6>,w_<8>)  

 // Statistical errors of Q-cumulants:       
 Double_t qc2Error = 0.;
 Double_t qc4Error = 0.;
 Double_t qc6Error = 0.;
 Double_t qc8Error = 0.; 
 // Squared statistical errors of Q-cumulants:       
 //Double_t qc2ErrorSquared = 0.;
 Double_t qc4ErrorSquared = 0.;
 Double_t qc6ErrorSquared = 0.;
 Double_t qc8ErrorSquared = 0.;        
 // Statistical error of QC{2}:              
 if(!(fDontGoBeyond < 2)){qc2Error = twoError;}                                                 
 // Statistical error of QC{4}:              
 qc4ErrorSquared = 16.*pow(two,2.)*pow(twoError,2)+pow(fourError,2.)
                 - 8.*two*wCov24;                     
 if(qc4ErrorSquared > 0. && !(fDontGoBeyond < 4))
 {
  qc4Error = pow(qc4ErrorSquared,0.5);
 } else{Warning(sMethodName.Data(),"Statistical error of QC{4} is imaginary !!!!"); fPropagateError = kFALSE;}
                                           
 // Statistical error of QC{6}:              
 qc6ErrorSquared = 81.*pow(4.*pow(two,2.)-four,2.)*pow(twoError,2.)
                 + 81.*pow(two,2.)*pow(fourError,2.)
                 + pow(sixError,2.)
                 - 162.*two*(4.*pow(two,2.)-four)*wCov24
                 + 18.*(4.*pow(two,2.)-four)*wCov26
                 - 18.*two*wCov46;                     
 if(qc6ErrorSquared > 0. && !(fDontGoBeyond < 6))
 {
  qc6Error = pow(qc6ErrorSquared,0.5);
 } else{Warning(sMethodName.Data(),"Statistical error of QC{6} is imaginary !!!!"); fPropagateError = kFALSE;}

 // Statistical error of QC{8}:              
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
 if(qc8ErrorSquared > 0. && !(fDontGoBeyond < 8))
 {
  qc8Error = pow(qc8ErrorSquared,0.5);
 } else{Warning(sMethodName.Data(),"Statistical error of QC{8} is imaginary !!!!"); fPropagateError = kFALSE;}

 if(!fPropagateError){fPropagateError = kTRUE; return;}

 // Store the statistical errors for Q-cumulants:
 if(TMath::Abs(qc2)>0.)
 {
  fQcumulantsHist->SetBinError(1,qc2Error);
 }
 if(TMath::Abs(qc4)>0.)
 {
  fQcumulantsHist->SetBinError(2,qc4Error);
 }
 if(TMath::Abs(qc6)>0.)
 {
  fQcumulantsHist->SetBinError(3,qc6Error);
 }
 if(TMath::Abs(qc8)>0.)
 {
  fQcumulantsHist->SetBinError(4,qc8Error);
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateQcumulants()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateReferenceFlow()
{
 // Calculate reference flow from Q-cumulants.

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::CalculateReferenceFlow()";

 Int_t n = fHarmonicQC;

 // Reference flow estimates:
 Double_t v2 = 0.; // v{2}  
 Double_t v4 = 0.; // v{4}  
 Double_t v6 = 0.; // v{6}  
 Double_t v8 = 0.; // v{8}

 // Reference flow's statistical errors:
 Double_t v2Error = 0.; // v{2} stat. error 
 Double_t v4Error = 0.; // v{4} stat. error
 Double_t v6Error = 0.; // v{6} stat. error
 Double_t v8Error = 0.; // v{8} stat. error
  
 // Q-cumulants:
 Double_t qc2 = fQcumulantsHist->GetBinContent(1); // QC{2}  
 Double_t qc4 = fQcumulantsHist->GetBinContent(2); // QC{4}  
 Double_t qc6 = fQcumulantsHist->GetBinContent(3); // QC{6}  
 Double_t qc8 = fQcumulantsHist->GetBinContent(4); // QC{8}

 // Q-cumulants's statistical errors: 
 Double_t qc2Error = fQcumulantsHist->GetBinError(1); // QC{2} stat. error  
 Double_t qc4Error = fQcumulantsHist->GetBinError(2); // QC{4} stat. error  
 Double_t qc6Error = fQcumulantsHist->GetBinError(3); // QC{6} stat. error  
 Double_t qc8Error = fQcumulantsHist->GetBinError(4); // QC{8} stat. error

 // Calculate reference flow estimates from Q-cumulants: 
 if(qc2>=0.){v2 = pow(qc2,0.5);} 
 if(qc4<=0.){v4 = pow(-1.*qc4,1./4.);} 
 if(qc6>=0.){v6 = pow((1./4.)*qc6,1./6.);}
 if(qc8<=0.){v8 = pow((-1./33.)*qc8,1./8.);}  

 // Calculate stat. error for reference flow estimates from stat. error of Q-cumulants:  
 if(qc2>0. && qc2Error>0.){v2Error = (1./2.)*pow(qc2,-0.5)*qc2Error;} 
 if(qc4<0. && qc4Error>0.){v4Error = (1./4.)*pow(-qc4,-3./4.)*qc4Error;} 
 if(qc6>0. && qc6Error>0.){v6Error = (1./6.)*pow(2.,-1./3.)*pow(qc6,-5./6.)*qc6Error;}   
 if(qc8<0. && qc8Error>0.){v8Error = (1./8.)*pow(33.,-1./8.)*pow(-qc8,-7./8.)*qc8Error;}   

 // Print warnings for the 'wrong sign' cumulants: 
 if(TMath::Abs(v2)<1.e-44){Warning(sMethodName.Data(),"Wrong sign QC{2}, couldn't calculate v{2} !!!!");}
 if(TMath::Abs(v4)<1.e-44){Warning(sMethodName.Data(),"Wrong sign QC{4}, couldn't calculate v{4} !!!!");} 
 if(TMath::Abs(v6)<1.e-44){Warning(sMethodName.Data(),"Wrong sign QC{6}, couldn't calculate v{6} !!!!");}
 if(TMath::Abs(v8)<1.e-44){Warning(sMethodName.Data(),"Wrong sign QC{8}, couldn't calculate v{8} !!!!");}                       

 // Store the results and statistical errors of reference flow estimates:
 for(Int_t qc=1;qc<=4;qc++) // [vn{2},vn{4},vn{6},vn{8}]
 {
  fReferenceFlowHist->GetXaxis()->SetBinLabel(qc,Form("v_{%d}{%d}",n,2*qc));
 } 
 fReferenceFlowHist->SetBinContent(1,v2);
 fReferenceFlowHist->SetBinError(1,v2Error);
 fReferenceFlowHist->SetBinContent(2,v4);
 fReferenceFlowHist->SetBinError(2,v4Error);
 fReferenceFlowHist->SetBinContent(3,v6);
 fReferenceFlowHist->SetBinError(3,v6Error);
 fReferenceFlowHist->SetBinContent(4,v8);
 fReferenceFlowHist->SetBinError(4,v8Error);  

 // Final printout:
 cout<<endl;
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<" flow estimates from Q-cumulants"<<endl;
 TString sWhichWeights = "no weights";
 if(fUseWeights[0][0]){sWhichWeights = "phi weights";} 
 else if(fUseWeights[0][1]){sWhichWeights = "pt weights";}
 else if(fUseWeights[0][2]){sWhichWeights = "eta weights";}
 cout<<Form("  (MPC class, RPs, %s)",sWhichWeights.Data())<<endl; 
 cout<<endl;
 for(Int_t co=0;co<4;co++) // cumulant order
 {
  cout<<Form("  v_%d{%d} = %.8f +/- %.8f",fHarmonicQC,2*(co+1),fReferenceFlowHist->GetBinContent(co+1),fReferenceFlowHist->GetBinError(co+1))<<endl;
 }
 cout<<endl;
 Int_t nEvts = 0;
 Double_t dAvM = 0.;
 if(fMultDistributionsHist[0])
 {
  nEvts = (Int_t)fMultDistributionsHist[0]->GetEntries();
  dAvM = fMultDistributionsHist[0]->GetMean();
 } else{Warning(sMethodName.Data(),"fMultDistributionsHist[0] is NULL !!!!");}
 cout<<Form("     nEvts = %d, <M> = %.2f",nEvts,dAvM)<<endl;
 cout<<"*************************************"<<endl;
 cout<<"*************************************"<<endl;
 cout<<endl;

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateReferenceFlow()

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::Covariance(const char *x, const char *y, TProfile2D *profile2D, Bool_t bUnbiasedEstimator)
{
 // Calculate covariance (multiplied by a weight dependent factor).

 // Remark: wCov = Cov(<x>,<y>) * (sum_{i=1}^{N} w_{<x>}_i w_{<y>}_i )/[(sum_{i=1}^{N} w_{<x>}_i) * (sum_{j=1}^{N} w_{<y>}_j)],  
 //         where Cov(<x>,<y>) is biased or unbiased estimator (specified via bUnbiasedEstimator) for the covariance.
 //         An unbiased estimator is given for instance in Eq. (C.12) on page 131 of my thesis. 

 Double_t wCov = 0.; // return value

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::Covariance(const char *x, const char *y, TProfile2D *profile2D, Bool_t bBiasedEstimator)";
 if(!profile2D){Fatal(sMethodName.Data(),"Sorry, 'profile2D' is on holidays.");}

 // Basic protection:
 if(!(TString(x).BeginsWith("Cos") || TString(x).BeginsWith("Sin")))
 {
  cout<<Form("And the fatal x is... '%s'. Congratulations!!",x)<<endl; 
  Fatal(sMethodName.Data(),"!(TString(x).BeginsWith(...");
 }
 if(!(TString(y).BeginsWith("Cos") || TString(y).BeginsWith("Sin")))
 {
  cout<<Form("And the fatal y is... '%s'. Congratulations!!",y)<<endl; 
  Fatal(sMethodName.Data(),"!(TString(y).BeginsWith(...");
 }

 // Determine 'cs' (cosine or sinus) indices for x and y:
 Int_t csx = 0; if(TString(x).BeginsWith("Sin")){csx = 1;}
 Int_t csy = 0; if(TString(y).BeginsWith("Sin")){csy = 1;}

 // Determine 'c' (order of correlator) indices for x and y:
 Int_t cx = -1;   
 for(Int_t t=0;t<=TString(x).Length();t++)
 {
  if(TString(x[t]).EqualTo(",") || TString(x[t]).EqualTo(")")) // TBI this is just ugly
  {
   cx++;
   if(cx>=8){Fatal(sMethodName.Data(),"cx>=8");} // not supporting corr. beyond 8p 
  } // if(TString(x[t]).EqualTo(",") || TString(x[t]).EqualTo(")")) // TBI this is just ugly
 } // for(Int_t t=0;t<=TString(x).Length();t++)
 Int_t cy = -1;   
 for(Int_t t=0;t<=TString(y).Length();t++)
 {
  if(TString(y[t]).EqualTo(",") || TString(y[t]).EqualTo(")")) // TBI this is just ugly
  {
   cy++;
   if(cy>=8){Fatal(sMethodName.Data(),"cy>=8");} // not supporting corr. beyond 8p 
  } // if(TString(y[t]).EqualTo(",") || TString(y[t]).EqualTo(")")) // TBI this is just ugly
 } // for(Int_t t=0;t<=TString(y).Length();t++)

 // Correlations corresponding to x and y:
 // x:
 Double_t dx = 0.; // <<x>>
 Double_t wx = 0.; // \sum w_x
 Int_t nbx = fCorrelationsPro[csx][cx]->GetNbinsX();
 for(Int_t b=1;b<=nbx;b++)
 {
  TString sBinLabel = fCorrelationsPro[csx][cx]->GetXaxis()->GetBinLabel(b);
  if(sBinLabel.EqualTo(x))
  {
   //cout<<Form("b = %d, binLabel = %s",b,sBinLabel.Data())<<endl;
   dx = fCorrelationsPro[csx][cx]->GetBinContent(b);
   wx = fCorrelationsPro[csx][cx]->GetBinEntries(b);
   break;
  } // if(sBinLabel.EqualTo(x))
  if(sBinLabel.EqualTo("")){break;}
 } // for(Int_t b=1;b<=nbx;b++)
 if(TMath::Abs(dx)<1.e-44)
 {
  Warning(sMethodName.Data(),"TMath::Abs(dx)<1.e-44, %s",x);
  fPropagateError = kFALSE;
  return wCov;
 } 
 if(TMath::Abs(wx)<1.e-44)
 {
  Warning(sMethodName.Data(),"TMath::Abs(wx)<1.e-44, %s",x);
  fPropagateError = kFALSE;
  return wCov;
 }

 // y:
 Double_t dy = 0.; // <<y>> 
 Double_t wy = 0.; // \sum w_y
 Int_t nby = fCorrelationsPro[csy][cy]->GetNbinsX();
 for(Int_t b=1;b<=nby;b++)
 {
  TString sBinLabel = fCorrelationsPro[csy][cy]->GetXaxis()->GetBinLabel(b);
  if(sBinLabel.EqualTo(y))
  {
   //cout<<Form("b = %d, binLabel = %s",b,sBinLabel.Data())<<endl;
   dy = fCorrelationsPro[csy][cy]->GetBinContent(b);
   wy = fCorrelationsPro[csy][cy]->GetBinEntries(b);
   break;
  } // if(sBinLabel.EqualTo(y))
  if(sBinLabel.EqualTo("")){break;}
 } // for(Int_t b=1;b<=nby;b++)
 if(TMath::Abs(dy)<1.e-44)
 {
  Warning(sMethodName.Data(),"TMath::Abs(dy)<1.e-44, %s",y);
  fPropagateError = kFALSE;
  return wCov; 
 }
 if(TMath::Abs(wy)<1.e-44)
 {
  Warning(sMethodName.Data(),"TMath::Abs(wy)<1.e-44, %s",y);
  fPropagateError = kFALSE;
  return wCov; 
 } 

 // Product: 
 Double_t dxy = 0.; // <<xy>>
 Double_t wxy = 0.; // \sum w_x*w_y
 // x:
 Int_t nBinsX = profile2D->GetNbinsX();
 Int_t gbx = 0; // generic bin for x 
 for(Int_t b=1;b<=nBinsX;b++)
 {
  TString sBinLabel = profile2D->GetXaxis()->GetBinLabel(b);
  if(sBinLabel.EqualTo(x))
  {
   gbx = b; break;
  } 
 } // for(Int_t b=1;b<=nBins2D;b++)
 if(0==gbx){Fatal(sMethodName.Data(),"0==gbx, %s",x);} 
 // y:
 Int_t nBinsY = profile2D->GetNbinsY();
 Int_t gby = 0; // generic bin for y 
 for(Int_t b=1;b<=nBinsY;b++)
 {
  TString sBinLabel = profile2D->GetYaxis()->GetBinLabel(b);
  if(sBinLabel.EqualTo(y))
  {
   gby = b; break;
  } 
 } // for(Int_t b=1;b<=nBinsY;b++)
 if(0==gby){Fatal(sMethodName.Data(),"0==gby, %s",y);} 

 if(gbx>gby)
 {
  dxy = profile2D->GetBinContent(profile2D->GetBin(gbx,gby));
  wxy = profile2D->GetBinEntries(profile2D->GetBin(gbx,gby));
 }
 else if(gbx<gby)
 {
  dxy = profile2D->GetBinContent(profile2D->GetBin(gby,gbx));
  wxy = profile2D->GetBinEntries(profile2D->GetBin(gby,gbx));
 } else{Fatal(sMethodName.Data(),"gbx==gby, %s, %s",x,y);}
 if(TMath::Abs(dxy)<1.e-44)
 {
  Warning(sMethodName.Data(),"TMath::Abs(dxy)<1.e-44, %s, %s",x,y);
  fPropagateError = kFALSE;
  return wCov; 
 } 
 if(TMath::Abs(wxy)<1.e-44)
 {
  Warning(sMethodName.Data(),"TMath::Abs(wxy)<1.e-44, %s, %s",x,y);
  fPropagateError = kFALSE;
  return wCov; 
 } 

 // Finally:
 if(bUnbiasedEstimator)
 {
  Double_t num = dxy-dx*dy; // numerator of Eq. (C.12) on page 131 of my thesis
  Double_t den = 1.-wxy/(wx*wy); // denominator of Eq. (C.12) on page 131 of my thesis 
  Double_t pre = wxy/(wx*wy); // prefactor
  if(TMath::Abs(den)<1.e-44)
  {
   Warning(sMethodName.Data(),"TMath::Abs(den)<1.e-44, %s, %s",x,y);
   fPropagateError = kFALSE;
   return wCov; 
  }
  wCov = pre*num/den;  
  if(TMath::Abs(wCov)<1.e-44)
  {
   Warning(sMethodName.Data(),"TMath::Abs(wCov)<1.e-44, %s, %s",x,y);
   fPropagateError = kFALSE;
   return wCov; 
  }
 } else
   {
    // TBI check if this is a correct formula for the biased estimator
    Double_t num = dxy-dx*dy; // numerator of Eq. (C.12) on page 131 of my thesis
    Double_t den = 1.; 
    Double_t pre = wxy/(wx*wy); // prefactor
    if(TMath::Abs(den)<1.e-44)
    {
     Warning(sMethodName.Data(),"TMath::Abs(den)<1.e-44, %s, %s",x,y);
     fPropagateError = kFALSE;
     return wCov; 
    }
    wCov = pre*num/den;  
    if(TMath::Abs(wCov)<1.e-44)
    {
     Warning(sMethodName.Data(),"TMath::Abs(wCov)<1.e-44, %s, %s",x,y);
     fPropagateError = kFALSE;
     return wCov; 
    }
   } 

 return wCov;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::Covariance(const char *x, const char *y, TProfile2D *profile2D, Bool_t bUnbiasedEstimator = kFALSE)

//=======================================================================================================================

/*
TComplex AliFlowAnalysisWithMultiparticleCorrelations::Recursion(Int_t n, Int_t* harmonic, Int_t* mult)
{
 // Calculate multi-particle correlators by using recursion originally developed by 
 // Kristjan Gulbrandsen (gulbrand@nbi.dk). 

 TComplex c = Q(harmonic[n-1], mult[n-1]);
 if (n == 1) return c;
 c *= Recursion(n-1, harmonic, mult);
 if (mult[n-1]>1) return c;
 for (Int_t i=0; i<(n-1); i++) {
    harmonic[i] += harmonic[n-1];
    mult[i]++;
    c -= (mult[i]-1.)*Recursion(n-1, harmonic, mult);
    mult[i]--;
    harmonic[i] -= harmonic[n-1];
  }

 return c;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Recursion(Int_t n, Int_t* harmonic, Int_t* mult)
*/

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Recursion(Int_t n, Int_t* harmonic, Int_t mult, Int_t skip) 
{
 // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by 
 // Kristjan Gulbrandsen (gulbrand@nbi.dk). 

  Int_t nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip) return c;

  Int_t multp1 = mult+1;
  Int_t nm2 = n-2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Recursion(Int_t n, Int_t* harmonic, Int_t mult, Int_t skip) 

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::OneDiff(Int_t n1)
{
 // Generic differential one-particle correlation <exp[i(n1*psi1)]>.
 // (psi labels POI, phi labels RPs)  

 TComplex one = p(n1,1);

 return one;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::OneDiff(Int_t n1)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::TwoDiff(Int_t n1, Int_t n2)
{
 // Generic differential two-particle correlation <exp[i(n1*psi1+n2*phi2)]>.
 // (psi labels POI, phi labels RPs)  

 TComplex two = p(n1,1)*Q(n2,1)-q(n1+n2,2);

 return two;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::TwoDiff(Int_t n1, Int_t n2)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::ThreeDiff(Int_t n1, Int_t n2, Int_t n3)
{
 // Generic differential three-particle correlation <exp[i(n1*psi1+n2*phi2+n3*phi3)]>.
 // (psi labels POI, phi labels RPs)  

 TComplex three = p(n1,1)*Q(n2,1)*Q(n3,1)-q(n1+n2,2)*Q(n3,1)-q(n1+n3,2)*Q(n2,1)
                - p(n1,1)*Q(n2+n3,2)+2.*q(n1+n2+n3,3); 

 return three;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::ThreeDiff(Int_t n1, Int_t n2, Int_t n3)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::FourDiff(Int_t n1, Int_t n2, Int_t n3, Int_t n4)
{
 // Generic differential four-particle correlation <exp[i(n1*psi1+n2*phi2+n3*phi3+n4*phi4)]>.
 // (psi labels POI, phi labels RPs)  

 TComplex four = p(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*q(n1+n3,2)*Q(n4,1)
               - p(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*q(n1+n4,2)
               + Q(n2+n3,2)*q(n1+n4,2)-p(n1,1)*Q(n3,1)*Q(n2+n4,2)+q(n1+n3,2)*Q(n2+n4,2)
               + 2.*Q(n3,1)*q(n1+n2+n4,3)-p(n1,1)*Q(n2,1)*Q(n3+n4,2)+q(n1+n2,2)*Q(n3+n4,2)
               + 2.*Q(n2,1)*q(n1+n3+n4,3)+2.*p(n1,1)*Q(n2+n3+n4,3)-6.*q(n1+n2+n3+n4,4);

 return four;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::FourDiff(Int_t n1, Int_t n2, Int_t n3, Int_t n4)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::SetDiffHarmonics(Int_t order, Int_t *harmonics)
{
 // Set harmonics for all differential correlators. 

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::SetDiffHarmonics(Int_t order, Int_t *harmonics)";

 // Basic protection:
 if(order<=0||order>4){Fatal(sMethodName.Data(),"order<=0||order>4");}
 if(!harmonics){Fatal(sMethodName.Data(),"!harmonics");}

 for(Int_t o=0;o<order;o++)
 {
  fDiffHarmonics[order-1][o] = harmonics[o];
 }
 
 fCalculateDiffCorrelations = kTRUE;

} // void AliFlowAnalysisWithMultiparticleCorrelations::SetDiffHarmonics(Int_t order, Int_t *harmonics)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::SetWeightsHist(TH1D* const hist, const char *type, const char *variable)
{
 // Pass histogram holding weights from an external file to the corresponding data member. 
 
 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::SetWeightsHist(TH1D* const hist, const char *type, const char *variable)";
 
 // Basic protection:
 if(!hist){Fatal(sMethodName.Data(),"hist");}
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI"))){Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);}
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))){Fatal(sMethodName.Data(),"!(TString(variable).EqualTo... variable = %s ",variable);}

 Int_t rp = 0; // [RP,POI]
 if(TString(type).EqualTo("POI")){rp=1;} 

 Int_t ppe = 0; // [phi,pt,eta]
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 // Finally:
 hist->SetDirectory(0);
 fWeightsHist[rp][ppe] = (TH1D*)hist->Clone();
 if(!fWeightsHist[rp][ppe]){Fatal(sMethodName.Data(),"fWeightsHist[%d][%d]",rp,ppe);}

 // Cosmetics:
 TString sType[2] = {"RP","POI"};
 TString sVariable[3] = {"phi","pt","eta"};
 fWeightsHist[rp][ppe]->SetName(Form("%s weights (%s)",sVariable[ppe].Data(),sType[rp].Data()));
 //fWeightsHist[rp][ppe]->SetTitle(Form("%s weights (%s)",sVariable[ppe].Data(),sType[rp].Data()));

 // Flag:
 fUseWeights[rp][ppe] = kTRUE; 

} // void AliFlowAnalysisWithMultiparticleCorrelations::SetWeightsHist(TH1D* const hwh, const char *type, const char *variable)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForWeights()
{
 // Initialize all arrays for weights. 

 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
  for(Int_t ppe=0;ppe<3;ppe++) // [phi,pt,eta]
  {
   fUseWeights[rp][ppe] = kFALSE;
   fWeightsHist[rp][ppe] = NULL; 
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForWeights()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::SetnBins(const char *type, const char *variable, Int_t nBins)
{
 // Set number of bins for histograms fKinematicsHist[2][3].

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::SetnBins(const char *type, const char *variable, Int_t nBins)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI")))
 {
  cout<<"Well, it would be better for you to use RP or POI here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta")))
 {
  cout<<"phi, pt or eta, please!"<<endl;
  Fatal(sMethodName.Data(),"!(TString(variable).EqualTo... variable = %s ",variable);
 }

 Int_t rp = 0; // [RP,POI]
 if(TString(type).EqualTo("POI")){rp=1;} 

 Int_t ppe = 0; // [phi,pt,eta]
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 fnBins[rp][ppe] = nBins;

} // void AliFlowAnalysisWithMultiparticleCorrelations::SetnBins(const char *type, const char *variable, Int_t nBins)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::SetMin(const char *type, const char *variable, Double_t min)
{
 // Set min bin range for histograms fKinematicsHist[2][3].

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::SetMin(const char *type, const char *variable, Double_t min)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI")))
 {
  cout<<"Well, it would be better for you to use RP or POI here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta")))
 {
  cout<<"phi, pt or eta, please!"<<endl;
  Fatal(sMethodName.Data(),"!(TString(variable).EqualTo... variable = %s ",variable);
 }

 Int_t rp = 0; // [RP,POI]
 if(TString(type).EqualTo("POI")){rp=1;} 

 Int_t ppe = 0; // [phi,pt,eta]
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 fMin[rp][ppe] = min;

} // void AliFlowAnalysisWithMultiparticleCorrelations::SetMin(const char *type, const char *variable, Double_t min)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::SetMax(const char *type, const char *variable, Double_t max)
{
 // Set max bin range for histograms fKinematicsHist[2][3].

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::SetMax(const char *type, const char *variable, Double_t max)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI")))
 {
  cout<<"Well, it would be better for you to use RP or POI here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta")))
 {
  cout<<"phi, pt or eta, please!"<<endl;
  Fatal(sMethodName.Data(),"!(TString(variable).EqualTo... variable = %s ",variable);
 }

 Int_t rp = 0; // [RP,POI]
 if(TString(type).EqualTo("POI")){rp=1;} 

 Int_t ppe = 0; // [phi,pt,eta]
 if(TString(variable).EqualTo("pt")){ppe=1;} 
 if(TString(variable).EqualTo("eta")){ppe=2;} 

 fMax[rp][ppe] = max;

} // void AliFlowAnalysisWithMultiparticleCorrelations::SetMax(const char *type, const char *variable, Double_t min)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::SetnBinsMult(const char *type, Int_t nBinsMult)
{
 // Set number of bins for histograms fMultDistributionsHist[3].

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::SetnBinsMult(const char *type, Int_t nBinsMult)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI") || TString(type).EqualTo("REF")))
 {
  cout<<"Well, it would be better for you to use RP, POI or REF here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }

 Int_t rpr = 0; // [RP,POI,REF]
 if(TString(type).EqualTo("POI")){rpr=1;} 
 else if(TString(type).EqualTo("REF")){rpr=2;} 

 fnBinsMult[rpr] = nBinsMult;

} // void AliFlowAnalysisWithMultiparticleCorrelations::SetnBinsMult(const char *type, Int_t nBinsMult)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::SetMinMult(const char *type, Double_t minMult)
{
 // Set min bin range for histograms fMultDistributionsHist[3].

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::SetMinMult(const char *type, Double_t minMult)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI") || TString(type).EqualTo("REF")))
 {
  cout<<"Well, it would be better for you to use RP, POI or REF here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }

 Int_t rpr = 0; // [RP,POI,REF]
 if(TString(type).EqualTo("POI")){rpr=1;} 
 else if(TString(type).EqualTo("REF")){rpr=2;} 

 fMinMult[rpr] = minMult;

} // void AliFlowAnalysisWithMultiparticleCorrelations::SetMinMult(const char *type, Double_t minMult)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::SetMaxMult(const char *type, Double_t maxMult)
{
 // Set max bin range for histograms fMultDistributionsHist[3].

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::SetMaxMult(const char *type, Double_t maxMult)";
 
 // Basic protection:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI") || TString(type).EqualTo("REF")))
 {
  cout<<"Well, it would be better for you to use RP, POI or REF here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(type).EqualTo... type = %s ",type);
 }

 Int_t rpr = 0; // [RP,POI,REF]
 if(TString(type).EqualTo("POI")){rpr=1;} 
 else if(TString(type).EqualTo("REF")){rpr=2;} 

 fMaxMult[rpr] = maxMult;

} // void AliFlowAnalysisWithMultiparticleCorrelations::SetMaxMult(const char *type, Double_t minMult)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::SetIntervalsToSkip(const char *ppe, Int_t nBoundaries, Double_t *boundaries)
{
 // Set all pt, phi and eta intervals to be skipped. 

 // Example usage in the steering macro (before Init()):
 //  Double_t skip[4] = {-0.1,0.2,0.8,0.9};
 //  mpc->SetIntervalsToSkip("Eta",4,skip);
 
 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::SetIntervalsToSkip(const char *ppe, Int_t n, Double_t *boundaries)";
 
 // Basic protection:
 if(!(TString(ppe).EqualTo("Phi") || TString(ppe).EqualTo("Pt") || TString(ppe).EqualTo("Eta")))
 {
  cout<<"Well, could you perhaps try to use only Phi, Pt or Eta here..."<<endl;
  Fatal(sMethodName.Data(),"!(TString(ppe).EqualTo... type = %s ",ppe);
 }
  
 if(nBoundaries>10)
 {
  cout<<"Maximum number of boundaries is hardwired to be 10 at the moment, sorry..."<<endl;
  Fatal(sMethodName.Data(),"nBoundaries = %d ",nBoundaries);
 }

 fSkipSomeIntervals = kTRUE;

 Int_t index = -44;
 if(TString(ppe).EqualTo("Phi"))
 {
  index = 0;
 } 
 else if(TString(ppe).EqualTo("Pt"))
 {
  index = 1;
 }
 else
 {
  index = 2;
 }

 for(Int_t b=0;b<nBoundaries;b++) // boundaries
 {
  fSkip[index][b] = boundaries[b];
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::SetIntervalsToSkip(const char *ppe, Int_t n, Double_t *boundaries)

//=======================================================================================================================


void AliFlowAnalysisWithMultiparticleCorrelations::DumpThePoints(AliFlowEventSimple *anEvent)
{
 // Dump the points into the external file. 
 
 // Dumping format: 
 // Event <eventNo> Multiplicity <multRP> 
 // phi pt eta

 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::DumpThePoints(AliFlowEventSimple *anEvent)";

 // Basic protection:
 if(!anEvent){Fatal(sMethodName.Data(),"if(!anEvent)");} 
 if(!fMultDistributionsHist[0]){Fatal(sMethodName.Data(),"if(!fMultDistributionsHist[0])");} 
 if(fMaxNoEventsPerFile<=0){Fatal(sMethodName.Data(),"if(fMaxNoEventsPerFile<=0)");} 

 // Determine event number and multiplicity:
 Int_t eventNo = (Int_t) fMultDistributionsHist[0]->GetEntries(); // TBI this is a little bit shaky...
 Int_t multRP = (Int_t) anEvent->GetNumberOfRPs(); // TBI shall I promote this variable into data member? 
 if(fSkipSomeIntervals){ multRP = multRP - fNumberOfSkippedRPParticles; } // TBI tmp gym

 // Determine external file name:
 Int_t fileCounter = (Int_t)((eventNo-1)/fMaxNoEventsPerFile);
 TString filename = Form("%s_%d.dat",fAnalysisTag.Data(),fileCounter);

 // Open external file and dump:
 ofstream myfile;
 myfile.open(filename.Data(),ios::app); 
 myfile << Form("Event %d Multiplicity %d\n",eventNo,multRP);   
 Int_t nTracks = (Int_t) anEvent->NumberOfTracks();
 Double_t dPhi = 0., dPt = 0., dEta = 0.;
 for(Int_t t=0;t<nTracks;t++) // loop over all tracks
 {
  AliFlowTrackSimple *pTrack = anEvent->GetTrack(t);
  if(!pTrack){printf("\n Error: pTrack is NULL in MPC::DumpThePoints(AliFlowEventSimple *anEvent) !!!!"); continue;}
  if(pTrack->InRPSelection()) 
  {
   dPhi = pTrack->Phi(); 
   dPt = pTrack->Pt();
   dEta = pTrack->Eta();
   myfile<<Form("%f %f %f\n",dPhi,dPt,dEta);
   //cout<<Form("%f %f %f",dPhi,dPt,dEta)<<endl;
  }
 } // for(Int_t t=0;t<nTracks;t++) // loop over all tracks
 myfile<<"\n";
 myfile.close();

} // void AliFlowAnalysisWithMultiparticleCorrelations::DumpThePoints(AliFlowEventSimple *anEvent)

//=======================================================================================================================

TH1D *AliFlowAnalysisWithMultiparticleCorrelations::GetHistogramWithWeights(const char *filePath, const char *listName, const char *type, const char *variable, const char *production)
{
 // Access from external ROOT file the desired histogram with weights. 

 // a) Return value; 
 // b) Method name; 
 // c) Basic protection for arguments; 
 // d) Check if the external ROOT file exists at specified path; 
 // e) Access the external ROOT file and fetch the desired histogram with weights;
 // f) Close the external ROOT file. 

 // a) Return value:
 TH1D *hist = NULL; 

 // b) Method name: 
 TString sMethodName = "Double_t AliFlowAnalysisWithMultiparticleCorrelations::GetHistogramWithWeights(const char *filePath, const char *listName, const char *type, const char *variable, const char *production)"; 

 // c) Basic protection for arguments:
 if(!(TString(type).EqualTo("RP") || TString(type).EqualTo("POI"))){Fatal(sMethodName.Data(),"!(TString(type).EqualTo...");}
 if(!(TString(variable).EqualTo("phi") || TString(variable).EqualTo("pt") || TString(variable).EqualTo("eta"))){Fatal(sMethodName.Data(),"!(TString(variable).EqualTo...");}
 if(!(TString(production).EqualTo("data") || TString(production).BeginsWith("LHC"))){Fatal(sMethodName.Data(),"!(TString(production).EqualTo...");}

 // d) Check if the external ROOT file exists at specified path:
 if(gSystem->AccessPathName(filePath,kFileExists))
 {
  Fatal(sMethodName.Data(),"if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s",filePath);
 }

 // e) Access the external ROOT file and fetch the desired histogram with weights:
 TFile *weightsFile = TFile::Open(filePath,"READ");
 TList *weightsFileLOK = weightsFile->GetListOfKeys(); 
 if(!weightsFileLOK || weightsFileLOK->GetEntries() != 1) // TBI get rid of the 2nd condition at some point...
 {
  //printf("\n => if(!weightsFileLOK || weightsFileLOK->GetEntries() != 1)\n\n"); 
  Fatal(sMethodName.Data(),"if(!weightsFileLOK || weightsFileLOK->GetEntries() != 1)");
 } 
 // Access TDirectoryFile "weightsMPCanalysis":
 TDirectoryFile *directoryFile = dynamic_cast<TDirectoryFile*>(weightsFile->Get("weightsMPCanalysis"));
 if(!directoryFile)
 {
  //printf("\n => if(!directoryFile)\n\n");   
  Fatal(sMethodName.Data(),"if(!directoryFile)");
 } 
 // Access the specified list:
 TList *list = dynamic_cast<TList*>(directoryFile->Get(listName));
 if(!list)
 {
  cout<<Form("listName = %s",listName)<<endl;
  Warning(sMethodName.Data(),"if(!list)"); 
  return NULL;
 }
 // Finally, access the desired histogram:
 hist = dynamic_cast<TH1D*>(list->FindObject(Form("%s,%s,%s",type,variable,production)));
 if(!hist)
 {
  //printf("\n => if(!hist)\n\n");   
  Warning(sMethodName.Data(),"if(!hist)");
  return NULL;
 } else { hist->SetDirectory(0); }

 // f) Close the external ROOT file: 
 weightsFile->Close(); delete weightsFile;

 return hist;

} // TH1D *AliFlowAnalysisWithMultiparticleCorrelations::GetHistogramWithWeights(const char *filePath, const char *listName, const char *type, const char *variable, const char *production)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateSymmetryPlanes(AliFlowEventSimple *anEvent)
{
 // Calculate symmetry plane correlations from Q-vector components.

 // a) Insanity checks;
 // b) Calculate and store symmetry plane correlations.

 // a) Insanity checks:
 TString sMethodName = "void AliFlowAnalysisWithMultiparticleCorrelations::CalculateSymmetryPlanes(AliFlowEventSimple *anEvent)";
 if(!anEvent){Fatal(sMethodName.Data(),"Sorry, anEvent is NULL.");}
 Double_t dMultRP = fSelectRandomlyRPs ? fnSelectedRandomlyRPs : anEvent->GetNumberOfRPs(); // TBI shall I promote this variable into data member?
 if(fSkipSomeIntervals){ dMultRP = dMultRP - fNumberOfSkippedRPParticles; } // TBI tmp gym
 if(dMultRP<0.){Fatal(sMethodName.Data(),"Sorry, dMultRP is negative.");}

 // b) Calculate and store symmetry plane correlations.
 for(Int_t o=0;o<4;o++) // o = optimizer TBI this clearly will be generalized further, also 4 is h.w.
 {
  fSymmetryPlanesPro[0][0]->Fill(0.5+o,CorrelationPsi2nPsi1n(1,o));
  fSymmetryPlanesPro[0][1]->Fill(0.5+o,CorrelationPsi2nPsi1n(2,o));
 } // for(Int_t o=0;o<4;o++) // o = optimizer TBI this clearly will be generalized further, also 4 is h.w.

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateSymmetryPlanes(AliFlowEventSimple *anEvent)

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::CorrelationPsi2nPsi1n(Int_t n, Int_t k)
{
 // TBI Comment the weather here eventually...

 // TBI the optimizer can have the non-trivial harmonic structure as well...
 //     now it's hardwired to n=2

 // a) Return value;
 // b) Method name;
 // c) Insanity checks;

 // a) Return value:
 Double_t ratio = -44.;

 // b) Method name:
 TString sMethodName = "Double_t AliFlowAnalysisWithMultiparticleCorrelations::CorrelationPsi2nPsi1n(Int_t n, Int_t k)";

 // c) Insanity checks:
 // TBI if n is this and that...
 // TBI if k is this and that...

 Int_t order = 6 + 2*k;
 TArrayI harmonics0 = TArrayI(order); harmonics0.Reset(0);
 TArrayI harmonics1 = TArrayI(order); harmonics1.Reset(0);
 TArrayI harmonics2 = TArrayI(order); harmonics2.Reset(0);

 // Harmonics for numerator:
 harmonics1.AddAt(2*n,0);
 harmonics1.AddAt(2*n,1);
 harmonics1.AddAt(-n,2);
 harmonics1.AddAt(-n,3);
 harmonics1.AddAt(-n,4);
 harmonics1.AddAt(-n,5);
 // Harmonics for denominator:
 harmonics2.AddAt(2*n,0);
 harmonics2.AddAt(-2*n,1);
 harmonics2.AddAt(n,2);
 harmonics2.AddAt(-n,3);
 harmonics2.AddAt(n,4);
 harmonics2.AddAt(-n,5);

 switch(k)
 {
  case 0:
   // TBI ???
  break;
  case 1:
   // Additional harmonics for numerator from the optimizer:
   harmonics1.AddAt(2,6);
   harmonics1.AddAt(-2,7);
   // Additional harmonics for denominator from the optimizer:
   harmonics2.AddAt(2,6);
   harmonics2.AddAt(-2,7);
  break;
  case 2:
   // Additional harmonics for numerator from the optimizer:
   harmonics1.AddAt(2,6);
   harmonics1.AddAt(-2,7);
   harmonics1.AddAt(2,8);
   harmonics1.AddAt(-2,9);
   // Additional harmonics for denominator from the optimizer:
   harmonics2.AddAt(2,6);
   harmonics2.AddAt(-2,7);
   harmonics2.AddAt(2,8);
   harmonics2.AddAt(-2,9);
  break;
  case 3:
   // Additional harmonics for numerator from the optimizer:
   harmonics1.AddAt(2,6);
   harmonics1.AddAt(-2,7);
   harmonics1.AddAt(2,8);
   harmonics1.AddAt(-2,9);
   harmonics1.AddAt(2,10);
   harmonics1.AddAt(-2,11);
   // Additional harmonics for denominator from the optimizer:
   harmonics2.AddAt(2,6);
   harmonics2.AddAt(-2,7);
   harmonics2.AddAt(2,8);
   harmonics2.AddAt(-2,9);
   harmonics2.AddAt(2,10);
   harmonics2.AddAt(-2,11);
  break;
  default:
   cout<<Form("And the fatal 'k' value is... %d. Congratulations!!",k)<<endl;
   Fatal(sMethodName.Data(),"switch(k)"); // TBI
 } // switch(k)

 // Calculate weight and correlators:
 Double_t dWeight = Recursion(order,harmonics0.GetArray()).Re(); // weight is 'number of combinations' by default
 TComplex cNum1 = Recursion(order,harmonics1.GetArray())/dWeight;
 TComplex cNum2 = Recursion(order,harmonics2.GetArray())/dWeight;
 ratio = cNum1.Re()/cNum2.Re();

 return ratio;

} // Double_t CorrelationPsi2nPsi1n(Int_t n, Int_t k)






