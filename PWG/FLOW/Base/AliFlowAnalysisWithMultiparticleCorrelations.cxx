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
 fFillControlHistograms(kFALSE),
 fFillKinematicsHist(kFALSE),
 fFillMultDistributionsHist(kFALSE),   
 fFillMultCorrelationsHist(kFALSE),
 // 2.) Q-vector:
 fQvectorList(NULL),       
 fQvectorFlagsPro(NULL),
 fCalculateQvector(kFALSE),
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
 // 4.) Cumulants:
 fCumulantsList(NULL),
 fCumulantsFlagsPro(NULL),
 f2pCumulantsPro(NULL),
 fCalculateCumulants(kFALSE),
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
 fStandardCandlesHist(NULL),
 fProductsPro2D(NULL)
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
 this->BookEverythingForQvector();
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
 // b) Cross-check all pointers used in this method;
 // c) Fill control histograms;
 // d) Fill Q-vector components;
 // e) Calculate multi-particle correlations from Q-vector components; 
 // f) Calculate products of multi-particle correlations (needed for error propagation of SCs);
 // g) Calculate e-b-e cumulants; 
 // h) Reset Q-vector components;
 // i) Cross-check results with nested loops.

 // a) Cross-check internal flags:
 if(fUseInternalFlags){if(!this->CrossCheckInternalFlags(anEvent)){return;}}

 // b) Cross-check all pointers used in this method:
 this->CrossCheckPointersUsedInMake(); 

 // c) Fill control histograms:
 if(fFillControlHistograms){this->FillControlHistograms(anEvent);}
 
 // d) Fill Q-vector components:
 if(fCalculateQvector){this->FillQvector(anEvent);}

 // e) Calculate multi-particle correlations from Q-vector components:
 if(fCalculateCorrelations){this->CalculateCorrelations(anEvent);}

 // f) Calculate products of multi-particle correlations (needed for error propagation of SCs):
 if(fCalculateStandardCandles){this->CalculateProductsOfCorrelations(anEvent);}

 // g) Calculate e-b-e cumulants: 
 if(fCalculateCumulants){this->CalculateCumulants();}

 // h) Reset Q-vector components:
 if(fCalculateQvector){this->ResetQvector();}

 // i) Cross-check results with nested loops:
 if(fCrossCheckWithNestedLoops){this->CrossCheckWithNestedLoops(anEvent);}

} // end of AliFlowAnalysisWithMultiparticleCorrelations::Make(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::Finish()
{
 // Closing the curtains. 

 // a) Cross-check pointers used in this method;
 // b) Calculate 'standard candles'.
 // ...
 
 // a) Cross-check pointers used in this method:
 this->CrossCheckPointersUsedInFinish();

 // b) Calculate 'standard candles':
 if(fCalculateStandardCandles){this->CalculateStandardCandles();}

 // ...
 printf("\n  ... Closing the curtains ... \n\n");

} // end of AliFlowAnalysisWithMultiparticleCorrelations::Finish()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInFinish()
{
 // Cross-check all pointers used in method Finish().

 // a) Correlations;
 // b) 'Standard candles'.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInFinish()"; 

 // a) Correlations:
 if(fCalculateCorrelations)
 {
  for(Int_t cs=0;cs<2;cs++) 
  {
   for(Int_t c=0;c<fMaxCorrelator;c++) 
   {
    if(!fCorrelationsPro[cs][c]){Fatal(sMethodName.Data(),"fCorrelationsPro[%d][%d] is NULL, for one reason or another...",cs,c);}
   }
  }
 } // if(fCalculateCorrelations)

 // b) 'Standard candles':
 if(fCalculateStandardCandles)
 {
  if(!fStandardCandlesHist){Fatal(sMethodName.Data(),"fStandardCandlesHist is NULL, for one reason or another...");}
  if(!fProductsPro2D){Fatal(sMethodName.Data(),"fProductsPro2D is NULL, for one reason or another...");}
 } // if(fCalculateStandardCandles)

} // void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInFinish()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInMake()
{
 // Cross-check all pointers used in method Make().

 // a) Correlations;
 // b) 'Standard candles'.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInMake()"; 

 // a) Correlations:
 if(fCalculateCorrelations)
 {
  for(Int_t cs=0;cs<2;cs++) 
  {
   for(Int_t c=0;c<fMaxCorrelator;c++) 
   {
    if(!fCorrelationsPro[cs][c]){Fatal(sMethodName.Data(),"fCorrelationsPro[%d][%d] is NULL, for one reason or another...",cs,c);}
   }
  }
 } // if(fCalculateCorrelations)

 // b) 'Standard candles':
 if(fCalculateStandardCandles)
 {
  if(!fProductsPro2D){Fatal(sMethodName.Data(),"fProductsPro2D is NULL, for one reason or another...");}
 } // if(fCalculateStandardCandles)

} // void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckPointersUsedInMake()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateStandardCandles()
{
 // Calculate 'standard candles'.

 // 'Standard candle' (SC) is defined in terms of average (all-event!) correlations as follows: 
 //   SC(-n1,-n2,n2,n1) = <<Cos(-n1,-n2,n2,n1)>> - <<Cos(-n1,n1)>>*<<Cos(-n2,n2)>>, n1 > n2.

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CalculateStandardCandles()"; 

 Int_t nBins2p = fCorrelationsPro[0][1]->GetXaxis()->GetNbins();
 Int_t nBins4p = fCorrelationsPro[0][3]->GetXaxis()->GetNbins();
 Int_t nBins2D = fProductsPro2D->GetXaxis()->GetNbins();
 Int_t nBinsSC = fStandardCandlesHist->GetXaxis()->GetNbins();
 Double_t dSCn1n2n2n1 = 0.; // SC(-n1,-n2,n2,n1)
 Double_t dSCn1n2n2n1Err = 0.; // stat. error of SC(-n1,-n2,n2,n1)
 Double_t dCovCosn2n2Cosn1n1 = 0.; // Cov[<Cos(-n2,n2)>*<Cos(-n1,n1)>]
 Double_t dCovCosn1n2n2n1Cosn1n1 = 0.; // Cov[<Cos(-n1,-n2,n2,n1)>]*<Cos(-n1,n1)>
 Double_t dCovCosn1n2n2n1Cosn2n2 = 0.; // Cov[<Cos(-n1,-n2,n2,n1)>*<Cos(-n2,n2)>]
 Int_t binNo = 1; // for fStandardCandlesHist

 for(Int_t n1=-fMaxHarmonic;n1<=-2;n1++)
 {
  Double_t dCosn1n1 = 0.; // <<Cos(-n1,n1)>>
  Double_t dCosn1n1Err = 0.; // stat. error of <<Cos(-n1,n1)>> 
  Double_t dCosn1n1SumW = 0.; // total sum of weights for <Cos(-n1,n1)> 
  TString sn1n1 = Form("(%d,%d)",n1,-1*n1); // pattern (-n1,n1) 
  for(Int_t bin2p=1;bin2p<=nBins2p;bin2p++)
  {
   TString binLabel2p = fCorrelationsPro[0][1]->GetXaxis()->GetBinLabel(bin2p);
   if(!binLabel2p.CompareTo("")){break;} // TBI this is a little bit shaky...
   if(binLabel2p.EndsWith(sn1n1.Data())) // TBI this is a little bit shaky as well...
   {
    dCosn1n1 = fCorrelationsPro[0][1]->GetBinContent(bin2p);
    dCosn1n1Err = fCorrelationsPro[0][1]->GetBinError(bin2p);
    dCosn1n1SumW = fCorrelationsPro[0][1]->GetBinEntries(bin2p);
   } // if(binLabel2p.EndsWith(sn1n1.Data()))
  } // for(Int_t bin2p=1;bin2p<=nBins2p;bin2p++)

  for(Int_t n2=n1+1;n2<=-1;n2++)
  {
   Double_t dCosn2n2 = 0.; // <<Cos(-n2,n2)>>
   Double_t dCosn2n2Err = 0.; // stat. error of <<Cos(-n2,n2)>> 
   Double_t dCosn2n2SumW = 0.; // total sum of weights for <Cos(-n2,n2)> 
   TString sn2n2 = Form("(%d,%d)",n2,-1*n2); // pattern (-n2,n2) 
   for(Int_t bin2p=1;bin2p<=nBins2p;bin2p++)
   {
    TString binLabel2p = fCorrelationsPro[0][1]->GetXaxis()->GetBinLabel(bin2p);
    if(!binLabel2p.CompareTo("")){break;} // TBI this is a little bit shaky...
    if(binLabel2p.EndsWith(sn2n2.Data())) // TBI this is a little bit shaky as well...
    {
     dCosn2n2 = fCorrelationsPro[0][1]->GetBinContent(bin2p);
     dCosn2n2Err = fCorrelationsPro[0][1]->GetBinError(bin2p);
     dCosn2n2SumW = fCorrelationsPro[0][1]->GetBinEntries(bin2p);
    } // if(binLabel2p.EndsWith(sn2n2.Data()))
   } // for(Int_t bin2p=1;bin2p<=nBins2p;bin2p++)
   Double_t dCosn1n2n2n1 = 0.; // <<Cos(-n1,-n2,n2,n1)>>
   Double_t dCosn1n2n2n1Err = 0.; // stat. error of <<Cos(-n1,-n2,n2,n1)>> 
   Double_t dCosn1n2n2n1SumW = 0.; // total sum of weights for <Cos(-n1,-n2,n2,n1)>
   TString sn1n2n2n1 = Form("(%d,%d,%d,%d)",n1,n2,-1*n2,-1*n1); // pattern (-n1,-n2,n2,n1) 
   for(Int_t bin4p=1;bin4p<=nBins4p;bin4p++)
   {
    TString binLabel4p = fCorrelationsPro[0][3]->GetXaxis()->GetBinLabel(bin4p);
    if(!binLabel4p.CompareTo("")){break;} // TBI this is a little bit shaky...
    if(binLabel4p.EndsWith(sn1n2n2n1.Data())) // TBI this is a little bit shaky as well...
    {
     dCosn1n2n2n1 = fCorrelationsPro[0][3]->GetBinContent(bin4p);
     dCosn1n2n2n1Err = fCorrelationsPro[0][3]->GetBinError(bin4p);
     dCosn1n2n2n1SumW = fCorrelationsPro[0][3]->GetBinEntries(bin4p);
    } // if(binLabel4p.EndsWith(sn1n2n2n1.Data()))
   } // for(Int_t bin4p=1;bin4p<=nBins4p;bin4p++)
 
   // Get the needed products, and calculate covariances for error propagation:
   for(Int_t bx=2;bx<=nBins2D;bx++) 
   { 
    if(!(TString(fProductsPro2D->GetXaxis()->GetBinLabel(bx)).Contains(sn1n1.Data()) 
         || TString(fProductsPro2D->GetXaxis()->GetBinLabel(bx)).Contains(sn2n2.Data()) 
         || TString(fProductsPro2D->GetXaxis()->GetBinLabel(bx)).Contains(sn1n2n2n1.Data()))){continue;}
    for(Int_t by=1;by<bx;by++) 
    {
     if(!(TString(fProductsPro2D->GetYaxis()->GetBinLabel(by)).Contains(sn1n1.Data()) 
          || TString(fProductsPro2D->GetYaxis()->GetBinLabel(by)).Contains(sn2n2.Data()) 
          || TString(fProductsPro2D->GetYaxis()->GetBinLabel(by)).Contains(sn1n2n2n1.Data()))){continue;}
     //  Cov[<Cos(-n2,n2)>*<Cos(-n1,n1)>]:
     if(TString(fProductsPro2D->GetXaxis()->GetBinLabel(bx)).Contains(sn2n2.Data()) 
        && TString(fProductsPro2D->GetYaxis()->GetBinLabel(by)).Contains(sn1n1.Data()))    
       { 
        if(dCosn2n2SumW > 0. && dCosn1n1SumW > 0.)
        {
         dCovCosn2n2Cosn1n1 = (fProductsPro2D->GetBinEntries(fProductsPro2D->GetBin(bx,by))/(dCosn2n2SumW*dCosn1n1SumW))
                            * (fProductsPro2D->GetBinContent(fProductsPro2D->GetBin(bx,by))-dCosn2n2*dCosn1n1);
        } else{Fatal(sMethodName.Data(),"dCosn2n2SumW > 0. && dCosn1n1SumW > 0.");}
       }
     // Cov[<Cos(-n1,-n2,n2,n1)>*<Cos(-n2,n2)>];
     else if(TString(fProductsPro2D->GetXaxis()->GetBinLabel(bx)).Contains(sn1n2n2n1.Data()) 
             && TString(fProductsPro2D->GetYaxis()->GetBinLabel(by)).Contains(sn2n2.Data()))    
       {
        if(dCosn1n2n2n1SumW > 0. && dCosn2n2SumW > 0.)
        {
         dCovCosn1n2n2n1Cosn2n2 = (fProductsPro2D->GetBinEntries(fProductsPro2D->GetBin(bx,by))/(dCosn2n2SumW*dCosn1n2n2n1SumW))
                                * (fProductsPro2D->GetBinContent(fProductsPro2D->GetBin(bx,by))-dCosn1n2n2n1*dCosn2n2);
        } else{Fatal(sMethodName.Data(),"dCosn1n2n2n1SumW > 0. && dCosn2n2SumW > 0.");}
       }
    // Cov[<Cos(-n1,-n2,n2,n1)>*<Cos(-n1,n1)>]:
    else if(TString(fProductsPro2D->GetXaxis()->GetBinLabel(bx)).Contains(sn1n2n2n1.Data()) 
            && TString(fProductsPro2D->GetYaxis()->GetBinLabel(by)).Contains(sn1n1.Data()))    
       {       
        if(dCosn1n2n2n1SumW > 0. && dCosn1n1SumW > 0.)
        {
         dCovCosn1n2n2n1Cosn1n1 = (fProductsPro2D->GetBinEntries(fProductsPro2D->GetBin(bx,by))/(dCosn1n1SumW*dCosn1n2n2n1SumW))
                                * (fProductsPro2D->GetBinContent(fProductsPro2D->GetBin(bx,by))-dCosn1n2n2n1*dCosn1n1);
        } else{Fatal(sMethodName.Data(),"if(dCosn1n2n2n1Cosn1n1SumW > 0. && dCosn1n1SumW > 0.)");}
       }
    } // for(Int_t by=1;by<bx;by++)
   } // for(Int_t bx=2;bx<=nBins2D;bx++)
   // Finally, calculate and store SCs:
   dSCn1n2n2n1 = dCosn1n2n2n1 - dCosn1n1*dCosn2n2;
   dSCn1n2n2n1Err = pow(dCosn1n1,2.)*pow(dCosn2n2Err,2.) + pow(dCosn2n2,2.)*pow(dCosn1n1Err,2.) + pow(dCosn1n2n2n1Err,2.)
                  + 2.*dCosn1n1*dCosn2n2*dCovCosn2n2Cosn1n1 
                  - 2.*dCosn1n1*dCovCosn1n2n2n1Cosn2n2 
                  - 2.*dCosn2n2*dCovCosn1n2n2n1Cosn1n1; // note that this is still error^2
   fStandardCandlesHist->SetBinContent(binNo,dSCn1n2n2n1); 
   if(dSCn1n2n2n1Err>0.){fStandardCandlesHist->SetBinError(binNo,pow(dSCn1n2n2n1Err,0.5));}
   else{Warning(sMethodName.Data(),"if(dSCn1n2n2n1Err>0.)");}
   binNo++;
  } // for(Int_t n2=n1+1;n2<=-1;n2++)
 } // for(Int_t n1=-fMaxHarmonic;n1<=-2;n1++)

 if(nBinsSC != binNo-1){Fatal(sMethodName.Data(),"nBinsSC != binNo-1");} // just in case...

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

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateCorrelations(AliFlowEventSimple *anEvent)
{
 // Calculate multi-particle correlations from Q-vector components.

 // Fill ... TBI this can be implemented better, most likely...

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CalculateCorrelations()"; 
 if(!anEvent){Fatal(sMethodName.Data(),"anEvent is apparently doing crazy things...");} // TBI

 Double_t dMultRP = anEvent->GetNumberOfRPs(); // TBI shall I promote this variable into data member? 

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
    fCorrelationsPro[0][0]->Fill(binNo[0]-.5,oneN.Re()/oneD,oneW);
    fCorrelationsPro[1][0]->Fill(binNo[0]++-.5,oneN.Im()/oneD,oneW);
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
    TComplex twoN = Two(n1,n2); // numerator
    Double_t twoD = Two(0,0).Re(); // denominator
    Double_t twoW = twoD; // weight TBI add other possibilities here for the weight
    if(twoD>0. && dMultRP>=2) 
    {
     fCorrelationsPro[0][1]->Fill(binNo[1]-.5,twoN.Re()/twoD,twoW);;
     fCorrelationsPro[1][1]->Fill(binNo[1]++-.5,twoN.Im()/twoD,twoW);;
    } else {Warning(sMethodName.Data(),"twoD>0. &&d MultRP>=2");} 
   } 
   if(2==fDontGoBeyond){continue;}
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
      fCorrelationsPro[0][2]->Fill(binNo[2]-.5,threeN.Re()/threeD,threeW);
      fCorrelationsPro[1][2]->Fill(binNo[2]++-.5,threeN.Im()/threeD,threeW);
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
       fCorrelationsPro[0][3]->Fill(binNo[3]-.5,fourN.Re()/fourD,fourW);
       fCorrelationsPro[1][3]->Fill(binNo[3]++-.5,fourN.Im()/fourD,fourW);
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
        fCorrelationsPro[0][4]->Fill(binNo[4]-.5,fiveN.Re()/fiveD,fiveW);
        fCorrelationsPro[1][4]->Fill(binNo[4]++-.5,fiveN.Im()/fiveD,fiveW);
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
         fCorrelationsPro[0][5]->Fill(binNo[5]-.5,sixN.Re()/sixD,sixW);
         fCorrelationsPro[1][5]->Fill(binNo[5]++-.5,sixN.Im()/sixD,sixW);
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
          fCorrelationsPro[0][6]->Fill(binNo[6]-.5,sevenN.Re()/sevenD,sevenW);
          fCorrelationsPro[1][6]->Fill(binNo[6]++-.5,sevenN.Im()/sevenD,sevenW);
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
           fCorrelationsPro[0][7]->Fill(binNo[7]-.5,eightN.Re()/eightD,eightW);
           fCorrelationsPro[1][7]->Fill(binNo[7]++-.5,eightN.Im()/eightD,eightW);
          }
         } 
        } // for(Int_t n8=n7;n8<=fMaxHarmonic;n8++)
       } // for(Int_t n7=n6;n7<=fMaxHarmonic;n7++) 
      } // for(Int_t n6=n5;n6<=fMaxHarmonic;n6++) 
     } // for(Int_t n5=n4;n5<=fMaxHarmonic;n5++) 
    } // for(Int_t n4=n3;n4<=fMaxHarmonic;n4++)   
   } // for(Int_t n3=n2;n3<=fMaxHarmonic;n3++) 
  } // for(Int_t n2=n1;n2<=fMaxHarmonic;n2++)
 } // for(Int_t n1=-fMaxHarmonic;n1<=fMaxHarmonic;n1++) 

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateCorrelations(AliFlowEventSimple *anEvent)

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::CastStringToCorrelation(const char *string, Bool_t numerator)
{
 // Cast string of the generic form Cos/Sin(-n_1,-n_2,...,n_{k-1},n_k) in the corresponding correlation value.
 // If you issue a call to this method with setting numerator = kFALSE, then you are getting back for free
 // the corresponding denumerator (a.k.a. weight 'number of combinations').

 // TBI:
 // a) add protection against cases a la:
 //     string = Cos(-3,-4,5,6,5,6-3)
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

void AliFlowAnalysisWithMultiparticleCorrelations::CalculateProductsOfCorrelations(AliFlowEventSimple *anEvent)
{
 // Calculate products of multi-particle correlations (needed for error propagation of SCs).

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::CalculateProductsOfCorrelations()"; 
 if(!anEvent){Fatal(sMethodName.Data(),"Sorry, 'anEvent' is on holidays.");} 

 Int_t nBins = fProductsPro2D->GetXaxis()->GetNbins();
 for(Int_t bx=2;bx<=nBins;bx++)
 {
  for(Int_t by=1;by<bx;by++)
  {
   const char *binLabelX = fProductsPro2D->GetXaxis()->GetBinLabel(bx);
   const char *binLabelY = fProductsPro2D->GetYaxis()->GetBinLabel(by);
   Double_t numX = this->CastStringToCorrelation(binLabelX,kTRUE); // numerator
   Double_t denX = this->CastStringToCorrelation(binLabelX,kFALSE); // denominator
   Double_t wX = denX; // weight TBI add support for other options
   Double_t numY = this->CastStringToCorrelation(binLabelY,kTRUE); // numerator
   Double_t denY = this->CastStringToCorrelation(binLabelY,kFALSE); // denominator
   Double_t wY = denY; // weight TBI add support for other options
   if(TMath::Abs(denX) > 0. && TMath::Abs(denY) > 0.)
   {
    fProductsPro2D->Fill(bx-0.5,by-0.5,(numX/denX)*(numY/denY),wX*wY);
   } 
   else{Fatal(sMethodName.Data(),"if(TMath::Abs(denX) > 0. && TMath::Abs(denY) > 0.)");}
  } // for(Int_t by=1;by<bx;by++)
 } // for(Int_t bx=2;bx<=nbins;bx++)

} // void AliFlowAnalysisWithMultiparticleCorrelations::CalculateProductsOfCorrelations(AliFlowEventSimple *anEvent)

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
   dPhi1 = aftsTrack->Phi(); 
   if(fUsePhiWeights){wPhi1 = PhiWeight(dPhi1);}
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
 if(dMultRP>=3)
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
 if(dMultRP>=4)
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
 if(dMultRP>=5)
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
 if(dMultRP>=6)
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
 if(dMultRP>=7)
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
 if(dMultRP>=8)
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
   for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)
   {
    for(Int_t p=0;p<fMaxCorrelator+1;p++)
    {
     if(fUsePhiWeights||fUsePtWeights||fUseEtaWeights){wToPowerP = pow(wPhi*wPt*wEta,p);} // TBI should I do something with the normalization of the product wPhi*wPt*wEta
     fQvector[h][p] += TComplex(wToPowerP*TMath::Cos(h*dPhi),wToPowerP*TMath::Sin(h*dPhi));
    } // for(Int_t p=0;p<fMaxCorrelator+1;p++)
   } // for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++)
  } // if(pTrack && pTrack->InRPSelection()) // fill Q-vector components only with reference particles
 } // for(Int_t t=0;t<nTracks;t++) // loop over all tracks

} // void AliFlowAnalysisWithMultiparticleCorrelations::FillQvector(AliFlowEventSimple *anEvent)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckSettings()
{
 // Cross-check all initial settings in this method. 
 
 // a) Few cross-checks for control histograms;
 // b) Few cross-checks for flags for correlations;
 // c) 'Standard candles'.

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

 // c) 'Standard candles':
 if(fCalculateStandardCandles && !fCalculateCorrelations)
 {
  Fatal(sMethodName.Data(),"fCalculateStandardCandles && !fCalculateCorrelations");
 }

} // end of void AliFlowAnalysisWithMultiparticleCorrelations::CrossCheckSettings()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for Q-vector;
 // c) Book and nest lists for correlations;
 // d) Book and nest lists for cumulants;
 // e) Book and nest lists for weights;
 // f) Book and nest lists for nested loops;
 // g) Book and nest lists for 'standard candles'.

 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("Control Histograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);

 // b) Book and nest lists for Q-vector:
 fQvectorList = new TList();
 fQvectorList->SetName("Q-vector");
 fQvectorList->SetOwner(kTRUE);
 fHistList->Add(fQvectorList);

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
 fControlHistogramsFlagsPro = new TProfile("fControlHistogramsFlagsPro","Flags and settings for control histograms",4,0,4);
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
 fControlHistogramsList->Add(fControlHistogramsFlagsPro);

 if(!fFillControlHistograms){return;} // TBI is this safe? Well, perhaps it is if I can't implement it better...

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
   if(fFillKinematicsHist){fControlHistogramsList->Add(fKinematicsHist[rp][ppe]);}
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
  if(fFillMultDistributionsHist){fControlHistogramsList->Add(fMultDistributionsHist[rprm]);}
 } // for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]

 //  b2) Book TH2D *fMultCorrelationsHist[3]: TBI too large objects to store in this way, perhaps,
 // ...
 fMultCorrelationsHist[0] = new TH2D("Multiplicity (RP vs. POI)","Multiplicity (RP vs. POI)",nBinsMult[0],minMult[0],maxMult[0],nBinsMult[1],minMult[1],maxMult[1]);
 fMultCorrelationsHist[0]->GetXaxis()->SetTitle(xAxisTitleMult[0].Data());
 fMultCorrelationsHist[0]->GetYaxis()->SetTitle(xAxisTitleMult[1].Data());
 if(fFillMultCorrelationsHist){fControlHistogramsList->Add(fMultCorrelationsHist[0]);}

 // ...
 fMultCorrelationsHist[1] = new TH2D("Multiplicity (RP vs. REF)","Multiplicity (RP vs. REF)",nBinsMult[0],minMult[0],maxMult[0],nBinsMult[2],minMult[2],maxMult[2]);
 fMultCorrelationsHist[1]->GetXaxis()->SetTitle(xAxisTitleMult[0].Data());
 fMultCorrelationsHist[1]->GetYaxis()->SetTitle(xAxisTitleMult[2].Data());
 if(fFillMultCorrelationsHist){fControlHistogramsList->Add(fMultCorrelationsHist[1]);}

 // ...
 fMultCorrelationsHist[2] = new TH2D("Multiplicity (POI vs. REF)","Multiplicity (POI vs. REF)",nBinsMult[1],minMult[1],maxMult[1],nBinsMult[2],minMult[2],maxMult[2]);
 fMultCorrelationsHist[2]->GetXaxis()->SetTitle(xAxisTitleMult[1].Data());
 fMultCorrelationsHist[2]->GetYaxis()->SetTitle(xAxisTitleMult[2].Data());
 if(fFillMultCorrelationsHist){fControlHistogramsList->Add(fMultCorrelationsHist[2]);}

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForControlHistograms()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::FillControlHistograms(AliFlowEventSimple *anEvent)
{
 // Fill control histograms. 
 // a) Fill TH1D *fKinematicsHist[2][3];
 // b) Fill TH1D *fMultDistributionsHist[3]; 
 // c) Fill TH2D *fMultCorrelationsHist[3].  

 // a) Fill TH1D *fKinematicsHist[2][3]:
 if(fFillKinematicsHist)
 {
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
 } // if(fFillKinematicsHist)

 // b) Fill TH1D *fMultDistributionsHist[3]: 
 Double_t dMultRP = anEvent->GetNumberOfRPs(); // TBI shall I promote these 3 variables into data members? 
 Double_t dMultPOI = anEvent->GetNumberOfPOIs();
 Double_t dMultREF = anEvent->GetReferenceMultiplicity();
 Double_t dMult[3] = {dMultRP,dMultPOI,dMultREF};
 for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]
 {
  if(fFillMultDistributionsHist){fMultDistributionsHist[rprm]->Fill(dMult[rprm]);}      
 } 

 // c) Fill TH2D *fMultCorrelationsHist[3]:  
 if(fFillMultCorrelationsHist)
 {
  fMultCorrelationsHist[0]->Fill((Int_t)dMultRP,(Int_t)dMultPOI); // RP vs. POI
  fMultCorrelationsHist[1]->Fill((Int_t)dMultRP,(Int_t)dMultREF); // RP vs. refMult
  fMultCorrelationsHist[2]->Fill((Int_t)dMultPOI,(Int_t)dMultREF); // POI vs. refMult
 } // if(fFillMultCorrelationsHist)

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

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForQvector()
{
 // Book all the stuff for Q-vector.

 // a) Book the profile holding all the flags for Q-vector;
 // ...

 // a) Book the profile holding all the flags for Q-vector:
 fQvectorFlagsPro = new TProfile("fQvectorFlagsPro","Flags for Q-vector",1,0,1);
 fQvectorFlagsPro->SetTickLength(-0.01,"Y");
 fQvectorFlagsPro->SetMarkerStyle(25);
 fQvectorFlagsPro->SetLabelSize(0.03);
 fQvectorFlagsPro->SetLabelOffset(0.02,"Y");
 fQvectorFlagsPro->SetStats(kFALSE);
 fQvectorFlagsPro->SetFillColor(kGray);
 fQvectorFlagsPro->SetLineColor(kBlack);
 fQvectorFlagsPro->GetXaxis()->SetBinLabel(1,"fCalculateQvector"); fQvectorFlagsPro->Fill(0.5,fCalculateQvector); 
 fQvectorList->Add(fQvectorFlagsPro);

 if(!fCalculateQvector){return;} // TBI is this safe enough? 

 // ...

} // void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForQvector()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForCorrelations()
{
 // Book all the stuff for correlations.

 // TBI this method can be implemented in a much more civilised way. 

 // a) Book the profile holding all the flags for correlations;
 // b) Book TProfile *fCorrelationsPro[2][8] ([0=cos,1=sin][1p,2p,...,8p]).

 // a) Book the profile holding all the flags for correlations:
 fCorrelationsFlagsPro = new TProfile("fCorrelationsFlagsPro","Flags for correlations",9,0,9);
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
 fCorrelationsList->Add(fCorrelationsFlagsPro);

 if(!fCalculateCorrelations){return;} // TBI is this safe enough? 

 // b) Book TProfile *fCorrelationsPro[2][8] ([0=cos,1=sin][1p,2p,...,8p]):
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
  if(c==fDontGoBeyond){nBins[c]=1;} // TBI a bit of spaghetti here...
 } // for(Int_t c=0;c<8;c++) // [1p,2p,...,8p]
 for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 {
  for(Int_t c=0;c<fMaxCorrelator;c++) // [1p,2p,...,8p]
  {
   fCorrelationsPro[cs][c] = new TProfile(Form("%dpCorrelations%s",c+1,sCosSin[cs].Data()),"",nBins[c],0.,1.*nBins[c]);
   fCorrelationsPro[cs][c]->Sumw2();
   fCorrelationsPro[cs][c]->SetStats(kFALSE);
   fCorrelationsPro[cs][c]->SetMarkerColor(markerColor[cs]);
   fCorrelationsPro[cs][c]->SetMarkerStyle(markerStyle[cs]);
   fCorrelationsList->Add(fCorrelationsPro[cs][c]);
  } // for(Int_t c=0;c<8;c++) // [1p,2p,...,8p]
 } // for(Int_t cs=0;cs<2;cs++) // [0=cos,1=sin]
 // Set all bin labels: TBI this can be implemented better, most likely...
 Int_t binNo[8]; for(Int_t c=0;c<fMaxCorrelator;c++){binNo[c]=1;} // TBI hardwired 8, shall be fMaxCorrelator
 for(Int_t n1=-fMaxHarmonic;n1<=fMaxHarmonic;n1++) 
 {
  if(fSkipZeroHarmonics && 0==n1){continue;}
  if(fCalculateAll)
  {
   fCorrelationsPro[0][0]->GetXaxis()->SetBinLabel(binNo[0],Form("Cos(%d)",n1));
   fCorrelationsPro[1][0]->GetXaxis()->SetBinLabel(binNo[0]++,Form("Sin(%d)",n1));
   nToBeFilled[0]++; 
  }
  if(1==fDontGoBeyond){continue;}
  for(Int_t n2=n1;n2<=fMaxHarmonic;n2++) 
  {
   if(fSkipZeroHarmonics && 0==n2){continue;}
   if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2) || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2)) 
      || (fCalculateSameIsotropic && 0==n1+n2 &&  TMath::Abs(n1)==TMath::Abs(n2)))
   {  
    fCorrelationsPro[0][1]->GetXaxis()->SetBinLabel(binNo[1],Form("Cos(%d,%d)",n1,n2));
    fCorrelationsPro[1][1]->GetXaxis()->SetBinLabel(binNo[1]++,Form("Sin(%d,%d)",n1,n2));
    nToBeFilled[1]++; 
   }
   if(2==fDontGoBeyond){continue;}
   for(Int_t n3=n2;n3<=fMaxHarmonic;n3++) 
   {
    if(fSkipZeroHarmonics && 0==n3){continue;}
    if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3) || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3))
       || (fCalculateSameIsotropic && 0==n1+n2+n3 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3)))
    {  
     fCorrelationsPro[0][2]->GetXaxis()->SetBinLabel(binNo[2],Form("Cos(%d,%d,%d)",n1,n2,n3));
     fCorrelationsPro[1][2]->GetXaxis()->SetBinLabel(binNo[2]++,Form("Sin(%d,%d,%d)",n1,n2,n3));
     nToBeFilled[2]++; 
    }
    if(3==fDontGoBeyond){continue;}
    for(Int_t n4=n3;n4<=fMaxHarmonic;n4++) 
    {
     if(fSkipZeroHarmonics && 0==n4){continue;}
     if(fCalculateAll || (fCalculateIsotropic && 0==n1+n2+n3+n4) || (fCalculateSame && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4))
       || (fCalculateSameIsotropic && 0==n1+n2+n3+n4 && TMath::Abs(n1)==TMath::Abs(n2) && TMath::Abs(n1)==TMath::Abs(n3) && TMath::Abs(n1)==TMath::Abs(n4)))
     {   
      fCorrelationsPro[0][3]->GetXaxis()->SetBinLabel(binNo[3],Form("Cos(%d,%d,%d,%d)",n1,n2,n3,n4));
      fCorrelationsPro[1][3]->GetXaxis()->SetBinLabel(binNo[3]++,Form("Sin(%d,%d,%d,%d)",n1,n2,n3,n4)); 
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
       fCorrelationsPro[0][4]->GetXaxis()->SetBinLabel(binNo[4],Form("Cos(%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5));
       fCorrelationsPro[1][4]->GetXaxis()->SetBinLabel(binNo[4]++,Form("Sin(%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5));
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
        fCorrelationsPro[0][5]->GetXaxis()->SetBinLabel(binNo[5],Form("Cos(%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6));
        fCorrelationsPro[1][5]->GetXaxis()->SetBinLabel(binNo[5]++,Form("Sin(%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6));
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
         fCorrelationsPro[0][6]->GetXaxis()->SetBinLabel(binNo[6],Form("Cos(%d,%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6,n7));
         fCorrelationsPro[1][6]->GetXaxis()->SetBinLabel(binNo[6]++,Form("Sin(%d,%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6,n7));
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
          fCorrelationsPro[0][7]->GetXaxis()->SetBinLabel(binNo[7],Form("Cos(%d,%d,%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6,n7,n8));
          fCorrelationsPro[1][7]->GetXaxis()->SetBinLabel(binNo[7]++,Form("Sin(%d,%d,%d,%d,%d,%d,%d,%d)",n1,n2,n3,n4,n5,n6,n7,n8));
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
   fCorrelationsPro[cs][c]->SetTitle(Form("%d-p correlations, %s terms, %d/%d in total",c+1,sCosSin[cs].Data(),nToBeFilled[c],nBinsTitle[c]));
   fCorrelationsPro[cs][c]->GetXaxis()->SetRangeUser(0.,fCorrelationsPro[cs][c]->GetBinLowEdge(nToBeFilled[c]+1));
  }
 } 

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

 // TBI this method is just ugly, who implemented it like this... 

 // a) Set default harmonic values; 
 // b) Book the profile holding all the flags for nested loops;
 // c) Book the profile holding all results for nested loops (cosine);
 // d) Book the profile holding all results for nested loops (sinus).

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
 // b) Book the histogram holding all results for 'standard candles';
 // c) Book 2D profile holding products of correlations (needed for error propagation).

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::BookEverythingForStandardCandles()";

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
 Int_t nBins = fMaxHarmonic*(fMaxHarmonic-1)/2;
 fStandardCandlesHist = new TH1D("fStandardCandlesHist","'Standard candles'",nBins,0.,1.*nBins); 
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

 // c) Book 2D profile holding products of correlations (needed for error propagation):
 Int_t nBins2p = fCorrelationsPro[0][1]->GetXaxis()->GetNbins(), nBins2pFilled = 0;
 Int_t nBins4p = fCorrelationsPro[0][3]->GetXaxis()->GetNbins(), nBins4pFilled = 0;
 for(Int_t b=1;b<=nBins2p;b++)
 {
  TString binLabel2p = fCorrelationsPro[0][1]->GetXaxis()->GetBinLabel(b);
  if(!binLabel2p.CompareTo("")){break;} 
  nBins2pFilled++;
 }
 for(Int_t b=1;b<=nBins4p;b++)
 {
  TString binLabel4p = fCorrelationsPro[0][3]->GetXaxis()->GetBinLabel(b);
  if(!binLabel4p.CompareTo("")){break;}
  nBins4pFilled++;
 }
 Int_t nBins2D = nBins2pFilled + nBins4pFilled;
 if(0==nBins2D){Fatal(sMethodName.Data(),"0==nBins2D");} // well, just in case...
 fProductsPro2D = new TProfile2D("fProductsPro2D","Products of correlations",nBins2D,0.,nBins2D,nBins2D,0.,nBins2D);
 fProductsPro2D->SetStats(kFALSE);
 fProductsPro2D->Sumw2();
 for(Int_t b=1;b<=nBins2D;b++)
 {
  TString binLabel = "";
  if(b>=1 && b<=nBins2pFilled)
  {
   binLabel = fCorrelationsPro[0][1]->GetXaxis()->GetBinLabel(b);
  }
  else if(b>nBins2pFilled && b<=nBins4pFilled+nBins2pFilled) // TBI landmine ?
  {
   binLabel = fCorrelationsPro[0][3]->GetXaxis()->GetBinLabel(b-nBins2pFilled); // TBI another landmine ?
  }
  fProductsPro2D->GetXaxis()->SetBinLabel(b,binLabel.Data());
  fProductsPro2D->GetYaxis()->SetBinLabel(b,binLabel.Data());
 } // for(Int_t b=1;b<=nBins2D;b++)
 fStandardCandlesList->Add(fProductsPro2D);

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
 // f) Get pointers for 'standard candles'.

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

 // c) Get pointers for control histograms:
 this->GetPointersForControlHistograms(); 

 // d) Get pointers for Q-vector:
 this->GetPointersForQvector(); 

 // e) Get pointers for correlations:
 this->GetPointersForCorrelations(); 

 // f) Get pointers for 'standard candles':
 this->GetPointersForStandardCandles(); 
  
} // void AliFlowAnalysisWithMultiparticleCorrelations::GetOutputHistograms(TList *histList)

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForQvector()
{
 // Get pointers for Q-vector objects.

 // a) Get pointer for fQvectorList; TBI
 // b) Get pointer for fQvectorFlagsPro; TBI
 // c) Set again all flags; TBI

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForQvector()";

 // a) Get pointer for fQvectorList: TBI
 fQvectorList = dynamic_cast<TList*>(fHistList->FindObject("Q-vector")); 
 if(!fQvectorList){Fatal(sMethodName.Data(),"fQvectorList");}

 // b) Get pointer for fQvectorFlagsPro: TBI
 fQvectorFlagsPro = dynamic_cast<TProfile*>(fQvectorList->FindObject("fQvectorFlagsPro"));
 if(!fQvectorFlagsPro){Fatal(sMethodName.Data(),"fQvectorFlagsPro");}

 // c) Set again all flags: TBI
 fCalculateQvector = (Bool_t)fQvectorFlagsPro->GetBinContent(1);

 if(!fCalculateQvector){return;} // TBI is this safe enough

} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForQvector()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForStandardCandles()
{
 // Get pointers for 'standard candles'.

 // a) Get pointer for fStandardCandlesList; TBI
 // b) Get pointer for fStandardCandlesFlagsPro; TBI
 // c) Set again all flags; TBI
 // d) Get pointer TH1D *fStandardCandlesHist; TBI
 // e) Get pointer TProfile2D *fProductsPro2D. TBI

 TString sMethodName = "AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForStandardCandles()";

 // a) Get pointer for fStandardCandlesList: TBI
 fStandardCandlesList = dynamic_cast<TList*>(fHistList->FindObject("Standard Candles"));
 if(!fStandardCandlesList){Fatal(sMethodName.Data(),"fStandardCandlesList");}

 // b) Get pointer for fStandardCandlesFlagsPro: TBI
 fStandardCandlesFlagsPro = dynamic_cast<TProfile*>(fStandardCandlesList->FindObject("fStandardCandlesFlagsPro"));
 if(!fStandardCandlesFlagsPro){Fatal(sMethodName.Data(),"fStandardCandlesFlagsPro");}

 // c) Set again all flags: TBI
 fCalculateStandardCandles = fStandardCandlesFlagsPro->GetBinContent(1);

 if(!fCalculateStandardCandles){return;} // TBI is this safe enough

 // d) Get pointer TH1D *fStandardCandlesHist: TBI
 fStandardCandlesHist = dynamic_cast<TH1D*>(fStandardCandlesList->FindObject("fStandardCandlesHist"));

 // e) Get pointer TProfile2D *fProductsPro2D: TBI
 fProductsPro2D = dynamic_cast<TProfile2D*>(fStandardCandlesList->FindObject("fProductsPro2D"));
 
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

 if(!fFillControlHistograms){return;} // TBI is this safe enough

 // d) Get pointers to fKinematicsHist[2][3]: TBI
 TString name[2][3] = {{"RP,phi","RP,pt","RP,eta"},{"POI,phi","POI,pt","POI,eta"}}; // [RP,POI][phi,pt,eta]
 for(Int_t rp=0;rp<2;rp++) // [RP,POI]
 {
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
  fMultDistributionsHist[rprm] = dynamic_cast<TH1D*>(fControlHistogramsList->FindObject(nameMult[rprm].Data()));
  if(!fMultDistributionsHist[rprm] && fFillMultDistributionsHist){Fatal(sMethodName.Data(),"%s",nameMult[rprm].Data());} // TBI 
 } // for(Int_t rprm=0;rprm<3;rprm++) // [RP,POI,reference multiplicity]

 // f) Get pointers to TH2D *fMultCorrelationsHist[3]: TBI automatize the things here...
 fMultCorrelationsHist[0] = dynamic_cast<TH2D*>(fControlHistogramsList->FindObject("Multiplicity (RP vs. POI)"));
 if(!fMultCorrelationsHist[0] && fFillMultCorrelationsHist){Fatal(sMethodName.Data(),"Multiplicity (RP vs. POI)");} // TBI 
 fMultCorrelationsHist[1] = dynamic_cast<TH2D*>(fControlHistogramsList->FindObject("Multiplicity (RP vs. REF)"));
 if(!fMultCorrelationsHist[1] && fFillMultCorrelationsHist){Fatal(sMethodName.Data(),"Multiplicity (RP vs. REF)");} // TBI 
 fMultCorrelationsHist[2] = dynamic_cast<TH2D*>(fControlHistogramsList->FindObject("Multiplicity (POI vs. REF)"));
 if(!fMultCorrelationsHist[2] && fFillMultCorrelationsHist){Fatal(sMethodName.Data(),"Multiplicity (POI vs. REF)");} // TBI 

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
 fCalculateIsotropic = fCorrelationsFlagsPro->GetBinContent(4);
 fCalculateSame = fCorrelationsFlagsPro->GetBinContent(5);
 fSkipZeroHarmonics = fCorrelationsFlagsPro->GetBinContent(6);
 fCalculateSameIsotropic = fCorrelationsFlagsPro->GetBinContent(7);
 fCalculateAll = fCorrelationsFlagsPro->GetBinContent(8);
 fDontGoBeyond = (Int_t)fCorrelationsFlagsPro->GetBinContent(9);

 if(!fCalculateCorrelations){return;} // TBI is this safe enough, that is the question...

 // d) Get pointers to TProfile *fCorrelationsPro[2][8]:
 TString sCosSin[2] = {"Cos","Sin"}; 
 for(Int_t cs=0;cs<2;cs++)
 {
  for(Int_t c=0;c<fMaxCorrelator;c++)
  {
   fCorrelationsPro[cs][c] = dynamic_cast<TProfile*>(fCorrelationsList->FindObject(Form("%dpCorrelations%s",c+1,sCosSin[cs].Data())));
   if(!fCorrelationsPro[cs][c]){Fatal(sMethodName.Data(),"%dpCorrelations%s",c+1,sCosSin[cs].Data());} 
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::GetPointersForCorrelations()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForQvector()
{
 // Initialize all arrays for Q-vector.

 for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++) 
 {
  for(Int_t p=0;p<fMaxCorrelator+1;p++) 
  {
   fQvector[h][p] = TComplex(0.,0.);
  }
 }

} // void AliFlowAnalysisWithMultiparticleCorrelations::InitializeArraysForQvector()

//=======================================================================================================================

void AliFlowAnalysisWithMultiparticleCorrelations::ResetQvector()
{
 // Reset all Q-vector components to zero before starting a new event. 

 for(Int_t h=0;h<fMaxHarmonic*fMaxCorrelator+1;h++) 
 {
  for(Int_t p=0;p<fMaxCorrelator+1;p++)
  {
   fQvector[h][p] = TComplex(0.,0.);
  } 
 } 

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

 TComplex seven = 0.*Q(n1-n2+n2-n3+n4-n5+n6-n7,1) + TComplex(1.,1.); // TBI implement the actual Eq. 

 return seven;

} // end of TComplex AliFlowAnalysisWithMultiparticleCorrelations::Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7)

//=======================================================================================================================

TComplex AliFlowAnalysisWithMultiparticleCorrelations::Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)
{
 // Generic eight-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8)]>.

 TComplex eight = 0.*Q(n1-n2+n2-n3+n4-n5+n6-n7+n8,1) + TComplex(1.,1.); // TBI implement the actual Eq. 

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

 if(!fPhiWeightsHist){Fatal("AliFlowAnalysisWithMultiparticleCorrelations::PhiWeight(const Double_t &dPhi)","fPhiWeightsHist");}

 Double_t wPhi = fPhiWeightsHist->GetBinContent(fPhiWeightsHist->FindBin(dPhi));

 return wPhi;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::PhiWeight(const Double_t &dPhi)

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::PtWeight(const Double_t &dPt)
{
 // Determine pt weight for a given pt. 

 if(!fPtWeightsHist){Fatal("AliFlowAnalysisWithMultiparticleCorrelations::PtWeight(const Double_t &dPt)","fPtWeightsHist");}

 Double_t wPt = fPtWeightsHist->GetBinContent(fPtWeightsHist->FindBin(dPt));

 return wPt;

} // Double_t AliFlowAnalysisWithMultiparticleCorrelations::PtWeight(const Double_t &dPt)

//=======================================================================================================================

Double_t AliFlowAnalysisWithMultiparticleCorrelations::EtaWeight(const Double_t &dEta)
{
 // Determine eta weight for a given eta. 

 if(!fEtaWeightsHist){Fatal("AliFlowAnalysisWithMultiparticleCorrelations::EtaWeight(const Double_t &dEta)","fEtaWeightsHist");}

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


