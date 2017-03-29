/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisMuMuFnorm.h"

#include "AliAnalysisMuMuGraphUtil.h"
#include "AliAnalysisMuMuResult.h"
#include "AliAnalysisTriggerScalers.h"
#include "AliCounterCollection.h"
#include "AliLog.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuMuConfig.h"
#include "Riostream.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH1.h"
#include "TList.h"
#include "TMap.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TPaveText.h"
#include "TStyle.h"
#include <cassert>
#include <numeric>

ClassImp(AliAnalysisMuMuFnorm)

using namespace std;

//_____________________________________________________________________________
AliAnalysisMuMuFnorm::AliAnalysisMuMuFnorm(AliCounterCollection& cc,AliAnalysisMuMuConfig& cf,
                                           AliAnalysisMuMuFnorm::ETriggerType refTriggerType,
                                           const char* ocdbpath,
                                           Bool_t compactGraphs) :
TObject(),
fCounterCollection(cc),
fConfig(cf),
fMergeableCollection(0x0),
fIsOwner(kTRUE),
fOCDBPath(ocdbpath),
fResult(0x0),
fIsCompactGraphs(compactGraphs),
fReferenceTriggerType(refTriggerType)
{
  // ctor
}

//_____________________________________________________________________________
AliAnalysisMuMuFnorm::~AliAnalysisMuMuFnorm()
{
  // dtor
  if ( fIsOwner )
  {
    delete fMergeableCollection;
  }
}

// //_____________________________________________________________________________
// void AliAnalysisMuMuFnorm::ComputeFnorm()
// {
//   /// Compute the REF to CINT ratio(s)
//   ///
//   /// Using offline method
//   ///   - in one go CINT/REF
//   ///   - in two steps CINT/CMSL and CMSL/REF
//   ///
//   /// Using scaler method
//   ///   - bare scaler values
//   ///   - scaler values corrected for pile-up
//   ///   - scaler values corrected for pile-up and physics selection

// //   const ETriggerType triggerTypes[] = { kMB, kMUL, kMSL, kMSH };
// //   const Bool_t trueFalse[] = { kTRUE, kFALSE };
// //   // Call ComputeNofEvent for every combination possible
// //   for ( Int_t i = 0; i < 4; ++i )
// //   {
// //     for ( Int_t pileup = 0; pileup < 2; ++pileup )
// //     {
// //       for ( Int_t ps = 0; ps < 2; ++ps )
// //       {
// //         ComputeNofEvents(triggerTypes[i],trueFalse[pileup],ps);
// //       }
// //     }
// //   }
// //   ComputeFnormOffline(1,kFALSE,0);
// //   ComputeFnormOffline(1,kFALSE,1);
// //   ComputeFnormOffline(1,kTRUE,1);

// //   ComputeFnormOffline(2,kFALSE,0);
// //   ComputeFnormOffline(2,kFALSE,1);
// //   ComputeFnormOffline(2,kTRUE,1);

// // //  ComputeFnormOffline(2,kFALSE,2);
// // //  ComputeFnormOffline(2,kTRUE,2);

// //   ComputeFnormScalers(kFALSE,0);
// //   ComputeFnormScalers(kTRUE,0);
// //   ComputeFnormScalers(kTRUE,1);
// // //  ComputeFnormScalers(kTRUE,2);

// //   WeightedMeanGraphs("Offline");
// //   WeightedMeanGraphs("Scalers");
// //   WeightedMeanGraphs("FnormOffline2PUPS,FnormOffline1PUPS","FnormOffline12PUPS");

// //   WeightedMeanGraphs("FnormOffline2PUPS,FnormScalersPUPS","FnormBest2");

// //   ComputeGraphRelDif("FnormOffline2PUPS","FnormScalersPUPS");

// //   ComputeGraphRelDif("FnormOffline2PUPS","FnormOffline2");
// //   ComputeGraphRelDif("FnormOffline2PUPS","FnormOffline2PS");

// //   ComputeGraphRelDif("CorrectionPSMB","CorrectionPSREF");

// // //  for ( Int_t i = 0; i < 4; ++i )
// // ///  {
// //     TString triggerEvents;

// // //  triggerEvents.Form("NofEvent%sPUPS",GetTriggerTypeName(triggerTypes[i]).Data());
// //   triggerEvents.Form("NofEvent%sPUPS",GetTriggerTypeName(fReferenceTriggerType).Data());

// //   MultiplyGraphs(triggerEvents.Data(),"FnormBest2","NMBeqBest2");

// //     MultiplyGraphs(triggerEvents.Data(),"FnormOffline2PUPS","NMBeqOffline2PUPS");
// //   MultiplyGraphs(triggerEvents.Data(),"FnormScalersPUPS","NMBeqScalersPUPS");
// // //  }

// // //  MultiplyGraphs(Form("NofEvent%sPUTS",GetTriggerTypeName(fReferenceTriggerType).Data()),"FnormOffline2PUTS","NMBeqOffline2PUTS");
// // //  MultiplyGraphs(Form("NofEvent%sPUTS",GetTriggerTypeName(fReferenceTriggerType).Data()),"FnormOffline2TS","NMBeqOffline2TS");

// //   ComputeResultsFromGraphs();

// //   AliAnalysisMuMuResult* result = GetResult("Fnorm");
// //   if (result)
// //   {
// //     result->Exclude("*");
// //     result->Include("FnormBest2");
// //   }
// }

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputeCorrectionFactors(Int_t eventSelectionCorrected)
{
  /// Compute individual graphs for the correction factors (PS_REF, PS_CINT,
  /// F_pile-up,PS_CINT/PS_REF) used in the computation of (some) Fnorm factors

  TString graphName(Form("CorrectionGlobal%s",GetEventSelectionName(eventSelectionCorrected).Data()));;

  if ( GetGraph(graphName) )
  {
    // insure we compute it only once
    return;
  }

  AliDebug(2,"");// ??


  vector<double> vx;
  vector<double> vxerr;
  vector<double> vy;
  vector<double> vyerr;

  vector<double> vyGlobal;
  vector<double> vyGlobalErr;

  const ETriggerType triggerTypes[] = { kMB, kMUL, kMSL, kMSH };

  for ( Int_t i = 0; i < 4; ++i )
  {
    ComputeEventSelectionGraph(triggerTypes[i],eventSelectionCorrected);// Problem from here
    ComputePileUpGraph(triggerTypes[i],eventSelectionCorrected);
  }

  TGraphErrors* gPSCINT = GetGraph(Form("Correction%s%s",GetEventSelectionName(eventSelectionCorrected).Data(),GetTriggerTypeName(AliAnalysisMuMuFnorm::kMB).Data()));

  TGraphErrors* gPSREF = GetGraph(Form("Correction%s%s", GetEventSelectionName(eventSelectionCorrected).Data(), GetTriggerTypeName(fReferenceTriggerType).Data()));

  TGraphErrors* gPU = GetGraph(Form("CorrectionPU%s%s",GetEventSelectionName(eventSelectionCorrected).Data(), GetTriggerTypeName(AliAnalysisMuMuFnorm::kMB).Data()));

  if ( !gPSCINT || !gPSREF || !gPU )
  {
    AliError("Did not get the relevant graphs. Cannot work");
    return;
  }

  for ( Int_t i = 0; i < gPSCINT->GetN(); ++i )
  {
    Double_t x,y,yerr,yGlobal,yGlobalErr;

    gPSCINT->GetPoint(i,x,y);

    if ( fIsCompactGraphs )
    {
      x = TString(gPSCINT->GetXaxis()->GetBinLabel(i)).Atoi();
    }

    yGlobal = gPSCINT->GetY()[i] * gPU->GetY()[i] / gPSREF->GetY()[i];

    yGlobalErr = yGlobal*AliAnalysisMuMuResult::ErrorABC(gPSCINT->GetY()[i],gPSCINT->GetEY()[i],
                                                         gPSREF->GetY()[i],gPSREF->GetEY()[i],
                                                         gPU->GetY()[i],gPU->GetEY()[i]);

    y = gPSCINT->GetY()[i] / gPSREF->GetY()[i];
    yerr = y * AliAnalysisMuMuResult::ErrorAB(gPSCINT->GetY()[i],gPSCINT->GetEY()[i],
                                              gPSREF->GetY()[i],gPSREF->GetEY()[i]);

    vx.push_back(x);
    vxerr.push_back(gPSCINT->GetEX()[i]);

    vyGlobal.push_back(yGlobal);
    vyGlobalErr.push_back(yGlobalErr);

    vy.push_back(y);
    vyerr.push_back(yerr);
  }

  TString name(Form("Correction%sRatio",GetEventSelectionName(eventSelectionCorrected).Data()));
  TString title(Form("%s_MB/%s_%s",GetEventSelectionName(eventSelectionCorrected).Data(),GetTriggerTypeName(fReferenceTriggerType).Data(),
                     GetEventSelectionName(eventSelectionCorrected).Data()));

  CreateAndAddGraph(name,title,vx,vxerr,vy,vyerr);

  title = TString::Format("%s_MB x Fpile-up / %s_%s ",GetEventSelectionName(eventSelectionCorrected).Data(),GetTriggerTypeName(fReferenceTriggerType).Data(),GetEventSelectionName(eventSelectionCorrected).Data());

  CreateAndAddGraph(graphName,title,vx,vxerr,vyGlobal,vyGlobalErr);
}


// //_____________________________________________________________________________
// void AliAnalysisMuMuFnorm::ComputeLuminosity(Int_t nstep, Bool_t pileUpCorrected, Int_t eventSelectionCorrected)
// {
//   /// Compute the run by run REF cross-section as sigma(ref) = sigma_(vdm)/FNorm and the integrated luminosity
//   /// This method ment to run originaly for PP system.
//   /// @argument nsteps : The FNorm reference graph, could be 1 (offline1), 2 (offline1) or 0 (Scaler)

//   TString centrality(fConfig.First(fConfig.CentralitySelectionKey()));

//   TString name("Luminosity");
//   TString title(Form("%s Cross-section",GetTriggerTypeName(fReferenceTriggerType).Data()));
//   TString nameFNorm("Fnorm");
//   TString nameDirectory;

//   // COnfigure name
//   if (nstep==1){
//     name          += "Offline1";
//     nameFNorm     += "Offline1";
//     title         += " with pile-up correction";
//     nameDirectory = "Offline"
//   }
//   else if (nstep==2){
//     name          += "Offline2";
//     nameFNorm     += "Offline2";
//     title         += " with pile-up correction";
//     nameDirectory = "Offline2"
//   }
//   else if (nstep==0){
//     name          += "Scaler";
//     nameFNorm     += "Scaler";
//     title         += " with pile-up correction";
//     nameDirectory = "Scaler"
//   }
//   if (pileUpCorrected){
//     name      += "PU";
//     nameFNorm += "PU";
//     title     += " with pile-up correction";
//   }
//   if (eventSelectionCorrected==1){
//     name      += "PS";
//     nameFNorm += "PS";
//     title     += " with (ps) purity corrections";
//   }
//   if (eventSelectionCorrected==2){
//     name      += "TS";
//     nameFNorm += "TS";
//     title     += " with (ts) purity corrections";
//   }

//   if (GetGraph(Form("Luminosity/%s",name.Data())))return;

//   AliDebug(2,name);

//   // Get the FNorm Graph
//   TGraphErrors* gFnorm = 0x0;
//   gFnorm = GetGraph(Form("%s",nameFNorm.Data());

//   if(!gFnorm) {
//    printf("cannot get  the FNorm histo %s\n",nameFNorm.Data());
//    return;
//   }


//   const set<int>& runs = RunNumbers();

//   vector<double> vx;
//   vector<double> vxerr;
//   vector<double> vy;
//   vector<double> vyerr;

//   Double_t purityREF(1.0);
//   Double_t purityMB(1.0);
//   Double_t purityREFerror(00);
//   Double_t purityMBerror(0.0);

//   // compute the per run values
//   for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
//   {
//     Int_t runNumber = *it;

//     TString mbTrigger = GetTriggerClassName(kMB,runNumber);
// //    TString refTrigger = GetTriggerClassName(fReferenceTriggerType,runNumber);

//     purityMB=purityREF=1.0;
//     purityMBerror=purityREFerror=0.0;

//     if (eventSelectionCorrected>0)
//     {
//       ComputePurityFactorForScalerGraph(kMB,eventSelectionCorrected,&centrality);
//       ComputePurityFactorForScalerGraph(fReferenceTriggerType,eventSelectionCorrected,&centrality);

//       TGraphErrors* gpsfactorMB  = GetGraph(Form("PurityFactorForScaler%s%s_%s",GetEventSelectionName(eventSelectionCorrected).Data(), GetTriggerTypeName(kMB).Data(),centrality.Data()));
//       TGraphErrors* gpsfactorREF = GetGraph(Form("PurityFactorForScaler%s%s_%s",GetEventSelectionName(eventSelectionCorrected).Data(), GetTriggerTypeName(fReferenceTriggerType).Data(),centrality.Data()));

//       GetValueAndErrorFromGraph(gpsfactorMB,runNumber,purityMB,purityMBerror);
//       GetValueAndErrorFromGraph(gpsfactorREF,runNumber,purityREF,purityREFerror);
//     }

//     if (purityMB<=0.0)
//     {
//       AliError(Form("Got purity=%e for MB for run %9d",purityMB,runNumber));
//       continue;
//     }

//     TGraphErrors* gl0bMB = GetGraph(Form("L0B%s",GetTriggerTypeName(kMB).Data()));
//     TGraphErrors* gl0bREF = GetGraph(Form("L0B%s",GetTriggerTypeName(fReferenceTriggerType).Data()));

//     Double_t L0bMB,L0bMBError;
//     Double_t L0bREF,L0bREFError;

//     GetValueAndErrorFromGraph(gl0bMB,runNumber,L0bMB,L0bMBError);
//     GetValueAndErrorFromGraph(gl0bREF,runNumber,L0bREF,L0bREFError);

//     Double_t pileUpFactor(1.0);
//     Double_t pileUpFactorError(0.0);

//     if (pileUpCorrected)
//     {
//       ComputePileUpGraph(kMB,eventSelectionCorrected);

//       TGraphErrors* gpu = GetGraph((Form("CorrectionPU%s%s",GetEventSelectionName(eventSelectionCorrected).Data(),GetTriggerTypeName(kMB).Data())));

//       GetValueAndErrorFromGraph(gpu,runNumber,pileUpFactor,pileUpFactorError);
//     }

//     Double_t value;
//     Double_t error;

//     ScalerFnorm(value,error,
//                 L0bREF,purityREF,purityREFerror,
//                 L0bMB,purityMB,purityMBerror,
//                 pileUpFactor,pileUpFactorError);

//     if ( value > 0.0 )
//     {
//       vx.push_back(1.0*runNumber);
//       vxerr.push_back(0.5);
//       vy.push_back(value);
//       vyerr.push_back(error);
//     }
//   }

//   CreateAndAddGraph(name,title,vx,vxerr,vy,vyerr);
// }

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputeFnormOffline(Int_t nstep, Bool_t pileUpCorrected, Int_t eventSelectionCorrected)
{
  /// Compute MB to REF ratio using offline method, either in 1 or 2 steps
  /// @argument nstep                    1 = one step method ; 2 = two steps method
  /// @argument pileUpCorrected          wheither or not we apply the pile-up correction factor
  /// @argument eventSelectionCorrected  wheither or not we apply the purity factor correction and how : 1 = ps ; 2= ts

  TString name("FnormOffline");
  TString title("Computed using offline information");
  TString refInput = Form("0%s",GetTriggerTypeName(fReferenceTriggerType).Data());
  printf("refInput = %s\n",refInput.Data() );

  // Select graph title
  if (nstep==1){
    name  += "1";
    title += Form(" in one step (CINT/CINT&%s)",refInput.Data());
  }
  else{
    name  += "2";
    title += Form(" in two steps (CMSL/CMSL&%s x CINT/CINT&0MSL)",refInput.Data());
  }
  if (pileUpCorrected){
    name += "PU";
    title += " with pile-up correction";
  }
  if (eventSelectionCorrected==1){
    name += "PS";
    title += " with (ps) purity corrections";
  }
  else if ( eventSelectionCorrected==2 ){
    name += "TS";
    title += " with (ts) purity corrections";
  }
  // insure we're computing it only once
  if ( GetGraph(name) )return;

  AliDebug(2,name);

  vector<double> vx;
  vector<double> vxerr;
  vector<double> vy;
  vector<double> vyerr;

  const set<int>& runs = RunNumbers();

  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it ){

    Int_t runNumber = *it;
    // Get trigger name
    TString mbTrigger   = GetTriggerClassName(kMB,runNumber);
    TString muonTrigger = GetTriggerClassName(kMSL,runNumber);
    if (!mbTrigger.Length()){
      AliError(Form("Cannot get MB trigger for run %d",runNumber));
      continue;
    }
    AliDebug(1,Form("mb trigger  : %s\n", mbTrigger.Data()));
    AliDebug(1,Form("muon trigger: %s\n", muonTrigger.Data()));

    Double_t nofMB = GetSum(mbTrigger.Data(),runNumber,eventSelectionCorrected);
    Double_t nofMSL(0.0);
    Double_t nofMSLw0REF(0.0);

    if ( nstep==2 ){
      nofMSL = GetSum(muonTrigger.Data(),runNumber,eventSelectionCorrected);
      TString counterName = muonTrigger;
      if ( fReferenceTriggerType != kMSL ) counterName += Form("&%s",refInput.Data());
      nofMSLw0REF = GetSum(counterName,runNumber,eventSelectionCorrected);
    }

    Double_t nofMBw0REF = GetSum(Form("%s&%s",mbTrigger.Data(),refInput.Data()),runNumber,eventSelectionCorrected);
    Double_t nofMBw0MSL = GetSum(Form("%s&0MSL",mbTrigger.Data()),runNumber,eventSelectionCorrected);
    if(mbTrigger.Contains("CINT7-B-NOPF-MUFAST&0TVX")) nofMBw0REF = GetSum(Form("CINT7-B-NOPF-MUFAST&%s",refInput.Data()),runNumber,eventSelectionCorrected); //Specific to LHC15n
    if(mbTrigger.Contains("CINT7-B-NOPF-MUFAST&0TVX")) nofMBw0MSL = GetSum("CINT7-B-NOPF-MUFAST&0MSL",runNumber,eventSelectionCorrected); //Specific to LHC15n

    // if ( !nofMBw0REF && nstep == 1) continue;
    // if ( !nofMBw0MSL && nstep == 2 ) continue;

    Double_t purityMB(1.0);
    Double_t purityMBerror(0.0);

    if ( eventSelectionCorrected > 0 ){
      ComputeEventSelectionGraph(kMB,eventSelectionCorrected);

      TGraphErrors* gps= 0x0;
      gps = GetGraph(Form("Correction%s%s",GetEventSelectionName(eventSelectionCorrected).Data(),GetTriggerTypeName(kMB).Data()));
      if (!gps)printf("Could not found Correction%s%s \n",GetEventSelectionName(eventSelectionCorrected).Data(),GetTriggerTypeName(kMB).Data());

      GetValueAndErrorFromGraph(gps,runNumber,purityMB,purityMBerror);
    }

    Double_t pileUpFactor(1.0);
    Double_t pileUpFactorError(0.0);

    if (pileUpCorrected){
      ComputePileUpGraph(kMB,eventSelectionCorrected);

      TGraphErrors* gpu = GetGraph(Form("CorrectionPU%s%s",GetEventSelectionName(eventSelectionCorrected).Data(),GetTriggerTypeName(kMB).Data()));

      GetValueAndErrorFromGraph(gpu,runNumber,pileUpFactor,pileUpFactorError);

      nofMB *= pileUpFactor;
    }

    double value = nofMBw0REF > 0.0 ? nofMB/nofMBw0REF : 0.0;
    double error = value*AliAnalysisMuMuResult::ErrorABC(nofMB,TMath::Sqrt(nofMB),
                                                              nofMBw0REF,TMath::Sqrt(nofMBw0REF),
                                                              pileUpFactor,pileUpFactorError);

    if ( nstep == 2 ){
      value =  nofMBw0MSL > 0.0 ? (nofMB/nofMSLw0REF)*(nofMSL/nofMBw0MSL) : 0.0;
      AliDebug(1,Form("RUN %09d %d-%d-%d value=%e nofMB %e nofMSLw%s %e nofMSL %e nofMBw0MSL %e",
                      runNumber,nstep,pileUpCorrected,eventSelectionCorrected,
                      value,nofMB,refInput.Data(),nofMSLw0REF,nofMSL,nofMBw0MSL));

      error = value*AliAnalysisMuMuResult::ErrorABCD(nofMB,TMath::Sqrt(nofMB),
                                                          nofMSLw0REF,TMath::Sqrt(nofMSLw0REF),
                                                          nofMSL,TMath::Sqrt(nofMSL),
                                                          nofMBw0MSL,TMath::Sqrt(nofMBw0MSL));
    }

    vx.push_back(1.0*runNumber);
    vxerr.push_back(0.5);
    vy.push_back(value);
    vyerr.push_back(error);
  }

  CreateAndAddGraph(name,title,vx,vxerr,vy,vyerr);
}


//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputeFnormScalers(Bool_t pileUpCorrected,
                                               Int_t eventSelectionCorrected)
{
  /// Compute the MB to REF ratio using the scalers method (from OCDB)
  ///
  /// i.e. Fnorm = L0B(MB) x PS(MB) x Fpile-up / ( L0B(REF) x PS(REF) )
  ///
  /// where
  ///
  ///   - MB is the minbias trigger
  ///
  ///
  ///   -REF is the fReferenceTrigger
  ///
  ///
  ///   -PS is the fraction of events selected by the physics selection
  ///
  /// The correction factor (the two PS and one Fpile-up) are
  /// taken from graphs computed in other methods
  ///

  TString centrality(fConfig.First(fConfig.CentralitySelectionKey()));


  TString name("FnormScalers");
  TString title("Computed using OCDB scalers");

  if (pileUpCorrected)
  {
    name += "PU";
    title += " with pile-up correction";
  }
  if (eventSelectionCorrected==1)
  {
    name += "PS";
    title += " with (ps) purity corrections";
  }
  if (eventSelectionCorrected==2)
  {
    name += "TS";
    title += " with (ts) purity corrections";
  }

  if ( GetGraph(name) )
  {
    // insure we're computing it only once
    return;
  }

  AliDebug(2,name);

  // insure we have all the graphs we need to work
  ComputeTriggerL0B(kMB);
  ComputeTriggerL0B(fReferenceTriggerType);

  const set<int>& runs = RunNumbers();

  vector<double> vx;
  vector<double> vxerr;
  vector<double> vy;
  vector<double> vyerr;

  Double_t purityREF(1.0);
  Double_t purityMB(1.0);
  Double_t purityREFerror(00);
  Double_t purityMBerror(0.0);

  // compute the per run values
  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    Int_t runNumber = *it;

    TString mbTrigger = GetTriggerClassName(kMB,runNumber);
//    TString refTrigger = GetTriggerClassName(fReferenceTriggerType,runNumber);

    purityMB=purityREF=1.0;
    purityMBerror=purityREFerror=0.0;

    if (eventSelectionCorrected>0)
    {
      ComputePurityFactorForScalerGraph(kMB,eventSelectionCorrected,&centrality);
      ComputePurityFactorForScalerGraph(fReferenceTriggerType,eventSelectionCorrected,&centrality);

      TGraphErrors* gpsfactorMB  = GetGraph(Form("PurityFactorForScaler%s%s_%s",GetEventSelectionName(eventSelectionCorrected).Data(), GetTriggerTypeName(kMB).Data(),centrality.Data()));
      TGraphErrors* gpsfactorREF = GetGraph(Form("PurityFactorForScaler%s%s_%s",GetEventSelectionName(eventSelectionCorrected).Data(), GetTriggerTypeName(fReferenceTriggerType).Data(),centrality.Data()));

      GetValueAndErrorFromGraph(gpsfactorMB,runNumber,purityMB,purityMBerror);
      GetValueAndErrorFromGraph(gpsfactorREF,runNumber,purityREF,purityREFerror);
    }

    if (purityMB<=0.0)
    {
      AliError(Form("Got purity=%e for MB for run %9d",purityMB,runNumber));
      continue;
    }

    TGraphErrors* gl0bMB = GetGraph(Form("L0B%s",GetTriggerTypeName(kMB).Data()));
    TGraphErrors* gl0bREF = GetGraph(Form("L0B%s",GetTriggerTypeName(fReferenceTriggerType).Data()));

    Double_t L0bMB,L0bMBError;
    Double_t L0bREF,L0bREFError;

    GetValueAndErrorFromGraph(gl0bMB,runNumber,L0bMB,L0bMBError);
    GetValueAndErrorFromGraph(gl0bREF,runNumber,L0bREF,L0bREFError);

    Double_t pileUpFactor(1.0);
    Double_t pileUpFactorError(0.0);

    if (pileUpCorrected)
    {
      ComputePileUpGraph(kMB,eventSelectionCorrected);

      TGraphErrors* gpu = GetGraph((Form("CorrectionPU%s%s",GetEventSelectionName(eventSelectionCorrected).Data(),GetTriggerTypeName(kMB).Data())));

      GetValueAndErrorFromGraph(gpu,runNumber,pileUpFactor,pileUpFactorError);
    }

    Double_t value;
    Double_t error;

    ScalerFnorm(value,error,
                L0bREF,purityREF,purityREFerror,
                L0bMB,purityMB,purityMBerror,
                pileUpFactor,pileUpFactorError);

    if ( value > 0.0 )
    {
      vx.push_back(1.0*runNumber);
      vxerr.push_back(0.5);
      vy.push_back(value);
      vyerr.push_back(error);
    }
  }

  CreateAndAddGraph(name,title,vx,vxerr,vy,vyerr);
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputeGraphRelDif(const char* a, const char* b) const
{
  // compute dispersion of b versus a
  //
  // computed differences graphs are put into the GRAPHS/ directory
  // computed differences results are put into the HISTOS/ directory

  TString name(Form("RelDif%svs%s",b,a));

  if ( GetGraph(name) )
  {
    // insure we compute it only once
    return;
  }

  AliDebug(2,name);

  TGraphErrors* ga = static_cast<TGraphErrors*>(MC()->GetObject(Form("/GRAPHS/%s",a)));
  TGraphErrors* gb = static_cast<TGraphErrors*>(MC()->GetObject(Form("/GRAPHS/%s",b)));

  if (!ga)
  {
    AliError(Form("Cannot get graph for %s",a));
    return;
  }

  if (!gb)
  {
    AliError(Form("Cannot get graph for %s",b));
    return;
  }

  if ( ga->GetN() != gb->GetN() )
  {
    AliError(Form("Cannot work with different number of points in the graphs : %d vs %d",
                  ga->GetN(),gb->GetN()));
    return;
  }

  TString title(Form("%s-%s (RelDif,%%)",b,a));

  vector<double> vx;
  vector<double> vxerr;
  vector<double> vy;
  vector<double> vyerr;

  for ( Int_t i = 0; i < ga->GetN(); ++i )
  {
    Double_t xa,xb,ya,yb;

    ga->GetPoint(i,xa,ya);
    gb->GetPoint(i,xb,yb);

    if ( xa != xb )
    {
      AliError(Form("Incompatible graphs : got xa=%e and xb=%e",xa,xb));
      return;
    }

    Double_t newvalue = 0.0;

    if ( TMath::Abs(xa) > 1E-12 )
    {
      newvalue = 100.0*( yb - ya ) / ya;
    }

    Double_t yerr = newvalue*AliAnalysisMuMuResult::ErrorAB(ya,ga->GetEY()[i],
                                                            yb,gb->GetEY()[i]);

    if ( fIsCompactGraphs )
    {
      xa = TString(ga->GetXaxis()->GetBinLabel(i+1)).Atoi()*1.0;
    }

    vx.push_back(xa);
    vxerr.push_back(0.5);
    vy.push_back(newvalue);
    vyerr.push_back(yerr);

  }

  CreateAndAddGraph(name,title,vx,vxerr,vy,vyerr);

  // FIXME : fill here an histogram from the graph to get the
  // weight of 1/e2 ?
  //  h->Fill(newvalue,1.0/(yerr*yerr));
  //  MC()->Adopt("/HISTOS/",h);

  //  AliAnalysisMuMuResult* r = GetRunIntegratedResult(*g,"FnormDispersion");
  //  if (r)
  //  {
  //    if (!dispersion)
  //    {
  //      dispersion = new AliAnalysisMuMuResult("FnormDispersion");
  //    }
  //    dispersion->AdoptSubResult(r);
  //    if ( !TString(g->GetName()).BeginsWith("Fnorm") )
  //    {
  //      dispersion->Exclude(r->Alias());
  //    }
  //  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputePileUpGraph(ETriggerType tt, Int_t eventSelectionCorrected)
{
  /// Compute the per-run graph of pile-up factor

  TString graphName(Form("CorrectionPU%s%s",GetEventSelectionName(eventSelectionCorrected).Data(),GetTriggerTypeName(tt).Data()));

  if ( GetGraph(graphName) )
  {
    // insure we're computing it only once
    return;
  }

  AliDebug(2,graphName);

  vector<double> vx;
  vector<double> vxerr;
  vector<double> vy;
  vector<double> vyerr;

  const set<int>& runs = RunNumbers();

  AliAnalysisTriggerScalers ts(runs,OCDBPath().Data());

  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    Int_t runNumber = *it;

    Double_t pileUpFactor(1.0);
    Double_t pileUpFactorError(0.0);
    Double_t purity(1.0);
    Double_t purityError(0.0);

    TString triggerClassName = GetTriggerClassName(tt,runNumber);

    if ( triggerClassName.Length()==0 )
    {
      AliError(Form("Unknown trigger type %d",tt));
      return;
    }

    if (eventSelectionCorrected)
    {
      GetPurity(triggerClassName.Data(),runNumber,purity,purityError,eventSelectionCorrected);
    }
    if(triggerClassName.Contains("CINT7-B-NOPF-MUFAST&0TVX")) ts.GetPileUpFactor(runNumber,"C0TVX-B-NOPF-CENTNOTRD",purity,pileUpFactor,pileUpFactorError); //specific to T0 in LHC15o.
    else ts.GetPileUpFactor(runNumber,triggerClassName.Data(),purity,pileUpFactor,pileUpFactorError);

    vx.push_back(runNumber);
    vxerr.push_back(0.5);
    vy.push_back(pileUpFactor);
    vyerr.push_back(pileUpFactorError);
  }

  TString title(Form("Pile-up correction factor for trigger %s",GetTriggerTypeName(tt).Data()));

  if (eventSelectionCorrected)
  {
    title += "( L0BRate corrected by event selection";
    title += GetEventSelectionName(eventSelectionCorrected);
    title += " accept fraction)";
  }

  CreateAndAddGraph(graphName,title,vx,vxerr,vy,vyerr);
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputeEventSelectionGraph(ETriggerType tt, Int_t eventSelectionCorrected)
{
  /// Compute the per-run graph of physics selection accept fraction
  /// for the given trigger

  TString graphName(Form("Correction%s%s",GetEventSelectionName(eventSelectionCorrected).Data(), GetTriggerTypeName(tt).Data()));

  if (GetGraph(graphName))
  {
    // insure we're computing it only once
    return;
  }

  AliDebug(2,graphName);

  vector<double> vx;
  vector<double> vxerr;
  vector<double> vy;
  vector<double> vyerr;

  const set<int>& runs = RunNumbers();

  // AliAnalysisTriggerScalers ts(runs,OCDBPath().Data()); // à quoi ça sert ??

  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    Int_t runNumber = *it;

    Double_t purity, purityError;

    TString triggerClassName = GetTriggerClassName(tt,runNumber);

    if ( triggerClassName.Length()==0 )
    {
      AliError(Form("Unknown trigger type %d",tt));
      return;
    }

    GetPurity(triggerClassName.Data(),runNumber,purity,purityError,eventSelectionCorrected);

    vx.push_back(runNumber);
    vxerr.push_back(0.5);
    vy.push_back(purity);
    vyerr.push_back(purityError);
  }

  TString title(Form("Fraction of events accepted by the event selection %s for trigger %s",GetTriggerTypeName(tt).Data(),
                     GetEventSelectionName(eventSelectionCorrected).Data()));

  CreateAndAddGraph(graphName,title,vx,vxerr,vy,vyerr);
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputePurityFactorForScalerGraph(ETriggerType tt, Int_t eventSelectionCorrected, TString* centrality)
{
  /// Compute the per-run graph of the scaler purity factor purity
  /// for the given trigger
  /// This factor is computed offline as FScal_purity = NTRIGG(PS,0-90%)/NTRIGG(All,All); Offline

  TString graphName(Form("PurityFactorForScaler%s%s_%s",GetEventSelectionName(eventSelectionCorrected).Data(), GetTriggerTypeName(tt).Data(),centrality->Data()));

  if (GetGraph(graphName))
  {
    // insure we're computing it only once
    return;
  }

  AliDebug(2,graphName);

  vector<double> vx;
  vector<double> vxerr;
  vector<double> vy;
  vector<double> vyerr;

  const set<int>& runs = RunNumbers();

  // AliAnalysisTriggerScalers ts(runs,OCDBPath().Data()); // à quoi ça sert ??

  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    Int_t runNumber = *it;

    Double_t purityfactor, purityfactorError;

    TString triggerClassName = GetTriggerClassName(tt,runNumber);

    if ( triggerClassName.Length()==0 )
    {
      AliError(Form("Unknown trigger type %d",tt));
      return;
    }

    GetPurityFactor(triggerClassName,runNumber,purityfactor,purityfactorError,eventSelectionCorrected,centrality);

    vx.push_back(runNumber);
    vxerr.push_back(0.5);
    vy.push_back(purityfactor);
    vyerr.push_back(purityfactorError);
  }

  TString title(Form("Fraction of %s events accepted for %s by all events for trigger %s",
    GetEventSelectionName(eventSelectionCorrected).Data(),centrality->Data(),GetTriggerTypeName(tt).Data()));

  CreateAndAddGraph(graphName,title,vx,vxerr,vy,vyerr);
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputeResultsFromGraphs()
{
  // Compute one single value for each graph, by weighting by the fraction
  // of events in each run
  // do this for certain groups of graphs

  TObjArray groups;
  groups.SetOwner(kTRUE);

  groups.Add(new TObjString("Fnorm"));
  groups.Add(new TObjString("NMBeq"));
  groups.Add(new TObjString("Correction"));
  groups.Add(new TObjString("RelDif"));

  TList* objectList = MC()->CreateListOfObjectNames("/GRAPHS/");

  TIter nextGroup(&groups);
  TObjString* grp;

  TIter next(objectList);
  TObjString* str(0x0);

  while ( ( grp = static_cast<TObjString*>(nextGroup()) ) )
  {
    TString oname(Form("/RESULTS/%s",grp->String().Data()));

    if ( MC()->GetObject(oname) )
    {
      // delete if we have it already so we can replace it
      MC()->Remove(oname);
    }

    AliAnalysisMuMuResult* result = new AliAnalysisMuMuResult(grp->String());

    MC()->Adopt("/RESULTS/",result);

    next.Reset();

    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      if ( ! ( str->String().BeginsWith(grp->String() ) ) ) continue;

      TGraphErrors* g = GetGraph(str->String());

      if (!g) continue;

      AliAnalysisMuMuResult* sub = GetRunIntegratedResult(*g);

      if ( !sub )
      {
        AliError(Form("Could not get result for graph %s",g->GetName()));
      }
      if ( sub )
      {
        result->AdoptSubResult(sub);
      }
    }
  }

  delete objectList;
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputeNofEvents(ETriggerType triggerType,
                                            Bool_t pileUpCorrected,
                                            Int_t eventSelectionCorrected)
{
  /// Compute trigger fractions

  TString graphName(Form("NofEvent%s%s%s",GetTriggerTypeName(triggerType).Data(),
                         pileUpCorrected ? "PU" : "",
                         GetEventSelectionName(eventSelectionCorrected).Data()));
  // Check if Compute has already been done
  if ( GetGraph(graphName) )
  {
    // compute it only once
    return;
  }

  ComputeCorrectionFactors(eventSelectionCorrected);

  TString gpsname(Form("Correction%s%s",GetEventSelectionName(eventSelectionCorrected).Data(),GetTriggerTypeName(triggerType).Data()));
  TString gpuname(Form("CorrectionPU%s%s",GetEventSelectionName(eventSelectionCorrected).Data(), GetTriggerTypeName(triggerType).Data()));

  TGraphErrors* gPS = GetGraph(gpsname);

  if (!gPS)
  {
    AliError(Form("Could not find %s",gpsname.Data()));
    return;
  }

  TGraphErrors* gPU = GetGraph(gpuname);

  if (!gPU)
  {
    AliError(Form("Could not find %s",gpuname.Data()));
    return;
  }

  AliDebug(2,graphName);

  vector<double> vx;
  vector<double> vxerr;
  vector<double> vy;
  vector<double> vyerr;

  const set<int>& runs = RunNumbers();

  Int_t i(0);

  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    Int_t runNumber = *it;

    TString triggerClassName = GetTriggerClassName(triggerType,runNumber);

    if ( triggerClassName.Length() )
    {
      Double_t n = GetSum(triggerClassName,runNumber,0);

      vx.push_back(runNumber);
      vxerr.push_back(0.5);

      assert(runNumber==TMath::Nint(gPU->GetX()[i]));

      Double_t y(n);
      Double_t y1(1.0);
      Double_t y2(1.0);
      Double_t e1(0);
      Double_t e2(0);

      if ( pileUpCorrected )
      {
        y1 = gPU->GetY()[i];
        e1 = gPU->GetEY()[i];
      }

      if ( eventSelectionCorrected > 0 )
      {
        y2 = gPS->GetY()[i];
        e2 = gPS->GetEY()[i];
      }

      y *= y1*y2;

      AliDebug(2,Form("RUN %09d n %e y1 %e y2 %e y% e",runNumber,n,y1,y2,y));

      Double_t yerr = y*AliAnalysisMuMuResult::ErrorABC( n, TMath::Sqrt(n),
                                                        y1, e1,
                                                        y2, e2);

      vy.push_back(y);
      vyerr.push_back(yerr);

      ++i;
    }
  }

  TString title(Form("Number of event of trigger %s",GetTriggerTypeName(triggerType).Data()));

  CreateAndAddGraph(graphName,title,vx,vxerr,vy,vyerr);
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputeTriggerFractions(ETriggerType triggerType,
                                                   Bool_t physicsSelectionCorrected)
{
  /// Compute trigger fractions

  TString graphName(Form("Fractions%s%s",GetTriggerTypeName(triggerType).Data(),physicsSelectionCorrected ? "PS" : ""));

  if ( GetGraph(graphName) )
  {
    // compute it only once
    return;
  }

  AliDebug(2,graphName);

  vector<double> vx;
  vector<double> vxerr;
  vector<double> vy;
  vector<double> vyerr;

  const set<int>& runs = RunNumbers();

  Double_t n(0.0);

  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    Int_t runNumber = *it;

    TString triggerClassName = GetTriggerClassName(triggerType,runNumber);

    if ( triggerClassName.Length() )
    {
      n += GetSum(triggerClassName,runNumber,physicsSelectionCorrected);
    }
  }

  if ( n <= 0.0 )
  {
    AliWarning(Form("Got zero trigger for %s",GetTriggerTypeName(triggerType).Data()));
    return;
  }

  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    Int_t runNumber = *it;

    TString triggerClassName = GetTriggerClassName(triggerType,runNumber);

    vx.push_back(runNumber);
    vxerr.push_back(0.5);
    Double_t y = GetSum(triggerClassName,runNumber,physicsSelectionCorrected);
    vy.push_back(y/n);
    vyerr.push_back( (y/n) * AliAnalysisMuMuResult::ErrorAB( y,TMath::Sqrt(y),
                                                        n, TMath::Sqrt(n)));
  }


  TString title(Form("Fraction of event of trigger %s",GetTriggerTypeName(triggerType).Data()));

  CreateAndAddGraph(graphName,title,vx,vxerr,vy,vyerr);
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ComputeTriggerL0B(ETriggerType triggerType)
{
  /// Compute trigger L0B

  vector<double> vx;
  vector<double> vxerr;
  vector<double> vy;
  vector<double> vyerr;

  TString graphName(Form("L0B%s",GetTriggerTypeName(triggerType).Data()));

  if ( GetGraph(graphName) )
  {
    // insure we're computing it only once
    return;
  }

  AliDebug(2,graphName);

  const set<int>& runs = RunNumbers();

  AliAnalysisTriggerScalers ts(runs,OCDBPath().Data());

  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    Int_t runNumber = *it;

    TString triggerClassName = GetTriggerClassName(triggerType,runNumber);

    AliAnalysisTriggerScalerItem* item = 0x0;

    if(triggerClassName.Contains("CINT7-B-NOPF-MUFAST&0TVX")) item=  ts.GetTriggerScaler(runNumber,"L0B","C0TVX-B-NOPF-CENTNOTRD"); //specific to T0 in LHC15o.
    else item= ts.GetTriggerScaler(runNumber,"L0B",triggerClassName);

    if (!item) continue;

    vx.push_back(runNumber);
    vxerr.push_back(0.5);

    Double_t y = item->Value();

    vy.push_back(y);

    vyerr.push_back( TMath::Sqrt(y) );
  }

  TString title(Form("L0B of trigger %s",GetTriggerTypeName(triggerType).Data()));

  CreateAndAddGraph(graphName,title,vx,vxerr,vy,vyerr);
}

////_____________________________________________________________________________
//void AliAnalysisMuMuFnorm::ComputeTSGraph(ETriggerType tt)
//{
//  /// Compute the per-run graph of physics selection accept fraction x track
//  /// accept fraction for the given trigger
//
//  TString graphName(Form("CorrectionTS%s",GetTriggerTypeName(tt).Data()));
//
//  if (GetGraph(graphName))
//  {
//    // insure we're computing it only once
//    return;
//  }
//
//  AliDebug(1,graphName);
//
//  vector<double> vx;
//  vector<double> vxerr;
//  vector<double> vy;
//  vector<double> vyerr;
//
//  const set<int>& runs = RunNumbers();
//
//  AliAnalysisTriggerScalers ts(runs,OCDBPath().Data());
//
//  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
//  {
//    Int_t runNumber = *it;
//
//    Double_t purity, purityError;
//
//    TString triggerClassName = GetTriggerClassName(tt,runNumber);
//
//    if ( triggerClassName.Length()==0 )
//    {
//      AliError(Form("Unknown trigger type %d",tt));
//      return;
//    }
//
//    GetPurity(triggerClassName.Data(),runNumber,purity,purityError,"OFFLINE1");
//
//    vx.push_back(runNumber);
//    vxerr.push_back(0.5);
//    vy.push_back(purity);
//    vyerr.push_back(purityError);
//  }
//
//  TString title(Form("Fraction of events accepted by the physics selection x track selection for trigger %s",GetTriggerTypeName(tt).Data()));
//
//  CreateAndAddGraph(graphName,title,vx,vxerr,vy,vyerr);
//}
//

//_____________________________________________________________________________
TGraphErrors* AliAnalysisMuMuFnorm::CreateAndAddGraph(const TString& name,
                                                      const TString& title,
                                                      const vector<double>& vx,
                                                      const vector<double>& vxerr,
                                                      const vector<double>& vy,
                                                      const vector<double>& vyerr) const
{
  /// Creates a graph and an histo and adds it to our mergeable collection

  TGraphErrors* g = new TGraphErrors(vx.size(),&vx[0],&vy[0],&vxerr[0],&vyerr[0]);
  g->SetName(name.Data());
  g->SetTitle(title.Data());
  if  (fIsCompactGraphs)
  {
    AliAnalysisMuMuGraphUtil::Compact(*g);
  }

//  g->GetXaxis()->SetTitle("Run number");

  TPaveText* text = new TPaveText(0.70,0.70,0.89,0.89,"NDC");
  text->SetBorderSize(0);
  text->SetFillColor(0);
  text->AddText(Form("Mean %e",g->GetMean(2)));
  text->AddText(Form("RMS  %e",g->GetRMS(2)));
  g->GetListOfFunctions()->Add(text);



  MC()->Adopt("/GRAPHS/",g);
  TH1* h = GetGraphAsHisto(name);
  MC()->Adopt("/GRAPHS/",h);
  return g;
}

//_____________________________________________________________________________
AliMergeableCollection* AliAnalysisMuMuFnorm::DetachMC()
{
  // let go the ownership of our mergeable collection
  fIsOwner = kFALSE;
  return fMergeableCollection;
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::DrawWith2Scales(const char* graphName1, const char* graphName2)
{
  TGraphErrors* g1 = static_cast<TGraphErrors*>(GetGraph(graphName1)->Clone());
  TGraphErrors* g2 = static_cast<TGraphErrors*>(GetGraph(graphName2)->Clone());

  if ( g1 && g2 )
  {
    gStyle->SetOptTitle(0);

    AliAnalysisMuMuGraphUtil gu;

    TCanvas* c = gu.DrawWith2Scales(*g1,*g2);
    c->Draw();

    for ( Int_t i = 0; i < 2; ++i )
    {
      c->cd(i);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.15);
    }

    c->Update();

    //    TGraphErrors* g = new TGraphErrors(g1->GetN());
    //
    //    Double_t check(0.0);
    //
    //    for ( Int_t i = 0; i < g->GetN(); ++i )
    //    {
    //      Double_t y = g1->GetY()[i]*g2->GetY()[i];
    //
    //      check += g2->GetY()[i];
    //
    //      g->SetPoint(i,g2->GetX()[i],y);
    //      g->SetPointError(i,g2->GetEX()[i],
    //                       y*AliAnalysisMuMuResult::ErrorAB(g1->GetY()[i],g1->GetEY()[i],
    //                                                        g2->GetY()[i],g2->GetEY()[i]));
    //    }
    //
    //    new TCanvas;
    //
    //    g->Draw("ap");
    //
    //    AliInfo(Form("check: %e g mean %e rms %e",check,g->GetMean(2),g->GetRMS(2)));

    /*

    // g1 vs g2

    TGraphErrors* g = new TGraphErrors(g1->GetN());

    for ( Int_t i = 0; i < g->GetN(); ++i )
    {
      g->SetPoint(i,g2->GetY()[i],g1->GetY()[i]);
      g->SetPointError(i,g2->GetEY()[i],g1->GetEY()[i]);
    }

    new TCanvas;

    g->Draw("ap");

    AliInfo(Form("g mean %e rms %e",g->GetMean(2),g->GetRMS(2)));

     */

  }

}

//_____________________________________________________________________________
TString AliAnalysisMuMuFnorm::GetEventSelectionName(Int_t eventSelectionCorrected) const
{
  if ( eventSelectionCorrected == 1 )
  {
    return "PS";
  }
  else if ( eventSelectionCorrected == 2 )
  {
    return "TS";
  }
  return "";
}

//_____________________________________________________________________________
TGraphErrors* AliAnalysisMuMuFnorm::GetGraph(const char* name) const
{
  /// shortcut method to give access to one graph

  TObject* o = MC()->GetObject(Form("/GRAPHS/%s",name));

  if (!o) return 0x0;

  if ( o->IsA() != TGraphErrors::Class() )
  {
    AliError(Form("Object %s is not of the expected type",o->GetName()));
    return 0x0;
  }

  return static_cast<TGraphErrors*>(o);
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMuFnorm::GetGraphAsHisto(const char* name) const
{
  /// shortcut method to give access to one graph and return it as an histo

  TObject* o = MC()->GetObject(Form("/GRAPHS/%s",name));

  if (!o) return 0x0;
  if ( o->IsA() != TGraphErrors::Class() )
  {
    AliError(Form("Object %s is not of the expected type",o->GetName()));
    return 0x0;
  }

  const set<int>& runs = RunNumbers();

  TGraphErrors *g = static_cast<TGraphErrors*>(o->Clone());
  TH1F * h =new TH1F(Form("%s_AsHisto",g->GetName()),Form("%s_AsHisto",g->GetName()),1,0.,1.);

  Double_t y = 0.;
  Double_t dy = 0.;

  //Fill
  for ( set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    Int_t runNumber = *it;

    GetValueAndErrorFromGraph(g,runNumber,y,dy);
    // if(y==0.)
    // {
    //   printf("Error : no values for run %d\n",runNumber );
    //   continue;
    // }

    Int_t bin = h->Fill(TString::Format("%d",runNumber).Data(),y);
    h->SetBinError(bin,dy);
  }

  // Set range
  int j=0;
  for (int i = 1; i < h->GetEntries()+1; i++)
  {
    if(TString(h->GetXaxis()->GetBinLabel(i)).IsNull()) continue;
    j++;
  }
  if(j == runs.size()) h->GetXaxis()->SetRange(1,j);

  delete g;
  return h;

}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::GetPurity(const char* triggerClassName, Int_t runNumber,
                                     Double_t& value, Double_t& error, Int_t eventSelectionCorrected) const
{
  /// Get the physics selection accept fraction for a given trigger

  value=error=0.0;

  TString runCondition;

  if (runNumber>0)
  {
    runCondition.Form("/run:%d",runNumber);
  }

  Double_t nall = fCounterCollection.GetSum(Form("trigger:%s/event:ALL%s",triggerClassName,runCondition.Data()));

  TString ename;

  if ( eventSelectionCorrected == 1 )
  {
    ename = fConfig.First(fConfig.EventSelectionKey(),kFALSE).Data();
  }
  else if ( eventSelectionCorrected == 2 )
  {
    ename = "OFFLINE1";
  }
  else
  {
    value = 1.0;
    return;

  }
  Double_t nps = fCounterCollection.GetSum(Form("trigger:%s/event:%s%s",
                                                triggerClassName,ename.Data(),runCondition.Data()));

  if ( nall <= 0.0 ) return;

  value = nps/nall;

  error = value*AliAnalysisMuMuResult::ErrorAB(nall,TMath::Sqrt(nall),nps,TMath::Sqrt(nps));
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::GetPurityFactor(TString triggerClassName, Int_t runNumber,
                                     Double_t& value, Double_t& error, Int_t eventSelectionCorrected, TString * centrality) const
{
  /// Get the physics selection accept fraction for a given trigger

  value=error=0.0;

  if (triggerClassName.Contains("CINT7-B-NOPF-MUFAST&0TVX")){
    value=1.;
    error=0.;
    return;
  }

  TString runCondition;

  if (runNumber>0)
  {
    runCondition.Form("/run:%d",runNumber);
  }

  TObjArray * centralityArray = centrality->Tokenize("_");
  TObjString * centralityType = static_cast<TObjString*>(centralityArray->At(0));

  TString ename;
  if ( eventSelectionCorrected == 1 ){
    if (triggerClassName.Contains("CMUL")) ename = "PSMUL";
    else ename = Form("%s",fConfig.First(fConfig.EventSelectionKey(),kFALSE).Data());
  }
  else if ( eventSelectionCorrected == 2 )ename = "OFFLINE1";
  else {
    value = 1.0;
    return;
  }

  Double_t b1 = fCounterCollection.GetSum(Form("trigger:%s/event:%s%s/centrality:%s",
                                                triggerClassName.Data(),ename.Data(),runCondition.Data(),centrality->Data()));
  Double_t b1sq = b1*b1;
  Double_t e1sq = b1;

  Double_t b2 = fCounterCollection.GetSum(Form("trigger:%s/event:ALL%s/centrality:%s",triggerClassName.Data(),runCondition.Data(),centralityType->String().Data()));
  Double_t b2sq = b2*b2;
  Double_t e2sq = b2;

  if ( b2 <= 0.0 ) return;

  value = b1/b2;

  //fully correlated bayasian
  if (b1 != b2) {
  error = TMath::Sqrt( TMath::Abs( ( (1. - 2.* b1 / b2) * e1sq  + b1sq * e2sq / b2sq ) / b2sq ) );
  }
  else error = 0;

}

//_____________________________________________________________________________
AliAnalysisMuMuResult* AliAnalysisMuMuFnorm::GetResult(const char* name) const
{
  /// shortcut method to give access to one result

  TObject* o = MC()->GetObject(Form("/RESULTS/%s",name));

  if (!o) return 0x0;

  if ( o->IsA() != AliAnalysisMuMuResult::Class() )
  {
    AliError(Form("Object %s is not of the expected type",o->GetName()));
    return 0x0;
  }

  return static_cast<AliAnalysisMuMuResult*>(o);
}

//______________________________________________________________________________
AliAnalysisMuMuResult* AliAnalysisMuMuFnorm::GetRunIntegratedResult(const TGraphErrors& g, const char* basename)
{
  /// get one value +- error for this graph (weighting the run-by-run values
  /// by the fraction of the number of triggers in that run)
  /// also the rms is computed.

  Bool_t physicsSelectionCorrected(kFALSE);

  if ( TString(g.GetName()).Contains("PS") )
  {
    physicsSelectionCorrected=kTRUE;
  }

  ComputeTriggerFractions(fReferenceTriggerType,physicsSelectionCorrected);

  TString fname(Form("Fractions%s%s",GetTriggerTypeName(fReferenceTriggerType).Data(),physicsSelectionCorrected ? "PS" : ""));

  TGraphErrors* gTriggerFractions = GetGraph(fname);

  StdoutToAliDebug(2,cout << g.GetName() << endl; g.Print(););

  if (!gTriggerFractions)
  {
    AliError(Form("Did not find graph for %s",fname.Data()));
    return 0x0;
  }

  Double_t mean = g.GetY()[0];
  Int_t n = g.GetN();

  Double_t errorOnMean = g.GetEY()[0];
  Double_t rms = 0.0;

  if ( n > 1 )
  {
    mean = errorOnMean = 0.0;
    Double_t d(0.0);
    Double_t v1(0.0);
    Double_t v2(0.0);

    for ( Int_t i = 0; i < n; ++i )
    {
      Double_t y = g.GetY()[i];

      Double_t weight = gTriggerFractions->GetY()[i]; // sum of weight should be 1.0

      AliDebug(2,Form("%s i %3d y %e weight %e",g.GetName(),i,y,weight));

      mean += y * weight;

      v1 += weight;
      v2 += weight*weight;
    }

    mean /= v1;

    for ( Int_t i = 0; i < n; ++i )
    {
      Double_t weight = gTriggerFractions->GetY()[i]; // sum of weight should be 1.0
      Double_t y = g.GetY()[i];

      d += (y-mean)*(y-mean)*weight;
    }

    AliDebug(2,Form("v1=%e v2=%e d=%e",v1,v2,d));

    if ( v1 <= 0) return 0x0;

    errorOnMean = TMath::Sqrt((1.0/v1)*(1.0/(n-1))*d);

    rms = TMath::Sqrt( (v1/(v1*v1-v2))*d );
  }

  AliAnalysisMuMuResult* result = new AliAnalysisMuMuResult(g.GetName(),g.GetTitle());

  result->Set(basename,mean,errorOnMean,rms);

  if ( TString(g.GetName()) == "FnormScalersPUPS" )
  {
    AliDebug(1,Form("mean %e errorOnMean %e rms %e",mean,errorOnMean,rms));
    StdoutToAliDebug(1,result->Print("full"));
  }

  return result;
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuFnorm::GetSum(const char* triggerClassName,
                                      Int_t runNumber,
                                      Int_t eventSelectionCorrected) const
{
  TString condition(Form("trigger:%s/run:%d",triggerClassName,runNumber));
  TString triggerName(triggerClassName);

  if (eventSelectionCorrected==1)
  {
    if(triggerName.Contains("CMSL") )condition += "/event:PSMSL";
    else condition += Form("/event:%s",fConfig.First(fConfig.EventSelectionKey(),kFALSE).Data());
  }
  else if ( eventSelectionCorrected == 2 )
  {
    condition += "/event:OFFLINE1";
  }
  else
  {
    condition += "/event:ALL";
  }

  condition += Form("/centrality:%s",fConfig.First(fConfig.CentralitySelectionKey(),kFALSE).Data());

  Double_t n = fCounterCollection.GetSum(condition.Data());

  printf("Sum of %s for %s = %f \n",triggerClassName,condition.Data(),n );

  if (n<=0)
  {
    AliError(Form("Got no count for %s for run %d (physicsSelected:%d)",triggerClassName,runNumber,eventSelectionCorrected));
    return 0;
  }

  return n;
}

//_____________________________________________________________________________
TString AliAnalysisMuMuFnorm::GetTriggerClassName(ETriggerType tt, Int_t runNumber) const
{
  // get the triggerclass to for a given trigger type and run number

  if ( tt == kMB )
  {
    return MBTriggerClassName(runNumber);
  }
  else if ( tt == kMUL )
  {
    return MULTriggerClassName(runNumber);
  }
  else if ( tt == kMSL)
  {
    return MSLTriggerClassName(runNumber);
  }
  else if ( tt == kMSH)
  {
    return MSHTriggerClassName(runNumber);
  }
  return "";
}

//_____________________________________________________________________________
TString AliAnalysisMuMuFnorm::GetTriggerTypeName(ETriggerType tt) const
{
  // get the name of the trigger type
  if ( tt == kMB )
  {
    return "MB";
  }
  else if ( tt == kMUL )
  {
    return "MUL";
  }
  else if ( tt == kMSL)
  {
    return "MSL";
  }
  else if ( tt == kMSH)
  {
    return "MSH";
  }
  return "";
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::GetValueAndErrorFromGraph(TGraphErrors* graph,
                                                     Int_t runNumber,
                                                     Double_t& value,
                                                     Double_t& error) const
{
  /// get (value,error) corresponding to run=runNumber.
  /// Works both for compact and non-compact graphs

  value = TMath::Limits<Double_t>::Max();
  error = 0.0;

  if (!graph) return;

  TAxis* axis = graph->GetXaxis();

  for ( Int_t i = 0; i < graph->GetN(); ++i )
  {
    Int_t rbin = TMath::Nint(graph->GetX()[i]);
    Int_t rlabel = TString(axis->GetBinLabel(i+1)).Atoi();
    if ( rbin == runNumber || rlabel == runNumber )
    {
      value = graph->GetY()[i];
      error = graph->GetEY()[i];
    }
  }
}

//_____________________________________________________________________________
AliMergeableCollection* AliAnalysisMuMuFnorm::MC() const
{
  // get our mergeable collection
  if (!fMergeableCollection)
  {
    fMergeableCollection = new AliMergeableCollection("Fnorm",Form("MB to %s trigger normalization results",GetTriggerTypeName(fReferenceTriggerType).Data()));
  }
  return fMergeableCollection;
}

//_____________________________________________________________________________
TString AliAnalysisMuMuFnorm::MBTriggerClassName(Int_t runNumber) const
{
  /// FIXME : find a better way ?

  // if ( TriggerClassnameTest("CINT7-B-NOPF-ALLNOTRD",runNumber) )
  // {
  //   return "CINT7-B-NOPF-ALLNOTRD";
  // }
  // else if ( TriggerClassnameTest("CPBI2_B1-B-NOPF-ALLNOTRD",runNumber) )
  // {
  //   return "CPBI2_B1-B-NOPF-ALLNOTRD";
  // }
  // return "";

  TString triggerType(fConfig.First(fConfig.MinbiasTriggerKey(),kFALSE).Data());
  return triggerType.Data();
}


//_____________________________________________________________________________
TString AliAnalysisMuMuFnorm::MSHTriggerClassName(Int_t runNumber) const
{
  /// FIXME : find a better way ?

  if ( TriggerClassnameTest("CMSH7-B-NOPF-ALLNOTRD",runNumber) )
  {
    return "CMSH7-B-NOPF-ALLNOTRD";
  }
  else if ( TriggerClassnameTest("CMSH7-B-NOPF-MUON",runNumber) )
  {
    return "CMSH7-B-NOPF-MUON";
  }
  return "";
}

//_____________________________________________________________________________
TString AliAnalysisMuMuFnorm::MSLTriggerClassName(Int_t runNumber) const
{
  /// Take the first key in the config. file.

 TString triggerType(fConfig.First(fConfig.MuonTriggerKey(),kFALSE).Data());
  return triggerType.Data();
}


//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::MultiplyGraphs(const char* g1name, const char* g2name, const char* name)
{
  /// Make a new graph = g1*g2
  vector<double> vx;
  vector<double> vy;
  vector<double> vxerr;
  vector<double> vyerr;

  TGraphErrors* g1 = GetGraph(g1name);
  TGraphErrors* g2 = GetGraph(g2name);

  if (!g1)
  {
    AliError(Form("Could not get graph %s",g1name));
    return;
  }

  if (!g2)
  {
    AliError(Form("Could not get graph %s",g2name));
    return;
  }

  if ( g1->GetN() != g2->GetN() )
  {
    AliError(Form("Could not multiply incompatible graphs %d pts vs %d pts",g1->GetN(),g2->GetN()));
    return;
  }

  for ( Int_t i = 0; i < g1->GetN(); ++i )
  {
    if ( g1->GetX()[i] != g2->GetX()[i] )
    {
      AliError(Form("Incompatible bin %d : %e vs %e",i,g1->GetX()[i],g2->GetX()[i]));
      return;
    }

    vx.push_back(g1->GetX()[i]);
    vxerr.push_back(g1->GetEX()[i]);

    Double_t y = g1->GetY()[i]*g2->GetY()[i];
    Double_t yerr = y*AliAnalysisMuMuResult::ErrorAB( g1->GetY()[i], g1->GetEY()[i],
                                                     g2->GetY()[i], g2->GetEY()[i]);

    vy.push_back(y);
    vyerr.push_back(yerr);
  }

  TString gname(name);

  if ( gname.Length() == 0 )
  {
    gname = g1->GetName();
    gname += "x";
    gname += g2->GetName();
  }

  TString title(Form("Product of %s by %s",g1->GetName(),g2->GetName()));

  CreateAndAddGraph(gname,title,vx,vxerr,vy,vyerr);
}

//_____________________________________________________________________________
TString AliAnalysisMuMuFnorm::MULTriggerClassName(Int_t runNumber) const
{
  /// Take the first key in the config. file.

  TString triggerType(fConfig.First(fConfig.DimuonTriggerKey(),kFALSE).Data());
    return triggerType.Data();


}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::Print(Option_t* opt) const
{
  if ( fMergeableCollection )
  {
    fMergeableCollection->Print(opt);
  }
}

//_____________________________________________________________________________
set<int> AliAnalysisMuMuFnorm::RunNumbers() const
{
  /// Extract the run numbers from our counter collection

  set<int> runset;

  TString sruns = fCounterCollection.GetKeyWords("run");
  TObjArray* runs = sruns.Tokenize(",");

  TIter next(runs);
  TObjString* s;

  while ( ( s = static_cast<TObjString*>(next())) )
  {
    runset.insert(s->String().Atoi());
  }
  delete runs;

  return runset;
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ScalerFnorm(Double_t& value, Double_t& error,
                                       Double_t L0bREF, Double_t purityREF, Double_t purityREFerror,
                                       Double_t L0bMB, Double_t purityMB, Double_t purityMBerror,
                                       Double_t pileUpFactor, Double_t pileUpFactorError)
{
  /// Compute the MB to REF ratio and its associated error

  value = error = 0.0;

  value = L0bREF*purityREF;

  if (!value) return;

  value = L0bMB*purityMB*pileUpFactor/value;

  error = value*AliAnalysisMuMuResult::ErrorABCDE(L0bREF,TMath::Sqrt(L0bREF),
                                                  purityREF,purityREFerror,
                                                  L0bMB,TMath::Sqrt(L0bMB),
                                                  purityMB,purityMBerror,
                                                  pileUpFactor,pileUpFactorError);
}

//_____________________________________________________________________________
void AliAnalysisMuMuFnorm::ShowFnorm(const TObjArray& a) const
{
  /// Print and plot Fnorm values
  TIter next(&a);
  TGraphErrors* g;

  while ( ( g = static_cast<TGraphErrors*>(next()) ) )
  {
    g->SetTitle(Form("Fnorm %s",g->GetTitle()));

    new TCanvas(Form("c%s",g->GetName()));

    if (fIsCompactGraphs)
    {
      AliAnalysisMuMuGraphUtil::Compact(*g);
    }
    else
    {
      g->GetXaxis()->SetNoExponent();
    }
    g->Draw("lpa");

    Double_t y(0.0);
    Double_t yerr(0.0);

    for ( Int_t i = 0; i < g->GetN(); ++i )
    {
      y += g->GetY()[i];
      Double_t e = ( y != 0.0  ? g->GetEY()[i]/y : 0.0);

      yerr += e*e;
    }

    y /= (g->GetN());
    yerr = TMath::Sqrt(yerr)*y;

    AliInfo(Form("%30s graph %e +- %e (%5.2f %%) RMS %e",g->GetName(),
                 y,yerr,yerr*100/y,g->GetRMS(2)));

  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFnorm::TriggerClassnameTest(const char* triggerClassName, Int_t runNumber) const
{
  /// Check if we have counts for that trigger,run combination

  TString runs = fCounterCollection.GetKeyWords("run");

  if ( !runs.Contains(Form("%d",runNumber)) ) return kFALSE;

  TString triggers = fCounterCollection.GetKeyWords("trigger");

  if (!triggers.Contains(triggerClassName)) return kFALSE;

  Double_t n = fCounterCollection.GetSum(Form("trigger:%s/run:%d",triggerClassName,runNumber));

  return ( n > 0.0 );
}

//_____________________________________________________________________________
void
AliAnalysisMuMuFnorm::WeightedMeanGraphs(const char* patternOrList, const char* graphName, AliMergeableCollection* oc)
{
  /// Sum the graphs which name matches pattern
  ///
  /// Sum is made using a weighted mean (each element is weighted by the inverse
  /// of its error squared)

  TString spattern(patternOrList);
  TObjArray* slist(0x0);

  if ( spattern.CountChar(',') )
  {
    // it's not a pattern but a list...
    slist = spattern.Tokenize(",");
    spattern = "";
  }

  TList* objectList = 0x0;
  if(oc) {
    printf("Adding %s to the list...\n", "/FNORM/Offline/GRAPHS/");
    objectList = oc->CreateListOfObjectNames("/FNORM/Offline/GRAPHS/");
    printf("Adding %s to the list...\n", "/FNORM/Scaler/GRAPHS/");
    objectList->Add(static_cast<TObject*>(oc->CreateListOfObjectNames("/FNORM/Scaler/GRAPHS/")));
  }
  else return;

  if (!objectList){
    printf("Cannot add list\n");
    return;
  }

  TIter next(objectList);
  TObjString* str(0x0);
  TObjArray selected;
  selected.SetOwner(kFALSE);

  printf("Selecting graphs ...\n");
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TGraphErrors* g = static_cast<TGraphErrors*>(oc->GetObject(Form("/FNORM/Offline/GRAPHS/%s",str->String().Data())));
    if (!g) g = static_cast<TGraphErrors*>(oc->GetObject(Form("/FNORM/Scaler/GRAPHS/%s",str->String().Data())));
    if (!g) continue;

    TString name(g->GetName());

    if ( spattern.Length() >0 && !name.Contains(spattern.Data()) ) continue;
    AliDebug(2,Form("name : %s !\n",name.Data() ));

    if ( slist && !slist->FindObject(name)) continue;

    AliDebug(2,Form("Selected for sum : %s",name.Data()));

    selected.Add(g->Clone());
  }

  delete slist;
  delete objectList;

  if (selected.GetLast()<0) return;

  vector<double> vx;
  vector<double> vy;
  vector<double> vxerr;
  vector<double> vyerr;

  Int_t npts = static_cast<TGraphErrors*>(selected.First())->GetN();
  if(npts==0 || npts==1) return;

  printf("Computing mean ...\n");

  for ( Int_t ipoint = 0; ipoint < npts; ++ipoint )
  {
    Double_t x(0.0),xref(0.0),xerr(0.0);
    Double_t sum(0.0);
    Double_t sume2(0.0);

    // Loop on graphs
    for ( Int_t igraph = 0; igraph <= selected.GetLast(); ++igraph )
    {
      TGraphErrors* g = static_cast<TGraphErrors*>(selected.At(igraph));

      if ( g->GetN() != npts )
      {
        AliError(Form("Graph %s does not have the expected %d points",g->GetName(),npts));
        continue;
      }
      Double_t runNumber;

      if ( fIsCompactGraphs )
      {
        runNumber = TString(g->GetXaxis()->GetBinLabel(ipoint+1)).Atoi()*1.0;
      }
      else
      {
        runNumber = g->GetX()[ipoint];
      }
      if ( igraph == 0 )
      {
        xref = g->GetX()[ipoint];
        x = runNumber;
        xerr = g->GetEX()[ipoint];
      }
      else
      {
        if ( xref != g->GetX()[ipoint] )
        {
          AliError(Form("Cannot sum graphs with different axis : get %e and expected %e : %s vs %s",xref,x,selected.At(0)->GetName(),g->GetName()));
          return;
        }
      }

      Double_t e2 = g->GetEY()[ipoint]*g->GetEY()[ipoint];

      if ( e2 != 0.0 )
      {
        sum += g->GetY()[ipoint]/e2;
        sume2 += 1.0/e2;
      }
    }

    if (sume2 != 0.0)
    {
      vx.push_back(x);
      vxerr.push_back(xerr);
      vy.push_back(sum/sume2);
      vyerr.push_back(TMath::Sqrt(1/sume2));
    }
  }

  printf("Changing titles ...\n");
  Int_t n = selected.GetEntries();

  TString name(graphName);
  TString title(Form("Weighted mean from %d individual graphs",n));

  if ( strlen(graphName) == 0)
  {
    name = TString::Format("WeightMeanFnorm%s",patternOrList);
    name.ReplaceAll(",","_");
    title = TString::Format("WeightMeanFnorm%s from %d individual graphs",patternOrList,n);
  }

  printf("Creating graph ...\n");
  CreateAndAddGraph(name,title,vx,vxerr,vy,vyerr);

}


