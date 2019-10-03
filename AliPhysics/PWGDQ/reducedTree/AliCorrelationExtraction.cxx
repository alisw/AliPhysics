/**************************************************************
 *                                                            *
 * class used for the extraction of J/psi-hadron correlations *
 *                                                            *
 * authors: Lucas Altenkamper (lucas.altenkamper@cern.ch)     *
 *          Antoine Lardeux   (antoine.lardeux@cern.ch)       *
 *                                                            *
 * 12/08/2018                                                 *
 *                                                            *
 **************************************************************/

#ifndef ALICORRELATIONEXTRACTION_H
#include "AliCorrelationExtraction.h"
#endif

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <memory>

#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THn.h>
#include <THnSparse.h>
#include <TMath.h>
#include <TRandom.h>
#include <TMinuit.h>

#include "AliReducedVarManager.h"

using std::cout;
using std::endl;

ClassImp(AliCorrelationExtraction)

Double_t EPSILON = 1.0e-6;

//_______________________________________________________________________________
AliCorrelationExtraction::AliCorrelationExtraction() :
  fSEOS(0x0),
  fSEOSSparse(0x0),
  fSEPP(0x0),
  fSEPPSparse(0x0),
  fSEMM(0x0),
  fSEMMSparse(0x0),
  fMEOS(0x0),
  fMEOSSparse(0x0),
  fMEOS2(0x0),
  fMEOS2Sparse(0x0),
  fMEPP(0x0),
  fMEPPSparse(0x0),
  fMEMM(0x0),
  fMEMMSparse(0x0),
  fSEPPPair(0x0),
  fSEMMPair(0x0),
  fMEOSPair(0x0),
  fMEPPPair(0x0),
  fMEMMPair(0x0),
  fHadronEff(0x0),
  fSEOSNorm(0x0),
  fSEOSNormBackgroundMassWindow(),
  fSEPPNormBackgroundMassWindow(),
  fSEMMNormBackgroundMassWindow(),
  fSEPPPairInvMass(0x0),
  fSEMMPairInvMass(0x0),
  fMEOSPairInvMass(0x0),
  fMEPPPairInvMass(0x0),
  fMEMMPairInvMass(0x0),
  fMEOSNorm(0x0),
  fMEOSNormBackgroundMassWindow(),
  fMEOS2NormBackgroundMassWindow(),
  fMEPPNormBackgroundMassWindow(),
  fMEMMNormBackgroundMassWindow(),
  fInclusiveCF1D(0x0),
  fInclusiveCF1DInvMass(),
  fInclusiveCF1DInvMassBackground(),
  fInclusiveCF1DBackgroundMassWindow(),
  fInclusiveCF2D(0x0),
  fInclusiveCF2DBackgroundMassWindow(),
  fInclusiveCF3D(0x0),
  fBackgroundCF1D(0x0),
  fBackgroundCF1DInvMassFit(),
  fBackgroundCF1DInvMass(),
  fBackgroundCF2D(0x0),
  fCombinatorialBackgroundCF1D(),
  fCombinatorialBackgroundCF2D(),
  fSignalCF1D(0x0),
  fSignalCF1DEffCorr(0x0),
  fSignalCF1DInvMass(),
  fSignalCF2D(0x0),
  fSignalCF2DEffCorr(0x0),
  fNVariables(0),
  fVariables(),
  fVarLimits(),
  fVarIndices(),
  fNMixingVariables(0),
  fMixingVariables(),
  fMixingVarIndices(),
  fNMixingVarBins(),
  fMixingVarBinLimits(),
  fMassVariable(AliReducedVarManager::kMass),
  fMassVariableIndex(-1),
  fMassVariableIndexPair(-1),
  fDeltaPhiVariable(AliReducedVarManager::kDeltaPhi),
  fDeltaPhiVariableIndex(-1),
  fDeltaEtaVariable(AliReducedVarManager::kDeltaEta),
  fDeltaEtaVariableIndex(-1),
  fEfficiencyVariable(AliReducedVarManager::kAssociatedPt),
  fEfficiencyVariableIndex(-1),
  fVerboseFlag(kFALSE),
  fUseMixingVars(kFALSE),
  fIntegrateDeltaEta(),
  fResonanceFits(0x0),
  fOptionBkgMethod(kBkgNone),
  fBkgFitFunction(0x0),
  fNBackgroundMassRanges(0),
  fBackgroundMassRanges(),
  fMassSignalRange(),
  fTrigValSig(),
  fTrigValBkg(),
  fProcessDone(kFALSE)
{
  //
  // default constructor
  //
  for (Int_t i=0; i<kNBackgroundMethods; ++i) fIntegrateDeltaEta[i] = kFALSE;
  fMassSignalRange[0] = -1;
  fMassSignalRange[1] = -1;
  for(Int_t i=0; i<kNMaxVariables; ++i) {
    fVariables[i]       = -1;
    fVarLimits[i][0]    = 0;
    fVarLimits[i][1]    = 0;
    fVarIndices[i]      = -1;
  }
  for (Int_t i=0; i<kNMaxMixingVariables; ++i) {
    fMixingVariables[i]     = -1;
    fMixingVarIndices[i]    = -1;
    fNMixingVarBins[i]      = 0;
    for (Int_t j=0; j<kNMaxMixingVarBins; ++j) {
      fMixingVarBinLimits[i][j][0] = 0;
      fMixingVarBinLimits[i][j][1] = 0;
    }
  }
  
  for (Int_t i=0; i<kNMaxBackgroundMassRanges; ++i) {
    fSEOSNormBackgroundMassWindow[i]      = 0x0;
    fSEPPNormBackgroundMassWindow[i]      = 0x0;
    fSEMMNormBackgroundMassWindow[i]      = 0x0;
    fMEOSNormBackgroundMassWindow[i]      = 0x0;
    fMEOS2NormBackgroundMassWindow[i]     = 0x0;
    fMEPPNormBackgroundMassWindow[i]      = 0x0;
    fMEMMNormBackgroundMassWindow[i]      = 0x0;
    fInclusiveCF1DBackgroundMassWindow[i] = 0x0;
    fInclusiveCF2DBackgroundMassWindow[i] = 0x0;
    fCombinatorialBackgroundCF1D[i]       = 0x0;
    fCombinatorialBackgroundCF2D[i]       = 0x0;
    fBackgroundMassRanges[i][0]           = -1;
    fBackgroundMassRanges[i][1]           = -1;
    for (Int_t j=0; j<kNTriggerValues; ++j)
      fTrigValBkg[i][j]                   = -9999.;
  }
  for (Int_t i=0; i<kNMaxDeltaPhiBins; ++i) {
    for (Int_t j=0; j<kNMaxDeltaEtaBins; ++j) {
      fInclusiveCF1DInvMass[i][j]           = 0x0;
      fInclusiveCF1DInvMassBackground[i][j] = 0x0;
      fBackgroundCF1DInvMassFit[i][j]       = 0x0;
      fBackgroundCF1DInvMass[i][j]          = 0x0;
      fSignalCF1DInvMass[i][j]              = 0x0;
    }
  }
  for (Int_t i=0; i<kNTriggerValues; ++i) fTrigValSig[i] = -9999.;
}

//_______________________________________________________________________________
AliCorrelationExtraction::~AliCorrelationExtraction()
{
  //
  // default destructor
  //
  if (fSEOS)            delete fSEOS;
  if (fSEOSSparse)      delete fSEOSSparse;
  if (fSEPP)            delete fSEPP;
  if (fSEPPSparse)      delete fSEPPSparse;
  if (fSEMM)            delete fSEMM;
  if (fSEMMSparse)      delete fSEMMSparse;
  if (fMEOS)            delete fMEOS;
  if (fMEOSSparse)      delete fMEOSSparse;
  if (fMEOS2)           delete fMEOS2;
  if (fMEOS2Sparse)     delete fMEOS2Sparse;
  if (fMEPP)            delete fMEPP;
  if (fMEPPSparse)      delete fMEPPSparse;
  if (fMEMM)            delete fMEMM;
  if (fMEMMSparse)      delete fMEMMSparse;
  if (fSEPPPair)        delete fSEPPPair;
  if (fSEMMPair)        delete fSEMMPair;
  if (fMEOSPair)        delete fMEOSPair;
  if (fMEPPPair)        delete fMEPPPair;
  if (fMEMMPair)        delete fMEMMPair;
  if (fHadronEff)       delete fHadronEff;
  if (fResonanceFits)   delete fResonanceFits;
  if (fBkgFitFunction)  delete fBkgFitFunction;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::SetHistograms(THnF* seos, THnF* sepp/*=0x0*/, THnF* semm/*=0x0*/,
                                             THnF* meos/*=0x0*/, THnF* mepp/*=0x0*/, THnF* memm/*=0x0*/) {
  //
  // set multi-dimensional histograms
  //
  fSEOS = seos;
  fSEPP = sepp;
  fSEMM = semm;
  fMEOS = meos;
  fMEPP = mepp;
  fMEMM = memm;
  fProcessDone = kFALSE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::SetHistograms(THnSparseF* seos, THnSparseF* sepp/*=0x0*/, THnSparseF* semm/*=0x0*/,
                                             THnSparseF* meos/*=0x0*/, THnSparseF* mepp/*=0x0*/, THnSparseF* memm/*=0x0*/) {
  //
  // set multi-dimensional histograms
  //
  fSEOSSparse = seos;
  fSEPPSparse = sepp;
  fSEMMSparse = semm;
  fMEOSSparse = meos;
  fMEPPSparse = mepp;
  fMEMMSparse = memm;
  fProcessDone = kFALSE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::SetPairHistograms(THnF* sepp/*=0x0*/, THnF* semm/*=0x0*/,
                                                 THnF* meos/*=0x0*/, THnF* mepp/*=0x0*/, THnF* memm/*=0x0*/) {
  //
  // set multi-dimensional histograms
  //
  fSEPPPair = sepp;
  fSEMMPair = semm;
  fMEOSPair = meos;
  fMEPPPair = mepp;
  fMEMMPair = memm;
  fProcessDone = kFALSE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::SetBackgroundMethod(Int_t method, Bool_t integrateDeltaEta/*=kFALSE*/) {
  fOptionBkgMethod = method;
  if (fOptionBkgMethod>-1 && integrateDeltaEta) fIntegrateDeltaEta[fOptionBkgMethod] = kTRUE;
  fProcessDone = kFALSE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::SetBackgroundMassWindows(Int_t n, Double_t* min, Double_t* max) {
  //
  // set mass windows for background determination
  //
  fNBackgroundMassRanges = n;
  for (Int_t i=0; i<n; ++i) {
    if (fNBackgroundMassRanges>=kNMaxBackgroundMassRanges) {
      cout << "AliCorrelationExtraction::SetBackgroundMassWindows() Fatal: Number of requested mass background windows " << n << " larger than allowed number " << kNMaxBackgroundMassRanges << "!" << endl;
      return;
    }
    fBackgroundMassRanges[i][0] = min[i];
    fBackgroundMassRanges[i][1] = max[i];
  }
  fProcessDone = kFALSE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::AddVariables(Int_t nVars, Int_t* vars, Int_t* indices) {
  //
  // initialize variable types and mapping of the dimensions in the THnF
  // vars     = variable enum pos. from AliReducedVarManager
  // indices  = dimensions in THnF
  //
  for(Int_t i=0; i<nVars; ++i) {
    if (fNVariables>=kNMaxVariables) {
      cout << "AliCorrelationExtraction::AddVariables() Fatal: Number of requested variables " << fNVariables << " larger than allowed number " << kNMaxVariables << "!" << endl;
      return;
    }
    fVariables[fNVariables]   = vars[i];
    fVarIndices[fNVariables]  = indices[i];
    if (fSEOS) {
      fVarLimits[fNVariables][0] = fSEOS->GetAxis(indices[i])->GetXmin()+EPSILON;
      fVarLimits[fNVariables][1] = fSEOS->GetAxis(indices[i])->GetXmax()-EPSILON;
    } else if (fSEOSSparse) {
      fVarLimits[fNVariables][0] = fSEOSSparse->GetAxis(indices[i])->GetXmin()+EPSILON;
      fVarLimits[fNVariables][1] = fSEOSSparse->GetAxis(indices[i])->GetXmax()-EPSILON;
    }
    fNVariables++;
  }
  fProcessDone = kFALSE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::AddVariable(Int_t var, Int_t index) {
  //
  // initialize a variable being used in the THnF
  // var    = variable enum pos. from AliReducedVarManager
  // index  = dimension in THnF
  //
  if (fNVariables>=kNMaxVariables) {
    cout << "AliCorrelationExtraction::AddVariable() Fatal: Number of requested variables " << fNVariables << " larger than allowed number " << kNMaxVariables << "!" << endl;
    return;
  }
  fVariables[fNVariables]   = var;
  fVarIndices[fNVariables]  = index;
  if(fSEOS) {
    fVarLimits[fNVariables][0] = fSEOS->GetAxis(index)->GetXmin()+EPSILON;
    fVarLimits[fNVariables][1] = fSEOS->GetAxis(index)->GetXmax()-EPSILON;
  } else if (fSEOSSparse) {
    fVarLimits[fNVariables][0] = fSEOSSparse->GetAxis(index)->GetXmin()+EPSILON;
    fVarLimits[fNVariables][1] = fSEOSSparse->GetAxis(index)->GetXmax()-EPSILON;
  }
  fNVariables++;
  fProcessDone = kFALSE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::SetVarRange(Int_t var, Double_t* lims) {
  //
  // set the user range for variable var
  // NOTE: Variable var is encoded using the AliReducedVarManager::Variables enum
  // NOTE: If the var is not found in fVariables, nothing happens
  // NOTE: The user provided limits are slightly modified to avoid bin edge problems
  //
  Int_t idx = -1;
  for (Int_t i=0; i<fNVariables; ++i) {
    if (fVariables[i]==var) idx = i;
  }
  if (idx==-1) return;
  fVarLimits[idx][0] = lims[0]+EPSILON;
  fVarLimits[idx][1] = lims[1]-EPSILON;
  fProcessDone = kFALSE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::SetVarRange(Int_t var, Double_t min, Double_t max) {
  //
  // set the user range for variable var
  // NOTE: Variable var is encoded using the AliReducedVarManager::Variables enum
  // NOTE: If the var is not found in fVariables, nothing happens
  // NOTE: The user provided limits are slightly modified to avoid bin edge problems
  //
  Int_t idx = -1;
  for (Int_t i=0; i<fNVariables; ++i) {
    if (fVariables[i]==var) idx = i;
  }
  if (idx==-1) return;
  fVarLimits[idx][0] = min+EPSILON;
  fVarLimits[idx][1] = max-EPSILON;
  fProcessDone = kFALSE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::AddMixingVariables(Int_t nVars, Int_t* vars, Int_t* indices) {
  //
  // initialize variable types and mapping of the dimensions in the THnF
  // vars     = variable enum pos. from AliReducedVarManager
  // indices  = dimensions in THnF
  //
  for(Int_t i=0; i<nVars; ++i) {
    if (fNMixingVariables>=kNMaxMixingVariables) {
      cout << "AliCorrelationExtraction::AddVariables() Fatal: Number of requested mixing variables " << fNMixingVariables << " larger than allowed number " << kNMaxMixingVariables << "!" << endl;
      return;
    }
    fMixingVariables[fNMixingVariables]   = vars[i];
    fMixingVarIndices[fNMixingVariables]  = indices[i];
    fNMixingVariables++;
  }
  fProcessDone = kFALSE;
  fUseMixingVars = kTRUE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::AddMixingVariable(Int_t var, Int_t index) {
  //
  // initialize a variable being used in the THnF
  // var    = variable enum pos. from AliReducedVarManager
  // index  = dimension in THnF
  //
  if (fNMixingVariables>=kNMaxMixingVariables) {
    cout << "AliCorrelationExtraction::AddVariable() Fatal: Number of requested mixing variables " << fNMixingVariables << " larger than allowed number " << kNMaxMixingVariables << "!" << endl;
    return;
  }
  fMixingVariables[fNMixingVariables]   = var;
  fMixingVarIndices[fNMixingVariables]  = index;
  fNMixingVariables++;
  fProcessDone = kFALSE;
  fUseMixingVars = kTRUE;
}

//_______________________________________________________________________________
void AliCorrelationExtraction::ApplyUserRanges(THnBase* h) {
  //
  // apply user ranges to the THnF histogram
  // NOTE:  user ranges are applied only to the specified fNVariables variables
  //        THnF may contain more than fNVariables, but those unspecified will be automatically integrated over
  //
  std::unique_ptr<THnBase> seosHist;
  if (fSEOS)            seosHist = std::unique_ptr<THnF>(       static_cast<THnF*>(       fSEOS->Clone(       "seosHist")));
  else if (fSEOSSparse) seosHist = std::unique_ptr<THnSparseF>( static_cast<THnSparseF*>( fSEOSSparse->Clone( "seosHist")));
  for (Int_t i=0; i<fNVariables; ++i) {
    if (TMath::Abs(fVarLimits[i][0]-fVarLimits[i][1])<EPSILON) {
      fVarLimits[i][0] = seosHist->GetAxis(fVarIndices[i])->GetXmin()+EPSILON;
      fVarLimits[i][1] = seosHist->GetAxis(fVarIndices[i])->GetXmax()-EPSILON;
    }
    if (fVarLimits[i][0] < seosHist->GetAxis(fVarIndices[i])->GetXmin()) fVarLimits[i][0] = seosHist->GetAxis(fVarIndices[i])->GetXmin()+EPSILON;
    if (fVarLimits[i][1] > seosHist->GetAxis(fVarIndices[i])->GetXmax()) fVarLimits[i][1] = seosHist->GetAxis(fVarIndices[i])->GetXmax()-EPSILON;
    if (fVarLimits[i][0]>fVarLimits[i][1]) {
      fVarLimits[i][0] = seosHist->GetAxis(fVarIndices[i])->GetXmin()+EPSILON;
      fVarLimits[i][1] = seosHist->GetAxis(fVarIndices[i])->GetXmax()-EPSILON;
    }
    h->GetAxis(fVarIndices[i])->SetRangeUser(fVarLimits[i][0], fVarLimits[i][1]);
  }
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::Initialize() {
  //
  // initialize correlation extraction
  //
  AliReducedVarManager::SetDefaultVarNames();

  // clean up output histograms
  if (fSEOSNorm)          {delete fSEOSNorm;          fSEOSNorm           = 0x0;}
  if (fMEOSNorm)          {delete fMEOSNorm;          fMEOSNorm           = 0x0;}
  if (fSEPPPairInvMass)   {delete fSEPPPairInvMass;   fSEPPPairInvMass    = 0x0;}
  if (fSEMMPairInvMass)   {delete fSEMMPairInvMass;   fSEMMPairInvMass    = 0x0;}
  if (fMEOSPairInvMass)   {delete fMEOSPairInvMass;   fMEOSPairInvMass    = 0x0;}
  if (fMEPPPairInvMass)   {delete fMEPPPairInvMass;   fMEPPPairInvMass    = 0x0;}
  if (fMEMMPairInvMass)   {delete fMEMMPairInvMass;   fMEMMPairInvMass    = 0x0;}
  if (fInclusiveCF1D)     {delete fInclusiveCF1D;     fInclusiveCF1D      = 0x0;}
  if (fInclusiveCF2D)     {delete fInclusiveCF2D;     fInclusiveCF2D      = 0x0;}
  if (fInclusiveCF3D)     {delete fInclusiveCF3D;     fInclusiveCF3D      = 0x0;}
  if (fBackgroundCF1D)    {delete fBackgroundCF1D;    fBackgroundCF1D     = 0x0;}
  if (fBackgroundCF2D)    {delete fBackgroundCF2D;    fBackgroundCF2D     = 0x0;}
  if (fSignalCF1D)        {delete fSignalCF1D;        fSignalCF1D         = 0x0;}
  if (fSignalCF1DEffCorr) {delete fSignalCF1DEffCorr; fSignalCF1DEffCorr  = 0x0;}
  if (fSignalCF2D)        {delete fSignalCF2D;        fSignalCF2D         = 0x0;}
  if (fSignalCF2DEffCorr) {delete fSignalCF2DEffCorr; fSignalCF2DEffCorr  = 0x0;}
  for (Int_t i=0; i<kNMaxBackgroundMassRanges; ++i) {
    if (fSEOSNormBackgroundMassWindow[i])       {delete fSEOSNormBackgroundMassWindow[i];       fSEOSNormBackgroundMassWindow[i]      = 0x0;}
    if (fSEPPNormBackgroundMassWindow[i])       {delete fSEPPNormBackgroundMassWindow[i];       fSEPPNormBackgroundMassWindow[i]      = 0x0;}
    if (fSEMMNormBackgroundMassWindow[i])       {delete fSEMMNormBackgroundMassWindow[i];       fSEMMNormBackgroundMassWindow[i]      = 0x0;}
    if (fMEOSNormBackgroundMassWindow[i])       {delete fMEOSNormBackgroundMassWindow[i];       fMEOSNormBackgroundMassWindow[i]      = 0x0;}
    if (fMEOS2NormBackgroundMassWindow[i])      {delete fMEOS2NormBackgroundMassWindow[i];      fMEOS2NormBackgroundMassWindow[i]     = 0x0;}
    if (fMEPPNormBackgroundMassWindow[i])       {delete fMEPPNormBackgroundMassWindow[i];       fMEPPNormBackgroundMassWindow[i]      = 0x0;}
    if (fMEMMNormBackgroundMassWindow[i])       {delete fMEMMNormBackgroundMassWindow[i];       fMEMMNormBackgroundMassWindow[i]      = 0x0;}
    if (fInclusiveCF1DBackgroundMassWindow[i])  {delete fInclusiveCF1DBackgroundMassWindow[i];  fInclusiveCF1DBackgroundMassWindow[i] = 0x0;}
    if (fInclusiveCF2DBackgroundMassWindow[i])  {delete fInclusiveCF2DBackgroundMassWindow[i];  fInclusiveCF2DBackgroundMassWindow[i] = 0x0;}
    if (fCombinatorialBackgroundCF1D[i])        {delete fCombinatorialBackgroundCF1D[i];        fCombinatorialBackgroundCF1D[i]       = 0x0;}
    if (fCombinatorialBackgroundCF2D[i])        {delete fCombinatorialBackgroundCF2D[i];        fCombinatorialBackgroundCF2D[i]       = 0x0;}
  }
  for (Int_t i=0; i<kNMaxDeltaPhiBins; ++i) {
    for (Int_t j=0; j<kNMaxDeltaEtaBins; ++j) {
      if (fInclusiveCF1DInvMass[i][j])            {delete fInclusiveCF1DInvMass[i][j];            fInclusiveCF1DInvMass[i][j]           = 0x0;}
      if (fInclusiveCF1DInvMassBackground[i][j])  {delete fInclusiveCF1DInvMassBackground[i][j];  fInclusiveCF1DInvMassBackground[i][j] = 0x0;}
      if (fBackgroundCF1DInvMassFit[i][j])        {delete fBackgroundCF1DInvMassFit[i][j];        fBackgroundCF1DInvMassFit[i][j]       = 0x0;}
      if (fBackgroundCF1DInvMass[i][j])           {delete fBackgroundCF1DInvMass[i][j];           fBackgroundCF1DInvMass[i][j]          = 0x0;}
      if (fSignalCF1DInvMass[i][j])               {delete fSignalCF1DInvMass[i][j];               fSignalCF1DInvMass[i][j]              = 0x0;}
    }
  }

  // check for required user provided histograms
  if (!fSEOS && !fSEOSSparse) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: No SE-OS histogram provided! This is always needed!" << endl;
    return kFALSE;
  }
  if (!fMEOS && !fMEOSSparse) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: No ME-OS histogram provided! This is always needed!" << endl;
    return kFALSE;
  }

  // check for required user options
  if(fNVariables<1) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: No dimensions have been defined. Use AddVariable(s)!" << endl;
    return kFALSE;
  }
  if (fOptionBkgMethod==kBkgNone) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: No background determination method specified!" << endl;
    return kFALSE;
  }
  if (!fResonanceFits) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: No resonance fits object provided!" << endl;
    return kFALSE;
  }
  if (fOptionBkgMethod==kBkgSideband && fNBackgroundMassRanges>1) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: More than 1 background mass range provided for sideband method!" << endl;
    return kFALSE;
  }
  if (fOptionBkgMethod==kBkgInterpolation && fNBackgroundMassRanges<2) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: Less than 2 background mass ranges provided for interpolation method!" << endl;
    return kFALSE;
  }
  if (fOptionBkgMethod==kBkgFitting) {
    if (fNBackgroundMassRanges<2) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: Less than 2 background mass ranges provided for fitting method!" << endl;
      return kFALSE;
    }
    if (!fBkgFitFunction) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: Fit function missing! This is required for fitting method!" << endl;
      return kFALSE;
    }
  }
  if (fOptionBkgMethod==kBkgLikeSign) {
    if (fNBackgroundMassRanges!=1) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: More/less than 1 background mass ranges provided for likesign method!" << endl;
      return kFALSE;
    }
    if ((!fSEPP || !fSEMM || !fMEPP || !fMEMM) && (!fSEPPSparse || !fSEMMSparse || !fMEPPSparse || !fMEMMSparse)) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: SE-LS and/or ME-LS histograms missing! These are required for like-sign method!" << endl;
      return kFALSE;
    }
    if (!fSEPPPair || !fSEMMPair) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: SE-LS pair histograms missing! These are required for like-sign method!" << endl;
      return kFALSE;
    }
    if (fMassVariableIndexPair<0) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: Mass variable index for LS pair histograms not defined! This is required for like-sign method!" << endl;
      return kFALSE;
    }
  }
  if (fOptionBkgMethod==kBkgSuperposition && fNBackgroundMassRanges!=2) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: More/less than 2 background mass ranges provided for superposition method!" << endl;
    return kFALSE;
  }
  if (fOptionBkgMethod==kBkgSuperpositionTwoComponent) {
    if (fNBackgroundMassRanges!=2) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: More/less than 2 background mass ranges provided for 2-component superposition method!" << endl;
      return kFALSE;
    }
    if (!fSEPPPair || !fSEMMPair) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: SE-LS pair histograms missing! These are required for 2-component superposition method!" << endl;
      return kFALSE;
    }
    if (!fMEOSPair || !fMEPPPair || !fMEMMPair) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: ME-OS or ME-LS pair histograms missing! These are required for 2-component superposition method!" << endl;
      return kFALSE;
    }
    if (!fMEOS2 && !fMEOS2Sparse) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: Second ME-OS histogram missing! This is required for 2-component superposition method!" << endl;
      return kFALSE;
    }
  }
  
  // check for inv. mass, delta phi and delta eta variables
  for(Int_t i=0; i<fNVariables; ++i) {
    if(fVariables[i]==fMassVariable)      fMassVariableIndex      = fVarIndices[i];
    if(fVariables[i]==fDeltaPhiVariable)  fDeltaPhiVariableIndex  = fVarIndices[i];
    if(fVariables[i]==fDeltaEtaVariable)  fDeltaEtaVariableIndex  = fVarIndices[i];
  }
  if (fMassVariableIndex<0) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: Mass dimension not found in list of user defined dimensions!" << endl;
    return kFALSE;
  }
  if (fDeltaPhiVariableIndex<0) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: Delta phi dimension not found in list of user defined dimensions!" << endl;
    return kFALSE;
  }
  if (fDeltaEtaVariableIndex<0) {
    cout << "AliCorrelationExtraction::Initialize() Fatal: Delta eta dimension not found in list of user defined dimensions!" << endl;
    return kFALSE;
  }

  // check for efficiency variable
  if (fHadronEff) {
    for(Int_t i=0; i<fNVariables; ++i) {
      if (fVariables[i]==fEfficiencyVariable) fEfficiencyVariableIndex = fVarIndices[i];
    }
    if (fEfficiencyVariableIndex<0) {
      cout << "AliCorrelationExtraction::Initialize() Fatal: Efficiency dimension not found in list of user defined dimensions!" << endl;
      return kFALSE;
    }
  }
  
  // apply user ranges
  // NOTE: user ranges for LS pair histograms (i.e. fSEPPPair, fSEMMPair, fMEOSPair, fMEPPPair, fMEMMPair) have to be set outside the class
  if (fSEOS)        ApplyUserRanges(fSEOS);
  if (fSEOSSparse)  ApplyUserRanges(fSEOSSparse);
  if (fMEOS)        ApplyUserRanges(fMEOS);
  if (fMEOSSparse)  ApplyUserRanges(fMEOSSparse);
  if (fMEOS2)       ApplyUserRanges(fMEOS2);
  if (fMEOS2Sparse) ApplyUserRanges(fMEOS2Sparse);
  if (fSEPP)        ApplyUserRanges(fSEPP);
  if (fSEPPSparse)  ApplyUserRanges(fSEPPSparse);
  if (fSEMM)        ApplyUserRanges(fSEMM);
  if (fSEMMSparse)  ApplyUserRanges(fSEMMSparse);
  if (fMEPP)        ApplyUserRanges(fMEPP);
  if (fMEPPSparse)  ApplyUserRanges(fMEPPSparse);
  if (fMEMM)        ApplyUserRanges(fMEMM);
  if (fMEMMSparse)  ApplyUserRanges(fMEMMSparse);

  // set mixing variable bin limits
  std::unique_ptr<TH1D> tmpHist;
  for (Int_t i=0; i<fNMixingVariables; ++i) {
    if (fSEOS)        tmpHist = std::unique_ptr<TH1D>(static_cast<TH1D*>(fSEOS->Projection(       fMixingVarIndices[i])));
    if (fSEOSSparse)  tmpHist = std::unique_ptr<TH1D>(static_cast<TH1D*>(fSEOSSparse->Projection( fMixingVarIndices[i])));
    tmpHist->SetName(Form("tmpHist_%.6f", gRandom->Rndm()));
    fNMixingVarBins[i] = tmpHist->GetNbinsX();
    for (Int_t bin=0; bin<fNMixingVarBins[i]; ++bin) {
      fMixingVarBinLimits[i][bin][0] = tmpHist->GetXaxis()->GetBinLowEdge(bin+1);
      fMixingVarBinLimits[i][bin][1] = tmpHist->GetXaxis()->GetBinUpEdge(bin+1);
    }
  }
  
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::NormalizeToNearSidePeak(TH2D* h) {
  //
  // normalize mixed-event (or any) 2D distribution to the near-side peak (delta phi, delta eta) = (0, 0)
  // NOTE: 4 bins around (0,0) are checked/used in case the bin boundaries are exactly at (0,0)
  // NOTE: no error propagation is done
  //
  if (!h) return kFALSE;
  Double_t signsA[4]  = {1.0, 1.0, -1.0, -1.0};
  Double_t signsB[4]  = {1.0, -1.0, -1.0, 1.0};
  std::vector<Int_t> binNumbers;
  for (Int_t i=0; i<4; i++) {
    Double_t tempX = 0.+signsA[i]*EPSILON;
    Double_t tempY = 0.+signsB[i]*EPSILON;
    if (h->GetXaxis()->GetXmin()>tempX || h->GetXaxis()->GetXmax()<tempX) continue;
    if (h->GetYaxis()->GetXmin()>tempY || h->GetYaxis()->GetXmax()<tempY) continue;
    Int_t tempBin = h->FindBin(tempX, tempY);
    if (std::find(binNumbers.begin(), binNumbers.end(), tempBin) != binNumbers.end()) continue;
    binNumbers.push_back(tempBin);
  }
  if (!binNumbers.size()) return kFALSE;
  Double_t norm = 0.;
  for (Int_t i=0; i<binNumbers.size(); ++i) norm += h->GetBinContent(binNumbers.at(i));
  norm /= binNumbers.size();
  binNumbers.clear();
  if (!norm) return kFALSE;
  h->Scale(1./norm);
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t  AliCorrelationExtraction::InBackgroundRange(Double_t min, Double_t max, Int_t& index) {
  //
  // test if [min, max] is within a given background mass range
  //
  for (Int_t i=0; i<fNBackgroundMassRanges; ++i) {
    if (fBackgroundMassRanges[i][0]<=min && max<=fBackgroundMassRanges[i][1]) {
      index = i;
      return kTRUE;
    }
  }
  return kFALSE;
}

//_______________________________________________________________________________
TH1D* AliCorrelationExtraction::ProjectToDeltaPhi(TH2D* hIn, TString name) {
  //
  // project 2D (delta eta, delta phi) distribution into 1D (delta phi)
  //
  if (!hIn) return NULL;
  Double_t deltaEtaRange = hIn->GetXaxis()->GetXmax() - hIn->GetXaxis()->GetXmin();
  std::unique_ptr<TH2D> hTmp = std::unique_ptr<TH2D>(static_cast<TH2D*>(hIn->Clone(Form("hTmp_%.6f", gRandom->Rndm()))));
  for (Int_t etaBin=1; etaBin<hTmp->GetNbinsX()+1; ++etaBin) {
    Double_t deltaEtaBinWidth = hTmp->GetXaxis()->GetBinWidth(etaBin);
    for (Int_t phiBin=1; phiBin<hTmp->GetNbinsY()+1; ++phiBin) {
      Double_t binContent = hTmp->GetBinContent( etaBin, phiBin);
      Double_t binError   = hTmp->GetBinError(   etaBin, phiBin);
      hTmp->SetBinContent( etaBin, phiBin, binContent*deltaEtaBinWidth);
      hTmp->SetBinError(   etaBin, phiBin, binError*deltaEtaBinWidth);
    }
  }
  TH1D* hOut = (TH1D*)hTmp->ProjectionY(name, 1, hTmp->GetNbinsX(), "e");
  hOut->Scale(1./deltaEtaRange);
  return hOut;
}

//_______________________________________________________________________________
Bool_t  AliCorrelationExtraction::CalculateInclusiveCorrelationInMixingBins(Int_t currentVar, Int_t& nCalls,
                                                                            THnBase* seos, THnBase* meos, TH2D* (&inclCF)) {
  //
  // calculate inclusive correlation in mixing bins
  //
  if (!seos || !meos) return kFALSE;
  if (currentVar>=fNMixingVariables) return kFALSE;
  if (nCalls && !inclCF) return kFALSE;
  
  // loop over mixinga variable bins
  std::unique_ptr<TH2D> seosTmp;
  std::unique_ptr<TH2D> meosTmp;
  std::unique_ptr<TH2D> inclCFTmp;
  for (Int_t bin=0; bin<fNMixingVarBins[currentVar]; ++bin) {
    if (currentVar<fNMixingVariables-1) {
      CalculateInclusiveCorrelationInMixingBins(currentVar+1, nCalls, seos, meos, inclCF);
    } else {
      // set mixing variable range
      seos->GetAxis(fMixingVarIndices[currentVar])->SetRangeUser(fMixingVarBinLimits[currentVar][bin][0]+EPSILON,
                                                                 fMixingVarBinLimits[currentVar][bin][1]-EPSILON);
      meos->GetAxis(fMixingVarIndices[currentVar])->SetRangeUser(fMixingVarBinLimits[currentVar][bin][0]+EPSILON,
                                                                 fMixingVarBinLimits[currentVar][bin][1]-EPSILON);

      // project
      seosTmp = std::unique_ptr<TH2D>(static_cast<TH2D*>(seos->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e")));
      meosTmp = std::unique_ptr<TH2D>(static_cast<TH2D*>(meos->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e")));
      seosTmp->SetName(Form("projSEOS_%.6f", gRandom->Rndm()));
      meosTmp->SetName(Form("projMEOS_%.6f", gRandom->Rndm()));

      // intercept empty histograms
      if (!seosTmp->GetEntries() || !meosTmp->GetEntries()) {
        cout << "AliCorrelationExtraction::CalculateInclusiveCorrelationInMixingBins() Warning: Empty histogram encountered, skipping mixing bin!" << endl;
        continue;
      }
      
      // normalize SE-OS
      if (fDeltaPhiVariable==AliReducedVarManager::kDeltaPhiSym) seosTmp->Scale(1./2.);  // normalize to symmetric delta phi range
      if (fDeltaEtaVariable==AliReducedVarManager::kDeltaEtaAbs) seosTmp->Scale(1./2.);  // normalize to absolute delta eta range, NOTE: do we need this?

      // normalize ME-OS
      if (!NormalizeToNearSidePeak(meosTmp.get())) {
        cout << "AliCorrelationExtraction::CalculateInclusiveCorrelationInMixingBins() Warning: ME-OS normalization failed, skipping mixing bin!" << endl;
        continue;
      }
      
      // calculate incl CF for current step in iteration
      inclCFTmp = std::unique_ptr<TH2D>(static_cast<TH2D*>(seosTmp->Clone(Form("inclCFTmp_%.6f", gRandom->Rndm()))));
      inclCFTmp->Divide(meosTmp.get());

      // sum bins
      if (!nCalls)  inclCF = (TH2D*)inclCFTmp->Clone(Form("inclCF_%.6f", gRandom->Rndm()));
      else          inclCF->Add(inclCFTmp.get());
      nCalls++;
    }
  }
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::CalculateInclusiveCorrelation(Double_t minMass, Double_t maxMass,
                                                               TH2D* (&seos), TH2D* (&meos),
                                                               TH2D* (&inclCF2D), TH1D* (&inclCF1D)) {
  //
  // calculate inclusive (signal+background) correlation in J/psi signal region
  // NOTE: error propagation for S+B error required in normalization step?
  //
  
  // calculate J/psi signal values in required mass window
  Double_t* vals = fResonanceFits->ComputeOutputValues(minMass, maxMass);
  if (fVerboseFlag) {
    cout << "AliCorrelationExtraction::CalculateInclusiveCorrelation() called for mass interval: " << minMass << " - " << maxMass << endl;
    fResonanceFits->PrintFitValues();
  }
  if (vals[AliResonanceFits::kSplusB]<=0) {
    cout << "AliCorrelationExtraction::CalculateInclusiveCorrelation() Fatal: AliResonanceFits::kSplusB <= 0!" << endl;
    seos      = NULL;
    meos      = NULL;
    inclCF2D  = NULL;
    inclCF1D  = NULL;
    return kFALSE;
  }
  
  // clone THnF and set mass range to signal region
  std::unique_ptr<THnBase> seostmp;
  std::unique_ptr<THnBase> meostmp;
  if (fSEOS)            seostmp = std::unique_ptr<THnF>(      static_cast<THnF*>(       fSEOS->Clone(       Form("seostmp_%.6f", gRandom->Rndm()))));
  else if (fSEOSSparse) seostmp = std::unique_ptr<THnSparseF>(static_cast<THnSparseF*>( fSEOSSparse->Clone( Form("seostmp_%.6f", gRandom->Rndm()))));
  if (fMEOS)            meostmp = std::unique_ptr<THnF>(      static_cast<THnF*>(       fMEOS->Clone(       Form("meostmp_%.6f", gRandom->Rndm()))));
  else if (fMEOSSparse) meostmp = std::unique_ptr<THnSparseF>(static_cast<THnSparseF*>( fMEOSSparse->Clone( Form("meostmp_%.6f", gRandom->Rndm()))));
  seostmp->GetAxis(fMassVariableIndex)->SetRangeUser(minMass+EPSILON, maxMass-EPSILON);
  meostmp->GetAxis(fMassVariableIndex)->SetRangeUser(minMass+EPSILON, maxMass-EPSILON);
  
  // project into 2D
  // NOTE: x = fDeltaEtaVariableIndex, y = fDeltaPhiVariableIndex
  seos = (TH2D*)seostmp->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e");
  meos = (TH2D*)meostmp->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e");
  seos->SetName(Form("projSEOS_%.6f", gRandom->Rndm()));
  meos->SetName(Form("projMEOS_%.6f", gRandom->Rndm()));
  
  // normalize same event
  seos->Scale(1./vals[AliResonanceFits::kSplusB]);                                // normalize to number of triggers
  if (fDeltaPhiVariable==AliReducedVarManager::kDeltaPhiSym) seos->Scale(1./2.);  // normalize to symmetric delta phi range
  if (fDeltaEtaVariable==AliReducedVarManager::kDeltaEtaAbs) seos->Scale(1./2.);  // normalize to absolute delta eta range, NOTE: do we need this?

  // normalize mixed event to near-side peak
  Bool_t meosNormFlag = NormalizeToNearSidePeak(meos);
  if (!meosNormFlag) {
    cout << "AliCorrelationExtraction::CalculateInclusiveCorrelation() Fatal: ME-OS normalization failed!" << endl;
    return kFALSE;
  }

  // calculate inclusive correlation
  if (!fUseMixingVars) {
    inclCF2D = (TH2D*)seos->Clone(Form("inclCF2D_%.6f", gRandom->Rndm()));
    inclCF2D->Divide(meos);
  } else {
    if (fVerboseFlag) cout << "AliCorrelationExtraction::CalculateInclusiveCorrelation() calculating inclusive correlation in mixing variable bins!" << endl;
    Int_t startVar  = 0;
    Int_t nCalls    = 0;
    CalculateInclusiveCorrelationInMixingBins(startVar, nCalls, seostmp.get(), meostmp.get(), inclCF2D);
    inclCF2D->SetName(Form("inclCF2D_%.6f", gRandom->Rndm()));
    inclCF2D->Scale(1./vals[AliResonanceFits::kSplusB]);
  }
  inclCF2D->Scale(1., "width"); // normalize to bin area

  // project 1D inclusive correlation from 2D
  inclCF1D = ProjectToDeltaPhi(inclCF2D, Form("inclCF1D_%.6f", gRandom->Rndm()));

  // fill trigger values
  Int_t massWindow = -1;
  if (InBackgroundRange(minMass, maxMass, massWindow)) {
    fTrigValBkg[massWindow][kSig]         = vals[AliResonanceFits::kSig];
    fTrigValBkg[massWindow][kSigErr]      = vals[AliResonanceFits::kSigErr];
    fTrigValBkg[massWindow][kBkg]         = vals[AliResonanceFits::kBkg];
    fTrigValBkg[massWindow][kBkgErr]      = vals[AliResonanceFits::kBkgErr];
    fTrigValBkg[massWindow][kSplusB]      = vals[AliResonanceFits::kSplusB];
    fTrigValBkg[massWindow][kSplusBErr]   = vals[AliResonanceFits::kSplusBerr];
    fTrigValBkg[massWindow][kSigFrac]     = fTrigValBkg[massWindow][kSig]/fTrigValBkg[massWindow][kSplusB];
    fTrigValBkg[massWindow][kSigFracErr]  = TMath::Sqrt((TMath::Power(fTrigValBkg[massWindow][kBkg]*fTrigValBkg[massWindow][kSigErr], 2) +
                                                         TMath::Power(fTrigValBkg[massWindow][kSig]*fTrigValBkg[massWindow][kBkgErr], 2)) /
                                                        TMath::Power(fTrigValBkg[massWindow][kSplusB], 4));
    fTrigValBkg[massWindow][kBkgFrac]     = fTrigValBkg[massWindow][kBkg]/fTrigValBkg[massWindow][kSplusB];
    fTrigValBkg[massWindow][kBkgFracErr]  = fTrigValBkg[massWindow][kSigFracErr];
  } else if (minMass==fMassSignalRange[0] && maxMass==fMassSignalRange[1]) {
    fTrigValSig[kSig]                     = vals[AliResonanceFits::kSig];
    fTrigValSig[kSigErr]                  = vals[AliResonanceFits::kSigErr];
    fTrigValSig[kBkg]                     = vals[AliResonanceFits::kBkg];
    fTrigValSig[kBkgErr]                  = vals[AliResonanceFits::kBkgErr];
    fTrigValSig[kSplusB]                  = vals[AliResonanceFits::kSplusB];
    fTrigValSig[kSplusBErr]               = vals[AliResonanceFits::kSplusBerr];
    fTrigValSig[kSigFrac]                 = fTrigValSig[kSig]/fTrigValSig[kSplusB];
    fTrigValSig[kSigFracErr]              = TMath::Sqrt((TMath::Power(fTrigValSig[kBkg]*fTrigValSig[kSigErr], 2) +
                                                         TMath::Power(fTrigValSig[kSig]*fTrigValSig[kBkgErr], 2)) /
                                                        TMath::Power(fTrigValSig[kSplusB], 4));
    fTrigValSig[kBkgFrac]                 = fTrigValSig[kBkg]/fTrigValSig[kSplusB];
    fTrigValSig[kBkgFracErr]              = fTrigValSig[kSigFracErr];
  }

  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::CalculateBackgroundCorrelationFitting() {
  //
  // calculate background and signal correlation from fitting method
  // TODO: add fit range for inv. mass
  //
  
  // calculate inclusive correlation in inv. mass slices
  if (fSEOS)            fInclusiveCF3D = (TH3D*)fSEOS->Projection(fMassVariableIndex, fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e");
  else if (fSEOSSparse) fInclusiveCF3D = (TH3D*)fSEOSSparse->Projection(fMassVariableIndex, fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e");
  Int_t nMassBins = fInclusiveCF3D->GetNbinsX();
  Int_t nPhiBins  = fInclusiveCF3D->GetNbinsY();
  if (fIntegrateDeltaEta[kBkgFitting]) fInclusiveCF3D->RebinZ(fInclusiveCF3D->GetNbinsZ());
  Int_t nEtaBins  = fInclusiveCF3D->GetNbinsZ();

  std::unique_ptr<THnBase> seosthnftmp;
  std::unique_ptr<THnBase> meosthnftmp;
  std::unique_ptr<TH2D> seostmp;
  std::unique_ptr<TH2D> meostmp;
  std::unique_ptr<TH2D> inclCF2Dtmp;
  for (Int_t massBin=1; massBin<=nMassBins; ++massBin) {
    // mass range
    Double_t minMass = fInclusiveCF3D->GetXaxis()->GetBinLowEdge( massBin);
    Double_t maxMass = fInclusiveCF3D->GetXaxis()->GetBinUpEdge(  massBin);
    
    if (fVerboseFlag) cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationFitting() Inclusive correlation slice in mass [" << minMass << ", " << maxMass << "]" << endl;
    
    // set range
    if (fSEOS)        seosthnftmp = std::unique_ptr<THnF>(      static_cast<THnF*>(       fSEOS->Clone(       Form("seosthnftmp_%d", massBin))));
    if (fSEOSSparse)  seosthnftmp = std::unique_ptr<THnSparseF>(static_cast<THnSparseF*>( fSEOSSparse->Clone( Form("seosthnftmp_%d", massBin))));
    if (fMEOS)        meosthnftmp = std::unique_ptr<THnF>(      static_cast<THnF*>(       fMEOS->Clone(       Form("meosthnftmp_%d", massBin))));
    if (fMEOSSparse)  meosthnftmp = std::unique_ptr<THnSparseF>(static_cast<THnSparseF*>( fMEOSSparse->Clone( Form("meosthnftmp_%d", massBin))));
    seosthnftmp->GetAxis(fMassVariableIndex)->SetRangeUser(minMass+EPSILON, maxMass-EPSILON);
    meosthnftmp->GetAxis(fMassVariableIndex)->SetRangeUser(minMass+EPSILON, maxMass-EPSILON);

    if (!fUseMixingVars) {
      // project into 2D
      // NOTE: x = fDeltaEtaVariableIndex, y = fDeltaPhiVariableIndex
      seostmp = std::unique_ptr<TH2D>(static_cast<TH2D*>(seosthnftmp->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e")));
      meostmp = std::unique_ptr<TH2D>(static_cast<TH2D*>(meosthnftmp->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e")));
      seostmp->SetName(Form("projSEOS_%df", massBin));
      meostmp->SetName(Form("projMEOS_%df", massBin));

      // normalize same event
      if (fDeltaPhiVariable==AliReducedVarManager::kDeltaPhiSym) seostmp->Scale(1./2.);  // normalize to symmetric delta phi range
      if (fDeltaEtaVariable==AliReducedVarManager::kDeltaEtaAbs) seostmp->Scale(1./2.);  // normalize to absolute delta eta range
      
      // normalize mixed event to near-side peak
      Bool_t meosNormFlag = NormalizeToNearSidePeak(meostmp.get());
      if (!meosNormFlag) {
        cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationFitting() Fatal: ME-OS normalization failed!" << endl;
        return kFALSE;
      }

      // calculate inclusive correlation
      inclCF2Dtmp = std::unique_ptr<TH2D>(static_cast<TH2D*>(seostmp->Clone(Form("inclCF2D_%d", massBin))));
      inclCF2Dtmp->Divide(meostmp.get());
    } else {
      if (fVerboseFlag) cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationFitting() calculating inclusive correlation in mixing variable bins!" << endl;
      Int_t startVar  = 0;
      Int_t nCalls    = 0;

      TH2D* inclCF2Dtmp_mixing = NULL;
      CalculateInclusiveCorrelationInMixingBins(startVar, nCalls, seosthnftmp.get(), meosthnftmp.get(), inclCF2Dtmp_mixing);
      inclCF2Dtmp = std::unique_ptr<TH2D>(static_cast<TH2D*>(inclCF2Dtmp_mixing->Clone(Form("inclCF2D_%d", massBin))));
      delete inclCF2Dtmp_mixing;
    }
    if (fIntegrateDeltaEta[kBkgFitting]) inclCF2Dtmp->RebinX(inclCF2Dtmp->GetNbinsX());

    // fill 3D inclusive correlation histogram
    for (Int_t phiBin=1; phiBin<=nPhiBins; ++phiBin) {
      for (Int_t etaBin=1; etaBin<=nEtaBins; ++etaBin) {
        Double_t binContent = inclCF2Dtmp->GetBinContent( etaBin, phiBin);
        Double_t binError   = inclCF2Dtmp->GetBinError(   etaBin, phiBin);
        fInclusiveCF3D->SetBinContent(massBin, phiBin, etaBin, binContent);
        fInclusiveCF3D->SetBinError(  massBin, phiBin, etaBin, binError);
      }
    }
  }
  
  // project 3D inclusive correlation into inv. mass dim. (in delta phi, delta eta bins)
  for (Int_t phiBin=0; phiBin<kNMaxDeltaPhiBins; ++phiBin) {
    for (Int_t etaBin=0; etaBin<kNMaxDeltaEtaBins; ++etaBin) {
      if (phiBin<fInclusiveCF3D->GetNbinsY() && etaBin<fInclusiveCF3D->GetNbinsZ())
        fInclusiveCF1DInvMass[phiBin][etaBin] = (TH1D*)fInclusiveCF3D->ProjectionX(Form("inclusiveCF_3D_projection_eta%d_phi%d", etaBin, phiBin),
                                                                                   phiBin+1, phiBin+1, etaBin+1, etaBin+1, "e");
      else fInclusiveCF1DInvMass[phiBin][etaBin] = NULL;
    }
  }
  
  // define signal and background correlation histogram
  fBackgroundCF2D = (TH2D*)fInclusiveCF3D->Project3D("yzoe");
  fBackgroundCF2D->SetName("backgroundCF_2D");
  fBackgroundCF2D->Reset("ICESM");
  fSignalCF2D = (TH2D*)fInclusiveCF3D->Project3D("yzoe");
  fSignalCF2D->SetName("signalCF_2D");
  fSignalCF2D->Reset("ICESM");

  // fit inv. mass using specified background ranges
  for (Int_t phiBin=0; phiBin<nPhiBins; ++phiBin) {
    for (Int_t etaBin=0; etaBin<nEtaBins; ++etaBin) {
      
      // exclude bins outside specified background range from fit
      fInclusiveCF1DInvMassBackground[phiBin][etaBin] = (TH1D*)fInclusiveCF1DInvMass[phiBin][etaBin]->Clone("inclusiveCF_3D_projection_bkgRegion_phi%d_eta%d");
      for (Int_t i=1; i<fInclusiveCF1DInvMassBackground[phiBin][etaBin]->GetNbinsX()+1; ++i) {
        Int_t     index   = 0;
        Double_t  binMin  = fInclusiveCF1DInvMassBackground[phiBin][etaBin]->GetXaxis()->GetBinLowEdge(i);
        Double_t  binMax  = fInclusiveCF1DInvMassBackground[phiBin][etaBin]->GetXaxis()->GetBinUpEdge(i);
        if (!InBackgroundRange(binMin, binMax, index)) {
          fInclusiveCF1DInvMassBackground[phiBin][etaBin]->SetBinContent( i, 0);
          fInclusiveCF1DInvMassBackground[phiBin][etaBin]->SetBinError(   i, 0);
        }
      }
      
      // fit distribution
      fInclusiveCF1DInvMassBackground[phiBin][etaBin]->Scale(1., "width"); // normalize by bin width prior to fitting
      fBackgroundCF1DInvMassFit[phiBin][etaBin] = (TF1*)fBkgFitFunction->Clone(Form("backgroundCF_1D_fit_eta%d_phi%d", etaBin, phiBin));
      fBackgroundCF1DInvMassFit[phiBin][etaBin]->SetNpx(10000.);
      TFitResultPtr fitResult = 0;
      if (fInclusiveCF1DInvMassBackground[phiBin][etaBin]->GetEntries()) {
        fitResult = fInclusiveCF1DInvMassBackground[phiBin][etaBin]->Fit(fBackgroundCF1DInvMassFit[phiBin][etaBin], "ISQMRNEF", "",
                                                                         fInclusiveCF1DInvMassBackground[phiBin][etaBin]->GetXaxis()->GetXmin(),
                                                                         fInclusiveCF1DInvMassBackground[phiBin][etaBin]->GetXaxis()->GetXmax());
      } else {
        fBackgroundCF1DInvMassFit[phiBin][etaBin] = NULL;
      }
      
      // write background correlation to histogram (inv. mass dist.)
      fBackgroundCF1DInvMass[phiBin][etaBin] = (TH1D*)fInclusiveCF1DInvMassBackground[phiBin][etaBin]->Clone(Form("backgroundCF_1D_eta%d_phi%d", etaBin, phiBin));
      for (Int_t i=1; i<fBackgroundCF1DInvMass[phiBin][etaBin]->GetNbinsX()+1; ++i) {
        Double_t minMass    = fBackgroundCF1DInvMass[phiBin][etaBin]->GetXaxis()->GetBinLowEdge(i);
        Double_t maxMass    = fBackgroundCF1DInvMass[phiBin][etaBin]->GetXaxis()->GetBinUpEdge(i);
        Double_t binContent = 0.;
        Double_t binError   = 0.;
        if (fBackgroundCF1DInvMassFit[phiBin][etaBin]) {
          binContent = fBackgroundCF1DInvMassFit[phiBin][etaBin]->Integral(minMass, maxMass);
          binError   = fBackgroundCF1DInvMassFit[phiBin][etaBin]->IntegralError(minMass, maxMass,
                                                                                fBackgroundCF1DInvMassFit[phiBin][etaBin]->GetParameters(),
                                                                                fitResult->GetCovarianceMatrix().GetMatrixArray());
        }
        fBackgroundCF1DInvMass[phiBin][etaBin]->SetBinContent( i, binContent);
        fBackgroundCF1DInvMass[phiBin][etaBin]->SetBinError(   i, binError);
      }
      
      // write signal correlation to histogram (inv. mass dist.)
      fSignalCF1DInvMass[phiBin][etaBin] = (TH1D*)fInclusiveCF1DInvMass[phiBin][etaBin]->Clone(Form("signalCF_1D_eta%d_phi%d", etaBin, phiBin));
      fSignalCF1DInvMass[phiBin][etaBin]->Add(fBackgroundCF1DInvMass[phiBin][etaBin], -1);
      
      // calculate background correlation
      Double_t bkg    = 0.;
      Double_t bkgErr = 0.;
      if (fBackgroundCF1DInvMassFit[phiBin][etaBin]) {
        bkg     = fBackgroundCF1DInvMassFit[phiBin][etaBin]->Integral(fMassSignalRange[0], fMassSignalRange[1]);
        bkgErr  = fBackgroundCF1DInvMassFit[phiBin][etaBin]->IntegralError(fMassSignalRange[0], fMassSignalRange[1],
                                                                           fBackgroundCF1DInvMassFit[phiBin][etaBin]->GetParameters(),
                                                                           fitResult->GetCovarianceMatrix().GetMatrixArray());
      }
      fBackgroundCF2D->SetBinContent( etaBin+1, phiBin+1, bkg);
      fBackgroundCF2D->SetBinError(   etaBin+1, phiBin+1, bkgErr);

      // calculate signal correlation
      Double_t sigPlusBkgErr  = 0.;
      Double_t sigPlusBkg     = fInclusiveCF1DInvMass[phiBin][etaBin]->IntegralAndError(fInclusiveCF1DInvMass[phiBin][etaBin]->FindBin(fMassSignalRange[0]+EPSILON),
                                                                                        fInclusiveCF1DInvMass[phiBin][etaBin]->FindBin(fMassSignalRange[1]-EPSILON),
                                                                                        sigPlusBkgErr);
      Double_t sig            = sigPlusBkg - bkg;
      Double_t sigErr         = TMath::Sqrt(sigPlusBkgErr*sigPlusBkgErr + bkgErr*bkgErr);
      // NOTE: signal correlation could also be taken directly from inv. mass dist. after bkg. sub
      //Double_t sigErr = 0.;
      //Double_t sig    = fSignalCF1DInvMass[phiBin][etaBin]->IntegralAndError(fSignalCF1DInvMass[phiBin][etaBin]->FindBin(fMassSignalRange[0]+EPSILON),
      //                                                                       fSignalCF1DInvMass[phiBin][etaBin]->FindBin(fMassSignalRange[1]-EPSILON),
      //                                                                       sigErr);
      fSignalCF2D->SetBinContent( etaBin+1, phiBin+1, sig);
      fSignalCF2D->SetBinError(   etaBin+1, phiBin+1, sigErr);
    }
  }
  
  // normalize background and signal correlation to number of triggers
  if (!fTrigValSig[kSig] || !fTrigValSig[kBkg]) {
    cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationFitting() Fatal: There was an issue with the J/psi signal extraction!" << endl;
    return kFALSE;
  }

  // scale by number of signal/background counts and propagate uncertainty
  for (Int_t xBin=1; xBin<=fSignalCF2D->GetNbinsX(); xBin++) {
    for (Int_t yBin=1; yBin<=fSignalCF2D->GetNbinsY(); yBin++) {
      Double_t bkgBinCont  = fBackgroundCF2D->GetBinContent( xBin, yBin);
      Double_t bkgBinErr   = fBackgroundCF2D->GetBinError(   xBin, yBin);
      Double_t sigBinCont  = fSignalCF2D->GetBinContent(     xBin, yBin);
      Double_t sigBinErr   = fSignalCF2D->GetBinError(       xBin, yBin);

      Double_t bkg    = bkgBinCont/fTrigValSig[kBkg];
      Double_t bkgErr = TMath::Sqrt(TMath::Power(bkgBinErr/fTrigValSig[kBkg], 2) + TMath::Power(bkgBinCont*fTrigValSig[kBkgErr]/fTrigValSig[kBkg]/fTrigValSig[kBkg], 2));
      Double_t sig    = sigBinCont/fTrigValSig[kSig];
      Double_t sigErr = TMath::Sqrt(TMath::Power(sigBinErr/fTrigValSig[kSig], 2) + TMath::Power(sigBinCont*fTrigValSig[kSigErr]/fTrigValSig[kSig]/fTrigValSig[kSig], 2));

      fBackgroundCF2D->SetBinContent( xBin, yBin, bkg);
      fBackgroundCF2D->SetBinError(   xBin, yBin, bkgErr);
      fSignalCF2D->SetBinContent(     xBin, yBin, sig);
      fSignalCF2D->SetBinError(       xBin, yBin, sigErr);
    }
  }
  
  // normalize background and signal correlation to bin width (delta eta, delta phi)
  fBackgroundCF2D->Scale( 1., "width");
  fSignalCF2D->Scale(     1., "width");

  // project 1D correlation from 2D
  if (!fIntegrateDeltaEta[kBkgFitting]) {
    fBackgroundCF1D = ProjectToDeltaPhi(fBackgroundCF2D,  "backgroundCF_1D");
    fSignalCF1D     = ProjectToDeltaPhi(fSignalCF2D,      "signalCF_1D");
  } else {
    fBackgroundCF1D = (TH1D*)fBackgroundCF2D->ProjectionY("backgroundCF_1D",  1, fBackgroundCF2D->GetNbinsX(),  "e");
    fSignalCF1D     = (TH1D*)fSignalCF2D->ProjectionY(    "signalCF_1D",      1, fSignalCF2D->GetNbinsX(),      "e");
  }
  
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::CalculateBackgroundCorrelationSideband() {
  //
  // calculate background correlation from inv. mass sideband
  //
  if (!CalculateInclusiveCorrelation(fBackgroundMassRanges[0][0], fBackgroundMassRanges[0][1],
                                     fSEOSNormBackgroundMassWindow[0], fMEOSNormBackgroundMassWindow[0],
                                     fInclusiveCF2DBackgroundMassWindow[0], fInclusiveCF1DBackgroundMassWindow[0])) return kFALSE;
  fSEOSNormBackgroundMassWindow[0]->SetName("SE-OS_backgroundRegion0_2D");
  fMEOSNormBackgroundMassWindow[0]->SetName("ME-OS_backgroundRegion0_2D");
  fInclusiveCF2DBackgroundMassWindow[0]->SetName("inclusiveCF_backgroundRegion0_2D");
  fInclusiveCF1DBackgroundMassWindow[0]->SetName("inclusiveCF_backgroundRegion0_1D");
  fBackgroundCF2D = (TH2D*)fInclusiveCF2DBackgroundMassWindow[0]->Clone("backgroundCF_2D");
  fBackgroundCF1D = (TH1D*)fInclusiveCF1DBackgroundMassWindow[0]->Clone("backgroundCF_1D");
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::CalculateBackgroundCorrelationLikeSign() {
  //
  // calculate background correlation from like-sign distributions
  // NOTE: difficult to set up nicely within this framework, all cross checks must be done outside the class
  // NOTE: future optimization welcome
  //
  
  // get number of triggers
  fSEPPPairInvMass = (TH1D*)fSEPPPair->Projection(fMassVariableIndexPair, "e");
  fSEPPPairInvMass->SetName("SE-PP_invMass_pair");
  fSEMMPairInvMass = (TH1D*)fSEMMPair->Projection(fMassVariableIndexPair, "e");
  fSEMMPairInvMass->SetName("SE-MM_invMass_pair");
  Double_t nTriggerPPErr  = 0.;
  Double_t nTriggerMMErr  = 0.;
  Double_t nTriggerPP     = fSEPPPairInvMass->IntegralAndError(fSEPPPairInvMass->FindBin(fBackgroundMassRanges[0][0]+EPSILON),
                                                               fSEPPPairInvMass->FindBin(fBackgroundMassRanges[0][1]-EPSILON),
                                                               nTriggerPPErr);
  Double_t nTriggerMM     = fSEMMPairInvMass->IntegralAndError(fSEMMPairInvMass->FindBin(fBackgroundMassRanges[0][0]+EPSILON),
                                                               fSEMMPairInvMass->FindBin(fBackgroundMassRanges[0][1]-EPSILON),
                                                               nTriggerMMErr);
  if (nTriggerPP<=0 || nTriggerMM<=0) {
    cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationLikeSign() Fatal: Number of LS triggers <= 0!" << endl;
    return kFALSE;
  }

  // fill trigger values
  // for LS method: fBackgroundMassRanges[0][0] = fSignalMassRange[0]
  //                fBackgroundMassRanges[0][1] = fSignalMassRange[1]
  fTrigValBkg[0][kNTrigSEPP]    = nTriggerPP;
  fTrigValBkg[0][kNTrigSEPPErr] = nTriggerPPErr;
  fTrigValBkg[0][kNTrigSEMM]    = nTriggerMM;
  fTrigValBkg[0][kNTrigSEMMErr] = nTriggerMMErr;
  fTrigValSig[kNTrigSEPP]       = nTriggerPP;
  fTrigValSig[kNTrigSEPPErr]    = nTriggerPPErr;
  fTrigValSig[kNTrigSEMM]       = nTriggerMM;
  fTrigValSig[kNTrigSEMMErr]    = nTriggerMMErr;

  // clone THnF (or THnSparseF) and set mass range to signal region
  std::unique_ptr<THnBase> sepptmp;
  std::unique_ptr<THnBase> semmtmp;
  std::unique_ptr<THnBase> mepptmp;
  std::unique_ptr<THnBase> memmtmp;
  if (fSEPP)        sepptmp = std::unique_ptr<THnF>(      static_cast<THnF*>(       fSEPP->Clone(       Form("sepptmp_%.6f", gRandom->Rndm()))));
  if (fSEPPSparse)  sepptmp = std::unique_ptr<THnSparseF>(static_cast<THnSparseF*>( fSEPPSparse->Clone( Form("sepptmp_%.6f", gRandom->Rndm()))));
  if (fSEMM)        semmtmp = std::unique_ptr<THnF>(      static_cast<THnF*>(       fSEMM->Clone(       Form("semmtmp_%.6f", gRandom->Rndm()))));
  if (fSEMMSparse)  semmtmp = std::unique_ptr<THnSparseF>(static_cast<THnSparseF*>( fSEMMSparse->Clone( Form("semmtmp_%.6f", gRandom->Rndm()))));
  if (fMEPP)        mepptmp = std::unique_ptr<THnF>(      static_cast<THnF*>(       fMEPP->Clone(       Form("mepptmp_%.6f", gRandom->Rndm()))));
  if (fMEPPSparse)  mepptmp = std::unique_ptr<THnSparseF>(static_cast<THnSparseF*>( fMEPPSparse->Clone( Form("mepptmp_%.6f", gRandom->Rndm()))));
  if (fMEMM)        memmtmp = std::unique_ptr<THnF>(      static_cast<THnF*>(       fMEMM->Clone(       Form("memmtmp_%.6f", gRandom->Rndm()))));
  if (fMEMMSparse)  memmtmp = std::unique_ptr<THnSparseF>(static_cast<THnSparseF*>( fMEMMSparse->Clone( Form("memmtmp_%.6f", gRandom->Rndm()))));
  sepptmp->GetAxis(fMassVariableIndex)->SetRangeUser(fBackgroundMassRanges[0][0]+EPSILON, fBackgroundMassRanges[0][1]-EPSILON);
  semmtmp->GetAxis(fMassVariableIndex)->SetRangeUser(fBackgroundMassRanges[0][0]+EPSILON, fBackgroundMassRanges[0][1]-EPSILON);
  mepptmp->GetAxis(fMassVariableIndex)->SetRangeUser(fBackgroundMassRanges[0][0]+EPSILON, fBackgroundMassRanges[0][1]-EPSILON);
  memmtmp->GetAxis(fMassVariableIndex)->SetRangeUser(fBackgroundMassRanges[0][0]+EPSILON, fBackgroundMassRanges[0][1]-EPSILON);

  // project into 2D
  // NOTE: x = fDeltaEtaVariableIndex, y = fDeltaPhiVariableIndex
  fSEPPNormBackgroundMassWindow[0] = (TH2D*)sepptmp->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e");
  fSEMMNormBackgroundMassWindow[0] = (TH2D*)semmtmp->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e");
  fMEPPNormBackgroundMassWindow[0] = (TH2D*)mepptmp->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e");
  fMEMMNormBackgroundMassWindow[0] = (TH2D*)memmtmp->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e");
  fSEPPNormBackgroundMassWindow[0]->SetName("SE-PP_backgroundRegion0_2D");
  fSEMMNormBackgroundMassWindow[0]->SetName("SE-MM_backgroundRegion0_2D");
  fMEPPNormBackgroundMassWindow[0]->SetName("ME-PP_backgroundRegion0_2D");
  fMEMMNormBackgroundMassWindow[0]->SetName("ME-MM_backgroundRegion0_2D");

  // normalize same event
  fSEPPNormBackgroundMassWindow[0]->Scale(1./nTriggerPP);       // normalize to number of triggers
  fSEMMNormBackgroundMassWindow[0]->Scale(1./nTriggerMM);       // normalize to number of triggers
  if (fDeltaPhiVariable==AliReducedVarManager::kDeltaPhiSym) {  // normalize to symmetric delta phi range
    fSEPPNormBackgroundMassWindow[0]->Scale(1./2.);
    fSEMMNormBackgroundMassWindow[0]->Scale(1./2.);
  }
  if (fDeltaEtaVariable==AliReducedVarManager::kDeltaEtaAbs) {  // normalize to absolute delta eta range, NOTE: do we need this?
    fSEPPNormBackgroundMassWindow[0]->Scale(1./2.);
    fSEMMNormBackgroundMassWindow[0]->Scale(1./2.);
  }

  // normalize mixed event to near-side peak
  Bool_t meosNormFlag = kFALSE;
  meosNormFlag = NormalizeToNearSidePeak(fMEPPNormBackgroundMassWindow[0]);
  meosNormFlag = NormalizeToNearSidePeak(fMEMMNormBackgroundMassWindow[0]);
  if (!meosNormFlag) {
    cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationLikeSign() Fatal: ME-OS normalization failed!" << endl;
    return kFALSE;
  }

  // calculate inclusive correlation - separately for both signs
  std::unique_ptr<TH2D> inclCF2DPP;
  std::unique_ptr<TH2D> inclCF2DMM;
  if (!fUseMixingVars) {
    inclCF2DPP = std::unique_ptr<TH2D>( static_cast<TH2D*>(fSEPPNormBackgroundMassWindow[0]->Clone(Form("inclCF2D_pp_%.6f", gRandom->Rndm()))));
    inclCF2DMM = std::unique_ptr<TH2D>( static_cast<TH2D*>(fSEMMNormBackgroundMassWindow[0]->Clone(Form("inclCF2D_mm_%.6f", gRandom->Rndm()))));
    inclCF2DPP->Divide(fMEPPNormBackgroundMassWindow[0]);
    inclCF2DMM->Divide(fMEMMNormBackgroundMassWindow[0]);
  } else {
    if (fVerboseFlag) cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationLikeSign() calculating inclusive correlation in mixing variable bins!" << endl;
    Int_t startVarPP  = 0;
    Int_t startVarMM  = 0;
    Int_t nCallsPP    = 0;
    Int_t nCallsMM    = 0;
    TH2D* inclCF2DPP_mixing = NULL;
    TH2D* inclCF2DMM_mixing = NULL;
    CalculateInclusiveCorrelationInMixingBins(startVarPP, nCallsPP, sepptmp.get(), mepptmp.get(), inclCF2DPP_mixing);
    CalculateInclusiveCorrelationInMixingBins(startVarMM, nCallsMM, semmtmp.get(), memmtmp.get(), inclCF2DMM_mixing);
    inclCF2DPP = std::unique_ptr<TH2D>( static_cast<TH2D*>(inclCF2DPP_mixing->Clone(Form("inclCF2D_pp_%.6f", gRandom->Rndm()))));
    inclCF2DMM = std::unique_ptr<TH2D>( static_cast<TH2D*>(inclCF2DMM_mixing->Clone(Form("inclCF2D_mm_%.6f", gRandom->Rndm()))));
    inclCF2DPP->Scale(1./nTriggerPP);
    inclCF2DMM->Scale(1./nTriggerMM);
    delete inclCF2DPP_mixing;
    delete inclCF2DMM_mixing;
  }
  inclCF2DPP->Scale(1., "width"); // normalize to bin area
  inclCF2DMM->Scale(1., "width"); // normalize to bin area

  // average ++ and -- contributions
  fInclusiveCF2DBackgroundMassWindow[0] = (TH2D*)inclCF2DPP->Clone("inclusiveCF_LS_backgroundRegion0_2D");
  fInclusiveCF2DBackgroundMassWindow[0]->Add(inclCF2DMM.get());
  fInclusiveCF2DBackgroundMassWindow[0]->Scale(0.5);
  
  // project inclusive LS correlation into 1D
  fInclusiveCF1DBackgroundMassWindow[0] = ProjectToDeltaPhi(fInclusiveCF2DBackgroundMassWindow[0], "inclusiveCF_LS_backgroundRegion0_1D");

  // set background correlation
  fBackgroundCF2D = (TH2D*)fInclusiveCF2DBackgroundMassWindow[0]->Clone("backgroundCF_2D");
  fBackgroundCF1D = (TH1D*)fInclusiveCF1DBackgroundMassWindow[0]->Clone("backgroundCF_1D");
  
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::CalculateBackgroundCorrelationInterpolation() {
  //
  // calculate background correlation from interpolation
  // NOTE: originally only designed for 2 background mass windows
  //
  
  // calculate inclusive correlation in background mass windows
  for (Int_t i=0; i<fNBackgroundMassRanges; ++i) {
    if (!CalculateInclusiveCorrelation(fBackgroundMassRanges[i][0], fBackgroundMassRanges[i][1],
                                       fSEOSNormBackgroundMassWindow[i], fMEOSNormBackgroundMassWindow[i],
                                       fInclusiveCF2DBackgroundMassWindow[i], fInclusiveCF1DBackgroundMassWindow[i])) return kFALSE;
  }
  
  // set histogram names
  for (Int_t i=0; i<fNBackgroundMassRanges; ++i) {
    fSEOSNormBackgroundMassWindow[i]->SetName(Form("SE-OS_backgroundRegion%d_2D", i));
    fMEOSNormBackgroundMassWindow[i]->SetName(Form("ME-OS_backgroundRegion%d_1D", i));
    fInclusiveCF2DBackgroundMassWindow[i]->SetName(Form("inclusiveCF_backgroundRegion%d_2D", i));
    fInclusiveCF1DBackgroundMassWindow[i]->SetName(Form("inclusiveCF_backgroundRegion%d_1D", i));
  }

  // calculate mass window centers
  Double_t centerSignal = (fMassSignalRange[0]+fMassSignalRange[1])/2.;
  Double_t centerBackground[kNMaxBackgroundMassRanges];
  for (Int_t i=0; i<kNMaxBackgroundMassRanges; ++i) {
    if (i<fNBackgroundMassRanges) centerBackground[i] = (fBackgroundMassRanges[i][0]+fBackgroundMassRanges[i][1])/2.;
    else                          centerBackground[i] = 0;
  }

  // calculate mass window widths
  Double_t widthBackground[kNMaxBackgroundMassRanges];
  for (Int_t i=0; i<kNMaxBackgroundMassRanges; ++i) {
    if (i<fNBackgroundMassRanges) widthBackground[i] = fBackgroundMassRanges[i][1] - fBackgroundMassRanges[i][0];
    else                          widthBackground[i] = 0;
  }

  // calculate weights
  Double_t weightBackground[kNMaxBackgroundMassRanges];
  for (Int_t i=0; i<kNMaxBackgroundMassRanges; ++i) {
    if (i<fNBackgroundMassRanges) weightBackground[i] = widthBackground[i]/TMath::Abs(centerSignal-centerBackground[i]);
    else                          weightBackground[i] = 0;
  }
  Double_t norm = 0.;
  for (Int_t i=0; i<kNMaxBackgroundMassRanges; ++i) {
    if (i<fNBackgroundMassRanges) norm += weightBackground[i];
  }
  if (!norm) return kFALSE;

  // calculate background correlation 2D
  fBackgroundCF2D = (TH2D*)fInclusiveCF2DBackgroundMassWindow[0]->Clone("backgroundCF_2D");
  fBackgroundCF2D->Scale(weightBackground[0]);
  for (Int_t i=1; i<fNBackgroundMassRanges; ++i) {
    fBackgroundCF2D->Add(fInclusiveCF2DBackgroundMassWindow[i], weightBackground[i]);
  }
  fBackgroundCF2D->Scale(1./norm);

  // project 1D background correlation from 2D
  fBackgroundCF1D = ProjectToDeltaPhi(fBackgroundCF2D, "backgroundCF_1D");
  
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::CalculateBackgroundCorrelationSuperposition() {
  //
  // calculate background and signal correlation from superposition principle
  //

  // calculate inclusive correlation in (background) mass windows
  for (Int_t i=0; i<fNBackgroundMassRanges; ++i) {
    if (!CalculateInclusiveCorrelation(fBackgroundMassRanges[i][0], fBackgroundMassRanges[i][1],
                                       fSEOSNormBackgroundMassWindow[i], fMEOSNormBackgroundMassWindow[i],
                                       fInclusiveCF2DBackgroundMassWindow[i], fInclusiveCF1DBackgroundMassWindow[i])) return kFALSE;
  }
  
  // set histogram names
  for (Int_t i=0; i<fNBackgroundMassRanges; ++i) {
    fSEOSNormBackgroundMassWindow[i]->SetName(Form("SE-OS_backgroundRegion%d_2D", i));
    fMEOSNormBackgroundMassWindow[i]->SetName(Form("ME-OS_backgroundRegion%d_2D", i));
    fInclusiveCF2DBackgroundMassWindow[i]->SetName(Form("inclusiveCF_backgroundRegion%d_2D", i));
    fInclusiveCF1DBackgroundMassWindow[i]->SetName(Form("inclusiveCF_backgroundRegion%d_1D", i));
  }

  //
  // NOTE:  from here on out it is more difficult to generalize to more than two windows
  //        but fNBackgroundMassRanges!=2 is currently intercepted in Initialize()
  //
  
  if (!(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac])) {
    cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationSuperposition() Fatal: There was an issue with the J/psi signal extraction!" << endl;
    return kFALSE;
  }

  // calculate background correlation 2D
  fBackgroundCF2D = (TH2D*)fInclusiveCF2DBackgroundMassWindow[0]->Clone("backgroundCF_2D");
  fBackgroundCF2D->Scale(-fTrigValBkg[1][kSigFrac]);
  fBackgroundCF2D->Add(fInclusiveCF2DBackgroundMassWindow[1], fTrigValBkg[0][kSigFrac]);
  fBackgroundCF2D->Scale(1./(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac]));

  // calculate signal correlation
  fSignalCF2D = (TH2D*)fInclusiveCF2DBackgroundMassWindow[0]->Clone("signalCF_2D");
  fSignalCF2D->Scale(1-fTrigValBkg[1][kSigFrac]);
  fSignalCF2D->Add(fInclusiveCF2DBackgroundMassWindow[1], -1+fTrigValBkg[0][kSigFrac]);
  fSignalCF2D->Scale(1./(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac]));
  
  // re-calculate uncertainties for signal and background correlation
  Double_t inclCFBinCont[2] = {0., 0.};
  Double_t inclCFBinErr[2]  = {0., 0.};
  for (Int_t xBin=1; xBin<=fBackgroundCF2D->GetNbinsX(); xBin++) {
    for (Int_t yBin=1; yBin<=fBackgroundCF2D->GetNbinsY(); yBin++) {
      // get bin contents and errors
      inclCFBinCont[0] = fInclusiveCF2DBackgroundMassWindow[0]->GetBinContent( xBin, yBin);
      inclCFBinCont[1] = fInclusiveCF2DBackgroundMassWindow[1]->GetBinContent( xBin, yBin);
      inclCFBinErr[0]  = fInclusiveCF2DBackgroundMassWindow[0]->GetBinError(   xBin, yBin);
      inclCFBinErr[1]  = fInclusiveCF2DBackgroundMassWindow[1]->GetBinError(   xBin, yBin);

      // calculate background uncertainties
      Double_t errA = fTrigValBkg[1][kSigFrac]/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac])*inclCFBinErr[0];
      Double_t errB = fTrigValBkg[0][kSigFrac]/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac])*inclCFBinErr[1];
      Double_t errC = fTrigValBkg[1][kSigFrac]/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac]) *
                      (inclCFBinCont[0]-inclCFBinCont[1])/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac]) *
                      fTrigValBkg[0][kSigFracErr];
      Double_t errD = fTrigValBkg[0][kSigFrac]/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac]) *
                      (inclCFBinCont[1]-inclCFBinCont[0])/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac]) *
                      fTrigValBkg[1][kSigFracErr];
      fBackgroundCF2D->SetBinError(xBin, yBin, TMath::Sqrt(errA*errA + errB*errB + errC*errC + errD*errD));

      // calculate signal uncertainties
      Double_t errE = (1-fTrigValBkg[1][kSigFrac])/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac])*inclCFBinErr[0];
      Double_t errF = (1-fTrigValBkg[0][kSigFrac])/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac])*inclCFBinErr[1];
      Double_t errG = (fTrigValBkg[1][kSigFrac]-1)/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac]) *
                      (inclCFBinCont[0]-inclCFBinCont[1])/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac]) *
                      fTrigValBkg[0][kSigFracErr];
      Double_t errH = (fTrigValBkg[0][kSigFrac]-1)/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac]) *
                      (inclCFBinCont[0]-inclCFBinCont[1])/(fTrigValBkg[0][kSigFrac]-fTrigValBkg[1][kSigFrac]) *
                      fTrigValBkg[1][kSigFracErr];
      fSignalCF2D->SetBinError(xBin, yBin, TMath::Sqrt(errE*errE + errF*errF + errG*errG + errH*errH));
    }
  }

  // project 1D correlation from 2D
  fBackgroundCF1D = ProjectToDeltaPhi(fBackgroundCF2D,  "backgroundCF_1D");
  fSignalCF1D     = ProjectToDeltaPhi(fSignalCF2D,      "signalCF_1D");
  
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::CalculateBackgroundCorrelationSuperpositionTwoComponent() {
  //
  // calculate background and signal correlation from superposition principle using a 2-component background
  //
  
  // calculate inclusive correlation in (background) mass windows
  for (Int_t i=0; i<fNBackgroundMassRanges; ++i) {
    if (!CalculateInclusiveCorrelation(fBackgroundMassRanges[i][0], fBackgroundMassRanges[i][1],
                                       fSEOSNormBackgroundMassWindow[i], fMEOSNormBackgroundMassWindow[i],
                                       fInclusiveCF2DBackgroundMassWindow[i], fInclusiveCF1DBackgroundMassWindow[i])) return kFALSE;
    
    // set histogram names
    fSEOSNormBackgroundMassWindow[i]->SetName(Form("SE-OS_backgroundRegion%d_2D", i));
    fMEOSNormBackgroundMassWindow[i]->SetName(Form("ME-OS_backgroundRegion%d_1D", i));
    fInclusiveCF2DBackgroundMassWindow[i]->SetName(Form("inclusiveCF_backgroundRegion%d_2D", i));
    fInclusiveCF1DBackgroundMassWindow[i]->SetName(Form("inclusiveCF_backgroundRegion%d_1D", i));
  }
  
  // get trigger inv. mass distributions
  fMEOSPairInvMass = (TH1D*)fMEOSPair->Projection(fMassVariableIndexPair, "e");
  fMEPPPairInvMass = (TH1D*)fMEPPPair->Projection(fMassVariableIndexPair, "e");
  fMEMMPairInvMass = (TH1D*)fMEMMPair->Projection(fMassVariableIndexPair, "e");
  fSEPPPairInvMass = (TH1D*)fSEPPPair->Projection(fMassVariableIndexPair, "e");
  fSEMMPairInvMass = (TH1D*)fSEMMPair->Projection(fMassVariableIndexPair, "e");
  fMEOSPairInvMass->SetName("ME-OS_invMass_pair");
  fMEPPPairInvMass->SetName("ME-PP_invMass_pair");
  fMEMMPairInvMass->SetName("ME-MM_invMass_pair");
  fSEPPPairInvMass->SetName("SE-PP_invMass_pair");
  fSEMMPairInvMass->SetName("SE-MM_invMass_pair");
  
  // calculate number of triggers
  Double_t  nTriggerSEPP[kNMaxBackgroundMassRanges];
  Double_t  nTriggerSEPPErr[kNMaxBackgroundMassRanges];
  Double_t  nTriggerSEMM[kNMaxBackgroundMassRanges];
  Double_t  nTriggerSEMMErr[kNMaxBackgroundMassRanges];
  Double_t  nTriggerMEOS[kNMaxBackgroundMassRanges];
  Double_t  nTriggerMEOSErr[kNMaxBackgroundMassRanges];
  Double_t  nTriggerMEPP[kNMaxBackgroundMassRanges];
  Double_t  nTriggerMEPPErr[kNMaxBackgroundMassRanges];
  Double_t  nTriggerMEMM[kNMaxBackgroundMassRanges];
  Double_t  nTriggerMEMMErr[kNMaxBackgroundMassRanges];
  for (Int_t i=0; i<kNMaxBackgroundMassRanges; ++i) {
    if (i<fNBackgroundMassRanges) {
      nTriggerMEOS[i] = fMEOSPairInvMass->IntegralAndError(fMEOSPairInvMass->FindBin(fBackgroundMassRanges[i][0]+EPSILON),
                                                           fMEOSPairInvMass->FindBin(fBackgroundMassRanges[i][1]-EPSILON),
                                                           nTriggerMEOSErr[i]);
      nTriggerMEPP[i] = fMEPPPairInvMass->IntegralAndError(fMEPPPairInvMass->FindBin(fBackgroundMassRanges[i][0]+EPSILON),
                                                           fMEPPPairInvMass->FindBin(fBackgroundMassRanges[i][1]-EPSILON),
                                                           nTriggerMEPPErr[i]);
      nTriggerMEMM[i] = fMEMMPairInvMass->IntegralAndError(fMEMMPairInvMass->FindBin(fBackgroundMassRanges[i][0]+EPSILON),
                                                           fMEMMPairInvMass->FindBin(fBackgroundMassRanges[i][1]-EPSILON),
                                                           nTriggerMEMMErr[i]);
      nTriggerSEPP[i] = fSEPPPairInvMass->IntegralAndError(fSEPPPairInvMass->FindBin(fBackgroundMassRanges[i][0]+EPSILON),
                                                           fSEPPPairInvMass->FindBin(fBackgroundMassRanges[i][1]-EPSILON),
                                                           nTriggerSEPPErr[i]);
      nTriggerSEMM[i] = fSEMMPairInvMass->IntegralAndError(fSEMMPairInvMass->FindBin(fBackgroundMassRanges[i][0]+EPSILON),
                                                           fSEMMPairInvMass->FindBin(fBackgroundMassRanges[i][1]-EPSILON),
                                                           nTriggerSEMMErr[i]);
      
      if (!nTriggerMEOS[i] || !nTriggerMEPP[i] || !nTriggerMEMM[i] || !nTriggerSEPP[i] || !nTriggerSEMM[i]) {
        cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationSuperpositionTwoComponent() Fatal: Number of triggers = 0 for inv. mass window " << fBackgroundMassRanges[i][0] << " - " << fBackgroundMassRanges[i][1] << "!" << endl;
        return kFALSE;
      }
    } else {
      nTriggerSEPP[i]     = 0.;
      nTriggerSEPPErr[i]  = 0.;
      nTriggerSEMM[i]     = 0.;
      nTriggerSEMMErr[i]  = 0.;
      nTriggerMEOS[i]     = 0.;
      nTriggerMEOSErr[i]  = 0.;
      nTriggerMEPP[i]     = 0.;
      nTriggerMEPPErr[i]  = 0.;
      nTriggerMEMM[i]     = 0.;
      nTriggerMEMMErr[i]  = 0.;
    }
  }
  
  // build 2D combinatorial background correlation
  std::unique_ptr<THnBase>  meos2tmp;
  std::unique_ptr<THnBase>  meostmp;
  std::unique_ptr<TH2D>     bkgtmp;
  for (Int_t i=0; i<fNBackgroundMassRanges; ++i) {
    if (i<fNBackgroundMassRanges) {
      if (fMEOS2)       meos2tmp  = std::unique_ptr<THnF>(      static_cast<THnF*>(       fMEOS2->Clone(      Form("meos2tmp_%.6f", gRandom->Rndm()))));
      if (fMEOS2Sparse) meos2tmp  = std::unique_ptr<THnSparseF>(static_cast<THnSparseF*>( fMEOS2Sparse->Clone(Form("meos2tmp_%.6f", gRandom->Rndm()))));
      if (fMEOS)        meostmp   = std::unique_ptr<THnF>(      static_cast<THnF*>(       fMEOS->Clone(       Form("meostmp_%.6f", gRandom->Rndm()))));
      if (fMEOSSparse)  meostmp   = std::unique_ptr<THnSparseF>(static_cast<THnSparseF*>( fMEOSSparse->Clone( Form("meostmp_%.6f", gRandom->Rndm()))));
      meos2tmp->GetAxis(fMassVariableIndex)->SetRangeUser(fBackgroundMassRanges[i][0]+EPSILON, fBackgroundMassRanges[i][1]-EPSILON);
      meostmp->GetAxis(fMassVariableIndex)->SetRangeUser(fBackgroundMassRanges[i][0]+EPSILON, fBackgroundMassRanges[i][1]-EPSILON);

      // projections
      fMEOS2NormBackgroundMassWindow[i] = (TH2D*)meos2tmp->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e");
      fMEOSNormBackgroundMassWindow[i]  = (TH2D*)meostmp->Projection(fDeltaPhiVariableIndex, fDeltaEtaVariableIndex, "e");
      fMEOS2NormBackgroundMassWindow[i]->SetName(Form("ME-OS-2_backgroundRegion%d_2D", i));
      fMEOSNormBackgroundMassWindow[i]->SetName(Form("ME-OS_backgroundRegion%d_2D", i));
        
      // normalize
      fMEOS2NormBackgroundMassWindow[i]->Scale(1./nTriggerMEOS[i]);                                               // normalize to number of triggers
      if (fDeltaPhiVariable==AliReducedVarManager::kDeltaPhiSym) fMEOS2NormBackgroundMassWindow[i]->Scale(1./2.); // normalize to symmetric delta phi range
      if (fDeltaEtaVariable==AliReducedVarManager::kDeltaEtaAbs) fMEOS2NormBackgroundMassWindow[i]->Scale(1./2.); // normalize to absolute delta eta range, NOTE: do we need this?
      Bool_t meosNormFlag = NormalizeToNearSidePeak(fMEOSNormBackgroundMassWindow[i]);                            // normalize mixed event to near-side peak
      if (!meosNormFlag) {
        cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationSuperpositionTwoComponent() Fatal: ME-OS normalization failed!" << endl;
        return kFALSE;
      }

      // calculate correlation
      if (!fUseMixingVars) {
        fCombinatorialBackgroundCF2D[i] = (TH2D*)fMEOS2NormBackgroundMassWindow[i]->Clone(Form("backgroundCF_combinatorial_backgroundRegion%d_2D", i));
        fCombinatorialBackgroundCF2D[i]->Divide(fMEOSNormBackgroundMassWindow[i]);
      } else {
        if (fVerboseFlag) cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationSuperpositionTwoComponent() calculating inclusive correlation in mixing variable bins!" << endl;
        Int_t startVar  = 0;
        Int_t nCalls    = 0;
        CalculateInclusiveCorrelationInMixingBins(startVar, nCalls, meos2tmp.get(), meostmp.get(), fCombinatorialBackgroundCF2D[i]);
        fCombinatorialBackgroundCF2D[i]->SetName(Form("backgroundCF_combinatorial_backgroundRegion%d_2D", i));
        fCombinatorialBackgroundCF2D[i]->Scale(1./nTriggerMEOS[i]);
      }
      fCombinatorialBackgroundCF2D[i]->Scale(1., "width"); // normalize to bin area

      // project 1D combinatorial background correlation from 2D
      Double_t deltaEtaRange = fCombinatorialBackgroundCF2D[i]->GetXaxis()->GetXmax() - fCombinatorialBackgroundCF2D[i]->GetXaxis()->GetXmin();
      bkgtmp = std::unique_ptr<TH2D>(static_cast<TH2D*>(fCombinatorialBackgroundCF2D[i]->Clone(Form("bkgCF2DTmp_%.6f", gRandom->Rndm()))));
      for (Int_t etaBin=1; etaBin<bkgtmp->GetNbinsX()+1; ++etaBin) {
        Double_t deltaEtaBinWidth = bkgtmp->GetXaxis()->GetBinWidth(etaBin);
        for (Int_t phiBin=1; phiBin<bkgtmp->GetNbinsY()+1; ++phiBin) {
          Double_t binContent = bkgtmp->GetBinContent( etaBin, phiBin);
          Double_t binError   = bkgtmp->GetBinError(   etaBin, phiBin);
          bkgtmp->SetBinContent( etaBin, phiBin, binContent*deltaEtaBinWidth);
          bkgtmp->SetBinError(   etaBin, phiBin, binError*deltaEtaBinWidth);
        }
      }
      fCombinatorialBackgroundCF1D[i] = (TH1D*)bkgtmp->ProjectionY(Form("backgroundCF_combinatorial_backgroundRegion%d_1D", i), 1, bkgtmp->GetNbinsX(), "e");
      fCombinatorialBackgroundCF1D[i]->Scale(1./deltaEtaRange);
    } else {
      fCombinatorialBackgroundCF2D[i] = NULL;
      fCombinatorialBackgroundCF1D[i] = NULL;
    }
  }
  
  // fill trigger values
  for (Int_t i=0; i<fNBackgroundMassRanges; ++i) {
    fTrigValBkg[i][kNTrigSEPP]      = nTriggerSEPP[i];
    fTrigValBkg[i][kNTrigSEPPErr]   = nTriggerSEPPErr[i];
    fTrigValBkg[i][kNTrigSEMM]      = nTriggerSEMM[i];
    fTrigValBkg[i][kNTrigSEMMErr]   = nTriggerSEMMErr[i];
    fTrigValBkg[i][kNTrigMEOS]      = nTriggerMEOS[i];
    fTrigValBkg[i][kNTrigMEOSErr]   = nTriggerMEOSErr[i];
    fTrigValBkg[i][kNTrigMEPP]      = nTriggerMEPP[i];
    fTrigValBkg[i][kNTrigMEPPErr]   = nTriggerMEPPErr[i];
    fTrigValBkg[i][kNTrigMEMM]      = nTriggerMEMM[i];
    fTrigValBkg[i][kNTrigMEMMErr]   = nTriggerMEMMErr[i];
    fTrigValBkg[i][kR]              = nTriggerMEOSErr[i]/(2*TMath::Sqrt(nTriggerMEPP[i]*nTriggerMEMM[i]));
    fTrigValBkg[i][kRErr]           = TMath::Sqrt(1/4/nTriggerMEPP[i]/nTriggerMEMM[i] *
                                                  (TMath::Power(nTriggerMEOSErr[i], 2) +
                                                   TMath::Power(nTriggerMEOS[i]*nTriggerMEPP[i]*nTriggerMEMMErr[i]/2/nTriggerMEPP[i]/nTriggerMEMM[i], 2) +
                                                   TMath::Power(nTriggerMEOS[i]*nTriggerMEMM[i]*nTriggerMEPPErr[i]/2/nTriggerMEPP[i]/nTriggerMEMM[i], 2)));
    fTrigValBkg[i][kBkgComb]        = 2*fTrigValBkg[i][kR]*TMath::Sqrt(nTriggerSEPP[i]*nTriggerSEMM[i]);
    fTrigValBkg[i][kBkgCombErr]     = TMath::Sqrt(4*nTriggerSEPP[i]*nTriggerSEMM[i]*TMath::Power(fTrigValBkg[i][kRErr], 2) +
                                                  nTriggerSEMM[i]*TMath::Power(fTrigValBkg[i][kR]*nTriggerSEPPErr[i], 2)/nTriggerSEPP[i] +
                                                  nTriggerSEPP[i]*TMath::Power(fTrigValBkg[i][kR]*nTriggerSEMMErr[i], 2)/nTriggerSEPP[i]);
    fTrigValBkg[i][kBkgCombFrac]    = fTrigValBkg[i][kBkgComb]/fTrigValBkg[i][kSplusB];
    fTrigValBkg[i][kBkgCombFracErr] = TMath::Sqrt(TMath::Power(fTrigValBkg[i][kBkgCombErr]/fTrigValBkg[i][kSplusB], 2) +
                                                  TMath::Power(fTrigValBkg[i][kBkgComb]*fTrigValBkg[i][kSplusBErr] /
                                                               TMath::Power(fTrigValBkg[i][kSplusB], 2), 2));
  }

  //
  // NOTE:  from here on out it is more difficult to generalize to more than two windows
  //        but fNBackgroundMassRanges!=2 is currently intercepted in Initialize()
  //
  
  // calculate normalization factor
  Double_t norm = fTrigValBkg[0][kSigFrac]*(1-fTrigValBkg[1][kBkgCombFrac])+fTrigValBkg[1][kSigFrac]*(-1+fTrigValBkg[0][kBkgCombFrac]);
  if (!norm) {
    cout << "AliCorrelationExtraction::CalculateBackgroundCorrelationSuperpositionTwoComponent() Fatal: There was an issue with the J/psi signal extraction!" << endl;
    return kFALSE;
  }
  
  // calculate residual background correlation
  //  C_B = (sigFrac[0] - sigFrac[0]*bkgFrac[1] - sigFrac[1] + sigFrac[1]*bkgFrac[0])^(-1) * (
  //          - sigFrac[1]*C_inc[0] + sigFrac[0]*C_inc[1]
  //          + sigFrac[1]*bkgFrac[0]*C_comb[0] - sigFrac[0]*bkgFrac[1]*C_comb[1] )
  fBackgroundCF2D = (TH2D*)fInclusiveCF2DBackgroundMassWindow[0]->Clone("backgroundCF_2D");
  fBackgroundCF2D->Scale(-fTrigValBkg[1][kSigFrac]);
  fBackgroundCF2D->Add(fInclusiveCF2DBackgroundMassWindow[1], fTrigValBkg[0][kSigFrac]);
  fBackgroundCF2D->Add(fCombinatorialBackgroundCF2D[0], fTrigValBkg[1][kSigFrac]*fTrigValBkg[0][kBkgCombFrac]);
  fBackgroundCF2D->Add(fCombinatorialBackgroundCF2D[1], -fTrigValBkg[0][kSigFrac]*fTrigValBkg[1][kBkgCombFrac]);
  fBackgroundCF2D->Scale(1./norm);
  
  // calculate signal correlation
  //  C_S = (sigFrac[0] - sigFrac[0]*bkgFrac[1] - sigFrac[1] + sigFrac[1]*bkgFrac[0])^(-1) * (
  //          + (1 - sigFrac[1] - bkgFrac[1])*C_inc[0] - (1 - sigFrac[0] - bkgFrac[0])*C_inc[1]
  //          - (1 - sigFrac[1] - bkgFrac[1])*bkgFrac[0]*C_comb[0] + (1 - sigFrac[0] - bkgFrac[0])*bkgFrac[1]*C_comb[1] )
  fSignalCF2D = (TH2D*)fInclusiveCF2DBackgroundMassWindow[0]->Clone("signalCF_2D");
  fSignalCF2D->Scale(1-fTrigValBkg[1][kSigFrac]-fTrigValBkg[1][kBkgCombFrac]);
  fSignalCF2D->Add(fInclusiveCF2DBackgroundMassWindow[1], -1+fTrigValBkg[0][kSigFrac]+fTrigValBkg[0][kBkgCombFrac]);
  fSignalCF2D->Add(fCombinatorialBackgroundCF2D[0], fTrigValBkg[0][kBkgCombFrac]*(-1+fTrigValBkg[1][kSigFrac]+fTrigValBkg[1][kBkgCombFrac]));
  fSignalCF2D->Add(fCombinatorialBackgroundCF2D[1], fTrigValBkg[1][kBkgCombFrac]*(1-fTrigValBkg[0][kSigFrac]-fTrigValBkg[0][kBkgCombFrac]));
  fSignalCF2D->Scale(1./norm);

  // re-calculate uncertainties for signal and background correlation
  Double_t inclCFBinCont[2] = {0., 0.};
  Double_t inclCFBinErr[2]  = {0., 0.};
  Double_t bkgCFBinCont[2]  = {0., 0.};
  Double_t bkgCFBinErr[2]   = {0., 0.};
  for (Int_t xBin=1; xBin<=fBackgroundCF2D->GetNbinsX(); xBin++) {
    for (Int_t yBin=1; yBin<=fBackgroundCF2D->GetNbinsY(); yBin++) {
      // get bin contents and errors
      inclCFBinCont[0]  = fInclusiveCF2DBackgroundMassWindow[0]->GetBinContent( xBin, yBin);
      inclCFBinCont[1]  = fInclusiveCF2DBackgroundMassWindow[1]->GetBinContent( xBin, yBin);
      inclCFBinErr[0]   = fInclusiveCF2DBackgroundMassWindow[0]->GetBinError(   xBin, yBin);
      inclCFBinErr[1]   = fInclusiveCF2DBackgroundMassWindow[1]->GetBinError(   xBin, yBin);
      bkgCFBinCont[0]   = fCombinatorialBackgroundCF2D[0]->GetBinContent(       xBin, yBin);
      bkgCFBinCont[1]   = fCombinatorialBackgroundCF2D[1]->GetBinContent(       xBin, yBin);
      bkgCFBinErr[0]    = fCombinatorialBackgroundCF2D[0]->GetBinError(         xBin, yBin);
      bkgCFBinErr[1]    = fCombinatorialBackgroundCF2D[1]->GetBinError(         xBin, yBin);

      // calculate background uncertainties
      Double_t errA = ((inclCFBinCont[1]-fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1])/norm -
                       ((1-fTrigValBkg[1][kBkgCombFrac])*(-fTrigValBkg[1][kSigFrac]*inclCFBinCont[0] +
                                                          fTrigValBkg[0][kSigFrac]*inclCFBinCont[1] +
                                                          fTrigValBkg[1][kSigFrac]*fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0] -
                                                          fTrigValBkg[0][kSigFrac]*fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1]))/TMath::Power(norm, 2)) *
                      fTrigValBkg[0][kSigFracErr];
      Double_t errB = ((-inclCFBinCont[0]+fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0])/norm -
                       ((-1+fTrigValBkg[0][kBkgCombFrac])*(-fTrigValBkg[1][kSigFrac]*inclCFBinCont[0] +
                                                           fTrigValBkg[0][kSigFrac]*inclCFBinCont[1] +
                                                           fTrigValBkg[1][kSigFrac]*fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0] -
                                                           fTrigValBkg[0][kSigFrac]*fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1]))/TMath::Power(norm, 2)) *
                      fTrigValBkg[1][kSigFracErr];
      Double_t errC = ((fTrigValBkg[1][kSigFrac]*bkgCFBinCont[0])/norm -
                       (fTrigValBkg[1][kSigFrac]*(-fTrigValBkg[1][kSigFrac]*inclCFBinCont[0] +
                                                  fTrigValBkg[0][kSigFrac]*inclCFBinCont[1] +
                                                  fTrigValBkg[1][kSigFrac]*fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0] -
                                                  fTrigValBkg[0][kSigFrac]*fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1]))/TMath::Power(norm, 2)) *
                      fTrigValBkg[0][kBkgCombFracErr];
      Double_t errD = (-(fTrigValBkg[0][kSigFrac]*bkgCFBinCont[1])/norm +
                       (fTrigValBkg[0][kSigFrac]*(-fTrigValBkg[1][kSigFrac]*inclCFBinCont[0] +
                                                  fTrigValBkg[0][kSigFrac]*inclCFBinCont[1] +
                                                  fTrigValBkg[1][kSigFrac]*fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0] -
                                                  fTrigValBkg[0][kSigFrac]*fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1]))/TMath::Power(norm, 2)) *
                      fTrigValBkg[1][kBkgCombFracErr];
      Double_t errE = fTrigValBkg[1][kSigFrac]/norm * inclCFBinErr[0];
      Double_t errF = fTrigValBkg[0][kSigFrac]/norm * inclCFBinErr[1];
      Double_t errG = fTrigValBkg[1][kSigFrac]*fTrigValBkg[0][kBkgCombFrac]/norm * bkgCFBinErr[0];
      Double_t errH = fTrigValBkg[0][kSigFrac]*fTrigValBkg[1][kBkgCombFrac]/norm * bkgCFBinErr[1];
      fBackgroundCF2D->SetBinError(xBin, yBin, TMath::Sqrt(errA*errA + errB*errB + errC*errC + errD*errD +
                                                           errE*errE + errF*errF + errG*errG + errH*errH));
      
      // calculate signal uncertainties
      Double_t errI = ((inclCFBinCont[1]-fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1])/norm -
                       ((1-fTrigValBkg[1][kBkgCombFrac]) * ((1-fTrigValBkg[1][kSigFrac]-fTrigValBkg[1][kBkgCombFrac]) *
                                                            (inclCFBinCont[0]-fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0]) -
                                                            (1-fTrigValBkg[0][kSigFrac]-fTrigValBkg[0][kBkgCombFrac]) *
                                                            (inclCFBinCont[1]-fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1])))/TMath::Power(norm, 2)) *
                      fTrigValBkg[0][kSigFracErr];
      Double_t errJ = ((-inclCFBinCont[0]+fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0])/norm -
                       ((-1+fTrigValBkg[0][kBkgCombFrac]) * ((1-fTrigValBkg[1][kSigFrac]-fTrigValBkg[1][kBkgCombFrac]) *
                                                             (inclCFBinCont[0]-fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0]) -
                                                             (1-fTrigValBkg[0][kSigFrac]-fTrigValBkg[0][kBkgCombFrac]) *
                                                             (inclCFBinCont[1]-fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1])))/TMath::Power(norm, 2)) *
                      fTrigValBkg[1][kSigFracErr];
      Double_t errK = ((-(1-fTrigValBkg[1][kSigFrac]-fTrigValBkg[1][kBkgCombFrac]) *
                        bkgCFBinCont[0]+inclCFBinCont[1]-fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1])/norm -
                       (fTrigValBkg[1][kSigFrac] * ((1-fTrigValBkg[1][kSigFrac]-fTrigValBkg[1][kBkgCombFrac]) *
                                                    (inclCFBinCont[0]-fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0]) -
                                                    (1-fTrigValBkg[0][kSigFrac]-fTrigValBkg[0][kBkgCombFrac]) *
                                                    (inclCFBinCont[1]-fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1])))/TMath::Power(norm, 2)) *
                      fTrigValBkg[0][kBkgCombFracErr];
      Double_t errL = (((1-fTrigValBkg[0][kSigFrac]-fTrigValBkg[0][kBkgCombFrac]) *
                        bkgCFBinCont[1]-inclCFBinCont[0]+fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0])/norm +
                       (fTrigValBkg[0][kSigFrac] * ((1-fTrigValBkg[1][kSigFrac]-fTrigValBkg[1][kBkgCombFrac]) *
                                                    (inclCFBinCont[0]-fTrigValBkg[0][kBkgCombFrac]*bkgCFBinCont[0]) -
                                                    (1-fTrigValBkg[0][kSigFrac]-fTrigValBkg[0][kBkgCombFrac]) *
                                                    (inclCFBinCont[1]-fTrigValBkg[1][kBkgCombFrac]*bkgCFBinCont[1])))/TMath::Power(norm, 2)) *
                      fTrigValBkg[1][kBkgCombFracErr];
      Double_t errO = (-1+fTrigValBkg[1][kSigFrac]+fTrigValBkg[1][kBkgCombFrac])/norm * inclCFBinErr[0];
      Double_t errP = (-1+fTrigValBkg[0][kSigFrac]+fTrigValBkg[1][kBkgCombFrac])/norm * inclCFBinErr[1];
      Double_t errQ = ((-1+fTrigValBkg[1][kSigFrac]+fTrigValBkg[1][kBkgCombFrac])*fTrigValBkg[0][kBkgCombFrac])/norm * bkgCFBinErr[0];
      Double_t errR = ((-1+fTrigValBkg[0][kSigFrac]+fTrigValBkg[0][kBkgCombFrac])*fTrigValBkg[1][kBkgCombFrac])/norm * bkgCFBinErr[1];
      fSignalCF2D->SetBinError(xBin, yBin, TMath::Sqrt(errI*errI + errJ*errJ + errK*errK + errL*errL +
                                                       errO*errO + errP*errP + errQ*errQ + errR*errR));
    }
  }
  
  // project 1D background correlation from 2D
  fBackgroundCF1D = ProjectToDeltaPhi(fBackgroundCF2D,  "backgroundCF_1D");
  fSignalCF1D     = ProjectToDeltaPhi(fSignalCF2D,      "signalCF_1D");
  
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::CalculateSignalCorrelation() {
  //
  // calculate signal correlation
  //
  if (!fInclusiveCF2D) {
    cout << "AliCorrelationExtraction::CalculateSignalCorrelation() Fatal: Inclusive correlation missing!" << endl;
    return kFALSE;
  }
  if (!fBackgroundCF2D) {
    cout << "AliCorrelationExtraction::CalculateSignalCorrelation() Fatal: Background correlation missing!" << endl;
    return kFALSE;
  }
  if (!fTrigValSig[kSig]) {
    cout << "AliCorrelationExtraction::CalculateSignalCorrelation() Fatal: There was an issue with the J/psi signal extraction!" << endl;
    return kFALSE;
  }
    
  // calculate signal correlation
  fSignalCF2D = (TH2D*)fInclusiveCF2D->Clone("signalCF_2D");
  fSignalCF2D->Scale(fTrigValSig[kSplusB]);
  fSignalCF2D->Add(fBackgroundCF2D, -fTrigValSig[kBkg]);
  fSignalCF2D->Scale(1./fTrigValSig[kSig]);

  // calculate uncertainty
  for (Int_t phiBin=1; phiBin<fSignalCF2D->GetNbinsY()+1; ++phiBin) {
    for (Int_t etaBin=1; etaBin<fSignalCF2D->GetNbinsX()+1; ++etaBin) {
      Double_t inclCFBinCont  = fInclusiveCF2D->GetBinContent(  etaBin, phiBin);
      Double_t inclCFBinErr   = fInclusiveCF2D->GetBinError(    etaBin, phiBin);
      Double_t bkgCFBinCont   = fBackgroundCF2D->GetBinContent( etaBin, phiBin);
      Double_t bkgCFBinErr    = fBackgroundCF2D->GetBinError(   etaBin, phiBin);
      Double_t errA = inclCFBinCont/fTrigValSig[kSig]*fTrigValSig[kSplusBErr];
      Double_t errB = bkgCFBinCont/fTrigValSig[kSig]*fTrigValSig[kBkgErr];
      Double_t errC = (fTrigValSig[kSplusB]*inclCFBinCont-fTrigValSig[kBkg]*bkgCFBinCont)/TMath::Power(fTrigValSig[kSig], 2)*fTrigValSig[kSigErr];
      Double_t errD = fTrigValSig[kSplusB]/fTrigValSig[kSig]*inclCFBinErr;
      Double_t errE = fTrigValSig[kBkg]/fTrigValSig[kSig]*bkgCFBinErr;
      Double_t signalCFBinErr = TMath::Sqrt(errA*errA + errB*errB + errC*errC + errD*errD + errE*errE);
      fSignalCF2D->SetBinError(etaBin, phiBin, signalCFBinErr);
    }
  }
  
  // project 1D signal correlation from 2D
  fSignalCF1D = ProjectToDeltaPhi(fSignalCF2D, "signalCF_1D");
  
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::EfficiencyCorrection() {
  //
  // correct signal correlation for hadron efficiency
  //
  std::unique_ptr<TH1D> normHist;
  if (fSEOS)        normHist = std::unique_ptr<TH1D>(static_cast<TH1D*>(fSEOS->Projection(      fEfficiencyVariableIndex, "e")));
  if (fSEOSSparse)  normHist = std::unique_ptr<TH1D>(static_cast<TH1D*>(fSEOSSparse->Projection(fEfficiencyVariableIndex, "e")));
  if (!normHist.get()) {
    cout << "AliCorrelationExtraction::EfficiencyCorrection() Fatal: There was an issue with the histogram needed for normalization! Efficiency correction can't be applied!" << endl;
    return kFALSE;
  }

  // compare normalization and efficiency bins and find common binning
  std::vector<Double_t> normBins;
  std::vector<Double_t> effBins;
  Double_t* normBinsArr = (Double_t*)(normHist->GetXaxis()->GetXbins())->GetArray();
  Double_t* effBinsArr  = (Double_t*)(fHadronEff->GetXaxis()->GetXbins())->GetArray();
  normBins.assign(normBinsArr, normBinsArr+(normHist->GetNbinsX()+1));
  effBins.assign(effBinsArr, effBinsArr+(fHadronEff->GetNbinsX()+1));
  if (fVerboseFlag) {
    cout << "AliCorrelationExtraction::EfficiencyCorrection() normBins   = {";
    for (Int_t i=0; i<normBins.size(); ++i) cout << normBins.at(i) << ", ";
    cout << "}" << endl;
    cout << "AliCorrelationExtraction::EfficiencyCorrection() effBins    = {";
    for (Int_t i=0; i<effBins.size(); ++i) cout << effBins.at(i) << ", ";
    cout << "}" << endl;
  }
  std::vector<Double_t> commonBins;
  std::set_intersection(normBins.begin(), normBins.end(), effBins.begin(), effBins.end(), std::back_inserter(commonBins));
  if (commonBins.size()<2) {
    cout << "AliCorrelationExtraction::EfficiencyCorrection() Fatal: No common binning in efficiency and normalization histograms found! Efficiency correction can't be applied!" << endl;
    return kFALSE;
  }
  if (fVerboseFlag) {
    cout << "AliCorrelationExtraction::EfficiencyCorrection() commonBins = {";
    for (Int_t i=0; i<commonBins.size(); ++i) cout << commonBins.at(i) << ", ";
    cout << "}" << endl;
  }

  // rebin norm and eff histograms to matching binning
  Double_t* newBins = new Double_t[commonBins.size()];
  for (Int_t i=0; i<commonBins.size(); ++i) newBins[i] = commonBins.at(i);
  std::unique_ptr<TH1D> normHistTmp = std::unique_ptr<TH1D>(static_cast<TH1D*>(normHist->Clone(   "normHistTmp")));
  std::unique_ptr<TH1D> effHistTmp  = std::unique_ptr<TH1D>(static_cast<TH1D*>(fHadronEff->Clone( "effHistTmp")));
  for (Int_t i=1; i<normHistTmp->GetNbinsX()+1; ++i) {
    normHistTmp->SetBinContent( i, normHistTmp->GetXaxis()->GetBinWidth(i)*normHistTmp->GetBinContent(i));
    normHistTmp->SetBinError(   i, normHistTmp->GetXaxis()->GetBinWidth(i)*normHistTmp->GetBinError(  i));
  }
  for (Int_t i=1; i<effHistTmp->GetNbinsX()+1; ++i) {
    effHistTmp->SetBinContent( i, effHistTmp->GetXaxis()->GetBinWidth(i)*effHistTmp->GetBinContent(i));
    effHistTmp->SetBinError(   i, effHistTmp->GetXaxis()->GetBinWidth(i)*effHistTmp->GetBinError(  i));
  }
  normHistTmp = std::unique_ptr<TH1D>(static_cast<TH1D*>(normHistTmp->Rebin(commonBins.size()-1, "normHistTmp", newBins)));
  effHistTmp  = std::unique_ptr<TH1D>(static_cast<TH1D*>(effHistTmp->Rebin( commonBins.size()-1, "effHistTmp",  newBins)));
  normHistTmp->Scale( 1., "width");
  effHistTmp->Scale(  1., "width");
  delete[] newBins;

  // calculate weighted average of efficiency:
  // 1/eff_{tot} = sum_{i=1}^{n} N_{bin i}/N_{tot}*1/eff_{bin i}
  Double_t normTot        = normHistTmp->Integral(1, normHistTmp->GetNbinsX());
  Double_t oneOverEff     = 0.;
  Double_t oneOverEffErr  = 0.;
  for (Int_t i=1; i<effHistTmp->GetNbinsX()+1; ++i) {
    Double_t normBin    = normHistTmp->GetBinContent(i);
    Double_t effBin     = effHistTmp->GetBinContent(i);
    Double_t effBinErr  = effHistTmp->GetBinError(i);

    oneOverEff    += normBin/normTot/effBin;
    oneOverEffErr += TMath::Power(normBin/normTot*effBinErr/TMath::Power(effBin, 2), 2);
  }
  oneOverEffErr   = TMath::Sqrt(oneOverEffErr);
  Double_t eff    = 1./oneOverEff;
  Double_t effErr = oneOverEffErr/TMath::Power(oneOverEff, 2);
  if (!eff || !effErr) {
    cout << "AliCorrelationExtraction::EfficiencyCorrection() Fatal: There was an issue with the efficiency calculation! Efficiency correction can't be applied!" << endl;
    return kFALSE;
  }

  // correct signal histogram (1D and 2D) for efficiency
  fSignalCF1DEffCorr = (TH1D*)fSignalCF1D->Clone("signalCF_1D_efficiencyCorrected");
  for (Int_t i=1; i<fSignalCF1DEffCorr->GetNbinsX()+1; ++i) {
    Double_t binCont    = fSignalCF1DEffCorr->GetBinContent( i);
    Double_t binErr     = fSignalCF1DEffCorr->GetBinError(   i);
    Double_t binContNew = binCont/eff;
    Double_t binErrNew  = TMath::Sqrt(TMath::Power(binErr/eff, 2) + TMath::Power(binCont*effErr/eff/eff, 2));
    fSignalCF1DEffCorr->SetBinContent(i, binContNew);
    fSignalCF1DEffCorr->SetBinError(  i, binErrNew);
  }
  fSignalCF2DEffCorr = (TH2D*)fSignalCF2D->Clone("signalCF_2D_efficiencyCorrected");
  for (Int_t i=1; i<fSignalCF2DEffCorr->GetNbinsX()+1; ++i) {
    for (Int_t j=1; j<fSignalCF2DEffCorr->GetNbinsY()+1; ++j) {
      Double_t binCont    = fSignalCF2DEffCorr->GetBinContent( i, j);
      Double_t binErr     = fSignalCF2DEffCorr->GetBinError(   i, j);
      Double_t binContNew = binCont/eff;
      Double_t binErrNew  = TMath::Sqrt(TMath::Power(binErr/eff, 2) + TMath::Power(binCont*effErr/eff/eff, 2));
      fSignalCF2DEffCorr->SetBinContent(i, j, binContNew);
      fSignalCF2DEffCorr->SetBinError(  i, j, binErrNew);
    }
  }

  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCorrelationExtraction::Process() {
  //
  // main function for correlation extraction
  //
  fProcessDone = kFALSE;
  
  // initialize AliCorrelationExtraction object
  if (!Initialize()) return kFALSE;

  // J/psi signal extraction
  if (fVerboseFlag) fResonanceFits->Print();
  fResonanceFits->Process();

  // print user options
  if (fVerboseFlag) PrintUserOptions();

  // calculate inclusive correlation function in J/psi signal region
  if (!CalculateInclusiveCorrelation(fMassSignalRange[0], fMassSignalRange[1],
                                     fSEOSNorm, fMEOSNorm,
                                     fInclusiveCF2D, fInclusiveCF1D)) return kFALSE;
  fSEOSNorm->SetName("SE-OS_2D");
  fMEOSNorm->SetName("ME-OS_2D");
  fInclusiveCF2D->SetName("inclusiveCF_2D");
  fInclusiveCF1D->SetName("inclusiveCF_1D");

  // calculate background (and signal) correlation according to user selected background method
  switch (fOptionBkgMethod) {
    case kBkgFitting:
      fProcessDone = CalculateBackgroundCorrelationFitting();
      break;
      
    case kBkgSideband:
      fProcessDone = CalculateBackgroundCorrelationSideband();
      break;

    case kBkgLikeSign:
      fProcessDone = CalculateBackgroundCorrelationLikeSign();
      break;

    case kBkgInterpolation:
      fProcessDone = CalculateBackgroundCorrelationInterpolation();
      break;

    case kBkgSuperposition:
      fProcessDone = CalculateBackgroundCorrelationSuperposition();
      break;

    case kBkgSuperpositionTwoComponent:
      fProcessDone = CalculateBackgroundCorrelationSuperpositionTwoComponent();
      break;

    default:
      cout << "AliCorrelationExtraction::Process() Fatal: fOptionBkgMethod " << fOptionBkgMethod << " not recognized!" << endl;
      fProcessDone = kFALSE;
      break;
  }

  // calculate signal correlation
  if (fProcessDone && fOptionBkgMethod!=kBkgSuperposition && fOptionBkgMethod!=kBkgSuperpositionTwoComponent && fOptionBkgMethod!=kBkgFitting) {
    fProcessDone = CalculateSignalCorrelation();
  }

  // correct for hadron efficiency
  if (fHadronEff && fProcessDone) fProcessDone = EfficiencyCorrection();

  return fProcessDone;
}
  

//_______________________________________________________________________________
void AliCorrelationExtraction::PrintUserOptions() {
  //
  // print user options
  //
  cout << endl;
  cout << "AliCorrelationExtraction summary of all user options ======================================" << endl;
  cout << "input histograms: +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if (fSEOS)        cout << "fSEOS        = " << fSEOS << endl;
  if (fSEOSSparse)  cout << "fSEOSSparse  = " << fSEOSSparse << endl;
  if (fMEOS)        cout << "fMEOS        = " << fMEOS << endl;
  if (fMEOSSparse)  cout << "fMEOSSparse  = " << fMEOSSparse << endl;
  if (fMEOS2)       cout << "fMEOS2       = " << fMEOS2 << endl;
  if (fMEOS2Sparse) cout << "fMEOS2Sparse = " << fMEOS2Sparse << endl;
  if (fSEPP)        cout << "fSEPP        = " << fSEPP << endl;
  if (fSEPPSparse)  cout << "fSEPPSparse  = " << fSEPPSparse << endl;
  if (fSEMM)        cout << "fSEMM        = " << fSEMM << endl;
  if (fSEMMSparse)  cout << "fSEMMSparse  = " << fSEMMSparse << endl;
  if (fMEPP)        cout << "fMEPP        = " << fMEPP << endl;
  if (fMEPPSparse)  cout << "fMEPPSparse  = " << fMEPPSparse << endl;
  if (fMEMM)        cout << "fMEMM        = " << fMEMM << endl;
  if (fMEMMSparse)  cout << "fMEMMSparse  = " << fMEMMSparse << endl;
  if (fSEPPPair)    cout << "fSEPPPair    = " << fSEPPPair << endl;
  if (fSEMMPair)    cout << "fSEMMPair    = " << fSEMMPair << endl;
  if (fMEOSPair)    cout << "fMEOSPair    = " << fMEOSPair << endl;
  if (fMEPPPair)    cout << "fMEPPPair    = " << fMEPPPair << endl;
  if (fMEMMPair)    cout << "fMEMMPair    = " << fMEMMPair << endl;
  if (fHadronEff)   cout << "fHadronEff   = " << fHadronEff << endl;
  cout << "variables: ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "fNVariables = " << fNVariables << endl;
  for(Int_t i=0; i<fNVariables; ++i) {
    cout << "idx = " << fVarIndices[i];
    cout << ", name = " << AliReducedVarManager::fgVariableNames[fVariables[i]];
    cout << ", range = [" << fVarLimits[i][0] << ", " << fVarLimits[i][1] << "]" << endl;
  }
  if (fUseMixingVars) {
    cout << "fNMixingVariables = " << fNMixingVariables << endl;
    for(Int_t i=0; i<fNMixingVariables; ++i) {
      cout << "idx = " << fMixingVarIndices[i];
      cout << ", name = " << AliReducedVarManager::fgVariableNames[fMixingVariables[i]];
      cout << ", nBins = " << fNMixingVarBins[i] << endl;
      for (Int_t j=0; j<fNMixingVarBins[i]; ++j) {
        cout << "\tbin " << j << " range = [" << fMixingVarBinLimits[i][j][0] << ", " << fMixingVarBinLimits[i][j][1] << "]" << endl;
      }
    }
  }
  cout << "fMassVariable       = " << AliReducedVarManager::fgVariableNames[fMassVariable] << ",\tindex = " << fMassVariableIndex << endl;
  cout << "fDeltaPhiVariable   = " << AliReducedVarManager::fgVariableNames[fDeltaPhiVariable] << ",\tindex = " << fDeltaPhiVariableIndex << endl;
  cout << "fDeltaEtaVariable   = " << AliReducedVarManager::fgVariableNames[fDeltaEtaVariable] << ",\tindex = " << fDeltaEtaVariableIndex << endl;
  if (fHadronEff)
    cout << "fEfficiencyVariable = " << AliReducedVarManager::fgVariableNames[fEfficiencyVariable] << ",\tindex = " << fEfficiencyVariableIndex << endl;

  if (fMassVariableIndexPair>=0) cout << "LS pair mass variable set: index = " << fMassVariableIndexPair << endl;
  cout << "general user options: +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "fVerboseFlag       = " << fVerboseFlag << endl;
  cout << "fUseMixingVars     = " << fUseMixingVars << endl;
  cout << "fIntegrateDeltaEta = ";
  for (Int_t i=0; i<kNBackgroundMethods; ++i) cout << fIntegrateDeltaEta[i] << " ";
  cout << endl;
  cout << "fResonanceFits     = " << fResonanceFits << " (" << fResonanceFits->GetName() << ")" << endl;
  cout << "fOptionBkgMethod   = " << fOptionBkgMethod << endl;
  if (fBkgFitFunction) cout << "fBkgFitFunction    = " << fBkgFitFunction << " (" << fBkgFitFunction->GetName() << ")" << endl;
  cout << "mass ranges: ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "signal mass range      = [" << fMassSignalRange[0] << ", " << fMassSignalRange[1] << "]" << endl;
  cout << "fNBackgroundMassRanges = " << fNBackgroundMassRanges << endl;
  for (Int_t i=0; i<fNBackgroundMassRanges; ++i) {
    cout << "background mass range " << i << "  = [" << fBackgroundMassRanges[i][0] << ", " << fBackgroundMassRanges[i][1] << "]" << endl;
  }
  cout << "===========================================================================================" << endl;
  cout << endl;
}
