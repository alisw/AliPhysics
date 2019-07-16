/**************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

#include <TTree.h>
#include <TFile.h>
#include "AliLog.h"
#include "AliProdInfo.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliCSAnalysisCutsBase.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMultiInputEventHandler.h"

/// \file AliCSAnalysisCutsBase.cxx
/// \brief Implementation of base analysis cuts class within the correlation studies analysis

TString                             AliCSAnalysisCutsBase::fgPeriodName = "";
AliCSAnalysisCutsBase::ProdPeriods  AliCSAnalysisCutsBase::fgDataPeriod = AliCSAnalysisCutsBase::kNoPeriod;
AliCSAnalysisCutsBase::ProdPeriods  AliCSAnalysisCutsBase::fgAnchorPeriod = AliCSAnalysisCutsBase::kNoPeriod;
AliCSAnalysisCutsBase::EnergyValue  AliCSAnalysisCutsBase::fgEnergy = AliCSAnalysisCutsBase::kUnset;
Bool_t                              AliCSAnalysisCutsBase::fgIsMC = kFALSE;
Bool_t                              AliCSAnalysisCutsBase::fgIsMConlyTruth = kFALSE;
AliInputEventHandler               *AliCSAnalysisCutsBase::fgInputHandler = NULL;
AliMCEventHandler                  *AliCSAnalysisCutsBase::fgMCHandler = NULL;
Bool_t                              AliCSAnalysisCutsBase::fgIsESD = kTRUE;
TClonesArray                       *AliCSAnalysisCutsBase::fgMCArray = NULL;

/// Default constructor for serialization
AliCSAnalysisCutsBase::AliCSAnalysisCutsBase() :
    TNamed(),
    fNParams(0),
    fCutsString(NULL),
    fNCuts(0),
    fQALevel(kQALevelLight),
    fParameters(NULL),
    fCutsEnabledMask(TBits()),
    fCutsActivatedMask(TBits()),
    fDataPeriod(kNoPeriod),
    fHistogramsList(NULL)
{
}

/// Constructor
/// Allocates the needed memory for the number of cuts to support
/// \param nCuts the number of cuts to support
/// \param nParams the number of cuts parameters
/// \param name name of the event cuts
/// \param title title of the event cuts
AliCSAnalysisCutsBase::AliCSAnalysisCutsBase(Int_t nCuts, Int_t nParams, const char *name, const char *title) :
    TNamed(name,title),
    fNParams(nParams),
    fNCuts(nCuts),
    fQALevel(kQALevelLight),
    fCutsEnabledMask(TBits(nCuts)),
    fCutsActivatedMask(TBits(nCuts)),
    fDataPeriod(kNoPeriod),
    fHistogramsList(NULL)
{
  fParameters = new Int_t[nParams];
  fCutsString = new Char_t[nParams+1];
  for (Int_t i=0;i<nParams;i++) {fParameters[i]=0;}
  for (Int_t i=0;i<nParams;i++) {fCutsString[i]='0';} fCutsString[nParams]='\0';
  fCutsEnabledMask.ResetAllBits();
  fCutsActivatedMask.ResetAllBits();
}

/// Destructor
AliCSAnalysisCutsBase::~AliCSAnalysisCutsBase()
{
  if (fParameters != NULL) delete [] fParameters;
  if (fCutsString != NULL) delete [] fCutsString;
}

/// The run to analyze has potentially changed
///
/// Stores the needed production information if required
void AliCSAnalysisCutsBase::NotifyRunGlobal() {

  TString szLHCPeriod = GetPeriodNameFromDataFilePath();

  /* if period has not changed do nothing */
  if (szLHCPeriod.EqualTo(fgPeriodName)) return;

  /* period name has changed */
  AliInfoClass(Form("Data period has changed. New data period: %s", szLHCPeriod.Data()));
  fgPeriodName = szLHCPeriod;
  fgDataPeriod = kNoPeriod;
  fgAnchorPeriod = kNoPeriod;
  fgIsMC = kFALSE;
  fgEnergy = kUnset;

  if (szLHCPeriod.CompareTo("LHC10b") == 0 ||
      szLHCPeriod.CompareTo("LHC10c") == 0 ||
      szLHCPeriod.CompareTo("LHC10d") == 0 ||
      szLHCPeriod.CompareTo("LHC10e") == 0 ||
      szLHCPeriod.CompareTo("LHC10f") == 0 ||
      szLHCPeriod.CompareTo("LHC10g") == 0)
  {
    fgDataPeriod = kLHC10bg;
    fgAnchorPeriod = kLHC10bg;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC10h") == 0) {
    fgDataPeriod = kLHC10h;
    fgAnchorPeriod = kLHC10h;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC11a") == 0) {
    fgDataPeriod = kLHC11a;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC11b") == 0) {
    fgDataPeriod = kLHC11b;
    fgAnchorPeriod = kLHC11b;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC11c") == 0 ||
      szLHCPeriod.CompareTo("LHC11d") == 0 ||
      szLHCPeriod.CompareTo("LHC11e") == 0 ||
      szLHCPeriod.CompareTo("LHC11f") == 0 ||
      szLHCPeriod.CompareTo("LHC11g") == 0)
  {
    fgDataPeriod = kLHC11cg;
    fgAnchorPeriod = kLHC11cg;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC11h") == 0) {
    fgDataPeriod = kLHC11h;
    fgAnchorPeriod = kLHC11h;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC12a") == 0 ||
      szLHCPeriod.CompareTo("LHC12b") == 0 ||
      szLHCPeriod.CompareTo("LHC12c") == 0 ||
      szLHCPeriod.CompareTo("LHC12d") == 0 ||
      szLHCPeriod.CompareTo("LHC12e") == 0 ||
      szLHCPeriod.CompareTo("LHC12f") == 0 ||
      szLHCPeriod.CompareTo("LHC12g") == 0 ||
      szLHCPeriod.CompareTo("LHC12h") == 0 ||
      szLHCPeriod.CompareTo("LHC12i") == 0)
  {
    fgDataPeriod = kLHC12;
    fgAnchorPeriod = kLHC12;
    fgEnergy = k8TeV;
  } else if (szLHCPeriod.CompareTo("LHC13b") == 0 ||
      szLHCPeriod.CompareTo("LHC13c") == 0 ){
    fgDataPeriod = kLHC13bc;
    fgAnchorPeriod = kLHC13bc;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC13d") == 0 ||
      szLHCPeriod.CompareTo("LHC13e") == 0 ){
    fgDataPeriod = kLHC13de;
    fgAnchorPeriod = kLHC13de;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC13f") == 0 ){
    fgDataPeriod = kLHC13f;
    fgAnchorPeriod = kLHC13f;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC13g") == 0 ){
    fgDataPeriod = kLHC13g;
    fgAnchorPeriod = kLHC13g;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC15f") == 0 ||
      szLHCPeriod.CompareTo("LHC15g") == 0 ||
      szLHCPeriod.CompareTo("LHC15h") == 0 ||
      szLHCPeriod.CompareTo("LHC15i") == 0 ||
      szLHCPeriod.CompareTo("LHC15j") == 0 ||
      szLHCPeriod.CompareTo("LHC15k") == 0 ||
      szLHCPeriod.CompareTo("LHC15l") == 0 ||
      szLHCPeriod.CompareTo("LHC15m") == 0) {
    fgDataPeriod = kLHC15fm;
    fgAnchorPeriod = kLHC15fm;
    fgEnergy = k13TeV;
  } else if (szLHCPeriod.CompareTo("LHC15n") == 0 ){
    fgDataPeriod = kLHC15n;
    fgAnchorPeriod = kLHC15n;
    fgEnergy = k5TeV;
  } else if (szLHCPeriod.CompareTo("LHC15o") == 0 ){
    /* we need to check if LIR or HIR                     */
    /* we do that by checking against the LIR run numbers */
    switch (GetCurrentRunNumber()) {
    case 244918:
    case 244975:
    case 244980:
    case 244982:
    case 244983:
    case 245061:
    case 245064:
    case 245066:
    case 245068:
    case 246390:
    case 246391:
    case 246392:
      fgDataPeriod = kLHC15oLIR;
      fgAnchorPeriod = kLHC15oLIR;
      break;
    default:
      fgDataPeriod = kLHC15oHIR;
      fgAnchorPeriod = kLHC15oHIR;
    }
    fgEnergy = kPbPb5TeV;
  // LHC16 pp periods
  } else if (szLHCPeriod.CompareTo("LHC16k") == 0) {
    fgDataPeriod = kLHC16k;
    fgAnchorPeriod = kLHC16k;
    fgEnergy = k13TeV;

  } else if (szLHCPeriod.CompareTo("LHC16l") == 0) {
    fgDataPeriod = kLHC16l;
    fgAnchorPeriod = kLHC16l;
    fgEnergy = k13TeV;

  // LHC17 XeXe periods
  } else if (szLHCPeriod.CompareTo("LHC17n") == 0) {
    fgDataPeriod = kLHC17n;
    fgAnchorPeriod = kLHC17n;
    fgEnergy = kXeXe5440GeV;

  // LHC10x anchored MCs
  } else if (szLHCPeriod.CompareTo("LHC10d1") == 0){
    fgDataPeriod = kLHC10d1;
    fgAnchorPeriod = kLHC10bg;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC10d2") == 0){
    fgDataPeriod = kLHC10d2;
    fgAnchorPeriod = kLHC10bg;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC10d4a") == 0){
    fgDataPeriod = kLHC10d4a;
    fgIsMC = kTRUE;
    fgAnchorPeriod = kLHC10bg;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC10d4") == 0){
    fgDataPeriod = kLHC10d4;
    fgAnchorPeriod = kLHC10bg;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC10e12") == 0){
    fgDataPeriod = kLHC10e12;
    fgAnchorPeriod = kLHC10bg;
    fgIsMC = kTRUE;
    fgEnergy = k900GeV;
  } else if (szLHCPeriod.CompareTo("LHC10e13") == 0){
    fgDataPeriod = kLHC10e13;
    fgAnchorPeriod = kLHC10bg;
    fgIsMC = kTRUE;
    fgEnergy = k900GeV;
  } else if (szLHCPeriod.CompareTo("LHC10e20") == 0){
    fgDataPeriod = kLHC10e20;
    fgAnchorPeriod = kLHC10bg;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC10e21") == 0){
    fgDataPeriod = kLHC10e21;
    fgAnchorPeriod = kLHC10bg;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC10f6a") == 0){
    fgDataPeriod = kLHC10f6a;
    fgAnchorPeriod = kLHC10bg;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC10f6") == 0){
    fgDataPeriod = kLHC10f6;
    fgAnchorPeriod = kLHC10bg;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.Contains("LHC14j4")){
    fgDataPeriod = kLHC14j4;
    fgAnchorPeriod = kLHC10bg;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC11a10a") == 0){
    fgDataPeriod = kLHC11a10a;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC11a10b") == 0){
    fgDataPeriod = kLHC11a10b;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC13d2") == 0){
    fgDataPeriod = kLHC13d2;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC13d2b") == 0){
    fgDataPeriod = kLHC13d2b;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC12a11a") == 0){
    fgDataPeriod = kLHC12a11a;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC12a11b") == 0){
    fgDataPeriod = kLHC12a11b;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC12a11c") == 0){
    fgDataPeriod = kLHC12a11c;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC12a11d") == 0){
    fgDataPeriod = kLHC12a11d;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC12a11e") == 0){
    fgDataPeriod = kLHC12a11e;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC12a11f") == 0){
    fgDataPeriod = kLHC12a11f;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  // LHC11x anchored MCs
  } else if (szLHCPeriod.CompareTo("LHC12a15c") == 0){
    fgDataPeriod = kLHC12a15c;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.Contains("LHC12f1a") ){
    fgDataPeriod = kLHC12f1a;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.Contains("LHC12f1b") ){
    fgDataPeriod = kLHC12f1b;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.Contains("LHC12i3") ){
    fgDataPeriod = kLHC12i3;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC15g1a") == 0){
    fgDataPeriod = kLHC15g1a;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC15g1b") == 0){
    fgDataPeriod = kLHC15g1b;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC13e4") == 0){
    fgDataPeriod = kLHC13e4;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC13e5") == 0){
    fgDataPeriod = kLHC13e5;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC14k1a") == 0){
    fgDataPeriod = kLHC14k1a;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC14k1b") == 0){
    fgDataPeriod = kLHC14k1b;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC12a15f") == 0){
    fgDataPeriod = kLHC12a15f;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC12a15g") == 0){
    fgDataPeriod = kLHC12a15g;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC12f2a") == 0){
    fgDataPeriod = kLHC12f2a;
    fgIsMC = kTRUE;
    fgEnergy = k7TeV;
  } else if (szLHCPeriod.CompareTo("LHC14a1a") == 0){
    fgDataPeriod = kLHC14a1a;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC14a1b") == 0){
    fgDataPeriod = kLHC14a1b;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC14a1c") == 0){
    fgDataPeriod = kLHC14a1c;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb2760GeV;
  // LHC12x anchored MCs
  } else if (szLHCPeriod.CompareTo("LHC14e2a") == 0){
    fgDataPeriod = kLHC14e2a;
    fgIsMC = kTRUE;
    fgEnergy = k8TeV;
  } else if (szLHCPeriod.CompareTo("LHC14e2b") == 0){
    fgDataPeriod = kLHC14e2b;
    fgIsMC = kTRUE;
    fgEnergy = k8TeV;
  } else if (szLHCPeriod.CompareTo("LHC14e2c") == 0){
    fgDataPeriod = kLHC14e2c;
    fgIsMC = kTRUE;
    fgEnergy = k8TeV;
  } else if (szLHCPeriod.Contains("LHC15h1")){
    fgDataPeriod = kLHC15h1;
    fgIsMC = kTRUE;
    fgEnergy = k8TeV;
  } else if (szLHCPeriod.Contains("LHC15h2")){
    fgDataPeriod = kLHC15h2;
    fgIsMC = kTRUE;
    fgEnergy = k8TeV;
  } else if (szLHCPeriod.CompareTo("LHC16c2") == 0){
    fgDataPeriod = kLHC16c2;
    fgIsMC = kTRUE;
    fgEnergy = k8TeV;
  // LHC13x anchored MCs
  } else if (szLHCPeriod.Contains("LHC13b2_efix")){
    fgDataPeriod = kLHC13b2_efix;
    fgIsMC = kTRUE;
    fgAnchorPeriod = kLHC13bc;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.Contains("LHC13b2")){
    fgDataPeriod = kLHC13b2;
    fgIsMC = kTRUE;
    fgAnchorPeriod = kLHC13bc;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC13e7") == 0){
    fgDataPeriod = kLHC13e7;
    fgIsMC = kTRUE;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC14b2") == 0){
    fgDataPeriod = kLHC14b2;
    fgIsMC = kTRUE;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC13b4_fix") == 0){
    fgDataPeriod = kLHC13b4_fix;
    fgIsMC = kTRUE;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC13b4_plus") == 0){
    fgDataPeriod = kLHC13b4_plus;
    fgIsMC = kTRUE;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC16c3a") == 0){
    fgDataPeriod = kLHC16c3a;
    fgIsMC = kTRUE;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC16c3b") == 0){
    fgDataPeriod = kLHC16c3b;
    fgIsMC = kTRUE;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC16c3c") == 0){
    fgDataPeriod = kLHC16c3c;
    fgIsMC = kTRUE;
    fgEnergy = kpPb5TeV;
  } else if (szLHCPeriod.CompareTo("LHC15g2") == 0){
    fgDataPeriod = kLHC15g2;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC15a3a") == 0){
    fgDataPeriod = kLHC15a3a;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC15a3a_plus") == 0){
    fgDataPeriod = kLHC15a3a_plus;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC15a3b") == 0){
    fgDataPeriod = kLHC15a3b;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC15d3a") == 0){
    fgDataPeriod = kLHC15d3a;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  } else if (szLHCPeriod.CompareTo("LHC15d3b") == 0){
    fgDataPeriod = kLHC15d3b;
    fgIsMC = kTRUE;
    fgEnergy = k2760GeV;
  // LHC15x anchored MCs
  } else if (szLHCPeriod.CompareTo("LHC15g3a3") == 0){
    fgDataPeriod = kLHC15g3a3;
    fgIsMC = kTRUE;
    fgEnergy = k13TeV;
  } else if (szLHCPeriod.CompareTo("LHC15g3a") == 0){
    fgDataPeriod = kLHC15g3a;
    fgIsMC = kTRUE;
    fgEnergy = k13TeV;
  } else if (szLHCPeriod.CompareTo("LHC15g3c2") == 0){
    fgDataPeriod = kLHC15g3c2;
    fgIsMC = kTRUE;
    fgEnergy = k13TeV;
  } else if (szLHCPeriod.CompareTo("LHC15g3c3") == 0){
    fgDataPeriod = kLHC15g3c3;
    fgIsMC = kTRUE;
    fgEnergy = k13TeV;
  } else if (szLHCPeriod.CompareTo("LHC15g3") == 0){
    fgDataPeriod = kLHC15g3;
    fgIsMC = kTRUE;
    fgEnergy = k13TeV;
  } else if (szLHCPeriod.CompareTo("LHC16a2a") == 0){
    fgDataPeriod = kLHC16a2a;
    fgIsMC = kTRUE;
    fgEnergy = k13TeV;
  } else if (szLHCPeriod.CompareTo("LHC16a2b") == 0){
    fgDataPeriod = kLHC16a2b;
    fgIsMC = kTRUE;
    fgEnergy = k13TeV;
  } else if (szLHCPeriod.CompareTo("LHC16a2c") == 0){
    fgDataPeriod = kLHC16a2c;
    fgIsMC = kTRUE;
    fgEnergy = k13TeV;
  } else if (szLHCPeriod.CompareTo("LHC15l1a2") == 0){
    fgDataPeriod = kLHC15l1a2;
    fgIsMC = kTRUE;
    fgEnergy = k5TeV;
  } else if (szLHCPeriod.CompareTo("LHC15l1b2") == 0){
    fgDataPeriod = kLHC15l1b2;
    fgIsMC = kTRUE;
    fgEnergy = k5TeV;
  } else if (szLHCPeriod.Contains("LHC15k1a1")){
    fgDataPeriod = kLHC15k1a1;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC15k1a2")){
    fgDataPeriod = kLHC15k1a2;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC15k1a3")){
    fgDataPeriod = kLHC15k1a3;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC16h4a")){
    fgDataPeriod = kLHC16h4a;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC16h4b")){
    fgDataPeriod = kLHC16h4b;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC16h4b2")){
    fgDataPeriod = kLHC16h4b2;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC16h4c")){
    fgDataPeriod = kLHC16h4c;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  // LHC15x anchored MCs
  } else if (szLHCPeriod.Contains("LHC16g1")){
    fgDataPeriod = kLHC16g1;
    fgAnchorPeriod = kLHC15oHIR;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC16g1a")){
    fgDataPeriod = kLHC16g1a;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC16g1b")){
    fgDataPeriod = kLHC16g1b;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC16g1c")){
    fgDataPeriod = kLHC16g1c;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC16h2a")){
    fgDataPeriod = kLHC16h2a;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC16h2b")){
    fgDataPeriod = kLHC16h2b;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("LHC16h2c")){
    fgDataPeriod = kLHC16h2c;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("kLHC16i3")){
    fgDataPeriod = kLHC16i3;
    fgAnchorPeriod = kLHC15oHIR;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;
  } else if (szLHCPeriod.Contains("kLHC17i2")){
    fgDataPeriod = kLHC17i2;
    fgAnchorPeriod = kLHC15oHIR;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;

  // LHC17x anchored MCs
  } else if (szLHCPeriod.Contains("LHC17j6")) {
    fgDataPeriod = kLHC17j6;
    fgAnchorPeriod = kLHC17n;
    fgIsMC = kTRUE;
    fgEnergy = kXeXe5440GeV;
  } else if (szLHCPeriod.Contains("LHC17j7")) {
    fgDataPeriod = kLHC17j7;
    fgAnchorPeriod = kLHC17n;
    fgIsMC = kTRUE;
    fgEnergy = kXeXe5440GeV;

  // MC upgrade
  } else if (szLHCPeriod.Contains("LHC13d19")){
    fgDataPeriod = kLHC13d19;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;

  // LHC18x PbPb productions
  } else if (szLHCPeriod.Contains("LHC18q")){
    fgDataPeriod = kLHC18q;
    fgAnchorPeriod = kLHC18q;
    fgIsMC = kFALSE;
    fgEnergy = kPbPb5TeV;

  } else if (szLHCPeriod.Contains("LHC18r")){
    fgDataPeriod = kLHC18r;
    fgAnchorPeriod = kLHC18r;
    fgIsMC = kFALSE;
    fgEnergy = kPbPb5TeV;

  // LHC18x anchored MCs
  } else if (szLHCPeriod.Contains("LHC18l8")) {
    fgDataPeriod = kLHC18l8;
    fgAnchorPeriod = kLHC18r;
    fgIsMC = kTRUE;
    fgEnergy = kPbPb5TeV;

  // fast MC productions
  } else if (szLHCPeriod.Contains("LHC13f3")) {
    fgDataPeriod = kLHC13f3;
    fgAnchorPeriod = kLHC10h;
    fgIsMC = kTRUE;
    fgIsMConlyTruth = kTRUE;
    fgEnergy = kPbPb2760GeV;

  } else {
    AliFatalClass(Form("Analysis period %s not supported. Please update the class!!!", szLHCPeriod.Data()));
    fgDataPeriod = kNoPeriod;
    fgEnergy = kUnset;
  }

  /* let's check the consistency of the MC flag and store the input and the MC handlers */
  fgInputHandler = (AliInputEventHandler *)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if (fgInputHandler != NULL && fgInputHandler->InheritsFrom("AliESDInputHandler")) {
    fgIsESD = kTRUE;
    AliMultiInputEventHandler *multiInputHandler = dynamic_cast<AliMultiInputEventHandler *>(fgInputHandler);
    if (multiInputHandler) {
      fgInputHandler = dynamic_cast<AliInputEventHandler *>(multiInputHandler->GetFirstInputEventHandler());
      fgMCHandler = dynamic_cast<AliMCEventHandler *>(multiInputHandler->GetFirstMCEventHandler());
    } else {
      fgMCHandler = dynamic_cast<AliMCEventHandler *>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    }

    if (fgIsMC && !(fgMCHandler != NULL))
      AliFatalClass("Data is from a MC production but there is no MC handler for it");
  }
  else {
    /* we know data format is AOD but still don't know whether it is MC or not */
    fgIsESD = kFALSE;
  }
}

/// Extract the period name from the data file path
/// \return the name of the period associated with the input data file
TString AliCSAnalysisCutsBase::GetPeriodNameFromDataFilePath()
{
  TString szPeriodName = "";
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();

  if(man) {
    AliVEventHandler* inputHandler = (AliVEventHandler*) (man->GetInputEventHandler());

    if (inputHandler){
      TString szFullFileName = inputHandler->GetTree()->GetCurrentFile()->GetName();

      AliInfoClass(Form("Detecting period name from filename: %s", szFullFileName.Data()));

      TObjArray *tokens = szFullFileName.Tokenize("/._");
      for (Int_t i = 0; i < tokens->GetEntriesFast();i++ ){
        TObjString* testObjString = (TObjString*)tokens->At(i);
        if (testObjString->GetString().BeginsWith("LHC")){
          szPeriodName = testObjString->GetString();
          break;
        }
      }
      delete tokens;
      if (szPeriodName.Length() == 0){
        tokens = szFullFileName.Tokenize("__");
        for (Int_t i = 0; i < tokens->GetEntriesFast();i++ ){
          TObjString* testObjString = (TObjString*) tokens->At(i);
          if (testObjString->GetString().BeginsWith("LHC")){
            szPeriodName = testObjString->GetString();
            break;
          }
        }
        delete tokens;
      }
    }
  }
  return szPeriodName;
}

/// Get the current run number being (or going to be) analyzed
/// \return the current run number
Int_t AliCSAnalysisCutsBase::GetCurrentRunNumber()
{
  Int_t runno = -1;
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();

  if(man != NULL) {
    AliVEventHandler* inputHandler = (AliVEventHandler*) (man->GetInputEventHandler());

    if (inputHandler != NULL){
      AliVEvent *event = inputHandler->GetEvent();

      if (event != NULL) {
        runno = event->GetRunNumber();
      }
    }
  }
  return runno;
}


/// Update the cut string (if it has been created yet)
/// \return kFALSE if inconsistent length or cut values
Bool_t AliCSAnalysisCutsBase::UpdateCutsString() {

  if (Int_t(strlen(fCutsString)) == fNParams) {
    /* before touching anything check consistency */
    for (Int_t i = 0; i < fNParams; i++) {
      if (fParameters[i] < 0 || 9 < fParameters[i]) {
        return kFALSE;
      }
    }
  } else {
    /* not consistent length */
    return kFALSE;
  }
  /* now it is safe to do it */
  for (Int_t i = 0; i < fNParams; i++) {
    fCutsString[i] = '0' + fParameters[i];
  }
  fCutsString[fNParams] = '\0';

  return kTRUE;
}

/// Gives values to the different event cuts from the passed string
/// \param cutsSelectionString string with the values for the different cuts
/// \return kTRUE if succeeded

Bool_t AliCSAnalysisCutsBase::InitializeCutsFromCutString(const TString cutsSelectionString ) {
  AliInfo(Form("Set Cuts String: %s",cutsSelectionString.Data()));
  if(cutsSelectionString.Length()!=fNParams) {
    AliError(Form("Cuts selection string has the wrong length! size is %d, number of cuts is %d", cutsSelectionString.Length(), fNParams));
    return kFALSE;
  }
  if(!cutsSelectionString.IsDigit()){
    AliError("Cut selection string contains non digits");
    return kFALSE;
  }

  // Set Individual Cuts
  for(Int_t i=0;i<fNParams;i++){
    if(!SetCutAndParams(i,cutsSelectionString[i] - '0'))return kFALSE;
  }
  /* produce some feedback */
  PrintCutsWithValues();
  return kTRUE;
}

/// Prints the cuts information and values in an usere friendly way
/// Ask to each cut to print its own information
void AliCSAnalysisCutsBase::PrintCutsWithValues() const {
  // Print out current Cut Selection with value
  printf("\n=========== %s information ===============\n", GetName());
  printf("Cuts string: %s\n",fCutsString);
  printf("Cuts values: ");
  for(Int_t i = 0; i < fNParams; i++) {
    printf("%d",fParameters[i]);
  }
  printf("\nIndividual cut information\n");

  for (Int_t i = 0; i < fNParams; i++) {
    PrintCutWithParams(i);
  }
  printf("=========== %s information end ===========\n\n", GetName());
}

/// \cond CLASSIMP
ClassImp(AliCSAnalysisCutsBase);
/// \endcond
