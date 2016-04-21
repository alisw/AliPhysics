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

/*************************************************
 * Charge-rapidity correlations with Q-cumulants *
 *                                               *
 * author: Jacopo Margutti                       *
 *         (margutti@nikhef.nl)                  *
 *************************************************/

#ifndef ALIFLOWANALYSISCRC_H
#define ALIFLOWANALYSISCRC_H

#include "TMatrixD.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "AliFlowCommonConstants.h"
#include "TNamed.h"

class TObjArray;
class TList;
class TFile;
class TGraph;
class TH1;
class TProfile;
class TProfile2D;
class TDirectoryFile;
class TRandom3;
class TNtuple;
class THnSparse;

class AliFlowEventSimple;
class AliFlowTrackSimple;
class AliFlowCommonConstants;
class AliFlowCommonHist;
class AliFlowCommonHistResults;
class AliFlowVector;

//==============================================================================================================

class AliFlowAnalysisCRC : public TNamed {
public:
 AliFlowAnalysisCRC(const char* name="AliFlowAnalysisCRC",
                    Int_t nCen=7,
                    Double_t CenWidth=10.);
 virtual ~AliFlowAnalysisCRC();
 
 enum CorrelationWeights { kMultiplicity,
  kUnit
 };
 
 enum DataSet { k2010,
  k2011,
  kAny
 };
 
 enum MagnetPol { kMAll,
  kMPos,
  kMNeg
 };
 
 // 0.) methods called in the constructor:
 virtual void InitializeArraysForIntFlow();
 virtual void InitializeArraysForDiffFlow();
 virtual void InitializeArraysForDistributions();
 virtual void InitializeArraysForVarious();
 virtual void InitializeArraysForParticleWeights();
 virtual void InitializeArraysForNestedLoops();
 virtual void InitializeArraysForMixedHarmonics();
 virtual void InitializeArraysForControlHistograms();
 virtual void InitializeArraysForBootstrap();
 virtual void InitializeCostantsForCRC();
 virtual void InitializeArraysForCRC();
 virtual void InitializeArraysForCRCVZ();
 virtual void InitializeArraysForCRCZDC();
 virtual void InitializeArraysForCRCPt();
 virtual void InitializeArraysForCME();
 virtual void InitializeArraysForCRC2();
 virtual void InitializeArraysForFlowEbE();
 virtual void InitializeArraysForFlowQC();
 virtual void InitializeArraysForFlowSPZDC();
 virtual void InitializeArraysForFlowSPVZ();
 
 // 1.) method Init() and methods called within Init():
 virtual void Init();
 virtual void InitializeArraysForQVec();
 virtual void CrossCheckSettings();
 virtual void CommonConstants(TString method);
 virtual void BookAndNestAllLists();
 virtual void BookAndNestListsForDifferentialFlow();
 virtual void BookCommonHistograms();
 virtual void BookAndFillWeightsHistograms();
 virtual void BookEverythingForIntegratedFlow();
 virtual void BookEverythingForDifferentialFlow();
 virtual void BookEverythingFor2DDifferentialFlow();
 virtual void BookEverythingForDistributions();
 virtual void BookEverythingForVarious();
 virtual void BookEverythingForNestedLoops();
 virtual void BookEverythingForMixedHarmonics();
 virtual void BookEverythingForControlHistograms();
 virtual void BookEverythingForBootstrap();
 virtual void SetRunList();
 virtual void BookEverythingForCRC();
 virtual void BookEverythingForCRCVZ();
 virtual void BookEverythingForCRCZDC();
 virtual void BookEverythingForCRCPt();
 virtual void BookEverythingForQVec();
 virtual void BookEverythingForCME();
 virtual void BookEverythingForCRC2();
 virtual void BookEverythingForFlowEbE();
 virtual void BookEverythingForFlowQC();
 virtual void BookEverythingForFlowSPZDC();
 virtual void BookEverythingForFlowSPVZ();
 virtual void StoreIntFlowFlags();
 virtual void StoreDiffFlowFlags();
 virtual void StoreFlagsForDistributions();
 virtual void StoreHarmonic();
 virtual void StoreMixedHarmonicsFlags();
 virtual void StoreControlHistogramsFlags();
 virtual void StoreBootstrapFlags();
 virtual void StoreCRCFlags();
 virtual void SetCentralityWeights();
 
 // 2.) method Make() and methods called within Make():
 virtual void Make(AliFlowEventSimple *anEvent);
 // 2a.) Common:
 virtual void CheckPointersUsedInMake();
 virtual void FillAverageMultiplicities(Int_t nRP);
 virtual void FillCommonControlHistograms(AliFlowEventSimple *anEvent);
 virtual void FillControlHistograms(AliFlowEventSimple *anEvent);
 virtual void ResetEventByEventQuantities();
 // 2b.) Reference flow:
 virtual void CalculateIntFlowCorrelations();
 virtual void CalculateIntFlowCorrelationsUsingParticleWeights();
 virtual void CalculateIntFlowProductOfCorrelations();
 virtual void CalculateIntFlowSumOfEventWeights();
 virtual void CalculateIntFlowSumOfProductOfEventWeights();
 virtual void CalculateIntFlowCorrectionsForNUACosTerms();
 virtual void CalculateIntFlowCorrectionsForNUACosTermsUsingParticleWeights();
 virtual void CalculateIntFlowCorrectionsForNUASinTerms();
 virtual void CalculateIntFlowCorrectionsForNUASinTermsUsingParticleWeights();
 virtual void CalculateIntFlowProductOfCorrectionTermsForNUA();
 virtual void CalculateIntFlowSumOfEventWeightsNUA();
 virtual void CalculateIntFlowSumOfProductOfEventWeightsNUA();
 virtual void CalculateMixedHarmonics();
 // 2c.) Cross-checking reference flow correlations with nested loops:
 virtual void EvaluateIntFlowNestedLoops(AliFlowEventSimple* const anEvent);
 virtual void EvaluateIntFlowCorrelationsWithNestedLoops(AliFlowEventSimple* const anEvent);
 virtual void EvaluateIntFlowCorrelationsWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* const anEvent);
 virtual void EvaluateIntFlowCorrectionsForNUAWithNestedLoops(AliFlowEventSimple* const anEvent);
 virtual void EvaluateIntFlowCorrectionsForNUAWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* const anEvent);
 virtual void EvaluateMixedHarmonicsWithNestedLoops(AliFlowEventSimple* const anEvent);
 // 2d.) Differential flow:
 virtual void CalculateDiffFlowCorrelations(TString type, TString ptOrEta); // type = RP or POI
 virtual void CalculateDiffFlowCorrelationsUsingParticleWeights(TString type, TString ptOrEta); // type = RP or POI
 virtual void CalculateDiffFlowProductOfCorrelations(TString type, TString ptOrEta); // type = RP or POI
 virtual void CalculateDiffFlowSumOfEventWeights(TString type, TString ptOrEta); // type = RP or POI
 virtual void CalculateDiffFlowSumOfProductOfEventWeights(TString type, TString ptOrEta); // type = RP or POI
 virtual void CalculateDiffFlowCorrectionsForNUACosTerms(TString type, TString ptOrEta);
 virtual void CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights(TString type, TString ptOrEta);
 virtual void CalculateDiffFlowCorrectionsForNUASinTerms(TString type, TString ptOrEta);
 virtual void CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights(TString type, TString ptOrEta);
 // 2e.) 2D differential flow:
 virtual void Calculate2DDiffFlowCorrelations(TString type); // type = RP or POI
 // 2f.) Other differential correlators (i.e. Teaney-Yan correlator):
 virtual void CalculateOtherDiffCorrelators(TString type, TString ptOrEta); // type = RP or POI
 // 2g.) Distributions of reference flow correlations:
 virtual void StoreDistributionsOfCorrelations();
 // 2h.) Cross-checking differential flow correlations with nested loops:
 virtual void EvaluateDiffFlowNestedLoops(AliFlowEventSimple* const anEvent);
 virtual void EvaluateDiffFlowCorrelationsWithNestedLoops(AliFlowEventSimple* const anEvent, TString type, TString ptOrEta);
 virtual void EvaluateDiffFlowCorrelationsWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* const anEvent, TString type, TString ptOrEta);
 virtual void EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(AliFlowEventSimple* const anEvent, TString type, TString ptOrEta);
 virtual void EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* const anEvent, TString type, TString ptOrEta);
 virtual void EvaluateOtherDiffCorrelatorsWithNestedLoops(AliFlowEventSimple* const anEvent, TString type, TString ptOrEta);
 // 2i.) Charge-Rapidity Correlations
 virtual void RecenterCRCQVec();
 virtual void RecenterCRCQVecZDC();
 virtual Bool_t PassQAZDCCuts();
 virtual void CalculateCRCCorr();
 virtual void CalculateCRCVZERO();
 virtual void CalculateCRCZDC();
 virtual void CalculateCRCPtCorr();
 virtual void CalculateCRCQVec();
 virtual void CalculateVZvsZDC();
 virtual void CalculateCMETPC();
 virtual void CalculateCMEZDC();
 virtual void CalculateCRC2Cor();
 virtual void CalculateFlowQC();
 virtual void CalculateFlowSPZDC();
 virtual void CalculateFlowSPVZ();
 // 2h.) Various
 virtual void FillVarious();
 
 // 3.) method Finish() and methods called within Finish():
 virtual void Finish();
 virtual void CheckPointersUsedInFinish();
 // 3a.) integrated flow:
 virtual void FinalizeCorrelationsIntFlow();
 virtual void FinalizeCorrectionTermsForNUAIntFlow();
 virtual void CalculateCovariancesIntFlow();
 virtual void CalculateCovariancesNUAIntFlow();
 virtual void CalculateCumulantsIntFlow();
 virtual void CalculateReferenceFlow();
 virtual void FillCommonHistResultsIntFlow();
 //  nua:
 virtual void CalculateQcumulantsCorrectedForNUAIntFlow();
 virtual void PrintFinalResultsForIntegratedFlow(TString type);
 virtual void CrossCheckIntFlowCorrelations();
 virtual void CrossCheckIntFlowExtraCorrelations(); // extra correlations which appear only when particle weights are used
 virtual void CrossCheckIntFlowCorrectionTermsForNUA();
 // 3b.) differential flow:
 virtual void FinalizeReducedCorrelations(TString type, TString ptOrEta);
 virtual void CalculateDiffFlowCovariances(TString type, TString ptOrEta);
 virtual void CalculateDiffFlowCumulants(TString type, TString ptOrEta);
 virtual void CalculateDiffFlow(TString type, TString ptOrEta);
 virtual void FinalizeCorrectionTermsForNUADiffFlow(TString type, TString ptOrEta);
 virtual void CalculateDiffFlowCumulantsCorrectedForNUA(TString type, TString ptOrEta);
 virtual void CalculateDiffFlowCorrectedForNUA(TString type, TString ptOrEta);
 virtual void CalculateFinalResultsForRPandPOIIntegratedFlow(TString type); // to be improved (add also possibility to integrate over eta yield)
 virtual void FillCommonHistResultsDiffFlow(TString type);
 virtual void CrossCheckDiffFlowCorrelations(TString type, TString ptOrEta);
 virtual void PrintNumberOfParticlesInSelectedBin();
 virtual void CrossCheckDiffFlowCorrectionTermsForNUA(TString type, TString ptOrEta);
 // 3c.) 2D:
 virtual void Calculate2DDiffFlowCumulants(TString type);
 virtual void Calculate2DDiffFlow(TString type);
 // 3d.) Other differential correlators:
 virtual void CrossCheckOtherDiffCorrelators(TString type, TString ptOrEta);
 // 3e.) Mixed harmonics:
 virtual void CalculateCumulantsMixedHarmonics();
 // 3f.) Bootstrap:
 virtual void CalculateCumulantsForBootstrap();
 // 3g.) CRC:
 virtual void FinalizeCRCCorr();
 virtual void FinalizeCRCVZERO();
 virtual void FinalizeCRCZDC();
 virtual void FinalizeCRCPtCorr();
 virtual void FinalizeCMETPC();
 virtual void FinalizeCMEZDC();
 virtual void FinalizeCRC2Cor();
 virtual void FinalizeFlowQC();
 virtual void FinalizeFlowSPZDC();
 virtual void FinalizeFlowSPVZ();
 virtual Bool_t CheckRunFullTPCFlow(Int_t RunNum);
 // 3h.) Various:
 virtual void FinalizeVarious();
 
 // 4.)  method GetOutputHistograms() and methods called within GetOutputHistograms():
 virtual void GetOutputHistograms(TList *outputListHistos);
 virtual void GetPointersForCommonHistograms();
 virtual void GetPointersForParticleWeightsHistograms();
 virtual void GetPointersForIntFlowHistograms();
 virtual void GetPointersForDiffFlowHistograms();
 virtual void GetPointersFor2DDiffFlowHistograms();
 virtual void GetPointersForOtherDiffCorrelators();
 virtual void GetPointersForNestedLoopsHistograms();
 virtual void GetPointersForMixedHarmonicsHistograms();
 virtual void GetPointersForControlHistograms();
 virtual void GetPointersForBootstrap();
 virtual void GetPointersForCRC();
 virtual void GetPointersForCRCVZ();
 virtual void GetPointersForCRCZDC();
 virtual void GetPointersForCRC2();
 virtual void GetPointersForCRCPt();
 virtual void GetPointersForQVec();
 virtual void GetPointersForCME();
 virtual void GetPointersForFlowQC();
 virtual void GetPointersForFlowSPZDC();
 virtual void GetPointersForFlowSPVZ();
 virtual void GetPointersForVarious();
 
 // 5.) other methods:
 TProfile* MakePtProjection(TProfile2D *profilePtEta) const;
 TProfile* MakeEtaProjection(TProfile2D *profilePtEta) const;
 virtual void WriteHistograms(TString outputFileName);
 virtual void WriteHistograms(TDirectoryFile *outputFileName);
 virtual Int_t GetCRCBin(Int_t c, Int_t y, Int_t c2, Int_t y2);
 virtual Int_t GetCRCVZBin(Int_t c, Int_t c2);
 virtual Int_t GetCRCQVecBin(Int_t c, Int_t y);
 virtual Int_t GetCRCRunBin(Int_t RunNum);
 virtual Int_t GetCRCCenBin(Double_t Centrality);
 virtual Int_t GetCRCPtBin(Double_t pt);
 
 virtual Double_t GetSPZDChar(Int_t har, Double_t QRe,Double_t QIm,Double_t ZARe,Double_t ZAIm,Double_t ZCRe,Double_t ZCIm);
 
 // **** SETTERS and GETTERS ****
 
 // 0.) base:
 void SetHistList(TList* const hlist) {this->fHistList = hlist;}
 TList* GetHistList() const {return this->fHistList;}
 
 // 1.) common:
 void SetBookOnlyBasicCCH(Bool_t const bobcch) {this->fBookOnlyBasicCCH = bobcch;};
 Bool_t GetBookOnlyBasicCCH() const {return this->fBookOnlyBasicCCH;};
 void SetCommonHists(AliFlowCommonHist* const ch) {this->fCommonHists = ch;};
 AliFlowCommonHist* GetCommonHists() const {return this->fCommonHists;};
 void SetCommonHists2nd(AliFlowCommonHist* const ch2nd) {this->fCommonHists2nd = ch2nd;};
 AliFlowCommonHist* GetCommonHists2nd() const {return this->fCommonHists2nd;};
 void SetCommonHists4th(AliFlowCommonHist* const ch4th) {this->fCommonHists4th = ch4th;};
 AliFlowCommonHist* GetCommonHists4th() const {return this->fCommonHists4th;};
 void SetCommonHists6th(AliFlowCommonHist* const ch6th) {this->fCommonHists6th = ch6th;};
 AliFlowCommonHist* GetCommonHists6th() const {return this->fCommonHists6th;};
 void SetCommonHists8th(AliFlowCommonHist* const ch8th) {this->fCommonHists8th = ch8th;};
 AliFlowCommonHist* GetCommonHists8th() const {return this->fCommonHists8th;};
 void SetCommonHistsResults2nd(AliFlowCommonHistResults* const chr2nd) {this->fCommonHistsResults2nd = chr2nd;};
 AliFlowCommonHistResults* GetCommonHistsResults2nd() const {return this->fCommonHistsResults2nd;};
 void SetCommonHistsResults4th(AliFlowCommonHistResults* const chr4th) {this->fCommonHistsResults4th = chr4th;};
 AliFlowCommonHistResults* GetCommonHistsResults4th() const {return this->fCommonHistsResults4th;};
 void SetCommonHistsResults6th(AliFlowCommonHistResults* const chr6th) {this->fCommonHistsResults6th = chr6th;};
 AliFlowCommonHistResults* GetCommonHistsResults6th() const {return this->fCommonHistsResults6th;};
 void SetCommonHistsResults8th(AliFlowCommonHistResults* const chr8th) {this->fCommonHistsResults8th = chr8th;};
 AliFlowCommonHistResults* GetCommonHistsResults8th() const {return this->fCommonHistsResults8th;};
 void SetCommonConstants(TProfile* const cc) {this->fCommonConstants = cc;};
 TProfile* GetCommonConstants() const {return this->fCommonConstants;};
 void SetFillMultipleControlHistograms(Bool_t const fmch) {this->fFillMultipleControlHistograms = fmch;};
 Bool_t GetFillMultipleControlHistograms() const {return this->fFillMultipleControlHistograms;};
 void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
 Int_t GetHarmonic() const {return this->fHarmonic;};
 void SetAnalysisLabel(const char *aLabel) {this->fAnalysisLabel->Append(*aLabel);}; // to be improved (Append(*aLabel) changed into Append(aLabel))
 TString *GetAnalysisLabel() const {return this->fAnalysisLabel;};
 void SetPrintFinalResults(Bool_t const printOrNot, Int_t const i) {this->fPrintFinalResults[i] = printOrNot;};
 Bool_t GetPrintFinalResults(Int_t i) const {return this->fPrintFinalResults[i];};
 
 // 2a.) particle weights:
 void SetWeightsList(TList* const wlist) {this->fWeightsList = wlist;}
 TList* GetWeightsList() const {return this->fWeightsList;}
 void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
 Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
 void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
 Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
 void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
 Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
 void SetUseTrackWeights(Bool_t const uTrackW) {this->fUseTrackWeights = uTrackW;};
 Bool_t GetUseTrackWeights() const {return this->fUseTrackWeights;};
 void SetUsePhiEtaWeights(Bool_t const uPhiEtaW) {this->fUsePhiEtaWeights = uPhiEtaW;};
 Bool_t GetUsePhiEtaWeights() const {return this->fUsePhiEtaWeights;};
 void SetUseParticleWeights(TProfile* const uPW) {this->fUseParticleWeights = uPW;};
 TProfile* GetUseParticleWeights() const {return this->fUseParticleWeights;};
 void SetPhiWeights(TH1F* const histPhiWeights) {this->fPhiWeightsRPs = histPhiWeights;};
 TH1F* GetPhiWeights() const {return this->fPhiWeightsRPs;};
 
 // 2b.) event weights:
 void SetMultiplicityWeight(const char *multiplicityWeight) {*this->fMultiplicityWeight = multiplicityWeight;};
 void SetMultiplicityIs(AliFlowCommonConstants::ERefMultSource mi) {this->fMultiplicityIs = mi;};
 
 // 3.) Reference flow:
 // Flags:
 void SetIntFlowFlags(TProfile* const intFlowFlags) {this->fIntFlowFlags = intFlowFlags;};
 TProfile* GetIntFlowFlags() const {return this->fIntFlowFlags;};
 void SetApplyCorrectionForNUA(Bool_t const applyCorrectionForNUA) {this->fApplyCorrectionForNUA = applyCorrectionForNUA;};
 Bool_t GetApplyCorrectionForNUA() const {return this->fApplyCorrectionForNUA;};
 void SetApplyCorrectionForNUAVsM(Bool_t const applyCorrectionForNUAVsM) {this->fApplyCorrectionForNUAVsM = applyCorrectionForNUAVsM;};
 Bool_t GetApplyCorrectionForNUAVsM() const {return this->fApplyCorrectionForNUAVsM;};
 void SetnBinsMult(Int_t const nbm) {this->fnBinsMult = nbm;};
 Int_t GetnBinsMult() const {return this->fnBinsMult;};
 void SetMinMult(Double_t const minm) {this->fMinMult = minm;};
 Double_t GetMinMult() const {return this->fMinMult;};
 void SetMaxMult(Double_t const maxm) {this->fMaxMult = maxm;};
 Double_t GetMaxMult() const {return this->fMaxMult;};
 void SetPropagateErrorAlsoFromNIT(Bool_t const peafNIT) {this->fPropagateErrorAlsoFromNIT = peafNIT;};
 Bool_t GetPropagateErrorAlsoFromNIT() const {return this->fPropagateErrorAlsoFromNIT;};
 void SetCalculateCumulantsVsM(Bool_t const ccvm) {this->fCalculateCumulantsVsM = ccvm;};
 Bool_t GetCalculateCumulantsVsM() const {return this->fCalculateCumulantsVsM;};
 void SetCalculateAllCorrelationsVsM(Bool_t const cacvm) {this->fCalculateAllCorrelationsVsM = cacvm;};
 Bool_t GetCalculateAllCorrelationsVsM() const {return this->fCalculateAllCorrelationsVsM;};
 void SetMinimumBiasReferenceFlow(Bool_t const mmrf) {this->fMinimumBiasReferenceFlow = mmrf;};
 Bool_t GetMinimumBiasReferenceFlow() const {return this->fMinimumBiasReferenceFlow;};
 void SetForgetAboutCovariances(Bool_t const fac) {this->fForgetAboutCovariances = fac;};
 Bool_t GetForgetAboutCovariances() const {return this->fForgetAboutCovariances;};
 void SetStoreVarious(Bool_t const spdfoe) {this->fStoreVarious = spdfoe;};
 Bool_t GetStoreVarious() const {return this->fStoreVarious;};
 void SetExactNoRPs(Int_t const enr) {this->fExactNoRPs = enr;};
 Int_t GetExactNoRPs() const {return this->fExactNoRPs;};
 void SetUse2DHistograms(Bool_t const u2dh){this->fUse2DHistograms = u2dh;if(u2dh){this->fStoreControlHistograms = kTRUE;}};
 Bool_t GetUse2DHistograms() const {return this->fUse2DHistograms;};
 void SetFillProfilesVsMUsingWeights(Bool_t const fpvmuw){this->fFillProfilesVsMUsingWeights = fpvmuw;};
 Bool_t GetFillProfilesVsMUsingWeights() const {return this->fFillProfilesVsMUsingWeights;};
 void SetUseQvectorTerms(Bool_t const uqvt){this->fUseQvectorTerms = uqvt;if(uqvt){this->fStoreControlHistograms = kTRUE;}};
 Bool_t GetUseQvectorTerms() const {return this->fUseQvectorTerms;};
 
 // Reference flow profiles:
 void SetAvMultiplicity(TProfile* const avMultiplicity) {this->fAvMultiplicity = avMultiplicity;};
 TProfile* GetAvMultiplicity() const {return this->fAvMultiplicity;};
 void SetIntFlowCorrelationsPro(TProfile* const intFlowCorrelationsPro) {this->fIntFlowCorrelationsPro = intFlowCorrelationsPro;};
 TProfile* GetIntFlowCorrelationsPro() const {return this->fIntFlowCorrelationsPro;};
 void SetIntFlowSquaredCorrelationsPro(TProfile* const ifscp) {this->fIntFlowSquaredCorrelationsPro = ifscp;};
 TProfile* GetIntFlowSquaredCorrelationsPro() const {return this->fIntFlowSquaredCorrelationsPro;};
 void SetIntFlowCorrelationsVsMPro(TProfile* const ifcvp, Int_t const ci) {this->fIntFlowCorrelationsVsMPro[ci] = ifcvp;};
 TProfile* GetIntFlowCorrelationsVsMPro(Int_t const ci) const {return this->fIntFlowCorrelationsVsMPro[ci];};
 void SetIntFlowSquaredCorrelationsVsMPro(TProfile* const ifscvp, Int_t const ci) {this->fIntFlowSquaredCorrelationsVsMPro[ci] = ifscvp;};
 TProfile* GetIntFlowSquaredCorrelationsVsMPro(Int_t const ci) const {return this->fIntFlowSquaredCorrelationsVsMPro[ci];};
 void SetIntFlowCorrelationsAllPro(TProfile* const intFlowCorrelationsAllPro) {this->fIntFlowCorrelationsAllPro = intFlowCorrelationsAllPro;};
 TProfile* GetIntFlowCorrelationsAllPro() const {return this->fIntFlowCorrelationsAllPro;};
 void SetIntFlowExtraCorrelationsPro(TProfile* const intFlowExtraCorrelationsPro) {this->fIntFlowExtraCorrelationsPro = intFlowExtraCorrelationsPro;};
 TProfile* GetIntFlowExtraCorrelationsPro() const {return this->fIntFlowExtraCorrelationsPro;};
 void SetIntFlowProductOfCorrelationsPro(TProfile* const intFlowProductOfCorrelationsPro) {this->fIntFlowProductOfCorrelationsPro = intFlowProductOfCorrelationsPro;};
 TProfile* GetIntFlowProductOfCorrelationsPro() const {return this->fIntFlowProductOfCorrelationsPro;};
 void SetIntFlowProductOfCorrelationsVsMPro(TProfile* const ifpocvm, Int_t const pi) {this->fIntFlowProductOfCorrelationsVsMPro[pi] = ifpocvm;};
 TProfile* GetIntFlowProductOfCorrelationsVsMPro(Int_t const pi) const {return this->fIntFlowProductOfCorrelationsVsMPro[pi];};
 void SetIntFlowProductOfCorrectionTermsForNUAPro(TProfile* const ifpoctfNUA) {this->fIntFlowProductOfCorrectionTermsForNUAPro = ifpoctfNUA;};
 TProfile* GetIntFlowProductOfCorrectionTermsForNUAPro() const {return this->fIntFlowProductOfCorrectionTermsForNUAPro;};
 void SetIntFlowCorrectionTermsForNUAPro(TProfile* const ifctfnp, Int_t const sc) {this->fIntFlowCorrectionTermsForNUAPro[sc] = ifctfnp;};
 TProfile* GetIntFlowCorrectionTermsForNUAPro(Int_t sc) const {return this->fIntFlowCorrectionTermsForNUAPro[sc];};
 void SetIntFlowCorrectionTermsForNUAVsMPro(TProfile* const ifctfnpvm, Int_t const sc, Int_t const ci) {this->fIntFlowCorrectionTermsForNUAVsMPro[sc][ci] = ifctfnpvm;};
 TProfile* GetIntFlowCorrectionTermsForNUAVsMPro(Int_t sc, Int_t ci) const {return this->fIntFlowCorrectionTermsForNUAVsMPro[sc][ci];};
 // integrated flow histograms holding all results:
 void SetIntFlowCorrelationsHist(TH1D* const intFlowCorrelationsHist) {this->fIntFlowCorrelationsHist = intFlowCorrelationsHist;};
 TH1D* GetIntFlowCorrelationsHist() const {return this->fIntFlowCorrelationsHist;};
 void SetIntFlowCorrelationsVsMHist(TH1D* const ifcvmh, Int_t const ci) {this->fIntFlowCorrelationsVsMHist[ci] = ifcvmh;};
 TH1D* GetIntFlowCorrelationsVsMHist(Int_t const ci) const {return this->fIntFlowCorrelationsVsMHist[ci];};
 void SetIntFlowCorrelationsAllHist(TH1D* const intFlowCorrelationsAllHist) {this->fIntFlowCorrelationsAllHist = intFlowCorrelationsAllHist;};
 TH1D* GetIntFlowCorrelationsAllHist() const {return this->fIntFlowCorrelationsAllHist;};
 void SetIntFlowCorrectionTermsForNUAHist(TH1D* const ifctfnh, Int_t const sc) {this->fIntFlowCorrectionTermsForNUAHist[sc] = ifctfnh;};
 TH1D* GetIntFlowCorrectionTermsForNUAHist(Int_t sc) const {return this->fIntFlowCorrectionTermsForNUAHist[sc];};
 void SetIntFlowCovariances(TH1D* const intFlowCovariances) {this->fIntFlowCovariances = intFlowCovariances;};
 TH1D* GetIntFlowCovariances() const {return this->fIntFlowCovariances;};
 void SetIntFlowSumOfEventWeights(TH1D* const intFlowSumOfEventWeights, Int_t const power) {this->fIntFlowSumOfEventWeights[power] = intFlowSumOfEventWeights;};
 TH1D* GetIntFlowSumOfEventWeights(Int_t power) const {return this->fIntFlowSumOfEventWeights[power];};
 void SetIntFlowSumOfProductOfEventWeights(TH1D* const intFlowSumOfProductOfEventWeights) {this->fIntFlowSumOfProductOfEventWeights = intFlowSumOfProductOfEventWeights;};
 TH1D* GetIntFlowSumOfProductOfEventWeights() const {return this->fIntFlowSumOfProductOfEventWeights;};
 void SetIntFlowCovariancesVsM(TH1D* const ifcvm, Int_t ci) {this->fIntFlowCovariancesVsM[ci] = ifcvm;};
 TH1D* GetIntFlowCovariancesVsM(Int_t ci) const {return this->fIntFlowCovariancesVsM[ci];};
 void SetIntFlowSumOfEventWeightsVsM(TH1D* const ifsoewvm, Int_t si, Int_t lc) {this->fIntFlowSumOfEventWeightsVsM[si][lc] = ifsoewvm;};
 TH1D* GetIntFlowSumOfEventWeightsVsM(Int_t si, Int_t lc) const {return this->fIntFlowSumOfEventWeightsVsM[si][lc];};
 void SetIntFlowSumOfProductOfEventWeightsVsM(TH1D* const ifsopoevm, Int_t si) {this->fIntFlowSumOfProductOfEventWeightsVsM[si] = ifsopoevm;};
 TH1D* GetIntFlowSumOfProductOfEventWeightsVsM(Int_t si) const {return this->fIntFlowSumOfProductOfEventWeightsVsM[si];};
 void SetIntFlowCovariancesNUA(TH1D* const intFlowCovariancesNUA) {this->fIntFlowCovariancesNUA = intFlowCovariancesNUA;};
 TH1D* GetIntFlowCovariancesNUA() const {return this->fIntFlowCovariancesNUA;};
 void SetIntFlowSumOfEventWeightsNUA(TH1D* const ifsoewNUA, Int_t const sc, Int_t const power) {this->fIntFlowSumOfEventWeightsNUA[sc][power] = ifsoewNUA;};
 TH1D* GetIntFlowSumOfEventWeightsNUA(Int_t sc, Int_t power) const {return this->fIntFlowSumOfEventWeightsNUA[sc][power];};
 void SetIntFlowSumOfProductOfEventWeightsNUA(TH1D* const ifsopoewNUA) {this->fIntFlowSumOfProductOfEventWeightsNUA = ifsopoewNUA;};
 TH1D* GetIntFlowSumOfProductOfEventWeightsNUA() const {return this->fIntFlowSumOfProductOfEventWeightsNUA;};
 void SetIntFlowQcumulants(TH1D* const intFlowQcumulants) {this->fIntFlowQcumulants = intFlowQcumulants;};
 TH1D* GetIntFlowQcumulants() const {return this->fIntFlowQcumulants;};
 void SetIntFlowQcumulantsVsM(TH1D* const intFlowQcumulantsVsM, Int_t co) {this->fIntFlowQcumulantsVsM[co] = intFlowQcumulantsVsM;};
 TH1D* GetIntFlowQcumulantsVsM(Int_t co) const {return this->fIntFlowQcumulantsVsM[co];};
 void SetIntFlowQcumulantsRebinnedInM(TH1D* const ifqcrim) {this->fIntFlowQcumulantsRebinnedInM = ifqcrim;};
 TH1D* GetIntFlowQcumulantsRebinnedInM() const {return this->fIntFlowQcumulantsRebinnedInM;};
 void SetIntFlowQcumulantsErrorSquaredRatio(TH1D* const ifqcesr) {this->fIntFlowQcumulantsErrorSquaredRatio = ifqcesr;};
 TH1D* GetIntFlowQcumulantsErrorSquaredRatio() const {return this->fIntFlowQcumulantsErrorSquaredRatio;};
 void SetIntFlow(TH1D* const intFlow) {this->fIntFlow = intFlow;};
 TH1D* GetIntFlow() const {return this->fIntFlow;};
 void SetIntFlowVsM(TH1D* const intFlowVsM, Int_t co) {this->fIntFlowVsM[co] = intFlowVsM;};
 TH1D* GetIntFlowVsM(Int_t co) const {return this->fIntFlowVsM[co];};
 void SetIntFlowRebinnedInM(TH1D* const ifrim) {this->fIntFlowRebinnedInM = ifrim;};
 TH1D* GetIntFlowRebinnedInM() const {return this->fIntFlowRebinnedInM;};
 void SetIntFlowDetectorBias(TH1D* const ifdb) {this->fIntFlowDetectorBias = ifdb;};
 TH1D* GetIntFlowDetectorBias() const {return this->fIntFlowDetectorBias;};
 void SetIntFlowDetectorBiasVsM(TH1D* const ifdbvm, Int_t ci) {this->fIntFlowDetectorBiasVsM[ci] = ifdbvm;};
 TH1D* GetIntFlowDetectorBiasVsM(Int_t ci) const {return this->fIntFlowDetectorBiasVsM[ci];};
 
 // 4.) Differential flow:
 //  Flags:
 void SetDiffFlowFlags(TProfile* const diffFlowFlags) {this->fDiffFlowFlags = diffFlowFlags;};
 TProfile* GetDiffFlowFlags() const {return this->fDiffFlowFlags;};
 void SetCalculateDiffFlow(Bool_t const cdf) {this->fCalculateDiffFlow = cdf;};
 Bool_t GetCalculateDiffFlow() const {return this->fCalculateDiffFlow;};
 void SetCalculate2DDiffFlow(Bool_t const c2ddf) {this->fCalculate2DDiffFlow = c2ddf;};
 Bool_t GetCalculate2DDiffFlow() const {return this->fCalculate2DDiffFlow;};
 void SetCalculateDiffFlowVsEta(Bool_t const cdfve) {this->fCalculateDiffFlowVsEta = cdfve;};
 Bool_t GetCalculateDiffFlowVsEta() const {return this->fCalculateDiffFlowVsEta;};
 //  Profiles:
 //   1D:
 void SetDiffFlowCorrelationsPro(TProfile* const diffFlowCorrelationsPro, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowCorrelationsPro[i][j][k] = diffFlowCorrelationsPro;};
 TProfile* GetDiffFlowCorrelationsPro(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowCorrelationsPro[i][j][k];};
 void SetDiffFlowSquaredCorrelationsPro(TProfile* const diffFlowSquaredCorrelationsPro, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowSquaredCorrelationsPro[i][j][k] = diffFlowSquaredCorrelationsPro;};
 TProfile* GetDiffFlowSquaredCorrelationsPro(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowSquaredCorrelationsPro[i][j][k];};
 void SetDiffFlowProductOfCorrelationsPro(TProfile* const dfpocp, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowProductOfCorrelationsPro[i][j][k][l] = dfpocp;};
 TProfile* GetDiffFlowProductOfCorrelationsPro(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowProductOfCorrelationsPro[i][j][k][l];};
 void SetDiffFlowCorrectionTermsForNUAPro(TProfile* const dfctfnp, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowCorrectionTermsForNUAPro[i][j][k][l] = dfctfnp;};
 TProfile* GetDiffFlowCorrectionTermsForNUAPro(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowCorrectionTermsForNUAPro[i][j][k][l];};
 //   2D:
 void Set2DDiffFlowCorrelationsPro(TProfile2D* const p2ddfcp, Int_t const i, Int_t const k) {this->f2DDiffFlowCorrelationsPro[i][k] = p2ddfcp;};
 TProfile2D* Get2DDiffFlowCorrelationsPro(Int_t i, Int_t k) const {return this->f2DDiffFlowCorrelationsPro[i][k];};
 //   Other differential correlators:
 void SetOtherDiffCorrelators(TProfile* const odc,Int_t const i,Int_t const j,Int_t const k,Int_t const l) {this->fOtherDiffCorrelators[i][j][k][l] = odc;};
 TProfile* GetOtherDiffCorrelators(Int_t i,Int_t j,Int_t k,Int_t l) const {return this->fOtherDiffCorrelators[i][j][k][l];};
 // histograms:
 void SetDiffFlowCorrelationsHist(TH1D* const diffFlowCorrelationsHist, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowCorrelationsHist[i][j][k] = diffFlowCorrelationsHist;};
 TH1D* GetDiffFlowCorrelationsHist(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowCorrelationsHist[i][j][k];};
 void SetDiffFlowCovariances(TH1D* const diffFlowCovariances, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowCovariances[i][j][k] = diffFlowCovariances;};
 TH1D* GetDiffFlowCovariances(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowCovariances[i][j][k];};
 void SetDiffFlowCumulants(TH1D* const diffFlowCumulants, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowCumulants[i][j][k] = diffFlowCumulants;};
 TH1D* GetDiffFlowCumulants(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowCumulants[i][j][k];};
 void SetDiffFlowDetectorBias(TH1D* const dfdb, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowDetectorBias[i][j][k] = dfdb;};
 TH1D* GetDiffFlowDetectorBias(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowDetectorBias[i][j][k];};
 void SetDiffFlow(TH1D* const diffFlow, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlow[i][j][k] = diffFlow;};
 TH1D* GetDiffFlow(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlow[i][j][k];};
 void SetDiffFlowSumOfEventWeights(TH1D* const dfsoew, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowSumOfEventWeights[i][j][k][l] = dfsoew;};
 TH1D* GetDiffFlowSumOfEventWeights(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowSumOfEventWeights[i][j][k][l];};
 void SetDiffFlowSumOfProductOfEventWeights(TH1D* const dfsopoew, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowSumOfProductOfEventWeights[i][j][k][l] = dfsopoew;};
 TH1D* GetDiffFlowSumOfProductOfEventWeights(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowSumOfProductOfEventWeights[i][j][k][l];};
 void SetDiffFlowCorrectionTermsForNUAHist(TH1D* const dfctfnh, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowCorrectionTermsForNUAHist[i][j][k][l] = dfctfnh;};
 TH1D* GetDiffFlowCorrectionTermsForNUAHist(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowCorrectionTermsForNUAHist[i][j][k][l];};
 //  2D:
 void Set2DDiffFlowCumulants(TH2D* const h2ddfc, Int_t const i, Int_t const j) {this->f2DDiffFlowCumulants[i][j] = h2ddfc;};
 TH2D* Get2DDiffFlowCumulants(Int_t i, Int_t j) const {return this->f2DDiffFlowCumulants[i][j];};
 void Set2DDiffFlow(TH2D* const h2ddf, Int_t const i, Int_t const j) {this->f2DDiffFlow[i][j] = h2ddf;};
 TH2D* Get2DDiffFlow(Int_t i, Int_t j) const {return this->f2DDiffFlow[i][j];};
 
 // 5.) distributions of correlations:
 // profile:
 void SetDistributionsFlags(TProfile* const distributionsFlags) {this->fDistributionsFlags = distributionsFlags;};
 TProfile* GetDistributionsFlags() const {return this->fDistributionsFlags;};
 // flags:
 void SetStoreDistributions(Bool_t const storeDistributions) {this->fStoreDistributions = storeDistributions;};
 Bool_t GetStoreDistributions() const {return this->fStoreDistributions;};
 // # of bins for correlation axis in fDistributions[4], fCorrelation2468VsMult[4] and fCorrelationProduct2468VsMult[1]:
 void SetnBinsForCorrelations(Int_t const nb) {this->fnBinsForCorrelations = nb;};
 Int_t GetnBinsForCorrelations() const {return this->fnBinsForCorrelations;};
 // histograms:
 void SetDistributions(TH1D* const distributions, Int_t const i) {this->fDistributions[i] = distributions;};
 TH1D* GetDistributions(Int_t i) const {return this->fDistributions[i];};
 // min and max values of correlations (ci is correlations index [0=<2>,1=<4>,2=<6>,3=<8>]):
 void SetMinValueOfCorrelation(Int_t const ci, Double_t const minValue) {this->fMinValueOfCorrelation[ci] = minValue;};
 Double_t GetMinValueOfCorrelation(Int_t ci) const {return this->fMinValueOfCorrelation[ci];};
 void SetMaxValueOfCorrelation(Int_t const ci, Double_t const maxValue) {this->fMaxValueOfCorrelation[ci] = maxValue;};
 Double_t GetMaxValueOfCorrelation(Int_t ci) const {return this->fMaxValueOfCorrelation[ci];};
 // min and max values of correlation products:
 void SetMinValueOfCorrelationProduct(Int_t const cpi, Double_t const minValue) {this->fMinValueOfCorrelationProduct[cpi] = minValue;};
 Double_t GetMinValueOfCorrelationProduct(Int_t cpi) const {return this->fMinValueOfCorrelationProduct[cpi];};
 void SetMaxValueOfCorrelationProduct(Int_t const cpi, Double_t const maxValue) {this->fMaxValueOfCorrelationProduct[cpi] = maxValue;};
 Double_t GetMaxValueOfCorrelationProduct(Int_t cpi) const {return this->fMaxValueOfCorrelationProduct[cpi];};
 // min and max values of QvectorTerms:
 void SetMinValueOfQvectorTerms(Int_t const qvti, Double_t const minValue) {this->fMinValueOfQvectorTerms[qvti] = minValue;};
 Double_t GetMinValueOfQvectorTerms(Int_t qvti) const {return this->fMinValueOfQvectorTerms[qvti];};
 void SetMaxValueOfQvectorTerms(Int_t const qvti, Double_t const maxValue) {this->fMaxValueOfQvectorTerms[qvti] = maxValue;};
 Double_t GetMaxValueOfQvectorTerms(Int_t qvti) const {return this->fMaxValueOfQvectorTerms[qvti];};
 
 // x.) debugging and cross-checking:
 void SetNestedLoopsList(TList* const nllist) {this->fNestedLoopsList = nllist;};
 TList* GetNestedLoopsList() const {return this->fNestedLoopsList;};
 void SetEvaluateIntFlowNestedLoops(Bool_t const eifnl) {this->fEvaluateIntFlowNestedLoops = eifnl;};
 Bool_t GetEvaluateIntFlowNestedLoops() const {return this->fEvaluateIntFlowNestedLoops;};
 void SetEvaluateDiffFlowNestedLoops(Bool_t const edfnl) {this->fEvaluateDiffFlowNestedLoops = edfnl;};
 Bool_t GetEvaluateDiffFlowNestedLoops() const {return this->fEvaluateDiffFlowNestedLoops;};
 void SetMaxAllowedMultiplicity(Int_t const maxAllowedMultiplicity) {this->fMaxAllowedMultiplicity = maxAllowedMultiplicity;};
 Int_t GetMaxAllowedMultiplicity() const {return this->fMaxAllowedMultiplicity;};
 void SetEvaluateNestedLoops(TProfile* const enl) {this->fEvaluateNestedLoops = enl;};
 TProfile* GetEvaluateNestedLoops() const {return this->fEvaluateNestedLoops;};
 void SetIntFlowDirectCorrelations(TProfile* const ifdc) {this->fIntFlowDirectCorrelations = ifdc;};
 TProfile* GetIntFlowDirectCorrelations() const {return this->fIntFlowDirectCorrelations;};
 void SetIntFlowExtraDirectCorrelations(TProfile* const ifedc) {this->fIntFlowExtraDirectCorrelations = ifedc;};
 TProfile* GetIntFlowExtraDirectCorrelations() const {return this->fIntFlowExtraDirectCorrelations;};
 void SetIntFlowDirectCorrectionTermsForNUA(TProfile* const ifdctfn, Int_t const sc) {this->fIntFlowDirectCorrectionTermsForNUA[sc] = ifdctfn;};
 TProfile* GetIntFlowDirectCorrectionTermsForNUA(Int_t sc) const {return this->fIntFlowDirectCorrectionTermsForNUA[sc];};
 void SetCrossCheckInPtBinNo(Int_t const crossCheckInPtBinNo) {this->fCrossCheckInPtBinNo = crossCheckInPtBinNo;};
 Int_t GetCrossCheckInPtBinNo() const {return this->fCrossCheckInPtBinNo;};
 void SetCrossCheckInEtaBinNo(Int_t const crossCheckInEtaBinNo) {this->fCrossCheckInEtaBinNo = crossCheckInEtaBinNo;};
 Int_t GetCrossCheckInEtaBinNo() const {return this->fCrossCheckInEtaBinNo;};
 void SetNoOfParticlesInBin(TH1D* const noOfParticlesInBin) {this->fNoOfParticlesInBin = noOfParticlesInBin;};
 TH1D* GetNoOfParticlesInBin() const {return this->fNoOfParticlesInBin;};
 void SetDiffFlowDirectCorrelations(TProfile* const diffFlowDirectCorrelations,Int_t const i,Int_t const j,Int_t const k){this->fDiffFlowDirectCorrelations[i][j][k]=diffFlowDirectCorrelations;};
 TProfile* GetDiffFlowDirectCorrelations(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowDirectCorrelations[i][j][k];};
 void SetDiffFlowDirectCorrectionTermsForNUA(TProfile* const dfdctfn, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowDirectCorrectionTermsForNUA[i][j][k][l] = dfdctfn;};
 TProfile* GetDiffFlowDirectCorrectionTermsForNUA(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowDirectCorrectionTermsForNUA[i][j][k][l];};
 void SetOtherDirectDiffCorrelators(TProfile* const oddc, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fOtherDirectDiffCorrelators[i][j][k][l] = oddc;};
 TProfile* GetOtherDirectDiffCorrelators(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fOtherDirectDiffCorrelators[i][j][k][l];};
 void SetMixedHarmonicsNestedLoops(TProfile* const mhnl) {this->fMixedHarmonicsNestedLoops = mhnl;};
 TProfile* GetMixedHarmonicsNestedLoops() const {return this->fMixedHarmonicsNestedLoops;};
 
 // 9.) Mixed harmonics:
 void SetMixedHarmonicsList(TList* const mhlist) {this->fMixedHarmonicsList = mhlist;};
 void SetMixedHarmonicsFlags(TProfile* const mhFlags) {this->fMixedHarmonicsFlags = mhFlags;};
 TProfile* GetMixedHarmonicsFlags() const {return this->fMixedHarmonicsFlags;};
 void SetCalculateMixedHarmonics(Bool_t const cmh) {this->fCalculateMixedHarmonics = cmh;};
 Bool_t GetCalculateMixedHarmonics() const {return this->fCalculateMixedHarmonics;};
 void SetCalculateMixedHarmonicsVsM(Bool_t const cmhvm) {this->fCalculateMixedHarmonicsVsM = cmhvm;};
 Bool_t GetCalculateMixedHarmonicsVsM() const {return this->fCalculateMixedHarmonicsVsM;};
 void Set2pCorrelations(TProfile* const p2pCorr) {this->f2pCorrelations = p2pCorr;};
 TProfile* Get2pCorrelations() const {return this->f2pCorrelations;};
 void Set3pCorrelations(TProfile* const p3pCorr) {this->f3pCorrelations = p3pCorr;};
 TProfile* Get3pCorrelations() const {return this->f3pCorrelations;};
 void Set4pCorrelations(TProfile* const p4pCorr) {this->f4pCorrelations = p4pCorr;};
 TProfile* Get4pCorrelations() const {return this->f4pCorrelations;};
 void Set5pCorrelations(TProfile* const p5pCorr) {this->f5pCorrelations = p5pCorr;};
 TProfile* Get5pCorrelations() const {return this->f5pCorrelations;};
 void Set6pCorrelations(TProfile* const p6pCorr) {this->f6pCorrelations = p6pCorr;};
 TProfile* Get6pCorrelations() const {return this->f6pCorrelations;};
 void Set7pCorrelations(TProfile* const p7pCorr) {this->f7pCorrelations = p7pCorr;};
 TProfile* Get7pCorrelations() const {return this->f7pCorrelations;};
 void Set8pCorrelations(TProfile* const p8pCorr) {this->f8pCorrelations = p8pCorr;};
 TProfile* Get8pCorrelations() const {return this->f8pCorrelations;};
 void Set2pCumulants(TH1D* const p2pC) {this->f2pCumulants = p2pC;};
 TH1D* Get2pCumulants() const {return this->f2pCumulants;};
 void Set3pCumulants(TH1D* const p3pC) {this->f3pCumulants = p3pC;};
 TH1D* Get3pCumulants() const {return this->f3pCumulants;};
 void Set4pCumulants(TH1D* const p4pC) {this->f4pCumulants = p4pC;};
 TH1D* Get4pCumulants() const {return this->f4pCumulants;};
 void Set5pCumulants(TH1D* const p5pC) {this->f5pCumulants = p5pC;};
 TH1D* Get5pCumulants() const {return this->f5pCumulants;};
 void Set6pCumulants(TH1D* const p6pC) {this->f6pCumulants = p6pC;};
 TH1D* Get6pCumulants() const {return this->f6pCumulants;};
 void Set7pCumulants(TH1D* const p7pC) {this->f7pCumulants = p7pC;};
 TH1D* Get7pCumulants() const {return this->f7pCumulants;};
 void Set8pCumulants(TH1D* const p8pC) {this->f8pCumulants = p8pC;};
 TH1D* Get8pCumulants() const {return this->f8pCumulants;};
 void SetMixedHarmonicEventWeights(TH1D* const mhew, Int_t const power) {this->fMixedHarmonicEventWeights[power] = mhew;};
 TH1D* GetMixedHarmonicEventWeights(Int_t power) const {return this->fMixedHarmonicEventWeights[power];};
 void SetMixedHarmonicProductOfEventWeights(TH2D* const mhpoew) {this->fMixedHarmonicProductOfEventWeights = mhpoew;};
 TH2D* GetMixedHarmonicProductOfEventWeights() const {return this->fMixedHarmonicProductOfEventWeights;};
 void SetMixedHarmonicProductOfCorrelations(TProfile2D* const mhpoc) {this->fMixedHarmonicProductOfCorrelations = mhpoc;};
 TProfile2D* GetMixedHarmonicProductOfCorrelations() const {return this->fMixedHarmonicProductOfCorrelations;};
 
 // 10.) Control histograms:
 void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
 void SetControlHistogramsFlags(TProfile* const chf) {this->fControlHistogramsFlags = chf;};
 TProfile* GetControlHistogramsFlags() const {return this->fControlHistogramsFlags;};
 void SetStoreControlHistograms(Bool_t const sch) {this->fStoreControlHistograms = sch;};
 Bool_t GetStoreControlHistograms() const {return this->fStoreControlHistograms;};
 void SetCorrelationNoRPsVsRefMult(TH2D* const cnrvrm) {this->fCorrelationNoRPsVsRefMult = cnrvrm;};
 TH2D* GetCorrelationNoRPsVsRefMult() const {return this->fCorrelationNoRPsVsRefMult;};
 void SetCorrelationNoPOIsVsRefMult(TH2D* const cnpvrm) {this->fCorrelationNoPOIsVsRefMult = cnpvrm;};
 TH2D* GetCorrelationNoPOIsVsRefMult() const {return this->fCorrelationNoPOIsVsRefMult;};
 void SetCorrelationNoRPsVsNoPOIs(TH2D* const cnrvnp) {this->fCorrelationNoRPsVsNoPOIs = cnrvnp;};
 TH2D* GetCorrelationNoRPsVsNoPOIs() const {return this->fCorrelationNoRPsVsNoPOIs;};
 void SetCorrelation2468VsMult(TH2D* const c2468vm, Int_t const ci) {this->fCorrelation2468VsMult[ci] = c2468vm;};
 TH2D* GetCorrelation2468VsMult(Int_t ci) const {return this->fCorrelation2468VsMult[ci];};
 void SetCorrelationProduct2468VsMult(TH2D* const cp2468vm, Int_t const ci) {this->fCorrelationProduct2468VsMult[ci] = cp2468vm;};
 TH2D* GetCorrelationProduct2468VsMult(Int_t ci) const {return this->fCorrelationProduct2468VsMult[ci];};
 void SetQvectorTermsVsMult(TH2D* const qvtvm, Int_t const qvti) {this->fQvectorTermsVsMult[qvti] = qvtvm;};
 TH2D* GetQvectorTermsVsMult(Int_t qvti) const {return this->fQvectorTermsVsMult[qvti];};
 
 // 11.) Bootstrap:
 void SetBootstrapList(TList* const bl) {this->fBootstrapList = bl;};
 void SetBootstrapProfilesList(TList* const bpl) {this->fBootstrapProfilesList = bpl;};
 void SetBootstrapResultsList(TList* const brl) {this->fBootstrapResultsList = brl;};
 void SetBootstrapFlags(TProfile* const bf) {this->fBootstrapFlags = bf;};
 TProfile* GetBootstrapFlags() const {return this->fBootstrapFlags;};
 void SetUseBootstrap(Bool_t const ub) {this->fUseBootstrap = ub;};
 Bool_t GetUseBootstrap() const {return this->fUseBootstrap;};
 void SetUseBootstrapVsM(Bool_t const ubVsM) {this->fUseBootstrapVsM = ubVsM;};
 Bool_t GetUseBootstrapVsM() const {return this->fUseBootstrapVsM;};
 void SetnSubsamples(Int_t const ns) {this->fnSubsamples = ns;};
 Int_t GetnSubsamples() const {return this->fnSubsamples;};
 void SetBootstrapCorrelations(TProfile2D* const bcp) {this->fBootstrapCorrelations = bcp;};
 TProfile2D* GetBootstrapCorrelations() const {return this->fBootstrapCorrelations;};
 void SetBootstrapCorrelationsVsM(TProfile2D* const bcpVsM, Int_t const qvti) {this->fBootstrapCorrelationsVsM[qvti] = bcpVsM;};
 TProfile2D* GetBootstrapCorrelationsVsM(Int_t qvti) const {return this->fBootstrapCorrelationsVsM[qvti];};
 void SetBootstrapCumulants(TH2D* const bc) {this->fBootstrapCumulants = bc;};
 TH2D* GetBootstrapCumulants() const {return this->fBootstrapCumulants;};
 void SetBootstrapCumulantsVsM(TH2D* const bcpVsM, Int_t const qvti) {this->fBootstrapCumulantsVsM[qvti] = bcpVsM;};
 TH2D* GetBootstrapCumulantsVsM(Int_t qvti) const {return this->fBootstrapCumulantsVsM[qvti];};
 
 // 12.) CRC
 void SetCRCList(TList* const CRCL) {this->fCRCList = CRCL;};
 void SetCRCIntList(TList* const CRCL) {this->fCRCIntList = CRCL;};
 void SetCRCIntRbRList(TList* const CRCL) {this->fCRCIntRbRList = CRCL;};
 void SetCRCIntRunsList(TList* const CRCL, Int_t r) {this->fCRCIntRunsList[r] = CRCL;};
 void SetCRCVZList(TList* const CRCL) {this->fCRCVZList = CRCL;};
 void SetCRCVZRbRList(TList* const CRCL) {this->fCRCVZRbRList = CRCL;};
 void SetCRCVZRunsList(TList* const CRCL, Int_t r) {this->fCRCVZRunsList[r] = CRCL;};
 void SetCRCZDCList(TList* const CRCL) {this->fCRCZDCList = CRCL;};
 void SetCRCZDCRbRList(TList* const CRCL) {this->fCRCZDCRbRList = CRCL;};
 void SetCRCZDCRunsList(TList* const CRCL, Int_t r) {this->fCRCZDCRunsList[r] = CRCL;};
 void SetCRCPtList(TList* const CRCL) {this->fCRCPtList = CRCL;};
 void SetCRC2List(TList* const CRCL) {this->fCRC2List = CRCL;};
 void SetCRC2RbRList(TList* const CRCL) {this->fCRC2RbRList = CRCL;};
 void SetCRC2RunsList(TList* const CRCL, Int_t r) {this->fCRC2RunsList[r] = CRCL;};
 void SetCRCQVecList(TList* const CRCL) {this->fCRCQVecList = CRCL;};
 void SetCRCQVecListRun(TList* const CRCL, Int_t r) {this->fCRCQVecListRun[r] = CRCL;};
 void SetCRCFlags(TProfile* const CRCF) {this->fCRCFlags = CRCF;};
 TProfile* GetCRCFlags() const {return this->fCRCFlags;};
 void SetCalculateCRC(Bool_t const cCRC) {this->fCalculateCRC = cCRC;};
 Bool_t GetCalculateCRC() const {return this->fCalculateCRC;};
 void SetCalculateCRCPt(Bool_t const cCRC) {this->fCalculateCRCPt = cCRC;};
 Bool_t GetCalculateCRCPt() const {return this->fCalculateCRCPt;};
 void SetCalculateCME(Bool_t const cCRC) {this->fCalculateCME = cCRC;};
 Bool_t GetCalculateCME() const {return this->fCalculateCME;};
 void SetCalculateFlowQC(Bool_t const cCRC) {this->fCalculateFlowQC = cCRC;};
 Bool_t GetCalculateFlowQC() const {return this->fCalculateFlowQC;};
 void SetCalculateFlowZDC(Bool_t const cCRC) {this->fCalculateFlowZDC = cCRC;};
 Bool_t GetCalculateFlowZDC() const {return this->fCalculateFlowZDC;};
 void SetCalculateFlowVZ(Bool_t const cCRC) {this->fCalculateFlowVZ = cCRC;};
 Bool_t GetCalculateFlowVZ() const {return this->fCalculateFlowVZ;};
 void SetCalculateCRC2(Bool_t const cCRC) {this->fCalculateCRC2 = cCRC;};
 Bool_t GetCalculateCRC2() const {return this->fCalculateCRC2;};
  void SetCalculateCRCVZ(Bool_t const cCRC) {this->fCalculateCRCVZ = cCRC;};
  Bool_t GetCalculateCRCVZ() const {return this->fCalculateCRCVZ;};
  void SetCalculateCRCZDC(Bool_t const cCRC) {this->fCalculateCRCZDC = cCRC;};
  Bool_t GetCalculateCRCZDC() const {return this->fCalculateCRCZDC;};
 void SetUseVZERO(Bool_t const cCRC) {this->fUseVZERO = cCRC;};
 Bool_t GetUseVZERO() const {return this->fUseVZERO;};
 void SetUseZDC(Bool_t const cCRC) {this->fUseZDC = cCRC;};
 Bool_t GetUseZDC() const {return this->fUseZDC;};
 void SetRecenterZDC(Bool_t const cCRC) {this->fRecenterZDC = cCRC;};
 Bool_t GetRecenterZDC() const {return this->fRecenterZDC;};
 void SetDivSigma(Bool_t const cCRC) {this->fDivSigma = cCRC;};
 Bool_t GetDivSigma() const {return this->fDivSigma;};
 void SetInvertZDC(Bool_t const cCRC) {this->fInvertZDC = cCRC;};
 Bool_t GetInvertZDC() const {return this->fInvertZDC;};
 void SetQAZDCCuts(Bool_t const cCRC) {this->fQAZDCCuts = cCRC;};
 Bool_t GetQAZDCCuts() const {return this->fQAZDCCuts;};
 void SetTestSin(Bool_t const cCRC) {this->fCRCTestSin = cCRC;};
 Bool_t GetTestSin() const {return this->fCRCTestSin;};
 void SetNUAforCRC(Bool_t const cCRC) {this->fNUAforCRC = cCRC;};
 Bool_t GetNUAforCRC() const {return this->fNUAforCRC;};
 void SetUseCRCRecenter(Bool_t const cCRC) {this->fUseCRCRecenter = cCRC;};
 Bool_t GetUseCRCRecenter() const {return this->fUseCRCRecenter;};
 void SetCRCEtaRange(Double_t const etamin, Double_t const etamax) {this->fCRCEtaMin = etamin; this->fCRCEtaMax = etamax;};
 void SetCRCQVecWeightsList(TList* const wlist) {this->fCRCQVecWeightsList = wlist;}
 TList* GetCRCQVecWeightsList() const {return this->fCRCQVecWeightsList;}
  void SetCRCZDCCalibList(TList* const wlist) {this->fCRCZDCCalibList = wlist;}
  TList* GetCRCZDCCalibList() const {return this->fCRCZDCCalibList;}
  void SetZDCESEList(TList* const kList) {this->fZDCESEList = kList;};
  TList* GetZDCESEList() const {return this->fZDCESEList;};
 // 12.a) EbE Corr:
 void SetCRCCorrProdTempHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCCorrProdTempHist[c][eg][h] = TH;};
 TH1D* GetCRCCorrProdTempHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCCorrProdTempHist[c][eg][h];};
 // 12.b) Final histo:
 void SetCRCCorrHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCCorrHist[c][eg][h] = TH;};
 TH1D* GetCRCCorrHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCCorrHist[c][eg][h];};
 void SetCRCCumHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCCumHist[c][eg][h] = TH;};
 TH1D* GetCRCCumHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCCumHist[c][eg][h];};
 void SetCRCCFunHist(TH1D* const TH, Int_t const eg, Int_t const h) {this->fCRCCFunHist[eg][h] = TH;};
 TH1D* GetCRCCFunHist(Int_t const eg, Int_t const h) const {return this->fCRCCFunHist[eg][h];};
 // 12.c) Covariances:
 void SetCRCCorrProd2p2pPro(TProfile* const TP, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCCorrProd2p2pPro[c][eg][h] = TP;};
 TProfile* GetCRCCorrProd2p2pPro(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCCorrProd2p2pPro[c][eg][h];};
 void SetCRCCovHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCCovHist[c][eg][h] = TH;};
 TH1D* GetCRCCovHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCCovHist[c][eg][h];};
 // 12.d) NUA corrections:
 void SetCRCNUATermsHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCNUATermsHist[c][eg][h] = TH;};
 TH1D* GetCRCNUATermsHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCNUATermsHist[c][eg][h];};
 
 void SetCRCCorrPro(TProfile* const TP, Int_t const r, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCCorrPro[r][c][eg][h] = TP;};
 TProfile* GetCRCCorrPro(Int_t const r, Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCCorrPro[r][c][eg][h];};
 void SetCRCSumWeigHist(TH1D* const TH, Int_t const r, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCSumWeigHist[r][c][eg][h] = TH;};
 TH1D* GetCRCSumWeigHist(Int_t const r, Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCSumWeigHist[r][c][eg][h];};
 void SetCRCNUATermsPro(TProfile* const TP, Int_t const r, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCNUATermsPro[r][c][eg][h] = TP;};
 TProfile* GetCRCNUATermsPro(Int_t const r, Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCNUATermsPro[r][c][eg][h];};
 
 // 12.e) Q Vectors:
 void SetCRCQ2ReHist(TProfile* const TH, Int_t const r) {this->fCRCQ2Re[r] = TH;};
 TProfile* GetCRCQ2ReHist(Int_t const r) const {return this->fCRCQ2Re[r];};
 void SetCRCQ2ImHist(TProfile* const TH, Int_t const r) {this->fCRCQ2Im[r] = TH;};
 TProfile* GetCRCQ2ImHist(Int_t const r) const {return this->fCRCQ2Im[r];};
 void SetCRCQ2ReCorrHist(TProfile* const TH, Int_t const r) {this->fCRCQ2ReCorr[r] = TH;};
 TProfile* GetCRCQ2ReCorrHist(Int_t const r) const {return this->fCRCQ2ReCorr[r];};
 void SetCRCQ2ImCorrHist(TProfile* const TH, Int_t const r) {this->fCRCQ2ImCorr[r] = TH;};
 TProfile* GetCRCQ2ImCorrHist(Int_t const r) const {return this->fCRCQ2ImCorr[r];};
 void SetCRCPhiHist(TH2D* const TH, Int_t const r, Int_t const c, Int_t const i) {this->fCRCPhiHist[r][c][i] = TH;};
 TH2D* GetCRCPhiHist(Int_t const r, Int_t const c, Int_t const i) const {return this->fCRCPhiHist[r][c][i];};
 
 void SetCRCVZEvPlA(TH1D* const TH, Int_t const r, Int_t const c, Int_t const h) {this->fCRCVZEvPlA[r][c][h] = TH;};
 TH1D* GetCRCVZEvPlA(Int_t const r, Int_t const c, Int_t const h) const {return this->fCRCVZEvPlA[r][c][h];};
 void SetCRCVZEvPlC(TH1D* const TH, Int_t const r, Int_t const c, Int_t const h) {this->fCRCVZEvPlC[r][c][h] = TH;};
 TH1D* GetCRCVZEvPlC(Int_t const r, Int_t const c, Int_t const h) const {return this->fCRCVZEvPlC[r][c][h];};
 void SetCRCVZQVecAHist(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCVZQVecA[r][c] = TH;};
 TProfile* GetCRCVZQVecAHist(Int_t const r, Int_t const c) const {return this->fCRCVZQVecA[r][c];};
 void SetCRCVZQVecCHist(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCVZQVecC[r][c] = TH;};
 TProfile* GetCRCVZQVecCHist(Int_t const r, Int_t const c) const {return this->fCRCVZQVecC[r][c];};
 void SetCRCVZQVecCov(TProfile* const TH, Int_t const r, Int_t const i) {this->fCRCVZQVecCov[r][i] = TH;};
 TProfile* GetCRCVZQVecCov(Int_t const r, Int_t const i) const {return this->fCRCVZQVecCov[r][i];};
 
 void SetCRCZDCEvPlA(TH1D* const TH, Int_t const r, Int_t const c) {this->fCRCZDCEvPlA[r][c] = TH;};
 TH1D* GetCRCZDCEvPlA(Int_t const r, Int_t const c) const {return this->fCRCZDCEvPlA[r][c];};
 void SetCRCZDCEvPlC(TH1D* const TH, Int_t const r, Int_t const c) {this->fCRCZDCEvPlC[r][c] = TH;};
 TH1D* GetCRCZDCEvPlC(Int_t const r, Int_t const c) const {return this->fCRCZDCEvPlC[r][c];};
 void SetCRCZDCQVecAHist(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCZDCQVecA[r][c] = TH;};
 TProfile* GetCRCZDCQVecAHist(Int_t const r, Int_t const c) const {return this->fCRCZDCQVecA[r][c];};
 void SetCRCZDCQVecCHist(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCZDCQVecC[r][c] = TH;};
 TProfile* GetCRCZDCQVecCHist(Int_t const r, Int_t const c) const {return this->fCRCZDCQVecC[r][c];};
 void SetCRCZDCQVecACorrHist(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCZDCQVecACorr[r][c] = TH;};
 TProfile* GetCRCZDCQVecACorrHist(Int_t const r, Int_t const c) const {return this->fCRCZDCQVecACorr[r][c];};
 void SetCRCZDCQVecCCorrHist(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCZDCQVecCCorr[r][c] = TH;};
 TProfile* GetCRCZDCQVecCCorrHist(Int_t const r, Int_t const c) const {return this->fCRCZDCQVecCCorr[r][c];};
 void SetCRCZDCQVecCov(TProfile* const TH, Int_t const r, Int_t const i) {this->fCRCZDCQVecCov[r][i] = TH;};
 TProfile* GetCRCZDCQVecCov(Int_t const r, Int_t const i) const {return this->fCRCZDCQVecCov[r][i];};
 
 void SetCRCVZvsZDCCov(TProfile* const TH, Int_t const r, Int_t const i) {this->fCRCVZvsZDCCov[r][i] = TH;};
 TProfile* GetCRCVZvsZDCCov(Int_t const r, Int_t const i) const {return this->fCRCVZvsZDCCov[r][i];};
 
 // CRC VZERO:
 // 12.a) EbE Corr:
 void SetCRCVZCorrPro(TProfile* const TP, Int_t const r, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCVZCorrPro[r][c][eg][h] = TP;};
 TProfile* GetCRCVZCorrPro(Int_t const r, Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCVZCorrPro[r][c][eg][h];};
 void SetCRCVZCorrProd2p2pHist(TProfile* const TH, Int_t const r, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCVZCorrProd2p2pHist[r][c][eg][h] = TH;};
 TProfile* GetCRCVZCorrProd2p2pHist(Int_t const r, Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCVZCorrProd2p2pHist[r][c][eg][h];};
 void SetCRCVZNUAPro(TProfile* const TP, Int_t const r, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCVZNUAPro[r][c][eg][h] = TP;};
 TProfile* GetCRCVZNUAPro(Int_t const r, Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCVZNUAPro[r][c][eg][h];};
 // 12.b) Final histo:
 void SetCRCVZCorrHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCVZCorrHist[c][eg][h] = TH;};
 TH1D* GetCRCVZCorrHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCVZCorrHist[c][eg][h];};
 void SetCRCVZCFunHist(TH1D* const TH, Int_t const eg, Int_t const h) {this->fCRCVZCFunHist[eg][h] = TH;};
 TH1D* GetCRCVZCFunHist(Int_t const eg, Int_t const h) const {return this->fCRCVZCFunHist[eg][h];};
 // 12.c) Covariances:
 void SetCRCVZCovHist(TH2D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCVZCovHist[c][eg][h] = TH;};
 TH2D* GetCRCVZCovHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCVZCovHist[c][eg][h];};
 void SetCRCVZCorrProdTempHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCVZCorrProdTempHist[c][eg][h] = TH;};
 TH1D* GetCRCVZCorrProdTempHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCVZCorrProdTempHist[c][eg][h];};
 
 // CRC ZDC:
 // 12.a) EbE Corr:
 void SetCRCZDCCorrPro(TProfile* const TP, Int_t const r, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCZDCCorrPro[r][c][eg][h] = TP;};
 TProfile* GetCRCZDCCorrPro(Int_t const r, Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCZDCCorrPro[r][c][eg][h];};
 void SetCRCZDCCorrProd2p2pHist(TProfile* const TH, Int_t const r, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCZDCCorrProd2p2pHist[r][c][eg][h] = TH;};
 TProfile* GetCRCZDCCorrProd2p2pHist(Int_t const r, Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCZDCCorrProd2p2pHist[r][c][eg][h];};
 void SetCRCZDCNUAPro(TProfile* const TP, Int_t const r, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCZDCNUAPro[r][c][eg][h] = TP;};
 TProfile* GetCRCZDCNUAPro(Int_t const r, Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCZDCNUAPro[r][c][eg][h];};
 // 12.b) Final histo:
 void SetCRCZDCCorrHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCZDCCorrHist[c][eg][h] = TH;};
 TH1D* GetCRCZDCCorrHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCZDCCorrHist[c][eg][h];};
 void SetCRCZDCCFunHist(TH1D* const TH, Int_t const eg, Int_t const h) {this->fCRCZDCCFunHist[eg][h] = TH;};
 TH1D* GetCRCZDCCFunHist(Int_t const eg, Int_t const h) const {return this->fCRCZDCCFunHist[eg][h];};
 void SetCRCZDCSpectra(TProfile* const TH, Int_t const eg, Int_t const h) {this->fCRCZDCSpectra[eg][h] = TH;};
 TProfile* GetCRCZDCSpectra(Int_t const eg, Int_t const h) const {return this->fCRCZDCSpectra[eg][h];};
 
 // 12.c) Covariances:
 void SetCRCZDCCovHist(TH2D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCZDCCovHist[c][eg][h] = TH;};
 TH2D* GetCRCZDCCovHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCZDCCovHist[c][eg][h];};
 void SetCRCZDCCorrProdTempHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fCRCZDCCorrProdTempHist[c][eg][h] = TH;};
 TH1D* GetCRCZDCCorrProdTempHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fCRCZDCCorrProdTempHist[c][eg][h];};
 
 // CRC2:
 void SetCRC2nEtaBins(Int_t NB) {this->fCRC2nEtaBins = NB;};
 Int_t GetCRC2nEtaBins() {return this->fCRC2nEtaBins;};
 void SetCRC2CorPro(TProfile* const TP, Int_t const r, Int_t const h, Int_t const c) {this->fCRC2CorPro[r][h][c] = TP;};
 TProfile* GetCRC2CorPro(Int_t const r, Int_t const h, Int_t const c) const {return this->fCRC2CorPro[r][h][c];};
 void SetCRC2NUAPro(TProfile* const TP, Int_t const r, Int_t const h, Int_t const c, Int_t const e) {this->fCRC2NUAPro[r][h][c][e] = TP;};
 TProfile* GetCRC2NUAPro(Int_t const r, Int_t const h, Int_t const c, Int_t const e) const {return this->fCRC2NUAPro[r][h][c][e];};
 void SetCRC2CorHist(TH1D* const TP, Int_t const h, Int_t const c, Int_t const e) {this->fCRC2CorHist[h][c][e] = TP;};
 TH1D* GetCRC2CorHist(Int_t const h, Int_t const c, Int_t const e) const {return this->fCRC2CorHist[h][c][e];};
 void SetCRC2NUAHist(TH1D* const TP, Int_t const h, Int_t const c, Int_t const e) {this->fCRC2NUAHist[h][c][e] = TP;};
 TH1D* GetCRC2NUAHist(Int_t const h, Int_t const c, Int_t const e) const {return this->fCRC2NUAHist[h][c][e];};
 void SetCRC2CovPro(TProfile* const TP, Int_t const r, Int_t const h, Int_t const c) {this->fCRC2CovPro[r][h][c] = TP;};
 TProfile* GetCRC2CovPro(Int_t const r, Int_t const h, Int_t const c) const {return this->fCRC2CovPro[r][h][c];};
 void SetCRC2CovHist(TH1D* const TP, Int_t const r, Int_t const h, Int_t const c) {this->fCRC2CovHist[r][h][c] = TP;};
 TH1D* GetCRC2CovHist(Int_t const r, Int_t const h, Int_t const c) const {return this->fCRC2CovHist[r][h][c];};
 
 //  13.) CRC Pt differential
 void SetCRCPtTPCTNt(TNtuple* const TH) {this->fCRCPtTPCTNt = TH;};
 TNtuple* GetCRCPtTPCTNt() const {return this->fCRCPtTPCTNt;};
 void SetCRCPtVZTNt(TNtuple* const TH) {this->fCRCPtVZTNt = TH;};
 TNtuple* GetCRCPtVZTNt() const {return this->fCRCPtVZTNt;};
 void SetCRCPtZDCTNt(TNtuple* const TH) {this->fCRCPtZDCTNt = TH;};
 TNtuple* GetCRCPtZDCTNt() const {return this->fCRCPtZDCTNt;};
 
 void SetCenWeightsHist(TH1D* const n) {this->fCenWeightsHist = n;};
 TH1D* GetCenWeightsHist() const {return this->fCenWeightsHist;};
 void SetPtWeightsHist(TH1D* const n, Int_t c) {this->fPtWeightsHist[c] = n;};
 TH1D* GetPtWeightsHist(Int_t c) const {return this->fPtWeightsHist[c];};
 void SetEtaWeightsHist(TH1D* const n, Int_t h, Int_t b, Int_t c) {this->fEtaWeightsHist[h][b][c] = n;};
 TH1D* GetEtaWeightsHist(Int_t h, Int_t b, Int_t c) const {return this->fEtaWeightsHist[h][b][c];};
 void SetNvsCenCut(TH1D* const n, Int_t c, Int_t h) {this->fNvsCenCut[c][h] = n;};
 TH1D* GetNvsCenCut(Int_t c, Int_t h) const {return this->fNvsCenCut[c][h];};
 void SetZNCentroid(TH2F* const n, Int_t const eg, Int_t const h) {this->fhZNCentroid[eg][h] = n;};
 TH2F* GetZNCentroid(Int_t const eg, Int_t const h) const {return this->fhZNCentroid[eg][h];};
 void SetZNSpectra(TH1F* const n, Int_t const eg, Int_t const h) {this->fhZNSpectra[eg][h] = n;};
 TH1F* GetZNSpectra(Int_t const eg, Int_t const h) const {return this->fhZNSpectra[eg][h];};
 void SetZNCvsZNA(TH2F* const n, Int_t const h) {this->fhZNCvsZNA[h] = n;};
 TH2F* GetZNCvsZNA(Int_t const h) const {return this->fhZNCvsZNA[h];};
 void SetZNvsCen(TH2F* const n, Int_t const h) {this->fhZNvsCen[h] = n;};
 TH2F* GetZNvsCen(Int_t const h) const {return this->fhZNvsCen[h];};
  void SetZNvsTCen(TH2F* const n, Int_t const h) {this->fhZNvsTCen[h] = n;};
  TH2F* GetZNvsTCen(Int_t const h) const {return this->fhZNvsTCen[h];};
  void SetCenvsMul(TH2F* const n, Int_t const h) {this->fhCenvsMul[h] = n;};
  TH2F* GetCenvsMul(Int_t const h) const {return this->fhCenvsMul[h];};
  void SetCenvsDif(TH2F* const n, Int_t const h) {this->fhCenvsDif[h] = n;};
  TH2F* GetCenvsDif(Int_t const h) const {return this->fhCenvsDif[h];};
  void SetZNvsMul(TH2F* const n, Int_t const h) {this->fhZNvsMul[h] = n;};
  TH2F* GetZNvsMul(Int_t const h) const {return this->fhZNvsMul[h];};
 void SetZNCenvsMul(TH2F* const n, Int_t const eg, Int_t const h) {this->fhZNCenvsMul[eg][h] = n;};
 TH2F* GetZNCenvsMul(Int_t const eg, Int_t const h) const {return this->fhZNCenvsMul[eg][h];};
 void SetZNResvsMul(TH2F* const n, Int_t const eg, Int_t const h) {this->fhZNResvsMul[eg][h] = n;};
 TH2F* GetZNResvsMul(Int_t const eg, Int_t const h) const {return this->fhZNResvsMul[eg][h];};
 void SetZNResvsCen(TH2F* const n, Int_t const eg, Int_t const h) {this->fhZNResvsCen[eg][h] = n;};
 TH2F* GetZNResvsCen(Int_t const eg, Int_t const h) const {return this->fhZNResvsCen[eg][h];};
 void SetZNQVecCov(TProfile* const n, Int_t const h) {this->fhZNQVecCov[h] = n;};
 TProfile* GetZNQVecCov(Int_t const h) const {return this->fhZNQVecCov[h];};
  void SetZDCESEHistEP(TH2D* const n, Int_t const h) {this->fZDCESEHistEP[h] = n;};
  TH2D* GetZDCESEHistEP(Int_t const h) const {return this->fZDCESEHistEP[h];};
  void SetZDCESEHistQV(TH2D* const n, Int_t const h) {this->fZDCESEHistQV[h] = n;};
  TH2D* GetZDCESEHistQV(Int_t const h) const {return this->fZDCESEHistQV[h];};
  
 void SetPtDiffNBins(Int_t nbins) {this->fPtDiffNBins=nbins;}
 void SetPtDiffRangePt(Double_t min, Double_t max) {this->fPtDiffMinPt=min; this->fPtDiffMaxPt=max;}
 
 // Flow QC
 void SetFlowQCList(TList* const TL) {this->fFlowQCList = TL;};
 void SetFlowQCCorPro(TProfile* const TP, Int_t const c, Int_t const eg, Int_t const h) {this->fFlowQCCorPro[c][eg][h] = TP;};
 TProfile* GetFlowQCCorPro(Int_t const c, Int_t const eg, Int_t const h) const {return this->fFlowQCCorPro[c][eg][h];};
 void SetFlowQCNUAPro(TProfile* const TP, Int_t const c, Int_t const eg, Int_t const k) {this->fFlowQCNUAPro[c][eg][k] = TP;};
 TProfile* GetFlowQCNUAPro(Int_t const c, Int_t const eg, Int_t const k) const {return this->fFlowQCNUAPro[c][eg][k];};
 void SetFlowQCCorHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fFlowQCCorHist[c][eg][h] = TH;};
 TH1D* GetFlowQCCorHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fFlowQCCorHist[c][eg][h];};
 TProfile* GetFlowQCIntCorPro(Int_t const eg, Int_t const h) const {return this->fFlowQCIntCorPro[eg][h];};
 void SetFlowQCIntCorPro(TProfile* const TP, Int_t const eg, Int_t const k) {this->fFlowQCIntCorPro[eg][k] = TP;};
 TH1D* GetFlowQCIntCorHist(Int_t const eg, Int_t const h) const {return this->fFlowQCIntCorHist[eg][h];};
 void SetFlowQCIntCorHist(TH1D* const TP, Int_t const eg, Int_t const k) {this->fFlowQCIntCorHist[eg][k] = TP;};
 
 void SetFlowQCFinalPtDifHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fFlowQCFinalPtDifHist[c][eg][h] = TH;};
 TH1D* GetFlowQCFinalPtDifHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fFlowQCFinalPtDifHist[c][eg][h];};
 void SetFlowQCFinalPtIntHist(TH1D* const TH, Int_t const c, Int_t const eg) {this->fFlowQCFinalPtIntHist[c][eg] = TH;};
 TH1D* GetFlowQCFinalPtIntHist(Int_t const c, Int_t const eg) const {return this->fFlowQCFinalPtIntHist[c][eg];};
 void SetFlowQCSpectra(TH1D* const TH, Int_t const c) {this->fFlowQCSpectra[c] = TH;};
 TH1D* GetFlowQCSpectra(Int_t const c) const {return this->fFlowQCSpectra[c];};
 
 // Flow SP ZDC
 void SetFlowSPZDCList(TList* const TL) {this->fFlowSPZDCList = TL;};
  void SetFlowSPZDCRbRList(TList* const TL) {this->fFlowSPZDCRbRList = TL;};
  void SetFlowSPZDCRunsList(TList* const TL, Int_t r) {this->fFlowSPZDCRunsList[r] = TL;};
  
 void SetFlowSPZDCCorPro(TProfile* const TP, Int_t const r, Int_t const c, Int_t const eg, Int_t const h) {this->fFlowSPZDCCorPro[r][c][eg][h] = TP;};
 TProfile* GetFlowSPZDCCorPro(Int_t const r, Int_t const c, Int_t const eg, Int_t const h) const {return this->fFlowSPZDCCorPro[r][c][eg][h];};
  void SetFlowSPZDCCorNUA(TProfile* const TP, Int_t const r, Int_t const c, Int_t const h) {this->fFlowSPZDCCorNUA[r][c][h] = TP;};
  TProfile* GetFlowSPZDCCorNUA(Int_t const r, Int_t const c, Int_t const h) const {return this->fFlowSPZDCCorNUA[r][c][h];};
 void SetFlowSPZDCCorHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fFlowSPZDCCorHist[c][eg][h] = TH;};
 TH1D* GetFlowSPZDCCorHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fFlowSPZDCCorHist[c][eg][h];};
 
 void SetFlowSPZDCFinalPtDifHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fFlowSPZDCFinalPtDifHist[c][eg][h] = TH;};
 TH1D* GetFlowSPZDCFinalPtDifHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fFlowSPZDCFinalPtDifHist[c][eg][h];};
 void SetFlowSPZDCIntHist(TH1D* const TH, Int_t const c, Int_t const eg) {this->fFlowSPZDCIntHist[c][eg] = TH;};
 TH1D* GetFlowSPZDCIntHist(Int_t const c, Int_t const eg) const {return this->fFlowSPZDCIntHist[c][eg];};
 void SetFlowSPZDCSpectra(TH2F* const TH) {this->fFlowSPZDCSpectra = TH;};
 TH2F* GetFlowSPZDCSpectra() const {return this->fFlowSPZDCSpectra;};
 void SetFlowSPZDCIntPro(TProfile* const TP, Int_t const r, Int_t const c, Int_t const eg) {this->fFlowSPZDCIntPro[r][c][eg] = TP;};
 TProfile* GetFlowSPZDCIntPro(Int_t const r, Int_t const c, Int_t const eg) const {return this->fFlowSPZDCIntPro[r][c][eg];};
  void SetFlowSPZDCIntNUA(TProfile* const TP, Int_t const r, Int_t const eg) {this->fFlowSPZDCIntNUA[r][eg] = TP;};
  TProfile* GetFlowSPZDCIntNUA(Int_t const r, Int_t const eg) const {return this->fFlowSPZDCIntNUA[r][eg];};

 
 // Flow SP VZ
 void SetFlowSPVZList(TList* const TL) {this->fFlowSPVZList = TL;};
 void SetFlowSPVZCorPro(TProfile* const TP, Int_t const c, Int_t const eg, Int_t const h) {this->fFlowSPVZCorPro[c][eg][h] = TP;};
 TProfile* GetFlowSPVZCorPro(Int_t const c, Int_t const eg, Int_t const h) const {return this->fFlowSPVZCorPro[c][eg][h];};
 void SetFlowSPVZNUAPro(TProfile* const TP, Int_t const c, Int_t const eg, Int_t const h, Int_t const k) {this->fFlowSPVZNUAPro[c][eg][h][k] = TP;};
 TProfile* GetFlowSPVZNUAPro(Int_t const c, Int_t const eg, Int_t const h, Int_t const k) const {return this->fFlowSPVZNUAPro[c][eg][h][k];};
 void SetFlowSPVZCorHist(TH1D* const TH, Int_t const c, Int_t const eg, Int_t const h) {this->fFlowSPVZCorHist[c][eg][h] = TH;};
 TH1D* GetFlowSPVZCorHist(Int_t const c, Int_t const eg, Int_t const h) const {return this->fFlowSPVZCorHist[c][eg][h];};
 
 Int_t GetnRun() const {return this->fCRCnRun;};
 Int_t GetCRCPtnCen() const {return this->fCRCPtnCenBin;};
 Double_t GetCRCPtwCen() const {return this->fCRCPtwCenBin;};
 Int_t GetnEG() const {return this->fCRCnEtaGap;};
 Int_t GetCRCPtnPt() const {return this->fCRCPtnPtBin;};
 Double_t GetCRCPtwPt() const {return this->fCRCPtwPtBin;};
 Double_t GetCRCPtMinPt() const {return this->fCRCPtMinPt;};
 Int_t GetnCR() const {return this->fCRCnCR;};
 Double_t* GetCRCPtvarPtBins() const {return this->fCRCPtvarPtBins;};
 Int_t GetCRCZDCnCR() const {return this->fCRCZDCnCR;};
 Int_t GetCRCZDCnEtaBin() const {return this->fCRCZDCnEtaBin;};
 
 // CME:
 void SetCMEList(TList* const TL) {this->fCMEList = TL;};
 void SetCMERbRList(TList* const TL) {this->fCMERbRList = TL;};
 void SetCMETPCList(TList* const TL) {this->fCMETPCList = TL;};
 void SetCMEZDCList(TList* const TL) {this->fCMEZDCList = TL;};
 void SetCMERunsList(TList* const TL, Int_t r) {this->fCMERunsList[r] = TL;};
 // CME TPC only:
 void SetCMETPCCorPro(TProfile* const TP, Int_t const r, Int_t const eg, Int_t const h) {this->fCMETPCCorPro[r][eg][h] = TP;};
 TProfile* GetCMETPCCorPro(Int_t const r, Int_t const eg, Int_t const h) const {return this->fCMETPCCorPro[r][eg][h];};
 void SetCMETPCCovPro(TProfile* const TP, Int_t const r, Int_t const eg, Int_t const h) {this->fCMETPCCovPro[r][eg][h] = TP;};
 TProfile* GetCMETPCCovPro(Int_t const r, Int_t const eg, Int_t const h) const {return this->fCMETPCCovPro[r][eg][h];};
 void SetCMETPCNUAPro(TProfile* const TP, Int_t const r, Int_t const eg, Int_t const h) {this->fCMETPCNUAPro[r][eg][h] = TP;};
 TProfile* GetCMETPCNUAPro(Int_t const r, Int_t const eg, Int_t const h) const {return this->fCMETPCNUAPro[r][eg][h];};
 void SetCMETPCCorHist(TH1D* const TH, Int_t const eg, Int_t const h) {this->fCMETPCCorHist[eg][h] = TH;};
 TH1D* GetCMETPCCorHist(Int_t const eg, Int_t const h) const {return this->fCMETPCCorHist[eg][h];};
 void SetCMETPCCovHist(TH2D* const TH, Int_t const eg, Int_t const h) {this->fCMETPCCovHist[eg][h] = TH;};
 TH2D* GetCMETPCCovHist(Int_t const eg, Int_t const h) const {return this->fCMETPCCovHist[eg][h];};
 void SetCMETPCDistHist(TH1D* const TH, Int_t const eg, Int_t const h, Int_t const k) {this->fCMETPCDistHist[eg][h][k] = TH;};
 TH1D* GetCMETPCDistHist(Int_t const eg, Int_t const h, Int_t const k) const {return this->fCMETPCDistHist[eg][h][k];};
 // CME TPC-ZDCs:
 void SetCMEZDCCorPro(TProfile* const TP, Int_t const r, Int_t const eg, Int_t const h) {this->fCMEZDCCorPro[r][eg][h] = TP;};
 TProfile* GetCMEZDCCorPro(Int_t const r, Int_t const eg, Int_t const h) const {return this->fCMEZDCCorPro[r][eg][h];};
 void SetCMEZDCCovPro(TProfile* const TP, Int_t const r, Int_t const eg, Int_t const h) {this->fCMEZDCCovPro[r][eg][h] = TP;};
 TProfile* GetCMEZDCCovPro(Int_t const r, Int_t const eg, Int_t const h) const {return this->fCMEZDCCovPro[r][eg][h];};
 void SetCMEZDCNUAPro(TProfile* const TP, Int_t const r, Int_t const eg, Int_t const h) {this->fCMEZDCNUAPro[r][eg][h] = TP;};
 TProfile* GetCMEZDCNUAPro(Int_t const r, Int_t const eg, Int_t const h) const {return this->fCMEZDCNUAPro[r][eg][h];};
 void SetCMEZDCCorHist(TH1D* const TH, Int_t const eg, Int_t const h) {this->fCMEZDCCorHist[eg][h] = TH;};
 TH1D* GetCMEZDCCorHist(Int_t const eg, Int_t const h) const {return this->fCMEZDCCorHist[eg][h];};
 void SetCMEZDCCovHist(TH2D* const TH, Int_t const eg, Int_t const h) {this->fCMEZDCCovHist[eg][h] = TH;};
 TH2D* GetCMEZDCCovHist(Int_t const eg, Int_t const h) const {return this->fCMEZDCCovHist[eg][h];};
 void SetCMEZDCDistHist(TH1D* const TH, Int_t const eg, Int_t const h, Int_t const k) {this->fCMEZDCDistHist[eg][h][k] = TH;};
 TH1D* GetCMEZDCDistHist(Int_t const eg, Int_t const h, Int_t const k) const {return this->fCMEZDCDistHist[eg][h][k];};
 
 // 15.) Various
 void SetVariousList(TList* const Various) {this->fVariousList = Various;};
 void SetMultHist(TH1D* const TH) {this->fMultHist = TH;};
 TH1D* GetMultHist() const {return this->fMultHist;}
 void SetCenHist(TH1D* const TH) {this->fCenHist = TH;};
 TH1D* GetCenHist() const {return this->fCenHist;}
 void SetRunNumber(Int_t const n) {this->fRunNum = n;};
 Int_t GetRunNumber() const {return this->fRunNum;}
 void SetDataSet(DataSet set) {this->fDataSet = set;};
 DataSet GetDataSet() const {return this->fDataSet;}
 void SetMPolSelec(MagnetPol set) {this->fMPolSelec = set;};
 MagnetPol GetMPolSelec() const {return this->fMPolSelec;}
 void SetCorrWeightTPC(CorrelationWeights weights) {this->fCorrWeightTPC = weights;};
 CorrelationWeights GetCorrWeightTPC() const {return this->fCorrWeightTPC;};
 void SetCorrWeightVZ(CorrelationWeights weights) {this->fCorrWeightVZ = weights;};
 CorrelationWeights GetCorrWeightVZ() const {return this->fCorrWeightVZ;};
 void SetCorrWeightZDC(CorrelationWeights weights) {this->fCorrWeightZDC = weights;};
 CorrelationWeights GetCorrWeightZDC() const {return this->fCorrWeightZDC;};
 void SetMinMulZN(Int_t weights) {this->fMinMulZN = weights;};
 Int_t GetMinMulZN() const {return this->fMinMulZN;};
 void SetMaxDevZN(Float_t weights) {this->fMaxDevZN = weights;};
 Float_t GetMaxDevZN() const {return this->fMaxDevZN;};
 
private:
 
 AliFlowAnalysisCRC(const AliFlowAnalysisCRC& afawQc);
 AliFlowAnalysisCRC& operator=(const AliFlowAnalysisCRC& afawQc);
 
 // 0.) base:
 TList* fHistList; //! base list to hold all output object
 
 // 1.) common:
 Bool_t fBookOnlyBasicCCH; // book only basis common control histrograms (by default book them all)
 AliFlowCommonHist *fCommonHists; //! common control histograms (taking into account ALL events)
 AliFlowCommonHist *fCommonHists2nd; //! common control histograms (taking into account only the events with 2 and more particles)
 AliFlowCommonHist *fCommonHists4th; //! common control histograms (taking into account only the events with 4 and more particles)
 AliFlowCommonHist *fCommonHists6th; //! common control histograms (taking into account only the events with 6 and more particles)
 AliFlowCommonHist *fCommonHists8th; //! common control histograms (taking into account only the events with 8 and more particles)
 AliFlowCommonHistResults *fCommonHistsResults2nd; //! final results for 2nd order int. and diff. flow for events with 2 and more particles
 AliFlowCommonHistResults *fCommonHistsResults4th; //! final results for 4th order int. and diff. flow for events with 4 and more particles
 AliFlowCommonHistResults *fCommonHistsResults6th; //! final results for 6th order int. and diff. flow for events with 6 and more particles
 AliFlowCommonHistResults *fCommonHistsResults8th; //! final results for 8th order int. and diff. flow for events with 8 and more particles
 Int_t fnBinsPhi; // number of phi bins
 Double_t fPhiMin; // minimum phi
 Double_t fPhiMax; // maximum phi
 Double_t fPhiBinWidth; // bin width for phi histograms
 Int_t fnBinsPt; // number of pt bins
 Double_t fPtMin; // minimum pt
 Double_t fPtMax; // maximum pt
 Double_t fPtBinWidth; // bin width for pt histograms
 Int_t fnBinsEta; // number of eta bins
 Double_t fEtaMin; // minimum eta
 Double_t fEtaMax; // maximum eta
 Double_t fEtaBinWidth; // bin width for eta histograms
 TProfile *fCommonConstants; //! profile to hold common constants
 Bool_t fFillMultipleControlHistograms; // fill separately control histos for events with >= 2, 4, 6 and 8 particles
 Int_t fHarmonic; // harmonic
 TString *fAnalysisLabel; //! analysis label (all histograms and output file will have this label)
 Bool_t fPrintFinalResults[4]; // print on the screen the final results (0=RF, 1=RP, 2=POI, 3=RF rebinned in M)
 
 // 2a.) particle weights:
 TList *fWeightsList; // list to hold all histograms with particle weights: fUseParticleWeights, fPhiWeights, fPtWeights and fEtaWeights
 Bool_t fUsePhiWeights; // use phi weights
 Bool_t fUsePtWeights; // use pt weights
 Bool_t fUseEtaWeights; // use eta weights
 Bool_t fUseTrackWeights; // use track weights (e.g. VZERO sector weights)
 Bool_t fUsePhiEtaWeights; // use phi,eta weights
 TProfile *fUseParticleWeights; //! profile with three bins to hold values of fUsePhiWeights, fUsePtWeights and fUseEtaWeights
 // TH1F *fPhiWeightsPOIs[2]; //! histogram holding phi weights
 // TH1D *fPtWeightsPOIs[2]; //! histogram holding pt weights
 // TH1D *fEtaWeightsPOIs[2]; //! histogram holding eta weights
 // TH2D *fPhiEtaWeightsPOIs[2]; //! histogram holding phi,eta weights
 TH1F *fPhiWeightsRPs; //!
 // TH1D *fPtWeightsRPs; //!
 // TH1D *fEtaWeightsRPs; //!
 // TH1F *fPhiDistrRefPOIs[2]; //! histogram holding phi weights
 // TH1D *fPtDistrRefPOIs[2]; //! histogram holding pt weights
 // TH1D *fEtaDistrRefPOIs[2]; //! histogram holding eta weights
 // TH2D *fPhiEtaDistrRefPOIs[2]; //! histogram holding phi,eta weights
 // TH1F *fPhiDistrRefRPs; //!
 // TH1D *fPtDistrRefRPs; //!
 // TH1D *fEtaDistrRefRPs; //!
 // TH2D *fPhiEtaDistrRefRPs; //!
 // TH1D *fPtWeights[2]; //! histogram holding pt weights
 TH2D *fPhiEtaWeights[2]; //!
 
 // 2b.) event weights:
 TString *fMultiplicityWeight; //! event-by-event weights for multiparticle correlations
 AliFlowCommonConstants::ERefMultSource fMultiplicityIs; // by default "kRP"
 
 // 3.) integrated flow
 //  3a.) lists:
 TList *fIntFlowList; //! list to hold all histograms and profiles relevant for integrated flow
 TList *fIntFlowProfiles; //! list to hold all profiles relevant for integrated flow
 TList *fIntFlowResults; //! list to hold all histograms with final results relevant for integrated flow
 TList *fIntFlowAllCorrelationsVsM; //! list to hold all profiles with correlations vs M
 //  3b.) flags:
 TProfile *fIntFlowFlags; //! profile to hold all flags for integrated flow
 Bool_t fApplyCorrectionForNUA; // apply correction for non-uniform acceptance
 Bool_t fApplyCorrectionForNUAVsM; // apply correction for non-uniform acceptance versus M
 Int_t fnBinsMult; // number of multiplicity bins for flow analysis versus multiplicity
 Double_t fMinMult; // minimal multiplicity for flow analysis versus multiplicity
 Double_t fMaxMult; // maximal multiplicity for flow analysis versus multiplicity
 Bool_t fPropagateErrorAlsoFromNIT; // propagate error by taking into account also non-isotropic terms (not sure if resulting error then is correct - to be improved)
 Bool_t fCalculateCumulantsVsM; // calculate cumulants versus multiplicity
 Bool_t fCalculateAllCorrelationsVsM; // calculate all correlations versus multiplicity
 Bool_t fMinimumBiasReferenceFlow; // store as reference flow in AliFlowCommonHistResults the minimum bias result (kFALSE by default)
 Bool_t fForgetAboutCovariances; // when propagating error forget about the covariances
 Bool_t fStoreVarious; // store phi distribution for one event to illustrate flow
 Int_t fExactNoRPs; // when shuffled, select only this number of RPs for the analysis
 Bool_t fUse2DHistograms; // use TH2D instead of TProfile to improve numerical stability in reference flow calculation
 Bool_t fFillProfilesVsMUsingWeights; // if the width of multiplicity bin is 1, weights are not needed
 Bool_t fUseQvectorTerms; // use TH2D with separate Q-vector terms instead of TProfile to improve numerical stability in reference flow calculation
 
 //  3c.) event-by-event quantities:
 TMatrixD *fReQ; //! fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
 TMatrixD *fImQ; //! fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
 TMatrixD *fSpk; //! fSM[p][k] = (sum_{i=1}^{M} w_{i}^{k})^{p+1}
 TH1D *fIntFlowCorrelationsEBE; //! 1st bin: <2>, 2nd bin: <4>, 3rd bin: <6>, 4th bin: <8>
 TH1D *fIntFlowEventWeightsForCorrelationsEBE; //! 1st bin: eW_<2>, 2nd bin: eW_<4>, 3rd bin: eW_<6>, 4th bin: eW_<8>
 TH1D *fIntFlowCorrelationsAllEBE; //! to be improved (add comment)
 TH1D *fIntFlowCorrectionTermsForNUAEBE[2]; //! [0=sin terms,1=cos terms], NUA = non-uniform acceptance
 TH1D *fIntFlowEventWeightForCorrectionTermsForNUAEBE[2]; //! [0=sin terms,1=cos terms], NUA = non-uniform acceptance
 Double_t fNumberOfRPsEBE; // # of Reference Particles
 Double_t fNumberOfPOIsEBE; // # of Particles of Interest
 Double_t fReferenceMultiplicityEBE; // reference multiplicity
 Double_t fCentralityEBE; // centrality percentile
 Double_t fCentralityVarEBE; // centrality (alternative estimation) percentile
 //  3d.) profiles:
 TProfile *fAvMultiplicity; //! profile to hold average multiplicities and number of events for events with nRP>=0, nRP>=1, ... , and nRP>=8
 TProfile *fIntFlowCorrelationsPro; //! average correlations <<2>>, <<4>>, <<6>> and <<8>> (with wrong errors!)
 TProfile *fIntFlowSquaredCorrelationsPro; //! average correlations squared <<2>^2>, <<4>^2>, <<6>^2> and <<8>^2>
 TProfile *fIntFlowCorrelationsVsMPro[4]; //! average correlations <<2>>, <<4>>, <<6>> and <<8>> versus multiplicity (error is wrong here!)
 TProfile *fIntFlowSquaredCorrelationsVsMPro[4]; //! average correlations <<2>^2>, <<4>^2>, <<6>^2> and <<8>^2> versus multiplicity
 TProfile *fIntFlowCorrelationsAllPro; //! average all correlations for integrated flow (with wrong errors!)
 TProfile *fIntFlowCorrelationsAllVsMPro[64]; //! average all correlations vs M (errors via Sumw2 - to me improved)
 TProfile *fIntFlowExtraCorrelationsPro; //! when particle weights are used some extra correlations appear
 TProfile *fIntFlowProductOfCorrelationsPro; //! average product of correlations <2>, <4>, <6> and <8>
 TProfile *fIntFlowProductOfCorrelationsVsMPro[6]; //! average product of correlations <2>, <4>, <6> and <8>
 // [0=<<2><4>>,1=<<2><6>>,2=<<2><8>>,3=<<4><6>>,4=<<4><8>>,5=<<6><8>>]
 TProfile *fIntFlowProductOfCorrectionTermsForNUAPro; //! average product of correction terms for NUA
 TProfile *fIntFlowCorrectionTermsForNUAPro[2]; //! average correction terms for non-uniform acceptance (with wrong errors!) [0=sin terms,1=cos terms]
 TProfile *fIntFlowCorrectionTermsForNUAVsMPro[2][4]; //! average correction terms for non-uniform acceptance (with wrong errors!) [0=sin terms,1=cos terms][correction term index] vs multiplicity
 
 //  3e.) histograms with final results:
 TH1D *fIntFlowCorrelationsHist; //! final results for average correlations <<2>>, <<4>>, <<6>> and <<8>> (with correct errors!)
 TH1D *fIntFlowCorrelationsVsMHist[4]; //! average correlations <<2>>, <<4>>, <<6>> and <<8>> versus multiplicity (error is correct here!)
 TH1D *fIntFlowCorrelationsAllHist; //! final results for all average correlations (with correct errors!)
 TH1D *fIntFlowCorrectionTermsForNUAHist[2]; //! final results for correction terms for non-uniform acceptance (with correct errors!) [0=sin terms,1=cos terms]
 TH1D *fIntFlowCovariances; //! final result for covariances of correlations (multiplied with weight dependent prefactor)
 TH1D *fIntFlowSumOfEventWeights[2]; //! sum of linear and quadratic event weights for <2>, <4>, <6> and <8>: [0=linear 1,1=quadratic]
 TH1D *fIntFlowSumOfProductOfEventWeights; //! sum of products of event weights for correlations <2>, <4>, <6> and <8>
 TH1D *fIntFlowCovariancesVsM[6]; //! final result for covariances of correlations (multiplied with weight dependent prefactor) versus M
 // [0=Cov(2,4),1=Cov(2,6),2=Cov(2,8),3=Cov(4,6),4=Cov(4,8),5=Cov(6,8)]
 TH1D *fIntFlowSumOfEventWeightsVsM[4][2]; //! sum of linear and quadratic event weights for <2>, <4>, <6> and <8> versum multiplicity
 // [0=sum{w_{<2>}},1=sum{w_{<4>}},2=sum{w_{<6>}},3=sum{w_{<8>}}][0=linear 1,1=quadratic]
 TH1D *fIntFlowSumOfProductOfEventWeightsVsM[6]; //! sum of products of event weights for correlations <2>, <4>, <6> and <8> vs M
 // [0=sum{w_{<2>}w_{<4>}},1=sum{w_{<2>}w_{<6>}},2=sum{w_{<2>}w_{<8>}},
 //  3=sum{w_{<4>}w_{<6>}},4=sum{w_{<4>}w_{<8>}},5=sum{w_{<6>}w_{<8>}}]
 TH1D *fIntFlowCovariancesNUA; //! final result for covariances of all terms needed for NUA (multiplied with weight dependent prefactor)
 TH1D *fIntFlowSumOfEventWeightsNUA[2][2]; //! sum of linear and quadratic event weights for NUA terms: [0=sin,1=cos][0=linear 1,1=quadratic]
 TH1D *fIntFlowSumOfProductOfEventWeightsNUA; //! sum of products of event weights for NUA terms
 TH1D *fIntFlowQcumulants; //! final results for integrated Q-cumulants QC{2}, QC{4}, QC{6} and QC{8}
 TH1D *fIntFlowQcumulantsVsM[4]; //! final results for integrated Q-cumulants QC{2}, QC{4}, QC{6} and QC{8} versus multiplicity
 TH1D *fIntFlowQcumulantsRebinnedInM; //! final results for reference Q-cumulants QC{2}, QC{4}, QC{6} and QC{8} rebinned in M
 TH1D *fIntFlowQcumulantsErrorSquaredRatio; //! ratio between error squared: with/without non-isotropic terms
 TH1D *fIntFlow; //! final results for integrated flow estimates v_n{2,QC}, v_n{4,QC}, v_n{6,QC} and v_n{8,QC}
 TH1D *fIntFlowVsM[4]; //! final results for integrated flow estimates v_n{2,QC}, v_n{4,QC}, v_n{6,QC} and v_n{8,QC} versus multiplicity
 TH1D *fIntFlowRebinnedInM; //! final results for ref. flow estimates v_n{2,QC}, v_n{4,QC}, v_n{6,QC} and v_n{8,QC} rebinned in M
 TH1D *fIntFlowDetectorBias; //! bias coming from detector inefficiencies to <<2>>, <<4>>, <<6>> and <<8>> (corrected/measured)
 TH1D *fIntFlowDetectorBiasVsM[4]; //! bias coming from detector inefficiencies to <<2>>, <<4>>, <<6>> and <<8>> vs M (corrected/measured)
 // 4.) differential flow
 //  4a.) lists:
 TList *fDiffFlowList; //! list to hold list with all histograms (fDiffFlowResults) and list with profiles (fDiffFlowProfiles) relevant for differential flow
 TList *fDiffFlowProfiles; //! list to hold all profiles relevant for differential flow
 TList *fDiffFlowResults; //! list to hold all histograms with final results relevant for differential flow
 TList *fDiffFlow2D; //! list to hold all objects relevant for 2D differential flow
 //    4aa.) nested list in list fDiffFlowProfiles:
 TList *fDiffFlowCorrelationsProList[2][2]; //! list to hold profiles with all correlations for differential flow [0=RP,1=POI][0=pt,1=eta]
 TList *fDiffFlowProductOfCorrelationsProList[2][2]; //! list to hold profiles with products of all correlations for differential flow [0=RP,1=POI][0=pt,1=eta]
 TList *fDiffFlowCorrectionsProList[2][2]; //! list to hold profiles with correction term for NUA for differential flow [0=RP,1=POI][0=pt,1=eta]
 TList *f2DDiffFlowCorrelationsProList[2]; //! list to hold profiles with all correlations for 2D differential flow [0=RP,1=POI]
 //    4ab.) nested list in list fDiffFlowResults:
 TList *fDiffFlowCorrelationsHistList[2][2]; //! list to hold histograms with all correlations for differential flow [0=RP,1=POI][0=pt,1=eta]
 TList *fDiffFlowSumOfEventWeightsHistList[2][2][2]; //! list to hold histograms with sum of linear/quadratic event weights [0=RP,1=POI][0=pt,1=eta][0=linear 1,1=quadratic]
 TList *fDiffFlowSumOfProductOfEventWeightsHistList[2][2]; //! list to hold histograms with sum of products of event weights [0=RP,1=POI][0=pt,1=eta]
 TList *fDiffFlowCorrectionsHistList[2][2]; //! list to hold histograms with correction term for NUA for differential flow [0=RP,1=POI][0=pt,1=eta]
 TList *fDiffFlowCovariancesHistList[2][2]; //! list to hold histograms with all covariances for differential flow [0=RP,1=POI][0=pt,1=eta]
 TList *fDiffFlowCumulantsHistList[2][2]; //! list to hold histograms with all cumulants for differential flow [0=RP,1=POI][0=pt,1=eta]
 TList *fDiffFlowDetectorBiasHistList[2][2]; //! list to hold histograms which quantify detector bias to differential cumulants [0=RP,1=POI][0=pt,1=eta]
 TList *fDiffFlowHistList[2][2]; //! list to hold histograms with final results for differential flow [0=RP,1=POI][0=pt,1=eta]
 //  4b.) flags:
 TProfile *fDiffFlowFlags; //! profile to hold all flags for differential flow
 Bool_t fCalculateDiffFlow; // if you set kFALSE only reference flow will be calculated
 Bool_t fCalculate2DDiffFlow; // calculate 2D differential flow vs (pt,eta) (Remark: this is expensive in terms of CPU time)
 Bool_t fCalculateDiffFlowVsEta; // if you set kFALSE only differential flow vs pt is calculated
 //  4c.) event-by-event quantities:
 //   1D:
 TProfile *fReRPQ1dEBE[3][2][4][9]; //! real part [0=r,1=p,2=q][0=pt,1=eta][m][k]
 TProfile *fImRPQ1dEBE[3][2][4][9]; //! imaginary part [0=r,1=p,2=q][0=pt,1=eta][m][k]
 TProfile *fs1dEBE[3][2][9]; //! [0=r,1=p,2=q][0=pt,1=eta][k] // to be improved
 TH1D *fDiffFlowCorrelationsEBE[2][2][4]; //! [0=RP,1=POI][0=pt,1=eta][reduced correlation index]
 TH1D *fDiffFlowEventWeightsForCorrelationsEBE[2][2][4]; //! [0=RP,1=POI][0=pt,1=eta][event weights for reduced correlation index]
 TH1D *fDiffFlowCorrectionTermsForNUAEBE[2][2][2][10]; //! [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correction term index]
 
 //   2D:
 TProfile2D *fReRPQ2dEBE[3][4][9]; //! real part of r_{m*n,k}(pt,eta), p_{m*n,k}(pt,eta) and q_{m*n,k}(pt,eta)
 TProfile2D *fImRPQ2dEBE[3][4][9]; //! imaginary part of r_{m*n,k}(pt,eta), p_{m*n,k}(pt,eta) and q_{m*n,k}(pt,eta)
 TProfile2D *fs2dEBE[3][9]; //! [t][k] // to be improved
 //  4d.) profiles:
 //   1D:
 TProfile *fDiffFlowCorrelationsPro[2][2][4]; //! [0=RP,1=POI][0=pt,1=eta][correlation index]
 TProfile *fDiffFlowSquaredCorrelationsPro[2][2][4]; //! [0=RP,1=POI][0=pt,1=eta][correlation index]
 TProfile *fDiffFlowProductOfCorrelationsPro[2][2][8][8]; //! [0=RP,1=POI][0=pt,1=eta] [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>] x
 //                          [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>]
 TProfile *fDiffFlowCorrectionTermsForNUAPro[2][2][2][10]; //! [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correction term index]
 //   2D:
 TProfile2D *f2DDiffFlowCorrelationsPro[2][4]; //! [0=RP,1=POI][correlation index]
 //   Other differential correlators:
 TList *fOtherDiffCorrelatorsList; //! list to hold profiles with other differential correlators
 TProfile *fOtherDiffCorrelators[2][2][2][1]; //! [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correlator index]
 //  4e.) histograms holding final results:
 //   1D:
 TH1D *fDiffFlowCorrelationsHist[2][2][4]; //! [0=RP,1=POI][0=pt,1=eta][correlation index]
 TH1D *fDiffFlowCovariances[2][2][5]; //! [0=RP,1=POI][0=pW not used,1=pW used][0=exact eW,1=non-exact eW][0=pt,1=eta][index of covariances]
 TH1D *fDiffFlowCumulants[2][2][4]; //! [0=RP,1=POI][0=pt,1=eta][0=QC{2'},1=QC{4'},2=QC{6'},3=QC{8'}]
 TH1D *fDiffFlowDetectorBias[2][2][4]; //! [0=RP,1=POI][0=pt,1=eta][0=gQC{2'}/QC{2'},1=gQC{4'}/QC{4'},2=gQC{6'}/QC{6'},3=gQC{8'}/QC{8'}]
 TH1D *fDiffFlow[2][2][4]; //! [0=RP,1=POI][0=pt,1=eta][0=v'{2},1=v'{4},2=v'{6},3=v'{8}]
 TH1D *fDiffFlowSumOfEventWeights[2][2][2][4]; //! [0=RP,1=POI][0=pt,1=eta][0=linear 1,1=quadratic][0=<2'>,1=<4'>,2=<6'>,3=<8'>]
 TH1D *fDiffFlowSumOfProductOfEventWeights[2][2][8][8]; //! [0=RP,1=POI][0=pt,1=eta]  [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>] x
 //                           [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>]
 TH1D *fDiffFlowCorrectionTermsForNUAHist[2][2][2][10]; //! [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correction term index]
 //   2D:
 TH2D *f2DDiffFlowCumulants[2][4]; //! 2D differential cumulants [0=RP,1=POI][cumulant order]
 TH2D *f2DDiffFlow[2][4]; //! 2D differential flow [0=RP,1=POI][cumulants order]
 // 6.) distributions:
 TList *fDistributionsList; //! list to hold all distributions of correlations
 TProfile *fDistributionsFlags; //! profile to hold all flags for distributions of correlations
 Bool_t fStoreDistributions; // store or not distributions of correlations
 TH1D *fDistributions[4]; //! [0=distribution of <2>,1=distribution of <4>,2=distribution of <6>,3=distribution of <8>]
 Int_t fnBinsForCorrelations; // # of bins for correlation axis in fDistributions[4], fCorrelation2468VsMult[4] and fCorrelationProduct2468VsMult[1]
 Double_t fMinValueOfCorrelation[4]; // min values of <2>, <4>, <6> and <8>
 Double_t fMaxValueOfCorrelation[4]; // max values of <2>, <4>, <6> and <8>
 Double_t fMinValueOfCorrelationProduct[1]; // min values of <2><4>, <2><6>, <2><8>, <4><6> etc. TBI add the other ones when needed first time
 Double_t fMaxValueOfCorrelationProduct[1]; // max values of <2><4>, <2><6>, <2><8>, <4><6> etc. TBI add the other ones when needed first time
 Double_t fMinValueOfQvectorTerms[4]; // MinValueOfQvectorTerms
 Double_t fMaxValueOfQvectorTerms[4]; // MaxValueOfQvectorTerms
 
 // 8.) debugging and cross-checking:
 TList *fNestedLoopsList; //! list to hold all profiles filled with nested loops
 Bool_t fEvaluateIntFlowNestedLoops; // evaluate nested loops relevant for integrated flow
 Bool_t fEvaluateDiffFlowNestedLoops; // evaluate nested loops relevant for differential flow
 Int_t fMaxAllowedMultiplicity; // nested loops will be evaluated only for events with multiplicity <= fMaxAllowedMultiplicity
 TProfile *fEvaluateNestedLoops; //! profile with four bins: fEvaluateIntFlowNestedLoops, fEvaluateDiffFlowNestedLoops, fCrossCheckInPtBinNo and fCrossCheckInEtaBinNo
 // integrated flow:
 TProfile *fIntFlowDirectCorrelations; //! multiparticle correlations relevant for int. flow calculated with nested loops
 TProfile *fIntFlowExtraDirectCorrelations; //! when particle weights are used some extra correlations appear
 TProfile *fIntFlowDirectCorrectionTermsForNUA[2]; //! average correction terms for non-uniform acceptance evaluated with nested loops [0=sin terms,1=cos terms]
 // differential flow:
 Int_t fCrossCheckInPtBinNo; // cross-check results for reduced correlations and corrections in this pt bin
 Int_t fCrossCheckInEtaBinNo; // cross-check results for reduced correlations and corrections in this eta bin
 TH1D *fNoOfParticlesInBin; // bin: 1 = # of RPs in pt bin, 2 = # of RPs in eta bin, 3 = # of POIs in pt bin, 4 = # of POIs in eta bin
 TProfile *fDiffFlowDirectCorrelations[2][2][4]; //! [0=RP,1=POI][0=pt,1=eta][correlation index]
 TProfile *fDiffFlowDirectCorrectionTermsForNUA[2][2][2][10]; //! [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correction term index]
 // other differential correlators:
 TProfile *fOtherDirectDiffCorrelators[2][2][2][1]; //! [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correlator index]
 // mixed harmonics:
 TProfile *fMixedHarmonicsNestedLoops; //! Cross-check mixed harmonics with nested loops.
 
 // 9.) mixed harmonics:
 //  9a.) lists:
 TList *fMixedHarmonicsList; //! list to hold all histograms and profiles for mixed harmonics
 TList *fMixedHarmonicsProfiles; //! list to hold all profiles for mixed harmonics
 TList *fMixedHarmonicsResults; //! list to hold all histograms with final results for mixed harmonics
 TList *fMixedHarmonicsErrorPropagation; //! list to hold all objects needed for statistical error propagation
 //TList *fIntFlowAllCorrelationsVsM; // list to hold all profiles with correlations vs M
 //  9b.) flags:
 TProfile *fMixedHarmonicsFlags; //! profile to hold all flags for mixed harmonics
 Bool_t fCalculateMixedHarmonics; // calculate or not mixed harmonics
 Bool_t fCalculateMixedHarmonicsVsM; // calculate or not mixed harmonics vs multiplicity
 //  9c.) profiles:
 TProfile *f2pCorrelations; //! profile to hold all 2-particle correlations
 TProfile *f3pCorrelations; //! profile to hold all 3-particle correlations
 TProfile *f4pCorrelations; //! profile to hold all 4-particle correlations
 TProfile *f5pCorrelations; //! profile to hold all 5-particle correlations
 TProfile *f6pCorrelations; //! profile to hold all 6-particle correlations
 TProfile *f7pCorrelations; //! profile to hold all 7-particle correlations
 TProfile *f8pCorrelations; //! profile to hold all 8-particle correlations
 //  9d.) results:
 TH1D *f2pCumulants; //! histogram to hold all 2-particle cumulants
 TH1D *f3pCumulants; //! histogram to hold all 3-particle cumulants
 TH1D *f4pCumulants; //! histogram to hold all 4-particle cumulants
 TH1D *f5pCumulants; //! histogram to hold all 5-particle cumulants
 TH1D *f6pCumulants; //! histogram to hold all 6-particle cumulants
 TH1D *f7pCumulants; //! histogram to hold all 7-particle cumulants
 TH1D *f8pCumulants; //! histogram to hold all 8-particle cumulants
 //  9e.) statistical error propagation:
 TH1D *fMixedHarmonicEventWeights[2]; //! sum of linear and quadratic event weights for mixed harmonics => [0=linear 1,1=quadratic]
 TH2D *fMixedHarmonicProductOfEventWeights; //! sum of products of event weights for mixed harmonics
 TProfile2D *fMixedHarmonicProductOfCorrelations; //! averages of products of mixed harmonics correlations
 // 10.) Control histograms:
 //  10a.) list:
 TList *fControlHistogramsList; // list to hold all control histograms
 //  10b.) flags:
 TProfile *fControlHistogramsFlags; // profile to hold all flags for control histograms
 Bool_t fStoreControlHistograms; // store or not control histograms
 //  10c.) histograms:
 TH2D *fCorrelationNoRPsVsRefMult; //! correlation between # RPs and ref. mult. determined centrally
 TH2D *fCorrelationNoPOIsVsRefMult; //! correlation between # POIs and ref. mult. determined centrally
 TH2D *fCorrelationNoRPsVsNoPOIs; //! correlation between # RPs and # POIs
 TH2D *fCorrelation2468VsMult[4]; //! <2>, <4>, <6> and <8> vs multiplicity (#RPs, #POIs or external)
 TH2D *fCorrelationProduct2468VsMult[1]; //! <2><4>, <2><6>, <2><8>, <4><6> etc. vs multiplicity (#RPs, #POIs or external)
 TH2D *fQvectorTermsVsMult[4]; //! |Qn|^2/M, |Q2n|^2/M, |Qn|^4/(M(2M-1)), Re[Q2nQn^*Qn^*]/M, ... vs multiplicity (#RPs, #POIs or external)
 // 11.) Bootstrap:
 //  11a) lists:
 TList *fBootstrapList; //! list to hold all output objects for bootstrap
 TList *fBootstrapProfilesList; //! list to hold all profiles for bootstrap
 TList *fBootstrapResultsList; //! list to hold all histograms for bootstrap
 //  11b) flags:
 TProfile *fBootstrapFlags; //! profile to hold all flags for mixed harmonics
 Bool_t fUseBootstrap; // use bootstrap to estimate statistical spread
 Bool_t fUseBootstrapVsM; // use bootstrap to estimate statistical spread for results vs M
 Int_t fnSubsamples; // number of subsamples (SS), by default 10
 TRandom3 *fRandom; //! local random generator
 //  11c) profiles:
 TProfile2D *fBootstrapCorrelations; //! x-axis => <2>, <4>, <6>, <8>; y-axis => subsample #
 TProfile2D *fBootstrapCorrelationsVsM[4]; //! index => <2>, <4>, <6>, <8>; x-axis => multiplicity; y-axis => subsample #
 //  11d) histograms:
 TH2D *fBootstrapCumulants; //! x-axis => QC{2}, QC{4}, QC{6}, QC{8}; y-axis => subsample #
 TH2D *fBootstrapCumulantsVsM[4]; //! index => QC{2}, QC{4}, QC{6}, QC{8}; x-axis => multiplicity; y-axis => subsample #
 
 // 12.) CRC
 
 TList *fCRCList; //! list to hold CRC histograms
 TList *fTempList; //! list to hold temp histograms
 TProfile *fCRCFlags; //! profile to hold all flags for CRC
 Bool_t fCalculateCRC; // calculate CRC
 Bool_t fCalculateCRCPt;
 Bool_t fCalculateCME;
 Bool_t fCalculateCRC2;
 Bool_t fCalculateCRCVZ;
 Bool_t fCalculateCRCZDC;
 Bool_t fCalculateFlowQC;
 Bool_t fCalculateFlowZDC;
 Bool_t fCalculateFlowVZ;
 Bool_t fUseVZERO;
 Bool_t fUseZDC;
 Bool_t fRecenterZDC;
 Bool_t fNUAforCRC;
 Bool_t fUseCRCRecenter;
 Bool_t fDivSigma;
 Bool_t fInvertZDC;
 Bool_t fCRCTestSin;
 Double_t fCRCEtaMin;
 Double_t fCRCEtaMax;
 Int_t fRunNum;
 Int_t fCachedRunNum;
 Int_t fRunBin;
 Int_t fCenBin;
 CorrelationWeights fCorrWeightTPC;
 CorrelationWeights fCorrWeightVZ;
 CorrelationWeights fCorrWeightZDC;
 
  TList *fCRCIntList; //! list to hold CRC histograms
  const static Int_t fCRCnCR = 16;
  const static Int_t fCRCnNUA = 14;
  const static Int_t fCRCnEtaGap = 7;
  const static Int_t fCRCnCorr = 4;
  const static Int_t fCRCnEtaBins = 32;
  const static Int_t fCRCMaxnCen = 10;
  Int_t fCRCnCen;
  Int_t fCRCCenBinWidth;
  const static Int_t fCRCnHar = 3;
  const static Int_t fCRCMaxnRun = 211;
  
  TList *fCRCIntRbRList; //! CRC list of histograms RbR
  TList *fCRCIntRunsList[fCRCMaxnRun]; //! list of runs
  TH1D *fCRCQRe[2][fCRCnHar]; //! real part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCRCQIm[2][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCRCMult[2][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][k]
  TProfile *fCRCCorrPro[fCRCMaxnRun][fCRCnCorr][fCRCnEtaGap][fCRCMaxnCen]; //! correlation profile, [CRCBin][eg]
  TH1D *fCRCSumWeigHist[fCRCMaxnRun][fCRCnCorr][fCRCnEtaGap][fCRCMaxnCen]; //! correlation weights histo, [CRCBin][eg]
  TProfile *fCRCNUATermsPro[fCRCMaxnRun][fCRCnNUA][fCRCnEtaGap][fCRCMaxnCen]; //! NUA terms profile
  
  TH1D *fCRCCorrProdTempHist[fCRCnCorr][fCRCnEtaGap][fCRCMaxnCen]; //! temporary correlation products for covariances, [CRCBin][eg]
  TH1D *fCRCCorrHist[fCRCnCorr][fCRCnEtaGap][fCRCMaxnCen]; //! <<2'>>, [CRCBin][eg]
  TH1D *fCRCCumHist[fCRCnCorr][fCRCnEtaGap][fCRCMaxnCen]; //! QC{2}, [CRCBin][eg]
  TProfile *fCRCCorrProd2p2pPro[fCRCnCorr][fCRCnEtaGap][fCRCMaxnCen]; //! correlation products
  TH1D *fCRCCovHist[fCRCnCorr][fCRCnEtaGap][fCRCMaxnCen]; //! covariances final histo
  TH1D *fCRCCFunHist[fCRCnEtaGap][fCRCMaxnCen]; //! correlation function histo, [CRCBin][eg]
  TH1D *fCRCNUATermsHist[fCRCnNUA][fCRCnEtaGap][fCRCMaxnCen]; //! NUA terms final histo
 
 // Q vectors
 const static Int_t fCRCQVecnCR = 64;
 Int_t fCRCnRun;
 DataSet fDataSet;
 MagnetPol fMPolSelec;
 TArrayI fRunList;    // Run list
 TList *fCRCQVecList; //! Q Vectors list
 TList *fCRCQVecListRun[fCRCMaxnRun]; //! Q Vectors list per run
 TList *fCRCQVecWeightsList; //! Weights for Q Vectors
 TList *fCRCZDCCalibList; //! ZDC calibration
 TList *fZDCESEList; //! ZDC ESE
 TProfile *fCRCQ2Re[fCRCMaxnRun]; //! Q Vectors Re
 TProfile *fCRCQ2Im[fCRCMaxnRun]; //! Q Vectors Im
 TProfile *fCRCQ2ReCorr[fCRCMaxnRun]; //! Q Vectors Re
 TProfile *fCRCQ2ImCorr[fCRCMaxnRun]; //! Q Vectors Im
 TH2D *fCRCPhiHist[fCRCMaxnRun][fCRCMaxnCen][2]; //! Phi Hist for weights
 AliFlowVector fTPCQ2Recenter;
 
 // temp
 TH2D *fRunPhiEtaHist[fCRCMaxnCen][2]; //! Run-by-run PhiEtaHist
 TProfile *fTPCQHist[2];  //! Run-by-run TPCQvecHist
 TProfile *fZDCQHist[4];               //! Run-by-run ZDCQvecHist
 TH1D *fZDCESEMinHist[2]; //!
 TH1D *fZDCESEMaxHist[2]; //!
 
 TH1D *fCRCVZEvPlA[fCRCMaxnRun][fCRCMaxnCen][fCRCnHar]; //! Ev Plane VZEROA
 TH1D *fCRCVZEvPlC[fCRCMaxnRun][fCRCMaxnCen][fCRCnHar]; //! Ev Plane VZEROC
 TProfile *fCRCVZQVecA[fCRCMaxnRun][2]; //! Q Vectors VZERO-A
 TProfile *fCRCVZQVecC[fCRCMaxnRun][2]; //! Q Vectors VZERO-C
 TProfile *fCRCVZQVecCov[fCRCMaxnRun][4]; //! VZs Q Vectors correlations
 
 TH1D *fCRCZDCEvPlA[fCRCMaxnRun][fCRCMaxnCen]; //! Ev Plane ZDCN-A
 TH1D *fCRCZDCEvPlC[fCRCMaxnRun][fCRCMaxnCen]; //! Ev Plane ZDCN-C
 TProfile *fCRCZDCQVecA[fCRCMaxnRun][2]; //! Q Vectors ZDCN-A
 TProfile *fCRCZDCQVecC[fCRCMaxnRun][2]; //! Q Vectors ZDCN-C
 TProfile *fCRCZDCQVecACorr[fCRCMaxnRun][2]; //! Q Vectors ZDCN-A
 TProfile *fCRCZDCQVecCCorr[fCRCMaxnRun][2]; //! Q Vectors ZDCN-C
 TProfile *fCRCZDCQVecCov[fCRCMaxnRun][4]; //! ZDCs Q Vectors correlations
 
 TProfile *fCRCVZvsZDCCov[fCRCMaxnRun][16]; //! ZDC vs VZ Q Vectors correlations
 
 // CRCVZERO
 TList *fCRCVZList; //! VZERO CRC List
 const static Int_t fCRCVZnCR = 13;
 const static Int_t fCRCVZnEtaBin = 5;
 TList *fCRCVZRbRList; //! CRC list of histograms RbR
 AliFlowVector fVZFlowVect[2][fCRCnHar];
 TList *fCRCVZRunsList[fCRCMaxnRun];  //! list of runs
 TProfile *fCRCVZCorrPro[fCRCMaxnRun][2][fCRCVZnEtaBin][fCRCMaxnCen]; //! correlation profile, [CRCBin][eg]
 TProfile *fCRCVZCorrProd2p2pHist[fCRCMaxnRun][2][fCRCVZnEtaBin][fCRCMaxnCen]; //! correlation weights histo, [CRCBin][eg]
 TProfile *fCRCVZNUAPro[fCRCMaxnRun][4][fCRCVZnEtaBin][fCRCMaxnCen]; //! correlation profile, [CRCBin][eg]
 TH1D *fCRCVZCorrProdTempHist[2][fCRCVZnEtaBin][fCRCMaxnCen]; //! temporary correlation products for covariances, [CRCBin][eg]
 TH1D *fCRCVZCorrHist[2][fCRCVZnEtaBin][fCRCMaxnCen]; //! <<2'>>, [CRCBin][eg]
 TH1D *fCRCVZCFunHist[fCRCVZnEtaBin][fCRCMaxnCen]; //! correlation function histo, [CRCBin][eg]
 TH2D *fCRCVZCovHist[2][fCRCVZnEtaBin][fCRCMaxnCen]; //! covariances final histo
 
 // CRCZDC
 TList *fCRCZDCList; //! ZDCERO CRC List
 const static Int_t fCRCZDCnCR = 13;
 const static Int_t fCRCZDCnEtaBin = 5;
 TList *fCRCZDCRbRList; //! CRC list of histograms RbR
 AliFlowVector fZDCFlowVect[2];
 TH1D *fCRCZDCQRe[4][fCRCnHar]; //! real part [0=pos,1=neg][0=back,1=forw][m]
 TH1D *fCRCZDCQIm[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][m]
 TH1D *fCRCZDCMult[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][k]
 TList *fCRCZDCRunsList[fCRCMaxnRun]; //! list of runs
 TProfile *fCRCZDCCorrPro[fCRCMaxnRun][2][fCRCZDCnEtaBin][fCRCMaxnCen]; //! correlation profile, [CRCBin][eg]
 TProfile *fCRCZDCCorrProd2p2pHist[fCRCMaxnRun][2][fCRCZDCnEtaBin][fCRCMaxnCen]; //! correlation weights histo, [CRCBin][eg]
 TProfile *fCRCZDCNUAPro[fCRCMaxnRun][4][fCRCZDCnEtaBin][fCRCMaxnCen]; //! correlation profile, [CRCBin][eg]
 TH1D *fCRCZDCCorrProdTempHist[2][fCRCZDCnEtaBin][fCRCMaxnCen]; //! temporary correlation products for covariances, [CRCBin][eg]
 TH1D *fCRCZDCCorrHist[2][fCRCZDCnEtaBin][fCRCMaxnCen]; //! <<2'>>, [CRCBin][eg]
 TH1D *fCRCZDCCFunHist[fCRCZDCnEtaBin][fCRCMaxnCen]; //! correlation function histo, [CRCBin][eg]
 TH2D *fCRCZDCCovHist[2][fCRCZDCnEtaBin][fCRCMaxnCen]; //! covariances final histo
 TProfile *fCRCZDCSpectra[fCRCZDCnEtaBin][fCRCMaxnCen]; //! spectra
 
 // CRC Pt differential
 TList *fCRCPtList; //! list to hold CRC histograms
 
 Int_t fCRCPtnPtBin;
 const static Int_t fCRCPtnPtBinMax = 48;
 Double_t fCRCPtMinPt;
 Double_t fCRCPtMaxPt;
 Double_t fCRCPtwPtBin;
 Double_t *fCRCPtvarPtBins; //!
 Int_t fCRCPtnCenBin;
 Double_t fCRCPtCenMin;
 Double_t fCRCPtCenMax;
 Double_t fCRCPtwCenBin;
 
 TNtuple *fCRCPtTPCTNt; //! TNtuple test
 TNtuple *fCRCPtVZTNt; //! TNtuple test
 TNtuple *fCRCPtZDCTNt; //! TNtuple test
 THnSparse *fCRCPtEbEQVec; //!
 
 // CME
 const static Int_t fCMETPCnCR = 20;
 const static Int_t fCMEnEtaBin = 2;
 const static Int_t fCMETPCnDist = 4;
 const static Int_t fCMEZDCnCR = 20;
 const static Int_t fCMEZDCnDist = 4;
 TList *fCMEList;    //! CME List
 TList *fCMERbRList; //! CME list of histograms RbR
 TList *fCMETPCList; //! CME list of histograms TPC only
 TList *fCMEZDCList; //! CME list of histograms TPC-ZDCs
 TList *fCMERunsList[fCRCMaxnRun]; //! list of runs
 TH1D *fCMEQRe[4][fCRCnHar]; //! real part [0=pos,1=neg][0=back,1=forw][m]
 TH1D *fCMEQIm[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][m]
 TH1D *fCMEMult[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][k]
 TH1D *fCMETPCCorHist[fCMEnEtaBin][fCRCMaxnCen]; //! <<2'>>, [CRCBin][eg]
 TH2D *fCMETPCCovHist[fCMEnEtaBin][fCRCMaxnCen]; //! correlation function histo, [CRCBin][eg]
 TH1D *fCMETPCDistHist[fCMEnEtaBin][fCRCMaxnCen][fCMETPCnDist]; //! <<2'>>, [CRCBin][eg]
 TH1D *fCMEZDCCorHist[fCMEnEtaBin][fCRCMaxnCen]; //! <<2'>>, [CRCBin][eg]
 TH2D *fCMEZDCCovHist[fCMEnEtaBin][fCRCMaxnCen]; //! correlation function histo, [CRCBin][eg]
 TH1D *fCMEZDCDistHist[fCMEnEtaBin][fCRCMaxnCen][fCMEZDCnDist]; //! <<2'>>, [CRCBin][eg]
 TProfile *fCMETPCCorPro[fCRCMaxnRun][fCMEnEtaBin][fCRCMaxnCen]; //! correlation profile, [CRCBin][eg]
 TProfile *fCMETPCCovPro[fCRCMaxnRun][fCMEnEtaBin][fCRCMaxnCen]; //! correlation weights histo, [CRCBin][eg]
 TProfile *fCMETPCNUAPro[fCRCMaxnRun][fCMEnEtaBin][fCRCMaxnCen]; //! correlation profile, [CRCBin][eg]
 TProfile *fCMEZDCCorPro[fCRCMaxnRun][fCMEnEtaBin][fCRCMaxnCen]; //! correlation profile, [CRCBin][eg]
 TProfile *fCMEZDCCovPro[fCRCMaxnRun][fCMEnEtaBin][fCRCMaxnCen]; //! correlation weights histo, [CRCBin][eg]
 TProfile *fCMEZDCNUAPro[fCRCMaxnRun][fCMEnEtaBin][fCRCMaxnCen]; //! correlation profile, [CRCBin][eg]
 
 // CRC2
  const static Int_t fkNCorCRC2 = 6;
 TList *fCRC2List; //! ZDCERO CRC List
 Int_t fCRC2nEtaBins; // CRC2 n eta bins
 TList *fCRC2RbRList; //! CRC list of histograms RbR
 TH1D *fCRC2QRe[3][fCRCnHar]; //! real part [0=pos,1=neg][0=back,1=forw][m]
 TH1D *fCRC2QIm[3][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][m]
 TH1D *fCRC2Mul[3][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][k]
 TList *fCRC2RunsList[fCRCMaxnRun]; //! list of runs
 TProfile *fCRC2CorPro[fCRCMaxnRun][fCRCMaxnCen][fkNCorCRC2]; //! correlation profile, [CRCBin][eg]
 TProfile *fCRC2NUAPro[fCRCMaxnRun][fCRCMaxnCen][fkNCorCRC2][4]; //! NUA terms profile, [CRCBin][eg]
 TH1D *fCRC2CorHist[fCRCMaxnCen][fkNCorCRC2][2]; //! <<2'>>, [centrality][pos,neg,all][corr,spectra]
 TH1D *fCRC2NUAHist[fCRCMaxnCen][fkNCorCRC2][4]; //! NUA hist
 TProfile *fCRC2CovPro[fCRCMaxnCen][fkNCorCRC2][fkNCorCRC2]; //! Cov pro
 TH1D *fCRC2CovHist[fCRCMaxnCen][fkNCorCRC2][fkNCorCRC2]; //! Cov hist
 
 // Flow all
 const static Int_t fFlowNHarm = 1;
 const static Int_t fFlowNHarmMax = 9;
 Int_t fPtDiffNBins;
 Int_t fEtaDiffNBins;
 Double_t fPtDiffMinPt;
 Double_t fPtDiffMaxPt;
 TH1D *fPtDiffQRe[2][fFlowNHarmMax]; //! real part [0=pos,1=neg][0=back,1=forw][m]
 TH1D *fPtDiffQIm[2][fFlowNHarmMax]; //! imaginary part [0=pos,1=neg][0=back,1=forw][m]
 TH1D *fPtDiffMul[2][fFlowNHarmMax]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][k]
  TH1D *fEtaDiffQRe[fFlowNHarmMax]; //! real part [0=pos,1=neg][0=back,1=forw][eta]
  TH1D *fEtaDiffQIm[fFlowNHarmMax]; //! imaginary part [0=pos,1=neg][0=back,1=forw][eta]
  TH1D *fEtaDiffMul[fFlowNHarmMax]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][eta]
 
 // Flow SP ZDC
 TList *fFlowSPZDCList;    //! SPZDC List
 TList *fFlowSPZDCRbRList; //! CRC list of histograms RbR
 TList *fFlowSPZDCRunsList[fCRCMaxnRun]; //! list of runs
 const static Int_t fFlowNPro = 14;
 const static Int_t fFlowNNUA = 4;
 TProfile *fFlowSPZDCCorPro[fCRCMaxnRun][fCRCMaxnCen][fFlowNHarm][fFlowNPro]; //! correlation profile, [CRCBin][eg]
 TProfile *fFlowSPZDCCorNUA[fCRCMaxnRun][fCRCMaxnCen][fFlowNPro]; //! NUA profile, [CRCBin][eg]
 TH1D *fFlowSPZDCCorHist[fCRCMaxnCen][fFlowNHarm][fFlowNPro]; //! <<2'>>, [CRCBin][eg]
 TProfile *fFlowSPZDCIntPro[fCRCMaxnRun][fFlowNHarm][fFlowNPro]; //! reference flow
 TH1D *fFlowSPZDCIntHist[fFlowNHarm][fFlowNPro]; //!
 TH1D *fFlowSPZDCFinalPtDifHist[fCRCMaxnCen][fFlowNHarm][fFlowNPro]; //!
 TH2F *fFlowSPZDCSpectra; //!
 TProfile *fFlowSPZDCIntNUA[fCRCMaxnRun][fFlowNNUA]; //!
 
 // Flow QC
 TList *fFlowQCList;    //! QC List
 TProfile *fFlowQCCorPro[fCRCMaxnCen][fFlowNHarm][5]; //! correlation profile, [CRCBin][eg]
 TProfile *fFlowQCNUAPro[fCRCMaxnCen][fFlowNHarm][6]; //! NUA profile, [CRCBin][eg]
 TH1D *fFlowQCCorHist[fCRCMaxnCen][fFlowNHarm][5]; //! <<2'>>, [CRCBin][eg]
 TProfile *fFlowQCIntCorPro[fFlowNHarm][3]; //!
 TH1D *fFlowQCIntCorHist[fFlowNHarm][3]; //!
 
 TH1D *fFlowQCFinalPtDifHist[fCRCMaxnCen][fFlowNHarm][5]; //!
 TH1D *fFlowQCFinalPtIntHist[fCRCMaxnCen][5]; //!
 TH1D *fFlowQCSpectra[fCRCMaxnCen]; //!
 
 // Flow SP VZ
 TList *fFlowSPVZList;    //! SPVZ List
 TProfile *fFlowSPVZCorPro[fCRCMaxnCen][fFlowNHarm][5]; //! correlation profile, [CRCBin][eg]
 TProfile *fFlowSPVZNUAPro[fCRCMaxnCen][fFlowNHarm][5][4]; //! NUA profile, [CRCBin][eg]
 TH1D *fFlowSPVZCorHist[fCRCMaxnCen][fFlowNHarm][5]; //! <<2'>>, [CRCBin][eg]
 
 // Various:
 TList *fVariousList; //! list to hold various unclassified objects
 TH1D *fMultHist; //! Multiplicity distribution
 TH1D *fCenHist; //! Centrality distribution
 TH1D* fCenWeightsHist; //! Centrality weights
 TH1D* fCenWeigCalHist; //! Centrality weights
 TH1D* fPtWeightsHist[10]; //! Pt weights
 TH1D* fEtaWeightsHist[10][21][2]; //! Eta weights
 TH1D* fNvsCenCut[2][2]; //! ZDC mult cuts
 Double_t *fCRCPtBins; //!
 Double_t *fCorrMap; //!
 Double_t fCenWeightEbE;
 TH2F* fhZNCentroid[fCRCMaxnCen][2]; //! Centroid position x-y
 TH1F* fhZNSpectra[fCRCMaxnCen][2]; //! ZN spectra
 TH2F* fhZNCvsZNA[fCRCMaxnCen]; //! ZNA-ZNC correlation
 TH2F* fhZNCenvsMul[fCRCMaxnCen][2]; //! rad vs mul
 TH2F* fhZNResvsMul[fCRCMaxnCen][2]; //! res vs mul
 TH2F* fhZNResvsCen[fCRCMaxnCen][2]; //! res vs rad
 TH2F* fhZNvsCen[2]; //! cen vs mul
  TH2F* fhZNvsTCen[2]; //! cen vs mul
  TH2F* fhCenvsMul[2]; //! cen vs mul
  TH2F* fhCenvsDif[2]; //! cen vs mul
  TH2F* fhZNvsMul[2]; //! cen vs mul
  
 TProfile* fhZNQVecCov[4]; //! Q-vec cov.
 TH2D *fZDCESEHistEP[fCRCMaxnCen]; //! Test ZDC ESE
 TH2D *fZDCESEHistQV[fCRCMaxnCen]; //! Test ZDC ESE
 Bool_t fQAZDCCuts;
 Bool_t fQAZDCCutsFlag;
 Int_t fMinMulZN;
 Float_t fMaxDevZN;
 
 ClassDef(AliFlowAnalysisCRC, 6);
 
};

//================================================================================================================

#endif





