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

/********************************
 * analysis task for CRC        *
 *                              *
 * author: Jacopo Margutti      *
 *         (margutti@nikhef.nl) *
 * ******************************/

#ifndef AliAnalysisTaskCRC_H
#define AliAnalysisTaskCRC_H

#include "AliAnalysisTaskSE.h"
#include "AliFlowCommonConstants.h"

class TString;
class TList;
class AliFlowEventSimple;
class AliFlowEvent;
class AliFlowAnalysisCRC;

//==============================================================================================================

class AliAnalysisTaskCRC : public AliAnalysisTaskSE{
public:
  AliAnalysisTaskCRC();
  AliAnalysisTaskCRC(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskCRC(){};
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  // Common:
  void SetBookOnlyBasicCCH(Bool_t const bobcch) {this->fBookOnlyBasicCCH = bobcch;};
  Bool_t GetBookOnlyBasicCCH() const {return this->fBookOnlyBasicCCH;};
  void SetFillMultipleControlHistograms(Bool_t const fmch) {this->fFillMultipleControlHistograms = fmch;};
  Bool_t GetFillMultipleControlHistograms() const {return this->fFillMultipleControlHistograms;};
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  Int_t GetHarmonic() const {return this->fHarmonic;};
  void SetApplyCorrectionForNUA(Bool_t const applyCorrectionForNUA) {this->fApplyCorrectionForNUA = applyCorrectionForNUA;};
  Bool_t GetApplyCorrectionForNUA() const {return this->fApplyCorrectionForNUA;};
  void SetApplyCorrectionForNUAVsM(Bool_t const applyCorrectionForNUAVsM) {this->fApplyCorrectionForNUAVsM = applyCorrectionForNUAVsM;};
  Bool_t GetApplyCorrectionForNUAVsM() const {return this->fApplyCorrectionForNUAVsM;};
  void SetPropagateErrorAlsoFromNIT(Bool_t const peafNIT) {this->fPropagateErrorAlsoFromNIT = peafNIT;};
  Bool_t GetPropagateErrorAlsoFromNIT() const {return this->fPropagateErrorAlsoFromNIT;};
  void SetCalculateDiffFlow(Bool_t const calculateDiffFlow) {this->fCalculateDiffFlow = calculateDiffFlow;};
  Bool_t GetCalculateDiffFlow() const {return this->fCalculateDiffFlow;};
  void SetCalculate2DDiffFlow(Bool_t const calculate2DDiffFlow) {this->fCalculate2DDiffFlow = calculate2DDiffFlow;};
  Bool_t GetCalculate2DDiffFlow() const {return this->fCalculate2DDiffFlow;};
  void SetCalculateDiffFlowVsEta(Bool_t const cdfve) {this->fCalculateDiffFlowVsEta = cdfve;};
  Bool_t GetCalculateDiffFlowVsEta() const {return this->fCalculateDiffFlowVsEta;};
  void SetStoreDistributions(Bool_t const storeDistributions) {this->fStoreDistributions = storeDistributions;};
  Bool_t GetStoreDistributions() const {return this->fStoreDistributions;};
  void SetCalculateCumulantsVsM(Bool_t const ccvm) {this->fCalculateCumulantsVsM = ccvm;};
  Bool_t GetCalculateCumulantsVsM() const {return this->fCalculateCumulantsVsM;};
  void SetCalculateAllCorrelationsVsM(Bool_t const cacvm) {this->fCalculateAllCorrelationsVsM = cacvm;};
  Bool_t GetCalculateAllCorrelationsVsM() const {return this->fCalculateAllCorrelationsVsM;};
  void SetCalculateMixedHarmonics(Bool_t const cmh) {this->fCalculateMixedHarmonics = cmh;};
  Bool_t GetCalculateMixedHarmonics() const {return this->fCalculateMixedHarmonics;};
  void SetCalculateMixedHarmonicsVsM(Bool_t const cmhvm) {this->fCalculateMixedHarmonicsVsM = cmhvm;};
  Bool_t GetCalculateMixedHarmonicsVsM() const {return this->fCalculateMixedHarmonicsVsM;};
  void SetStoreControlHistograms(Bool_t const sch) {this->fStoreControlHistograms = sch;};
  Bool_t GetStoreControlHistograms() const {return this->fStoreControlHistograms;};
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
  void SetWeightsList(TList* const kList) {this->fWeightsList = (TList*)kList->Clone();};
  TList* GetWeightsList() const {return this->fWeightsList;};
  
  // Multiparticle correlations vs multiplicity:
  void SetnBinsMult(Int_t const nbm) {this->fnBinsMult = nbm;};
  Int_t GetnBinsMult() const {return this->fnBinsMult;};
  void SetMinMult(Double_t const minm) {this->fMinMult = minm;};
  Double_t GetMinMult() const {return this->fMinMult;};
  void SetMaxMult(Double_t const maxm) {this->fMaxMult = maxm;};
  Double_t GetMaxMult() const {return this->fMaxMult;};
  
  // Particle weights:
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
  void SetUsePhiEtaWeightsChDep(Bool_t const uPhiEtaW) {this->fUsePhiEtaWeightsChDep = uPhiEtaW;};
  Bool_t GetUsePhiEtaWeightsChDep() const {return this->fUsePhiEtaWeightsChDep;};
  void SetUsePhiEtaWeightsVtxDep(Bool_t const uPhiEtaW) {this->fUsePhiEtaWeightsVtxDep = uPhiEtaW;};
  Bool_t GetUsePhiEtaWeightsVtxDep() const {return this->fUsePhiEtaWeightsVtxDep;};
  void SetUsePhiEtaCuts(Bool_t const uPhiEtaW) {this->fUsePhiEtaCuts = uPhiEtaW;};
  Bool_t GetUsePhiEtaCuts() const {return this->fUsePhiEtaCuts;};
  void SetUseZDCESEMulWeights(Bool_t const uPhiEtaW) {this->fUseZDCESEMulWeights = uPhiEtaW;};
  Bool_t GetUseZDCESEMulWeights() const {return this->fUseZDCESEMulWeights;};
  void SetUseZDCESESpecWeights(Bool_t const uPhiEtaW) {this->fUseZDCESESpecWeights = uPhiEtaW;};
  Bool_t GetUseZDCESESpecWeights() const {return this->fUseZDCESESpecWeights;};
  
  // Event weights:
  void SetMultiplicityWeight(const char *multiplicityWeight) {*this->fMultiplicityWeight = multiplicityWeight;};
  void SetMultiplicityIs(AliFlowCommonConstants::ERefMultSource mi) {this->fMultiplicityIs = mi;};
  // # of bins for correlation axis in fDistributions[4], fCorrelation2468VsMult[4] and fCorrelationProduct2468VsMult[1]
  void SetnBinsForCorrelations(Int_t const nb) {this->fnBinsForCorrelations = nb;};
  Int_t GetnBinsForCorrelations() const {return this->fnBinsForCorrelations;};
  
  // Boundaries for distributions of correlations:
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
  
  // bootstrap:
  void SetUseBootstrap(Bool_t const ub) {this->fUseBootstrap = ub;};
  Bool_t GetUseBootstrap() const {return this->fUseBootstrap;};
  void SetUseBootstrapVsM(Bool_t const ubVsM) {this->fUseBootstrapVsM = ubVsM;};
  Bool_t GetUseBootstrapVsM() const {return this->fUseBootstrapVsM;};
  void SetnSubsamples(Int_t const ns) {this->fnSubsamples = ns;};
  Int_t GetnSubsamples() const {return this->fnSubsamples;};
  
  // Charge-Rapidity Correlations:
  void SetCalculateCRC(Bool_t const cCRC) {this->fCalculateCRC = cCRC;};
  Bool_t GetCalculateCRC() const {return this->fCalculateCRC;};
  void SetCalculateCRCPt(Bool_t const cCRC) {this->fCalculateCRCPt = cCRC;};
  Bool_t GetCalculateCRCPt() const {return this->fCalculateCRCPt;};
  void SetCalculateCME(Bool_t const cCRC) {this->fCalculateCME = cCRC;};
  Bool_t GetCalculateCME() const {return this->fCalculateCME;};
  void SetCalculateCRCInt(Bool_t const cCRC) {this->fCalculateCRCInt = cCRC;};
  Bool_t GetCalculateCRCInt() const {return this->fCalculateCRCInt;};
  void SetCalculateCRC2(Bool_t const cCRC) {this->fCalculateCRC2 = cCRC;};
  Bool_t GetCalculateCRC2() const {return this->fCalculateCRC2;};
  void SetCalculateCRCVZ(Bool_t const cCRC) {this->fCalculateCRCVZ = cCRC;};
  Bool_t GetCalculateCRCVZ() const {return this->fCalculateCRCVZ;};
  void SetCalculateCRCZDC(Bool_t const cCRC) {this->fCalculateCRCZDC = cCRC;};
  Bool_t GetCalculateCRCZDC() const {return this->fCalculateCRCZDC;};
  void SetCalculateEbEFlow(Bool_t const cCRC) {this->fCalculateEbEFlow = cCRC;};
  Bool_t GetCalculateEbEFlow() const {return this->fCalculateEbEFlow;};
  void SetStoreZDCQVecVtxPos(Bool_t const cCRC) {this->fStoreZDCQVecVtxPos = cCRC;};
  Bool_t GetStoreZDCQVecVtxPos() const {return this->fStoreZDCQVecVtxPos;};
  void SetCRC2nEtaBins(Int_t NB) {this->fCRC2nEtaBins = NB;};
  Int_t GetCRC2nEtaBins() {return this->fCRC2nEtaBins;};
  void SetCalculateFlowQC(Bool_t const cCRC) {this->fCalculateFlowQC = cCRC;};
  Bool_t GetCalculateFlowQC() const {return this->fCalculateFlowQC;};
  void SetCalculateFlowZDC(Bool_t const cCRC) {this->fCalculateFlowZDC = cCRC;};
  Bool_t GetCalculateFlowZDC() const {return this->fCalculateFlowZDC;};
  void SetCalculateFlowVZ(Bool_t const cCRC) {this->fCalculateFlowVZ = cCRC;};
  Bool_t GetCalculateFlowVZ() const {return this->fCalculateFlowVZ;};
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
  void SetTestSin(Bool_t const cCRC) {this->fCRCTestSin = cCRC;};
  Bool_t GetTestSin() const {return this->fCRCTestSin;};
  void SetNUAforCRC(Bool_t const cCRC) {this->fUseNUAforCRC = cCRC;};
  Bool_t GetNUAforCRC() const {return this->fUseNUAforCRC;};
  void SetUseCRCRecenter(Bool_t const cCRC) {this->fUseCRCRecenter = cCRC;};
  Bool_t GetUseCRCRecenter() const {return this->fUseCRCRecenter;};
  void SetCRCEtaRange(Double_t const etamin, Double_t const etamax) {this->fCRCEtaMin = etamin; this->fCRCEtaMax = etamax;};
  void SetQVecList(TList* const kList) {this->fQVecList = (TList*)kList->Clone();};
  TList* GetQVecList() const {return this->fQVecList;};
  void SetZDCESEList(TList* const kList) {this->fZDCESEList = (TList*)kList->Clone();};
  TList* GetZDCESEList() const {return this->fZDCESEList;};
  void SetCRCZDCCalibList(TList* const wlist) {this->fCRCZDCCalibList = (TList*)wlist->Clone();}
  TList* GetCRCZDCCalibList() const {return this->fCRCZDCCalibList;}
  void SetCRCVZEROCalibList(TList* const wlist) {this->fCRCVZEROCalibList = (TList*)wlist->Clone();}
  TList* GetCRCVZEROCalibList() const {return this->fCRCVZEROCalibList;}
  void SetCRCZDCResList(TList* const wlist) {this->fCRCZDCResList = (TList*)wlist->Clone();}
  TList* GetCRCZDCResList() const {return this->fCRCZDCResList;}
  void SetnCenBin(Int_t const n) {this->fnCenBin = n;};
  Int_t GetnCenBin() const {return this->fnCenBin;};
  void SetFlowQCCenBin(Int_t const TL) {this->fFlowQCCenBin = TL;};
  Int_t GetFlowQCCenBin() const {return this->fFlowQCCenBin;};
  void SetFlowQCDeltaEta(Double_t const TL) {this->fFlowQCDeltaEta = TL;};
  Double_t GetFlowQCDeltaEta() const {return this->fFlowQCDeltaEta;};
  void SetCenBinWidth(Double_t const n) {this->fCenBinWidth = n;};
  Double_t GetCenBinWidth() const {return this->fCenBinWidth;};
  void SetDataSet(TString const n) {this->fDataSet = n;};
  TString GetDataSet() const {return this->fDataSet;};
  void SetInteractionRate(TString const n) {this->fInteractionRate = n;};
  TString GetInteractionRate() const {return this->fInteractionRate;}
  void SetSelectCharge(TString const n) {this->fSelectCharge = n;};
  TString GetSelectCharge() const {return this->fSelectCharge;}
  void SetPOIExtraWeights(TString const n) {this->fPOIExtraWeights = n;};
  TString GetPOIExtraWeights() const {return this->fPOIExtraWeights;}
  void SetCorrWeight(TString const n) {this->fCorrWeight = n;};
  TString GetCorrWeight() const {return this->fCorrWeight;};
  void SetCenWeightsHist(TH1D* const n) {this->fCenWeightsHist = n;};
  TH1D* GetCenWeightsHist() const {return this->fCenWeightsHist;};
  void SetPtWeightsHist(TH1D* const n, Int_t c) {this->fPtWeightsHist[c] = n;};
  TH1D* GetPtWeightsHist(Int_t c) const {return this->fPtWeightsHist[c];};
  void SetEtaWeightsHist(TH1D* const n, Int_t h, Int_t b, Int_t c) {this->fEtaWeightsHist[h][b][c] = n;};
  TH1D* GetEtaWeightsHist(Int_t h, Int_t b, Int_t c) const {return this->fEtaWeightsHist[h][b][c];};
  void SetZDCESEMultWeightsHist(TH2F* const n, Int_t h) {this->fZDCESEMultWeightsHist[h] = n;};
  TH2F* GetZDCESEMultWeightsHist(Int_t h) const {return this->fZDCESEMultWeightsHist[h];};
  void SetZDCESESpecWeightsHist(TH2F* const n, Int_t h) {this->fZDCESESpecWeightsHist[h] = n;};
  TH2F* GetZDCESESpecWeightsHist(Int_t h) const {return this->fZDCESESpecWeightsHist[h];};
  void SetNvsCenCut(TH1D* const n, Int_t c, Int_t h) {this->fNvsCenCut[c][h] = n;};
  TH1D* GetNvsCenCut(Int_t c, Int_t h) const {return this->fNvsCenCut[c][h];};
  void SetQAZDCCuts(Bool_t const cCRC) {this->fQAZDCCuts = cCRC;};
  Bool_t GetQAZDCCuts() const {return this->fQAZDCCuts;};
  void SetMinMulZN(Int_t weights) {this->fMinMulZN = weights;};
  Int_t GetMinMulZN() const {return this->fMinMulZN;};
  void SetMaxDevZN(Float_t weights) {this->fMaxDevZN = weights;};
  Float_t GetMaxDevZN() const {return this->fMaxDevZN;};
  void SetZDCGainAlpha( Float_t a ) { fZDCGainAlpha = a; }
  
private:
  AliAnalysisTaskCRC(const AliAnalysisTaskCRC& aatqc);
  AliAnalysisTaskCRC& operator=(const AliAnalysisTaskCRC& aatqc);
  
  AliFlowEvent *fEvent;         // the input event
  AliFlowAnalysisCRC *fQC;            // CRC object
  TList *fListHistos;                 // collection of output
  // Common:
  Bool_t fBookOnlyBasicCCH;              // book only basis common control histrograms (by default book them all)
  Bool_t fFillMultipleControlHistograms; // fill separately control histos for events with >= 2, 4, 6 and 8 particles
  Int_t fHarmonic;                       // harmonic
  Bool_t fApplyCorrectionForNUA;         // apply correction for non-uniform acceptance
  Bool_t fApplyCorrectionForNUAVsM;      // apply correction for non-uniform acceptance versus M
  Bool_t fPropagateErrorAlsoFromNIT;     // propagate error by taking into account also non-isotrpic terms
  Bool_t fCalculateDiffFlow;             // calculate differential flow in pt or eta
  Bool_t fCalculate2DDiffFlow;           // calculate differential flow in (pt,eta) (Remark: this is very expensive in terms of CPU time)
  Bool_t fCalculateDiffFlowVsEta;        // if you set kFALSE only differential flow vs pt is calculated
  Bool_t fStoreDistributions;            // store or not distributions of correlations
  Bool_t fCalculateCumulantsVsM;         // calculate cumulants versus multiplicity
  Bool_t fCalculateAllCorrelationsVsM;   // calculate all correlations versus multiplicity
  Bool_t fCalculateMixedHarmonics;       // calculate all mixed harmonics correlations
  Bool_t fCalculateMixedHarmonicsVsM;    // calculate all mixed harmonics correlations versus multiplicity
  Bool_t fStoreControlHistograms;        // store or not control histograms
  Bool_t fMinimumBiasReferenceFlow;      // store as reference flow in AliFlowCommonHistResults the minimum bias result (kFALSE by default)
  Bool_t fForgetAboutCovariances;        // when propagating error forget about the covariances
  Bool_t fStoreVarious; // store phi distribution for one event to illustrate flow
  Int_t fExactNoRPs;                     // when shuffled, select only this number of RPs for the analysis
  Bool_t fUse2DHistograms;               // use TH2D instead of TProfile to improve numerical stability in reference flow calculation
  Bool_t fFillProfilesVsMUsingWeights;   // if the width of multiplicity bin is 1, weights are not needed
  Bool_t fUseQvectorTerms; // use TH2D with separate Q-vector terms instead of TProfile to improve numerical stability in reference flow calculation
  // Multiparticle correlations vs multiplicity:
  Int_t fnBinsMult;                   // number of multiplicity bins for flow analysis versus multiplicity
  Double_t fMinMult;                  // minimal multiplicity for flow analysis versus multiplicity
  Double_t fMaxMult;                  // maximal multiplicity for flow analysis versus multiplicity
  // Particle weights:
  Bool_t fUseParticleWeights;         // use any particle weights
  Bool_t fUsePhiWeights;              // use phi weights
  Bool_t fUsePhiEtaCuts;              // use phi,eta cuts (for NUA)
  Bool_t fUsePtWeights;               // use pt weights
  Bool_t fUseEtaWeights;              // use eta weights
  Bool_t fUseTrackWeights;            // use track weights (e.g. VZERO sector weights)
  Bool_t fUsePhiEtaWeights;           // use phi,eta weights
  Bool_t fUsePhiEtaWeightsChDep;      // use phi,eta weights ch dep
  Bool_t fUsePhiEtaWeightsVtxDep;      // use phi,eta weights vtx dep
  Bool_t fUseZDCESEMulWeights;        // use ZDC-ESE mult. weights
  Bool_t fUseZDCESESpecWeights;       // use ZDC-ESE mult. weights
  TList *fWeightsList;                // list with weights
  // Event weights:
  TString *fMultiplicityWeight;       // event-by-event weights for multiparticle correlations ("combinations","unit" or "multiplicity")
  AliFlowCommonConstants::ERefMultSource fMultiplicityIs;           // by default "#RPs", other supported options are "RefMultFromESD" = ref. mult. from ESD, and "#POIs"
  Int_t fnBinsForCorrelations; // # of bins for correlation axis in fDistributions[4], fCorrelation2468VsMult[4] and fCorrelationProduct2468VsMult[1]
  // Boundaries for distributions of correlations:
  Double_t fMinValueOfCorrelation[4]; // min values of <2>, <4>, <6> and <8>
  Double_t fMaxValueOfCorrelation[4]; // max values of <2>, <4>, <6> and <8>
  Double_t fMinValueOfCorrelationProduct[1]; // min values of <2><4>, <2><6>, <2><8>, <4><6> etc. TBI add the other ones when needed first time
  Double_t fMaxValueOfCorrelationProduct[1]; // max values of <2><4>, <2><6>, <2><8>, <4><6> etc. TBI add the other ones when needed first time
  Double_t fMinValueOfQvectorTerms[4]; // min value of Q-vector terms
  Double_t fMaxValueOfQvectorTerms[4]; // max value of Q-vector terms
  // Bootstrap:
  Bool_t fUseBootstrap; // use bootstrap to estimate statistical spread
  Bool_t fUseBootstrapVsM; // use bootstrap to estimate statistical spread for results vs M
  Int_t fnSubsamples; // number of subsamples (SS), by default 10
  // Charge-Eta Asymmetry
  Bool_t fCalculateCRC; // calculate CRC quantities
  Bool_t fCalculateCRCPt;
  Bool_t fCalculateCME;
  Bool_t fCalculateCRCInt;
  Bool_t fCalculateCRC2;
  Bool_t fCalculateCRCVZ;
  Bool_t fCalculateCRCZDC;
  Bool_t fCalculateEbEFlow;
  Bool_t fStoreZDCQVecVtxPos;
  Int_t fCRC2nEtaBins; // CRC2 n eta bins
  Bool_t fCalculateFlowQC;
  Bool_t fCalculateFlowZDC;
  Bool_t fCalculateFlowVZ;
  Bool_t fUseVZERO;
  Bool_t fUseZDC;
  Bool_t fRecenterZDC;
  Bool_t fDivSigma;
  Bool_t fInvertZDC;
  Bool_t fCRCTestSin;
  Bool_t fUseNUAforCRC;
  Bool_t fUseCRCRecenter;
  Double_t fCRCEtaMin;
  Double_t fCRCEtaMax;
  Int_t fnCenBin;
  Int_t fFlowQCCenBin;
  Double_t fFlowQCDeltaEta;
  Double_t fCenBinWidth;
  TString fDataSet;
  TString fInteractionRate;
  TString fSelectCharge;
  TString fPOIExtraWeights;
  TString fCorrWeight;
  TList *fQVecList;       // list with weights
  TList *fCRCZDCCalibList; // ZDC calibration
  TList *fCRCVZEROCalibList; // ZDC calibration
  TList *fCRCZDCResList; // ZDC rescaling
  TList *fZDCESEList;       // list with weights
  TH1D* fCenWeightsHist;
  TH1D* fPtWeightsHist[10];
  TH1D* fEtaWeightsHist[10][21][2];
  TH1D* fNvsCenCut[2][2]; //! ZDC mult cuts
  TH2F* fZDCESEMultWeightsHist[5];
  TH2F* fZDCESESpecWeightsHist[5];
  Bool_t fQAZDCCuts;
  Int_t fMinMulZN;
  Float_t fMaxDevZN;
  Float_t fZDCGainAlpha;
  
  ClassDef(AliAnalysisTaskCRC,11);
};

//================================================================================================================

#endif











