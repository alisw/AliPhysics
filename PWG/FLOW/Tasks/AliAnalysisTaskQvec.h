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

#ifndef AliAnalysisTaskQvec_H
#define AliAnalysisTaskQvec_H

#include "AliAnalysisTaskSE.h"
#include "AliFlowCommonConstants.h"
#include "TGrid.h"

class TString;
class TList;
class TProfile2D;
class AliFlowEventSimple;
class AliFlowEvent;
class AliFlowAnalysisQvec;

//==============================================================================================================

class AliAnalysisTaskQvec : public AliAnalysisTaskSE{
public:
  AliAnalysisTaskQvec();
  AliAnalysisTaskQvec(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskQvec(){};

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  virtual void NotifyRun();
  
  // Common:
  void SetExactNoRPs(Int_t const enr) {this->fExactNoRPs = enr;};
  Int_t GetExactNoRPs() const {return this->fExactNoRPs;};
  
  void SetWeightsList(TList* const kList) {this->fWeightsList = (TList*)kList->Clone();}; // useful
  TList* GetWeightsList() const {return this->fWeightsList;}; // useful

  // Multiparticle correlations vs multiplicity:
  
  // Particle weights:
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  void SetUsePhiEtaCuts(Bool_t const uPhiEtaW) {this->fUsePhiEtaCuts = uPhiEtaW;};
  Bool_t GetUsePhiEtaCuts() const {return this->fUsePhiEtaCuts;};
  void SetUseZDCESEMulWeights(Bool_t const uPhiEtaW) {this->fUseZDCESEMulWeights = uPhiEtaW;};
  Bool_t GetUseZDCESEMulWeights() const {return this->fUseZDCESEMulWeights;};
  void SetUseZDCESESpecWeights(Bool_t const uPhiEtaW) {this->fUseZDCESESpecWeights = uPhiEtaW;};
  Bool_t GetUseZDCESESpecWeights() const {return this->fUseZDCESESpecWeights;};
  void SetCutMultiplicityOutliers(Bool_t const uPhiEtaW) {this->fCutMultiplicityOutliers = uPhiEtaW;};
  Bool_t GetCutMultiplicityOutliers() const {return this->fCutMultiplicityOutliers;};

  // Event weights:
  void SetMultiplicityWeight(const char *multiplicityWeight) {*this->fMultiplicityWeight = multiplicityWeight;};

  // # of bins for correlation axis in fDistributions[4], fCorrelation2468VsMult[4] and fCorrelationProduct2468VsMult[1]

  // Boundaries for distributions of correlations:

  // min and max values of correlation products:

  // min and max values of QvectorTerms:
  
  // bootstrap:

  // Charge-Rapidity Correlations:
  void SetCalculateCRC(Bool_t const cCRC) {this->fCalculateCRC = cCRC;};
  Bool_t GetCalculateCRC() const {return this->fCalculateCRC;};
  void SetCalculateCME(Bool_t const cCRC) {this->fCalculateCME = cCRC;};
  Bool_t GetCalculateCME() const {return this->fCalculateCME;};
  void SetStoreZDCQVecVtxPos(Bool_t const cCRC) {this->fStoreZDCQVecVtxPos = cCRC;};
  Bool_t GetStoreZDCQVecVtxPos() const {return this->fStoreZDCQVecVtxPos;};
  void SetUseVZERO(Bool_t const cCRC) {this->fUseVZERO = cCRC;};
  Bool_t GetUseVZERO() const {return this->fUseVZERO;};
  void SetUseZDC(Bool_t const cCRC) {this->fUseZDC = cCRC;};
  Bool_t GetUseZDC() const {return this->fUseZDC;};
  void SetRemoveSplitMergedTracks(Bool_t const uPhiEtaW) {this->fRemoveSplitMergedTracks = uPhiEtaW;};
  Bool_t GetRemoveSplitMergedTracks() const {return this->fRemoveSplitMergedTracks;};
  void SetRecenterZDC(Bool_t const cCRC) {this->fRecenterZDC = cCRC;};
  Bool_t GetRecenterZDC() const {return this->fRecenterZDC;};
  void SetDivSigma(Bool_t const cCRC) {this->fDivSigma = cCRC;};
  Bool_t GetDivSigma() const {return this->fDivSigma;};
  void SetInvertZDC(Bool_t const cCRC) {this->fInvertZDC = cCRC;};
  Bool_t GetInvertZDC() const {return this->fInvertZDC;};
  void SetTestSin(Bool_t const cCRC) {this->fCRCTestSin = cCRC;};
  Bool_t GetTestSin() const {return this->fCRCTestSin;};
  void SetZDCESEList(TList* const kList) {this->fZDCESEList = (TList*)kList->Clone();};
  TList* GetZDCESEList() const {return this->fZDCESEList;};
  void SetCRCZDCCalibList(TList* const wlist) {this->fCRCZDCCalibList = (TList*)wlist->Clone();}
  TList* GetCRCZDCCalibList() const {return this->fCRCZDCCalibList;}
  void SetCRCZDC2DCutList(TList* const wlist) {this->fCRCZDC2DCutList = (TList*)wlist->Clone();}
  void SetCRCVZEROCalibList(TList* const wlist) {this->fCRCVZEROCalibList = (TList*)wlist->Clone();}
  TList* GetCRCVZEROCalibList() const {return this->fCRCVZEROCalibList;}
  void SetCRCZDCResList(TList* const wlist) {this->fCRCZDCResList = (TList*)wlist->Clone();}
  TList* GetCRCZDCResList() const {return this->fCRCZDCResList;}
  void SetnCenBin(Int_t const n) {this->fnCenBin = n;};
  Int_t GetnCenBin() const {return this->fnCenBin;};
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
  void SetCenWeightsHist(TH1D* const n) {this->fCenWeightsHist = n;};
  TH1D* GetCenWeightsHist() const {return this->fCenWeightsHist;};
  void SetRefMultRbRPro(TProfile2D* const n) {this->fRefMultRbRPro = n;};
  void SetAvEZDCRbRPro(TProfile2D* const A, TProfile2D* const B) {this->fAvEZDCCRbRPro = A; this->fAvEZDCARbRPro = B;};
  void SetPtWeightsHist(TH1D* const n, Int_t c) {this->fPtWeightsHist[c] = n;};
  TH1D* GetPtWeightsHist(Int_t c) const {return this->fPtWeightsHist[c];};
  void SetZDCESEMultWeightsHist(TH2F* const n, Int_t h) {this->fZDCESEMultWeightsHist[h] = n;};
  TH2F* GetZDCESEMultWeightsHist(Int_t h) const {return this->fZDCESEMultWeightsHist[h];};
  void SetZDCESESpecWeightsHist(TH2F* const n, Int_t h) {this->fZDCESESpecWeightsHist[h] = n;};
  TH2F* GetZDCESESpecWeightsHist(Int_t h) const {return this->fZDCESESpecWeightsHist[h];};
  void SetQAZDCCuts(Bool_t const cCRC) {this->fQAZDCCuts = cCRC;};
  Bool_t GetQAZDCCuts() const {return this->fQAZDCCuts;};
  void SetMinMulZN(Int_t weights) {this->fMinMulZN = weights;};
  Int_t GetMinMulZN() const {return this->fMinMulZN;};
  void SetUseTracklets(Bool_t const cCRC) {this->fUseTracklets = cCRC;};
  void SetCRCEtaRange(Double_t const etamin, Double_t const etamax) {this->fCRCEtaMin = etamin; this->fCRCEtaMax = etamax;};
  
  //@Shi set histogram for recentering
  void SetZDCCalibListFinalCommonPart(TList* const kList) {this->fZDCCalibListFinalCommonPart = (TList*)kList->Clone();};
  TList* GetZDCCalibListFinalCommonPart() const {return this->fZDCCalibListFinalCommonPart;};
  
private:
  AliAnalysisTaskQvec(const AliAnalysisTaskQvec& aatqc);
  AliAnalysisTaskQvec& operator=(const AliAnalysisTaskQvec& aatqc);

  AliFlowEvent *fEvent;         // the input event
  AliFlowAnalysisQvec *fQC;            // CRC object
  TList *fListHistos;                 // collection of output
  // Common:
  Int_t fExactNoRPs;                     // when shuffled, select only this number of RPs for the analysis
  
  // Particle weights:
  Bool_t fUseParticleWeights;         // use any particle weights
  Bool_t fUsePtWeights;               // use pt weights
  Bool_t fUsePhiEtaCuts;              // use phi,eta cuts (for NUA)
  Bool_t fUseZDCESEMulWeights;        // use ZDC-ESE mult. weights
  Bool_t fUseZDCESESpecWeights;       // use ZDC-ESE mult. weights
  Bool_t fCutMultiplicityOutliers;    // cut on reference multiplicity
  
  TList *fWeightsList;                // list with weights
  // Event weights:
  TString *fMultiplicityWeight;       // event-by-event weights for multiparticle correlations ("combinations","unit" or "multiplicity")
  
  // Boundaries for distributions of correlations:

  // Bootstrap:
  
  // Charge-Eta Asymmetry
  Bool_t fCalculateCRC; // calculate CRC quantities
  Bool_t fCalculateCME;
  Bool_t fStoreZDCQVecVtxPos;
  Bool_t fUseVZERO;
  Bool_t fUseZDC;
  Double_t fCenBinWidth;
  TString fDataSet;
  TString fInteractionRate;
  TString fSelectCharge;
  TString fPOIExtraWeights;
  TH1D* fCenWeightsHist;
  TProfile2D *fRefMultRbRPro;
  TProfile2D *fAvEZDCCRbRPro;
  TProfile2D *fAvEZDCARbRPro;
  TH1D* fPtWeightsHist[10];
  TH2F* fZDCESEMultWeightsHist[5];
  TH2F* fZDCESESpecWeightsHist[5];
  Double_t fCRCEtaMin;
  Double_t fCRCEtaMax;
  
  Bool_t fRemoveSplitMergedTracks;
  Bool_t fRecenterZDC;
  Bool_t fDivSigma;
  Bool_t fInvertZDC;
  Bool_t fCRCTestSin;
  Int_t fnCenBin;
  TList *fCRCZDCCalibList; // ZDC calibration
  TList *fCRCZDC2DCutList; // ZDC calibration
  TList *fCRCVZEROCalibList; // ZDC calibration
  TList *fCRCZDCResList; // ZDC rescaling
  TList *fZDCESEList;       // list with weights
  Bool_t fQAZDCCuts;
  Bool_t fUseTracklets;
  Int_t fMinMulZN;

  //@Shi ZDC calib recenter TList
  TList *fZDCCalibListFinalCommonPart; //
  
  ClassDef(AliAnalysisTaskQvec,1);
};

//================================================================================================================

#endif
