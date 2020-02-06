/*
 * AliFemtoDreamCollConfig.h
 *
 *  Created on: Sep 13, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMCOLLCONFIG_H_
#define ALIFEMTODREAMCOLLCONFIG_H_
#include <vector>
#include "Rtypes.h"
#include "TNamed.h"
#include "TNtuple.h"
#include "AliFemtoDreamEvent.h"

class AliFemtoDreamCollConfig : public TNamed {
 public:
  enum UncorrelatedMode {
    kNone = 0,
    kPhiSpin = 1,
    kStravinsky = 2,
    kCorrelatedPhi = 3
  };
  AliFemtoDreamCollConfig();
  AliFemtoDreamCollConfig(const AliFemtoDreamCollConfig& config);
  AliFemtoDreamCollConfig(const char *name, const char *title, bool QACouts =
                              false);
  AliFemtoDreamCollConfig &operator=(const AliFemtoDreamCollConfig& config);
  virtual ~AliFemtoDreamCollConfig();
  void SetMultBinning(bool doIt) {
    fMultBinning = doIt;
  }
  ;
  void SetCentBinning(bool doIt) {
    fCentBinning = doIt;
  }
  ;
  void SetkTBinning(bool doIt) {
    fkTBinning = doIt;
  }
  ;
  void SetmTBinning(bool doIt) {
    fmTBinning = doIt;
  }
  ;
  void SetkTandMultBinning(bool doIt) {
    fkTandMultBinning = doIt;
  }
  ;
  void SetPtQA(bool doIt) {
    fPtQA = doIt;
  }
  void SetMassQA(bool doIt) {
    fMassQA = doIt;
  }
  void SetMomentumResolution(bool doIt) {
    fMomentumResolution = doIt;
  }
  ;
  void SetkTCentralityBinning(bool doIt) {
    fkTCentrality = doIt;
  }
  ;
  void SetPhiEtaBinnign(bool doIt) {
    fPhiEtaBinning = doIt;
    fNumberRadii = 9;
  }
  ;
  void SetdPhidEtaPlots(bool doIt) {
    fdPhidEtaPlots = doIt;
  }
  void SetdPhidEtaPlotsSmallK(bool doIt) {
    fdPhidEtaPlotsSmallK = doIt;
  }
  ;
  void SetUseEventMixing(bool use) {
    fMixedEventStatistics = use;
  }
  ;
  void SetUsePhiSpinning(bool use) {
    fGetTheControlSampel = use;
  }
  ;
  void SetControlMethod(AliFemtoDreamCollConfig::UncorrelatedMode mode) {
    fMode = mode;
  }
  ;
    void SetAncestors(bool doIt) {
    fAncestors = doIt;
  }
  void SetZBins(std::vector<float> ZBins);
  void SetMultBins(std::vector<int> MultBins);
  void SetPDGCodes(std::vector<int> PDGCodes);
  void SetNBinsHist(std::vector<int> NBins);
  void SetMinKRel(std::vector<float> minKRel);
  void SetMaxKRel(std::vector<float> maxKRel);
  void SetCentBins(std::vector<int> CentBins);
  void SetmTdEtadPhiBins(std::vector<float> mTBins);
  //TODO: should be renamed since besides the QA it also specifies the
  // number of tracks to compare when doing the CPR cut
  void SetExtendedQAPairs(std::vector<int> whichPairs);
  void SetClosePairRejection(std::vector<bool> whichPairs);
  void SetMixingDepth(int MixingDepth) {
    fMixingDepth = MixingDepth;
  }
  ;
  void SetSpinningDepth(int SpinningDepth) {
    fSpinningDepth = SpinningDepth;
  }
  ;
  void SetCorrelationRange(float CorrRange) {
    fCorrelationRange = CorrRange;
  }
  void SetMinimalBookingME(bool doIt) {
    fMinimalBookingME = doIt;
  }
  ;
  void SetMinimalBookingSample(bool doIt) {
    fMinimalBookingSample = doIt;
  }
  ;
  void SetMultiplicityEstimator(AliFemtoDreamEvent::MultEstimator est) {
    fEst = est;
  }

  bool GetDoMultBinning() {
    return fMultBinning;
  }
  ;
  bool GetDoCentBinning() {
    return fCentBinning;
  }
  ;
  bool GetDokTBinning() {
    return fkTBinning;
  }
  ;
  bool GetDomTBinning() {
    return fmTBinning;
  }
  ;
  bool GetDokTandMultBinning() {
    return fkTandMultBinning;
  }
  ;
  bool GetDoPtQA() {
    return fPtQA;
  }
  bool GetDoMassQA() {
    return fMassQA;
  }
  bool GetDoMomResolution() {
    return fMomentumResolution;
  }
  ;
  bool GetDoPhiEtaBinning() {
    return fPhiEtaBinning;
  }
  ;
  bool GetDokTCentralityBinning() {
    return fkTCentrality;
  }
  ;
  bool GetUseEventMixing() {
    return fMixedEventStatistics;
  }
  ;
  bool GetUsePhiSpinning() {
    return fGetTheControlSampel;
  }
  ;
  bool GetDoAncestorsPlots() {
    return fAncestors;
  }
  ;
  AliFemtoDreamCollConfig::UncorrelatedMode GetControlMode() {
    return fMode;
  }
  ;
  bool GetdPhidEtaPlots() {
    return fdPhidEtaPlots;
  }
  bool GetdPhidEtaPlotsSmallK() {
    return fdPhidEtaPlotsSmallK;
  }
  ;
  bool GetdPhidEtamTPlots() {
    return (fdPhidEtaPlots && fmTdEtadPhi);
  }
  bool GetMinimalBookingME() {
    return fMinimalBookingME;
  }
  ;
  bool GetMinimalBookingSample() {
    return fMinimalBookingSample;
  }
  ;
  AliFemtoDreamEvent::MultEstimator GetMultiplicityEstimator() {
    return fEst;
  }
  int GetNRadii() {
    return fNumberRadii;
  }
  ;
  std::vector<float> GetZVtxBins();
  int GetNZVtxBins() {
    return ((int) fZVtxBins.size() - 1);
  }
  ;
  std::vector<int> GetMultBins();
  int GetNMultBins() {
    return (int) fMultBins.size();
  }
  ;
  std::vector<int> GetPDGCodes();
  int GetNParticles() {
//    return fPDGParticleSpecies->GetEntries();
    return (int) fPDGParticleSpecies.size();
  }
  ;
  int GetNParticleCombinations();
  std::vector<int> GetNBinsHist();
  std::vector<float> GetMinKRel();
  std::vector<float> GetMaxKRel();
  std::vector<int> GetCentBins();
  std::vector<float> GetmTBins();
  std::vector<unsigned int> GetWhichPairs();
  std::vector<bool> GetClosePairRej();
  std::vector<float> GetStandardmTBins();
  std::vector<int> GetStandardPairs();
  std::vector<bool> GetStandardPairRejection();
  std::vector<bool> GetAllPairRejection();
  int GetMixingDepth() {
    return fMixingDepth;
  }
  ;
  int GetSpinningDepth() {
    return fSpinningDepth;
  }
  ;
  float GetCorrelationRange() {
    return fCorrelationRange;
  }
  ;
  void SetDeltaEtaMax(float delta) {
    fDoDeltaEtaDeltaPhiCut = true;
    fDeltaEtaMax = delta;
  }
  float GetDeltaEtaMax() const {
    return fDeltaEtaMax;
  }

  void SetDeltaPhiMax(float delta) {
    fDoDeltaEtaDeltaPhiCut = true;
    fDeltaPhiMax = delta;
  }
  float GetDeltaPhiMax() const {
    return fDeltaPhiMax;
  }
  float GetSqDeltaPhiEtaMax() const {
    return fDeltaEtaMax * fDeltaEtaMax + fDeltaPhiMax * fDeltaPhiMax;
  }
  ;
  void DoDeltaEtaDeltaPhiCut(bool doIt) {
    fDoDeltaEtaDeltaPhiCut = doIt;
  }
  bool GetDoDeltaEtaDeltaPhiCut() const {
    return fDoDeltaEtaDeltaPhiCut;
  }

 private:
  bool fMultBinning;            //
  bool fCentBinning;            //
  bool fkTBinning;              //
  bool fmTBinning;              //
  bool fkTandMultBinning;	//
  bool fPtQA;                   //
  bool fMassQA;                 //
  bool fMomentumResolution;     //
  bool fPhiEtaBinning;          //
  bool fdPhidEtaPlots;          //
  bool fdPhidEtaPlotsSmallK;    //
  bool fMixedEventStatistics;   //
  bool fGetTheControlSampel;    //
  bool fAncestors;              //
  AliFemtoDreamCollConfig::UncorrelatedMode fMode;  //
  bool fMinimalBookingME;       //
  bool fMinimalBookingSample;   //
  int fNumberRadii;             //
  std::vector<float> fZVtxBins;           //
  std::vector<int> fMultBins;           //
  std::vector<int> fPDGParticleSpecies;  //
  std::vector<int> fNBinsHists;         //
  std::vector<float> fMinK_rel;           //
  std::vector<float> fMaxK_rel;           //
  std::vector<int> fCentBins;           //
  std::vector<float> fmTBins;             //
  std::vector<unsigned int> fWhichQAPairs;       //
  std::vector<bool> fClosePairRej;       //
  int fMixingDepth;             //
  int fSpinningDepth;			      //
  float fCorrelationRange;	      //
  bool fkTCentrality;           //
  bool fmTdEtadPhi;             //
  AliFemtoDreamEvent::MultEstimator fEst;  //
  float fDeltaEtaMax;           //
  float fDeltaPhiMax;           //
  bool fDoDeltaEtaDeltaPhiCut;  //
  bool fCoutVariables;
  ClassDef(AliFemtoDreamCollConfig,16);
};

#endif /* ALIFEMTODREAMCOLLCONFIG_H_ */
