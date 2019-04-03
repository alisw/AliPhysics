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
  AliFemtoDreamCollConfig();
  AliFemtoDreamCollConfig(const AliFemtoDreamCollConfig& config);
  AliFemtoDreamCollConfig(const char *name, const char *title);
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
    fdPhidEtaPlotsSmallK = doIt;
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
  void SetUseStravinskyMethod(bool use) {
    fStravinsky = use;
  }
  ;
  void SetZBins(std::vector<float> ZBins);
  void SetMultBins(std::vector<int> MultBins);
  void SetPDGCodes(std::vector<int> PDGCodes);
  void SetNBinsHist(std::vector<int> NBins);
  void SetMinKRel(std::vector<float> minKRel);
  void SetMaxKRel(std::vector<float> maxKRel);
  void SetCentBins(std::vector<float> CentBins);
  void SetmTdEtadPhiBins(std::vector<float> mTBins);
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
  void SetInvMassPairs(bool doIt) {
    fInvMassPairs = doIt;
  }
  ;
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
  bool GetDoStravinsky() {
    return fStravinsky;
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
    return (fdPhidEtaPlots&&fmTdEtadPhi);
  }
  bool GetInvMassPairs() {
    return fInvMassPairs;
  }
  ;
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
    return (fZVtxBins->GetEntries() - 1);
  }
  ;
  std::vector<int> GetMultBins();
  int GetNMultBins() {
    return fMultBins->GetEntries();
  }
  ;
  std::vector<int> GetPDGCodes();
  int GetNParticles() {
    return fPDGParticleSpecies->GetEntries();
  }
  ;
  int GetNParticleCombinations();
  std::vector<int> GetNBinsHist();
  std::vector<float> GetMinKRel();
  std::vector<float> GetMaxKRel();
  std::vector<float> GetCentBins();
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

  void SetDeltaEtaMax(float delta) {
    fDoDeltaEtaDeltaPhiCut = true;
    fDeltaEtaMax = delta;
  }
  float GetDeltaEtaMax() const { return fDeltaEtaMax; }

  void SetDeltaPhiMax(float delta) {
    fDoDeltaEtaDeltaPhiCut = true;
    fDeltaPhiMax = delta;
  }
  float GetDeltaPhiMax() const { return fDeltaPhiMax; }

  void DoDeltaEtaDeltaPhiCut(bool doIt) { fDoDeltaEtaDeltaPhiCut = doIt; }
  float GetDoDeltaEtaDeltaPhiCut() const { return fDoDeltaEtaDeltaPhiCut; }

 private:
  bool fMultBinning;            //
  bool fCentBinning;            //
  bool fkTBinning;              //
  bool fmTBinning;              //
  bool fMomentumResolution;     //
  bool fPhiEtaBinning;          //
  bool fdPhidEtaPlots;          //
  bool fdPhidEtaPlotsSmallK;    //
  bool fMixedEventStatistics;   //
  bool fGetTheControlSampel;    //
  bool fStravinsky;             //
  bool fInvMassPairs;           //
  bool fMinimalBookingME;       //
  bool fMinimalBookingSample;   //
  int fNumberRadii;             //
  TNtuple *fZVtxBins;           //
  TNtuple *fMultBins;           //
  TNtuple *fPDGParticleSpecies;  //
  TNtuple *fNBinsHists;         //
  TNtuple *fMinK_rel;           //
  TNtuple *fMaxK_rel;           //
  TNtuple *fCentBins;           //
  TNtuple *fmTBins;             //
  TNtuple *fWhichQAPairs;       //
  TNtuple *fClosePairRej;       //
  int fMixingDepth;             //
  int fSpinningDepth;			      //
  bool fkTCentrality;           //
  bool fmTdEtadPhi;             //
  AliFemtoDreamEvent::MultEstimator fEst;  //
  float fDeltaEtaMax;           //
  float fDeltaPhiMax;           //
  bool fDoDeltaEtaDeltaPhiCut;  //

ClassDef(AliFemtoDreamCollConfig,10)
  ;
};

#endif /* ALIFEMTODREAMCOLLCONFIG_H_ */
