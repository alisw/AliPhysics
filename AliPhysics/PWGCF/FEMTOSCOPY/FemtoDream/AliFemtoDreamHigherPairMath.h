/*
 * AliFemtoDreamHigherPairMath.h
 *
 *  Created on: Jul 3, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIFEMTODREAMHIGHERPAIRMATH_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIFEMTODREAMHIGHERPAIRMATH_H_
#include "Rtypes.h"
#include "AliLog.h"
#include "TRandom3.h"
#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamCorrHists.h"
#include <vector>
class AliFemtoDreamHigherPairMath {
 public:
  AliFemtoDreamHigherPairMath(AliFemtoDreamCollConfig *conf);
  virtual ~AliFemtoDreamHigherPairMath();
  AliFemtoDreamHigherPairMath(const AliFemtoDreamHigherPairMath& samp);
  AliFemtoDreamHigherPairMath& operator=(const AliFemtoDreamHigherPairMath& math);
  void SetBField(float bielefield) {
    fBField = bielefield;
  }
  int RandomNumber() {
    AliFatal("Why would you use this? \n");
    //decided by a fair xkcd dice roll.
    return 4;
  }
  ;
  bool PassesPairSelection(AliFemtoDreamBasePart* part1,
                           AliFemtoDreamBasePart* part2, bool Recalculate);
  void RecalculatePhiStar(AliFemtoDreamBasePart *part);
  float FillSameEvent(int iHC, int Mult, float cent, TVector3 Part1Momentum,
                      int PDGPart1, TVector3 Part2Momentum, int PDGPart2);
  void MassQA(int iHC, float RelK, float massOne, float massTwo);
  void SEDetaDPhiPlots(int iHC, AliFemtoDreamBasePart& part1, int PDGPart1,
                       AliFemtoDreamBasePart& part2, int PDGPart2,
                       float RelativeK, bool recalculate) {
    SEDetaDPhiPlots(iHC, &part1, PDGPart1, &part2, PDGPart2, RelativeK,
                    recalculate);
  }
  void SEDetaDPhiPlots(int iHC, AliFemtoDreamBasePart* part1, int PDGPart1,
                       AliFemtoDreamBasePart* part2, int PDGPart2,
                       float RelativeK, bool recalculate);
  void FillPairCounterSE(int iHC, unsigned int sizePartOne,
                         unsigned int sizePartTwo) {
    fHists->FillPartnersSE(iHC, sizePartOne, sizePartTwo);
  }
  float FillMixedEvent(int iHC, int Mult, float cent, TVector3 Part1Momentum,
                       int PDGPart1, TVector3 Part2Momentum, int PDGPart2,
                       AliFemtoDreamCollConfig::UncorrelatedMode mode);
  void MEMomentumResolution(int iHC, AliFemtoDreamBasePart* part1, int PDGPart1,
                            AliFemtoDreamBasePart* part2, int PDGPart2,
                            float RelativeK);
  void MEDetaDPhiPlots(int iHC, AliFemtoDreamBasePart& part1, int PDGPart1,
                       AliFemtoDreamBasePart& part2, int PDGPart2,
                       float RelativeK, bool recalculate) {
    MEDetaDPhiPlots(iHC, &part1, PDGPart1, &part2, PDGPart2, RelativeK,
                    recalculate);
  }
  void MEDetaDPhiPlots(int iHC, AliFemtoDreamBasePart* part1, int PDGPart1,
                       AliFemtoDreamBasePart* part2, int PDGPart2,
                       float RelativeK, bool recalculate);
  void FillPairCounterME(int iHC, unsigned int sizePartOne,
                         unsigned int sizePartTwo) {
    fHists->FillPartnersME(iHC, sizePartOne, sizePartTwo);
  }
  TList* GetHistList() {
    return fHists->GetHistList();
  }
  TList* GetQAHists() {
    return fHists->GetQAHists();
  }
  TString ClassName() {
    return "HighMaths";
  }
  ;
 private:
  float RelativePairMomentum(TLorentzVector &PartOne, TLorentzVector &PartTwo);
  float RelativePairkT(TLorentzVector &PartOne, TLorentzVector &PartTwo);
  float RelativePairmT(TLorentzVector &PartOne, TLorentzVector &PartTwo);
  void DeltaEtaDeltaPhi(int Hist, AliFemtoDreamBasePart *part1,
                        AliFemtoDreamBasePart *part2, bool SEorME, float relk,
                        bool recalculate);
  AliFemtoDreamCorrHists *fHists;
  std::vector<unsigned int> fWhichPairs;
  float fBField;
  float fDeltaPhiEtaMax;
  TRandom3 fRandom;
  double fPi;

};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIFEMTODREAMHIGHERPAIRMATH_H_ */
