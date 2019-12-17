#ifndef ALIFEMTODREAMCONTROLSAMPLE_H_
#define ALIFEMTODREAMCONTROLSAMPLE_H_

#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamCorrHists.h"
#include "AliFemtoDreamHigherPairMath.h"
#include "AliFemtoDreamEvent.h"
#include "Rtypes.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "vector"

class AliFemtoDreamControlSample {
 public:
  AliFemtoDreamControlSample();
  AliFemtoDreamControlSample(const AliFemtoDreamControlSample& samp);
  AliFemtoDreamControlSample(AliFemtoDreamCollConfig *conf);
  AliFemtoDreamControlSample& operator=(const AliFemtoDreamControlSample& samp);
  virtual ~AliFemtoDreamControlSample();

  void SetEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
                AliFemtoDreamEvent *evt);
  TString ClassName() {
    return "Control sample leaking";
  }
  ;
  TList* GetHistList() {
    return fHigherMath->GetHistList();
  }
  TList* GetQAList() {
    return fHigherMath->GetQAHists();
  }

 private:
  void CorrelatedSample(std::vector<AliFemtoDreamBasePart> &part1,
                        int &PDGPart1,
                        std::vector<AliFemtoDreamBasePart> &part2,
                        int &PDGPart2, bool SameParticle, int HistCounter);
  void UncorrelatedSample(std::vector<AliFemtoDreamBasePart> &part1,
                          int &PDGPart1,
                          std::vector<AliFemtoDreamBasePart> &part2,
                          int &PDGPart2, bool SameParticle, int HistCounter);
  void PhiSpinning(std::vector<AliFemtoDreamBasePart> &part1, int PDGPart1,
                   std::vector<AliFemtoDreamBasePart> &part2, int PDGPart2,
                   bool SameParticle, int HistCounter);
  void LimitedPhiSpinning(std::vector<AliFemtoDreamBasePart> &part1,
                          int PDGPart1,
                          std::vector<AliFemtoDreamBasePart> &part2,
                          int PDGPart2, bool SameParticle, int HistCounter);
  void Randomizer(std::vector<AliFemtoDreamBasePart*> &part);
  int FindBin(float Multiplicity);
  AliFemtoDreamHigherPairMath* fHigherMath;
  std::vector<int> fPDGParticleSpecies;
  std::vector<bool> fRejPairs;
  std::vector<int> fMultBins;
  TRandom3 fRandom;
  double fPi;
  AliFemtoDreamCollConfig::UncorrelatedMode fmode;
  int fSpinningDepth;
  double fCorrelationRange;
  float fDeltaEtaMax;
  float fDeltaPhiMax;
  bool fDoDeltaEtaDeltaPhiCut;
  int fMult;
  float fCent;ClassDef(AliFemtoDreamControlSample, 4)
};

#endif
