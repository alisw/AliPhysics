#ifndef ALIFEMTODREAMCONTROLSAMPLE_H_
#define ALIFEMTODREAMCONTROLSAMPLE_H_

#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamCorrHists.h"
#include "Rtypes.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "vector"

class AliFemtoDreamControlSample {
 public:
  AliFemtoDreamControlSample();
  AliFemtoDreamControlSample(const AliFemtoDreamControlSample& samp);
  AliFemtoDreamControlSample(AliFemtoDreamCollConfig *conf,
                             bool minimalBooking = false);
  AliFemtoDreamControlSample& operator=(const AliFemtoDreamControlSample& samp);
  virtual ~AliFemtoDreamControlSample();

  void SetEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
                float mult);

  float RelativePairMomentum(TVector3 Part1Momentum, int PDGPart1,
                             TVector3 Part2Momentum, int PDGPart2, bool random =
                                 false);
  float ComputeDeltaEta(AliFemtoDreamBasePart &part1,
                        AliFemtoDreamBasePart &part2);
  float ComputeDeltaPhi(AliFemtoDreamBasePart &part1,
                        AliFemtoDreamBasePart &part2);
  int FindBin(float Multiplicity);
  TString ClassName() {
    return "Control sample leaking";
  }
  ;

  TList* GetHistList() {
    return fHists->GetHistList();
  }
  TList* GetQAList() {
    return fHists->GetQAHists();
  }

 private:

  AliFemtoDreamCorrHists *fHists;
  std::vector<int> fPDGParticleSpecies;
  std::vector<int> fMultBins;
  TRandom3 fRandom;
  double fPi;
  int fSpinningDepth;
  bool fStravinsky;
  float fDeltaEtaMax;
  float fDeltaPhiMax;

ClassDef(AliFemtoDreamControlSample, 2)
};

#endif
