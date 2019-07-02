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
  void CorrelatedSample(std::vector<AliFemtoDreamBasePart> &part1,
                        int &PDGPart1,
                        std::vector<AliFemtoDreamBasePart> &part2,
                        int &PDGPart2, bool SameParticle, int Mult,
                        int HistCounter);
  void UncorrelatedSample(std::vector<AliFemtoDreamBasePart> &part1,
                          int &PDGPart1,
                          std::vector<AliFemtoDreamBasePart> &part2,
                          int &PDGPart2, bool SameParticle, int Mult,
                          int HistCounter);
  void PhiSpinning(std::vector<AliFemtoDreamBasePart> &part1, int PDGPart1,
                   std::vector<AliFemtoDreamBasePart> &part2, int PDGPart2,
                   bool SameParticle, int Mult, int HistCounter);
  void LimitedPhiSpinning(std::vector<AliFemtoDreamBasePart> &part1,
                          int PDGPart1,
                          std::vector<AliFemtoDreamBasePart> &part2,
                          int PDGPart2, bool SameParticle, int Mult,
                          int HistCounter);
  void Randomizer(std::vector<AliFemtoDreamBasePart*> &part);
  float RelativePairMomentum(TVector3 Part1Momentum, int PDGPart1,
                             TVector3 Part2Momentum, int PDGPart2,
                             AliFemtoDreamCollConfig::UncorrelatedMode mode = AliFemtoDreamCollConfig::kNone);
  float ComputeDeltaEta(AliFemtoDreamBasePart *part1,
                        AliFemtoDreamBasePart *part2);
  float ComputeDeltaPhi(AliFemtoDreamBasePart *part1,
                        AliFemtoDreamBasePart *part2);
  float ComputeDeltaEta(AliFemtoDreamBasePart &part1,
                        AliFemtoDreamBasePart &part2);
  float ComputeDeltaPhi(AliFemtoDreamBasePart &part1,
                        AliFemtoDreamBasePart &part2);
  int FindBin(float Multiplicity);
  AliFemtoDreamCorrHists *fHists;
  std::vector<int> fPDGParticleSpecies;
  std::vector<unsigned int> fWhichPairs;
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

ClassDef(AliFemtoDreamControlSample, 3)
};

#endif
