/*
 * AliFemtoDreamPairCleaner.h
 *
 *  Created on: 8 Nov 2017
 *      Author: bernhardhohlweger
 */

#ifndef ALIFEMTODREAMPAIRCLEANER_H_
#define ALIFEMTODREAMPAIRCLEANER_H_
#include <vector>
#include "Rtypes.h"
#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamPairCleanerHists.h"
#include "TDatabasePDG.h"
#include "TVector3.h"

class AliFemtoDreamPairCleaner {
 public:
  AliFemtoDreamPairCleaner();
  AliFemtoDreamPairCleaner(const AliFemtoDreamPairCleaner& cleaner);
  AliFemtoDreamPairCleaner(int nTrackDecayChecks, int nDecayDecayChecks,
                           bool MinimalBooking);
  AliFemtoDreamPairCleaner& operator=(const AliFemtoDreamPairCleaner& cleaner);
  virtual ~AliFemtoDreamPairCleaner();
  void CleanTrackAndDecay(std::vector<AliFemtoDreamBasePart> *Tracks,
                          std::vector<AliFemtoDreamBasePart> *Decay,
                          int histnumber);
  void CleanDecay(std::vector<AliFemtoDreamBasePart> *Decay, int histnumber);
  void CleanDecayAndDecay(std::vector<AliFemtoDreamBasePart> *Decay1,
                          std::vector<AliFemtoDreamBasePart> *Decay2,
                          int histnumber);
  void FillInvMassPair(std::vector<AliFemtoDreamBasePart> &Part1, int PDGCode1,
                       std::vector<AliFemtoDreamBasePart> &Part2, int PDGCode2,
                       int histnumber);
  void StoreParticle(std::vector<AliFemtoDreamBasePart> Particles);
  TList* GetHistList() {
    return fHists->GetHistList();
  }
  ;
  std::vector<std::vector<AliFemtoDreamBasePart>>& GetCleanParticles() {
    return fParticles;
  }
  ;
  void ResetArray();
  float RelativePairMomentum(TVector3 Part1Momentum, int PDGPart1,
                             TVector3 Part2Momentum, int PDGPart2);
  int GetCounter() const {return fCounter;};
 private:
  double InvMassPair(TVector3 Part1, int PDG1, TVector3 Part2, int PDG2);
  double E2(int pdgCode, double Ptot2);
  bool fMinimalBooking;
  int fCounter;
  std::vector<std::vector<AliFemtoDreamBasePart>> fParticles;
  AliFemtoDreamPairCleanerHists *fHists;ClassDef(AliFemtoDreamPairCleaner,3)
};

inline double AliFemtoDreamPairCleaner::E2(int pdgCode, double Ptot2) {
  double mass = TDatabasePDG::Instance()->GetParticle(pdgCode)->Mass();
  return mass * mass + Ptot2;
}

inline double AliFemtoDreamPairCleaner::InvMassPair(TVector3 Part1, int PDG1,
                                                    TVector3 Part2, int PDG2) {
  double EPart1 = ::sqrt(E2(PDG1, Part1.Mag2()));
  double EPart2 = ::sqrt(E2(PDG2, Part2.Mag2()));
  return ::sqrt((EPart1 + EPart2) * (EPart1 + EPart2) - (Part1 + Part2).Mag2());
}

#endif /* ALIFEMTODREAMPAIRCLEANER_H_ */
