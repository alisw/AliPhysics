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
class AliFemtoDreamPairCleaner {
 public:
  AliFemtoDreamPairCleaner();
  AliFemtoDreamPairCleaner(
      int nTrackDecayChecks, int nDecayDecayChecks,bool MinimalBooking);
  virtual ~AliFemtoDreamPairCleaner();
  void CleanTrackAndDecay(std::vector<AliFemtoDreamBasePart> *Tracks,
                          std::vector<AliFemtoDreamBasePart> *Decay,
                          int histnumber);
  void CleanDecay(std::vector<AliFemtoDreamBasePart> *Decay,int histnumber);
  void CleanDecayAndDecay(std::vector<AliFemtoDreamBasePart> *Decay1,
                          std::vector<AliFemtoDreamBasePart> *Decay2,
                          int histnumber);
  void StoreParticle(std::vector<AliFemtoDreamBasePart> Particles);
  TList* GetHistList(){return fHists->GetHistList();};
  std::vector<std::vector<AliFemtoDreamBasePart>>& GetCleanParticles()
      {return fParticles;};
  void ResetArray();
 private:
  bool fMinimalBooking;
  std::vector<std::vector<AliFemtoDreamBasePart>> fParticles;
  AliFemtoDreamPairCleanerHists *fHists;
  ClassDef(AliFemtoDreamPairCleaner,2)
};

#endif /* ALIFEMTODREAMPAIRCLEANER_H_ */
