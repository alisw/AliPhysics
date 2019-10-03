/*
 * AliFemtoPairCleanerHists.h
 *
 *  Created on: 22 Dec 2017
 *      Author: bernhardhohlweger
 */

#ifndef ALIFEMTODREAMPAIRCLEANERHISTS_H_
#define ALIFEMTODREAMPAIRCLEANERHISTS_H_
#include "Rtypes.h"
#include "TH1F.h"
#include "TList.h"
#include "TNtuple.h"

class AliFemtoDreamPairCleanerHists {
 public:
  AliFemtoDreamPairCleanerHists();
  AliFemtoDreamPairCleanerHists(const AliFemtoDreamPairCleanerHists& hists);
  AliFemtoDreamPairCleanerHists(int nTrackDecays, int nDecayDecays);
  AliFemtoDreamPairCleanerHists& operator=(
      const AliFemtoDreamPairCleanerHists& hists);
  virtual ~AliFemtoDreamPairCleanerHists();
  void FillDaughtersSharedTrack(int Hist, int counter) {
    fTrackDecays[Hist]->Fill(counter);
  }
  ;
  void FillDaughtersSharedDaughter(int Hist, int counter) {
    fDecayDecays[Hist]->Fill(counter);
  }
  ;
  void FillPairInvMass(int Hist, float mass) {
    fPairInvMass[Hist]->Fill(mass);
  }
  TList *GetHistList() {
    return fOutput;
  }
  ;
  TH1F **fTrackDecays;
  TH1F **fDecayDecays;
  TH1F **fPairInvMass;
  TList *fOutput;
  ClassDef(AliFemtoDreamPairCleanerHists,2)
};

#endif /* ALIFEMTODREAMPAIRCLEANERHISTS_H_ */
