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
class AliFemtoDreamPairCleanerHists {
 public:
  AliFemtoDreamPairCleanerHists();
  AliFemtoDreamPairCleanerHists(int nTrackDecays,int nDecayDecays);
  virtual ~AliFemtoDreamPairCleanerHists();
  void FillDaughtersSharedTrack(int Hist,int counter) {
    fTrackDecays[Hist]->Fill(counter);
  };
  void FillDaughtersSharedDaughter(int Hist,int counter) {
    fDecayDecays[Hist]->Fill(counter);
  };
  TList *GetHistList(){return fOutput;};
  TH1F **fTrackDecays;
  TH1F **fDecayDecays;
  TList *fOutput;
  ClassDef(AliFemtoDreamPairCleanerHists,2)
};

#endif /* ALIFEMTODREAMPAIRCLEANERHISTS_H_ */
