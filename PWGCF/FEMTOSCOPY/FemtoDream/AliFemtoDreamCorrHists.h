/*
 * AliFemtoDreamCorrHists.h
 *
 *  Created on: Sep 12, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMCORRHISTS_H_
#define ALIFEMTODREAMCORRHISTS_H_
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"

#include "AliFemtoDreamCollConfig.h"

class AliFemtoDreamCorrHists {
 public:
  AliFemtoDreamCorrHists();
  AliFemtoDreamCorrHists(AliFemtoDreamCollConfig *conf);
  virtual ~AliFemtoDreamCorrHists();
  void FillSameEventDist(int i,double RelK){fSameEventDist[i]->Fill(RelK);};
  void FillMixedEventDist(int i,double RelK){fMixedEventDist[i]->Fill(RelK);};
  void FillPartnersSE(int hist,int nPart1,int nPart2){
    fPairCounterSE[hist]->Fill(nPart1,nPart2);
  }
  void FillPartnersME(int hist,int nPart1,int nPart2){
    fPairCounterME[hist]->Fill(nPart1,nPart2);
  }
  TList* GetHistList(){return fResults;};
  TList* GetQAHists(){return fQA;};
 private:
  TList         *fQA;
  TList         *fResults;
  TList         **fPairs;
  TList         **fPairQA;
  TH1F          **fSameEventDist;
  TH2F          **fPairCounterSE;
  TH1F          **fMixedEventDist;
  TH2F          **fPairCounterME;

  ClassDef(AliFemtoDreamCorrHists,1);
};

#endif /* ALIFEMTODREAMCORRHISTS_H_ */
