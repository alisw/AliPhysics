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
  AliFemtoDreamCorrHists(AliFemtoDreamCollConfig *conf,bool MinimalBooking);
  virtual ~AliFemtoDreamCorrHists();
  bool GetDoMultBinning(){return fDoMultBinning;};
  void FillSameEventDist(int i,double RelK){fSameEventDist[i]->Fill(RelK);};
  void FillSameEventMultDist(int i,int iMult,double RelK){
    if (fSameEventMultDist[i])fSameEventMultDist[i]->Fill(RelK,iMult);
  }
  void FillMixedEventDist(int i,double RelK){fMixedEventDist[i]->Fill(RelK);};
  void FillMixedEventMultDist(int i,int iMult,double RelK){
     if (fMixedEventMultDist[i])fMixedEventMultDist[i]->Fill(RelK,iMult);
   }
  void FillPartnersSE(int hist,int nPart1,int nPart2){
    if (!fMinimalBooking)fPairCounterSE[hist]->Fill(nPart1,nPart2);
  }
  void FillPartnersME(int hist,int nPart1,int nPart2){
    if (!fMinimalBooking)fPairCounterME[hist]->Fill(nPart1,nPart2);
  }
  TList* GetHistList(){return fResults;};
  TList* GetQAHists(){return fQA;};
 private:
  TList         *fQA;
  TList         *fResults;
  TList         **fPairs;
  TList         **fPairQA;
  bool          fMinimalBooking;
  TH1F          **fSameEventDist;
  TH2F          **fSameEventMultDist;
  TH2F          **fPairCounterSE;
  TH1F          **fMixedEventDist;
  TH2F          **fMixedEventMultDist;
  TH2F          **fPairCounterME;
  bool          fDoMultBinning;

  ClassDef(AliFemtoDreamCorrHists,1);
};

#endif /* ALIFEMTODREAMCORRHISTS_H_ */
