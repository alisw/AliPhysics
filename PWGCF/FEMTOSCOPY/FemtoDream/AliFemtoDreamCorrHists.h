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
  bool GetObtainMomentumResolution() {return fMomentumResolution;};
  bool GetEtaPhiPlots() {return fPhiEtaPlots;};
  void FillSameEventDist(int i,float RelK){fSameEventDist[i]->Fill(RelK);};
  void FillSameEventMultDist(int i,int iMult,float RelK){
    if (fSameEventMultDist[i])fSameEventMultDist[i]->Fill(RelK,iMult);
  }
  void FillMixedEventDist(int i,float RelK){fMixedEventDist[i]->Fill(RelK);};
  void FillMixedEventMultDist(int i,int iMult,float RelK){
     if (fMixedEventMultDist[i])fMixedEventMultDist[i]->Fill(RelK,iMult);
   }
  void FillPartnersSE(int hist,int nPart1,int nPart2){
    if (!fMinimalBooking)fPairCounterSE[hist]->Fill(nPart1,nPart2);
  }
  void FillPartnersME(int hist,int nPart1,int nPart2){
    if (!fMinimalBooking)fPairCounterME[hist]->Fill(nPart1,nPart2);
  }
  void FillMomentumResolution(int hist,float RelKTrue,float RelKReco) {
    if (!fMinimalBooking)fMomResolution[hist]->Fill(RelKTrue,RelKReco);
  }
  void FillEtaPhiAtRadiiSE(int hist,int iDaug,int iRad,float dPhi,float dEta){
    if (!fMinimalBooking&&fPhiEtaPlots)fRadiiEtaPhiSE[hist][iDaug][iRad]->Fill(dEta,dPhi);
  }
  void FillEtaPhiAtRadiiME(int hist,int iDaug,int iRad,float dPhi,float dEta){
    if (!fMinimalBooking&&fPhiEtaPlots)fRadiiEtaPhiME[hist][iDaug][iRad]->Fill(dEta,dPhi);
  }
  TList* GetHistList(){return fResults;};
  TList* GetQAHists(){return fQA;};
 private:
  TList         *fQA;
  TList         *fResults;
  TList         **fPairs;
  TList         **fPairQA;
  bool          fMinimalBooking;
  bool          fMomentumResolution;
  bool          fPhiEtaPlots;
  TH1F          **fSameEventDist;
  TH2F          **fSameEventMultDist;
  TH2F          **fPairCounterSE;
  TH1F          **fMixedEventDist;
  TH2F          **fMixedEventMultDist;
  TH2F          **fPairCounterME;
  TH2F          **fMomResolution;
  TH2F          ****fRadiiEtaPhiSE;
  TH2F          ****fRadiiEtaPhiME;
  bool          fDoMultBinning;

  ClassDef(AliFemtoDreamCorrHists,2);
};

#endif /* ALIFEMTODREAMCORRHISTS_H_ */
