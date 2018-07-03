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
  bool GetDoCentBinning(){return fDoCentBinning;};
  bool GetDokTBinning(){return fDokTBinning;};
  bool GetDomTBinning(){return fDomTBinning;};
  bool GetObtainMomentumResolution() {return fMomentumResolution;};
  bool GetEtaPhiPlots() {return fPhiEtaPlots;};
  bool GetDoMCCommonAncest() {return fDoMCCommonAncest;};
  bool GetDodPhidEtaPlots() {return fdPhidEtaPlots;};
  void FillSameEventDist(int i,float RelK){fSameEventDist[i]->Fill(RelK);};
  void FillSameEventCommonAncestDist(int i,float RelK){
    fSameEventCommonAncestDist[i]->Fill(RelK);};
  void FillSameEventNonCommonAncestDist(int i,float RelK){
    fSameEventNonCommonAncestDist[i]->Fill(RelK);};
  void FillSameEventMultDist(int i,int iMult,float RelK){
    if (fSameEventMultDist[i])fSameEventMultDist[i]->Fill(RelK,iMult);
  }
  void FillSameEventCentDist(int i,int iCent,float RelK){
    if (fSameEventCentDist[i])fSameEventCentDist[i]->Fill(RelK,iCent);
  }
  void FillSameEventkTDist(int i,float kT,float RelK,float cent){
    if (fSameEventkTDist[i])fSameEventkTDist[i]->Fill(RelK,kT);
    if (fDokTCentralityBins) FillSameEventkTCentDist(i,kT,RelK,cent);
  }
  void FillSameEventkTCentDist(int i,float kT,float RelK,float cent);
  void FillSameEventmTDist(int i,float mT, float RelK){
    if (fSameEventmTDist[i])fSameEventmTDist[i]->Fill(RelK,mT);
  }
  void FillMixedEventDist(int i,float RelK){fMixedEventDist[i]->Fill(RelK);};
  void FillMixedEventMultDist(int i,int iMult,float RelK){
    if (fMixedEventMultDist[i])fMixedEventMultDist[i]->Fill(RelK,iMult);
  }
  void FillMixedEventCentDist(int i,int iCent,float RelK){
    if (fMixedEventCentDist[i])fMixedEventCentDist[i]->Fill(RelK,iCent);
  }
  void FillMixedEventkTDist(int i,float kT,float RelK,float cent){
    if (fMixedEventkTDist[i])fMixedEventkTDist[i]->Fill(RelK,kT);
    if (fDokTCentralityBins) FillMixedEventkTCentDist(i,kT,RelK,cent);
  }
  void FillMixedEventkTCentDist(int i,float kT,float RelK,float cent);

  void FillMixedEventmTDist(int i,float mT,float RelK){
    if (fMixedEventmTDist[i])fMixedEventmTDist[i]->Fill(RelK,mT);
  }
  void FillPartnersSE(int hist,int nPart1,int nPart2){
    if (!fMinimalBooking)fPairCounterSE[hist]->Fill(nPart1,nPart2);
  }
  void FillPartnersME(int hist,int nPart1,int nPart2){
    if (!fMinimalBooking)fPairCounterME[hist]->Fill(nPart1,nPart2);
  }
  void FillMomentumResolution(int hist,float RelKTrue,float RelKReco) {
    if (!fMinimalBooking)fMomResolution[hist]->Fill(RelKTrue,RelKReco);
    if (!fMinimalBooking)fMomResolutionDist[hist]->Fill(RelKReco-RelKTrue,RelKTrue);
  }
  void FillEtaPhiAtRadiiSE(int hist,int iDaug,int iRad,float dPhi,float dEta){
    if (!fMinimalBooking&&fPhiEtaPlots)fRadiiEtaPhiSE[hist][iDaug][iRad]->Fill(dEta,dPhi);
  }
  void FillEtaPhiAtRadiiME(int hist,int iDaug,int iRad,float dPhi,float dEta){
    if (!fMinimalBooking&&fPhiEtaPlots)fRadiiEtaPhiME[hist][iDaug][iRad]->Fill(dEta,dPhi);
  }
  void FilldPhidEtaSE(int iHist,float dPhi,float dEta) {
    if (fdPhidEtaPlots) fdEtadPhiSE[iHist]->Fill(dEta,dPhi);
  }
  void FilldPhidEtaME(int iHist,float dPhi,float dEta) {
    if (fdPhidEtaPlots) fdEtadPhiME[iHist]->Fill(dEta,dPhi);
  }
  void FillEffectiveMixingDepth(int iHist,int iDepth) {
    if (!fMinimalBooking)fEffMixingDepth[iHist]->Fill(iDepth);
  }
  TList* GetHistList(){return fResults;};
  TList* GetQAHists(){return fQA;};
  TString ClassName(){return "AliFemtoDreamCorrHists";}
 private:
  TList         *fQA;
  TList         *fResults;
  TList         **fPairs;
  TList         **fPairQA;
  bool          fMinimalBooking;
  bool          fMomentumResolution;
  bool          fPhiEtaPlots;
  TH1F          **fSameEventDist;
  TH1F          **fSameEventCommonAncestDist;
  TH1F          **fSameEventNonCommonAncestDist;
  TH2F          **fSameEventMultDist;
  TH2F          **fSameEventCentDist;
  TH2F          **fSameEventmTDist;
  TH2F          **fSameEventkTDist;
  TH2F          ***fSameEventkTCentDist;
  TH2F          **fPairCounterSE;
  TH1F          **fMixedEventDist;
  TH2F          **fMixedEventMultDist;
  TH2F          **fMixedEventCentDist;
  TH2F          **fMixedEventmTDist;
  TH2F          **fMixedEventkTDist;
  TH2F          ***fMixedEventkTCentDist;
  TH2F          **fPairCounterME;
  TH2F          **fMomResolution;
  TH2F          **fMomResolutionDist;
  TH2F          ****fRadiiEtaPhiSE;
  TH2F          ****fRadiiEtaPhiME;
  TH2F          **fdEtadPhiSE;
  TH2F          **fdEtadPhiME;
  TH1F          **fEffMixingDepth;
  bool          fDoMultBinning;
  bool          fDoCentBinning;
  bool          fDokTBinning;
  bool          fDomTBinning;
  bool          fDokTCentralityBins;
  bool          fDoMCCommonAncest;
  bool          fdPhidEtaPlots;
  std::vector<float> fCentBins;
  ClassDef(AliFemtoDreamCorrHists,3);
};

#endif /* ALIFEMTODREAMCORRHISTS_H_ */
