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
#include "AliFemtoDreamBasePart.h"

class AliFemtoDreamCorrHists {
 public:
  AliFemtoDreamCorrHists();
  AliFemtoDreamCorrHists(const AliFemtoDreamCorrHists& hists);
  AliFemtoDreamCorrHists(AliFemtoDreamCollConfig *conf, bool MinimalBooking);
  AliFemtoDreamCorrHists &operator=(const AliFemtoDreamCorrHists& hists);
  virtual ~AliFemtoDreamCorrHists();
  bool GetDoMultBinning() {
    return fDoMultBinning;
  }
  ;
  bool GetDoCentBinning() {
    return fDoCentBinning;
  }
  ;
  bool GetDokTBinning() {
    return fDokTBinning;
  }
  ;
  bool GetDomTBinning() {
    return fDomTBinning;
  }
  ;
  bool GetDoPtQA() const {
    return fPtQA;
  }
  bool GetDoMassQA() const {
    return fMassQA;
  }
  bool GetObtainMomentumResolution() {
    return fMomentumResolution;
  }
  ;
  bool GetEtaPhiPlots() {
    return fPhiEtaPlots;
  }
  ;
  bool GetDodPhidEtaPlots() {
    return fdPhidEtaPlots;
  }
  ;
  bool GetDodPhidEtaPlotsSmallK() {
    return fPhiEtaPlotsSmallK;
  }
  bool GetDodPhidEtamTPlots() {
    return fmTDetaDPhi;
  }
  ;
  void FillSameEventDist(int i, float RelK) {
    fSameEventDist[i]->Fill(RelK);
  }
  ;
  void FillSameEventMultDist(int i, int iMult, float RelK) {
    if (fSameEventMultDist[i])
      fSameEventMultDist[i]->Fill(RelK, iMult);
  }
  void FillSameEventCentDist(int i, float iCent, float RelK) {
    if (fSameEventCentDist[i])
      fSameEventCentDist[i]->Fill(RelK, iCent);
  }
  void FillSameEventkTDist(int i, float kT, float RelK, float cent) {
    if (fSameEventkTDist[i])
      fSameEventkTDist[i]->Fill(RelK, kT);
    if (fDokTCentralityBins)
      FillSameEventkTCentDist(i, kT, RelK, cent);
  }
  void FillSameEventkTCentDist(int i, float kT, float RelK, float cent);
  void FillSameEventmTDist(int i, float mT, float RelK) {
    if (fSameEventmTDist[i])
      fSameEventmTDist[i]->Fill(RelK, mT);
  }
  void FillMixedEventDist(int i, float RelK) {
    fMixedEventDist[i]->Fill(RelK);
  }
  ;
  void FillPtQADist(int i, float kstar, float pt1, float pt2) {
    // TODO for the moment the threshold is hardcoded to 200 MeV/c
    if(fPtQADist[i] && kstar < 0.2) {
      fPtQADist[i]->Fill(pt1, pt2);
    }
  }
  void FillMassQADist(int i, float kstar, float invMass1, float invMass2) {
    if(fMassQADistPart1[i] && fMassQADistPart2[i]) {
      fMassQADistPart1[i]->Fill(invMass1, kstar);
      fMassQADistPart2[i]->Fill(invMass2, kstar);      
    }

}

  void FillPairInvMassQAD(int i, AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2) {
    if(fPairInvMassQAD[i]) {
 	TVector3 momPart1 = part1.GetMomentum();
        TVector3 momPart2 = part2.GetMomentum();
        TLorentzVector trackPos, trackNeg;
          trackPos.SetXYZM(momPart1.Px(), momPart1.Py(), momPart1.Pz(), part1.GetInvMass());
          trackNeg.SetXYZM(momPart2.Px(), momPart2.Py(), momPart2.Pz(), part2.GetInvMass());
         TLorentzVector trackSum = trackPos + trackNeg;
    
       fPairInvMassQAD[i]->Fill(trackSum.M());
       }
 }
  void FillMixedEventMultDist(int i, int iMult, float RelK) {
    if (fMixedEventMultDist[i])
      fMixedEventMultDist[i]->Fill(RelK, iMult);
  }
  void FillMixedEventCentDist(int i, float iCent, float RelK) {
    if (fMixedEventCentDist[i])
      fMixedEventCentDist[i]->Fill(RelK, iCent);
  }
  void FillMixedEventkTDist(int i, float kT, float RelK, float cent) {
    if (fMixedEventkTDist[i])
      fMixedEventkTDist[i]->Fill(RelK, kT);
    if (fDokTCentralityBins)
      FillMixedEventkTCentDist(i, kT, RelK, cent);
  }
  void FillMixedEventkTCentDist(int i, float kT, float RelK, float cent);

  void FillMixedEventmTDist(int i, float mT, float RelK) {
    if (fMixedEventmTDist[i])
      fMixedEventmTDist[i]->Fill(RelK, mT);
  }
  void FillPartnersSE(int hist, int nPart1, int nPart2) {
    if (!fMinimalBooking)
      fPairCounterSE[hist]->Fill(nPart1, nPart2);
  }
  void FillPartnersME(int hist, int nPart1, int nPart2) {
    if (!fMinimalBooking)
      fPairCounterME[hist]->Fill(nPart1, nPart2);
  }
  void FillMomentumResolution(int hist, float RelKTrue, float RelKReco) {
    if (!fMinimalBooking)
      fMomResolution[hist]->Fill(RelKTrue, RelKReco);
    if (!fMinimalBooking)
      fMomResolutionDist[hist]->Fill(RelKReco - RelKTrue, RelKTrue);
  }
  void FillEtaPhiAtRadiiSE(int hist, int iDaug, int iRad, float dPhi,
                           float dEta, float relk) {
    if (!fMinimalBooking && fPhiEtaPlots) {
      fRadiiEtaPhiSE[hist][iDaug][iRad]->Fill(dEta, dPhi);
      if (relk < fRelKThreshold && fPhiEtaPlotsSmallK) {
        fRadiiEtaPhiSEsmallK[hist][iDaug][iRad]->Fill(dEta, dPhi);
      }
    }
  }
  void FillEtaPhiAtRadiiME(int hist, int iDaug, int iRad, float dPhi,
                           float dEta, float relk) {
    if (!fMinimalBooking && fPhiEtaPlots) {
      fRadiiEtaPhiME[hist][iDaug][iRad]->Fill(dEta, dPhi);
      if (relk < fRelKThreshold && fPhiEtaPlotsSmallK) {
        fRadiiEtaPhiMEsmallK[hist][iDaug][iRad]->Fill(dEta, dPhi);
      }
    }
  }
  void FillEtaPhiAverageSE(int hist, int iDaug, float dPhi, float dEta, bool BeforeOrAfter) {
    if (!fMinimalBooking && fPhiEtaPlots) {
      if (BeforeOrAfter) {
        fIntRadiiQAEtaPhiSEBefore[hist][iDaug]->Fill(dEta, dPhi);
      } else {
        fIntRadiiQAEtaPhiSEAfter[hist][iDaug]->Fill(dEta, dPhi);
      }
    }
  }
  void FillEtaPhiAverageME(int hist, int iDaug, float dPhi, float dEta, bool BeforeOrAfter) {
    if (!fMinimalBooking && fPhiEtaPlots) {
      if (BeforeOrAfter) {
        fIntRadiiQAEtaPhiMEBefore[hist][iDaug]->Fill(dEta, dPhi);
      } else {
        fIntRadiiQAEtaPhiMEAfter[hist][iDaug]->Fill(dEta, dPhi);
      }
    }
  }
  void FilldPhidEtaSE(int iHist, float dPhi, float dEta, float mT);
  void FilldPhidEtaME(int iHist, float dPhi, float dEta, float mT);
  void FillEffectiveMixingDepth(int iHist, int iDepth) {
    if (!fMinimalBooking)
      fEffMixingDepth[iHist]->Fill(iDepth);
  }
  TList* GetHistList() {
    return fResults;
  }
  ;
  TList* GetQAHists() {
    return fQA;
  }
  ;
  TString ClassName() {
    return "AliFemtoDreamCorrHists";
  }
 private:
  TList *fQA;
  TList *fResults;
  TList **fPairs;
  TList **fPairQA;
  bool fMinimalBooking;
  bool fMomentumResolution;
  bool fPhiEtaPlots;
  float fRelKThreshold;
  TH1F **fSameEventDist;
  TH2F **fSameEventMultDist;
  TH2F **fSameEventCentDist;
  TH2F **fSameEventmTDist;
  TH2F **fSameEventkTDist;
  TH2F ***fSameEventkTCentDist;
  TH2F **fPtQADist;
  TH2F **fMassQADistPart1;
  TH2F **fMassQADistPart2;
  TH1F **fPairInvMassQAD;
  TH2F **fPairCounterSE;
  TH1F **fMixedEventDist;
  TH2F **fMixedEventMultDist;
  TH2F **fMixedEventCentDist;
  TH2F **fMixedEventmTDist;
  TH2F **fMixedEventkTDist;
  TH2F ***fMixedEventkTCentDist;
  TH2F **fPairCounterME;
  TH2F **fMomResolution;
  TH2F **fMomResolutionDist;
  TH2F ****fRadiiEtaPhiSE;
  TH2F ****fRadiiEtaPhiME;
  TH2F ***fIntRadiiQAEtaPhiSEBefore;
  TH2F ***fIntRadiiQAEtaPhiMEBefore;
  TH2F ***fIntRadiiQAEtaPhiSEAfter;
  TH2F ***fIntRadiiQAEtaPhiMEAfter;
  TH2F ****fRadiiEtaPhiSEsmallK;
  TH2F ****fRadiiEtaPhiMEsmallK;
  TH2F **fdEtadPhiSE;
  TH2F **fdEtadPhiME;
  TH2F ***fdEtadPhiSEmT;
  TH2F ***fdEtadPhiMEmT;
  TH1F **fEffMixingDepth;
  bool fDoMultBinning;
  bool fDoCentBinning;
  bool fDokTBinning;
  bool fDomTBinning;
  bool fPtQA;
  bool fMassQA;
  bool fDokTCentralityBins;
  bool fdPhidEtaPlots;
  bool fPhiEtaPlotsSmallK;
  bool fmTDetaDPhi;
  std::vector<int> fPDGCode;
  std::vector<float> fmTdEtadPhiBins;
  std::vector<unsigned int> fWhichPairs;
  std::vector<int> fCentBins;
  ClassDef(AliFemtoDreamCorrHists,8);
};

#endif /* ALIFEMTODREAMCORRHISTS_H_ */
