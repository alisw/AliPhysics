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
  bool GetDokTandMultBinning() {
    return fDokTandMultBinning;
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
  bool GetDomTMultPlots() {
    return fmTMultPlots;
  }
  bool GetDoAncestorsPlots() {
    return fAncestors;
  }
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

  void FillSameEventkTandMultDist(int i, float kT, float RelK, int multBin) {
    fSameEventkTandMultDist[i][multBin]->Fill(RelK, kT);
  }
  void FillMixedEventkTandMultDist(int i, float kT, float RelK, int multBin) {
    fMixedEventkTandMultDist[i][multBin]->Fill(RelK, kT);
  }
  void FillSameEventmTDist(int i, float mT, float RelK) {
    if (fSameEventmTDist[i])
      fSameEventmTDist[i]->Fill(RelK, mT);
  }
  void FillSameEventmTMultDist(int i, float mT, int iMult, float RelK);
  void FillMixedEventDist(int i, float RelK) {
    fMixedEventDist[i]->Fill(RelK);
  }
  ;
  void FillPtQADist(int i, float kstar, float pt1, float pt2) {
    if (!fMinimalBooking) {
      // TODO for the moment the threshold is hardcoded to 200 MeV/c
      if (fPtQADist[i] && kstar < 0.2) {
        fPtQADist[i]->Fill(pt1, pt2);
      }
    }
  }
  void FillPtSEOneQADist(int i, float pt, int mult) {
    if (!fMinimalBooking) {
      fPtQADistSEPartOne[i]->Fill(pt, mult);
    }
  }
  void FillPtSETwoQADist(int i, float pt, int mult) {
    if (!fMinimalBooking) {
      fPtQADistSEPartTwo[i]->Fill(pt, mult);
    }
  }
  void FillPtMEOneQADist(int i, float pt, int mult) {
    if (!fMinimalBooking) {
      fPtQADistMEPartOne[i]->Fill(pt, mult);
    }
  }
  void FillPtMETwoQADist(int i, float pt, int mult) {
    if (!fMinimalBooking) {
      fPtQADistMEPartTwo[i]->Fill(pt, mult);
    }
  }
  void FillMassQADist(int i, float kstar, float invMass1, float invMass2) {
    if (!fMinimalBooking) {
      if (fMassQADistPart1[i] && fMassQADistPart2[i]) {
        fMassQADistPart1[i]->Fill(invMass1, kstar);
        fMassQADistPart2[i]->Fill(invMass2, kstar);
      }
    }
  }

  void FillPairInvMassQAD(int i, AliFemtoDreamBasePart &part1,
                          AliFemtoDreamBasePart &part2) {
    if (!fMinimalBooking) {
      if (fPairInvMassQAD[i]) {
        TVector3 momPart1 = part1.GetMomentum();
        TVector3 momPart2 = part2.GetMomentum();
        TLorentzVector trackPos, trackNeg;
        trackPos.SetXYZM(momPart1.Px(), momPart1.Py(), momPart1.Pz(),
                         part1.GetInvMass());
        trackNeg.SetXYZM(momPart2.Px(), momPart2.Py(), momPart2.Pz(),
                         part2.GetInvMass());
        TLorentzVector trackSum = trackPos + trackNeg;

        fPairInvMassQAD[i]->Fill(trackSum.M());
      }
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
  void FillMixedEventmTMultDist(int i, float mT, int iMult, float RelK);
  void FillPartnersSE(int hist, int nPart1, int nPart2) {
    if (!fMinimalBooking)
      fPairCounterSE[hist]->Fill(nPart1, nPart2);
  }
  void FillPartnersME(int hist, int nPart1, int nPart2) {
    if (!fMinimalBooking)
      fPairCounterME[hist]->Fill(nPart1, nPart2);
  }
  void FillMomentumResolutionSE(int hist, float RelKTrue, float RelKReco) {
    if (!fMinimalBooking)
      fMomResolutionSE[hist]->Fill(RelKTrue, RelKReco);
  }
  void FillMomentumResolutionSEAll(int hist, float RelKTrue, float RelKReco) {
    if (!fMinimalBooking)
      fMomResolutionSEAll[hist]->Fill(RelKTrue, RelKReco);
  }
  void FillMomentumResolutionME(int hist, float RelKTrue, float RelKReco) {
    if (!fMinimalBooking)
      fMomResolutionME[hist]->Fill(RelKTrue, RelKReco);
    if (!fMinimalBooking)
      fMomResolutionDist[hist]->Fill(RelKReco - RelKTrue, RelKTrue);
  }
  void FillMomentumResolutionMEAll(int hist, float RelKTrue, float RelKReco) {
    if (!fMinimalBooking)
      fMomResolutionMEAll[hist]->Fill(RelKTrue, RelKReco);
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
  void FillEtaPhiAverageSE(int hist, int iDaug, float dPhi, float dEta,
                           bool BeforeOrAfter) {
    if (!fMinimalBooking && fPhiEtaPlots) {
      if (BeforeOrAfter) {
        fIntRadiiQAEtaPhiSEBefore[hist][iDaug]->Fill(dEta, dPhi);
      } else {
        fIntRadiiQAEtaPhiSEAfter[hist][iDaug]->Fill(dEta, dPhi);
      }
    }
  }
  void FillEtaPhiAverageME(int hist, int iDaug, float dPhi, float dEta,
                           bool BeforeOrAfter) {
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

  void FillSameEventDistCommon(int i, float RelK) {
    if (!fMinimalBooking) {
      if (fAncestors)
	fSameEventDistCommon[i]->Fill(RelK);
    }
  }
  ;
  void FillSameEventDistNonCommon(int i, float RelK) {
    if (!fMinimalBooking) {
      if (fAncestors)
	fSameEventDistNonCommon[i]->Fill(RelK);
    }
  }
  ;
  void FillSameEventMultDistCommon(int i, int iMult, float RelK) {
    if (!fMinimalBooking) {
      if (fSameEventMultDistCommon[i])
	fSameEventMultDistCommon[i]->Fill(RelK, iMult);
    }
  }
  void FillSameEventMultDistNonCommon(int i, int iMult, float RelK) {
    if (!fMinimalBooking) {
      if (fSameEventMultDistNonCommon[i])
	fSameEventMultDistNonCommon[i]->Fill(RelK, iMult);
    }
  }
  void FillSameEventmTDistCommon(int i, float mT, float RelK) {
    if (fSameEventmTDistCommon[i])
      fSameEventmTDistCommon[i]->Fill(RelK, mT);
  }
  void FillSameEventmTDistNonCommon(int i, float mT, float RelK) {
    if (fSameEventmTDistNonCommon[i])
      fSameEventmTDistNonCommon[i]->Fill(RelK, mT);
  }
  void FilldPhidEtaSECommon(int iHist, float dPhi, float dEta, float mT);
  void FilldPhidEtaSENonCommon(int iHist, float dPhi, float dEta, float mT);


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
  TH2F **fSameEventmTvsMultDist;
  TH2F **fSameEventkTDist;
  TH2F ***fSameEventkTandMultDist;
  TH2F ***fSameEventkTCentDist;
  TH2F ***fSameEventmTMultDist;
  TH2F **fPtQADist;
  TH2F **fPtQADistSEPartOne;
  TH2F **fPtQADistSEPartTwo;
  TH2F **fPtQADistMEPartOne;
  TH2F **fPtQADistMEPartTwo;
  TH2F **fMassQADistPart1;
  TH2F **fMassQADistPart2;
  TH1F **fPairInvMassQAD;
  TH2F **fPairCounterSE;
  TH1F **fMixedEventDist;
  TH2F **fMixedEventMultDist;
  TH2F **fMixedEventCentDist;
  TH2F **fMixedEventmTDist;
  TH2F **fMixedEventmTvsMultDist;
  TH2F **fMixedEventkTDist;
  TH2F ***fMixedEventkTandMultDist;
  TH2F ***fMixedEventkTCentDist;
  TH2F ***fMixedEventmTMultDist;
  TH2F **fPairCounterME;
  TH2F **fMomResolutionSE;
  TH2F **fMomResolutionSEAll;
  TH2F **fMomResolutionME;
  TH2F **fMomResolutionMEAll;
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
  TH1F **fSameEventDistCommon;
  TH1F **fSameEventDistNonCommon;   
  TH2F **fdEtadPhiSECommon;
  TH2F **fdEtadPhiSENonCommon;
  TH2F **fSameEventMultDistCommon;
  TH2F **fSameEventMultDistNonCommon;
  TH2F **fSameEventmTDistCommon;
  TH2F **fSameEventmTDistNonCommon;




  bool fDoMultBinning;
  bool fDoCentBinning;
  bool fDokTandMultBinning;
  bool fDokTBinning;
  bool fDomTBinning;
  bool fmTMultPlots;
  bool fPtQA;
  bool fMassQA;
  bool fDokTCentralityBins;
  bool fdPhidEtaPlots;
  bool fPhiEtaPlotsSmallK;
  bool fmTDetaDPhi;
  bool fAncestors;
  std::vector<int> fPDGCode;
  std::vector<float> fmTBins;
  std::vector<unsigned int> fWhichPairs;
  std::vector<int> fCentBins;
  ClassDef(AliFemtoDreamCorrHists,10);
};

#endif /* ALIFEMTODREAMCORRHISTS_H_ */
