/*
 * AliFemtoDreamv0Cuts.cxx
 *
 *  Created on: Dec 12, 2017
 *      Author: gu74req
 */

#include "AliFemtoDreamv0Cuts.h"
#include "TDatabasePDG.h"
#include "AliLog.h"
ClassImp(AliFemtoDreamv0Cuts)
AliFemtoDreamv0Cuts::AliFemtoDreamv0Cuts()
    : fHistList(),
      fMCHistList(),
      fMCHist(nullptr),
      fHist(nullptr),
      fPosCuts(0),
      fNegCuts(0),
      fCutDaughters(false),
      fMinimalBooking(false),
      fMCData(false),
      fCPAPlots(false),
      fContribSplitting(false),
      fDoMultBinning(false),
      fCheckMother(false),
      fRunNumberQA(false),
      fMinRunNumber(0),
      fMaxRunNumber(0),
      fCutOnFlyStatus(false),
      fOnFlyStatus(false),
      fCutCharge(false),
      fCharge(0),
      fDoCombinedTimingCut(false),
      fCombinedTiming(BothDaughtersCombined),
      fDoArmenterosCut(false),
      fArmenterosQtLow(0),
      fArmenterosQtUp(1),
      fArmenterosAlphaLow(-1),
      fArmenterosAlphaUp(1),
      fCutPt(false),
      fpTmin(0),
      fpTmax(0),
      fKaonRejection(false),
      fKaonRejLow(0),
      fKaonRejUp(0),
      fCutDecayVtxXYZ(false),
      fMaxDecayVtxXYZ(0),
      fCutTransRadius(false),
      fMinTransRadius(0),
      fMaxTransRadius(0),
      fCutMinDCADaugPrimVtx(false),
      fMinDCADaugToPrimVtx(0),
      fCutMaxDCADaugToDecayVtx(false),
      fMaxDCADaugToDecayVtx(0),
      fCutCPA(false),
      fMinCPA(0),
      fCutInvMass(false),
      fInvMassCutWidth(0),
      fCutInvMassSidebands(false),
      fInvMassCutSBdown(0),
      fInvMassCutSBup(0),
      fAxisMinMass(0),
      fAxisMaxMass(1),
      fNumberXBins(1),
      fPDGv0(0),
      fPDGDaugP(0),
      fPDGDaugN(0) {
}

AliFemtoDreamv0Cuts::AliFemtoDreamv0Cuts(const AliFemtoDreamv0Cuts& cuts)
    : fHistList(cuts.fHistList),
      fMCHistList(cuts.fMCHistList),
      fMCHist(cuts.fMCHist),
      fHist(cuts.fHist),
      fPosCuts(cuts.fPosCuts),
      fNegCuts(cuts.fNegCuts),
      fCutDaughters(cuts.fCutDaughters),
      fMinimalBooking(cuts.fMinimalBooking),
      fMCData(cuts.fMCData),
      fCPAPlots(cuts.fCPAPlots),
      fContribSplitting(cuts.fContribSplitting),
      fDoMultBinning(cuts.fDoMultBinning),
      fCheckMother(cuts.fCheckMother),
      fRunNumberQA(cuts.fRunNumberQA),
      fMinRunNumber(cuts.fMinRunNumber),
      fMaxRunNumber(cuts.fMaxRunNumber),
      fCutOnFlyStatus(cuts.fCutOnFlyStatus),
      fOnFlyStatus(cuts.fOnFlyStatus),
      fCutCharge(cuts.fCutCharge),
      fCharge(cuts.fCharge),
      fDoCombinedTimingCut(cuts.fDoCombinedTimingCut),
      fCombinedTiming(cuts.fCombinedTiming),
      fDoArmenterosCut(cuts.fDoArmenterosCut),
      fArmenterosQtLow(cuts.fArmenterosQtLow),
      fArmenterosQtUp(cuts.fArmenterosQtUp),
      fArmenterosAlphaLow(cuts.fArmenterosAlphaLow),
      fArmenterosAlphaUp(cuts.fArmenterosAlphaUp),
      fCutPt(cuts.fCutPt),
      fpTmin(cuts.fpTmin),
      fpTmax(cuts.fpTmax),
      fKaonRejection(cuts.fKaonRejection),
      fKaonRejLow(cuts.fKaonRejLow),
      fKaonRejUp(cuts.fKaonRejUp),
      fCutDecayVtxXYZ(cuts.fCutDecayVtxXYZ),
      fMaxDecayVtxXYZ(cuts.fMaxDecayVtxXYZ),
      fCutTransRadius(cuts.fCutTransRadius),
      fMinTransRadius(cuts.fMinTransRadius),
      fMaxTransRadius(cuts.fMaxTransRadius),
      fCutMinDCADaugPrimVtx(cuts.fCutMinDCADaugPrimVtx),
      fMinDCADaugToPrimVtx(cuts.fMinDCADaugToPrimVtx),
      fCutMaxDCADaugToDecayVtx(cuts.fCutMaxDCADaugToDecayVtx),
      fMaxDCADaugToDecayVtx(cuts.fMaxDCADaugToDecayVtx),
      fCutCPA(cuts.fCutCPA),
      fMinCPA(cuts.fMinCPA),
      fCutInvMass(cuts.fCutInvMass),
      fInvMassCutWidth(cuts.fInvMassCutWidth),
      fCutInvMassSidebands(cuts.fCutInvMassSidebands),
      fInvMassCutSBdown(cuts.fInvMassCutSBdown),
      fInvMassCutSBup(cuts.fInvMassCutSBup),
      fAxisMinMass(cuts.fAxisMinMass),
      fAxisMaxMass(cuts.fAxisMaxMass),
      fNumberXBins(cuts.fNumberXBins),
      fPDGv0(cuts.fPDGv0),
      fPDGDaugP(cuts.fPDGDaugP),
      fPDGDaugN(cuts.fPDGDaugN) {
}

AliFemtoDreamv0Cuts& AliFemtoDreamv0Cuts::operator=(
    const AliFemtoDreamv0Cuts& cuts) {
  if (this != &cuts) {
    this->fHistList = cuts.fHistList;
    this->fMCHistList = cuts.fMCHistList;
    this->fMCHist = cuts.fMCHist;
    this->fHist = cuts.fHist;
    this->fPosCuts = cuts.fPosCuts;
    this->fNegCuts = cuts.fNegCuts;
    this->fCutDaughters = cuts.fCutDaughters;
    this->fMinimalBooking = cuts.fMinimalBooking;
    this->fMCData = cuts.fMCData;
    this->fCPAPlots = cuts.fCPAPlots;
    this->fContribSplitting = cuts.fContribSplitting;
    this->fDoMultBinning = cuts.fDoMultBinning;
    this->fCheckMother = cuts.fCheckMother;
    this->fRunNumberQA = cuts.fRunNumberQA;
    this->fMinRunNumber = cuts.fMinRunNumber;
    this->fMaxRunNumber = cuts.fMaxRunNumber;
    this->fCutOnFlyStatus = cuts.fCutOnFlyStatus;
    this->fOnFlyStatus = cuts.fOnFlyStatus;
    this->fCutCharge = cuts.fCutCharge;
    this->fCharge = cuts.fCharge;
    this->fDoCombinedTimingCut = cuts.fDoCombinedTimingCut;
    this->fCombinedTiming = cuts.fCombinedTiming;
    this->fDoArmenterosCut = cuts.fDoArmenterosCut;
    this->fArmenterosQtLow = cuts.fArmenterosQtLow;
    this->fArmenterosQtUp = cuts.fArmenterosQtUp;
    this->fArmenterosAlphaLow = cuts.fArmenterosAlphaLow;
    this->fArmenterosAlphaUp = cuts.fArmenterosAlphaUp;
    this->fCutPt = cuts.fCutPt;
    this->fpTmin = cuts.fpTmin;
    this->fpTmax = cuts.fpTmax;
    this->fKaonRejection = cuts.fKaonRejection;
    this->fKaonRejLow = cuts.fKaonRejLow;
    this->fKaonRejUp = cuts.fKaonRejUp;
    this->fCutDecayVtxXYZ = cuts.fCutDecayVtxXYZ;
    this->fMaxDecayVtxXYZ = cuts.fMaxDecayVtxXYZ;
    this->fCutTransRadius = cuts.fCutTransRadius;
    this->fMinTransRadius = cuts.fMinTransRadius;
    this->fMaxTransRadius = cuts.fMaxTransRadius;
    this->fCutMinDCADaugPrimVtx = cuts.fCutMinDCADaugPrimVtx;
    this->fMinDCADaugToPrimVtx = cuts.fMinDCADaugToPrimVtx;
    this->fCutMaxDCADaugToDecayVtx = cuts.fCutMaxDCADaugToDecayVtx;
    this->fMaxDCADaugToDecayVtx = cuts.fMaxDCADaugToDecayVtx;
    this->fCutCPA = cuts.fCutCPA;
    this->fMinCPA = cuts.fMinCPA;
    this->fCutInvMass = cuts.fCutInvMass;
    this->fInvMassCutWidth = cuts.fInvMassCutWidth;
    this->fCutInvMassSidebands = cuts.fCutInvMassSidebands;
    this->fInvMassCutSBdown = cuts.fInvMassCutSBdown;
    this->fInvMassCutSBup = cuts.fInvMassCutSBup;
    this->fAxisMinMass = cuts.fAxisMinMass;
    this->fAxisMaxMass = cuts.fAxisMaxMass;
    this->fNumberXBins = cuts.fNumberXBins;
    this->fPDGv0 = cuts.fPDGv0;
    this->fPDGDaugP = cuts.fPDGDaugP;
    this->fPDGDaugN = cuts.fPDGDaugN;
  }
  return *this;
}

AliFemtoDreamv0Cuts::~AliFemtoDreamv0Cuts() {
  if (fHist) {
    delete fHist;
  }
  //  if (fMCHist) {
  //    delete fMCHist;
  //  }
}
AliFemtoDreamv0Cuts* AliFemtoDreamv0Cuts::LambdaCuts(bool isMC, bool CPAPlots,
                                                     bool SplitContrib) {
  AliFemtoDreamv0Cuts *LambdaCuts = new AliFemtoDreamv0Cuts();
  LambdaCuts->SetIsMonteCarlo(isMC);
  LambdaCuts->SetPlotCPADist(CPAPlots);
  LambdaCuts->SetPlotContrib(SplitContrib);

  LambdaCuts->SetCheckOnFlyStatus(false);  //online = kTRUE, offline = kFALSE
  LambdaCuts->SetCutCharge(0);
  LambdaCuts->SetPtRange(0.3, 999.);
  LambdaCuts->SetKaonRejection(0.48, 0.515);
  LambdaCuts->SetCutMaxDecayVtx(100);
  LambdaCuts->SetCutTransverseRadius(0.2, 100);
  LambdaCuts->SetCutDCADaugToPrimVtx(0.05);
  LambdaCuts->SetCutDCADaugTov0Vtx(1.5);
  LambdaCuts->SetCutCPA(0.99);
  LambdaCuts->SetCutInvMass(0.004);
  LambdaCuts->SetAxisInvMassPlots(400, 1.0, 1.2);

  return LambdaCuts;
}

AliFemtoDreamv0Cuts *AliFemtoDreamv0Cuts::LambdaSigma0Cuts(bool isMC,
                                                           bool CPAPlots,
                                                           bool SplitContrib) {
  AliFemtoDreamv0Cuts *LambdaCuts = new AliFemtoDreamv0Cuts();
  LambdaCuts->SetIsMonteCarlo(isMC);
  LambdaCuts->SetPlotCPADist(CPAPlots);
  LambdaCuts->SetPlotContrib(SplitContrib);

  // Changes w.r.t the default Lambda cuts
  // 1. Kaon rejection only 1sigma of the peak
  // 2. CPA > 0.999
  // 3. Timing cut only for one daughter

  LambdaCuts->SetCheckOnFlyStatus(false);  //online = kTRUE, offline = kFALSE
  LambdaCuts->SetCutCharge(0);
  LambdaCuts->SetPtRange(0.3, 999.);
  LambdaCuts->SetKaonRejection(0.493, 0.504);
  LambdaCuts->SetCutMaxDecayVtx(100);
  LambdaCuts->SetCutTransverseRadius(0.2, 100);
  LambdaCuts->SetCutDCADaugToPrimVtx(0.05);
  LambdaCuts->SetCutDCADaugTov0Vtx(1.5);
  LambdaCuts->SetCutCPA(0.999);
  LambdaCuts->SetCutInvMass(0.006);
  LambdaCuts->SetAxisInvMassPlots(400, 1.0, 1.2);
  LambdaCuts->SetDaughterTimingCut(OneDaughterCombined);

  return LambdaCuts;
}

bool AliFemtoDreamv0Cuts::isSelected(AliFemtoDreamv0 *v0) {
  bool pass = true;
  if (!v0->IsSet()) {
    pass = false;
  } else {
    if (!fMinimalBooking)
      fHist->FillTrackCounter(0);
  }
  if (pass) {
    if (fCutDaughters) {
      if (!DaughtersPassCuts(v0)) {
        pass = false;
      }
    }
  }
  if (pass) {
    if (!MotherPassCuts(v0)) {
      pass = false;
    }
  }
  if (pass) {
    if (fDoArmenterosCut && !ArmenterosSelection(v0)) {
      pass = false;
    }
  }
  if (pass) {
    if (fKaonRejection && !RejectAsKaon(v0)) {
      pass = false;
    }
  }
  if (pass) {
    if (!CPAandMassCuts(v0)) {
      pass = false;
    }
  }
  v0->SetUse(pass);
  BookQA(v0);
  if (!fMinimalBooking) {
    if (fMCData) {
      BookMC(v0);
    }
  }
  return pass;
}

bool AliFemtoDreamv0Cuts::DaughtersPassCuts(AliFemtoDreamv0 *v0) {
  bool pass = true;
  bool passD1 = true;
  bool passD2 = true;
  if (!v0->GetHasDaughters()) {
    pass = false;
  } else {
    if (!fMinimalBooking)
      fHist->FillTrackCounter(1);
  }
  if (pass) {
    if (v0->GetCharge().at(1) < 0 && v0->GetCharge().at(2) > 0) {
      //at 1: Negative daughter, at 2: Positive Daughter should be the way it
      //was set, but sometimes it it stored in the wrong way
      if (!fMinimalBooking)
        fHist->FillTrackCounter(15);
      v0->Setv0Mass(CalculateInvMass(v0, fPDGDaugP, fPDGDaugN));
      passD1 = fNegCuts->isSelected(v0->GetNegDaughter());
      passD2 = fPosCuts->isSelected(v0->GetPosDaughter());
      if (passD1 && passD2) {
        if (!fMinimalBooking)
          fHist->FillTrackCounter(16);
      }
    } else {
      if (!fMinimalBooking)
        fHist->FillTrackCounter(17);
      v0->Setv0Mass(CalculateInvMass(v0, fPDGDaugN, fPDGDaugP));
      passD1 = fPosCuts->isSelected(v0->GetNegDaughter());
      passD2 = fNegCuts->isSelected(v0->GetPosDaughter());
      if (passD1 && passD2) {
        if (!fMinimalBooking)
          fHist->FillTrackCounter(18);
      }
    }
    if (!(passD1 && passD2)) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillTrackCounter(2);
    }
  }
  if (pass)
    pass = passD1 && passD2;
  return pass;
}
bool AliFemtoDreamv0Cuts::MotherPassCuts(AliFemtoDreamv0 *v0) {
  //all topological and kinematic cuts on the mother except for the CPA, this
  //will be filled later, in case also the CPA distributions are required.
  //Special CPA checks impelemented in the RejKaons later, to still ensure
  //proper Invariant Mass Plots.
  bool pass = true;
  if (fCutOnFlyStatus) {
    if (!(v0->GetOnlinev0() == fOnFlyStatus)) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillTrackCounter(3);
    }
  }
  if (pass && fCutCharge) {
    if (v0->GetCharge().at(0) != fCharge) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillTrackCounter(4);
    }
  }
  if (pass && fCutPt) {
    if ((v0->GetPt() < fpTmin) || (v0->GetPt() > fpTmax)) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillTrackCounter(5);
    }
  }
  if (pass && fCutDecayVtxXYZ) {
    if ((v0->GetDCAv0Vtx(0) > fMaxDecayVtxXYZ)
        || (v0->GetDCAv0Vtx(1) > fMaxDecayVtxXYZ)
        || (v0->GetDCAv0Vtx(2) > fMaxDecayVtxXYZ)) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillTrackCounter(6);
    }
  }
  if (pass && fCutTransRadius) {
    if ((v0->GetTransverseRadius() < fMinTransRadius)
        || (v0->GetTransverseRadius() > fMaxTransRadius)) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillTrackCounter(7);
    }
  }
  if (pass && fCutMinDCADaugPrimVtx) {
    if ((v0->GetDCADaugPosVtx() < fMinDCADaugToPrimVtx)
        || (v0->GetDCADaugNegVtx() < fMinDCADaugToPrimVtx)) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillTrackCounter(8);
    }
  }
  if (pass && fCutMaxDCADaugToDecayVtx) {
    if (v0->GetDaugDCA() > fMaxDCADaugToDecayVtx) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillTrackCounter(9);
    }
  }
  if (pass && fDoCombinedTimingCut) {
    if (!TimingCut(v0)) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillTrackCounter(10);
    }
  }
  return pass;
}

bool AliFemtoDreamv0Cuts::TimingCut(AliFemtoDreamv0 *v0) {
  auto *pos = v0->GetPosDaughter();
  auto *neg = v0->GetNegDaughter();
  bool posSPD = pos->GetHasSPDHit();
  bool posSSD = pos->GetITSHit(4) || pos->GetITSHit(5);
  bool posTOF = pos->GetTOFTimingReuqirement();
  bool posITS = (posSPD || posSSD);
  bool posComb = (posITS || posTOF);

  bool negSPD = neg->GetHasSPDHit();
  bool negSSD = neg->GetITSHit(4) || pos->GetITSHit(5);
  bool negTOF = neg->GetTOFTimingReuqirement();
  bool negITS = (negSPD || negSSD);
  bool negComb = (negITS || negTOF);

  switch (fCombinedTiming) {
    case (BothDaughtersCombined):
      return (posComb && negComb);
      break;
    case (OneDaughterCombined):
      return (posComb || negComb);
      break;
    case (BothDaughtersITSonly):
      return (posITS && negITS);
      break;
    case (BothDaughtersTOFonly):
      return (posTOF && negTOF);
      break;
    case (OneDaughterITSonly):
      return (posITS || negITS);
      break;
    case (OneDaughterTOFonly):
      return (posTOF || negTOF);
      break;
    case (BothDaughersSPDonly):
      return (posSPD && negSPD);
      break;
    case (OneDaugherSPDonly):
      return (posSPD || negSPD);
      break;
    case (BothDaughersSPDTOF):
      return ((posSPD || posTOF) && (negSPD || negTOF));
      break;
    case (OneDaugherSPDTOF):
      return ((posSPD || posTOF) || (negSPD || negTOF));
      break;
  }
  return false;
}

bool AliFemtoDreamv0Cuts::RejectAsKaon(AliFemtoDreamv0 *v0) {
  bool pass = true;
  bool cpaPass = true;
  if (v0->GetCPA() < fMinCPA) {
    cpaPass = false;
  }
  float massKaon = CalculateInvMass(v0, 211, 211);
  if (cpaPass) {
    if (!fMinimalBooking)
      fHist->FillInvMassBefKaonRej(v0->Getv0Mass());
    if (!fMinimalBooking)
      fHist->FillInvMassKaon(massKaon);
  }
  if (fKaonRejLow < massKaon && massKaon < fKaonRejUp) {
    pass = false;
  } else {
    if (!fMinimalBooking)
      fHist->FillTrackCounter(12);
  }
  return pass;
}

bool AliFemtoDreamv0Cuts::ArmenterosSelection(AliFemtoDreamv0 *v0) {
  bool pass = true;
  float armQt = v0->GetArmenterosQt();
  float armAlpha = v0->GetArmenterosAlpha();
  if (fDoArmenterosCut) {
    if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
      pass = false;
    }
    float prefactorAlpha = (fPDGv0 > 0) ? 1.f : -1.f;  // for anti-particles
                                                       // negative to compensate
                                                       // for sign change of qt
    if (armAlpha * prefactorAlpha > fArmenterosAlphaUp
        || armAlpha * prefactorAlpha < fArmenterosAlphaLow)
      pass = false;
  } else {
    if (!fMinimalBooking) {
      fHist->FillTrackCounter(11);
    }
  }
  return pass;
}

bool AliFemtoDreamv0Cuts::CPAandMassCuts(AliFemtoDreamv0 *v0) {
  //here we cut mass and cpa to fill the cpa distribution properly
  bool cpaPass = true;
  bool massPass = true;
  if (fCutCPA && (v0->GetCPA() < fMinCPA)) {
    cpaPass = false;
  }
  if (fCutInvMass) {
    float massv0 = TDatabasePDG::Instance()->GetParticle(fPDGv0)->Mass();
    if ((v0->Getv0Mass() < massv0 - fInvMassCutWidth)
        || (massv0 + fInvMassCutWidth < v0->Getv0Mass())) {
      massPass = false;
    }
  } else if (fCutInvMassSidebands) {
    if ((v0->Getv0Mass() < fInvMassCutSBdown)
        || (fInvMassCutSBup < v0->Getv0Mass())) {
      massPass = false;
    }
  }
  //now with this information fill the histograms
  if (cpaPass) {
    fHist->FillInvMassPtBins(v0->GetPt(), v0->Getv0Mass());
    if (!fMinimalBooking)
      fHist->Fillv0MassDist(v0->Getv0Mass());
    if (fRunNumberQA) {
      fHist->FillInvMassPerRunNumber(v0->GetEvtNumber(), v0->Getv0Mass());
    }
  }
  if (massPass && fCPAPlots && !fMinimalBooking) {
    fHist->FillCPAPtBins(v0->GetPt(), v0->GetCPA(), v0->GetEventMultiplicity());
    if (fMCData) {
      fMCHist->FillMCCPAPtBins(v0->GetParticleOrigin(), v0->GetPt(),
                               v0->GetCPA(), v0->GetEventMultiplicity());
    }
  }
  if (massPass) {
    if (!fMinimalBooking)
      fHist->FillTrackCounter(13);
  }
  if (massPass && cpaPass) {
    if (!fMinimalBooking)
      fHist->FillTrackCounter(14);
  }
  bool pass = massPass && cpaPass;
  return pass;
}
void AliFemtoDreamv0Cuts::Init() {
  //Cant be set externally cause else the lists don't exist. Needs to be changed in case
  //it is needed
  if (fCutDaughters) {
    fPosCuts->SetMinimalBooking(fMinimalBooking);
    fPosCuts->Init();
    fNegCuts->SetMinimalBooking(fMinimalBooking);
    fNegCuts->Init();
  }
  if (!fMinimalBooking) {
    fHist = new AliFemtoDreamv0Hist(fNumberXBins, fAxisMinMass, fAxisMaxMass,
                                    fCPAPlots, fRunNumberQA, fMinRunNumber,
                                    fMaxRunNumber);
    BookTrackCuts();
    fHistList = new TList();
    fHistList->SetOwner();
    fHistList->SetName("v0Cuts");
    fHistList->Add(fHist->GetHistList());
    if (fCutDaughters) {
      fPosCuts->SetName("PosCuts");
      fHistList->Add(fPosCuts->GetQAHists());
      fNegCuts->SetName("NegCuts");
      fHistList->Add(fNegCuts->GetQAHists());
    }

    if (fMCData) {
      fMCHist = new AliFemtoDreamv0MCHist(fNumberXBins, fAxisMinMass,
                                          fAxisMaxMass, fContribSplitting,
                                          fCPAPlots, fDoMultBinning,
                                          fCheckMother);
      fMCHistList = new TList();
      fMCHistList->SetOwner();
      fMCHistList->SetName("v0MCCuts");
      fMCHistList->Add(fMCHist->GetHistList());
      if (fCutDaughters) {
        fPosCuts->SetMCName("PosCuts");
        fMCHistList->Add(fPosCuts->GetMCQAHists());
        fNegCuts->SetMCName("NegCuts");
        fMCHistList->Add(fNegCuts->GetMCQAHists());
      }
    }
  } else {
    fHist = new AliFemtoDreamv0Hist("MinimalBooking", fNumberXBins,
                                    fAxisMinMass, fAxisMaxMass);
    fHistList = new TList();
    fHistList->SetOwner();
    fHistList->SetName("v0Cuts");
    fHistList->Add(fHist->GetHistList());
    if (fCutDaughters) {
      fPosCuts->SetName("PosCuts");
      fHistList->Add(fPosCuts->GetQAHists());
      fNegCuts->SetName("NegCuts");
      fHistList->Add(fNegCuts->GetQAHists());
    }
  }
}

void AliFemtoDreamv0Cuts::BookQA(AliFemtoDreamv0 *v0) {
  if (!fMinimalBooking) {
    for (int i = 0; i < 2; ++i) {
      if (i == 0 || (i == 1 && v0->UseParticle())) {
        if (!v0->GetOnlinev0()) {
          fHist->FillOnFlyStatus(i, 1);
        } else if (v0->GetOnlinev0()) {
          fHist->FillOnFlyStatus(i, 0);
        }
        fHist->FillpTCut(i, v0->GetPt());
        fHist->FillPhi(i, v0->GetPhi().at(0));
        fHist->FillEtaCut(i, v0->GetEta().at(0));
        fHist->Fillv0DecayVtxXCut(i, v0->GetDCAv0Vtx(0));
        fHist->Fillv0DecayVtxYCut(i, v0->GetDCAv0Vtx(1));
        fHist->Fillv0DecayVtxZCut(i, v0->GetDCAv0Vtx(2));
        fHist->FillTransverRadiusCut(i, v0->GetTransverseRadius());
        fHist->FillDCAPosDaugToPrimVtxCut(i, v0->GetDCADaugPosVtx());
        fHist->FillDCANegDaugToPrimVtxCut(i, v0->GetDCADaugNegVtx());
        fHist->FillDCADaugTov0VtxCut(i, v0->GetDaugDCA());
        fHist->FillCPACut(i, v0->GetCPA());
        fHist->FillArmenterosPodolandski(i, v0->GetArmenterosAlpha(),
                                         v0->GetArmenterosQt());
        fHist->FillInvMass(i, v0->Getv0Mass());
      }
    }
  }
  v0->GetPosDaughter()->SetUse(v0->UseParticle());
  v0->GetNegDaughter()->SetUse(v0->UseParticle());

  if (v0->IsSet() && fCutDaughters) {
    fPosCuts->BookQA(v0->GetPosDaughter());
    fNegCuts->BookQA(v0->GetNegDaughter());
  }
}

void AliFemtoDreamv0Cuts::BookMC(AliFemtoDreamv0 *v0) {
  if (!fMinimalBooking) {
    float pT = v0->GetPt();
    if (v0->GetHasDaughters()) {
      float etaNegDaug = v0->GetEta().at(1);
      float etaPosDaug = v0->GetEta().at(2);
      if (v0->GetMCPDGCode() == fPDGv0) {
        if (fpTmin < pT && pT < fpTmax) {
          if (fPosCuts->GetEtaMin() < etaPosDaug
              && etaPosDaug < fPosCuts->GetEtaMax()) {
            if (fNegCuts->GetEtaMin() < etaNegDaug
                && etaNegDaug < fNegCuts->GetEtaMax()) {
              fMCHist->FillMCGen(pT);
            }
          }
        }
      }
    }
    if (v0->UseParticle()) {
      fMCHist->FillMCIdent(pT);
      AliFemtoDreamBasePart::PartOrigin tmpOrg = v0->GetParticleOrigin();
      if (v0->GetMCPDGCode() == fPDGv0) {
        fMCHist->FillMCCorr(pT);
        if (tmpOrg == AliFemtoDreamBasePart::kPhysPrimary) {
          fMCHist->FillMCPtResolution(v0->GetMCPt(), v0->GetPt());
          fMCHist->FillMCThetaResolution(v0->GetMCTheta().at(0),
                                         v0->GetTheta().at(0), v0->GetMCPt());
          fMCHist->FillMCPhiResolution(v0->GetMCPhi().at(0), v0->GetPhi().at(0),
                                       v0->GetMCPt());
        }
      } else {
        v0->SetParticleOrigin(AliFemtoDreamBasePart::kContamination);
      }
      if (fContribSplitting) {
        FillMCContributions(v0);
      }
      v0->GetPosDaughter()->SetParticleOrigin(v0->GetParticleOrigin());
      v0->GetNegDaughter()->SetParticleOrigin(v0->GetParticleOrigin());
      fPosCuts->BookMC(v0->GetPosDaughter());
      fNegCuts->BookMC(v0->GetNegDaughter());
      v0->SetParticleOrigin(tmpOrg);
      if (fCheckMother)
        fMCHist->FillMCMother(v0->GetPt(), v0->GetMotherPDG());
    }
  }
}

void AliFemtoDreamv0Cuts::FillMCContributions(AliFemtoDreamv0 *v0) {
  if (!fMinimalBooking) {
    Double_t pT = v0->GetPt();
    Int_t iFill = -1;
    switch (v0->GetParticleOrigin()) {
      case AliFemtoDreamBasePart::kPhysPrimary:
        fMCHist->FillMCPrimary(pT);
        iFill = 0;
        break;
      case AliFemtoDreamBasePart::kWeak:
        fMCHist->FillMCFeeddown(pT, TMath::Abs(v0->GetMotherWeak()));
        iFill = 1;
        break;
      case AliFemtoDreamBasePart::kMaterial:
        fMCHist->FillMCMaterial(pT);
        iFill = 2;
        break;
      case AliFemtoDreamBasePart::kContamination:
        fMCHist->FillMCCont(pT);
        iFill = 3;
        break;
      case AliFemtoDreamBasePart::kFake:
        fMCHist->FillMCCont(pT);
        iFill = 3;
        break;
      default:
        AliFatal("Type Not implemented");
        break;
    }
    if (iFill != -1) {
      fMCHist->FillMCpT(iFill, pT);
      fMCHist->FillMCEta(iFill, v0->GetEta().at(0));
      fMCHist->FillMCPhi(iFill, v0->GetPhi().at(0));
      fMCHist->FillMCDCAVtxX(iFill, pT, v0->GetDCAv0Vtx(0));
      fMCHist->FillMCDCAVtxY(iFill, pT, v0->GetDCAv0Vtx(1));
      fMCHist->FillMCDCAVtxZ(iFill, pT, v0->GetDCAv0Vtx(2));
      fMCHist->FillMCTransverseRadius(iFill, pT, v0->GetTransverseRadius());
      fMCHist->FillMCDCAPosDaugPrimVtx(iFill, pT, v0->GetDCADaugPosVtx());
      fMCHist->FillMCDCANegDaugPrimVtx(iFill, pT, v0->GetDCADaugNegVtx());
      fMCHist->FillMCDCADaugVtx(iFill, pT, v0->GetDaugDCA());
      fMCHist->FillMCCosPoint(iFill, pT, v0->GetCPA());
      fMCHist->FillMCInvMass(iFill, v0->Getv0Mass());
    } else {
      std::cout << "this should not happen \n";
    }
  }
}
void AliFemtoDreamv0Cuts::BookTrackCuts() {
  if (!fMinimalBooking) {
    if (!fHist) {
      AliFatal("No Histograms available");
    }
    if (fCutOnFlyStatus) {
      fHist->FillConfig(0, 1);
    }
    if (fCutCharge) {
      fHist->FillConfig(1, fCharge);
    }
    if (fCutPt) {
      fHist->FillConfig(2, fpTmin);
      fHist->FillConfig(3, fpTmax);
    }
    if (fKaonRejection) {
      fHist->FillConfig(4, 1);
    }
    if (fCutDecayVtxXYZ) {
      fHist->FillConfig(5, fMaxDecayVtxXYZ);
    }
    if (fCutTransRadius) {
      fHist->FillConfig(6, fMinTransRadius);
      fHist->FillConfig(7, fMaxTransRadius);
    }
    if (fCutMinDCADaugPrimVtx) {
      fHist->FillConfig(8, fMinDCADaugToPrimVtx);
    }
    if (fCutMaxDCADaugToDecayVtx) {
      fHist->FillConfig(9, fMaxDCADaugToDecayVtx);
    }
    if (fCutInvMass) {
      float massv0 = TDatabasePDG::Instance()->GetParticle(fPDGv0)->Mass();
      fHist->FillConfig(10, massv0 - fInvMassCutWidth);
      fHist->FillConfig(11, massv0 + fInvMassCutWidth);
    } else if (fCutInvMassSidebands) {
      fHist->FillConfig(10, fInvMassCutSBdown);
      fHist->FillConfig(11, fInvMassCutSBup);
    }
    if (fCutCPA) {
      fHist->FillConfig(12, fMinCPA);
    }
    if (fDoArmenterosCut) {
      fHist->FillConfig(13, fArmenterosQtLow);
      fHist->FillConfig(14, fArmenterosQtUp);
      fHist->FillConfig(15, fArmenterosAlphaLow);
      fHist->FillConfig(16, fArmenterosAlphaUp);
    }
    if (fDoCombinedTimingCut) {
      fHist->FillConfig(17, fCombinedTiming);
    }
  }
}

float AliFemtoDreamv0Cuts::CalculateInvMass(AliFemtoDreamv0 *v0, int PDGPosDaug,
                                            int PDGNegDaug) {
  Double_t invMass = 0;
  float massDP = TDatabasePDG::Instance()->GetParticle(PDGPosDaug)->Mass();
  float massDN = TDatabasePDG::Instance()->GetParticle(PDGNegDaug)->Mass();

  float EDaugP = TMath::Sqrt(
      massDP * massDP
          + v0->GetPosDaughter()->GetMomentum().X()
              * v0->GetPosDaughter()->GetMomentum().X()
          + v0->GetPosDaughter()->GetMomentum().Y()
              * v0->GetPosDaughter()->GetMomentum().Y()
          + v0->GetPosDaughter()->GetMomentum().Z()
              * v0->GetPosDaughter()->GetMomentum().Z());
  float EDaugN = TMath::Sqrt(
      massDN * massDN
          + v0->GetNegDaughter()->GetMomentum().X()
              * v0->GetNegDaughter()->GetMomentum().X()
          + v0->GetNegDaughter()->GetMomentum().Y()
              * v0->GetNegDaughter()->GetMomentum().Y()
          + v0->GetNegDaughter()->GetMomentum().Z()
              * v0->GetNegDaughter()->GetMomentum().Z());

  float energysum = EDaugP + EDaugN;
  float pSum2 = (v0->GetNegDaughter()->GetMomentum().X()
      + v0->GetPosDaughter()->GetMomentum().X())
      * (v0->GetNegDaughter()->GetMomentum().X()
          + v0->GetPosDaughter()->GetMomentum().X())
      +

      (v0->GetNegDaughter()->GetMomentum().Y()
          + v0->GetPosDaughter()->GetMomentum().Y())
          * (v0->GetNegDaughter()->GetMomentum().Y()
              + v0->GetPosDaughter()->GetMomentum().Y())
      +

      (v0->GetNegDaughter()->GetMomentum().Z()
          + v0->GetPosDaughter()->GetMomentum().Z())
          * (v0->GetNegDaughter()->GetMomentum().Z()
              + v0->GetPosDaughter()->GetMomentum().Z());
  invMass = TMath::Sqrt(energysum * energysum - pSum2);
  return invMass;
}
