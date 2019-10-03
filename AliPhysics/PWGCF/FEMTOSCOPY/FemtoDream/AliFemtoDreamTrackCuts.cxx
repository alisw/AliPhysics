/*
 * AliFemtoDreamTrackCuts.cxx
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#include "AliFemtoDreamTrackCuts.h"
#include "TMath.h"
#include "AliLog.h"
#include <iostream>
ClassImp(AliFemtoDreamTrackCuts)

AliFemtoDreamTrackCuts::AliFemtoDreamTrackCuts()
    : fMCHists(nullptr),
      fHists(nullptr),
      fMinimalBooking(false),
      fMCData(false),
      fDCAPlots(false),
      fTOFM(false),
      fDoMultBinning(false),
      fCheckMother(false),
      fCombSigma(false),
      fContribSplitting(false),
      fFillQALater(false),
      fCheckFilterBit(false),
      fCheckESDFiltering(false),
      fCheckPileUpITS(false),
      fCheckPileUpSPDTOF(false),
      fCheckPileUpTOF(false),
      fCheckPileUp(false),
      fFilterBit(0),
      fpTmin(0.),
      fpTmax(0.),
      fcutPt(false),
      fetamin(0.),
      fetamax(0.),
      fcutEta(false),
      fcutCharge(false),
      fCharge(0),
      fnTPCCls(0),
      fcutnTPCCls(false),
      fMaxSharedClsTPC(160),
      fCutSharedClsTPC(false),
      fCutChi2(false),
      fMinCutChi2(0.f),
      fMaxCutChi2(0.f),
      fDCAProp(false),
      fDCAToVertexXY(0),
      fCutDCAToVtxXY(false),
      fDCAToVertexZ(0),
      fCutDCAToVtxZ(false),
      fCutSharedCls(false),
      fCheckTPCRefit(false),
      fCutTPCCrossedRows(false),
      fCrossedRows(70),
      fRatioCrossedRows(0.83),
      fCutPID(false),
      fAllowITSonly(false),
      fCutHighPtSig(false),
      fParticleID(AliPID::kUnknown),
      fNSigValue(3.),
      fNSigValueITS(3.),
      fPIDPTPCThreshold(0),
      fRejectPions(false) {
}

AliFemtoDreamTrackCuts::AliFemtoDreamTrackCuts(
    const AliFemtoDreamTrackCuts& cuts)
    : fMCHists(cuts.fMCHists),
      fHists(cuts.fHists),
      fMinimalBooking(cuts.fMinimalBooking),
      fMCData(cuts.fMCData),
      fDCAPlots(cuts.fDCAPlots),
      fTOFM(cuts.fTOFM),
      fDoMultBinning(cuts.fDoMultBinning),
      fCheckMother(cuts.fCheckMother),
      fCombSigma(cuts.fCombSigma),
      fContribSplitting(cuts.fContribSplitting),
      fFillQALater(cuts.fFillQALater),
      fCheckFilterBit(cuts.fCheckFilterBit),
      fCheckESDFiltering(cuts.fCheckESDFiltering),
      fCheckPileUpITS(cuts.fCheckPileUpITS),
      fCheckPileUpSPDTOF(cuts.fCheckPileUpSPDTOF),
      fCheckPileUpTOF(cuts.fCheckPileUpTOF),
      fCheckPileUp(cuts.fCheckPileUp),
      fFilterBit(cuts.fFilterBit),
      fpTmin(cuts.fpTmin),
      fpTmax(cuts.fpTmax),
      fcutPt(cuts.fcutPt),
      fetamin(cuts.fetamin),
      fetamax(cuts.fetamax),
      fcutEta(cuts.fcutEta),
      fcutCharge(cuts.fcutCharge),
      fCharge(cuts.fCharge),
      fnTPCCls(cuts.fnTPCCls),
      fcutnTPCCls(cuts.fcutnTPCCls),
      fMaxSharedClsTPC(cuts.fMaxSharedClsTPC),
      fCutSharedClsTPC(cuts.fCutSharedClsTPC),
      fCutChi2(cuts.fCutChi2),
      fMinCutChi2(cuts.fMinCutChi2),
      fMaxCutChi2(cuts.fMaxCutChi2),
      fDCAProp(cuts.fDCAProp),
      fDCAToVertexXY(cuts.fDCAToVertexXY),
      fCutDCAToVtxXY(cuts.fCutDCAToVtxXY),
      fDCAToVertexZ(cuts.fDCAToVertexZ),
      fCutDCAToVtxZ(cuts.fCutDCAToVtxZ),
      fCutSharedCls(cuts.fCutSharedCls),
      fCheckTPCRefit(cuts.fCheckTPCRefit),
      fCutTPCCrossedRows(cuts.fCutTPCCrossedRows),
      fCrossedRows(cuts.fCrossedRows),
      fRatioCrossedRows(cuts.fRatioCrossedRows),
      fCutPID(cuts.fCutPID),
      fAllowITSonly(false),
      fCutHighPtSig(cuts.fCutHighPtSig),
      fParticleID(cuts.fParticleID),
      fNSigValue(cuts.fNSigValue),
      fNSigValueITS(cuts.fNSigValueITS),
      fPIDPTPCThreshold(cuts.fPIDPTPCThreshold),
      fRejectPions(cuts.fRejectPions) {
}

AliFemtoDreamTrackCuts &AliFemtoDreamTrackCuts::operator =(
    const AliFemtoDreamTrackCuts& cuts) {
  if (this == &cuts) {
    return *this;
  }
  this->fMCHists = cuts.fMCHists;
  this->fHists = cuts.fHists;
  this->fMinimalBooking = cuts.fMinimalBooking;
  this->fMCData = cuts.fMCData;
  this->fDCAPlots = cuts.fDCAPlots;
  this->fTOFM= cuts.fTOFM;
  this->fDoMultBinning = cuts.fDoMultBinning;
  this->fCheckMother = cuts.fCheckMother;
  this->fCombSigma = cuts.fCombSigma;
  this->fContribSplitting = cuts.fContribSplitting;
  this->fFillQALater = cuts.fFillQALater;
  this->fCheckFilterBit = cuts.fCheckFilterBit;
  this->fCheckESDFiltering = cuts.fCheckESDFiltering;
  this->fCheckPileUpITS = cuts.fCheckPileUpITS;
  this->fCheckPileUpSPDTOF = cuts.fCheckPileUpSPDTOF;
  this->fCheckPileUpTOF = cuts.fCheckPileUpTOF;
  this->fCheckPileUp = cuts.fCheckPileUp;
  this->fFilterBit = cuts.fFilterBit;
  this->fpTmin = cuts.fpTmin;
  this->fpTmax = cuts.fpTmax;
  this->fcutPt = cuts.fcutPt;
  this->fetamin = cuts.fetamin;
  this->fetamax = cuts.fetamax;
  this->fcutEta = cuts.fcutEta;
  this->fcutCharge = cuts.fcutCharge;
  this->fCharge = cuts.fCharge;
  this->fnTPCCls = cuts.fnTPCCls;
  this->fcutnTPCCls = cuts.fcutnTPCCls;
  this->fMaxSharedClsTPC = cuts.fMaxSharedClsTPC;
  this->fCutSharedClsTPC = cuts.fCutSharedClsTPC;
  this->fCutChi2 = cuts.fCutChi2;
  this->fMinCutChi2 = cuts.fMinCutChi2;
  this->fMaxCutChi2 = cuts.fMaxCutChi2;
  this->fDCAProp = cuts.fDCAProp;
  this->fDCAToVertexXY = cuts.fDCAToVertexXY;
  this->fCutDCAToVtxXY = cuts.fCutDCAToVtxXY;
  this->fDCAToVertexZ = cuts.fDCAToVertexZ;
  this->fCutDCAToVtxZ = cuts.fCutDCAToVtxZ;
  this->fCutSharedCls = cuts.fCutSharedCls;
  this->fCheckTPCRefit = cuts.fCheckTPCRefit;
  this->fCutTPCCrossedRows = cuts.fCutTPCCrossedRows;
  this->fCrossedRows = cuts.fCrossedRows;
  this->fRatioCrossedRows = cuts.fRatioCrossedRows;
  this->fCutPID = cuts.fCutPID;
  this->fAllowITSonly = cuts.fAllowITSonly;
  this->fCutHighPtSig = cuts.fCutHighPtSig;
  this->fParticleID = cuts.fParticleID;
  this->fNSigValue = cuts.fNSigValue;
  this->fNSigValueITS = cuts.fNSigValueITS;
  this->fPIDPTPCThreshold = cuts.fPIDPTPCThreshold;
  this->fRejectPions = cuts.fRejectPions;
  return *this;
}

AliFemtoDreamTrackCuts::~AliFemtoDreamTrackCuts() {
  if (fMCHists) {
    delete fMCHists;
  }
  if (fHists) {
    delete fHists;
  }
}

bool AliFemtoDreamTrackCuts::isSelected(AliFemtoDreamTrack *Track) {
  if (!Track) {
    AliFatal("No Input Track received");
  }
  bool pass = true;
  if (!Track->IsSet()) {
    pass = false;
  } else {
    if (!fMinimalBooking)
      fHists->FillTrackCounter(0);
  }
  if (pass) {
    if (!TrackingCuts(Track)) {
      pass = false;
    }
  }
  if (pass && fCutPID) {
    if (!PIDCuts(Track)) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(22);
    }
  }
  if (pass) {
    if (!DCACuts(Track)) {
      pass = false;
    }
  }
  Track->SetUse(pass);
  if (!fFillQALater && Track->IsSet()) {
    BookQA(Track);
    if (fMCData && !fMinimalBooking) {
      BookMC(Track);
    }
  }
  return pass;
}

bool AliFemtoDreamTrackCuts::TrackingCuts(AliFemtoDreamTrack *Track) {
  bool pass = true;
  std::vector<float> eta = Track->GetEta();
  std::vector<int> charge = Track->GetCharge();
  if (fCheckFilterBit) {
    if (!Track->TestFilterBit(fFilterBit)) {
      //if there is Filterbit -1 .. see Prong cuts! Daughters don't check for FB
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(1);
    }
  }
  if (fCheckESDFiltering) {
    if (!Track->PassESDFiltering()) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(1);
    }
  }
  if (pass && fcutPt) {
    if (Track->GetPt() < fpTmin || Track->GetPt() > fpTmax) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(2);
    }
  }
  if (pass && fcutEta) {
    if (eta[0] < fetamin || eta[0] > fetamax) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(3);
    }
  }
  if (pass && fcutCharge) {
    if (!(charge[0] == fCharge)) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(4);
    }
  }
  if (pass && fCheckPileUpITS) {
    if (!Track->GetHasITSHit()) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(5);
    }
  }
  if (pass && fCheckPileUpSPDTOF) {
    if (!(Track->GetTOFTimingReuqirement() || Track->GetHasSPDHit())) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(6);
    }
  }
  if (pass && fCheckPileUpTOF) {
    if (!Track->GetTOFTimingReuqirement()) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(7);
    }
  }
  if (pass && fCheckPileUp) {
    if (!(Track->GetTOFTimingReuqirement() || Track->GetHasITSHit())) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(8);
    }
  }
  if (pass && fcutnTPCCls) {
    if (Track->GetNClsTPC() < fnTPCCls) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(9);
    }
  }
  if (pass && fCutSharedClsTPC) {
    if (Track->GetTPCClsC() > fMaxSharedClsTPC) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(10);
    }
  }
  if (pass && fCutSharedCls) {
    if (!Track->isnoSharedClst()) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(11);
    }
  }
  if (pass && fCheckTPCRefit) {
    if (!Track->GetHasTPCRefit()) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(12);
    }
  }
  if (pass && fCutTPCCrossedRows) {
    if (Track->GetTPCCrossedRows() < fCrossedRows) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(13);
    }
    if (pass) {
      if (Track->GetRatioCr() < fRatioCrossedRows) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHists->FillTrackCounter(14);
      }
    }
  }
  if (pass && fCutChi2) {
    if (Track->GetChiSquare() > fMaxCutChi2
        || Track->GetChiSquare() < fMinCutChi2) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(15);
    }
  }

  return pass;
}

bool AliFemtoDreamTrackCuts::PIDCuts(AliFemtoDreamTrack *Track) {
  bool pass = true;
  //PID Method with an nSigma cut, just use the TPC below threshold,
  //and above TPC and TOF Combined
  //there must be TPC & TOF signal (TOF for P>0.75 GeV/c)
  bool ITSisthere = false;
  bool TPCisthere = false;
  bool TOFisthere = false;

  if (Track->GetstatusITS() == AliPIDResponse::kDetPidOk) {
    ITSisthere = true;
  }
  if (Track->GetstatusTPC() == AliPIDResponse::kDetPidOk) {
    TPCisthere = true;
  }
  if (Track->GetstatusTOF() == AliPIDResponse::kDetPidOk) {
    TOFisthere = true;
  }
  //Below a threshold where the bands are well seperated in the TPC use only
  //TPC for PID, since the TOF has only limited matching efficiency. Above
  //threshold use both detectors and perform a purity check, if another
  //particle species doesn't have a smaller sigma value

  if (Track->GetMomTPC() < fPIDPTPCThreshold) {
    if (!fAllowITSonly) {
      if (!TPCisthere) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHists->FillTrackCounter(16);
        if (fRejectPions && TOFisthere) {
          float nSigTOF = (Track->GetnSigmaTOF((int) (AliPID::kPion)));
          if (TMath::Abs(nSigTOF) < fNSigValue) {
            if (fParticleID == AliPID::kPion) {
              AliWarning(
                  "Sure you want to use this method? Propably want to set"
                  " SetRejLowPtPionsTOF(kFALSE), since you are selecting Pions");
            }
            //if the particle is a Pion according to the TOF, reject it!
            pass = false;
          } else {
            if (!fMinimalBooking)
              fHists->FillTrackCounter(17);
          }
        }
        if (pass) {
          float nSigTPC = (Track->GetnSigmaTPC((int) (fParticleID)));
          if (!(TMath::Abs(nSigTPC) < fNSigValue)) {
            pass = false;
          } else {
            if (!fMinimalBooking)
              fHists->FillTrackCounter(18);
          }
        }
      }
    } else {  //exception for omega bachelor: enable use of ITS pid
      if (!ITSisthere && !TPCisthere) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHists->FillTrackCounter(16);
        if (fRejectPions && TOFisthere) {
          float nSigTOF = (Track->GetnSigmaTOF((int) (AliPID::kPion)));
          if (TMath::Abs(nSigTOF) < fNSigValue) {
            if (fParticleID == AliPID::kPion) {
              AliWarning(
                  "Sure you want to use this method? Propably want to set"
                  " SetRejLowPtPionsTOF(kFALSE), since you are selecting Pions");
            }
            //if the particle is a Pion according to the TOF, reject it!
            pass = false;
          } else {
            if (!fMinimalBooking)
              fHists->FillTrackCounter(17);
          }
        }
        if (pass) {
          if (TPCisthere) {  //tpc prevails
            float nSigTPC = (Track->GetnSigmaTPC((int) (fParticleID)));
            if (!(TMath::Abs(nSigTPC) < fNSigValue)) {
              pass = false;
            } else {
              if (!fMinimalBooking)
                fHists->FillTrackCounter(18);
            }
          } else {
            if (ITSisthere) {  //if there is no tpc, check its
              float nSigITS = (Track->GetnSigmaITS((int) (fParticleID)));
              if (!(TMath::Abs(nSigITS) < fNSigValueITS)) {
                pass = false;
              } else {
                if (!fMinimalBooking)
                  fHists->FillTrackCounter(19);
              }
            }
          }
        }
      }
    }  // if fallowitsonly
  } else {
    if (!(TPCisthere && TOFisthere)) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHists->FillTrackCounter(20);
      float nSigTPC = (Track->GetnSigmaTPC((int) (fParticleID)));
      float nSigTOF = (Track->GetnSigmaTOF((int) (fParticleID)));
      float nSigComb = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);
      if (!(nSigComb < fNSigValue)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHists->FillTrackCounter(21);
        if (fCutHighPtSig) {
          if (!SmallestNSig(Track)) {
            pass = false;
          } else {
            if (!fMinimalBooking)
              fHists->FillTrackCounter(22);
          }
        }
      }
    }
  }
  return pass;
}

bool AliFemtoDreamTrackCuts::SmallestNSig(AliFemtoDreamTrack *Track) {
  bool pass = true;
  //check before if TPC and TOF PID are available
  //This should just be for PID of high pT particles
  AliPID::EParticleType type[6] = { AliPID::kElectron, AliPID::kMuon,
      AliPID::kPion, AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron };
  float nSigmaComb[6];
  //Form the combination:
  for (int i = 0; i < 6; ++i) {
    nSigmaComb[i] = TMath::Sqrt(
        pow((Track->GetnSigmaTPC(i)), 2.) + pow((Track->GetnSigmaTOF(i)), 2.));
  }
  int index = 0;
  for (int i = 0; i < 6; ++i) {
    if (nSigmaComb[index] > nSigmaComb[i]) {
      index = i;
    }
  }
  if (!(type[index] == fParticleID)) {
    pass = false;
  }
  return pass;
}
bool AliFemtoDreamTrackCuts::DCACuts(AliFemtoDreamTrack *Track) {
  bool pass = true;

  if (fCutDCAToVtxZ) {
    if (fDCAProp) {
      if (!(TMath::Abs(Track->GetDCAZProp()) < fDCAToVertexZ)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHists->FillTrackCounter(24);
      }
    } else {
      if (!(TMath::Abs(Track->GetDCAZ()) < fDCAToVertexZ)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHists->FillTrackCounter(24);
      }
    }
  }
  if (pass && fDCAPlots && !fMinimalBooking) {
    if (fDCAProp) {
      fHists->FillDCAXYPtBins(Track->GetPt(), Track->GetDCAXYProp(),
                              Track->GetEventMultiplicity());
    } else {
      fHists->FillDCAXYPtBins(Track->GetPt(), Track->GetDCAXY(),
                              Track->GetEventMultiplicity());
    }
    if (fMCData) {
      if (fDCAPlots) {
        if (fDCAProp) {
          fMCHists->FillMCDCAXYPtBins(Track->GetParticleOrigin(),
                                      Track->GetMotherWeak(), Track->GetPt(),
                                      Track->GetDCAXYProp(),
                                      Track->GetEventMultiplicity());
        } else {
          fMCHists->FillMCDCAXYPtBins(Track->GetParticleOrigin(),
                                      Track->GetMotherWeak(), Track->GetPt(),
                                      Track->GetDCAXY(),
                                      Track->GetEventMultiplicity());
        }
      }
    }
  }
  if (pass && fCutDCAToVtxXY) {
    if (fDCAProp) {
      if (!(TMath::Abs(Track->GetDCAXYProp()) < fDCAToVertexXY)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHists->FillTrackCounter(25);
      }
    } else {
      if (!(TMath::Abs(Track->GetDCAXY()) < fDCAToVertexXY)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHists->FillTrackCounter(25);
      }
    }
  }
  return pass;
}

void AliFemtoDreamTrackCuts::Init(TString name) {
  if (!fMinimalBooking) {
    fHists = new AliFemtoDreamTrackHist(fDCAPlots, fCombSigma, fTOFM);
    if (fMCData) {
      fMCHists = new AliFemtoDreamTrackMCHist(fContribSplitting, fDCAPlots,
                                              fDoMultBinning, fCheckMother);
    }
    BookTrackCuts();
  } else {
    fHists = new AliFemtoDreamTrackHist(name);
  }
}

void AliFemtoDreamTrackCuts::BookQA(AliFemtoDreamTrack *Track) {
  if (!fMinimalBooking) {
    std::vector<float> eta = Track->GetEta();
    std::vector<float> phi = Track->GetPhi();
    float pT = Track->GetPt();
    float p = Track->GetMomTPC();
    float Pprim = Track->GetP();
    for (int i = 0; i < 2; ++i) {
      if (i == 0 || (i == 1 && Track->UseParticle())) {
        fHists->FilletaCut(i, eta.at(0));
        fHists->FillphiCut(i, phi.at(0));
        fHists->FillpTCut(i, pT);
        fHists->FillpTPCCut(i, p);
        fHists->FillTPCclsCut(i, Track->GetNClsTPC());
        fHists->FillTrackChiSquare(i, pT, Track->GetChiSquare());
        fHists->FillDCAxyCut(i, pT, Track->GetDCAXY());
        fHists->FillDCAzCut(i, pT, Track->GetDCAZ());
        fHists->FillDCAxyPropCut(i, pT, Track->GetDCAXYProp());
        fHists->FillDCAzPropCut(i, pT, Track->GetDCAZProp());
        fHists->FillTPCCrossedRowCut(i, Track->GetTPCCrossedRows());
        fHists->FillTPCRatioCut(i, Track->GetRatioCr());
        fHists->FillTPCClsS(i, Track->GetTPCClsC());
        for (size_t j = 0; j < Track->GetITSHits().size(); ++j) {
          if (Track->GetITSHit(j)) {
            fHists->FillTPCClsCPileUp(i, j, Track->GetTPCClsC());
          } else if (Track->GetHasITSHit()
              || Track->GetTOFTimingReuqirement()) {
            fHists->FillTPCClsCPileUp(i, j + 7, Track->GetTPCClsC());
          }
        }
        if (Track->GetTOFTimingReuqirement()) {
          fHists->FillTPCClsCPileUp(i, 6, Track->GetTPCClsC());
        } else if (Track->GetHasITSHit()) {
          fHists->FillTPCClsCPileUp(i, 13, Track->GetTPCClsC());
        } else {
          fHists->FillTPCClsCPileUp(i, 14, Track->GetTPCClsC());
        }

        for (size_t j = 0; j < Track->GetSharedClusterITS().size(); ++j) {
          if (Track->GetSharedClusterITS(j)) {
            fHists->FillHasSharedClsITS(i, j + 1, 0);
          } else {
            fHists->FillHasSharedClsITS(i, j + 1, 1);
          }
          for (int k = 0; k < 6; ++k) {
            if (Track->GetITSHit(k)) {
              if (Track->GetSharedClusterITS(j)) {
                fHists->FillITSSharedPileUp(i, k, j);
              } else {
                fHists->FillITSSharedPileUp(i, k, j + 6);
              }
            } else if (Track->GetHasITSHit()
                || Track->GetTOFTimingReuqirement()) {
              if (Track->GetSharedClusterITS(j)) {
                fHists->FillITSSharedPileUp(i, k + 7, j);
              } else {
                fHists->FillITSSharedPileUp(i, k + 7, j + 6);
              }
            }
          }
          if (Track->GetTOFTimingReuqirement()) {
            if (Track->GetSharedClusterITS(j)) {
              fHists->FillITSSharedPileUp(i, 6, j);
            } else {
              fHists->FillITSSharedPileUp(i, 6, j + 6);
            }
          } else if (Track->GetHasITSHit()) {
            if (Track->GetSharedClusterITS(j)) {
              fHists->FillITSSharedPileUp(i, 13, j);
            } else {
              fHists->FillITSSharedPileUp(i, 13, j + 6);
            }
          } else {
            if (Track->GetSharedClusterITS(j)) {
              fHists->FillITSSharedPileUp(i, 14, j);
            } else {
              fHists->FillITSSharedPileUp(i, 14, j + 6);
            }
          }
        }
        fHists->FillTPCdedx(i, p, Track->GetdEdxTPC());
        fHists->FillTOFbeta(i, p, Track->GetbetaTOF());

        fHists->FillNSigTPC(i, p, (Track->GetnSigmaTPC(fParticleID)));
        fHists->FillNSigTPCMod(i, p, (Track->GetnSigmaTPC(fParticleID)));
        fHists->FillNSigTOF(i, p, (Track->GetnSigmaTOF(fParticleID)));
        fHists->FillTPCStatus(i, Track->GetstatusTPC());
        fHists->FillTOFStatus(i, Track->GetstatusTOF());
        fHists->FillNSigComTPCTOF(i, p, Track->GetnSigmaTPC(fParticleID),
                                  Track->GetnSigmaTOF(fParticleID));
        //Fill These Before
        if (i == 0 && fCombSigma) {
          fHists->FillNSigComb(pT, Track->GetnSigmaTPC(fParticleID),
                               Track->GetnSigmaTOF(fParticleID));
        }
        //Fill These After
        if (i == 1 && fTOFM) {
            fHists->FillTOFMass(Pprim, Track->GetbetaTOF());
        }
      }
    }
  } else {
    if (Track->UseParticle())
      fHists->FillpTCut(1, Track->GetPt());
  }
  return;
}

void AliFemtoDreamTrackCuts::BookMC(AliFemtoDreamTrack *Track) {
  if (!fMinimalBooking) {
    Int_t PDGcode[6] = { 11, 13, 211, 321, 2212, 1000010020 };
    //this is not the correct way to do it, since there might be float counting
    if (fpTmin < Track->GetPt() && Track->GetPt() < fpTmax) {
      if (fetamin < Track->GetEta().at(0) && Track->GetEta().at(0) < fetamax) {
        if (!fcutCharge) {
          if (TMath::Abs(Track->GetMCPDGCode()) == PDGcode[fParticleID]) {
            fMCHists->FillMCGen(Track->GetPt());
          }
        } else {
          Int_t sign = (Int_t) fCharge / TMath::Abs(fCharge);
          if (Track->GetMCPDGCode() == sign * PDGcode[fParticleID]) {
            fMCHists->FillMCGen(Track->GetPt());
          }
        }
      }
    }
    if (Track->UseParticle()) {
      float pT = Track->GetPt();
      int PDGcode[6] = { 11, 13, 211, 321, 2212, 1000010020};
//Fill Identified
      fMCHists->FillMCIdent(pT);
      AliFemtoDreamBasePart::PartOrigin tmpOrg = Track->GetParticleOrigin();
      if (!fcutCharge) {
        if (TMath::Abs(Track->GetMCPDGCode())
            == TMath::Abs(PDGcode[fParticleID])) {
          fMCHists->FillMCCorr(pT);
          if (tmpOrg == AliFemtoDreamBasePart::kPhysPrimary) {
            fMCHists->FillMCPtResolution(Track->GetMCPt(), Track->GetPt());
            float phi = Track->GetMomentum().Phi();
            if (phi < 0) {
              //Root handles phi from - pi to + pi, while AliROOT handles it from 0 to 2 pi
              phi += 2 * TMath::Pi();
            }
            fMCHists->FillMCThetaResolution(Track->GetMCTheta().at(0),
                                            Track->GetMomentum().Theta(),
                                            Track->GetMCPt());
            fMCHists->FillMCPhiResolution(Track->GetMCPhi().at(0), phi,
                                          Track->GetMCPt());
          }
        } else {
          Track->SetParticleOrigin(AliFemtoDreamBasePart::kContamination);
        }
      } else {
        Int_t sign = (Int_t) fCharge / TMath::Abs(fCharge);
        if (Track->GetMCPDGCode() == sign * PDGcode[fParticleID]) {
          fMCHists->FillMCCorr(pT);
          if (tmpOrg == AliFemtoDreamBasePart::kPhysPrimary) {
            fMCHists->FillMCPtResolution(Track->GetMCPt(), Track->GetPt());
            fMCHists->FillMCThetaResolution(Track->GetMCTheta().at(0),
                                            Track->GetMomentum().Theta(),
                                            Track->GetMCPt());
            float phi = Track->GetMomentum().Phi();
//            std::cout << phi << '\t' << Track->GetMCPhi().at(0) << std::endl;
            if (phi < 0) {
              //Root handles phi from - pi to + pi, while AliROOT handles it from 0 to 2 pi
              phi += 2 * TMath::Pi();
            }
//            std::cout << phi << '\t' << Track->GetMCPhi().at(0) << std::endl;
            fMCHists->FillMCPhiResolution(Track->GetMCPhi().at(0), phi,
                                          Track->GetMCPt());
          }
        } else {
          Track->SetParticleOrigin(AliFemtoDreamBasePart::kContamination);
        }
      }
      if (fContribSplitting) {
        FillMCContributions(Track);
      }
      Track->SetParticleOrigin(tmpOrg);
      if (fCheckMother)
        fMCHists->FillMCMother(Track->GetPt(), Track->GetMotherPDG());
    }
  }
}

void AliFemtoDreamTrackCuts::FillMCContributions(AliFemtoDreamTrack *Track) {
  if (!fMinimalBooking) {
    float pT = Track->GetPt();
    AliFemtoDreamBasePart::PartOrigin org = Track->GetParticleOrigin();
    Int_t iFill = -1;
    switch (org) {
      case AliFemtoDreamBasePart::kPhysPrimary:
        fMCHists->FillMCPrimary(pT);
        iFill = 0;
        break;
      case AliFemtoDreamBasePart::kWeak:
        fMCHists->FillMCFeeddown(pT, TMath::Abs(Track->GetMotherWeak()));
        iFill = 1;
        break;
      case AliFemtoDreamBasePart::kMaterial:
        fMCHists->FillMCMaterial(pT);
        iFill = 2;
        break;
      case AliFemtoDreamBasePart::kContamination:
        fMCHists->FillMCCont(pT);
        iFill = 3;
        break;
      default:
        AliFatal("Type Not implemented");
        break;
    }
    if (iFill >= 0 && iFill < 4) {
      std::vector<float> eta = Track->GetEta();
      std::vector<float> phi = Track->GetPhi();
      fMCHists->FillMCpTPCCut(iFill, Track->GetMomTPC());
      fMCHists->FillMCetaCut(iFill, eta[0]);
      fMCHists->FillMCphiCut(iFill, phi[0]);
      fMCHists->FillMCTPCclsCut(iFill, pT, Track->GetNClsTPC());
      if (fDCAProp) {
        fMCHists->FillMCDCAxyCut(iFill, pT, Track->GetDCAXYProp());
        fMCHists->FillMCDCAzCut(iFill, pT, Track->GetDCAZProp());
      } else {
        fMCHists->FillMCDCAxyCut(iFill, pT, Track->GetDCAXY());
        fMCHists->FillMCDCAzCut(iFill, pT, Track->GetDCAZ());
      }
      fMCHists->FillMCTPCCrossedRowCut(iFill, pT, Track->GetTPCCrossedRows());
      fMCHists->FillMCTPCRatioCut(iFill, pT, Track->GetRatioCr());
      fMCHists->FillMCTPCdedx(iFill, pT, Track->GetdEdxTPC());
      fMCHists->FillMCTOFbeta(iFill, pT, Track->GetbetaTOF());
      fMCHists->FillMCNSigTPC(iFill, pT, Track->GetnSigmaTPC(fParticleID));
      fMCHists->FillMCNSigTOF(iFill, pT, Track->GetnSigmaTOF(fParticleID));
    } else {
      TString errMSG = Form("iFill = %d", iFill);
      AliFatal(errMSG.Data());
    }
  }
  return;
}

void AliFemtoDreamTrackCuts::BookTrackCuts() {
  if (!fMinimalBooking) {
    if (!fHists) {
      AliFatal(
          "AliFemtoPPbpbLamSpTrackCuts::BookTrackCuts No Histograms to work with");
    }
    if (fcutPt) {
      fHists->FillConfig(0, fpTmin);
      fHists->FillConfig(1, fpTmax);
    }
    if (fcutEta) {
      fHists->FillConfig(2, fetamin);
      fHists->FillConfig(3, fetamax);
    }
    if (fcutCharge) {
      fHists->FillConfig(4, fCharge);
    }

    if (fcutnTPCCls) {
      fHists->FillConfig(5, fnTPCCls);
    }

    if (fCheckFilterBit) {
      fHists->FillConfig(6, fFilterBit);
    } else {
      fHists->FillConfig(6, -1);
    }
    if (fMCData) {
      fHists->FillConfig(7, 1);
    }
    if (fCutDCAToVtxXY) {
      fHists->FillConfig(8, fDCAToVertexXY);
    }
    if (fCutDCAToVtxZ) {
      fHists->FillConfig(9, fDCAToVertexZ);
    }
    if (fCutSharedClsTPC) {
      fHists->FillConfig(10, fMaxSharedClsTPC);
    }
    if (fCutSharedCls) {
      fHists->FillConfig(11, 1);
    }
    if (fCutTPCCrossedRows) {
      fHists->FillConfig(12, fCrossedRows);
      fHists->FillConfig(13, fRatioCrossedRows);
    } else {
      fHists->FillConfig(12, 0);
      fHists->FillConfig(13, 0);
    }
    if (fCutPID) {
      fHists->FillConfig(14, fPIDPTPCThreshold);
      fHists->FillConfig(15, fNSigValue);
      if (fAllowITSonly) {
        fHists->FillConfig(16, fNSigValueITS);
      } else {
        fHists->FillConfig(16, 0);
      }
      if (fRejectPions) {
        fHists->FillConfig(17, 1);
      } else {
        fHists->FillConfig(17, 0);
      }
      if (fCutHighPtSig) {
        fHists->FillConfig(18, 1);
      } else {
        fHists->FillConfig(18, 0);
      }
    } else {
      fHists->FillConfig(14, 0);
      fHists->FillConfig(15, 0);
      fHists->FillConfig(16, 0);
      fHists->FillConfig(17, 0);
      fHists->FillConfig(18, 0);
    }
    if (fCheckPileUpITS) {
      fHists->FillConfig(19, 1);
    }
    if (fCheckPileUpSPDTOF) {
      fHists->FillConfig(20, 1);
    }
    if (fCheckPileUpTOF) {
      fHists->FillConfig(21, 1);
    }
    if (fCheckPileUp) {
      fHists->FillConfig(22, 1);
    }
    if (fCheckTPCRefit) {
      fHists->FillConfig(23, 1);
    }
    if (fCutChi2) {
      fHists->FillConfig(24, fMinCutChi2);
      fHists->FillConfig(25, fMaxCutChi2);
    }
    if (fCheckESDFiltering) {
      fHists->FillConfig(26, 1);
    }
  }
}

AliFemtoDreamTrackCuts* AliFemtoDreamTrackCuts::PrimProtonCuts(
    bool isMC, bool DCAPlots, bool CombSigma, bool ContribSplitting) {
  AliFemtoDreamTrackCuts *trackCuts = new AliFemtoDreamTrackCuts();
  //you can leave DCA cut active, this will still be filled
  //over the whole DCA_xy range
  trackCuts->SetPlotDCADist(DCAPlots);
  trackCuts->SetPlotCombSigma(CombSigma);
  trackCuts->SetPlotContrib(ContribSplitting);
  trackCuts->SetIsMonteCarlo(isMC);

  trackCuts->SetFilterBit(128);
  trackCuts->SetPtRange(0.5, 4.05);
  trackCuts->SetEtaRange(-0.8, 0.8);
  trackCuts->SetNClsTPC(80);
  trackCuts->SetDCAReCalculation(true);  //Get the dca from the PropagateToVetex
  trackCuts->SetDCAVtxZ(0.2);
  trackCuts->SetDCAVtxXY(0.1);
  trackCuts->SetCutSharedCls(true);
  trackCuts->SetCutTPCCrossedRows(true, 70, 0.83);
  trackCuts->SetPID(AliPID::kProton, 0.75);
  trackCuts->SetRejLowPtPionsTOF(true);
  trackCuts->SetCutSmallestSig(true);

  return trackCuts;
}

AliFemtoDreamTrackCuts *AliFemtoDreamTrackCuts::PrimKaonCuts(
    bool isMC, bool DCAPlots, bool CombSigma, bool ContribSplitting) {
  AliFemtoDreamTrackCuts *trackCuts = new AliFemtoDreamTrackCuts();
  trackCuts->SetPlotDCADist(DCAPlots);
  trackCuts->SetPlotCombSigma(CombSigma);
  trackCuts->SetPlotContrib(ContribSplitting);
  trackCuts->SetIsMonteCarlo(isMC);

  trackCuts->SetFilterBit(96);
  trackCuts->SetPtRange(0.15, 999);
  trackCuts->SetEtaRange(-0.8, 0.8);
  trackCuts->SetNClsTPC(80);
  trackCuts->SetDCAReCalculation(true);
  trackCuts->SetDCAVtxZ(0.2);
  trackCuts->SetDCAVtxXY(0.1);
  trackCuts->SetCutSharedCls(true);
  trackCuts->SetCutTPCCrossedRows(true, 70, 0.80);
  trackCuts->SetPID(AliPID::kKaon, 0.4, 5);
  // trackCuts->SetRejLowPtPionsTOF(false);
  trackCuts->SetCutSmallestSig(true);

  return trackCuts;
}

AliFemtoDreamTrackCuts* AliFemtoDreamTrackCuts::PrimDeuteronCuts(
    bool isMC, bool DCAPlots, bool CombSigma, bool ContribSplitting) {
  AliFemtoDreamTrackCuts *trackCuts = new AliFemtoDreamTrackCuts();
  //you can leave DCA cut active, this will still be filled
  //over the whole DCA_xy range
  trackCuts->SetPlotDCADist(DCAPlots);
  trackCuts->SetPlotCombSigma(CombSigma);
  trackCuts->SetPlotContrib(ContribSplitting);
  trackCuts->SetIsMonteCarlo(isMC);

  trackCuts->SetFilterBit(128);
  trackCuts->SetPtRange(0.4, 4.);
  trackCuts->SetEtaRange(-0.8, 0.8);
  trackCuts->SetNClsTPC(80);
  trackCuts->SetDCAReCalculation(true);  //Get the dca from the PropagateToVetex
  trackCuts->SetDCAVtxZ(0.2);
  trackCuts->SetDCAVtxXY(0.1);
  trackCuts->SetCutSharedCls(true);
  trackCuts->SetCutTPCCrossedRows(true, 70, 0.83);
  trackCuts->SetPID(AliPID::kDeuteron, 1.4);
  trackCuts->SetRejLowPtPionsTOF(true);
  trackCuts->SetCutSmallestSig(true);

  return trackCuts;
}

AliFemtoDreamTrackCuts* AliFemtoDreamTrackCuts::DecayProtonCuts(
    bool isMC, bool PileUpRej, bool ContribSplitting) {
  AliFemtoDreamTrackCuts *trackCuts = new AliFemtoDreamTrackCuts();
  trackCuts->SetPlotDCADist(false);
  trackCuts->SetPlotCombSigma(false);
  trackCuts->SetPlotContrib(ContribSplitting);
  trackCuts->SetIsMonteCarlo(isMC);
  trackCuts->SetFillQALater(true);

  trackCuts->SetCheckPileUp(PileUpRej);
  trackCuts->SetCheckFilterBit(false);
  trackCuts->SetEtaRange(-0.8, 0.8);
  trackCuts->SetNClsTPC(70);
  trackCuts->SetDCAReCalculation(true);
  trackCuts->SetCutCharge(1);
  trackCuts->SetPID(AliPID::kProton, 999., 5);
  return trackCuts;
}

AliFemtoDreamTrackCuts* AliFemtoDreamTrackCuts::DecayPionCuts(
    bool isMC, bool PileUpRej, bool ContribSplitting) {
  AliFemtoDreamTrackCuts *trackCuts = new AliFemtoDreamTrackCuts();
  trackCuts->SetPlotDCADist(false);
  trackCuts->SetPlotCombSigma(false);
  trackCuts->SetPlotContrib(ContribSplitting);
  trackCuts->SetIsMonteCarlo(isMC);
  trackCuts->SetFillQALater(true);

  trackCuts->SetCheckPileUp(PileUpRej);
  trackCuts->SetCheckFilterBit(kFALSE);
  trackCuts->SetEtaRange(-0.8, 0.8);
  trackCuts->SetNClsTPC(70);
  trackCuts->SetDCAReCalculation(kTRUE);
  trackCuts->SetCutCharge(-1);
  trackCuts->SetPID(AliPID::kPion, 999., 5);
  return trackCuts;
}

AliFemtoDreamTrackCuts* AliFemtoDreamTrackCuts::Xiv0PionCuts(
    bool isMC, bool PileUpRej, bool ContribSplitting) {
  AliFemtoDreamTrackCuts *trackCuts = new AliFemtoDreamTrackCuts();
  trackCuts->SetPlotDCADist(false);
  trackCuts->SetPlotCombSigma(false);
  trackCuts->SetPlotContrib(ContribSplitting);
  trackCuts->SetIsMonteCarlo(isMC);
  trackCuts->SetFillQALater(true);

  trackCuts->SetCheckPileUp(PileUpRej);
  trackCuts->SetCheckFilterBit(kFALSE);
  trackCuts->SetEtaRange(-0.8, 0.8);
//  trackCuts->SetPtRange(0.3,999);
  trackCuts->SetCutTPCCrossedRows(true, 70, 0.83);
  trackCuts->SetDCAReCalculation(kTRUE);
  trackCuts->SetCutCharge(-1);
  trackCuts->SetCheckTPCRefit(true);
  trackCuts->SetPID(AliPID::kPion, 999., 4);
  return trackCuts;
}

AliFemtoDreamTrackCuts* AliFemtoDreamTrackCuts::Xiv0ProtonCuts(
    bool isMC, bool PileUpRej, bool ContribSplitting) {
  AliFemtoDreamTrackCuts *trackCuts = new AliFemtoDreamTrackCuts();
  trackCuts->SetPlotDCADist(false);
  trackCuts->SetPlotCombSigma(false);
  trackCuts->SetPlotContrib(ContribSplitting);
  trackCuts->SetIsMonteCarlo(isMC);
  trackCuts->SetFillQALater(true);

  trackCuts->SetCheckPileUp(PileUpRej);
  trackCuts->SetCheckFilterBit(kFALSE);
  trackCuts->SetEtaRange(-0.8, 0.8);
//  trackCuts->SetPtRange(0.3,999);
  trackCuts->SetCutTPCCrossedRows(true, 70, 0.83);
  trackCuts->SetDCAReCalculation(kTRUE);
  trackCuts->SetCutCharge(1);
  trackCuts->SetCheckTPCRefit(true);
  trackCuts->SetPID(AliPID::kProton, 999., 4);
  return trackCuts;
}

AliFemtoDreamTrackCuts* AliFemtoDreamTrackCuts::XiBachPionCuts(
    bool isMC, bool PileUpRej, bool ContribSplitting) {
  AliFemtoDreamTrackCuts *trackCuts = new AliFemtoDreamTrackCuts();
  trackCuts->SetPlotDCADist(false);
  trackCuts->SetPlotCombSigma(false);
  trackCuts->SetPlotContrib(ContribSplitting);
  trackCuts->SetIsMonteCarlo(isMC);
  trackCuts->SetFillQALater(true);

  trackCuts->SetCheckPileUp(PileUpRej);
  trackCuts->SetCheckFilterBit(kFALSE);
  trackCuts->SetEtaRange(-0.8, 0.8);
//  trackCuts->SetPtRange(0.3,999);
  trackCuts->SetCutTPCCrossedRows(true, 70, 0.83);
  trackCuts->SetDCAReCalculation(kTRUE);
  trackCuts->SetCutCharge(-1);
  trackCuts->SetCheckTPCRefit(true);
  trackCuts->SetPID(AliPID::kPion, 999., 4);
  return trackCuts;
}

AliFemtoDreamTrackCuts* AliFemtoDreamTrackCuts::OmegaBachKaonCuts(
    bool isMC, bool PileUpRej, bool ContribSplitting) {
  AliFemtoDreamTrackCuts *trackCuts = new AliFemtoDreamTrackCuts();
  trackCuts->SetPlotDCADist(false);
  trackCuts->SetPlotCombSigma(false);
  trackCuts->SetPlotContrib(ContribSplitting);
  trackCuts->SetIsMonteCarlo(isMC);
  trackCuts->SetFillQALater(true);

  trackCuts->SetCheckPileUp(PileUpRej);
  trackCuts->SetCheckFilterBit(kFALSE);
  trackCuts->SetEtaRange(-0.8, 0.8);
//  trackCuts->SetPtRange(0.3,999);
  trackCuts->SetCutTPCCrossedRows(true, 70, 0.83);
  trackCuts->SetDCAReCalculation(kTRUE);
  trackCuts->SetCutCharge(-1);
  trackCuts->SetCheckTPCRefit(true);
  trackCuts->SetPID(AliPID::kKaon, 999., 4);
  return trackCuts;
}

int AliFemtoDreamTrackCuts::GetPDGCode() {
  int PDGcode[6] = { 11, 13, 211, 321, 2212, 1000010020 };
  if (fParticleID < 6) {
    return fCharge * PDGcode[fParticleID];
  } else {
    return 0;
  }
}
;
