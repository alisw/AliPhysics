/*
 * AliFemtoDreamTrack.cxx
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */
#include "AliAnalysisManager.h"
#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliInputEventHandler.h"
#include "AliFemtoDreamTrack.h"
#include "AliLog.h"
#include "TClonesArray.h"
#include "TMath.h"
#include <iostream>
ClassImp(AliFemtoDreamTrack)
AliFemtoDreamTrack::AliFemtoDreamTrack()
    : AliFemtoDreamBasePart(1),
      fPIDResponse(0),
      fstatusITS(AliPIDResponse::kDetNoParams),
      fstatusTPC(AliPIDResponse::kDetNoParams),
      fstatusTOF(AliPIDResponse::kDetNoParams),
      fFilterMap(0),
      fPassFiltering(false),
      fdcaXY(-99),
      fdcaZ(-99),
      fdcaXYProp(-99),
      fdcaZProp(-99),
      fNClsTPC(0),
      fTPCCrossedRows(0),
      fRatioCR(0),
      fnoSharedClst(0),
      fTPCClsS(0),
      fTPCClsSRatio(0),
      fChi2(0),
      fChi2TPC(0),
      fChi2ITS(0),
      fSharedClsITSLayer(0),
      fHasSharedClsITSLayer(false),
      fdEdxITS(0),
      fdEdxTPC(0),
      fbetaTOF(0),
      fHasITSHit(false),
      fITSHit(0),
      fTOFTiming(false),
      fTPCRefit(false),
      fESDStatus(0),
      fESDnClusterITS(0),
      fESDnClusterTPC(0),
      fESDTrack(0),
      fESDTPCOnlyTrack(0),
      fESDTrackCuts(AliESDtrackCuts::GetStandardTPCOnlyTrackCuts()),
      fVTrack(nullptr),
      fVGlobalTrack(nullptr),
      fAODTrack(0),
      fAODGlobalTrack(0) {
  for (int i = 0; i < 5; ++i) {
    fnSigmaITS[i] = 0;
    fnSigmaTPC[i] = 0;
    fnSigmaTOF[i] = 0;
  }
}

AliFemtoDreamTrack::~AliFemtoDreamTrack() {
  if (fESDTrackCuts) {
    delete fESDTrackCuts;
  }
  if (fESDTPCOnlyTrack) {
    delete fESDTPCOnlyTrack;
  }
}

void AliFemtoDreamTrack::SetTrack(AliAODTrack *track, const int multiplicity) {
  this->Reset();
  SetEventMultiplicity(multiplicity);
  fAODTrack = track;
  int trackID = fAODTrack->GetID();
  if (trackID < 0) {
    if (!fGTI) {
      AliFatal("AliFemtoSPTrack::SetTrack No fGTI Set");
      fAODGlobalTrack = NULL;
    } else if (-trackID - 1 >= fTrackBufferSize) {
      //      AliFatal("Buffer Size too small");
      this->fIsSet = false;
      fIsReset = false;
      fAODGlobalTrack = NULL;
    } else if (!CheckGlobalTrack(trackID)) {
      fAODGlobalTrack = NULL;
    } else {
      fAODGlobalTrack = fGTI[-trackID - 1];
    }
  } else {
    fAODGlobalTrack = track;
  }
  fIsReset = false;
  if (fAODGlobalTrack && fAODTrack) {
    this->SetIDTracks(fAODGlobalTrack->GetID());
    this->SetAODTrackingInformation();
    this->SetAODPIDInformation();
    this->SetEvtNumber(fAODTrack->GetAODEvent()->GetRunNumber());
    if (fIsMC) {
      this->SetMCInformation();
    }
  } else {
    this->fIsSet = false;
  }
}

void AliFemtoDreamTrack::SetTrack(AliVTrack *track, AliVEvent *event,
                                  const int multiplicity) {
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man
        ->GetInputEventHandler());
    if (inputHandler) {
      fPIDResponse = inputHandler->GetPIDResponse();
    }
  }
  AliNanoAODTrack* nanoTrack = dynamic_cast<AliNanoAODTrack*>(track);
  this->Reset();
  SetEventMultiplicity(multiplicity);
  fVTrack = track;
  int trackID = nanoTrack->GetID();
  if (trackID < 0) {
    if (!fVGTI) {
      AliFatal("AliFemtoSPTrack::SetTrack No fVGTI Set");
      fVGlobalTrack = NULL;
    } else if (-trackID - 1 >= fTrackBufferSize) {
      //      AliFatal("Buffer Size too small");
      this->fIsSet = false;
      fIsReset = false;
      fVGlobalTrack = NULL;
    } else if (!CheckGlobalVTrack(trackID)) {
      fVGlobalTrack = NULL;
    } else {
      fVGlobalTrack = fVGTI[-trackID - 1];
    }
  } else {
    fVGlobalTrack = track;
  }
  fIsReset = false;
  if (fVGlobalTrack && fVTrack) {
    SetVInformation(event);
    this->SetIDTracks(dynamic_cast<AliNanoAODTrack *>(fVGlobalTrack)->GetID());
    this->SetEvtNumber(event->GetRunNumber());
    if (fIsMC) {
      this->SetMCInformation(event);
    }
  } else {
    this->fIsSet = false;
  }
}

void AliFemtoDreamTrack::SetTrack(AliESDtrack *track, AliMCEvent *mcEvent,
                                  const int multiplicity,
                                  const bool TPCOnlyTrack,
                                  const bool IsOmegaTrack) {
  this->Reset();
  SetEventMultiplicity(multiplicity);
  fESDTrack = track;
  if (fESDTrack) {
    this->fIsReset = false;
    int trackID = fESDTrack->GetID();
    if (trackID < 0) {
      AliWarning("Negative ID for ESD tracks");
    }
    this->SetIDTracks(trackID);
    if(!IsOmegaTrack){
     this->SetESDTrackingInformation(TPCOnlyTrack);
    }else{
     this->SetESDTrackingInformationOmega();
    }
    this->SetEvtNumber(fESDTrack->GetESDEvent()->GetRunNumber());
    if (this->fIsSet) {
      this->SetESDPIDInformation();
      if (fIsMC) {
        this->SetMCInformation(mcEvent);
      }
    }
    if (TPCOnlyTrack && fESDTPCOnlyTrack) {
      delete fESDTPCOnlyTrack;
    }
  }
}

void AliFemtoDreamTrack::ApplyESDtoAODFilter(const bool TPCOnlyTrack) {
  if (fESDTrackCuts->IsSelected(fESDTrack)) {
    fPassFiltering = true;
  } else {
    fPassFiltering = false;
  }
  if (TPCOnlyTrack) {
    fESDTPCOnlyTrack = AliESDtrackCuts::GetTPCOnlyTrack(
        fESDTrack->GetESDEvent(), fESDTrack->GetID());
  } else {
    fESDTPCOnlyTrack = fESDTrack;
    fFilterMap = 96; // to mimic the filterbit of global tracks
  }
}
void AliFemtoDreamTrack::SetESDTrackingInformation(const bool TPCOnlyTrack) {
  fESDStatus = fESDTrack->GetStatus();
  this->ApplyESDtoAODFilter(TPCOnlyTrack);
  if (!fESDTPCOnlyTrack) {
    this->fIsSet = false;
    return;
  }
  double pos[3] = { 0. };
  double covTr[21] = { 0. };
  //  Double_t pid[10]={0.};
  double p[3] = { 0. };
  double pDCA[3] = { 0. };  // momentum at DCA
  double rDCA[3] = { 0. };  // position at DCA
  float dDCA[2] = { 0. };    // DCA to the vertex d and z
  float cDCA[3] = { 0. };    // covariance of impact parameters
  // Loop over the ESD trcks and pick out the tracks passing TPC only cuts
  const AliESDVertex *vtxSPD = fESDTrack->GetESDEvent()->GetPrimaryVertexSPD();
//  const AliESDVertex *vtx = fESDTrack->GetESDEvent()->GetPrimaryVertex();

  if (fESDTPCOnlyTrack->Pt() > 0) {
    AliExternalTrackParam exParam;
    Bool_t relate = false;
    relate = fESDTPCOnlyTrack->RelateToVertexTPC(
        vtxSPD, fESDTrack->GetESDEvent()->GetMagneticField(), 1e30, &exParam);
    if (!relate) {
      this->fIsSet = false;
      return;
    }
    // fetch the track parameters at the DCA (unconstraint)
    if (fESDTPCOnlyTrack->GetTPCInnerParam()) {
      fESDTPCOnlyTrack->GetTPCInnerParam()->GetPxPyPz(pDCA);
      fESDTPCOnlyTrack->GetTPCInnerParam()->GetXYZ(rDCA);
    }
    // get the DCA to the vertex:
    fESDTrack->GetImpactParameters(dDCA, cDCA);
    // set the constrained parameters to the track
    fESDTPCOnlyTrack->Set(exParam.GetX(), exParam.GetAlpha(),
                          exParam.GetParameter(), exParam.GetCovariance());
    fESDTPCOnlyTrack->GetPxPyPz(p);
    fESDTPCOnlyTrack->GetXYZ(pos);
    fESDTPCOnlyTrack->GetCovarianceXYZPxPyPz(covTr);
  } else {
    this->fIsSet = false;
    return;
  }
  this->SetEta(fESDTPCOnlyTrack->Eta());
  this->SetPhi(fESDTPCOnlyTrack->Phi());
  this->SetTheta(fESDTPCOnlyTrack->Theta());
  this->SetCharge(fESDTPCOnlyTrack->Charge());
  this->SetMomTPC(fESDTrack->GetTPCmomentum());
  this->SetMomentum(0, p[0], p[1], p[2]);
  this->SetPt(fESDTPCOnlyTrack->Pt());

  this->fdcaXY = dDCA[0];
  this->fdcaZ = dDCA[1];
  fESDnClusterITS = fESDTrack->GetITSclusters(0);
  if (fESDnClusterITS != 0) {
    fChi2ITS = fESDTrack->GetITSchi2() / Float_t(fESDnClusterITS);
  }

  fESDnClusterTPC = fESDTPCOnlyTrack->GetTPCclusters(0);
  fTPCCrossedRows = fESDTPCOnlyTrack->GetTPCCrossedRows();
  if (fESDTPCOnlyTrack->GetTPCNclsF() > 0) {
    fRatioCR = fTPCCrossedRows / fESDTPCOnlyTrack->GetTPCNclsF();
  }
  fTPCClsS = fESDTPCOnlyTrack->GetTPCnclsS();
  if (fESDnClusterTPC != 0) {
    fChi2TPC = fESDTPCOnlyTrack->GetTPCchi2() / Float_t(fESDnClusterTPC);
    fTPCClsSRatio = Float_t(fTPCClsS) / Float_t(fESDnClusterTPC);
  }
  fNClsTPC = fESDTPCOnlyTrack->GetTPCNcls();
  if (fNClsTPC > 5) {
    fChi2 = fESDTPCOnlyTrack->GetTPCchi2() / Float_t(fNClsTPC - 5);
  } else {
    fChi2 = -1.;
  }
  //loop over the 6 ITS Layrs and check for a hit!
  for (int i = 0; i < 6; ++i) {
    fITSHit.push_back(fESDTrack->HasPointOnITSLayer(i));
    if (i == 2 || i == 3) continue;
    if (fESDTrack->HasPointOnITSLayer(i)) {
      this->fHasITSHit = true;
    }
  }
  if (fESDTrack->IsOn(AliESDtrack::kTPCrefit)) {
    fTPCRefit = true;
  }
  if (fESDTrack->GetTOFBunchCrossing() == 0) {
    this->fTOFTiming = true;
  } else {
    this->fTOFTiming = false;
  }
  const TBits sharedMap = fESDTrack->GetTPCSharedMap();
  if ((sharedMap.CountBits()) >= 1) {
    // Bad Track, has too many shared clusters!
    this->fnoSharedClst = false;
  } else {
    this->fnoSharedClst = true;
  }
  for (int i = 0; i < 6; ++i) {
    fSharedClsITSLayer.push_back(fESDTrack->HasSharedPointOnITSLayer(i));
    if (fSharedClsITSLayer.at(i)) {
      fHasSharedClsITSLayer = true;
    }
  }
  SetPhiAtRadii(fESDTrack->GetESDEvent()->GetMagneticField());
}


//_____________________________________________________________________________________
void AliFemtoDreamTrack::SetESDTrackingInformationOmega() {
  fESDStatus = fESDTrack->GetStatus();
    //Get primary vertex
    Double_t PrimVtx[3];
    fESDTrack->GetESDEvent()->GetPrimaryVertex()->GetXYZ(PrimVtx);
    //get tpc mom
    this->SetMomTPC(fESDTrack->GetTPCmomentum());
    //REQUIRE EVERYTHING TO HAVE TPCMOM > 50MeV
    if(fESDTrack->GetTPCmomentum()<.050){
      this->fIsSet = false;
      return;
    }

    //fill dca  //CHECK!!!
    float cDCA[3] = { 0. };    // covariance of impact parameters
    float dDCA[2] = { 0. };    // DCA to the vertex d and z
    fESDTrack->GetImpactParameters(dDCA, cDCA);
    this->fdcaXY = dDCA[0];
    this->fdcaZ = dDCA[1];

    //fill momentum. This will be overwritten later in the cascade setting.
    double p[3] = { 0. };
    fESDTrack->GetPxPyPz(p);
    this->SetMomentum(0, p[0], p[1], p[2]);
    this->SetPt(fESDTrack->Pt());

    //fill Eta etc:
    this->SetEta(fESDTrack->Eta());
    this->SetPhi(fESDTrack->Phi());
    this->SetTheta(fESDTrack->Theta());
    this->SetCharge(fESDTrack->Charge());

    //fill #TPCclusters
    fESDnClusterITS = fESDTrack->GetITSclusters(0);
    if (fESDnClusterITS != 0) {
      fChi2ITS = fESDTrack->GetITSchi2() / Float_t(fESDnClusterITS);
    }
    fESDnClusterTPC = fESDTrack->GetTPCclusters(0);
    fTPCCrossedRows = fESDTrack->GetTPCCrossedRows();
    if (fESDTrack->GetTPCNclsF() > 0) {
      fRatioCR = fTPCCrossedRows / fESDTrack->GetTPCNclsF();
    }
    fTPCClsS = fESDTrack->GetTPCnclsS();
    if (fESDnClusterTPC != 0) {
      fChi2TPC = fESDTrack->GetTPCchi2() / Float_t(fESDnClusterTPC);
      fTPCClsSRatio = Float_t(fTPCClsS) / Float_t(fESDnClusterTPC);
    }
    fNClsTPC = fESDTrack->GetTPCNcls();
    if (fNClsTPC > 5) {
      fChi2 = fESDTrack->GetTPCchi2() / Float_t(fNClsTPC - 5);
    } else {
      fChi2 = -1.;
    }
    //loop over the 6 ITS Layrs and check for a hit!
    for (int i = 0; i < 6; ++i) {
      fITSHit.push_back(fESDTrack->HasPointOnITSLayer(i));
      if (i == 2 || i == 3) continue;
      if (fESDTrack->HasPointOnITSLayer(i)) {
        this->fHasITSHit = true;
      }
    }
    if (fESDTrack->IsOn(AliESDtrack::kTPCrefit)) {
      fTPCRefit = true;
    }
    if (fESDTrack->GetTOFBunchCrossing() == 0) {
      this->fTOFTiming = true;
    } else {
      this->fTOFTiming = false;
    }
    const TBits sharedMap = fESDTrack->GetTPCSharedMap();
    if ((sharedMap.CountBits()) >= 1) {
      // Bad Track, has too many shared clusters!
      this->fnoSharedClst = false;
    } else {
      this->fnoSharedClst = true;
    }
    for (int i = 0; i < 6; ++i) {
      fSharedClsITSLayer.push_back(fESDTrack->HasSharedPointOnITSLayer(i));
      if (fSharedClsITSLayer.at(i)) {
        fHasSharedClsITSLayer = true;
      }
    }
    SetPhiAtRadii(fESDTrack->GetESDEvent()->GetMagneticField());

    //not to be filled so far:
    //fFilterMap = 0;
    //fMCP.SetXYZ(0, 0, 0);
    //fOrigin = AliFemtoDreamBasePart::kUnknown;
    //fMCPDGCode = 0;
    //fPDGMotherWeak = 0;
    //fdcaXYProp = -99;
    //fdcaZProp = -99;
    //fMCTheta.clear();
    //fPhiAtRadius.clear();
    //fXYZAtRadius.clear();
    //fMCPhi.clear();
    //fCPA = 0;
}


//_______________________________________________
void AliFemtoDreamTrack::SetESDPIDInformation() {
  AliPID::EParticleType particleID[6] = { AliPID::kElectron, AliPID::kMuon,
      AliPID::kPion, AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron };
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man
        ->GetInputEventHandler());
    if (inputHandler) {
      fPIDResponse = inputHandler->GetPIDResponse();
      if (!fPIDResponse) {
        AliFatal("No PID Response, did you run your PID Task?");
      }
    } else {
      AliFatal("No Input Handler");
    }
  } else {
    AliFatal("No PID Response");
  }
  AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(
      AliPIDResponse::kITS, fESDTrack);
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(
      AliPIDResponse::kTPC, fESDTrack);
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(
      AliPIDResponse::kTOF, fESDTrack);
  this->fstatusITS = statusITS;
  this->fstatusTPC = statusTPC;
  this->fstatusTOF = statusTOF;
  this->fdEdxITS = fESDTrack->GetITSsignal();
  this->fdEdxTPC = fESDTrack->GetTPCsignal();
  this->fbetaTOF = GetBeta(fESDTrack);
  for (int i = 0; i < 6; ++i) {
    if (statusITS == AliPIDResponse::kDetPidOk) {
      (this->fnSigmaITS)[i] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS,
                                                           fESDTrack,
                                                           particleID[i]);
    } else {
      (this->fnSigmaITS)[i] = -999.;
    }
    if (statusTPC == AliPIDResponse::kDetPidOk) {
      (this->fnSigmaTPC)[i] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                           fESDTrack,
                                                           particleID[i]);
    } else {
      (this->fnSigmaTPC)[i] = -999.;
    }
    if (statusTOF == AliPIDResponse::kDetPidOk) {
      (this->fnSigmaTOF)[i] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF,
                                                           fESDTrack,
                                                           particleID[i]);
    } else {
      (this->fnSigmaTOF)[i] = -999.;
    }
  }
}

//_______________________________________________
void AliFemtoDreamTrack::SetVInformation(AliVEvent *event) {
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(fVTrack);
  AliNanoAODTrack *globalNanoTrack =
      dynamic_cast<AliNanoAODTrack *>(fVGlobalTrack);
  this->SetEta(fVTrack->Eta());
  this->SetPhi(fVTrack->Phi());
  this->SetTheta(fVTrack->Theta());
  this->SetCharge(fVTrack->Charge());
  this->SetMomentum(0, fVTrack->Px(), fVTrack->Py(), fVTrack->Pz());
  this->SetPt(fVTrack->Pt());

  // loop over the 6 ITS Layrs and check for a hit!
  for (int i = 0; i < 6; ++i) {
    fITSHit.push_back(fVGlobalTrack->HasPointOnITSLayer(i));
    if (i == 2 || i == 3) continue;
    if (fVGlobalTrack->HasPointOnITSLayer(i)) {
      this->fHasITSHit = true;
    }
  }
//  if (fVTrack->IsOn(AliVTrack::kTPCrefit)) {
//    //doesn't seem to work for nanos.
//    fTPCRefit = true;
//  }
  if (fVGlobalTrack->GetTOFBunchCrossing() == 0) {
    this->fTOFTiming = true;
  } else {
    this->fTOFTiming = false;
  }
  this->fNClsTPC = fVTrack->GetTPCNcls();
  this->fTPCCrossedRows = nanoTrack->GetTPCNCrossedRows();
  const float findable = nanoTrack->GetTPCNclsF();
  this->fRatioCR = (findable > 0) ? fTPCCrossedRows / findable : 0;
  this->fTPCClsS = nanoTrack->GetTPCnclsS();
  this->fChi2 = nanoTrack->Chi2perNDF();
  this->fFilterMap = nanoTrack->GetFilterMap();
  this->SetMomTPC(globalNanoTrack->GetTPCmomentum());
  this->fdcaXY = nanoTrack->DCA();
  this->fdcaZ = nanoTrack->ZAtDCA();
  globalNanoTrack->GetImpactParameters(fdcaXYProp, fdcaZProp);
  SetPhiAtRadii(event->GetMagneticField());

  // PID stuff
  static Bool_t bPIDAvailable = AliNanoAODTrack::InitPIDIndex();
  if (!bPIDAvailable) {
    this->fIsSet = false;
    return;
  }

//  this->fdEdxTPC = nanoTrack->GetTPCsignal();
  this->fdEdxTPC =0;
  AliPID::EParticleType particleID[6] = {AliPID::kElectron, AliPID::kMuon,
                                         AliPID::kPion, AliPID::kKaon,
                                         AliPID::kProton, AliPID::kDeuteron};

  this->fstatusTPC =
      (std::abs(globalNanoTrack->GetVar(AliNanoAODTrack::GetPIDIndex(
                    AliNanoAODTrack::kSigmaTPC, AliPID::kPion)) +
                999.f) > 1E-6)
          ? AliPIDResponse::kDetPidOk
          : AliPIDResponse::kDetNoSignal;
  this->fstatusTOF =
      (std::abs(globalNanoTrack->GetVar(AliNanoAODTrack::GetPIDIndex(
                    AliNanoAODTrack::kSigmaTOF, AliPID::kPion)) +
                999.f) > 1E-6)
          ? AliPIDResponse::kDetPidOk
          : AliPIDResponse::kDetNoSignal;

  for (int i = 0; i < 6; ++i) {
    this->fnSigmaTPC[i] = globalNanoTrack->GetVar(AliNanoAODTrack::GetPIDIndex(
        AliNanoAODTrack::kSigmaTPC, particleID[i]));
    this->fnSigmaTOF[i] = globalNanoTrack->GetVar(AliNanoAODTrack::GetPIDIndex(
        AliNanoAODTrack::kSigmaTOF, particleID[i]));
  }

  if (fTPCClsS > 0) {
    // Bad Track, has too many shared clusters!
    this->fnoSharedClst = false;
  } else {
    this->fnoSharedClst = true;
  }

  if (fPIDResponse) {
    this->fbetaTOF = GetBeta(globalNanoTrack);
  }

  // TODO
  //   For the moment we don't need ITS PID
  //  this->fstatusITS = statusITS;
  //          this->fnSigmaITS)[i] =
  //    AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaITS,
  //    particleID[i]);
  //  for (int i = 0; i < 6; ++i) {
  //    fSharedClsITSLayer.push_back(fAODTrack->HasSharedPointOnITSLayer(i));
  //    if (track->HasSharedPointOnITSLayer(i)) {
  //      fHasSharedClsITSLayer = true;
  //    }
  //  }
  return;
}


void AliFemtoDreamTrack::SetAODTrackingInformation() {
  this->fFilterMap = fAODTrack->GetFilterMap();
  this->SetEta(fAODTrack->Eta());
  this->SetPhi(fAODTrack->Phi());
  this->SetTheta(fAODTrack->Theta());
  this->SetCharge(fAODTrack->Charge());
  this->SetMomentum(0, fAODTrack->Px(), fAODTrack->Py(), fAODTrack->Pz());
  this->SetMomTPC(fAODGlobalTrack->GetTPCmomentum());
  this->SetPt(fAODTrack->Pt());
  this->fdcaXY = fAODTrack->DCA();
  this->fdcaZ = fAODTrack->ZAtDCA();
  this->fChi2 = fAODTrack->Chi2perNDF();
  fAODGlobalTrack->GetImpactParameters(fdcaXYProp,fdcaZProp);
//  AliAODTrack copy(*fAODGlobalTrack);
//  fAODGlobalTrack->GetPosition(pos);
//  if (pos[0] * pos[0] + pos[1] * pos[1] <= 3. * 3.
//      && copy.PropagateToDCA(copy.GetAODEvent()->GetPrimaryVertex(),
//                             copy.GetAODEvent()->GetMagneticField(), 10,
//                             dcaVals, covar)) {
//    this->fdcaXYProp = dcaVals[0];
//    this->fdcaZProp = dcaVals[1];
//  } else {
//    this->fdcaXYProp = -99;
//    this->fdcaZProp = -99;
//  }
  //loop over the 6 ITS Layrs and check for a hit!
  for (int i = 0; i < 6; ++i) {
    fITSHit.push_back(fAODGlobalTrack->HasPointOnITSLayer(i));
    if (fAODGlobalTrack->HasPointOnITSLayer(i)) {
      if (i == 2 || i == 3) continue;
      this->fHasITSHit = true;
    }
  }
  if (fAODTrack->IsOn(AliAODTrack::kTPCrefit)) {
    fTPCRefit = true;
  } else {
    fTPCRefit = false;
  }
  if (fAODGlobalTrack->GetTOFBunchCrossing() == 0) {
    this->fTOFTiming = true;
  } else {
    this->fTOFTiming = false;
  }

  this->fNClsTPC = fAODTrack->GetTPCNcls();
  //This method was inherited from H. Beck analysis
  // In the documents
  // https://alisoft.cern.ch/AliRoot/trunk/TPC/doc/Definitions/Definitions.pdf
  // TPC people describe the cut strategy for the TPC. It is explicitly
  // stated that a cut on the number of crossed rows and a cut on the
  // number of crossed rows over findable clusters is recommended to
  // remove fakes. In the pdf a cut value of .83 on the ratio
  // is stated, no value for the number of crossed rows. Looking at the
  // AliESDTrackCuts.cxx one sees that exactly this cut is used with
  // 0.8 on the ratio and 70 on the crossed rows.

  // Checked the filter task and AliAODTrack and AliESDtrack and
  // AliESDtrackCuts and the Definitions.pdf:
  // The function to get the findable clusters is GetTPCNclsF()

  // For the number fo crossed rows for ESD tracks, the function
  // GetTPCCrossedRows() usually is used. Looking at the AliESDtrack.cxx
  // one sees that it's just an alias (with additional caching) for
  // GetTPCClusterInfo(2, 1); The identical function exists in the
  // AliAODTrack.cxx
  this->fTPCCrossedRows = fAODTrack->GetTPCClusterInfo(2, 1);
  if (!(fAODTrack->GetTPCNclsF() > 0)) {
    this->fRatioCR = 0.;
  } else {
    this->fRatioCR = fTPCCrossedRows / float(fAODTrack->GetTPCNclsF());
  }
  const TBits sharedMap = fAODTrack->GetTPCSharedMap();
  if ((sharedMap.CountBits()) >= 1) {
    // Bad Track, has too many shared clusters!
    this->fnoSharedClst = false;
  } else {
    this->fnoSharedClst = true;
  }
  for (int i = 0; i < 6; ++i) {
    fSharedClsITSLayer.push_back(fAODTrack->HasSharedPointOnITSLayer(i));
    if (fAODTrack->HasSharedPointOnITSLayer(i)) {
      fHasSharedClsITSLayer = true;
    }
  }
  this->fTPCClsS = fAODTrack->GetTPCnclsS();
  SetPhiAtRadii(fAODTrack->GetAODEvent()->GetMagneticField());
}
void AliFemtoDreamTrack::SetPhiAtRadii(const float bfield) {
  float TPCradii[9] = { 85., 105., 125., 145., 165., 185., 205., 225., 245. };
  float phi0 = GetPhi().at(0);
  float pt = GetPt();
  float chg = GetCharge().at(0);
  std::vector<float> phiatRadius;
  for (int radius = 0; radius < 9; radius++) {
    phiatRadius.push_back(
        phi0
            - TMath::ASin(
                0.1 * chg * bfield * 0.3 * TPCradii[radius] * 0.01
                    / (2. * pt)));
  }
  fPhiAtRadius.push_back(phiatRadius);
  return;
}

void AliFemtoDreamTrack::SetGlobalCoordAtRadii(const float bfield) {
  float TPCradii[9] = { 85., 105., 125., 145., 165., 185., 205., 225., 245. };
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(fAODGlobalTrack);
  for (int iRad = 0; iRad < 9; ++iRad) {
    double posBuffer[3] = { 0., 0., 0. };
    bool good = etp.GetXYZatR(TPCradii[iRad], bfield, posBuffer, nullptr);
    if (good) {
      fXYZAtRadius.push_back(TVector3(posBuffer));
    } else {
      fXYZAtRadius.push_back(TVector3(999, 999, 999));
    }
  }

}

void AliFemtoDreamTrack::SetAODPIDInformation() {
  AliPID::EParticleType particleID[6] = { AliPID::kElectron, AliPID::kMuon,
      AliPID::kPion, AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron };
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man
        ->GetInputEventHandler());
    if (inputHandler) {
      fPIDResponse = inputHandler->GetPIDResponse();
      if (!fPIDResponse) {
        AliFatal("No PID Response, did you run your PID Task?");
      }
    } else {
      AliFatal("No Input Handler");
    }
  } else {
    AliFatal("No PID Response");
  }
  AliPIDResponse::EDetPidStatus statusITS = fPIDResponse->CheckPIDStatus(
      AliPIDResponse::kITS, fAODGlobalTrack);
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(
      AliPIDResponse::kTPC, fAODGlobalTrack);
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(
      AliPIDResponse::kTOF, fAODGlobalTrack);
  this->fstatusITS = statusITS;
  this->fstatusTPC = statusTPC;
  this->fstatusTOF = statusTOF;
  this->fdEdxITS= fAODGlobalTrack->GetITSsignal();
  this->fdEdxTPC = fAODGlobalTrack->GetTPCsignal();
  this->fbetaTOF = GetBeta(fAODGlobalTrack);
  for (int i = 0; i < 6; ++i) {
    if (statusITS == AliPIDResponse::kDetPidOk) {
      (this->fnSigmaITS)[i] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS,
                                                           fAODGlobalTrack,
                                                           particleID[i]);
    } else {
      (this->fnSigmaITS)[i] = -999.;
    }
    if (statusTPC == AliPIDResponse::kDetPidOk) {
      (this->fnSigmaTPC)[i] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                           fAODGlobalTrack,
                                                           particleID[i]);
    } else {
      (this->fnSigmaTPC)[i] = -999.;
    }
    if (statusTOF == AliPIDResponse::kDetPidOk) {
      (this->fnSigmaTOF)[i] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF,
                                                           fAODGlobalTrack,
                                                           particleID[i]);
    } else {
      (this->fnSigmaTOF)[i] = -999.;
    }
  }
}

void AliFemtoDreamTrack::SetMCInformation() {
  //Set the phi at radii at the TPC information
  //      SetPhiStar(track,fphiAtRadius);
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(fAODGlobalTrack
      ->GetAODEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcarray) {
    AliError("SPTrack: MC Array not found");
  }
  SetMCInformation(mcarray, fAODGlobalTrack->GetLabel());
}

void AliFemtoDreamTrack::SetMCInformation(AliVEvent *event) {
  TClonesArray *mcarray = dynamic_cast<TClonesArray *>(event->FindListObject(
      "mcparticles"));
  if (!mcarray) {
    AliError("SPTrack: MC Array not found");
  }
  SetMCInformation(mcarray, fVGlobalTrack->GetLabel());
}

void AliFemtoDreamTrack::SetMCInformation(TClonesArray* mcarray, int label) {
  if (label > 0) {
    AliAODMCParticle * mcPart = (AliAODMCParticle*) mcarray->At(label);
    this->SetID(label);
    if (!(mcPart)) {
      this->fIsSet = false;
    } else {
      this->SetMCPhi(mcPart->Phi());
      this->SetMCTheta(mcPart->Theta());
      this->SetMCPDGCode(mcPart->PdgCode());
      this->SetMCPt(mcPart->Pt());
      this->SetMCMomentum(mcPart->Px(), mcPart->Py(), mcPart->Pz());

      //check for secondary and set origin and mother
      if (mcPart->IsPhysicalPrimary() && !mcPart->IsSecondaryFromWeakDecay()) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);
      } else if (mcPart->IsSecondaryFromWeakDecay()
          && !mcPart->IsSecondaryFromMaterial()) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kWeak);
        this->SetPDGMotherWeak(
            ((AliAODMCParticle*) mcarray->At(mcPart->GetMother()))->PdgCode());
      } else if (mcPart->IsSecondaryFromMaterial()) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kMaterial);
      } else {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kUnknown);
      }
      int motherID = mcPart->GetMother();
      int lastMother = motherID;
      AliAODMCParticle *mcMother = nullptr;
      while (motherID != -1) {
        lastMother = motherID;
        mcMother = (AliAODMCParticle *) mcarray->At(motherID);
        motherID = mcMother->GetMother();
      }
      if (lastMother != -1) {
        mcMother = (AliAODMCParticle *) mcarray->At(lastMother);
      }
      if (mcMother) {
        this->SetMotherPDG(mcMother->GetPdgCode());
        this->SetMotherID(lastMother);
      }
    }
  } else {
    this->fIsSet = false;  //if we don't have MC Information, don't use that track
  }
}

void AliFemtoDreamTrack::SetMCInformation(AliMCEvent *mcEvent) {
  if (!mcEvent) {
    AliError("SPTrack: MC Event not found");
  }
  if (fESDTrack->GetLabel() > 0) {
    AliMCParticle * mcPart = (AliMCParticle*) mcEvent->GetTrack(
        fESDTrack->GetLabel());
    if (!(mcPart)) {
      this->fIsSet = false;
    } else {
      this->SetMCPhi(mcPart->Phi());
      this->SetMCTheta(mcPart->Theta());
      this->SetMCPDGCode(mcPart->PdgCode());
      this->SetMCPt(mcPart->Pt());
      this->SetMCMomentum(mcPart->Px(), mcPart->Py(), mcPart->Pz());

      //check for secondary and set origin and mother
      if (mcPart->IsPhysicalPrimary() && !mcPart->IsSecondaryFromWeakDecay()) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);
      } else if (mcPart->IsSecondaryFromWeakDecay()
          && !mcPart->IsSecondaryFromMaterial()) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kWeak);
        this->SetPDGMotherWeak(
            ((AliMCParticle*) mcEvent->GetTrack(mcPart->GetMother()))->PdgCode());
      } else if (mcPart->IsSecondaryFromMaterial()) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kMaterial);
      } else {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kUnknown);
      }
    }
  } else {
    this->fIsSet = false;  //if we don't have MC Information, don't use that track
  }
}

float AliFemtoDreamTrack::GetBeta(AliNanoAODTrack *track) const {
  static float c = 2.99792457999999984e-02;
  if (!fPIDResponse) {
    return -0.075f;
  }

  const float len = track->GetIntegratedLength();
  if (!(track->HasTOFpid() && (len > 350.))) {
    return -0.05f;
  }
  const float tim = track->GetTOFsignal()
      - fPIDResponse->GetTOFResponse().GetStartTime(track->GetTPCmomentum());
  return len / (tim * c);
}

float AliFemtoDreamTrack::GetBeta(AliAODTrack *track) const {
  static float c = 2.99792457999999984e-02;
  const float len = track->GetIntegratedLength();
  if (!(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track)
      && (len > 350.))) {
    return -0.05f;
  }
  const float tim = track->GetTOFsignal()
      - fPIDResponse->GetTOFResponse().GetStartTime(track->GetTPCmomentum());
  return len / (tim * c);
}

float AliFemtoDreamTrack::GetBeta(AliESDtrack *track) const {
  static float c = 2.99792457999999984e-02;
  const float len = track->GetIntegratedLength();
  if (!(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track)
      && (len > 350.))) {
    return -0.05f;
  }
  const float tim = track->GetTOFsignal()
      - fPIDResponse->GetTOFResponse().GetStartTime(track->GetTPCmomentum());
  return len / (tim * c);
}

bool AliFemtoDreamTrack::CheckGlobalTrack(const Int_t TrackID) {
  //This method was inherited from H. Beck analysis
  //Checks if to the corresponding track a global track exists
  //This is especially useful if one has TPC only tracks and needs the PID information
  bool isGlobal = true;
  if (TMath::Abs(TrackID) < fTrackBufferSize) {
    if (!(fGTI[-TrackID - 1])) {
      isGlobal = false;
    }
  }
  return isGlobal;
}
bool AliFemtoDreamTrack::CheckGlobalVTrack(const Int_t TrackID) {
  bool isGlobal = true;
  if (TMath::Abs(TrackID) < fTrackBufferSize) {
    if (!(fVGTI[-TrackID - 1])) {
      isGlobal = false;
    }
  }
  return isGlobal;
}
void AliFemtoDreamTrack::Reset() {
  if (!fIsReset) {
    fstatusITS = AliPIDResponse::kDetNoParams;
    fstatusTPC = AliPIDResponse::kDetNoParams;
    fstatusTOF = AliPIDResponse::kDetNoParams;
    fFilterMap = 0;
    fdcaXY = -99;
    fdcaZ = -99;
    fdcaXYProp = -99;
    fdcaZProp = -99;
    fNClsTPC = 0;
    fTPCCrossedRows = 0;
    fRatioCR = 0;
    fnoSharedClst = 0;
    fTPCClsS = 0;
    fTPCClsSRatio = 0;
    fChi2 = 0;
    fChi2TPC = 0;
    fChi2ITS = 0;
    fSharedClsITSLayer.clear();
    fHasSharedClsITSLayer = false;
    fdEdxITS = -999;
    fdEdxTPC = -999;
    fbetaTOF = 1.1;
    fHasITSHit = false;
    fITSHit.clear();
    fTOFTiming = false;
    fTPCRefit = false;
    for (int i = 0; i < 5; ++i) {
      fnSigmaITS[i] = 99;
      fnSigmaTPC[i] = 99;
      fnSigmaTOF[i] = 99;
    }
    fESDStatus = 0;
    fESDnClusterITS = 0;
    fESDnClusterTPC = 0;
    GetMomentum(0).SetXYZ(0, 0, 0);
    fMCP.SetXYZ(0, 0, 0);
    fPt = 0;
    fMCPt = 0;
    fMCPt = 0;
    fP_TPC = 0;
    fEta.clear();
    fTheta.clear();
    fMCTheta.clear();
    fPhi.clear();
    fPhiAtRadius.clear();
    fXYZAtRadius.clear();
    fMCPhi.clear();
    fIDTracks.clear();
    fCharge.clear();
    fCPA = 0;
    fOrigin = AliFemtoDreamBasePart::kUnknown;
    //we don't want to reset the fPDGCode
    fMCPDGCode = 0;
    fPDGMotherWeak = 0;
    //we don't want to reset isMC
    fUse = false;
    fIsSet = true;
    fIsReset = true;
  }
}
