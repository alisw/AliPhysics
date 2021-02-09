/*
 * AliFemtoDreamEvent.cxx
 *
 *  Created on: 22 Nov 2017
 *      Author: bernhardhohlweger
 */

#include "AliFemtoDreamEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliVEvent.h"
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "AliNanoAODTrack.h"

ClassImp(AliFemtoDreamEvent)
AliFemtoDreamEvent::AliFemtoDreamEvent()
    : fUtils(nullptr),
      fEvtCuts(nullptr),
      fuseAliEvtCuts(false),
      fEvtCutList(nullptr),
      fxVtx(0),
      fyVtx(0),
      fzVtx(0),
      fzVtxTracks(0),
      fzVtxSPD(0),
      fBField(-99),
      fSPDMult(0),
      fNSPDClusterLy0(0),
      fNSPDClusterLy1(0),
      fRefMult08(0),
      fV0AMult(0),
      fV0CMult(0),
      fV0ATime(0),
      fV0CTime(0),
      fV0MCentrality(0),
      fnContrib(0),
      fPassAliEvtSelection(false),
      fisPileUp(false),
      fHasVertex(false),
      fHasMagField(false),
      fisSelected(false),
      fEstimator(AliFemtoDreamEvent::kRef08),
      fspher(0),
      fLowPtSpherCalc(),
      fsphero(0),
      fcalcsphero(false){

}

AliFemtoDreamEvent::AliFemtoDreamEvent(bool mvPileUp, bool EvtCutQA,
                                       UInt_t trigger, bool useEvtCuts, float LowPtSpherCalc)
    : fUtils(new AliAnalysisUtils()),
      fEvtCuts(nullptr),
      fuseAliEvtCuts(useEvtCuts),
      fEvtCutList(nullptr),
      fxVtx(0),
      fyVtx(0),
      fzVtx(0),
      fzVtxTracks(0),
      fzVtxSPD(0),
      fBField(-99),
      fSPDMult(0),
      fNSPDClusterLy0(0),
      fNSPDClusterLy1(0),
      fRefMult08(0),
      fV0AMult(0),
      fV0CMult(0),
      fV0ATime(0),
      fV0CTime(0),
      fV0MCentrality(0),
      fnContrib(0),
      fPassAliEvtSelection(false),
      fisPileUp(false),
      fHasVertex(false),
      fHasMagField(false),
      fisSelected(false),
      fEstimator(kRef08),
      fspher(0),
      fLowPtSpherCalc(LowPtSpherCalc),
      fsphero(0),
      fcalcsphero(false) {
  if (fuseAliEvtCuts) {
    fEvtCuts = new AliEventCuts();
    std::cout << "Setting up Event Cuts correspondingly for pp trigger: "
        << trigger << std::endl;
    if (trigger != AliVEvent::kINT7) {
      fEvtCuts->OverrideAutomaticTriggerSelection(trigger);
    }
  }
  if (mvPileUp) {
    //For pPb this is necessary according to DPG Processing status news
    //(week 29 April - 5 May 2017)
    fUtils->SetUseMVPlpSelection(true);
  } else {
    //Following the analysis in pp Run1 of O.Arnold
    fUtils->SetMinPlpContribSPD(3);
  }

  if (EvtCutQA&&fuseAliEvtCuts) {
    fEvtCutList = new TList();
    fEvtCutList->SetName("AliEventCuts");
    fEvtCutList->SetOwner(true);
    fEvtCuts->AddQAplotsToList(fEvtCutList);
  } else {
    fEvtCutList = nullptr;
  }
}


AliFemtoDreamEvent::~AliFemtoDreamEvent() {
  if (fEvtCuts) {
    delete fEvtCuts;
  }
//  if (fUtils) {
//    delete fUtils;
//  }
}

AliFemtoDreamEvent &AliFemtoDreamEvent::operator=(
    const AliFemtoDreamEvent &obj) {
  if (this == &obj) {
    return *this;
  }
  fEvtCutList = obj.fEvtCutList;
  fxVtx = obj.fxVtx;
  fyVtx = obj.fyVtx;
  fzVtx = obj.fzVtx;
  fzVtxTracks = obj.fzVtxTracks;
  fzVtxSPD = obj.fzVtxSPD;
  fBField = obj.fBField;
  fSPDMult = obj.fSPDMult;
  fNSPDClusterLy0 = obj.fNSPDClusterLy0;
  fNSPDClusterLy1 = obj.fNSPDClusterLy1;
  fRefMult08 = obj.fRefMult08;
  fV0AMult = obj.fV0AMult;
  fV0CMult = obj.fV0CMult;
  fV0ATime = obj.fV0ATime;
  fV0CTime = obj.fV0CTime;
  fV0MCentrality = obj.fV0MCentrality;
  fnContrib = obj.fnContrib;
  fPassAliEvtSelection = obj.fPassAliEvtSelection;
  fisPileUp = obj.fisPileUp;
  fHasVertex = obj.fHasVertex;
  fHasMagField = obj.fHasMagField;
  fisSelected = obj.fisSelected;
  fEstimator = obj.fEstimator;
  fspher = obj.fspher;
  fLowPtSpherCalc = obj.fLowPtSpherCalc;
  fsphero = obj.fsphero;
  fcalcsphero= obj.fcalcsphero;
  return (*this);
}

void AliFemtoDreamEvent::SetEvent(AliAODEvent *evt) {
  AliAODHeader *header = dynamic_cast<AliAODHeader*>(evt->GetHeader());

  AliAODVertex *vtx = evt->GetPrimaryVertex();
  AliAODVZERO *vZERO = evt->GetVZEROData();
  if (!vtx) {
    this->fHasVertex = false;
  } else {
    this->fHasVertex = true;
  }
  if (TMath::Abs(evt->GetMagneticField()) < 0.001) {
    this->fHasMagField = false;
  } else {
    this->fHasMagField = true;
    this->fBField = evt->GetMagneticField();
  }
  this->fnContrib = vtx->GetNContributors();
  this->fxVtx = vtx->GetX();
  this->fyVtx = vtx->GetY();
  this->fzVtx = vtx->GetZ();
  if (fUtils->IsPileUpEvent(evt)) {
    this->fisPileUp = true;
  } else {
    this->fisPileUp = false;
  }
  if (fuseAliEvtCuts) {
    if (fEvtCuts->AcceptEvent(evt)) {
      this->fPassAliEvtSelection = true;
    } else {
      this->fPassAliEvtSelection = false;
    }
  }
  this->fSPDMult = CalculateITSMultiplicity(evt);
  this->fNSPDClusterLy0 = evt->GetNumberOfITSClusters(0);
  this->fNSPDClusterLy1 = evt->GetNumberOfITSClusters(1);
  this->fV0AMult = vZERO->GetMTotV0A();
  this->fV0CMult = vZERO->GetMTotV0C();
  this->fV0ATime = vZERO->GetV0ATime();
  this->fV0CTime = vZERO->GetV0CTime();
  this->fspher = CalculateSphericityEvent(evt, this->fLowPtSpherCalc);
  if (fcalcsphero){
       this->fsphero = CalculateSpherocityEvent(evt);
  }
  fRefMult08 = header->GetRefMultiplicityComb08();
  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection *) evt->FindListObject("MultSelection");
  if (!MultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  } else {
    fV0MCentrality = MultSelection->GetMultiplicityPercentile("V0M");
  }
  return;
}

void AliFemtoDreamEvent::SetEvent(AliVEvent *evt) {
  AliNanoAODHeader* nanoHeader = dynamic_cast<AliNanoAODHeader*>(evt->GetHeader());
  const AliVVertex *vtx = evt->GetPrimaryVertex();
  if (!vtx) {
    this->fHasVertex = false;
  } else {
    this->fHasVertex = true;
    this->fnContrib = vtx->GetNContributors();
    this->fxVtx = vtx->GetX();
    this->fyVtx = vtx->GetY();
    this->fzVtx = vtx->GetZ();
  }
  if (TMath::Abs(evt->GetMagneticField()) < 0.001) {
    this->fHasMagField = false;
  } else {
    this->fHasMagField = true;
    this->fBField = evt->GetMagneticField();
  }

  // This we already know since the NanoAOD filtering was done
  this->fisPileUp = false;
  this->fPassAliEvtSelection = true;

  static const Int_t kRefMult =
      nanoHeader->GetVarIndex("MultSelection.RefMult08.Value");
  static const Int_t kSPDCluster =
      nanoHeader->GetVarIndex("MultSelection.SPDClusters.Value");
  static const Int_t kSPDTracklets =
      nanoHeader->GetVarIndex("MultSelection.SPDTracklets.Value");
  static const Int_t kV0A =
      nanoHeader->GetVarIndex("MultSelection.V0A.Value");
  static const Int_t kV0C =
      nanoHeader->GetVarIndex("MultSelection.V0C.Value");
  if (kSPDCluster != -1) {
    //dirty trick to have the combined plot properly
    this->fNSPDClusterLy0 = nanoHeader->GetVar(kSPDCluster);
    this->fNSPDClusterLy1 = fNSPDClusterLy0;
  }
  if (kRefMult != -1) {
    this->fRefMult08 = nanoHeader->GetVar(kRefMult);
  }
  this->fV0MCentrality = nanoHeader->GetCentr("V0M");
  if (kSPDTracklets!=-1) {
    this->fSPDMult = nanoHeader->GetVar(kSPDTracklets);
  }
  if (kV0A!=-1) {
    this->fV0AMult = nanoHeader->GetVar(kV0A);
  }
  if (kV0C!=-1) {
    this->fV0CMult = nanoHeader->GetVar(kV0C);
  }
  this->fspher = CalculateSphericityEvent(evt);

  if (fcalcsphero){
  this->fsphero = CalculateSpherocityEvent(evt);
  }
}

void AliFemtoDreamEvent::SetEvent(AliESDEvent *evt) {
  const AliVVertex *vtx = evt->GetPrimaryVertex();
  const AliVVertex* vtxSPD = evt->GetPrimaryVertexSPD();
  fzVtxSPD = (vtxSPD) ? vtxSPD->GetZ() : vtx->GetZ();
  const AliESDVertex* vtxTrk = evt->GetPrimaryVertexTracks();
  fzVtxTracks = (vtxTrk) ? vtxTrk->GetZ() : vtx->GetZ();
  AliESDVZERO *vZERO = evt->GetVZEROData();
  if (!vtx) {
    this->fHasVertex = false;
  } else {
    this->fHasVertex = true;
  }
  if (TMath::Abs(evt->GetMagneticField()) < 0.001) {
    this->fHasMagField = false;
  } else {
    this->fHasMagField = true;
    this->fBField = evt->GetMagneticField();
  }
  this->fnContrib = vtx->GetNContributors();
  this->fxVtx = vtx->GetX();
  this->fyVtx = vtx->GetY();
  this->fzVtx = vtx->GetZ();
  if (fUtils->IsPileUpEvent(evt)) {
    this->fisPileUp = true;
  } else {
    this->fisPileUp = false;
  }
  if (fEvtCuts->AcceptEvent(evt)) {
    this->fPassAliEvtSelection = true;
  } else {
    this->fPassAliEvtSelection = false;
  }
  //!to do: Check event multiplicity estimation!
  if (evt->GetMultiplicity()) {
    this->fSPDMult = evt->GetMultiplicity()->GetNumberOfTracklets();
    this->fNSPDClusterLy0 = evt->GetMultiplicity()->GetNumberOfITSClusters(0,
                                                                           0);
    this->fNSPDClusterLy1 = evt->GetMultiplicity()->GetNumberOfITSClusters(1,
                                                                           1);
  } else {
    this->fSPDMult = evt->GetNumberOfITSClusters(1);
    this->fNSPDClusterLy0 = 0;
    this->fNSPDClusterLy1 = 0;
  }
  this->fRefMult08 = AliESDtrackCuts::GetReferenceMultiplicity(
      evt, AliESDtrackCuts::kTrackletsITSTPC, 0.8, 0);
  this->fV0AMult = vZERO->GetMTotV0A();
  this->fV0CMult = vZERO->GetMTotV0C();
  this->fV0ATime = vZERO->GetV0ATime();
  this->fV0CTime = vZERO->GetV0CTime();
  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection *) evt->FindListObject("MultSelection");
  if (!MultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  } else {
    fV0MCentrality = MultSelection->GetMultiplicityPercentile("V0M");
  }
  return;
}
//rebuild
int AliFemtoDreamEvent::CalculateITSMultiplicity(AliAODEvent *evt) {
  AliAODTracklets* tracklets = evt->GetTracklets();
  int nTr = tracklets->GetNumberOfTracklets();
  int count = 0;
  for (int iTr = 0; iTr < nTr; iTr++) {
    float theta = tracklets->GetTheta(iTr);
    float eta = -TMath::Log(TMath::Tan(theta / 2.));
    if (TMath::Abs(eta) < 0.8) {
      count++;
    }
  }
  return count;
}

int AliFemtoDreamEvent::GetMultiplicity() {
  int mult = 0;
  switch (fEstimator) {
    case AliFemtoDreamEvent::kRef08:
      mult = GetRefMult08();
      break;
    case AliFemtoDreamEvent::kSPD:
      mult = GetSPDMult();
      break;
    case AliFemtoDreamEvent::kV0M:
      mult = GetV0MMult();
      break;
    case AliFemtoDreamEvent::kV0A:
      mult = GetV0AMult();
      break;
    case AliFemtoDreamEvent::kV0C:
      mult = GetV0CMult();
      break;
    default:
      AliFatal("Type Not implemented");
      break;
  }

  return mult;
}

double AliFemtoDreamEvent::CalculateSphericityEvent(AliAODEvent *evt, float lowerPtbound) {
//Initializing
  double ptTot = 0.;
  double s00 = 0.;  //elements of the sphericity matrix taken form EPJC72:2124
  double s01 = 0.;
  double s10 = 0.;
  double s11 = 0.;

  int numOfTracks = evt->GetNumberOfTracks();
  if (numOfTracks < 3)
    return -9999.;

  for (int iTrack = 0; iTrack < numOfTracks; iTrack++) {

    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(evt->GetTrack(iTrack));
    double pt = aodtrack->Pt();
    double eta = aodtrack->Eta();
    double px = aodtrack->Px();
    double py = aodtrack->Py();
    if(!aodtrack->TestFilterBit(96)) continue;
    if (TMath::Abs(pt) < lowerPtbound || TMath::Abs(eta) > 0.8) {
      continue;
    }

    ptTot += pt;

    s00 += px * px / pt;
    s01 += px * py / pt;
    s10 = s01;
    s11 += py * py / pt;
  }

  //normalize to total Pt to obtain a linear form:
  if (ptTot == 0.)
    return -9999.;
  s00 /= ptTot;
  s11 /= ptTot;
  s10 /= ptTot;

  //Calculate the trace of the sphericity matrix:
  double T = s00 + s11;
  //Calculate the determinant of the sphericity matrix:
  double D = s00 * s11 - s10 * s10;  //S10 = S01

  //Calculate the eigenvalues of the sphericity matrix:
  double lambda1 = 0.5 * (T + std::sqrt(T * T - 4. * D));
  double lambda2 = 0.5 * (T - std::sqrt(T * T - 4. * D));

  if ((lambda1 + lambda2) == 0.)
    return -9999.;

  double spt = -1.;

  if (lambda2 > lambda1) {
    spt = 2. * lambda1 / (lambda1 + lambda2);
  } else {
    spt = 2. * lambda2 / (lambda1 + lambda2);
  }

  return spt;
}


double AliFemtoDreamEvent::CalculateSphericityEvent(AliVEvent *evt) {
//Initializing
  double ptTot = 0.;
  double s00 = 0.;  //elements of the sphericity matrix taken form EPJC72:2124
  double s01 = 0.;
  double s10 = 0.;
  double s11 = 0.;

  int numOfTracks = evt->GetNumberOfTracks();
  if (numOfTracks < 3)
    return -9999.;

  for (int iTrack = 0; iTrack < numOfTracks; iTrack++) {

    AliNanoAODTrack *track = dynamic_cast<AliNanoAODTrack*>(evt->GetTrack(iTrack));

    double pt = track->Pt();
    double eta = track->Eta();
    double px = track->Px();
    double py = track->Py();
    if(!track->TestFilterBit(96)) continue;
    if (TMath::Abs(pt) < 0.5 || TMath::Abs(eta) > 0.8) {
      continue;
    }

    ptTot += pt;

    s00 += px * px / pt;
    s01 += px * py / pt;
    s10 = s01;
    s11 += py * py / pt;
  }

  //normalize to total Pt to obtain a linear form:
  if (ptTot == 0.)
    return -9999.;
  s00 /= ptTot;
  s11 /= ptTot;
  s10 /= ptTot;

  //Calculate the trace of the sphericity matrix:
  double T = s00 + s11;
  //Calculate the determinant of the sphericity matrix:
  double D = s00 * s11 - s10 * s10;  //S10 = S01

  //Calculate the eigenvalues of the sphericity matrix:
  double lambda1 = 0.5 * (T + std::sqrt(T * T - 4. * D));
  double lambda2 = 0.5 * (T - std::sqrt(T * T - 4. * D));

  if ((lambda1 + lambda2) == 0.)
    return -9999.;

  double spt = -1.;

  if (lambda2 > lambda1) {
    spt = 2. * lambda1 / (lambda1 + lambda2);
  } else {
    spt = 2. * lambda2 / (lambda1 + lambda2);
  }

  return spt;
}

double AliFemtoDreamEvent::CalculateSpherocityEvent(AliAODEvent *evt) {
  float pFull = 0.f;
  float Spherocity = 2.f;

  const float pi = TMath::Pi();

  float pTtot = 0.f;

  std::vector<float> pXVec;
  std::vector<float> pYVec;

  int numOfTracks = evt->GetNumberOfTracks();
  if (numOfTracks < 3)
    return -9999.;

  for (int iTrack = 0; iTrack < numOfTracks; iTrack++) {
    AliAODTrack *track = dynamic_cast<AliAODTrack *>(evt->GetTrack(iTrack));
    if (!track->TestFilterBit(96))
      continue;
    double pt = track->Pt();
    if (TMath::Abs(pt) < 0.5 || TMath::Abs(track->Eta()) > 0.8) {
      continue;
    }
    pTtot += pt;
    pXVec.push_back(track->Px());
    pYVec.push_back(track->Py());
  }

  if (pTtot == 0.f)
    return -9999.;

  const float OneOverPtTotal = 1.f / pTtot;

  float numerator = 0.f;
  float phiparam = 0.f;
  float nx = 0.f;
  float ny = 0.f;

  for (int i = 0; i < 360 / 0.1; ++i) {
    numerator = 0.f;
    phiparam = (pi * i * 0.1 / 180);  // parametrization of the angle
    nx = TMath::Cos(phiparam);  // x component of an unitary vector n
    ny = TMath::Sin(phiparam);  // y component of an unitary vector n

    for (size_t itTrack = 0; itTrack < pXVec.size(); ++itTrack) {
      numerator += TMath::Abs(ny * pXVec[itTrack] - nx * pYVec[itTrack]);  // product between p
      // proyection in XY plane and
      // the unitary vector
    }
    pFull = std::pow((numerator * OneOverPtTotal), 2);

    if (pFull < Spherocity)  // maximization of pFull
        {
      Spherocity = pFull;
    }
  }
  return ((Spherocity) * pi * pi) / 4.0;
}

double AliFemtoDreamEvent::CalculateSpherocityEvent(AliVEvent *evt) {
  float pFull = 0.f;
  float Spherocity = 2.f;

  const float pi = TMath::Pi();

  float pTtot = 0.f;

  int numOfTracks = evt->GetNumberOfTracks();
  if (numOfTracks < 3)
    return -9999.;

  std::vector<float> pXVec;
  std::vector<float> pYVec;

  for (int iTrack = 0; iTrack < numOfTracks; iTrack++) {
    AliNanoAODTrack *track = dynamic_cast<AliNanoAODTrack *>(evt->GetTrack(
        iTrack));
    if (!track->TestFilterBit(96))
      continue;
    double pt = track->Pt();
    if (TMath::Abs(pt) < 0.5 || TMath::Abs(track->Eta()) > 0.8) {
      continue;
    }
    pTtot += pt;
    pXVec.push_back(track->Px());
    pYVec.push_back(track->Py());
  }

  if (pTtot == 0.f) {
    return -9999.;
  }

  const float OneOverPtTotal = 1.f / pTtot;

  float numerator = 0.f;
  float phiparam = 0.f;
  float nx = 0.f;
  float ny = 0.f;

  for (int i = 0; i < 360 / 0.1; ++i) {
    numerator = 0.f;
    phiparam = (pi * i * 0.1 / 180);  // parametrization of the angle
    nx = std::cos(phiparam);  // x component of an unitary vector n
    ny = std::sin(phiparam);  // y component of an unitary vector n

    for (size_t itTrack = 0; itTrack < pXVec.size(); ++itTrack) {
      numerator += TMath::Abs(ny * pXVec[itTrack] - nx * pYVec[itTrack]);  // product between p
      // proyection in XY plane and
      // the unitary vector
    }
    pFull = std::pow((numerator * OneOverPtTotal), 2);

    if (pFull < Spherocity)  // maximization of pFull
        {
      Spherocity = pFull;
    }
  }
  return ((Spherocity) * pi * pi) / 4.0;
}
