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
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"

ClassImp(AliFemtoDreamEvent)
AliFemtoDreamEvent::AliFemtoDreamEvent()
    : fUtils(nullptr),
      fEvtCuts(nullptr),
      fEvtCutList(nullptr),
      fxVtx(0),
      fyVtx(0),
      fzVtx(0),
      fzVtxTracks(0),
      fzVtxSPD(0),
      fBField(-99),
      fSPDMult(0),
      fNSPDCluster(0),
      fRefMult08(0),
      fV0AMult(0),
      fV0CMult(0),
      fV0MCentrality(0),
      fnContrib(0),
      fPassAliEvtSelection(false),
      fisPileUp(false),
      fHasVertex(false),
      fHasMagField(false),
      fisSelected(false),
      fEstimator(AliFemtoDreamEvent::kRef08),
      fspher(0) {

}

AliFemtoDreamEvent::AliFemtoDreamEvent(bool mvPileUp, bool EvtCutQA,
                                       UInt_t trigger)
    : fUtils(new AliAnalysisUtils()),
      fEvtCuts(new AliEventCuts()),
      fEvtCutList(nullptr),
      fxVtx(0),
      fyVtx(0),
      fzVtx(0),
      fzVtxTracks(0),
      fzVtxSPD(0),
      fBField(-99),
      fSPDMult(0),
      fNSPDCluster(0),
      fRefMult08(0),
      fV0AMult(0),
      fV0CMult(0),
      fV0MCentrality(0),
      fnContrib(0),
      fPassAliEvtSelection(false),
      fisPileUp(false),
      fHasVertex(false),
      fHasMagField(false),
      fisSelected(false),
      fEstimator(kRef08),
      fspher(0) {
  if (mvPileUp) {
    //For pPb this is necessary according to DPG Processing status news
    //(week 29 April - 5 May 2017)
    fUtils->SetUseMVPlpSelection(true);
  } else {
    //Following the analysis in pp Run1 of O.Arnold
    fUtils->SetMinPlpContribSPD(3);
  }

  if (trigger != AliVEvent::kINT7) {
    fEvtCuts->SetManualMode();
    fEvtCuts->SetupRun2pp();
    std::cout << "Setting up Track Cuts correspondingly for pp trigger: "
              << trigger << std::endl;
    fEvtCuts->fTriggerMask = trigger;
  }
  if (EvtCutQA) {
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
  fNSPDCluster = obj.fNSPDCluster;
  fRefMult08 = obj.fRefMult08;
  fV0AMult = obj.fV0AMult;
  fV0CMult = obj.fV0CMult;
  fV0MCentrality = obj.fV0MCentrality;
  fnContrib = obj.fnContrib;
  fPassAliEvtSelection = obj.fPassAliEvtSelection;
  fisPileUp = obj.fisPileUp;
  fHasVertex = obj.fHasVertex;
  fHasMagField = obj.fHasMagField;
  fisSelected = obj.fisSelected;
  fEstimator = obj.fEstimator;
  fspher = obj.fspher;
  return (*this);
}

void AliFemtoDreamEvent::SetEvent(AliAODEvent *evt) {
  AliAODVertex *vtx = evt->GetPrimaryVertex();
  AliAODVZERO *vZERO = evt->GetVZEROData();
  AliAODHeader *header = dynamic_cast<AliAODHeader*>(evt->GetHeader());
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
  this->fSPDMult = CalculateITSMultiplicity(evt);
  this->fNSPDCluster = evt->GetMultiplicity()->GetNumberOfSPDClusters();
  this->fV0AMult = vZERO->GetMTotV0A();
  this->fV0CMult = vZERO->GetMTotV0C();
  this->fRefMult08 = header->GetRefMultiplicityComb08();
  this->fspher = CalculateSphericityEvent(evt);

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
    this->fNSPDCluster = evt->GetMultiplicity()->GetNumberOfITSClusters(0, 1);
  } else {
    this->fSPDMult = evt->GetNumberOfITSClusters(1);
    this->fNSPDCluster = 0;
  }
  this->fRefMult08 = AliESDtrackCuts::GetReferenceMultiplicity(
      evt, AliESDtrackCuts::kTrackletsITSTPC, 0.8, 0);
  this->fV0AMult = vZERO->GetMTotV0A();
  this->fV0CMult = vZERO->GetMTotV0C();
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

double AliFemtoDreamEvent::CalculateSphericityEvent(AliAODEvent *evt) {
//Initializing
  double ptTot = 0.;
  double s00 = 0.; //elements of the sphericity matrix taken form EPJC72:2124
  double s01 = 0.;
  double s10 = 0.;
  double s11 = 0.;

  int  numOfTracks = evt->GetNumberOfTracks();
  if(numOfTracks<3) return -9999.;

  int nTracks = 0;
  for (int iTrack = 0; iTrack < numOfTracks; iTrack++) {

    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(evt->GetTrack(iTrack));

    double pt = aodtrack->Pt();
    double px = aodtrack->Px();
    double py = aodtrack->Py();
    double pz = aodtrack->Pz();
    double eta = aodtrack->Eta();

    ptTot +=pt;

    s00 +=px*px/pt;
    s01 +=px*py/pt;
    s10  =s01;
    s11 +=py*py/pt;
    nTracks++;
  }

  //normalize to total Pt to obtain a linear form:
  if(ptTot == 0.) return -9999.;
  s00 /= ptTot;
  s11 /= ptTot;
  s10 /= ptTot;

  //Calculate the trace of the sphericity matrix:
  double T = s00+s11;
  //Calculate the determinant of the sphericity matrix:
  double D = s00*s11 - s10*s10;//S10 = S01

  //Calculate the eigenvalues of the sphericity matrix:
  double lambda1 = 0.5*(T + std::sqrt(T*T - 4.*D));
  double lambda2 = 0.5*(T - std::sqrt(T*T - 4.*D));

  if((lambda1 + lambda2) == 0.) return -9999.;

  double spt = -1.;

  if(lambda2>lambda1)
    {
      spt = 2.*lambda1/(lambda1+lambda2);
    }
  else
    {
      spt = 2.*lambda2/(lambda1+lambda2);
    }

  return spt;
}
