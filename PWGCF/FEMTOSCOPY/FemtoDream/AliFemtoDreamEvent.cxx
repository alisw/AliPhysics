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

ClassImp(AliFemtoDreamEvent)
AliFemtoDreamEvent::AliFemtoDreamEvent()
:fUtils(new AliAnalysisUtils())
,fEvtCuts(new AliEventCuts())
,fEvtCutList()
,fxVtx(0)
,fyVtx(0)
,fzVtx(0)
,fSPDMult(0)
,fRefMult08(0)
,fV0AMult(0)
,fV0CMult(0)
,fV0MCentrality(0)
,fnContrib(0)
,fPassAliEvtSelection(false)
,fisPileUp(false)
,fHasVertex(false)
,fHasMagField(false)
,fisSelected(false)
,fEstimator(AliFemtoDreamEvent::kRef08)
{

}

AliFemtoDreamEvent::AliFemtoDreamEvent(
    bool mvPileUp,bool EvtCutQA, UInt_t trigger)
:fUtils(new AliAnalysisUtils())
,fEvtCuts(new AliEventCuts())
,fxVtx(0)
,fyVtx(0)
,fzVtx(0)
,fSPDMult(0)
,fRefMult08(0)
,fV0AMult(0)
,fV0CMult(0)
,fV0MCentrality(0)
,fnContrib(0)
,fPassAliEvtSelection(false)
,fisPileUp(false)
,fHasVertex(false)
,fHasMagField(false)
,fisSelected(false)
,fEstimator(kRef08)
{
  if (mvPileUp) {
    //For pPb this is necessary according to DPG Processing status news
    //(week 29 April - 5 May 2017)
    fUtils->SetUseMVPlpSelection(true);
  } else {
    //Following the analysis in pp Run1 of O.Arnold
    fUtils->SetMinPlpContribSPD(3);
  }

  if(trigger != AliVEvent::kINT7) {
    fEvtCuts->SetManualMode();
    fEvtCuts->SetupRun2pp();
    std::cout << "Setting up Track Cuts correspondingly for pp trigger: " <<
        trigger << std::endl;
    fEvtCuts->fTriggerMask = trigger;
  }
  if (EvtCutQA) {
    fEvtCutList=new TList();
    fEvtCutList->SetName("AliEventCuts");
    fEvtCutList->SetOwner(true);
    fEvtCuts->AddQAplotsToList(fEvtCutList);
  } else {
    fEvtCutList=nullptr;
  }
}

AliFemtoDreamEvent::~AliFemtoDreamEvent() {
  if (fEvtCutList) {
    delete fEvtCutList;
  }
  if (fEvtCutList) {
    delete fEvtCutList;
  }
  if (fEvtCuts) {
    delete fEvtCuts;
  }
//  if (fUtils) {
//    delete fUtils;
//  }
}

void AliFemtoDreamEvent::SetEvent(AliAODEvent *evt) {
  AliAODVertex *vtx=evt->GetPrimaryVertex();
  AliAODVZERO *vZERO = evt->GetVZEROData();
  AliAODHeader *header = dynamic_cast<AliAODHeader*>(evt->GetHeader());
  if (!vtx) {
    this->fHasVertex=false;
  } else {
    this->fHasVertex=true;
  }
  if (TMath::Abs(evt->GetMagneticField())< 0.001) {
    this->fHasMagField=false;
  } else {
    this->fHasMagField=true;
  }
  this->fnContrib=vtx->GetNContributors();
  this->fxVtx=vtx->GetX();
  this->fyVtx=vtx->GetY();
  this->fzVtx=vtx->GetZ();
  if (fUtils->IsPileUpEvent(evt)) {
    this->fisPileUp=true;
  } else {
    this->fisPileUp=false;
  }
  if (fEvtCuts->AcceptEvent(evt)) {
    this->fPassAliEvtSelection=true;
  } else {
    this->fPassAliEvtSelection=false;
  }
  this->fSPDMult=CalculateITSMultiplicity(evt);
  this->fV0AMult=vZERO->GetMTotV0A();
  this->fV0CMult=vZERO->GetMTotV0C();
  this->fRefMult08=header->GetRefMultiplicityComb08();
  float lPercentile = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection * ) evt->FindListObject("MultSelection");
  if( !MultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  }else{
    fV0MCentrality= MultSelection->GetMultiplicityPercentile("V0M");
  }
  return;
}

int AliFemtoDreamEvent::CalculateITSMultiplicity(AliAODEvent *evt) {
  AliAODTracklets* tracklets=evt->GetTracklets();
  int nTr=tracklets->GetNumberOfTracklets();
  int count=0;
  for(int iTr=0; iTr<nTr; iTr++){
    float theta=tracklets->GetTheta(iTr);
    float eta=-TMath::Log(TMath::Tan(theta/2.));
    if(TMath::Abs(eta) < 0.8){
      count++;
    }
  }
  return count;
}

int AliFemtoDreamEvent::GetMultiplicity() {
  int mult=0;
  switch(fEstimator) {
    case AliFemtoDreamEvent::kRef08:
      mult=GetRefMult08();
      break;
    case AliFemtoDreamEvent::kSPD:
      mult=GetSPDMult();
      break;
    case AliFemtoDreamEvent::kV0M:
      mult=GetV0MMult();
      break;
    case AliFemtoDreamEvent::kV0A:
      mult=GetV0AMult();
      break;
    case AliFemtoDreamEvent::kV0C:
      mult=GetV0CMult();
      break;
    default:
      AliFatal("Type Not implemented");
      break;
  }

  return mult;
}
