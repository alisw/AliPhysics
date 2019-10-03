/*
 * AliLightNEvent.cxx
 *
 *  Created on: 22 Nov 2017
 *      Author: bernhardhohlweger
 */

#include "AliLightNEvent.h"
ClassImp(AliLightNEvent)
AliLightNEvent::AliLightNEvent()
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
,fnContrib(0)
,fPassAliEvtSelection(false)
,fisPileUp(false)
,fHasVertex(false)
,fHasMagField(false)
,fisSelected(false)
{
    
}

AliLightNEvent::AliLightNEvent(bool mvPileUp,bool EvtCutQA)
:fUtils(new AliAnalysisUtils())
,fEvtCuts(new AliEventCuts())
,fxVtx(0)
,fyVtx(0)
,fzVtx(0)
,fSPDMult(0)
,fRefMult08(0)
,fV0AMult(0)
,fV0CMult(0)
,fnContrib(0)
,fPassAliEvtSelection(false)
,fisPileUp(false)
,fHasVertex(false)
,fHasMagField(false)
,fisSelected(false)
{
    //  if (mvPileUp) {
    //    //For pPb this is necessary according to DPG Processing status news
    //    //(week 29 April - 5 May 2017)
    //    fUtils->SetUseMVPlpSelection(true);
    //  } else {
    //    //Following the analysis in pp Run1 of O.Arnold
    //    fUtils->SetMinPlpContribSPD(3);
    //  }
    if (EvtCutQA) {
        fEvtCutList=new TList();
        fEvtCutList->SetName("AliEventCuts");
        fEvtCutList->SetOwner(true);
        fEvtCuts->AddQAplotsToList(fEvtCutList);
    } else {
        fEvtCutList=nullptr;
    }
}

AliLightNEvent::~AliLightNEvent() {
    if (fEvtCutList) {
        delete fEvtCutList;
    }
    if (fEvtCutList) {
        delete fEvtCutList;
    }
    if (fEvtCuts) {
        delete fEvtCuts;
    }
    if (fUtils) {
        delete fUtils;
    }
}

void AliLightNEvent::SetEvent(AliAODEvent *evt) {
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
    }else{
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
    return;
}

int AliLightNEvent::CalculateITSMultiplicity(AliAODEvent *evt) {
    AliAODTracklets* tracklets=evt->GetTracklets();
    int nTr=tracklets->GetNumberOfTracklets();
    int count=0;
    for(int iTr=0; iTr<nTr; iTr++){
        double theta=tracklets->GetTheta(iTr);
        double eta=-TMath::Log(TMath::Tan(theta/2.));
        if(TMath::Abs(eta) < 0.8){
            count++;
        }
    }
    return count;
}
