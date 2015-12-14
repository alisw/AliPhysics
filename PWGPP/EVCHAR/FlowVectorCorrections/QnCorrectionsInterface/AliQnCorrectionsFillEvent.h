#ifndef ALIQNCORRECTIONS_FILLEVENT_H
#define ALIQNCORRECTIONS_FILLEVENT_H

#include <AliQnCorrectionsManager.h>
#include "AliQnCorrectionsHistos.h"
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskSE.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>

class AliQnCorrectionsFillEvent : public TNamed{
  public:


    AliQnCorrectionsFillEvent();
    ~AliQnCorrectionsFillEvent();


    void Process(AliAnalysisTaskSE* task, AliVEvent* event, Float_t* values);

    void FillDetectors(AliAnalysisTaskSE* task, Float_t* values);
    void FillTPC(Float_t* values);
    void FillEsdTPC(Float_t* values);
    void FillAodTPC(Float_t* values);
    void FillVZERO();
    void FillTZERO();
    void FillZDC();
    void FillFMD(AliAnalysisTaskSE* task);
    void FillRawFMD(Float_t* values);
    void FillSPDTracklets(Float_t* values);

    void FillEventInfo(Float_t* values); 
    void FillTrackInfo(AliESDtrack* p, Float_t* values);
    void FillTrackInfo(AliVParticle* p, Float_t* values);

    void SetUseTPCStandaloneTracks(Bool_t b=kTRUE) {fUseTPCStandaloneTracks=b;}

    void SetEvent(AliVEvent* ev) {fEvent=ev;}
    void SetEventPlaneManager(AliQnCorrectionsManager* EPmanager) {fEventPlaneManager=EPmanager; SetDetectors();}
    void SetEventPlaneHistos(AliQnCorrectionsHistos* EPmanager) {fQAhistos=EPmanager;}
    void SetDetectors();


  private:

    AliQnCorrectionsFillEvent(const AliQnCorrectionsFillEvent &c);
    AliQnCorrectionsFillEvent& operator= (const AliQnCorrectionsFillEvent &c);

    AliVEvent* fEvent;
    AliQnCorrectionsManager* fEventPlaneManager;
    AliQnCorrectionsHistos* fQAhistos;

    Bool_t fUseTPCStandaloneTracks;

    Bool_t fFillVZERO;
    Bool_t fFillTPC;
    Bool_t fFillZDC;
    Bool_t fFillTZERO;
    Bool_t fFillFMD;
    Bool_t fFillRawFMD;
    Bool_t fFillSPD;

    Bool_t fIsAOD;
    Bool_t fIsESD;

    ClassDef(AliQnCorrectionsFillEvent, 1);

};

#endif
