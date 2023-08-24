#ifndef AliGFWFilterTask_H
#define AliGFWFilterTask_H

#include <vector>
#include <complex>

#include "TChain.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TFile.h"
#include "TProfile.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "AliGFWFilter.h"

class AliGFWFilterTask : public AliAnalysisTaskSE
{
    public:
        AliGFWFilterTask();
        AliGFWFilterTask(const char *name);
        virtual ~AliGFWFilterTask();
        virtual void UserCreateOutputObjects();
        virtual void UserExec(Option_t* option);
        virtual void Terminate(Option_t* option);
        virtual void NotifyRun();
        void EmbedEventSelection(Bool_t newVal, Bool_t addQA, UInt_t lTrigger=0) {fEmbedES=newVal; fAddQA=addQA; fTrigger=lTrigger; };
        void SetExtentV0MAcceptance(Bool_t newVal) { fExtendV0MAcceptance=newVal; };
        void AddCustomCuts(UInt_t lEvCuts, UInt_t lTrCuts) {fCustomCuts.push_back(make_pair(lEvCuts,lTrCuts)); };
        void SetDisableDefaultCuts(Bool_t newval) {fDisableDefaultCuts=newval;};
        void SetupLWCuts(Bool_t fb96=kTRUE, Bool_t fb768=kTRUE, Bool_t use25Chi=kTRUE); //Run only with FB96 &/ FB768 cuts
        void SetDefaultChi2Cut(GFWFlags::kLocalTrackFlags newval) { fStandardChi2Cut = newval; };
        void SetPt(Double_t lPtMin, Double_t lPtMax) {fPtMin=lPtMin; fPtMax=lPtMax; };
        void SetEta(Double_t lEtaMin, Double_t lEtaMax) {fEtaMin=lEtaMin; fEtaMax=lEtaMax; };
    private:
        AliAODEvent *fAOD;
        TList *fOutList;
        AliEventCuts *fEventCuts; //!
        AliGFWFilter *fFilter; //!
        Bool_t fEmbedES;
        Bool_t fAddQA;
        Bool_t fTrigger;
        Bool_t fExtendV0MAcceptance;
        vector< pair<UInt_t, UInt_t> > fCustomCuts;
        Bool_t fDisableDefaultCuts;
        GFWFlags::kLocalTrackFlags fStandardChi2Cut;
        Double_t fPtMin;
        Double_t fPtMax;
        Double_t fEtaMin;
        Double_t fEtaMax;
        ClassDef(AliGFWFilterTask, 2);
};

#endif
