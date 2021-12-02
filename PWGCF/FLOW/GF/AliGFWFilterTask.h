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
        ClassDef(AliGFWFilterTask, 1);
};

#endif
