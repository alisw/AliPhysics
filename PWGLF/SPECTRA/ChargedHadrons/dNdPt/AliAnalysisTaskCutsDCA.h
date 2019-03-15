#ifndef AliAnalysisTaskCutsDCA_H
#define AliAnalysisTaskCutsDCA_H

class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskCutsDCA : public AliAnalysisTaskSE, AlidNdPtTools
{
    public:
                                AliAnalysisTaskCutsDCA();
                                AliAnalysisTaskCutsDCA(const char *name);
        virtual                 ~AliAnalysisTaskCutsDCA();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        virtual void            FinishTaskOutput();
        
        virtual void            SetESDtrackCuts(AliESDtrackCuts* cut) { fESDtrackCuts0 = cut; }
        virtual void            SetESDtrackCutsTPC(AliESDtrackCuts* cut) { fESDtrackCutsTPC = cut; }        
    
    private:    
        AliVEvent*              fEvent;     //! input event       
        AliESDEvent*            fESD;       //! input ESD event
        AliAODEvent*            fAOD;       //! input AOD event
        AliMCEvent*             fMCEvent;   //! MC event
        AliStack*               fMCStack;   //! MC stack
        AliHeader*              fMCHeader;      //! MC header
        AliGenEventHeader*      fMCGenHeader;   //! MC gen header
        
        Double_t                fEventMCzv;   //!
        Double_t                fEventCent;   //! 
        Int_t                   fEventNtracks;   //! 
        Int_t                   fEventNtracksTPC;   //! 
        Int_t                   fEventMultMB;   //! 
        Double_t                fEventMultV0M;   //!
        Double_t                fEventMCb;       //!
        Int_t                   fEventMCnPrim;   //!
        Int_t                   fEventMCnPrimV0M;  //!
        Bool_t                  fEventIsTrigger;   //!
        Bool_t                  fEventHasVertex;   //!

        AliESDtrackCuts*        fESDtrackCuts0;        
        AliESDtrackCuts*        fESDtrackCutsTPC;      
                
        TList*                  fOutputList;    //! output list
        TH1D*                   fLogHist;       //! histogram with errors and logs
        TH1D*                   fLogEvent;      //! event count        
        TH1D*                   fRunHist;       //! run numbers
        THnSparseD*             fHistDCAxy;     //! DCAr TPC-ITS
        THnSparseD*             fHistDCAxyTPC;  //! DCAr TPC
        THnSparseD*             fHistDCAxyMC;   //! DCAr MC
        THnSparseD*             fHistDCAxyMC;   //! DCAr MC
        
      
        void            Initialize();        
        void            Log(const char* name) { Log(fLogHist,name); }        
    
        AliAnalysisTaskCutsDCA(const AliAnalysisTaskCutsDCA&); // not implemented
        AliAnalysisTaskCutsDCA& operator=(const AliAnalysisTaskCutsDCA&); // not implemented

        ClassDef(AliAnalysisTaskCutsDCA, 1);
};

#endif
