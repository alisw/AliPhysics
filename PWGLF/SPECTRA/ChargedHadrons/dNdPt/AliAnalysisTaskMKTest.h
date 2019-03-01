#ifndef AliAnalysisTaskMKTest_H
#define AliAnalysisTaskMKTest_H

class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"
#include "AlidNdPtTools.h"


class AliAnalysisTaskMKTest : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskMKTest();
                                AliAnalysisTaskMKTest(const char *name);
        virtual                 ~AliAnalysisTaskMKTest();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        virtual void            FinishTaskOutput();
        
        virtual void            SetESDtrackCuts(AliESDtrackCuts* cut) { fESDtrackCuts0 = cut; }
    
    private:    
        AliVEvent*              fEvent;     //! input event           
        AliESDEvent*            fESD;       //! input ESD event
        AliMCEvent*             fMCEvent;   //! MC event
        AliStack*               fMCStack;   //! MC stack
        AliHeader*              fMCHeader;      //! MC header
        AliGenEventHeader*      fMCGenHeader;   //! MC gen header
        AliESDtrackCuts*        fESDtrackCuts0; //track cuts 0
        TList*                  fOutputList;    //! output list
        TH1D*                   fLogHist;       //! histogram with errors and logs
        THnSparseD*             fHistDCA;     //! DCA

        void                    Log(const char* name);
    
        AliAnalysisTaskMKTest(const AliAnalysisTaskMKTest&); // not implemented
        AliAnalysisTaskMKTest& operator=(const AliAnalysisTaskMKTest&); // not implemented

        ClassDef(AliAnalysisTaskMKTest, 1);
};

#endif
