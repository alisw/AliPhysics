#ifndef AliAnalysisTaskFlowExample_H
#define AliAnalysisTaskFlowExample_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFlowExample : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskFlowExample();
                                AliAnalysisTaskFlowExample(const char *name);
        virtual                 ~AliAnalysisTaskFlowExample();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TH2F*                   fHistPhiEta;    //!

        AliAnalysisTaskFlowExample(const AliAnalysisTaskFlowExample&); // not implemented
        AliAnalysisTaskFlowExample& operator=(const AliAnalysisTaskFlowExample&); // not implemented

        ClassDef(AliAnalysisTaskFlowExample, 1);
};

#endif
