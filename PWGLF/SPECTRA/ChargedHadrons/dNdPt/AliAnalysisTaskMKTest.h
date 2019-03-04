#ifndef AliAnalysisTaskMKTest_H
#define AliAnalysisTaskMKTest_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMKTest : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskMKTest();
                                AliAnalysisTaskMKTest(const char *name);
        virtual                 ~AliAnalysisTaskMKTest();
    
    private:    
    
//         AliAnalysisTaskMKTest(const AliAnalysisTaskMKTest&); // not implemented
//         AliAnalysisTaskMKTest& operator=(const AliAnalysisTaskMKTest&); // not implemented

        ClassDef(AliAnalysisTaskMKTest, 1);
};

#endif
