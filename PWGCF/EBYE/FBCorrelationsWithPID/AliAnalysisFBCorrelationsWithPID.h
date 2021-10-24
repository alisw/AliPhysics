/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisFBCorrelationsWithPID_H
#define AliAnalysisFBCorrelationsWithPID_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisFBCorrelationsWithPID : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisFBCorrelationsWithPID();
                                AliAnalysisFBCorrelationsWithPID(const char *name);
        virtual                 ~AliAnalysisFBCorrelationsWithPID();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TH1F*                   fHistPt;        //! dummy histogram
        TH2F*                   fHistTPCSignPLog;  //!
        
        AliAnalysisFBCorrelationsWithPID(const AliAnalysisFBCorrelationsWithPID&); // not implemented
        AliAnalysisFBCorrelationsWithPID& operator=(const AliAnalysisFBCorrelationsWithPID&); // not implemented

        ClassDef(AliAnalysisFBCorrelationsWithPID, 1);
};

#endif
