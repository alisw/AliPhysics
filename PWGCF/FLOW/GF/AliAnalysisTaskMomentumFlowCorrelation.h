/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskMomentumFlowCorrelation_H
#define AliAnalysisTaskMomentumFlowCorrelation_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMomentumFlowCorrelation : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskMomentumFlowCorrelation();
                                AliAnalysisTaskMomentumFlowCorrelation(const char *name);
        virtual                 ~AliAnalysisTaskMomentumFlowCorrelation();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TProfile*               fProfileRho;    //! 
        TProfile*               fProfileVar;    //! 
        TProfile*               fProfileCov;    //! 
        TProfile*               fProfileC2;    //! 

        TProfile*               fProfileMeanPt;    //! 
        TProfile*               fProfileTwoParPtCorr;    //! 

        TProfile*               fProfileOneParPtCorrCent;    //! 
        TProfile*               fProfileTwoParPtCorrCent;    //! 

        TProfile*               fProfileTwoParCorr;         //! 
        TProfile*               fProfileThreeParCorr;       //! 
        TProfile*               fProfileFourParCorr;        //! 


        AliAnalysisTaskMomentumFlowCorrelation(const AliAnalysisTaskMomentumFlowCorrelation&); // not implemented
        AliAnalysisTaskMomentumFlowCorrelation& operator=(const AliAnalysisTaskMomentumFlowCorrelation&); // not implemented

        ClassDef(AliAnalysisTaskMomentumFlowCorrelation, 1);
};

#endif
