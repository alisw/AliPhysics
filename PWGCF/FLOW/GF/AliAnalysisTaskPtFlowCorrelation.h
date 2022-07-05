/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskPtFlowCorrelation_H
#define AliAnalysisTaskPtFlowCorrelation_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPtFlowCorrelation : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskPtFlowCorrelation();
                                AliAnalysisTaskPtFlowCorrelation(const char *name);
        virtual                 ~AliAnalysisTaskPtFlowCorrelation();

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

        TFile*                  inputFile;          //!
        TProfile*               inputProfileMeanPt; //!
        TAxis*                  inputAxisMeanPt;    //!

        TH1D*                   inputHistPhi;       //!
        TAxis*                  inputAxisPhi;       //!
        Int_t                   N_max;              //


        Double_t GetMeanPt(Double_t centrality);
        Double_t GetPhiWeight(Double_t phi);

        AliAnalysisTaskPtFlowCorrelation(const AliAnalysisTaskPtFlowCorrelation&); // not implemented
        AliAnalysisTaskPtFlowCorrelation& operator=(const AliAnalysisTaskPtFlowCorrelation&); // not implemented

        ClassDef(AliAnalysisTaskPtFlowCorrelation, 1);
};

#endif
