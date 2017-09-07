/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFQA_H
#define AliAnalysisTaskCorPIDTOFQA_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"


class AliAODTrack;

//namespace BSchaefer_devel{
    
class AliAnalysisTaskCorPIDTOFQA : public AliAnalysisTaskSE  
{
    public:
                              AliAnalysisTaskCorPIDTOFQA();
                              AliAnalysisTaskCorPIDTOFQA(const char *name);
        virtual              ~AliAnalysisTaskCorPIDTOFQA();

        virtual void          UserCreateOutputObjects();
        virtual void          UserExec(Option_t* option);
        virtual void          Terminate(Option_t* option);
	virtual Double_t      Beta(AliAODTrack *track);
	virtual Double_t      tof_minus_tpion(AliAODTrack *track);
	virtual Double_t      get_mass_squared(AliAODTrack *track);


	Double_t deut_curves[2][2][3];  /* [charge][mean,sigma][par]  */
	TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x)",  1.0, 4.4);

	Double_t cut_width  = 4.0;
	
    private:

        AliAODEvent*          fAOD;               //! input event
        TList*                fOutputList;        //! output list
	AliPIDResponse*       fPIDResponse;


	TTree*                htree;
	AliAODEvent*          helperAOD;                   //! input event
	
	TH1F*                 fHistPt;                     //  1

	TH1I*                 deut_per_event;              // 15
	
        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

//}  //// namespace
#endif
