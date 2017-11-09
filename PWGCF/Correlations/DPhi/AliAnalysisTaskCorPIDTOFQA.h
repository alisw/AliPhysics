/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFQA_H
#define AliAnalysisTaskCorPIDTOFQA_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"


//#define AliAnalysisTaskCorPIDTOFQA AliAnalysisTaskCorPIDTOFQA2

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


//	Double_t deut_curves[2][2][3];  /* [charge][mean,sigma][par]  */
//	TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x)",  1.0, 4.4);

//	Double_t cut_width  = 3.0;
	
    private:

        AliAODEvent*          fAOD;               //! input event
        TList*                fOutputList;        //! output list
	AliPIDResponse*       fPIDResponse;
	AliAnalysisUtils*     fAnalysisUtils;
	
	TH1F*                 primary_vertex_z;            //  1
	TH1F*                 primary_vertex_z_cut;        //  2
	TH1I*                 deut_per_event;              //  3

	TH1F*                 fHistPt;                     //  4

	TH2F*                 m2_pt_pos;                   //  5
	TH2F*                 m2_pt_neg;                   //  6

	TH2F*                 m2_pt_pos_cut;               //  7
	TH2F*                 m2_pt_neg_cut;               //  8

	



	
	
        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

//}  //// namespace
#endif
