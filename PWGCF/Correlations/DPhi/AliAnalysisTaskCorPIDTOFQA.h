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

	Double_t cut_width  = 2.0;
	
    private:

        AliAODEvent*          fAOD;               //! input event
        TList*                fOutputList;        //! output list
	AliPIDResponse*       fPIDResponse;
	
	TH1F*                 fHistPt;                     //  1
	TH2F*                 cent_ntracks;                //  2
	TH2F*                 m2_pt_pos;                   //  3
	TH2F*                 m2_pt_neg;                   //  4
	TH2F*                 beta_p_pos;                  //  5
	TH2F*                 beta_p_neg;                  //  6
	TH2F*	              deltat_pt_pos;               //  7
	TH2F*	              deltat_pt_neg;               //  8

	TH2F*                 m2_pt_pos_cut;               //  9
	TH2F*                 m2_pt_neg_cut;	           // 10
	TH2F*                 beta_p_pos_cut;              // 11
	TH2F*                 beta_p_neg_cut;              // 12
	TH2F*	              deltat_pt_pos_cut;           // 13
	TH2F*	              deltat_pt_neg_cut;           // 14

	TH1I*                 deut_per_event;              // 15
//	TH1I*                 deut_per_event_pos;          // 16
//	TH1I*                 deut_per_event_neg;  	   // 17
	TH2F*                 m2_pt_pos_cut_T;             // 18
	TH2F*                 m2_pt_neg_cut_T;             // 19

	TH2F*                 deut_phi_pt_pos;             // 20
	TH2F*                 deut_phi_pt_neg;             // 21
	
//	TH2F*                 deut_q2_pt_pos_pos;          // 22
//	TH2F*                 deut_q2_pt_pos_neg;          // 23
//	TH2F*                 deut_q2_pt_neg_neg;          // 24

	TH2F*                 deut_phi_pt;                 // 25
	TH2F*                 trig_05_phi_pt;              // 26a
	TH2F*                 trig_05_phi_pt_pos;          // 26b
	TH2F*                 trig_05_phi_pt_neg;          // 26c
	
	TH2F*                 trig_08_phi_pt;              // 27a
	TH2F*                 trig_08_phi_pt_pos;          // 27b
	TH2F*                 trig_08_phi_pt_neg;          // 27c
	
	TH2F*                 deut_dphi_pt_pos_pos_05;     // 28
	TH2F*                 deut_dphi_pt_pos_neg_05;     // 29
	TH2F*                 deut_dphi_pt_neg_neg_05;     // 30

	
	TH2F*                 deut_dphi_pt_pos_pos_08;     // 31
	TH2F*                 deut_dphi_pt_pos_neg_08;     // 32
	TH2F*                 deut_dphi_pt_neg_neg_08;     // 33

//	TH2F*                 track_cor_radius_pt;         // 34
//	TH2F*                 track_cor_radius_pt_cut;     // 35
	
        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

//}  //// namespace
#endif
