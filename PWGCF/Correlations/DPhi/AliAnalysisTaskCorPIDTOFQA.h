/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFQA_H
#define AliAnalysisTaskCorPIDTOFQA_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"


class AliAODTrack;
//class AliEmcalTrackSelection;

//namespace BSchaefer_devel{
    
class AliAnalysisTaskCorPIDTOFQA : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskCorPIDTOFQA();
                                AliAnalysisTaskCorPIDTOFQA(const char *name);
        virtual                 ~AliAnalysisTaskCorPIDTOFQA();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
	virtual Double_t        Beta(AliAODTrack *track);
	virtual Double_t        tof_minus_tpion(AliAODTrack *track);
	virtual Double_t        get_mass_squared(AliAODTrack *track);


	Double_t deut_curves[2][2][3];  /* [charge][mean,sigma][par]  */
	TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x)",  1.1, 4.4);

    private:

        AliAODEvent*            fAOD;               //! input event
        TList*                  fOutputList;        //! output list
	AliPIDResponse*         fPIDResponse;
	
	TH1F*                   fHistPt; 
	TH2F*                   cent_ntracks;

	TH2F*                   m_squared_pos;
	TH2F*                   m_squared_neg;
	TH2F*                   beta_vs_mom_pos;
	TH2F*                   beta_vs_mom_neg;
	TH2F*	                deltat_vs_mom_pos;
	TH2F*	                deltat_vs_mom_neg;
		
	TH2F*                   dedx_vs_mom_pos;
	TH2F*                   dedx_vs_mom_neg;
	TH2F*                   dedx_vs_deltat_pos;
	TH2F*                   dedx_vs_deltat_neg;
	TH3F*                   dedx_mom_deltat_pos;
	TH3F*                   dedx_mom_deltat_neg;
	
	TH2F*                   m_squared_pos_cut;
	TH2F*                   m_squared_neg_cut;	
	TH2F*                   beta_vs_mom_pos_cut;
	TH2F*                   beta_vs_mom_neg_cut;
	TH2F*	                deltat_vs_mom_pos_cut;
	TH2F*	                deltat_vs_mom_neg_cut;
	
	TH2F*                   dedx_vs_mom_pos_cut;
	TH2F*                   dedx_vs_mom_neg_cut;
	TH2F*                   dedx_vs_deltat_pos_cut;
	TH2F*                   dedx_vs_deltat_neg_cut;
	TH3F*                   dedx_mom_deltat_pos_cut;
	TH3F*                   dedx_mom_deltat_neg_cut;

	TH2F*                   dphi_ket_deut_T;
	TH2F*                   dphi_ket_deut_A;
	TH2F*                   dphi_ket_deut_B;

	
	TH2F*                   m_squared_pos_cut_T;
	TH2F*                   m_squared_pos_cut_A;
	TH2F*                   m_squared_pos_cut_B;
	TH2F*                   m_squared_neg_cut_T;
	TH2F*                   m_squared_neg_cut_A;
	TH2F*                   m_squared_neg_cut_B;
	TH2F*                   deut_dphi_T;
	TH2F*                   deut_dphi_A;
	TH2F*                   deut_dphi_B;
	TH2F*                   deut_dphi_pos_T;
	TH2F*                   deut_dphi_pos_A;
	TH2F*                   deut_dphi_pos_B;
	TH2F*                   deut_dphi_neg_T;
	TH2F*                   deut_dphi_neg_A;
	TH2F*                   deut_dphi_neg_B;
	
	TH1I*                   deuterons_per_event;
	TH1I*                   deuterons_per_event_pos;
	TH1I*                   deuterons_per_event_neg;
	
//	TH2F*                   mom_compare_pos;
//	TH2F*                   mom_compare_neg;
//	TH1F*                   path_length;
//	TH1F*                   ttof;
//	TH1F*                   alpha;
//	TH1F*                   theta;
//	TH2F*                   dca;
//	TH1F*                   generic;
	
        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

//}  //// namespace
#endif
