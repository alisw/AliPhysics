/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFQA_H
#define AliAnalysisTaskCorPIDTOFQA_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"


class AliAODTrack;


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
	TH2F*                   m_squared_pos_cut_T;
	TH2F*                   m_squared_pos_cut_A;
	TH2F*                   m_squared_pos_cut_B;
	TH2F*                   m_squared_pos_cut_V;
	TH2F*                   m_squared_neg;
	TH2F*                   m_squared_neg_cut_T;
	TH2F*                   m_squared_neg_cut_A;
	TH2F*                   m_squared_neg_cut_B;
	TH2F*                   m_squared_neg_cut_V;
	TH2F*                   plength_vs_mom_pos;
	TH2F*                   plength_vs_mom_neg;
	TH2F*                   ttof_vs_mom_pos;
	TH2F*                   ttof_vs_mom_neg;
	TH2F*                   beta_vs_mom_pos;
	TH2F*                   beta_vs_mom_neg;
	TH2F*                   deltat_vs_mom_pos;
	TH2F*                   deltat_vs_mom_neg;
	TH2F*                   deut_dphi_T;
	TH2F*                   deut_dphi_A;
	TH2F*                   deut_dphi_B;
	TH2F*                   deut_dphi_V;
	TH2F*                   deut_dphi_pos_T;
	TH2F*                   deut_dphi_pos_A;
	TH2F*                   deut_dphi_pos_B;
	TH2F*                   deut_dphi_pos_V;
	TH2F*                   deut_dphi_neg_T;
	TH2F*                   deut_dphi_neg_A;
	TH2F*                   deut_dphi_neg_B;
	TH2F*                   deut_dphi_neg_V;
	TH1I*                   deuterons_per_event;
	TH1I*                   deuterons_per_event_pos;
	TH1I*                   deuterons_per_event_neg;
	TH1F*                   track_phi;
	TH1F*                   track_phi_hybrid;
	TH1F*                   track_eta;
	TH1F*                   track_eta_hybrid;	
	TH2F*                   deut_dphi_deta_p0510;
	TH2F*                   deut_dphi_deta_p1020;
	TH2F*                   deut_dphi_deta_p2030;
	TH2F*                   deut_dphi_deta_p3040;
	TH2F*                   deut_dphi_deta_p4050;
	
        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented
        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

#endif
