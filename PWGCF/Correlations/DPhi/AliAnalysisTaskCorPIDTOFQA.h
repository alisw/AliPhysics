/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFQA_H
#define AliAnalysisTaskCorPIDTOFQA_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"


class AliAODTrack;
//class AliEmcalTrackSelection;

////namespace BSchaefer_devel{
    
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
	TH2F*                   m_squared_pos_cut;
	TH2F*                   m_squared_neg;
	TH2F*                   m_squared_neg_cut;
	
	TH2F*                   plength_vs_mom_pos;
	TH2F*                   plength_vs_mom_neg;
	TH2F*                   ttof_vs_mom_pos;
	TH2F*                   ttof_vs_mom_neg;
	TH2F*                   beta_vs_mom_pos;
	TH2F*                   beta_vs_mom_neg;
	TH2F*	                deltat_vs_mom_pos;
	TH2F*	                deltat_vs_mom_neg;
	
	TH2F*                   plength_vs_mom_pos_cut;
	TH2F*                   plength_vs_mom_neg_cut;
	TH2F*                   ttof_vs_mom_pos_cut;
	TH2F*                   ttof_vs_mom_neg_cut;
	TH2F*                   beta_vs_mom_pos_cut;
	TH2F*                   beta_vs_mom_neg_cut;
	TH2F*	                deltat_vs_mom_pos_cut;
	TH2F*	                deltat_vs_mom_neg_cut;

	TH1F*                   h_dedx_exp_deut;
	TH1F*                   h_dedx_exp_prot;
	TH1F*                   h_dedx;
	TH1F*                   h_TPC_resolution_deut;
	TH1F*                   h_TPC_resolution_prot;
	
        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

////}  //// namespace
#endif
