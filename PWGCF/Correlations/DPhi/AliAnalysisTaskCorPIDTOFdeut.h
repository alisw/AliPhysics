/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFdeut_H
#define AliAnalysisTaskCorPIDTOFdeut_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"


class AliAODTrack;
//class AliEmcalTrackSelection;

//namespace BSchaefer_devel{
    
class AliAnalysisTaskCorPIDTOFdeut : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskCorPIDTOFdeut();
                                AliAnalysisTaskCorPIDTOFdeut(const char *name);
        virtual                 ~AliAnalysisTaskCorPIDTOFdeut();

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
	
	TH1F*                   fHistPt;                //  1
	TH2F*                   cent_ntracks;           //  2
	
	TH2F*                   m2_pos;                 //  3
	TH2F*                   m2_neg;                 //  4
	TH2F*                   beta_p_pos;             //  5
	TH2F*                   beta_p_neg;             //  6
	TH2F*	                deltat_p_pos;           //  7
	TH2F*	                deltat_p_neg;           //  8
	TH2F*                   m2_pos_cut;             //  9
	TH2F*                   m2_neg_cut;	        // 10
	TH2F*                   beta_p_pos_cut;         // 11
	TH2F*                   beta_p_neg_cut;         // 12
	TH2F*	                deltat_p_pos_cut;       // 13
	TH2F*	                deltat_p_neg_cut;       // 14
	TH2F*                   m2_pos_cut_T;           // 15
	TH2F*                   m2_neg_cut_T;           // 16
	TH1I*                   deut_per_event;         // 17	
	
	TH2F*                   dphi_pt_deut_C;         // 18
	TH2F*                   dphi_kt_deut_C;         // 19
	TH2F*                   dphi_et_deut_C;         // 20
	TH2F*                   dphi_en_deut_C;         // 21

	TH1F*                   deut_pt_count_C;        // 22
	TH1F*                   deut_kt_count_C;        // 23
	TH1F*                   deut_et_count_C;        // 24
	TH1F*                   deut_en_count_C;        // 25

	TH2F*                   dphi_pt_deut_P;         // 26
	TH2F*                   dphi_kt_deut_P;         // 27
	TH2F*                   dphi_et_deut_P;         // 28
	TH2F*                   dphi_en_deut_P;         // 29

	TH1F*                   deut_pt_count_P;        // 30
	TH1F*                   deut_kt_count_P;        // 31
	TH1F*                   deut_et_count_P;        // 32
	TH1F*                   deut_en_count_P;        // 33
	
	TH1F*                   centrality_dist;        // 34
	TH1F*                   multiplicity_dist;      // 35
	TH1F*                   centrality_raw;         // 36
	TH1F*                   multiplicity_raw;       // 37
	
	TH2F*                   deut_dphi_deta_p1_0510; // 38
	TH2F*                   deut_dphi_deta_p1_1020; // 39
	TH2F*                   deut_dphi_deta_p1_2030; // 40
	TH2F*                   deut_dphi_deta_p1_3040; // 41
	TH2F*                   deut_dphi_deta_p1_4050; // 42

	TH2F*                   deut_dphi_deta_p2_0510; // 43
	TH2F*                   deut_dphi_deta_p2_1020; // 44
	TH2F*                   deut_dphi_deta_p2_2030; // 45
	TH2F*                   deut_dphi_deta_p2_3040; // 46
	TH2F*                   deut_dphi_deta_p2_4050; // 47
	
	TH2F*                   deut_dphi_deta_p3_0510; // 48
	TH2F*                   deut_dphi_deta_p3_1020; // 49
	TH2F*                   deut_dphi_deta_p3_2030; // 50
	TH2F*                   deut_dphi_deta_p3_3040; // 51
	TH2F*                   deut_dphi_deta_p3_4050; // 52
	
	TH2F*                   deut_cor_radius;        // 53

	
        AliAnalysisTaskCorPIDTOFdeut(const AliAnalysisTaskCorPIDTOFdeut&);                        // not implemented
        AliAnalysisTaskCorPIDTOFdeut& operator=(const AliAnalysisTaskCorPIDTOFdeut&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFdeut, 1);
};

//}  //// namespace
#endif
