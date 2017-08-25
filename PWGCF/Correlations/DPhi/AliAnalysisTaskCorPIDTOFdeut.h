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
	TH2F*                   dedx_p_pos;             //  9
	TH2F*                   dedx_p_neg;             // 10
	TH2F*                   dedx_deltat_pos;        // 11
	TH2F*                   dedx_deltat_neg;        // 12
	TH3F*                   dedx_p_deltat_pos;      // 13
	TH3F*                   dedx_p_deltat_neg;      // 14
	TH2F*                   m2_pos_cut;             // 15
	TH2F*                   m2_neg_cut;	        // 16
	TH2F*                   beta_p_pos_cut;         // 17
	TH2F*                   beta_p_neg_cut;         // 18
	TH2F*	                deltat_p_pos_cut;       // 19
	TH2F*	                deltat_p_neg_cut;       // 20
	TH2F*                   dedx_p_pos_cut;         // 21
	TH2F*                   dedx_p_neg_cut;         // 22
	TH2F*                   dedx_deltat_pos_cut;    // 23
	TH2F*                   dedx_deltat_neg_cut;    // 24
	TH3F*                   dedx_p_deltat_pos_cut;  // 25
	TH3F*                   dedx_p_deltat_neg_cut;  // 26
	TH2F*                   deut_dphi_T;            // 27
	TH2F*                   deut_dphi_pos_T;        // 28
	TH2F*                   deut_dphi_neg_T;        // 29
//	TH2F*                   deut_dphi_A;            // 30
//	TH2F*                   deut_dphi_pos_A;        // 31
//	TH2F*                   deut_dphi_neg_A;	// 32
//	TH2F*                   deut_dphi_B;            // 33
//	TH2F*                   deut_dphi_pos_B;        // 34
//	TH2F*                   deut_dphi_neg_B;        // 35
	TH1I*                   deut_per_event;         // 36
	TH1I*                   deut_per_event_pos;     // 37
	TH1I*                   deut_per_event_neg;	// 38
	TH2F*                   m2_pos_cut_T;           // 39
//	TH2F*                   m2_pos_cut_A;           // 40
//	TH2F*                   m2_pos_cut_B;           // 41
	TH2F*                   m2_neg_cut_T;           // 42
//	TH2F*                   m2_neg_cut_A;           // 43
//	TH2F*                   m2_neg_cut_B;           // 44
	TH2F*                   dphi_et_deut_T;         // 45
//	TH2F*                   dphi_et_deut_A;         // 46
//	TH2F*                   dphi_et_deut_B;         // 47
	TProfile*               deltat_channel;         // 48

	TH1F*                   deut_et_count;          // 49
	TH1F*                   deut_e_count;           // 50
	TH2F*                   dphi_e_deut_T;          // 51
	TH1F*                   deut_pt_count;          // 52
	TH1F*                   asso_phi_hist;          // 53
	TH1F*                   deut_phi_hist;          // 54

	
        AliAnalysisTaskCorPIDTOFdeut(const AliAnalysisTaskCorPIDTOFdeut&);                        // not implemented
        AliAnalysisTaskCorPIDTOFdeut& operator=(const AliAnalysisTaskCorPIDTOFdeut&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFdeut, 1);
};

//}  //// namespace
#endif
