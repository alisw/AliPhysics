/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFprot_H
#define AliAnalysisTaskCorPIDTOFprot_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"


class AliAODTrack;
//class AliEmcalTrackSelection;

//namespace BSchaefer_devel{
    
class AliAnalysisTaskCorPIDTOFprot : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskCorPIDTOFprot();
                                AliAnalysisTaskCorPIDTOFprot(const char *name);
        virtual                 ~AliAnalysisTaskCorPIDTOFprot();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
	virtual Double_t        Beta(AliAODTrack *track);
	virtual Double_t        tof_minus_tpion(AliAODTrack *track);
	virtual Double_t        get_mass_squared(AliAODTrack *track);


	Double_t prot_curves[2][2][4];  /* [charge][mean,sigma][par]  */
	TF1 *fit_prot_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) + [3]*x^5",  1.0, 4.0);

	Double_t cut_width      = 1.0;
	
    private:

        AliAODEvent*            fAOD;               //! input event
        TList*                  fOutputList;        //! output list
	AliPIDResponse*         fPIDResponse;
	
	TH1F*                   fHistPt;                    //  1
	TH2F*                   cent_ntracks;               //  2
	TH2F*                   m2_pos;                     //  3
	TH2F*                   m2_neg;                     //  4
	TH2F*                   beta_p_pos;                 //  5
	TH2F*                   beta_p_neg;                 //  6
//	TH2F*	                deltat_p_pos;               //  7
//	TH2F*	                deltat_p_neg;               //  8

	TH2F*                   m2_pos_cut;                 //  9
	TH2F*                   m2_neg_cut;	            // 10
	TH2F*                   beta_p_pos_cut;             // 11
	TH2F*                   beta_p_neg_cut;             // 12
//	TH2F*	                deltat_p_pos_cut;           // 13
//	TH2F*	                deltat_p_neg_cut;           // 14

	TH1I*                   prot_per_event;             // 15
	TH1I*                   prot_per_event_pos;         // 16
	TH1I*                   prot_per_event_neg;	    // 17
	TH2F*                   m2_pos_cut_T;               // 18
	TH2F*                   m2_neg_cut_T;               // 19

	TH1F*                   prot_phi_hist;              // 20

	TH2F*                   pt_mom;                     // 21
	TH2F*                   proton_pt_q_pos_pos;        // 22
	TH2F*                   proton_pt_q_pos_neg;        // 23
	TH2F*                   proton_pt_q_neg_neg;        // 24

	TH1F*                   di_proton_phi;              // 25
	TH1F*                   hi_05_phi;                  // 26
	TH1F*                   hi_08_phi;                  // 27

	TH2F*                   di_prot_q_dphi_pos_pos_05;  // 28
	TH2F*                   di_prot_q_dphi_pos_neg_05;  // 29
	TH2F*                   di_prot_q_dphi_neg_neg_05;  // 30

	
	TH2F*                   di_prot_q_dphi_pos_pos_08;  // 31
	TH2F*                   di_prot_q_dphi_pos_neg_08;  // 32
	TH2F*                   di_prot_q_dphi_neg_neg_08;  // 33

	
        AliAnalysisTaskCorPIDTOFprot(const AliAnalysisTaskCorPIDTOFprot&);                        // not implemented
        AliAnalysisTaskCorPIDTOFprot& operator=(const AliAnalysisTaskCorPIDTOFprot&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFprot, 1);
};

//}  //// namespace
#endif
