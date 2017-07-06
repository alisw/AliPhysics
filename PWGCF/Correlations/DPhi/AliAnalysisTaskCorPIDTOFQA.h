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
        virtual                ~AliAnalysisTaskCorPIDTOFQA();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
	virtual Double_t        Beta(AliAODTrack *track);
	virtual Double_t        tof_minus_tpion(AliAODTrack *track);
	virtual Double_t        get_mass_squared(AliAODTrack *track);


	Double_t deut_curves[2][2][3];  /* [charge][mean,sigma][par]  */
	Double_t prot_curves[2][2][4];  /* [charge][mean,sigma][par]  */
	
	TF1 *fit_deut_curve = new TF1("fit_d_mean",   "[0] + [1]*x + [2]/sqrt(x)",            1.1, 4.4);
	TF1 *fit_prot_curve = new TF1("fit_p_mean",   "[0] + [1]*x + [2]/sqrt(x) + [3]*x^5",  1.0, 4.4);
	
    private:

        AliAODEvent*            fAOD;               //! input event
        TList*                  fOutputList;        //! output list
	AliPIDResponse*         fPIDResponse;
	

	TH1F*	    fHistPt;             //  1
	TH2F*       cent_ntracks;        //  2
    
	TH1I*       high_per_event_05;   //  3a
	TH1I*       high_per_event_10;   //  3b
	TH1I*       high_per_event_pos;  //  4
	TH1I*       high_per_event_neg;  //  5

	TH2F*       dphi_pt_deut_05;     //  6
	TH2F*       dphi_et_deut_05;     //  7
	TH2F*       dphi_en_deut_05;     //  8
	TH2F*       dphi_kt_deut_05;     //  9
	TH2F*       dphi_pt_prot_05;     // 10
	TH2F*       dphi_et_prot_05;     // 11
	TH2F*       dphi_en_prot_05;     // 12
	TH2F*       dphi_kt_prot_05;     // 13
	TH2F*       dphi_pt_hadr_05;     // 14
	TH2F*       dphi_et_hadr_05;     // 15
	TH2F*       dphi_en_hadr_05;     // 16
	TH2F*       dphi_kt_hadr_05;     // 17
    
//	TH1F*       high_pt_count_05;    // 18
//	TH1F*       high_et_count_05;    // 19
//	TH1F*       high_en_count_05;    // 20
//	TH1F*       high_kt_count_05;    // 21
    
	TH1F*       deut_phi_hist_05;    // 22
	TH1F*       prot_phi_hist_05;    // 23
	TH1F*       hadr_phi_hist_05;    // 24
	TH1F*       high_phi_hist_05;    // 25


	TH2F*       dphi_pt_deut_10;     // 26
	TH2F*       dphi_et_deut_10;     // 27
	TH2F*       dphi_en_deut_10;     // 28
	TH2F*       dphi_kt_deut_10;     // 29    
	TH2F*       dphi_pt_prot_10;     // 30
	TH2F*       dphi_et_prot_10;     // 31
	TH2F*       dphi_en_prot_10;     // 32
	TH2F*       dphi_kt_prot_10;     // 33
	TH2F*       dphi_pt_hadr_10;     // 34
	TH2F*       dphi_et_hadr_10;     // 35
	TH2F*       dphi_en_hadr_10;     // 36
	TH2F*       dphi_kt_hadr_10;     // 37
    
//	TH1F*       high_pt_count_10;    // 38
//	TH1F*       high_et_count_10;    // 39
//	TH1F*       high_en_count_10;    // 40
//	TH1F*       high_kt_count_10;    // 41
    
	TH1F*       deut_phi_hist_10;    // 42
	TH1F*       prot_phi_hist_10;    // 43
	TH1F*       hadr_phi_hist_10;    // 44
	TH1F*       high_phi_hist_10;    // 45

	TH2F*       m2_cut;              // 46

	
        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

//}  //// namespace
#endif
