/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFQA_H
#define AliAnalysisTaskCorPIDTOFQA_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"

class AliAODTrack;
class AliPIDResponse;


//namespace BSchaefer_devel{
    
class AliAnalysisTaskCorPIDTOFQA : public AliAnalysisTaskSE  
{
    public:
                              AliAnalysisTaskCorPIDTOFQA();
                              AliAnalysisTaskCorPIDTOFQA(const char *name);
        virtual              ~AliAnalysisTaskCorPIDTOFQA();

        virtual void          UserCreateOutputObjects();
        virtual void          UserExec( Option_t* option);
        virtual void          Terminate(Option_t* option);
	virtual Double_t      Beta(            AliAODTrack *track);
	virtual Double_t      tof_minus_tpion( AliAODTrack *track);
	virtual Double_t      get_mass_squared(AliAODTrack *track);
	virtual Double_t      get_deut_tof_pt( AliAODTrack *track);
	virtual Double_t      get_deut_tof_p(  AliAODTrack *track);

	Double_t deut_curves[2][2][6];  /* [charge][mean,sigma][par]  */
//	TF1 *fit_deut_curve   = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x)",  1.0, 4.4);
	TF1 *fit_deut_curve   = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/x^2 + [3]*x^2 + [4]*x^3 + [5]/x^4",      0.8, 4.5);
	
	Double_t cut_width    = 2.0;
	short    do_lead_only = 0;     // 0 = all,  1 = only leading
	short    run_mode     = 0;

	Float_t pio2   = TMath::PiOver2();
	Float_t twopi  = TMath::TwoPi();

    private:    

        AliAODEvent*          fAOD;               //! input event
        TList*                fOutputList;        //! output list
	AliPIDResponse*       fPIDResponse;
	AliAnalysisUtils*     fAnalysisUtils;
	
	TH1F*                 fHistPt;                     //  1
	TH2F*                 cent_ntracks;                //  2
	
	TH2F*                 m2_pt_pos;                   //  3
	TH2F*                 m2_pt_neg;                   //  4
	TH2F*                 m2_pt_pos_TPC;               //  5
	TH2F*                 m2_pt_neg_TPC;	           //  6
	TH2F*                 m2_pt_pos_cut_T;             //  7
	TH2F*                 m2_pt_neg_cut_T;             //  8
	TH2F*                 m2_pt_pos_cut_G;             //  9
	TH2F*                 m2_pt_neg_cut_G;             // 10
	TH2F*                 m2_pt_pos_cut_A;             // 11
	TH2F*                 m2_pt_neg_cut_A;             // 12
	TH2F*                 m2_pt_pos_cut_B;             // 13
	TH2F*                 m2_pt_neg_cut_B;             // 14

	TH2F*                 m2_pt_pos_cut_with_trig_04;  // 15
	TH2F*                 m2_pt_neg_cut_with_trig_04;  // 16
	TH2F*                 m2_pt_pos_cut_with_trig_05;  // 17
	TH2F*                 m2_pt_neg_cut_with_trig_05;  // 18
	TH2F*                 m2_pt_pos_cut_with_trig_06;  // 19
	TH2F*                 m2_pt_neg_cut_with_trig_06;  // 20
	
	TH2F*                 deut_phi_pt_pos_T;           // 21
	TH2F*                 deut_phi_pt_neg_T;           // 22
	TH2F*                 deut_phi_pt_pos_A;           // 23
	TH2F*                 deut_phi_pt_neg_A;           // 24
	TH2F*                 deut_phi_pt_pos_B;           // 25
	TH2F*                 deut_phi_pt_neg_B;           // 26

	TH1I*                 deut_per_event;              // 27
	TH1I*                 trig_04_per_event;           // 28
	TH1I*                 trig_05_per_event;           // 29
	TH1I*                 trig_06_per_event;           // 30
		
 	TH2F*                 trig_04_phi_pt_pos;          // 31
	TH2F*                 trig_04_phi_pt_neg;          // 32
 	TH2F*                 trig_05_phi_pt_pos;          // 33
	TH2F*                 trig_05_phi_pt_neg;          // 34
	TH2F*                 trig_06_phi_pt_pos;          // 35
	TH2F*                 trig_06_phi_pt_neg;          // 36

	TH2F*                 tof_phi_eta_pos;             // 37
	TH2F*                 tof_phi_eta_neg;             // 38
	TH2F*                 tof_phi_eta_pos_deut;        // 39
	TH2F*                 tof_phi_eta_neg_deut;        // 40
	
	TH1F*                 deut_pt_compare_pos;         // 41
	TH1F*                 deut_pt_compare_neg;         // 42
	TH1F*                 tpc_sector_fraction;         // 43
	TH1F*                 primary_vertex_z;            // 44
	TH1F*                 primary_vertex_z_cut1;       // 45
	TH1F*                 primary_vertex_z_cut2;       // 46
	
	TH2F*                 deut_dphi_pt_pos_pos_04_T;   // 47
	TH2F*                 deut_dphi_pt_pos_neg_04_T;   // 48
	TH2F*                 deut_dphi_pt_neg_neg_04_T;   // 49
	TH2F*                 deut_dphi_pt_pos_pos_05_T;   // 50
	TH2F*                 deut_dphi_pt_pos_neg_05_T;   // 51
	TH2F*                 deut_dphi_pt_neg_neg_05_T;   // 52
	TH2F*                 deut_dphi_pt_pos_pos_06_T;   // 53
	TH2F*                 deut_dphi_pt_pos_neg_06_T;   // 54
	TH2F*                 deut_dphi_pt_neg_neg_06_T;   // 55

	TH2F*                 deut_dphi_pt_pos_pos_04_A;   // 56
	TH2F*                 deut_dphi_pt_pos_neg_04_A;   // 57
	TH2F*                 deut_dphi_pt_neg_neg_04_A;   // 58
	TH2F*                 deut_dphi_pt_pos_pos_05_A;   // 59
	TH2F*                 deut_dphi_pt_pos_neg_05_A;   // 60
	TH2F*                 deut_dphi_pt_neg_neg_05_A;   // 61
	TH2F*                 deut_dphi_pt_pos_pos_06_A;   // 62
	TH2F*                 deut_dphi_pt_pos_neg_06_A;   // 63
	TH2F*                 deut_dphi_pt_neg_neg_06_A;   // 64

	TH2F*                 deut_dphi_pt_pos_pos_04_B;   // 65
	TH2F*                 deut_dphi_pt_pos_neg_04_B;   // 66
	TH2F*                 deut_dphi_pt_neg_neg_04_B;   // 67
	TH2F*                 deut_dphi_pt_pos_pos_05_B;   // 68
	TH2F*                 deut_dphi_pt_pos_neg_05_B;   // 69
	TH2F*                 deut_dphi_pt_neg_neg_05_B;   // 70
	TH2F*                 deut_dphi_pt_pos_pos_06_B;   // 71
	TH2F*                 deut_dphi_pt_pos_neg_06_B;   // 72
	TH2F*                 deut_dphi_pt_neg_neg_06_B;   // 73

	TH1F*                 DCAxy_pos;                   // 74
	TH1F*                 DCAz_pos;                    // 75
	TH1F*                 DCAxy_neg;                   // 76
	TH1F*                 DCAz_neg;                    // 77

	TH2F*                 m2_pt_pos_fine;              // 78
	TH2F*                 m2_pt_neg_fine;              // 79
	TH2F*                 m2_pt_pos_TPC_fine;          // 80
	TH2F*                 m2_pt_neg_TPC_fine;          // 81
	TH2F*                 m2_pt_pos_cut_T_fine;        // 82
	TH2F*                 m2_pt_neg_cut_T_fine;        // 83
	TH2F*                 m2_pt_pos_cut_A_fine;        // 84
	TH2F*                 m2_pt_neg_cut_A_fine;        // 85
	TH2F*                 m2_pt_pos_cut_B_fine;        // 86
	TH2F*                 m2_pt_neg_cut_B_fine;        // 87

	TH2F*                 deut_m2_dphi_04_pt1;         // 88
	TH2F*                 deut_m2_dphi_04_pt2;         // 89
	TH2F*                 deut_m2_dphi_04_pt3;         // 90
	TH2F*                 deut_m2_dphi_04_pt4;         // 91
	
        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

//}  //// namespace
#endif
