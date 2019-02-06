/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFdeut_H
#define AliAnalysisTaskCorPIDTOFdeut_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"

class AliAODTrack;

//namespace BSchaefer_devel{
    
class AliAnalysisTaskCorPIDTOFdeut : public AliAnalysisTaskSE  
{
    public:
                              AliAnalysisTaskCorPIDTOFdeut();
                              AliAnalysisTaskCorPIDTOFdeut(const char *name);
        virtual              ~AliAnalysisTaskCorPIDTOFdeut();

        virtual void          UserCreateOutputObjects();
        virtual void          UserExec( Option_t* option);
        virtual void          Terminate(Option_t* option);
	virtual Double_t      Beta(            AliAODTrack *track);
	virtual Double_t      tof_minus_tpion( AliAODTrack *track);
	virtual Double_t      get_mass_squared(AliAODTrack *track);
	virtual Double_t      get_deut_tof_pt( AliAODTrack *track);
	virtual Double_t      get_deut_tof_p(  AliAODTrack *track);

	Double_t deut_curves[2][2][3];  /* [charge][mean,sigma][par]  */
	TF1 *fit_deut_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x)",  1.0, 4.4);

	Double_t cut_width  = 2.0;
	short    do_lead_only = 1;     // 0 = all,  1 = only leading

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



	TH2F*                 m2_pt_pos_cut_T;             // 15
	TH2F*                 m2_pt_neg_cut_T;             // 16
	TH2F*                 m2_pt_pos_cut_G;             // 17
	TH2F*                 m2_pt_neg_cut_G;             // 18
	TH2F*                 m2_pt_pos_cut_A;             // 19
	TH2F*                 m2_pt_neg_cut_A;             // 20
	TH2F*                 m2_pt_pos_cut_B;             // 21
	TH2F*                 m2_pt_neg_cut_B;             // 22
	
	TH2F*                 deut_phi_pt_pos_T;           // 23
	TH2F*                 deut_phi_pt_neg_T;           // 24
	TH2F*                 deut_phi_pt_pos_A;           // 25
	TH2F*                 deut_phi_pt_neg_A;           // 26
	TH2F*                 deut_phi_pt_pos_B;           // 27
	TH2F*                 deut_phi_pt_neg_B;           // 28

	TH1I*                 deut_per_event;              // 29
	TH1I*                 trig_03_per_event;           // 30
	TH1I*                 trig_05_per_event;           // 31
	TH1I*                 trig_08_per_event;           // 32
		
 	TH2F*                 trig_03_phi_pt_pos;          // 33
	TH2F*                 trig_03_phi_pt_neg;          // 34
 	TH2F*                 trig_05_phi_pt_pos;          // 35
	TH2F*                 trig_05_phi_pt_neg;          // 36
	TH2F*                 trig_08_phi_pt_pos;          // 37
	TH2F*                 trig_08_phi_pt_neg;          // 38

	TH2F*                 tof_phi_eta_pos;             // 39
	TH2F*                 tof_phi_eta_neg;             // 40
	TH2F*                 tof_phi_eta_pos_deut;        // 41
	TH2F*                 tof_phi_eta_neg_deut;        // 42
	
	TH1F*                 deut_pt_compare_pos;         // 43
	TH1F*                 deut_pt_compare_neg;         // 44
	TH1F*                 tpc_sector_fraction;         // 45
	TH1F*                 primary_vertex_z;            // 46
	TH1F*                 primary_vertex_z_cut;        // 47
	
	TH2F*                 deut_dphi_pt_pos_pos_03_T;   // 48
	TH2F*                 deut_dphi_pt_pos_neg_03_T;   // 49
	TH2F*                 deut_dphi_pt_neg_neg_03_T;   // 50
	TH2F*                 deut_dphi_pt_pos_pos_05_T;   // 51
	TH2F*                 deut_dphi_pt_pos_neg_05_T;   // 52
	TH2F*                 deut_dphi_pt_neg_neg_05_T;   // 53
	TH2F*                 deut_dphi_pt_pos_pos_08_T;   // 54
	TH2F*                 deut_dphi_pt_pos_neg_08_T;   // 55
	TH2F*                 deut_dphi_pt_neg_neg_08_T;   // 56

	TH2F*                 deut_dphi_pt_pos_pos_03_A;   // 57
	TH2F*                 deut_dphi_pt_pos_neg_03_A;   // 58
	TH2F*                 deut_dphi_pt_neg_neg_03_A;   // 59
	TH2F*                 deut_dphi_pt_pos_pos_05_A;   // 60
	TH2F*                 deut_dphi_pt_pos_neg_05_A;   // 61
	TH2F*                 deut_dphi_pt_neg_neg_05_A;   // 62
	TH2F*                 deut_dphi_pt_pos_pos_08_A;   // 63
	TH2F*                 deut_dphi_pt_pos_neg_08_A;   // 64
	TH2F*                 deut_dphi_pt_neg_neg_08_A;   // 65

	TH2F*                 deut_dphi_pt_pos_pos_03_B;   // 66
	TH2F*                 deut_dphi_pt_pos_neg_03_B;   // 67
	TH2F*                 deut_dphi_pt_neg_neg_03_B;   // 68
	TH2F*                 deut_dphi_pt_pos_pos_05_B;   // 69
	TH2F*                 deut_dphi_pt_pos_neg_05_B;   // 70
	TH2F*                 deut_dphi_pt_neg_neg_05_B;   // 71
	TH2F*                 deut_dphi_pt_pos_pos_08_B;   // 72
	TH2F*                 deut_dphi_pt_pos_neg_08_B;   // 73
	TH2F*                 deut_dphi_pt_neg_neg_08_B;   // 74

	TH1F*                 DCAxy_pos;                   // 75
	TH1F*                 DCAz_pos;                    // 76
	TH1F*                 DCAxy_neg;                   // 77
	TH1F*                 DCAz_neg;                    // 78

	TH2F*                 m2_pt_pos_fine;              // 79
	TH2F*                 m2_pt_neg_fine;              // 80
	TH2F*                 m2_pt_pos_cut_fine;          // 81
	TH2F*                 m2_pt_neg_cut_fine;          // 82
	TH2F*                 m2_pt_pos_cut_T_fine;        // 83
	TH2F*                 m2_pt_neg_cut_T_fine;        // 84
	
        AliAnalysisTaskCorPIDTOFdeut(const AliAnalysisTaskCorPIDTOFdeut&);                        // not implemented
        AliAnalysisTaskCorPIDTOFdeut& operator=(const AliAnalysisTaskCorPIDTOFdeut&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFdeut, 1);
};

//}  //// namespace
#endif
