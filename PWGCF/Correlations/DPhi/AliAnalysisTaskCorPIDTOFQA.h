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
	Double_t prot_curves[2][2][4];  /* [charge][mean,sigma][par]  */
	TF1 *fit_prot_curve   = new TF1("fit_m_prot",   "[0] + [1]*x + [2]/sqrt(x) + [3]*x^5", 1.0, 4.0);
	TF1 *fit_deut_curve   = new TF1("fit_m_deut",   "[0] + [1]*x + [2]/x^2 + [3]*x^2 + [4]*x^3 + [5]/x^4",      0.8, 4.5);
	

	short    run_mode     = 0;

	Float_t pio2   = TMath::PiOver2();
	Float_t twopi  = TMath::TwoPi();

	int      eta_limiting = 1;
	Double_t cut_width    = 2.0;
	short    do_lead_only = 0;     // 0 = all,  1 = only leading

	
    private:    

        AliAODEvent*          fAOD;               //! input event
        TList*                fOutputList;        //! output list
	AliPIDResponse*       fPIDResponse;
	AliAnalysisUtils*     fAnalysisUtils;

	TH1F* deut_eta_pos_pt1;
	TH1F* deut_eta_neg_pt1;
	TH1F* deut_eta_pos_pt2;
	TH1F* deut_eta_neg_pt2;
	TH1F* deut_eta_pos_pt3;
	TH1F* deut_eta_neg_pt3;
	TH1F* deut_eta_pos_pt4;
	TH1F* deut_eta_neg_pt4;
	TH1F* deut_eta_pos_pt5;
	TH1F* deut_eta_neg_pt5;
	
	TH1F* hadr_eta_pos;
	TH1F* hadr_eta_neg;

	TH2F* deut_hadr_dphi36_deta20_pt1;
	TH2F* deut_hadr_dphi36_deta20_pt2;
	TH2F* deut_hadr_dphi36_deta20_pt3;
	TH2F* deut_hadr_dphi36_deta20_pt4;
	TH2F* deut_hadr_dphi36_deta20_pt5;

	
	TH2F* deut_hadr_dphi36_deta21_pt1;
	TH2F* deut_hadr_dphi36_deta21_pt2;
	TH2F* deut_hadr_dphi36_deta21_pt3;
	TH2F* deut_hadr_dphi36_deta21_pt4;
	TH2F* deut_hadr_dphi36_deta21_pt5;
	
	TH1I* multiplicity_w_deut_1;                       //  n1
	TH1I* multiplicity_w_deut_2;                       //  n2
	TH1I* multiplicity_w_deut_3;                       //  n3
	TH1I* multiplicity_w_deut_4;                       //  n4
	TH1I* multiplicity_w_deut_5;                       //  n5

	TH1I* multiplicity_w_trig_and_deut_1;              //  n6
	TH1I* multiplicity_w_trig_and_deut_2;              //  n7
	TH1I* multiplicity_w_trig_and_deut_3;              //  n8
	TH1I* multiplicity_w_trig_and_deut_4;              //  n9
	TH1I* multiplicity_w_trig_and_deut_5;              //  n10
    
	
	TH1F*                 fHistPt;                     //  1
	TH2F*                 cent_ntracks;                //  2
	
//	TH2F*                 m2_pt_pos;                   //  3
//	TH2F*                 m2_pt_neg;                   //  4
	TH2F*                 m2_pt_pos_TPC;               //  5
	TH2F*                 m2_pt_neg_TPC;	           //  6
	TH2F*                 m2_pt_pos_cut_T;             //  7
	TH2F*                 m2_pt_neg_cut_T;             //  8
	TH2F*                 m2_pt_pos_cut_G;             //  9
	TH2F*                 m2_pt_neg_cut_G;             // 10

	TH2F*                 m2_pt_pos_cut_sa0;           //
	TH2F*                 m2_pt_neg_cut_sa0;           //
	TH2F*                 m2_pt_pos_cut_sa1;           //
	TH2F*                 m2_pt_neg_cut_sa1;           //
	TH2F*                 m2_pt_pos_cut_sa2;           //
	TH2F*                 m2_pt_neg_cut_sa2;           //
	TH2F*                 m2_pt_pos_cut_sa3;           //
	TH2F*                 m2_pt_neg_cut_sa3;           //
	TH2F*                 m2_pt_pos_cut_sa4;           //
	TH2F*                 m2_pt_neg_cut_sa4;           //
	TH2F*                 m2_pt_pos_cut_sa5;           //
	TH2F*                 m2_pt_neg_cut_sa5;           //
	TH2F*                 m2_pt_pos_cut_sa6;           //
	TH2F*                 m2_pt_neg_cut_sa6;           //
	TH2F*                 m2_pt_pos_cut_sa7;           //
	TH2F*                 m2_pt_neg_cut_sa7;           //
	
	TH2F*                 m2_pt_pos_cut_sb0;           //
	TH2F*                 m2_pt_neg_cut_sb0;           //
	TH2F*                 m2_pt_pos_cut_sb1;           //
	TH2F*                 m2_pt_neg_cut_sb1;           //
	TH2F*                 m2_pt_pos_cut_sb2;           //
	TH2F*                 m2_pt_neg_cut_sb2;           //
	TH2F*                 m2_pt_pos_cut_sb3;           //
	TH2F*                 m2_pt_neg_cut_sb3;           //
	TH2F*                 m2_pt_pos_cut_sb4;           //
	TH2F*                 m2_pt_neg_cut_sb4;           //
	TH2F*                 m2_pt_pos_cut_sb5;           //
	TH2F*                 m2_pt_neg_cut_sb5;           //
	TH2F*                 m2_pt_pos_cut_sb6;           //
	TH2F*                 m2_pt_neg_cut_sb6;           //
	TH2F*                 m2_pt_pos_cut_sb7;           //
	TH2F*                 m2_pt_neg_cut_sb7;           //

	TH2F*                 m2_pt_pos_cut_T_prot;        // 19
	TH2F*                 m2_pt_neg_cut_T_prot;        // 20
	
	TH2F*                 m2_pt_pos_cut_with_trig_05;  // 21
	TH2F*                 m2_pt_neg_cut_with_trig_05;  // 22
	
	TH2F*                 deut_phi_pt_pos_T;           // 23
	TH2F*                 deut_phi_pt_neg_T;           // 24
	
	TH2F*                 prot_phi_pt_pos_T;           // 25
	TH2F*                 prot_phi_pt_neg_T;           // 26

	TH2F*                 deut_phi_pt_pos_sa0;         //
	TH2F*                 deut_phi_pt_neg_sa0;         //
	TH2F*                 deut_phi_pt_pos_sa1;         //
	TH2F*                 deut_phi_pt_neg_sa1;         //
	TH2F*                 deut_phi_pt_pos_sa2;         //
	TH2F*                 deut_phi_pt_neg_sa2;         //
	TH2F*                 deut_phi_pt_pos_sa3;         //
	TH2F*                 deut_phi_pt_neg_sa3;         //
	TH2F*                 deut_phi_pt_pos_sa4;         //
	TH2F*                 deut_phi_pt_neg_sa4;         //
	TH2F*                 deut_phi_pt_pos_sa5;         //
	TH2F*                 deut_phi_pt_neg_sa5;         //
	TH2F*                 deut_phi_pt_pos_sa6;         //
	TH2F*                 deut_phi_pt_neg_sa6;         //
	TH2F*                 deut_phi_pt_pos_sa7;         //
	TH2F*                 deut_phi_pt_neg_sa7;         //
	
	TH2F*                 deut_phi_pt_pos_sb0;         //
	TH2F*                 deut_phi_pt_neg_sb0;         //
	TH2F*                 deut_phi_pt_pos_sb1;         //
	TH2F*                 deut_phi_pt_neg_sb1;         //
	TH2F*                 deut_phi_pt_pos_sb2;         //
	TH2F*                 deut_phi_pt_neg_sb2;         //
	TH2F*                 deut_phi_pt_pos_sb3;         //
	TH2F*                 deut_phi_pt_neg_sb3;         //
	TH2F*                 deut_phi_pt_pos_sb4;         //
	TH2F*                 deut_phi_pt_neg_sb4;         //
	TH2F*                 deut_phi_pt_pos_sb5;         //
	TH2F*                 deut_phi_pt_neg_sb5;         //
	TH2F*                 deut_phi_pt_pos_sb6;         //
	TH2F*                 deut_phi_pt_neg_sb6;         //
	TH2F*                 deut_phi_pt_pos_sb7;         //
	TH2F*                 deut_phi_pt_neg_sb7;         //

	TH1I*                 deut_per_event;              // 35
	TH1I*                 prot_per_event;              // 36
	TH1I*                 deut_per_event_with_trig;    // 35
	TH1I*                 trig_05_per_event;           // 37

 	TH2F*                 trig_05_phi_pt_pos;          // 38
	TH2F*                 trig_05_phi_pt_neg;          // 39

	TH2F*                 tof_phi_eta_pos;             // 40
	TH2F*                 tof_phi_eta_neg;             // 41
	TH2F*                 tof_phi_eta_pos_deut;        // 42
	TH2F*                 tof_phi_eta_neg_deut;        // 43
	
	TH1F*                 deut_pt_compare_pos;         // 44
	TH1F*                 deut_pt_compare_neg;         // 45
	TH1F*                 tpc_sector_fraction;         // 46
	TH1F*                 primary_vertex_z;            // 47
	TH1F*                 primary_vertex_z_cut1;       // 48
	TH1F*                 primary_vertex_z_cut2;       // 49
	
	TH2F*                 deut_dphi_pt_pos_pos_05_T;   // 50
	TH2F*                 deut_dphi_pt_pos_neg_05_T;   // 51
	TH2F*                 deut_dphi_pt_neg_neg_05_T;   // 52

	TH2F*                 prot_dphi_pt_pos_pos_05_T;   // 53
	TH2F*                 prot_dphi_pt_pos_neg_05_T;   // 54
	TH2F*                 prot_dphi_pt_neg_neg_05_T;   // 55

	
	TH2F*                 deut_dphi_pt_pos_pos_05_sa0; // 62
	TH2F*                 deut_dphi_pt_pos_neg_05_sa0; // 63
	TH2F*                 deut_dphi_pt_neg_neg_05_sa0; // 64
	
	TH2F*                 deut_dphi_pt_pos_pos_05_sa1; // 65
	TH2F*                 deut_dphi_pt_pos_neg_05_sa1; // 66
	TH2F*                 deut_dphi_pt_neg_neg_05_sa1; // 67

	TH2F*                 deut_dphi_pt_pos_pos_05_sa2; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sa2; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sa2; //

	TH2F*                 deut_dphi_pt_pos_pos_05_sa3; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sa3; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sa3; //
	
	TH2F*                 deut_dphi_pt_pos_pos_05_sa4; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sa4; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sa4; //
	
	TH2F*                 deut_dphi_pt_pos_pos_05_sa5; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sa5; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sa5; //

	TH2F*                 deut_dphi_pt_pos_pos_05_sa6; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sa6; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sa6; //
	
	TH2F*                 deut_dphi_pt_pos_pos_05_sa7; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sa7; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sa7; //
	
	
	TH2F*                 deut_dphi_pt_pos_pos_05_sb0; // 62
	TH2F*                 deut_dphi_pt_pos_neg_05_sb0; // 63
	TH2F*                 deut_dphi_pt_neg_neg_05_sb0; // 64
	
	TH2F*                 deut_dphi_pt_pos_pos_05_sb1; // 65
	TH2F*                 deut_dphi_pt_pos_neg_05_sb1; // 66
	TH2F*                 deut_dphi_pt_neg_neg_05_sb1; // 67

	TH2F*                 deut_dphi_pt_pos_pos_05_sb2; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sb2; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sb2; //

	TH2F*                 deut_dphi_pt_pos_pos_05_sb3; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sb3; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sb3; //
	
	TH2F*                 deut_dphi_pt_pos_pos_05_sb4; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sb4; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sb4; //
	
	TH2F*                 deut_dphi_pt_pos_pos_05_sb5; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sb5; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sb5; //

	TH2F*                 deut_dphi_pt_pos_pos_05_sb6; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sb6; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sb6; //
	
	TH2F*                 deut_dphi_pt_pos_pos_05_sb7; //
	TH2F*                 deut_dphi_pt_pos_neg_05_sb7; //
	TH2F*                 deut_dphi_pt_neg_neg_05_sb7; //
	
	TH1F*                 DCAxy_pos;                   // 68
	TH1F*                 DCAz_pos;                    // 69
	TH1F*                 DCAxy_neg;                   // 70
	TH1F*                 DCAz_neg;                    // 71

	TH2F*                 m2_pt_pos_fine;              // 72
	TH2F*                 m2_pt_neg_fine;              // 73
	TH2F*                 m2_pt_pos_TPC_fine;          // 74
	TH2F*                 m2_pt_neg_TPC_fine;          // 75
	TH2F*                 m2_pt_pos_cut_T_fine;        // 76
	TH2F*                 m2_pt_neg_cut_T_fine;        // 77

	TH2F*                 m2_pt_pos_cut_sa0_fine;      // 82
	TH2F*                 m2_pt_neg_cut_sa0_fine;      // 83
	TH2F*                 m2_pt_pos_cut_sa1_fine;      // 84
	TH2F*                 m2_pt_neg_cut_sa1_fine;      // 85
	TH2F*                 m2_pt_pos_cut_sa2_fine;      //
	TH2F*                 m2_pt_neg_cut_sa2_fine;      //
	TH2F*                 m2_pt_pos_cut_sa3_fine;      //
	TH2F*                 m2_pt_neg_cut_sa3_fine;      //
	TH2F*                 m2_pt_pos_cut_sa4_fine;      //
	TH2F*                 m2_pt_neg_cut_sa4_fine;      //
	TH2F*                 m2_pt_pos_cut_sa5_fine;      //
	TH2F*                 m2_pt_neg_cut_sa5_fine;      //
	TH2F*                 m2_pt_pos_cut_sa6_fine;      //
	TH2F*                 m2_pt_neg_cut_sa6_fine;      //
	TH2F*                 m2_pt_pos_cut_sa7_fine;      //
	TH2F*                 m2_pt_neg_cut_sa7_fine;      //
	
	TH2F*                 m2_pt_pos_cut_sb0_fine;      // 82
	TH2F*                 m2_pt_neg_cut_sb0_fine;      // 83
	TH2F*                 m2_pt_pos_cut_sb1_fine;      // 84
	TH2F*                 m2_pt_neg_cut_sb1_fine;      // 85
	TH2F*                 m2_pt_pos_cut_sb2_fine;      //
	TH2F*                 m2_pt_neg_cut_sb2_fine;      //
	TH2F*                 m2_pt_pos_cut_sb3_fine;      //
	TH2F*                 m2_pt_neg_cut_sb3_fine;      //
	TH2F*                 m2_pt_pos_cut_sb4_fine;      //
	TH2F*                 m2_pt_neg_cut_sb4_fine;      //
	TH2F*                 m2_pt_pos_cut_sb5_fine;      //
	TH2F*                 m2_pt_neg_cut_sb5_fine;      //
	TH2F*                 m2_pt_pos_cut_sb6_fine;      //
	TH2F*                 m2_pt_neg_cut_sb6_fine;      //
	TH2F*                 m2_pt_pos_cut_sb7_fine;      //
	TH2F*                 m2_pt_neg_cut_sb7_fine;      //

	TH2F*                 m2_pt_pos_TPC_prot_fine;     // 86
	TH2F*                 m2_pt_neg_TPC_prot_fine;     // 87
	TH2F*                 m2_pt_pos_cut_T_prot_fine;   // 88
	TH2F*                 m2_pt_neg_cut_T_prot_fine;   // 89

	TH2F*                 deut_deut_dphi_deta_same;
	TH2F*                 deut_deut_dphi_deta_diff;
	TH2F*                 deut_deut_dpt_dcharge;
	TH2F*                 deut_deut_ptA_ptB_same;
	TH2F*                 deut_deut_ptA_ptB_diff;

	TH1F*                 deut_deut_inv_mass;

        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

//}  //// namespace
#endif
