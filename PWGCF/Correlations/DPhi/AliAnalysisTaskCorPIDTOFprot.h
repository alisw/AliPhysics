/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFprot_H
#define AliAnalysisTaskCorPIDTOFprot_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"

class AliAODTrack;
class AliPIDResponse;


//namespace BSchaefer_devel{
    
class AliAnalysisTaskCorPIDTOFprot : public AliAnalysisTaskSE  
{
    public:
                              AliAnalysisTaskCorPIDTOFprot();
                              AliAnalysisTaskCorPIDTOFprot(const char *name);
        virtual              ~AliAnalysisTaskCorPIDTOFprot();

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
	
	TH1F*                 fHistPt;                     //  1
	TH2F*                 cent_ntracks;                //  2

	TH2F*                 m2_pt_pos_cut_T_prot;        // 19
	TH2F*                 m2_pt_neg_cut_T_prot;        // 20

	TH2F*                 prot_phi_pt_pos_T;           // 25
	TH2F*                 prot_phi_pt_neg_T;           // 26

	TH1I*                 prot_per_event;              // 36

	TH2F*                 tof_phi_eta_pos;             // 40
	TH2F*                 tof_phi_eta_neg;             // 41

	TH1F*                 tpc_sector_fraction;         // 46
	TH1F*                 primary_vertex_z;            // 47
	TH1F*                 primary_vertex_z_cut1;       // 48
	TH1F*                 primary_vertex_z_cut2;       // 49

	TH2F*                 m2_pt_pos_fine;              // 72
	TH2F*                 m2_pt_neg_fine;              // 73

	TH2F*                 m2_pt_pos_TPC_prot_fine;     // 86
	TH2F*                 m2_pt_neg_TPC_prot_fine;     // 87
	TH2F*                 m2_pt_pos_cut_T_prot_fine;   // 88
	TH2F*                 m2_pt_neg_cut_T_prot_fine;   // 89

	TH1F*                 prot_prot_0509_same;
	TH1F*                 prot_prot_0509_diff;
	TH1F*                 prot_prot_1018_same;
	TH1F*                 prot_prot_1018_diff;


        AliAnalysisTaskCorPIDTOFprot(const AliAnalysisTaskCorPIDTOFprot&);                        // not implemented
        AliAnalysisTaskCorPIDTOFprot& operator=(const AliAnalysisTaskCorPIDTOFprot&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFprot, 1);
};

//}  //// namespace
#endif
