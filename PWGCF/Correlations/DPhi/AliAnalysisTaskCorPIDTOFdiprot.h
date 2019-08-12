/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFdiprot_H
#define AliAnalysisTaskCorPIDTOFdiprot_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"

class AliAODTrack;

//namespace BSchaefer_devel{
    
class AliAnalysisTaskCorPIDTOFdiprot : public AliAnalysisTaskSE  
{
    public:
                              AliAnalysisTaskCorPIDTOFdiprot();
                              AliAnalysisTaskCorPIDTOFdiprot(const char *name);
        virtual              ~AliAnalysisTaskCorPIDTOFdiprot();

        virtual void          UserCreateOutputObjects();
        virtual void          UserExec(Option_t* option);
        virtual void          Terminate(Option_t* option);
	virtual Double_t      Beta(AliAODTrack *track);
	virtual Double_t      tof_minus_tpion(AliAODTrack *track);
	virtual Double_t      get_mass_squared(AliAODTrack *track);


	Double_t prot_curves[2][2][4];  /* [charge][mean,sigma][par]  */
	TF1 *fit_prot_curve = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x) + [3]*x^5",  1.0, 4.0);

	Double_t cut_width  = 1.0;
	
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

	TH1I*                 prot_per_event;              // 15
//	TH1I*                 prot_per_event_pos;          // 16
//	TH1I*                 prot_per_event_neg;  	   // 17
	TH2F*                 m2_pt_pos_cut_T;             // 18
	TH2F*                 m2_pt_neg_cut_T;             // 19

	TH2F*                 prot_phi_pt_pos;             // 20
	TH2F*                 prot_phi_pt_neg;             // 21
	
	TH2F*                 prot_p0_pt_pos_pos;          // 22
	TH2F*                 prot_p0_pt_pos_neg;          // 23
	TH2F*                 prot_p0_pt_neg_neg;          // 24

	TH2F*                 di_prot_phi_pt;              // 25
	TH2F*                 trig_05_phi_pt;              // 26
	TH2F*                 trig_08_phi_pt;              // 27

	TH2F*                 di_prot_dphi_p0_pos_pos_05;  // 28
	TH2F*                 di_prot_dphi_p0_pos_neg_05;  // 29
	TH2F*                 di_prot_dphi_p0_neg_neg_05;  // 30

	TH2F*                 di_prot_dphi_p0_pos_pos_08;  // 31
	TH2F*                 di_prot_dphi_p0_pos_neg_08;  // 32
	TH2F*                 di_prot_dphi_p0_neg_neg_08;  // 33

	TH1F*                 primary_vertex_z;            // 47
	TH1F*                 primary_vertex_z_cut1;       // 48
	TH1F*                 primary_vertex_z_cut2;       // 49
	
	TH1F*  di_prot_dphi_r02_pos_pos_05; // 28
	TH1F*  di_prot_dphi_r02_pos_neg_05; // 29
	TH1F*  di_prot_dphi_r02_neg_neg_05; // 30

	TH1F*  di_prot_dphi_r02_pos_pos_08; // 31
	TH1F*  di_prot_dphi_r02_pos_neg_08; // 32
	TH1F*  di_prot_dphi_r02_neg_neg_08; // 33

	TH1F*  di_prot_dphi_r04_pos_pos_05; // 28
	TH1F*  di_prot_dphi_r04_pos_neg_05; // 29
	TH1F*  di_prot_dphi_r04_neg_neg_05; // 30

	TH1F*  di_prot_dphi_r04_pos_pos_08; // 31
	TH1F*  di_prot_dphi_r04_pos_neg_08; // 32
	TH1F*  di_prot_dphi_r04_neg_neg_08; // 33
	
//	TH2F*                 track_cor_radius_pt;         // 34
//	TH2F*                 track_cor_radius_pt_cut;     // 35
	
        AliAnalysisTaskCorPIDTOFdiprot(const AliAnalysisTaskCorPIDTOFdiprot&);                        // not implemented
        AliAnalysisTaskCorPIDTOFdiprot& operator=(const AliAnalysisTaskCorPIDTOFdiprot&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFdiprot, 1);
};

//}  //// namespace
#endif
