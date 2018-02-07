/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFQA_H
#define AliAnalysisTaskCorPIDTOFQA_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"
#include <iosfwd>

//#define AliAnalysisTaskCorPIDTOFQA AliAnalysisTaskCorPIDTOFQA2

class AliAODTrack;

//namespace BSchaefer_devel{

//namespace std{
    
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



	Double_t deut_curves[2][2][3];  /* [charge][mean,sigma][par]  */
	TF1 *fit_deut_curve   = new TF1("fit_m_mean",   "[0] + [1]*x + [2]/sqrt(x)",  1.0, 4.4);

	Double_t cut_width    = 2.0;
	short    do_lead_only = 0;     // 0 = all,  1 = only leading
	short    do_aod_skim  = 0;
	
	Float_t pio2          = TMath::PiOver2();
	Float_t twopi         = TMath::TwoPi();

	static const int trig_bins    = 3;
	static const int pt_bins      = 4;
	static const int zvert_bins   = 2;
	static const int mult_bins    = 5;
	static const int physics_bins = 4;
	
	TH2F deut_dphi_deta[trig_bins][pt_bins][zvert_bins][mult_bins][physics_bins];

//	ofstream file("output.txt");;

    private:

        AliAODEvent*          fAOD;               //! input event
        TList*                fOutputList;        //! output list
	AliPIDResponse*       fPIDResponse;
	AliAnalysisUtils*     fAnalysisUtils;
	
	TH1F*                 primary_vertex_z;            //  1
	TH1F*                 primary_vertex_z_cut;        //  2
	
	TH1I*                 multiplicity_hybrid;
	TH1I*                 multiplicity_global;

	TH1I*                 deut_per_event;              //  3
	TH1F*                 fHistPt;                     //  4

	TH2F*                 m2_pt_pos;                   //  5
	TH2F*                 m2_pt_neg;                   //  6
	TH2F*                 m2_pt_pos_cut;               //  7
	TH2F*                 m2_pt_neg_cut;               //  8
	
	TH2F*                 m2_pt_pos_fine;              //  9
	TH2F*                 m2_pt_neg_fine;              // 10
	TH2F*                 m2_pt_pos_cut_fine;          // 11
	TH2F*                 m2_pt_neg_cut_fine;          // 12

	TH2F*                 trig_03_phi_pt_pos;          // 13
	TH2F*                 trig_03_phi_pt_neg;          // 14
	TH2F*                 trig_05_phi_pt_pos;          // 15
	TH2F*                 trig_05_phi_pt_neg;          // 16
	TH2F*                 trig_08_phi_pt_pos;          // 17
	TH2F*                 trig_08_phi_pt_neg;          // 18

	std::ofstream*             associates;
	std::ofstream*             triggers;
	
        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

//}  //// namespace
#endif
