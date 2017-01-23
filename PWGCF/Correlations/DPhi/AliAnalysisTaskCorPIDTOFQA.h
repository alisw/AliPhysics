/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorPIDTOFQA_H
#define AliAnalysisTaskCorPIDTOFQA_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"

class AliAODTrack;

class AliAnalysisTaskCorPIDTOFQA : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskCorPIDTOFQA();
                                AliAnalysisTaskCorPIDTOFQA(const char *name);
        virtual                 ~AliAnalysisTaskCorPIDTOFQA();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
	virtual Double_t        Beta(AliAODTrack *track);
	virtual Double_t        tof_minus_tpion(AliAODTrack *track);
	virtual Double_t        get_mass_squared(AliAODTrack *track);
	
    private:

        AliAODEvent*            fAOD;               //! input event
        TList*                  fOutputList;        //! output list
	AliPIDResponse*         fPIDResponse;       //! PID response object          //// added by Brennan
	TH1F*                   fHistPt;            //! pt histogram

	TH2F*                   cent_ntracks;       //!                              //// added by Brennan
	TH2F*                   m_squared_pos_raw;  //!                              //// added by Brennan
	TH2F*                   m_squared_pos_cut;  //!                              //// added by Brennan
	TH2F*                   m_squared_pos;      //!                              //// added by Brennan
	TH2F*                   m_squared_neg_raw;  //!                              //// added by Brennan
	TH2F*                   m_squared_neg_cut;  //!                              //// added by Brennan
	TH2F*                   m_squared_neg;      //!                              //// added by Brennan
//	TH1F*                   species_num;        //!                              //// added by Brennan
	TH2F*                   beta_vs_mom;        //!                              //// added by Brennan
//	TH2F*                   sigma_vs_p;         //!                              //// added by Brennan
	TH2F*                   dtof_dEdx;          //!                              //// added by Brennan
	TH1F*                   mapped_ttof;        //!                              //// added by Brennan
	TH2F*                   m_squared_pos_deut; //!                              //// added by Brennan
	TH2F*                   m_squared_neg_deut; //!                              //// added by Brennan
	TH2F*                   plength_vs_mom;     //!                              //// added by Brennan
	TH2F*                   ttof_vs_mom;        //!                              //// added by Brennan
	TH2F*                   deltat_vs_mom;      //!                              //// added by Brennan
	
        AliAnalysisTaskCorPIDTOFQA(const AliAnalysisTaskCorPIDTOFQA&);                        // not implemented
        AliAnalysisTaskCorPIDTOFQA& operator=(const AliAnalysisTaskCorPIDTOFQA&);             // not implemented

        ClassDef(AliAnalysisTaskCorPIDTOFQA, 1);
};

#endif
