#ifndef ALIANALYSISTASKEMCALJETENERGYFLOW_H
#define ALIANALYSISTASKEMCALJETENERGYFLOW_H
/**
 * \file AliAnalysisTaskEmcalJetEnergyFlow.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetEnergyFlow
 *
 * In this header file the class AliAnalysisTaskEmcalJetEnergyFlow is declared
 * This is an analysis task which calculates the energy flow between jets of increasing jet radii, defined as the difference of pt in spatially matched jets.
 * The goal of this task is to define an observable for the effect of the recoil of medium particles from the jet equivalent to the feature modeled to Jewel generator.
 *
 * \author Christos Pliatskas <c.pliatskas.stylianidis@nikhef.nl>, Nikhef
 * \date Oct 26, 2021
 */

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved.*
 * See cxx source for full copyright notice */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskEmcalJetEnergyFlow
 * \brief Implementation of analysis task which calculates the energy flow between jets of increasing jet radii. It derives from AliAnalysisTaskEmcalJet.
 *  After performing a geometrical matching for the jets of the different radii, it calculates the pt difference between subsequent pairs of them. Additionally some basic analysis on track and jet spectra is performed.
 */

class AliAnalysisTaskEmcalJetEnergyFlow: public AliAnalysisTaskEmcalJet {
	public:
        enum AnalysisType{
                kppData         = 0,
                kppMC           = 1,
                kPbPbData       = 2,
                kEmbedded       = 3 
                };
	
	AliAnalysisTaskEmcalJetEnergyFlow()			;
	AliAnalysisTaskEmcalJetEnergyFlow(const char* name)	;
	virtual ~AliAnalysisTaskEmcalJetEnergyFlow()		;
	
	void 	UserCreateOutputObjects()				;
	void	Terminate(Option_t* option)				;
	
	static  AliAnalysisTaskEmcalJetEnergyFlow* AddTaskEmcalJetEnergyFlow(
		const char *ntracks		= "usedefault",
		const char *nclusters		= "usedefault",
		const char *ncells		= "usedefault",
		Double_t Rstep_EF               = 0.1,              
                Double_t Max_match_dr           = 0.2,
                Double_t Lead_pt_cut            = 0.0,
                AnalysisType fAnType            =kppData,
                const char *suffix              = "" );
	
	protected:


	void			JetMatcher(const TList *genJetsList,const Int_t &kGenJets,
   					   const TList *recJetsList,const Int_t &kRecJets,
   					   TArrayI &iGenIndex,TArrayI &iRecIndex,
   					   Int_t iDebug = 0,Float_t maxDist = 0.3,Float_t max_eta = 0.5);	///< This is a re-implementation of the AliAnalysisHelperJetTasks::GetClosestJets adjusted for AliEmcalJet instead of AliAODjets for the function's inputs.
	void			ExecOnce()				;
	Bool_t			FillHistograms()			;
	Bool_t			Run()					;
        void                    SetAnalysisType(AnalysisType a){fAnalysisType = a;}
        AnalysisType            GetAnalysisType(){return fAnalysisType;}
        void                    SetMaxMatchDR(Double_t dr){Max_match_dist = dr;}
        Double_t                SetMaxMatchDR(){return Max_match_dist;}
	void			AllocateJetHistograms()			;
	void			AllocateTrackHistograms()		; ///<Same as Sample task
	void                    AllocateClusterHistograms()             ; ///<May remove later
  	void                    AllocateCellHistograms()                ; ///<May remove later
	void			AllocateEnergyflowHistograms()		;

	void                    DoJetLoop()                             ;
	void                    DoTrackLoop()                           ;
	void			FillEFHistograms()			;
	void                    DoClusterLoop()                         ; ///<May remove later	
	void                    DoCellLoop()                            ; ///<May remove later

        Double_t                LeadPtCut                               ;///<Pt cut on the jet's leading track
	Double_t                R_jet_step				;///<Radial step for the dpt calculation
        Double_t                Max_match_dist                          ;///<Maximum distance for the matching between Rjet
	THistManager            fHistManager                            ;///<Hist manager
//	TList*			fOutput					;///!<! Output list
 	private:
        AnalysisType            fAnalysisType                           ;///<Flag for type of analysis
 	 AliAnalysisTaskEmcalJetEnergyFlow(const AliAnalysisTaskEmcalJetEnergyFlow&); // not implemented
  	 AliAnalysisTaskEmcalJetEnergyFlow &operator=(const AliAnalysisTaskEmcalJetEnergyFlow&); // not implemented

  	  ClassDef(AliAnalysisTaskEmcalJetEnergyFlow,20);
	/// \endcond
};
#endif
