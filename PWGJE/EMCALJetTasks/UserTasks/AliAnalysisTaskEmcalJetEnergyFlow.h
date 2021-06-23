#ifndef ALIANALYSISTASKEMCALJETENERGYFLOW_H
#define ALIANALYSISTASKEMCALJETENERGYFLOW_H
/**
 * \file AliAnalysisTaskEmcalJetEnergyFlow.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetEnergyFlow
 *
 * In this header file the class AliAnalysisTaskEmcalJetEnergyFlow is declared
 * This is an analysis task which calculates the energy flow between jets of increasing jet radii, defined as the difference of pt in matched jets.
 * The goal of this task is to define an observable for the effect of the recoil of medium particles from the jet equivalent to the feature modeled to Jewel generator.
 *
 * \author Christos Pliatskas <c.pliatskas.stylianidis@nikhef.nl>, Nikhef
 * \date Jan 04, 2021
 */

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved.*
 * See cxx source for full copyright notice */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskEmcalJetEnergyFlow
 * \brief Implementation of analysis task which calculates the energy flow between jets of increasing
 *  jet radii. It derives from AliAnalysisTaskEmcalJet.
 *  After performing a geometrical matching for the jets of the different radii, it calculates the pt
 *  difference between subsequent pairs of them. Additionally some basic analysis on track and jet
 *  spectra is performed.
 */

class AliAnalysisTaskEmcalJetEnergyFlow: public AliAnalysisTaskEmcalJet {
	public:
	
	AliAnalysisTaskEmcalJetEnergyFlow()				;
	AliAnalysisTaskEmcalJetEnergyFlow(const char* name)		;
	virtual ~AliAnalysisTaskEmcalJetEnergyFlow()			;
	
	void 	UserCreateOutputObjects()				;
	void	Terminate(Option_t* option)				;
	
	static  AliAnalysisTaskEmcalJetEnergyFlow* AddTaskEmcalJetEnergyFlow(
		const char *ntracks		= "usedefault",
		const char *nclusters		= "usedefault",
		const char *ncells		= "usedefault",
		Bool_t SetMCprod                 = kTRUE,
                const char *suffix              = "" );
	
	protected:


	void			JetMatcher(const TList *genJetsList,const Int_t &kGenJets,
   					   const TList *recJetsList,const Int_t &kRecJets,
   					   TArrayI &iGenIndex,TArrayI &iRecIndex,
   					   Int_t iDebug = 0,Float_t maxDist = 0.3,Float_t max_eta = 0.5);	// /< This is a re-implementation of the AliAnalysisHelperJetTasks::GetClosestJets adjusted for AliEmcalJet instead of AliAODjets for the function's inputs.
	void			ExecOnce()				;
	Bool_t			FillHistograms()			;
	Bool_t			Run()					;

	void			AllocateJetHistograms()			;
	void			AllocateTrackHistograms()		; // /<Same as Sample task
	void                    AllocateClusterHistograms()             ; // /<May remove later
  	void                    AllocateCellHistograms()                ; // /<May remove later
	void			AllocateEnergyflowHistograms()		;

	void                    DoJetLoop()                             ;
	void                    DoTrackLoop()                           ;
	void			FillEFHistograms()			;
	void                    DoClusterLoop()                         ; // /<May remove later	
	void                    DoCellLoop()                            ; // /<May remove later

        Bool_t                  IsMCprod                                ;// /<Flag for MC productions
	THistManager            fHistManager                            ;// /< Histogram manager
//	TList*			fOutput					;// /!<! Output list
 	private:
 	 AliAnalysisTaskEmcalJetEnergyFlow(const AliAnalysisTaskEmcalJetEnergyFlow&); // not implemented
  	 AliAnalysisTaskEmcalJetEnergyFlow &operator=(const AliAnalysisTaskEmcalJetEnergyFlow&); // not implemented

  	/// \cond CLASSIMP
  	  ClassDef(AliAnalysisTaskEmcalJetEnergyFlow,10);
	/// \endcond
};
#endif
