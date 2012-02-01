/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskLukeV0.h 45956 2010-12-10 12:55:37Z agheata $ */
/* AliAnalysisTaskLukeV0.h
 *
 * Header file for task AliAnalysisTaskLukeV0.cxx
 * The task is designed to extract Lambda, Antilambda & K0s spectra from PbPb collision data/MC 
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 * Further adapted by Luke Hanratty
 */

#ifndef ALIANALYSISTASKLukeV0_H
#define ALIANALYSISTASKLukeV0_H

class TH1F;
class TH2F;
class TList;
class AliESDtrackCuts;
class AliPIDResponse;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskLukeV0 : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskLukeV0();
    AliAnalysisTaskLukeV0(const char *name);
    virtual ~AliAnalysisTaskLukeV0();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    
 private:
    TList           *fOutputList;			// Output list
    AliESDtrackCuts *fTrackCuts;			// Track cuts
	AliPIDResponse	*fPIDResponse;			// PID
    TH1F			*fHistCosPA;			// Cosine of pointing angle
	TH1F			*fHistDCAV0Daughters;	//	DCA between the v0 daughters
	TH1F			*fHistDecayL;			//	3d decay length of V0
	TH1F			*fHistImpactxyN;		//	Impact paramater of negative daughter in xy plane
	TH1F			*fHistImpactzN;			//  Impact paramater of negative daughter in z direction
	TH1F			*fHistImpactxyP;		//  Impact paramater of positive daughter in xy plane
	TH1F			*fHistImpactzP;			//  Impact paramater of positive daughter in z direction
	TH1F			*fHistMCLambdacTau;		// cTau distribution of MC lambdas
	TH1F			*fHistMCLambdaNotcTau;	// cTau distribution of MC lambda background
	TH1F			*fHistMCLambdaDecayL;	// Decay Length distribution of MC lambdas
	TH1F			*fHistMCLambdaNotDecayL;	// Decay Length distribution of MC lambda background
	TH1F			*fHistMcNLambdaPrimary;	//	Number of primary lambdas in MC event
	TH1F			*fHistMcNLambda;		//	Number of lambdas in MC event
	TH1F			*fHistMcNAntilambda;	//	Number of antilambdas in MC event
	TH1F			*fHistMcNKshort;		//	Number of K0s in MC event
	TH1F			*fHistMK0;				//	Reconstructed K0 mass for all V0s
	TH1F			*fHistMLa;				//	Reconstructed Lambda mass for all V0s
	TH1F			*fHistMLb;				//	Reconstructed Antilambda mass for all V0s
	TH1F			*fHistNLambda;			//	Number of lambda candidates in peak region per event
	TH1F			*fHistNV0;				//	Number of V0s per event
	TH1F			*fHistPtV0;				//	Transverse momentum distribution of V0s
	TH1F			*fHistPVZ;				//	Distribution of Z coordinate of primary vertex
	TH1F			*fHistTauLa;			//	Distribution of 'lambda' lifetime*c
	TH1F			*fHistV0Z;				//	Z coordinate of V0 decay
	TH1F			*fHistZ;				//	Distance in Z between primary vertex and v0 decay
	TH2F			*fHistBetheBlochTPCNeg; //	TPC response for negative daughters vs momentum
	TH2F			*fHistBetheBlochTPCPos; //	TPC response for positive daughters vs momentum
	TH2F			*fHistImpactxyImpactz;	//	Impact paramater in xy vs z
	TH2F			*fHistMcPMK0Pt;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC
	TH2F			*fHistMcPMLaPt;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC
	TH2F			*fHistMcPMLbPt;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC
	TH2F			*fHistMcV0MK0Pt;		//	Reconstructed K0 mass vs transverse momentum of V0s passing cuts in MC
	TH2F			*fHistMcV0MLaPt;		//	Reconstructed Lambda mass vs transverse momentum of V0s passing cuts in MC
	TH2F			*fHistMcV0MLbPt;		//	Reconstructed Antilambda mass vs transverse momentum of V0s passing cuts in MC
	TH2F			*fHistMK0Pt;			//	Mass of 'K0' vs transverse momentum
	TH2F			*fHistMLaPt;			//	Mass of 'Lambda' vs transverse momentum
	TH2F			*fHistMLbPt;			//	Mass of 'Antilambda' vs transverse momentum
	TH2F			*fHistPtArm;			//	Podolanski-Armenteros plot of V0s
	TH2F			*fHistPtV0Z;			//	Transverse momentum vs Z coordinate of v0 decay
	TH2F			*fHistRZ;				//	Radius vs z of V0 decay
	TH2F			*fHistXZ;				//	X vs z of V0 decay
	TH2F			*fHistYZ;				//	Y vs Z of V0 decay
	
    // NEW HISTO to be declared here
    
    AliAnalysisTaskLukeV0(const AliAnalysisTaskLukeV0&); // not implemented
    AliAnalysisTaskLukeV0& operator=(const AliAnalysisTaskLukeV0&); // not implemented
    
    ClassDef(AliAnalysisTaskLukeV0, 1); // example of analysis
};

#endif

