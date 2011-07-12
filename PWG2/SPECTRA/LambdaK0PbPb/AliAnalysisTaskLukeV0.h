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
    TH1F            *fHistPt;				// Pt spectrum
    TH1F            *fHistEta;				// pseudorapidity spectrum
	TH1F			*fHistLuke;				// test added histogram
	TH1F			*fHistBetaV0;			// beta of V0s
	TH1F			*fHistCosPA;			// Cosine of pointing angle
	TH1F			*fHistCosTheta;			// Cosine of decay angle in 'lambda' rest frame
	TH1F			*fHistCosThetaPi2;		//	cosine of 'pion' decay angle in lab frame
	TH1F			*fHistCosThetaProton2;	//	cosine of 'proton' decay angle in lab frame
	TH1F			*fHistDCAV0Daughters;	//	DCA between the v0 daughters
	TH1F			*fHistDecayL;			//	3d decay length of V0
	TH1F			*fHistDecayLxy;			//	2d decay length of V0 in xy plane
	TH1F			*fHistDeltaTheta;		//	Difference between theta of 'pion' & 'proton' in 'lambda' rest frame
	TH1F			*fHistdNV0sdT2;			//	Ratio of number of v0s to number of tracks squared
	TH1F			*fHistEtaNTracks;		//	Eta of negative tracks
	TH1F			*fHistEtaPTracks;		//	Eta of positive tracks
	TH1F			*fHistEtaV0s;			//	Eta of v0s
	TH1F			*fHistImpactxyN;		//	Impact paramater of negative daughter in xy plane
	TH1F			*fHistImpactzN;			//  Impact paramater of negative daughter in z direction
	TH1F			*fHistImpactxyP;		//  Impact paramater of positive daughter in xy plane
	TH1F			*fHistImpactzP;			//  Impact paramater of positive daughter in z direction
	TH1F			*fHistKinkIndexFalse;	//	Lambda mass for non-kink candidates
	TH1F			*fHistKinkIndexTrue;	//	Lambda mass for kink candidates
	TH1F			*fHistLambdaBgRapidity;	//	Rapididity of V0s in sidebands of lambda mass
	TH1F			*fHistLambdaBgEta;		//	Eta of V0s in sidebands of lambda mass
	TH1F			*fHistLambdaBgPt;		//	Transverse momentum of V0s in sidebands of lambda mass
	TH1F			*fHistLambdaSigRapidity;//	Rapididity of V0s in peak of lambda mass
	TH1F			*fHistLambdaSigEta;		//	Eta of V0s in peak of lambda mass
	TH1F			*fHistLambdaSigPt;		//	Transverse momentum of V0s in peak of lambda mass
	TH1F			*fHistMagneticField;	//	Magnetic field of events
	TH1F			*fHistMcLog;			//	Log of monte-carlo related information. Code in histogram.
	TH1F			*fHistMcNLambdaPrimary;	//	Number of primary lambdas in MC event
	TH1F			*fHistMcNLambda;		//	Number of lambdas in MC event
	TH1F			*fHistMcNAntilambda;	//	Number of antilambdas in MC event
	TH1F			*fHistMcNKshort;		//	Number of K0s in MC event
	TH1F			*fHistMcPtV0La;			//	Transverse Momentum distribution of lambdas in MC data
	TH1F			*fHistMcPtV0Lb;			//	Transverse Momentum distribution of antilambdas in MC data
	TH1F			*fHistMcPtV0K0;			//	Transverse Momentum distribution of K0s in MC data
	TH1F			*fHistMcSigmaPProton;	//	TPC response sigma of being a proton for MC protons
	TH1F			*fHistMcSigmaPNProton;	//	TPC response sigma of being a proton for MC non-protons
	TH1F			*fHistMcSLambdaEta;		//	Eta distribution of secondary lambdas in MC
	TH1F			*fHistMcSLambdaDaughter;//	Daughter PDG code of Secondary lambdas in MC
	TH1F			*fHistMcSLambdaPt;		//	Transverse momentum distribution of secondary lambdas in MC
	TH1F			*fHistMcTPCTrackLength;	//	TPC track length in MC - DOESNT WORK
	TH1F			*fHistMK0;				//	Reconstructed K0 mass for all V0s
	TH1F			*fHistMK00;				//	Reconstructed K0 mass satisfying chosen conditions, typically momentum range 
	TH1F			*fHistMK01; 			//	"
	TH1F			*fHistMK02; 			//	"
	TH1F			*fHistMK03; 			//	"
	TH1F			*fHistMK04; 			//	"
	TH1F			*fHistMK05; 			//	"
	TH1F			*fHistMK06; 			//	"
	TH1F			*fHistMK07; 			//	"
	TH1F			*fHistMK08; 			//	"
	TH1F			*fHistMK09; 			//	"
	TH1F			*fHistMLa;				//	Reconstructed Lambda mass for all V0s
	TH1F			*fHistMLa0; 			//	Reconstructed Lambda mass satisfying chosen conditions, typically momentum range
	TH1F			*fHistMLa1; 			//	"
	TH1F			*fHistMLa2; 			//	"
	TH1F			*fHistMLa3; 			//	"
	TH1F			*fHistMLa4; 			//	"
	TH1F			*fHistMLa5; 			//	"
	TH1F			*fHistMLa6; 			//	"
	TH1F			*fHistMLa7; 			//	"
	TH1F			*fHistMLa8; 			//	"
	TH1F			*fHistMLa9; 			//	"
	TH1F			*fHistMLb;				//	Reconstructed Antilambda mass for all V0s
	TH1F			*fHistMLb0; 			//	Reconstructed Antilambda mass satisfying chosen conditions, typically momentum range
	TH1F			*fHistMLb1; 			//	"
	TH1F			*fHistMLb2; 			//	"
	TH1F			*fHistMLb3; 			//	"
	TH1F			*fHistMLb4; 			//	"
	TH1F			*fHistMLb5; 			//	"
	TH1F			*fHistMLb6; 			//	"
	TH1F			*fHistMLb7; 			//	"
	TH1F			*fHistMLb8; 			//	"
	TH1F			*fHistMLb9; 			//	"
	TH1F			*fHistNLambda;			//	Number of lambda candidates in peak region per event
	TH1F			*fHistNV0;				//	Number of V0s per event
	TH1F			*fHistNTracks;			//	Number of tracks per event
	TH1F			*fHistPtV0;				//	Transverse momentum distribution of V0s
	TH1F			*fHistPVZ;				//	Distribution of Z coordinate of primary vertex
	TH1F			*fHistTauLa;			//	Distribution of 'lambda' lifetime*c
	TH1F			*fHistTheta;			//	Angle of decay in 'lambda' rest frame
	TH1F			*fHistThetaPi2;			//	Angle of decay of 'pion' in rest frame
	TH1F			*fHistThetaProton2;		//	Angle of decay of 'proton' in rest frame
	TH1F			*fHistV0Z;				//	Z coordinate of V0 decay
	TH1F			*fHistZ;				//	Distance in Z between primary vertex and v0 decay
	TH2F			*fHistBetaERatio;		//	Beta of V0 vs ratio of energies between 'pion' & 'proton'
	TH2F			*fHistBetaPRatio;		//	Beta of V0 vs ratio of momenta between 'pion' & 'proton'
	TH2F			*fHistBetaPXRatio;		//	Beta of V0 vs ratio of x-momenta between 'pion' & 'proton'
	TH2F			*fHistBetheBlochITSNeg; //	ITS response for negative daughters vs momentum
	TH2F			*fHistBetheBlochITSPos; //	ITS response for positive daughters vs momentum
	TH2F			*fHistBetheBlochTPCNeg; //	TPC response for negative daughters vs momentum
	TH2F			*fHistBetheBlochTPCPos; //	TPC response for positive daughters vs momentum
	TH2F			*fHistCosPAMLa;			//	Cosine of pointing angle vs reconstructed lambda mass
	TH2F			*fHistDCAPtPSig;		//	DCA vs momenum for 'proton' from 'lambda' for V0s in lambda mass peak
	TH2F			*fHistDCAPtPBg;			//	DCA vs momenum for 'proton' from 'lambda' for V0s in lambda mass sidebands
	TH2F			*fHistDCAPtPbarSig;		//	DCA vs momenum for 'antiproton' from 'antilambda' for V0s in antilambda mass peak
	TH2F			*fHistDCAPtPbarBg;		//	DCA vs momenum for 'antiproton' from 'antilambda' for V0s in antilambda mass sidebands
	TH2F			*fHistDecayLDCA;		//	Decay length of V0 vs DCA
	TH2F			*fHistDecayLMLa;		//	Decay length of v0 vs reconstructed 'lambda' mass
	TH2F			*fHistDeltaThetaMLa;	//	Difference between theta of 'pion' & 'proton' in 'lambda' rest frame vs reconstructed 'lambda' mass
	TH2F			*fHistImpactxyImpactz;	//	Impact paramater in xy vs z
	TH2F			*fHistImpactxyMLa;		//	Impact paramater in xy vs reconstructed lambda mass
	TH2F			*fHistImpactzMLa;		//  Impact paramater in z vs reconstructed lambda mass
	TH2F			*fHistMcLamdaPProductionVertex;		//	Production vertex of primary lambdas in MC 
	TH2F			*fHistMcLamdaSProductionVertex;		//	Production vertex of secondary lambdas in MC 
	TH2F			*fHistMcLamdaSDecayVertex;			//	Decay vertex of secondary lambdas in MC 
	TH2F			*fHistMcPMK0Pt;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC
	TH2F			*fHistMcPMLaPt;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC
	TH2F			*fHistMcPMLbPt;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC
	TH2F			*fHistMcSLambdaDaughterPairs;		//	PDG codes of daughters of Secondary lambdas in MC
	TH2F			*fHistMcV0MK0Pt;		//	Reconstructed K0 mass vs transverse momentum of V0s passing cuts in MC
	TH2F			*fHistMcV0MLaPt;		//	Reconstructed Lambda mass vs transverse momentum of V0s passing cuts in MC
	TH2F			*fHistMcV0MLbPt;		//	Reconstructed Antilambda mass vs transverse momentum of V0s passing cuts in MC
	TH2F			*fHistMcV0LamdaSProductionVertex;	//	Production vertex of V0 secondary Lambdas in RZ plane in MC
	TH2F			*fHistMcV0IDLamdaSProductionVertex;	//  Production vertex of identified V0 secondary Lambdas in RZ plane in MC
	TH2F			*fHistMK0PtArm;			//	Mass of 'K0' vs Armenteros transverse momentum
	TH2F			*fHistMK0Pt;			//	Mass of 'K0' vs transverse momentum
	TH2F			*fHistMK0R;				//	Mass of 'K0' vs radius
	TH2F			*fHistMLaPtArm;			//	Mass of 'Lambda' vs Armenteros transverse momentum
	TH2F			*fHistMLaPt;			//	Mass of 'Lambda' vs transverse momentum
	TH2F			*fHistMLaR;				//	Mass of 'Lambda' vs radius
	TH2F			*fHistMLbPtArm;			//	Mass of 'Antilambda' vs Armenteros transverse momentum
	TH2F			*fHistMLbPt;			//	Mass of 'Antilambda' vs transverse momentum
	TH2F			*fHistMLbR;				//	Mass of 'Antilambda' vs radius  
	TH2F			*fHistNegChi2PerClusterMLa;			//	Chi2 per cluster of negative daughter vs reconstructed 'lambda' mass
	TH2F			*fHistNegTPCRefitMLa;	//	TPC refit t/f of negative daughter vs reconstructed 'lambda' mass
	TH2F			*fHistNTPCNegClustersMLa;			//	Number of TPC clusters  of negative daughter vs reconstructed 'lambda' mass
	TH2F			*fHistNTPCPosClustersMLa;			// Number of TPC clusters of positive daughter vs reconstructed 'lambda' mass
	TH2F			*fHistNV0sNTracks;		//	Number of V0s vs number of tracks per event
	TH2F			*fHistPosChi2PerClusterMLa;			//	Chi2 per cluster of positive daughter vs reconstructed 'lambda' mass
	TH2F			*fHistPosTPCRefitMLa;	//	TPC refit t/f of positive daughter vs reconstructed 'lambda' mass
	TH2F			*fHistPtArm;			//	Podolanski-Armenteros plot of V0s
	TH2F			*fHistPtArmR;			//	PtArm vs radius of v0
	TH2F			*fHistPtV0Z;			//	Transverse momentum vs Z coordinate of v0 decay
	TH2F			*fHistRZ;				//	Radius vs z of V0 decay
	TH2F			*fHistTauLaMLa;			//	'lambda' lifetime vs reconstructed 'lambda' mass
	TH2F			*fHistXZ;				//	X vs z of V0 decay
	TH2F			*fHistYZ;				//	Y vs Z of V0 decay
	
    // NEW HISTO to be declared here
    
    AliAnalysisTaskLukeV0(const AliAnalysisTaskLukeV0&); // not implemented
    AliAnalysisTaskLukeV0& operator=(const AliAnalysisTaskLukeV0&); // not implemented
    
    ClassDef(AliAnalysisTaskLukeV0, 1); // example of analysis
};

#endif

