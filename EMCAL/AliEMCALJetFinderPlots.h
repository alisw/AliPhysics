#ifndef ALIEMCALJETFINDERPLOTS_H
#define ALIEMCALJETFINDERPLOTS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  *  * See cxx source for full Copyright notice     */


/* $Id$ */

//_________________________________________________________________________
//  Class for Filling jetfinder plots
//
//*-- Author: Mark Horner (LBL/UCT)
//
//



#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"

#include "AliEMCALJetFinderOutput.h"

class AliEMCALJetFinderPlots : public TObject
{
	public:	
	AliEMCALJetFinderPlots();
	~AliEMCALJetFinderPlots();
	void SetConeRadius(Float_t coneradius){fConeRadius = coneradius;}
	void SetNominalEnergy(Float_t energy){fNominalEnergy = energy;}
	void SetDebug(Int_t debug){fDebug = debug;}
	void FillFromOutput(AliEMCALJetFinderOutput* output);
	//========================== CASE 1 ========================
	// Only consider events with only 1 jet
	TH1F* GetFragmFcn(){return hFragmFcn;}	
	TH1F* GetPartonFragmFcn(){return hPartonFragmFcn;}	
	TH1F* GetJetJT(){return hJetJT;}	
	TH1F* GetPartonJT(){return hPartonJT;}	
	TH1F* GetJetPL(){return hJetPL;}
	TH1F* GetPartonPL(){return hPartonPL;}
	TH1F* GetJetEt(){return hJetEt;}
	TH1F* GetJetEta(){return hJetEta;}
	TH1F* GetPartonEta(){return hPartonEta;}
	TH1F* GetPartonPhi(){return hPartonPhi;}
	TH1F* GetJetPhi(){return hJetPhi;}
	TH1F* GetEtaDiff(){return hEtaDiff;}
	TH1F* GetPhiDiff(){return hPhiDiff;}
	TH2F* GetEtaPhiSpread(){return hEtaPhiSpread;}
	TH1F* GetNJets(){return hNJets;}

	//========================== CASE 2 ========================
	// Only consider events with at least 2 jets
	TH1F* GetFragmFcn2(){return hFragmFcn2;}	
	TH1F* GetPartonFragmFcn2(){return hPartonFragmFcn2;}	
	TH1F* GetJetJT2(){return hJetJT2;}	
	TH1F* GetPartonJT2(){return hPartonJT2;}	
	TH1F* GetJetPL2(){return hJetPL2;}
	TH1F* GetPartonPL2(){return hPartonPL2;}
	TH1F* GetJetEt2(){return hJetEt2;}
	TH1F* GetJetEta2(){return hJetEta2;}
	TH1F* GetPartonEta2(){return hPartonEta2;}
	TH1F* GetPartonPhi2(){return hPartonPhi2;}
	TH1F* GetJetPhi2(){return hJetPhi2;}
	TH1F* GetEtaDiff2(){return hEtaDiff2;}
	TH1F* GetPhiDiff2(){return hPhiDiff2;}
	TH2F* GetEtaPhiSpread2(){return hEtaPhiSpread2;}
	TH1F* GetNJets2(){return hNJets2;}
	TH1F* GetJetEtSecond2(){return hJetEtSecond2;}
	TH1F* GetJetEtRatio2(){return hJetEtRatio2;}
	TH1F* GetEtaPhiDist2(){return hEtaPhiDist2;}

	private:
	void InitPlots();
	Int_t				fDebug;		// Debug value
	Float_t 			fConeRadius; 	// Cone radius to be used in filling
	Float_t				fNominalEnergy; // Force a nominal energy - specifically for 80+20 jets
	AliEMCALJetFinderOutput* 	fOutput;     	// Output object to be analysed
	//===================== CASE 1 ===========================================
	TH1F				*hFragmFcn;	// ("hFragmFcn","Fragmentation Function",100,0,1);
	TH1F				*hPartonFragmFcn;// ("hFragmFcn","Parton Fragmentation Function",100,0,1);
	TH1F				*hPartonJT;	// ("hPartonJT","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
	TH1F				*hPartonPL;	// ("hPartonPL","Track Momentum Parallel to Parton Axis ",100,0.,100.);
	TH1F				*hJetJT;	// ("hJetJT","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
	TH1F				*hJetPL;	// ("hJetPL","Track Momentum Parallel to Jet Axis ",100,0.,100.);
	TH1F				*hJetEt;	// ("hJetEt","E_{T}^{reco}",250,0.,250.);
	TH1F				*hJetEta;	// ("hJetEta","#eta_{jet}^{reco}",180,-0.9,0.9);
	TH1F				*hJetPhi;	// ("hJetPhi","#phi_{jet}^{reco}",62,0.,3.1);
	TH1F				*hPartonEta;	// ("hPartonEta","#eta_{Parton}",180,-0.9,0.9);
	TH1F				*hPartonPhi;	// ("hPartonPhi","#phi_{Parton}",62,0.,3.1);
	TH1F				*hEtaDiff;	// ("hEtaDiff","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
	TH1F				*hPhiDiff;	// ("hPhiDiff","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
	TH2F				*hEtaPhiSpread;	// ("hEtaPhiSpread","#eta - #phi Distribution 
							//of Reconstructed Jets",192,-0.7,0.7,288,pi/3,pi);
	TH1F				*hNJets;	// ("hNJets","N Reconstructed jets",11,-0.5,10.5);
  
	//============================== CASE 2 ============================================

	TH1F				*hFragmFcn2;	// ("hFragmFcn2","Fragmentation Function",100,0,1);
	TH1F				*hPartonFragmFcn2;// ("hFragmFcn2","Parton Fragmentation Function",100,0,1);
	TH1F				*hPartonJT2;	// ("hPartonJT2","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
	TH1F				*hPartonPL2;	// ("hPartonPL2","Track Momentum Parallel to Parton Axis ",100,0.,100.);
	TH1F				*hJetJT2;	// ("hJetJT2","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
	TH1F				*hJetPL2;	// ("hJetPL2","Track Momentum Parallel to Jet Axis ",100,0.,100.);
	TH1F				*hJetEt2;	// ("hJetEt2","E_{T}^{reco}",250,0.,250.);
	TH1F				*hJetEta2;	// ("hJetEta2","#eta_{jet}^{reco}",180,-0.9,0.9);
	TH1F				*hJetPhi2;	// ("hJetPhi2","#phi_{jet}^{reco}",62,0.,3.1);
	TH1F				*hPartonEta2;	// ("hPartonEta2","#eta_{Parton}",180,-0.9,0.9);
	TH1F				*hPartonPhi2;	// ("hPartonPhi2","#phi_{Parton}",62,0.,3.1);
	TH1F				*hEtaDiff2;	// ("hEtaDiff2","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
	TH1F				*hPhiDiff2;	// ("hPhiDiff2","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
	TH2F				*hEtaPhiSpread2;	// ("hEtaPhiSpread2","#eta - #phi Distribution 
							//of Reconstructed Jets",192,-0.7,0.7,288,pi/3,pi);
	TH1F				*hNJets2;	// ("hNJets2","N Reconstructed jets",11,-0.5,10.5);
	TH1F				*hJetEtSecond2; //("hJetEtSecond2","E_{T}^{reco}",250,0.,250.); 
	TH1F				*hJetEtRatio2;  //("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);
	TH1F 				*hEtaPhiDist2;  //("hEtaPhiDist2","Angular Distance Between First and Second",100,0,3);



	
	Bool_t				fInitialised;
		
	
	ClassDef(AliEMCALJetFinderPlots,3)
	
};
#endif

