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
	
	private:
	void InitPlots();
	Int_t				fDebug;		// Debug value
	Float_t 			fConeRadius; 	// Cone radius to be used in filling
	Float_t				fNominalEnergy; // Force a nominal energy - specifically for 80+20 jets
	AliEMCALJetFinderOutput* 	fOutput;     	// Output object to be analysed
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
	Bool_t				fInitialised;
		
	
	ClassDef(AliEMCALJetFinderPlots,2)
	
};
#endif

