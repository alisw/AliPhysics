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

	AliEMCALJetFinderPlots (const AliEMCALJetFinderPlots&);
	AliEMCALJetFinderPlots & operator = (AliEMCALJetFinderPlots &) {
	  Fatal("operator =", "not implemented") ;
	  return *this ;
	}

	void SetConeRadius(Float_t coneradius){fConeRadius = coneradius;}
	void SetNominalEnergy(Float_t energy){fNominalEnergy = energy;}
	void SetDebug(Int_t debug){fDebug = debug;}
	void SetBackHisto(TH1F* histo){fhBackHisto=histo;}
	void FillFromOutput(AliEMCALJetFinderOutput* output,Float_t weight=1.0);
	//========================== CASE 1 ========================
	// Only consider events with only 1 jet
	TH1F* GetFragmFcn() const {return fhFragmFcn;}	
	TH1F* GetPartonFragmFcn() const {return fhPartonFragmFcn;}	
	TH1F* GetPT() const {return fhJetPT;}	
	TH1F* GetPartonPT() const {return fhPartonPT;}	
	TH1F* GetJetJT() const {return fhJetJT;}	
	TH1F* GetPartonJT() const {return fhPartonJT;}	
	TH1F* GetJetPL() const {return fhJetPL;}
	TH1F* GetPartonPL() const {return fhPartonPL;}
	TH1F* GetJetEt() const {return fhJetEt;}
	TH1F* GetJetEtDiff() const {return fhJetEtDiff;}
	TH1F* GetJetEta() const {return fhJetEta;}
	TH1F* GetPartonEta() const {return fhPartonEta;}
	TH1F* GetPartonPhi() const {return fhPartonPhi;}
	TH1F* GetJetPhi() const {return fhJetPhi;}
	TH1F* GetEtaDiff() const {return fhEtaDiff;}
	TH1F* GetPhiDiff()const {return fhPhiDiff;}
	TH2F* GetEtaPhiSpread() const {return fhEtaPhiSpread;}
	TH1F* GetNJets() const {return fhNJets;}

	//========================== CASE 2 ========================
	// Only consider events with at least 2 jets
	TH1F* GetFragmFcn2() const {return fhFragmFcn2;}	
	TH1F* GetPartonFragmFcn2() const {return fhPartonFragmFcn2;}	
	TH1F* GetPT2() const {return fhJetPT2;}	
	TH1F* GetPartonPT2() const {return fhPartonPT2;}	
	TH1F* GetJetJT2() const {return fhJetJT2;}	
	TH1F* GetPartonJT2() const {return fhPartonJT2;}	
	TH1F* GetJetPL2() const {return fhJetPL2;}
	TH1F* GetPartonPL2() const {return fhPartonPL2;}
	TH1F* GetJetEt2() const {return fhJetEt2;}
	TH1F* GetJetEtDiff2() const {return fhJetEtDiff2;}
	TH1F* GetJetEta2() const {return fhJetEta2;}
	TH1F* GetPartonEta2() const {return fhPartonEta2;}
	TH1F* GetPartonPhi2() const {return fhPartonPhi2;}
	TH1F* GetJetPhi2() const {return fhJetPhi2;}
	TH1F* GetEtaDiff2() const {return fhEtaDiff2;}
	TH1F* GetPhiDiff2() const {return fhPhiDiff2;}
	TH2F* GetEtaPhiSpread2() const {return fhEtaPhiSpread2;}
	TH1F* GetNJets2() const {return fhNJets2;}
	TH1F* GetJetEtSecond2() const {return fhJetEtSecond2;}
	TH1F* GetJetEtRatio2() const {return fhJetEtRatio2;}
	TH1F* GetEtaPhiDist2() const {return fhEtaPhiDist2;}

	
	TH1F* GetJetPt() const {return	fhJetPT ;}
	TH1F* GetPartonPt() const {return fhPartonPT ;}
	TH1F* GetJetPt2() const {return 	fhJetPT2;}
	TH1F* GetPartonPt2() const {return fhPartonPT2;}
	TH1F* GetRecoBinFragmFcn() const {return fhRecoBinFragmFcn;}
	TH1F* GetRecoBinFragmFcnNoBg() const {return fhRecoBinFragmFcnNoBg;}
	TH1F* GetRecoBinPartonFragmFcn() const {return	fhRecoBinPartonFragmFcn;}
	
	//============================== ALL CASES ============================================
	
	TH2F* GetInputOutput() const {return fhInputOutput;}
	
	//============================== Reconstruction Bin Comparison  ============================================
	
	TH1F* GetRecoBinPt() const {return fhRecoBinPt;}	           // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
	TH1F* GetRecoBinPtNoBg() const {return fhRecoBinPtNoBg;}	           // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
	TH1F* GetRecoBinPartonPt() const {return fhRecoBinPartonPt;}      // ("fhRecoBinPartonPt","Input Pt Distribution",100,0,1);
	TH1F* GetRecoBinJetEt() const {return fhRecoBinJetEt;}            // ("fhRecoJetEt","E_{T}^{reco}",250,0.,250.);
	TH1F* GetRecoBinInputJetEt() const {return fhRecoBinInputJetEt;}  // ("fhRecoInputJetEt","E_{T}^{reco}",250,0.,250.);

	private:
	void InitPlots();
	Int_t				fDebug;		// Debug value
	Float_t 			fConeRadius; 	// Cone radius to be used in filling
	Float_t				fNominalEnergy; // Force a nominal energy - specifically for 80+20 jets
	AliEMCALJetFinderOutput* 	fOutput;     	// Output object to be analysed
	//===================== CASE 1 ===========================================
	TH1F				*fhFragmFcn;	// ("hFragmFcn","Fragmentation Function",100,0,1);
	TH1F				*fhPartonFragmFcn;// ("hFragmFcn","Parton Fragmentation Function",100,0,1);
	TH1F				*fhPartonJT;	// ("hPartonJT","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
	TH1F				*fhPartonPL;	// ("hPartonPL","Track Momentum Parallel to Parton Axis ",100,0.,100.);
	TH1F				*fhJetJT;	// ("hJetJT","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
	TH1F				*fhJetPL;	// ("hJetPL","Track Momentum Parallel to Jet Axis ",100,0.,100.);
	TH1F				*fhJetEt;	// ("hJetEt","E_{T}^{reco}",250,0.,250.);
	TH1F				*fhJetEtDiff;	// ("hJetEt","E_{T}^{reco}",250,0.,250.);
	TH1F				*fhJetEta;	// ("hJetEta","#eta_{jet}^{reco}",180,-0.9,0.9);
	TH1F				*fhJetPhi;	// ("hJetPhi","#phi_{jet}^{reco}",62,0.,3.1);
	TH1F				*fhPartonEta;	// ("hPartonEta","#eta_{Parton}",180,-0.9,0.9);
	TH1F				*fhPartonPhi;	// ("hPartonPhi","#phi_{Parton}",62,0.,3.1);
	TH1F				*fhEtaDiff;	// ("hEtaDiff","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
	TH1F				*fhPhiDiff;	// ("hPhiDiff","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
	TH2F				*fhEtaPhiSpread;	// ("hEtaPhiSpread","#eta - #phi Distribution 
							//of Reconstructed Jets",192,-0.7,0.7,288,pi/3,pi);
	TH1F				*fhNJets;	// ("hNJets","N Reconstructed jets",11,-0.5,10.5);
  
	//============================== CASE 2 ============================================

	TH1F				*fhFragmFcn2;	// ("hFragmFcn2","Fragmentation Function",100,0,1);
	TH1F				*fhPartonFragmFcn2;// ("hFragmFcn2","Parton Fragmentation Function",100,0,1);
	TH1F				*fhPartonJT2;	// ("hPartonJT2","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
	TH1F				*fhPartonPL2;	// ("hPartonPL2","Track Momentum Parallel to Parton Axis ",100,0.,100.);
	TH1F				*fhJetJT2;	// ("hJetJT2","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
	TH1F				*fhJetPL2;	// ("hJetPL2","Track Momentum Parallel to Jet Axis ",100,0.,100.);
	TH1F				*fhJetEt2;	// ("hJetEt2","E_{T}^{reco}",250,0.,250.);
	TH1F				*fhJetEtDiff2;	// ("hJetEt","E_{T}^{reco}",250,0.,250.);
	TH1F				*fhJetEta2;	// ("hJetEta2","#eta_{jet}^{reco}",180,-0.9,0.9);
	TH1F				*fhJetPhi2;	// ("hJetPhi2","#phi_{jet}^{reco}",62,0.,3.1);
	TH1F				*fhPartonEta2;	// ("hPartonEta2","#eta_{Parton}",180,-0.9,0.9);
	TH1F				*fhPartonPhi2;	// ("hPartonPhi2","#phi_{Parton}",62,0.,3.1);
	TH1F				*fhEtaDiff2;	// ("hEtaDiff2","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
	TH1F				*fhPhiDiff2;	// ("hPhiDiff2","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
	TH2F				*fhEtaPhiSpread2;	// ("hEtaPhiSpread2","#eta - #phi Distribution 
							//of Reconstructed Jets",192,-0.7,0.7,288,pi/3,pi);
	TH1F				*fhNJets2;	// ("hNJets2","N Reconstructed jets",11,-0.5,10.5);
	TH1F				*fhJetEtSecond2; //("hJetEtSecond2","E_{T}^{reco}",250,0.,250.); 
	TH1F				*fhJetEtRatio2;  //("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);
	TH1F 				*fhEtaPhiDist2;  //("hEtaPhiDist2","Angular Distance Between First and Second",100,0,3);

	//============================== ALL CASES ============================================

	TH2F				*fhInputOutput;  //("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);

	//============================== Reconstruction Bin Comparison  ============================================
	
	TH1F				*fhRecoBinPt;	       // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
	TH1F				*fhRecoBinPtNoBg;	       // ("fhRecoBinPt","Reconstructed Pt Distribution",100,0,1);
	TH1F				*fhRecoBinPartonPt;    // ("fhRecoBinPartonPt","Input Pt Distribution",100,0,1);
	TH1F				*fhRecoBinJetEt;       // ("fhRecoJetEt","E_{T}^{reco}",250,0.,250.);
	TH1F				*fhRecoBinInputJetEt;  // ("fhRecoInputJetEt","E_{T}^{reco}",250,0.,250.);
    TH1F* 				fhJetPT ;// new TH1F("hJetPT","P_{T} Distribution",200,0,200);
    TH1F* 				fhPartonPT ;// new TH1F("hPartonPT","Parton P_{T} Distribution",200,0,1);
    TH1F* 				fhJetPT2 ;// new TH1F("hJetPT","P_{T} Distribution",200,0,200);
    TH1F* 				fhPartonPT2 ;// new TH1F("hPartonPT","Parton P_{T} Distribution",200,0,1);
    TH1F* 				fhRecoBinFragmFcn;//new TH1F("fhRecoBinFragmFcn","Reconstructed Frag. Fcn",100,0,1);
    TH1F* 				fhRecoBinFragmFcnNoBg;//new TH1F("fhRecoBinFragmFcn","Reconstructed Frag. Fcn",100,0,1);
    TH1F* 				fhRecoBinPartonFragmFcn;// new TH1F("fhRecoBinPartonFragmFcn","Input Bin Fragm Fcn Distribution",100,0,1);
        
	TH1F* 				fhJetInvE;// new TH1F("fhJetInvE","#frac{1}{E_{R}}",100,0,1);
	TH1F* 				fhJetInvE2;// new TH1F("fhJetInvE2","#frac{1}{E_{R}}",100,0,1);
 
	TH1F*				fhBackHisto;
	Float_t				fScaleFactor; //Scaling to get back to correct energy
	Bool_t				fInitialised; // have histograms been initialised
		
	
	ClassDef(AliEMCALJetFinderPlots,6)
	
};
#endif

