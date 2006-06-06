/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnAnalysis
//             Reconstruction and analysis of K* Rsn
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef AliRsnANALYSIS_H
#define AliRsnANALYSIS_H

#include <TNamed.h>

#include "AliPID.h"

class TH1D;
class TObjArray;
class AliRsnEvent;
class AliRsnDaughter;
class AliRsnDaughterCut;

class AliRsnAnalysis : public TObject
{
public:

	         AliRsnAnalysis();
	virtual ~AliRsnAnalysis() {Clear();}
	
	void     AddCutPair(AliRsnDaughterCut *cut);
	void     AddCutSingle(AliPID::EParticleType type, AliRsnDaughterCut *cut);
	void     AddMixPairDef(AliPID::EParticleType p1, Char_t s1, AliPID::EParticleType p2, Char_t s2);
	void     AddPairDef(AliPID::EParticleType p1, Char_t s1, AliPID::EParticleType p2, Char_t s2, Bool_t onlyTrue = kFALSE);
	void     Clear(Option_t *option = "");	
	Stat_t   EventMix(Int_t nmix = 5, Int_t multDiffMax = 10, Double_t vzDiffMax = 0.01, Bool_t compareTotalMult = kFALSE);
	Stat_t   Process();
	void     SetBins(Int_t nbins, Double_t min, Double_t max) {fNBins=nbins;fHistoMin=min;fHistoMax=max;}
	void     SetEventsTree(TTree *tree)                       {fEventsTree = tree;}
	void     SetRejectFakes(Bool_t doit=kTRUE)                {fRejectFakes = doit;}		
	void     SetTrueMotherPDG(Int_t pdg)                      {fTrueMotherPDG = pdg;}
	void     WriteHistograms();
	
private:

	class AliPairDef : public TNamed
	{
	public:
	
		AliPairDef(AliPID::EParticleType p1, Char_t s1, AliPID::EParticleType p2, Char_t s2, Int_t pdg, Bool_t onlyTrue = kFALSE);
				   
		virtual ~AliPairDef()                 { }
		
		Char_t                 GetSign1()     {return fSign1;}
		Char_t                 GetSign2()     {return fSign2;}
		Bool_t                 GetOnlyTrue()  {return fOnlyTrue;}
		AliPID::EParticleType  GetParticle1() {return fParticle1;}
		AliPID::EParticleType  GetParticle2() {return fParticle2;}
		Double_t               GetMass1()     {return fMass1;}
		Double_t               GetMass2()     {return fMass2;}
				
	    void SetSign1(Char_t value)                 {fSign1 = value;}
		void SetSign2(Char_t value)                 {fSign2 = value;}
		void SetOnlyTrue(Bool_t value = kTRUE)      {fOnlyTrue = value;}
		void SetTrueMotherPDG(Int_t pdg)            {fTrueMotherPDG = pdg;}
		void SetParticle1(AliPID::EParticleType p)  {fParticle1 = p;}
		void SetParticle2(AliPID::EParticleType p)  {fParticle2 = p;}
		
		Text_t* ParticleName(AliPID::EParticleType part);
		
	private:
	
		Bool_t                 fOnlyTrue;	   // flag to be used for spectra of true pairs
		Int_t                  fTrueMotherPDG; // PDG code of true mother (if requested)
		
		Double_t               fMass1;     //
		Char_t                 fSign1;     // info about particle 1
		AliPID::EParticleType  fParticle1; //
	
		Double_t               fMass2;     //
		Char_t                 fSign2;     // info about particle 2
		AliPID::EParticleType  fParticle2; //
	};
	
	Stat_t     Compute(AliPairDef *pd, TH1D* &h, AliRsnEvent *ev1, AliRsnEvent *ev2);
	Bool_t     SingleCutCheck(Int_t itype, AliRsnDaughter *track);
	Bool_t     PairCutCheck(AliRsnDaughter *track1, AliRsnDaughter *track2);
	
	Bool_t     fRejectFakes;             // reject particles labeled as fake
	
	Int_t      fNBins;                   // number of histogram bins
	Double_t   fHistoMin;                // minimum of the histograms
	Double_t   fHistoMax;                // maximum of the histograms
	
	Int_t      fTrueMotherPDG;           // PDG code of true mother (used to create 'true' histos)
	
	TObjArray *fMixPairDefs;             //  list of pair definitions for histograms (event mixing)
	TObjArray *fMixHistograms;           //! list of invmass histograms created (event mixing)
	
	TObjArray *fPairDefs;                //  list of pair definitions for histograms
	TObjArray *fHistograms;              //! list of invmass histograms created

	TObjArray *fCuts[AliPID::kSPECIES];  //! list of single particle cuts for each particle type
	TObjArray *fPairCuts;                //! list of pair cuts
	TTree     *fEventsTree;              //! TTree of events (can not be created here, must be passed)
	
	// Rsn analysis implementation
	ClassDef(AliRsnAnalysis,1)
};

#endif
