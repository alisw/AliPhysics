/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnAnalysis
//             Reconstruction and analysis of K* Rsn
// ........................................
// ........................................
// ........................................
// ........................................
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNPAIR_H
#define ALIRSNPAIR_H

#include <TNamed.h>
#include "AliRsnPID.h"

class TH1D;
class TRefArray;
class TObjArray;
class AliRsnEvent;
class AliRsnDaughter;
class AliRsnDaughterCut;
class AliRsnDaughterCutPair;

class AliRsnPair : public TNamed
{

public:
	
	AliRsnPair();
	AliRsnPair(const char *name, const char *title, 
	           Int_t nbins, Double_t min, Double_t max,
	           Double_t ptmin = 0., Double_t ptmax = 0., Double_t dmax = 0.);
	AliRsnPair(const AliRsnPair &copy);
	const AliRsnPair& operator=(const AliRsnPair &copy);
	virtual ~AliRsnPair() {Clear();}
	virtual void Clear(Option_t *option = "");
	
	/* getters */
	TH1D*             GetHistogram() {return fHistogram;}
	Char_t            GetCharge(Int_t i) const {if (i>=0&&i<2) return fCharge[i]; else return 0;}
	AliRsnPID::EType  GetParticle(Int_t i) const {if (i>=0&&i<2) return fType[i]; else return AliRsnPID::kUnknown;}
	Double_t          GetMass(Int_t i) const {if (i>=0&&i<2) return fMass[i]; else return 0.0;}
	Bool_t            StoreOnlyTruePairs() const {return fStoreOnlyTrue;}
	Bool_t            IsForMixing() const {return fForMixing;}
	
	/* setters */
	void SetMass(Int_t i, Double_t value) {if (i>=0&&i<2) fMass[i] = value;}
	void SetTrueMotherPDG(Int_t pdg) {fTrueMotherPDG = pdg;}
	void SetPair(Char_t charge1, AliRsnPID::EType pid1, Char_t charge2, AliRsnPID::EType pid2);
    void SetPtBin(Double_t min, Double_t max) {fPtMin = min; fPtMax = max;}
    void SetImpactMax(Double_t max) {fVtMax = max;}
    void SetStoreOnlyTrue(Bool_t doit = kTRUE) {fStoreOnlyTrue = doit;}
    void SetForMixing(Bool_t doit = kTRUE) {fForMixing = doit;}
	
	/* working parameters */
	void   AddCutPair(AliRsnDaughterCutPair *cut);
	void   AddCutSingle(Int_t i, AliRsnDaughterCut *cut);
	Stat_t Process(AliRsnEvent *event1, AliRsnEvent *event2 = 0, Bool_t usePID = kTRUE);
	
private:

	/* private functions */
	void   InitHistogram(Int_t nbins, Double_t min, Double_t max);
	Bool_t SingleCutCheck(Int_t ipart, AliRsnDaughter *track) const;
	Bool_t PairCutCheck(AliRsnDaughter *track1, AliRsnDaughter *track2) const;
    Stat_t Fill(TRefArray *list1, TRefArray *list2, Bool_t skipSameIndex = kTRUE);

	/* flags */
	Bool_t               fForMixing;       // flag is true for objects created for event mixing
	Bool_t               fStoreOnlyTrue;   // output = only spectra of true pairs
	
	/* parameters */
	Int_t                fTrueMotherPDG;   // PDG code of true mother (if known)
	Double_t             fMass[2];         // nominal mass of particles
	Char_t               fCharge[2];       // charge of particles
	AliRsnPID::EType     fType[2];         // particles types
    
    /* basic cuts */
    Double_t             fPtMin;           // minimum allowed pt for the pair
    Double_t             fPtMax;           // maximum allowed pt for the pair
    Double_t             fVtMax;           // maximum transverse impact parameter for each track
	
	/* cuts */
	TObjArray           *fCutsSingle[2];   // single-particle cuts
	TObjArray           *fCutsPair;        // pair cuts
	
	/* output */
	TH1D                *fHistogram;       // invariant mass distribution
	
	/* ROOT dictionary */
	ClassDef(AliRsnPair, 1) 
};

#endif
