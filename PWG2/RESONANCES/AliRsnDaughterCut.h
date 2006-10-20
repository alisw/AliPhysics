/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnDaughterCut
//  
//           Implementation of various cut which can be applied
//           during resonance analysis. 
//           First is the virtual base class.
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNDAUGHTERCUT_H
#define ALIRSNDAUGHTERCUT_H

class AliRsnDaughter;

class AliRsnDaughterCut : public TObject
{
public:
			       AliRsnDaughterCut() {  }
	virtual       ~AliRsnDaughterCut() {  }
	
	        Bool_t IsPairCut() const {return fPairCut;}
	virtual Bool_t Pass(AliRsnDaughter *track1, AliRsnDaughter *track2 = 0) const;
	
protected:

	Bool_t  fPairCut;  // this flag is TRUE for all pair cuts
	
	ClassDef(AliRsnDaughterCut,1)
};

//-------------------------------------------------------------------------

class AliRsnDaughterCutPtSingle : public AliRsnDaughterCut
{
public:
			       AliRsnDaughterCutPtSingle(Double_t min, Double_t max) {fPairCut=kFALSE;fPtMin=min;fPtMax=max;}
	virtual       ~AliRsnDaughterCutPtSingle()                           { }
	
	virtual Bool_t Pass(AliRsnDaughter *track1, AliRsnDaughter *track2 = 0) const;

protected:

	Double_t fPtMin; // smallest allowed Pt
	Double_t fPtMax; // largest allowed Pt
	
	ClassDef(AliRsnDaughterCutPtSingle,1)
};
	
//-------------------------------------------------------------------------

class AliRsnDaughterCutPtPair : public AliRsnDaughterCut
{
public:
			       AliRsnDaughterCutPtPair(Double_t min, Double_t max) {fPairCut=kTRUE;fPtMin=min;fPtMax=max;}
	virtual       ~AliRsnDaughterCutPtPair()                           { }
	
	virtual Bool_t Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const;

protected:

	Double_t fPtMin; // smallest allowed Pt
	Double_t fPtMax; // largest allowed Pt
	
	ClassDef(AliRsnDaughterCutPtPair,1)
};

//-------------------------------------------------------------------------

#endif
