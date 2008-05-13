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
	virtual Bool_t Pass(AliRsnDaughter *track) const;
private:
	ClassDef(AliRsnDaughterCut,1)
};
//
//-------------------------------------------------------------------------
//
class AliRsnDaughterCutPt : public AliRsnDaughterCut
{
public:
			       AliRsnDaughterCutPt(Double_t min, Double_t max) : fPtMin(min),fPtMax(max) {}
	virtual       ~AliRsnDaughterCutPt() {}
	virtual Bool_t Pass(AliRsnDaughter *track) const;
private:
	Double_t fPtMin; // smallest allowed Pt
	Double_t fPtMax; // largest allowed Pt
	
	ClassDef(AliRsnDaughterCutPt,1)
};
//
//-------------------------------------------------------------------------
//
class AliRsnDaughterCutVt : public AliRsnDaughterCut
{
public:
			       AliRsnDaughterCutVt(Double_t max) : fVtMax(max) {}
	virtual       ~AliRsnDaughterCutVt() {}
	virtual Bool_t Pass(AliRsnDaughter *track) const;
private:
	Double_t fVtMax; // largest allowed transverse impact parameter
	
	ClassDef(AliRsnDaughterCutVt,1)
};

#endif
