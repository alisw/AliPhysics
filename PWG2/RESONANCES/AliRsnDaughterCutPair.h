/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnDaughterCutPair
//  
//           Implementation of various cut which can be applied
//           during resonance analysis. 
//           First is the virtual base class.
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNDAUGHTERCUTPAIR_H
#define ALIRSNDAUGHTERCUTPAIR_H

class AliRsnDaughter;

class AliRsnDaughterCutPair : public TObject
{
public:
   AliRsnDaughterCutPair() {}
   virtual        ~AliRsnDaughterCutPair() {}
   virtual Bool_t  Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const;
private:
   ClassDef(AliRsnDaughterCutPair,1)
};
//	
//-------------------------------------------------------------------------
//
class AliRsnDaughterCutPairPt : public AliRsnDaughterCutPair
{
public:
   AliRsnDaughterCutPairPt(Double_t min, Double_t max) : fPtMin(min),fPtMax(max) {}
   virtual       ~AliRsnDaughterCutPairPt() {}
   virtual Bool_t Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const;
protected:
   Double_t fPtMin; // smallest allowed Pt
   Double_t fPtMax; // largest allowed Pt
	
   ClassDef(AliRsnDaughterCutPairPt,1)
};
//	
//-------------------------------------------------------------------------
//
class AliRsnDaughterCutPairAngle : public AliRsnDaughterCutPair
{
public:
			       AliRsnDaughterCutPairAngle(Double_t min, Double_t max) : fAngleMin(min),fAngleMax(max) {}
	virtual       ~AliRsnDaughterCutPairAngle() {}
	virtual Bool_t Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const;
protected:
	Double_t fAngleMin; // smallest allowed relative angle
	Double_t fAngleMax; // largest allowed relative angle
	
	ClassDef(AliRsnDaughterCutPairAngle,1)
};
//
//-------------------------------------------------------------------------
//
class AliRsnDaughterCutPairArmenteros : public AliRsnDaughterCutPair
{
public:
   AliRsnDaughterCutPairArmenteros(Double_t qMin, Double_t qMax, Double_t aMin, Double_t aMax) :
     fQtMin(qMin),fQtMax(qMax),fAlphaMin(aMin),fAlphaMax(aMax) {}
   virtual       ~AliRsnDaughterCutPairArmenteros() {}
   virtual Bool_t Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const;
   void   Compute(AliRsnDaughter *track1, AliRsnDaughter *track2, Double_t &qt, Double_t &alpha) const;
private:
   Double_t fQtMin;    // minimum Qt
   Double_t fQtMax;    // minimum Qt
   Double_t fAlphaMin; // minimum alpha
   Double_t fAlphaMax; // minimum alpha
	
   ClassDef(AliRsnDaughterCutPairArmenteros,1)
};

#endif
