#ifndef ALITRACKERBASE_H
#define ALITRACKERBASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTrackerBase.h 38069 2009-12-24 16:56:18Z belikov $ */

//-------------------------------------------------------------------------
//                          class AliTrackerBase
//        This is the base tracker class, independent on STEER 
//       Origin: Iouri Belikov, IPHC, iouri.belikov@in2p3.fr 
//    Most of the functions implemented by Marian.Ivanov@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TGeoGlobalMagField.h>

#include "AliMagF.h"

class AliExternalTrackParam;

class AliTrackerBase : public TObject {
public:
  AliTrackerBase();
  virtual ~AliTrackerBase(){}

  void SetVertex(const Double_t *xyz, const Double_t *ers=0) { 
     fX=xyz[0]; fY=xyz[1]; fZ=xyz[2];
     if (ers) { fSigmaX=ers[0]; fSigmaY=ers[1]; fSigmaZ=ers[2]; } 
  }
  Double_t GetX() const {return fX;}
  Double_t GetY() const {return fY;}
  Double_t GetZ() const {return fZ;}
  Double_t GetSigmaX() const {return fSigmaX;}
  Double_t GetSigmaY() const {return fSigmaY;}
  Double_t GetSigmaZ() const {return fSigmaZ;}

  static Double_t GetTrackPredictedChi2(AliExternalTrackParam *track,
                                        Double_t mass, Double_t step, 
			          const AliExternalTrackParam *backup);
  static 
  Double_t MeanMaterialBudget(const Double_t *start, const Double_t *end, 
  Double_t *mparam);
  static
  Bool_t PropagateTrackTo(AliExternalTrackParam *track, Double_t x, Double_t m,
			  Double_t maxStep, Bool_t rotateTo=kTRUE, Double_t maxSnp=0.8, Double_t sign=1.);  
  static Bool_t PropagateTrackToBxByBz(AliExternalTrackParam *track, Double_t x, 
         Double_t m,
				Double_t maxStep, Bool_t rotateTo=kTRUE, Double_t maxSnp=0.8,Double_t sign=1.);  
  //
  static Double_t GetBz(const Double_t *r);
  static void GetBxByBz(const Double_t r[3], Double_t b[3]);
  static Double_t GetBz();
  static Bool_t UniformField(); 
  //

protected:
  AliTrackerBase(const AliTrackerBase &atr);
private:
  AliTrackerBase & operator=(const AliTrackerBase & atr);

  Double_t fX;  //X-coordinate of the primary vertex
  Double_t fY;  //Y-coordinate of the primary vertex
  Double_t fZ;  //Z-coordinate of the primary vertex
 
  Double_t fSigmaX; // error of the primary vertex position in X
  Double_t fSigmaY; // error of the primary vertex position in Y
  Double_t fSigmaZ; // error of the primary vertex position in Z
  
  ClassDef(AliTrackerBase,1) //base tracker
};

//__________________________________________________________________________
inline Bool_t AliTrackerBase::UniformField()
{
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  return fld ? fld->IsUniform():kTRUE;
}

#endif
