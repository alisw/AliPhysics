#ifndef ALIESDKINK_H
#define ALIESDKINK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          ESD V0 Vertex Class
//          This class is part of the Event Summary Data set of classes
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>
#include "AliExternalTrackParam.h"
#include <TPDGCode.h>

class AliESDtrack;

class AliESDkink : public TObject {
public:
  AliESDkink();             //constructor
  //
  void SetID(Int_t id){fID=id;}
  Int_t GetID(){return fID;}
  void SetMother(const AliExternalTrackParam & pmother); 
  void SetDaughter(const AliExternalTrackParam & pdaughter);
  void Update();            //update
  Float_t GetTPCDensityFactor() const;
  Float_t GetQt() const;    
  //  
  Int_t          fID;       // kink ID
  AliExternalTrackParam fParamDaughter;
  AliExternalTrackParam fParamMother;
  Double_t       fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t       fDist2;    //info about closest distance parabolic DCA
  //
  Double_t       fPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t       fXr[3];     //rec. position according helix
  //
  Double_t       fPm[3];    //momentum at the vertex mother
  Double_t       fAngle[3]; //three angles
  Double_t       fRr;       // rec position of the vertex 
  Int_t          fLab[2];   //MC label of the partecle
  Int_t          fIndex[2]; //reconstructed labels of the tracks
  Int_t          fStatus;       //status 
  Float_t        fTPCdensity[2][2];  //tpc cluster density before and after kink
  Float_t        fTPCdensity2[2][2];  //tpc cluster density before and after kink - after second iteration
  Float_t        fShapeFactor;       // tpc clusters shape factor
  Int_t          fRow0;              // critical pad row number
  Int_t          fMultiple[2];       //how many times the track's were used
  Float_t        fZm[2];                // z at the middle of chamber
  Float_t        fFi[2];
  ClassDef(AliESDkink,1)      // ESD V0 vertex
};

#endif


