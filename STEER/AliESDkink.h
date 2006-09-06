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
  Float_t GetTPCDensityFactor() const;
  Float_t GetQt() const;    
  //
  Float_t GetR() const {return fRr;}
  Float_t GetDistance() const {return fDist2;}
  Int_t   GetTPCRow0() const {return fRow0;}
  Float_t GetAngle(Int_t i) const {return fAngle[i];}
  const Double_t *GetPosition() const   {return fXr;}
  const Double_t *GetMotherP()  const   {return fPm;}
  const Double_t *GetDaughterP()  const {return fPdr;}
  void SetTPCRow0(Int_t row0){fRow0 = row0;}
  Int_t GetLabel(Int_t i) const {return fLab[i];}
  void SetLabel(Int_t label, Int_t pos) {fLab[pos]=label;}
  Int_t GetIndex(Int_t i) const {return fIndex[i];}
  void SetIndex(Int_t index, Int_t pos){fIndex[pos]=index;}
  void SetStatus(Int_t status, Int_t pos){fStatus[pos]=status;}
  Int_t GetStatus(Int_t pos) const {return fStatus[pos];}
  void SetTPCncls(Int_t ncls,Int_t pos) {fTPCncls[pos]=ncls;}
  const Int_t *GetTPCncls() const {return fTPCncls;} 
  void  SetTPCDensity(Float_t dens, Int_t pos0,Int_t pos1){fTPCdensity[pos0][pos1]=dens;}
  Float_t GetTPCDensity(Int_t pos0,Int_t pos1) const {return fTPCdensity[pos0][pos1];}
  void    SetTPCDensity2(Float_t dens, Int_t pos0,Int_t pos1){fTPCdensity[pos0][pos1]=dens;}
  Float_t GetTPCDensity2(Int_t pos0,Int_t pos1) const {return fTPCdensity[pos0][pos1];}
  Float_t GetShapeFactor() const {return fShapeFactor;}
  void    SetShapeFactor(Float_t factor){fShapeFactor = factor;}
  void  SetMultiple(Int_t mult,Int_t pos){fMultiple[pos]=mult;}
  const Int_t * GetMultiple() const {return fMultiple;}
  //  
  const AliExternalTrackParam& RefParamDaughter() {return fParamDaughter;}
  const AliExternalTrackParam& RefParamMother()   {return fParamMother;}
 protected:
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
  Char_t         fStatus[12];       //status of kink - first 4 mother (ITS,TPC,TRD,TOF)  other daughter
  Float_t        fTPCdensity[2][2];  //tpc cluster density before and after kink
  Float_t        fTPCdensity2[2][2];  //tpc cluster density before and after kink - after second iteration
  Float_t        fShapeFactor;       // tpc clusters shape factor
  Int_t          fRow0;              // critical pad row number
  Int_t          fMultiple[2];       //how many times the track's were used
  Int_t          fTPCncls[2];     //number of clusters for mother particle
  ClassDef(AliESDkink,2)      // ESD V0 vertex
};

#endif


