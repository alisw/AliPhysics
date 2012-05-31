#ifndef ALIITSDRIFTSPEEDSDD_H
#define ALIITSDRIFTSPEEDSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for SDD drift speed extracted from injector triggers    //
// drift speed dependence on anode number is parametrized via    //
// a polynomial function (3rd degree by default, max. 5th degree)//
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include<TObject.h>
#include<TMath.h>

class AliITSDriftSpeedSDD : public TObject {
 public:
  AliITSDriftSpeedSDD();
  AliITSDriftSpeedSDD(Int_t ev, UInt_t timest, Int_t deg, Double_t *coeff);
  AliITSDriftSpeedSDD(const AliITSDriftSpeedSDD& drSpeed);
  AliITSDriftSpeedSDD& operator=(const AliITSDriftSpeedSDD& drSpeed); 
  virtual ~AliITSDriftSpeedSDD(){};

  static Float_t DefaultDriftSpeed() {return fgkDriftSpeedDefault;}
  virtual Bool_t IsEqual(const TObject *obj) const 
    {return fEvNum == ((AliITSDriftSpeedSDD*)obj)->fEvNum;}
  virtual Bool_t      IsSortable() const { return kTRUE; }
  virtual Int_t       Compare(const TObject *obj) const 
    {if(fEvNum<((AliITSDriftSpeedSDD*)obj)->fEvNum) return -1;
    else if(fEvNum>((AliITSDriftSpeedSDD*)obj)->fEvNum) return 1;
    else return 0; }

  void PrintDriftSpeedParameters() const;
  void SetDegreeofPoly(Int_t deg) {fPolDeg = deg>fgkMaxPolDeg ? fgkMaxPolDeg : deg;}
  Int_t GetDegreeofPoly() const {return fPolDeg;}
  Int_t GetEventNumber() const {return fEvNum;}
  UInt_t GetEventTimestamp() const {return fTimestamp;}
  Float_t GetDriftSpeedParameter(Int_t i) const {return fDriftSpeedParam[i];}
  void    SetDriftSpeedParameter(Int_t i,Float_t par)  {if (i<=fPolDeg) fDriftSpeedParam[i] = par;}
  Double_t GetDriftSpeedAtAnode(Double_t nAnode) const{
    Double_t drSpeed=fDriftSpeedParam[fgkMaxPolDeg];
    for(Int_t i=fgkMaxPolDeg-1; i>=0; --i) drSpeed=fDriftSpeedParam[i]+nAnode*drSpeed;
    return drSpeed;
  }
  static UShort_t GetMaxPolDeg() {return fgkMaxPolDeg;}

 protected:
  static const Float_t fgkDriftSpeedDefault; // default for drift speed
  static const UShort_t fgkMaxPolDeg=5; // max. degree of the poly fit

  Int_t fEvNum;  // event number of injector event
  UInt_t fTimestamp; // event timestamp
  Char_t fPolDeg;    // degree of the ploy fit to drift speed vs. anode
                     // saved as char since 1 byte is enough
  Float_t fDriftSpeedParam[fgkMaxPolDeg+1];  // coefficients of the poly fit
  ClassDef(AliITSDriftSpeedSDD,4);
};
#endif
