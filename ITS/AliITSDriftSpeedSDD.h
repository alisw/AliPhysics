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
  virtual ~AliITSDriftSpeedSDD(){};

  virtual Bool_t IsEqual(const TObject *obj) const 
    {return fEvNum == ((AliITSDriftSpeedSDD*)obj)->fEvNum;}
  virtual Bool_t      IsSortable() const { return kTRUE; }
  virtual Int_t       Compare(const TObject *obj) const 
    {if(fEvNum<((AliITSDriftSpeedSDD*)obj)->fEvNum) return -1;
    else if(fEvNum>((AliITSDriftSpeedSDD*)obj)->fEvNum) return 1;
    else return 0; }

  void PrintDriftSpeedParameters() const;

  Int_t GetEventNumber() const {return fEvNum;}
  UInt_t GetEventTimestamp() const {return fTimestamp;}
  Double_t GetDriftSpeedParameter(Int_t i) const {return fDriftSpeedParam[i];}
  Double_t GetDriftSpeedAtAnode(Double_t nAnode) const{
    Double_t drSpeed=fDriftSpeedParam[fgkMaxPolDeg];
    for(Int_t i=fgkMaxPolDeg-1; i>=0; --i) drSpeed=fDriftSpeedParam[i]+nAnode*drSpeed;
    return drSpeed;
  }

 protected:
  static const Int_t fgkMaxPolDeg=5; // max. degree of the poly fit
  Int_t fEvNum;  // event number of injector event
  UInt_t fTimestamp; // event timestamp
  Int_t fPolDeg;    // degree of the ploy fit to drift speed vs. anode
  Double32_t fDriftSpeedParam[fgkMaxPolDeg+1];  // coefficients of the poly fit
  ClassDef(AliITSDriftSpeedSDD,2);
};
#endif
