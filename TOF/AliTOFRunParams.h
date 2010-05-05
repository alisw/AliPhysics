#ifndef ALITOFRUNPARAMS_H
#define ALITOFRUNPARAMS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

// *
// *
// *
// * this class defines the TOF object to be stored
// * in OCDB on a run-by-run basis in order to have the measurement
// * of the time evolution of T0 and of TOF resolution including
// * average T0 uncertainty
// *
// *
// *

#include "TObject.h"

class AliTOFRunParams :
public TObject
{

 public:

  AliTOFRunParams(); // default constructor
  AliTOFRunParams(Int_t nPoints); // standard constructor
  virtual ~AliTOFRunParams(); // default destructor
  AliTOFRunParams(const AliTOFRunParams &source); // copy constructor
  AliTOFRunParams &operator=(const AliTOFRunParams &source); // operator=

  Int_t GetNPoints() const {return fNPoints;}; // getter
  UInt_t GetTimestamp(Int_t i) const {return fTimestamp && i < fNPoints ? fTimestamp[i] : 0;}; // getter
  Float_t GetT0(Int_t i) const {return fT0 && i < fNPoints ? fT0[i] : 0.;}; // getter
  Float_t GetTOFResolution(Int_t i) const {return fTOFResolution && i < fNPoints ? fTOFResolution[i] : 0.;}; // getter
  Float_t GetT0Spread(Int_t i) const {return fT0Spread && i < fNPoints ? fT0Spread[i] : 0.;}; // getter
  
  void SetTimestamp(UInt_t *value) {if (fTimestamp) for (Int_t i = 0; i < fNPoints; i++) fTimestamp[i] = value[i];}; // setter
  void SetT0(Float_t *value) {if (fT0) for (Int_t i = 0; i < fNPoints; i++) fT0[i] = value[i];}; // setter
  void SetTOFResolution(Float_t *value) {if (fTOFResolution) for (Int_t i = 0; i < fNPoints; i++) fTOFResolution[i] = value[i];}; // setter
  void SetT0Spread(Float_t *value) {if (fT0Spread) for (Int_t i = 0; i < fNPoints; i++) fT0Spread[i] = value[i];}; // setter

  Float_t EvalT0(UInt_t timestamp); // eval T0
  Float_t EvalTOFResolution(UInt_t timestamp); // eval TOF resolution
  Float_t EvalT0Spread(UInt_t timestamp); // eval T0 spread

 private:

  Int_t fNPoints;
  UInt_t *fTimestamp; //[fNPoints] time stamp
  Float_t *fT0; //[fNPoints] average T0 (ps)
  Float_t *fTOFResolution; //[fNPoints] average TOF resolution (T0 uncertainty included) (ps)
  Float_t *fT0Spread; //[fNPoints] estimated T0 spread (from vertex spread z) (ps)

  ClassDef(AliTOFRunParams, 1);
};

#endif /* ALITOFRUNPARAMS_H */
