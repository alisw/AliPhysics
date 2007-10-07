// -*- mode: C++ -*- 
#ifndef ALIESDTZERO_H
#define ALIESDTZERO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//                          Class AliESDTZERO
//   This is a class that summarizes the TZERO data for the ESD   
//   Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------



#include <TObject.h>

class AliESDTZERO: public TObject {
public:
  AliESDTZERO();
  AliESDTZERO(const AliESDTZERO& tzero);
  AliESDTZERO& operator=(const AliESDTZERO& tzero);

  Double_t GetT0zVertex() const {return fT0zVertex;}
  void SetT0zVertex(Double_t z) {fT0zVertex=z;}
  Double_t GetT0() const {return fT0timeStart;}
  void SetT0(Double_t timeStart) {fT0timeStart = timeStart;}
  const Double_t * GetT0time() const {return fT0time;}
  void SetT0time(Float_t time[24]) {
    for (Int_t i=0; i<24; i++) fT0time[i] = time[i];
  }
  const Double_t * GetT0amplitude() const {return fT0amplitude;}
  void SetT0amplitude(Float_t amp[24]) {
    for (Int_t i=0; i<24; i++) fT0amplitude[i] = amp[i];
  }

  void    Reset();
  void    Print(const Option_t *opt=0) const;

private:

  Double32_t      fT0zVertex;       // vertex z position estimated by the T0
  Double32_t      fT0timeStart;     // interaction time estimated by the T0
  Double32_t      fT0time[24];      // best TOF on each T0 PMT
  Double32_t      fT0amplitude[24]; // number of particles(MIPs) on each T0 PMT

  ClassDef(AliESDTZERO,2)
};


#endif
