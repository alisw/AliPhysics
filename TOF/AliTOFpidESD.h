#ifndef ALITOFPIDESD_H
#define ALITOFPIDESD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TOF PID class
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include "TObject.h"

class AliESD;

class AliTOFGeometry;

class AliTOFpidESD : public TObject {
enum {kMaxCluster=77777}; //maximal number of the TOF clusters
public:
  AliTOFpidESD(){fN=0; fEventN=0;}
  AliTOFpidESD(Double_t *param);
  ~AliTOFpidESD(){}

  Int_t MakePID(AliESD *event);
  void  SetEventNumber(Int_t n) {fEventN=n;}
  Int_t GetEventNumber() const {return fEventN;}

private:
 
  Int_t fN;               // number of the TOF clusters
  Int_t fEventN;          // event number
  Double_t fSigma;        // intrinsic TOF resolution
  Double_t fRange;        // one particle type PID range (in sigmas)

  ClassDef(AliTOFpidESD,1)   // TOF PID class
};

#endif
