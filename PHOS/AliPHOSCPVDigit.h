#ifndef ALIPHOSCPVDIGIT_H
#define ALIPHOSCPVDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Digit class for CPV                       //
//                                            //
//  Author: Yuri Kharlov, IHEP, Protvino      //
//  e-mail: Yuri.Kharlov@cern.ch              //
//  Last modified: 23 March 2000              //
////////////////////////////////////////////////
 
// --- ROOT system ---
#include <TObject.h> 

class AliPHOSCPVDigit : public TObject {
  
public:
  virtual ~AliPHOSCPVDigit() {}
           AliPHOSCPVDigit() {}
           AliPHOSCPVDigit(Int_t x, Int_t y, Float_t q);
  
  void     SetQpad(Float_t q) { fQpad = q;     }
  Int_t    GetXpad()          { return  fXpad; }
  Int_t    GetYpad()          { return  fYpad; }
  Float_t  GetQpad()          { return  fQpad; }

private:
  Int_t    fXpad;       // Digit's pad number in Phi
  Int_t    fYpad;       // Digit's pad number in Z
  Float_t  fQpad;       // Digit's pad amplitude
  
  ClassDef(AliPHOSCPVDigit,1)  // Digit object in one CPV pad
};
 
#endif
