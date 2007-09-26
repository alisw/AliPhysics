#ifndef ALIITSPIDESD1_H
#define ALIITSPIDESD1_H
/* Copyright(c) 2005-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
// PID method # 1
//                    ITS PID class
// A very naive design... Should be made better by the detector experts...
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------
#include "AliITSpidESD.h"
#include "AliLog.h"

class AliITSpidESD1 : public AliITSpidESD {
public:
  AliITSpidESD1();
  AliITSpidESD1(Double_t *param);
  virtual ~AliITSpidESD1() {}
  virtual Int_t MakePID(AliESDEvent *event);
  // This method is here because of the AliITSpidESD2 class...
  virtual Int_t MakePID(TTree */*clustersTree*/, AliESDEvent */*event*/)
  {AliError("Method should not be used!"); return 0;}

private:
  Double_t fRes;          // relative dEdx resolution
  Double_t fRange;        // one particle type PID range (in sigmas)
  ClassDef(AliITSpidESD1,2)   // ITS PID class
};

#endif

