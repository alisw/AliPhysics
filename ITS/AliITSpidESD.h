#ifndef ALIITSpIDESD_H
#define ALIITSpIDESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                    ITS PID class
// Base class:
// See the implementations AliITSpidESD1 and AliITSpidESD2
//-------------------------------------------------------
//#include <Rtypes.h>
#include <TObject.h>

class TTree;
class AliESDEvent;

class AliITSpidESD : public TObject {
public:
  AliITSpidESD();
  virtual ~AliITSpidESD() {}
  virtual Int_t MakePID(AliESDEvent *event) = 0;
  // This method is here because of the AliITSpidESD2 class...
  virtual Int_t MakePID(TTree *clustersTree, AliESDEvent *event) = 0;
  static Double_t Bethe(Double_t p,Double_t mass);
private:
  ClassDef(AliITSpidESD,1)   // ITS PID class
};

#endif


