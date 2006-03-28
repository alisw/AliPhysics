#ifndef ALITOFCALSTRIP_H
#define ALITOFCALSTRIP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//              class for TOF calibration:: strips              //
//////////////////////////////////////////////////////////////////

//_____________________________________________________________

#include "TObject.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TClass.h"
#include "AliTOFGeometry.h"
#include "AliTOFChannel.h"

class AliTOFCalStrip: public TObject 
{
 public:
  AliTOFCalStrip(AliTOFGeometry *geom);
  AliTOFCalStrip(AliTOFGeometry *geom,AliTOFChannel *ch);
  AliTOFCalStrip();
  AliTOFCalStrip(AliTOFChannel *ch);
  AliTOFCalStrip(const AliTOFCalStrip& strip);
  virtual ~AliTOFCalStrip();
  Int_t NpadZ()const {return fNpadZ;}
  Int_t NpadX()const {return fNpadX;}
  void Browse(TBrowser *b);
  Bool_t IsFolder() const{return kTRUE;}
private:
  Int_t fNpadZ;    // number of TOF pads Z
  Int_t fNpadX;    // number of TOF pads X

  AliTOFGeometry *fGeom;    // AliTOFgeometry pointer
  AliTOFChannel *fCh; //array of AliTOFChannel storing calib parameters
  ClassDef(AliTOFCalStrip,1)
};

#endif


