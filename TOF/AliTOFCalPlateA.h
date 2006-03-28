#ifndef ALITOFCALPLATEA_H
#define ALITOFCALPLATEA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//              class for TOF calibration:: PlateC              //
//////////////////////////////////////////////////////////////////

//_____________________________________________________________

#include "TObject.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TClass.h"
#include "AliTOFGeometry.h"
#include "AliTOFChannel.h"


class AliTOFCalPlateA: public TObject 
{
 public:
  AliTOFCalPlateA();
  AliTOFCalPlateA(AliTOFGeometry *geom);
  AliTOFCalPlateA(AliTOFChannel *ch);
  AliTOFCalPlateA(AliTOFGeometry *geom,AliTOFChannel *ch);
  AliTOFCalPlateA(const AliTOFCalPlateA& pl);
  virtual ~AliTOFCalPlateA();
  Int_t NStripA()const {return fNStripA;}
  Int_t NpadZ()const {return fNpadZ;}
  Int_t NpadX()const {return fNpadX;}
  void Browse(TBrowser *b);
  Bool_t IsFolder() const{return kTRUE;}
private:
  Int_t fNStripA;  // number of TOF strips A
  Int_t fNpadZ;    // number of TOF pads Z
  Int_t fNpadX;    // number of TOF pads X

  AliTOFGeometry *fGeom;    // AliTOFgeometry pointer
  AliTOFChannel *fCh;  //array of AliTOFChannel storing calib parameters
  ClassDef(AliTOFCalPlateA,1)
};

#endif


