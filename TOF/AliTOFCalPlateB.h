#ifndef ALITOFCALPLATEB_H
#define ALITOFCALPLATEB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//              class for TOF calibration:: PlateB              //
//////////////////////////////////////////////////////////////////

//_____________________________________________________________

#include "TObject.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TClass.h"
#include "AliTOFGeometry.h"
#include "AliTOFChannel.h"


class AliTOFCalPlateB: public TObject 
{
 public:
  AliTOFCalPlateB();
  AliTOFCalPlateB(AliTOFGeometry *geom);
  AliTOFCalPlateB(AliTOFChannel *ch);
  AliTOFCalPlateB(AliTOFGeometry *geom,AliTOFChannel *ch);
  AliTOFCalPlateB(const AliTOFCalPlateB& pl);
  AliTOFCalPlateB& operator=(const AliTOFCalPlateB &source); // ass. op.
  virtual ~AliTOFCalPlateB();
  Int_t NStripB()const {return fNStripB;}
  Int_t NpadZ()const {return fNpadZ;}
  Int_t NpadX()const {return fNpadX;}
  void Browse(TBrowser *b);
  Bool_t IsFolder() const{return kTRUE;}
private:
  Int_t fNStripB;  // number of TOF strips B
  Int_t fNpadZ;    // number of TOF pads Z
  Int_t fNpadX;    // number of TOF pads X

  AliTOFGeometry *fGeom;    // AliTOFgeometry pointer
  AliTOFChannel *fCh;  //array of AliTOFChannel storing calib parameters
  ClassDef(AliTOFCalPlateB,1)
};

#endif


