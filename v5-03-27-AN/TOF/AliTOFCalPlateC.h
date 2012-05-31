#ifndef ALITOFCALPLATEC_H
#define ALITOFCALPLATEC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//              class for TOF calibration:: PlateC              //
//////////////////////////////////////////////////////////////////

//_____________________________________________________________

#include "TObject.h"

class TBrowser;

class AliTOFChannel;
class AliTOFGeometry;

class AliTOFCalPlateC: public TObject 
{
 public:
  AliTOFCalPlateC();
  AliTOFCalPlateC(AliTOFGeometry *geom);
  AliTOFCalPlateC(AliTOFChannel *ch);
  AliTOFCalPlateC(AliTOFGeometry *geom,AliTOFChannel *ch);
  AliTOFCalPlateC(const AliTOFCalPlateC& pl);
  AliTOFCalPlateC& operator=(const AliTOFCalPlateC &source); // ass. op.
  virtual ~AliTOFCalPlateC();
  Int_t NStripC()const {return fNStripC;}
  Int_t NpadZ()const {return fNpadZ;}
  Int_t NpadX()const {return fNpadX;}
  void Browse(TBrowser *b);
  Bool_t IsFolder() const{return kTRUE;}
private:
  Int_t fNStripC;  // number of TOF strips C
  Int_t fNpadZ;    // number of TOF pads Z
  Int_t fNpadX;    // number of TOF pads X

  AliTOFGeometry *fGeom;    // AliTOFgeometry pointer
  AliTOFChannel *fCh;  //array of AliTOFChannel storing calib parameters
  ClassDef(AliTOFCalPlateC,1)
};

#endif


