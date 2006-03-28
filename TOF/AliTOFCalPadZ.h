#ifndef ALITOFCALPADZ_H
#define ALITOFCALPADZ_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//              class for TOF calibration:: PadZ                //
//////////////////////////////////////////////////////////////////

//_____________________________________________________________

#include "TObject.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TClass.h"
#include "AliTOFGeometry.h"
#include "AliTOFChannel.h"


class AliTOFCalPadZ: public TObject 
{
 public:
  AliTOFCalPadZ();
  AliTOFCalPadZ(AliTOFChannel *ch);
  AliTOFCalPadZ(AliTOFGeometry *geom);
  AliTOFCalPadZ(AliTOFGeometry *geom,AliTOFChannel *ch);
  virtual ~AliTOFCalPadZ();
  Int_t NpadX()const {return fNpadX;}
  void Browse(TBrowser *b);
  Bool_t IsFolder() const{return kTRUE;}
private:
  Int_t fNpadX;    // number of TOF pads X

  AliTOFGeometry *fGeom;    // AliTOFgeometry pointer
  AliTOFChannel *fCh; //array of AliTOFChannel storing calib parameters
  ClassDef(AliTOFCalPadZ,1)
};

#endif //AliTOFPadZ


