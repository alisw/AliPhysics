#ifndef ALITOFCALPLATEC_H
#define ALITOFCALPLATEC_H

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
#include "AliTOFGeometryV4.h"
#include "AliTOFChannel.h"

class AliTOFCalPlateC: public TObject 
{
 public:
  AliTOFCalPlateC();
  AliTOFCalPlateC(AliTOFChannel *ch);
  AliTOFCalPlateC(const AliTOFCalPlateC& pl);
  virtual ~AliTOFCalPlateC();
  Int_t NSector()const {return fNSector;}
  Int_t NPlate()const {return fNPlate;}
  Int_t NStripA()const {return fNStripA;}
  Int_t NStripB()const {return fNStripB;}
  Int_t NStripC()const {return fNStripC;}
  Int_t NpadZ()const {return fNpadZ;}
  Int_t NpadX()const {return fNpadX;}
  void Browse(TBrowser *b);
  Bool_t IsFolder() const{return kTRUE;}
private:
  Int_t fNSector;  // number of TOF sectors
  Int_t fNPlate;   // number of TOF plates
  Int_t fNStripA;  // number of TOF strips A
  Int_t fNStripB;  // number of TOF strips B
  Int_t fNStripC;  // number of TOF strips C
  Int_t fNpadZ;    // number of TOF pads Z
  Int_t fNpadX;    // number of TOF pads X

  AliTOFChannel *fCh; //array of AliTOFChannel storing calib parameters
  ClassDef(AliTOFCalPlateC,1)
};

#endif


