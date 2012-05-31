#ifndef ALITOFCALONLINE_H
#define ALITOFCALONLINE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////////
//     class for TOF online calibration:: array of AliTOFChannelOnline      //
/////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________

#include "TObject.h"
#include "AliTOFChannelOnline.h"

class TBrowser;
class AliTOFGeometry;

class AliTOFCalOnline: public TObject 
{
 public:
  AliTOFCalOnline();
  AliTOFCalOnline(AliTOFGeometry *geom);
  AliTOFCalOnline(const AliTOFCalOnline& cal);
  AliTOFCalOnline& operator=(const AliTOFCalOnline &source); // ass. op.
  virtual ~AliTOFCalOnline();
  Int_t NSector()const {return fNSector;}
  Int_t NPlate()const {return fNPlate;}
  Int_t NStripA()const {return fNStripA;}
  Int_t NStripB()const {return fNStripB;}
  Int_t NStripC()const {return fNStripC;}
  Int_t NpadZ()const {return fNpadZ;}
  Int_t NpadX()const {return fNpadX;}
  Int_t NPads()const {return fnpad;}
  AliTOFChannelOnline* GetChannel(Int_t i) {return &fPads[i];}
  AliTOFChannelOnline* GetArray() {return fPads;}
  void CreateArray();
private:

  Int_t fNSector;  // number of TOF sectors
  Int_t fNPlate;   // number of TOF platess
  Int_t fNStripA;  // number of TOF strips A
  Int_t fNStripB;  // number of TOF strips B
  Int_t fNStripC;  // number of TOF strips C
  Int_t fNpadZ;    // number of TOF pads Z
  Int_t fNpadX;    // number of TOF pads X
  Int_t fnpad;     // number of TOF channels

  AliTOFGeometry *fGeom;    // AliTOFgeometry pointer
  AliTOFChannelOnline* fPads;  //[fnpad]  
                               // array of AliTOFChannels storing 
                               // the calib parameters
  ClassDef(AliTOFCalOnline,1)
};

#endif


