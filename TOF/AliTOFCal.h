#ifndef ALITOFCAL_H
#define ALITOFCAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//     class for TOF calibration:: array of AliTOFChannels      //
//////////////////////////////////////////////////////////////////

//_____________________________________________________________

#include "TObject.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TClass.h"
#include "AliTOFGeometryV4.h"
#include "AliTOFChannel.h"


class AliTOFCal: public TObject 
{
 public:
  AliTOFCal();
  AliTOFCal(const AliTOFCal& cal);
  virtual ~AliTOFCal();
  void Browse(TBrowser *b);
  Bool_t IsFolder() const{return kTRUE;}
  Int_t NSector()const {return fNSector;}
  Int_t NPlate()const {return fNPlate;}
  Int_t NStripA()const {return fNStripA;}
  Int_t NStripB()const {return fNStripB;}
  Int_t NStripC()const {return fNStripC;}
  Int_t NpadZ()const {return fNpadZ;}
  Int_t NpadX()const {return fNpadX;}
  Int_t NPads()const {return fnpad;}
 //void PlotPad(Int_t n) {} // *MENU*
  AliTOFChannel* GetChannel(Int_t i) {return &fPads[i];}
  AliTOFChannel* GetArray() {return fPads;}
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

  AliTOFChannel* fPads;  //[fnpad]  
                         // array of AliTOFChannels storing the calib parameters
  ClassDef(AliTOFCal,1)
};

#endif


