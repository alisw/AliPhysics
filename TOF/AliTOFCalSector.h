#ifndef ALITOFCALSECTOR_H
#define ALITOFCALSECTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//              class for TOF calibration:: Sectors             //
//////////////////////////////////////////////////////////////////

//_____________________________________________________________

#include "TObject.h"

class TBrowser;
class AliTOFGeometry;
class AliTOFChannel;

class AliTOFCalSector: public TObject 
{
public:
  AliTOFCalSector();
  AliTOFCalSector(AliTOFChannel *ch);
  AliTOFCalSector(AliTOFGeometry *geom);
  AliTOFCalSector(AliTOFGeometry *geom, AliTOFChannel *ch);
  AliTOFCalSector(const AliTOFCalSector& sec);
  AliTOFCalSector& operator=(const AliTOFCalSector &source); // ass. op.
  virtual ~AliTOFCalSector();
  Int_t NPlate()const {return fNPlate;}
  Int_t NStripA()const {return fNStripA;}
  Int_t NStripB()const {return fNStripB;}
  Int_t NStripC()const {return fNStripC;}
  Int_t NpadZ()const {return fNpadZ;}
  Int_t NpadX()const {return fNpadX;}
  void Browse(TBrowser *b);
  Bool_t IsFolder() const{return kTRUE;}
private:
  Int_t fNPlate;   // number of TOF plates
  Int_t fNStripA;  // number of TOF strips A
  Int_t fNStripB;  // number of TOF strips B
  Int_t fNStripC;  // number of TOF strips C
  Int_t fNpadZ;    // number of TOF pads Z
  Int_t fNpadX;    // number of TOF pads X

  AliTOFGeometry *fGeom;    // AliTOFgeometry pointer
  AliTOFChannel *fCh; //array of AliTOFChannel storing calib parameters
  ClassDef(AliTOFCalSector,1)
};

#endif
