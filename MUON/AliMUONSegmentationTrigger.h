#ifndef ALIMUONSEGMENTATIONTRIGGER_H
#define ALIMUONSEGMENTATIONTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMUONSegmentationV0.h"
class AliMUONChamber;
//----------------------------------------------
//
// Chamber segmentation virtual base class
//
class AliMUONSegmentationTrigger :
public AliMUONSegmentationV0 {
 public:
  AliMUONSegmentationTrigger(){};
  virtual ~AliMUONSegmentationTrigger(){}  
  virtual void Init(AliMUONChamber* chamber);         // Initialization
  Int_t ModuleNumber(Int_t imodule);  // returns module number of ModuleId
  // Set pad position -> in SegRes X & Y
  //       virtual void     SetPad(Int_t, Int_t);
  // Set hit position
  virtual void     SetHit(Float_t xhit, Float_t yhit);
  virtual void     SetHit(Float_t xhit, Float_t yhit, Float_t zhit)
      {SetHit(xhit, yhit);}
  
  ClassDef(AliMUONSegmentationTrigger,1) //Segmentation class for trigger
    protected:
  Float_t StripSizeX(Int_t imodule);
  Float_t StripSizeY(Int_t imodule);    
 protected:
  Float_t fYcmin[126];       // y min position of modules
  Float_t fYcmax[126];       // y max position of modules
  Float_t fZscale;            // scaling factor (Zx/Z1, x=1,2,3,4)

// Current pad during integration (cursor for disintegration)
  Int_t fix;  // pad coord.  x 
  Int_t fiy;  // pad coord.  y 
  Float_t fx; // real coord. x
  Float_t fy; // real ccord. y
  
  Float_t fxhit;  // x-position of hit
  Float_t fyhit;  // y-position of hit
  Int_t   fSector;// Segmentation Sector
  
};

#endif













