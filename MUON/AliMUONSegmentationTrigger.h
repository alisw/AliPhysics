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
  
  // Current Pad during Integration
  // x-coordinate
  //    virtual Int_t  Ix();
  // y-coordinate
  //    virtual Int_t  Iy();
  
  ClassDef(AliMUONSegmentationTrigger,1) //Segmentation class for trigger
    protected:
  //  Returns x-strip size for given module imodule
  Float_t StripSizeX(Int_t imodule);
  //  Returns y-strip size for given module imodule
  Float_t StripSizeY(Int_t imodule);    
 protected:
// Geometry Parameters

  Int_t fgNum[126];           // circuit Id. 
  Int_t fgNmodule;        // total number of modules
  Int_t fgNstripx[126];       // number of X strip / module
  Int_t fgNstripy[126];       // number of Y strip / module
  Float_t fgXcmin[126];       // x min position of modules
  Float_t fgXcmax[126];       // x max position of modules
  Float_t fgYcmin[126];       // y min position of modules
  Float_t fgYcmax[126];       // y max position of modules    
  Float_t    fZscale;            // scaling factor (Zx/Z1, x=1,2,3,4)

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













