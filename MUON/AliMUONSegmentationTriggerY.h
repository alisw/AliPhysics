#ifndef ALIMUONSEGMENTATIONTRIGGERY_H
#define ALIMUONSEGMENTATIONTRIGGERY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliMUONSegmentationTrigger.h"

class AliMUONChamber;
//----------------------------------------------
//
// Chamber segmentation virtual base class
//
class AliMUONSegmentationTriggerY :
public AliMUONSegmentationTrigger {
 public:
  AliMUONSegmentationTriggerY(){}
  virtual ~AliMUONSegmentationTriggerY(){}
  // Transform from pad to real coordinates
  virtual void    GetPadI(Float_t x,Float_t y,Int_t &ix,Int_t &iy);
  virtual void    GetPadI(Float_t x, Float_t y, Float_t z, Int_t &ix, Int_t &iy) 
      {GetPadI(x, y, ix, iy);}
  // Transform from real to pad coordinates
  virtual void    GetPadC(Int_t ix,Int_t iy,Float_t &x,Float_t &y);
  virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z) 
      {z=0; GetPadC(ix, iy, x , y);}
  // Pad size Dx*Dy 
  virtual void SetPadSize(Float_t dp1, Float_t dp2);
  // Strip size by Module
  virtual Float_t Dpx(Int_t imodule) const;
  virtual Float_t Dpy(Int_t imodule) const;
  // Set pad position
  virtual void     SetPad(Int_t ix, Int_t iy);
  // Set hit position
  virtual void     SetHit(Float_t xhit , Float_t yhit);
  virtual void     SetHit(Float_t xhit, Float_t yhit, Float_t zhit)
       {SetHit(xhit, yhit);}
  // Current integration parameters
  virtual void IntegrationLimits(Float_t& x1, Float_t& x2, Float_t& x3, Float_t& width);
  // Current Pad during Integration
  // x-coordinate
  virtual Int_t  Ix();
  // y-coordinate
  virtual Int_t  Iy();
  // Sector
  virtual Int_t ISector();
  // calculate sector from pad coordinates
  virtual Int_t Sector(Int_t ix, Int_t iy);

  
  // Get next neighbours 
  virtual void Neighbours
    (Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[2], Int_t Ylist[2]);
  //
  // Initialisation
  virtual void Init(Int_t chamber);    
  
  ClassDef(AliMUONSegmentationTriggerY,1) //Segmentation class for trigger X
    protected:    
  void  IntegrationParam(Float_t& x1, Float_t& x2, Float_t& y1);
  
// Geometry Parameters
  float fXofysmin[126][16]; // x-min   
  float fXofysmax[126][16]; // x-max
  float fYofysmin[126][16]; // y-min
  float fYofysmax[126][16]; // y-max

// Current pad during integration (cursor for disintegration)
  Int_t fIx;  // pad coord.  x 
  Int_t fIy;  // pad coord.  y 
  Float_t fX; // real coord. x
  Float_t fY; // real ccord. y

};
#endif



