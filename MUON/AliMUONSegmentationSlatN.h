#ifndef ALIMUONSEGMENTATIONSLATN_H
#define ALIMUONSEGMENTATIONSLATN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation classes for slat modules          //
//  to be used with AluMUONSegmentationSlat        //
/////////////////////////////////////////////////////

#include "AliMUONSegmentationSlat.h"

class TArrayF;
class TList;
class AliMUONSegmentationSlatModuleN;


class  AliMUONSegmentationSlatN :
public AliMUONSegmentationSlat {
 public:
    AliMUONSegmentationSlatN();
    virtual ~AliMUONSegmentationSlatN(){}
    //    
    // Set Chamber Segmentation Parameters
    //
    // Transform from pad to real coordinates and vice versa
    virtual void GetPadI(Float_t x, Float_t y , Float_t z, Int_t &ix, Int_t &iy);
    //
    // Get member data
    //
    // Pad size in x by Sector
    virtual Float_t Dpx(Int_t isec) const;
    // Pad size in y by Sector
    virtual Float_t Dpy(Int_t isec) const;
    //
    //
    // Class specific methods
    virtual void GlobalToLocal(
	Int_t ix, Int_t iy, Int_t &islat, Int_t &ixlocal, Int_t &iylocal);
    virtual void LocalToGlobal(
	 Int_t islat, Int_t ixlocal, Int_t iylocal, Int_t &ix, Int_t &iy);
    // Factory method for associated slat module class
     AliMUONSegmentationSlatModule* CreateSlatModule();
 private:
    ClassDef(AliMUONSegmentationSlatN,1) // Segmentation for Muon Chamber built from Slat Modules
	
};
	

#endif






