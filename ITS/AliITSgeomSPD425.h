#ifndef ALIITSGEOMSPD425_H
#define ALIITSGEOMSPD425_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */



#include "AliITSgeomSPD.h"
#include "TBRIK.h"


// temporary - this will migrate into the segmentation class


class AliITSgeomSPD425 : public AliITSgeomSPD {

 public:
    AliITSgeomSPD425();
    AliITSgeomSPD425(Float_t cy,Int_t Nx,Float_t *bx,Int_t Nz,Float_t *bz);
    AliITSgeomSPD425(AliITSgeomSPD425 &source);
    virtual ~AliITSgeomSPD425();
    AliITSgeomSPD425& operator=(AliITSgeomSPD425 &source);
    TShape *GetShape() {return fShapeSPD;} // returns shape
    Float_t GetDx() {return fdx;} // returns half detector width cm
    Float_t GetDy() {return fdy;} // returns half detector thickness cm
    Float_t GetDz() {return fdz;} // returns half detector lenght cm
    Int_t   GetNbinsX() {return fNbinx;} // returns number of x pixels
    Float_t GetBinSizeX(Int_t i) {return fBinSizeX[i];} // returns x pixel size
    Int_t   GetNbinsZ() {return fNbinz;} // returns number of z pixels
    Float_t GetBinSizeZ(Int_t i) {return fBinSizeZ[i];} // returns z pixel size
    void    ReSetBins(Float_t dy,Int_t Nx, Float_t *bx,Int_t Nz,Float_t *bz);

 private:
    // define shape of active area using ROOT shapes so that they can
    // be easly plotted. Inputs to TBRIK are
    // Shape name (what ever that is)
    // Shape title (what ever that means)
    // name of material (something I took from ITSgeometry.tme file
    // fdx => 1/2 thickness of wafer's active volume (cm)
    // fdy => 1/2 r*phi size of active volume (cm)
    // fdz => 1/2 size of active volume (cm)
    TBRIK *fShapeSPD;     // Shape of sensitive volume.
    Float_t fdx;          // Brick half width cm
    Float_t fdy;          // Brick half thickness cm
    Float_t fdz;          // Brick half length cm
    Float_t *fBinSizeX;  // Pixel size in X, cm
    Int_t   fNbinx;      // Number of pixels in x
    Float_t *fBinSizeZ;  // Pixel size in Z, cm
    Int_t   fNbinz;      // Number of pixels in z

    ClassDef(AliITSgeomSPD425,2) // ITS SPD detector geometry class for 425X50 micron pixel size.

};
#endif
