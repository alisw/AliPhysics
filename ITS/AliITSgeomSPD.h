#ifndef ALIITSGEOMSPD_H
#define ALIITSGEOMSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "TShape.h"
#include "TBRIK.h"

#include "TObject.h"


class AliITSgeomSPD:public TObject {

 public:
    virtual ~AliITSgeomSPD(){}; // destructor
    virtual TShape *GetShape(){return (TShape*)0;} // get TShape
    virtual Float_t GetDx(){return 0.0;}; // get width
    virtual Float_t GetDy(){return 0.0;}; // get thickness
    virtual Float_t GetDz(){return 0.0;}; // get length
    virtual Int_t   GetNbinsX(){return 0;}; // get number of x pixels
    virtual Float_t GetBinSizeX(Int_t i){return 0.0;}; // get x pixel size
    virtual Int_t   GetNbinsZ(){return 0;}; // get number of z pixels
    virtual Float_t GetBinSizeZ(Int_t i){return 0.0;}; // get z pixel size
    virtual void    ReSetBins(Float_t dy,Int_t Nx, Float_t *bx,
			      Int_t Nz,Float_t *bz){;}; // change pixel sizes

    ClassDef(AliITSgeomSPD,2) // ITS SPD detector geometry class

};
#endif
