#ifndef ALIITSGEOMSDD_H
#define ALIITSGEOMSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include "TBRIK.h"
#include "AliITSgeom.h"

// temporary 

class AliITSgeomSDD: public TObject {
 public:
    AliITSgeomSDD();
    AliITSgeomSDD(AliITSgeomSDD &source);
    virtual ~AliITSgeomSDD(){};
    AliITSgeomSDD& operator=(AliITSgeomSDD &source);
    TBRIK *GetShape() const {return fShapeSDD;}
    Float_t GetDx() {return fDx;}
    Float_t GetDy() {return fDy;}
    Float_t GetDz() {return fDz;}
    // or what other or different information that is needed.

 private:
    // define shape of active area using ROOT shapes so that they can
    // be easly plotted. Inputs to TBRIK are
    // Shape name (what ever that is)
    // Shape title (what ever that means)
    // name of material (something I took from ITSgeometry.tme file
    // dx => 1/2 thickness of wafer's active volume (cm)
    // dy => 1/2 r*phi size of active volume (cm)
    // dz => 1/2 size of active volume (cm)
    Float_t fDx;          // Brick half width cm
    Float_t fDy;          // Brick half thickness cm
    Float_t fDz;          // Brick half length cm
    TBRIK *fShapeSDD;     // shape of sensitive volume

    ClassDef(AliITSgeomSDD,1) // ITS SDD detector geometry class
};
#endif
