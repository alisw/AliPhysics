#ifndef ITSgeomSPD_H
#define ITSgeomSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TShape.h"
#include "TBRIK.h"

class AliITSgeomSPD: public TObject {
 private:
    // define shape of active area using ROOT shapes so that they can
    // be easly plotted. Inputs to TBRIK are
    // Shape name (what ever that is)
    // Shape title (what ever that means)
    // name of material (something I took from ITSgeometry.tme file
    // dx => 1/2 thickness of wafer's active volume (cm)
    // dy => 1/2 r*phi size of active volume (cm)
    // dz => 1/2 size of active volume (cm)
    TBRIK *fShapeSPD;
    // Other infomation like.
    // Float_t fPitchZ;     // cm
    // Float_t fPitchY;     // cm
    // Float_t fCellZ;      // cm
    // Float_t fCellY;      // cm
    // or what other or different information that is needed.
 public:
    AliITSgeomSPD();
    virtual ~AliITSgeomSPD(){};
    TBRIK *GetShape() const {return fShapeSPD;}

    ClassDef(AliITSgeomSPD,1)	
};
#endif
