#ifndef ALIITSGEOMSSD_H
#define ALIITSGEOMSSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include "TBRIK.h"
//#include "AliITSgeom.h"


class AliITSgeom;

class AliITSgeomSSD: public TObject {

 public:
    AliITSgeomSSD();
    virtual ~AliITSgeomSSD() {
      // destructor
    }
    AliITSgeomSSD(const AliITSgeomSSD &source); 
    AliITSgeomSSD& operator=(const AliITSgeomSSD &source); 
    
    TBRIK *GetShape() const {
      // get shape
      return fShapeSSD;
    }
    Float_t GetDx(){
      // get Dx
      return fShapeSSD->GetDx();
    }
    Float_t GetDy(){
      // get Dx
      return fShapeSSD->GetDy();
    }
    Float_t GetDz(){
      // get Dx
      return fShapeSSD->GetDz();
    }
    
 private:
    // define shape of active area using ROOT shapes so that they can
    // be easly plotted. Inputs to TBRIK are
    // Shape name (what ever that is)
    // Shape title (what ever that means)
    // name of material (something I took from ITSgeometry.tme file
    // dx => 1/2 thickness of wafer's active volume (cm)
    // dy => 1/2 r*phi size of active volume (cm)
    // dz => 1/2 size of active volume (cm)
    TBRIK *fShapeSSD; // comment
    // Other infomation like.
    // Float_t fTopPitch;      // cm
    // Float_t fTopWidth;      // cm
    // Float_t fTopLength;     // cm
    // Float_t fTopAngle;      // cm
    // Float_t fBottomPitch;   // cm
    // Float_t fBottomWidth;   // cm
    // Float_t fBottomLength;  // cm
    // Float_t fBottomAngle;   // cm
    // or what other or different information that is needed.
    
    ClassDef(AliITSgeomSSD,1) // ITS SSD detector geometry class
      };

#endif
      











