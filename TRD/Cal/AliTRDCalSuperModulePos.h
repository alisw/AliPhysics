#ifndef AliTRDCALSUPERMODULEPOS_H
#define AliTRDCALSUPERMODULEPOS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for position parameters of the supermodules        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTRDCalSuperModulePos : public TNamed {
  public:
    enum { kNdet = 540, kNstacks = 90, kNcham = 5, kNsect = 18 };
  
    AliTRDCalSuperModulePos();
    AliTRDCalSuperModulePos(const Text_t* name, const Text_t* title);
  
    const Float_t* GetPos(Int_t sm) const { return fPos[sm]; };
    const Float_t* GetRot(Int_t sm) const { return fRot[sm]; };

    inline void SetPos(Int_t sm, Float_t x, Float_t y, Float_t z);
    void SetPos(Int_t sm, Float_t* xyz) { SetPos(sm, xyz[0], xyz[1], xyz[2]); };
  
    inline void SetRot(Int_t sm, Float_t x, Float_t y, Float_t z);
    void SetRot(Int_t sm, Float_t* xyz) { SetRot(sm, xyz[0], xyz[1], xyz[2]); };
  
  protected:
    Float_t fPos[kNsect][3];                    //  Deviations of the positions of the supermodules from the ideal position
    Float_t fRot[kNsect][3];                    //  Rotation of the supermodules in respect to the ideal position
    
    ClassDef(AliTRDCalSuperModulePos,1)                     
};
    
void AliTRDCalSuperModulePos::SetPos(Int_t sm, Float_t x, Float_t y, Float_t z) 
{ 
  fPos[sm][0] = x; 
  fPos[sm][1] = y;
  fPos[sm][2] = z; 
}

void AliTRDCalSuperModulePos::SetRot(Int_t sm, Float_t x, Float_t y, Float_t z) 
{ 
  fRot[sm][0] = x;
  fRot[sm][1] = y;
  fRot[sm][2] = z;
}

#endif
