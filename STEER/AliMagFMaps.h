#ifndef ALIMAGFMAPS_H
#define ALIMAGFMAPS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Author: Andreas Morsch <andreas.morsch@cern.ch>
//

#include "AliMagFC.h"
class AliFieldMap;

class AliMagFMaps : public AliMagFC
{
  //Alice Magnetic Field with constant mesh

public:
    enum constants {k2kG, k4kG, k5kG};
    AliMagFMaps();
    AliMagFMaps(const char *name, const char *title, const Int_t integ,
		const Float_t factor, const Float_t fmax, const Int_t map = k2kG,
		const Int_t l3 = 1);
    AliMagFMaps(const AliMagFMaps &mag);
    virtual ~AliMagFMaps();
    virtual void    Field(Float_t *x, Float_t *b);
    AliFieldMap* FieldMap(Int_t i) {return fFieldMap[i];}
    virtual void ReadField();
    virtual Float_t SolenoidField() const;
    virtual void    SetL3ConstField(Int_t flag = 0) {fL3Option = flag;}
    virtual void    SetL3ConstField(Float_t bsol, Int_t flag = 0)
	{fL3Option = flag; fSolenoidUser = bsol;}
    
    virtual AliMagFMaps & operator=(const AliMagFMaps &magf)
      {magf.Copy(*this); return *this;}

protected:
    void Copy(TObject &magf) const;

    AliFieldMap* fFieldMap[3];     // Field maps
    Float_t      fSolenoid;        // Solenoid field setting
    Float_t      fSolenoidUser;    // User set solenoid field setting  
    Int_t        fL3Option;        // Option for field inside L3
    Int_t        fFieldRead;       // Field has been read in
    ClassDef(AliMagFMaps,3)        // Class for all Alice MagField using three Maps with Constant Mesh
};

#endif
