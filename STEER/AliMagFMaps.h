#ifndef ALIMAGFMAPS_H
#define ALIMAGFMAPS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Magnetic field composed by 3 maps: the L3 magnet, extended region, and
// dipole magnet
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
    AliMagFMaps(const char *name, const char *title, Int_t integ,
		Float_t factor, Float_t fmax, Int_t map = k2kG,
		Int_t l3 = 1);
    AliMagFMaps(const AliMagFMaps &mag);
    virtual ~AliMagFMaps();
    virtual void    Field(const float  *x, float  *b) const;
    virtual void    Field(const double *x, double *b) const;
    AliFieldMap* FieldMap(Int_t i) {return fFieldMap[i];}
    virtual void ReadField();
    virtual Float_t SolenoidField() const;
    virtual void    SetL3ConstField(Int_t flag = 0) {fL3Option = flag;}
    virtual void    SetL3ConstField(Float_t bsol, Int_t flag = 0)
	{fL3Option = flag; fSolenoidUser = bsol;}
    
    AliMagFMaps & operator=(const AliMagFMaps &magf)
      {magf.Copy(*this); return *this;}

protected:
    void Copy(TObject &magf) const;

    AliFieldMap* fFieldMap[3];     // Field maps
    Float_t      fSolenoid;        // Solenoid field setting
    Float_t      fSolenoidUser;    // User set solenoid field setting  
    Int_t        fL3Option;        // Option for field inside L3
    ClassDef(AliMagFMaps,4)        // Class for all Alice MagField using three Maps with Constant Mesh
};

#endif
