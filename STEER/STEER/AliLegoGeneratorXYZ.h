#ifndef ALILEGOGENERATORXYZ_H
#define ALILEGOGENERATORXYZ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//    Utility class to compute and draw Radiation Length Map                 //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliLegoGenerator.h"

class AliLegoGeneratorXYZ : public AliLegoGenerator
{

 public:
    AliLegoGeneratorXYZ();
    AliLegoGeneratorXYZ(char* axis);
    AliLegoGeneratorXYZ(Int_t nc1, Float_t c1min, Float_t c1max,
			Int_t nc2, Float_t c2min, Float_t c2max,
			Float_t rmin, Float_t rmax, Float_t zmax);
    virtual ~AliLegoGeneratorXYZ() {}
    virtual void    Generate();
 protected:
    Float_t fDir1[3];  // 1st unit vector spanning the scanning plane
    Float_t fDir2[3];  // 2nd unit vector spanning the scanning plane
    Float_t fDir3[3];  // Direction of flight for geantinos
    
    ClassDef(AliLegoGeneratorXYZ,1) //Lego GeneratorXYZ
};

#endif








