#ifndef ALIG3MATERIAL_H
#define ALIG3MATERIAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TMaterial.h"

class AliG3Material : public TMaterial 
{
public:
    AliG3Material(){}
    AliG3Material(char* name, char* title,
		   Float_t a, Float_t z, Float_t dens, Float_t radl, Float_t intl);
    
    virtual ~AliG3Material(){}
    // Dump material parameters
    virtual void  Dump();
    // Get material id
    virtual Int_t Id()    {return fId;}
    virtual void  SetId(Int_t id) {fId = id;}
    
private:
    Int_t   fId;          // Id number of the material
    AliG3Material(const AliG3Material &) {}
    AliG3Material &operator=(const AliG3Material &) {return *this;}

    ClassDef(AliG3Material,1) // Material Object for GUI 
};

#endif








