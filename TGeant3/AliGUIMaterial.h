#ifndef ALIGUIMATERIAL_H
#define ALIGUIMATERIAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"

class AliGUIMaterial : public TObject 
{
public:
    AliGUIMaterial();
    AliGUIMaterial(Int_t imat, char* name, Float_t a, Float_t z,
		   Float_t dens, Float_t radl, Float_t absl);
    virtual ~AliGUIMaterial(){}
    // Dump material parameters
    virtual void  Dump();
    // Get material id
    virtual Int_t Id();
    // Get material name
    virtual char* Name();
    // Get mass number 
    virtual Float_t A();
    // Get charge number 
    virtual Float_t Z();
    // Get density
    virtual Float_t Density();
    // Get radiation length
    virtual Float_t RadiationLength();
    // Get absorption lenth
    virtual Float_t AbsorptionLength();
    // Plot
    virtual void  Plot();
    // Set and get link to widget entry
    virtual Int_t ItemId() {return fItem;}
    virtual void  SetItemId(Int_t id) {fItem=id;}
private:
    Int_t   fId;          // Id number of the material
    char*   fName;        // name of the material 
    Float_t fA;           // mass number of the material
    Float_t fZ;           // charge number of the material
    Float_t fDensity;     // density of the material
    Float_t fRadl;        // radiation length of the material
    Float_t fAbsl;        // absorption length
    //
    Int_t   fItem;            // Link to Widget Entry

  AliGUIMaterial(const AliGUIMaterial &) {}
  AliGUIMaterial &operator=(const AliGUIMaterial &) {return *this;}

    ClassDef(AliGUIMaterial,1) // Material Object for GUI 
};

#endif








