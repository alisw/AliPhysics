#ifndef ALIDRAWVOLUME_H
#define ALIDRAWVOLUME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TGListTree.h>
#include "TROOT.h"

#include "THIGZ.h"
#include "TGeant3.h"

class AliDrawVolume : public TObject 
{
public:
    AliDrawVolume(char* name);
    virtual ~AliDrawVolume(){;}
    // Draw the volume
    virtual void    Draw(Option_t * option =0);
    // Draw volume specs
    virtual void    DrawSpec();
    // Return volume name
    virtual char*   Name();
    // Set volume parameter i
    virtual void    SetParam(Int_t i, Float_t param);
    // Get volume parameters i
    virtual Float_t GetParam(Int_t i);
    // Set volume id
    virtual void  SetIdVolume(Int_t id) {fIdVolume = id;}
    // Set volume copy number
    virtual void  SetIdCopy(Int_t id)   {fIdCopy = id;}
    // Set volume medium number
    virtual void  SetIdMedium(Int_t id)   {fIdMedium = id;}
    // Set volume material number
    virtual void  SetIdMaterial(Int_t id) {fIdMaterial = id;}
    // Get volume id
    virtual Int_t GetIdVolume()         {return fIdVolume;}
    // Get copy number
    virtual Int_t GetIdCopy()           {return fIdCopy;}
    // Get medium number
    virtual Int_t Medium()   {return fIdMedium;}
    // Get material number
    virtual Int_t Material() {return fIdMaterial;}
    // Increase copy number by one
    virtual void  AddCopy()             {fIdCopy ++;}
    // Set link to ListTree Item 
    virtual void  SetItem(TGListTreeItem *item) {fItem = item;}
    // Get link to ListTree Item
    virtual TGListTreeItem* GetItem() {return fItem;}
	    
private:
    char*   fName;        // name of the volume 
    Float_t fTheta;       // theta-angle for drawing
    Float_t fPhi;         // phi-angle   for drawing
    Float_t fPsi;         // psi-angle   for drawing 
    Float_t fU;           // u-position
    Float_t fV;           // v-position
    Float_t fUscale;      // u-scaling factor
    Float_t fVscale;      // v-scaling factor
    Bool_t  fHide;        // hide flag
    Bool_t  fShadow;      // shadow flag
    Int_t   fFill;        // fill option 1-6
    Int_t   fSeen;        // seen option -2 - 1
    Bool_t  fClip;        // clipping flag
    Float_t fClipXmin;    // clip box range xmin
    Float_t fClipXmax;    // clip box range xmax
    Float_t fClipYmin;    // clip box range ymin
    Float_t fClipYmax;    // clip box range ymax
    Float_t fClipZmin;    // clip box range zmin
    Float_t fClipZmax;    // clip box range zmax
    Int_t   fIdVolume;    // geant volume id
    Int_t   fIdMedium;    // geant medium id
    Int_t   fIdMaterial;  // geant material id    
    Int_t   fIdCopy;      // copy flag
    TGListTreeItem        *fItem; // current item

  AliDrawVolume(const AliDrawVolume&) {}
  AliDrawVolume & operator=(const AliDrawVolume&) {return *this;}

    ClassDef(AliDrawVolume,1) // Volume Object for Drawing 
};

//
// Drawing parameter tags
enum AliDrawParamId {
   kTheta,
   kPhi,
   kPsi,
   kU,
   kV,
   kUscale,
   kVscale,
   kShadow,
   kHide,
   kFill,
   kSeen,
   kClip,
   kClipXmin,
   kClipXmax,
   kClipYmin,
   kClipYmax,
   kClipZmin,
   kClipZmax
};


#endif
