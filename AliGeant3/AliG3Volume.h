#ifndef ALIG3VOLUME_H
#define ALIG3VOLUME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TGListTree.h>
#include "TROOT.h"

#include "THIGZ.h"
#include "TGeant3.h"
#include "TArrayF.h"
#include "TNamed.h"
#include "TList.h"

class TShape;
class TMaterial;

class AliG3Volume : public TNamed 
{
public:
    AliG3Volume() {;}
    AliG3Volume(const char* name);
    virtual ~AliG3Volume(){;}
    // G3 the volume
    virtual void    Draw(Option_t * option =0);
    // G3 volume specs
    virtual void    DrawSpec();
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
    virtual Int_t GetIdVolume() const {return fIdVolume;}
    // Get copy number
    virtual Int_t GetIdCopy() const {return fIdCopy;}
    // Get medium number
    virtual Int_t Medium() const  {return fIdMedium;}
    // Get material number
    virtual Int_t Material() const {return fIdMaterial;}
    // Increase copy number by one
    virtual void  AddCopy() {fIdCopy ++;}
    // Set link to ListTree Item 
    virtual void  SetItem(TObject *item) {fItem = item;}
    // Get link to ListTree Item
    virtual void  SetPosition(Float_t x, Float_t y, Float_t z);
    virtual TArrayF Position(Int_t i) const;
    
    virtual void  SetRotMatrix(Int_t irot) {fRotMatrix = irot;}
    virtual Int_t RotMatrix() const {return fRotMatrix;}
    virtual void  SetShape(Int_t shape) {fShape = shape;}
    virtual Int_t Shape() const {return fShape;}
    virtual void  SetParameters(Int_t np, Float_t* param);
    virtual Int_t NParam() const {return fNParam;} 
    virtual void  Parameters(Int_t i, TArrayF& param) const;
    virtual TList* Copies() const {return fCopies;}
    virtual void  AddCopy(AliG3Volume* volume);
    virtual AliG3Volume* Copy(Int_t i);
    
    virtual Int_t  NCopies() const {return fNCopies;}
    virtual Bool_t Posp() const {return fPosp;}
    virtual void   SetPosp(Bool_t flag) {fPosp = flag;}
    virtual void   CreateTShape(char* nameV, TMaterial* mat);
    virtual void   SetDivision(Int_t ndiv, Int_t axis, Float_t start, Float_t step);
    virtual void   Division(Int_t& ndiv, Int_t& axis, Float_t& start, Float_t& step) const;
    virtual Int_t   Axis()   {return fAxis;}
    virtual Int_t   Ndiv()   {return fNdiv;}
    virtual Float_t Step()   {return fStep;}
    virtual Float_t StartC() {return fStartC;}
    
	    
	    
    virtual TObject* GetItem() {return fItem;}

    AliG3Volume(const AliG3Volume&);
    

private:
    
    TArrayF  fPosition;     // position with respect to mother volume
    TArrayF  fParameters;   // volume parameters
    TList*   fCopies;       // volume copies
    Bool_t   fPosp;         // flag for G3 POSP
    Int_t    fNCopies;      // number of copies
    Int_t    fRotMatrix;    // rotation with respect to mother volume
    Int_t    fNParam;       // number of volume parameters
    Int_t    fAxis;         // division axis
    Int_t    fNdiv;         // number of divisions
    Float_t  fStep;         // number of steps
    Float_t  fStartC;       // start coordinate
    Int_t    fShape;       // G3 volume shape
    Float_t  fTheta;       // theta-angle for drawing
    Float_t  fPhi;         // phi-angle   for drawing
    Float_t  fPsi;         // psi-angle   for drawing 
    Float_t  fU;           // u-position
    Float_t  fV;           // v-position
    Float_t  fUscale;      // u-scaling factor
    Float_t  fVscale;      // v-scaling factor
    Bool_t   fHide;        // hide flag
    Bool_t   fShadow;      // shadow flag
    Int_t    fFill;        // fill option 1-6
    Int_t    fSeen;        // seen option -2 - 1
    Bool_t   fClip;        // clipping flag
    Float_t  fClipXmin;    // clip box range xmin
    Float_t  fClipXmax;    // clip box range xmax
    Float_t  fClipYmin;    // clip box range ymin
    Float_t  fClipYmax;    // clip box range ymax
    Float_t  fClipZmin;    // clip box range zmin
    Float_t  fClipZmax;    // clip box range zmax
    Int_t    fIdVolume;    // geant volume id
    Int_t    fIdMedium;    // geant medium id
    Int_t    fIdMaterial;  // geant material id    
    Int_t    fIdCopy;      // copy flag
    TObject* fItem;        //!current item
    AliG3Volume & operator=(const AliG3Volume&) {return *this;}

    ClassDef(AliG3Volume,1) // Volume Object for Drawing 
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
