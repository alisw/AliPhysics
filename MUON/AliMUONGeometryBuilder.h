#ifndef ALI_MUON_GEOMETRY_BUILDER_H
#define ALI_MUON_GEOMETRY_BUILDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
// Class AliMUONGeometryBuilder
// ----------------------------
// MUON manager class for geometry construction,
// separated form AliMUONv1
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <TObject.h>

class TGeoCombiTrans;
class TGeoHMatrix;
class TObjArray;

class AliMUON;
class AliMUONVGeometryBuilder;

class AliMUONGeometryBuilder : public TObject 
{
  public:
    AliMUONGeometryBuilder();
    AliMUONGeometryBuilder(AliMUON* muon);
    virtual  ~AliMUONGeometryBuilder();

    void  CreateGeometry();
    void  CreateMaterials();
    void  InitGeometry();
    void  WriteTransformations();
    void  WriteSVMaps(Bool_t rebuild = true);

    // Alignement
    virtual Bool_t  GetAlign() const;
    virtual void    SetAlign(Bool_t align);
 
    void  AddBuilder(AliMUONVGeometryBuilder* geomBuilder);
   
  protected:
    AliMUONGeometryBuilder(const AliMUONGeometryBuilder& right);
    AliMUONGeometryBuilder&  operator = (const AliMUONGeometryBuilder& right);
 
  private:
    // method
    void PlaceVolume(const TString& name, const TString& mName, Int_t copyNo, 
             const TGeoHMatrix& matrix, Int_t npar, Double_t* param,
	     const char* only) const;

    // data members
    AliMUON*        fMUON;                // MUON detector
    Bool_t          fAlign;               // option to read transformations 
                                          // from a file
    TGeoCombiTrans* fGlobalTransformation;// global transformation 
                                          // applied to the whole geometry 
    TObjArray*      fGeometryBuilders;    // list of Geometry Builders

  ClassDef(AliMUONGeometryBuilder,2)  // MUON Detector class Version 1
};

// inline functions

inline Bool_t  AliMUONGeometryBuilder::GetAlign() const
{ return fAlign; }

#endif //ALI_MUON_GEOMETRY_BUILDER_H







