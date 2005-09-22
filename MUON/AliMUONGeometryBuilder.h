/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONGeometryBuilder
/// \brief Manager class for geometry construction via geometry builders.
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_BUILDER_H
#define ALI_MUON_GEOMETRY_BUILDER_H

#include <TObject.h>
#include <TGeoMatrix.h>

class TObjArray;

class AliModule;
class AliMUONVGeometryBuilder;

class AliMUONGeometryBuilder : public TObject 
{
  public:
    AliMUONGeometryBuilder(AliModule* detector);
    AliMUONGeometryBuilder();
    virtual  ~AliMUONGeometryBuilder();
    
    // static methods
    static TGeoHMatrix Multiply(const TGeoMatrix& m1, const TGeoMatrix& m2); 
    static TGeoHMatrix Multiply(const TGeoMatrix& m1, const TGeoMatrix& m2,
                                const TGeoMatrix& m3); 
    static TGeoHMatrix Multiply(const TGeoMatrix& m1, const TGeoMatrix& m2,
                                const TGeoMatrix& m3, const TGeoMatrix& m4); 

    // methods
    void  AddBuilder(AliMUONVGeometryBuilder* geomBuilder);
    void  CreateGeometry();
    void  CreateMaterials();
    void  InitGeometry();
    void  WriteTransformations();
    void  WriteSVMaps(Bool_t rebuild = true);
    void  SetGlobalTransformation(const TGeoCombiTrans& transform);

    // Alignement
    virtual Bool_t  GetAlign() const;
    virtual void    SetAlign(Bool_t align);
 
  protected:
    AliMUONGeometryBuilder(const AliMUONGeometryBuilder& right);
    AliMUONGeometryBuilder&  operator = (const AliMUONGeometryBuilder& right);
 
  private:
    // method
    void PlaceVolume(const TString& name, const TString& mName, Int_t copyNo, 
             const TGeoHMatrix& matrix, Int_t npar, Double_t* param,
	     const char* only) const;
    void FillGlobalTransformations(AliMUONVGeometryBuilder* builder);	     

    // data members
    AliModule*      fModule;              // the AliRoot module
    Bool_t          fAlign;               // option to read transformations 
                                          // from a file
    TGeoCombiTrans  fGlobalTransformation;// global transformation 
                                          // applied to the whole geometry 
    TObjArray*      fGeometryBuilders;    // list of Geometry Builders

  ClassDef(AliMUONGeometryBuilder,4)  // Geometry builder
};

// inline functions

inline Bool_t  AliMUONGeometryBuilder::GetAlign() const
{ return fAlign; }

#endif //ALI_MUON_GEOMETRY_BUILDER_H







