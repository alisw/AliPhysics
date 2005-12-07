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

#include "AliMUONGeometry.h"

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
    //
    void  AddBuilder(AliMUONVGeometryBuilder* geomBuilder);
    void  CreateGeometry();
    void  CreateMaterials();

    void  InitGeometry();
    void  InitGeometry(const TString& svmapFileName);

    void  ReadTransformations();
    void  ReadTransformations(const TString& fileName);

    void  WriteTransformations();
    void  WriteTransformations(const TString& fileName);

    void  WriteSVMaps();
    void  WriteSVMaps(const TString& fileName, Bool_t rebuild = true);
    
    // Geometry parametrisation
    const AliMUONGeometry*            GetGeometry() const;
    const AliMUONGeometryTransformer* GetTransformer() const;

    // Alignement
    virtual Bool_t  GetAlign() const;
    virtual void    SetAlign(Bool_t align = true);
    virtual void    SetAlign(const TString& fileName, Bool_t align = true);
 
  protected:
    AliMUONGeometryBuilder(const AliMUONGeometryBuilder& right);
    AliMUONGeometryBuilder&  operator = (const AliMUONGeometryBuilder& right);
 
  private:
    // method
    void PlaceVolume(const TString& name, const TString& mName, Int_t copyNo, 
             const TGeoHMatrix& matrix, Int_t npar, Double_t* param,
	     const char* only) const;
    void SetAlign(AliMUONVGeometryBuilder* builder);	     

    // static data members
    static const TString  fgkDefaultTransformFileName; // default transformations file name					   
    static const TString  fgkDefaultSVMapFileName;     // default svmaps file name					   
    static const TString  fgkOutFileNameExtension;     // default output file name extension					   

    // data members
    AliModule*       fModule;              // the AliRoot module
    Bool_t           fAlign;               // option to read transformations 
                                           // from a file
    TString          fTransformFileName;   // transformations file name					   
    TString          fSVMapFileName;       // svmaps file name					   
    TGeoCombiTrans   fGlobalTransformation;// global transformation 
                                           // applied to the whole geometry 
    TObjArray*       fGeometryBuilders;    // list of Geometry Builders
    AliMUONGeometry* fGeometry;            // geometry parametrisation

  ClassDef(AliMUONGeometryBuilder,4)  // Geometry builder
};

// inline functions

inline void  AliMUONGeometryBuilder::InitGeometry()
{ InitGeometry(fSVMapFileName); }

inline void  AliMUONGeometryBuilder::ReadTransformations()
{ ReadTransformations(fTransformFileName); }

inline void  AliMUONGeometryBuilder::WriteTransformations()
{ WriteTransformations(fTransformFileName + fgkOutFileNameExtension); }

inline void  AliMUONGeometryBuilder::WriteSVMaps()
{ WriteSVMaps(fSVMapFileName + fgkOutFileNameExtension); }

inline 
const AliMUONGeometry* AliMUONGeometryBuilder::GetGeometry() const
{ return fGeometry; }

inline 
const AliMUONGeometryTransformer* AliMUONGeometryBuilder::GetTransformer() const
{ return fGeometry->GetTransformer(); }

inline Bool_t  AliMUONGeometryBuilder::GetAlign() const
{ return fAlign; }

#endif //ALI_MUON_GEOMETRY_BUILDER_H







