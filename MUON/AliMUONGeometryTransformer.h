/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONGeometryTransformer
/// \brief Top container class for geometry transformations
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_TRANSFORMER_H
#define ALI_MUON_GEOMETRY_TRANSFORMER_H

#include <TObject.h>
#include <TGeoMatrix.h>

class TObjArray;

class AliMUONGeometryModuleTransformer;

class AliMUONGeometryTransformer : public TObject
{
  public:
    AliMUONGeometryTransformer(Bool_t isOwner);
    AliMUONGeometryTransformer();
    virtual  ~AliMUONGeometryTransformer();
    
    // methods
    void  AddModuleTransformer(AliMUONGeometryModuleTransformer* transformer);

    // IO
    Bool_t  ReadTransformations(const TString& fileName);
    Bool_t  WriteTransformations(const TString& fileName) const;

    // transformation methods 
    void Global2Local(Int_t detElemId,
                 Float_t xg, Float_t yg, Float_t zg, 
                 Float_t& xl, Float_t& yl, Float_t& zl) const;
    void Global2Local(Int_t detElemId,
                 Double_t xg, Double_t yg, Double_t zg, 
                 Double_t& xl, Double_t& yl, Double_t& zl) const;

    void Local2Global(Int_t detElemId,
                 Float_t xl, Float_t yl, Float_t zl, 
                 Float_t& xg, Float_t& yg, Float_t& zg) const;
    void Local2Global(Int_t detElemId,
                 Double_t xl, Double_t yl, Double_t zl, 
                 Double_t& xg, Double_t& yg, Double_t& zg) const;

    // get methods
    const AliMUONGeometryModuleTransformer* GetModuleTransformer(
                               Int_t index, Bool_t warn = true) const;

    const AliMUONGeometryModuleTransformer* GetModuleTransformerByDEId(
                               Int_t detElemId, Bool_t warn = true) const;

  protected:
    AliMUONGeometryTransformer(const AliMUONGeometryTransformer& right);
    AliMUONGeometryTransformer&  operator = (const AliMUONGeometryTransformer& right);
 
  private:
    // methods
    AliMUONGeometryModuleTransformer* GetModuleTransformerNonConst(
                                    Int_t index, Bool_t warn = true) const;
    TString  ComposePath(const TString& volName, Int_t copyNo) const; 

    TGeoHMatrix GetTransform(
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3, 
 		  Double_t a4, Double_t a5, Double_t a6) const;
    void FillData(Int_t moduleId, 
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3, 
 		  Double_t a4, Double_t a5, Double_t a6); 
    void FillData(Int_t id, const TString& volName, Int_t copyNo,
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3, 
 		  Double_t a4, Double_t a5, Double_t a6);

    TString ReadData1(ifstream& in);
    TString ReadData2(ifstream& in);

    void WriteTransform(ofstream& out, const TGeoCombiTrans* transform) const;
    void WriteData1(ofstream& out) const;
    void WriteData2(ofstream& out) const;

    // data members
    TObjArray*  fModuleTransformers; // list of module transformers

  ClassDef(AliMUONGeometryTransformer,1)  // Geometry parametrisation
};

#endif //ALI_MUON_GEOMETRY_TRANSFORMER_H







