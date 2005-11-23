/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup geometry
/// \class AliMUONGeometryModuleTransformer
/// \brief Geometry transformationer for detector module
///
/// Class for definition of the detector module parameters
/// (the transformations of detection elements, mapping between
///  sensitive volumes and detection elements).
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_MODULE_TRANSFORMS_H
#define ALI_MUON_GEOMETRY_MODULE_TRANSFORMS_H

#include <TObject.h>
#include <TString.h>

class TGeoTranslation;
class TGeoRotation;
class TGeoCombiTrans;
class TObjArray;
class TArrayI;

class AliMUONGeometryDetElement;
class AliMUONGeometryStore;

class AliMUONGeometryModuleTransformer : public TObject
{
  public:
    AliMUONGeometryModuleTransformer(Int_t moduleId);
    AliMUONGeometryModuleTransformer();
    virtual ~AliMUONGeometryModuleTransformer();

    // methods
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

    // set methods
    void  SetTransformation(const TGeoCombiTrans& transform);
 
    // get methods
    Int_t                  GetModuleId() const;
    const TGeoCombiTrans*  GetTransformation() const;    

    AliMUONGeometryStore*       GetDetElementStore() const;
    AliMUONGeometryDetElement*  GetDetElement(
                                   Int_t detElemId, Bool_t warn = true) const;    

  protected:
    AliMUONGeometryModuleTransformer(const AliMUONGeometryModuleTransformer& rhs);
    // operators  
    AliMUONGeometryModuleTransformer& operator = (const AliMUONGeometryModuleTransformer& rhs);

  private:
    // data members
    Int_t                 fModuleId;      // the module Id
    TGeoCombiTrans*       fTransformation;// the module transformation wrt to top
                                          // volume
    AliMUONGeometryStore* fDetElements;   // detection elements
 
  ClassDef(AliMUONGeometryModuleTransformer,2) // MUON geometry module class
};

// inline functions

inline Int_t  AliMUONGeometryModuleTransformer::GetModuleId() const
{ return fModuleId; }

inline const TGeoCombiTrans* AliMUONGeometryModuleTransformer::GetTransformation() const 
{ return fTransformation; }

inline  AliMUONGeometryStore* AliMUONGeometryModuleTransformer::GetDetElementStore() const
{ return fDetElements; }

#endif //ALI_MUON_GEOMETRY_MODULE_PARAM_H
