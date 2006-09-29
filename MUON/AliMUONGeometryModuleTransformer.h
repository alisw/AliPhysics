/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup geometry
/// \class AliMUONGeometryModuleTransformer
/// \brief Geometry transformer for a detector module
///
/// Class for definition of the trasformation for adetector module
/// and its detection elements
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_MODULE_TRANSFORMER_H
#define ALI_MUON_GEOMETRY_MODULE_TRANSFORMER_H

#include <TObject.h>
#include <TString.h>

class AliMUONGeometryDetElement;
class AliMUONGeometryStore;

class TGeoTranslation;
class TGeoRotation;
class TGeoHMatrix;
class TObjArray;
class TArrayI;

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
    void  SetTransformation(const TGeoHMatrix& transform);
    void  SetVolumePath(const TString& volumePath);
 
    // get methods
    Int_t    GetModuleId() const;
    TString  GetModuleName() const;
    TString  GetVolumePath() const;
    TString  GetVolumeName() const;
    TString  GetMotherVolumeName() const;

    const TGeoHMatrix*  GetTransformation() const;    

    AliMUONGeometryStore*       GetDetElementStore() const;
    AliMUONGeometryDetElement*  GetDetElement(
                                   Int_t detElemId, Bool_t warn = true) const;    

  protected:
    AliMUONGeometryModuleTransformer(const AliMUONGeometryModuleTransformer& rhs);
    AliMUONGeometryModuleTransformer& 
      operator = (const AliMUONGeometryModuleTransformer& rhs);

  private:
    // static data members
    static const TString  fgkModuleNamePrefix; /// < Geometry module name prefix

    // data members
    Int_t                 fModuleId;   ///< the module Id
    TString               fModuleName; ///< the module name
    TString               fVolumePath; ///< \brief the full path of aligned module volume
                                       /// or envelope in geometry
    TGeoHMatrix*          fTransformation;///< \brief the module transformation wrt to top
                                          /// volume (world)
    AliMUONGeometryStore* fDetElements;   ///< detection elements
 
  ClassDef(AliMUONGeometryModuleTransformer,3) // MUON geometry module class
};

// inline functions

/// Set the full path of aligned module volume or envelope in geometry
inline void 
AliMUONGeometryModuleTransformer::SetVolumePath(const TString& volumePath)
{ fVolumePath = volumePath; }

/// Return module ID
inline Int_t  
AliMUONGeometryModuleTransformer::GetModuleId() const
{ return fModuleId; }

/// Return module name
inline TString
AliMUONGeometryModuleTransformer::GetModuleName() const
{ return fModuleName; }

/// Return the full path of aligned module volume or envelope in geometry
inline TString 
AliMUONGeometryModuleTransformer::GetVolumePath() const
{ return fVolumePath; }

/// Return the module transformation wrt to the top volume (world)
inline const TGeoHMatrix* 
AliMUONGeometryModuleTransformer::GetTransformation() const 
{ return fTransformation; }

/// Return detection elements associated with this module
inline  AliMUONGeometryStore* 
AliMUONGeometryModuleTransformer::GetDetElementStore() const
{ return fDetElements; }

#endif //ALI_MUON_GEOMETRY_MODULE_TRANSFORMER_H
