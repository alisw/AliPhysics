/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONGeometryDetElement
/// \brief Class for storing detection element transformations 
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_DET_ELEMENT_H
#define ALI_MUON_GEOMETRY_DET_ELEMENT_H

#include <TObject.h>
#include <TString.h>

class TGeoHMatrix;

class AliMUONGeometryDetElement : public TObject
{
  public:
    AliMUONGeometryDetElement(Int_t detElemId,
                              const TString& volumePath);
    AliMUONGeometryDetElement();
    virtual ~AliMUONGeometryDetElement();

    // static methods
    static TString GetDENamePrefix();

    // methods
    void Global2Local(
                 Float_t xg, Float_t yg, Float_t zg, 
                 Float_t& xl, Float_t& yl, Float_t& zl) const;
    void Global2Local(
                 Double_t xg, Double_t yg, Double_t zg, 
                 Double_t& xl, Double_t& yl, Double_t& zl) const;

    void Local2Global(
                 Float_t xl, Float_t yl, Float_t zl, 
                 Float_t& xg, Float_t& yg, Float_t& zg) const;
    void Local2Global(
                 Double_t xl, Double_t yl, Double_t zl, 
                 Double_t& xg, Double_t& yg, Double_t& zg) const;
    void PrintLocalTransform() const;
    void PrintGlobalTransform() const;

    // set methods
    void SetLocalTransformation(const TGeoHMatrix& transform);
    void SetGlobalTransformation(const TGeoHMatrix& transform);
    void SetVolumePath(const TString& volumePath);
    
    // get methods
    Int_t    GetId() const;
    TString  GetDEName() const;
    TString  GetVolumePath() const;
    TString  GetVolumeName() const;
    Int_t    GetVolumeCopyNo() const;
    const TGeoHMatrix*  GetLocalTransformation() const;
    const TGeoHMatrix*  GetGlobalTransformation() const;

  protected:
    AliMUONGeometryDetElement(const AliMUONGeometryDetElement& rhs);
    AliMUONGeometryDetElement& operator = (const AliMUONGeometryDetElement& rhs);
  
  private:
    // methods
    void PrintTransform(const TGeoHMatrix* transform) const;
 
     // static data members
    static const TString  fgkDENamePrefix; /// < Geometry module name prefix
 
    // data members
    TString       fDEName;     ///< detection element name
    TString       fVolumePath; ///< \brief the full path of aligned volume
                               ///  or envelope in geometry
    TGeoHMatrix*  fLocalTransformation;  ///< the transformation wrt module
    TGeoHMatrix*  fGlobalTransformation; ///< the transformation wrt world

  ClassDef(AliMUONGeometryDetElement,2) // MUON det element transformations
};

// inline functions

/// Return module name prefix
inline TString AliMUONGeometryDetElement::GetDENamePrefix()
{ return fgkDENamePrefix; }

/// Set the full path of the aligned volume or envelope in geometry
inline void AliMUONGeometryDetElement::SetVolumePath(const TString& volumePath)
{ fVolumePath = volumePath; }

/// Return detection element ID
inline Int_t AliMUONGeometryDetElement::GetId() const
{ return GetUniqueID(); }

/// Return detection element ID
inline TString AliMUONGeometryDetElement::GetDEName() const
{ return fDEName; }

/// Return the full path of the aligned volume or envelope in geometry
inline TString AliMUONGeometryDetElement::GetVolumePath() const
{ return fVolumePath; }

/// Return the detection element transformation wrt module
inline const TGeoHMatrix* 
AliMUONGeometryDetElement::GetLocalTransformation() const
{ return fLocalTransformation; }

/// Return the detection element transformation wrt world
inline const TGeoHMatrix* 
AliMUONGeometryDetElement::GetGlobalTransformation() const
{ return fGlobalTransformation; }

#endif //ALI_MUON_GEOMETRY_DET_ELEMENT_H
