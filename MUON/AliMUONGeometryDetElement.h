/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
// Class AliMUONGeometryDetElement
// --------------------------------
// The class defines the detection element.
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_DET_ELEMENT_H
#define ALI_MUON_GEOMETRY_DET_ELEMENT_H

#include <TObject.h>

class TGeoCombiTrans;

class AliMUONGeometryDetElement : public TObject
{
  public:
    AliMUONGeometryDetElement(Int_t detElemId,
                              const TString& alignedVolume, 
			      const TGeoCombiTrans& relTransform);
    AliMUONGeometryDetElement();
    virtual ~AliMUONGeometryDetElement();

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
    void SetGlobalTransformation(const TGeoCombiTrans& transform);
    
    // get methods
    Int_t  GetId() const;
    const TString&        GetAlignedVolume() const;
    const TGeoCombiTrans* GetLocalTransformation() const;
    const TGeoCombiTrans* GetGlobalTransformation() const;

  protected:
    AliMUONGeometryDetElement(const AliMUONGeometryDetElement& rhs);

    // operators  
    AliMUONGeometryDetElement& operator = (const AliMUONGeometryDetElement& rhs);
  
  private:
    // methods
    void PrintTransform(const TGeoCombiTrans* transform) const;
  
    // data members
    TString          fAlignedVolume; // the name of aligned volume or envelope
                                      // representing this detection element
    TGeoCombiTrans*  fLocalTransformation;  // the transformation wrt module
    TGeoCombiTrans*  fGlobalTransformation; // the transformation wrt world

  ClassDef(AliMUONGeometryDetElement,1) // MUON transformations store
};

// inline functions

inline Int_t AliMUONGeometryDetElement::GetId() const
{ return GetUniqueID(); }

inline const TString& AliMUONGeometryDetElement::GetAlignedVolume() const
{ return fAlignedVolume; }

inline const TGeoCombiTrans* 
AliMUONGeometryDetElement::GetLocalTransformation() const
{ return fLocalTransformation; }

inline const TGeoCombiTrans* 
AliMUONGeometryDetElement::GetGlobalTransformation() const
{ return fGlobalTransformation; }

#endif //ALI_MUON_GEOMETRY_DET_ELEMENT_H
