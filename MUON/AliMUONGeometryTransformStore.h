// $Id$
//
// Class AliMUONGeometryTransformStore
// -----------------------------------
// The class contains the array of transformations fDETransforms
// of the detection elements from a defined reference frame (MUON chamber ) 
// to the detection element frame.
// The detection elements numbering:
//    DetElemId = chamberId*100 + [50] + detElemNum
//                where  chamberId  = 1, 2, ..., 14
//                       detElemNum = 0, 1, ...
// The number 50 is added to distinguish detector elements 
// in the left and the right hemispheres.
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_TRANSFORM_STORE_H
#define ALI_MUON_GEOMETRY_TRANSFORM_STORE_H

#include <Riostream.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TArrayI.h>

#include "AliMUONGeometrySVMap.h"

class TGeoCombiTrans;
class TGeoTranslation;

class AliMUONGeometryTransformStore : public TObject
{
  public:
    AliMUONGeometryTransformStore(
           Int_t firstDetElemId, Int_t  nofDetElems,
           AliMUONGeometrySVMap* svMap);
    AliMUONGeometryTransformStore();
    virtual ~AliMUONGeometryTransformStore();

    // methods
    void Add(Int_t detElemId,
             const TString& alignedVolume, 
             const TGeoCombiTrans& transformation);  

    virtual void  Print(Option_t* /*option*/) const;
    
    // get methods
    const TGeoCombiTrans* Get(Int_t detElemId) const;
    const TGeoCombiTrans* FindBySensitiveVolume(const TString& volumePath) const;

  protected:
    AliMUONGeometryTransformStore(const AliMUONGeometryTransformStore& rhs);

    // operators  
    AliMUONGeometryTransformStore& operator 
      = (const AliMUONGeometryTransformStore& rhs);
  
  private:
    // methods
    Int_t GetDetElementIndex(Int_t detElemId) const;
    Int_t GetDetElementId(Int_t detElemIndex) const;

    // data members
    static const Int_t   fgkHemisphere;   // The constant to distinguish 
                                          // the left/right hemispere            

    Int_t      fFirstDetElemId; // The first detection element Id 
    Int_t      fNofDetElems;    // Number of detection elements
    TObjArray  fDETransforms;   // The array of transformations
    AliMUONGeometrySVMap*  fSVMap; // The map of sensitive volumes 
                                   // and detector element id 

  ClassDef(AliMUONGeometryTransformStore,1) // MUON transformations store
};

#endif //ALI_MUON_GEOMETRY_TRANSFORM_STORE_H
