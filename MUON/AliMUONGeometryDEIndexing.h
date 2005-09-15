/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
// Class AliMUONGeometryDEIndexing
// -------------------------------
// The class that provides conversion between the detection element Id
// and the index in the array.
// Used in storing DE transformations and segmentations.
// The detection elements numbering:
//    DetElemId = chamberId*100 + detElemNum
//                where  chamberId  = 1, 2, ..., 14
//                       detElemNum = 0, 1, ...
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_DE_INDEXING_H
#define ALI_MUON_GEOMETRY_DE_INDEXING_H

#include "AliMUONVGeometryDEIndexing.h"

class AliMUONGeometryDEIndexing : public AliMUONVGeometryDEIndexing
{
  public:
    AliMUONGeometryDEIndexing(Int_t moduleId, Int_t nofDetElements);
    AliMUONGeometryDEIndexing();
    virtual ~AliMUONGeometryDEIndexing();
    
    // methods for conversion between det element Id and index
    virtual Int_t GetDetElementIndex(Int_t detElemId) const;
    virtual Int_t GetDetElementId(Int_t detElemIndex) const;

    // set method
    virtual Int_t GetNofDetElements() const ;
    virtual void  SetNofDetElements(Int_t  nofDetElements);  

  private:
    Int_t GetFirstDetElemId() const;

    // data members
    Int_t  fModuleId;       // module Id					           
    Int_t  fNofDetElements; // number of detection elements in the module					           

  ClassDef(AliMUONGeometryDEIndexing,1) // MUON transformations store
};

// inline functions

inline Int_t AliMUONGeometryDEIndexing::GetNofDetElements() const
{ return fNofDetElements; }

inline void AliMUONGeometryDEIndexing::SetNofDetElements(Int_t nofDetElements)
{ fNofDetElements =  nofDetElements; }

#endif //ALI_MUON_GEOMETRY_DE_INDEXING_H
