/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
// Class AliMUONVGeometryDEIndexing
// --------------------------------
// The abstract singleton base class for definition of
// the conversion between the detection element Ids and 
// the indexing in a simple array.
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_V_GEOMETRY_DE_INDEXING_H
#define ALI_MUON_V_GEOMETRY_DE_INDEXING_H

#include <TObject.h>

class AliMUONVGeometryDEIndexing;

class AliMUONVGeometryDEIndexing : public TObject
{
  public:
    AliMUONVGeometryDEIndexing();
    virtual ~AliMUONVGeometryDEIndexing();

    // static method
    static  Int_t GetModuleId(Int_t detElemId);
            
    // methods
    virtual Int_t GetDetElementIndex(Int_t detElemId) const = 0;
    virtual Int_t GetDetElementId(Int_t detElemIndex) const = 0;

    virtual Int_t GetNofDetElements() const = 0;
    virtual void  SetNofDetElements(Int_t  nofDetElements) = 0;  

  ClassDef(AliMUONVGeometryDEIndexing, 1) // MUON transformations store
};

#endif //ALI_MUON_DE_INDEXING_H
