/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
// Class AliMUONGeometryStore
// --------------------------
// The class contains the array of the detection elements,
// which are sorted using the AliMUONVGeometryDEIndexing class.
// The class provides fast access to detection element via detElemId.
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_STORE_H
#define ALI_MUON_GEOMETRY_STORE_H

#include <TObject.h>
#include <TObjArray.h>

class AliMUONVGeometryDEIndexing;

class AliMUONGeometryStore : public TObject
{
  public:
    AliMUONGeometryStore(AliMUONVGeometryDEIndexing* indexing,
                         Bool_t isOwner = true);
    AliMUONGeometryStore();
    virtual ~AliMUONGeometryStore();

    // methods
    void Add(Int_t objectId, TObject* object);  

    // get methods
    TObject* Get(Int_t objectId, Bool_t warn = true) const;

    // methods for looping
    Int_t     GetNofEntries() const;
    TObject*  GetEntry(Int_t index) const;

  protected:
    AliMUONGeometryStore(const AliMUONGeometryStore& rhs);

    // operators  
    AliMUONGeometryStore& operator = (const AliMUONGeometryStore& rhs);
  
  private:
    // static data members
    static const Int_t fgkInitSize; // Initial size of array of objects

    // data members
    TObjArray                   fObjects;// The array of detection elements
    AliMUONVGeometryDEIndexing* fDEIndexing; // The DE indexing definitions

  ClassDef(AliMUONGeometryStore,1) // MUON geometry store
};

#endif //ALI_MUON_GEOMETRY_STORE_H
