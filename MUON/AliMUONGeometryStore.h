/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONGeometryStore
/// \brief Array of objects sorted using the AliMUONVGeometryDEIndexing class
///
/// The class contains the array of objects (derived from TObject),
/// which are sorted using the AliMUONVGeometryDEIndexing class.
/// The class provides fast access to detection element via detElemId.
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_STORE_H
#define ALI_MUON_GEOMETRY_STORE_H

#include <TObject.h>
#include <TObjArray.h>

class AliMUONGeometryStore : public TObject
{
  public:
    AliMUONGeometryStore(Bool_t isOwner);
    AliMUONGeometryStore();
    virtual ~AliMUONGeometryStore();

    // static method
    static  Int_t GetModuleId(Int_t detElemId);

    // methods
    void Add(Int_t objectId, TObject* object);  

    // get methods
    TObject* Get(Int_t objectId, Bool_t warn = true) const;

    // methods for looping
    Int_t     GetNofEntries() const;
    TObject*  GetEntry(Int_t index) const;

  protected:
    AliMUONGeometryStore(const AliMUONGeometryStore& rhs);
    AliMUONGeometryStore& operator = (const AliMUONGeometryStore& rhs);
  
  private:
    // static data members
    static const Int_t fgkInitSize;    ///< Initial size of array of objects
    static const Int_t fgkCoefficient; ///< Coefficient used in DE Id <-> Module Id

    // methods
    Int_t GetDEIndex(Int_t detElemId) const;

    // data members
    TObjArray  fObjects; ///< The array of detection elements

  ClassDef(AliMUONGeometryStore,1) // MUON geometry store
};

#endif //ALI_MUON_GEOMETRY_STORE_H
