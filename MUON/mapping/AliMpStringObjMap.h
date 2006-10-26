/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup basic
/// \class AliMpStringObjMap
/// \brief Substitutes map <string, TObject> which ALICE does not allow to use 
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MP_STRING_OBJ_MAP_H
#define ALI_MP_STRING_OBJ_MAP_H

#include <TObject.h>
#include <TObjArray.h>
#include <TArrayI.h>

class AliMpStringObjMap : public TObject
{
  public:
    AliMpStringObjMap(Bool_t isOwner = false);
    virtual ~AliMpStringObjMap();
    
    // methods
    Bool_t    Add(const TString& first, TObject* second);
    TObject*  Get(const TString& first) const;
    Int_t     GetNofItems() const;
    virtual void Clear(Option_t* /*option*/ ="");
    virtual void Print(const char* /*option*/ = "") const;
    void Print(const TString& key, ofstream& out) const;
    
  private:
    AliMpStringObjMap(const AliMpStringObjMap& rhs);
    AliMpStringObjMap& operator = (const AliMpStringObjMap& rhs);

    // data members
    Int_t      fNofItems;    ///<  number of items
    TObjArray  fFirstArray;  ///<  first item array
    TObjArray  fSecondArray; ///<  second item array
 
  ClassDef(AliMpStringObjMap,1)  // motif map
};    

#endif //ALI_MP_STRING_OBJ_MAP_H
