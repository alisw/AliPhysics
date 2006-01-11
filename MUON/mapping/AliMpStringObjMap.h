/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup basic
/// \class AliMpStringObjMap
/// \brief Substitutes map <string, TObject> which ALICE does not allow to use 
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MP_STRING_OBJ_MAP_H
#define ALI_MP_STRING_OBJ_MAP_H

#include <TObject.h>
#include <TObjArray.h>
#include <TArrayI.h>

class AliMpStringObjMap : public TObject
{
  public:
    AliMpStringObjMap();
    virtual ~AliMpStringObjMap();
    
    // methods
    Bool_t    Add(const TString& first, TObject* second);
    TObject*  Get(const TString& first) const;
    Int_t     GetNofItems() const;
    virtual void Clear(Option_t* /*option*/ ="");
    virtual void Print(const char* /*option*/ = "") const;
    void Print(const TString& key, ofstream& out) const;
    
  protected:
    AliMpStringObjMap(const AliMpStringObjMap& rhs);

    // operators  
    AliMpStringObjMap& operator = (const AliMpStringObjMap& rhs);
 
  private:
    // data members
    Int_t      fNofItems;    // number of items
    TObjArray  fFirstArray;  // first item array
    TObjArray  fSecondArray; // second item array
 
  ClassDef(AliMpStringObjMap,1)  // motif map
};    

#endif //ALI_MP_STRING_OBJ_MAP_H
