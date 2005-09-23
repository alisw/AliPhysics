/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONStringIntMap
/// \brief Substitutes map <string, int> which ALICE does not allow to use 
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_STRING_INT_MAP_H
#define ALI_MUON_STRING_INT_MAP_H

#include <TObject.h>
#include <TObjArray.h>
#include <TArrayI.h>

class TGeoCombiTrans;
class TGeoTranslation;

class AliMUONStringIntMap : public TObject
{
  public:
    AliMUONStringIntMap();
    virtual ~AliMUONStringIntMap();
    
    // methods
    Bool_t  Add(const TString& first, Int_t second);
    Int_t   Get(const TString& first) const;
    Int_t   GetNofItems() const;
    virtual void Clear(Option_t* /*option*/ ="");
    virtual void Print(const char* /*option*/ = "") const;
    void Print(const TString& key, ofstream& out) const;
    
  protected:
    AliMUONStringIntMap(const AliMUONStringIntMap& rhs);

    // operators  
    AliMUONStringIntMap& operator = (const AliMUONStringIntMap& rhs);
 
  private:
    // data members
    Int_t      fNofItems;    // number of items
    TObjArray  fFirstArray;  // first item array
    TArrayI    fSecondArray; // second item array
 
  ClassDef(AliMUONStringIntMap,1)  // motif map
};    

#endif //ALI_MUON_STRING_INT_MAP_H
