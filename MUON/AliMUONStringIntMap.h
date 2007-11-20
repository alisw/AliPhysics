/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONStringIntMap
/// \brief Substitutes map <string, int> which ALICE does not allow to use 
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_STRING_INT_MAP_H
#define ALI_MUON_STRING_INT_MAP_H

#include <TObject.h>
#include <TObjArray.h>
#include <TArrayI.h>

class AliMUONStringIntMap : public TObject
{
  public:
    AliMUONStringIntMap();
    virtual ~AliMUONStringIntMap();
    
    // methods
    Bool_t  Add(const TString& first, Int_t second);
    Bool_t  Set(const TString& first, Int_t second);
    Int_t Contains(const TString& first) const;
    
    Int_t   Get(const TString& first) const;
    Int_t   GetNofItems() const;
    virtual void Clear(Option_t* /*option*/ ="");
    virtual void Print(const char* /*option*/ = "") const;
    void Print(const TString& key, ofstream& out) const;

    // Methods for iterating over all elements    
    Bool_t  Next(TString& first, Int_t& second);
    void    ResetItr();

  protected:
    /// Not implemented
    AliMUONStringIntMap(const AliMUONStringIntMap& rhs);
    /// Not implemented
    AliMUONStringIntMap& operator = (const AliMUONStringIntMap& rhs);
 
  private:
    // data members
    Int_t      fNofItems;    ///< number of items
    TObjArray  fFirstArray;  ///< first item array
    TArrayI    fSecondArray; ///< second item array
    Int_t      fCurrentIndex;///< current index

  ClassDef(AliMUONStringIntMap,2)  // motif map
};    

#endif //ALI_MUON_STRING_INT_MAP_H
