#ifndef ALIMUONCHECKITEM_H
#define ALIMUONCHECKITEM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMUONCheckItem
/// \brief A structure used to gather information at different levels (ch,manu,de,chamber)
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif

class AliMpExMap;
class TIterator;

class AliMUONCheckItem : public TNamed
{
public:

  AliMUONCheckItem(Int_t id, Int_t maxNumber, const char* name);
  virtual ~AliMUONCheckItem();
  
  /// Return the identifier of this item
  Int_t GetID() const { return fID; }
  
  TObject* GetItem(Int_t id) const;
  Bool_t AddItem(Int_t id, TObject* item);
  
  Bool_t IsFull() const;
  Bool_t IsDead() const;
  
  void Print(Option_t* opt="") const;
  
    TIterator* CreateIterator() const;
    
private:
  /// Not implemented
  AliMUONCheckItem(const AliMUONCheckItem&);
  /// Not implemented
  AliMUONCheckItem& operator=(const AliMUONCheckItem&);

  void ComputeDead() const;
  
private:
  Int_t fID; //!<! identifier of this item
  mutable Int_t fDead; //!<! whether this object is completely dead
  Int_t fMaximum; //!<! maximum number of sub-items possible within this item
  AliMpExMap* fMissing; //!<! pointers to the sub-items
  
  ClassDef(AliMUONCheckItem,1) // A composite object
};

#endif
