#ifndef ALIMUONCHECKITEM_H
#define ALIMUONCHECKITEM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONCheckItem
/// \brief A structure used to gather information at different levels (ch,manu,de,chamber)
/// 
/// \author Laurent Aphecetche

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif

class AliMpExMap;
class TExMapIter;
class AliMUONCheckItemIterator;

class AliMUONCheckItem : public TNamed
{
  friend class AliMUONCheckItemIterator;

public:

  AliMUONCheckItem(Int_t id, Int_t maxNumber, const char* name);
  virtual ~AliMUONCheckItem();
  
  Int_t GetID() const { return fID; }
  
  TObject* GetItem(Int_t id) const;
  Bool_t AddItem(Int_t id, TObject* item);
  
  Bool_t IsFull() const;
  Bool_t IsDead() const;
  
  void Print(Option_t* opt="") const;
  
private:
  AliMUONCheckItem(const AliMUONCheckItem&);
  AliMUONCheckItem& operator=(const AliMUONCheckItem&);
  void ComputeDead() const;
  
private:
  Int_t fID; //!< identifier of this item
  mutable Int_t fDead; //!< whether this object is completely dead
  Int_t fMaximum; //!< maximum number of sub-items possible within this item
  AliMpExMap* fMissing; //!< pointers to the sub-items
  
  ClassDef(AliMUONCheckItem,1) // A composite object
};

#endif
