#ifndef ALIMUONOBJECTPAIR_H
#define ALIMUONOBJECTPAIR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONObjectPair
/// \brief The equivalent of a std::pair<TObject*,TObject*> ;-)
/// 
// Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONObjectPair : public TObject
{
public:
  AliMUONObjectPair();
  AliMUONObjectPair(TObject* first, 
                    TObject* second,
                    Bool_t isOwnerOfFirst=kTRUE,
                    Bool_t isOwnerOfSecond=kFALSE);
  AliMUONObjectPair(const AliMUONObjectPair& other);
  AliMUONObjectPair& operator=(const AliMUONObjectPair& other);
  
  virtual ~AliMUONObjectPair();

  TObject* First() const { return fFirst; }
  TObject* Second() const { return fSecond; }

  TObject* Key() const { return fFirst; }
  TObject* Value() const { return fSecond; }

  virtual void Copy(TObject& other) const;
  
private:

  TObject* fFirst; ///< first element of the pair
  TObject* fSecond; ///< second element of the pair
  Bool_t fIsOwnerOfFirst; ///< whether we own the first element
  Bool_t fIsOwnerOfSecond; ///<whether we own the second element
  
  ClassDef(AliMUONObjectPair,1) // A pair of TObject*
};

#endif
