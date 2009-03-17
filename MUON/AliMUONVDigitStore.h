#ifndef ALIMUONVDIGITSTORE_H
#define ALIMUONVDIGITSTORE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONVDigitStore
/// \brief Interface for a digit container
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVSTORE_H
#  include "AliMUONVStore.h"
#endif

#ifndef ALIMUONVDIGIT_H
#  include "AliMUONVDigit.h" // must be there for covariant return type of FindObjet methods
#endif

class AliMUONVDigitStore : public AliMUONVStore
{  
public:
  
  /// Replacement policy : what to do when adding a digit to the store
  enum EReplacePolicy { kAllow, kDeny, kMerge, kIgnore };
  
public:
  AliMUONVDigitStore();
  virtual ~AliMUONVDigitStore();
  
  /// Add an object, if it is of the right class
  virtual Bool_t Add(TObject* object);
  
  /// Create an (empty) object of the same concrete class as *this
  virtual AliMUONVDigitStore* Create() const = 0;
  
  static AliMUONVDigitStore* Create(TTree& tree);
  
  static AliMUONVDigitStore* Create(const char* classname);
  
  /// Create a digit
  virtual AliMUONVDigit* CreateDigit(Int_t detElemId, Int_t manuId,
                                     Int_t manuChannel, Int_t cathode) const = 0;

  /// Add a digit and return the newly created digit
  virtual AliMUONVDigit* Add(Int_t detElemId, 
                             Int_t manuId,
                             Int_t manuChannel,
                             Int_t cathode,
                             EReplacePolicy replace);
  
  /** Add a (s)digit. Digit is adopted. 
    @param digit the digit to be added
    @param replace specify what to do if the digit is already there.
    kAllow means replacement is allowed, kDeny means it is forbidden (in which
    case we return 0x0), and kMerge means both digits will be merged).
  Finally, kIgnore means no check is done at all. This is the most rapid option,
  but also the more dangerous ;-)
  */
  virtual AliMUONVDigit* Add(const AliMUONVDigit& digit, EReplacePolicy replace) = 0;

  /// Create an iterator to loop over all our digits.
  virtual TIterator* CreateIterator() const = 0;
  
  /** Create an iterator to loop over all digits of a group of detection elements,
    and a given cathode (if cathode != -1)
    */
  virtual TIterator* CreateIterator(Int_t firstDetElemId, 
                                    Int_t lastDetElemId,
                                    Int_t cathode=2) const = 0;
  
  /// Create an iterator to loop over tracker digits only
  virtual TIterator* CreateTrackerIterator() const = 0;

  /// Create an iterator to loop over trigger digits only
  virtual TIterator* CreateTriggerIterator() const = 0;

  using AliMUONVStore::FindObject;

  /// Find an object (default is to forward to FindObject(object->GetUniqueID())
  virtual AliMUONVDigit* FindObject(const TObject* object) const;
  
  /// Find an object by its uniqueID
  virtual AliMUONVDigit* FindObject(UInt_t uniqueID) const;
  
  /// Find a digit by the quadruplet (de,manu,channel,cathode)
  virtual AliMUONVDigit* FindObject(Int_t detElemId, Int_t manuId, Int_t manuChannel, Int_t cathode) const = 0;
  
  /// Number of digits we store
  virtual Int_t GetSize() const = 0;
  
  /// Remove an element
  virtual AliMUONVDigit* Remove(AliMUONVDigit& digit) = 0;
  
  /// Number of digits in a given detection element
  virtual Int_t GetSize(Int_t detElemId) const { return GetSize(detElemId,2); }
  
  /// Number of digits in a given detection element and a given cathode (2 for both cathodes)
  virtual Int_t GetSize(Int_t detElemId, Int_t cathode) const;
  
  /// Whether we have any MC related information (e.g. at least one simulated digit)
  virtual Bool_t HasMCInformation() const = 0;
  
  ClassDef(AliMUONVDigitStore,1) // Digit container interface
};

#endif
