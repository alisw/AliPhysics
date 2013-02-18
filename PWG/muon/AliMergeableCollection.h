#ifndef ALIMERGEABLECOLLECTION_H
#define ALIMERGEABLECOLLECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id: AliMergeableCollection.h 50593 2011-07-14 17:42:28Z martinez $

///////////////////////////////////////////////////////////////////////////////
///
/// AliMergeableCollection
///
/// Collection of mergeable objects, indexed by key-tuples
///
/// Important point is that AliMergeableCollection is *always* the
/// owner of the objects it holds. This is why you should not
/// use the (inherited from TCollection) Add() method but the Adopt() methods
///
/// \author Diego Stocco

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif
#ifndef ROOT_TCollection
#  include "TCollection.h"
#endif
#include "Riostream.h"
#include <map>
#include <string>

class TMap;
class AliMergeableCollectionIterator;
class TH1;

class AliMergeableCollection : public TNamed
{
  friend class AliMergeableCollectionIterator; // our iterator class

public:

  AliMergeableCollection(const char* name="", const char* title="");
  virtual ~AliMergeableCollection();

  virtual AliMergeableCollection* Clone(const char* name="") const;
  
  Bool_t Adopt(TObject* obj);
  Bool_t Adopt(const char* identifier, TObject* obj);
    
  virtual void Clear(Option_t *option="") { Delete(option); }
  
  virtual TObject* FindObject(const char* fullIdentifier) const;

  virtual TObject* FindObject(const TObject* object) const;

  virtual void Delete(Option_t *option="");
  
  virtual Int_t NumberOfObjects() const;

  virtual Int_t NumberOfKeys() const;

  TObject* GetObject(const char* fullIdentifier) const;
  TObject* GetObject(const char* identifier, const char* objectName) const;

  TH1* Histo(const char* fullIdentifier) const;
  TH1* Histo(const char* identifier, const char* objectName) const;

  virtual TIterator* CreateIterator(Bool_t dir = kIterForward) const;
  
  virtual TList* CreateListOfKeys(Int_t index) const;
  
  virtual TList* CreateListOfObjectNames(const char* identifier) const;
  
  virtual TObject* Remove(const char* fullIdentifier);
  
  Int_t RemoveByType(const char* typeName);
  
  TString GetKey(const char* identifier, Int_t index, Bool_t idContainsObjName = kFALSE) const;
  TString GetIdentifier(const char* fullIdentifier) const;
  TString GetObjectName(const char* fullIdentifier) const;
  
  void Print(Option_t *option="") const;
  
  void ClearMessages();
  void PrintMessages(const char* prefix="") const;
  
  Long64_t Merge(TCollection* list);
  
  AliMergeableCollection* Project(const char* identifier) const;
  
  UInt_t EstimateSize(Bool_t show=kFALSE) const;
  
  /// Turn on the display of empty objects for the Print method
  void ShowEmptyObjects(Bool_t show=kTRUE) {
    fMustShowEmptyObject = show;
  }
  
  void PruneEmptyObjects();
  
  static Bool_t MergeObject(TObject* baseObject, TObject* objToAdd);
  
  TObject* GetSum(const char* idPattern);
  
  Bool_t IsEmptyObject(TObject* obj) const;
  
private:
  
  AliMergeableCollection(const AliMergeableCollection& rhs);
  AliMergeableCollection& operator=(const AliMergeableCollection& rhs);
  
  TH1* HistoWithAction(const char* identifier, TObject* o, const TString& action) const;

  Bool_t InternalAdopt(const char* identifier, TObject* obj);
  
  TString InternalDecode(const char* fullIdentifier, Int_t index) const;
  
  TObject* InternalObject(const char* identifier, const char* objectName) const;
  
public:
  TObjArray* SortAllIdentifiers() const;

  TString NormalizeName(const char* identifier, const char* action) const;
  
  TMap* Map() const;

private:
  
  mutable TMap* fMap; /// map of TMap of THashList* of TObject*...
  Bool_t fMustShowEmptyObject; /// Whether or not to show empty objects with the Print method
  mutable Int_t fMapVersion; /// internal version of map (to avoid custom streamer...)
  mutable std::map<std::string,int> fMessages; //! log messages
  
  ClassDef(AliMergeableCollection,1) /// A collection of mergeable objects
};

class AliMergeableCollectionIterator : public TIterator
{  
public:
  virtual ~AliMergeableCollectionIterator();
  
  AliMergeableCollectionIterator(const AliMergeableCollection* hcol, Bool_t direction=kIterForward);
  AliMergeableCollectionIterator& operator=(const TIterator &rhs);
  
  const TCollection *GetCollection() const { return 0x0; }

  TObject* Next();
  
  void Reset();
  
private:
  const AliMergeableCollection* fkMergeableCollection; // Mergeable objects collection being iterated
  TIterator* fMapIterator; // Iterator for the internal map
  TIterator* fHashListIterator; // Iterator for the current hash list
  Bool_t fDirection; // forward or reverse
  
  AliMergeableCollectionIterator() : fkMergeableCollection(0x0), fMapIterator(0x0), fHashListIterator(0x0), fDirection(kIterForward) {}
  
  /// not implemented
  AliMergeableCollectionIterator& operator=(const AliMergeableCollectionIterator &rhs);
  /// not implemented
  AliMergeableCollectionIterator(const AliMergeableCollectionIterator &iter);
    
  ClassDef(AliMergeableCollectionIterator,0)  // Mergeable object collection iterator
};

#endif
