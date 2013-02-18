#ifndef ALIHISTOGRAMCOLLECTION_H
#define ALIHISTOGRAMCOLLECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

///////////////////////////////////////////////////////////////////////////////
///
/// AliHistogramCollection
///
/// Collection of histograms, indexed by key-tuples
///
/// Important point is that AliHistogramCollection is *always* the
/// owner of the histograms it holds. This is why you should not
/// use the (inherited from TCollection) Add() method but the Adopt() methods
///
/// \author Laurent Aphecetche

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

class TH1;
class TMap;
class AliHistogramCollectionIterator;
class AliMergeableCollection;

class AliHistogramCollection : public TNamed
{
  friend class AliHistogramCollectionIterator; // our iterator class

public:

  AliHistogramCollection(const char* name="", const char* title="");
  virtual ~AliHistogramCollection();

  virtual AliHistogramCollection* Clone(const char* name="") const;
  
  AliMergeableCollection* Convert() const;
  
  Bool_t Adopt(TH1* histo);
  Bool_t Adopt(const char* keyA, TH1* histo);
  Bool_t Adopt(const char* keyA, const char* keyB, TH1* histo);
  Bool_t Adopt(const char* keyA, const char* keyB, const char* keyC, TH1* histo);
  Bool_t Adopt(const char* keyA, const char* keyB, const char* keyC, const char* keyD, TH1* histo);
    
  virtual void Clear(Option_t *option="") { Delete(option); }
  
  virtual TObject* FindObject(const char* identifier) const;

  virtual TObject* FindObject(const TObject* key) const;

  virtual void Delete(Option_t *option="");
  
  virtual Int_t NumberOfHistograms() const;

  virtual Int_t NumberOfKeys() const;

  TH1* Histo(const char* identifier) const;
  TH1* Histo(const char* keyA, const char* histoname) const;
  TH1* Histo(const char* keyA, const char* keyB, const char* histoname) const;
  TH1* Histo(const char* keyA, const char* keyB, const char* keyC, const char* histoname) const;
  TH1* Histo(const char* keyA, const char* keyB, const char* keyC, const char* keyD, const char* histoname) const;
  
  virtual TIterator* CreateIterator(Bool_t dir = kIterForward) const;
  
  virtual TList* CreateListOfKeysA() const;
  virtual TList* CreateListOfKeysB() const;
  virtual TList* CreateListOfKeysC() const;
  virtual TList* CreateListOfKeysD() const;
  
  virtual TList* CreateListOfHistogramNames(const char* keyA, const char* keyB="", const char* keyC="", const char* keyD="") const;
  
  virtual TObject* Remove(TObject *obj);

  TString KeyA(const char* identifier) const;
  TString KeyB(const char* identifier) const;
  TString KeyC(const char* identifier) const;
  TString KeyD(const char* identifier) const;
  TString HistoName(const char* identifier) const;
  
  void Print(Option_t *option="") const;
  
  void ClearMessages();
  void PrintMessages(const char* prefix="") const;
  
  Long64_t Merge(TCollection* list);
  
  AliHistogramCollection* Project(const char* keyA, const char* keyB="", const char* keyC="", const char* keyD="") const;
  
  UInt_t EstimateSize(Bool_t show=kFALSE) const;
  
  /// Turn on the display of empty histograms for the Print method
  void ShowEmptyHistograms(Bool_t show=kTRUE) {
    fMustShowEmptyHistogram = show;
  }
  
  void PruneEmptyHistograms();
  
private:
  
  AliHistogramCollection(const AliHistogramCollection& rhs);
  AliHistogramCollection& operator=(const AliHistogramCollection& rhs);

  TList* CreateListOfKeys(Int_t index) const;

  Bool_t InternalAdopt(const char* identifier, TH1* histo);
  
  Bool_t HistoSameAxis(TH1 *h0, TH1 *h1) const;

  TString InternalDecode(const char* identifier, Int_t index) const;
  
  TH1* InternalHisto(const char* identifier, const char* histoname) const;  
  TObjArray* SortAllIdentifiers() const;
  
  TString NormalizeName(const char* identifier, const char* action) const;
  
  TMap* Map() const;

private:
  
  mutable TMap* fMap; // map of TMap of THashList* of TH1*...
  Bool_t fMustShowEmptyHistogram; // Whether or not to show empty histograms with the Print method
  mutable Int_t fMapVersion; // internal version of map (to avoid custom streamer...)
  mutable std::map<std::string,int> fMessages; //! log messages
  
  ClassDef(AliHistogramCollection,7) // A collection of histograms
};

class AliHistogramCollectionIterator : public TIterator
{  
public:
  virtual ~AliHistogramCollectionIterator();
  
  AliHistogramCollectionIterator(const AliHistogramCollection* hcol, Bool_t direction=kIterForward);
  AliHistogramCollectionIterator& operator=(const TIterator &rhs);
  
  const TCollection *GetCollection() const { return 0x0; }

  TObject* Next();
  
  void Reset();
  
private:
  const AliHistogramCollection* fkHistogramCollection; // histogram collection being iterated
  TIterator* fMapIterator; // Iterator for the internal map
  TIterator* fHashListIterator; // Iterator for the current hash list
  Bool_t fDirection; // forward or reverse
  
  AliHistogramCollectionIterator() : fkHistogramCollection(0x0), fMapIterator(0x0), fHashListIterator(0x0), fDirection(kIterForward) {}
  
  /// not implemented
  AliHistogramCollectionIterator& operator=(const AliHistogramCollectionIterator &rhs);
  /// not implemented
  AliHistogramCollectionIterator(const AliHistogramCollectionIterator &iter);
    
  ClassDef(AliHistogramCollectionIterator,0)  // Histogram collection iterator
};

#endif
