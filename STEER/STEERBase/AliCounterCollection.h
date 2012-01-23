#ifndef ALICOUNTERCOLLECTION_H
#define ALICOUNTERCOLLECTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup PWG3muon
/// \class AliCounterCollection
/// \brief generic class to handle a collection of counters
// Author: Philippe Pillot

#include <TNamed.h>

class TString;
class TObjArray;
class THnSparse;
class THashList;
class TArrayI;
class TH1D;
class TH2D;
class TCollection;

class AliCounterCollection : public TNamed {
public:
  
  AliCounterCollection(const char* name = "counters");
  virtual ~AliCounterCollection();
  
  virtual void Clear(Option_t* = "");
  
  // Add a new rubric with the complete list of related key words separated by "/"
  void AddRubric(TString name, TString listOfKeyWords);
  // Add a new rubric containing at maximum maxNKeyWords key words
  void AddRubric(TString name, Int_t maxNKeyWords);
  // Initialize the internal counters from the added rubrics
  void Init(Bool_t weightedCounters = kFALSE);
  
  // return the list of key words for the given rubric
  TString GetKeyWords(TString rubric) const;
  
  // Add "value" to the counter referenced by "externalKey"
  void Count(TString externalKey, Int_t value = 1);
  void Count(TString externalKey, Double_t value);
  
  // Get the overall statistics for the given selection (result is integrated over not specified rubrics)
  Double_t GetSum(TString selections = "", Bool_t* longCounters = 0x0);
  // Get counters of the rubric "rubric1" for the given "selection"
  TH1D* Get(TString rubric1, TString selections);
  // Get counters of the "rubric1" vs "rubric2" for the given "selection"
  TH2D* Get(TString rubric1, TString rubric2, TString selections);
  
  // Print every individual counters if opt=="" or call "Print(opt, "")".
  virtual void Print(const Option_t* opt = "") const;
  // Print the full list of key words
  void PrintKeyWords() const;
  // Print value of selected counter
  void PrintValue(TString selections);
  // Print desired rubrics for the given selection
  void Print(TString rubrics, TString selections, Bool_t removeEmpty = kFALSE);
  // Print the overall statistics for the given selection (result is integrated over not specified rubrics)
  void PrintSum(TString selections = "");
  
  /// Overload TObject::Draw(Option_t*): Call "Draw(TString rubric1=opt, TString selections="")"
  virtual void Draw(Option_t* opt = "") {Draw(opt, "");}
  // Draw counters of the rubric "rubric1" for the given "selection"
  TH1D* Draw(TString rubric1, TString selections);
  // Draw counters of the "rubric1" vs "rubric2" for the given "selection"
  TH2D* Draw(TString rubric1, TString rubric2, TString selections);
  
  // Add the given AliCounterCollections to this
  void Add(const AliCounterCollection* counter);
  
  // Merge this with a list of AliCounterCollections
  Long64_t Merge(TCollection* list);
  
  // Sort rubrics defined without a list of authorized key words or all rubrics if opt=="all"
  void Sort(Option_t* opt = "", Bool_t asInt = kFALSE);
  /// Sort only that rubric. If asInt=kTRUE, key words are ordered as interger instead of alphabetically
  void SortRubric(TString rubric, Bool_t asInt = kFALSE);
  
private:
  
  /// Not implemented
  AliCounterCollection(const AliCounterCollection& rhs);
  /// Not implemented
  AliCounterCollection& operator = (const AliCounterCollection& rhs);
  
  // return the number of labels in that rubric
  Int_t GetNActiveBins(Int_t dim);
  // return kTRUE if that rubric contains the keyWord "ANY"
  Bool_t ContainsAny(Int_t dim);
    
  // Return the corresponding bins ordered by rubric or 0x0 if externalKey is not valid
  const Int_t* FindBins(const TString& externalKey, Bool_t allocate, Int_t& nEmptySlots);
  // Return the dimension corresponding to that rubric (or -1)
  Int_t FindDim(const TString& rubricName) const;
  // Return the bin number corresponding to that key word (or -1)
  Int_t FindBin(Int_t dim, const TString& keyWord, Bool_t allocate);
  
  // Tag the selected keywords in each rubric (-1=subtract; 0=discard; 1=add)
  Short_t** DecodeSelection(const TString& selections, const TObjArray& displayedRubrics);
  // Tag the selected keywords (separated by ',') in that rubric (-1=subtract; 0=discard; 1=add)
  Bool_t Select(Bool_t include, const TString& rubric, const TString& keywords, Bool_t displayed, Short_t* selectBins[]);
    
  // Make sure all strings appear only once in this list
  void CleanListOfStrings(TObjArray* list);
  
  // Add "value" to the counter referenced by "externalKey"
  void CountAsDouble(TString externalKey, Double_t value);
  
  // Print the content of 1D histogram as a list
  void PrintList(const TH1D* hist, Bool_t removeEmpty, Bool_t longCounters) const;
  // Print the content of 2D histogram as an array
  void PrintArray(const TH2D* hist, Bool_t removeEmpty, Bool_t longCounters) const;
  // Print the content of nD histogram as a list of arrays
  void PrintListOfArrays(const THnSparse* hist, Bool_t removeEmpty, Bool_t longCounters) const;
  
  // Return the number of characters of the longest label
  Int_t GetMaxLabelSize(THashList* labels) const;
  
  // Return desired "data" for the given "selection" stored in a new histogram or 0x0
  TObject* Projection(const TObjArray& data, const TString& selections, Bool_t& longCounters);
  
  // Consistency check of the two counter collections
  Int_t* CheckConsistency(const AliCounterCollection* c);
  
  // Sort labels (alphabetically or as integer) in each rubric flagged in "rubricsToSort"
  void Sort(const Bool_t* rubricsToSort, Bool_t asInt);
  // Return a list (not owner) of labels sorted assuming they are integers
  THashList* SortAsInt(const THashList* labels);
  
  // Convert the given THnSparse to a THnSparseL (able to handle numbers >= 2^31)
  void ConvertToTHnSparseL(THnSparse* &h);
  
private:
  
  THashList* fRubrics;        ///< list of rubrics with associated key words
  TArrayI*   fRubricsSize;    ///< maximum number of key words in the corresponding rubric
  THnSparse* fCounters;       ///< histogram of nRubrics dimensions used as n-dimensional counter
  Bool_t fWeightedCounters;   ///< use THnSparseF instead of THnSparseI
  Bool_t fLongCounters;       ///< use THnSparseL instead of THnSparseI
  
  ClassDef(AliCounterCollection, 3); // collection of mergeable counters
};

#endif

