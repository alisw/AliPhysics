/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliCounterCollection
/// 
/// generic class to handle a collection of counters
///
/// \author Philippe Pillot
//-----------------------------------------------------------------------------

#include "AliCounterCollection.h"

#include <limits.h>

#include <AliLog.h>

#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <THnSparse.h>
#include <THashList.h>
#include <TArrayI.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCollection.h>

ClassImp(AliCounterCollection)

//-----------------------------------------------------------------------
AliCounterCollection::AliCounterCollection(const char* name) :
TNamed(name,name),
fRubrics(new THashList(10)),
fRubricsSize(new TArrayI(10)),
fCounters(0x0),
fWeightedCounters(kFALSE),
fLongCounters(kFALSE)
{
  /// Constructor
  fRubrics->SetOwner();
}

//-----------------------------------------------------------------------
AliCounterCollection::~AliCounterCollection()
{
  /// Destructor
  delete fRubrics;
  delete fRubricsSize;
  delete fCounters;
}

//-----------------------------------------------------------------------
void AliCounterCollection::Clear(Option_t*)
{
  /// Clear counters
  fRubrics->Clear();
  fRubricsSize->Reset();
  delete fCounters; fCounters = 0x0;
  fWeightedCounters = kFALSE;
  fLongCounters = kFALSE;
}

//-----------------------------------------------------------------------
void AliCounterCollection::AddRubric(TString name, TString listOfKeyWords)
{
  /// Add a new rubric with the complete list of related key words separated by "/".
  /// If the key word "any" is not defined, the overall statistics is
  /// assumed to be the sum of the statistics under each key word.
  
  name.ToUpper();
  listOfKeyWords.ToUpper();
  
  if (fRubrics->Contains(name.Data())) {
    AliError(Form("rubric named %s already exist",name.Data()));
    return;
  }
  
  // add the list of autorized key words
  TObjArray* rubric = listOfKeyWords.Tokenize("/");
  CleanListOfStrings(rubric);
  rubric->SetName(name.Data());
  Int_t nRubrics = fRubrics->GetSize();
  rubric->SetUniqueID(nRubrics);
  fRubrics->AddLast(rubric);
  
  // save the number of autorized key words (expand the array if needed)
  if (nRubrics+1 > fRubricsSize->GetSize()) fRubricsSize->Set(2*fRubricsSize->GetSize());
  (*fRubricsSize)[nRubrics] = rubric->GetEntriesFast();
}

//-----------------------------------------------------------------------
void AliCounterCollection::AddRubric(TString name, Int_t maxNKeyWords)
{
  /// Add a new rubric containing at maximum maxNKeyWords key words.
  /// Key words will be added as the counters get filled until the maximum is reached.
  /// If the key word "any" is never defined, the overall statistics is
  /// assumed to be the sum of the statistics under each key word.
  
  name.ToUpper();
  
  if (fRubrics->Contains(name.Data())) {
    AliError(Form("rubric named %s already exist",name.Data()));
    return;
  }
  
  // create the empty rubric
  TObjString* rubric = new TObjString(name.Data());
  Int_t nRubrics = fRubrics->GetSize();
  rubric->SetUniqueID(nRubrics);
  fRubrics->AddLast(rubric);
  
  // save the maximum number of autorized key words
  if (nRubrics+1 > fRubricsSize->GetSize()) fRubricsSize->Set(2*fRubricsSize->GetSize());
  (*fRubricsSize)[nRubrics] = maxNKeyWords;
}

//-----------------------------------------------------------------------
void AliCounterCollection::Init(Bool_t weightedCounters)
{
  /// Initialize the internal counters from the added rubrics.
  
  // create the counters
  delete fCounters;
  fWeightedCounters = weightedCounters;
  if (fWeightedCounters)
    fCounters = new THnSparseT<TArrayF>("hCounters", "hCounters", fRubrics->GetSize(), fRubricsSize->GetArray(), 0x0, 0x0);
  else
    fCounters = new THnSparseT<TArrayI>("hCounters", "hCounters", fRubrics->GetSize(), fRubricsSize->GetArray(), 0x0, 0x0);
  
  // loop over axis
  TObject* rubric = 0x0;
  TIter nextRubric(fRubrics);
  while ((rubric = nextRubric())) {
    TAxis* axis = fCounters->GetAxis((Int_t)rubric->GetUniqueID());
    
    // set axis name
    axis->SetName(rubric->GetName());
    
    // set labels if already known
    TObjArray* keyWords = dynamic_cast<TObjArray*>(rubric);
    if (keyWords) {
      TObjString* label = 0x0;
      Int_t bin = 1;
      TIter nextLabel(keyWords);
      while ((label = static_cast<TObjString*>(nextLabel()))) axis->SetBinLabel(bin++, label->String().Data());
    }
    
  }
  
}

//-----------------------------------------------------------------------
Int_t AliCounterCollection::GetNActiveBins(Int_t dim)
{
  /// return the number of labels in that rubric.
  THashList* labels = fCounters->GetAxis(dim)->GetLabels();
  return (labels) ? labels->GetSize() : 0;
}

//-----------------------------------------------------------------------
Bool_t AliCounterCollection::ContainsAny(Int_t dim)
{
  /// return kTRUE if that rubric contains the keyWord "ANY".
  THashList* labels = fCounters->GetAxis(dim)->GetLabels();
  return (labels && labels->Contains("ANY"));
}

//-----------------------------------------------------------------------
const Int_t* AliCounterCollection::FindBins(const TString& externalKey, Bool_t allocate, Int_t& nEmptySlots)
{
  /// Return the corresponding bins ordered by rubric or 0x0 if externalKey is not valid.
  /// The externalKey format must be rubric:keyWord/rubric:keyWord/rubric:keyWord/...
  /// If allocate = kTRUE, new key words are added to the corresponding rubric if possible.
  /// If a rubric is not filled in, the coresponding slot contain -1 in the array.
  /// It is the responsability of the user to delete the returned array.
  
  // produce an empty array of keys
  Int_t nRubrics = fRubrics->GetSize();
  Int_t* bins = new Int_t[nRubrics];
  for (Int_t i=0; i<nRubrics; i++) bins[i] = -1;
  nEmptySlots = nRubrics;
  Bool_t isValid = kTRUE;
  
  // get the list of rubric:keyWord pairs
  TObjArray* rubricKeyPairs = externalKey.Tokenize("/");
  
  // loop over each rubric:keyWord pair
  TObjString* pair = 0x0;
  TIter next(rubricKeyPairs);
  while ((pair = static_cast<TObjString*>(next()))) {
    
    // get both rubric and associated key word
    TObjArray* rubricKeyPair = pair->String().Tokenize(":");
    
    // check the format of the pair
    if (rubricKeyPair->GetEntriesFast() != 2) {
      AliError("invalid key format");
      isValid = kFALSE;
      delete rubricKeyPair;
      break;
    }
    
    // get the axis corresponding to that rubric
    Int_t dim = FindDim(static_cast<TObjString*>(rubricKeyPair->UncheckedAt(0))->String());
    if (dim < 0) {
      isValid = kFALSE;
      delete rubricKeyPair;
      break;
    }
    
    // find the bin corresponding to that key word
    Int_t bin = FindBin(dim, static_cast<TObjString*>(rubricKeyPair->UncheckedAt(1))->String(), allocate);
    if (bin < 0) {
      isValid = kFALSE;
      delete rubricKeyPair;
      break;
    }
    
    // check if the array of keys already contains something for that rubric
    if (bins[dim] >= 0) {
      AliWarning("key already given for that rubric --> ignored");
      delete rubricKeyPair;
      continue;
    }
    
    // store the corresponding bin for that slot
    bins[dim] = bin;
    nEmptySlots--;
    
    // clean memory
    delete rubricKeyPair;
  }
  
  // delete the array in case of problem
  if (!isValid) {
    delete[] bins;
    bins = 0x0;
    nEmptySlots = nRubrics;
  }
  
  // clean memory
  delete rubricKeyPairs;
  
  return bins;
}

//-----------------------------------------------------------------------
Int_t AliCounterCollection::FindDim(const TString& rubricName) const
{
  /// Return the dimension corresponding to that rubric (or -1 in case of failure).
  TObject* rubric = fRubrics->FindObject(rubricName.Data());
  if (!rubric) {
    AliError(Form("invalid rubric: %s",rubricName.Data()));
    return -1;
  }
  return (Int_t) rubric->GetUniqueID();
}

//-----------------------------------------------------------------------
Int_t AliCounterCollection::FindBin(Int_t dim, const TString& keyWord, Bool_t allocate)
{
  /// Return the bin number corresponding to that key word (or -1 in case of failure).
  /// If allocate = kTRUE, try to add the key word if possible.
  
  TAxis* axis = fCounters->GetAxis(dim);
  
  // look for the bin corresponding to keyWord
  THashList* labels = axis->GetLabels();
  TObjString* label = (labels) ? static_cast<TObjString*>(labels->FindObject(keyWord.Data())) : 0x0;
  Int_t bin = (label) ? (Int_t)label->GetUniqueID() : -1;
  
  // in case the keyWord does not exist, try to add it if required
  if (bin<0 && allocate) {
    Int_t nLabels = (labels) ? labels->GetSize() : 0;
    if (nLabels < axis->GetNbins()) {
      bin = nLabels+1;
      axis->SetBinLabel(bin, keyWord.Data());
    }
  }
  
  if (bin<0) AliError(Form("invalid key word: %s:%s",axis->GetName(),keyWord.Data()));
  
  return bin;
}

//-----------------------------------------------------------------------
Short_t** AliCounterCollection::DecodeSelection(const TString& selections, const TObjArray& displayedRubrics)
{
  /// Tag the selected keywords in each rubric (-1=subtract; 0=discard; 1=add). Format:
  /// "rubric:[any-]keyWord,keyWord,../rubric:[any-]keyWord,.." (order does not matter).
  /// It is the responsability of the user to delete the returned array.
  
  // produce an empty array of selected keys
  Int_t nRubrics = fCounters->GetNdimensions();
  Short_t** selects = new Short_t*[nRubrics];
  for (Int_t i=0; i<nRubrics; i++) selects[i] = 0x0;
  
  // get the list of rubric:LisOfKeyWord pairs
  TObjArray* rubricKeyPairs = selections.Tokenize("/");
  
  // loop over each rubric:keyWord pair
  TObjString* pair = 0x0;
  TIter next(rubricKeyPairs);
  while ((pair = static_cast<TObjString*>(next()))) {
    
    // get both rubric and associated list of key words
    TObjArray* rubricKeyPair = pair->String().Tokenize(":");
    
    // check the format of the pair
    if (rubricKeyPair->GetEntriesFast() != 2) {
      AliError("invalid key format");
      delete rubricKeyPair;
      delete rubricKeyPairs;
      for (Int_t i=0; i<nRubrics; i++) if (selects[i]) delete[] selects[i];
      delete[] selects;
      return 0x0;
    }
    
    // check wether to select or to discard the keyWords
    Int_t include = kTRUE;
    TString ListOfKeyWords(static_cast<TObjString*>(rubricKeyPair->UncheckedAt(1))->String());
    if (ListOfKeyWords.BeginsWith("ANY-")) {
      ListOfKeyWords.Remove(0,4);
      include = kFALSE;
    }
    
    // select the key words
    const TString& rubric = static_cast<TObjString*>(rubricKeyPair->UncheckedAt(0))->String();
    if (!Select(include, rubric, ListOfKeyWords, displayedRubrics.Contains(rubric.Data()), selects)) {
      delete rubricKeyPair;
      delete rubricKeyPairs;
      for (Int_t i=0; i<nRubrics; i++) if (selects[i]) delete[] selects[i];
      delete[] selects;
      return 0x0;
    }
    
    // clean memory
    delete rubricKeyPair;
  }
  
  // clean memory
  delete rubricKeyPairs;
  
  // complete the selection of other rubrics
  for (Int_t i=0; i<nRubrics; i++) {
    
    // skip already processed rubrics
    if (selects[i]) continue;
    
    // create the list of bins
    Int_t nBins = GetNActiveBins(i) + 1;
    selects[i] = new Short_t[nBins];
    
    // select all key words or only the key work "ANY"
    if (ContainsAny(i) && !displayedRubrics.Contains(fCounters->GetAxis(i)->GetName())) {
      memset(selects[i], 0, sizeof(Short_t) * nBins);
      selects[i][FindBin(i, "ANY", kFALSE)] = 1;
    } else for (Int_t j=0; j<nBins; j++) selects[i][j] = 1;
    
  }
  
  return selects;
}

//-----------------------------------------------------------------------
Bool_t AliCounterCollection::Select(Bool_t include, const TString& rubric, const TString& keywords,
				    Bool_t displayed, Short_t* selectBins[])
{
  /// Tag the selected keywords (separated by ',') in that rubric (-1=subtract; 0=discard; 1=add).
  
  Int_t dim = FindDim(rubric);
  if (dim < 0) return kFALSE;
  
  if (selectBins[dim]) {
    AliWarning(Form("selection already made for rubric %s --> ignored",rubric.Data()));
    return kTRUE;
  }
  
  // get list of key words to select
  TObjArray* keys = keywords.Tokenize(",");
  if (keys->GetEntriesFast() == 0) {
    AliError(Form("no key word specified for rubric %s",rubric.Data()));
    delete keys;
    return kFALSE;
  }
  
  // create the list of bins
  Int_t nBins = GetNActiveBins(dim) + 1;
  selectBins[dim] = new Short_t[nBins];
  
  // select/unselect all bins
  Bool_t containsAny = ContainsAny(dim);
  if (include || (containsAny && !displayed)) {
    memset(selectBins[dim], 0, sizeof(Short_t) * nBins);
    if (!include) selectBins[dim][FindBin(dim, "ANY", kFALSE)] = 1;
  } else for (Int_t j=0; j<nBins; j++) selectBins[dim][j] = 1;
  
  // select/unselect specific key words
  TObjString* key = 0x0;
  TIter nextKey(keys);
  while ((key = static_cast<TObjString*>(nextKey()))) {
    
    // special case of key word "ANY"
    if (key->String() == "ANY") {
      
      if (containsAny) {
	
	Int_t binAny = FindBin(dim, "ANY", kFALSE);
	if (include) selectBins[dim][binAny] = 1;
	else selectBins[dim][binAny] = 0;
	
      } else {
	
	if (include) for (Int_t j=0; j<nBins; j++) selectBins[dim][j] = 1;
	else memset(selectBins[dim], 0, sizeof(Short_t) * nBins);
	
      }
      
    } else { // other cases
      
      // find the corresponding bin
      Int_t bin = FindBin(dim, key->String().Data(), kFALSE);
      if (bin < 0) {
	delete keys;
	return kFALSE;
      }
      
      // select/unselect it
      if (include) selectBins[dim][bin] = 1;
      else if (containsAny && !displayed) selectBins[dim][bin] = -1;
      else selectBins[dim][bin] = 0;
      
    }
    
  }
  
  // clean memory
  delete keys;
  
  return kTRUE;
}

//-----------------------------------------------------------------------
void AliCounterCollection::CleanListOfStrings(TObjArray* list)
{
  /// Make sure all strings appear only once in this list
  
  // remove multiple-occurrence
  Int_t nEntries = list->GetEntriesFast();
  for (Int_t i = 0; i < nEntries; i++) {
    TObjString* entry1 = static_cast<TObjString*>(list->UncheckedAt(i));
    if (!entry1) continue;
    for (Int_t j = i+1; j < nEntries; j++) {
      TObjString* entry2 = static_cast<TObjString*>(list->UncheckedAt(j));
      if (entry2 && entry2->IsEqual(entry1)) {
	AliWarning(Form("multiple-occurence of string \"%s\" --> removed",entry2->String().Data()));
	list->RemoveAt(j);
      }
    }
  }
  
  // remove empty slots
  list->Compress();
}

//-----------------------------------------------------------------------
void AliCounterCollection::Count(TString externalKey, Int_t value)
{
  /// Add "value" to the counter referenced by "externalKey".
  /// The externalKey format must be rubric:keyWord/rubric:keyWord/rubric:keyWord/...
  if (value > 0) CountAsDouble(externalKey, (Double_t)value);
  else if (value < 0) AliError("cannot count negative values");
}

//-----------------------------------------------------------------------
void AliCounterCollection::Count(TString externalKey, Double_t value)
{
  /// Add "value" to the counter referenced by "externalKey".
  /// The externalKey format must be rubric:keyWord/rubric:keyWord/rubric:keyWord/...
  if (fWeightedCounters) CountAsDouble(externalKey, value);
  else AliError("non-weighted counters can only be filled with intergers");
}

//-----------------------------------------------------------------------
void AliCounterCollection::CountAsDouble(TString externalKey, Double_t value)
{
  /// Add "value" to the counter referenced by "externalKey".
  /// The externalKey format must be rubric:keyWord/rubric:keyWord/rubric:keyWord/...
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return;
  }
  
  externalKey.ToUpper();
  
  // convert external to internal key
  Int_t nEmptySlots = 0;
  const Int_t* bins = FindBins(externalKey, kTRUE, nEmptySlots);
  if (!bins) return;
  
  // check for empty slots
  if (nEmptySlots > 0) {
    AliError("incomplete key");
    delete[] bins;
    return;
  }
  
  if (fWeightedCounters || fLongCounters) {
    
    // increment the corresponding counter
    fCounters->AddBinContent(bins, value);
    
  } else {
    
    // switch to long counters if needed before incrementing
    Long64_t linBin = fCounters->GetBin(bins, kTRUE);
    Double_t currentValue = fCounters->GetBinContent(linBin);
    if (currentValue+value > INT_MAX) {
      ConvertToTHnSparseL(fCounters);
      fLongCounters = kTRUE;
      fCounters->AddBinContent(bins, value);
    } else fCounters->AddBinContent(linBin, value);
    
  }
  
  // clean memory
  delete[] bins;
}

//-----------------------------------------------------------------------
void AliCounterCollection::Print(const Option_t* opt) const
{
  /// Print every individual counters if opt=="" or call "Print(opt, "")".
  
  if (strcmp(opt,"")) {
    const_cast<AliCounterCollection*>(this)->Print(opt, "");
    return;
  }
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return;
  }
  
  if (fCounters->GetNbins() == 0) {
    printf("\nall counters are empty\n\n");
    return;
  }
  
  Int_t nRubrics = fCounters->GetNdimensions();
  Int_t* bins = new Int_t[nRubrics];
  
  // loop over every filled counters
  for (Long64_t i=0; i<fCounters->GetNbins(); ++i) {
    
    // get the content of the bin
    Double_t value = fCounters->GetBinContent(i, bins);
    
    // build the corresponding counter name
    TString counter;
    for (Int_t j=0; j<nRubrics; j++) counter += Form("/%s",fCounters->GetAxis(j)->GetBinLabel(bins[j]));
    counter += "/";
    
    // print value
    if (fWeightedCounters) printf("\n%s   %g", counter.Data(), value);
    else if (fLongCounters) printf("\n%s   %ld", counter.Data(), (Long_t)value);
    else printf("\n%s   %d", counter.Data(), (Int_t)value);
  }
  printf("\n\n");
  
  // clean memory
  delete[] bins;
}

//-----------------------------------------------------------------------
TString AliCounterCollection::GetKeyWords(TString rubric) const
{
  /// return the list of key words for the given rubric.
  
  TString keyWords = "";
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return keyWords;
  }
  
  rubric.ToUpper();
  
  // get the dimension corresponding to that rubric
  Int_t dim = FindDim(rubric);
  if (dim < 0) return keyWords;
  
  // build list of key words
  TObjString* label = 0x0;
  TIter nextLabel(fCounters->GetAxis(dim)->GetLabels());
  while ((label = static_cast<TObjString*>(nextLabel()))) keyWords += Form("%s,",label->String().Data());
  keyWords.Remove(TString::kTrailing, ',');
  
  return keyWords;
}

//-----------------------------------------------------------------------
Double_t AliCounterCollection::GetSum(TString selections, Bool_t* longCounters)
{
  /// Get the overall statistics for the given selection (result is integrated over not specified rubrics):
  /// - format of "selections" is rubric:[any-]keyWord,keyWord,../rubric:[any-]keyWord,.. (order does not matter).
  /// - The flag "longCounters" tells whether the histogram contains value(s) larger than INT_MAX.
  
  // set the flag "longCounters" according to the content of current counters
  if (longCounters) *longCounters = fLongCounters;
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return 0.;
  }
  
  selections.ToUpper();
  
  // decode the selections
  Short_t** select = DecodeSelection(selections, TObjArray());
  if (!select) return 0.;
  
  // loop over every filled counters and compute integral
  Double_t sum = 0.;
  Int_t nDims = fCounters->GetNdimensions();
  Int_t* coord = new Int_t[nDims];
  for (Long64_t i=0; i<fCounters->GetNbins(); ++i) {
    
    // get the content of the counter
    Double_t value = fCounters->GetBinContent(i, coord);
    
    // discard not selected counters and compute the selection factor
    Int_t selectionFactor = 1;
    for (Int_t dim = 0; dim < nDims && selectionFactor != 0; dim++) selectionFactor *= select[dim][coord[dim]];
    if (selectionFactor == 0) continue;
    
    // compute integral
    sum += selectionFactor * value;
  }
  
  // check if the sum exceed INT_MAX (only in case of integer counters)
  if (longCounters && !fWeightedCounters && sum > INT_MAX) *longCounters = kTRUE;
  
  // clean memory
  for (Int_t iDim=0; iDim<nDims; iDim++) delete[] select[iDim];
  delete[] select;
  delete[] coord;
  
  return sum;
}

//-----------------------------------------------------------------------
void AliCounterCollection::PrintKeyWords() const
{
  /// Print the full list of key words.
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return;
  }
  
  // loop over rubrics
  Int_t nRubrics = fCounters->GetNdimensions();
  for (Int_t iDim=0; iDim<nRubrics; iDim++) {
    TAxis* axis = fCounters->GetAxis(iDim);
    
    // print rubric's name
    printf("\n%s:", axis->GetName());
    
    // loop over key words
    Bool_t first = kTRUE;
    TObjString* label = 0x0;
    TIter nextLabel(axis->GetLabels());
    while ((label = static_cast<TObjString*>(nextLabel()))) {
      
      //print key word's name
      if (first) {
	printf("%s", label->String().Data());
	first = kFALSE;
      } else printf(",%s", label->String().Data());
      
    }
  }
  printf("\n\n");
}

//-----------------------------------------------------------------------
void AliCounterCollection::PrintValue(TString selections)
{
  /// Print value of selected counter.
  /// format of "selections" is rubric:keyWord/rubric:keyWord/rubric:keyWord/...
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return;
  }
  
  selections.ToUpper();
  
  // convert external to internal key
  Int_t nEmptySlots = 0;
  const Int_t* selectedBins = FindBins(selections, kFALSE, nEmptySlots);
  if (!selectedBins) return;
  
  // check for empty slots
  if (nEmptySlots > 0) {
    AliError("incomplete key");
    delete[] selectedBins;
    return;
  }
  
  // print value
  if (fWeightedCounters) printf("\n%g\n\n", fCounters->GetBinContent(selectedBins));
  else if (fLongCounters) printf("\n%ld\n\n", (Long_t) fCounters->GetBinContent(selectedBins));
  else printf("\n%d\n\n", (Int_t) fCounters->GetBinContent(selectedBins));
  
  // clean memory
  delete[] selectedBins;
}

//-----------------------------------------------------------------------
void AliCounterCollection::Print(TString rubrics, TString selections, Bool_t removeEmpty)
{
  /// Print desired rubrics for the given selection:
  /// - format of "rubrics" is rubric1/rubric2/.. (order matters only for output).
  /// - format of "selections" is rubric:[any-]keyWord,keyWord,../rubric:[any-]keyWord,.. (order does not matter).
  /// If "data" contains 1 rubric, the output will be one counter for each element of that rubric.
  /// If "data" contains 2 rubrics, the output will be an array of counters, rubric1 vs rubric2.
  /// If "data" contains 3 rubrics, the output will be an array rubric1 vs rubric2 for each element in rubric3.
  /// ...
  /// Results are integrated over rubrics not specified neither in "rubrics" nor in "selections".
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return;
  }
  
  rubrics.ToUpper();
  selections.ToUpper();
  
  // get the rubrics to print
  TObjArray* rubricsToPrint = rubrics.Tokenize("/");
  if (rubricsToPrint->GetEntriesFast() == 0) {
    delete rubricsToPrint;
    return;
  }
  
  // remove rubrics called twice
  CleanListOfStrings(rubricsToPrint);
  
  // project counters in the rubrics to print according to the selections
  Bool_t longCounters = kFALSE;
  TObject* hist = Projection(*rubricsToPrint, selections, longCounters);
  if (!hist) {
    delete rubricsToPrint;
    return;
  }
  
  // print counters
  Int_t nRubricsToPrint = rubricsToPrint->GetEntriesFast();
  if (nRubricsToPrint == 1 && (static_cast<TH1D*>(hist))->GetEntries() > 0)
    PrintList(static_cast<TH1D*>(hist), removeEmpty, longCounters);
  else if (nRubricsToPrint == 2 && (static_cast<TH2D*>(hist))->GetEntries() > 0)
    PrintArray(static_cast<TH2D*>(hist), removeEmpty, longCounters);
  else if (nRubricsToPrint > 2 && (static_cast<THnSparse*>(hist))->GetEntries() > 0)
    PrintListOfArrays(static_cast<THnSparse*>(hist), removeEmpty, longCounters);
  else
    printf("\nselected counters are empty\n\n");
  
  // clean memory
  delete rubricsToPrint;
  delete hist;
}

//-----------------------------------------------------------------------
void AliCounterCollection::PrintSum(TString selections)
{
  /// Print the overall statistics for the given selection (result is integrated over not specified rubrics):
  /// - format of "selections" is rubric:[any-]keyWord,keyWord,../rubric:[any-]keyWord,.. (order does not matter).
  Bool_t longCounters = kFALSE;
  Double_t sum = GetSum(selections, &longCounters);
  if (fWeightedCounters) printf("\n%g\n\n", sum);
  else if (longCounters) printf("\n%ld\n\n", (Long_t)sum);
  else printf("\n%d\n\n", (Int_t)sum);
}

//-----------------------------------------------------------------------
void AliCounterCollection::PrintList(const TH1D* hist, Bool_t removeEmpty, Bool_t longCounters) const
{
  /// Print the content of 1D histogram as a list.
  
  // set the format to print labels
  THashList* labels = hist->GetXaxis()->GetLabels();
  TString format = "";
  if (fWeightedCounters) format = Form("\n%%%ds %%10g",GetMaxLabelSize(labels));
  else if (longCounters) format = Form("\n%%%ds %%16ld",GetMaxLabelSize(labels));
  else format = Form("\n%%%ds %%10d",GetMaxLabelSize(labels));
  
  // print value for each label
  TObjString* label = 0x0;
  TIter nextLabel(labels);
  while ((label = static_cast<TObjString*>(nextLabel()))) {
    Int_t bin = (Int_t) label->GetUniqueID();
    Double_t value = hist->GetBinContent(bin);
    if (removeEmpty && value == 0.) continue;
    if (fWeightedCounters) printf(format.Data(), label->String().Data(), value);
    else if (longCounters) printf(format.Data(), label->String().Data(), (Long_t) value);
    else printf(format.Data(), label->String().Data(), (Int_t) value);
  }
  printf("\n\n");
}

//-----------------------------------------------------------------------
void AliCounterCollection::PrintArray(const TH2D* hist, Bool_t removeEmpty, Bool_t longCounters) const
{
  /// Print the content of 2D histogram as an array.
  
  // set the format to print labels in X direction
  THashList* labelsX = hist->GetXaxis()->GetLabels();
  TString formatX(Form("\n%%%ds ",GetMaxLabelSize(labelsX)));
  
  // set the format to print labels in Y direction
  THashList* labelsY = hist->GetYaxis()->GetLabels();
  Int_t maxLabelSizeY = GetMaxLabelSize(labelsY);
  TString formatYv = "";
  if (fWeightedCounters) {
    maxLabelSizeY = TMath::Max(10, maxLabelSizeY);
    formatYv = Form("%%%dg ",maxLabelSizeY);
  } else if (longCounters) {
    maxLabelSizeY = TMath::Max(16, maxLabelSizeY);
    formatYv = Form("%%%dld ",maxLabelSizeY);
  } else {
    maxLabelSizeY = TMath::Max(10, maxLabelSizeY);
    formatYv = Form("%%%dd ",maxLabelSizeY);
  }
  TString formatYs = Form("%%%ds ",maxLabelSizeY);
  
  // if required, set the list of labels for which all counters are not empty
  Bool_t *useLabelX = 0x0, *useLabelY = 0x0;
  TObjString *labelX = 0x0, *labelY = 0x0;
  TIter nextLabelX(labelsX);
  TIter nextLabelY(labelsY);
  if (removeEmpty) {
    
    // create label flags and set them as unused
    useLabelX = new Bool_t[labelsX->GetSize()+1];
    memset(useLabelX, kFALSE, sizeof(Bool_t) * (labelsX->GetSize()+1));
    useLabelY = new Bool_t[labelsY->GetSize()+1];
    memset(useLabelY, kFALSE, sizeof(Bool_t) * (labelsY->GetSize()+1));
    
    // loop over labels in X direction
    while ((labelX = static_cast<TObjString*>(nextLabelX()))) {
      Int_t binX = (Int_t) labelX->GetUniqueID();
      
      // loop over labels in Y direction
      nextLabelY.Reset();
      while ((labelY = static_cast<TObjString*>(nextLabelY()))) {
	Int_t binY = (Int_t) labelY->GetUniqueID();
	
	// both labels already set as used
	if (useLabelX[binX] && useLabelY[binY]) continue;
	
	// skip empty bins
	if (hist->GetBinContent(binX, binY) == 0.) continue;
	
	// set label as used
	useLabelX[binX] = kTRUE;
	useLabelY[binY] = kTRUE;
	
      }
      
    }
    
  }
  
  // print labels in Y axis
  printf(formatX.Data()," ");
  nextLabelY.Reset();
  while ((labelY = static_cast<TObjString*>(nextLabelY()))) {
    if (removeEmpty && !useLabelY[labelY->GetUniqueID()]) continue;
    printf(formatYs.Data(), labelY->String().Data());
  }
  
  // fill array for each label in X axis
  nextLabelX.Reset();
  while ((labelX = static_cast<TObjString*>(nextLabelX()))) {
    Int_t binX = (Int_t) labelX->GetUniqueID();
    
    if (removeEmpty && !useLabelX[binX]) continue;
    
    // print label X
    printf(formatX.Data(), labelX->String().Data());
    
    // print value for each label in Y axis
    nextLabelY.Reset();
    while ((labelY = static_cast<TObjString*>(nextLabelY()))) {
      Int_t binY = (Int_t) labelY->GetUniqueID();
      if (removeEmpty && !useLabelY[binY]) continue;
      if (fWeightedCounters) printf(formatYv.Data(), hist->GetBinContent(binX, binY));
      else if (longCounters) printf(formatYv.Data(), (Long_t) hist->GetBinContent(binX, binY));
      else printf(formatYv.Data(), (Int_t) hist->GetBinContent(binX, binY));
    }
  }
  printf("\n\n");
  
  // clean memory
  if (removeEmpty) {
    delete[] useLabelX;
    delete[] useLabelY;
  }
}

//-----------------------------------------------------------------------
void AliCounterCollection::PrintListOfArrays(const THnSparse* hist, Bool_t removeEmpty, Bool_t longCounters) const
{
  /// Print the content of nD histogram as a list of arrays.
  
  // set the format to print labels in X direction
  THashList* labelsX = hist->GetAxis(0)->GetLabels();
  TString formatX(Form("\n%%%ds ",GetMaxLabelSize(labelsX)));
  
  // set the format to print labels in Y direction
  THashList* labelsY = hist->GetAxis(1)->GetLabels();
  Int_t maxLabelSizeY = GetMaxLabelSize(labelsY);
  TString formatYv = "";
  if (fWeightedCounters) {
    maxLabelSizeY = TMath::Max(10, maxLabelSizeY);
    formatYv = Form("%%%dg ",maxLabelSizeY);
  } else if (longCounters) {
    maxLabelSizeY = TMath::Max(16, maxLabelSizeY);
    formatYv = Form("%%%dld ",maxLabelSizeY);
  } else {
    maxLabelSizeY = TMath::Max(10, maxLabelSizeY);
    formatYv = Form("%%%dd ",maxLabelSizeY);
  }
  TString formatYs(Form("%%%ds ",maxLabelSizeY));
  
  // create a list containing each combination of labels refering the arrays to be printout
  TList listOfCombis;
  listOfCombis.SetOwner();
  
  // add a first empty combination
  Int_t nDim = hist->GetNdimensions();
  listOfCombis.AddLast(new TObjArray(nDim-2));
  
  // loop over the nDim-2 other rubrics
  for (Int_t i=2; i<nDim; i++) {
    
    // save the last label of that rubic
    THashList* labels = hist->GetAxis(i)->GetLabels();
    TObjString* lastLabel = (labels) ? static_cast<TObjString*>(labels->Last()) : 0x0;
    if (!lastLabel) return;
    
    // prepare iteration over the list of labels
    TIter nextLabel(labels);
    
    // loop over existing combinations
    TObjLink* lnk = listOfCombis.FirstLink();
    while (lnk) {
      
      // get the current combination
      TObjArray* currentCombi = static_cast<TObjArray*>(lnk->GetObject());
      
      // loop over labels in the current rubric
      nextLabel.Reset();
      TObjString* label = 0x0;
      while ((label = static_cast<TObjString*>(nextLabel()))) {
	
	// stop at the last one
	if (label == lastLabel) break;
	
	// copy the current combination, add the current label to it and add it to the list of combinations
	TObjArray* combi = new TObjArray(*currentCombi);
	combi->AddLast(label);
	listOfCombis.AddBefore(lnk, combi);
      }
      
      // add the last label to the current combination
      currentCombi->AddLast(lastLabel);
      
      lnk = lnk->Next();
    }
    
  }
  
  // create bin coordinates to access individual counters
  Int_t* bins = new Int_t[nDim];
  
  // create label flags
  Bool_t *useLabelX = 0x0, *useLabelY = 0x0;
  if (removeEmpty) {
    useLabelX = new Bool_t[labelsX->GetSize()+1];
    useLabelY = new Bool_t[labelsY->GetSize()+1];
  }
  
  // loop over each combination of labels
  TObjArray* combi = 0x0;
  TIter nextCombi(&listOfCombis);
  while ((combi = static_cast<TObjArray*>(nextCombi()))) {
    
    // make the name of the combination and fill the corresponding bin coordinates
    TString combiName = "/";
    for (Int_t i=2; i<nDim; i++) {
      TObjString* label = static_cast<TObjString*>(combi->UncheckedAt(i-2));
      combiName += Form("%s/",label->String().Data());
      bins[i] = (Int_t)label->GetUniqueID();
    }
    
    // reset the list of labels for which all counters are not empty
    if (removeEmpty) {
      memset(useLabelX, kFALSE, sizeof(Bool_t) * (labelsX->GetSize()+1));
      memset(useLabelY, kFALSE, sizeof(Bool_t) * (labelsY->GetSize()+1));
    }
    
    Bool_t empty = kTRUE;
    TObjString* labelX = 0x0;
    TObjString* labelY = 0x0;
    TIter nextLabelX(labelsX);
    TIter nextLabelY(labelsY);
    // loop over labels in X direction
    while ((labelX = static_cast<TObjString*>(nextLabelX()))) {
      bins[0] = (Int_t) labelX->GetUniqueID();
      
      // loop over labels in Y direction
      nextLabelY.Reset();
      while ((labelY = static_cast<TObjString*>(nextLabelY()))) {
	bins[1] = (Int_t) labelY->GetUniqueID();
	
	// both labels already set as used
	if (removeEmpty && useLabelX[bins[0]] && useLabelY[bins[1]]) continue;
	
	// skip empty bins
	if (hist->GetBinContent(bins) == 0.) continue;
	
	// set label as used and array as not empty
	empty = kFALSE;
	if (removeEmpty) {
	  useLabelX[bins[0]] = kTRUE;
	  useLabelY[bins[1]] = kTRUE;
	} else break;
	
      }
      
      if (!removeEmpty && !empty) break;
      
    }
    
    // skip empty arrays
    if (empty) continue;
    
    // print the name of the combination of labels refering the incoming array
    printf("\n%s:\n",combiName.Data());
    
    // print labels in Y axis
    printf(formatX.Data()," ");
    nextLabelY.Reset();
    while ((labelY = static_cast<TObjString*>(nextLabelY()))) {
      if (removeEmpty && !useLabelY[labelY->GetUniqueID()]) continue;
      printf(formatYs.Data(), labelY->String().Data());
    }
    
    // fill array for each label in X axis
    nextLabelX.Reset();
    while ((labelX = static_cast<TObjString*>(nextLabelX()))) {
      bins[0] = (Int_t) labelX->GetUniqueID();
      
      if (removeEmpty && !useLabelX[bins[0]]) continue;
      
      // print label X
      printf(formatX.Data(), labelX->String().Data());
      
      // print value for each label in Y axis
      nextLabelY.Reset();
      while ((labelY = static_cast<TObjString*>(nextLabelY()))) {
	bins[1] = (Int_t) labelY->GetUniqueID();
	if (removeEmpty && !useLabelY[bins[1]]) continue;
	if (fWeightedCounters) printf(formatYv.Data(), hist->GetBinContent(bins));
	else if (longCounters) printf(formatYv.Data(), (Long_t) hist->GetBinContent(bins));
	else printf(formatYv.Data(), (Int_t) hist->GetBinContent(bins));
      }
    }
    printf("\n\n");
  }
  
  // clean memory
  delete[] bins;
  if (removeEmpty) {
    delete[] useLabelX;
    delete[] useLabelY;
  }
}

//-----------------------------------------------------------------------
Int_t AliCounterCollection::GetMaxLabelSize(THashList* labels) const
{
  /// Return the number of characters of the longest label.
  Int_t maxLabelSize = 0;
  TObjString* label = 0x0;
  TIter nextLabel(labels);
  while ((label = static_cast<TObjString*>(nextLabel())))
    maxLabelSize = TMath::Max(maxLabelSize, label->String().Length());
  return maxLabelSize;
}

//-----------------------------------------------------------------------
TH1D* AliCounterCollection::Get(TString rubric, TString selections)
{
  /// Get counters of the rubric "rubric" for the given "selection".
  /// Format of "selections" is rubric:[any-]keyWord,keyWord,../rubric:[any-]keyWord,.. (order does not matter).
  /// Results are integrated over rubrics not specified neither in "rubric1" nor in "selections".
  /// It is the responsability of the user to delete the returned histogram.
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return 0x0;
  }
  
  rubric.ToUpper();
  selections.ToUpper();
  
  // fill the rubrics to print
  TObjArray rubricsToPrint(1);
  rubricsToPrint.SetOwner();
  rubricsToPrint.AddLast(new TObjString(rubric.Data()));
  
  // project counters in the rubrics to print according to the selections
  Bool_t longCounters = kFALSE;
  TH1D* hist = static_cast<TH1D*>(Projection(rubricsToPrint, selections, longCounters));
  
  // make it ready to display
  if (hist) {
    
    // remove statistic box
    hist->SetStats(kFALSE);
    
    // prepare X axis
    TAxis* axis = hist->GetXaxis();
    THashList* labels = axis->GetLabels();
    Int_t nLabels = (labels) ? labels->GetSize() : 1;
    axis->SetRange(1,nLabels);
    axis->SetNdivisions(1,kFALSE);
    axis->SetTitle(rubric.Data());
    
    // prepare Y axis
    hist->GetYaxis()->SetTitle("Counts");
  }
  
  return hist;
}

//-----------------------------------------------------------------------
TH2D* AliCounterCollection::Get(TString rubric1, TString rubric2, TString selections)
{
  /// Get counters of the "rubric1" vs "rubric2" for the given "selection".
  /// Format of "selections" is rubric:[any-]keyWord,keyWord,../rubric:[any-]keyWord,.. (order does not matter).
  /// Results are integrated over rubrics not specified neither in "rubric1", "rubric2" nor in "selections".
  /// It is the responsability of the user to delete the returned histogram.
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return 0x0;
  }
  
  rubric1.ToUpper();
  rubric2.ToUpper();
  selections.ToUpper();
  
  // fill the rubrics to print
  TObjArray rubricsToPrint(2);
  rubricsToPrint.SetOwner();
  rubricsToPrint.AddLast(new TObjString(rubric2.Data()));
  rubricsToPrint.AddLast(new TObjString(rubric1.Data()));
  
  // project counters in the rubrics to print according to the selections
  Bool_t longCounters = kFALSE;
  TH2D* hist = static_cast<TH2D*>(Projection(rubricsToPrint, selections, longCounters));
  
  // draw counters
  if (hist) {
    
    // remove statistic box
    hist->SetStats(kFALSE);
    
    // prepare X axis
    TAxis* axisX = hist->GetXaxis();
    THashList* labelsX = axisX->GetLabels();
    Int_t nLabelsX = (labelsX) ? labelsX->GetSize() : 1;
    axisX->SetRange(1,nLabelsX);
    axisX->SetNdivisions(1,kFALSE);
    axisX->SetTitle(rubric2.Data());
    
    // prepare Y axis
    TAxis* axisY = hist->GetYaxis();
    THashList* labelsY = axisY->GetLabels();
    Int_t nLabelsY = (labelsY) ? labelsY->GetSize() : 1;
    axisY->SetRange(1,nLabelsY);
    axisY->SetNdivisions(1,kFALSE);
    axisY->SetTitle(rubric1.Data());
  }
  
  return hist;
}

//-----------------------------------------------------------------------
TH1D* AliCounterCollection::Draw(TString rubric, TString selections)
{
  /// Draw counters of the rubric "rubric" for the given "selection".
  /// Format of "selections" is rubric:[any-]keyWord,keyWord,../rubric:[any-]keyWord,.. (order does not matter).
  /// Results are integrated over rubrics not specified neither in "rubric1" nor in "selections".
  /// It is the responsability of the user to delete the returned histogram.
  TH1D* hist = Get(rubric, selections);
  if (hist) hist->Draw("htext");
  return hist;
}

//-----------------------------------------------------------------------
TH2D* AliCounterCollection::Draw(TString rubric1, TString rubric2, TString selections)
{
  /// Draw counters of the "rubric1" vs "rubric2" for the given "selection".
  /// Format of "selections" is rubric:[any-]keyWord,keyWord,../rubric:[any-]keyWord,.. (order does not matter).
  /// Results are integrated over rubrics not specified neither in "rubric1", "rubric2" nor in "selections".
  /// It is the responsability of the user to delete the returned histogram.
  TH2D* hist = Get(rubric1, rubric2, selections);
  if (hist) hist->Draw("text");
  return hist;
}

//-----------------------------------------------------------------------
TObject* AliCounterCollection::Projection(const TObjArray& data, const TString& selections, Bool_t& longCounters)
{
  /// Return desired "data" for the given "selection" stored in a new histogram or 0x0 in case of failure.
  /// The type of the histogram (TH1D, TH2D or THnSparse) depend on the number of data.
  /// The flag "longCounters" tells whether the histogram contains value(s) larger than INT_MAX.
  /// It is the responsability of the user to delete the returned histogram.
  
  // set the flag "longCounters" according to the content of current counters
  longCounters = fLongCounters;
  
  // decode the selections
  Short_t** select = DecodeSelection(selections, data);
  if (!select) return 0x0;
  
  // define name and dimensions of projection histo
  TString name(fCounters->GetName());
  Int_t nDims = fCounters->GetNdimensions();
  Int_t nTargetDims = data.GetEntriesFast();
  TArrayI targetDims(nTargetDims);
  TArrayI nNewBins(nTargetDims);
  TArrayI* OldToNewCoord = new TArrayI[nTargetDims];
  for (Int_t i=0; i<nTargetDims; i++) {
    
    // histo name
    name += Form("_%s",static_cast<TObjString*>(data.UncheckedAt(i))->String().Data());
    
    // find target dims
    targetDims[i] = FindDim(static_cast<TObjString*>(data.UncheckedAt(i))->String());
    
    // set number of selected bins in the target dims and make the correspondence between old and new coordinates
    nNewBins[i] = 0;
    if (targetDims[i] > -1) {
      Int_t nBins = GetNActiveBins(targetDims[i]) + 1;
      OldToNewCoord[i].Set(nBins);
      for (Int_t j=1; j<nBins; j++) if (select[targetDims[i]][j] > 0) OldToNewCoord[i][j] = ++nNewBins[i];
    }
    
    // clean memory and return 0x0 in case of problem
    if (nNewBins[i] == 0) {
      for (Int_t iDim=0; iDim<nDims; iDim++) delete[] select[iDim];
      delete[] select;
      delete[] OldToNewCoord;
      return 0x0;
    }
    
  }
  
  // define title of projection histo
  TString title = "Selections:  ";
  TString selectionString(selections);
  selectionString.Remove(TString::kBoth, '/');
  if (selectionString.Length() > 0) title += Form("%s/", selectionString.Data());
  TObject* rub = 0x0;
  TIter nextRubric(fRubrics);
  while ((rub = nextRubric())) {
    if (selectionString.Contains(Form("%s:",rub->GetName()))) continue;
    if (data.Contains(rub->GetName())) continue;
    title += Form("%s:ANY/", rub->GetName());
  }
  title.ReplaceAll("/", "  ");
  
  // Create new histograms
  TObject* hist;
  if (nTargetDims == 1) hist = new TH1D(name.Data(), title.Data(), nNewBins[0], 0., 1.);
  else if (nTargetDims == 2) hist = new TH2D(name.Data(), title.Data(), nNewBins[0], 0., 1., nNewBins[1], 0., 1.);
  else if (fWeightedCounters) hist = new THnSparseT<TArrayF>(name.Data(), title.Data(), nTargetDims, nNewBins.GetArray(), 0x0, 0x0);
  else if (fLongCounters) hist = new THnSparseT<TArrayL>(name.Data(), title.Data(), nTargetDims, nNewBins.GetArray(), 0x0, 0x0);
  else hist = new THnSparseT<TArrayI>(name.Data(), title.Data(), nTargetDims, nNewBins.GetArray(), 0x0, 0x0);
  
  // Set new axis labels
  TObjString* label;
  if (nTargetDims < 3) {
    
    // X axis
    TIter nextLabelX(fCounters->GetAxis(targetDims[0])->GetLabels());
    while ((label = static_cast<TObjString*>(nextLabelX()))) {
      if (select[targetDims[0]][label->GetUniqueID()] > 0) {
	static_cast<TH1*>(hist)->GetXaxis()->SetBinLabel(OldToNewCoord[0][label->GetUniqueID()], label->String().Data());
      }
    }
    
    // Y axis if any
    if (nTargetDims == 2) {
      TIter nextLabelY(fCounters->GetAxis(targetDims[1])->GetLabels());
      while ((label = static_cast<TObjString*>(nextLabelY()))) {
	if (select[targetDims[1]][label->GetUniqueID()] > 0) {
	  static_cast<TH1*>(hist)->GetYaxis()->SetBinLabel(OldToNewCoord[1][label->GetUniqueID()], label->String().Data());
	}
      }
    }
    
  } else {
    
    // all axes
    for (Int_t i=0; i<nTargetDims; i++) {
      TIter nextLabel(fCounters->GetAxis(targetDims[i])->GetLabels());
      while ((label = static_cast<TObjString*>(nextLabel()))) {
	if (select[targetDims[i]][label->GetUniqueID()] > 0) {
	  static_cast<THnSparse*>(hist)->GetAxis(i)->SetBinLabel(OldToNewCoord[i][label->GetUniqueID()], label->String().Data());
	}
      }
    }
    
  }
  
  // loop over every filled counters
  Int_t* coord = new Int_t[nDims];
  Int_t* newCoord = new Int_t[nTargetDims];
  Int_t nEntries = 0;
  for (Long64_t i=0; i<fCounters->GetNbins(); ++i) {
    
    // get the content of the counter
    Double_t value = fCounters->GetBinContent(i, coord);
    
    // discard not selected counters and compute the selection factor
    Int_t selectionFactor = 1;
    for (Int_t dim = 0; dim < nDims && selectionFactor != 0; dim++) selectionFactor *= select[dim][coord[dim]];
    if (selectionFactor == 0) continue;
    
    // find new coordinates in the projection histo
    for (Int_t d = 0; d < nTargetDims; ++d) newCoord[d] = OldToNewCoord[d][coord[targetDims[d]]];
    
    // fill projection histo
    if (nTargetDims < 3) {
      
      Int_t linBin = (nTargetDims == 1) ? newCoord[0] : static_cast<TH1*>(hist)->GetBin(newCoord[0], newCoord[1]);
      static_cast<TH1*>(hist)->AddBinContent(linBin, selectionFactor*value);
      
      // check if the new value exceed INT_MAX (only in case of integer counters)
      if (!fWeightedCounters && !longCounters && static_cast<TH1*>(hist)->GetBinContent(linBin) > INT_MAX)
	longCounters = kTRUE;
      
    } else if (fWeightedCounters || longCounters) {
      
      static_cast<THnSparse*>(hist)->AddBinContent(newCoord, selectionFactor*value);
      
    } else {
      
      // switch to long counters if needed before filling
      Long64_t linBin = static_cast<THnSparse*>(hist)->GetBin(newCoord, kTRUE);
      Double_t currentValue = static_cast<THnSparse*>(hist)->GetBinContent(linBin);
      if (currentValue+selectionFactor*value > INT_MAX) {
	THnSparse* h = static_cast<THnSparse*>(hist);
	ConvertToTHnSparseL(h);
	hist = h;
	longCounters = kTRUE;
	static_cast<THnSparse*>(hist)->AddBinContent(newCoord, selectionFactor*value);
      } else static_cast<THnSparse*>(hist)->AddBinContent(linBin, selectionFactor*value);
      
    }
    
    nEntries++;
  }
  
  // update the number of entries
  if (nTargetDims < 3) static_cast<TH1*>(hist)->SetEntries(nEntries);
  else static_cast<THnSparse*>(hist)->SetEntries(nEntries);
  
  // clean memory
  for (Int_t iDim=0; iDim<nDims; iDim++) delete[] select[iDim];
  delete[] select;
  delete[] coord;
  delete[] newCoord;
  delete[] OldToNewCoord;
  
  return hist;
}

//-----------------------------------------------------------------------
Int_t* AliCounterCollection::CheckConsistency(const AliCounterCollection* c)
{
  /// Consistency check of the two counter collections. To be consistent, both counters
  /// must have the same rubrics with the same list of authorized key words if any.
  /// Return the correspondence between the local rubric ordering and the one of the other counter,
  /// or 0x0 in case of problem. It is the responsability of the user to delete the returned array.
  
  if (!fCounters || !c->fCounters) {
    AliError("counters are not initialized");
    return 0x0;
  }
  
  // check if both counters are weighted or not
  if (c->fWeightedCounters != fWeightedCounters) AliWarning("merging non-weighted with weigthed counters");
  
  // check if the number of rubrics is the same
  Int_t nRubrics = fRubrics->GetSize();
  if (c->fRubrics->GetSize() != nRubrics) {
    AliError("both counters do not contain the same number of rubrics");
    return 0x0;
  }
  
  Int_t* otherDims = new Int_t[nRubrics];
  
  // loop over local rubrics
  TObject* rubric1 = 0x0;
  TIter nextRubric(fRubrics);
  while ((rubric1 = nextRubric())) {
    
    // find that rubric in the other counter
    TObject* rubric2 = c->fRubrics->FindObject(rubric1->GetName());
    if (!rubric2) {
      AliError(Form("the other counter does not contain the rubric %s", rubric1->GetName()));
      delete[] otherDims;
      return 0x0;
    }
    
    // check the list of authorized key words if any
    TObjArray* keyWords1 = dynamic_cast<TObjArray*>(rubric1);
    TObjArray* keyWords2 = dynamic_cast<TObjArray*>(rubric2);
    if (keyWords1 && keyWords2) {
      
      // check if the number of key words is the same
      if (keyWords1->GetEntriesFast() != keyWords2->GetEntriesFast()) {
	AliError("that rubric does not contain the same number of authorized key words in both counters");
	delete[] otherDims;
	return 0x0;
      }
      
      // loop over local key words
      TObjString* keyWord = 0x0;
      TIter nextKeyWord(keyWords1);
      while ((keyWord = static_cast<TObjString*>(nextKeyWord()))) {
	
	// find that key word in the corresponding rubric of the other counter
	if (!keyWords2->FindObject(keyWord->String().Data())) {
	  AliError(Form("rubric %s does not contain the key word %s in the other counter", rubric1->GetName(), keyWord->String().Data()));
	  delete[] otherDims;
	  return 0x0;
	}
	
      }
      
    } else if (keyWords1 || keyWords2) {
      
      // that rubric has not been initialized the same way in both counter
      if (keyWords1) {
	AliError(Form("rubric %s of the other counter does not contain a list of authorized key words while this does", rubric1->GetName()));
      } else {
	AliError(Form("rubric %s of this counter does not contain a list of authorized key words while the other does", rubric1->GetName()));
      }
      delete[] otherDims;
      return 0x0;
      
    }
    
    // save the correspondence of rubric IDs in both counters
    otherDims[rubric1->GetUniqueID()] = rubric2->GetUniqueID();
    
  }
  
  return otherDims;
}

//-----------------------------------------------------------------------
void AliCounterCollection::Add(const AliCounterCollection* counter)
{
  /// Add the given AliCounterCollections to this. They must have the
  /// same rubrics with the same list of authorized key words if any.
  
  // check the consistency between the other counter and this and get the correspondences between rubric IDs.
  Int_t* otherDims = CheckConsistency(counter);
  if (!otherDims) return;
  
  // switch to long counters if the given counter collection is of that type
  if (counter->fLongCounters && !fLongCounters) {
    ConvertToTHnSparseL(fCounters);
    fLongCounters = kTRUE;
  }
  
  Int_t nRubrics = fCounters->GetNdimensions();
  Int_t* thisBins = new Int_t[nRubrics];
  Int_t* otherBins = new Int_t[nRubrics];
  
  // loop over every filled bins inside the other counter
  for (Long64_t i = 0; i < counter->fCounters->GetNbins(); i++) {
    
    // get the content of the bin
    Double_t value = counter->fCounters->GetBinContent(i, otherBins);
    
    // convert "other" bin coordinates to "this" bin coordinates
    Bool_t ok = kTRUE;
    for (Int_t dim = 0; dim < nRubrics; dim++) {
      TString label = counter->fCounters->GetAxis(otherDims[dim])->GetBinLabel(otherBins[otherDims[dim]]);
      thisBins[dim] = FindBin(dim, label, kTRUE);
      if (thisBins[dim] < 0) {
	AliError("this counter is full, unable to add that key word");
	ok = kFALSE;
	break;
      }
    }
    if (!ok) continue;
    
    if (fWeightedCounters || fLongCounters) {
      
      // increment the corresponding local counter
      fCounters->AddBinContent(thisBins, value);
      
    } else {
      
      // switch to long counters if needed before incrementing
      Long64_t linBin = fCounters->GetBin(thisBins, kTRUE);
      Double_t currentValue = fCounters->GetBinContent(linBin);
      if (currentValue+value > INT_MAX) {
	ConvertToTHnSparseL(fCounters);
	fLongCounters = kTRUE;
	fCounters->AddBinContent(thisBins, value);
      } else fCounters->AddBinContent(linBin, value);
      
    }
    
  }
  
  // clean memory
  delete[] otherDims;
  delete[] thisBins;
  delete[] otherBins;
}

//-----------------------------------------------------------------------
Long64_t AliCounterCollection::Merge(TCollection* list)
{
  /// Merge this with a list of AliCounterCollections. All AliCounterCollections provided
  /// must have the same rubrics with the same list of authorized key words if any.
  
  if (!list || !fCounters) return 0;
  if (list->IsEmpty()) return (Long64_t)fCounters->GetEntries();
  
  TIter next(list);
  const TObject* obj = 0x0;
  while ((obj = next())) {
    
    // check that "obj" is an object of the class AliCounterCollection
    const AliCounterCollection* counter = dynamic_cast<const AliCounterCollection*>(obj);
    if (!counter) {
      AliFatal(Form("object named \"%s\" is a %s instead of an AliCounterCollection!", obj->GetName(), obj->ClassName()));
      continue;
    }
    
    // merge counter to this one
    Add(counter);
    
  }
  
  return (Long64_t)fCounters->GetEntries();
}

//-----------------------------------------------------------------------
void AliCounterCollection::Sort(Option_t* opt, Bool_t asInt)
{
  /// Sort rubrics defined without a list of authorized key words or all rubrics if opt=="all".
  /// If asInt=kTRUE, key words are ordered as interger instead of alphabetically.
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return;
  }
  
  Bool_t all = (!strcasecmp(opt, "all"));
  
  Bool_t somethingToSort = kFALSE;
  Int_t nRubrics = fRubrics->GetSize();
  Bool_t* rubricsToSort = new Bool_t[nRubrics];
  memset(rubricsToSort, kFALSE, sizeof(Bool_t) * nRubrics);
  
  // choose rubrics to sort
  TObject* rubric = 0x0;
  TIter nextRubric(fRubrics);
  while ((rubric = nextRubric())) {
    
    if (all || dynamic_cast<TObjString*>(rubric)) {
      
      // check if something to sort
      THashList* labels = fCounters->GetAxis((Int_t)rubric->GetUniqueID())->GetLabels();
      if (!labels || labels->GetSize() < 2) continue;
      
      // select that rubric
      rubricsToSort[(Int_t)rubric->GetUniqueID()] = kTRUE;
      somethingToSort = kTRUE;
      
    }
    
  }
  
  // sort selected rubrics if any
  if (somethingToSort) Sort(rubricsToSort, asInt);
  
  // clean memory
  delete[] rubricsToSort;
}

//-----------------------------------------------------------------------
void AliCounterCollection::SortRubric(TString rubric, Bool_t asInt)
{
  /// Sort only that rubric. If asInt=kTRUE, key words are ordered as interger instead of alphabetically.
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return;
  }
  
  rubric.ToUpper();
  
  // find the rubric to sort
  Int_t dim = FindDim(rubric);
  if (dim < 0) return;
  
  // check if something to sort
  THashList* labels = fCounters->GetAxis(dim)->GetLabels();
  if (!labels || labels->GetSize() < 2) return;
  
  // select that rubric
  Int_t nRubrics = fRubrics->GetSize();
  Bool_t* rubricsToSort = new Bool_t[nRubrics];
  memset(rubricsToSort, kFALSE, sizeof(Bool_t) * nRubrics);
  rubricsToSort[dim] = kTRUE;
  
  // sort it
  Sort(rubricsToSort, asInt);
  
  // clean memory
  delete[] rubricsToSort;
}

//-----------------------------------------------------------------------
void AliCounterCollection::Sort(const Bool_t* rubricsToSort, Bool_t asInt)
{
  /// Sort labels (alphabetically or as integer) in each rubric flagged in "rubricsToSort".
  
  // create a new counter
  THnSparse* oldCounters = fCounters;
  Int_t nRubrics = fRubrics->GetSize();
  if (fWeightedCounters)
    fCounters = new THnSparseT<TArrayF>("hCounters", "hCounters", nRubrics, fRubricsSize->GetArray(), 0x0, 0x0);
  else if (fLongCounters)
    fCounters = new THnSparseT<TArrayL>("hCounters", "hCounters", nRubrics, fRubricsSize->GetArray(), 0x0, 0x0);
  else 
    fCounters = new THnSparseT<TArrayI>("hCounters", "hCounters", nRubrics, fRubricsSize->GetArray(), 0x0, 0x0);
  Int_t** newBins = new Int_t*[nRubrics];
  Bool_t newBinsFilled = kTRUE;
  
  // define the new axes
  for (Int_t i=0; i<nRubrics; i++) {
    TAxis* oldAxis = oldCounters->GetAxis(i);
    TAxis* newAxis = fCounters->GetAxis(i);
    
    // set the name of the new axis
    newAxis->SetName(oldAxis->GetName());
    
    // get old labels
    THashList* oldLabels = oldAxis->GetLabels();
    if (!oldLabels) {
      newBins[i] = 0x0;
      newBinsFilled = kFALSE;
      continue;
    }
    
    // sort them if required
    if (rubricsToSort[i]) {
      if (asInt) { oldLabels = SortAsInt(oldLabels); }
      else { oldLabels->Sort(); }
    }
    
    // set labels in the new axis and save the correspondence between new and old bins
    newBins[i] = new Int_t[oldLabels->GetSize()+1];
    TObjString* label = 0x0;
    Int_t bin = 1;
    TIter nextLabel(oldLabels);
    while ((label = static_cast<TObjString*>(nextLabel()))) {
      newAxis->SetBinLabel(bin, label->String().Data());
      newBins[i][(Int_t)label->GetUniqueID()] = bin;
      bin++;
    }
    
    // clean memory
    if (rubricsToSort[i] && asInt) delete oldLabels;
  }
  
  // fill the new fCounters only if all axes have label(s) defined (otherwise it is empty)
  if (newBinsFilled) {
    
    // fill the new counters
    Int_t* oldCoor = new Int_t[nRubrics];
    Int_t* newCoor = new Int_t[nRubrics];
    for (Long64_t i = 0; i < oldCounters->GetNbins(); i++) {
      Double_t value = oldCounters->GetBinContent(i, oldCoor);
      for (Int_t dim = 0; dim < nRubrics; dim++) newCoor[dim] = newBins[dim][oldCoor[dim]];    
      fCounters->AddBinContent(newCoor, value);
    }
    
    // clean memory
    delete[] oldCoor;
    delete[] newCoor;
  }
  
  // clean memory
  for (Int_t i=0; i<nRubrics; i++) delete[] newBins[i];
  delete[] newBins;
  delete oldCounters;
}

//-----------------------------------------------------------------------
THashList* AliCounterCollection::SortAsInt(const THashList* labels)
{
  /// Return a list (not owner) of labels sorted assuming they are integers.
  /// It is the responsability of user to delete the returned list.
  
  THashList* sortedLabels = new THashList(labels->GetSize());
  TIter nextSortedLabel(sortedLabels);
  
  // loop over labels
  TObjString* label = 0x0;
  TIter nextLabel(labels);
  while ((label = static_cast<TObjString*>(nextLabel()))) {
    
    // find where to add it
    TObjString* sortedLabel = 0x0;
    nextSortedLabel.Reset();
    while ((sortedLabel = static_cast<TObjString*>(nextSortedLabel())) &&
	   (sortedLabel->String().Atoi() <= label->String().Atoi())) {}
    
    // add it
    if (sortedLabel) sortedLabels->AddBefore(sortedLabel, label);
    else sortedLabels->AddLast(label);
  }
  
  return sortedLabels;
}

//-----------------------------------------------------------------------
void AliCounterCollection::ConvertToTHnSparseL(THnSparse* &h)
{
  /// Convert the given THnSparse to a THnSparseL (able to handle numbers >= 2^31)
  
  // create the new THnSparse
  Int_t nDims = h->GetNdimensions();
  Int_t* nBins = new Int_t[nDims];
  for (Int_t i=0; i<nDims; i++) nBins[i] = h->GetAxis(i)->GetNbins();
  THnSparse* hNew = new THnSparseT<TArrayL>("new", "new", nDims, nBins, 0x0, 0x0);
  delete[] nBins;
  
  // transfer the axes
  for (Int_t i=0; i<nDims; i++) {
    TAxis* oldAxis = h->GetAxis(i);
    TAxis* newAxis = hNew->GetAxis(i);
    
    // transfer the name
    newAxis->SetName(oldAxis->GetName());
    
    // transfer labels
    TObjString* label = 0x0;
    TIter nextLabel(oldAxis->GetLabels());
    while ((label = static_cast<TObjString*>(nextLabel())))
      newAxis->SetBinLabel(label->GetUniqueID(), label->String().Data());
  }
  
  // fill the new THnSparse
  Int_t* coor = new Int_t[nDims];
  for (Long64_t i = 0; i < h->GetNbins(); i++) {
    Double_t value = h->GetBinContent(i, coor);
    hNew->AddBinContent(coor, value);
  }
  delete[] coor;
  
  // transfer the number of entries
  hNew->SetEntries(h->GetEntries());
  
  // remove old THnSparse and transfer its name and title to the new one
  TString name(h->GetName());
  TString title(h->GetTitle());
  delete h;
  h = hNew;
  h->SetNameTitle(name.Data(), title.Data());
}

