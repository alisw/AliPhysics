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

//-----------------------------------------------------------------------------
/// \class AliCounterCollection
/// 
/// generic class to handle a collection of counters
///
/// \author Philippe Pillot
//-----------------------------------------------------------------------------

#include "AliCounterCollection.h"

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
fCounters(0x0)
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
void AliCounterCollection::Init()
{
  /// Initialize the internal counters from the added rubrics.
  
  // create the counters
  delete fCounters;
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
  memset(bins, -1, sizeof(Int_t) * nRubrics);
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
  
  if (value < 1) return;
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
  
  // increment the corresponding counter
  fCounters->AddBinContent(bins, (Double_t)value);
  
  // clean memory
  delete[] bins;
}

//-----------------------------------------------------------------------
void AliCounterCollection::Print(const Option_t* opt) const
{
  /// Print every individual counters if opt=="", else call "Print(TString rubrics=opt, TString selections="")".
  
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
    Int_t value = (Int_t) fCounters->GetBinContent(i, bins);
    
    // build the corresponding counter name
    TString counter;
    for (Int_t j=0; j<nRubrics; j++) counter += Form("/%s",fCounters->GetAxis(j)->GetBinLabel(bins[j]));
    counter += "/";
    
    // print value
    printf("\n%s   %d", counter.Data(), value);
  }
  printf("\n\n");
  
  // clean memory
  delete[] bins;
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
  printf("\n%d\n\n", (Int_t) fCounters->GetBinContent(selectedBins));
  
  // clean memory
  delete[] selectedBins;
}

//-----------------------------------------------------------------------
void AliCounterCollection::Print(TString rubrics, TString selections)
{
  /// Print desired rubrics for the given selection:
  /// - format of "rubrics" is rubric1/rubric2/.. (order matters only for output).
  /// - format of "selections" is rubric:keyWord/rubric:keyWord/.. (order does not matter).
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
  TObject* hist = Projection(*rubricsToPrint, selections);
  if (!hist) {
    delete rubricsToPrint;
    return;
  }
  
  // print counters
  Int_t nRubricsToPrint = rubricsToPrint->GetEntriesFast();
  if (nRubricsToPrint == 1 && (static_cast<TH1D*>(hist))->Integral() > 0.)
    PrintList(static_cast<TH1D*>(hist));
  else if (nRubricsToPrint == 2 && (static_cast<TH2D*>(hist))->Integral() > 0.)
    PrintArray(static_cast<TH2D*>(hist));
  else if (nRubricsToPrint > 2 && (static_cast<THnSparse*>(hist))->GetNbins() > 0)
    PrintListOfArrays(static_cast<THnSparse*>(hist));
  else
    printf("\nselected counters are empty\n\n");
  
  // clean memory
  delete rubricsToPrint;
  delete hist;
}

//-----------------------------------------------------------------------
void AliCounterCollection::PrintSum(TString rubric, TString selections)
{
  /// Print the overall statistics under the given rubric for the given selection:
  /// - format of "selections" is rubric:keyWord/rubric:keyWord/.. (order does not matter).
  /// Result is integrated over rubrics not specified neither in "rubric" nor in "selections".
  
  if (!fCounters) {
    AliError("counters are not initialized");
    return;
  }
  
  rubric.ToUpper();
  selections.ToUpper();
  
  // fill the rubric to sum
  TObjArray rubricsToSum(1);
  rubricsToSum.SetOwner();
  rubricsToSum.AddLast(new TObjString(rubric.Data()));
  
  // project counters in the rubric to sum according to the selections
  TH1D* hist = static_cast<TH1D*>(Projection(rubricsToSum, selections));
  if (!hist) return;
  
  // check for empty rubric
  THashList* labels = hist->GetXaxis()->GetLabels();
  if (!labels) {
    printf("\n0\n\n");
    return;
  }
  
  // print the sum of counters under that rubric
  TObjString* any = static_cast<TObjString*>(labels->FindObject("ANY"));
  if (any) printf("\n%d\n\n", (Int_t) hist->GetBinContent((Int_t)any->GetUniqueID()));
  else printf("\n%d\n\n", (Int_t) hist->Integral());
  
  // clean memory
  delete hist;
}

//-----------------------------------------------------------------------
void AliCounterCollection::PrintList(const TH1D* hist) const
{
  /// Print the content of 1D histogram as a list.
  
  // set the format to print labels
  THashList* labels = hist->GetXaxis()->GetLabels();
  TString format(Form("\n%%%ds %%9d",GetMaxLabelSize(labels)));
  
  // print value for each label
  TObjString* label = 0x0;
  TIter nextLabel(labels);
  while ((label = static_cast<TObjString*>(nextLabel()))) {
    Int_t bin = (Int_t) label->GetUniqueID();
    printf(format.Data(), label->String().Data(), (Int_t) hist->GetBinContent(bin));
  }
  printf("\n\n");
}

//-----------------------------------------------------------------------
void AliCounterCollection::PrintArray(const TH2D* hist) const
{
  /// Print the content of 2D histogram as an array.
  
  // set the format to print labels in X direction
  THashList* labelsX = hist->GetXaxis()->GetLabels();
  TString formatX(Form("\n%%%ds ",GetMaxLabelSize(labelsX)));
  
  // set the format to print labels in Y direction and values
  THashList* labelsY = hist->GetYaxis()->GetLabels();
  Int_t maxLabelSizeY = TMath::Max(9, GetMaxLabelSize(labelsY));
  TString formatYs(Form("%%%ds ",maxLabelSizeY));
  TString formatYd(Form("%%%dd ",maxLabelSizeY));
  
  // print labels in Y axis
  printf(formatX.Data()," ");
  TObjString* labelY = 0x0;
  TIter nextLabelY(labelsY);
  while ((labelY = static_cast<TObjString*>(nextLabelY())))
    printf(formatYs.Data(), labelY->String().Data());
  
  // fill array for each label in X axis
  TObjString* labelX = 0x0;
  TIter nextLabelX(labelsX);
  while ((labelX = static_cast<TObjString*>(nextLabelX()))) {
    Int_t binX = (Int_t) labelX->GetUniqueID();
    
    // print label X
    printf(formatX.Data(), labelX->String().Data());
    
    // print value for each label in Y axis
    nextLabelY.Reset();
    while ((labelY = static_cast<TObjString*>(nextLabelY()))) {
      Int_t binY = (Int_t) labelY->GetUniqueID();
      printf(formatYd.Data(), (Int_t) hist->GetBinContent(binX, binY));
    }
  }
  printf("\n\n");
}

//-----------------------------------------------------------------------
void AliCounterCollection::PrintListOfArrays(const THnSparse* hist) const
{
  /// Print the content of nD histogram as a list of arrays.
  
  // set the format to print labels in X direction
  THashList* labelsX = hist->GetAxis(0)->GetLabels();
  TString formatX(Form("\n%%%ds ",GetMaxLabelSize(labelsX)));
  
  // set the format to print labels in Y direction and values
  THashList* labelsY = hist->GetAxis(1)->GetLabels();
  Int_t maxLabelSizeY = TMath::Max(9, GetMaxLabelSize(labelsY));
  TString formatYs(Form("%%%ds ",maxLabelSizeY));
  TString formatYd(Form("%%%dd ",maxLabelSizeY));
  
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
    
    // skip empty array
    Bool_t empty = kTRUE;
    TObjString* labelX = 0x0;
    TObjString* labelY = 0x0;
    TIter nextLabelX(labelsX);
    TIter nextLabelY(labelsY);
    while ((labelX = static_cast<TObjString*>(nextLabelX()))) {
      bins[0] = (Int_t) labelX->GetUniqueID();
      nextLabelY.Reset();
      while ((labelY = static_cast<TObjString*>(nextLabelY()))) {
	bins[1] = (Int_t) labelY->GetUniqueID();
	if (((Int_t) hist->GetBinContent(bins)) > 0) {
	  empty = kFALSE;
	  break;
	}
      }
      if (!empty) break;
    }
    if (empty) continue;
    
    // print the name of the combination of labels refering the incoming array
    printf("\n%s:\n",combiName.Data());
    
    // print labels in Y axis
    printf(formatX.Data()," ");
    nextLabelY.Reset();
    while ((labelY = static_cast<TObjString*>(nextLabelY())))
      printf(formatYs.Data(), labelY->String().Data());
    
    // fill array for each label in X axis
    nextLabelX.Reset();
    while ((labelX = static_cast<TObjString*>(nextLabelX()))) {
      bins[0] = (Int_t) labelX->GetUniqueID();
      
      // print label X
      printf(formatX.Data(), labelX->String().Data());
      
      // print value for each label in Y axis
      nextLabelY.Reset();
      while ((labelY = static_cast<TObjString*>(nextLabelY()))) {
	bins[1] = (Int_t) labelY->GetUniqueID();
	printf(formatYd.Data(), (Int_t) hist->GetBinContent(bins));
      }
    }
    printf("\n\n");
  }
  
  // clean memory
  delete[] bins;
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
TH1D* AliCounterCollection::Draw(TString rubric, TString selections)
{
  /// Draw counters of the rubric "rubric" for the given "selection".
  /// Format of "selections" is rubric:keyWord/rubric:keyWord/.. (order does not matter).
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
  TH1D* hist = static_cast<TH1D*>(Projection(rubricsToPrint, selections));
  
  // draw counters
  if (hist) {
    
    // draw histogram
    hist->Draw("htext");
    hist->SetStats(kFALSE);
    
    // set title
    TString title = "Selections:  ";
    selections.Remove(TString::kBoth, '/');
    if (selections.Length() > 0) title += Form("%s/", selections.Data());
    TObject* rub = 0x0;
    TIter nextRubric(fRubrics);
    while ((rub = nextRubric())) {
      if (selections.Contains(Form("%s:",rub->GetName()))) continue;
      if (rubricsToPrint.Contains(rub->GetName())) continue;
      title += Form("%s:ANY/", rub->GetName());
    }
    title.ReplaceAll("/", "  ");
    hist->SetTitle(title.Data());
    
    // draw X axis
    TAxis* axis = hist->GetXaxis();
    THashList* labels = axis->GetLabels();
    Int_t nLabels = (labels) ? labels->GetSize() : 1;
    axis->SetRange(1,nLabels);
    axis->SetNdivisions(1,kFALSE);
    axis->SetTitle(rubric.Data());
    
    // draw Y axis
    hist->GetYaxis()->SetTitle("Counts");
  }
  
  return hist;
}

//-----------------------------------------------------------------------
TH2D* AliCounterCollection::Draw(TString rubric1, TString rubric2, TString selections)
{
  /// Draw counters of the "rubric1" vs "rubric2" for the given "selection".
  /// Format of "selections" is rubric:keyWord/rubric:keyWord/.. (order does not matter).
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
  TH2D* hist = static_cast<TH2D*>(Projection(rubricsToPrint, selections));
  
  // draw counters
  if (hist) {
    
    // draw histogram
    hist->Draw("text");
    hist->SetStats(kFALSE);
    
    // set title
    TString title = "Selections:  ";
    selections.Remove(TString::kBoth, '/');
    if (selections.Length() > 0) title += Form("%s/", selections.Data());
    TObject* rub = 0x0;
    TIter nextRubric(fRubrics);
    while ((rub = nextRubric())) {
      if (selections.Contains(Form("%s:",rub->GetName()))) continue;
      if (rubricsToPrint.Contains(rub->GetName())) continue;
      title += Form("%s:ANY/", rub->GetName());
    }
    title.ReplaceAll("/", "  ");
    hist->SetTitle(title.Data());
    
    // draw X axis
    TAxis* axisX = hist->GetXaxis();
    THashList* labelsX = axisX->GetLabels();
    Int_t nLabelsX = (labelsX) ? labelsX->GetSize() : 1;
    axisX->SetRange(1,nLabelsX);
    axisX->SetNdivisions(1,kFALSE);
    axisX->SetTitle(rubric2.Data());
    
    // draw Y axis
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
TObject* AliCounterCollection::Projection(const TObjArray& data, const TString& selections)
{
  /// Return desired "data" for the given "selection" stored in a new histogram or 0x0 in case of failure.
  /// The type of the histogram (TH1D, TH2D or THnSparse) depend on the number of data.
  /// It is the responsability of the user to delete the returned histogram.
  
  // get the corresponding dimensions
  Int_t nTargetDim = data.GetEntriesFast();
  Int_t* targetDims = new Int_t[nTargetDim];
  for (Int_t i=0; i<nTargetDim; i++) {
    targetDims[i] = FindDim(static_cast<TObjString*>(data.UncheckedAt(i))->String());
    if (targetDims[i] < 0) {
      delete[] targetDims;
      return 0x0;
    }
  }
  
  // find bins to select
  Int_t nEmptySlots = 0;
  const Int_t* selectedBins = FindBins(selections, kFALSE, nEmptySlots);
  if (!selectedBins) {
    delete[] targetDims;
    return 0x0;
  }
  
  // apply selection for each rubric
  Int_t nRubrics = fCounters->GetNdimensions();
  for (Int_t iDim=0; iDim<nRubrics; iDim++) {
    TAxis* axis = fCounters->GetAxis(iDim);
    
    // select the desired key word
    if (selectedBins[iDim] >= 0) axis->SetRange(selectedBins[iDim], selectedBins[iDim]);
    
    // or select all key words
    else if (data.Contains(axis->GetName())) axis->SetRange();
    
    // or integrate over all cases
    else {
      THashList* labels = axis->GetLabels();
      TObjString* label = (labels) ? static_cast<TObjString*>(labels->FindObject("ANY")) : 0x0;
      Int_t binAny = (label) ? (Int_t)label->GetUniqueID() : -1;
      if (binAny >= 0) axis->SetRange(binAny, binAny);
      else axis->SetRange();
    }
  }
  
  // do projection
  TObject* hist = 0x0;
  if (nTargetDim == 1) {
    
    // project counters to TH1D
    hist = fCounters->Projection(targetDims[0]);
    
    // reset bin labels lost when producing TH1D
    if (selectedBins[targetDims[0]] >= 0)
      static_cast<TH1D*>(hist)->GetXaxis()->SetBinLabel(1, fCounters->GetAxis(targetDims[0])->GetBinLabel(selectedBins[targetDims[0]]));
    else {
      TObjString* label;
      TIter nextLabel(fCounters->GetAxis(targetDims[0])->GetLabels());
      while ((label = static_cast<TObjString*>(nextLabel())))
	static_cast<TH1D*>(hist)->GetXaxis()->SetBinLabel((Int_t)label->GetUniqueID(), label->String().Data());
    }
    
  } else if (nTargetDim == 2) {
    
    // project counters to TH2D (warning X and Y inverted in THnSparse::Projection(X,Y))
    hist = fCounters->Projection(targetDims[1], targetDims[0]);
    
    // reset bin labels in X axis lost when producing TH2D
    if (selectedBins[targetDims[0]] >= 0)
      static_cast<TH2D*>(hist)->GetXaxis()->SetBinLabel(1, fCounters->GetAxis(targetDims[0])->GetBinLabel(selectedBins[targetDims[0]]));
    else {
      TObjString* label;
      TIter nextLabel(fCounters->GetAxis(targetDims[0])->GetLabels());
      while ((label = static_cast<TObjString*>(nextLabel())))
	static_cast<TH2D*>(hist)->GetXaxis()->SetBinLabel((Int_t)label->GetUniqueID(), label->String().Data());
    }
    
    // reset bin labels in Y axis lost when producing TH2D
    if (selectedBins[targetDims[1]] >= 0)
      static_cast<TH2D*>(hist)->GetYaxis()->SetBinLabel(1, fCounters->GetAxis(targetDims[1])->GetBinLabel(selectedBins[targetDims[1]]));
    else {
      TObjString* label;
      TIter nextLabel(fCounters->GetAxis(targetDims[1])->GetLabels());
      while ((label = static_cast<TObjString*>(nextLabel())))
	static_cast<TH2D*>(hist)->GetYaxis()->SetBinLabel((Int_t)label->GetUniqueID(), label->String().Data());
    }
    
  } else {
    
    // project counters to THnSparse (labels are not lost in that case)
    hist = fCounters->Projection(nTargetDim, targetDims);
    
    // reset bin labels in case only one bin has been selected
    for (Int_t i=0; i<nTargetDim; i++) {
      if (selectedBins[targetDims[i]] >= 0) {
	TAxis* axis = static_cast<THnSparse*>(hist)->GetAxis(i);
	axis->GetLabels()->Clear();
	axis->SetBinLabel(1, fCounters->GetAxis(targetDims[i])->GetBinLabel(selectedBins[targetDims[i]]));
      }
    }
    
  }
  
  // clean memory
  delete[] targetDims;
  delete[] selectedBins;
  
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
    
    // increment the corresponding local counter
    fCounters->AddBinContent(thisBins, value);
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
      AliError(Form("object named %s is not AliCounterCollection! Skipping it.", counter->GetName()));
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
  fCounters = new THnSparseT<TArrayI>("hCounters", "hCounters", nRubrics, fRubricsSize->GetArray(), 0x0, 0x0);
  Int_t** newBins = new Int_t*[nRubrics];
  
  // define the new axes
  for (Int_t i=0; i<nRubrics; i++) {
    TAxis* oldAxis = oldCounters->GetAxis(i);
    TAxis* newAxis = fCounters->GetAxis(i);
    
    // set the name of the new axis
    newAxis->SetName(oldAxis->GetName());
    
    // get old labels
    THashList* oldLabels = oldAxis->GetLabels();
    if (!oldLabels) continue;
    
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
  
  // fill the new counters
  Int_t* oldCoor = new Int_t[nRubrics];
  Int_t* newCoor = new Int_t[nRubrics];
  for (Long64_t i = 0; i < oldCounters->GetNbins(); i++) {
    Double_t value = oldCounters->GetBinContent(i, oldCoor);
    for (Int_t dim = 0; dim < nRubrics; dim++) newCoor[dim] = newBins[dim][oldCoor[dim]];    
    fCounters->AddBinContent(newCoor, value);
  }
  
  // clean memory
  for (Int_t i=0; i<nRubrics; i++) delete[] newBins[i];
  delete[] newBins;
  delete[] oldCoor;
  delete[] newCoor;
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

