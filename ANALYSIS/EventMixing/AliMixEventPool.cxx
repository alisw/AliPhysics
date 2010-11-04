//
// Class AliMixEventPool
//
// AliMixEventPool is used to find 
// similar events
//
// author:
//        Martin Vala (martin.vala@cern.ch)
//

#include <TEntryList.h>

#include "AliLog.h"
#include "AliMixEventCutObj.h"

#include "AliMixEventPool.h"

ClassImp(AliMixEventPool)

//_________________________________________________________________________________________________
AliMixEventPool::AliMixEventPool(const char* name, const char* title) : TNamed(name, title),
    fListOfEntryList(),
    fListOfEventCuts(),
    fBinNumber(0) {
  //
  // Default constructor.
  //

  AliDebug(AliLog::kDebug + 5, "<-");
  AliDebug(AliLog::kDebug + 5, "->");
}
//_________________________________________________________________________________________________
AliMixEventPool::AliMixEventPool(const AliMixEventPool& obj) : TNamed(obj),
    fListOfEntryList(obj.fListOfEntryList),
    fListOfEventCuts(obj.fListOfEventCuts),
    fBinNumber(obj.fBinNumber) {
  //
  // Copy constructor
  //
  AliDebug(AliLog::kDebug + 5, "<-");
  AliDebug(AliLog::kDebug + 5, "->");
}
//_________________________________________________________________________________________________
AliMixEventPool::~AliMixEventPool() {
  //
  // Destructor
  //
  AliDebug(AliLog::kDebug + 5, "<-");
  AliDebug(AliLog::kDebug + 5, "->");
}
//_________________________________________________________________________________________________
void AliMixEventPool::AddCut(AliMixEventCutObj* cut) {
  //
  // Adds cut
  //

  if (cut) fListOfEventCuts.Add(cut);
}
//_________________________________________________________________________________________________
void AliMixEventPool::Print(const Option_t* option) const {
  //
  // Prints usefull information
  //

  TObjArrayIter next(&fListOfEventCuts);
  //   Int_t c=0;
  AliMixEventCutObj *cut;
  while ((cut = (AliMixEventCutObj *) next())) {
    cut->Print(option);
  }

  AliInfo(Form("NumOfEntryList %d", fListOfEntryList.GetEntries()));

  TEntryList *el;
  for (Int_t i = 0;i < fListOfEntryList.GetEntries();i++) {
    el = (TEntryList*) fListOfEntryList.At(i);

    AliInfo(Form("EntryList[%d] %lld", i, el->GetN()));
  }
}
//_________________________________________________________________________________________________
Int_t AliMixEventPool::Init() {
  //
  // Init event pool
  //
  AliDebug(AliLog::kDebug+5,"<-");
  CreateEntryListsRecursivly(fListOfEventCuts.GetEntries() - 1);

  fBinNumber++;
  AliDebug(AliLog::kDebug, Form("fBinnumber = %d", fBinNumber));

  AddEntryList();
  AliDebug(AliLog::kDebug+5,"->");
  return 0;
}

//_________________________________________________________________________________________________
void AliMixEventPool::CreateEntryListsRecursivly(Int_t index) {
  //
  // Helper function which create entrylist recursivly
  //
  AliDebug(AliLog::kDebug+5,"<-");
  AliMixEventCutObj *cut;
  if (index >= 0) {
    AliDebug(AliLog::kDebug, Form("index = %d", index));
    cut = dynamic_cast<AliMixEventCutObj*>(fListOfEventCuts.At(index));
    cut->Reset();

    while (cut->HasMore()) {
      cut->AddStep();
      CreateEntryListsRecursivly(index - 1);
      if (cut->HasMore()) {
        fBinNumber++;
        AliDebug(AliLog::kDebug, Form("fBinnumber = %d", fBinNumber));
        AddEntryList();
        //                 PrintCurrentCutIntervals();
      }
    }

  }
  AliDebug(AliLog::kDebug+5,"->");
}

//_________________________________________________________________________________________________
TEntryList* AliMixEventPool::AddEntryList() {
  //
  // Adds endtry list
  //
  
  AliDebug(AliLog::kDebug+5,"<-");
  
  TObjArrayIter next(&fListOfEventCuts);
  AliMixEventCutObj *cut;
  while ((cut = (AliMixEventCutObj*) next())) {
    if (cut) cut->PrintCurrentInterval();
  }

  TEntryList *el = new TEntryList;
  fListOfEntryList.Add(el);

  AliDebug(AliLog::kDebug, Form("Number in Entry list -> %lld", el->GetN()));
  AliDebug(AliLog::kDebug+5,"->");
  return el;
}

//_________________________________________________________________________________________________
Bool_t AliMixEventPool::AddEntry(Long64_t entry, AliVEvent* ev) {
  //
  // Adds entry to correct entry list
  //
  
  AliDebug(AliLog::kDebug+5,"<-");
  AliDebug(AliLog::kDebug + 5, Form("AddEntry(%lld,%p)", entry, ev));
  if (entry < 0) {
    AliDebug(AliLog::kDebug, Form("Entry %lld was NOT added !!!", entry));
    return kFALSE;

  }

  TEntryList *el =  FindEntryList(ev);
  if (el) {
    el->Enter(entry);
    AliDebug(AliLog::kDebug, Form("Entry %lld was added !!!", entry));
    return kTRUE;
  }

  AliDebug(AliLog::kDebug, Form("Entry %lld was NOT added !!!", entry));
  AliDebug(AliLog::kDebug+5,"->");
  return kFALSE;
}

//_________________________________________________________________________________________________
TEntryList* AliMixEventPool::FindEntryList(AliVEvent* ev) {
  //
  // Find entrlist in list of entrlist
  //
  
  AliDebug(AliLog::kDebug+5,"<-");
  const Int_t num = fListOfEventCuts.GetEntries();
  if (num <= 0) return 0;

  Int_t indexes[num];
  Int_t lenght[num];
  Int_t i = 0;
  TObjArrayIter next(&fListOfEventCuts);
  AliMixEventCutObj *cut;
  while ((cut = (AliMixEventCutObj*) next())) {
    indexes[i] = cut->GetIndex(ev);
    if (indexes[i] < 0) {
      AliDebug(AliLog::kDebug, Form("retIndex %d", -1));
      return 0;
    }
    lenght[i] = cut->GetNumberOfBins();
    AliDebug(AliLog::kDebug + 1, Form("indexes[%d] %d", i, indexes[i]));
    i++;
  }

  Int_t retIndex = 0;
  SearchIndexRecursive(fListOfEventCuts.GetEntries() - 1, &indexes[0], &lenght[0], retIndex);
  AliDebug(AliLog::kDebug, Form("retIndex %d", retIndex - 1));
  // index which start with 0 (retIndex-1)
  AliDebug(AliLog::kDebug+5,"->");
  return (TEntryList*) fListOfEntryList.At(retIndex - 1);
}

//_________________________________________________________________________________________________
void AliMixEventPool::SearchIndexRecursive(Int_t num, Int_t* i, Int_t* d, Int_t& index) {
  //
  // Search for index of entrylist
  //
  
  AliDebug(AliLog::kDebug+5,"<-");
  if (num > 0) {
    index += (i[num] - 1) * d[num-1];
    SearchIndexRecursive(num - 1, i, d, index);
  } else {
    index += i[num];
  }
  AliDebug(AliLog::kDebug+5,"->");
}
