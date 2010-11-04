//
// Class AliMixInputHandlerInfo
//
// AliMixInputHandlerInfo is interface with mixed 
// input handlers
//
// author: 
//        Martin Vala (martin.vala@cern.ch)
//
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TChainElement.h>

#include "AliLog.h"
#include "AliInputEventHandler.h"

#include "AliMixInputHandlerInfo.h"

ClassImp(AliMixInputHandlerInfo)

//_____________________________________________________________________________
AliMixInputHandlerInfo::AliMixInputHandlerInfo(const char* name, const char* title): TNamed(name, title),
    fChain(0),
    fChainEntriesArray(),
    fZeroEntryNumber(0),
    fNeedNotify(kFALSE) {
  //
  // Default constructor.
  //
}
//_____________________________________________________________________________
AliMixInputHandlerInfo::~AliMixInputHandlerInfo() {
  //
  // Destructor
  //
  if (fChain) delete fChain;
}

//_____________________________________________________________________________
TChain* AliMixInputHandlerInfo::GetChain() {
  //
  // Returns curren chain. When chain is null it will create it
  //
  if (!fChain) fChain = new TChain(GetName());
  return fChain;
}

//_____________________________________________________________________________
void AliMixInputHandlerInfo::AddChain(TChain* chain) {
  //
  // Add chain
  //
  AliDebug(AliLog::kDebug, "<-");

  if (!chain) return;

  if (fChain) delete fChain;
  fChain = new TChain(GetName());
  fChain->Add(chain);

  AliDebug(AliLog::kDebug, "->");
}

//_____________________________________________________________________________
void AliMixInputHandlerInfo::AddTreeToChain(TTree *tree) {
  //
  // Adds tree to chain
  //
  AliDebug(AliLog::kDebug, Form("%s %lld", tree->GetCurrentFile()->GetName(), tree->GetEntries()));

  GetChain();
  fChain->AddFile(tree->GetCurrentFile()->GetName());


  fChainEntriesArray.Set(fChain->GetListOfFiles()->GetEntries());

  AliDebug(AliLog::kDebug, Form("Adding %lld to id %d", tree->GetEntries(), fChain->GetListOfFiles()->GetEntries() - 1));
  fChainEntriesArray.AddAt(tree->GetEntries(), fChain->GetListOfFiles()->GetEntries() - 1);

}

//_____________________________________________________________________________
TChainElement* AliMixInputHandlerInfo::GetEntryInTree(Long64_t& entry) {
  //
  // Get entry in current tree
  //
  fZeroEntryNumber = 0;
  if (entry < fZeroEntryNumber) {
    AliError(Form("Num %lld is less then ZeroEntryNumber(%lld)", entry, fZeroEntryNumber));
    entry = -1;
    return 0;
  }

  Long64_t sumTree = fZeroEntryNumber;
  for (Int_t i = 0;i < fChainEntriesArray.GetSize() ;i++) {
    sumTree += fChainEntriesArray.At(i);
    if (sumTree > entry) {
      sumTree = entry - sumTree + fChainEntriesArray.At(i);
      AliDebug(AliLog::kDebug, Form("Entry in current tree num is %lld with i=%d", sumTree, i));

      entry = sumTree;
      TChainElement *chEl = (TChainElement*) fChain->GetListOfFiles()->At(i);
      AliDebug(AliLog::kDebug, Form("Real filename is %s %s", chEl->GetName(), chEl->GetTitle()));

      AliDebug(AliLog::kDebug, Form("And filename is %s %lld", fChain->GetTree()->GetCurrentFile()->GetName(), fChain->GetEntries()));
      return chEl;
    }
  }

  entry = -1;
  return 0;
}

//_____________________________________________________________________________
void AliMixInputHandlerInfo::PrepareEntry(TChainElement *te, Long64_t entry, AliInputEventHandler *eh) {
  //
  // Prepare Entry 
  //
  if (!te) return;

  if (te) {
    if (entry < 0) {
      AliDebug(AliLog::kDebug, Form("We are creating new chain from file %s ...", te->GetTitle()));
      if (!fChain) {
        fChain = new TChain(te->GetName());
        fChain->AddFile(te->GetTitle());
        fChain->GetEntry(0);
        eh->Init(fChain->GetTree(), "proof");
//       eh->Notify(te->GetTitle());
      }
      fNeedNotify = kTRUE;
      return;
    }

  }

  if (fChain) {
    AliDebug(AliLog::kDebug, Form("Filename is %s", fChain->GetCurrentFile()->GetName()));
    TString fn = fChain->GetCurrentFile()->GetName();
    if (fn.CompareTo(te->GetTitle())) {
      AliDebug(AliLog::kDebug, Form("Filename %s is NOT same ...", te->GetTitle()));
      AliDebug(AliLog::kDebug, Form("We are changing to file %s ...", te->GetTitle()));
      // change file
      delete fChain;
      fChain = new TChain(te->GetName());
      fChain->AddFile(te->GetTitle());
      fChain->GetEntry(0);
      eh->Init(fChain->GetTree(), "proof");

      eh->Notify(te->GetTitle());
      eh->BeginEvent(entry);
      fChain->GetEntry(entry);
      fNeedNotify = kFALSE;
    } else {
      AliDebug(AliLog::kDebug, Form("We are reusing file %s ...", te->GetTitle()));
      if (fNeedNotify) eh->Notify(te->GetTitle());
      fNeedNotify = kFALSE;
      AliDebug(AliLog::kDebug, Form("Entry is %lld + GetEntries %lld ...", entry, fChain->GetEntries()));
      eh->BeginEvent(entry);
      fChain->GetEntry(entry);
      // file is in tree fChain already
    }
  }

  AliDebug(AliLog::kDebug, Form("We are USING file %s ...", te->GetTitle()));
  AliDebug(AliLog::kDebug, Form("We are USING file from fChain->GetTree() %s ...", fChain->GetTree()->GetCurrentFile()->GetName()));

  // here we have correct chain with 1 tree only
  AliDebug(AliLog::kDebug, Form("Entry is %lld ...", entry));

}
