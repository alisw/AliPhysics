#include "AliHistoListWrapper.h"
#include "TH2F.h"
#include "TList.h"
#include "AliLog.h"
#include "TString.h"
#include "AliESDInputHandlerRP.h"
#include "AliAnalysisManager.h"
#include "TTree.h"
#include "AliMultiplicity.h"

ClassImp(AliHistoListWrapper)

AliHistoListWrapper::AliHistoListWrapper():
  TNamed(), fList(0)
{
  
  // constructor
  fList = new TList();
  fList->SetOwner();

}

AliHistoListWrapper::AliHistoListWrapper(const char* name, const char* title):
  TNamed(name,title), fList(0)
{
  // constructor

  fList = new TList();
  fList->SetOwner();

}

AliHistoListWrapper::AliHistoListWrapper(const AliHistoListWrapper& obj) : 
  TNamed(obj),
  fList(0)
{
  // Copy ctor
  fList  = obj.fList;
}

AliHistoListWrapper::~AliHistoListWrapper() {
  // Destructor
  if(fList) {
    delete fList;
    fList = 0;
  }

}

Long64_t AliHistoListWrapper::Merge(TCollection* list)
{
  // Merge a list of AliHistoListWrapper objects with this.
  // Returns the number of merged objects (including this).

  // We have to make sure that all the list contain the same histos in
  // the same order. We thus also have to sort the list (sorting is
  // done by name in TList).

  //  AliInfo("Merging");

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  TList collections;

  Int_t count = 0;

  while ((obj = iter->Next())) {
    Bool_t foundDiffinThisIterStep = kFALSE;

    //    Printf("%d - %s",count, obj->GetName());
    AliHistoListWrapper* entry = dynamic_cast<AliHistoListWrapper*> (obj);
    if (entry == 0) 
      continue;

    TList * hlist = entry->fList;

    // Check if all histos in this fList are also in the one from entry and viceversa
    // Use getters to automatically book non defined histos    

    Bool_t areListsDifferent=kTRUE;
    Int_t iloop = 0;
    Int_t max_loops = hlist->GetSize() + fList->GetSize(); // In the worst case all of the histos will be different...
    while(areListsDifferent) {
      if(iloop>max_loops) AliFatal("Infinite Loop?");
      iloop++;
      // sort
      hlist->Sort();
      fList->Sort();
      // loop over the largest 
      TObject * hist =0;
      TIterator * iterlist = 0;
      TList * thislist  = 0; // the list over which I'm iterating
      TList * otherlist = 0; // the other

      if (hlist->GetSize() >= fList->GetSize()) { 
	thislist  = hlist;
	otherlist = fList;
      }
      else{
	thislist  = fList;
	otherlist = hlist;	
      }
      iterlist = thislist->MakeIterator();

      while ((hist= iterlist->Next())){ 
	if(!otherlist->FindObject(hist->GetName())){
	  AliInfo(Form("Adding object %s",hist->GetName()));	  
	  TH1 * hclone =  (TH1*) hist->Clone();
	  if (!hclone->InheritsFrom("TH1")) AliFatal(Form("Found a %s. This class only supports objects inheriting from TH1",hclone->ClassName()));
	  hclone->Reset();
	  otherlist->Add(hclone);
	  foundDiffinThisIterStep=kTRUE;
	}
      }

      // re-sort before checking
      hlist->Sort();
      fList->Sort();

      // check if everything is fine    
      areListsDifferent=kFALSE;
      if (hlist->GetSize() == fList->GetSize()) {	
	Int_t nhist =  fList->GetSize();
	for(Int_t ihist = 0; ihist < nhist; ihist++){
	  if(strcmp(fList->At(ihist)->GetName(),hlist->At(ihist)->GetName())) areListsDifferent = kTRUE;
	}
      } else {
	areListsDifferent=kTRUE;
      }
    }

    // last check: if something is not ok die loudly 
    if (hlist->GetSize() != fList->GetSize()) {
      AliFatal("Mismatching size!");
    }
    Int_t nhist =  fList->GetSize();
    for(Int_t ihist = 0; ihist < nhist; ihist++){
      if(strcmp(fList->At(ihist)->GetName(),hlist->At(ihist)->GetName())){
	AliFatal(Form("Mismatching histos: %s -> %s", fList->At(ihist)->GetName(),hlist->At(ihist)->GetName()));
      }
    }
    
    if (foundDiffinThisIterStep){
      iter->Reset(); // We found a difference: previous lists could
		     // also be affected... We start from scratch
      collections.Clear();
      count = 0;
    }
    else {
      
      collections.Add(hlist);
      
      count++;
    }
  }

  fList->Merge(&collections);
  
  delete iter;

  return count+1;
}


AliHistoListWrapper& AliHistoListWrapper::operator=(const AliHistoListWrapper& wrap) {

  // Assignment operator
  if(this!=&wrap) {
    
    fList = new TList();
    fList->SetOwner();
    TIterator* iter = wrap.fList->MakeIterator();
    TObject* obj;

    while ((obj = iter->Next())) {
      fList->Add(obj->Clone());
    }

  }
  return *this;
}

