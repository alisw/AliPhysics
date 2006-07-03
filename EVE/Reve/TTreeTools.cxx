// $Header$

//__________________________________________________________________________
// TTreeTools
//
// Collection of classes for TTree interaction.

#include "TTreeTools.h"
#include <TTree.h>
#include <TTreeFormula.h>

/**************************************************************************/
/**************************************************************************/

ClassImp(TSelectorToEventList)

TSelectorToEventList::TSelectorToEventList(TEventList* evl, const Text_t* sel) :
  TSelectorDraw(), fEvList(evl)
{
  fInput.Add(new TNamed("varexp", ""));
  fInput.Add(new TNamed("selection", sel));
  SetInputList(&fInput);
}

Bool_t TSelectorToEventList::Process(Long64_t entry)
{
  if(GetSelect()->EvalInstance(0) != 0)
    fEvList->Enter(entry);
  return kTRUE;
}

/**************************************************************************/
/**************************************************************************/

ClassImp(TTreeQuery)

Int_t TTreeQuery::Select(TTree* t, const Text_t* selection)
{
  TSelectorToEventList sel(this, selection);
  t->Process(&sel, "goff");
  return GetN();
}

/**************************************************************************/
// TPointSelectorConsumer, TPointSelector
/**************************************************************************/

ClassImp(TPointSelectorConsumer)
ClassImp(TPointSelector)

TPointSelector::TPointSelector(TTree* t,
			       TPointSelectorConsumer* c,
			       const Text_t* vexp, const Text_t* sel) :
  TSelectorDraw(),

  fTree      (t),
  fConsumer  (c),
  fVarexp    (vexp),
  fSelection (sel)
{
  SetInputList(&fInput);
}

Long64_t TPointSelector::Select(const Text_t* selection)
{
  if(selection != 0)
    fSelection = selection;

  fInput.Delete();
  fInput.Add(new TNamed("varexp",    fVarexp.Data()));
  fInput.Add(new TNamed("selection", fSelection.Data()));

  if(fTree)
    fTree->Process(this, "goff");
  return fSelectedRows;
}

Long64_t TPointSelector::Select(TTree* t, const Text_t* selection)
{
  fTree = t;
  return Select(selection);
}

void TPointSelector::TakeAction()
{
  fSelectedRows += fNfill;
  // printf("TPointSelector::TakeAction nfill=%d, nall=%lld\n", fNfill, fSelectedRows);
  if(fConsumer) {
    fConsumer->TakeAction(this);
  }
}
