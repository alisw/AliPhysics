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
  fSelection (sel),
  fSubIdExp  (),
  fSubIdNum  (0)
{
  SetInputList(&fInput);
}

Long64_t TPointSelector::Select(const Text_t* selection)
{
  TString var(fVarexp);
  if (fSubIdExp.IsNull()) {
    fSubIdNum = 0;
  } else {
    fSubIdNum = fSubIdExp.CountChar(':') + 1;
    var += ":" + fSubIdExp;
  }

  TString sel;
  if (selection != 0)
    sel = selection;
  else
    sel = fSelection;

  fInput.Delete();
  fInput.Add(new TNamed("varexp",    var.Data()));
  fInput.Add(new TNamed("selection", sel.Data()));

  if (fConsumer)
    fConsumer->InitFill(fSubIdNum);

  // 'para' option -> hack allowing arbitrary dimensions.
  if(fTree)
    fTree->Process(this, "goff para");

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
  if (fConsumer) {
    fConsumer->TakeAction(this);
  }
}
