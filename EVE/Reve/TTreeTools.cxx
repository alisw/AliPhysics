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
  if(ProcessCut(entry)) { ProcessFill(entry); return true; }
  return false;
}

Bool_t TSelectorToEventList::ProcessCut(Long64_t )
{
  return GetSelect()->EvalInstance(0) != 0;
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
