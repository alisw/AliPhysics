//------------------------------------------------------------------------------
// Implementation of the AliRecInfoCuts class. It keeps selection cuts for 
// reconstructed tracks. 
//
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

#include "AliRecInfoCuts.h"

ClassImp(AliRecInfoCuts)

//_____________________________________________________________________________
AliRecInfoCuts::AliRecInfoCuts(const Char_t* name,const Char_t *title) : AliESDtrackCuts(name, title)
, fMinTPCsignalN(0)
, fMaxAbsTanTheta(0)
{
  // init data members with defaults
  Init();
}

//_____________________________________________________________________________
void AliRecInfoCuts::Init()
{
  // set default values 
  SetMinTPCsignalN();
  SetMaxAbsTanTheta();
}


//_____________________________________________________________________________
/*
Long64_t AliRecInfoCuts::Merge(TCollection* list) const 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
  AliRecInfoCuts* entry = dynamic_cast<AliRecInfoCuts*>(obj);
  if (entry == 0) 
   continue;

  count++;
  }

return count;
}
*/

