// $Header$

#include "AliEVEHOMERSourceList.h"

//______________________________________________________________________
// AliEVEHOMERSourceList
//

ClassImp(AliEVEHOMERSourceList)

AliEVEHOMERSourceList::AliEVEHOMERSourceList(const Text_t* n, const Text_t* t) :
  TEveElementList(n, t)
{

}

void AliEVEHOMERSourceList::SelectAll()
{
  EnableListElements(kTRUE, kTRUE);
}

void AliEVEHOMERSourceList::DeselectAll()
{
  DisableListElements (kFALSE, kFALSE);
}
