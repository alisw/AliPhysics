// $Header$

#include "AliEVEHOMERSourceList.h"

using namespace Reve;
//using namespace Alieve;

//______________________________________________________________________
// AliEVEHOMERSourceList
//

ClassImp(AliEVEHOMERSourceList)

AliEVEHOMERSourceList::AliEVEHOMERSourceList(const Text_t* n, const Text_t* t) :
  Reve::RenderElementList(n, t)
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
