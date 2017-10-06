#include "AliHLTexampleMergeableContainer.h"
#include "TCollection.h"
#include "TList.h"
ClassImp(AliHLTexampleMergeableContainer)

Long64_t AliHLTexampleMergeableContainer::Merge(TCollection *list)
{
  if (!fContainer) return -1;
  TList* tmplist = new TList();
  TIter next(list);
  while (TObject* tmp = next()) {
    AliHLTexampleMergeableContainer* tmpC = dynamic_cast<AliHLTexampleMergeableContainer*>(tmp);
    if (tmpC) {
      tmplist->Add(tmpC->fContainer);
    }
  }
  return fContainer->Merge(tmplist);
}

TCollection* AliHLTexampleMergeableContainer::GetListOfDrawableObjects()
{
  //move the container out
  AliHLTObjArray* tmp = fContainer;
  fContainer=NULL;
  return tmp;
}

void AliHLTexampleMergeableContainer::Add(TObject* object)
{
  if (fContainer) fContainer->AddLast(object);
}
