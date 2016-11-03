#ifndef AliHLTexampleMergeableContainer_h
#define AliHLTexampleMergeableContainer_h

#include "AliMergeable.h"
#include "AliHLTObjArray.h"
#include "TObject.h"

class TRootIoCtor;

class AliHLTexampleMergeableContainer : public TObject, public AliMergeable {

private:
  AliHLTObjArray* fContainer;

public:
  AliHLTexampleMergeableContainer(TRootIoCtor*): TObject(), AliMergeable(), fContainer(NULL) {
  }

  AliHLTexampleMergeableContainer(): TObject(), AliMergeable(),
    fContainer(new AliHLTObjArray(1))
  {
    //printf("AliHLTexampleMergeableContainer ctor, %p\n}
    fContainer->SetOwner(kTRUE);
  }

  AliHLTexampleMergeableContainer(const char* name): TObject(), AliMergeable(),
    fContainer(new AliHLTObjArray(1))
  {
    //printf("AliHLTexampleMergeableContainer ctor, %p\n}
    fContainer->SetName(name);
    fContainer->SetOwner(kTRUE);
  }

  virtual ~AliHLTexampleMergeableContainer() {
    //printf("AliHLTexampleMergeableContainer dtor, %p\n",this);
    if (fContainer) fContainer->Delete();
    delete fContainer; fContainer=NULL;
  }

void Add(TObject*);

virtual Long64_t Merge(TCollection *list);
  virtual TCollection* GetListOfDrawableObjects();

  ClassDef(AliHLTexampleMergeableContainer,1);
};
#endif
