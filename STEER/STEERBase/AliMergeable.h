#ifndef ALIMERGEABLE_H
#define ALIMERGEABLE_H
//pure virtual interface class to allow online merging
//Author: Mikolaj Krzewicki, mkrzewic@cern.ch
#include "TCollection.h"

class AliMergeable {
public:
  virtual ~AliMergeable() {}
  virtual Long64_t Merge(TCollection *list) = 0;
  virtual TCollection* GetListOfDrawableObjects() = 0;
  ClassDef(AliMergeable,0);
};

#endif
