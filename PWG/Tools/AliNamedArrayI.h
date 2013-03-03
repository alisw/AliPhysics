#ifndef ALINAMEDARRAYI_H
#define ALINAMEDARRAYI_H

// $Id: AliNamedArrayI.h  $

#include <TArrayI.h>
#include <TNamed.h>

class AliNamedArrayI : public TNamed, public TArrayI {
 public: 
  AliNamedArrayI();
  AliNamedArrayI(const char *name, Int_t n);
  AliNamedArrayI(const char *name, Int_t n, const Int_t* array);
  AliNamedArrayI(const char *name, const TArrayI& array);

  void Clear(Option_t *option="");

private:
  AliNamedArrayI(const AliNamedArrayI&);             // not implemented
  AliNamedArrayI& operator=(const AliNamedArrayI&);  // not implemented
  
  ClassDef(AliNamedArrayI, 1); // Named integer array
};
#endif
