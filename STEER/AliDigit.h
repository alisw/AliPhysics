#ifndef AliDigit_H
#define AliDigit_H
////////////////////////////////////////////////
//  Base class for Alice Digits               //
////////////////////////////////////////////////

#include "TObject.h"

class AliDigit : public TObject {
public:
  Int_t     fTracks[3];   //tracks number making this digit (up to 3)

public:
  AliDigit();
  AliDigit(Int_t *track);
  ~AliDigit() {;}
  inline virtual int *GetTracks() {return &fTracks[0];}
  
  ClassDef(AliDigit,1)  //Base class for all Alice digits
};
#endif
