#ifndef ALIHLTCALOTRIGGERDATASTRUCT_H
#define ALIHLTCALOTRIGGERDATASTRUCT_H

#include "Rtypes.h"

struct AliHLTCaloTriggerDataStruct {
  UChar_t             fRow;
  UChar_t             fCol;
  Float_t             fAmplitude;
  Float_t             fTime;
  UChar_t             fNL0Times;
  UChar_t             fL0Times[10];
  UInt_t              fL1TimeSum;
  Int_t               fTriggerBits;
};

#endif
