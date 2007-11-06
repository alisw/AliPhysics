#ifndef ALIHLTALTRODATA_H
#define ALIHLTALTRODATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "Rtypes.h"
#include "AliHLTAltroBunch.h"

class AliHLTAltroData
{
public:
  AliHLTAltroData();
  ~AliHLTAltroData();
  bool NextBunch(AliHLTAltroBunch *altrobunch);

  int GetChannel();
  int GetChip();
  int GetCard();
  int GetBranch();
  void Reset();

  UInt_t *fData; //comment
  UInt_t *fBunchData; //comment
  int fDataSize; //comment
  int fWc; //comment
  int fHadd; //comment
  int fBunchCounter; //comment
  bool fIsComplete; //comment

};


#endif

