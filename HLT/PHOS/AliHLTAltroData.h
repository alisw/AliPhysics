#ifndef ALIHLTALTRODATA_H
#define ALIHLTALTRODATA_H

#include "Rtypes.h"
#include "AliHLTAltroBunch.h"

class AliHLTAltroData
{
public:
  AliHLTAltroData();
  ~ AliHLTAltroData();
  bool NextBunch(AliHLTAltroBunch &altrobunch);

  int GetChannel();
  int GetChip();
  int GetCard();
  int GetBranch();
  void Reset();

  UInt_t *fData;
  UInt_t *fBunchData;
  int fDataSize;
  int fWc;
  int fHadd;
  int fBunchCounter;
  bool fIsComplete;

};


#endif

