#ifndef ALIHLTALTROBUNCH_H
#define ALIHLTALTROBUNCH_H


#include "Rtypes.h"
#include <iostream>

using namespace std;


class AliHLTAltroBunch
{
public:
  AliHLTAltroBunch();
  ~ AliHLTAltroBunch();
  //  unsigned int fData[1024];
  //  unsigned int fBunchSize;
  //  UInt_t fData[1024];
  UInt_t *fData;
  //  UInt_t fBunchSize;
 Int_t fBunchSize;

  unsigned int fEndTimeBin;

 private:
  unsigned int fStartTimeBin;
  //  int fBunchPos;
};




#endif

