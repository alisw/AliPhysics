#include "AliHBTWeights.h"

ClassImp(AliHBTWeights)

AliHBTWeights* AliHBTWeights::fgWeights = 0x0;

AliHBTWeights::~AliHBTWeights()
 {
   delete fgWeights;
   fgWeights = 0x0;
 }
