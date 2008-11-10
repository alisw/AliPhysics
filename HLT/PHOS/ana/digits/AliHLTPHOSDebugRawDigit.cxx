// $Id$

//insert copyright

#include "AliHLTPHOSDebugRawDigit.h"

ClassImp(AliHLTPHOSDebugRawDigit);

AliHLTPHOSDebugRawDigit::AliHLTPHOSDebugRawDigit() :
  TObject(),
  fX(-1),
  fZ(-1),
  fAmplitude(-1),
  fTime(-1),
  fEnergy(-1),
  fGain(-1)
  
    
{
}

AliHLTPHOSDebugRawDigit::~AliHLTPHOSDebugRawDigit()
{
}

void 
AliHLTPHOSDebugRawDigit::SetRawData(UInt_t *dataPtr)
{
  for(Int_t i = 0; i < 70; i++)
    {
      fData[i] = dataPtr[i];
    }
}
 
