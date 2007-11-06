//insert copyright

#include "AliHLTPHOSDigit.h"
#include "AliHLTPHOSAltroConfig.h"

ClassImp(AliHLTPHOSDigit);

AliHLTPHOSDigit::AliHLTPHOSDigit() :
  TObject(),
  AliHLTPHOSBase(),
  fX(-1),
  fZ(-1),
  fAmplitude(-1),
  fTime(-1),
  fEnergy(-1),
  fGain(-1),
  fSamples(55),
  fPreSamples(15),
  fTotalSamples(70),
  fDebugVar(-1)
{
  //added by PT
  fSamples = fNSamples;
  fPreSamples = fNPresamples;
  fTotalSamples = fNTotalSamples;
  //   fData = new Int_t[fNTotalSamples];
  fData = new Int_t[fNTotalSamples];

}

AliHLTPHOSDigit::~AliHLTPHOSDigit()
{
  //comment
}

void 
AliHLTPHOSDigit::SetRawData(Int_t *dataPtr)
{
  //modified by PT
  //  for(Int_t i = 0; i < 70; i++)
  //    {
  //     fData[i] = dataPtr[i];
  //    }
  for(Int_t i = 0; i < fNTotalSamples; i++)
    {
      fData[i] = dataPtr[i];
    }
}
 

void 
AliHLTPHOSDigit::ResetDigit()
  {
    fZ = -1;
    fX = -1;
    fAmplitude = -1;
    fTime = -1;
    fEnergy =-1;
    fGain = -1;
    fSamples = 55;
    fPreSamples =15;
    fTotalSamples =70;
    fDebugVar = -1;
  }
