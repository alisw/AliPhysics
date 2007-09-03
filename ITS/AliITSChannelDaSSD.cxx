
#include "AliITSChannelDaSSD.h"

ClassImp(AliITSChannelDaSSD)

using namespace std;

const Float_t  AliITSChannelDaSSD::fgkUndefinedValue  = 32639.0f;  // = 0x7F7F

AliITSChannelDaSSD::AliITSChannelDaSSD() :
  fStripId(0),
  fEventsNumber(0),
  fSignal(NULL),
  fPedestal(fgkUndefinedValue),
  fNoise(fgkUndefinedValue),
  fZsThresholdFactor(0.0f)
{
}


AliITSChannelDaSSD::AliITSChannelDaSSD(const UShort_t stripID) :
  fStripId(stripID),
  fEventsNumber(0),
  fSignal(NULL),
  fPedestal(fgkUndefinedValue),
  fNoise(fgkUndefinedValue),
  fZsThresholdFactor(0.0f)
{
}


AliITSChannelDaSSD::AliITSChannelDaSSD(const UShort_t stripID, const Long_t eventsnumber) :
  fStripId(stripID),
  fEventsNumber(0),
  fSignal(NULL),
  fPedestal(fgkUndefinedValue),
  fNoise(fgkUndefinedValue),
  fZsThresholdFactor(0.0f)
{
  if (stripID > fgkMaxStripId)
    Warning("AliITSChannelDaSSD", "Wrong StripID: %i", stripID);
  try
  {
     fSignal = new Short_t [eventsnumber];
     fEventsNumber = eventsnumber;
     memset(fSignal, fgkDefaultSignal, (eventsnumber * sizeof(Short_t)));
  }
  catch (bad_alloc&) 
  {
     cout << "Error allocating memory for " << (long) eventsnumber << " Short_t objects in AliITSChannelDaSSD constructor!" << endl;
     fSignal = NULL;
     fEventsNumber = 0;
  }  
}



AliITSChannelDaSSD::AliITSChannelDaSSD(const AliITSChannelDaSSD& strip) :
  TObject(strip),
  fStripId(strip.fStripId),
  fEventsNumber(strip.fEventsNumber),
  fSignal(strip.fSignal),
  fPedestal(strip.fPedestal),
  fNoise(strip.fNoise),
  fZsThresholdFactor(strip.fZsThresholdFactor)
{
  // copy constructor

  Fatal("AliITSChannelDaSSD", "copy constructor not implemented");
}

AliITSChannelDaSSD& AliITSChannelDaSSD::operator = (const AliITSChannelDaSSD& strip)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}


AliITSChannelDaSSD::~AliITSChannelDaSSD()
{
  if (fSignal) 
  {
     delete [] fSignal;
  }
}


Bool_t AliITSChannelDaSSD::SetEvenetsNumber(const Long_t eventsnumber)
{
  try
  {
     fSignal = new Short_t[eventsnumber];
     fEventsNumber = eventsnumber;
     memset(fSignal, fgkDefaultSignal, (eventsnumber * sizeof(Short_t)));
     return kTRUE;
  }
  catch (bad_alloc&) 
  {
     cout << "Error allocating memory for " << (long) eventsnumber << " Short_t objects!" << endl;
     fSignal = NULL;
     fEventsNumber = 0;
     return kFALSE;
  }  
}



Bool_t AliITSChannelDaSSD::SetSignal(const Long_t eventnumber, const Short_t signal)
{
  if (eventnumber < fEventsNumber && fSignal)
  {
     fSignal[eventnumber] = signal;
     return kTRUE;
  }
  return kFALSE;
}
