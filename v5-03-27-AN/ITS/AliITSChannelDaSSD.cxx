/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$  */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides storage container ITS SSD channel callibration data
/// used by DA. 
///
///////////////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include "AliITSChannelDaSSD.h"
#include "TString.h"
#include "AliLog.h"

ClassImp(AliITSChannelDaSSD)

using namespace std;

const Short_t AliITSChannelDaSSD::fgkMinStripId = 0;               // minimum strip id
const Short_t AliITSChannelDaSSD::fgkMaxStripId = 1535;            // maximum strip id

const Short_t  AliITSChannelDaSSD::fgkSignalOverflow  =  2047;      // ADC overflow value
const Short_t  AliITSChannelDaSSD::fgkSignalUnderflow = -2048;      // ADC underflow value
const UShort_t AliITSChannelDaSSD::fgkDefaultSignal   =  0x7F;      // initialization value for fNoise, fPedestal, fSignal[i]
const Float_t  AliITSChannelDaSSD::fgkUndefinedValue  =  32639.0f;  // = 0x7F7F


//______________________________________________________________________________
AliITSChannelDaSSD::AliITSChannelDaSSD() :
  fStripId(0),
  fEventsNumber(0),
  fSignal(NULL),
  fPedestal(fgkUndefinedValue),
  fNoise(fgkUndefinedValue),
  fNoiseCM(fgkUndefinedValue),
  fNOverflowEv(0)
{
// Default costructor
}


//______________________________________________________________________________
AliITSChannelDaSSD::AliITSChannelDaSSD(const UShort_t stripID) :
  fStripId(stripID),
  fEventsNumber(0),
  fSignal(NULL),
  fPedestal(fgkUndefinedValue),
  fNoise(fgkUndefinedValue),
  fNoiseCM(fgkUndefinedValue),
  fNOverflowEv(0)
{
// Costructor, initialize channal id
}


//______________________________________________________________________________
AliITSChannelDaSSD::AliITSChannelDaSSD(const UShort_t stripID, const Long_t eventsnumber) :
  fStripId(stripID),
  fEventsNumber(0),
  fSignal(NULL),
  fPedestal(fgkUndefinedValue),
  fNoise(fgkUndefinedValue),
  fNoiseCM(fgkUndefinedValue),
  fNOverflowEv(0)
{
// Costructor, initialize channal id and allocate array for events data
  if (stripID > fgkMaxStripId)
    AliWarning(Form("AliITSChannelDaSSD: Wrong StripID: %i", stripID));
  fSignal = new (nothrow) Short_t[eventsnumber];
  if (fSignal) {
    fEventsNumber = eventsnumber;
    memset(fSignal, fgkDefaultSignal, (eventsnumber * sizeof(Short_t)));
  } else {
    AliError(Form("AliITSChannelDaSSD: Error allocating memory for %ld Short_t objects!", eventsnumber));
    fSignal = NULL;
    fEventsNumber = 0;
  }
}



//______________________________________________________________________________
AliITSChannelDaSSD::AliITSChannelDaSSD(const AliITSChannelDaSSD& strip) :
  TObject(strip),
  fStripId(strip.fStripId),
  fEventsNumber(strip.fEventsNumber),
  fSignal(NULL),
  fPedestal(strip.fPedestal),
  fNoise(strip.fNoise),
  fNoiseCM(strip.fNoiseCM),
  fNOverflowEv(strip.fNOverflowEv)
{
  // copy constructor
  if ((strip.fEventsNumber > 0) && (strip.fSignal)) {
    fSignal = new (nothrow) Short_t[strip.fEventsNumber];
    if (fSignal) {
      memcpy(fSignal, strip.fSignal, (strip.fEventsNumber * sizeof(Short_t)));
    } else {
      AliError(Form("AliITSChannelDaSSD: Error allocating memory for %ld Short_t objects!", strip.fEventsNumber));
      fSignal = NULL;
      fEventsNumber = 0;
    }
  }  
}



//______________________________________________________________________________
AliITSChannelDaSSD& AliITSChannelDaSSD::operator = (const AliITSChannelDaSSD& strip)
{
// assignment operator
  if (this == &strip)  return *this;
  TObject::operator=(strip);  
  if (fSignal) { delete [] fSignal; fSignal = NULL; }
  fStripId = strip.fStripId;
  fEventsNumber = strip.fEventsNumber;
  fPedestal = strip.fPedestal;
  fNoise = strip.fNoise;
  fNoiseCM = strip.fNoiseCM;
  fNOverflowEv = strip.fNOverflowEv;
  if ((strip.fEventsNumber > 0) && (strip.fSignal)) fSignal = new (nothrow) Short_t[strip.fEventsNumber];
  else return *this;
  if (fSignal) {
    memcpy(fSignal, strip.fSignal, (strip.fEventsNumber * sizeof(Short_t)));
  } else {
    AliError(Form("AliITSChannelDaSSD: Error allocating memory for %ld Short_t objects!", strip.fEventsNumber));
    fSignal = NULL;
    fEventsNumber = 0;
  }
  return *this;
}


//______________________________________________________________________________
AliITSChannelDaSSD::~AliITSChannelDaSSD()
{
// Destructor
  if (fSignal) delete [] fSignal;
}


//______________________________________________________________________________
Bool_t AliITSChannelDaSSD::SetEvenetsNumber(const Long_t eventsnumber)
{
// Allocate array for events data
  if (fSignal) {delete [] fSignal; fSignal = NULL; }
  fSignal = new (nothrow) Short_t[eventsnumber];
  if (fSignal) {
    fEventsNumber = eventsnumber;
    memset(fSignal, fgkDefaultSignal, (eventsnumber * sizeof(Short_t)));
    return kTRUE;
  } else {
    AliError(Form("AliITSChannelDaSSD: Error allocating memory for %ld Short_t objects!", eventsnumber));
    fSignal = NULL;
    fEventsNumber = 0;
    return kFALSE;
  }
}


//______________________________________________________________________________
Bool_t AliITSChannelDaSSD::SetSignal(const Long_t eventnumber, const Short_t signal)
{
// put signal value to array 
  if (eventnumber < fEventsNumber && fSignal)
  {
     fSignal[eventnumber] = signal;
     return kTRUE;
  }
  return kFALSE;
}
