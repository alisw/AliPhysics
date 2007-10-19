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
// Default costructor
}


AliITSChannelDaSSD::AliITSChannelDaSSD(const UShort_t stripID) :
  fStripId(stripID),
  fEventsNumber(0),
  fSignal(NULL),
  fPedestal(fgkUndefinedValue),
  fNoise(fgkUndefinedValue),
  fZsThresholdFactor(0.0f)
{
// Costructor, initialize channal id
}


AliITSChannelDaSSD::AliITSChannelDaSSD(const UShort_t stripID, const Long_t eventsnumber) :
  fStripId(stripID),
  fEventsNumber(0),
  fSignal(NULL),
  fPedestal(fgkUndefinedValue),
  fNoise(fgkUndefinedValue),
  fZsThresholdFactor(0.0f)
{
// Costructor, initialize channal id and allocate array for events data
  if (stripID > fgkMaxStripId)
    Warning("AliITSChannelDaSSD", "Wrong StripID: %i", stripID);
  fSignal = new (nothrow) Short_t[eventsnumber];
  if (fSignal) {
    fEventsNumber = eventsnumber;
    memset(fSignal, fgkDefaultSignal, (eventsnumber * sizeof(Short_t)));
  } else {
    Error("AliITSChannelDaSSD", "Error allocating memory for %i Short_t objects!", eventsnumber);
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
// Destructor
  if (fSignal) 
  {
     delete [] fSignal;
  }
}


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
    Error("AliITSChannelDaSSD", "Error allocating memory for %i Short_t objects!", eventsnumber);
    fSignal = NULL;
    fEventsNumber = 0;
    return kFALSE;
  }
}



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
