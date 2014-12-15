#ifndef ALITOFCHANNELONLINEARRAY_H
#define ALITOFCHANNELONLINEARRAY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

////////////////////////////////////////////////
//  class for TOF Online calibration          //
//  to define the delay of the channels.      //
//  New object created, to use an array       //
//  instead of a TObjArray.                   // 
////////////////////////////////////////////////

#include "TObject.h"

class AliTOFChannelOnlineArray: public TObject {

public:

  enum{
    kTOFOnlineUnknown=0x0, kTOFOnlineOk=0x15, kTOFOnlineBad=0x2a, 
	    kTOFHWUnknown=0x0, kTOFHWOk=0x1, kTOFHWBad=0x2, kTOFHWReset=0xfc, kTOFHW=0x3,
	    kTOFPulserUnknown=0x0, kTOFPulserOk=0x4, kTOFPulserBad=0x8, kTOFPulserReset=0xf3, kTOFPulser=0xc,
	    kTOFNoiseUnknown=0x0, kTOFNoiseOk=0x10, kTOFNoiseBad=0x20, kTOFNoiseReset=0xcf, kTOFNoise=0x30
      };

  AliTOFChannelOnlineArray();
  AliTOFChannelOnlineArray(Int_t size);
  AliTOFChannelOnlineArray(const AliTOFChannelOnlineArray &source);
  AliTOFChannelOnlineArray& operator= (const AliTOFChannelOnlineArray &source);
  ~AliTOFChannelOnlineArray();
  Int_t GetSize() const {return fSize;}
  void  SetDelay(Int_t pos, Float_t parr);
  Float_t  GetDelay(Int_t pos) const;

private:
  Int_t fSize;      // Size of the array of UChar_t
  Float_t * fArray; //[fSize]

  ClassDef(AliTOFChannelOnlineArray,1)    // TOF Sensor Online Calibration object
};

#endif
