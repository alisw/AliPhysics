#ifndef ALITOFCHANNELONLINESTATUSARRAY_H
#define ALITOFCHANNELONLINESTATUSARRAY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

////////////////////////////////////////////////
//  class for TOF Online calibration          //
//  to define the status of the channels      //
//  New object created, to use an array       //
//  instead of a TObjArray.                   //
//  Storing all the info coming from          // 
//  HW FEE map, pulser runs, and noise runs   //
//  in a single object (char).                // 
////////////////////////////////////////////////

#include "TObject.h"

class AliTOFChannelOnlineStatusArray: public TObject {

public:

  enum{
    kTOFOnlineUnknown=0x0, kTOFOnlineOk=0x15, kTOFOnlineBad=0x2a, 
	    kTOFHWUnknown=0x0, kTOFHWOk=0x1, kTOFHWBad=0x2, kTOFHWReset=0xfc, kTOFHW=0x3,
	    kTOFPulserUnknown=0x0, kTOFPulserOk=0x4, kTOFPulserBad=0x8, kTOFPulserReset=0xf3, kTOFPulser=0xc,
	    kTOFNoiseUnknown=0x0, kTOFNoiseOk=0x10, kTOFNoiseBad=0x20, kTOFNoiseReset=0xcf, kTOFNoise=0x30
      };

  AliTOFChannelOnlineStatusArray();
  AliTOFChannelOnlineStatusArray(Int_t size);
  AliTOFChannelOnlineStatusArray(const AliTOFChannelOnlineStatusArray &source);
  AliTOFChannelOnlineStatusArray& operator= (const AliTOFChannelOnlineStatusArray &source);
  ~AliTOFChannelOnlineStatusArray();
  Int_t GetSize() const {return fSize;}
  void  SetStatus(Int_t pos, UChar_t parr);
  void  SetHWStatus(Int_t pos, UChar_t parr);
  void  SetPulserStatus(Int_t pos, UChar_t parr);
  void  SetNoiseStatus(Int_t pos, UChar_t parr);
  void SetLatencyWindow(Int_t pos, Int_t parr);
  UChar_t  GetStatus(Int_t pos) const;
  UChar_t  GetHWStatus(Int_t pos) const;
  UChar_t  GetPulserStatus(Int_t pos) const;
  UChar_t  GetNoiseStatus(Int_t pos) const;
  Int_t GetLatencyWindow(Int_t pos) const;
  Bool_t HasLatencyWindow() {return fLatencyWindow != NULL;};

private:
  Int_t fSize;      // Size of the array of UChar_t
  UChar_t * fArray; //[fSize]
  Int_t *fLatencyWindow; //[fSize]

  ClassDef(AliTOFChannelOnlineStatusArray,2)    // TOF Sensor Online Calibration object
};

#endif
