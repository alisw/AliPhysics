#ifndef ALITOFCHANNELONLINESTATUS_H
#define ALITOFCHANNELONLINESTATUS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////
//  class for TOF Online calibration for noise run  //
//////////////////////////////////////////////////////

#include "TObject.h"

class AliTOFChannelOnlineStatus: public TObject {

public:
  AliTOFChannelOnlineStatus();
  AliTOFChannelOnlineStatus(UChar_t status);
  AliTOFChannelOnlineStatus(const AliTOFChannelOnlineStatus &chan);
  AliTOFChannelOnlineStatus& operator= (const AliTOFChannelOnlineStatus &chan);
  virtual ~AliTOFChannelOnlineStatus(){};
  UChar_t  GetStatus()   	const {return fStatus;}
  void	   SetStatus(UChar_t status) {fStatus=status;}

  enum{
    kTOFOnlineUnknown=0x0, kTOFOnlineOk=0x15, kTOFOnlineBad=0x2a, 
      kTOFHWOk=0x1, kTOFHWBad=0x2, kTOFHVUnknown=0x0,
      kTOFPulserOk=0x4, kTOFPulserBad=0x8, kTOFPulserUnknown=0x0,
      kTOFNoiseOk=0x10, kTOFNoiseBad=0x20, kTOFNoiseUnknown=0x0
      };
  
  enum{
    kRightShiftHW=0,kRightShiftPulser=2,kRightShiftNoise=4
      };

private:
  UChar_t  fStatus;     // Status of the channel

  ClassDef(AliTOFChannelOnlineStatus,1) // TOF Online Calibration object 
                                        // setting status
};

#endif
