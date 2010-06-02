#ifndef ALITOFCALIBESD_H
#define ALITOFCALIBESD_H


////////////////////////////////////////////////////////////////////
//                                                                //
// recalculate the TOF signal from the TOF raw signal             //
// using the updates in the OCDB                                  //
//                                                                //
////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliTOFChannelOnlineStatusArray;
class TObjArray;
class AliTOFDeltaBCOffset;
class AliTOFCTPLatency;
class AliTOFRunParams;
class AliESDEvent;

class AliTOFcalibESD :
public TObject
{

 public:

  AliTOFcalibESD(); // default constructor
  virtual ~AliTOFcalibESD(); // default destructor

  Bool_t Init(Int_t run = -1); // init
  void CalibrateESD(AliESDEvent *event); // calibrate ESD

  Float_t GetTimeZero() const {return fTimeZero;}
  Float_t GetTOFResolution() const {return fTOFResolution;}
 private:

  AliTOFcalibESD(const AliTOFcalibESD &); // copy constructor
  AliTOFcalibESD &operator=(const AliTOFcalibESD &); // operator=

  Bool_t fInitFlag; //! init flag
  AliTOFChannelOnlineStatusArray *fChannelStatusArray; //! channel status array
  TObjArray *fParOfflineArray;            //! par offline array
  AliTOFDeltaBCOffset *fDeltaBCOffsetObj; //! deltaBC offset object
  AliTOFCTPLatency *fCTPLatencyObj;       //! CTP latency object
  AliTOFRunParams *fRunParamsObj;         //! run params object
  Float_t fTimeZero;                      //! TOF Time0
  Float_t fTOFResolution;                 //! TOF Time0
  
  ClassDef(AliTOFcalibESD, 1);
};

#endif /* ALITOFCALIBESD_H */
