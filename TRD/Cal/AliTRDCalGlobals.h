#ifndef AliTRDCALGLOBALS_H
#define AliTRDCALGLOBALS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for global TRD parameters                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalGlobals : public TNamed {

 public:

  AliTRDCalGlobals();
  AliTRDCalGlobals(const Text_t *name, const Text_t *title);
  virtual ~AliTRDCalGlobals() { };
    
  void    SetNumberOfTimeBins(Int_t value)   { fNumberOfTimeBins    = value; }
  void    SetTailCancelationTau1(Int_t tau1) { fTailCancelationTau1 = tau1;  }
  void    SetTailCancelationTau2(Int_t tau2) { fTailCancelationTau2 = tau2;  }
  void    SetTailCancelationAmp(Int_t amp)   { fTailCancelationAmp  = amp;   }
  void    SetPedestal(Int_t ped)             { fPedestal            = ped;   }
  void    SetADCClockphase(Float_t cp)       { fADCClockphase       = cp;    }
  void    SetT0Offset(Float_t t0)            { fT0Offset            = t0;    }
  void    SetConfigID(TString id)            { fConfigID            = id;    }
  void    SetGainTableID(TString id)         { fGainTableID         = id;    }
  void    SetPretriggerConf(TString conf)    { fPretriggerConf      = conf;  }

  Int_t   GetNumberOfTimeBins() const        { return fNumberOfTimeBins;     }
  Int_t   GetTailCancelationTau1() const     { return fTailCancelationTau1;  }
  Int_t   GetTailCancelationTau2() const     { return fTailCancelationTau2;  }
  Int_t   GetTailCancelationAmp() const      { return fTailCancelationAmp;   }
  Int_t   GetPedestal() const                { return fPedestal;             }
  Float_t GetADCClockphase() const           { return fADCClockphase;        }
  Float_t GetT0Offset() const                { return fT0Offset;             }
  TString GetConfigID() const                { return fConfigID;             }
  TString GetGainTableID() const             { return fGainTableID;          }
  TString GetPretriggerConf() const          { return fPretriggerConf;       }

 protected:

  Int_t   fNumberOfTimeBins;       //  Number of timebins  

  Int_t   fTailCancelationTau1;    //  Tau1 of tail cancelation
  Int_t   fTailCancelationTau2;    //  Tau2 of tail cancelation
  Int_t   fTailCancelationAmp;     //  Amplitude of tail cancelation

  Int_t   fPedestal;               //  Pedestal

  Float_t fADCClockphase;          //  ADC clockphase in respect to TTC

  Float_t fT0Offset;               //  Global offset on t0

  TString fConfigID;               //  Configuration ID
  TString fGainTableID;            //  Gain table ID
  TString fPretriggerConf;         //  Pretrigger configuration

  ClassDef(AliTRDCalGlobals,3)     //  TRD calibration class for global TRD parameters

};
#endif
