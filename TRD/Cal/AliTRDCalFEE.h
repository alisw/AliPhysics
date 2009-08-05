#ifndef ALITRDCALFEE_H
#define ALITRDCALFEE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD FEE parameters                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalFEE : public TNamed {

 public:

  AliTRDCalFEE();
  AliTRDCalFEE(const Text_t *name, const Text_t *title);
  virtual ~AliTRDCalFEE() { };
    
  void    SetNumberOfTimeBins(Int_t value)   { fNumberOfTimeBins    = value; }
  void    SetTailCancelationTau1(Int_t tau1) { fTailCancelationTau1 = tau1;  }
  void    SetTailCancelationTau2(Int_t tau2) { fTailCancelationTau2 = tau2;  }
  void    SetTailCancelationAmp(Int_t amp)   { fTailCancelationAmp  = amp;   }
  void    SetPedestal(Int_t ped)             { fPedestal            = ped;   }
  void    SetConfigID(TString id)            { fConfigID            = id;    }
  void    SetGainTableID(TString id)         { fGainTableID         = id;    }

  Int_t   GetNumberOfTimeBins() const        { return fNumberOfTimeBins;     }
  Int_t   GetTailCancelationTau1() const     { return fTailCancelationTau1;  }
  Int_t   GetTailCancelationTau2() const     { return fTailCancelationTau2;  }
  Int_t   GetTailCancelationAmp() const      { return fTailCancelationAmp;   }
  Int_t   GetPedestal() const                { return fPedestal;             }
  TString GetConfigID() const                { return fConfigID;             }
  TString GetGainTableID() const             { return fGainTableID;          }

 protected:

  Int_t   fNumberOfTimeBins;       //  Number of timebins  

  Int_t   fTailCancelationTau1;    //  Tau1 of tail cancelation
  Int_t   fTailCancelationTau2;    //  Tau2 of tail cancelation
  Int_t   fTailCancelationAmp;     //  Amplitude of tail cancelation

  Int_t   fPedestal;               //  Pedestal

  TString fConfigID;               //  Configuration ID
  TString fGainTableID;            //  Gain table ID

  ClassDef(AliTRDCalFEE,1)         //  TRD calibration class for TRD FEE parameters

};
#endif
