#ifndef AliTRDCALDCSFEE_H
#define AliTRDCALDCSFEE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSFEE.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD FEE configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSFEE : public TNamed {

 public:

  AliTRDCalDCSFEE();
  AliTRDCalDCSFEE(const char *name, const char *title);
  virtual ~AliTRDCalDCSFEE() { };

  void    SetDCSid(Int_t dcsid)              { fDCSID               = dcsid; }  
  void    SetSM(Int_t smid)                  { fSM                  = smid;  }
  void    SetStack(Int_t stid)               { fStack               = stid;  }
  void    SetLayer(Int_t lyid)               { fLayer               = lyid;  }
  void    SetNumberOfTimeBins(Int_t value)   { fNumberOfTimeBins    = value; }
  void    SetTailCancelationTau1(Int_t tau1) { fTailCancelationTau1 = tau1;  }
  void    SetTailCancelationTau2(Int_t tau2) { fTailCancelationTau2 = tau2;  }
  void    SetTailCancelationAmp(Int_t amp)   { fTailCancelationAmp  = amp;   }
  void    SetPedestal(Int_t ped)             { fPedestal            = ped;   }
  void    SetConfigID(TString id)            { fConfigID            = id;    }
  void    SetGainTableID(TString id)         { fGainTableID         = id;    }

  Int_t   GetDCSid() const                   { return fDCSID;                }
  Int_t   GetSM() const                      { return fSM;                   }
  Int_t   GetStack() const                   { return fStack;                }
  Int_t   GetLayer() const                   { return fLayer;                }
  Int_t   GetNumberOfTimeBins() const        { return fNumberOfTimeBins;     }
  Int_t   GetTailCancelationTau1() const     { return fTailCancelationTau1;  }
  Int_t   GetTailCancelationTau2() const     { return fTailCancelationTau2;  }
  Int_t   GetTailCancelationAmp() const      { return fTailCancelationAmp;   }
  Int_t   GetPedestal() const                { return fPedestal;             }
  TString GetConfigID() const                { return fConfigID;             }
  TString GetGainTableID() const             { return fGainTableID;          }

 protected:
  
  Int_t   fDCSID;                  //  ID of the DCS-Board
  Int_t   fSM;                     //  the number of the supermode 0..17
  Int_t   fStack;                  //  the number of the stack 0..4
  Int_t   fLayer;                  //  the number of the layer 0..5
  Int_t   fNumberOfTimeBins;       //  Number of timebins  

  Int_t   fTailCancelationTau1;    //  Tau1 of tail cancelation
  Int_t   fTailCancelationTau2;    //  Tau2 of tail cancelation
  Int_t   fTailCancelationAmp;     //  Amplitude of tail cancelation

  Int_t   fPedestal;               //  Pedestal

  TString fConfigID;               //  Configuration ID
  TString fGainTableID;            //  Gain table ID

  ClassDef(AliTRDCalDCSFEE,1)         //  TRD calibration class for TRD FEE parameters

};
#endif
