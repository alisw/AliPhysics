#ifndef AliTRDCALDCS_H
#define AliTRDCALDCS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCS.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS parameters                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TObjArray.h"

class TString;

class AliTRDCalDCSFEE;
class AliTRDCalDCSPTR;
class AliTRDCalDCSGTU;

class AliTRDCalDCS : public TNamed {

 public:

  AliTRDCalDCS();
  AliTRDCalDCS(const Text_t *name, const Text_t *title);
  AliTRDCalDCS(const AliTRDCalDCS &cd);
  AliTRDCalDCS &operator=(const AliTRDCalDCS &cd);
  virtual ~AliTRDCalDCS() { };
    
  void    SetNumberOfTimeBins(Int_t value)    { fNumberOfTimeBins    = value; }
  void    SetTailCancelationTau1(Int_t tau1)  { fTailCancelationTau1 = tau1;  }
  void    SetTailCancelationTau2(Int_t tau2)  { fTailCancelationTau2 = tau2;  }
  void    SetTailCancelationAmp(Int_t amp)    { fTailCancelationAmp  = amp;   }
  void    SetPedestal(Int_t ped)              { fPedestal            = ped;   }
  void    SetConfigID(TString id)             { fConfigID            = id;    }
  void    SetGainTableID(TString id)          { fGainTableID         = id;    }
  void    SetFEEArr(TObjArray *fa)            { fFEEArr              = fa;    }
  void    SetPTRArr(TObjArray *pa)            { fPTRArr              = pa;    }
  void    SetGTUArr(TObjArray *ga)            { fGTUArr              = ga;    }

  Int_t   GetNumberOfTimeBins() const         { return fNumberOfTimeBins;     }
  Int_t   GetTailCancelationTau1() const      { return fTailCancelationTau1;  }
  Int_t   GetTailCancelationTau2() const      { return fTailCancelationTau2;  }
  Int_t   GetTailCancelationAmp() const       { return fTailCancelationAmp;   }
  Int_t   GetPedestal() const                 { return fPedestal;             }
  TString GetConfigID() const                 { return fConfigID;             }
  TString GetGainTableID() const              { return fGainTableID;          }
  TObjArray*       GetFEEArr() const          { return fFEEArr;               }
  TObjArray*       GetPTRArr() const          { return fPTRArr;               }
  TObjArray*       GetGTUArr() const          { return fGTUArr;               }
  AliTRDCalDCSFEE* GetCalDCSFEEObj(Int_t det) 
  		  	    { return (AliTRDCalDCSFEE*)fFEEArr->At(det);      }
  AliTRDCalDCSPTR* GetCalDCSPTRObj(Int_t det) 
  			    { return (AliTRDCalDCSPTR*)fPTRArr->At(det);      }
  AliTRDCalDCSGTU* GetCalDCSGTUObj(Int_t det) 
           		    { return (AliTRDCalDCSGTU*)fGTUArr->At(det);      }

 protected:

  // global configuration parameters
  Int_t   fNumberOfTimeBins;       //  Number of timebins  
  Int_t   fTailCancelationTau1;    //  Tau1 of tail cancelation
  Int_t   fTailCancelationTau2;    //  Tau2 of tail cancelation
  Int_t   fTailCancelationAmp;     //  Amplitude of tail cancelation
  Int_t   fPedestal;               //  Pedestal
  TString fConfigID;               //  Configuration ID
  TString fGainTableID;            //  Gain table ID

  //individual configuration parameters
  TObjArray *fFEEArr; // config param of the individual chambers
  TObjArray *fPTRArr; // config param of the pretrigger
  TObjArray *fGTUArr; // config param of the GTU

  ClassDef(AliTRDCalDCS,1)         //  TRD calibration class for TRD DCS parameters

};
#endif
