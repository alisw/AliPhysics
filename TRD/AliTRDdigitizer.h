#ifndef TRDdigitizer_h
#define TRDdigitizer_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>
#include <TFile.h>
#include <TF1.h>

#include "AliHit.h" 
#include "AliTRDdigit.h"
#include "AliTRDconst.h"
#include "AliTRDgeometry.h"

class AliTRDdigitsManager;

///////////////////////////////////////////////////////
//  Produces digits from the hits information        //
///////////////////////////////////////////////////////

class AliTRDdigitizer : public TNamed {

 public:

  AliTRDdigitizer();
  AliTRDdigitizer(const Text_t* name, const Text_t* title);
  ~AliTRDdigitizer();

  virtual void         Init();
  virtual Bool_t       Open(const Char_t *name, Int_t nEvent = 0);
  virtual Bool_t       MakeDigits();
  virtual Bool_t       WriteDigits();

  virtual void         SetGasGain(Float_t gasgain)     { fGasGain       = gasgain;  };
  virtual void         SetNoise(Float_t noise)         { fNoise         = noise;    };
  virtual void         SetChipGain(Float_t chipgain)   { fChipGain      = chipgain; };
  virtual void         SetADCoutRange(Float_t range)   { fADCoutRange   = range;    };
  virtual void         SetADCinRange(Float_t range)    { fADCinRange    = range;    };
  virtual void         SetADCthreshold(Int_t thresh)   { fADCthreshold  = thresh;   };
  virtual void         SetDiffusion(Int_t diff_on = 1) { fDiffusionOn   = diff_on;  };
  virtual void         SetDiffusionT(Float_t diff)     { fDiffusionT    = diff;     };
  virtual void         SetDiffusionL(Float_t diff)     { fDiffusionL    = diff;     };
  virtual void         SetElAttach(Int_t el_on = 1)    { fElAttachOn    = el_on;    };
  virtual void         SetElAttachProp(Float_t prop)   { fElAttachProp  = prop;     };
  virtual void         SetExB(Int_t exb_on = 1)        { fExBOn         = exb_on;   };
  virtual void         SetLorentzAngle(Float_t angle)  { fLorentzAngle  = angle;    };
  virtual void         SetPadResponse(TF1 *PRF)        { if (fPRF) delete fPRF;
                                                         fPRF           = PRF;      };

  AliTRDdigitsManager *Digits()                        { return fDigits;        };

  virtual Float_t      GetGasGain()                    { return fGasGain;       };
  virtual Float_t      GetNoise()                      { return fNoise;         };
  virtual Float_t      GetChipGain()                   { return fChipGain;      };
  virtual Float_t      GetADCoutRange()                { return fADCoutRange;   };
  virtual Float_t      GetADCinRange()                 { return fADCinRange;    };
  virtual Int_t        GetADCthreshold()               { return fADCthreshold;  };
  virtual Float_t      GetDiffusionT()                 { return fDiffusionT;    };
  virtual Float_t      GetDiffusionL()                 { return fDiffusionL;    };
  virtual Float_t      GetElAttachProp()               { return fElAttachProp;  };
  virtual Float_t      GetLorentzAngle()               { return fLorentzAngle;  };
  virtual TF1         *GetPadResponse()                { return fPRF;           };

 protected:

  TFile               *fInputFile;       //! ALIROOT-filename
  AliTRDdigitsManager *fDigits;          //! TRD digits manager
  AliTRD              *fTRD;             //! TRD detector class
  AliTRDgeometry      *fGeo;             //! TRD geometry
  
  Int_t                fEvent;           //! Event number

  Float_t              fGasGain;         // Gas gain
  Float_t              fNoise;           // Electronics noise
  Float_t              fChipGain;        // Electronics gain
  Float_t              fADCoutRange;     // ADC output range (number of channels)
  Float_t              fADCinRange;      // ADC input range (input charge)
  Int_t                fADCthreshold;    // ADC threshold in ADC channel
  Int_t                fDiffusionOn;     // Switch for the diffusion
  Float_t              fDiffusionT;      // Diffusion in transverse direction
  Float_t              fDiffusionL;      // Diffusion in longitudinal direction
  Int_t                fElAttachOn;      // Switch for the electron attachment
  Float_t              fElAttachProp;    // Propability for electron attachment (for 1m)
  Int_t                fExBOn;           // Switch for the ExB effects
  Float_t              fLorentzAngle;    // Lorentz angle 
  Float_t              fLorentzFactor;   // Factor due to Lorentz force
  TF1                 *fPRF;             // Pad response function

 private:

  virtual Int_t       Diffusion(Float_t driftlength, Float_t *xyz);
  virtual Int_t       ExB(Float_t driftlength, Float_t *xyz);  
  
  ClassDef(AliTRDdigitizer,1)            // Produces TRD-Digits

};

#endif
