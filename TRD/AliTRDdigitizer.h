#ifndef ALITRDDIGITIZER_H
#define ALITRDDIGITIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>

class TFile;
class TF1;

class AliTRD;
class AliTRDdigitsManager;
class AliTRDgeometry;

///////////////////////////////////////////////////////
//  Produces digits from the hits information        //
///////////////////////////////////////////////////////

class AliTRDdigitizer : public TNamed {

 public:

  AliTRDdigitizer();
  AliTRDdigitizer(const Text_t* name, const Text_t* title);
  AliTRDdigitizer(const AliTRDdigitizer &d);
  virtual ~AliTRDdigitizer();
  AliTRDdigitizer &operator=(const AliTRDdigitizer &d);

  virtual void         Copy(TObject &d);
  virtual void         Init();
  virtual Bool_t       Open(const Char_t *name, Int_t nEvent = 0);
  virtual Bool_t       MakeDigits();
  virtual void         ReInit();
  virtual Bool_t       SumSDigits();
  virtual Bool_t       WriteDigits();
  virtual Bool_t       InitDetector();

  virtual void         SetGasGain(Float_t gasgain)      { fGasGain        = gasgain;  };
  virtual void         SetNoise(Float_t noise)          { fNoise          = noise;    };
  virtual void         SetChipGain(Float_t chipgain)    { fChipGain       = chipgain; };
  virtual void         SetADCoutRange(Float_t range)    { fADCoutRange    = range;    };
  virtual void         SetADCinRange(Float_t range)     { fADCinRange     = range;    };
  virtual void         SetADCthreshold(Int_t thresh)    { fADCthreshold   = thresh;   };
  virtual void         SetDiffusion(Int_t diffOn = 1)   { fDiffusionOn    = diffOn;   };
  virtual void         SetElAttach(Int_t elOn = 1)      { fElAttachOn     = elOn;     };
  virtual void         SetElAttachProp(Float_t prop)    { fElAttachProp   = prop;     };
  virtual void         SetExB(Int_t exbOn = 1)          { fExBOn          = exbOn;    };
  virtual void         SetPadResponse(Int_t prfOn = 1)  { fPRFOn          = prfOn;    };
  virtual void         SetTimeResponse(Int_t trfOn = 1) { fTRFOn          = trfOn;   
                                                          ReInit();                   };
  virtual void         SetDriftVelocity(Float_t v)      { fDriftVelocity  = v;       
                                                          ReInit();                   };
  virtual void         SetPadCoupling(Float_t v)        { fPadCoupling    = v;        };
  virtual void         SetTimeCoupling(Float_t v)       { fTimeCoupling   = v;        };
  virtual void         SetCompress(Int_t c = 1)         { fCompress       = c;        };
  virtual void         SetVerbose(Int_t v = 1)          { fVerbose        = v;        };
  virtual void         SetSDigits(Int_t v = 1)          { fSDigits        = v;        };
  virtual void         SetEvent(Int_t v = 0)            { fEvent          = v;        };

  AliTRDdigitsManager *Digits() const                   { return fDigits;             };

          Float_t      GetGasGain() const               { return fGasGain;            };
          Float_t      GetNoise() const                 { return fNoise;              };
          Float_t      GetChipGain() const              { return fChipGain;           };
          Float_t      GetADCoutRange() const           { return fADCoutRange;        };
          Float_t      GetADCinRange() const            { return fADCinRange;         };
          Int_t        GetADCthreshold() const          { return fADCthreshold;       };
          Float_t      GetDiffusionT() const            { return fDiffusionT;         };
          Float_t      GetDiffusionL() const            { return fDiffusionL;         };
          Float_t      GetElAttachProp() const          { return fElAttachProp;       };
          Float_t      GetOmegaTau() const              { return fOmegaTau;           };
          Float_t      GetDriftVelocity() const         { return fDriftVelocity;      };
          Float_t      GetPadCoupling() const           { return fPadCoupling;        };
          Float_t      GetTimeCoupling() const          { return fTimeCoupling;       };
          Bool_t       GetCompress() const              { return fCompress;           };
          Bool_t       GetSDigits() const               { return fSDigits;            };
          Float_t      GetTimeBinWidth() const          { return fTimeBinWidth;       };
  virtual Float_t      GetDiffusionL(Float_t vd, Float_t b);
  virtual Float_t      GetDiffusionT(Float_t vd, Float_t b);
  virtual Float_t      GetOmegaTau(Float_t vd, Float_t b);

 protected:

  TFile               *fInputFile;       //! ALIROOT-filename
  AliTRDdigitsManager *fDigits;          //! TRD digits manager
  AliTRD              *fTRD;             //! TRD detector class
  AliTRDgeometry      *fGeo;             //! TRD geometry
  
  Int_t                fEvent;           //! Event number

  Float_t              fField;           // Magnetic field
  Float_t              fGasGain;         // Gas gain
  Float_t              fNoise;           // Electronics noise
  Float_t              fChipGain;        // Electronics gain
  Float_t              fADCoutRange;     // ADC output range (number of channels)
  Float_t              fADCinRange;      // ADC input range (input charge)
  Float_t              fSinRange;        // Input range for summable digits 
  Float_t              fSoutRange;       // Output range for summable digits 
  Int_t                fADCthreshold;    // ADC threshold in ADC channel
  Int_t                fDiffusionOn;     // Switch for the diffusion
  Float_t              fDiffusionT;      // Diffusion in transverse direction
  Float_t              fDiffusionL;      // Diffusion in longitudinal direction
  Int_t                fElAttachOn;      // Switch for the electron attachment
  Float_t              fElAttachProp;    // Propability for electron attachment (for 1m)
  Int_t                fExBOn;           // Switch for the ExB effects
  Float_t              fOmegaTau;        // Tangens of the Lorentz angle 
  Float_t              fLorentzFactor;   // Factor due to Lorentz force
  Int_t                fPRFOn;           // Switch for the pad response
  Float_t             *fPRFsmp;          //!Sampled pad response
  Int_t                fPRFbin;          // Number of bins for the PRF
  Float_t              fPRFlo;           // Lower boundary of the PRF
  Float_t              fPRFhi;           // Higher boundary of the PRF
  Float_t              fPRFwid;          // Bin width of the sampled PRF
  Int_t                fPRFpad;          // Distance to next pad in PRF
  Int_t                fTRFOn;           // Switch for the time response
  Float_t             *fTRFsmp;          //!Integrated time response
  Int_t                fTRFbin;          // Number of bins for the TRF
  Float_t              fTRFlo;           // Lower boundary of the TRF
  Float_t              fTRFhi;           // Higher boundary of the TRF
  Float_t              fTRFwid;          // Bin width of the integrated TRF
  Float_t              fDriftVelocity;   // Drift velocity (cm / mus)
  Float_t              fTimeBinWidth;    // Time bin width in ns
  Float_t              fPadCoupling;     // Pad coupling factor
  Float_t              fTimeCoupling;    // Time coupling factor (image charge of moving ions)
  Bool_t               fCompress;        // Switch to keep only compressed data in memory
  Int_t                fVerbose;         // Sets the verbose level
  Bool_t               fSDigits;         // Switch for the summable digits

 private:

  virtual Int_t        Diffusion(Float_t driftlength, Float_t *xyz);
  virtual Int_t        ExB(Float_t driftlength, Float_t *xyz);  
  virtual Int_t        PadResponse(Float_t signal, Float_t dist, Float_t *pad);
  virtual Float_t      TimeResponse(Float_t time);  
  virtual Bool_t       CheckDetector(Int_t plane, Int_t chamber, Int_t sector);
  virtual void         SamplePRF();
  virtual void         SampleTRF();

  ClassDef(AliTRDdigitizer,4)            // Produces TRD-Digits

};

#endif
