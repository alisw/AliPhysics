#ifndef ALITRDDIGITIZER_H
#define ALITRDDIGITIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigitizer.h"

class TFile;
class TF1;

class AliRunDigitizer;

class AliTRD;
class AliTRDdigitsManager;
class AliTRDgeometry;

///////////////////////////////////////////////////////
//  Produces digits from the hits information        //
///////////////////////////////////////////////////////

class AliTRDdigitizer : public AliDigitizer {

 public:

  AliTRDdigitizer();
  AliTRDdigitizer(const Text_t* name, const Text_t* title);
  AliTRDdigitizer(AliRunDigitizer *manager, const Text_t* name, const Text_t* title);
  AliTRDdigitizer(AliRunDigitizer *manager);
  AliTRDdigitizer(const AliTRDdigitizer &d);
  virtual ~AliTRDdigitizer();
  AliTRDdigitizer &operator=(const AliTRDdigitizer &d);

  virtual void         Copy(TObject &d);
  virtual Bool_t       Init();
  virtual Bool_t       InitDetector();
  virtual Bool_t       ReInit();
  virtual void         Exec(Option_t* option = 0);  
  virtual Bool_t       Open(const Char_t *file, Int_t nEvent = 0);
  virtual Bool_t       MakeBranch(const Char_t *file = 0);
  virtual Bool_t       MakeDigits();
  virtual void         AddSDigitsManager(AliTRDdigitsManager *manager);
  virtual void         DeleteSDigitsManager();
  virtual Bool_t       ConvertSDigits();
  virtual Bool_t       MergeSDigits();
  virtual Bool_t       SDigits2Digits();
  virtual Bool_t       WriteDigits();

  virtual void         SetGasGain(Float_t gasgain)          { fGasGain        = gasgain;  };
  virtual void         SetNoise(Float_t noise)              { fNoise          = noise;    };
  virtual void         SetChipGain(Float_t chipgain)        { fChipGain       = chipgain; };
  virtual void         SetADCoutRange(Float_t range)        { fADCoutRange    = range;    };
  virtual void         SetADCinRange(Float_t range)         { fADCinRange     = range;    };
  virtual void         SetADCthreshold(Int_t thresh)        { fADCthreshold   = thresh;   };
  virtual void         SetDiffusion(Int_t diffOn = 1)       { fDiffusionOn    = diffOn;   };
  virtual void         SetElAttach(Int_t elOn = 1)          { fElAttachOn     = elOn;     };
  virtual void         SetElAttachProp(Float_t prop)        { fElAttachProp   = prop;     };
  virtual void         SetExB(Int_t exbOn = 1)              { fExBOn          = exbOn;    };
  virtual void         SetPadResponse(Int_t prfOn = 1)      { fPRFOn          = prfOn;    };
  virtual void         SetTimeResponse(Int_t trfOn = 1)     { fTRFOn          = trfOn;   
                                                              ReInit();                   };
  virtual void         SetCrossTalk(Int_t ctOn = 1)         { fCTOn           = ctOn;   
                                                              ReInit();                   };
  virtual void         SetTailCancelation(Int_t tcOn = 1)   { fTCOn           = tcOn;     };
  virtual void         SetNexponential(Int_t nexp)          { fTCnexp         = nexp;     };
  virtual void         SetDriftVelocity(Float_t v)          { fDriftVelocity  = v;       
                                                              ReInit();                   };
  virtual void         SetPadCoupling(Float_t v)            { fPadCoupling    = v;        };
  virtual void         SetTimeCoupling(Float_t v)           { fTimeCoupling   = v;        };
  virtual void         SetTiltingAngle(Float_t v);
  virtual void         SetCompress(Int_t c = 1)             { fCompress       = c;        };
  virtual void         SetDebug(Int_t v = 1)                { fDebug          = v;        };
  virtual void         SetSDigits(Int_t v = 1)              { fSDigits        = v;        };
  virtual void         SetSDigitsScale(Float_t s)           { fSDigitsScale   = s;        };
  virtual void         SetEvent(Int_t v = 0)                { fEvent          = v;        };
  virtual void         SetManager(AliTRDdigitsManager *man) { fDigitsManager  = man;      };    

  AliTRDdigitsManager *Digits() const                       { return fDigitsManager;      };

          Float_t      GetGasGain() const                   { return fGasGain;            };
          Float_t      GetNoise() const                     { return fNoise;              };
          Float_t      GetChipGain() const                  { return fChipGain;           };
          Float_t      GetADCoutRange() const               { return fADCoutRange;        };
          Float_t      GetADCinRange() const                { return fADCinRange;         };
          Int_t        GetADCthreshold() const              { return fADCthreshold;       };
          Float_t      GetDiffusionT() const                { return fDiffusionT;         };
          Float_t      GetDiffusionL() const                { return fDiffusionL;         };
          Float_t      GetElAttachProp() const              { return fElAttachProp;       };
          Int_t        GetExB() const                       { return fExBOn;              };
          Float_t      GetOmegaTau() const                  { return fOmegaTau;           };
          Float_t      GetDriftVelocity() const             { return fDriftVelocity;      };
          Float_t      GetPadCoupling() const               { return fPadCoupling;        };
          Float_t      GetTimeCoupling() const              { return fTimeCoupling;       };
          Bool_t       GetCompress() const                  { return fCompress;           };
          Bool_t       GetSDigits() const                   { return fSDigits;            };
          Float_t      GetSDigitsScale() const              { return fSDigitsScale;       };
          Float_t      GetTimeBinWidth() const              { return fTimeBinWidth;       };
          Float_t      GetTiltingAngle() const;
  virtual Float_t      GetDiffusionL(Float_t vd, Float_t b);
  virtual Float_t      GetDiffusionT(Float_t vd, Float_t b);
  virtual Float_t      GetOmegaTau(Float_t vd, Float_t b);

 protected:

  TFile               *fInputFile;          //! ALIROOT-file
  AliTRDdigitsManager *fDigitsManager;      //! Manager for the output digits
  AliTRDdigitsManager *fSDigitsManager;     //! Manager for the summed input s-digits
  TList               *fSDigitsManagerList; //! List of managers of input s-digits
  AliTRD              *fTRD;                //! TRD detector class
  AliTRDgeometry      *fGeo;                //! TRD geometry
  
  Int_t                fEvent;              //! Event number

  Int_t               *fMasks;              //! Masks for the merging

  Float_t              fField;              //  Magnetic field
  Float_t              fGasGain;            //  Gas gain
  Float_t              fNoise;              //  Electronics noise
  Float_t              fChipGain;           //  Electronics gain
  Float_t              fADCoutRange;        //  ADC output range (number of channels)
  Float_t              fADCinRange;         //  ADC input range (input charge)
  Int_t                fADCthreshold;       //  ADC threshold in ADC channel
  Int_t                fDiffusionOn;        //  Switch for the diffusion
  Float_t              fDiffusionT;         //  Diffusion in transverse direction
  Float_t              fDiffusionL;         //  Diffusion in longitudinal direction
  Int_t                fElAttachOn;         //  Switch for the electron attachment
  Float_t              fElAttachProp;       //  Propability for electron attachment (for 1m)
  Int_t                fExBOn;              //  Switch for the ExB effects
  Float_t              fOmegaTau;           //  Tangens of the Lorentz angle 
  Float_t              fLorentzFactor;      //  Factor due to Lorentz force
  Int_t                fPRFOn;              //  Switch for the pad response
  Float_t             *fPRFsmp;             //! Sampled pad response
  Int_t                fPRFbin;             //  Number of bins for the PRF
  Float_t              fPRFlo;              //  Lower boundary of the PRF
  Float_t              fPRFhi;              //  Higher boundary of the PRF
  Float_t              fPRFwid;             //  Bin width of the sampled PRF
  Int_t                fPRFpad;             //  Distance to next pad in PRF
  Int_t                fTRFOn;              //  Switch for the time response
  Float_t             *fTRFsmp;             //! Integrated time response
  Int_t                fTRFbin;             //  Number of bins for the TRF
  Float_t              fTRFlo;              //  Lower boundary of the TRF
  Float_t              fTRFhi;              //  Higher boundary of the TRF
  Float_t              fTRFwid;             //  Bin width of the integrated TRF
  Int_t                fCTOn;               //  Switch for cross talk
  Float_t             *fCTsmp;              //! Integrated cross talk
  Int_t                fTCOn;               //  Switch for the tail cancelation
  Int_t                fTCnexp;             //  Number of exponential of the digital filter
  Float_t              fDriftVelocity;      //  Drift velocity (cm / mus)
  Float_t              fTimeBinWidth;       //  Time bin width in ns
  Float_t              fPadCoupling;        //  Pad coupling factor
  Float_t              fTimeCoupling;       //  Time coupling factor (image charge of moving ions)
  Float_t              fTiltingAngle;       //  Tilting angle of the readout pads
  Bool_t               fCompress;           //  Switch to keep only compressed data in memory
  Int_t                fDebug;              //  Sets the debug level
  Bool_t               fSDigits;            //  Switch for the summable digits
  Float_t              fSDigitsScale;       //  Scale factor for the summable digits 

 private:

  virtual Float_t      Col0Tilted(Float_t col0, Float_t rowOffset, Int_t plane);
  virtual Float_t      CrossTalk(Float_t time); 
  virtual Int_t        Diffusion(Float_t driftlength, Float_t *xyz);
  virtual Int_t        ExB(Float_t driftlength, Float_t *xyz);  
  virtual Int_t        PadResponse(Float_t signal, Float_t dist, Int_t plane, Float_t *pad);
  virtual Float_t      TimeResponse(Float_t time);  
  virtual void         DeConvExp(Double_t *source, Double_t *target, Int_t n, Int_t nexp);
  virtual Bool_t       CheckDetector(Int_t plane, Int_t chamber, Int_t sector);
  virtual void         SamplePRF();
  virtual void         SampleTRF();

  ClassDef(AliTRDdigitizer,6)               //  Produces TRD-Digits

};

#endif
