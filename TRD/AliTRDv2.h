#ifndef TRDv2_H
#define TRDv2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 2    //
////////////////////////////////////////////////////////
 
#include <TF1.h> 
#include "AliTRD.h"

// Energy spectrum of the delta-rays 
Double_t Ermilova(Double_t *x, Double_t *par);

class AliTRDv2 : public AliTRD {

public:
  AliTRDv2() {}
  AliTRDv2(const char *name, const char *title);
  virtual        ~AliTRDv2();
  virtual void    CreateGeometry();
  virtual void    CreateMaterials();
  virtual Int_t   IsVersion() const {return 2;}
  virtual void    MakeBranch(Option_t* option);
  virtual void    StepManager();
  virtual void    SetSensPlane(Int_t iplane = 0);
  virtual void    SetSensChamber(Int_t ichamber = 0);
  virtual void    SetSensSector(Int_t isector = 0);
  virtual void    Init();
  virtual void    Hits2Digits(); 

  virtual void    SetRowPadSize(Float_t size)         { fRowPadSize   = size;     };
  virtual void    SetColPadSize(Float_t size)         { fColPadSize   = size;     };
  virtual void    SetTimeBinSize(Float_t size)        { fTimeBinSize  = size;     };

  virtual void    SetGasGain(Float_t gasgain)         { fGasGain      = gasgain;  };
  virtual void    SetNoise(Float_t noise)             { fNoise        = noise;    };
  virtual void    SetChipGain(Float_t chipgain)       { fChipGain     = chipgain; };
  virtual void    SetADCoutRange(Float_t range)       { fADCoutRange  = range;    };
  virtual void    SetADCinRange(Float_t range)        { fADCinRange   = range;    };
  virtual void    SetADCthreshold(Int_t thresh)       { fADCthreshold = thresh;   };
  virtual void    SetDiffusionT(Float_t diff)         { fDiffusionT   = diff;     };
  virtual void    SetDiffusionL(Float_t diff)         { fDiffusionL   = diff;     };

  virtual Float_t GetRowPadSize()                     { return fRowPadSize;   };
  virtual Float_t GetColPadSize()                     { return fColPadSize;   };
  virtual Float_t GetTimeBinSize()                    { return fTimeBinSize;  };

  virtual Float_t GetGasGain()                        { return fGasGain;      };
  virtual Float_t GetNoise()                          { return fNoise;        };
  virtual Float_t GetChipGain()                       { return fChipGain;     };
  virtual Float_t GetADCoutRange()                    { return fADCoutRange;  };
  virtual Float_t GetADCinRange()                     { return fADCinRange;   };
  virtual Int_t   GetADCthreshold()                   { return fADCthreshold; };
  virtual Float_t GetDiffusionT()                     { return fDiffusionT;   };
  virtual Float_t GetDiffusionL()                     { return fDiffusionL;   };

  virtual Int_t   GetRowMax(Int_t iplan, Int_t icham) { return fRowMax[iplan-1][icham-1]; };
  virtual Int_t   GetColMax(Int_t iplan)              { return fColMax[iplan-1];          };
  virtual Int_t   GetTimeMax()                        { return fTimeMax;                  };

protected:
  Int_t        fIdSens;                 // Sensitive volume identifier

  Int_t        fIdSpace1;               // Spaceframe volume identifier
  Int_t        fIdSpace2;               // 
  Int_t        fIdSpace3;               // 

  Int_t        fIdChamber1;             // Driftchamber volume identifier
  Int_t        fIdChamber2;             // 
  Int_t        fIdChamber3;             // 

  Int_t        fSensSelect;             // Switch to select only parts of the detector
  Int_t        fSensPlane;              // Sensitive detector plane
  Int_t        fSensChamber;            // Sensitive detector chamber
  Int_t        fSensSector;             // Sensitive detector sector

  Int_t        fRowMax[kNplan][kNcham]; // Number of pad-rows
  Int_t        fColMax[kNplan];         // Number of pad-columns
  Int_t        fTimeMax;                // Number of time buckets

  Float_t      fRowPadSize;             // Pad size in z-direction
  Float_t      fColPadSize;             // Pad size in rphi-direction
  Float_t      fTimeBinSize;            // Size of the time buckets

  Float_t      fGasGain;                // Gas gain
  Float_t      fNoise;                  // Electronics noise
  Float_t      fChipGain;               // Electronics gain
  Float_t      fADCoutRange;            // ADC output range (number of channels)
  Float_t      fADCinRange;             // ADC input range (input charge)
  Int_t        fADCthreshold;           // ADC threshold in ADC channel
  Float_t      fDiffusionT;             // Diffusion in transverse direction
  Float_t      fDiffusionL;             // Diffusion in logitudinal direction

private:
  virtual Double_t BetheBloch(Double_t bg);
  virtual void     Diffusion(Float_t driftlength, Float_t *xyz);
  virtual Float_t  PadResponse(Float_t x);

  TF1         *fDeltaE;                 // Energy distribution of the delta-electrons
   
  ClassDef(AliTRDv2,1)                  // Transition Radiation Detector version 2 (slow simulator)

};

#endif
