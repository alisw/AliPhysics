#ifndef TRDv1_H
#define TRDv1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 1    //
////////////////////////////////////////////////////////
 
#include <TF1.h> 
#include "AliTRD.h"

// Energy spectrum of the delta-rays 
Double_t Ermilova(Double_t *x, Double_t *par);

class AliTRDv1 : public AliTRD {

public:
  AliTRDv1() {}
  AliTRDv1(const char *name, const char *title);
  virtual        ~AliTRDv1();
  virtual void    CreateGeometry();
  virtual void    CreateMaterials();
  virtual Int_t   IsVersion() const { return 1; };
  virtual void    StepManager();
  virtual void    SetSensPlane(Int_t iplane = 0);
  virtual void    SetSensChamber(Int_t ichamber = 0);
  virtual void    SetSensSector(Int_t isector = 0);
  virtual void    Init();
  virtual void    Hits2Digits(); 
  virtual void    Digits2Clusters();

  virtual void    SetGasGain(Float_t gasgain)         { fGasGain       = gasgain;  };
  virtual void    SetNoise(Float_t noise)             { fNoise         = noise;    };
  virtual void    SetChipGain(Float_t chipgain)       { fChipGain      = chipgain; };
  virtual void    SetADCoutRange(Float_t range)       { fADCoutRange   = range;    };
  virtual void    SetADCinRange(Float_t range)        { fADCinRange    = range;    };
  virtual void    SetADCthreshold(Int_t thresh)       { fADCthreshold  = thresh;   };
  virtual void    SetDiffusionT(Float_t diff)         { fDiffusionT    = diff;     };
  virtual void    SetDiffusionL(Float_t diff)         { fDiffusionL    = diff;     };

  virtual void    SetClusMaxThresh(Float_t thresh)    { fClusMaxThresh = thresh;   };
  virtual void    SetClusSigThresh(Float_t thresh)    { fClusSigThresh = thresh;   };
  virtual void    SetClusMethod(Int_t meth)           { fClusMethod    = meth;     };

  virtual Float_t GetGasGain()                        { return fGasGain;       };
  virtual Float_t GetNoise()                          { return fNoise;         };
  virtual Float_t GetChipGain()                       { return fChipGain;      };
  virtual Float_t GetADCoutRange()                    { return fADCoutRange;   };
  virtual Float_t GetADCinRange()                     { return fADCinRange;    };
  virtual Int_t   GetADCthreshold()                   { return fADCthreshold;  };
  virtual Float_t GetDiffusionT()                     { return fDiffusionT;    };
  virtual Float_t GetDiffusionL()                     { return fDiffusionL;    };

  virtual Float_t GetClusMaxThresh()                  { return fClusMaxThresh; };
  virtual Float_t GetClusSigThresh()                  { return fClusSigThresh; };
  virtual Int_t   GetClusMethod()                     { return fClusMethod;    };

protected:
  Int_t        fIdSens;                 // Sensitive volume identifier

  Int_t        fIdChamber1;             // Driftchamber volume identifier
  Int_t        fIdChamber2;             // 
  Int_t        fIdChamber3;             // 

  Int_t        fSensSelect;             // Switch to select only parts of the detector
  Int_t        fSensPlane;              // Sensitive detector plane
  Int_t        fSensChamber;            // Sensitive detector chamber
  Int_t        fSensSector;             // Sensitive detector sector

  Float_t      fGasGain;                // Gas gain
  Float_t      fNoise;                  // Electronics noise
  Float_t      fChipGain;               // Electronics gain
  Float_t      fADCoutRange;            // ADC output range (number of channels)
  Float_t      fADCinRange;             // ADC input range (input charge)
  Int_t        fADCthreshold;           // ADC threshold in ADC channel
  Float_t      fDiffusionT;             // Diffusion in transverse direction
  Float_t      fDiffusionL;             // Diffusion in logitudinal direction

  Float_t      fClusMaxThresh;          // Threshold value for cluster maximum
  Float_t      fClusSigThresh;          // Threshold value for cluster signal
  Int_t        fClusMethod;             // Clustering method

private:
  virtual Double_t BetheBloch(Double_t bg);
  virtual void     Diffusion(Float_t driftlength, Float_t *xyz);
  virtual Float_t  PadResponse(Float_t x);
  virtual void     Pads2XYZ(Float_t *pads, Float_t *pos);
  virtual Float_t  Unfold(Float_t eps, Float_t *padSignal);

  TF1         *fDeltaE;                 // Energy distribution of the delta-electrons
   
  ClassDef(AliTRDv1,1)                  // Transition Radiation Detector version 1 (slow simulator)

};

#endif
