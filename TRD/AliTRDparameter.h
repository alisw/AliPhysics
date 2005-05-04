#ifndef ALITRDPARAMETER_H
#define ALITRDPARAMETER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD parameter class                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TObjArray;
class AliTRDgeometry;
class AliTRDpadPlane;

class AliTRDparameter : public TNamed {

 public:

  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };
 
  AliTRDparameter();
  AliTRDparameter(const Text_t* name, const Text_t* title);
  AliTRDparameter(const AliTRDparameter &p);   
  virtual ~AliTRDparameter();
  AliTRDparameter &operator=(const AliTRDparameter &p); 

  virtual void     Copy(TObject &p) const;
  virtual void     Init();
  virtual void     ReInit();
 
  virtual void     SetSamplingFrequency(Float_t freq)             { fSamplingFrequency = freq;
                                                                    ReInit();                   };
  virtual void     SetDriftVelocity(Float_t vd)                   { fDriftVelocity     = vd;
                                                                    SampleTimeStruct();         };
  virtual void     SetExpandTimeBin(Int_t nbefore, Int_t nafter)
                                                                  { fTimeBefore = nbefore;
                                                                    fTimeAfter  = nafter;       };

  virtual void     SetGasGain(Float_t gasgain)                    { fGasGain        = gasgain;  };
  virtual void     SetNoise(Float_t noise)                        { fNoise          = noise;    };
  virtual void     SetChipGain(Float_t chipgain)                  { fChipGain       = chipgain; };
  virtual void     SetADCoutRange(Float_t range)                  { fADCoutRange    = range;    };
  virtual void     SetADCinRange(Float_t range)                   { fADCinRange     = range;    };
  virtual void     SetADCthreshold(Int_t thresh)                  { fADCthreshold   = thresh;   };
  virtual void     SetADCbaseline(Int_t basel)                    { fADCbaseline    = basel;    };   
  virtual void     SetDiffusion(Int_t diffOn = 1)                 { fDiffusionOn    = diffOn;   };
  virtual void     SetElAttach(Int_t elOn = 1)                    { fElAttachOn     = elOn;     };
  virtual void     SetElAttachProp(Float_t prop)                  { fElAttachProp   = prop;     };
  virtual void     SetExB(Int_t exbOn = 1)                        { fExBOn          = exbOn;    };
  virtual void     SetPadResponse(Int_t prfOn = 1)                { fPRFOn          = prfOn;    };
  virtual void     SetTimeResponse(Int_t trfOn = 1)               { fTRFOn          = trfOn;   
                                                                    ReInit();                   };
  virtual void     SetCrossTalk(Int_t ctOn = 1)                   { fCTOn           = ctOn;   
                                                                    ReInit();                   };
  virtual void     SetTimeStruct(Bool_t tsOn = 1)                 { fTimeStructOn   = tsOn;     };
  virtual void     SetTailCancelation(Int_t tcOn = 1)             { fTCOn           = tcOn;     };
  virtual void     SetNexponential(Int_t nexp)                    { fTCnexp         = nexp;     };
  virtual void     SetPadCoupling(Float_t v)                      { fPadCoupling    = v;        };
  virtual void     SetTimeCoupling(Float_t v)                     { fTimeCoupling   = v;        };

  virtual void     SetLUT(Int_t lutOn = 1)                        { fLUTOn          = lutOn;    };
  virtual void     SetClusMaxThresh(Int_t thresh)                 { fClusMaxThresh  = thresh;   };
  virtual void     SetClusSigThresh(Int_t thresh)                 { fClusSigThresh  = thresh;   };

  virtual void     SetAnodeWireOffset(Float_t offset = 0.25)      { fAnodeWireOffset = offset;};
  
          Int_t    GetTimeMax()                             const { return fTimeMax;           };
          Int_t    GetTimeBefore()                          const { return fTimeBefore;        }; 
          Int_t    GetTimeAfter()                           const { return fTimeAfter;         }; 
          Int_t    GetTimeTotal()                           const { return fTimeMax 
                                                                         + fTimeBefore 
                                                                         + fTimeAfter;         };
          Float_t  GetTime0(Int_t p)                        const { return fTime0[p];          };

          Float_t  GetGasGain()                             const { return fGasGain;           };
          Float_t  GetNoise()                               const { return fNoise;             };
          Float_t  GetChipGain()                            const { return fChipGain;          };
          Float_t  GetADCoutRange()                         const { return fADCoutRange;       };
          Float_t  GetADCinRange()                          const { return fADCinRange;        };
          Int_t    GetADCthreshold()                        const { return fADCthreshold;      };
          Int_t    GetADCbaseline()                         const { return fADCbaseline;       };
          Float_t  GetDiffusionT()                          const { return fDiffusionT;        };
          Float_t  GetDiffusionL()                          const { return fDiffusionL;        };
          Float_t  GetElAttachProp()                        const { return fElAttachProp;      };
          Float_t  GetOmegaTau()                            const { return fOmegaTau;          };
	  Float_t  GetSamplingFrequency()                   const { return fSamplingFrequency; };
          Float_t  GetDriftVelocity()                       const { return fDriftVelocity;     };
          Float_t  GetPadCoupling()                         const { return fPadCoupling;       };
          Float_t  GetTimeCoupling()                        const { return fTimeCoupling;      };
          Float_t  GetTRFlo()                               const { return fTRFlo;             };
          Float_t  GetTRFhi()                               const { return fTRFhi;             };
          Float_t  GetLorentzFactor()                       const { return fLorentzFactor;     };
          Float_t  GetAnodeWireOffset()                     const { return fAnodeWireOffset;   };
          Int_t    GetTCnexp()                              const { return fTCnexp;            };
  virtual Float_t  GetDiffusionL(Float_t vd, Float_t b);
  virtual Float_t  GetDiffusionT(Float_t vd, Float_t b);
  virtual Float_t  GetOmegaTau(Float_t vd, Float_t b);

  virtual Int_t    GetClusMaxThresh()                       const { return fClusMaxThresh; };
  virtual Int_t    GetClusSigThresh()                       const { return fClusSigThresh; };

  virtual AliTRDpadPlane *GetPadPlane(Int_t p, Int_t c) const;
          Int_t    GetRowMax(Int_t p, Int_t c, Int_t /*s*/) const;
          Int_t    GetColMax(Int_t p) const;
          Double_t GetRow0(Int_t p, Int_t c, Int_t /*s*/) const;
          Double_t GetCol0(Int_t p) const;

          void     PrintDriftVelocity();

          Bool_t   TimeStructOn()                           const { return fTimeStructOn;  };
          Bool_t   ExBOn()                                  const { return fExBOn;         };
          Bool_t   PRFOn()                                  const { return fPRFOn;         };
          Bool_t   TRFOn()                                  const { return fTRFOn;         };
          Bool_t   ElAttachOn()                             const { return fElAttachOn;    }; 
          Bool_t   DiffusionOn()                            const { return fDiffusionOn;   };
          Bool_t   CTOn()                                   const { return fCTOn;          };
          Bool_t   TCOn()                                   const { return fTCOn;          };
          Bool_t   LUTOn()                                  const { return fLUTOn;         };

  virtual Int_t     Diffusion(Double_t driftlength, Double_t *xyz);
  virtual Int_t     ExB(Double_t driftlength, Double_t *xyz) const;  
  virtual Int_t     PadResponse(Double_t signal, Double_t dist, Int_t plane, Double_t *pad) const;
  virtual Double_t  CrossTalk(Double_t time) const; 
  virtual Double_t  TimeResponse(Double_t time) const;  
  virtual Double_t  TimeStruct(Double_t time, Double_t z) const;  
  virtual Double_t  LUTposition(Int_t iplane, Double_t ampL, Double_t ampC, Double_t ampR) const;

 protected:

  AliTRDgeometry      *fGeo;                                //! TRD geometry       
  TObjArray           *fPadPlaneArray;                      //  Array of pad plane objects

  Int_t                fTimeMax;                            //  Number of timebins in the drift region
  Int_t                fTimeBefore;                         //  Number of timebins before the drift region
  Float_t              fTime0[kNplan];                      //  Time-position of pad 0
  Int_t                fTimeAfter;                          //  Number of timebins after the drift region

  // Digitization parameter
  Float_t              fField;                              //  Magnetic field
  Float_t              fGasGain;                            //  Gas gain
  Float_t              fNoise;                              //  Electronics noise
  Float_t              fChipGain;                           //  Electronics gain
  Float_t              fADCoutRange;                        //  ADC output range (number of channels)
  Float_t              fADCinRange;                         //  ADC input range (input charge)
  Int_t                fADCthreshold;                       //  ADC threshold in ADC channel
  Int_t                fADCbaseline;                        //  ADC baseline in ADC chann
  Int_t                fDiffusionOn;                        //  Switch for the diffusion
  Float_t              fDiffusionT;                         //  Diffusion in transverse direction
  Float_t              fDiffusionL;                         //  Diffusion in longitudinal direction
  Int_t                fElAttachOn;                         //  Switch for the electron attachment
  Float_t              fElAttachProp;                       //  Propability for electron attachment (for 1m)
  Int_t                fExBOn;                              //  Switch for the ExB effects
  Float_t              fOmegaTau;                           //  Tangens of the Lorentz angle 
  Float_t              fLorentzFactor;                      //  Factor due to Lorentz force
  Int_t                fPRFOn;                              //  Switch for the pad response
  Float_t             *fPRFsmp;                             //! Sampled pad response
  Int_t                fPRFbin;                             //  Number of bins for the PRF
  Float_t              fPRFlo;                              //  Lower boundary of the PRF
  Float_t              fPRFhi;                              //  Higher boundary of the PRF
  Float_t              fPRFwid;                             //  Bin width of the sampled PRF
  Int_t                fPRFpad;                             //  Distance to next pad in PRF
  Int_t                fTRFOn;                              //  Switch for the time response
  Float_t             *fTRFsmp;                             //! Integrated time response
  Int_t                fTRFbin;                             //  Number of bins for the TRF
  Float_t              fTRFlo;                              //  Lower boundary of the TRF
  Float_t              fTRFhi;                              //  Higher boundary of the TRF
  Float_t              fTRFwid;                             //  Bin width of the integrated TRF
  Int_t                fCTOn;                               //  Switch for cross talk
  Float_t             *fCTsmp;                              //! Integrated cross talk
  Int_t                fTCOn;                               //  Switch for the tail cancelation
  Int_t                fTCnexp;                             //  Number of exponential of the digital filter
  Float_t             *fTimeStruct1;                        //! Time Structure of Drift Cells
  Float_t             *fTimeStruct2;                        //! Time Structure of Drift Cells
  Int_t                fTimeStructOn;                       //  Switch for cell time structure
  Float_t              fAnodeWireOffset;                    //  Distance of first anode wire from pad edge

  Float_t              fDriftVelocity;                      //  Drift velocity (cm / mus)
  Float_t              fSamplingFrequency;                  //  Sampling Frequency in MHz
  Float_t              fPadCoupling;                        //  Pad coupling factor
  Float_t              fTimeCoupling;                       //  Time coupling factor (image charge of moving ions)

  // Clusterization parameter
  Int_t                fClusMaxThresh;                      //  Threshold value for cluster maximum
  Int_t                fClusSigThresh;                      //  Threshold value for cluster signal
  Int_t                fLUTOn;                              //  Switch for the lookup table method
  Int_t                fLUTbin;                             //  Number of bins of the LUT
  Float_t             *fLUT;                                //! The lookup table

  Float_t              fVDlo;                               //  Lower drift velocity, for interpolation
  Float_t              fVDhi;                               //  Higher drift velocity, for interpolation

 private:

  virtual void         SamplePRF();
  virtual void         SampleTRF();
  virtual void         FillLUT();
  virtual void         SampleTimeStruct();

  ClassDef(AliTRDparameter,6)                               //  TRD parameter class

};

#endif
