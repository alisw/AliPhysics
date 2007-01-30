#ifndef ALITRDSIMPARAM_H
#define ALITRDSIMPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Class containing constant simulation parameters                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliTRDSimParam : public TObject {
  
 public:
  
          enum { kNplan =   6
               , kNcham =   5
               , kNsect =  18
               , kNdet  = 540 };

  static  AliTRDSimParam *Instance();
  static  void     Terminate();
  
  AliTRDSimParam(const AliTRDSimParam &p);   
  AliTRDSimParam &operator=(const AliTRDSimParam &p); 

  virtual void     Copy(TObject &p) const;
  
          void     SetGasGain(Float_t gasgain)               { fGasGain           = gasgain;         }
          void     SetNoise(Float_t noise)                   { fNoise             = noise;           }
          void     SetChipGain(Float_t chipgain)             { fChipGain          = chipgain;        }
          void     SetADCoutRange(Float_t range)             { fADCoutRange       = range;           }
          void     SetADCinRange(Float_t range)              { fADCinRange        = range;           }
          void     SetADCthreshold(Int_t thresh)             { fADCthreshold      = thresh;          }
          void     SetADCbaseline(Int_t basel)               { fADCbaseline       = basel;           }   
          void     SetDiffusion(Int_t diffOn = 1)            { fDiffusionOn       = diffOn;          }
          void     SetElAttach(Int_t elOn = 1)               { fElAttachOn        = elOn;            }
          void     SetElAttachProp(Float_t prop)             { fElAttachProp      = prop;            }
          void     SetTimeResponse(Int_t trfOn = 1)          { fTRFOn             = trfOn; ReInit(); }  
          void     SetCrossTalk(Int_t ctOn = 1)              { fCTOn              = ctOn; ReInit();  }
          void     SetPadCoupling(Float_t v)                 { fPadCoupling       = v;               }
          void     SetTimeCoupling(Float_t v)                { fTimeCoupling      = v;               }
          void     SetAnodeWireOffset(Float_t offset = 0.25) { fAnodeWireOffset   = offset;          }
          void     SetTimeStruct(Bool_t tsOn = 1)            { fTimeStructOn      = tsOn;            }
          void     SetPadResponse(Int_t prfOn = 1)           { fPRFOn             = prfOn;           }
    
          Float_t  GetGasGain() const                        { return fGasGain;                      }
          Float_t  GetNoise() const                          { return fNoise;                        }
          Float_t  GetChipGain() const                       { return fChipGain;                     }
          Float_t  GetADCoutRange() const                    { return fADCoutRange;                  }
          Float_t  GetADCinRange() const                     { return fADCinRange;                   }
          Int_t    GetADCthreshold() const                   { return fADCthreshold;                 }
          Int_t    GetADCbaseline() const                    { return fADCbaseline;                  }
          Float_t  GetTRFlo() const                          { return fTRFlo;                        }
          Float_t  GetTRFhi() const                          { return fTRFhi;                        }
          Float_t  GetPadCoupling() const                    { return fPadCoupling;                  }
          Float_t  GetTimeCoupling() const                   { return fTimeCoupling;                 }
          Float_t  GetAnodeWireOffset() const                { return fAnodeWireOffset;              }

          Bool_t   DiffusionOn() const                       { return fDiffusionOn;                  }
          Bool_t   ElAttachOn() const                        { return fElAttachOn;                   } 
          Float_t  GetElAttachProp() const                   { return fElAttachProp;                 }
          Bool_t   TRFOn() const                             { return fTRFOn;                        }
          Bool_t   CTOn() const                              { return fCTOn;                         }
          Bool_t   TimeStructOn() const                      { return fTimeStructOn;                 }
          Bool_t   PRFOn() const                             { return fPRFOn;                        }

          Double_t TimeResponse(Double_t time) const;  
          Double_t CrossTalk(Double_t time) const; 
  
protected:

  static AliTRDSimParam* fgInstance;   //  Instance of this class (singleton implementation)
  static  Bool_t   fgTerminated;       //  Defines if this class has already been terminated and
                                       //  therefore does not return instances in GetInstance anymore
  
          // Digitization parameter
          Float_t  fGasGain;           //  Gas gain
          Float_t  fNoise;             //  Electronics noise
          Float_t  fChipGain;          //  Electronics gain
  
          Float_t  fADCoutRange;       //  ADC output range (number of channels)
          Float_t  fADCinRange;        //  ADC input range (input charge)
          Int_t    fADCthreshold;      //  ADC threshold in ADC channel
          Int_t    fADCbaseline;       //  ADC baseline in ADC chann
  
          Int_t    fDiffusionOn;       //  Switch for the diffusion
  
          Int_t    fElAttachOn;        //  Switch for the electron attachment
          Float_t  fElAttachProp;      //  Propability for electron attachment (for 1m)
  
          Int_t    fTRFOn;             //  Switch for the time response
          Float_t *fTRFsmp;            //! Integrated time response
          Int_t    fTRFbin;            //  Number of bins for the TRF
          Float_t  fTRFlo;             //  Lower boundary of the TRF
          Float_t  fTRFhi;             //  Higher boundary of the TRF
          Float_t  fTRFwid;            //  Bin width of the integrated TRF
  
          Int_t    fCTOn;              //  Switch for cross talk
          Float_t *fCTsmp;             //! Integrated cross talk
  
          Float_t  fAnodeWireOffset;   //  Distance of first anode wire from pad edge
          Float_t  fPadCoupling;       //  Pad coupling factor
          Float_t  fTimeCoupling;      //  Time coupling factor (image charge of moving ions)
          Int_t    fTimeStructOn;      //  Switch for cell time structure
  
          Int_t    fPRFOn;             //  Switch for the pad response
  
 private:

  // This is a singleton, constructor is private!  
  AliTRDSimParam();
  virtual ~AliTRDSimParam();

          void Init();
          void ReInit();
          void SampleTRF();
  
  ClassDef(AliTRDSimParam,1)          // The TRD simulation parameters

};

#endif
