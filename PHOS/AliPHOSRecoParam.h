#ifndef ALIPHOSRECOPARAM_H
#define ALIPHOSRECOPARAM_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */
                                              
// Base class for the PHOS reconstruction parameters.
// Do not use in the reconstruction; use derivative classes instead.

#include "TNamed.h"


class AliPHOSRecoParam : public TNamed {

public:

  AliPHOSRecoParam();
  AliPHOSRecoParam(const AliPHOSRecoParam& recoParam);
  AliPHOSRecoParam& operator = (const AliPHOSRecoParam& recoParam);
  virtual ~AliPHOSRecoParam() {}

  Float_t GetClusteringThreshold() const { return fClusteringThreshold;  }
  Float_t GetLocalMaxCut()         const { return fLocMaxCut;            }
  Float_t GetMinE()                const { return fMinE;                 }
  Float_t GetLogWeight()           const { return fW0;                   }
  Float_t GetSampleQualityCut()    const { return fSampleQualityCut;     }
  Float_t GetEcoreRadius()         const { return fEcoreRadius;          }
  Bool_t  Ecore2ESD()              const { return fEcore2ESD;            }
  Bool_t  SubtractPedestals()      const { return fSubtractPedestals;    }
  Bool_t  ToUnfold()               const { return fUnfold;               }
  const char* DecoderVersion()     const { return fDecoderVersion.Data();}

  void SetClusteringThreshold(Float_t cluth)      { fClusteringThreshold=cluth;   }
  void SetLocalMaxCut(Float_t cut)                { fLocMaxCut          =cut;     }
  void SetMinE(Float_t minE)                      { fMinE               =minE;    }
  void SetLogWeight(Float_t w)                    { fW0                 =w;       }
  void SetSampleQualityCut(Float_t qu)            { fSampleQualityCut   =qu;      }
  void SetEcoreRadius(Float_t rCore)              { fEcoreRadius        =rCore;   }
  void SetEcore2ESD(Bool_t ecore)                 { fEcore2ESD          =ecore;   }
  void SetSubtractPedestals(Bool_t subtract)      { fSubtractPedestals  =subtract;} 
  void SetDecoderVersion(const char* version="v1"){fDecoderVersion      =version ;}
  void SetUnfolding(Bool_t toUnfold=kFALSE)       {fUnfold              =toUnfold;}

protected:

  Float_t fClusteringThreshold; // Min.digit energy to start a new cluster, in GeV
  Float_t fLocMaxCut;           // Min.energy difference between two local maxima, in GeV
  Float_t fMinE;                // Min.E in the digits list associated with rec.point, in GeV
  Float_t fW0;                  // Log.weight to evaluate a local coordinate of rec.point
  Float_t fSampleQualityCut;    // Cut on pusle shape fit quality
  Float_t fEcoreRadius;         // Radius within which the core energy is calculated, in cm
  Bool_t  fEcore2ESD;           // true if Ecore is stored in ESD instead of Etot
  Bool_t  fSubtractPedestals;   // true if pedestal should be subtracted (in non-ZS)
  Bool_t  fUnfold;              // true if overlapped clusters should be unfolded
  TString fDecoderVersion ;     // AliPHOSRawDecoder version

  ClassDef(AliPHOSRecoParam,4)
};

#endif
