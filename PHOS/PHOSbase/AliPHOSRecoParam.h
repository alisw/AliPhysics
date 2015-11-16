#ifndef ALIPHOSRECOPARAM_H
#define ALIPHOSRECOPARAM_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */
                                              
// Base class for the PHOS reconstruction parameters.
// Do not use in the reconstruction; use derivative classes instead.

#include "TArrayF.h"
#include "AliDetectorRecoParam.h"

class AliPHOSRecoParam : public AliDetectorRecoParam {

public:

  AliPHOSRecoParam();
  AliPHOSRecoParam(const AliPHOSRecoParam& recoParam);
  AliPHOSRecoParam& operator = (const AliPHOSRecoParam& recoParam);
  virtual ~AliPHOSRecoParam() {}

  Float_t GetEMCClusteringThreshold() const { return fEMCClusteringThreshold;  }
  Float_t GetEMCLocalMaxCut()         const { return fEMCLocMaxCut;            }
  Float_t GetEMCRawDigitThreshold()   const { return fEMCRawDigitThreshold;  }
  Float_t GetEMCMinE()                const { return fEMCMinE;                 }
  Float_t GetEMCLogWeight()           const { return fEMCW0;                   }
  Float_t GetEMCSampleQualityCut()    const { return fEMCSampleQualityCut;     }
  Float_t GetEMCEcoreRadius()         const { return fEMCEcoreRadius;          }
  const Float_t * GetNonlinearityParams()   const { return fNonLinearityParams.GetArray() ;}
  Bool_t  EMCEcore2ESD()              const { return fEMCEcore2ESD;            }
  Bool_t  EMCSubtractPedestals()      const { return fEMCSubtractPedestals;    }
  Bool_t  EMCToUnfold()               const { return fEMCUnfold;               }
  const char* EMCFitterVersion()      const { return fEMCFitterVersion.Data(); }
  Bool_t  GetEMCEnergyCorrectionOn()  const { return fEMCEnergyCorrectionOn;   }
  Int_t   GetGlobalAltroOffset()      const { return fGlobalAltroOffset ;      }
  Int_t   GetGlobalAltroThreshold()   const { return fGlobalAltroThreshold ;   }
  Float_t GetTimeGateAmpThresh()      const { return fTimeGateAmpThresh ;      }
  Float_t GetTimeGateLow()            const { return fTimeGateLow ;            }
  Float_t GetTimeGateHigh()           const { return fTimeGateHigh ;           }
  Bool_t  GetSubtractL1phase()        const { return fSubtractL1phase;         }

  const char* GetNonlinearityCorrectionVersion()const{return fNonlinearityCorrVersion.Data();}

  Float_t GetCPVClusteringThreshold() const { return fCPVClusteringThreshold;  }
  Float_t GetCPVLocalMaxCut()         const { return fCPVLocMaxCut;            }
  Float_t GetCPVMinE()                const { return fCPVMinE;                 }
  Float_t GetCPVLogWeight()           const { return fCPVW0;                   }
  Bool_t  CPVToUnfold()               const { return fCPVUnfold;               }

  void SetEMCClusteringThreshold(Float_t cluth)      { fEMCClusteringThreshold=cluth;   }
  void SetEMCLocalMaxCut(Float_t cut)                { fEMCLocMaxCut          =cut;     }
  void SetEMCRawDigitThreshold(Float_t rawDigTh)     { fEMCRawDigitThreshold  =rawDigTh;}
  void SetEMCMinE(Float_t minE)                      { fEMCMinE               =minE;    }
  void SetEMCLogWeight(Float_t w)                    { fEMCW0                 =w;       }
  void SetEMCSampleQualityCut(Float_t qu)            { fEMCSampleQualityCut   =qu;      }
  void SetEMCEcoreRadius(Float_t rCore)              { fEMCEcoreRadius        =rCore;   }
  void SetEMCEcore2ESD(Bool_t ecore)                 { fEMCEcore2ESD          =ecore;   }
  void SetEMCSubtractPedestals(Bool_t subtract)      { fEMCSubtractPedestals  =subtract;} 
  void SetEMCFitterVersion(const char* version="v1") { fEMCFitterVersion     =version ; }
  void SetEMCUnfolding(Bool_t toUnfold=kFALSE)       { fEMCUnfold             =toUnfold;}
  void SetEMCEnergyCorrectionOn(Bool_t on=kTRUE)     { fEMCEnergyCorrectionOn =on;      }
  //Make sure to set first version and then parameters, otherwise parameters will be overwritten by default ones.
  void SetNonlinearityParams(Int_t n, Float_t * params){fNonLinearityParams.Set(n,params);}
  void SetGlobalAltroOffset(Int_t offset=5)          { fGlobalAltroOffset     =offset ; }
  void SetGlobalAltroThreshold(Int_t ZSth=5)         { fGlobalAltroThreshold  =ZSth;    }
  void SetTimeGateAmpThresh(Float_t thrs=10)         { fTimeGateAmpThresh     = thrs ;  }
  void SetTimeGateLow(Float_t gate=1.e-7)            { fTimeGateLow = gate ;            }
  void SetTimeGateHigh(Float_t gate=1.e-8)           { fTimeGateHigh = gate;            }
  void SetSubtractL1phase(Bool_t a=kTRUE)            { fSubtractL1phase = a ;          }
  
  void SetNonlinearityCorrectionVersion(const char * ver="Gustavo2005");

  void SetCPVClusteringThreshold(Float_t cluth)      { fCPVClusteringThreshold=cluth;   }
  void SetCPVLocalMaxCut(Float_t cut)                { fCPVLocMaxCut          =cut;     }
  void SetCPVMinE(Float_t minE)                      { fCPVMinE               =minE;    }
  void SetCPVLogWeight(Float_t w)                    { fCPVW0                 =w;       }
  void SetCPVUnfolding(Bool_t toUnfold=kFALSE)       { fCPVUnfold            =toUnfold;}

  virtual void Print(const Option_t *option="RecoParam") const;

  static AliPHOSRecoParam* GetDefaultParameters();
  static const  TObjArray* GetMappings();
  static const  TObject*   GetTriggerParameters();

protected:

  TArrayF fNonLinearityParams;     // EMC: Array of non-linearity correction parameters
  Float_t fEMCClusteringThreshold; // EMC: Min.digit energy to start a new cluster, in GeV
  Float_t fEMCLocMaxCut;           // EMC: Min.energy difference between two local maxima, in GeV
  Float_t fEMCRawDigitThreshold;   // EMC: Min.amplitude of a digit produced from raw data in ADC
  Float_t fEMCMinE;                // EMC: Min.E in the digits list associated with rec.point, in GeV
  Float_t fEMCW0;                  // EMC: Log.weight to evaluate a local coordinate of rec.point
  Float_t fEMCSampleQualityCut;    // EMC: Cut on pulse shape fit quality
  Float_t fTimeGateAmpThresh ;     // EMC: Threshold for good/bad time calculation
  Float_t fTimeGateLow ;           // EMC: Time difference between cells in PHOS cluster (bad time estimate)
  Float_t fTimeGateHigh ;          // EMC: Time difference between cells in PHOS cluster (good time estimate)
  Float_t fEMCEcoreRadius;         // EMC: Radius within which the core energy is calculated, in cm
  Bool_t  fEMCEcore2ESD;           // EMC: true if Ecore is stored in ESD instead of Etot
  Bool_t  fEMCSubtractPedestals;   // EMC: true if pedestal should be subtracted (in non-ZS)
  Bool_t  fEMCUnfold;              // EMC: true if overlapped clusters should be unfolded
  Bool_t  fEMCEnergyCorrectionOn;  // EMC: if true do non-linear correction of cluster energy
  Bool_t  fSubtractL1phase;        // EMC: should one sabtract L1 phase in time calculation
  TString fEMCFitterVersion ;      // EMC: AliPHOSRawFitter version
  TString fNonlinearityCorrVersion ;// EMC: choose which version of nenlinearity correction
  Int_t   fGlobalAltroOffset ;     // Offset used in ALTRO chips in SZ runs
  Int_t   fGlobalAltroThreshold ;  // Threshold used in ALTRO chips in SZ runs

  Float_t fCPVClusteringThreshold; // CPV: Min.digit energy to start a new cluster, in GeV
  Float_t fCPVLocMaxCut;           // CPV: Min.energy difference between two local maxima, in GeV
  Float_t fCPVMinE;                // CPV: Min.E in the digits list associated with rec.point, in GeV
  Float_t fCPVW0;                  // CPV: Log.weight to evaluate a local coordinate of rec.point
  Bool_t  fCPVUnfold;              // CPV: true if overlapped clusters should be unfolded

  static TObjArray* fgkMaps;       // ALTRO mappings for RCU0..RCU3
  static TObject* fgkTrigParams;   // PHOS trigger parameters

  ClassDef(AliPHOSRecoParam,13)
};

#endif
