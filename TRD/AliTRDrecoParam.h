#ifndef ALITRDRECOPARAM_H
#define ALITRDRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Parameter class for the TRD reconstruction                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIDETECTORRECOPARAM_H
#include "AliDetectorRecoParam.h"
#endif

class AliTRDrecoParam : public AliDetectorRecoParam
{
public:
  enum AliTRDpidMethod{
    kLQPID = 0,
    kNNPID = 1
  };
  
  AliTRDrecoParam();
  AliTRDrecoParam(const AliTRDrecoParam &rec);
  ~AliTRDrecoParam() { }

  Int_t    GetADCbaseline() const           { return fADCbaseline;     };
  Double_t GetChi2Y() const                 { return fkChi2Y;    }
  Double_t GetChi2Z() const                 { return fkChi2Z;    }
  Double_t GetFindableClusters() const      { return fkFindable; }
  Double_t GetMaxTheta() const              { return fkMaxTheta; }
  Double_t GetMaxPhi() const                { return fkMaxPhi;   }
  Int_t    GetNdEdxSlices() const           { return fkPIDMethod == kNNPID ? kNNslices : kLQslices;}
  AliTRDpidMethod    GetPIDMethod() const   { return fkPIDMethod;}
  Double_t GetPlaneQualityThreshold() const { return fkPlaneQualityThreshold; }
  Double_t GetRoad0y() const                { return fkRoad0y;   }
  Double_t GetRoad0z() const                { return fkRoad0z;   }
  Double_t GetRoad1y() const                { return fkRoad1y;   }
  Double_t GetRoad1z() const                { return fkRoad1z;   }
  Double_t GetRoad2y() const                { return fkRoad2y;   }
  Double_t GetRoad2z() const                { return fkRoad2z;   }
  Double_t GetTrackLikelihood() const       { return fkTrackLikelihood;       }
  Int_t    GetStreamLevel() const           { return fkStreamLevel;           }
  inline void GetSysCovMatrix(Double_t *sys);  
  Double_t GetMinMaxCutSigma() const        { return fMinMaxCutSigma;     };
  Double_t GetMinLeftRightCutSigma() const  { return fMinLeftRightCutSigma;  };
  Double_t GetClusMaxThresh() const         { return fClusMaxThresh;   };
  Double_t GetClusSigThresh() const         { return fClusSigThresh;   };
  Int_t    GetTCnexp() const                { return fTCnexp;          };

        
  static   AliTRDrecoParam *GetLowFluxParam();
  static   AliTRDrecoParam *GetHighFluxParam();
  static   AliTRDrecoParam *GetCosmicTestParam();

  Bool_t   IsClusterSharing() const         { return TestBit(kClusterSharing);}
  Bool_t   IsLUT() const                    { return TestBit(kLUT);}
  Bool_t   IsSeeding() const                { return TestBit(kSeeding); }
  Bool_t   IsTailCancelation() const        { return TestBit(kTC);}
  Bool_t   IsVertexConstrained() const      { return TestBit(kVertexConstrained); }


  void     SetClusterSharing(Bool_t share = kTRUE)            { SetBit(kClusterSharing, share);  };
  void     SetPIDMethod(AliTRDpidMethod pid)                  { fkPIDMethod = pid; };
  void     SetSeeding(Bool_t so = kTRUE)             { SetBit(kSeeding, so); }
  void     SetVertexConstrained(Bool_t vc = kTRUE) { SetBit(kVertexConstrained, vc); }
  void     SetStreamLevel(Int_t streamLevel= 1)               { fkStreamLevel = streamLevel; }
  void     SetLUT(Bool_t lut = kTRUE)                            { SetBit(kLUT, lut);};
  void     SetMinMaxCutSigma(Float_t minMaxCutSigma)          { fMinMaxCutSigma   = minMaxCutSigma; };
  void     SetMinLeftRightCutSigma(Float_t minLeftRightCutSigma) { fMinLeftRightCutSigma   = minLeftRightCutSigma; };
  void     SetClusMaxThresh(Float_t thresh)                   { fClusMaxThresh   = thresh; };
  void     SetClusSigThresh(Float_t thresh)                   { fClusSigThresh   = thresh; };
  void     SetTailCancelation(Bool_t tc = kTRUE)                 { SetBit(kTC, tc);  };
  void     SetNexponential(Int_t nexp)                        { fTCnexp          = nexp;   };
  void     SetADCbaseline(Int_t base)                         { fADCbaseline     = base;   };
  inline void SetSysCovMatrix(Double_t *sys);

private:
  enum{
    kNNslices = 8
   ,kLQslices = 3
  };
  
  enum{
    kClusterSharing    = 1 // Toggle cluster sharing
   ,kSeeding           = 2 // Do stand alone tracking in the TRD
   ,kVertexConstrained = 3 // Perform vertex constrained fit
   ,kLUT               = 4 // 
   ,kTC                = 5 // tail cancelation
  };

  AliTRDpidMethod fkPIDMethod;       // PID method selector 0(LQ) 1(NN)
  Double_t  fkMaxTheta;              // Maximum theta
  Double_t  fkMaxPhi;                // Maximum phi

  Double_t  fkRoad0y;                // Road for middle cluster
  Double_t  fkRoad0z;                // Road for middle cluster

  Double_t  fkRoad1y;                // Road in y for seeded cluster
  Double_t  fkRoad1z;                // Road in z for seeded cluster

  Double_t  fkRoad2y;                // Road in y for extrapolated cluster
  Double_t  fkRoad2z;                // Road in z for extrapolated cluster
  
  Double_t  fkPlaneQualityThreshold; // Quality threshold
  Double_t  fkFindable;              // Ratio of clusters from a track in one chamber which are at minimum supposed to be found.
  Double_t  fkChi2Z;                 // Max chi2 on the z direction for seeding clusters fit
  Double_t  fkChi2Y;                 // Max chi2 on the y direction for seeding clusters Rieman fit
  Double_t  fkTrackLikelihood;       // Track likelihood for tracklets Rieman fit
  Int_t     fkStreamLevel;					 // Streaming Level in TRD Reconstruction
  
  Double_t  fSysCovMatrix[5];        // Systematic uncertainty from calibration and alignment for each tracklet

  // Clusterization parameter
  Double_t  fMinMaxCutSigma;         // Threshold sigma noise pad middle
  Double_t  fMinLeftRightCutSigma;   // Threshold sigma noise sum pad
  Double_t  fClusMaxThresh;          // Threshold value for cluster maximum
  Double_t  fClusSigThresh;          // Threshold value for cluster signal
  Int_t     fTCnexp;                 // Number of exponentials, digital filter
  
  // ADC parameter
  Int_t     fADCbaseline;            // ADC baseline to be subtracted

  ClassDef(AliTRDrecoParam, 4)       // Reconstruction parameters for TRD detector

};

//___________________________________________________
inline void AliTRDrecoParam::GetSysCovMatrix(Double_t *sys)
{
  if(!sys) return;
  memcpy(sys, fSysCovMatrix, 5*sizeof(Double_t));
}

//___________________________________________________
inline void AliTRDrecoParam::SetSysCovMatrix(Double_t *sys)
{
  if(!sys) return;
  memcpy(fSysCovMatrix, sys, 5*sizeof(Double_t));
}



#endif
