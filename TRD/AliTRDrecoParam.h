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
#ifndef ALITRDCALPID_H
#include "Cal/AliTRDCalPID.h"
#endif

class AliTRDrecoParam : public AliDetectorRecoParam
{
public:
  AliTRDrecoParam();
  AliTRDrecoParam(const AliTRDrecoParam &rec);
  ~AliTRDrecoParam() { }

  Double_t GetChi2Y() const                 { return fkChi2Y;    }
  Double_t GetChi2Z() const                 { return fkChi2Z;    }
  Double_t GetChi2YCut() const              { return fkChi2YCut; }
  Double_t GetChi2ZCut() const              { return fkChi2ZCut; }
  Double_t GetPhiCut() const                { return fkPhiCut;   }
  Float_t  GetNClusters() const;
  Double_t GetNMeanClusters() const         { return fkNMeanClusters; }
  Double_t GetNSigmaClusters() const        { return fkNSigmaClusters; }
  Double_t GetFindableClusters() const      { return fkFindable; }
  Double_t GetMaxTheta() const              { return fkMaxTheta; }
  Double_t GetMaxPhi() const                { return fkMaxPhi;   }
  Double_t GetPlaneQualityThreshold() const { return fkPlaneQualityThreshold; }
  Double_t GetPIDThreshold(Float_t /*p*/) const { return 0.;}
  Double_t GetRoad0y() const                { return fkRoad0y;   }
  Double_t GetRoad0z() const                { return fkRoad0z;   }
  Double_t GetRoad1y() const                { return fkRoad1y;   }
  Double_t GetRoad1z() const                { return fkRoad1z;   }
  Double_t GetRoad2y() const                { return fkRoad2y;   }
  Double_t GetRoad2z() const                { return fkRoad2z;   }
  Double_t GetTrackLikelihood() const       { return fkTrackLikelihood;       }
  inline void GetSysCovMatrix(Double_t *sys) const;  
  Double_t GetMinMaxCutSigma() const        { return fMinMaxCutSigma;     };
  Double_t GetMinLeftRightCutSigma() const  { return fMinLeftRightCutSigma;  };
  Double_t GetClusMaxThresh() const         { return fClusMaxThresh;   };
  Double_t GetClusSigThresh() const         { return fClusSigThresh;   };
  Int_t    GetTCnexp() const                { return fTCnexp;          };
  Int_t     GetNumberOfPresamples()  const {return fNumberOfPresamples;}
  Int_t    GetNumberOfPostsamples() const {return fNumberOfPostsamples;}

        
  static   AliTRDrecoParam *GetLowFluxParam();
  static   AliTRDrecoParam *GetHighFluxParam();
  static   AliTRDrecoParam *GetCosmicTestParam();

  Bool_t   IsClusterSharing() const         { return TestBit(kClusterSharing);}
  Bool_t   IsLUT() const                    { return TestBit(kLUT);}
  Bool_t   IsTailCancelation() const        { return TestBit(kTC);}
  Bool_t   IsVertexConstrained() const      { return TestBit(kVertexConstrained); }
  Bool_t   HasImproveTracklets() const       { return TestBit(kImproveTracklet); }

  void     SetMaxTheta(Double_t maxTheta) {fkMaxTheta = maxTheta;}
  void     SetMaxPhi(Double_t maxPhi) {fkMaxPhi = maxPhi;}
  void     SetFindableClusters(Double_t r) {fkFindable = r;}
  void     SetChi2Y(Double_t chi2) {fkChi2Y = chi2;}
  void     SetChi2Z(Double_t chi2) {fkChi2Z = chi2;}
  void     SetChi2YCut(Double_t chi2YCut) {fkChi2YCut = chi2YCut;}
  void     SetChi2ZCut(Double_t chi2ZCut) {fkChi2ZCut = chi2ZCut;}
  void     SetPhiCut(Double_t phiCut) {fkPhiCut = phiCut;}
  void     SetNMeanClusters(Double_t meanNclusters) {fkNMeanClusters = meanNclusters;}
  void     SetNSigmaClusters(Double_t sigmaNclusters) {fkNSigmaClusters = sigmaNclusters;} 
  void     SetClusterSharing(Bool_t share = kTRUE)            { SetBit(kClusterSharing, share);}
  void     SetImproveTracklets(Bool_t improve = kTRUE)         { SetBit(kImproveTracklet, improve);}
  void     SetVertexConstrained(Bool_t vc = kTRUE)            { SetBit(kVertexConstrained, vc); }
  void     SetLUT(Bool_t lut = kTRUE)                         { SetBit(kLUT, lut);};
  void     SetMinMaxCutSigma(Float_t minMaxCutSigma)          { fMinMaxCutSigma   = minMaxCutSigma; }
  void     SetMinLeftRightCutSigma(Float_t minLeftRightCutSigma) { fMinLeftRightCutSigma   = minLeftRightCutSigma; };
  void     SetClusMaxThresh(Float_t thresh)                   { fClusMaxThresh   = thresh; };
  void     SetClusSigThresh(Float_t thresh)                   { fClusSigThresh   = thresh; };
  void     SetTailCancelation(Bool_t tc = kTRUE)              { SetBit(kTC, tc);  };
  inline void SetPIDThreshold(Double_t *pid);
  void     SetNexponential(Int_t nexp)                        { fTCnexp          = nexp;   };
  inline void SetSysCovMatrix(Double_t *sys);
  void     SetNumberOfPresamples(Int_t n)                     { fNumberOfPresamples = n;}
  void     SetNumberOfPostsamples(Int_t n)                    { fNumberOfPostsamples = n;}

private:
  enum{
    kTC                = BIT(14) // tail cancelation
   ,kLUT               = BIT(15) // look up table for cluster position determination 
   ,kClusterSharing    = BIT(16) // Toggle cluster sharing
   ,kVertexConstrained = BIT(17) // Perform vertex constrained fit
   ,kImproveTracklet   = BIT(18) // Improve tracklet in the SA TRD track finder 
  };
  // Physics reference values for TRD
  Double_t  fkdNchdy;                // dNch/dy
  Double_t  fkMaxTheta;              // Maximum theta
  Double_t  fkMaxPhi;                // Maximum phi - momentum cut

  Double_t  fkRoad0y;                // Road for middle cluster
  Double_t  fkRoad0z;                // Road for middle cluster

  Double_t  fkRoad1y;                // Road in y for seeded cluster
  Double_t  fkRoad1z;                // Road in z for seeded cluster

  Double_t  fkRoad2y;                // Road in y for extrapolated cluster
  Double_t  fkRoad2z;                // Road in z for extrapolated cluster
  
  Double_t  fkPlaneQualityThreshold; // Quality threshold
  Double_t  fkFindable;              // minimum ratio of clusters per tracklet supposed to be attached.
  Double_t  fkChi2Z;                 // Max chi2 on the z direction for seeding clusters fit
  Double_t  fkChi2Y;                 // Max chi2 on the y direction for seeding clusters Rieman fit
  Double_t  fkChi2YCut;              // Cut on the Chi2 in y-direction in the likelihood filter
  Double_t  fkChi2ZCut;              // Cut on the Chi2 in z-direction in the likelihood filter
  Double_t  fkPhiCut;                // Cut on the deviation of the phi angles between tracklet and track fit (lik. filter)
  Double_t  fkNMeanClusters;         // Mean number of clusters per tracklet
  Double_t  fkNSigmaClusters;        // Sigma of the number of clusters per tracklet
  Double_t  fkNClusterNoise;         // ratio of noisy clusters to the true one
  Double_t  fkNMeanTracklets;        // Mean number of tracklets per track
  Double_t  fkTrackLikelihood;       // Track likelihood for tracklets Rieman fit
  
  Double_t  fSysCovMatrix[5];        // Systematic uncertainty from calibration and alignment for each tracklet
  Double_t  fPIDThreshold[AliTRDCalPID::kNMom];

  // Clusterization parameter
  Double_t  fMinMaxCutSigma;         // Threshold sigma noise pad middle
  Double_t  fMinLeftRightCutSigma;   // Threshold sigma noise sum pad
  Double_t  fClusMaxThresh;          // Threshold value for cluster maximum
  Double_t  fClusSigThresh;          // Threshold value for cluster signal
  Int_t     fTCnexp;                 // Number of exponentials, digital filter
  
  // ADC parameter
  Int_t     fNumberOfPresamples;     // number of presamples 
  Int_t     fNumberOfPostsamples;     // number of postsamples 

  ClassDef(AliTRDrecoParam, 6)       // Reconstruction parameters for TRD detector

};

//___________________________________________________
inline void AliTRDrecoParam::GetSysCovMatrix(Double_t *sys) const
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

//___________________________________________________
inline void AliTRDrecoParam::SetPIDThreshold(Double_t *pid)
{
  if(!pid) return;
  memcpy(fPIDThreshold, pid, AliTRDCalPID::kNMom*sizeof(Double_t));
}


#endif
