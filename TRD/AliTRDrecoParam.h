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
  Double_t GetChi2YSlope() const            { return fkChi2YSlope; }
  Double_t GetChi2ZSlope() const            { return fkChi2ZSlope; }
	Double_t GetChi2YCut() const              { return fkChi2YCut; }
  Double_t GetPhiSlope() const              { return fkPhiSlope;   }
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

  void     SetMaxTheta(Double_t maxTheta) {fkMaxTheta = maxTheta;}
  void     SetMaxPhi(Double_t maxPhi) {fkMaxPhi = maxPhi;}
  void     SetFindableClusters(Double_t r) {fkFindable = r;}
  void     SetChi2Y(Double_t chi2) {fkChi2Y = chi2;}
  void     SetChi2Z(Double_t chi2) {fkChi2Z = chi2;}
  void     SetChi2YSlope(Double_t chi2YSlope) {fkChi2YSlope = chi2YSlope;}
  void     SetChi2ZSlope(Double_t chi2ZSlope) {fkChi2ZSlope = chi2ZSlope;}
	void	   SetChi2YCut(Double_t chi2Cut) {fkChi2YCut = chi2Cut; }
  void     SetPhiSlope(Double_t phiSlope) {fkPhiSlope = phiSlope;}
  void     SetNMeanClusters(Double_t meanNclusters) {fkNMeanClusters = meanNclusters;}
  void     SetNSigmaClusters(Double_t sigmaNclusters) {fkNSigmaClusters = sigmaNclusters;} 
  void     SetMinMaxCutSigma(Float_t minMaxCutSigma)          { fMinMaxCutSigma   = minMaxCutSigma; }
  void     SetMinLeftRightCutSigma(Float_t minLeftRightCutSigma) { fMinLeftRightCutSigma   = minLeftRightCutSigma; };
  void     SetClusMaxThresh(Float_t thresh)                   { fClusMaxThresh   = thresh; };
  void     SetClusSigThresh(Float_t thresh)                   { fClusSigThresh   = thresh; };
  inline void SetPIDThreshold(Double_t *pid);
  void     SetNexponential(Int_t nexp)                        { fTCnexp          = nexp;   };
  inline void SetSysCovMatrix(Double_t *sys);
  void     SetNumberOfPresamples(Int_t n)                     { fNumberOfPresamples = n;}
  void     SetNumberOfPostsamples(Int_t n)                    { fNumberOfPostsamples = n;}

private:
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
  Double_t  fkChi2YSlope;            // Slope of the chi2-distribution in y-direction
  Double_t  fkChi2ZSlope;            // Slope of the chi2-distribution in z-direction
	Double_t  fkChi2YCut;							 // Cut on the Chi2 in y-direction in the likelihood filter
  Double_t  fkPhiSlope;              // Slope of the distribution of the deviation between track angle and tracklet angle
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

  ClassDef(AliTRDrecoParam, 7)       // Reconstruction parameters for TRD detector

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
