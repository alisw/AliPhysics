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

#ifndef ALITRDPIDRESPONSE_H
#include "AliTRDPIDResponse.h"
#endif

class TString;

class AliTRDrecoParam : public AliDetectorRecoParam
{
public:
  enum ETRDReconstructionTask{
    kClusterizer = 0,
    kTracker = 1,
    kPID = 2,
    kTRDreconstructionTasks = 3
  };
  enum ETRDflags {
    kDriftGas
    ,kVertexConstraint
    ,kTailCancelation
    ,kImproveTracklet
    ,kLUT
    ,kGAUS
    ,kClusterSharing
    ,kSteerPID
    ,kEightSlices
    ,kCheckTimeConsistency
    ,kLQ2D
  };
  AliTRDrecoParam();
  AliTRDrecoParam(const AliTRDrecoParam &rec);
  const AliTRDrecoParam& operator=(const AliTRDrecoParam &rec);
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
  inline Int_t    GetPIDLQslices() const;
  inline AliTRDPIDResponse::ETRDPIDMethod GetPIDmethod() const;
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
  Double_t GetRoadzMultiplicator() const    { return fkRoadzMultiplicator; }
  Double_t GetTrackLikelihood() const       { return fkTrackLikelihood;       }
  inline void GetSysCovMatrix(Double_t *sys) const;  
  inline void GetTCParams(Double_t *par) const;
  inline Int_t GetStreamLevel(ETRDReconstructionTask task) const;
  const TString *GetRawStreamVersion() const{ return &fRawStreamVersion; };
  Double_t GetMinMaxCutSigma() const        { return fMinMaxCutSigma;     };
  Double_t GetMinLeftRightCutSigma() const  { return fMinLeftRightCutSigma;  };
  Double_t GetClusMaxThresh() const         { return fClusMaxThresh;   };
  Double_t GetClusSigThresh() const         { return fClusSigThresh;   };
  Int_t    GetTCnexp() const                { return fTCnexp;          };
  Int_t    GetNumberOfPresamples()  const   { return fNumberOfPresamples;}
  Int_t    GetNumberOfPostsamples() const   { return fNumberOfPostsamples;}
  Int_t    GetNumberOfSeedConfigs() const   { return fNumberOfConfigs;}
  Int_t    GetRecEveryNTB() const           { return fRecEveryNTB; }
  Bool_t   IsArgon() const                  { return TESTBIT(fFlags, kDriftGas); }
  Bool_t   IsCheckTimeConsistency() const   { return kCheckTimeConsistency;}
  Bool_t   IsOverPtThreshold(Double_t pt) const {return Bool_t(pt>fkPtThreshold);}
  Bool_t   IsXenon() const                  { return !TESTBIT(fFlags, kDriftGas); }
  Bool_t   IsPIDNeuralNetwork() const       { return TESTBIT(fFlags, kSteerPID);}
  Bool_t   IsVertexConstrained() const      { return TESTBIT(fFlags, kVertexConstraint); }
  Bool_t   IsEightSlices() const            { return TESTBIT(fFlags, kEightSlices);}
  Bool_t   HasImproveTracklets() const      { return TESTBIT(fFlags, kImproveTracklet);}
  Bool_t   UseClusterSharing() const        { return TESTBIT(fFlags, kClusterSharing);}
  Bool_t   UseLUT() const                   { return TESTBIT(fFlags, kLUT);}
  Bool_t   UseGAUS() const                  { return TESTBIT(fFlags, kGAUS);}
  Bool_t   UseTailCancelation() const       { return TESTBIT(fFlags, kTailCancelation); }
        
  static   AliTRDrecoParam *GetLowFluxParam();
  static   AliTRDrecoParam *GetLowFluxHLTParam();
  static   AliTRDrecoParam *GetHighFluxParam();
  static   AliTRDrecoParam *GetHighFluxHLTParam();
  static   AliTRDrecoParam *GetCosmicTestParam();

  void     SetArgon(Bool_t b = kTRUE)                         {if(b) SETBIT(fFlags, kDriftGas); else CLRBIT(fFlags, kDriftGas);}
  void     SetCheckTimeConsistency(Bool_t b = kTRUE)          {if(b) SETBIT(fFlags, kCheckTimeConsistency); else CLRBIT(fFlags, kCheckTimeConsistency);}
  void     SetClusterSharing(Bool_t b = kTRUE)                {if(b) SETBIT(fFlags, kClusterSharing); else CLRBIT(fFlags, kClusterSharing);}
  void     SetEightSlices(Bool_t b = kTRUE)                   {if(b) SETBIT(fFlags, kEightSlices); else CLRBIT(fFlags, kEightSlices);}
  void     SetImproveTracklets(Bool_t b = kTRUE)              {if(b) SETBIT(fFlags, kImproveTracklet); else CLRBIT(fFlags, kImproveTracklet);}
  void     SetLUT(Bool_t b=kTRUE)                             {if(b) SETBIT(fFlags, kLUT); else CLRBIT(fFlags, kLUT);}
  void     SetGAUS(Bool_t b=kTRUE)                            {if(b) SETBIT(fFlags, kGAUS); else CLRBIT(fFlags, kGAUS);}
  void     SetPIDNeuralNetwork(Bool_t b=kTRUE)                {if(b) SETBIT(fFlags, kSteerPID); else CLRBIT(fFlags, kSteerPID);}
  inline void  SetPIDmethod(AliTRDPIDResponse::ETRDPIDMethod method);
  void     SetPIDLQslices(Int_t s);
  void     SetTailCancelation(Bool_t b=kTRUE)                 {if(b) SETBIT(fFlags, kTailCancelation); else CLRBIT(fFlags, kTailCancelation);}
  void     SetXenon(Bool_t b = kTRUE)                         {if(b) CLRBIT(fFlags, kDriftGas); else SETBIT(fFlags, kDriftGas);}
  void     SetVertexConstrained()                             {SETBIT(fFlags, kVertexConstraint);}
  void     SetMaxTheta(Double_t maxTheta)                     {fkMaxTheta = maxTheta;}
  void     SetMaxPhi(Double_t maxPhi)                         {fkMaxPhi = maxPhi;}
  void     SetFindableClusters(Double_t r)                    {fkFindable = r;}
  void     SetChi2Y(Double_t chi2)                            {fkChi2Y = chi2;}
  void     SetChi2Z(Double_t chi2)                            {fkChi2Z = chi2;}
  void     SetChi2YSlope(Double_t chi2YSlope)                 {fkChi2YSlope = chi2YSlope;}
  void     SetChi2ZSlope(Double_t chi2ZSlope)                 {fkChi2ZSlope = chi2ZSlope;}
	void	   SetChi2YCut(Double_t chi2Cut)                      {fkChi2YCut = chi2Cut; }
  void     SetPhiSlope(Double_t phiSlope)                     {fkPhiSlope = phiSlope;}
  void     SetNMeanClusters(Double_t meanNclusters)           {fkNMeanClusters = meanNclusters;}
  void     SetNSigmaClusters(Double_t sigmaNclusters)         {fkNSigmaClusters = sigmaNclusters;} 
  void     SetRawStreamVersion(const Char_t *version)         {fRawStreamVersion = version; }
  void     SetRoadzMultiplicator(Double_t mult)               {fkRoadzMultiplicator = mult; } 
  void     SetMinMaxCutSigma(Float_t minMaxCutSigma)          { fMinMaxCutSigma   = minMaxCutSigma; }
  void     SetMinLeftRightCutSigma(Float_t minLeftRightCutSigma) { fMinLeftRightCutSigma   = minLeftRightCutSigma; };
  void     SetClusMaxThresh(Float_t thresh)                   { fClusMaxThresh   = thresh; };
  void     SetClusSigThresh(Float_t thresh)                   { fClusSigThresh   = thresh; };
  inline void SetPIDThreshold(Double_t *pid);
  void     SetPtThreshold(Double_t pt) {fkPtThreshold = pt;}
  void     SetNexponential(Int_t nexp)                        { fTCnexp          = nexp;   };
  inline void SetTCParams(Double_t *par);
  inline void SetStreamLevel(ETRDReconstructionTask task, Int_t level);
  inline void SetSysCovMatrix(Double_t *sys);
  void     SetNumberOfPresamples(Int_t n)                     { fNumberOfPresamples = n;}
  void     SetNumberOfPostsamples(Int_t n)                    { fNumberOfPostsamples = n;}
  void     SetRecEveryTwoTB()                                 { fRecEveryNTB = 2; fkNMeanClusters = 10; }

private:
  // Physics reference values for TRD
  Double_t  fkdNchdy;                // dNch/dy
  Double_t  fkMaxTheta;              // Maximum theta
  Double_t  fkMaxPhi;                // Maximum phi - momentum cut
  // Tracker params 
  Double_t  fkRoad0y;                // Road for middle cluster
  Double_t  fkRoad0z;                // Road for middle cluster

  Double_t  fkRoad1y;                // Road in y for seeded cluster
  Double_t  fkRoad1z;                // Road in z for seeded cluster

  Double_t  fkRoad2y;                // Road in y for extrapolated cluster
  Double_t  fkRoad2z;                // Road in z for extrapolated cluster
  Double_t  fkPtThreshold;           // pt threshold for using TRD points for updating Kalaman track
  Double_t  fkPlaneQualityThreshold; // Quality threshold
  Double_t  fkRoadzMultiplicator;    // Multiplicator for the Roads in z 
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
  Double_t  fPIDThreshold[AliTRDCalPID::kNMom];   // PID Thresholds for Electron candidate decision
  Int_t     fNumberOfConfigs;        // Used number of seed configurations

  // Reconstruction Options for TRD reconstruction
  Int_t     fStreamLevel[kTRDreconstructionTasks]; // Stream Level
  Long64_t  fFlags;                  // option Flags

  // Raw Reader Params
  TString   fRawStreamVersion;       // Raw Reader version

  // Clusterization parameter
  Double_t  fMinMaxCutSigma;         // Threshold sigma noise pad middle
  Double_t  fMinLeftRightCutSigma;   // Threshold sigma noise sum pad
  Double_t  fClusMaxThresh;          // Threshold value for cluster maximum
  Double_t  fClusSigThresh;          // Threshold value for cluster signal
  Int_t     fTCnexp;                 // Number of exponentials, digital filter
  Double_t  fTCParams[8];            // Tail Cancellation parameters for drift gases 
  Int_t     fRecEveryNTB;            // Reconstruct each nth timebin

  // ADC parameter
  Int_t     fNumberOfPresamples;     // number of presamples 
  Int_t     fNumberOfPostsamples;     // number of postsamples 

  ClassDef(AliTRDrecoParam, 12)       // Reconstruction parameters for TRD detector

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

//___________________________________________________
inline void AliTRDrecoParam::SetStreamLevel(ETRDReconstructionTask task, Int_t level){
  if(task >= kTRDreconstructionTasks) return;
  fStreamLevel[static_cast<Int_t>(task)] = level;
}

//___________________________________________________
inline Int_t AliTRDrecoParam::GetStreamLevel(ETRDReconstructionTask task) const{
  if(task >= kTRDreconstructionTasks) return 0;
  return fStreamLevel[static_cast<Int_t>(task)];
}

//___________________________________________________
inline void AliTRDrecoParam::GetTCParams(Double_t *par) const
{
  if(!par) return;
  if(IsArgon()) memcpy(par, &fTCParams[4], 4*sizeof(Double_t));
  else memcpy(par, &fTCParams[0], 4*sizeof(Double_t));
}

//___________________________________________________
inline void AliTRDrecoParam::SetTCParams(Double_t *par)
{
  if(!par) return;
  memcpy(fTCParams, par, 8*sizeof(Double_t));
}

//___________________________________________________
inline Int_t AliTRDrecoParam::GetPIDLQslices() const
{
  if(IsPIDNeuralNetwork()) return -1;
  return TESTBIT(fFlags, kLQ2D) ? 2 : 1;
}

//___________________________________________________
inline AliTRDPIDResponse::ETRDPIDMethod AliTRDrecoParam::GetPIDmethod() const
{
  AliTRDPIDResponse::ETRDPIDMethod method = AliTRDPIDResponse::kLQ1D;
  if(IsPIDNeuralNetwork()) method = AliTRDPIDResponse::kNN;
  else if(TESTBIT(fFlags, kLQ2D)) method = AliTRDPIDResponse::kLQ2D;
  return method;
}

//___________________________________________________
inline void  AliTRDrecoParam::SetPIDmethod(AliTRDPIDResponse::ETRDPIDMethod method)
{
  switch(method){
  case AliTRDPIDResponse::kLQ2D:
    CLRBIT(fFlags, kSteerPID); 
    SETBIT(fFlags, kLQ2D);
    break;
  case AliTRDPIDResponse::kNN:
    SETBIT(fFlags, kSteerPID); 
    break;
  case AliTRDPIDResponse::kLQ1D:
  default:
    CLRBIT(fFlags, kSteerPID); 
    CLRBIT(fFlags, kLQ2D);
    break;
  }
}

#endif
