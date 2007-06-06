#ifndef ALIITSRECOPARAM_H
#define ALIITSRECOPARAM_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with ITS reconstruction parameters                                  //
// Origin: andrea.dainese@lnl.infn.it                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TObject.h"

//--------------- move from AliITSrecoV2.h ---------------------------    
const Int_t kMaxLayer  = 6;

const Int_t kLayersNotToSkip[6]={0,0,0,0,0,0};
const Int_t kLastLayerToTrackTo=0;

const Int_t kMaxClusterPerLayer=7000*10;
const Int_t kMaxClusterPerLayer5=7000*10*2/5;
const Int_t kMaxClusterPerLayer10=7000*10*2/10;
const Int_t kMaxClusterPerLayer20=7000*10*2/20;
const Int_t kMaxDetectorPerLayer=1000;
//------------- end of move from AliITSrecoV2.h --------------------


class AliITSRecoParam : public TObject
{
 public: 
  AliITSRecoParam();
  virtual ~AliITSRecoParam();

  static AliITSRecoParam *GetLowFluxParam();// make reco parameters for low flux env.
  static AliITSRecoParam *GetHighFluxParam();// make reco parameters for high flux env. 
  static AliITSRecoParam *GetCosmicTestParam();// special setting for cosmic  

  
  Double_t GetSigmaY2(Int_t i) const { return fSigmaY2[i]; }
  Double_t GetSigmaZ2(Int_t i) const { return fSigmaZ2[i]; }

  Double_t GetMaxSnp() const { return fMaxSnp; }

  Double_t GetNSigmaYLayerForRoadY() const { return fNSigmaYLayerForRoadY; }
  Double_t GetNSigmaRoadY() const { return fNSigmaRoadY; }
  Double_t GetNSigmaZLayerForRoadZ() const { return fNSigmaZLayerForRoadZ; }
  Double_t GetNSigmaRoadZ() const { return fNSigmaRoadZ; }
  Double_t GetNSigma2RoadYC() const { return fNSigma2RoadYC; }
  Double_t GetNSigma2RoadZC() const { return fNSigma2RoadZC; }
  Double_t GetNSigma2RoadYNonC() const { return fNSigma2RoadYNonC; }
  Double_t GetNSigma2RoadZNonC() const { return fNSigma2RoadZNonC; }

  Double_t GetChi2PerCluster() const { return fChi2PerCluster; }
  Double_t GetMaxChi2PerCluster(Int_t i) const { return fMaxChi2PerCluster[i]; }
  Double_t GetMaxNormChi2NonC(Int_t i) const { return fMaxNormChi2NonC[i]; }
  Double_t GetMaxNormChi2C(Int_t i) const { return fMaxNormChi2C[i]; }
  Double_t GetMaxChi2() const { return fMaxChi2; }
  Double_t GetMaxChi2s(Int_t i) const { return fMaxChi2s[i]; }
  Double_t GetMaxChi2sR(Int_t i) const { return fMaxChi2sR[i]; }
  Double_t GetMaxChi2In() const { return fMaxChi2In; }
  Double_t GetVertexCut() const { return fVertexCut; }
  Double_t GetMaxRoad() const { return fMaxRoad; }

  Double_t GetXVdef() const { return fXV; }
  Double_t GetYVdef() const { return fYV; }
  Double_t GetZVdef() const { return fZV; }
  Double_t GetSigmaXVdef() const { return fSigmaXV; }
  Double_t GetSigmaYVdef() const { return fSigmaYV; }
  Double_t GetSigmaZVdef() const { return fSigmaZV; }
  
  void SetLayersParameters();
  //
 protected:
  //
  // spatial resolutions of the detectors
  Double_t fSigmaY2[kMaxLayer];
  Double_t fSigmaZ2[kMaxLayer];
  //
  Double_t fMaxSnp; // maximum of sin(phi)  (MI)
  //
  // search road (MI)
  Double_t fNSigmaYLayerForRoadY;
  Double_t fNSigmaRoadY;
  Double_t fNSigmaZLayerForRoadZ;
  Double_t fNSigmaRoadZ;
  Double_t fNSigma2RoadZC;
  Double_t fNSigma2RoadYC;
  Double_t fNSigma2RoadZNonC;
  Double_t fNSigma2RoadYNonC;
  //
  // chi2 cuts
  Double_t fMaxChi2PerCluster[kMaxLayer-1]; // max chi2 for MIP (MI)
  Double_t fMaxNormChi2NonC[kMaxLayer]; //max norm chi2 for non constrained tracks (MI)
  Double_t fMaxNormChi2C[kMaxLayer];  //max norm chi2 for constrained tracks (MI)
  Double_t fMaxChi2; // used to initialize variables needed to find minimum chi2 (MI,V2)
  Double_t fMaxChi2s[kMaxLayer];   // max predicted chi2 (cluster & track prol.) (MI)
  //
  Double_t fMaxRoad;   // (V2)
  //
  Double_t fMaxChi2In; // (NOT USED)
  Double_t fMaxChi2sR[kMaxLayer];  // (NOT USED) 
  Double_t fChi2PerCluster; // (NOT USED)
  //
  // default primary vertex (MI,V2)
  Double_t fXV; 
  Double_t fYV;
  Double_t fZV;
  Double_t fSigmaXV;
  Double_t fSigmaYV;
  Double_t fSigmaZV;
  Double_t fVertexCut; // (V2)
  //
  ClassDef(AliITSRecoParam,1) // ITS reco parameters
};

#endif
