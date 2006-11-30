// @(#) $Id$
// Original: AliHLTModelTrack.h,v 1.13 2004/06/15 10:26:57 hristov 

#ifndef AliHLTTPCModel_Track
#define AliHLTTPCModel_Track

#include "AliHLTTPCTrack.h"
#include "AliHLTTPCModels.h"

class AliHLTTPCModelTrack : public AliHLTTPCTrack {

 public:
  AliHLTTPCModelTrack();
  virtual ~AliHLTTPCModelTrack();
  
  void Init(Int_t slice,Int_t patch);
  void CalculateClusterWidths(Int_t row,Bool_t parametrize=kFALSE);
  void SetCluster(Int_t row,Float_t dpad,Float_t dtime,Float_t charge,Float_t sigmaY2,Float_t sigmaZ2,Int_t npads);
  void FillModel();
  void FillTrack();
  void Print(Bool_t everything=kTRUE);
  void Set(AliHLTTPCTrack *tpt);

  void SetPadHit(Int_t row,Float_t f);
  void SetTimeHit(Int_t row,Float_t f);
  void SetCrossingAngleLUT(Int_t row,Float_t angle);
  void SetOverlap(Int_t row,Int_t id);
  void SetClusterLabel(Int_t row,Int_t *trackID);
  void SetNClusters(Int_t i) {fNClusters = i;}
  
  Int_t GetNPresentClusters();
  Bool_t IsPresent(Int_t row);
  Bool_t IsSet(Int_t row);
  
  AliHLTTPCClusterModel *GetClusters() {return fClusters;}
  AliHLTTPCTrackModel *GetModel() {return fTrackModel;}
  AliHLTTPCClusterModel *GetClusterModel(Int_t row);
  Int_t *GetOverlaps(Int_t row);
  Int_t GetNOverlaps(Int_t row);
  Int_t GetNPads(Int_t row);
  Int_t GetSlice(Int_t row);
  Float_t GetPadHit(Int_t row);
  Float_t GetTimeHit(Int_t row);
  Float_t GetCrossingAngleLUT(Int_t row);
  Float_t GetParSigmaY2(Int_t row);
  Float_t GetParSigmaZ2(Int_t row);
  Bool_t GetPad(Int_t row,Float_t &pad);
  Bool_t GetTime(Int_t row,Float_t &time);
  Bool_t GetClusterCharge(Int_t row,Int_t &charge);
  Bool_t GetSigmaY2(Int_t row,Float_t &sigma2);
  Bool_t GetSigmaZ2(Int_t row,Float_t &sigma2);
  Bool_t GetPadResidual(Int_t row,Float_t &res);
  Bool_t GetTimeResidual(Int_t row,Float_t &res);
  Bool_t GetSigmaYResidual(Int_t row,Float_t &res);
  Bool_t GetSigmaZResidual(Int_t row,Float_t &res);
  Int_t GetNClusters() const {return fNClusters;}
  void GetClusterLabel(Int_t row,Int_t *trackID);
    
 private:
  
  Short_t fClusterCharge; //Average cluster charge
  AliHLTTPCClusterModel *fClusters; //! Clusters
  AliHLTTPCTrackModel *fTrackModel; //! Track model
  Short_t fNClusters; // Number of clusters
  Int_t fMaxOverlaps; // Max overlaps (?)
  Int_t *fNoverlaps; //! Number of overlaps
  Int_t **fOverlap; //! Table of overlaps(?)
  Float_t *fParSigmaY2;    //! Parameter SigmaY2 (?)
  Float_t *fParSigmaZ2;    //! Parameter SigmaZ2 (?)
  Float_t *fCrossingAngle; //! Crossing angle
  Int_t fPatch; // Current patch
  Bool_t fArraysCreated; // Flag if arrays were created

  //Crossing points with padrows
  Float_t *fPad;  //! Current pad
  Float_t *fTime;  //! Current time
  
  void DeleteArrays();
  
  ClassDef(AliHLTTPCModelTrack,1)

};

#endif
