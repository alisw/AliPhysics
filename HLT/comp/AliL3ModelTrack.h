#ifndef AliL3Model_Track
#define AliL3Model_Track

#include "AliL3Track.h"
#include "AliL3Models.h"

class AliL3ModelTrack : public AliL3Track {

 private:
  
  Short_t fClusterCharge; //Average cluster charge
  AliL3ClusterModel *fClusters; //!
  AliL3TrackModel *fTrackModel; //!
  Short_t fNClusters;
  Int_t fOverlap;
  Float_t fXYResidualQ; //Quantization steps.
  Float_t fZResidualQ;
  
  //Crossing points with padrows
  Float_t *fPad; //!
  Float_t *fTime; //!
  
 public:
  AliL3ModelTrack();
  virtual ~AliL3ModelTrack();
  
  void Init(Int_t slice,Int_t patch);
  void SetCluster(Float_t dpad,Float_t dtime,Float_t charge,Float_t sigmaY2,Float_t sigmaZ2);
  void FillModel();
  void Print();
  
  void SetPadHit(Int_t row,Float_t f) {fPad[row]=f;}
  void SetTimeHit(Int_t row,Float_t f) {fTime[row]=f;}
  void SetOverlap(Int_t i) {fOverlap=i;}
  
  AliL3ClusterModel *GetClusters() {return fClusters;}
  AliL3ClusterModel *GetClusterModel(Int_t i) {if(!fClusters) return 0; return &fClusters[i];}
  AliL3TrackModel *GetModel() {return fTrackModel;}
  Int_t GetOverlap() {return fOverlap;}
  Float_t GetPadHit(Int_t row) {return fPad[row];}
  Float_t GetTimeHit(Int_t row) {return fTime[row];}
  Bool_t GetPadResidual(Int_t row,Float_t &res);
  Bool_t GetTimeResidual(Int_t row,Float_t &res);
  Int_t GetNClusters() {return fNClusters;}
  
  Double_t GetParSigmaY2(Double_t r);
  Double_t GetParSigmaZ2(Double_t r);
  
  ClassDef(AliL3ModelTrack,1)

};

#endif
