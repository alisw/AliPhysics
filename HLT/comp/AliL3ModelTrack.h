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
  Int_t *fOverlap; //!
  Float_t fXYResidualQ; //Quantization steps.
  Float_t fZResidualQ;
  Float_t fXYResolution;
  Float_t fZResolution;
  Float_t fXYWidthQ;
  Float_t fZWidthQ;
  Int_t fSlice;
  Int_t fPatch;
  Int_t fLabel;
  
  //Crossing points with padrows
  Float_t *fPad; //!
  Float_t *fTime; //!
  
 public:
  AliL3ModelTrack();
  virtual ~AliL3ModelTrack();
  
  void Init(Int_t slice,Int_t patch);
  void SetCluster(Int_t row,Float_t dpad,Float_t dtime,Float_t charge,Float_t sigmaY2,Float_t sigmaZ2,Int_t npads);
  void FillModel();
  void FillTrack();
  void Print();
  void AssignTrackID(Float_t wrong=0.10);
  
  void SetTrackID(Int_t row,Int_t *trackID);
  void SetPadHit(Int_t row,Float_t f);
  void SetTimeHit(Int_t row,Float_t f);
  void SetOverlap(Int_t row,Int_t id);
  void SetXYResolution(Float_t f) {fXYResolution=f;}
  void SetZResolution(Float_t f) {fZResolution=f;}
  void SetLabel(Int_t i) {fLabel = i;}
  Int_t CheckClustersQuality(UInt_t npads=3);
  

  Int_t GetTrackID(Int_t row,Int_t idindex);
  AliL3ClusterModel *GetClusters() {return fClusters;}
  AliL3TrackModel *GetModel() {return fTrackModel;}
  AliL3ClusterModel *GetClusterModel(Int_t row);
  Int_t GetOverlap(Int_t row);
  Int_t GetNPads(Int_t row);
  Float_t GetPadHit(Int_t row);
  Float_t GetTimeHit(Int_t row);
  Bool_t GetPad(Int_t row,Float_t &pad);
  Bool_t GetTime(Int_t row,Float_t &time);
  Bool_t GetClusterCharge(Int_t row,Int_t &charge);
  Bool_t GetXYWidth(Int_t row,Float_t &width);
  Bool_t GetZWidth(Int_t row,Float_t &width);
  Bool_t GetPadResidual(Int_t row,Float_t &res);
  Bool_t GetTimeResidual(Int_t row,Float_t &res);
  Bool_t GetXYWidthResidual(Int_t row,Float_t &res);
  Bool_t GetZWidthResidual(Int_t row,Float_t &res);
  Int_t GetNClusters() {return fNClusters;}
  Int_t GetLabel() {return fLabel;}

  Double_t GetParSigmaY2(Int_t row);
  Double_t GetParSigmaZ2(Int_t row);
  
  ClassDef(AliL3ModelTrack,1)

};

#endif
