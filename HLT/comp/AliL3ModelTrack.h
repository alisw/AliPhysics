// @(#) $Id$

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
  Int_t fMaxOverlaps;
  Int_t *fNoverlaps; //!
  Int_t **fOverlap; //!
  Float_t *fParSigmaY2;    //!
  Float_t *fParSigmaZ2;    //!
  Float_t *fCrossingAngle; //!
  Int_t fPatch;

  //Crossing points with padrows
  Float_t *fPad; //!
  Float_t *fTime; //!
  
  Float_t QuantizePad(Int_t row,Float_t pad);
  Float_t QuantizeTime(Int_t row,Float_t time);
  Float_t QuantizeSigmaY2(Int_t row,Float_t dsigmaY2);
  Float_t QuantizeSigmaZ2(Int_t row,Float_t dsigmaZ2);
  Float_t RetrievePad(Int_t row,Float_t dpad);
  Float_t RetrieveTime(Int_t row,Float_t time);
  Float_t RetrieveSigmaY2(Int_t row,Float_t dsigmaY2);
  Float_t RetrieveSigmaZ2(Int_t row,Float_t dsigmaZ2);
  
 public:
  AliL3ModelTrack();
  virtual ~AliL3ModelTrack();
  
  void Init(Int_t slice,Int_t patch);
  void CalculateClusterWidths(Int_t row,Bool_t parametrize=kFALSE);
  void SetCluster(Int_t row,Float_t dpad,Float_t dtime,Float_t charge,Float_t sigmaY2,Float_t sigmaZ2,Int_t npads);
  void FillModel();
  void FillTrack();
  void Print(Bool_t everything=kTRUE);
  void Set(AliL3Track *tpt);

  void SetPadHit(Int_t row,Float_t f);
  void SetTimeHit(Int_t row,Float_t f);
  void SetCrossingAngleLUT(Int_t row,Float_t angle);
  void SetOverlap(Int_t row,Int_t id);
  void SetClusterLabel(Int_t row,Int_t *trackID);
  void SetNClusters(Int_t i) {fNClusters = i;}
  
  Int_t GetNPresentClusters();
  Bool_t IsPresent(Int_t row);
  Bool_t IsSet(Int_t row);
  
  AliL3ClusterModel *GetClusters() {return fClusters;}
  AliL3TrackModel *GetModel() {return fTrackModel;}
  AliL3ClusterModel *GetClusterModel(Int_t row);
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
  Bool_t GetXYWidth(Int_t row,Float_t &width);
  Bool_t GetZWidth(Int_t row,Float_t &width);
  Bool_t GetPadResidual(Int_t row,Float_t &res);
  Bool_t GetTimeResidual(Int_t row,Float_t &res);
  Bool_t GetXYWidthResidual(Int_t row,Float_t &res);
  Bool_t GetZWidthResidual(Int_t row,Float_t &res);
  Int_t GetNClusters() {return fNClusters;}
  void GetClusterLabel(Int_t row,Int_t *trackID);
    
  ClassDef(AliL3ModelTrack,1)

};

#endif
