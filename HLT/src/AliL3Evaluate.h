// @(#) $Id$

#ifndef ALIL3_Evaluate
#define ALIL3_Evaluate


class TClonesArray;

#include <TObject.h>
#include <TH1.h>
#include <TTree.h>
#include <TNtuple.h>
#include "AliSimDigits.h"


struct GoodTrack 
{
  Int_t label;
  Double_t eta;
  Int_t code;
  Double_t px,py,pz;
  Double_t x,y,z;
  Int_t nhits;
  Int_t sector;
};
typedef struct GoodTrack GoodTrack;


class AliL3Track;
class TClonesArray;
class TFile;
class AliL3TrackArray;
class AliL3SpacePointData;
class TH1F;
class AliTPCParam;
class TTree;
class AliSimDigits;
class TObjArray;
class TParticle;
class AliL3FileHandler;

class AliL3Evaluate {

 private:

  AliL3TrackArray *fTracks; //!
  //AliTPCParam *fParam;
  AliL3SpacePointData *fClusters[36][6]; //!
  AliL3FileHandler *fClustersFile[36][6]; //!
  //TTree *fDigitsTree;
  //AliSimDigits *fDigits;
  Char_t fPath[1024];
  Int_t fMinSlice;
  Int_t fMaxSlice;
  UInt_t fNcl[36][6];
  Int_t fRowid[36][176];
  Int_t fMinPointsOnTrack;  //Minimum points on track to be considered.
  Int_t fMinHitsFromParticle;
  GoodTrack *fGoodTracks; //!
  Float_t fMaxFalseClusters;
  
  Int_t fNFastPoints;
  UInt_t *fMcIndex;//!
  Int_t *fMcId;//!
  Int_t fGoodFound;
  Int_t fGoodGen;
  Double_t fMinGoodPt;
  Double_t fMaxGoodPt;
  
  //Histograms
  TNtuple *fNtuppel;
  TH1F *fPtRes;
  TH1F *fNGoodTracksPt;
  TH1F *fNFoundTracksPt;
  TH1F *fNFakeTracksPt;
  TH1F *fTrackEffPt;
  TH1F *fFakeTrackEffPt;
  TH1F *fNGoodTracksEta;
  TH1F *fNFoundTracksEta;
  TH1F *fNFakeTracksEta;
  TH1F *fTrackEffEta;
  TH1F *fFakeTrackEffEta;
  TNtuple *fNtupleRes;
  Bool_t fStandardComparison;
  
  void Clear();
 public:
  AliL3Evaluate();
  AliL3Evaluate(Char_t *path,Int_t min_clusters,Int_t minhits,Double_t minpt=0.1,Double_t maxpt=4.,Int_t *slice=0);
  virtual ~AliL3Evaluate();
  
  void LoadData(Int_t event=-1,Bool_t sp=kFALSE);
  void CreateHistos(Int_t nbin=20,Float_t xlow=0,Float_t xup=4);
  void Write2File(Char_t *outputfile);
  void FillEffHistos();
  void FillEffHistosNAIVE();
  void CalcEffHistos();
  void AssignPIDs();
  void AssignIDs();
  void GetGoodParticles(Char_t *particle_file,Int_t event=-1,Int_t *padrowrange=0);
  void GetFastClusterIDs(Char_t *path);
  void GetCFeff(Char_t *path,Char_t *outfile,Int_t nevent=0,Bool_t sp=kFALSE);
  Int_t GetMCTrackLabel(AliL3Track *track);
  Float_t GetTrackPID(AliL3Track *track);
  void CalculateResiduals();
  void EvaluatePoints(Char_t *rootfile,Char_t *exactfile,Char_t *tofile,Int_t nevent=1,Bool_t offline=kFALSE,Bool_t sp=kFALSE);
  Float_t GetCrossingAngle(TParticle *part,Int_t slice,Int_t padrow,Float_t *xyz);
  Int_t FindPrimaries(Int_t nparticles);
  void SetMinPoints(Int_t f) {fMinPointsOnTrack = f;}
  void SetMinGoodPt(Double_t f) {fMinGoodPt = f;}
  void SetMaxFalseClusters(Float_t f) {fMaxFalseClusters = f;}
  void SetStandardComparison(Bool_t b) {fStandardComparison = b;}  
  TNtuple *GetNtuple();
  Int_t GetNGoodTracks() {return fGoodGen;}
  Int_t GetNFoundTracks() {return fGoodFound;}
  TH1F *GetTrackEffPt() {return fTrackEffPt;}
  TH1F *GetTrackEffEta() {return fTrackEffEta;}
  TH1F *GetFakeEffEta() {return fFakeTrackEffEta;}
  TH1F *GetFakeEffPt() {return fFakeTrackEffPt;}
  TH1F *GetPtRes() {return fPtRes;}
  AliL3TrackArray *GetTracks() {return fTracks;}
  
  ClassDef(AliL3Evaluate,1) //Tracking evaluation class
};

#endif
