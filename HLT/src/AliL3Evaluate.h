#ifndef ALIL3_Evaluate
#define ALIL3_Evaluate

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

class AliL3Evaluate {

 private:

  AliL3TrackArray *fTracks; //!
  AliTPCParam *fParam;
  AliL3SpacePointData *fClusters[36][6]; //!
  AliSimDigits *fDigits;
  TTree *fDigitsTree;
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
  
  
  
 public:
  AliL3Evaluate();
  AliL3Evaluate(Char_t *path,Int_t min_clusters,Int_t minhits,Double_t minpt=0.1,Int_t *slice=0);
  virtual ~AliL3Evaluate();

  void CreateHistos(Int_t nbin=20,Float_t xlow=0,Float_t xup=4);
  void Write2File(Char_t *outputfile);
  void FillEffHistos();
  void FillEffHistosNAIVE();
  void CalcEffHistos();
  void AssignIDs();
  void GetGoodParticles(Char_t *particle_file);
  void GetFastClusterIDs(Char_t *path);
  void GetCFeff(Char_t *outfile);
  Int_t GetMCTrackLabel(AliL3Track *track);
  TNtuple *CalculateResiduals();
  TNtuple *EvaluatePoints(Char_t *rootfile);
  
  void SetMinPoints(Int_t f) {fMinPointsOnTrack = f;}
  void SetMinGoodPt(Double_t f) {fMinGoodPt = f;}
  void SetMaxFalseClusters(Float_t f) {fMaxFalseClusters = f;}

  Int_t GetNGoodTracks() {return fGoodGen;}
  Int_t GetNFoundTracks() {return fGoodFound;}
  TH1F *GetTrackEffPt() {return fTrackEffPt;}
  TH1F *GetTrackEffEta() {return fTrackEffEta;}
  TH1F *GetFakeEffEta() {return fFakeTrackEffEta;}
  TH1F *GetFakeEffPt() {return fFakeTrackEffPt;}
  TH1F *GetPtRes() {return fPtRes;}

  
  ClassDef(AliL3Evaluate,1) //Tracking evaluation class
};

#endif
