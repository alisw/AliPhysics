#ifndef ALIL3_Evaluate
#define ALIL3_Evaluate

#include <TObject.h>
#include <TH1.h>
#include <TTree.h>
#include "AliSimDigits.h"

class AliL3Track;
class TClonesArray;
class TFile;
class AliL3TrackArray;
class AliL3SpacePointData;
class TH1F;
class AliL3Transform;
class AliTPCParam;
class TTree;
class AliSimDigits;
class TObjArray;

class AliL3Evaluate : public TObject {

 private:

  TFile *fMCFile;
  TFile *fMCclusterfile;  //If you run the fast simulator.
  TClonesArray *fParticles;
  AliL3TrackArray *fTracks; //!
  AliTPCParam *fParam;
  AliL3SpacePointData *fClusters[36][5]; //!
  AliL3Transform *fTransform; //!
  AliSimDigits *fDigits;
  TTree *fDigitsTree;
  Int_t fMinSlice;
  Int_t fMaxSlice;
  UInt_t fNcl[36][5];
  Int_t fRowid[36][174];
  Int_t fMinPointsOnTrack;  //Minimum points on track to be considered.
  Bool_t fIsSlow;
  Int_t fNFastPoints;
  UInt_t fMcIndex;
  Int_t fMcId;
  
  //Histograms
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
  
  void FillEffHistos(TObjArray *good_particles,Int_t *particle_id);
  void AssignIDs();
  Bool_t SetDigitsTree();
  Bool_t SetMCParticleArray();
  TObjArray *DefineGoodTracks(Int_t slice,Int_t *padrow,Int_t good_number,Int_t *particle_id);
  Int_t GetMCTrackLabel(AliL3Track *track,UInt_t *index=0,Int_t *pID=0,Int_t npoints=0);
  Int_t **GetClusterIDs(AliL3Track *track,UInt_t *index=0,Int_t *pID=0,Int_t npoints=0);
  Int_t *GetFastIDs(UInt_t &tmp_ind,Int_t &npoints);
 
 public:
  AliL3Evaluate();
  AliL3Evaluate(Char_t *mcfile,Int_t *slice);

  virtual ~AliL3Evaluate();

  void Setup(Char_t *trackfile,Char_t *clustfile);
  void SetupFast(Char_t *trackfile,Char_t *clustfile,Char_t *mcClusterfile);
  void CreateHistos(Int_t nbin=20,Int_t xlow=0,Int_t xup=4);
  void EvaluatePatch(Int_t slice,Int_t patch,Int_t min_points,Int_t good_number);
  void EvaluateSlice(Int_t slice,Int_t min_points,Int_t good_number);
  void EvaluateGlobal();
  void Write2File(Char_t *outputfile);
    
  TH1F *GetTrackEffPt() {return fTrackEffPt;}
  TH1F *GetTrackEffEta() {return fTrackEffEta;}
  TH1F *GetPtRes() {return fPtRes;}

  void SetMinPoints(Int_t f) {fMinPointsOnTrack = f;}

  ClassDef(AliL3Evaluate,1) 
};

#endif
