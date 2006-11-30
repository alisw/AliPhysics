// @(#) $Id$

#ifndef ALIL3_Evaluate
#define ALIL3_Evaluate


class AliHLTTrack;
class TClonesArray;
class TFile;
class AliHLTTrackArray;
class AliHLTSpacePointData;
class TH1F;
class AliTPCParam;
class TTree;
class AliSimDigits;
class TObjArray;
class TParticle;
class AliHLTFileHandler;
class TClonesArray;
class TNtuple;
class TH1F;
//#include <TNtuple.h>
//#include <TH1F.h>

class AliHLTEvaluate {

struct AliGoodTrack 
{
  Int_t flabel;  //label
  Double_t feta; //eta
  Int_t fcode;   //pcode
  Double_t fpx,fpy,fpz; //momentum
  Double_t fx,fy,fz;    //pos at entrance 
  Int_t fnhits;  //nhits
  Int_t fsector; //sector number
};
typedef struct AliGoodTrack AliGoodTrack;

struct AliS 
{
  Int_t flab; //lab
  Int_t fmax; //max
};
typedef struct AliS AliS;

 private:

  AliHLTTrackArray *fTracks; //!
  AliHLTSpacePointData *fClusters[36][6]; //!
  AliHLTFileHandler *fClustersFile[36][6]; //!
  Char_t fPath[1024];      //path
  Int_t fMinSlice;         //min slice
  Int_t fMaxSlice;         //max slice
  UInt_t fNcl[36][6];      //cluster numbers
  Int_t fRowid[36][176];   //row ids
  Int_t fMinPointsOnTrack;    //minimum points on track to be considered.
  Int_t fMinHitsFromParticle; //minimums hits a particle has to create
  AliGoodTrack *fGoodTracks; //!
  Float_t fMaxFalseClusters;  //maximum number of false assigned clusters
  
  Int_t fNFastPoints; //fast access to points
  UInt_t *fMcIndex;//!
  Int_t *fMcId;//!
  Int_t fGoodFound;    //good found
  Int_t fGoodGen;      //good generated found
  Double_t fMinGoodPt; //min pt
  Double_t fMaxGoodPt; //max pt
  
  //Histograms
  TNtuple *fNtuppel;//!
  TH1F *fPtRes;//!
  TH1F *fNGoodTracksPt;//!
  TH1F *fNFoundTracksPt;//!
  TH1F *fNFakeTracksPt;//!
  TH1F *fTrackEffPt;//!
  TH1F *fFakeTrackEffPt;//!
  TH1F *fNGoodTracksEta;//!
  TH1F *fNFoundTracksEta;//!
  TH1F *fNFakeTracksEta;//!
  TH1F *fTrackEffEta;//!
  TH1F *fFakeTrackEffEta;//!
  TNtuple *fNtupleRes;//!
  Bool_t fStandardComparison; // take standard macro
  
  void Clear();
 public:
  AliHLTEvaluate();
  AliHLTEvaluate(Char_t *path,Int_t min_clusters,Int_t minhits,Double_t minpt=0.1,Double_t maxpt=4.,Int_t *slice=0);
  virtual ~AliHLTEvaluate();
  
  void LoadData(Int_t event=-1,Bool_t sp=kFALSE);
  void CreateHistos(Int_t nbin=20,Float_t xlow=0,Float_t xup=4);
  void Write2File(Char_t *outputfile);
  void FillEffHistos();
  void FillEffHistosNAIVE();
  void CalcEffHistos();
  void AssignPIDs();
  void AssignIDs();
  void GetGoodParticles(Char_t *particlefile,Int_t event=-1,Int_t *padrowrange=0);
  void GetFastClusterIDs(Char_t *path);
  void GetCFeff(Char_t *path,Char_t *outfile,Int_t nevent=0,Bool_t sp=kFALSE);
  Int_t GetMCTrackLabel(AliHLTTrack *track);
  Float_t GetTrackPID(AliHLTTrack *track);
  void CalculateResiduals();
  void EvaluatePoints(Char_t *rootfile,Char_t *exactfile,Char_t *tofile,Int_t nevent=1,Bool_t offline=kFALSE,Bool_t sp=kFALSE);
  Float_t GetCrossingAngle(TParticle *part,Int_t slice,Int_t padrow,Float_t *xyz);
  Int_t FindPrimaries(Int_t nparticles);
  void SetMinPoints(Int_t f) {fMinPointsOnTrack = f;}
  void SetMinGoodPt(Double_t f) {fMinGoodPt = f;}
  void SetMaxFalseClusters(Float_t f) {fMaxFalseClusters = f;}
  void SetStandardComparison(Bool_t b) {fStandardComparison = b;}  
  TNtuple *GetNtuple();
  Int_t GetNGoodTracks() const {return fGoodGen;}
  Int_t GetNFoundTracks() const {return fGoodFound;}
  TH1F *GetTrackEffPt() {return fTrackEffPt;}
  TH1F *GetTrackEffEta() {return fTrackEffEta;}
  TH1F *GetFakeEffEta() {return fFakeTrackEffEta;}
  TH1F *GetFakeEffPt() {return fFakeTrackEffPt;}
  TH1F *GetPtRes() {return fPtRes;}
  AliHLTTrackArray *GetTracks() {return fTracks;}
  
  ClassDef(AliHLTEvaluate,1) //Tracking evaluation class
};

typedef AliHLTEvaluate AliL3Evaluate; // for backward compatibility

#endif
