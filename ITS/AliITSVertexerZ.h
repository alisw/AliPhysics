#ifndef ALIITSVERTEXERZ_H
#define ALIITSVERTEXERZ_H

#include<AliITSVertexer.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for primary vertex finding                              //
//                                                               //
///////////////////////////////////////////////////////////////////

class TFile;
class TString;
class AliESDVertex;
class TH1F;

class AliITSVertexerZ : public AliITSVertexer {

 public:

  AliITSVertexerZ();
  AliITSVertexerZ(Float_t x0, Float_t y0);
  virtual ~AliITSVertexerZ();
  // The number of contributors set in the AliESDVertex object is the
  // number of tracklets used to determine the vertex position
  // If this number is <1, the procedure could not find a vertex position
  // and by default the Z coordinate is set to 0
  // Number of contributors = -1  --> No tracklets 
  // Number of contributors = -2  --> No SPD recpoints
  virtual AliESDVertex* FindVertexForCurrentEvent(TTree *itsClusterTree);
  virtual void PrintStatus() const;
  void SetDiffPhiMax(Float_t pm = 0.01){fDiffPhiMax = pm;}
  void ConfigIterations(Int_t noiter=4,Float_t *ptr=0);
  void SetFirstLayerModules(Int_t m1 = 0, Int_t m2 = 79){fFirstL1 = m1; fLastL1 = m2;}
  void SetSecondLayerModules(Int_t m1 = 80, Int_t m2 = 239){fFirstL2 = m1; fLastL2 = m2;}
  void SetLowLimit(Float_t lim=-40.){fLowLim = lim;}
  void SetHighLimit(Float_t lim=40.){fHighLim = lim;}
  Float_t GetLowLimit() const {return fLowLim;}
  Float_t GetHighLimit() const {return fHighLim;}
  void SetBinWidthCoarse(Float_t bw=0.01){fStepCoarse = bw;}
  void SetPPsetting(Float_t cl2=250., Float_t coarsebin=0.02){fPPsetting[0]=cl2; fPPsetting[1]=coarsebin;}
  static Int_t GetPeakRegion(TH1F* h, Int_t &binmin, Int_t &binmax);
  static Int_t FindSecondPeak(TH1F* h, Int_t binmin,Int_t binmax, Float_t& secPeakPos);
  Float_t GetBinWidthCoarse() const {return fStepCoarse;}
  void SetTolerance(Float_t tol = 20./10000.){fTolerance = tol;}
  void SetWindowWidth(Float_t ww=0.2){fWindowWidth=ww;}
  Float_t GetTolerance() const {return fTolerance;}
  //  virtual void MakeTracklet(Double_t * /* pA */, Double_t * /*pB */, Int_t & /* nolines */) {} // implemented in a derived class

  void SetSearchForPileup(Bool_t opt){fSearchForPileup=opt;}
  Bool_t IsSearchForPileupActive() const { return fSearchForPileup;}

 protected:
  void ResetHistograms();
  void VertexZFinder(TTree *itsClusterTree);
  Float_t GetPhiMaxIter(Int_t i) const {return fPhiDiffIter[i];}


  Int_t fFirstL1;          // first module of the first pixel layer used
  Int_t fLastL1;           // last module of the first pixel layer used
  Int_t fFirstL2;          // first module of the second pixel layer used
  Int_t fLastL2;           // last module of the second pixel layer used
  Float_t fDiffPhiMax;     // Maximum delta phi allowed among corr. pixels
  Float_t fZFound;         //! found value for the current event
  Float_t fZsig;           //! RMS of Z
  TH1F *fZCombc;           //! histogram with coarse z distribution
  Float_t fLowLim;         // low limit for fZComb histograms
  Float_t fHighLim;        // high limit for fZComb histograms
  Float_t fStepCoarse;     // bin width for fZCombc
  Float_t fTolerance;      // tolerance on the symmetry of the Z interval 
  Float_t fPPsetting[2];   // [0] is the max. number of clusters on L2 to use [1] as fStepCoarse
  Int_t fMaxIter;            // Maximum number of iterations (<=5)
  Float_t fPhiDiffIter[5];   // Delta phi used in iterations
  Float_t fWindowWidth;      // Z window width for symmetrization
  Bool_t  fSearchForPileup;  // flag to switch pileup off/on

 private:
  AliITSVertexerZ(const AliITSVertexerZ& vtxr);
  AliITSVertexerZ& operator=(const AliITSVertexerZ& vtxr );

  ClassDef(AliITSVertexerZ,11);
};

#endif
