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
  AliITSVertexerZ(TString filename,Float_t x0=0., Float_t y0=0.);
  virtual ~AliITSVertexerZ();
  // The number of contributors set in the AliESDVertex object is the
  // number of tracklets used to determine the vertex position
  // If this number is <1, the procedure could not find a vertex position
  // and by default the Z coordinate is set to 0
  // Number of contributors = -1  --> No tracklets 
  // Number of contributors = -2  --> No SPD recpoints
  virtual AliESDVertex* FindVertexForCurrentEvent(Int_t evnumb);
  virtual void FindVertices();
  virtual void PrintStatus() const;
  void SetDiffPhiMax(Float_t pm = 0.01){fDiffPhiMax = pm;}
  void ConfigIterations(Int_t noiter=3,Float_t *ptr=0);
  void SetFirstLayerModules(Int_t m1 = 0, Int_t m2 = 79){fFirstL1 = m1; fLastL1 = m2;}
  void SetSecondLayerModules(Int_t m1 = 80, Int_t m2 = 239){fFirstL2 = m1; fLastL2 = m2;}
  void SetLowLimit(Float_t lim=-20.){fLowLim = lim;}
  void SetHighLimit(Float_t lim=20.){fHighLim = lim;}
  Float_t GetLowLimit() const {return fLowLim;}
  Float_t GetHighLimit() const {return fHighLim;}
  void SetBinWidthCoarse(Float_t bw=0.01){fStepCoarse = bw;}
  void SetBinWidthFine(Float_t bw=0.0005){fStepFine = bw;}
  void SetPPsetting(Float_t cl2=250., Float_t coarsebin=0.02){fPPsetting[0]=cl2; fPPsetting[1]=coarsebin;}
  Float_t GetBinWidthCoarse() const {return fStepCoarse;}
  Float_t GetBinWidthFine() const {return fStepFine;}
  void SetTolerance(Float_t tol = 20./10000.){fTolerance = tol;}
  Float_t GetTolerance() const {return fTolerance;}
  //  virtual void MakeTracklet(Double_t * /* pA */, Double_t * /*pB */, Int_t & /* nolines */) {} // implemented in a derived class

 protected:
  AliITSVertexerZ(const AliITSVertexerZ& vtxr);
  AliITSVertexerZ& operator=(const AliITSVertexerZ& vtxr );
  void ResetHistograms();
  void VertexZFinder(Int_t evnumber);
  Float_t GetPhiMaxIter(Int_t i) const {return fPhiDiffIter[i];}


  Int_t fFirstL1;          // first module of the first pixel layer used
  Int_t fLastL1;           // last module of the first pixel layer used
  Int_t fFirstL2;          // first module of the second pixel layer used
  Int_t fLastL2;           // last module of the second pixel layer used
  Float_t fDiffPhiMax;     // Maximum delta phi allowed among corr. pixels
  Float_t fX0;             // Nominal x coordinate of the vertex
  Float_t fY0;             // Nominal y coordinate of the vertex
  Float_t fZFound;         //! found value for the current event
  Float_t fZsig;           //! RMS of Z
  TH1F *fZCombc;           //! histogram with coarse z distribution
  TH1F *fZCombv;           //! histogram with very coarse z (step=3coarse)distribution
  TH1F *fZCombf;           //! histogram with fine z distribution
  Float_t fLowLim;         // low limit for fZComb histograms
  Float_t fHighLim;        // high limit for fZComb histograms
  Float_t fStepCoarse;     // bin width for fZCombc
  Float_t fStepFine;       // bin width for fZCombf
  Float_t fTolerance;      // tolerance on the symmetry of the Z interval 
  Float_t fPPsetting[2];   // [0] is the max. number of clusters on L2 to use [1] as fStepCoarse
  Int_t fMaxIter;            // Maximum number of iterations (<=5)
  Float_t fPhiDiffIter[5];   // Delta phi used in iterations

  ClassDef(AliITSVertexerZ,4);
};

#endif
