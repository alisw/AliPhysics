#ifndef AliAnalysisTaskC2_cxx
#define AliAnalysisTaskC2_cxx

#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "THn.h"
#include "TRandom3.h"

#include "AliAnalysisTaskC2Base.h"

class AliVEvent;
class AliAODTrack;

class AliAnalysisTaskC2 : public AliAnalysisTaskC2Base {
 public:
  AliAnalysisTaskC2();
  AliAnalysisTaskC2(const char *name);
  virtual ~AliAnalysisTaskC2() {};

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

 private:
  // Histograms to construct C2
  THn *fEventCounter;
  THn *fSingles;
  THn *fPairs;

  // QA Histograms
  TH1F *fmultDistribution;

  // struct to buffer the track informations. Makes it easy to replace tracks with random data
  struct cNano_track {
    Double_t eta;
    Double_t phi;
    Double_t pt;
  };
  struct cEventCounterDims {
    enum type {kMult, kZvtx, kNdimensions};
  };
  struct cPairsDims {
    enum type {kEta1, kEta2, kPhi1, kPhi2, kPtPair, kMult, kZvtx, kNdimensions};
  };
  struct cSinglesDims {
    enum type {kEta, kPhi, kPt, kMult, kZvtx, kNdimensions};
  };


  TRandom3 *fRndmGenerator;

  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskC2(const AliAnalysisTaskC2&); // not implemented
  AliAnalysisTaskC2& operator=(const AliAnalysisTaskC2&); // not implemented

  ClassDef(AliAnalysisTaskC2, 1); // example of analysis
};

#endif
