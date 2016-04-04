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
  AliAnalysisTaskC2(const char *name, Int_t mode);
  virtual ~AliAnalysisTaskC2() {};

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

 private:
  // Histograms to construct C2
  THnF *fEventCounter;
  THnF *fSingle1;
  THnF *fSingle2;
  THnS *fPairs;

  // Histograms to compute efficiency on the fly
  THnF *fSingleParticleForEfficiency;

  // QA Histograms
  THnS *fPairsConventional;
  TH1F *fdNdeta;
  TH1F *fmultDistribution;

  // struct to buffer the track informations. Makes it easy to replace tracks with random data
  struct cNano_track {
    Double_t eta;
    Double_t phi;
    Double_t pt;
  };

  TRandom3 *fRndmGenerator;

  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskC2(const AliAnalysisTaskC2&); // not implemented
  AliAnalysisTaskC2& operator=(const AliAnalysisTaskC2&); // not implemented

  ClassDef(AliAnalysisTaskC2, 1); // example of analysis
};

#endif
