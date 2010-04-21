#ifndef AliAnalysisTaskdNdetaMC_H
#define AliAnalysisTaskdNdetaMC_H

//
// Task used to analize simulations at generation level (i.e. only
// needs galice.root and Kinematics.root).
// 

class TH1F;
class TH1I;
class TGraphErrors;

enum {kHistINEL,kHistNSD,kHistND,kHistSiD,kHistHL,kNHist};
enum {kMult05,kMult10,kMult14,kNMultHist};// 
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskdNdetaMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskdNdetaMC();
  AliAnalysisTaskdNdetaMC(const char *name );
  AliAnalysisTaskdNdetaMC(const char *name, const char *file );
  virtual ~AliAnalysisTaskdNdetaMC();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  TH1F* BookHetaHist(const char * name, const char * title);
  TH1F* BookHptHist(const char * name, const char * title);
  TH1F* BookMultHisto(const char * name, const char * title) ;
  void SetEtaMax(Float_t eta) {fEtaMax = eta;}

  void SkipNormalization(Bool_t flag = kTRUE) { fSkipNormalization = flag; }
  void Finalize();
  TList * GetList() const { return fMyOut;} 
 private:
  TH1F         *fHistEta[kNHist]; //Eta spectrum 
  TH1F         *fHistPt[kNHist]; //Eta spectrum  , |eta| < 0.8
  TGraphErrors *fNchDens; // <dN/deta>
  TList * fMyOut; // list of output histos
  TH1I * fHistIev; // number of events per class
  TH1I * fHistNParticlesAtMidRapidity;  // number of particles at midrapidity per class
  static Float_t fEtaMax; // max eta
  Bool_t fSkipNormalization; // Use this when you are running the job on the grid, so that you can normalize dNdeta after merging

  Float_t  fEtaBins[kNMultHist];    // array of eta_max values
  TH1F * fHistMult[kNHist][kNMultHist];   // array of multiplicity histos in the different eta ranges values, for the different event classes

  AliAnalysisTaskdNdetaMC(const AliAnalysisTaskdNdetaMC&); // not implemented
  AliAnalysisTaskdNdetaMC& operator=(const AliAnalysisTaskdNdetaMC&); // not implemented
  
  ClassDef(AliAnalysisTaskdNdetaMC, 2); 
};

#endif
