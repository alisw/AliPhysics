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
enum {kEta05,kEta10,kEta14,kNEtaHist};// 
enum {kPionPos, kProtonPos, kKaonPos, kElectronPos, kMuonPos,
      kPionNeg, kProtonNeg, kKaonNeg, kElectronNeg, kMuonNeg,
      kLambda, kLambdaBar, kLambdaInclusive, kLambdaBarInclusive,
      kNPart}; //Particles used for identified particles pt spectra
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
  TH1F         *fHistPt[kNHist+1]; // Pt spectrum  , |eta| < 0.8
  TH1F         *fHistPtID[kNHist][kNPart+1]; //Pt identified particles, |y| < 0.5 
  
  TGraphErrors *fNchDens; // <dN/deta>
  TList * fMyOut; // list of output histos
  TH1I * fHistIev; // number of events per class
  TH1I * fHistNParticlesAtMidRapidity;  // number of particles at midrapidity per class
  static Float_t fEtaMax; // max eta
  Bool_t fSkipNormalization; // Use this when you are running the job on the grid, so that you can normalize dNdeta after merging

  Float_t  fEtaBins[kNEtaHist];    // array of eta_max values
  TH1F * fHistMult[kNHist][kNEtaHist];   // array of multiplicity histos in the different eta ranges values, for the different event classes

  AliAnalysisTaskdNdetaMC(const AliAnalysisTaskdNdetaMC&); // not implemented
  AliAnalysisTaskdNdetaMC& operator=(const AliAnalysisTaskdNdetaMC&); // not implemented
  
  static Int_t fPDGCodes[kNPart+1]; // array of PDG codes of particles for ID Spectra plots
  static const char *  fPartNames[kNPart+1]; // array of particles names for ID Spectra plots

  TH1F * fHistSpecies; // Histogram of particle species contributing to primaries

  ClassDef(AliAnalysisTaskdNdetaMC, 2); 
};

#endif
