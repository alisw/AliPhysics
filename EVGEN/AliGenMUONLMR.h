#ifndef AliGenMUONLMR_h
#define AliGenMUONLMR_h

#include <TH1F.h> 
#include <TH1D.h> 
#include <TF1.h> 
#include <TParticle.h> 
#include <TLorentzVector.h> 
#include "AliGenMC.h" 
 
class AliGenMUONLMR : public AliGenMC { 
 public:
  enum parttype_t {kPionLMR, kKaonLMR, kEtaLMR, kRhoLMR, kOmegaLMR, kPhiLMR, kEtaPrimeLMR};
  AliGenMUONLMR(); 
  ~AliGenMUONLMR(); 
  static Double_t PtDistr(const Double_t *x, const Double_t *par); 
  static Double_t YDistr(const Double_t *x, const Double_t *par); 
  virtual void Decay2Body(const TParticle *mother);
  virtual void DalitzDecay(const TParticle *mother);
  virtual void DecayPiK(TParticle *mother, Bool_t &hadDecayed);
  virtual Double_t FormFactor(Double_t q2, Int_t decay); 
  virtual void Generate(); 
  virtual TParticle* GetMuon(Int_t i) {return fMu[i];} 
  virtual void SetNMuMin(Int_t nmin) {fNMuMin = nmin; }
  virtual void GenerateSingleProcess(Int_t whichproc) { fGenSingleProc = whichproc;}
  virtual void SetScaleMultiplicity(Int_t ipart, Double_t scale) { fScaleMult[ipart] = scale; } 
  static Double_t RhoLineShapeNew(const Double_t *, const Double_t *); 
  virtual void FinishRun(); 
 private:
  AliGenMUONLMR(const AliGenMUONLMR &lmr);
  AliGenMUONLMR & operator=(const AliGenMUONLMR & rhs);

 private: 
  static const Int_t fgkNpart = 7; // number of particles to be generated 
  Int_t fNMuMin;                   // min. number of muons to accept the event for writing
  Int_t fGenSingleProc;            // flag to generate a single process (1) or the whole cocktail (0)
  Int_t fPDG[7];                   // pdg code of particle to be generated 
  Double_t fScaleMult[7];          // multiplicity scaling factor (w.r.t. pythia@7TeV)
  TF1 *fPt[7];                     // pt distribution
  TF1 *fY[7];                      // rapidity distribution
  TF1 *fMult[7];                   // multiplicity distribution 
  TF1 *fDecay[2];                  // fDecay[0] = pion, fDecay[1] = kaon
  TH1F *fDalitz[3];                // Dalitz decay form factor for eta, omega, etaprime
  TF1 *fCosTheta;                  // function for polarized theta distributions
  TF1 *fRhoLineShape;              // rho line shape 
  TParticle* fParticle[7];         // TPaticle object for the particles to be generated
  TParticle* fMu[2];               // fMu[0] = mu+    fMu[1] = mu-
  TH1D *fHMultMu;                  // muon multiplicity 
  TH1D *fHNProc;                   // number of events generated per process
  ClassDef(AliGenMUONLMR, 1)  // example of analysis
}; 

#endif
