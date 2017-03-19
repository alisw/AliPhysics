#ifndef ALIANALYSISTASKDIJET_H
#define ALIANALYSISTASKDIJET_H
class TH1;
class TH1D;
class TH2;
class TH3;
class TF1;
class TRandom3;
class TGraphErrors;
class TProfile;
class THistManager;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliEmcalTrackSelection;

#include "AliAnalysisTaskEmcalJet.h"
#include <vector>

class AliAnalysisTaskDijet : public AliAnalysisTaskEmcalJet {
public:
  typedef std::vector<TLorentzVector>   TLorentzVector1D;
  typedef std::vector<TLorentzVector1D> TLorentzVector2D;
  typedef std::vector<Double_t> Double1D; 
  enum { kBTypeBegin=1,kBData=1, kBMixingSimple, kBTypeEnd };
  enum { kBDiJetSelBegin=1, kBLS=1, kBLS_PI2, kBL_PI2_S, kBDiJetSelEnd }; 
  AliAnalysisTaskDijet();
  AliAnalysisTaskDijet(const char *name);
  AliAnalysisTaskDijet(const char *name, const char *option);
  AliAnalysisTaskDijet(const AliAnalysisTaskDijet& ap);
  AliAnalysisTaskDijet& operator =(const AliAnalysisTaskDijet& ap);
  ~AliAnalysisTaskDijet();
  virtual void UserCreateOutputObjects();
  void SetIsAA(bool is=true){ fIsAA = is; }
  void SetLeadingParticlePtMin(Double_t m){ fLeadingParticlePtMin=m; };
  void SetIsMC(Bool_t ismc){fIsMC = ismc;};
  Double1D& GetDijetPtPair(){return fDijetPtPair;};
  Double1D& GetDijetInvM(){return fDijetInvM;};
  Double1D& GetDijetapt(){return fDijetapt;};
protected:
  Bool_t       Run();
  TAxis AxisFix( TString name, int nbin, Double_t xmin, Double_t xmax);
  TAxis AxisVar( TString name, std::vector<Double_t> bin );
  TAxis AxisLog( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0);
  TAxis AxisStr( TString name, std::vector<TString> bin );
  THnSparse * CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt="");
  THnSparse * CreateTHnSparse(TString name, TString title, TString templ, Option_t * opt="");
  Long64_t FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w=1. );
  Long64_t FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w=1. );

private:
  TList*                          fOutput; //!
  TString                         fOption=""; 
  AliJetContainer*                fJetsCont; //!
  UInt_t                          fFilterBit=0;
  Bool_t                          IsFirstEvent=1;
  TH1D*                           fhJetPt; //!
  THistManager*                   fHistos=nullptr; //!
  Bool_t                          fIsAA=0; 
  Int_t                           fNDiJetSelection=0; 
  Int_t                           fNType=0;
  Bool_t                          fIsMC=kFALSE; 
  Double_t                        fLeadingParticlePtMin=0; 
  TRandom3*                       fBSRandom=nullptr; //!
  Double1D                        fDijetPtPair;
  Double1D                        fDijetInvM;
  Double1D                        fDijetapt;
  TF1                             *invmscale=nullptr;//!
  ClassDef(AliAnalysisTaskDijet, 1)
};
#endif
