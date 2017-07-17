#ifndef ALIANLYSISTASKHADRONICCOCKTAILMC_cxx
#define ALIANLYSISTASKHADRONICCOCKTAILMC_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliGenEMCocktailV2.h"
#include "AliMCGenHandler.h"
#include "AliGenerator.h"

class AliAnalysisTaskHadronicCocktailMC : public AliAnalysisTaskSE {
  public:
  
    AliAnalysisTaskHadronicCocktailMC();
    AliAnalysisTaskHadronicCocktailMC(const char *name);
    virtual ~AliAnalysisTaskHadronicCocktailMC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void ProcessMCParticles();

    // setters
    void SetMaxY(Double_t maxy){fMaxY = maxy;}
    void SetLightOutput(Bool_t flag) {fDoLightOutput = flag;}
    void SetAnalyzedParticle(Int_t flag);
    void SetHasMother(UInt_t selectedMothers);
    TH1* SetHist1D(TH1* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Bool_t optSumw2);
    TH2* SetHist2D(TH2* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t yMin, Double_t yMax, Bool_t optSumw2);
    TH2* SetHist2D(TH2* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t* binsY, Bool_t optSumw2);
    void SetLogBinningXTH1(TH1* histoRebin);
    void SetLogBinningXTH2(TH2* histoRebin);

    // getters
    Int_t   GetParticlePosLocal(Int_t pdg);
    Float_t GetDecayChannel(AliMCEvent* mcEvent, TParticle* part);
    void    GetAndSetPtParametrizations(AliGenEMCocktailV2* mcCocktailGen);
    void    GetAndSetPtYDistributions(AliGenEMCocktailV2* mcCocktailGen);

    // additional functions
    void InitializeDecayChannelHist(TH1F* hist, Int_t np);
    void FillPythiaBranchingRatio(TH1F* histo, Int_t np);

  protected:
    AliVEvent*                  fInputEvent;                      // current event
    AliMCEvent*                 fMCEvent;                         // corresponding MC event
    AliMCGenHandler*            fMCGenHandler;
    const AliGenerator*         fMCGenerator;
    AliGenEMCocktailV2*         fMCCocktailGen;

    TList*                      fUserInfo;
    TTree*                      fOutputTree;
    TList*                      fOutputContainer;                 // Output container
  
    Int_t*                      fParticleList;                    // array with particle Pdg values
    TString*                    fParticleListNames;               // array with particle names
  
    Int_t                       fAnalyzedMeson;                   // switch for analyzing pi0 (0), eta (1), pi+-(2)
    Bool_t                      fAnalyzeNeutralPi;                // switch for pi0 analysis
    Bool_t                      fAnalyzeChargedPi;                // switch for pi+- analysis
    Bool_t                      fDoLightOutput;                   // switch for running light
    Bool_t                      fHasMother[24];                   // mother i produced
  
    // nEvent histogram
    TH1F*                       fHistNEvents;                     // number of events histo
  
    // histograms mesons
    TH2F**                      fHistPtYInput;                    //! histo for pi0/eta from input particles
    TH2F**                      fHistPtYDaughterSource;           //! histo for input particles
    TH1F**                      fHistDecayChannelsInput;          //! histo for input particle decay channels
    TH1F**                      fHistPythiaBR;                    //! histo for input particle BR from pythia

    Int_t                       fIsMC;                            // MC flag
    Double_t                    fMaxY;                            // Max y
  
    TH2F**                      fHistPtPhiDaughterSource;         //! histo for phi of pi0/eta from input particles
    TH2F**                      fHistPtPhiInput;                  //! histo for phi of input particles

    TH2F**                      fHistPtYGammaFromXFromInput;      //! gammas from X from k0s, k0l, lambda
    TH2F**                      fHistPtPhiGammaFromXFromInput;    //!
    TH2F**                      fHistPtYGammaFromPi0FromInput;    //! gammas from pi0 from k0s, k0l, lambda
    TH2F**                      fHistPtPhiGammaFromPi0FromInput;  //!
  
    TH2F**                      fHistPtDaughterPtSourceInput;     //! histo for pt correlation of gammas from input particles to source
    TH2F**                      fHistPhiDaughterPhiSourceInput;   //! histo for phi correlation of gammas from input particles to source
  
    TH1I*                       fHistPdgInputRest;                //! histo for rest
    TH1I*                       fHistPdgDaughterSourceRest;       //! histo for gamma from rest
  
    // generator settings
    TF1*                        fPtParametrization[24];           //!
    TF1*                        fPtParametrizationProton;         //!
    TF1*                        fPtParametrizationPi0;            //!
    TObjString*                 fCocktailSettings[12];            //!
    TH1D*                       fMtScalingFactors;                //!
    TH2F*                       fPtYDistributions[24];            //!
  
  private:
    AliAnalysisTaskHadronicCocktailMC(const AliAnalysisTaskHadronicCocktailMC&);              // Prevent copy-construction
    AliAnalysisTaskHadronicCocktailMC &operator=(const AliAnalysisTaskHadronicCocktailMC&);   // Prevent assignment
  
    ClassDef(AliAnalysisTaskHadronicCocktailMC, 8);
};

#endif






  
  
  
  
  
  
