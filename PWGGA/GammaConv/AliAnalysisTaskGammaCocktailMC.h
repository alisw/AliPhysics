#ifndef ALIANLYSISTASKGAMMACOCKTAILMC_cxx
#define ALIANLYSISTASKGAMMACOCKTAILMC_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliGenEMCocktailV2.h"
#include "AliMCGenHandler.h"
#include "AliGenerator.h"

class AliAnalysisTaskGammaCocktailMC : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskGammaCocktailMC();
    AliAnalysisTaskGammaCocktailMC(const char *name);
    virtual ~AliAnalysisTaskGammaCocktailMC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void ProcessMCParticles();
    
    // additional functions
    void SetLogBinningXTH1(TH1* histoRebin);
    void SetLogBinningXTH2(TH2* histoRebin);
    void SetMaxY(Double_t maxy){fMaxY = maxy;}
    void SetLightOutput(Bool_t flag) {fDoLightOutput = flag;}
    void InitializeDecayChannelHist(TH1F* hist, Int_t np);
    void FillPythiaBranchingRatio(TH1F* histo, Int_t np);
    void GetAndSetPtParametrizations(AliGenEMCocktailV2* mcCocktailGen);
    void GetAndSetPtYDistributions(AliGenEMCocktailV2* mcCocktailGen);
    void SetHasMother(UInt_t selectedMothers);
    Int_t GetParticlePosLocal(Int_t pdg);
    TH1* SetHist1D(TH1* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Bool_t optSumw2);
    TH2* SetHist2D(TH2* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t yMin, Double_t yMax, Bool_t optSumw2);
    TH2* SetHist2D(TH2* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t* binsY, Bool_t optSumw2);
    Float_t GetDecayChannel(AliStack* stack, TParticle* part);
    
  protected:
    AliVEvent*                  fInputEvent;                    // current event
    AliMCEvent*                 fMCEvent;                       // corresponding MC event
    AliStack*                   fMCStack;                       // stack belonging to MC event
    AliMCGenHandler*            fMCGenHandler;                  // MC gen handler
    const AliGenerator*         fMCGenerator;                   //
    AliGenEMCocktailV2*         fMCCocktailGen;                 // cocktail generator
  
    TList*                      fUserInfo;
    TTree*                      fOutputTree;
    TList*                      fOutputContainer;               // Output container
  
    Int_t*                      fParticleList;                  // array with particle Pdg values
    TString*                    fParticleListNames;             // array with particle names
    
    Bool_t                      fDoLightOutput;                 // switch for running light
    Bool_t                      fHasMother[17];                 // mother i produced
  
    // histograms events
    TH1F*                       fHistNEvents;                   // number of events histo

    // histograms mesons
    TH2F*                       fHistPtYGamma;                  //! histo for gammas
    TH2F**                      fHistPtYInput;                  //! histo for gammas from input particles
    TH2F**                      fHistPtYGammaSource;            //! histo for input particles
    TH2F**                      fHistPtAlphaInput;              //! histo for asymmetry
    TH2F**                      fHistPtDeltaPhiInput;           //! histo for asymmetry
    TH1F**                      fHistDecayChannelsInput;        //! histo for input particle decay channels
    TH1F**                      fHistPythiaBR;                  //! histo for input particle BR from pythia

    Int_t                       fIsMC;                          // MC flag
    Double_t                    fMaxY;                          // Max y
    
    TH2F*                       fHistPtPhiGamma;                //! histo for phi of gammas
    TH2F**                      fHistPtPhiGammaSource;          //! histo for phi of gammas from input particles
    TH2F**                      fHistPtPhiInput;                //! histo for phi of input particles

    TH2F**                      fHistPtGammaSourcePtInput;      //! histo for pt correlation of gammas from input particles to source
    TH2F**                      fHistPhiGammaSourcePhiInput;    //! histo for phi correlation of gammas from input particles to source
    
    TH1I*                       fHistPdgInputRest;              //! histo for rest
    TH1I*                       fHistPdgGammaSourceRest;        //! histo for gamma from rest
  
    TF1*                        fPtParametrization[17];         //!
    TF1*                        fPtParametrizationProton;       //!
    TObjString*                 fCocktailSettings[12];          //!
    TH1D*                       fMtScalingFactors;              //!
    TH2F*                       fPtYDistributions[17];          //!

  private:
    AliAnalysisTaskGammaCocktailMC(const AliAnalysisTaskGammaCocktailMC&);            // Prevent copy-construction
    AliAnalysisTaskGammaCocktailMC &operator=(const AliAnalysisTaskGammaCocktailMC&); // Prevent assignment

    ClassDef(AliAnalysisTaskGammaCocktailMC, 4);
};

#endif
