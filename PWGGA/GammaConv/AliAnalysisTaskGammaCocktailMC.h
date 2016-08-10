#ifndef ALIANLYSISTASKGAMMACOCKTAILMC_cxx
#define ALIANLYSISTASKGAMMACOCKTAILMC_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

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
    void InitializeDecayChannelHist(Int_t inputParticle);
    Float_t GetDecayChannel(AliStack* stack, TParticle* part);
    Bool_t IsMotherInList(TParticle* mother);
    
  protected:
    AliVEvent*            fInputEvent;                // current event
    AliMCEvent*           fMCEvent;                   // corresponding MC event
    AliStack*             fMCStack;                   // stack belonging to MC event
    TList*                fOutputContainer;           // Output container
    
    Int_t*                fParticleList;              // array with particle Pdg values
    TString*              fParticleListNames;         // array with particle names
    
    Bool_t                fDoLightOutput;             // switch for running light
    
    // histograms events
    TH1F*                 fHistNEvents;               // number of events histo

    // histograms mesons
    TH2F*                 fHistPtYGamma;              //! histo for gammas
    TH2F**                fHistPtYInput;              //! histo for gammas from input particles
    TH2F**                fHistPtYGammaSource;        //! histo for input particles
    TH2F**                fHistPtAlphaInput;          //! histo for asymmetry
    TH2F**                fHistPtDeltaPhiInput;       //! histo for asymmetry
    TH1F**                fHistDecayChannelsInput;    //! histo for input particle decay channels
    Int_t                 fIsMC;                      // MC flag
    Double_t              fMaxY;                      // Max y
    
    TH2F*                 fHistPtPhiGamma;            //! histo for phi of gammas
    TH2F**                fHistPtPhiGammaSource;      //! histo for phi of gammas from input particles
    TH2F**                fHistPtPhiInput;            //! histo for phi of input particles
    
    TH2F**                fHistPtGammaSourceInput;    //! histo for pt correlation of gammas from input particles to source
    TH2F**                fHistPhiGammaSourceInput;   //! histo for phi correlation of gammas from input particles to source
    
    TH1I*                 fHistPtYInputRest;          //! histo for rest
    TH1I*                 fHistPtYGammaSourceRest;    //! histo for gamma from rest

  private:
    AliAnalysisTaskGammaCocktailMC(const AliAnalysisTaskGammaCocktailMC&); // Prevent copy-construction
    AliAnalysisTaskGammaCocktailMC &operator=(const AliAnalysisTaskGammaCocktailMC&); // Prevent assignment

    ClassDef(AliAnalysisTaskGammaCocktailMC, 1);
};

#endif
