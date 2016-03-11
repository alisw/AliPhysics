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
    Bool_t IsMotherInList(TParticle* mother);
    
  protected:
    AliVEvent*            fInputEvent;                // current event
    AliMCEvent*           fMCEvent;                   // corresponding MC event
    AliStack*             fMCStack;                   // stack belonging to MC event
    TList*                fOutputContainer;           // Output container
    
    Int_t*                fParticleList;              // array with particle Pdg values
    TString*              fParticleListNames;         // array with particle names
    
    // histograms events
    TH1F*                 fHistNEvents;               // number of events histo

    // histograms mesons
    TH2F*                 fHistPtYGamma;              //! histo for Gammas
    TH2F**                fHistPtYInput;              //! histo for Input particles
    TH2F**                fHistPtYGammaSource;        //! histo for Input particles
    Int_t                 fIsMC;                      // MC flag
    Double_t              fMaxY;                      // Max y
    
    TH2F*                 fHistPhiGamma;              //! histo for phi of gamma
    TH2F**                fHistPhiInput;              //! histo for phi of input particles
    
    TH1I*                 fHistPtYInputRest;          //! histo for rest
    TH1I*                 fHistPtYGammaSourceRest;    //! histo for gamma from rest

  private:
    AliAnalysisTaskGammaCocktailMC(const AliAnalysisTaskGammaCocktailMC&); // Prevent copy-construction
    AliAnalysisTaskGammaCocktailMC &operator=(const AliAnalysisTaskGammaCocktailMC&); // Prevent assignment

    ClassDef(AliAnalysisTaskGammaCocktailMC, 1);
};

#endif
