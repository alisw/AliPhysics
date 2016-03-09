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
    
  protected:
    AliVEvent*            fInputEvent;                // current event
    AliMCEvent*           fMCEvent;                   // corresponding MC event
    AliStack*             fMCStack;                   // stack belonging to MC event
    TList*                fOutputContainer;           // Output container
    // histograms events
    TH1F*                 fHistNEvents;               // number of events histo

    // histograms mesons
    TH2F*                 fHistPtYPi0;                //! histo for Pi0s
    TH2F*                 fHistPtYPiPl;               //! histo for Pi+s 
    TH2F*                 fHistPtYPiMi;               //! histo for Pi-s
    TH2F*                 fHistPtYEta;                //! histo for Etas
    TH2F*                 fHistPtYEtaPrim;            //! histo for EtaPrims
    TH2F*                 fHistPtYOmega;              //! histo for Omegas
    Int_t                 fIsMC;                      // MC flag

  private:
    AliAnalysisTaskGammaCocktailMC(const AliAnalysisTaskGammaCocktailMC&); // Prevent copy-construction
    AliAnalysisTaskGammaCocktailMC &operator=(const AliAnalysisTaskGammaCocktailMC&); // Prevent assignment

    ClassDef(AliAnalysisTaskGammaCocktailMC, 1);
};

#endif
