#ifndef ALIANLYSISTASKLMEEPUREMC_cxx
#define ALIANLYSISTASKLMEEPUREMC_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

class AliAnalysisTaskLMeePureMC : public AliAnalysisTaskSE {
  public:

     /**
     * @enum SupportedPdg_t
     * @brief Definition of constants for PDG codes used within the task
     */
    enum SupportedPdg_t {
      kPdgPi0 = 111,           //!< kPdgPi0
      kPdgRho0 = 113,          //!< kPdgRho0
      kPdgK0Long = 130,        //!< kPdgK0Long
      kPdgPiPlus = 211,        //!< kPdgPiPlus
      kPdgPiMinus = -211,      //!< kPdgPiMinus
      kPdgRhoPlus = 213,       //!< kPdgRhoPlus
      kPdgRhoMinus = -213,     //!< kPdgRhoMinus
      kPdgEta = 221,           //!< kPdgEta
      kPdgOmega = 223,         //!< kPdgOmega
      kPdgK0Short = 310,       //!< kPdgK0Short
      kPdgKStar = 313,         //!< kPdgKStar
      kPdgKPlus = 321,         //!< kPdgKPlus
      kPdgKMinus = -321,       //!< kPdgKMinus
      kPdgEtaPrime = 331,      //!< kPdgEtaPrime
      kPdgPhi = 333,           //!< kPdgPhi
      kPdgJPsi = 443,          //!< kPdgJPsi
      kPdgDeltaMinus = 1114,   //!< kPdgDeltaMinus
      kPdgDelta0 = 2114,       //!< kPdgDelta0s
      kPdgDeltaPlus = 2214,    //!< kPdgDeltaPlus
      kPdgDeltaPlusPlus = 2224,//!< kPdgDeltaPlusPlus
      kPdgSigmaMinus = 3112,   //!< kPdgSigmaMinus
      kPdgSigma0 = 3212,       //!< kPdgSigma0
      kPdgLambda = 3122,       //!< kPdgLambda
      kPdgSigmaPlus = 3222,    //!< kPdgSigmaPlus
      kPdgXiMinus = 3312,      //!< kPdgXiMinus
      kPdgXi0 = 3322           //!< kPdgXi0
    };

    AliAnalysisTaskLMeePureMC();
    AliAnalysisTaskLMeePureMC(const char *name);
    virtual ~AliAnalysisTaskLMeePureMC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    // MC functions
    void ProcessMCParticles();

  protected:
    TList*                fOutputContainer;           //! Output container
    // histograms events
    TH1F*                 fHistNEvents;               //! number of events histo
    TH1D*                 fHistXSection;              //! xSection
    TH1F*                 fHistPtHard;                //! ptHard
    // histograms mesons
    TH2F*                 fHistPtYPi0;                //! histo for Pi0s
    TH2F*                 fHistPtYPiPl;               //! histo for Pi+s
    TH2F*                 fHistPtYPiMi;               //! histo for Pi-s
    TH2F*                 fHistPtYEta;                //! histo for Etas
    TH2F*                 fHistPtYEtaPrime;           //! histo for EtaPrims
    TH2F*                 fHistPtYOmega;              //! histo for Omegas
    TH2F*                 fHistPtYRho0;               //! histo for rho0
    TH2F*                 fHistPtYRhoPl;              //! histo for rho+
    TH2F*                 fHistPtYRhoMi;              //! histo for rho-
    TH2F*                 fHistPtYPhi;                //! histo for phi
    TH2F*                 fHistPtYJPsi;               //! histo for J/psi
    TH2F*                 fHistPtYSigma0;             //! histo for Sigma0
    TH2F*                 fHistPtYK0s;                //! histo for K0s
    TH2F*                 fHistPtYK0l;                //! histo for K0l
    TH2F*                 fHistPtYK0star;             //! histo for K0*
    TH2F*                 fHistPtYDeltaPlPl;          //! histo for Delta++
    TH2F*                 fHistPtYDeltaPl;            //! histo for Delta+
    TH2F*                 fHistPtYDeltaMi;            //! histo for Delta-
    TH2F*                 fHistPtYDelta0;             //! histo for Delta0
    TH2F*                 fHistPtYLambda;             //! histo for Lambda
    TH2F*                 fHistPtYKPl;                //! histo for K+s
    TH2F*                 fHistPtYKMi;                //! histo for K-s

    TH2F*                 fHistPtYPi0FromEta;         //! histo for Pi0s from Etas
    TH2F*                 fHistPtYPi0FromLambda;      //! histo for Pi0s from Lambdas
    TH2F*                 fHistPtYPi0FromK;           //! histo for Pi0s from Ks
    TH2F*                 fHistPtYPiPlFromK;          //! histo for Pi+s form Ks
    TH2F*                 fHistPtYPiPlFromEta;        //! histo for Pi+s form Etas
    TH2F*                 fHistPtYPiMiFromK;          //! histo for Pi-s from Ks
    TH2F*                 fHistPtYPiMiFromEta;        //! histo for Pi-s from Etas
    Int_t                 fIsMC;                      // MC flag

  private:
    AliAnalysisTaskLMeePureMC(const AliAnalysisTaskLMeePureMC&); // Prevent copy-construction
    AliAnalysisTaskLMeePureMC &operator=(const AliAnalysisTaskLMeePureMC&); // Prevent assignment

    ClassDef(AliAnalysisTaskLMeePureMC, 1);
};

#endif
