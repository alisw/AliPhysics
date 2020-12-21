#ifndef ALIANLYSISTASKOMEGAMCSTUDIES_cxx
#define ALIANLYSISTASKOMEGAMCSTUDIES_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

class AliAnalysisTaskOmegaMCStudies : public AliAnalysisTaskSE {
  public:
    /**
     * @enum AcceptanceType_t
     * @brief Enumeration for acceptance type
     */
    enum AcceptanceType_t {
      kPCMAcceptance = 1,       //!< PCM Acceptance
      kEMCALAcceptance = 2,     //!< EMCAL Acceptance
      kPHOSAcceptance = 3       //!< PHOS Acceptance
    };

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
      kPdgXi0 = 3322,          //!< kPdgXi0
      kPdgPhoton = 22,         //!< kPdgPhoton
      kPdgElectron = 11,       //!< kPdgElectron
      kPdgGluon = 21,          //!< kPdgGluon
    };

    // definition of hard initial processes
    enum ProcessCodes_t {
      kPromptPhotonCompton          = 201,           //!< kPromptPhotonCompton
      kPromptPhotonAnnihilation     = 202,           //!< kPromptPhotonAnnihilation
      kPromptPhotongg2ggamma        = 203,           //!< kPromptPhotongg2ggamma 
      kPromptPhotonffbar2gammagamma = 204,           //!< kPromptPhotonffbar2gammagamma 
      kPromptPhotongg2gammagamma    = 205              //!< kPromptPhotongg2gammagamma 
    };

    AliAnalysisTaskOmegaMCStudies();
    AliAnalysisTaskOmegaMCStudies(const char *name);
    virtual ~AliAnalysisTaskOmegaMCStudies();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void ProcessMCParticles();
    bool IsPiPlusPiMinusPiZeroDecay(TParticle* part) const;
    bool IsPiPlusPiMinusEtaDecay(TParticle* part) const;
    Int_t ReturnPi0FromOmega(TParticle* part);
    Int_t ReturnEtaFromEtaPrime(TParticle* part);
    bool IsInPCMAcceptance(TParticle* part) const;
    bool IsInPHOSAcceptance(TParticle* part) const;
    bool IsInEMCalAcceptance(TParticle* part) const;
    bool IsInFOCALAcceptance(TParticle* part) const;
    bool IsInLHCbAcceptance(TParticle* part) const;
    bool IsInMidAcceptance(TParticle* part) const;

    // additional functions
    void SetLogBinningXTH1(TH1* histoRebin);
    void SetLogBinningXTH2(TH2* histoRebin);
    void SetLogBinningYTH2(TH2* histoRebin);
    void SetMaxPt(Double_t pTmax){fMaxpT = pTmax;}
    void SetNEvents(Double_t nevents){fNTotEvents = nevents;}



  protected:
    TList*                fOutputContainer;           //! Output container
    // histograms events
    TH1F*                 fHistNEvents;               //! number of events histo
    TH1D*                 fHistXSection;              //! xSection
    TH1F*                 fHistPtHard;                //! ptHard
    // histograms mesons
    TH2F*                 fHistPtYOmega;              //! histo for omega pT distribution
    TH2F*                 fHistPtYOmegaPiPiPi;        //! distribution of omegas that decay to pi+pi-pi0
    TH2F*                 fHistPtYPi0;                //! histo for pi0   in general
    TH2F*                 fHistPtYEtaPrime;           //! histo for eta' pT distribution
    TH2F*                 fHistPtYEtaPrimeEtaPiPi;    //! distribution of eta' that decay to pi+pi-eta
    TH2F*                 fHistOmegaPtPi0Pt;          //! pT of omega vs pT pi0
    TH2F*                 fHistEtaPrimePtEtaPt;       //! pT of eta' vs pT eta
    
    Int_t                 fIsMC;                      // MC flag
    Double_t              fMaxpT;                     // pt hat
    Int_t                 fProcessCode;               // process of initial scattering
    Float_t               fWeight;                    // cross section/ ntrials
    Int_t                 fEventCounter;              // counter for events
    Int_t                 fNTotEvents;                // total number of events that will be generated ( needed for book keeping)
    Int_t                 fNRejectEvents;             // number of events rejected at the beginning to get a stable xsec

  private:
    AliAnalysisTaskOmegaMCStudies(const AliAnalysisTaskOmegaMCStudies&); // Prevent copy-construction
    AliAnalysisTaskOmegaMCStudies &operator=(const AliAnalysisTaskOmegaMCStudies&); // Prevent assignment

    ClassDef(AliAnalysisTaskOmegaMCStudies, 1);
};

#endif
