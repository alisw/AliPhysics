#ifndef ALIANLYSISTASKGAMMAMCSTUDIES_cxx
#define ALIANLYSISTASKGAMMAMCSTUDIES_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

class AliAnalysisTaskGammaMCStudies : public AliAnalysisTaskSE {
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

    AliAnalysisTaskGammaMCStudies();
    AliAnalysisTaskGammaMCStudies(const char *name);
    virtual ~AliAnalysisTaskGammaMCStudies();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void ProcessMCParticles();
    void ProcessPhoton(TParticle* part);
    TParticle* GetInitialPhoton(TParticle* part) const;
    Int_t IsDirectPhoton(TParticle* photon) const;
    bool IsInPCMAcceptance(TParticle* part) const;
    bool IsInPHOSAcceptance(TParticle* part) const;
    bool IsInEMCalAcceptance(TParticle* part) const;
    bool IsInFOCALAcceptance(TParticle* part) const;
    bool IsInLHCbAcceptance(TParticle* part) const;
    bool IsInMidAcceptance(TParticle* part) const;

    Float_t CalculateEt(TParticle* part, Float_t rmax);
    Float_t CalculateX1LO(Float_t y, Float_t q, Float_t cme);
    Float_t CalculateX2LO(Float_t y, Float_t q, Float_t cme);
    
    void CalculatePhotonCorrelations(TParticle* part);

    // additional functions
    void SetLogBinningXTH1(TH1* histoRebin);
    void SetLogBinningXTH2(TH2* histoRebin);
    void SetLogBinningYTH2(TH2* histoRebin);
    void SetMaxPt(Double_t pTmax){fMaxpT = pTmax;}
    void SetNEvents(Double_t nevents){fNTotEvents = nevents;}
    void SetNRejectEvents(Double_t nrejectevents){fNRejectEvents = nrejectevents;}



  protected:
    TList*                fOutputContainer;           //! Output container
    // histograms events
    TH1F*                 fHistNEvents;               //! number of events histo
    TH1D*                 fHistXSection;              //! xSection
    TH1F*                 fHistPtHard;                //! ptHard
    // histograms mesons
    TH2F*                 fHistPtYPi0;                //! histo for Pi0s
    TH2F*                 fHistIsoYPi0;                //! histo for Pi0s

    // histograms photons:
    TH2F*                 fHistPtYAllPhotons;              //! Et vs y of all photons found
    TH2F*                 fHistPtYDecayPhotons;            //! Et vs y of all decay photons found
    TH2F*                 fHistPtYDirectPhotons;           //! Et vs y of all direct photons found
    TH2F*                 fHistPtYComptonPhotons;          //! Et vs y of all q g->gamma g photons found
    TH2F*                 fHistPtYAnnihilationPhotons;     //! Et vs y of all q anti q->gamma g photons found
    TH2F*                 fHistPtYgg2ggammaPhotons;        //! Et vs y of all g g-> g gamme found
    TH2F*                 fHistPtYffbar2gammagammaPhotons; //! Et vs y of all q anti q->gamma gamma photons found
    TH2F*                 fHistPtYgg2gammagammaPhotons;    //! Et vs y of all g g->gamma gamma photons found
    TH2F*                 fHistPtYFragmentationPhotons;    //! Et vs y of all frag photons found

    TH2F*                 fHistIsoYAllPhotons;              //! Et vs y of all photons found
    TH2F*                 fHistIsoYDecayPhotons;            //! Et vs y of all decay photons found
    TH2F*                 fHistIsoYDirectPhotons;           //! Et vs y of all direct photons found
    TH2F*                 fHistIsoYComptonPhotons;          //! Et vs y of all q g->gamma g photons found
    TH2F*                 fHistIsoYAnnihilationPhotons;     //! Et vs y of all q anti q->gamma g photons found
    TH2F*                 fHistIsoYgg2ggammaPhotons;        //! Et vs y of all g g-> g gamme found
    TH2F*                 fHistIsoYffbar2gammagammaPhotons; //! Et vs y of all q anti q->gamma gamma photons found
    TH2F*                 fHistIsoYgg2gammagammaPhotons;    //! Et vs y of all g g->gamma gamma photons found
    TH2F*                 fHistIsoYFragmentationPhotons;    //! Et vs y of all frag photons found

    TH2F*                 fHistXGluonQComptonPhotons;      //! X of gluon vs Q quared for compton photons
    TH2F*                 fHistXGluonQComptonPhotonsFOCAL; //! X of gluon vs Q quared for compton photons
    TH2F*                 fHistXGluonQComptonPhotonsLHCb;  //! X of gluon vs Q quared for compton photons
    TH2F*                 fHistXGluonQComptonPhotonsMid;   //! X of gluon vs Q quared for compton photons

    TH2F*                 fHistXGluonPtComptonPhotons;      //! X of gluon vs pT of compton photons
    TH2F*                 fHistXGluonPtComptonPhotonsFOCAL; //! X of gluon vs pT of photons
    TH2F*                 fHistXGluonPtComptonPhotonsLHCb;  //! X of gluon vs pT of compton photons
    TH2F*                 fHistXGluonPtComptonPhotonsMid;   //! X of gluon vs pT of compton photons

    TH2F*                 fHistXGluonXQuark;            // x gluon vs x quark for conmpton processes
    TH2F*                 fHistXGluonXQuarkMid;            // x gluon vs x quark for conmpton processes
    TH2F*                 fHistXGluonXQuarkFocal;            // x gluon vs x quark for conmpton processes
    TH2F*                 fHistXGluonXQuarkLHCb;            // x gluon vs x quark for conmpton processes

    TH2F*                 fHistEtaPhotonEtaQuark;            // x gluon vs x quark for conmpton processes
    TH2F*                 fHistEtaPhotonEtaQuarkMid;            // x gluon vs x quark for conmpton processes
    TH2F*                 fHistEtaPhotonEtaQuarkFocal;            // x gluon vs x quark for conmpton processes
    TH2F*                 fHistEtaPhotonEtaQuarkLHCb;            // x gluon vs x quark for conmpton processes

    TH2F*                 fHistDeltaPhiPhotonQuark;            // x gluon vs x quark for conmpton processes
    
    TH1F*                 fHistDeltaPhiPhotonPionMid;            // x gluon vs x quark for conmpton processes
  
    TH2F*                 fHistX1vsX1Calc;            // pythia x1 vs calculated x1
    TH2F*                 fHistX2vsX2Calc;            // pythia x1 vs calculated x1
    
    TH1F*                 fHistXSecVsEvent;            // xsec vs n events

  

    TParticle*            fOutPart1;            // outgoing particle
    TParticle*            fOutPart2;            // outgoing particle


    Int_t                 fIsMC;                      // MC flag
    Double_t              fMaxpT;                     // pt hat
    Double_t              fMaxIso;                    // max iso for plot
    Float_t               fMinX;                    // max iso for plot
    Int_t                 fProcessCode;               // process of initial scattering
    Float_t               fX1;                        // x fractions of the two partons for which parton density values are defined
    Float_t               fX2;                        // x fractions of the two partons for which parton density values are defined
    Float_t               fQFrac;                     // the Q actorization scale at which the densities were evaluated
    Int_t                 fId1;                       // id of parton enetering scattering
    Int_t                 fId2;                       // id of parton enetering scattering
    Int_t                 fStatus1;                   // status of parton entering scattering
    Int_t                 fStatus2;                   // status of parton enetering scattering
    Float_t               fWeight;                    // cross section/ ntrials
    Int_t                 fEventCounter;              // counter for events
    Int_t                 fNTotEvents;                // total number of events that will be generated ( needed for book keeping)
    Int_t                 fNRejectEvents;             // number of events rejected at the beginning to get a stable xsec

  private:
    AliAnalysisTaskGammaMCStudies(const AliAnalysisTaskGammaMCStudies&); // Prevent copy-construction
    AliAnalysisTaskGammaMCStudies &operator=(const AliAnalysisTaskGammaMCStudies&); // Prevent assignment

    ClassDef(AliAnalysisTaskGammaMCStudies, 2);
};

#endif
