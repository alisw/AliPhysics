#ifndef ALIANLYSISTASKGAMMAPYTHIA_cxx
#define ALIANLYSISTASKGAMMAPYTHIA_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

class AliAnalysisTaskGammaPythia : public AliAnalysisTaskSE {
  public:
    /**
     * @enum SupportedPdg_t
     * @brief Definition of constants for PDG codes used within the task
     */
    enum SupportedPdg_t {
      kPdgGamma = 22,          //!< kPdgGamma
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

    AliAnalysisTaskGammaPythia();
    AliAnalysisTaskGammaPythia(const char *name);
    virtual ~AliAnalysisTaskGammaPythia();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void ProcessMultiplicity();
    bool IsInV0Acceptance(AliVParticle* part) const;

    // additional functions
    void SetMaxY(Double_t Ymax){fMaxY = Ymax;}
    void SetMinNch(Double_t Nch_min){fNch_min = Nch_min;}
    void SetMaxNch(Double_t Nch_max){fNch_max = Nch_max;}
    void SetMaxPt(Double_t pTmax){fMaxpT = pTmax;}
    void SetDoMultStudies(Int_t tmp){fDoMultStudies = tmp;}

  protected:
    TList*                fOutputContainer;           //! Output container
    // histograms events
    TH1F*                 fHistNEvents;               //! number of events histo
    TH1F*                 fHistNEventsHM;             //! number of events histo
    TH1D*                 fHistXSection;              //! xSection
    TH1F*                 fHistPtHard;                //! ptHard

    TH1I*                 fHistParticlePDG;           //! histo for Gammas
    TH1I*                 fHistMotherParticlePDG;     //! histo for Gammas
    TH1I*                 fHistGrandMotherParticlePDG;//! histo for Gammas
    TH2I*                 fHistParticlevsMother;      //! histo for Gammas

  // histograms gammas
    TH1F*                 fHistPtGamma;                //! histo for Gammas
    TH2F*                 fHistPtYGamma;               //! histo for Gammas
    TH2F*                 fHistPtMultGamma;            //! histo for Gammas

    TH1D*                 fHistMult;                   //! histo for Pi0 pt vs multiplicity
    TH1D*                 fHistV0Mult;                 //! histo for Pi0 pt vs V0 multiplicity
    TH1D*                 fHistV0MultHM;                 //! histo for Pi0 pt vs V0 multiplicity

    Int_t                 fIsMC;            // MC flag
    Double_t              fMaxY;            // Max y flag
    Int_t                 fNch_min;          // Min Nch flag
    Int_t                 fNch_max;          // Max Nch flag
    Double_t              fMaxpT;           // Max pT flag
    Int_t                 fDoMultStudies;   // enable multiplicity dependent studies (0 -> No mult studies, 1 -> Mult estimation with V0, 2 -> Mult estimation with V0 and INEL>0 criterium for multiplicity)
    Int_t                 fNTracks;  // number of tracks
    Int_t                 fNTracksInV0Acc;  // number of tracks in V0A+C acceptance for multiplicity studies
    Int_t                 fNTracksInV0AccHM;  // number of tracks in V0A+C acceptance for multiplicity studies
    Bool_t                fIsEvtINELgtZERO; // flag if event is INEL>0

  private:
    AliAnalysisTaskGammaPythia(const AliAnalysisTaskGammaPythia&); // Prevent copy-construction
    AliAnalysisTaskGammaPythia &operator=(const AliAnalysisTaskGammaPythia&); // Prevent assignment

    ClassDef(AliAnalysisTaskGammaPythia, 3);
};

#endif
