#ifndef ALIANLYSISTASKGAMMAJETSPUREMC_cxx
#define ALIANLYSISTASKGAMMAJETSPUREMC_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "TRandom3.h"
#include "TH3.h"
#include "TF1.h"
#include "TObjString.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequenceArea.hh>


class AliAnalysisTaskGammaJetsPureMC : public AliAnalysisTaskSE {
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
      kPdgXi0 = 3322,          //!< kPdgXi0
      kPdgd = 1,               //!< kPdgd
      kPdgu = 2,               //!< kPdgu
      kPdgs = 3,               //!< kPdgs
      kPdgc = 4,               //!< kPdgc
      kPdgb = 5,               //!< kPdgb
      kPdgt = 6,               //!< kPdgt
      kPdgg = 21,              //!< kPdgg
      kGamma = 22,              //!< kGamma
      kNeutron = 2112,          //!< neutron
    };

    AliAnalysisTaskGammaJetsPureMC();
    AliAnalysisTaskGammaJetsPureMC(const char *name);
    virtual ~AliAnalysisTaskGammaJetsPureMC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    // MC functions
    void ProcessJets();

    // additional functions
    void SetLogBinningXTH1(TH1* histoRebin);
    void SetLogBinningXTH2(TH2* histoRebin);

    // setters
    void SetEfficiency(TString strNeutral = "0.7", TString strCharged = "0.8");
    void SetParticlesNonMeas(TString strPart);
    bool IsNonMeasureable(int pdgCode, int charge) const;
    void SetJetEnergyShift(double shift) { fDoJetEnergyShift = true; fJetEnergyShift = shift; }
    

  protected:
    bool AcceptParticle(AliVParticle* particle);
    void FillResponseMatrixAndEffi(std::vector<fastjet::PseudoJet> vecTrueJet, std::vector<fastjet::PseudoJet> vecRecJet, TH2F* hResp, TH1D* hUnMatched, TH1D* hMultiMatched, TH1D* hRecUnMatched, double energyShiftRec = 1. );
    int GetParticleIndex(const int pdgcode, const int motherpdg) const;

    TList*                fOutputContainer;           //! Output container
    // histograms events
    TH1F*                 fHistNEvents;               //! number of events histo
    TH1D*                 fHistXSection;              //! xSection
    TH1F*                 fHistPtHard;                //! ptHard
    
    // jet studies
    TH2F*                 fHistJetPtY_Std;  //! histo for jet pt
    TH1D*                 fHistJetEta_Std;  //! histo for jet eta
    TH1D*                 fHistJetPhi_Std;  //! histo for jet phi

    TH2F*                 fHistJetPtY_StdNN;  //! histo for jet pt
    TH1D*                 fHistJetEta_StdNN;  //! histo for jet eta
    TH1D*                 fHistJetPhi_StdNN;  //! histo for jet phi

    TH2F*                 fHistJetPtY_Det;  //! histo for jet pt
    TH1D*                 fHistJetEta_Det;  //! histo for jet eta
    TH1D*                 fHistJetPhi_Det;  //! histo for jet phi

    TH2F*                 fHistJetPtY_DetNN;  //! histo for jet pt
    TH1D*                 fHistJetEta_DetNN;  //! histo for jet eta
    TH1D*                 fHistJetPhi_DetNN;  //! histo for jet phi

    TH2F*                 fHistResponse_Std_Det;  //!
    TH1D*                 fHistUnMatched_Std_Det; //!
    TH1D*                 fHistRecUnMatched_Std_Det; //!
    TH1D*                 fHistMultiMatched_Std_Det; //!

    TH2F*                 fHistResponse_Std_DetNN;  //!
    TH1D*                 fHistUnMatched_Std_DetNN; //!
    TH1D*                 fHistRecUnMatched_Std_DetNN; //!
    TH1D*                 fHistMultiMatched_Std_DetNN; //!

    TH3F*                 fHistJetPtPartPtVsPart; //!
    TH3F*                 fHistJetPtPartFragVsPart; //!

    TH3F*                 fHistJetPtPartPtVsPartLead; //!
    TH3F*                 fHistJetPtPartFragVsPartLead; //!

    TH1D*                 fHistPrimaryParticles; //!

    TH3F*                 fHistEnergyFracParticleCat; //!

    // jet finding
    double fJetRadius;                                      // jet radius parameter
    double fJetMinE;                                        // minimum jet energy
    double fJetAccEta;                                      // eta acceptance of jet
    double fJetParticleAcc;                                 // acceptance of particles that contribute to jet for pt spectra
    double fJetParticleAccFF;                               // acceptance of particles that contribute to jet for FF
    fastjet::JetAlgorithm   	    fJetAlgorithm;            // jet algorithm
    fastjet::Strategy 		        fJetStrategy;             // jet strategy parameter
    fastjet::AreaType 		        fJetAreaType;             // jet rea type parameter
    fastjet::RecombinationScheme  fJetRecombScheme;         // jet recomb. scheme parameter
    double fJetGhostArea;                                   // jet ghost area
    double fGhostEtaMax;                                    // maximum eta of ghost particles
    int fActiveAreaRepeats;                                 // jet active area
    fastjet::AreaType 		  fAreaType;                      // jet area type
    TF1* fEffiNeutral;                                      // efficiency function for neutral particles
    TF1* fEffiCharged;                                      // efficiency function for charged particles
    std::vector<int> fVecNonMeasureable;                    // vector holding the non-measureable particles

    bool fDoJetEnergyShift;                                 // Boolean for shift jet energy for response studies
    double fJetEnergyShift;                                 // Shift jet energy for response studies

    TRandom3 fRand;
     

  private:
    AliAnalysisTaskGammaJetsPureMC(const AliAnalysisTaskGammaJetsPureMC&); // Prevent copy-construction
    AliAnalysisTaskGammaJetsPureMC &operator=(const AliAnalysisTaskGammaJetsPureMC&); // Prevent assignment

    ClassDef(AliAnalysisTaskGammaJetsPureMC, 3);
};

#endif
