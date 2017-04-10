#ifndef ALIANLYSISTASKGAMMAPUREMC_cxx
#define ALIANLYSISTASKGAMMAPUREMC_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

class AliAnalysisTaskGammaPureMC : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskGammaPureMC();
    AliAnalysisTaskGammaPureMC(const char *name);
    virtual ~AliAnalysisTaskGammaPureMC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void ProcessMCParticles();
    void IsInPCMAcceptance(TParticle* part, Int_t& acceptance);
    void IsInPHOSAcceptance(TParticle* part, Int_t& acceptance);
    void IsInEMCalAcceptance(TParticle* part, Int_t& acceptance);
    
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
    TH1D*                 fHistXSection;              // xSection
    TH1F*                 fHistPtHard;                // ptHard 
    // histograms mesons
    TH2F*                 fHistPtYPi0;                //! histo for Pi0s
    TH2F*                 fHistPtYPiPl;               //! histo for Pi+s 
    TH2F*                 fHistPtYPiMi;               //! histo for Pi-s
    TH2F*                 fHistPtYEta;                //! histo for Etas
    TH2F*                 fHistPtYEtaPrim;            //! histo for EtaPrims
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

    TH2F*                 fHistPtYPi0FromEta;         //! histo for Pi0s from Etas
    TH2F*                 fHistPtYPi0FromLambda;      //! histo for Pi0s from Lambdas
    TH2F*                 fHistPtYPi0FromK;           //! histo for Pi0s from Ks
    TH2F*                 fHistPtYPiPlFromK;          //! histo for Pi+s form Ks
    TH2F*                 fHistPtYPiMiFromK;          //! histo for Pi-s from Ks
    TH2F*                 fHistPtYPi0GG;              //! histo for Pi0s gamma gamma channel 
    TH2F*                 fHistPtYPi0GGPCMAcc;        //! histo for Pi0s gamma gamma channel, PCM acceptance 
    TH2F*                 fHistPtYPi0GGEMCAcc;        //! histo for Pi0s gamma gamma channel, EMCal acceptance 
    TH2F*                 fHistPtYPi0GGPHOAcc;        //! histo for Pi0s gamma gamma channel, PHOS acceptance 
    TH2F*                 fHistPtYPi0GGPCMEMCAcc;     //! histo for Pi0s gamma gamma channel, PCM-EMCal acceptance 
    TH2F*                 fHistPtYPi0GGPCMPHOAcc;     //! histo for Pi0s gamma gamma channel, PCM-PHOS acceptance
    TH2F*                 fHistPtAlphaPi0GGPCMAcc;    //! histo for Pi0s gamma gamma channel, PCM acceptance 
    TH2F*                 fHistPtAlphaPi0GGEMCAcc;    //! histo for Pi0s gamma gamma channel, EMCal acceptance 
    TH2F*                 fHistPtAlphaPi0GGPHOAcc;    //! histo for Pi0s gamma gamma channel, PHOS acceptance 
    TH2F*                 fHistPtAlphaPi0GGPCMEMCAcc; //! histo for Pi0s gamma gamma channel, PCM-EMCal acceptance 
    TH2F*                 fHistPtAlphaPi0GGPCMPHOAcc; //! histo for Pi0s gamma gamma channel, PCM-PHOS acceptance

    TH2F*                 fHistPtYEtaGG;              //! histo for Etas gamma gamma channel 
    TH2F*                 fHistPtYEtaGGPCMAcc;        //! histo for Etas gamma gamma channel, PCM acceptance 
    TH2F*                 fHistPtYEtaGGEMCAcc;        //! histo for Etas gamma gamma channel, EMCal acceptance 
    TH2F*                 fHistPtYEtaGGPHOAcc;        //! histo for Etas gamma gamma channel, PHOS acceptance 
    TH2F*                 fHistPtYEtaGGPCMEMCAcc;     //! histo for Etas gamma gamma channel, PCM-EMCal acceptance 
    TH2F*                 fHistPtYEtaGGPCMPHOAcc;     //! histo for Etas gamma gamma channel, PCM-PHOS acceptance
    TH2F*                 fHistPtAlphaEtaGGPCMAcc;    //! histo for Etas gamma gamma channel, PCM acceptance 
    TH2F*                 fHistPtAlphaEtaGGEMCAcc;    //! histo for Etas gamma gamma channel, EMCal acceptance 
    TH2F*                 fHistPtAlphaEtaGGPHOAcc;    //! histo for Etas gamma gamma channel, PHOS acceptance 
    TH2F*                 fHistPtAlphaEtaGGPCMEMCAcc; //! histo for Etas gamma gamma channel, PCM-EMCal acceptance 
    TH2F*                 fHistPtAlphaEtaGGPCMPHOAcc; //! histo for Etas gamma gamma channel, PCM-PHOS acceptance

    TH2F*                 fHistPtYEtaPrimGG;          //! histo for Etas gamma gamma channel 
    TH2F*                 fHistPtYEtaPrimGGPCMAcc;    //! histo for EtaPrims gamma gamma channel, PCM acceptance 
    TH2F*                 fHistPtYEtaPrimGGEMCAcc;    //! histo for EtaPrims gamma gamma channel, EMCal acceptance 
    TH2F*                 fHistPtYEtaPrimGGPHOAcc;    //! histo for EtaPrims gamma gamma channel, PHOS acceptance 
    TH2F*                 fHistPtYEtaPrimGGPCMEMCAcc; //! histo for EtaPrims gamma gamma channel, PCM-EMCal acceptance 
    TH2F*                 fHistPtYEtaPrimGGPCMPHOAcc; //! histo for Pi0s gamma gamma channel, PCM-PHOS acceptance

    Int_t                 fIsMC;                      // MC flag

    
  private:
    AliAnalysisTaskGammaPureMC(const AliAnalysisTaskGammaPureMC&); // Prevent copy-construction
    AliAnalysisTaskGammaPureMC &operator=(const AliAnalysisTaskGammaPureMC&); // Prevent assignment

    ClassDef(AliAnalysisTaskGammaPureMC, 1);
};

#endif
