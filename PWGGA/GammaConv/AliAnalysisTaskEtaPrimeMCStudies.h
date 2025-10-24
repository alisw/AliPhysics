#ifndef ALIANALYSISTASKETAPRIMEMCSTUDIES_cxx
#define ALIANALYSISTASKETAPRIMEMCSTUDIES_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

// heavy meson decaying into 3 particles: eta' -> pi+ pi- eta

class AliAnalysisTaskEtaPrimeMCStudies : public AliAnalysisTaskSE{
    public:
        AliAnalysisTaskEtaPrimeMCStudies();
        AliAnalysisTaskEtaPrimeMCStudies(const char* name);
        virtual ~AliAnalysisTaskEtaPrimeMCStudies();

        virtual void UserCreateOutputObjects();
        virtual void UserExec(Option_t *);
        virtual void Terminate(const Option_t *);

        // MC functions
        void ProcessMCParticles();
        bool IsPhotonInPCMAcceptance(AliMCParticle* part) const;
        bool IsPhotonInEMCalAcceptance(AliMCParticle* part) const;
        bool IsEtaSelected(AliMCParticle* part) const;
        bool IsPionSelected(AliMCParticle* part) const;
        // bool IsInV0Acceptance(AliVParticle* part) const;



    protected:
        TList*              fOutputContainer;                               //! Output container
        AliVEvent*          fInputEvent;                                    //!<! current event
        AliMCEvent*         fMCEvent;                                       //!<! current MC event

        // histogram events
        TH1F*               fHistNEvents;                                   //! number of events histo
        TH1D*               fHistXSection;                                  //! xSection
        TH1F*               fHistPtHard;                                    //! pT-hard

        // histograms heavy meson and its decay products
        TH1F*               fHistEtaPrime_Pt;                               //! histo for eta prime decaying into pi+ pi- eta
        TH1F*               fHistEtaPrimeAll_Pt;                            //! histo for all eta prime
        TH1F*               fHistEtaFromEtaPrime_Pt;                        //! histo for eta from eta prime decay
        TH1F*               fHistPiPl_Pt;                                   //! histo for pi+
        TH1F*               fHistPiMi_Pt;                                   //! histo for pi-
        TH1F*               fHistGamma1_Pt;                                 //! histo for gamma1 from eta decay
        TH1F*               fHistGamma2_Pt;                                 //! histo for gamma2 from eta decay
        TH1F*               fHistPi0_Pt;                                    //! histo for pi0 for validation
        TH1F*               fHistEta_Pt;                                    //! histo for eta for validation

        // histograms of eta prime within acceptance ?
        TH1F*               fHistEtaPrime_InAcc_EMC;                        //! histo for eta' within EMC acceptance
        TH1F*               fHistEtaPrime_InAcc_PCMECM;                     //! histo for eta' within PCM-EMC acceptance

        // tree for EtaPrime and its decay products
        TTree*              fTreeEtaPrimeKinematics;                             //!<! tree for EtaPrime and its decay products

        // gamma branches
        UShort_t            fBuffer_Gamma1_pt;
        Short_t             fBuffer_Gamma1_eta;
        UShort_t            fBuffer_Gamma1_phi;
        Short_t             fBuffer_Gamma1_y;
        UShort_t            fBuffer_Gamma1_E;

        UShort_t            fBuffer_Gamma2_pt;
        Short_t             fBuffer_Gamma2_eta;
        UShort_t            fBuffer_Gamma2_phi;
        Short_t             fBuffer_Gamma2_y;
        UShort_t            fBuffer_Gamma2_E;

        // charged pions
        UShort_t            fBuffer_PiPl_pt;
        UShort_t            fBuffer_PiPl_E;
        Short_t             fBuffer_PiPl_eta;
        UShort_t            fBuffer_PiPl_phi;
        Short_t             fBuffer_PiPl_y;

        UShort_t            fBuffer_PiMi_pt;
        UShort_t            fBuffer_PiMi_E;
        Short_t             fBuffer_PiMi_eta;
        UShort_t            fBuffer_PiMi_phi;
        Short_t             fBuffer_PiMi_y;

        // Eta
        UShort_t            fBuffer_Eta_pt;
        UShort_t            fBuffer_Eta_E;
        Short_t             fBuffer_Eta_eta;
        UShort_t            fBuffer_Eta_phi;
        Short_t             fBuffer_Eta_y;

        // EtaPrime
        UShort_t            fBuffer_EtaPrime_pt;
        UShort_t            fBuffer_EtaPrime_E;
        Short_t             fBuffer_EtaPrime_eta;
        UShort_t            fBuffer_EtaPrime_phi;
        Short_t             fBuffer_EtaPrime_y;


    private:
        AliAnalysisTaskEtaPrimeMCStudies(const AliAnalysisTaskEtaPrimeMCStudies&);                // prevent copy-construction
        AliAnalysisTaskEtaPrimeMCStudies &operator=(const AliAnalysisTaskEtaPrimeMCStudies&);     // prevent assignment

        ClassDef(AliAnalysisTaskEtaPrimeMCStudies, 1);
};


#endif
