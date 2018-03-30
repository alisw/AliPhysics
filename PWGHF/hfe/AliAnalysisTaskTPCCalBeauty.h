//
//  AliAnalysisTaskTPCCalBeauty.h
//  
//
//  Created by Erin Gauger
//
//

#ifndef AliAnalysisTaskTPCCalBeauty_cxx
#define AliAnalysisTaskTPCCalBeauty_cxx

#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
//#include "AliSelectNonHFE.h"
#include "AliAODMCParticle.h"

class THnSparse;
class AliMultSelection;
class AliAODMCParticle;
class AliAODMCHeader;

class AliAnalysisTaskTPCCalBeauty : public AliAnalysisTaskSE
{
    public:
        //two class constructors
                        AliAnalysisTaskTPCCalBeauty();
                        AliAnalysisTaskTPCCalBeauty(const char *name);
        //class destructor
        virtual         ~AliAnalysisTaskTPCCalBeauty();
        //called once at beginning of runtime
        virtual void    UserCreateOutputObjects();
        //called for each event
        virtual void    UserExec(Option_t *option);
        //called at end of analysis
        virtual void    Terminate(Option_t *option);
    
        Bool_t          GetTenderSwitch() {return fUseTender;};
        void            SetTenderSwitch(Bool_t usetender){fUseTender = usetender;};
    
        Bool_t          GetEMCalTriggerEG1() { return fEMCEG1; };
        void            SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; };
        Bool_t          GetEMCalTriggerDG1() { return fDCalDG1; };
        void            SetEMCalTriggerDG1(Bool_t flagTr1) { fDCalDG1=flagTr1; };
    
        void            SetClusterTypeEMC(Bool_t flagClsEMC) {fFlagClsTypeEMC = flagClsEMC;};
        void            SetClusterTypeDCAL(Bool_t flagClsDCAL) {fFlagClsTypeDCAL = flagClsDCAL;};
    
        void            SetCentralitySelection(Double_t centMin, Double_t centMax){fCentralityMin = centMin; fCentralityMax = centMax;};
        Double_t        CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass);
    
        Bool_t          GetNMCPartProduced();
        void            GetPi0EtaWeight(THnSparse *SparseWeight);

        void            GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
        void            FindMother(AliAODMCParticle* part, Int_t &fpidSort, Bool_t &kEmbEta, Bool_t &kEmbPi0, Bool_t &kHijing, Double_t &momPt);
        void            InvMassCheckData(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign);
        void            InvMassCheckMC(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign, Bool_t kHijing, Bool_t kEmbEta, Bool_t kEmbPi0, Bool_t &kFlagReco, Double_t fWeight, Int_t fpidSort);
    
    private:
        AliAODEvent         *fAOD;           //! input event
        AliAODMCHeader      *fMCHeader;      //!
        TClonesArray        *fMCarray;       //! MC array
        AliAODMCParticle    *fMCparticle;    //! MC particle
        AliPIDResponse      *fpidResponse;   //! pid response
        TList               *fOutputList;    //! output list
    
        AliMultSelection    *fMultSelection; //! need to get centrality
        Double_t            fCentrality;     //!
        Double_t            fCentralityMin;  // set event centrality min
        Double_t            fCentralityMax;  // set event centrality max
    
        Bool_t              fEMCEG1;         // EMCal Threshold EG1
        Bool_t              fDCalDG1;        // DCal Threshold DG1
        Bool_t              fFlagClsTypeEMC; // switch to select EMC clusters
        Bool_t              fFlagClsTypeDCAL;// switch to select DCAL clusters
        Bool_t              fUseTender;      // switch to add tender
        Bool_t              fFlagULS;        // flag ULS
        Bool_t              fFlagLS;         // flag LS
    
        TH1F                *fNevents;       //! no of events
        TH1F                *fVtX;           //! vertex x
        TH1F                *fVtY;           //! vertex y
        TH1F                *fVtZ;           //! vertex z
    
        TH1F                *fTrkPtB4TC;     //! track pT before track cuts
    
        TH1F                *fTrkPt;         //! track pT
        TH1F                *fTrkP;          //! track p
        TH1F                *fTrkClsPhi;     //! track and cluster delta phi
        TH1F                *fTrkClsEta;     //! track and cluster delta eta
        TH1F                *fClsPhi;         //! cluster phi
        TH1F                *fClsEta;         //! cluster eta
        TH1F                *fClsEamDCal;    //! cluster energy after matching to DCal
        TH1F                *fClsEamEMCal;   //! cluster energy after matching to EMCal
        TH1F                *fClsEAll;   //! cluster energy of all track-matched clusters
        TH1F                *fClsEamElecEMC;   //! cluster energy of e- after matching to EMCal
        TH1F                *fClsEamElecDC;   //! cluster energy of e- after matching to EMCal
        TH1F                *fTrkPhi;        //! track phi after track matching
        TH1F                *fTrkEta;        //! track eta after track matching
        TH1F                *fdEdx;          //! track dEdx
        TH1F                *fCentCheck;     //! event centrality
        TH1F                *fTrigCheck;     //! checking trigger used
        TH2F                *fEMCTrkMatch;   //! plots distance of cluster from closest track
    
        TH1F                *fInvmassLS;     //! Plots LS mass dist
        TH1F                *fInvmassULS;    //! Plots ULS mass dist
    
        TH1F                *fInvmassLSWeightEnhEta;   //! Plots LS mass dist
        TH1F                *fInvmassULSWeightEnhEta;  //!
        TH1F                *fInvmassLSWeightEnhPi0;   //!
        TH1F                *fInvmassULSWeightEnhPi0;  //!
        TH1F                *fInvmassLSHijingEta;      //!
        TH1F                *fInvmassULSHijingEta;     //!
        TH1F                *fInvmassLSHijingPi0;      //!
        TH1F                *fInvmassULSHijingPi0;     //!
        TH1F                *fInvmassLSHijingPhoton;   //!
        TH1F                *fInvmassULSHijingPhoton;  //!
        TH1F                *fInvmassLSEnhPhoton;         //!
        TH1F                *fInvmassULSEnhPhoton;        //!
    
        TH2F                *fULSdcaBelow;   //! ULS electron DCA vs. pT, m<0.1
        TH2F                *fLSdcaBelow;    //! LS electron DCA vs. pT, m<0.1
    
        TH1F                *fLSWeightEnhEta;     //! LS for Weighted enhanced eta
        TH1F                *fULSWeightEnhEta; //! ULS for Weighted enhanced eta
        TH1F                *fLSWeightEnhPi0;  //! LS for Weighted enhanced pi0
        TH1F                *fULSWeightEnhPi0; //! ULS for Weighted enhanced pi0
        TH1F                *fLSHijingEta;     //! LS for hijing eta
        TH1F                *fULSHijingEta;    //! ULS for hijing eta
        TH1F                *fLSHijingPi0;     //! LS for hijing pi0
        TH1F                *fULSHijingPi0;    //! ULS for hijing pi0
        TH1F                *fLSHijingPhoton;  //! LS for hijing photon
        TH1F                *fULSHijingPhoton; //! ULS for hijing photon
        TH1F                *fLSEnhPhoton;     //! LS for all photon e
        TH1F                *fULSEnhPhoton;    //! ULS for all photon e
    
        TH2F                *fPhotonicDCA;   //! Photonic DCA using MC PID
        TH2F                *fInclElecDCA;   //! Inclusive electron DCA vs. pT
        TH2F                *fInclElecEoP;   //! Inclusive electron EoP vs. pT
        TH2F                *fHadronEoP;     //! Hadron EoP vs. pT
        TH2F                *fHadronDCA;     //! Hadron DCA vs. pT
    
        TF1                 *fPi0Weight;    //! Function to weight enhanced pi0
        TF1                 *fEtaWeight;    //! Function to weight enhanced eta
        //TF1                 *fPi0EtaWeight; //! Function to weight enhanced eta+pi0
        TH1F                *fDWeight; //!
    
    
        Double_t            fWeight;        //!
    
        TH2F                *fPi0DCA;           //! Pi0 DCA vs. pT
        TH2F                *fEtaDCA;           //! Eta DCA vs. pT
    
        TH2F                *fEnhEtaDCA;           //! 
        TH1F                *fEnhEtaWeightedPt;    //!
        TH2F                *fEnhPi0DCA;           //!
        TH1F                *fEnhPi0WeightedPt;    //!
        TH2F                *fEtaHijingDCA;        //!
        TH1F                *fEtaHijingPt;         //!
        TH2F                *fPi0HijingDCA;        //!
        TH1F                *fPi0HijingPt;         //!
        TH2F                *fPhotonHijingDCA;     //!
        TH1F                *fPhotonHijingPt;      //!
        TH2F                *fEnhPhotonDCA;        //!
        TH1F                *fEnhPhotonWeightedPt; //!
    
        TH1F                *ComboNumWeight;       //!
        TH1F                *ComboNumNoWeight;     //!
        TH1F                *ComboDenomWeight;     //!
        TH1F                *ComboDenomNoWeight;   //!
    
        TH1F                *DMesonPDG; //! plots abs(pdg) of D mesons in the stack
        TH1F                *fD0MesonPt;  //!
        TH1F                *fD0MesonFromDStarPt; //!
        TH1F                *fDPlusMesonPt; //!
        TH1F                *fDsMesonPt;//!
        TH1F                *fDStarMesonPt;//!
        TH1F                *fAllDMesonPt; //!
    
        TH1F                *fLambdaCPt; //!
        TH1F                *fEtaCPt; //!
        TH1F                *fCBaryonPt; //!
    
        TH2F    *fPromptD0DCAWeight; //!
        TH2F    *fD0FromDStarDCAWeight; //!
    TH2F    *fPromptD0DCANoWeight; //!
    TH2F    *fD0FromDStarDCANoWeight; //!

    
        Int_t               fNtotMCpart;     //! N of total MC particles produced by generator
        Int_t               fNpureMC;        //! N of particles from main generator (Hijing/Pythia)
        Int_t               fNembMCpi0;      //! N > fNembMCpi0 = particles from pi0 generator
        Int_t               fNembMCeta;      //! N > fNembMCeta = particles from eta generator
    
        THnSparse           *fSprsPi0EtaWeightCal;  //! Sparse for pi0,eta weight calc
        THnSparse           *fSprsTemplatesNoWeight;  //! Sparse for templates
        THnSparse           *fSprsTemplatesWeight;  //! Sparse for templates
    
        TH2F                *fDTemplateWeight; //!
        TH2F                *fDTemplateNoWeight; //!
    
        //THnSparse           *fElectronSprs;  //! Sparse with electron cut parameters
        //Double_t            *fvalueElectron; //! Electron info
    
    AliAnalysisTaskTPCCalBeauty(const AliAnalysisTaskTPCCalBeauty&); // not implemented???
    AliAnalysisTaskTPCCalBeauty& operator=(const AliAnalysisTaskTPCCalBeauty&); // not implemented???
    ClassDef(AliAnalysisTaskTPCCalBeauty, 1);
};

#endif
