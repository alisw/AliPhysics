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
//#include "AliCentrality.h"
//#include "AliSelectNonHFE.h"

class THnSparse;
class AliMultSelection;
class AliAODMCParticle;

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
    
        void            GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
        void            FindMother(AliAODMCParticle* part, int &ilabelM, int &pidM, Int_t &fpidSort);
        void            FindMother2(AliAODMCParticle* part, int &ilabelM, int &pidM, Int_t &fpidSort);
        void            InvMassCheck(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign);
    
    private:
        AliAODEvent         *fAOD;           //! input event
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
        TH1F                *fTrkPhi;        //! track phi after track matching
        TH1F                *fTrkEta;        //! track eta after track matching
        TH1F                *fdEdx;          //! track dEdx
        TH1F                *fCentCheck;     //! event centrality
        TH1F                *fTrigCheck;     //! checking trigger used
        TH2F                *fEMCTrkMatch;   //! plots distance of cluster from closest track
        TH1F                *fInvmassLS;     //! Plots LS mass dist
        TH1F                *fInvmassULS;    //! Plots ULS mass dist
        TH1F                *fPhotonicElecYield; //! Photonic electron yield vs. pT
        TH2F                *fULSdcaBelow;   //! ULS electron DCA vs. pT, m<0.1
        TH2F                *fLSdcaBelow;    //! LS electron DCA vs. pT, m<0.1
        THnSparse           *fElectronSprs;  //! Sparse with electron cut parameters
        //Double_t            *fvalueElectron; //! Electron info
    
    AliAnalysisTaskTPCCalBeauty(const AliAnalysisTaskTPCCalBeauty&); // not implemented???
    AliAnalysisTaskTPCCalBeauty& operator=(const AliAnalysisTaskTPCCalBeauty&); // not implemented???
    ClassDef(AliAnalysisTaskTPCCalBeauty, 1);
};

#endif
