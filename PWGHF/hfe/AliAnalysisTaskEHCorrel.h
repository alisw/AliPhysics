#ifndef AliAnalysisTaskEHCorrel_h
#define AliAnalysisTaskEHCorrel_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//Task for Heavy Flavour Electron-Hadron DeltaPhi Correlation in Run 2//
//                                                                    //
//  Author: Deepa Thomas (University of Texas at Austin)              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class THnSparse;
class TH2F;
class TLorentzVector;

class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODTrack;
class AliCFManager;
class AliEventPoolManager;
class AliMultSelection;

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliSelectNonHFE.h"

class AliAnalysisTaskEHCorrel : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskEHCorrel();
    AliAnalysisTaskEHCorrel(const char *name);
    virtual ~AliAnalysisTaskEHCorrel();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    Bool_t  PassEventSelect(AliVEvent *fVevent);
    void SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;};
    void SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;};
    void    CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass);
    Bool_t  PassTrackCuts(AliAODTrack *atrack);
    void    GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
    Bool_t  PassEIDCuts(AliVTrack *track, AliVCluster *clust);
    Bool_t  PassHadronCuts(AliAODTrack *HadTrack);
    void    HadronInfo(Int_t itrack);
    void    ElectronHadCorrel(Int_t itrack, AliVTrack *track, THnSparse *SparseEHCorrl);
    void    EMCalClusterInfo();
    void    SelectNonHFElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec);
    
    
private:
    AliVEvent 		    *fVevent;		//!V event object
    AliAODEvent 		*fAOD;			//!AOD object
    AliPIDResponse      *fpidResponse; //!pid response
    AliMultSelection    *fMultSelection;//!

    
    Double_t            fCentrality;//!
    Double_t            fCentralityMin;//
    Double_t            fCentralityMax;//
    Bool_t              fEMCEG1;//
    Bool_t              fEMCEG2;//
    Bool_t              fFlagClsTypeEMC;//switch to select EMC clusters
    Bool_t              fFlagClsTypeDCAL;//switch to select DCAL clusters
    Double_t            fTPCnSigma;//!
    Double_t            fTPCnSigmaMin;//
    Double_t            fTPCnSigmaMax;//
    Double_t            fM02Min;//
    Double_t            fM02Max;//
    Double_t            fM20Min;//
    Double_t            fM20Max;//
    Double_t            fEovPMin;//
    Double_t            fEovPMax;//
    Int_t               fTPCNClsHad;// Had track TPC NClusters
    Double_t            fInvmassCut;//

    TList       	   	*fOutputList;		//!output list
    TH1F                *fNevents;		//!no of events
    TH1F                *fVtxZ;//!
    TH1F                *fVtxX;//!
    TH1F                *fVtxY;//!
    TH1F                *fCentralityNoPass;//!
    TH1F                *fCentralityPass;//!
    TH1F                *fHistClustE;//!
    TH2F                *fEMCClsEtaPhi;//!
    
    TH1F                *fNegTrkIDPt;//!
    TH1F                *fTrkPt;//!
    TH1F                *fTrketa;//!
    TH1F                *fTrkphi;//!
    TH2F                *fdEdx;//!
    TH2F                *fTPCnsig;//!
    TH1F                *fHistPtMatch;//!
    TH2F                *fEMCTrkMatch;//!
    TH1F                *fEMCTrkPt;//!
    TH1F                *fEMCTrketa;//!
    TH1F                *fEMCTrkphi;//!
    TH2F                *fEMCTPCnsig;//!
    TH1F                *fClsEAftMatch;//!
    TH2F                *fClsEtaPhiAftMatch;//!
    TH2F                *fHistNsigEop;//!
    TH2F                *fM20EovP;//!
    TH2F                *fM02EovP;//!
    TH2F                *fHistEop;//!
    TH2F                *fM20;//!
    TH2F                *fM02;//!
    TH2F                *fHistEop_AftEID;//!
    TH1F                *InclsElecPt;//!
    
    TH2F                *fHadronPhiPt;//!
    TH1F                *fHadronPhi;//!
    TH1F                *fHadronPhiTPChalf;//!
    TH1F                *fHadronPt;//!
    
    TH1F                *fInvmassLS;//!
    TH1F                *fInvmassULS;//!
    
    THnSparse           *fSprsInclusiveEHCorrl;//!
    
    
    AliAnalysisTaskEHCorrel(const AliAnalysisTaskEHCorrel&); // not implemented
    AliAnalysisTaskEHCorrel& operator=(const AliAnalysisTaskEHCorrel&); // not implemented
    
    ClassDef(AliAnalysisTaskEHCorrel, 2); //!example of analysis
};
#endif