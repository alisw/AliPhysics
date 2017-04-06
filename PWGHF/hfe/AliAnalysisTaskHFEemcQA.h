#ifndef AliAnalysisTaskHFEemcQA_cxx
#define AliAnalysisTaskHFEemcQA_cxx

//QA task for EMCAL electron analysis

class TH1F;
class THnSparse;
class AliESDEvent;
class AliAODEvent;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;
class AliAODMCHeader;
class AliAODMCParticle; // sample
class AliEMCALTriggerPatchInfo;
class AliMultSelection;
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHFEemcQA : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskHFEemcQA();
    AliAnalysisTaskHFEemcQA(const char *name);
    virtual ~AliAnalysisTaskHFEemcQA();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    
    Bool_t GetElecIDsparse() {return fFlagSparse;};
    void SetElecIDsparse(Bool_t flagelecIDsparse){fFlagSparse = flagelecIDsparse;};
    
    Bool_t GetTenderSwitch() {return fUseTender;};
    void SetTenderSwitch(Bool_t usetender){fUseTender = usetender;};
    
    Bool_t GetEMCalTriggerEG1() { return fEMCEG1; };
    Bool_t GetEMCalTriggerEG2() { return fEMCEG2; };
    void SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;};
    void SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;};
    Bool_t GetEMCalTriggerDG1() { return fDCalDG1; };
    Bool_t GetEMCalTriggerDG2() { return fDCalDG2; };
    void SetEMCalTriggerDG1(Bool_t flagTr1) { fDCalDG1=flagTr1; fDCalDG2=kFALSE;};
    void SetEMCalTriggerDG2(Bool_t flagTr2) { fDCalDG2=flagTr2; fDCalDG1=kFALSE;};
    
    void SetClusterTypeEMC(Bool_t flagClsEMC) {fFlagClsTypeEMC = flagClsEMC;};
    void SetClusterTypeDCAL(Bool_t flagClsDCAL) {fFlagClsTypeDCAL = flagClsDCAL;};
    
    void SetCentralityMim(Int_t centMim) {fcentMim = centMim;};
    void SetCentralityMax(Int_t centMax) {fcentMax = centMax;};
    void SetCentralityEstimator(const char *estimator) { fCentralityEstimator = estimator; }
    
    void CheckMCgen(AliAODMCHeader* mcHeader);
    void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC);
    void SetThresholdEG2(Int_t threshold) { fThresholdEG2=threshold; };
    void SetThresholdEG1(Int_t threshold) { fThresholdEG1=threshold; };
    void FindPatches(Bool_t &hasfiredEG1,Bool_t &hasfiredEG2,Double_t emceta, Double_t emcphi);
    void FindMother(AliAODMCParticle* part, Int_t &label, Int_t &pid);
    void GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
private:
    enum{
        kAODanalysis = BIT(20),
    };
    
    AliVEvent   *fVevent;  //!event object
    AliESDEvent *fESD;    //!ESD object
    AliAODEvent *fAOD;    //!AOD object
    AliAODMCHeader *fMCheader;
    AliPIDResponse *fpidResponse; //!pid response
    AliEMCALGeometry *fEMCALGeo;
    
    Bool_t      fFlagSparse;// switch to THnspare
    Bool_t       fUseTender;// switch to add tender
    
    Bool_t    fEMCEG1;//EMcal Threshold EG1
    Bool_t    fEMCEG2;//EMcal Threshold EG2
    Bool_t    fDCalDG1;//DCal Threshold DG1
    Bool_t    fDCalDG2;//DCal Threshold DG2
    
    TClonesArray  *fTracks_tender;//Tender tracks
    TClonesArray  *fCaloClusters_tender;//Tender cluster
    
    AliAODMCParticle  *fMCparticle;//! MC particle
    TClonesArray  *fMCarray;//! MC array
    
    AliMultSelection *fMultSelection;
    
    TClonesArray *fTriggersInfo;//TClonesArray to access container from EMCalTriggerMaker
    Int_t fThresholdEG2;// Threshold for EG2 trigger in ADC for trigger patches
    Int_t fThresholdEG1;// Threshold for EG1 trigger in ADC for trigger patches
    
    Bool_t fFlagClsTypeEMC;//switch to select EMC clusters
    Bool_t fFlagClsTypeDCAL;//switch to select DCAL clusters
    
    Int_t fcentMim; // mim. centrality
    Int_t fcentMax; // max. centrality
    
    TString fCentralityEstimator;         // Centrality Estimator
    
    TList       *fOutputList; //!Output list
    TH1F        *fNevents;//! no of events
    TH1F        *fCent;//! centrality
    TH2F        *fMult;//! track multiplicity vs centrality
    TH2F        *fEvPlaneV0;//! V0 event plane
    TH2F        *fEvPlaneV0A;//! V0A event plane
    TH2F        *fEvPlaneV0C;//! V0C event plane
    TH2F        *fEvPlaneTPC;//! TPC event plane
    TH1F        *fVtxZ;//!Vertex z
    TH1F        *fVtxX;//!Vertex x
    TH1F        *fVtxY;//!Vertex y
    TH2F        *fTrigMulti;//!trigger multiplicity
    TH1F        *fHistClustE;//!cluster energy
    TH1F        *fHistNonLinClustE;//!Nonlinear corrected cluster energy
    TH2F        *fHistClustEcent;//!cluster energy
    TH2F        *fEMCClsEtaPhi;//! EMC cluster eta and phi
    TH1F        *fHistClustEEG1;//! Cluster Energy, Trigger patch > ThresholdEG1
    TH2F        *fHistClustEEG1cent;//! Cluster Energy, Trigger patch > ThresholdEG1
    TH1F        *fHistClustEEG2;//! Cluster Energy, Trigger patch > ThresholdEG1
    TH2F        *fHistClustEEG2cent;//! Cluster Energy, Trigger patch > ThresholdEG1
    TH2F        *fEMCClsEtaPhiEG1;//! EMC cluster eta and phi, Trigger patch > ThresholdEG1
    TH2F        *fEMCClsEtaPhiEG2;//! EMC cluster eta and phi, Trigger patch > ThresholdEG2
    TH1F        *fHistoNCls;//! No of clusters per event
    TH1F        *fHistoNClsE1;//! No of clusters per event
    TH1F        *fHistoNClsE2;//! No of clusters per event
    TH1F        *fHistoNClsE3;//! No of clusters per event
    //TH1F        *fHistoNCells;//! No of cells per cluster
    TH2F        *fHistoNCells;//! No of cells per cluster
    TH2F        *fHistoEperCell;
    TH2F        *fHistoCalCell;//! No of cells per cluster
    TH2F        *fHistoTimeEMC;//! No of cells per cluster
    THnSparse   *fHistoTimeEMCcorr;//! No of cells per cluster
    TH1F        *fNegTrkIDPt;//!neg track ID
    TH1F        *fTrkPt;//!track pt
    TH1F        *fTrketa;//!track eta
    TH1F        *fTrkphi;//!track phi
    TH2F        *fdEdx;//!dedx vs pt
    TH2F        *fTPCNpts;//!TPC Npoints used for dedx
    TH2F        *fTPCnsig;//!TPC Nsigma
    TH2F        *fTPCnsigMcEle;//!TPC Nsigma
    TH2F        *fTPCnsigMcHad;//!TPC Nsigma
    TH2F        *fTPCnsig_Pi;//!TPC Nsigma wrt pion
    TH2F        *fTPCnsigEta0;//!TPC Nsigma
    TH2F        *fTPCnsigEta1;//!TPC Nsigma
    TH2F        *fTPCnsigEta2;//!TPC Nsigma
    TH1F        *fHistPtMatch;//!tracks matched to EMCAL
    TH2F        *fEMCTrkMatch;//!Distance of EMC cluster to closest track in phi and z
    TH1F        *fEMCTrkPt;//!tracks with EMCAL cluster
    TH1F        *fEMCTrketa;//!EMC trk eta
    TH1F        *fEMCTrkphi;//!EMC trk phi
    TH2F        *fEMCdEdx;//!EMC trk dedx
    TH2F        *fEMCTPCnsig;//! EMC trk nsig
    TH2F        *fEMCTPCNpts;//!EMC Npoints used for dedx
    TH1F        *fClsEAftMatch;//!EMC Cluster energy after track matching
    TH1F        *fNonLinClsEAftMatch;//!Nonlinear corrected EMC Cluster energy after track matching
    TH2F        *fClsEtaPhiAftMatch;//!EMC Cluster eta phi distribution after track matching
    TH2F        *fClsEtaPhiAftMatchEMCin;//!EMC Cluster eta phi distribution after track matching inside EMC phi acceptance
    TH2F        *fClsEtaPhiAftMatchEMCout;//!EMC Cluster eta phi distribution after track matching outside EMC phi acceptance
    TH2F        *fHistdEdxEop;//!E/p vs dedx
    TH2F        *fHistNsigEop;//!E/p vs dedx
    TH2F        *fHistNsigEop_Most;//!pt vs E/p
    TH2F        *fHistNsigEop_Semi;//!pt vs E/p
    TH2F        *fHistNsigEop_Peri;//!pt vs E/p
    TH2F        *fHistEop;//!pt vs E/p
    TH2F        *fHistMcEopEle;//!pt vs E/p
    TH2F        *fHistMcEopHad;//!pt vs E/p
    TH2F        *fM20;//!M20 vs pt
    TH2F        *fM02;//!M20 vs pt
    TH2F        *fM20EovP;//!M20 vs E/p
    TH2F        *fM02EovP;//!M20 vs E/p
    TH2F        *fEleCanTPCNpts;//!ele cand TPC Npoints used for dedx
    TH2F        *fEleCanTPCNCls;//!ele cand TPC N clusters
    TH2F        *fEleCanITSNCls;//!ele cand ITS N clusters
    TH1F        *fEleCanITShit;//!ele cand ITS hit map
    TH2F        *fEleCanSPD1;//!ele cand hit SPD layer 1
    TH2F        *fEleCanSPD2;//!ele cand hit SPD layer 2
    TH2F        *fEleCanSPDBoth;//!ele cand SPD both layer
    TH2F        *fEleCanSPDOr;//!ele cand SPD or
    TH2F        *fITShitPhi;//!ele cand SPD or
    TH1F        *fInvmassULS;//!Invmass of ULS
    TH1F        *fInvmassLS;//!Invmass of LS
    TH2F        *fInvmassULS_MCtrue;//!Invmass of ULS
    THnSparse   *fInvmassPi0Dalitz;//!Invmass of ULS
    TH2F        *fMCcheckMother;
    TH2F        *fMCneutral;
    
    THnSparse  *fSparseElectron;//!Electron info
    Double_t *fvalueElectron;//!Electron info
    
    AliAnalysisTaskHFEemcQA(const AliAnalysisTaskHFEemcQA&); // not implemented
    AliAnalysisTaskHFEemcQA& operator=(const AliAnalysisTaskHFEemcQA&); // not implemented
    
    ClassDef(AliAnalysisTaskHFEemcQA, 1); // example of analysis
};

#endif
