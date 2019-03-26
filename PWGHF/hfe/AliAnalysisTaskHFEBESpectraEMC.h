#ifndef AliAnalysisTaskHFEBESpectraEMC_cxx
#define AliAnalysisTaskHFEBESpectraEMC_cxx

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

class AliAnalysisTaskHFEBESpectraEMC : public AliAnalysisTaskSE {
public:
    enum HijingOrNot {kHijing,kElse};
    enum pi0etaType {kNoMother, kNoFeedDown, kNotIsPrimary, kLightMesons, kBeauty, kCharm};
    
    AliAnalysisTaskHFEBESpectraEMC();
    AliAnalysisTaskHFEBESpectraEMC(const char *name);
    virtual ~AliAnalysisTaskHFEBESpectraEMC();
    
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
    
    void SetClusterTypeEMC(Bool_t flagClsEMC) {fFlagClsTypeEMC = flagClsEMC;};
    void SetClusterTypeDCAL(Bool_t flagClsDCAL) {fFlagClsTypeDCAL = flagClsDCAL;};
    
    void SetCentralityMim(Int_t centMim) {fcentMim = centMim;};
    void SetCentralityMax(Int_t centMax) {fcentMax = centMax;};
    void SetCentralityEstimator(const char *estimator) { fCentralityEstimator = estimator; }
    
    void GetEMCalClusterInfo();
    void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC);
    
    void FindMother(AliAODMCParticle* part, Int_t &label, Int_t &pid);
    void GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
    
    Bool_t  PassEIDCuts(AliVTrack *track, AliVCluster *clust, Bool_t &Hadtrack);
    
    Bool_t  GetNMCPartProduced();
    void    GetPi0EtaWeight(THnSparse *SparseWeight);
    Int_t   GetPi0EtaType(AliAODMCParticle *part);
    Bool_t  GetNonHFEEffiDenom(AliVTrack *track);
    Bool_t  IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromHijing, Int_t &type, Int_t &iMom, Int_t &MomPDG, Double_t &MomPt);
    Bool_t  GetNonHFEEffiRecoTag(AliVTrack *track);
    Bool_t  GetNonHFEEffiULSLS(AliVTrack *track, AliVTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass);
    void    SwitchPi0EtaWeightCalc(Bool_t fSwitch) {fCalculateWeight = fSwitch;};
    void    SetNonHFEEffi(Bool_t fSwitch) {fCalculateNonHFEEffi = fSwitch;};

private:
    enum{
        kAODanalysis = BIT(20),
    };
    
    AliVEvent   *fVevent;  //!event object
    AliESDEvent *fESD;    //!ESD object
    AliAODEvent *fAOD;    //!AOD object
    AliAODMCHeader *fMCHeader;
    AliPIDResponse *fpidResponse; //!pid response
    AliEMCALGeometry *fEMCALGeo;
    
    Bool_t      fFlagSparse;// switch to THnspare
    Bool_t       fUseTender;// switch to add tender
    
    TClonesArray  *fTracks_tender;//Tender tracks
    TClonesArray  *fCaloClusters_tender;//Tender cluster
    
    AliAODMCParticle  *fMCparticle;//! MC particle
    TClonesArray  *fMCArray;//! MC array
    
    AliMultSelection *fMultSelection;
    
    Bool_t fFlagClsTypeEMC;//switch to select EMC clusters
    Bool_t fFlagClsTypeDCAL;//switch to select DCAL clusters
    
    Int_t   fcentMim; // mim. centrality
    Int_t   fcentMax; // max. centrality
    TString fCentralityEstimator;         // Centrality Estimator
    
    Double_t            fTPCnSigma;//!
    Double_t            fTPCnSigmaMin;//!
    Double_t            fTPCnSigmaMax;//!
    Double_t            fM02Min;//!
    Double_t            fM02Max;//!
    Double_t            fM20Min;//!
    Double_t            fM20Max;//!
    Double_t            fEovPMin;//!
    Double_t            fEovPMax;//!
    Int_t               fNEle;//!
    Double_t            fTPCnSigmaHadMin;//
    Double_t            fTPCnSigmaHadMax;//
    Double_t            fInvmassCut;//
    
    Bool_t              fCalculateWeight;//
    Bool_t              fCalculateNonHFEEffi;//
    Int_t               fNTotMCpart; //! N of total MC particles produced by generator
    Int_t               fNpureMC;//! N of particles from main generator (Hijing/Pythia)
    Int_t               fNembMCpi0; //! N > fNembMCpi0 = particles from pi0 generator
    Int_t               fNembMCeta; //! N > fNembMCeta = particles from eta generator
    Bool_t              fIsFrmEmbPi0;//!
    Bool_t              fIsFrmEmbEta;//!
    Int_t               ftype;//!
    Double_t            fWeight;//!
    
    TF1                 *fPi0Weight;//!
    TF1                 *fEtaWeight;//!
    
    TList       *fOutputList; //!Output list
    TH1F        *fNevents;//! no of events
    TH1F        *fCent;//! centrality
    TH2F        *fMult;//! track multiplicity vs centrality
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
    TH2F        *fHistoNCells;//! No of cells per cluster
    TH2F        *fHistoEperCell;//! Energy per cell
    TH2F        *fHistoCalCell;//! Cells address and amp
    TH2F        *fHistoTimeEMC;//! cluster time
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
    TH2F        *fHistEop;//!pt vs E/p
    TH2F        *fHistMcEopEle;//!pt vs E/p
    TH2F        *fHistMcEopHad;//!pt vs E/p
    TH2F        *fM20;//!M20 vs pt
    TH2F        *fM02;//!M20 vs pt
    TH2F        *fM20EovP;//!M20 vs E/p
    TH2F        *fM02EovP;//!M20 vs E/p
    TH1F        *fEleCanITShit;//!ele cand ITS hit map
    TH2F        *fInvmassULSPt;//!Invmass of ULS
    TH2F        *fInvmassLSPt;//!Invmass of LS
    TH2F        *fHFmomCorr;//!
    TH2F        *fEMCTrkMatch_Phi;//!
    TH2F        *fEMCTrkMatch_Eta;//!
    
    TH1F                *fInclsElecPt;//!
    TH1F                *fHadPt_AftEID;//!
    TH2F                *fHistEop_AftEID;//!
    TH1F                *fNElecInEvt;//!
    TH1F                *fULSElecPt;//!
    TH1F                *fLSElecPt;//!
    
    TH1F                *fRealInclsElecPt;//!
    TH1F                *fNonHFeTrkPt;//!
    TH1F                *fMissingEmbEtaEleTrkPt;//!
    TH1F                *fNonHFeEmbAllTypeTrkPt;//!
    TH1F                *fNonHFeEmbTrkPt;//!
    TH1F                *fNonHFeEmbWeightTrkPt;//!
    TH1F                *fPi0eEmbWeightTrkPt;//!
    TH1F                *fEtaeEmbWeightTrkPt;//!
    TH1F                *fRecoNonHFeTrkPt;//!
    TH1F                *fRecoNonHFeEmbTrkPt;//!
    TH1F                *fRecoNonHFeEmbWeightTrkPt;//!
    TH1F                *fRecoPi0eEmbWeightTrkPt;//!
    TH1F                *fRecoEtaeEmbWeightTrkPt;//!
    TH1F                *fNonHFePairInvmassLS;//!
    TH1F                *fNonHFePairInvmassULS;//!
    TH1F                *fNonHFeEmbInvmassLS;//!
    TH1F                *fNonHFeEmbInvmassULS;//!
    TH1F                *fNonHFeEmbWeightInvmassLS;//!
    TH1F                *fNonHFeEmbWeightInvmassULS;//!
    TH1F                *fPi0EmbInvmassLS;//!
    TH1F                *fPi0EmbInvmassULS;//!
    TH1F                *fPi0EmbWeightInvmassLS;//!
    TH1F                *fPi0EmbWeightInvmassULS;//!
    TH1F                *fEtaEmbInvmassLS;//!
    TH1F                *fEtaEmbInvmassULS;//!
    TH1F                *fEtaEmbWeightInvmassLS;//!
    TH1F                *fEtaEmbWeightInvmassULS;//!
    TH1F                *fRecoLSeEmbTrkPt;//!
    TH1F                *fRecoLSeEmbWeightTrkPt;//!
    TH1F                *fRecoPi0LSeEmbWeightTrkPt;//!
    TH1F                *fRecoEtaLSeEmbWeightTrkPt;//!
    TH1F                *fRecoULSeEmbTrkPt;//!
    TH1F                *fRecoULSeEmbWeightTrkPt;//!
    TH1F                *fRecoPi0ULSeEmbWeightTrkPt;//!
    TH1F                *fRecoEtaULSeEmbWeightTrkPt;//!
    
    THnSparse  *fSparseElectron;//!Electron info
    Double_t *fvalueElectron;//!Electron info
    
    THnSparse           *fSprsPi0EtaWeightCal;//!
    
    AliAnalysisTaskHFEBESpectraEMC(const AliAnalysisTaskHFEBESpectraEMC&); // not implemented
    AliAnalysisTaskHFEBESpectraEMC& operator=(const AliAnalysisTaskHFEBESpectraEMC&); // not implemented
    
    ClassDef(AliAnalysisTaskHFEBESpectraEMC, 1); // example of analysis
};

#endif

