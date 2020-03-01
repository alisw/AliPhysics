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
    enum EnhanceSigOrNot {kMB,kEnhance};
    enum pi0etaType {kNoMother, kNoFeedDown, kNotIsPrimary, kLightMesons, kBeauty, kCharm};
    
    AliAnalysisTaskHFEBESpectraEMC();
    AliAnalysisTaskHFEBESpectraEMC(const char *name);
    virtual ~AliAnalysisTaskHFEBESpectraEMC();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    void IsAnalysispp(Bool_t isPP) {fIsAnapp = isPP;};
    
    Bool_t GetEMCalTriggerEG1() { return fEMCEG1; };
    Bool_t GetEMCalTriggerEG2() { return fEMCEG2; };
    void SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;};
    void SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;};
    Bool_t GetEMCalTriggerDG1() { return fDCalDG1; };
    Bool_t GetEMCalTriggerDG2() { return fDCalDG2; };
    void SetEMCalTriggerDG1(Bool_t flagTr1) { fDCalDG1=flagTr1; fDCalDG2=kFALSE;};
    void SetEMCalTriggerDG2(Bool_t flagTr2) { fDCalDG2=flagTr2; fDCalDG1=kFALSE;};
    
    void FindPatches(Bool_t &hasfiredEG1,Bool_t &hasfiredEG2,Double_t emceta, Double_t emcphi);
    void SetThresholdEG2(Int_t threshold) { fThresholdEG2=threshold; };
    void SetThresholdEG1(Int_t threshold) { fThresholdEG1=threshold; };

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
    
    void    SetTrackMatchPar(Double_t deltaEta, Double_t deltaPhi){fDeltaEta = deltaEta; fDeltaPhi = deltaPhi;};
    void    SetM02Cut(Double_t m02Min, Double_t m02Max1, Double_t m02Max2) {fM02Min = m02Min; fM02Max1 = m02Max1; fM02Max2 = m02Max2;};
    void    SetM20Cut(Double_t m20Min, Double_t m20Max) {fM20Min = m20Min; fM20Max = m20Max;};
    void    SetEovPCut(Double_t eovpMin, Double_t eovpMax) {fEovPMin = eovpMin; fEovPMax = eovpMax;};

    
    Bool_t  PassEIDCuts(AliVTrack *track, AliVCluster *clust, Bool_t &Hadtrack);
    
    Bool_t  GetNMCPartProduced();
    void    GetPi0EtaWeight(THnSparse *SparseWeight);
    Int_t   GetPi0EtaType(AliAODMCParticle *part);
    
    Bool_t  IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromMB, Int_t &type, Int_t &iMom, Int_t &MomPDG, Double_t &MomPt);
    
    Bool_t  GetNonHFEEffiDenom(AliVTrack *track);
    Bool_t  GetNonHFEEffiRecoTag(AliVTrack *track);
    Bool_t  GetNonHFEEffiULSLS(AliVTrack *track, AliVTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass);
    
    void    SwitchPi0EtaWeightCalc(Bool_t fSwitch) {fCalculateWeight = fSwitch;};
    void    SetNonHFEEffi(Bool_t fSwitch) {fCalculateNonHFEEffi = fSwitch;};
    void    SetElecRecoEffi(Bool_t fSwitch) {fCalculateElecRecoEffi = fSwitch;};
    void    SwitchMCTemplateWeightCalc(Bool_t fSwitch) {fCalculateMCTemplWeightCalc = fSwitch;};
    void    SwitchFillMCTemplate(Bool_t fSwitch) {fFillMCTemplates = fSwitch;};

    void    GetElectronFromStack();
    void    GetTrackHFStatus(AliVTrack *track, Bool_t &IsMCEle, Bool_t &IsMCHFEle, Bool_t &IsMCBEle, Bool_t &IsMCDEle);
    void    GetEIDRecoEffi(AliVTrack *track, AliVCluster *clust, Bool_t IsMCEle, Bool_t IsMCHFEle, Bool_t IsMCBEle, Bool_t IsMCDEle);

    void    GetMCTemplateWeight();
    Bool_t  GetMCDCATemplates(AliVTrack *track, Double_t TrkDCA);
    
   // void    InputWeightCorrectionMaps();
    void SetDmesonWeightHist(TH1 *D1, TH1 *D2, TH1 *D3);
    void SetBmesonWeightHist(TH1 *B1, TH1 *B2, TH1 *B3);
    void    GetBWeight(AliAODMCParticle *Part, Double_t &BCentWeight, Double_t &BMinWeight, Double_t &BMaxWeight);
    void    GetDWeight(AliAODMCParticle *Part, Double_t &DCentWeight, Double_t &DMinWeight, Double_t &DMaxWeight);
    
    void SetDmesonWeightHistPbPb(TH1 *D0, TH1 *DPlus, TH1 *Ds, TH1 *Lc);
    void SetBmesonWeightHistPbPb(TH1 *B);
    void GetDWeightPbPb(AliAODMCParticle *Part, Int_t PDG, Double_t &DCentWeight);
    void GetBWeightPbPb(AliAODMCParticle *Part, Double_t &BCentWeight);

    void    SwitchRecalImpPar(Bool_t fSwitch) {fRecalIP = fSwitch;};
    void    RecalImpactParam(const AliVTrack * const track, Double_t dz[2], Double_t covar[3]);
    AliAODVertex*   RemoveDaughtersFromPrimaryVtx(const AliVTrack * const track);
    
    
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
    
    Bool_t    fEMCEG1;//EMcal Threshold EG1
    Bool_t    fEMCEG2;//EMcal Threshold EG2
    Bool_t    fDCalDG1;//DCal Threshold DG1
    Bool_t    fDCalDG2;//DCal Threshold DG2
    
    TClonesArray *fTriggersInfo;//TClonesArray to access container from EMCalTriggerMaker
    Int_t fThresholdEG2;// Threshold for EG2 trigger in ADC for trigger patches
    Int_t fThresholdEG1;// Threshold for EG1 trigger in ADC for trigger patches
    
    AliAODMCParticle  *fMCparticle;//! MC particle
    TClonesArray  *fMCArray;//! MC array
    
    AliMultSelection *fMultSelection;
    Bool_t  fIsAnapp;// Is analysis pp

    Bool_t fFlagClsTypeEMC;//switch to select EMC clusters
    Bool_t fFlagClsTypeDCAL;//switch to select DCAL clusters
    
    Int_t   fcentMim; // mim. centrality
    Int_t   fcentMax; // max. centrality
    TString fCentralityEstimator;         // Centrality Estimator
    
    Bool_t              fRecalIP;//
    
    Double_t            fDeltaEta;//
    Double_t            fDeltaPhi;//
    Double_t            fTPCnSigma;//
    Double_t            fTPCnSigmaMin;//
    Double_t            fTPCnSigmaMax;//
    Double_t            fM02Min;//
    Double_t            fM02Max1;//
    Double_t            fM02Max2;//
    Double_t            fM20Min;//
    Double_t            fM20Max;//
    Double_t            fEovPMin;//
    Double_t            fEovPMax;//
    Int_t               fNEle;//!
    Double_t            fTPCnSigmaHadMin;//
    Double_t            fTPCnSigmaHadMax;//
    Double_t            fInvmassCut;//
    
    Bool_t              fCalculateWeight;//
    Bool_t              fCalculateNonHFEEffi;//
    Bool_t              fCalculateElecRecoEffi;//
    Bool_t              fCalculateMCTemplWeightCalc;//
    Bool_t              fFillMCTemplates;//
    Int_t               fNTotMCpart; //! N of total MC particles produced by generator
    Int_t               fNpureMC;//! N of particles from main generator (MB/Enhanced)
    Int_t               fNembMCpi0; //! N > fNembMCpi0 = particles from pi0 generator
    Int_t               fNembMCeta; //! N > fNembMCeta = particles from eta generator
    Bool_t              fIsFrmEmbPi0;//!
    Bool_t              fIsFrmEmbEta;//!
    Int_t               ftype;//!
    Double_t            fWeight;//!
    Double_t            fWeightPi0;//!
    Double_t            fWeightEta;//!

    TF1                 *fPi0Weight;//!
    TF1                 *fEtaWeight;//!
    Int_t               fnBinsDCAHisto;//!
    Double_t            fTrkDCA;//!
    
    TH1F                *fDcent;//
    TH1F                *fDUp;//
    TH1F                *fDDown;//
    TH1F                *fBcent;//
    TH1F                *fBMin;//
    TH1F                *fBMax;//
    TH1F                *fD0;//
    TH1F                *fDPlus;//
    TH1F                *fDs;//
    TH1F                *fLc;//
    TH1D                *fB;//
    Double_t            fWeightB;//!
    Double_t            fWeightBMin;//!
    Double_t            fWeightBMax;//!
    Double_t            fWeightD;//!
    Double_t            fWeightDUp;//!
    Double_t            fWeightDDown;//!

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
    TH2F                *fHadEovp_AftEID;//!
    TH2F                *fHadEovpNL_AftEID;//!
    TH2F                *fEop_AftEID;//!
    TH2F                *fEopNL_AftEID;//!
    TH1F                *fNElecInEvt;//!
    TH1F                *fULSElecPt;//!
    TH1F                *fLSElecPt;//!
    TH2F                *fHadDCA;//!
    TH2F                *fInclElecDCA;//!
    TH2F                *fULSElecDCA;//!
    TH2F                *fLSElecDCA;//!

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
    
    TH1F                *fInclElePhysPriAll;//!
    TH1F                *fHFEPhysPriAll;//!
    TH1F                *fBEPhysPriAll;//!
    TH1F                *fDEPhysPriAll;//!
    TH1F                *fInclElePhysPriTrkCuts;//!
    TH1F                *fHFEPhysPriTrkCuts;//!
    TH1F                *fBEPhysPriTrkCuts;//!
    TH1F                *fDEPhysPriTrkCuts;//!
    TH1F                *fInclElePhysPriEMCMatch;//!
    TH1F                *fHFEPhysPriEMCMatch;//!
    TH1F                *fBEPhysPriEMCMatch;//!
    TH1F                *fDEPhysPriEMCMatch;//!
    TH1F                *fInclElePhysPriTPCnsig;//!
    TH1F                *fHFEPhysPriTPCnsig;//!
    TH1F                *fBEPhysPriTPCnsig;//!
    TH1F                *fDEPhysPriTPCnsig;//!
    TH1F                *fInclElePhysPriEovPBfrSS;//!
    TH1F                *fHFEPhysPriEovPBfrSS;//!
    TH1F                *fBEPhysPriEovPBfrSS;//!
    TH1F                *fDEPhysPriEovPBfrSS;//!
    TH1F                *fInclElePhysPriSS;//!
    TH1F                *fHFEPhysPriSS;//!
    TH1F                *fBEPhysPriSS;//!
    TH1F                *fDEPhysPriSS;//!
    TH1F                *fInclElePhysPriEovP;//!
    TH1F                *fHFEPhysPriEovP;//!
    TH1F                *fBEPhysPriEovP;//!
    TH1F                *fDEPhysPriEovP;//!
    
    TH1F                *fBHadpT;//!
    TH1F                *fBMesonpT;//!
    TH1F                *fBDHadpT;//!
    TH1F                *fDHadpT;//!
    TH1F                *fDMesonpT;//!
    TH1F                *fD0pT;//!
    TH1F                *fDPluspT;//!
    TH1F                *fDspT;//!
    TH1F                *fLambdaCpT;//!
    
    TH2F                *fDElecDCA;//!
    TH2F                *fBElecDCA;//!
    TH2F                *fBHadElecDCA;//!
    TH2F                *fBMesonElecDCA;//!
    TH2F                *fBBaryonElecDCA;//!
    TH2F                *fDHadElecDCA;//!
    TH2F                *fDMesonElecDCA;//!
    TH2F                *fDBaryonElecDCA;//!
    TH2F                *fLambdaCElecDCA;//!
    TH2F                *fD0ElecDCA;//!

    THnSparse  *fSparseElectron;//!Electron info
    Double_t *fvalueElectron;//!Electron info
    
    THnSparse           *fSprsPi0EtaWeightCal;//!
    THnSparse           *fSprsTemplatesNoWeight;//!
    THnSparse           *fSprsTemplatesWeight;//!
    THnSparse           *fSprsTemplatesWeightVar1;//!
    THnSparse           *fSprsTemplatesWeightVar2;//!

    AliAnalysisTaskHFEBESpectraEMC(const AliAnalysisTaskHFEBESpectraEMC&); // not implemented
    AliAnalysisTaskHFEBESpectraEMC& operator=(const AliAnalysisTaskHFEBESpectraEMC&); // not implemented
    
    ClassDef(AliAnalysisTaskHFEBESpectraEMC, 1); // example of analysis
};

#endif

