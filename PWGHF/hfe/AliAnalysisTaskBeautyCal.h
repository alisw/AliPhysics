#ifndef AliAnalysisTaskBeautyCal_cxx
#define AliAnalysisTaskBeautyCal_cxx

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

class AliAnalysisTaskBeautyCal : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskBeautyCal();
    AliAnalysisTaskBeautyCal(const char *name);
    virtual ~AliAnalysisTaskBeautyCal();
    
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
    
    void SetClusterTypeEMC(Bool_t flagClsEMC) {fFlagClsTypeEMC = flagClsEMC;};
    void SetClusterTypeDCAL(Bool_t flagClsDCAL) {fFlagClsTypeDCAL = flagClsDCAL;};
    
    void SetCentralityMim(Int_t centMim) {fcentMim = centMim;};
    void SetCentralityMax(Int_t centMax) {fcentMax = centMax;};

    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
    void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagULSElec, Bool_t &fFlagLSElec);
    void ElectronAway(Int_t itrack, AliVTrack *track);
    void SetThresholdEG2(Int_t threshold) { fThresholdEG2=threshold; };
    void SetThresholdEG1(Int_t threshold) { fThresholdEG1=threshold; };
    void FindPatches(Bool_t &hasfiredEG1,Bool_t &hasfiredEG2,Double_t emceta, Double_t emcphi);
    void FindMother(AliAODMCParticle* part, int &label, int &pid);
    Bool_t IsDdecay(int mpid);
    Bool_t IsBdecay(int mpid);
    Bool_t IsPdecay(int mpid);

    void SetHFECuts(AliHFEcuts * const hfecuts) {fhfeCuts = hfecuts;};

private:
    enum{
        kAODanalysis = BIT(20),
    };
    
    AliVEvent   *fVevent;  //!event object
    AliESDEvent *fESD;    //!ESD object
    AliAODEvent *fAOD;    //!AOD object
    AliPIDResponse *fpidResponse; //!pid response
    AliCFManager 	   	*fCFM;                  //!Correction Framework Manager
    
    Bool_t      fFlagSparse;// switch to THnspare
    Bool_t       fUseTender;// switch to add tender
    
    Bool_t	 fEMCEG1;//EMcal Threshold EG1
    Bool_t	 fEMCEG2;//EMcal Threshold EG2
    
    TClonesArray  *fTracks_tender;//Tender tracks
    TClonesArray  *fCaloClusters_tender;//Tender cluster
    
    AliAODMCParticle 	*fMCparticle;//! MC particle
    AliAODMCParticle 	*fMCparticleAss;//! MC particle
    TClonesArray 	*fMCarray;//! MC array
 
    AliMultSelection *fMultSelection;
   
    TClonesArray *fTriggersInfo;//TClonesArray to access container from EMCalTriggerMaker
    Int_t fThresholdEG2;// Threshold for EG2 trigger in ADC for trigger patches
    Int_t fThresholdEG1;// Threshold for EG1 trigger in ADC for trigger patches
    
    Bool_t fFlagClsTypeEMC;//switch to select EMC clusters
    Bool_t fFlagClsTypeDCAL;//switch to select DCAL clusters
    
    Int_t fcentMim; // mim. centrality
    Int_t fcentMax; // max. centrality

    TList       *fOutputList; //!Output list
    TH1F        *fNevents;//! no of events
    TH1F        *fCent;//! centrality
    TH1F        *fVtxZ;//!Vertex z
    TH1F        *fHistClustE;//!cluster energy
    TH2F        *fHistClustEcent;//!cluster energy
    TH2F        *fEMCClsEtaPhi;//! EMC cluster eta and phi
    TH1F        *fHistClustEEG1;//! Cluster Energy, Trigger patch > ThresholdEG1
    TH2F        *fHistClustEEG1cent;//! Cluster Energy, Trigger patch > ThresholdEG1
    TH1F        *fHistClustEEG2;//! Cluster Energy, Trigger patch > ThresholdEG1
    TH2F        *fHistClustEEG2cent;//! Cluster Energy, Trigger patch > ThresholdEG1
    TH2F        *fEMCClsEtaPhiEG1;//! EMC cluster eta and phi, Trigger patch > ThresholdEG1
    TH2F        *fEMCClsEtaPhiEG2;//! EMC cluster eta and phi, Trigger patch > ThresholdEG2
    TH1F        *fNegTrkIDPt;//!neg track ID
    TH1F        *fTrkPt;//!track pt
    TH1F        *fTrketa;//!track eta
    TH1F        *fTrkphi;//!track phi
    TH2F        *fdEdx;//!dedx vs pt
    TH2F        *fTPCNpts;//!TPC Npoints used for dedx
    TH2F        *fTPCnsig;//!TPC Nsigma
    TH1F        *fHistPtMatch;//!tracks matched to EMCAL
    TH1F        *fClsEAftMatch;//!EMC Cluster energy after track matching
    TH2F        *fClsEtaPhiAftMatch;//!EMC Cluster eta phi distribution after track matching
    TH2F        *fClsEtaPhiAftMatchEMCin;//!EMC Cluster eta phi distribution after track matching inside EMC phi acceptance
    TH2F        *fClsEtaPhiAftMatchEMCout;//!EMC Cluster eta phi distribution after track matching outside EMC phi acceptance
    TH2F        *fHistdEdxEop;//!E/p vs dedx
    TH2F        *fHistNsigEop;//!E/p vs dedx
    TH2F        *fHistEop;//!pt vs E/p
    TH2F        *fM20;//!M20 vs pt
    TH2F        *fM02;//!M20 vs pt
    TH2F        *fM20EovP;//!M20 vs E/p
    TH2F        *fM02EovP;//!M20 vs E/p
    TH1F        *fInvmassULS;//!Invmass of ULS
    TH1F        *fInvmassLS;//!Invmass of LS
    TH1F        *fMCcheckMother;
    
    THnSparse  *fSparseElectron;//!Electron info
    Double_t *fvalueElectron;//!Electron info
    
    TH1D        *fHistPhoReco0;//!ele cand SPD or
    TH1D        *fHistPhoReco1;//!ele cand SPD or
    TH1D        *fHistPhoReco2;//!ele cand SPD or

    TH2D        *fHistDCAinc;//!ele cand SPD or
    TH2D        *fHistDCApho;//!ele cand SPD or
    TH2D        *fHistDCAcomb;//!ele cand SPD or
    TH2D        *fHistDCAhfe;//!ele cand SPD or
    TH2D        *fHistDCAhad;//!ele cand SPD or

    TH2D        *fHistDCAde;//!ele cand SPD or
    TH2D        *fHistDCAbe;//!ele cand SPD or
    TH2D        *fHistDCApe;//!ele cand SPD or
    TH2D        *fHistDCAdeSemi;//!ele cand SPD or
    TH2D        *fHistDCAbeSemi;//!ele cand SPD or
    TH2D        *fHistDCApeSemi;//!ele cand SPD or
    TH2D        *fHistDCAhaSemi;//!ele cand SPD or
    TH1D        *fHistOtherEle;//!ele cand SPD or

    TH1D        *fHistHFEcorr;//!ele cand SPD or

    AliHFEcuts  *fhfeCuts;

    AliAnalysisTaskBeautyCal(const AliAnalysisTaskBeautyCal&); // not implemented
    AliAnalysisTaskBeautyCal& operator=(const AliAnalysisTaskBeautyCal&); // not implemented
    
    ClassDef(AliAnalysisTaskBeautyCal, 1); // example of analysis
};

#endif
