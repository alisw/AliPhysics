#ifndef AliAnalysisTaskHFEEfficiency_cxx
#define AliAnalysisTaskHFEEfficiency_cxx

class TH1F;
class TH2F;
class THnSparse;
class TLorentzVector;

class AliEMCALTrack;
class AliESDEvent;
class AliEMCALGeometry;

class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;

#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"



class AliAnalysisTaskHFEEfficiency : public AliAnalysisTaskSE{
public:
    enum{
        kHasMCdata = BIT(19),
    };
    enum pi0etaType {kNoMother, kNoFeedDown, kNoIsPrimary, kLightMesons, kBeauty, kCharm};
    enum ESourceType {kNoMotherE, kPi0NoFeedDown, kEtaNoFeedDown, kGPi0NoFeedDown, kGEtaNoFeedDown, kDirectGamma, kOthersE};
    enum HijingOr {kHijing,kElse};
    enum EInvSourceType {kPi0,kPi0Both,kEta,kEtaBoth,kGamma,kGammaBoth};
    enum WhichHeavyMeson {kelectronCharm,kelectronBeauty};
    
    AliAnalysisTaskHFEEfficiency();
    AliAnalysisTaskHFEEfficiency(const char *name);
    virtual ~AliAnalysisTaskHFEEfficiency();
    
    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); };
    
    void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
    void SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod); //select centrality
    void CheckCentrality(AliVEvent *event,Bool_t &centralitypass); //to use only events with the correct centrality....
    void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
    //void SelectPhotonicElectron(Int_t itrack, AliESDtrack *track, Bool_t &fFlagPhotonicElec70,Bool_t &fFlagPhotonicElec100,Bool_t &fFlagPhotonicElec500);
    //  void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec70,Bool_t &fFlagPhotonicElec100,Bool_t &fFlagPhotonicElec500);
    void SelectPhotonicElectronR(Int_t itrack, AliVTrack *track, Int_t hijing, Double_t weight, Int_t motherindex, Int_t pdg, Int_t source);
    void SelectPhotonicElectronR_ULSLS(Int_t itrack, AliVTrack *track, Int_t hijing, Double_t weight, Int_t pdg, Int_t source);
    void SetEtaRange(Double_t etaminimum, Double_t etamaximum);

    
    Double_t GetMCweightPi0(Double_t mcPi0pT);
    Double_t GetMCweightPi0tiltUp(Double_t mcPi0pT);
    Double_t GetMCweightPi0tiltDw(Double_t mcPi0pT);
    Double_t GetMCweightEta(Double_t mcEtapT);
    Double_t GetMCweightEtatiltUp(Double_t mcEtapT);
    Double_t GetMCweightEtatiltDw(Double_t mcEtapT);
    
    Double_t GetMCweightPi02040(Double_t mcPi0pT);
    Double_t GetMCweightPi0tiltUp2040(Double_t mcPi0pT);
    Double_t GetMCweightPi0tiltDw2040(Double_t mcPi0pT);
    Double_t GetMCweightEta2040(Double_t mcEtapT);
    Double_t GetMCweightEtatiltUp2040(Double_t mcEtapT);
    Double_t GetMCweightEtatiltDw2040(Double_t mcEtapT);
    
    
    Double_t GiveHFWeight(Double_t HFpt);
    
    
    Bool_t IsFromBGEventAOD(Int_t Index);
    AliHFEpid *GetPID() const { return fPID; }
    
    Bool_t HasMCData() const { return TestBit(kHasMCdata); }
    void SetHasMCData(Bool_t hasMC = kTRUE) { SetBit(kHasMCdata, hasMC); };
    template <typename T> void           PlotVZeroMultiplcities(const T* event) const;
    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    void                                 SetAssoTPCCluster(Int_t tpc_clust) {fAssoTPCCluster = tpc_clust;};
    void                                 SetAssoITSRefit(Bool_t itsref) {fAssoITSRefit = itsref;};
    void                                 SetAssopTmin(Double_t ptminassotrack) {fAssopTMin = ptminassotrack;};
    void                                 SetMassCut(Double_t masscutpair) {fMassCut = masscutpair;};
    void                                 SetStackLoop(Bool_t stack) {fStackloop = stack;};
    Int_t                                GetPrimary(Int_t id, TClonesArray *mcArray);
    Int_t                                GetPi0EtaType(AliAODMCParticle *pi0eta, TClonesArray *mcArray);
    Int_t                                GetElecSourceType(AliAODMCParticle *electron, TClonesArray *mcArray,Double_t &ptm);
    void                                 SetTiltUpWeights(Bool_t UPWeights){tiltup = UPWeights;};
    void                                 SetTiltDwWeights(Bool_t DWWeights){tiltdw = DWWeights;};
    
    void                                 SetCentralWeights(Bool_t mostcentral010){fcentral = mostcentral010;};
    void                                 SetSemicentralWeights(Bool_t semicentral2040){fsemicentral = semicentral2040;};
    
    void                                 SetWeights(Bool_t Weights){WeightsForEnhanced = Weights;};
    Bool_t                               FillNumerator(TClonesArray *mcArray,AliVTrack *track, Double_t ForHFEReco[]);
    void                                 SetIDCuts(Double_t minTOFnSigma, Double_t maxTOFnSigma, Double_t minITSnsigmalowpt, Double_t maxITSnsigmalowpt,Double_t minITSnsigmahighpt, Double_t maxITSnsigmahighpt, Double_t minTPCnsigmalowpt, Double_t maxTPCnsigmalowpt,Double_t minTPCnsigmahighpt, Double_t maxTPCnsigmahighpt );
    void SetTPCS(Int_t sig) {fTPCS = sig;};
    
    void                                 SetWeightsHF(Bool_t Weights){WeightsForHF = Weights;};
    
    
    
private:
    
    enum{
        kAODanalysis = BIT(20),
    };
    
    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
    
    AliESDEvent	*fESD;			//ESD object
    AliAODEvent	*fAOD;			//ESD object
    AliVEvent	   *fVevent;			//ESD object
    
    TList       	*fOutputList;		//output list
    AliESDtrackCuts	*fTrackCuts;		//ESD track cuts
    
    Int_t 		fNoTrks;		//No of tracks
    AliHFEcuts           *fCuts;                 //Cut Collection
    AliCFManager         *fCFM;                  //!Correction Framework Manager
    AliHFEpid            *fPID;                  //PID
    AliHFEpidQAmanager   *fPIDqa;    //! PID QA manager
    Bool_t               fRejectKinkMother;      //Reject Kink Mother
    TH1F		*fNoEvents;		//!no of events
    TH1F 	    *fMCPdgmom;		//!MC pdg of elec mother
    TH1F		*fMCElecPtPhoto;	//!MC electron Pt photonic
    TH1F		*fMCElecPtIncl;		//!MC electron Pt inclusive
    TH1F		*fElecNos;		//!no of electrons
    TH1F		*fInclusiveElecPt;		//!Photonic elec pt
    TH1F		*fPhotoElecPt;		//!Photonic elec pt
    TH1F		*fTruePhotoElecPt;	//!True photonic elec pt
    TH1F		*fMCtagPhotoElecPtAll;	//!All MC tagged photonic elec pt
    TH1F		*fMCtagTruePhotoElecPtAll;//!All True MC tagged photonic elec pt
    TH1F        *fMCtagGammaElecPtAll; //!
    TH1F   	    *fMCtagPi0ElecPtAll;//!
    TH1F   		*fMCtagEtaElecPtAll;//!
    TH1F		*fSingleElecPt;		//!Single elec pt
    TH1F        *fTrueSingleElecPt;	//!True single elec pt
    TH1F    	*fInvmassLS;		//!Inv mass of LS (e,e)
    TH1F		*fInvmassULS;		//!Inv mass of ULS (e,e)
    TH1F       	*fInvmassULSNHFE;    //!Inv mass of NHF (e,e) pairs
    TH1F		*fOpeningAngleLS;	//!opening angle for LS pairs
    TH1F		*fOpeningAngleULS;	//!opening angle for ULS pairs
    TH1F        	 *fTrackPtBefTrkCuts; //!track pt before track cuts
    TH1F        	 *fTrackPtAftTrkCuts; //! track pt after track cuts
    TH2F        	 *fTPCnsigma; //!TPC nsigma
    TH2F        	 *fITSnsigma; //!ITS nsigma
    TH2F        	 *fITSnsigmaElectron; //!ITS nsigma elect
    TH2F        	 *fTOFnsigma; //!TOF nsigma
    TH2F        	 *fTrkEovPBef; //!Track E/p before
    TH2F       		 *fdEdxBef; //!track dEdx before
    TH1F      	 	 *fPhotoElecCandiPdgMom; //! pdg code of Photo elec candi mother
    TH1F       		 *fInvmassULSNHFEAllPartn;  //!Inv mass of photo ele ULS with all pair
    TH1F        	 *fInvmassLSNHFEAllPartn;  //!Inv mass of photo ele LS with all pair
    TH1F                 *fCentralityPass; // ! QA histogram of events that pass centrality cut
    TH1F                 *fCentralityNoPass; //! QA histogram of events that do not pass centrality cut
    THnSparseF           *fPi0EtaSpectra;// ! QA for weights
    THnSparseF           *fHFEDenominator;// ! EFFfor EHFE
    THnSparseF           *fHFENumeratorTRKCUTS_Filt;// ! EFFfor EHFE
    THnSparseF           *fHFENumeratorTRKCUTS_ITSTPC;// ! EFFfor EHFE
    //  THnSparseF           *fHFENumeratorTRKCUTS_ITStrkCut;// ! EFFfor EHFE
    THnSparseF           *fHFENumeratorTRKCUTS;// ! EFFfor EHFE
    THnSparseF           *fHFENumeratorTOF;// ! EFFfor EHFE
    THnSparseF           *fHFENumeratorITS;// ! EFFfor EHFE
    THnSparseF           *fHFENumerator;// ! EFFfor EHFE
    THnSparseF           *fInvMULS;// ! QA for weights
    THnSparseF           *fULS;// ! QA for weights
    THnSparseF           *fLS;// ! QA for weights
    THnSparseF           *fIncl;// ! QA for weights
    Double_t              fCentrality; // event centrality for QA
    Double_t              fCentralityMin; // lower bound of cenrality bin
    Double_t              fCentralityMax; // upper bound of centrality bin
    const char           *fkCentralityMethod; // method used to determine centrality (V0 by default)
    TH1F                 *fVZEROA;//!vezo A multi
    TH1F                 *fVZEROC;//!vezo A multi
    Int_t                fAssoTPCCluster;//asso tpc cluster
    Bool_t               fAssoITSRefit;//asso its refit
    Double_t             fAssopTMin;//assoptmin
    Bool_t               fStackloop;//wvwe
    TH2F                 *fTPCnsigmaAft;//!TPC n sigma vs p after HFE pid
    TH2F                 *fTPCnsigmaVSptAft;   //!lodviow
    TH2F                 *fTPCnsigmaAftITSTOF;//!TPC n sigma vs p
    TH2F                 *Rconv_pT;//!FAKE ITS HITS
    Double_t              fMassCut;//mass cut on pair
    THnSparseF           *fSparseMassULS;//!ssss
    THnSparseF           *fSparseMassLS;//!ssssss
    TH1F                 *fakepT;//!ve
    TH1F                 *fakepTgraterFive;//!ee
    Bool_t                WeightsForEnhanced;//set the weigh
    TH1F                 *fNoEventsStackHFE;		//!no of events
    Double_t              fminITSnsigmaLowpT;  //ID cuts its
    Double_t              fmaxITSnsigmaLowpT;  //ID cuts its
    Double_t              fminITSnsigmaHighpT;  //ID cuts its
    Double_t              fmaxITSnsigmaHighpT;  //ID cuts its
    Double_t              fminTPCnsigmaLowpT;  //ID cuts tpc
    Double_t              fmaxTPCnsigmaLowpT;  //ID cuts tpc
    Double_t              fminTPCnsigmaHighpT;  //ID cuts tpc
    Double_t              fmaxTPCnsigmaHighpT;  //ID cuts tpc
    Double_t              fminTOFnSigma;  //ID cuts tof
    Double_t              fmaxTOFnSigma;//ID cuts tof
    Int_t                 fTPCS;//tpc signal cluster
    //    TH1F                 *fhHFStackBeauty;//!TPC n sigma vs p after HFE pid
    //    TH1F                 *fhHFStackCharm;//!TPC n sigma vs p after HFE pid
    Bool_t                tiltup;//tilting for weights
    Bool_t                tiltdw;//tilting for weights
    Bool_t                fcentral;//select weights for 010
    Bool_t                fsemicentral;//select weights for 2040
    Bool_t                WeightsForHF;//set the weigh
    Double_t             fmineta;  //eta cuts
    Double_t             fmaxeta;//eta cuts
    
    
    AliAnalysisTaskHFEEfficiency(const AliAnalysisTaskHFEEfficiency&); // not implemented
    AliAnalysisTaskHFEEfficiency& operator=(const AliAnalysisTaskHFEEfficiency&); // not implemented
    
    ClassDef(AliAnalysisTaskHFEEfficiency,4);//!example of analysis
    
};

#endif
