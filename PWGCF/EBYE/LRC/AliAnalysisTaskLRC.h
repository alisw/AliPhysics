#ifndef AliAnalysisTaskLRC_H
#define AliAnalysisTaskLRC_H

// Analysis task for Long Range Correlation (LRC) analysis using TPC data
// This includes a TList of AliLRCProcess objects that are processing LRC analysis
// for a given eta-phi windows

// Author : Andrey Ivanov, Igor Altsybeev, St.Peterburg State University
// Email: Igor.Altsybeev.ch


#include <AliPIDResponse.h>
#include <AliPIDCombined.h>
#include <AliAnalysisTaskSE.h>

//#define kNumberOfParentParticleClassesInMC 8


//class AliLRCProcess;
class AliLRCBase;
class AliESDtrackCuts;
class TH1D;
class TH2D;
class TH1I;
class TRandom3;
class TParticle;
//class AliSimpleEvent;
class TTree;

class TStopwatch;
//enum en_AnalysisType
//{
//    en_AnalysisType_ESD = 0,
//    en_AnalysisType_AOD

//};

class AliAnalysisTaskLRC : public AliAnalysisTaskSE {

public:
    //Constructors
    AliAnalysisTaskLRC( const char *name = "AliAnalysisTaskLRC", Bool_t runKine = kFALSE );
    virtual ~AliAnalysisTaskLRC() {}

    //AliAnalysisTaskSE overloading

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
//    virtual void   UserExecLoop( Double_t phiAdditional = 0 );//Option_t *option);
    virtual void   Terminate(Option_t *);
    //----------------------------------

    void AddLRCProcess(AliLRCBase *newProc); //Adds new AliLRCProcess to analysis task

    // Setters
    void SetMinNumberOfSPDtracklets( Int_t MinSPDtracklets );   //Sets  min number of SPD tracklets
    void SetMaxPtLimit(Double_t MaxPtLimit);   //Sets  Max Pt filter
    void SetMinPtLimit(Double_t MinPtLimit);   //Sets  Min Pt filter
    void SetCheckForkVtx(Bool_t CheckForkVtx){fCheckForkVtx=CheckForkVtx;} // Accept only events with veretex
    void SetCheckForVtxPosition(Bool_t CheckForVtxPosition ){fCheckForVtxPosition=CheckForVtxPosition;} //Accept only events with veretex in slected range
    void SetTrackCuts(AliESDtrackCuts* const cuts)  { fEsdTrackCuts = cuts; }
    void SetAODtrackCutBit(Int_t bit){ fAODtrackCutBit = bit;  } //AOD track cut bit
    void SetShowEventStats(Bool_t ShowEventStats)  {fShowEventStats= ShowEventStats;}
    void SetShowPerTrackStats(Bool_t ShowPerTrackStats) {fShowPerTrackStats=ShowPerTrackStats;}
    void SetVtxDiamond(Double_t Vx, Double_t Vy, Double_t Vz) {fVxMax = Vx;fVyMax =Vy;fVzMax = Vz;}
    void SetNchCuts(Int_t minNch, Int_t maxNch){fMinAcceptedTracksCut=minNch; fMaxAcceptedTracksCut=maxNch;}

    void SetNumberOfPhiSectors(Int_t nSectors){ fNumberOfPhiSectors = nSectors; }//fNeedToRotateSector = kTRUE; }
    void SetCentralityClass(Float_t minCentralityClass, Float_t maxCentralityClass ){ fMinCentralityClass = minCentralityClass; fMaxCentralityClass = maxCentralityClass; }
    
    void SetIonsAnalysis(Bool_t isIonsFlag ){ fIsIonsAnalysis = isIonsFlag; }
    void SetEtAnalysis(Bool_t isEtAnalysisFlag ){ fEtInsteadOfPt = isEtAnalysisFlag; }
    void SetArtificialInefficiencyCoeff( Double_t artificialInefficiencyCoeff ) { fArtificialInefficiency = artificialInefficiencyCoeff; }   //Sets coeff for artificial inefficiency

    //void SetNumberOfPhiSectorsByHand( Int_t numberOfPhiSectorsByHand ) { fNumberOfPhiSectorsByHand = numberOfPhiSectorsByHand; }

    // Getters
    TList* GetListOfProcessors() { return &fLRCproc;} // Returns list of included
    AliESDtrackCuts* GetTrackCuts() const                         { return fEsdTrackCuts; }
    AliLRCBase * Proc(Int_t index);// Get Processor i

    void SetParticleTypeForTask( char* strF, char* strB );
//    void SetMCparticleClassForFillingLRC( TString strParticleType ) { fStrMCparticleClassForFillingLRC = strParticleType; }
//    void SetEtaCutsForSpecMCanalysis( double etaMin, double etaMax  ) { fEtaMCanalysisCutMin = etaMin; fEtaMCanalysisCutMax = etaMax; }

    void SetV0ACMultThreshold( int minMult ) { fThresholdOnV0mult = minMult; }

    //    void SetIncludeEventTreeInOutput( Bool_t flag ) { fSetIncludeEventTreeInOutput = flag; }
    //    Bool_t GetIncludeEventTreeInOutput() { return fSetIncludeEventTreeInOutput; }

    //    void SetAnalysisType( en_AnalysisType analysisType ) { fAnalysisType = analysisType; }
    //    en_AnalysisType GetAnalysisType() { return fAnalysisType; }

    Double_t GetEventPlane(AliVEvent *event);


    //track cuts array stuff
    void AddTrackCutForBits(AliESDtrackCuts* const cuts, TString cutsName );
    void SetUseListCuts( Bool_t const useCutsList )  { fSwitchToListingCuts = useCutsList; }
    Int_t GetNumberOfTrackCutForBits() const                         { return fArrTrackCuts.GetEntries();/*fNumberOfCutsToRemember*/; }

    void SetAnalysisLevel(const char* analysisLevel) {  fAnalysisLevel = analysisLevel;}
    const char* GetAnalysisLevel() {return fAnalysisLevel.Data();}

    void SetFlagWatchZDC( Bool_t flagWatchZDC) {  fFlagWatchZDC = flagWatchZDC;}
    void SetFlagWatchV0 ( Bool_t flagWatchV0 ) {  fFlagWatchV0  = flagWatchV0 ;}
    void SetFlagWatchFMD( Bool_t flagWatchFMD) {  fFlagWatchFMD = flagWatchFMD;}
    Bool_t GetFlagWatchZDC() {  return fFlagWatchZDC; }
    Bool_t GetFlagWatchV0 () {  return fFlagWatchV0 ; }
    Bool_t GetFlagWatchFMD() {  return fFlagWatchFMD; }

    enum enTaskObjectParameters { kMaxParticlesNumber = 10000, kMaxLRCprocArrayPointers = 1000 }; // default TPC & TOF pid (via GetTPCpid & GetTOFpid)

protected:
    void SetParticleTypeToProcessors( int windowId, char* strPid );
    
    // Track cuts
    TString fAnalysisLevel; //ESD, AOD or MC

    AliESDtrackCuts *fEsdTrackCuts;               // esd track cuts
    Int_t fAODtrackCutBit;//track cut bit from track selection (only used for AODs)

    Int_t fNumberOfPhiSectors; // n of phi rotations
    AliLRCBase *fLRCprocArrayPointers[kMaxLRCprocArrayPointers];

    //arrays with data for LRC processors
    float fArrayTracksPt[kMaxParticlesNumber];
    float fArrayTracksEta[kMaxParticlesNumber];
    float fArrayTracksPhi[kMaxParticlesNumber];
    Short_t fArrayTracksCharge[kMaxParticlesNumber];
    Int_t fArrayTracksPID[kMaxParticlesNumber];


    // Array with different track cuts to remember in simple event Tree
    TObjArray fArrTrackCuts;    //AliESDtrackCuts*  [100];     // Arr with different track cuts
    TString  fArrCutsNames[100];     // Arr with names of different track cuts
    TH1I    *fHistCutsNamesBins;        //!  tracks passed different cut sets in histogram bins
    Bool_t fSwitchToListingCuts;      // switch to remember cuts desicions and not to drop track by fTrackCuts


    //    en_AnalysisType fAnalysisType; // type of analysis

    //SPD tracklets cut
    Int_t fMinNumberOfSPDtracklets;   //Minimum number of SPD tracklets in ESD event

    // Acceptance cuts
    Double_t fMaxPtLimit;  //Max Pt filter
    Double_t fMinPtLimit;  // Min Pt filter

    // Nch cuts
    Int_t fMinAcceptedTracksCut;   //Minimum number of accepted tracks in event
    Int_t fMaxAcceptedTracksCut;   //Maximum number of accepted tracks in event

    // Vtx cuts
    Bool_t fCheckForkVtx;		// Check for vertex
    Bool_t fCheckForVtxPosition;  // Check if vertex position in range
    Double_t fVxMax;	// X vrtx max
    Double_t fVyMax;	// Y vrtx max
    Double_t fVzMax;	// Z vrtx max


    TList fLRCproc;       //  AliLRCProcess objects list
    TList* fOutList;      //! Task Output data container

    Bool_t fRunKine;      // ESD/AOD  - KINE switch
    Bool_t fShowEventStats; //  Allows per event debug output (trigger Nch, cuts etc)
    Bool_t fShowPerTrackStats; // Allows per track debug output



    // QA histos

    TH1I *fHistEventCutStats;  //! Event cut statistics
    TH1I *fHistTrackCutStats;  //! Track cut statistics
    TH1I *fHistAODTrackStats;  //! AOD track bits statistics


    TH1D *fHistVx;  //!Vx hist
    TH1D *fHistVy;  //!Vy hist
    TH1D *fHistVz;  //!Vz hist

    TH1I *fHistVertexNconributors;  //!vertex contributors number

    TH2F *fHistEventPlane; //event plane distribution


    TH1F *fHistPt; //! Overal Pt spectrum
    TH1F *fHistEta; //! Overal Eta spectrum
    TH1F *fHistPhi; //! Overal Phi spectrum
    TH2D *fHistEtaPhi;       //! 2D plot for checking acceptance
    
    TH2D *fHistEtaVsZvCoverage; //! Statistics on tracks Zv and Eta for all tracks
    TH2D *fHistEtaVsZvCoverageAccepted; //!  Statistics on tracks Zv and Eta for accepted tracks

    TH1D *fHistMultBeforeCuts;   //! Histo: Number of tracks before applying cuts
    TH1D *fHistAcceptedMult;   //! Number of accepted tracks histo
    TH1D *fHistAcceptedTracks;   //! Number of tracks accepted for filling LRC processors, histo
    TH1D *fHistMultiplicityInEtaRegion; //! Number of tracks in |eta|<1
    TH1D *fHistAcceptedTracksAfterPtCuts;   //! Number of tracks accepted for filling LRC processors, histo
    TH1D *fHistAcceptedTPCtracks;   //! Number of accepted tracks with TPC inner param
    TH1D *fHistClustersTPC;   //! Number of TPC clusters distribution
    TH1D *fHistClustersTPCafterCuts;   //! Number of TPC clusters distribution after cuts
    TH1D *fHistCrossedRowsTPC;   //! Number of TPC crossed rows
    TH1D *fHistCrossedRowsTPCafterCuts;   //! Number of TPC crossed rows after cuts
    

    TH1D *fHistClustersITS;   //! Number of ITS clusters distribution
    TH1D *fHistTrackletsITS;   //! Number of ITS tracklets distribution
    TH2D *fHist2DClustersTPCvsPt;   //! Number of TPC clusters vs Pt distribution (to see the dependence!)
    TH2D *fHist2DClustersTPCvsEta;   //! Number of TPC clusters vs Eta distribution (to see the dependence!)

    TH2D *fHist2DAcceptedTracksPtvsEta;   //! rejected tracks pt vs eta

    TH1D *fHistMClabels;   //! MC labels
    TH1D *fHistRejectedTracksCharge;   //! Charge of rejected tracks
    TH1D *fHistTracksCharge;   //! Charge of accepted tracks (zero is filled only for MC truth)

    AliAnalysisTaskLRC(const AliAnalysisTaskLRC&); // not implemented
    AliAnalysisTaskLRC& operator=(const AliAnalysisTaskLRC&); // not implemented

    TH1D *fHistProbabilitiesPID;  //!hist of esd pid prob-s
    //Double_t *fProbabilitiesPID;	//! array of esd pid prob-s
    TH1D *fHistESDtrackMass;  //!hist of esd particle masses
    TH1D *fHistProbabilityPion;  //!hist of pion probability
    TH1D *fHistProbabilityKaon;  //!hist of kaon probability
    TH1D *fHistProbabilityProton;  //!hist of proton probability
    TH1D *fHistParticlesDistr;  //!hist of particles distr
    TH1D *fHistParticlesDistrBeforeCuts;  //!hist of particles distr



    TH1D *fHistCentralityPercentile;        //! centrality class
    TH1D *fHistCentralityClass10;           //! centrality class by 10
    TH1D *fHistCentralityClass5;            //! centrality class by 5


    //ZDC stuff
    TH1D *fHistZDCenergy[5];     //! ZDC energy for diff mult conditions
    TH1D *fHistZDCparticipants;             //! ZDC participants

    //V0 stuff
    TH1D *fHistV0multiplicity;     //! V0 mult
    TH1D *fHistV0Amultiplicity;     //! V0 A mult
    TH1D *fHistV0Cmultiplicity;     //! V0 C mult
    TH2D *fHist2DV0ACmultiplicity;     //! V0 A-C mult
    //TH1D *fHistV0spectra;     //! V0 particle masses
    TH2D *fHist2DTracksAcceptedVsV0multiplicity;     //! V0 mult Vs N tracks Accepted

    TH1D *fHistV0AmultiplicityRing[4];     //! V0 A mult in rings
    TH1D *fHistV0CmultiplicityRing[4];     //! V0 C mult in rings
    TH2D *fHist2DV0ACmultiplicityRing[4];     //! V0 A-C mult in rings
    TH2D *fHist2DTracksAcceptedVsV0AmultiplicityRing[4];     //! V0A mult Rings Vs N tracks Accepted
    TH2D *fHist2DTracksAcceptedVsV0CmultiplicityRing[4];     //! V0C mult Rings Vs N tracks Accepted

    TH1D *fHistV0cells         ;     //! V0 cells
    TH1D *fHistV0Acells        ;     //! V0 A cells
    TH1D *fHistV0Ccells        ;     //! V0 C cells
    TH2D *fHist2DV0ACcells     ;     //! V0 A-C cells

    Int_t fThresholdOnV0mult; //min V0AC mult to analyse this event (default is 0)


    Float_t fMinCentralityClass;    // min bound on centrality percentile
    Float_t fMaxCentralityClass;    // max bound on centrality percentile


    Bool_t fIsIonsAnalysis; //Ions analysis flag
    Bool_t fEtInsteadOfPt; //pass the Et instead of Pt to LRC processors

    Int_t fTmpCounter; //! TMP

    const AliPIDResponse *fPIDResponse;     //! PID response object
    AliPIDCombined       *fPIDCombined;     //! combined PID object
    TH1F *fPriors[AliPID::kSPECIES];           //! priors
    TH2D *fPriorsUsed[AliPID::kSPECIES];       //! priors used
    TH2D *fProbTPCTOF[AliPID::kSPECIES];       //! combined probabilities vs mom TPC-TOF
    TH2D *fProbAllDets[AliPID::kSPECIES];       //! combined probabilities ALL dets vs mom

    TH1D *fHistPidMaxProbability;  //!hist of max probabilities for arrays PID species
    TH1D *fHistPidPureMaxProbability;  //!hist of max probabilities for arrays PID species (when detId is TPC+TOF)

    char fStrPIDforFwd[20];		//PID name for FWD win
    char fStrPIDforBwd[20];		//PID name for BWD win
    Bool_t fPIDsensingFlag; 		//flag that we sense PID in processors


    int fMultForZDCstudy[5]; //! threshold multiplicities for ZDC study

    //artificial inefficiency (27.09.12)
    Double_t fArtificialInefficiency;	// inefficiency by hand in [0,1], default is 0
    TH2D *fHistNumberOfDroppedByHandTracks;   //! Number of tracks which were dropped by hand vs N of accepted tracks
    TRandom3 *fRand; //random generator for some uses

    //phi artificial gaps
    Double_t fPhiArtificialGapBegin;    // inefficiency in phi - gap position edge begins
    Double_t fPhiArtificialGapEnd;      // inefficiency in phi - gap position edge ends

    //flags for inclusion of detectors info:
    Bool_t fFlagWatchZDC;   //study ZDC issues
    Bool_t fFlagWatchV0;    //study V0 issues
    Bool_t fFlagWatchFMD;   //study FMD issues


    TStopwatch *fAnalysisTimer;

    //test MC particles
//    TH1D *fHistMCvertexRdeltaFromParent;  //!MC R hist
//    TH1F *fHistMCparentsStat;  //! MC parent ratios for different partile classes
//    TH1F *fHistMCparentsEta[kNumberOfParentParticleClassesInMC];  //! MC parents eta distributions for different particle classes
//    TH1F *fHistMCchildsEta[kNumberOfParentParticleClassesInMC];  //! MC childs eta distributions for different partile classes
//    TH1F *fHistMCdeltaEtaChildParent[kNumberOfParentParticleClassesInMC];  //! MC delta eta b/n parent and child
//    TH1F *fHistMCdeltaPhiChildParent[kNumberOfParentParticleClassesInMC];  //! MC delta phi b/n parent and child
//    TH2D *fHist2DMCchildrenPhiChildParent[kNumberOfParentParticleClassesInMC];  //! MC delta b/n parent and child in eta-phi

//    TH1F *fHistMCNumberOfChildren[kNumberOfParentParticleClassesInMC]; //! Number of children for fathers in MC (for different father classes)
//    TH1D *fHistMCchildrenEtaDeviationsFromAverage[kNumberOfParentParticleClassesInMC];  //! MC delta b/n av. eta and each child's eta for each father (for different father classes)
//    TH1D *fHistMCchildrenPhiDeviationsFromAverage[kNumberOfParentParticleClassesInMC];  //! MC delta b/n av. phi and each child's phi for each father (for different father classes)
//    TH2D *fHist2DMCchildrenPhiDeviationsFromAverage[kNumberOfParentParticleClassesInMC];  //! MC delta b/n av. eta-phi and each child's eta and phi for each father (for different father classes)

//    TString fStrMCparticleClassForFillingLRC; // name of particle class for LRC filling (default is All)
//    Double_t fEtaMCanalysisCutMin; // spec MC analysis: cut on eta
//    Double_t fEtaMCanalysisCutMax; // spec MC analysis: cut on eta

//    TH1F *fHistMCparentDeepness;  //! MC deepness of parent tree from "physical primary" children
//    TH1F *fHistMCparentsInitialStat;  //! MC initial papa particle class distr

//    TH1F *fHistMCEtaInitialQuark;  //! MC eta distr of "initial" quarks
//    TH1F *fHistMCEtaInitialGluon;  //! MC eta distr of "initial" gluons
//    TH1F *fHistMCEtaInitialProton;  //! MC eta distr of "initial" protons

//    TH1F *fHistMCnumberInitialQuarksInEvent;  //! MC initial quarks number distr
//    TH1F *fHistMCnumberInitialGluonsInEvent;  //! MC initial gluons number distr
//    TH1F *fHistMCnumberInitialProtonInEvent;  //! MC initial proton number distr (check that there are 2)

//    TH1F *fHistMCnumberChildrenFromInitialQuarksInEvent;  //! MC children number from initial quarks distr
//    TH1F *fHistMCnumberChildrenFromInitialGluonsInEvent;  //! MC children number from initial gluons distr
//    TH1F *fHistMCnumberChildrenFromInitialProtonInEvent;  //! MC children number from initial proton distr (check that there are 2)


    // 4.01.2012: MyTree stuff
    //    AliSimpleEvent *fSimpleEvent;   // instance of simple event to be filled in analysis loop
    //    Int_t fNsimpleEvents;
    //    TTree *fEventTree;              //! event tree to write into output file
    //    Bool_t fSetIncludeEventTreeInOutput;    // flag to use event tree or not

    ClassDef(AliAnalysisTaskLRC, 10 );
};

#endif
