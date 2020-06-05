/// \class AliAnalysisTaskMKBase
/// \brief Base class providing functionality for derived AnalaysisTasks
///
/// This class does some basic event analysis and provides a framework and
/// functionality to derived analysis tasks. The idea is not duplicate code
/// which is needed by many tasks and keep the derived tasks very simple and
/// clean.
///
/// In addition, this class provides (static) functionality to easily create and fill THnSparseF.
///
/// A task derived from this one should override the following functions:
/// AddOutput()        - define output histograms
/// IsEventSelected()  - define events selection (return kTRUE if events hould be selected)
/// AnaEvent()         - called for selected events (MC and DATA)
/// AnaEventMC()       - called for selected MC events only
/// AnaEventDATA()     - called for selected DATA events only
///
/// For single particle/track loops the above functions can make use of
/// LoopOverAllParticles(Int_t flag) and LoopOverAllTracks(Int_t flag)
/// where the flag can be used to distinguish multiple loops
/// default is flag=0
///
/// Inside the loops the following functions are called
/// (these functions should be overwritten by the derived task)
/// AnaParticleMC(flag) for LoopOverAllParticles (only for MC events)
/// AnaTrack(flag)      for LoopOverAllTracks
/// AnaTrackMC(flag)    for LoopOverAllTracks    (only for MC events)
/// AnaEventDATA(flag)  for LoopOverAllTracks    (only for DATA events)
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Mar 8, 2019

#ifndef AliAnalysisTaskMKBase_H
#define AliAnalysisTaskMKBase_H

#include "THnSparse.h"
#include "THn.h"
#include "AliAnalysisTaskSE.h"
#include "AlidNdPtTools.h"
#include "AliEventCuts.h"

#include <vector>

class AliExternalTrackParam;
class AliESDtrackCuts;
class AliInputEventHandler;
class AliAnalysisManager;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class AliESDtrack;
class AliEventCuts;
class AliMCParticle;
class AliMultSelection;
class AliGenPythiaEventHeader;
class AliGenDPMjetEventHeader;

class AliAnalysisTaskMKBase : public AliAnalysisTaskSE
{
  public:
    AliAnalysisTaskMKBase();
    AliAnalysisTaskMKBase(const char *name);
    virtual                 ~AliAnalysisTaskMKBase();

    virtual void           UserCreateOutputObjects();
    virtual void           UserExec(Option_t* option);
    virtual void           Terminate(Option_t* option);

    virtual void           SetESDtrackCutsM(AliESDtrackCuts* cut) { fESDtrackCutsM = cut; }
    virtual Bool_t         SetESDtrackCuts(Int_t i, AliESDtrackCuts* cut) { if (i<10) fESDtrackCuts[i] = cut; return (i<10); }
    virtual void           SetTriggerMaskRequired(UInt_t trigger) { fTriggerMaskRequired = trigger; }
    virtual void           SetTriggerMaskRejected(UInt_t trigger) { fTriggerMaskRejected = trigger; }
    virtual void           SetUseEventCuts(Bool_t use = kTRUE) { fUseEventCuts = use; }

    AliESDtrackCuts*       GetESDtrackCutsM() { return fESDtrackCutsM; }
    AliESDtrackCuts*       GetESDtrackCuts(Int_t i) { return (i < 10) ? fESDtrackCuts[i] : 0; }

    static Long64_t        FillHist(THnSparseF* s, Double_t x1, Double_t x2=0, Double_t x3=0, Double_t x4=0, Double_t x5=0, Double_t x6=0, Double_t x7 =0, Double_t x8 =0, Double_t x9 =0, Double_t x10 =0, Double_t x11 =0, Double_t x12 =0) { return AlidNdPtTools::FillHist(s, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12); }
    static Long64_t        FillHist(THnF* s, Double_t x1, Double_t x2=0, Double_t x3=0, Double_t x4=0, Double_t x5=0, Double_t x6=0, Double_t x7 =0, Double_t x8 =0, Double_t x9 =0, Double_t x10 =0, Double_t x11 =0, Double_t x12 =0) { return AlidNdPtTools::FillHist(s, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12); }
    static Long64_t        FillHistWeighted(THnSparseF* s, std::vector<double> const& val, double weight) {return AlidNdPtTools::FillHistWeighted(s, val, weight);}
    static Int_t           AddAxis(const char* label, Int_t nbins, Double_t xmin, Double_t xmax, const char* option = 0) { return AlidNdPtTools::AddAxis(label, nbins, xmin, xmax, option); }
    static Int_t           AddAxis(const char* label, const char* title, Int_t nbins, Double_t xmin, Double_t xmax, const char* option = 0) { return AlidNdPtTools::AddAxis(label, title, nbins, xmin, xmax, option); }
    static Int_t           AddAxis(const char* label, Int_t nbins, Double_t* xbins, const char* option = 0) { return AlidNdPtTools::AddAxis(label, nbins, xbins, option); }
    static Int_t           AddAxis(const char* label, const char* title, Int_t nbins, Double_t* xbins, const char* option = 0) { return AlidNdPtTools::AddAxis(label, title, nbins, xbins, option = 0); }
    static Int_t           AddAxis(const char* label, const char* title, const char* option) { return AlidNdPtTools::AddAxis(label, title, option); }
    static Int_t           AddAxis(const char* label, const char* option) { return AlidNdPtTools::AddAxis(label, option); }
    static Int_t           AddAxis(const char* option) { return AlidNdPtTools::AddAxis(option); }
    static THnSparseF*     CreateHist(const char* name) { return AlidNdPtTools::CreateHist(name); }
    static void            ResetHist() { AlidNdPtTools::ResetHist(); }
    static TH1D*           CreateLogHist(const char* name, const char* title) { return AlidNdPtTools::CreateLogHist(name, title); }
    static TH1D*           CreateLogHist(const char* name) { return AlidNdPtTools::CreateLogHist(name); }

    static AliAnalysisTaskMKBase* AddTaskMKBase(const char* name = "TaskMKBase", const char* outfile = 0);
    enum CentralityEstimator {kV0M=0, kCL0, kCL1, kV0Mplus05, kV0Mplus10, kV0Mminus05, kV0Mminus10, kSPDClustersCorr, kSPDTracklets};
    void                   SetCentralityEstimator(CentralityEstimator const _est) {fCentralityEstimator=_est;}
    CentralityEstimator    GetCentralityEstimator() const {return fCentralityEstimator;}

    // setters for switches
    void                   SetUseBaseOutput(Bool_t _use=kTRUE){fUseBaseOutput=_use;}
    void                   SetNeedEventVertex(Bool_t _use=kTRUE){fNeedEventVertex=_use;}
    void                   SetNeedEventCent(Bool_t _use=kTRUE){fNeedEventCent=_use;}
    void                   SetNeedEventMult(Bool_t _use=kTRUE){fNeedEventMult=_use;}
    void                   SetNeedEventVZERO(Bool_t _use=kTRUE){fNeedEventVZERO=_use;}
    void                   SetNeedTrackIP(Bool_t _use=kTRUE){fNeedTrackIP=_use;}
    void                   SetNeedTrackTPC(Bool_t _use=kTRUE){fNeedTrackTPC=_use;}
    void                   SetNeedTrackPID(Bool_t _use=kTRUE){fNeedTrackPID=_use;}
  protected:

    virtual void          Log(const char* name) {if(fUseBaseOutput) Log(fLogHist,name); }
    virtual void          Err(const char* name) {if(fUseBaseOutput) Log(fLogErr,name); }
    virtual void          LogEvent(const char* name) {if(fUseBaseOutput) Log(fLogEvent,name); }

    virtual void          Log(const char* name, Int_t n)    {if(fUseBaseOutput) Log(fLogHist,name,n); }
    virtual void          Log(const char* name, Double_t n) {if(fUseBaseOutput) Log(fLogHist,name,n); }

    virtual void          Log(TH1D* h, const char* name) { if (h) h->Fill(name,1); }
    virtual void          Log(TH1D* h, const char* name, Int_t n)    { TString s(name); s+=n; Log(h,s.Data()); }
    virtual void          Log(TH1D* h, const char* name, Double_t n) { TString s(name); s+=n; Log(h,s.Data()); }

    virtual Bool_t        ReadEvent(); // read the event info
    virtual Bool_t        ReadMCEvent(); // read the mc event info
    virtual void          FillDefaultHistograms(Int_t step); // fill default histograms of the base class, step=0(before),1(after) user evt selection

    virtual void          AddOutput() {}; // add histograms to fOutputList

    virtual void          AnaEvent() {}; //called for every event, both data and mc (to be implemented in derived class)
    virtual void          AnaEventMC() {}; //called for every mc event (to be implemented in derived class)
    virtual void          AnaEventDATA() {}; //called only every data event (to be implemented in derived class)
    virtual void          AnaTrack(Int_t flag = 0) {}; //called for every track, both data and mc (to be implemented in derived class)
    virtual void          AnaTrackMC(Int_t flag = 0) {}; //called for every track, both data and mc (to be implemented in derived class)
    virtual void          AnaTrackDATA(Int_t flag = 0) {}; //called for every track with MC info (to be implemented in derived class)
    virtual void          AnaParticleMC(Int_t flag = 0) {}; //called for every MC Particle (to be implemented in derived class)
    virtual Bool_t        IsEventSelected() { return kTRUE; }; //user defined event selection, default is all events are accepted

    virtual void          BaseAnaTrack(Int_t flag = 0);      // wraps AnaTracK, to be used for mult counting
    virtual void          BaseAnaParticleMC(Int_t flag = 0); // wraps AnaParticleMC, to be used for mult counting
    virtual void          BaseAddOutput();  //used to

    //
    virtual Bool_t          InitEvent();   // loads event-related properties
    virtual Bool_t          InitEventMult();   //initialize multiplicity specific variables, requires corresponding task
    virtual Bool_t          InitEventCent();   //initialize (old) centrality specific variables, requires corresponding task
    virtual Bool_t          InitEventVZERO(); // initalize the VZERO information (to be completed)
    virtual Bool_t          InitEventVertex(); // initalize the event vertex information

    virtual Bool_t          InitMCEvent(); //load mc event-related properties
    virtual Bool_t          InitMCEventType(); // load information about mc event type sd,nd,dd etc. works for dpmjet and pythia only

    virtual Bool_t          InitTrack();   //initializes track related quantities
    virtual Bool_t          InitTrackQA();   //initializes track related quantities needed for QA
    virtual Bool_t          InitTrackCuts(); //check all track cuts and set corresponding variables
    virtual Bool_t          InitTrackIP();  //initialize inner params
    virtual Bool_t          InitTrackTPC();  //initialize inner params tpc
    virtual Bool_t          InitTrackPID(); //initilaize PID related quantities

    virtual Bool_t          InitMCTrack();

    virtual Bool_t          InitMCParticle();

    virtual void            LoopOverAllTracks(Int_t flag = 0);    // loops over all tracks in the event, calls AnaTrack(), AnaTrackMC() and AnaTrackDATA() for each track
    virtual void            LoopOverAllParticles(Int_t flag = 0); // loops over all MC particles in the event, calls AnaParticleMC() for each particle

    //         virtual void            FillTrigInfo(TH1D* h);   // fill the trigger histogram
    virtual void            FillTrigHist(TH1D* h);   // fill the trigger histogram
    //         virtual void            FillRunHist(TH1D* h);   // fill the trigger histogram

    virtual void            InitEventChecks();       // do the event checks

    // event related properties
    AliAnalysisManager*             fAnalysisManager;           //!<!  analysis manager                                                  --ReadEvent()
    AliInputEventHandler*           fInputEventHandler;         //!<!  input event handler                                               --ReadEvent()
    UInt_t                          fEventSelected;             //!<!  AliInputEventHandler::IsEventSelected()                           --ReadEvent()
    AliVEvent*                      fEvent;                     //!<!  input event                                                       --ReadEvent()
    AliESDEvent*                    fESD;                       //!<!  input ESD event                                                   --ReadEvent()
    AliAODEvent*                    fAOD;                       //!<!  input AOD event                                                   --ReadEvent()
    AliMCEvent*                     fMC;                        //!<!  MC event                                                          --ReadMCEvent()
    AliESDVZERO*                    fVZERO;                     //!<!  VZERO information                                                 --InitEventVZERO()
    AliStack*                       fMCStack;                   //!<!  MC stack                                                          --ReadMCEvent()
    AliHeader*                      fMCHeader;                  //!<!  MC header                                                         --ReadMCEvent()
    AliGenEventHeader*              fMCGenHeader;               //!<!  MC gen event header                                               --ReadMCEvent()
    AliGenPythiaEventHeader*        fMCGenHeaderPythia;         //!<!  Pythia event header                                               --InitMCEventType()
    AliGenDPMjetEventHeader*        fMCGenHeaderDPMjet;         //!<!  DPMjet gen event header                                           --InitMCEventType()
    AlidNdPtTools::EventType        fMCEventType;               //!<!  event type (sd,nd,cd,dd,elastic)                                  --InitMCEventType()
    Int_t                           fMCProcessTypeFlag;         //!<!  event type flag form mc (depende on generator)                    --InitMCEventType()
    AliMultSelection*               fMultSelection;             //!<!  AliMultSelection                                                  --InitEventMult()
    AliCentrality*                  fCentrality;                //!<!  AliCentrality                                                     --InitEventCent()
    Int_t                           fNTracksESD;                //!<!  number of esd tracks in the event                                 --InitEvent()
    Int_t                           fNTracksAcc;                //!<!  number of accepted trackswith trackcuts                           --InitEvent()
    Bool_t                          fIsMC;                      //!<!  do we have an MC event?                                           --ReadMCEvent()
    Int_t                           fRunNumber;                 //!<!  run number                                                             --InitEvent()
    TString                         fRunNumberString;           //!<!  run number as string                                              --InitEvent()
    UInt_t                          fTimeStamp;                 //!<!  event time stamp                                                  --InitEvent()
    Int_t                           fEventNumberInFile;         //!<!  event number in file                                              --InitEvent()
    TString                         fFiredTriggerClasses;       //!<!  all trigger classes as string                                     --ReadEvent()
    UInt_t                          fEventSpecie;               //!<!  event specie                                                      --InitEvent()
    Double_t                        fOldCentPercentileV0M;      //!<!  centrality percentile from old framework                          --InitEventCent()
    Double_t                        fMultPercentileV0M;         //!<!  centrality/multiplicity percentile from new framework             --InitEventMult()
    Bool_t                          fIsAcceptedAliEventCuts;    //!<!  accepted by AliEventCuts?                                         --InitEventChecks()

    const AliESDVertex*             fVtx;                       //!<!  best available vertex                                             --InitEventVertex()
    Double_t                        fXv;                        //!<!  x vertex position                                                 --InitEventVertex()
    Double_t                        fYv;                        //!<!  y vertex position                                                 --InitEventVertex()
    Double_t                        fZv;                        //!<!  z vertex position                                                 --InitEventVertex()
    Double_t                        fXvRes;                     //!<!  x vertex resolution                                               --InitEventVertex()
    Double_t                        fYvRes;                     //!<!  y vertex resolution                                               --InitEventVertex()
    Double_t                        fZvRes;                     //!<!  z vertex resolution                                               --InitEventVertex()
    Int_t                           fVtxNContrib;               //!<!  vertex N Contributers                                             --InitEventVertex()
    Bool_t                          fVtxStatus;                 //!<!  vertex status                                                     --InitEventVertex()
    Double_t                        fVtxDispersion;             //!<!  vertex dipserion                                                  --InitEventVertex()
    Bool_t                          fUsedVtxTRK;                //!<!  is vertex tracks                                                  --InitEventVertex()
    Bool_t                          fUsedVtxSPD;                //!<!  is vertex spd                                                     --InitEventVertex()
    Bool_t                          fUsedVtxTPC;                //!<!  is vertex tpc                                                     --InitEventVertex()

    const AliESDVertex*             fVtxTRK;                    //!<!  track vertex                                                      --InitEventVertex()
    Double_t                        fXvTRK;                     //!<!  x vertex position                                                 --InitEventVertex()
    Double_t                        fYvTRK;                     //!<!  y vertex position                                                 --InitEventVertex()
    Double_t                        fZvTRK;                     //!<!  z vertex position                                                 --InitEventVertex()
    Double_t                        fXvResTRK;                  //!<!  x vertex resolution                                               --InitEventVertex()
    Double_t                        fYvResTRK;                  //!<!  y vertex resolution                                               --InitEventVertex()
    Double_t                        fZvResTRK;                  //!<!  z vertex resolution                                               --InitEventVertex()
    Int_t                           fVtxNContribTRK;            //!<!  z vertex N Contributers                                           --InitEventVertex()
    Bool_t                          fVtxStatusTRK;              //!<!  z vertex status                                                   --InitEventVertex()
    Double_t                        fVtxDispersionTRK;          //!<!  z vertex dipserion                                                --InitEventVertex()

    const AliESDVertex*             fVtxSPD;                    //!<!  SPD vertex                                                        --InitEventVertex()
    Double_t                        fXvSPD;                     //!<!  x vertex position                                                 --InitEventVertex()
    Double_t                        fYvSPD;                     //!<!  y vertex position                                                 --InitEventVertex()
    Double_t                        fZvSPD;                     //!<!  z vertex position                                                 --InitEventVertex()
    Double_t                        fXvResSPD;                  //!<!  x vertex resolution                                               --InitEventVertex()
    Double_t                        fYvResSPD;                  //!<!  y vertex resolution                                               --InitEventVertex()
    Double_t                        fZvResSPD;                  //!<!  z vertex resolution                                               --InitEventVertex()
    Int_t                           fVtxNContribSPD;            //!<!  z vertex N Contributers                                           --InitEventVertex()
    Bool_t                          fVtxStatusSPD;              //!<!  z vertex status                                                   --InitEventVertex()
    Double_t                        fVtxDispersionSPD;          //!<!  z vertex dipserion                                                --InitEventVertex()

    const AliESDVertex*             fVtxTPC;                    //!<!  TPC vertex                                                        --InitEventVertex()
    Double_t                        fXvTPC;                     //!<!  x vertex position                                                 --InitEventVertex()
    Double_t                        fYvTPC;                     //!<!  y vertex position                                                 --InitEventVertex()
    Double_t                        fZvTPC;                     //!<!  z vertex position                                                 --InitEventVertex()
    Double_t                        fXvResTPC;                  //!<!  x vertex resolution                                               --InitEventVertex()
    Double_t                        fYvResTPC;                  //!<!  y vertex resolution                                               --InitEventVertex()
    Double_t                        fZvResTPC;                  //!<!  z vertex resolution                                               --InitEventVertex()
    Int_t                           fVtxNContribTPC;            //!<!  z vertex N Contributers                                           --InitEventVertex()
    Bool_t                          fVtxStatusTPC;              //!<!  z vertex status                                                   --InitEventVertex()
    Double_t                        fVtxDispersionTPC;          //!<!  z vertex dipserion                                                --InitEventVertex()

    Double_t                        fMCxv;                      //!<!  mc truth x vertex position                                        --InitMCEvent()
    Double_t                        fMCyv;                      //!<!  mc truth y vertex position                                        --InitMCEvent()
    Double_t                        fMCzv;                      //!<!  mc truth z vertex position                                        --InitMCEvent()
    Int_t                           fMultMB;                    //!<!  MinBias Multiplicity (no of contributers to vertex)               --InitEventVertex()
    Double_t                        fMultV0A;                   //!<!  v0a amplitude                                                     --InitEventVZERO()
    Double_t                        fMultV0C;                   //!<!  v0c amplitude                                                     --InitEventVZERO()
    Double_t                        fMultV0M;                   //!<!  v0a + v0c ampliude                                                --InitEventVZERO()
    Double_t                        fMultV0MmultSelection;      //!<!  v0a + v0c ampliude as used in AliMultSelection                    --
    Double_t                        fMCb;                       //!<!  impact parameter in MC                                            --
    Int_t                           fMCnPrimPtCut;              //!<!  ch. prim. particles according to mario def (eta<0.8, pt>150 MeV)  --InitMCEvent()
    Int_t                           fMCnPrim10;                 //!<!  ch. prim. particles according to mc in eta<1.0                    --InitMCEvent()
    Int_t                           fMCnPrim08;                 //!<!  ch. prim. particles according to mc in eta<0.8                    --InitMCEvent()
    Int_t                           fMCnPrim05;                 //!<!  ch. prim. particles according to mc in eta<0.5                    --InitMCEvent()
    Int_t                           fMCnPrimV0M;                //!<!  ch. prim. particles in the v0 acceptance                          --InitMCEvent()
    Int_t                           fMCnTracks;                 //!<!  number of "tracks" i.e. particles in MCevent                      --InitMCEvent()
    Int_t                           fMCnPrim;                   //!<!  number of charged physical primaries in MCEvent                   --InitMCEvent()
    Bool_t                          fIsTrigger;                 //!<!  is event triggered?                                               --
    Bool_t                          fHasVertex;                 //!<!  has the event a vertex?                                           --InitEventVertex()
    Bool_t                          fIsIncompleteDAQ;           //!<!  incomplete daq event                                              --InitEventChecks()
    Bool_t                          fIsSPDClusterVsTrackletBG;  //!<!  spd cluster vs tracklet bg                                        --InitEventChecks()
    Bool_t                          fIsFirstEventInChunk;       //!<!  first event in chunk                                              --InitEventChecks()
    Bool_t                          fIsPileUpMV;                //!<!  is pileup mv                                                      --InitEventChecks()
    Bool_t                          fIsOutOfBunchPileUp;        //!<!  is out of bunch pileup                                            --InitEventChecks()
    Bool_t                          fIsPileUpEvent;             //!<!  is pileup event                                                   --InitEventChecks()
    Bool_t                          fIsPileUpSPD;               //!<!  is pileup from spd                                                --InitEventChecks()
    Bool_t                          fIsVertexRejected2013pA;    //!<!  vertex 2013 pA rejected                                           --InitEventChecks()
    Bool_t                          fIsPileupFromSPD508;        //!<!  is pileup from spd(5,0,8)                                         --InitEventChecks()
    Bool_t                          fIsEventAccepted;           //!<!  is Event selected according to user definition by user            --

    // track related properties
    AliESDtrack*                    fESDTrack;                  //!<!  current esd track                                                 --
    Double_t                        fPt;                        //!<!  track pT                                                          --InitTrack()
    Double_t                        fP;                         //!<!  track p                                                           --InitTrack()
    Double_t                        fEta;                       //!<!  track Eta                                                         --InitTrack()
    Double_t                        fPhi;                       //!<!  track Phi                                                         --InitTrack()
    Float_t                         fDCA[2];                    //!<!  impact parameter (DCA)                                            --InitTrack()
    Float_t                         fDCACov[3];                 //!<!  impat parameter (DCA) covariance                                  --InitTrack()
    Double_t                        fDCAr;                      //!<!  impact parameter (DCA) in xy-direction                            --InitTrack()
    Double_t                        fDCAz;                      //!<!  impact parameter (DCA) in z-direction                             --InitTrack()
    Double_t                        fSigma1Pt2;                 //!<!  sigma(1/pT)**2                                                    --InitTrack()
    Double_t                        fSigma1Pt;                  //!<!  sigma(1/pT)                                                       --InitTrack()
    Double_t                        fSigned1Pt;                 //!<!  signed 1/pT                                                       --InitTrack()
    Double_t                        f1Pt;                       //!<!  1/pT                                                              --InitTrack()
    Short_t                         fChargeSign;                //!<!  sign of the track charge                                          --InitTrack()
    UShort_t                        fTPCSignalN;                //!<!  number of clusters for PID                            --InitTrack()

//BEGIN NEW
    Double_t                        fX;                         //!<! x at dca (radial distance to vertex)
    Double_t                        fY;                         //!<! local Y-coordinate of track at dca  (cm)
    Double_t                        fZ;                         //!<! local Z-coordinate of track at dca  (cm)
    Double_t                        fAlpha;                     //!<! local to global angle

    Double_t                        fSnp;                       //!<! local sine of the track momentum azimuthal angle
    Double_t                        fTgl;                       //!<! tangent of the track momentum dip angle
    ULong_t                         fFlags;                     //!<! flags assigned to the track
  
    Double_t                        fITSFoundClusters;          //!<! found clusters ITS
    Double_t                        fITSChi2PerCluster;         //!<! chi2 per cluster ITS
    UChar_t                         fITSClusterMap;             //!<! hitmap ITS

    Double_t                        fTPCFindableClusters;       //!<! findable clusters TPC
    Double_t                        fTPCFoundClusters;          //!<! found clusters TPC
    Double_t                        fTPCSharedClusters;         //!<! shared clusters TPC
    Double_t                        fTPCFractionSharedClusters; //!<! fraction of shared clusters TPC
    Double_t                        fTPCCrossedRows;            //!<! crossed rows in TPC
    Double_t                        fTPCCrossedRowsOverFindableClusters;            //!<! crossed rows over findable clusters in TPC
    Double_t                        fTPCChi2PerCluster;            //!<! chi2 per cluster TPC
//END NEW

    AliMCParticle*                  fMCParticle;                //!<!  mc particle                                                       --
    Int_t                           fMCLabel;                   //!<!  mc label                                                          --
    Double_t                        fMCPt;                      //!<!  mc pt                                                             --InitMCParticle()
    Double_t                        fMCEta;                     //!<!  mc eta                                                            --InitMCParticle()
    Double_t                        fMCPhi;                     //!<!  mc phi                                                            --InitMCParticle()
    Bool_t                          fMCisPrim;                  //!<!  is physical primary?                                              --InitMCParticle()
    Bool_t                          fMCisSec;                   //!<!  is secondary?                                                     --InitMCParticle()
    Bool_t                          fMCisSecDecay;              //!<!  is secondary from decay?                                          --InitMCParticle()
    Bool_t                          fMCisSecMat;                //!<!  is secondary from material?                                       --InitMCParticle()
    Int_t                           fMCPrimSec;                 //!<!  status of mc track: 0=prim, 1=decay 2=material                    --InitMCParticle()
    AlidNdPtTools::ParticleType     fMCParticleType;            //!<!  which particle is it                                              --InitMCParticle()
    AlidNdPtTools::ProductionType   fMCProdcutionType;          //!<!  production mechanism (prim,material,decay)                        --InitMCParticle()
    Int_t                           fMCPDGCode;                 //!<!  PDG code                                                          --InitMCParticle()
    Short_t                         fMCCharge;                  //!<!  charge in units of 1/3e                                           --InitMCParticle()
    Double_t                        fMCQ;                       //!<!  charge in units of e                                              --InitMCParticle()
    Bool_t                          fMCIsCharged;               //!<!  charged particle                                                  --InitMCParticle()
    Short_t                         fMCChargeSign;              //!<!  Sign of the charge                                                --InitMCParticle()

    const AliExternalTrackParam*    fInnerP;                    //!<!  innerparams                                                       --InitTrackIP()
    const AliExternalTrackParam*    fTPCinnerP;                 //!<!  TPC inner params                                                  --InitTrackTPC()
    Double_t                        fPtInner;                   //!<!  inner param pt                                                    --InitTrackIP()
    Double_t                        fEtaInner;                  //!<!  inner param eta                                                   --InitTrackIP()
    Double_t                        fPhiInner;                  //!<!  inner param phi                                                   --InitTrackIP()
    Double_t                        fZInner;                    //!<!  inner param Z                                                   --InitTrackIP()
    Double_t                        fPtInnerTPC;                //!<!  TPC inner pt                                                      --InitTrackTPC()
    Double_t                        fEtaInnerTPC;               //!<!  TPC inner eta                                                     --InitTrackTPC()
    Double_t                        fPhiInnerTPC;               //!<!  TPC inner phi                                                     --InitTrackTPC()
    Double_t                        fZInnerTPC;                 //!<!  TPC inner parameter of track                              --InitTrackTPC
    Float_t                         fDCATPC[2];                 //!<!  TPC impact parameter (DCA)                                        --InitTrackTPC()
    Float_t                         fDCACovTPC[3];              //!<!  TPC impat parameter (DCA) covariance                              --InitTrackTPC()
    Double_t                        fDCArTPC;                   //!<!  TPC impact parameter (DCA) in xy-direction                        --InitTrackTPC()
    Double_t                        fDCAzTPC;                   //!<!  TPC impact parameter (DCA) in z-direction                         --InitTrackTPC()

    AliEventCuts                    fEventCuts;                 //!<!  event cuts                                                          --
    Bool_t                          fUseEventCuts;              ///<  use event cuts?                                                     --
    AliESDtrackCuts*                fESDtrackCutsM;             //-> trackcuts used for mult estimate                                   --
    Bool_t                          fAcceptTrackM;              ///<  is track accepted by fESDtrackCutsM                                --InitTrackCuts()
    AliESDtrackCuts*                fESDtrackCuts[10];          //-> several track cuts that can be used in the analysis                --
    Bool_t                          fAcceptTrack[10];           ///<  is track accepted by fESDtrackCuts[10]                             --InitTrackCuts()
    TString                         fMultEstimator;             ///<  mult/cent estimator                                                  --
    TString                         fCentEstimator;             ///<  old cent estimator                                                   --
    UInt_t                          fTriggerMaskRequired;       ///<  only events with this trigger mask are accepted                      --
    UInt_t                          fTriggerMaskRejected;       ///<  reject events with this trigger mask                                 --

    Bool_t                          fInternalLoop;              ///<  used to flag an internal particle/track loop for AliAnalysisTaskMKBase

    // output list and control histograms
    TList*                          fOutputList;            //!<!  output list
    TH1D*                           fLogHist;               //!<!  generic log histogram use Log() to fill
    TH1D*                           fLogErr;                //!<!  histogram with errors that should not occur use Err() to fill
    TH1D*                           fLogEvent;              //!<!  event related logs use LogEvent() to fill
    TH1D*                           fRunHist;               //!<!  distribution of events in runs before event selection
    TH1D*                           fRunHistSelected;       //!<!  distribution of events in runs after event selection
    TH1D*                           fTrigInfo;              //!<!  all the trigger strings
    TH1D*                           fTrigInfoSelected;      //!<!  all the trigger strings of selected events
    TH1D*                           fTrigHist;              //!<!  AliVEvent trigger classes
    TH1D*                           fTrigHistSelected;      //!<!  AliVEvent trigger classes of selected events

    CentralityEstimator             fCentralityEstimator;   ///<

    // switches to enable or disable certain loops
    Bool_t                          fUseBaseOutput;         ///<
    // event level
    Bool_t                          fNeedEventVertex;       ///<
    Bool_t                          fNeedEventCent;         ///<
    Bool_t                          fNeedEventMult;         ///<
    Bool_t                          fNeedEventVZERO;        ///<
    // track level
    Bool_t                          fNeedTrackIP;           ///<
    Bool_t                          fNeedTrackTPC;          ///<
    Bool_t                          fNeedTrackPID;          ///<

  private:
    AliAnalysisTaskMKBase(const AliAnalysisTaskMKBase&); // not implemented
    AliAnalysisTaskMKBase& operator=(const AliAnalysisTaskMKBase&); // not implemented

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskMKBase, 7);
    /// \endcond

};

#endif
