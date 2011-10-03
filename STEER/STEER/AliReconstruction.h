#ifndef ALIRECONSTRUCTION_H
#define ALIRECONSTRUCTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for running the reconstruction                                      //
// Clusters and tracks are created for all detectors and all events by       //
// typing:                                                                   //
//                                                                           //
//   AliReconstruction rec;                                                  //
//   rec.Run();                                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TSelector.h>
#include <TString.h>
#include <TObjArray.h>

class AliReconstructor;
class AliRunLoader;
class AliRawReader;
class AliLoader;
class AliTracker;
class AliMagF;
class AliVertexer;
class AliTrackleter;
class AliESDVertex;
class AliESDEvent;
class AliESDfriend;
class AliESDpid;
class AliVertexerTracks;
class TFile;
class TTree;
class TList;
class AliQAManager; 
class TMap;
class AliRecoParam;
class AliDetectorRecoParam;
class AliRunInfo;
class AliGRPObject;
class THashTable;
class AliTriggerRunScalers;
class AliCTPTimeParams;
class AliAnalysisManager;
class AliRecoInputHandler;

#include "AliQAv1.h"
#include "AliEventInfo.h"
#include "AliRecoParam.h"

class AliReconstruction: public TSelector {
public:
  AliReconstruction(const char* gAliceFilename = "galice.root");
  virtual ~AliReconstruction();

  void           SetGAliceFile(const char* fileName);
  void           SetInput(const char* input);
  void           SetOutput(const char* output);

  void           SetEquipmentIdMap(const char *mapFile) {fEquipIdMap = mapFile;};
  void           SetEventRange(Int_t firstEvent = 0, Int_t lastEvent = -1) 
    {fFirstEvent = firstEvent; fLastEvent = lastEvent;};
  void           SetNumberOfEventsPerFile(UInt_t nEvents)
    {fNumberOfEventsPerFile = nEvents;};
  void           SetFractionFriends(Double32_t frac = 0.04)
    {fFractionFriends = frac;};
  void           SetOption(const char* detector, const char* option);
  void           SetRecoParam(const char* detector, AliDetectorRecoParam *par);

  void           SetRunLocalReconstruction(const char* detectors) {
    fRunLocalReconstruction = detectors;};
  void           SetRunTracking(const char* detectors) {
    fRunTracking = detectors;};
  void           SetFillESD(const char* detectors) {fFillESD = detectors;};
  void           SetRunReconstruction(const char* detectors) {
    SetRunLocalReconstruction(detectors); 
    SetRunTracking(detectors);
    SetFillESD(detectors);};
  void           SetUseTrackingErrorsForAlignment(const char* detectors) 
    {fUseTrackingErrorsForAlignment = detectors;};
  void           SetLoadAlignFromCDB(Bool_t load)  {fLoadAlignFromCDB = load;};
  void           SetLoadAlignData(const char* detectors) 
    {fLoadAlignData = detectors;};

  //*** Global reconstruction flag setters
  void SetRunMultFinder(Bool_t flag=kTRUE) {fRunMultFinder=flag;};
  void SetRunVertexFinder(Bool_t flag=kTRUE) {fRunVertexFinder=flag;};
  void SetRunVertexFinderTracks(Bool_t flag=kTRUE) {fRunVertexFinderTracks=flag;};
  void SetRunHLTTracking(Bool_t flag=kTRUE) {fRunHLTTracking=flag;};
  void SetRunV0Finder(Bool_t flag=kTRUE) {fRunV0Finder=flag;};
  void SetRunCascadeFinder(Bool_t flag=kTRUE) {fRunCascadeFinder=flag;};
  void SetStopOnError(Bool_t flag=kTRUE) {fStopOnError=flag;}
  void SetWriteAlignmentData(Bool_t flag=kTRUE){fWriteAlignmentData=flag;}
  void SetWriteESDfriend(Bool_t flag=kTRUE){fWriteESDfriend=flag;}
  void SetFillTriggerESD(Bool_t flag=kTRUE){fFillTriggerESD=flag;}
  void SetDiamondProfileSPD(AliESDVertex *dp) {fDiamondProfileSPD=dp;}
  void SetDiamondProfile(AliESDVertex *dp) {fDiamondProfile=dp;}
  void SetDiamondProfileTPC(AliESDVertex *dp) {fDiamondProfileTPC=dp;}
		   
  void SetCleanESD(Bool_t flag=kTRUE){fCleanESD=flag;}
  void SetUseHLTData(const char* detectors){fUseHLTData=detectors;}
  void SetV0DCAmax(Float_t d) {fV0DCAmax=d;}
  void SetV0CsPmin(Float_t d) {fV0CsPmin=d;}
  void SetDmax(Float_t d) {fDmax=d;}
  void SetZmax(Float_t z) {fZmax=z;}
  Float_t GetV0DCAmax() const {return fV0DCAmax;}
  Float_t GetV0CsPmin() const {return fV0CsPmin;}
  Float_t GetDmax() const {return fDmax;}
  Float_t GetZmax() const {return fZmax;}
  //
  Bool_t  IsRunMultFinder()   const {return fRunMultFinder;}
  
  // CDB storage activation
  void SetDefaultStorage(const char* uri);
  void SetSpecificStorage(const char* calibType, const char* uri);
  void SetFromCDBSnapshot(const char* snapshotFileName) {fFromCDBSnapshot = kTRUE; fSnapshotFileName=snapshotFileName;}
  void UnSetFromCDBSnapshot() {fFromCDBSnapshot = kFALSE;}

  Bool_t MisalignGeometry(const TString& detectors);

  void           SetAlignObjArray(TObjArray *array)
                   {fAlignObjArray = array;
		   fLoadAlignFromCDB = kFALSE;}

  virtual Int_t  Version() const {return 2;}
  virtual void   Begin(TTree*);
  virtual void   SlaveBegin(TTree*);
  virtual void   Init(TTree *tree);
  virtual Bool_t Process(Long64_t entry);
  virtual Bool_t ProcessEvent(Int_t iEvent);
  virtual void   SlaveTerminate();
  virtual void   Terminate();
  virtual Bool_t Run(const char* input = NULL);
  void           Abort(const char *method, EAbort what);
  virtual void	 SetOption(const char* option) {
    TSelector::SetOption(option);
  }

  // Trackers
  AliTracker* GetTracker(Int_t idx) const { return fTracker[idx]; }
  Bool_t      CreateTrackers(const TString& detectors);
  void        ImportRunLoader(AliRunLoader* rl) { fRunLoader = rl; }

  // Quality Assurance 
  void    SetQACycles(AliQAv1::DETECTORINDEX_t det, Int_t cycles) { fQACycles[det] = cycles ; }
  void    SetQAWriteExpert(AliQAv1::DETECTORINDEX_t det) { fQAWriteExpert[det] = kTRUE ; }
  Bool_t  SetRunQA(TString detAndAction="ALL:ALL") ; 
  void    SetRunGlobalQA(Bool_t flag=kTRUE){fRunGlobalQA = flag;}
  void    SetQARefDefaultStorage(const char* uri);
  void    InitQA();
  void    MergeQA(const char *fileName = NULL); 
  void    SetWriteQAExpert() { fWriteQAExpertData = kTRUE ; }
  Bool_t  IsRunQA() const { return fRunQA ; } 
  Bool_t  IsRunGlobalQA() const { return fRunGlobalQA ; }
  Bool_t  IsInTasks(AliQAv1::TASKINDEX_t tsk) { return fQATasks.Contains(Form("%d", tsk)) ; }
  
  // Plane Efficiency Evaluation
  void    SetRunPlaneEff(Bool_t flag=kFALSE)  {fRunPlaneEff = flag;}

  enum {
    kNDetectors = 15   // number of detectors
  };
  static Int_t   GetDetIndex(const char * detector);

  // Upgrade
 void SetUpgradeModule(const char* detectors)  {fUpgradeModule = detectors; MatchUpgradeDetector() ; }
 void MatchUpgradeDetector();

  // Connection to analysis manager. 
  // If a macro is set it should contain a function returning an analysis 
  // manager object. The macro is ignored in case the analysis manager is 
  // set via SetAnalysisManager
  void        SetAnalysisMacro(const char *macropath) {fAnalysisMacro = macropath;}
  const char *GetAnalysisMacro() const {return fAnalysisMacro.Data();}
  void        SetAnalysisManager(AliAnalysisManager *mgr) {fAnalysis = mgr;}
  AliAnalysisManager *GetAnalysisManager() const {return fAnalysis;}

protected:
  virtual Bool_t ProcessEvent(void* event);
  void           InitRun(const char* input);
  // Functions needed to select events for which we store the ESD friends
  Bool_t IsHighPt() const;
  Bool_t IsCosmicOrCalibSpecie() const;
  void WriteESDfriend();

private:
  AliReconstruction(const AliReconstruction& rec);
  AliReconstruction& operator = (const AliReconstruction& rec);

  void           InitRawReader(const char* input);
  void           InitCDB();
  Bool_t         InitGRP();
  void           SetCDBLock();
  Bool_t         SetRunNumberFromData();
  Bool_t         LoadCDB();
  Bool_t         LoadTriggerScalersCDB();
  Bool_t         LoadCTPTimeParamsCDB();
  Bool_t         ReadIntensityInfoCDB();
  Bool_t         RunLocalEventReconstruction(const TString& detectors);
  Bool_t         RunVertexFinder(AliESDEvent*& esd);
  Bool_t         RunHLTTracking(AliESDEvent*& esd);
  Bool_t         RunMuonTracking(AliESDEvent*& esd);
  Bool_t         RunSPDTrackleting(AliESDEvent*& esd);
  Bool_t         RunMultFinder(AliESDEvent*& esd);
  Bool_t         RunTracking(AliESDEvent*& esd, AliESDpid &PID);
  Bool_t         CleanESD(AliESDEvent *esd);
  Bool_t         FillESD(AliESDEvent*& esd, const TString& detectors);
  Bool_t         FillTriggerESD(AliESDEvent*& esd);
  Bool_t         FillTriggerScalers(AliESDEvent*& esd);
  Bool_t         FillRawEventHeaderESD(AliESDEvent*& esd);

  Bool_t         IsSelected(TString detName, TString& detectors) const;
  Bool_t         InitRunLoader();
  AliReconstructor* GetReconstructor(Int_t iDet);
  AliVertexer*   CreateVertexer();
  AliTrackleter* CreateMultFinder();
  void           CleanUp();

  Bool_t         ParseOutput();

  //==========================================//
  void           WriteAlignmentData(AliESDEvent* esd);

  void           FillRawDataErrorLog(Int_t iEvent, AliESDEvent* esd);

  //Quality Assurance
  //  void                 CheckQA() ;

  // Plane Efficiency evaluation
  Bool_t  FinishPlaneEff(); //ultimate tasks related to Plane Eff. evaluation 
  Bool_t  InitPlaneEff();   // initialize what is needed for Plane Eff. evaluation

  Bool_t               InitAliEVE();
  void                 RunAliEVE();

  Bool_t         InitRecoParams(); // init the array with the reconstruciton parameters
  Bool_t         GetEventInfo();   // fill the event info inside the event loop

  const char    *MatchDetectorList(const char *detectorList, UInt_t detectorMask);

  //*** Global reconstruction flags *******************
  Bool_t         fRunVertexFinder;    // run the vertex finder
  Bool_t         fRunVertexFinderTracks;    // run the vertex finder with tracks
  Bool_t         fRunHLTTracking;     // run the HLT tracking
  Bool_t         fRunMuonTracking;    // run the HLT tracking
  Bool_t         fRunV0Finder;        // run the ESD V0 finder
  Bool_t         fRunCascadeFinder;   // run the ESD cascade finder
  Bool_t         fRunMultFinder;      // run the trackleter for ITS clusters
  Bool_t         fStopOnError;        // stop or continue on errors
  Bool_t         fWriteAlignmentData; // write track space-points flag
  Bool_t         fWriteESDfriend;     // write ESD friend flag
  Bool_t         fFillTriggerESD;     // fill trigger info into ESD

  //*** Clean ESD flag and parameters *******************
  Bool_t         fCleanESD;      // clean ESD flag
  Float_t        fV0DCAmax;      // max. allowed DCA between V0 daugthers 
  Float_t        fV0CsPmin;      // min. allowed cosine of V0 pointing angle 
  Float_t        fDmax;          // max. allowed transverse impact parameter 
  Float_t        fZmax;          // max. allowed longitudinal impact parameter 

  TString        fRunLocalReconstruction; // run the local reconstruction for these detectors
  TString        fRunTracking;        // run the tracking for these detectors
  TString        fFillESD;            // fill ESD for these detectors
  TString        fLoadCDB;            // prefetch CDB entries and init reco-params for these detectors
  TString        fUseTrackingErrorsForAlignment; // for these detectors
  TString        fGAliceFileName;     // name of the galice file
  TString        fRawInput;           // name of input raw-data file or directory
  TString        fESDOutput;          //! name of output ESD file
  TString        fProofOutputFileName;//! name of output file in PROOF mode
  TString        fProofOutputLocation;//! output location in PROOF mode
  Bool_t         fProofOutputDataset; //! keep output file in PROOF mode in dataset or merge
  TString        fProofOutputArchive; //! list of file to be kept in output zip archive in PROOF mode
  TString        fEquipIdMap;         // name of file with equipment id map
  Int_t          fFirstEvent;         // index of first event to be reconstr.
  Int_t          fLastEvent;          // index of last event to be reconstr.
  UInt_t         fNumberOfEventsPerFile; // number of events per file in case of raw-data reconstruction
  Double32_t     fFractionFriends; // fraction of ESD friends to be stored
  TObjArray      fOptions;            // options for reconstructor objects
  Bool_t         fLoadAlignFromCDB;   // Load alignment data from CDB and apply it to geometry or not
  TString        fLoadAlignData;      // Load alignment data from CDB for these detectors
  TString        fUseHLTData;        // Detectors for which the HLT data is used as input
  AliRunInfo*    fRunInfo;            // an object which contains essential global conditions information
  AliEventInfo   fEventInfo;          // an object which contains essential event information
  AliTriggerRunScalers* fRunScalers;  // an object which contains the trigger counters for xsection
  AliCTPTimeParams* fCTPTimeParams;   // an object which contains the trigger delays
  AliCTPTimeParams* fCTPTimeAlign;    // an object which contains the CTP (mis)alignment data
  AliRunLoader*  fRunLoader;          //! current run loader object
  AliRawReader*  fRawReader;          //! current raw data reader
  AliRawReader*  fParentRawReader;    //! parent raw data reader in case of AliRawReaderHLT

  static const char* fgkDetectorName[kNDetectors]; //! names of detectors
  AliReconstructor*  fReconstructor[kNDetectors];  //! array of reconstructor objects
  AliRecoParam   fRecoParam;                      // container for the reco-param objects for detectors
  AliLoader*     fLoader[kNDetectors];   //! detector loaders
  AliTracker*    fTracker[kNDetectors];  //! trackers
  AliTracker*    fSPDTrackleter;           //! trackleter [for SPD Plane Efficiency purposes]
  AliESDVertex*  fDiamondProfileSPD;       // (x,y) diamond profile from SPD for AliITSVertexer3D(Z)
  AliESDVertex*  fDiamondProfile;          // (x,y) diamond profile for AliVertexerTracks (ITS+TPC)
  AliESDVertex*  fDiamondProfileTPC;       // (x,y) diamond profile from TPC for AliVertexerTracks
  THashTable*    fListOfCosmicTriggers;    // list of cosmic triggers as defined by trigger coordination

  AliGRPObject*  fGRPData;              // Data from the GRP/GRP/Data CDB folder

  TObjArray* 	 fAlignObjArray;      //! array with the alignment objects to be applied to the geometry

  TString	 fCDBUri;	      //! Uri of the default CDB storage
  TString	 fQARefUri;	    //! Uri of the default QA reference storage
  TObjArray      fSpecCDBUri;         //! Array with detector specific CDB storages
  Bool_t 	 fInitCDBCalled;               //! flag to check if CDB storages are already initialized
  Bool_t 	 fFromCDBSnapshot;             //! flag to check if we are loading the CDB from a snapshot
  TString        fSnapshotFileName;	       //! string for the file containing the CDB snapshot
  Bool_t 	 fSetRunNumberFromDataCalled;  //! flag to check if run number is already loaded from run loader

  //Quality Assurance
  Int_t  fQACycles[     AliQAv1::kNDET];  // # events over which QA data are accumulated
  Bool_t fQAWriteExpert[AliQAv1::kNDET];  // Flag to save or not expert QA data
  TString               fQADetectors ;  // list of detectors to be QA'ed 	
  TString               fQATasks ;      // list of QA tasks to be performed	
  Bool_t                fRunQA ;        // Run QA flag
  Bool_t                fRunGlobalQA;   // Run global QA flag
  Bool_t                fSameQACycle;   //! open a new QA data file or not
  Bool_t                fInitQACalled;  //! control of QA intialisation
  Bool_t                fWriteQAExpertData ; //! decides wheter or not to write experts QA data; true by default

  // Plane Efficiency Evaluation
  Bool_t         fRunPlaneEff ;      // Evaluate Plane Efficiency

  // New members needed in order to split Run method
  // into InitRun,RunEvent,FinishRun methods
  AliESDEvent*         fesd;        //! Pointer to the ESD event object
  AliESDEvent*         fhltesd;     //! Pointer to the HLT ESD event object
  AliESDfriend*        fesdf;       //! Pointer to the ESD friend object
  TFile*               ffile;       //! Pointer to the ESD file
  TFile*               ffileF;      //! Pointer to the ESD friend file
  TTree*               ftree;       //! Pointer to the ESD tree
  TTree*               ftreeF;      //! Pointer to the ESD friend tree
  TTree*               fhlttree;    //! Pointer to the HLT ESD tree
  AliVertexerTracks*   ftVertexer;  //! Pointer to the vertexer based on ESD tracks
  Bool_t               fIsNewRunLoader; // galice.root created from scratch (real raw data case)
  Bool_t               fRunAliEVE;  // Run AliEVE or not
  Double_t             fBeamInt[2][2]; //! mean beam intensities from AliLHCData
  
  TTree*              fChain;      //! The raw-data chain in case of AliRawReaderChain

  // Counters used in WriteESDfriend()
  Int_t fNall;    //! Number of reconstructed events
  Int_t fNspecie; //! Number of events, selected by IsCosmicOrCalibSpecie
  Int_t fSspecie; //! Number of events, sampled from fNspecie
  Int_t fNhighPt; //! Number of events, selected by IsHighPt 
  Int_t fShighPt; //! Number of events, sampled from fNhighPt
 
  // Upgrade detector reconstruction
  TString fUpgradeModule;
  Bool_t  fUpgradeMask[kNDetectors];
  // Connect an analysis manager in the event loop
  TString              fAnalysisMacro; // Full path to a macro creating an analysis manager train
  AliAnalysisManager  *fAnalysis;      //! Analysis manager
  AliRecoInputHandler *fRecoHandler;   //! Input handler adapted for reconstruction
  ClassDef(AliReconstruction, 41)      // class for running the reconstruction
};

#endif
