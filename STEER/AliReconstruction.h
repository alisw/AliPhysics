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


#include <TNamed.h>
#include <TString.h>
#include <TObjArray.h>

class AliReconstructor;
class AliRunLoader;
class AliRawReader;
class AliLoader;
class AliTracker;
class AliVertexer;
class AliESDVertex;
class AliESDEvent;
class TFile;
class TTree;
class TList;
class AliQADataMaker;

class AliReconstruction: public TNamed {
public:
  AliReconstruction(const char* gAliceFilename = "galice.root",
  		    const char* cdbUri = "local://$ALICE_ROOT",
		    const char* name = "AliReconstruction", 
		    const char* title = "reconstruction");
  AliReconstruction(const AliReconstruction& rec);
  AliReconstruction& operator = (const AliReconstruction& rec);
  virtual ~AliReconstruction();

  void           SetGAliceFile(const char* fileName);
  void           SetInput(const char* input) {fInput = input;};
  void           SetEquipmentIdMap(const char *mapFile) {fEquipIdMap = mapFile;};
  void           SetEventRange(Int_t firstEvent = 0, Int_t lastEvent = -1) 
    {fFirstEvent = firstEvent; fLastEvent = lastEvent;};
  void           SetNumberOfEventsPerFile(UInt_t nEvents)
    {fNumberOfEventsPerFile = nEvents;};
  void           SetOption(const char* detector, const char* option);

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
  void           SetESDParLocation(const char *c){fESDPar = c;}
  TNamed *CopyFileToTNamed(TString fPath,TString fName);
  void  TNamedToFile(TTree* fTree, TString fName);

  //*** Global reconstruction flag setters
  void SetUniformFieldTracking(Bool_t flag=kTRUE){fUniformField=flag;} 
  void SetRunVertexFinder(Bool_t flag=kTRUE) {fRunVertexFinder=flag;};
  void SetRunHLTTracking(Bool_t flag=kTRUE) {fRunHLTTracking=flag;};
  void SetRunV0Finder(Bool_t flag=kTRUE) {fRunV0Finder=flag;};
  void SetRunCascadeFinder(Bool_t flag=kTRUE) {fRunCascadeFinder=flag;};
  void SetStopOnError(Bool_t flag=kTRUE) {fStopOnError=flag;}
  void SetWriteAlignmentData(Bool_t flag=kTRUE){fWriteAlignmentData=flag;}
  void SetWriteESDfriend(Bool_t flag=kTRUE){fWriteESDfriend=flag;}
  void SetWriteAOD(Bool_t flag=kTRUE){fWriteAOD=flag;}
  void SetFillTriggerESD(Bool_t flag=kTRUE){fFillTriggerESD=flag;}
  void SetDiamondProfile(AliESDVertex *dp) {fDiamondProfile=dp;}
		   
  void SetCleanESD(Bool_t flag=kTRUE){fCleanESD=flag;}
  void SetV0DCAmax(Float_t d) {fV0DCAmax=d;}
  void SetV0CsPmin(Float_t d) {fV0CsPmin=d;}
  void SetDmax(Float_t d) {fDmax=d;}
  void SetZmax(Float_t z) {fZmax=z;}
  Float_t GetV0DCAmax() const {return fV0DCAmax;}
  Float_t GetV0CsPmin() const {return fV0CsPmin;}
  Float_t GetDmax() const {return fDmax;}
  Float_t GetZmax() const {return fZmax;}

  void           SetCheckPointLevel(Int_t checkPointLevel)
    {fCheckPointLevel = checkPointLevel;}
  // CDB storage activation
  void InitCDBStorage();
  void SetDefaultStorage(const char* uri);
  void SetRemoteStorage(const char* uri);
  void SetSpecificStorage(const char* calibType, const char* uri);

  Bool_t SetRunNumber();

  Bool_t MisalignGeometry(const TString& detectors);

  void           SetAlignObjArray(TObjArray *array)
                   {fAlignObjArray = array;
		   fLoadAlignFromCDB = kFALSE;}

  virtual Bool_t Run(const char* input = NULL);

  // AOD creation
  void ESDFile2AODFile(TFile* esdFile, TFile* aodFile);

  // Quality Assurance 
//  virtual Bool_t RunQA(const char* detectors, AliESDEvent *& esd);
//  void    SetQACycles(const char * detector, const Int_t cycles) { fQACycles[GetDetIndex(detector)] = cycles ; }

private:
  Bool_t         RunLocalReconstruction(const TString& detectors);
  Bool_t         RunLocalEventReconstruction(const TString& detectors);
  Bool_t         RunVertexFinder(AliESDEvent*& esd);
  Bool_t         RunHLTTracking(AliESDEvent*& esd);
  Bool_t         RunMuonTracking(AliESDEvent*& esd);
  Bool_t         RunTracking(AliESDEvent*& esd);
  Bool_t         CleanESD(AliESDEvent *esd);
  Bool_t         FillESD(AliESDEvent*& esd, const TString& detectors);
  Bool_t         FillTriggerESD(AliESDEvent*& esd);
  Bool_t         FillRawEventHeaderESD(AliESDEvent*& esd);

  Bool_t         IsSelected(TString detName, TString& detectors) const;
  Bool_t         InitRunLoader();
  AliReconstructor* GetReconstructor(Int_t iDet);
  Bool_t         CreateVertexer();
  Bool_t         CreateTrackers(const TString& detectors);
  void           CleanUp(TFile* file = NULL, TFile* fileOld = NULL);

  Bool_t         ReadESD(AliESDEvent*& esd, const char* recStep) const;
  void           WriteESD(AliESDEvent* esd, const char* recStep) const;

  //==========================================//
  void           WriteAlignmentData(AliESDEvent* esd);

  void           FillRawDataErrorLog(Int_t iEvent, AliESDEvent* esd);

  //Quality Assurance
  Int_t                GetDetIndex(const char * detector);
//  AliQADataMaker*      GetQADataMaker(Int_t iDet);
//  const Int_t          GetQACycles(const char * detector) { return fQACycles[GetDetIndex(detector)] ; }
  void                 CheckQA() ;
  void        SetQA(const Bool_t val) { fRunQA = val ; } 

  //*** Global reconstruction flags *******************
  Bool_t         fUniformField;       // uniform field tracking flag
  Bool_t         fRunVertexFinder;    // run the vertex finder
  Bool_t         fRunHLTTracking;     // run the HLT tracking
  Bool_t         fRunMuonTracking;    // run the HLT tracking
  Bool_t         fRunV0Finder;        // run the ESD V0 finder
  Bool_t         fRunCascadeFinder;   // run the ESD cascade finder
  Bool_t         fStopOnError;        // stop or continue on errors
  Bool_t         fWriteAlignmentData; // write track space-points flag
  Bool_t         fWriteESDfriend;     // write ESD friend flag
  Bool_t         fWriteAOD;           // write AOD flag
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
  TString        fUseTrackingErrorsForAlignment; // for these detectors
  TString        fGAliceFileName;     // name of the galice file
  TString        fInput;              // name of input file or directory
  TString        fEquipIdMap;         // name of file with equipment id map
  Int_t          fFirstEvent;         // index of first event to be reconstr.
  Int_t          fLastEvent;          // index of last event to be reconstr.
  UInt_t         fNumberOfEventsPerFile; // number of events per file in case of raw-data reconstruction
  Int_t          fCheckPointLevel;    // level of ESD check points
  TObjArray      fOptions;            // options for reconstructor objects
  Bool_t         fLoadAlignFromCDB;   // Load alignment data from CDB and apply it to geometry or not
  TString        fLoadAlignData;      // Load alignment data from CDB for these detectors
  TString        fESDPar;             // String where the esd.par is stored, will be attached to the tree         


  AliRunLoader*  fRunLoader;          //! current run loader object
  AliRawReader*  fRawReader;          //! current raw data reader

  static const Int_t fgkNDetectors = 15;   //! number of detectors
  static const char* fgkDetectorName[fgkNDetectors]; //! names of detectors
  AliReconstructor*  fReconstructor[fgkNDetectors];  //! array of reconstructor objects
  AliLoader*     fLoader[fgkNDetectors];   //! detector loaders
  AliVertexer*   fVertexer;                //! vertexer for ITS
  AliTracker*    fTracker[fgkNDetectors];  //! trackers
  AliESDVertex*  fDiamondProfile;          // (x,y) diamond profile for AliVertexerTracks

  TList*         fGRPList;            // TList from the GRP/GRP/Data CDB folder

  TObjArray* 	 fAlignObjArray;      // array with the alignment objects to be applied to the geometry

  TString	 fCDBUri;	      // Uri of the default CDB storage
  TString	 fRemoteCDBUri;	      // Uri of the remote CDB storage
  TObjArray      fSpecCDBUri;         // Array with detector specific CDB storages

  //Quality Assurance
//  AliQADataMaker * fQADataMaker[fgkNDetectors];  //! array of QA data maker objects
//  Int_t                 fQACycles[fgkNDetectors] ;         // cycle length (# events) over which QA data are accumulated
  Bool_t              fRunQA ;                        // Runs the QA at the end of simulation

  ClassDef(AliReconstruction, 17)      // class for running the reconstruction
};

#endif
