#ifndef ALISIMULATION_H
#define ALISIMULATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliSimulation.h 63170 2013-06-24 14:20:10Z morsch $ */

//
// class for running generation, simulation and digitization
// Hits, sdigits and digits are created for all detectors by typing:
//   AliSimulation sim;
//   sim.Run();
//

#include <TNamed.h>
#include <TString.h>
#include <TObjArray.h>
#include "AliQAv1.h"
#include "AliQAManager.h"
#include <time.h>
#include <algorithm>

class AliCDBId;
class AliCDBParam;
class AliRunLoader;
class AliLegoGenerator;
class AliLego;
class AliMagF;
class AliHLTSimulation;

class AliSimulation: public TNamed {
public:
  AliSimulation(const char* configFileName = "Config.C",
		const char* name = "AliSimulation", 
		const char* title = "generation, simulation and digitization");
  virtual ~AliSimulation();

  static AliSimulation *Instance() {if(!fgInstance) fgInstance=new AliSimulation(); return fgInstance;}

  void           SetNumberOfEvents(Int_t nEvents);
  void           SetConfigFile(const char* fileName);
  void           SetGAliceFile(const char* fileName);
  void           SetEventsPerFile(const char* detector, const char* type, 
				  Int_t nEvents);

  void           SetRunGeneratorOnly(Bool_t val) {fRunGeneratorOnly = val;};
  void           SetRunGeneration(Bool_t run) {fRunGeneration = run;};
  void           SetRunSimulation(Bool_t run) {fRunSimulation = run;};
  void           SetLoadAlignFromCDB(Bool_t load)  {fLoadAlignFromCDB = load;};
  void           SetLoadAlignData(const char* detectors) 
                   {fLoadAlObjsListOfDets = detectors;};
  void           SetMakeSDigits(const char* detectors) 
                   {fMakeSDigits = detectors;};
  void           MergeWith(const char* fileName, Int_t nSignalPerBkgrd = 0);
  void           EmbedInto(const char* fileName, Int_t nSignalPerBkgrd = 0);
  void           SetUseBkgrdVertex(Bool_t useBkgrdVertex)
                   {fUseBkgrdVertex = useBkgrdVertex;};
  void           SetRegionOfInterest(Bool_t flag) {fRegionOfInterest = flag;};
  void           SetMakeDigits(const char* detectors)
                   {fMakeDigits = detectors;};
  void           SetMakeDigitsFromHits(const char* detectors)
                   {fMakeDigitsFromHits = detectors;};
  void           SetWriteRawData(const char* detectors, 
				 const char* fileName = NULL,
				 Bool_t deleteIntermediateFiles = kFALSE)
                   {fWriteRawData = detectors; fRawDataFileName = fileName;
		   fDeleteIntermediateFiles = deleteIntermediateFiles;};
  void           SetWriteSelRawData(Bool_t sel = kTRUE)
                   {fWriteSelRawData = sel;}
  void           SetTriggerConfig(TString conf) {fTriggerConfig=conf;}
  const Char_t*  GetTriggerConfig() const {return fTriggerConfig.Data();}
  void           SetAlignObjArray(TObjArray *array)
                   {fAlignObjArray = array;
		   fLoadAlignFromCDB = kFALSE;}
  void           SetUseMonitoring(Bool_t flag=kTRUE) {fUseMonitoring = flag;}

  Bool_t         MisalignGeometry(AliRunLoader *runLoader = NULL);

  void           SetRunNumber(Int_t run);
  void           SetSeed(Int_t seed);
    
  void 		 ProcessEnvironmentVars();
		   
  // CDB storage activation
  void SetDefaultStorage(const char* uri);
  void SetSpecificStorage(const char* calibType, const char* uri);
  void SetCDBSnapshotMode(const char* snapshotFileName) {AliCDBManager::Instance()->SetSnapshotMode(snapshotFileName);}

  virtual Bool_t Run(Int_t nEvents = 0);
  virtual Bool_t RunLego(const char *setup="Config.C",Int_t nc1=60,Float_t c1min=2,Float_t c1max=178,
			 Int_t nc2=60,Float_t c2min=0,Float_t c2max=360,Float_t rmin=0,
			 Float_t rmax=430,Float_t zmax=10000, AliLegoGenerator* gener=NULL, Int_t nev = -1);

  virtual Bool_t RunSimulation(Int_t nEvents = 0);
  virtual Bool_t RunGeneratorOnly();
  virtual Bool_t RunSDigitization(const char* detectors = "ALL");
  virtual Bool_t RunTrigger(const char* descriptors ="", const char* detectors = "ALL");
  virtual Bool_t WriteTriggerRawData();
  virtual Bool_t RunDigitization(const char* detectors = "ALL",
				 const char* excludeDetectors = "");
  virtual Bool_t RunHitsDigitization(const char* detectors = "ALL");
  virtual Bool_t WriteRawData(const char* detectors = "ALL",
			      const char* fileName = NULL,
			      Bool_t deleteIntermediateFiles = kFALSE,
			      Bool_t selrawdata = kFALSE);
  virtual Bool_t WriteRawFiles(const char* detectors = "ALL");
  virtual Bool_t ConvertRawFilesToDate(const char* dateFileName = "raw.date",
				       const char* rootFileName = "");
  virtual Bool_t ConvertDateToRoot(const char* dateFileName = "raw.date",
				   const char* rootFileName = "raw.root");
  virtual Int_t ConvertRaw2SDigits(const char* rawDirectory, const char* esdFile = "", Int_t N=-1, Int_t nSkip=0);

  // Sets the name of the file from which the geometry is loaded
  virtual void SetGeometryFile(const Char_t* filename) {fGeometryFile=filename;}
  virtual const Char_t* GetGeometryFile() const {return fGeometryFile.Data();}
  virtual Bool_t IsGeometryFromFile() const {return !fGeometryFile.IsNull();}

  void  DeactivateDetectorsAbsentInGRP(TObjArray* detArr);

  // HLT
  void SetRunHLT(const char* options) {fRunHLT=options;}
  virtual Bool_t CreateHLT();
  virtual Bool_t RunHLT();
  virtual  Bool_t IsLegoRun() const {return (fLego!=0);}
  AliLego* Lego() const {return fLego;}
  virtual  void  FinishRun();
  //
  void StoreUsedCDBMaps() const; 

  //Quality Assurance
  Int_t       GetDetIndex(const char * detector);
  void        SetQACycles(AliQAv1::DETECTORINDEX_t det, const Int_t cycles) {  AliQAManager::QAManager()->SetCycleLength(det, cycles) ; }
  Bool_t      RunQA() ;
  Bool_t      SetRunQA(TString detAndAction="ALL:ALL") ; 
  void        SetQAWriteExpert(AliQAv1::DETECTORINDEX_t det) { AliQAManager::QAManager()->SetWriteExpert(det) ; }  
  void        SetQARefDefaultStorage(const char* uri);
  void        InitQA();
  void        SetEventSpecie(AliRecoParam::EventSpecie_t es) { fEventSpecie = es ; }
  void        SetWriteQAExpert() { fWriteQAExpertData = kTRUE ; }

  void SetWriteGRPEntry(Bool_t flag = kTRUE) { fWriteGRPEntry = flag; }
  void WriteGRPEntry();
  void UseVertexFromCDB()   {fUseVertexFromCDB   = kTRUE;}
  void UseMagFieldFromGRP() {fUseMagFieldFromGRP = kTRUE;} 
  void SetGRPWriteLocation(char* loc) {fGRPWriteLocation = loc;}

  Bool_t GetUseTimeStampFromCDB()  const {return fUseTimeStampFromCDB;}
  void UseTimeStampFromCDB(Double_t decayTimeHours=-1);
  time_t GenerateTimeStamp() const;
  //
  Bool_t          GetUseDetectorsFromGRP()               const {return fUseDetectorsFromGRP;}
  void            SetUseDetectorsFromGRP(Bool_t v=kTRUE)       {fUseDetectorsFromGRP = v;}
  //
private:

  AliSimulation(const AliSimulation&); // Not implemented
  AliSimulation& operator = (const AliSimulation&); // Not implemented

  void 		 InitCDB();
  void 		 InitRunNumber();
  void 		 SetCDBLock();
  Bool_t         SetRunNumberFromData();
  AliRunLoader*  LoadRun(const char* mode = "UPDATE") const;
  Int_t          GetNSignalPerBkgrd(Int_t nEvents = 0) const;
  Bool_t         IsSelected(TString detName, TString& detectors) const;

  static AliSimulation *fgInstance;    // Static pointer to object

  Bool_t         fRunGeneratorOnly;   // run code for a generator only production
  Bool_t         fRunGeneration;      // generate prim. particles or not
  Bool_t         fRunSimulation;      // simulate detectors (hits) or not
  Bool_t         fLoadAlignFromCDB;   // Load alignment data from CDB and apply it to geometry or not
  TString        fLoadAlObjsListOfDets;   // Load alignment data from CDB for these detectors
  TString        fMakeSDigits;        // create sdigits for these detectors
  TString        fMakeDigits;         // create digits for these detectors
  TString        fTriggerConfig;      // run trigger for these descriptors
  TString        fMakeDigitsFromHits; // create digits from hits for these detectors
  TString        fWriteRawData;       // write raw data for these detectors
  TString        fRawDataFileName;    // file name for the raw data file
  Bool_t         fDeleteIntermediateFiles; // delete intermediate raw data files
  Bool_t         fWriteSelRawData;    // write detectors raw data in a separate file accoring to the trigger cluster
  Bool_t         fStopOnError;        // stop or continue on errors
  Bool_t         fUseMonitoring;      // monitor simulation timing per volume

  Int_t          fNEvents;            // number of events
  TString        fConfigFileName;     // name of the config file
  TString        fGAliceFileName;     // name of the galice file
  TObjArray      fEventsPerFile;      // number of events per file for given detectors and data types

  TObjArray*     fBkgrdFileNames;     // names of background files for merging
  TObjArray*     fAlignObjArray;      // array with the alignment objects to be applied to the geometry
  Bool_t         fUseBkgrdVertex;     // use vertex from background in case of merging
  Bool_t         fRegionOfInterest;   // digitization in region of interest

  TString 	 fCDBUri;	                     //! Uri of the default CDB storage
  TString 	 fQARefUri;	                     //! Uri of the default QA reference storage
  TObjArray      fSpecCDBUri;                        //! Array with detector specific CDB storages
  Int_t 	   fRun; 		                     //! Run number, will be passed to CDB and gAlice!!
  Int_t 	   fSeed;                        //! Seed for random number generator 
  Bool_t 	   fInitCDBCalled;               //! flag to check if CDB storages are already initialized
  Bool_t 	   fInitRunNumberCalled;         //! flag to check if run number is already initialized
  Bool_t 	   fSetRunNumberFromDataCalled;  //! flag to check if run number is already loaded from run loader
  
  Bool_t         fEmbeddingFlag;       // Flag for embedding
  AliLego       *fLego;                //! Pointer to aliLego object if it exists
  // OCDB
  ULong64_t       fKey;                //! current CDB key
  Bool_t          fUseVertexFromCDB;   // Flag to use Vertex from CDB
  Bool_t          fUseMagFieldFromGRP; // Use magnetic field settings from GRP
  TString         fGRPWriteLocation;   // Location to write the GRP entry from simulation
  
  Bool_t          fUseDetectorsFromGRP; // do not simulate detectors absent in the GRP

  Bool_t          fUseTimeStampFromCDB;// Flag to generate event time-stamps according to SOR/EOR from GRP
  time_t          fTimeStart;          // SOR time-stamp
  time_t          fTimeEnd;            // EOR time-stamp
  Float_t         fLumiDecayH;         // luminosity decay time in hours
  std::vector<time_t> fOrderedTimeStamps;  //! optional ordered time stamps
  
  //QA stuff
  static const Int_t   fgkNDetectors = 18 ;             // number of detectors   +FIT // alla
  static const char *  fgkDetectorName[fgkNDetectors] ; // names of detectors
  TString              fQADetectors ;                   // list of detectors to be QA'ed 	
  TString              fQATasks ;                       // list of QA tasks to be performed	
  Bool_t               fRunQA ;                         // Runs the QA at the end of simulation
  AliRecoParam::EventSpecie_t fEventSpecie ;            // type of event (see AliRecoParam::EventSpecie_t)
  Bool_t               fWriteQAExpertData ;             //! decides wheter or not to write experts QA data; true by default

  TString              fGeometryFile;                   // Geometry file

  //HLT
  TString              fRunHLT;       //! HLT options, HLT is disabled if empty, default='default'
  AliHLTSimulation*    fpHLT;         //! The instance of HLT simulation

  Bool_t         fWriteGRPEntry;      // Write or not GRP entry corresponding to the settings in Config.C

  static const Char_t *fgkRunHLTAuto;         // flag for automatic HLT mode detection
  static const Char_t *fgkHLTDefConf;         // default configuration to run HLT

  ClassDef(AliSimulation, 15)  // class for running generation, simulation and digitization
};

#endif
