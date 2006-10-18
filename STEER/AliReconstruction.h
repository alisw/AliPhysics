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
class AliESD;
class TFile;

class AliRunTag;
class AliLHCTag;
class AliDetectorTag;
class AliEventTag;


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
  void           SetLoadAlignFromCDB(Bool_t load)  {fLoadAlignFromCDB = load;};
  void           SetLoadAlignData(const char* detectors) 
    {fLoadAlignData = detectors;};


  //*** Global reconstruction flag setters
  void SetUniformFieldTracking(Bool_t flag=kTRUE){fUniformField=flag;} 
  void SetRunVertexFinder(Bool_t flag=kTRUE) {fRunVertexFinder=flag;};
  void SetRunHLTTracking(Bool_t flag=kTRUE) {fRunHLTTracking=flag;};
  void SetStopOnError(Bool_t flag=kTRUE) {fStopOnError=flag;}
  void SetWriteAlignmentData(Bool_t flag=kTRUE){fWriteAlignmentData=flag;}
  void SetWriteESDfriend(Bool_t flag=kTRUE){fWriteESDfriend=flag;}
  void SetFillTriggerESD(Bool_t flag=kTRUE){fFillTriggerESD=flag;}
  void SetDiamondProfile(AliESDVertex *dp) {fDiamondProfile=dp;}
		   
  void           SetCheckPointLevel(Int_t checkPointLevel)
    {fCheckPointLevel = checkPointLevel;}
  // CDB storage activation
  void InitCDBStorage();
  void SetDefaultStorage(const char* uri);
  void SetSpecificStorage(const char* calibType, const char* uri);

  Bool_t SetRunNumber();

  Bool_t SetAlignObjArraySingleDet(const char* detName);
  Bool_t MisalignGeometry(const TString& detectors);

  void           SetAlignObjArray(TObjArray *array)
                   {fAlignObjArray = array;
		   fLoadAlignFromCDB = kFALSE;}
  Bool_t         ApplyAlignObjsToGeom(TObjArray* alObjArray);

  virtual Bool_t Run(const char* input, 
		     Int_t firstEvent, Int_t lastEvent = -1);
  Bool_t         Run(const char* input = NULL)
    {return Run(input, fFirstEvent, fLastEvent);};
  Bool_t         Run(Int_t firstEvent, Int_t lastEvent = -1)
    {return Run(NULL, firstEvent, lastEvent);};


private:
  Bool_t         RunLocalReconstruction(const TString& detectors);
  Bool_t         RunLocalEventReconstruction(const TString& detectors);
  Bool_t         RunVertexFinder(AliESD*& esd);
  Bool_t         RunHLTTracking(AliESD*& esd);
  Bool_t         RunTracking(AliESD*& esd);
  Bool_t         FillESD(AliESD*& esd, const TString& detectors);
  Bool_t         FillTriggerESD(AliESD*& esd);

  Bool_t         IsSelected(TString detName, TString& detectors) const;
  Bool_t         InitRunLoader();
  AliReconstructor* GetReconstructor(Int_t iDet);
  Bool_t         CreateVertexer();
  Bool_t         CreateTrackers(const TString& detectors);
  void           CleanUp(TFile* file = NULL, TFile* fileOld = NULL);

  Bool_t         ReadESD(AliESD*& esd, const char* recStep) const;
  void           WriteESD(AliESD* esd, const char* recStep) const;

 
  //===========================================//
  void           CreateTag(TFile* file);
  //==========================================//

  void           WriteAlignmentData(AliESD* esd);


  //*** Global reconstruction flags *******************
  Bool_t         fUniformField;       // uniform field tracking flag
  Bool_t         fRunVertexFinder;    // run the vertex finder
  Bool_t         fRunHLTTracking;     // run the HLT tracking
  Bool_t         fStopOnError;        // stop or continue on errors
  Bool_t         fWriteAlignmentData; // write track space-points flag
  Bool_t         fWriteESDfriend;     // write ESD friend flag
  Bool_t         fFillTriggerESD;     // fill trigger info into ESD

  TString        fRunLocalReconstruction; // run the local reconstruction for these detectors
  TString        fRunTracking;        // run the tracking for these detectors
  TString        fFillESD;            // fill ESD for these detectors
  TString        fGAliceFileName;     // name of the galice file
  TString        fInput;              // name of input file or directory
  TString        fEquipIdMap;         // name of file with equipment id map
  Int_t          fFirstEvent;         // index of first event to be reconstr.
  Int_t          fLastEvent;          // index of last event to be reconstr.
  Int_t          fCheckPointLevel;    // level of ESD check points
  TObjArray      fOptions;            // options for reconstructor objects
  Bool_t         fLoadAlignFromCDB;   // Load alignment data from CDB and apply it to geometry or not
  TString        fLoadAlignData;      // Load alignment data from CDB for these detectors

  AliRunLoader*  fRunLoader;          //! current run loader object
  AliRawReader*  fRawReader;          //! current raw data reader

  static const Int_t fgkNDetectors = 15;   //! number of detectors
  static const char* fgkDetectorName[fgkNDetectors]; //! names of detectors
  AliReconstructor*  fReconstructor[fgkNDetectors];  //! array of reconstructor objects
  AliLoader*     fLoader[fgkNDetectors];   //! detector loaders
  AliVertexer*   fVertexer;                //! vertexer for ITS
  AliTracker*    fTracker[fgkNDetectors];  //! trackers
  AliESDVertex*  fDiamondProfile;          // (x,y) diamond profile for AliVertexerTracks

  TObjArray* 	 fAlignObjArray;      // array with the alignment objects to be applied to the geometry

  TString	 fCDBUri;	      // Uri of the default CDB storage
  TObjArray      fSpecCDBUri;         // Array with detector specific CDB storages

  ClassDef(AliReconstruction, 9)      // class for running the reconstruction
};

#endif
