#ifndef ALIMUONREALIGNTASK_H
#define ALIMUONREALIGNTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calign
/// \class AliMUONReAlignTask
/// \brief Task to refit ESD tracks with relaigned geometry
///
//  Author Javier Castillo, CEA/Saclay - Irfu/SPhN

class TTree;
class TString;
class AliESDEvent;
class AliMUONClusterInfo;
class AliMUONESDInterface;
class AliMUONRefitter;
class AliMUONRecoParam;
class AliMUONGeometryTransformer;
class AliMUONVStore;
class AliMUONTrack;

#include "AliAnalysisTask.h"

class AliMUONReAlignTask : public AliAnalysisTask {
 public:
  AliMUONReAlignTask(const char *name = "AliMUONReAlignTask", const char *geofilename = "geometry.root", const char *defaultocdb = "local://$ALICE_ROOT/OCDB", const char *misalignocdb = "local://ReAlignOCDB");
  AliMUONReAlignTask(const AliMUONReAlignTask& obj);
  AliMUONReAlignTask& operator=(const AliMUONReAlignTask& other); 
  virtual ~AliMUONReAlignTask();
  
  virtual void   LocalInit();
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(const Option_t*);

  /// Set geoemetry file name
  void SetGeoFilename(const char* geoFilename) {fGeoFilename = geoFilename;}
  /// Set mis align ocdb
  void SetMisAlignOCDB(const char* misalignOCDB) {fMisAlignOCDB = misalignOCDB;}
  /// Set default ocdb
  void SetDefaultOCDB(const char* defaultOCDB) {fDefaultOCDB = defaultOCDB;}
  void Prepare(const char* geoFilename, const char* defaultOCDB, const char* misalignOCDB);
  UInt_t BuildClusterMap(AliMUONTrack &track);

 private:
  AliESDEvent *fESD;                  ///< ESD object
  TTree *fClusterInfoTree;            ///< ClusterInfo tree 
  AliMUONClusterInfo *fClusterInfo;   ///< ClusterInfo object
  AliMUONESDInterface *fESDInterface; //!< MUONESDInterface
  AliMUONRefitter  *fRefitter;        //!< The refitter class 
  AliMUONRecoParam *fRecoParam;       //!< Parameters for reconstruction
  TString fGeoFilename;               ///< Geometry file name
  TString fMisAlignOCDB;              ///< OCDB with misalignment file
  TString fDefaultOCDB;               ///< Default OCDB
  AliMUONGeometryTransformer *fGeoTransformer;    //!< Original geometry
  AliMUONGeometryTransformer *fNewGeoTransformer; //!< Aligned geometry
  AliMUONVStore *fGainStore; ///< Store for gains
  AliMUONVStore *fPedStore;  ///< Store for pedestals
  Int_t fPrintLevel;         //!< Print information
  Int_t fLastRun;            //!< Last run number

  ClassDef(AliMUONReAlignTask, 1) // example of analysis
};

#endif

