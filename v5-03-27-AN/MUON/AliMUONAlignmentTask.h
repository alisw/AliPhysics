#ifndef ALIMUONALIGNMENTTASK_H
#define ALIMUONALIGNMENTTASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calign
/// \class AliMUONAlignmentTask
/// \brief Task to align the muon spectrometer
///
//  Author Javier Castillo, CEA/Saclay - Irfu/SPhN

class TList;
class TString;
class TGraphErrors;
class TClonesArray;
class AliMUONAlignment;
class AliMUONGeometryTransformer;

#include "AliAnalysisTaskSE.h"

class AliMUONAlignmentTask : public AliAnalysisTaskSE
{

  public:


  /// constructor
  AliMUONAlignmentTask(const char *name = "AliMUONAlignmentTask", const char *newalignocdb = "local://ReAlignOCDB", const char *oldalignocdb = "none", const char *defaultocdb = "raw://", const char *geofilename = "geometry.root");

  /// copy constructor
  AliMUONAlignmentTask(const AliMUONAlignmentTask& obj);

  /// asignment operator
  AliMUONAlignmentTask& operator=(const AliMUONAlignmentTask& other);

  /// destructor
  virtual ~AliMUONAlignmentTask();

  virtual void   LocalInit();

  ///@name flags
  //@{

  /// perform alignment from alignment records is true. Use ESD tracks otherwise
  void SetReadRecords( Bool_t value = kTRUE )
  { fReadRecords = value; }

  /// write alignment records to AOD if true
  void SetWriteRecords( Bool_t value = kTRUE )
  { fWriteRecords = value; }

  /// perform alignment (from either tracks or records depending on fReadRecords)
  void SetDoAlignment( Bool_t value )
  { fDoAlignment = value; }

  //@}

  /// output data
  virtual void UserCreateOutputObjects();

  /// per-event method
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();

  /// termination cleanup
  virtual void Terminate(const Option_t*);

  /// end of task execution
  virtual void FinishTaskOutput();

  /// Set geoemetry file name
  void SetGeoFilename(const char* geoFilename)
  {fGeoFilename = geoFilename;}

  /// Set default ocdb
  void SetDefaultStorage(const char* defaultOCDB)
  {fDefaultStorage = defaultOCDB;}

  /// Set mis align ocdb
  void SetOldAlignStorage(const char* oldalignOCDB)
  {fOldAlignStorage = oldalignOCDB;}

  /// Set mis align ocdb
  void SetNewAlignStorage(const char* newalignOCDB)
  {fNewAlignStorage = newalignOCDB;}

  /// Flag to set OCDB once at first notify
  void SetLoadOCDBOnce(Bool_t loadOCDBOnce = kTRUE)
  {fLoadOCDBOnce = loadOCDBOnce;}

  private:

  ///@name flags
  //@{

  /// perform alignment from alignment records is true. Use ESD tracks otherwise
  Bool_t fReadRecords;

  /// write alignment records to AOD if true
  Bool_t fWriteRecords;

  /// perform alignment (from either tracks or records depending on fReadRecords)
  Bool_t fDoAlignment;

  //@}

  /// The MUON alignment object
  AliMUONAlignment *fAlign;

  /// Geometry file name
  TString fGeoFilename;

  /// location of the default OCDB storage
  TString fDefaultStorage;

  /// location of the OCDB storage where to find old MUON/Align/Data (use the default one if empty)
  TString fOldAlignStorage;

  /// location of the OCDB storage where to put new MUON/Align/Data (use the default one if empty)
  TString fNewAlignStorage;

  /// geometry transformer used to recontruct the present data
  AliMUONGeometryTransformer* fOldGeoTransformer;

  /// new geometry transformer containing the new alignment to be applied
  AliMUONGeometryTransformer* fNewGeoTransformer;

  Bool_t fLoadOCDBOnce;
  Bool_t fOCDBLoaded;

  /// Total number of track read
  Int_t fTrackTot;

  /// Number of tracks used for alignment
  Int_t fTrackOk;

  /// Last run number
  Int_t fLastRunNumber;

  /// Array of alignment parameters
  Double_t fParameters[4*156];

  /// Array of alignment parameters errors
  Double_t fErrors[4*156];

  /// Array of alignment parameters pulls
  Double_t fPulls[4*156];

  /// Graph of translations along x
  TGraphErrors *fMSDEx ;

  /// Graph of translations along y
  TGraphErrors *fMSDEy ;

  /// Graph of translations along z
  TGraphErrors *fMSDEz ;

  /// Graph of rotation about z
  TGraphErrors *fMSDEp;

  /// list of graphs
  TList *fList;

  /// list of track records
  TClonesArray *fRecords;

  /// number of records
  Int_t fRecordCount;

  ClassDef(AliMUONAlignmentTask, 3)

};

#endif
