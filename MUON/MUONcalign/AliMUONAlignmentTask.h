#ifndef ALIMUONALIGNMENTTASK_H
#define ALIMUONALIGNMENTTASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calign
/// \class AliMUONAlignmentTask
/// \brief Task to align the muon spectrometer
///
/// \author Javier Castillo, CEA/Saclay - Irfu/SPhN
/// \author Hugo Pereira Da Costa, CEA/Saclay - Irfu/SPhN

#include "AliAnalysisTaskSE.h"
#include "AliMUONAlignment.h"

class TString;
class TClonesArray;
class AliMUONGeometryTransformer;

class AliMUONAlignmentTrackParam: public TObject
{
  public:

  /// constructor
  AliMUONAlignmentTrackParam(void):
    fPx(0), fPy(0), fPz(0), fHitPattern(0)
    {}

  /// constructor
  AliMUONAlignmentTrackParam( Double_t px, Double_t py, Double_t pz ):
    fPx(px), fPy(py), fPz(pz), fHitPattern(0)
    {}

  /// object name
  const char* GetName( void ) const
  { return "AliMUONAlignmentTrackParam"; }

  /// return number of chambers hit
  Int_t GetNChambers( void ) const;

  /// return number of stations hit
  Int_t GetNStations( void ) const;

  ///@name momentum
  //@{
  Double_t fPx;
  Double_t fPy;
  Double_t fPz;
  //@}

  /// hit pattern
  UInt_t fHitPattern;

  ClassDef(AliMUONAlignmentTrackParam,2)

};

class AliMUONAlignmentTask : public AliAnalysisTaskSE
{

  public:

  /// constructor
  AliMUONAlignmentTask( const char *name = "AliMUONAlignmentTask" );

  /// destructor
  virtual ~AliMUONAlignmentTask();

  /// get pointer to alignment class
  AliMUONAlignment* alignment( void ) const
  { return fAlign; }

  virtual void LocalInit();

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

  /// merge old and new Align CDBs into the new one.
  void SetMergeAlignmentCDBs( Bool_t value )
  { fMergeAlignmentCDBs = value; }

  /// field on alignment
  void SetBFieldOn( Bool_t value )
  {
    fForceBField = kTRUE;
    fBFieldOn = value;
  }

  /// refit straight tracks
  void SetRefitStraightTracks( Bool_t value )
  { fRefitStraightTracks = value; }

  /// do evaluation
  void SetDoEvaluation( Bool_t value )
  { fDoEvaluation = value; }

  /// min Pt
  void SetMinPt( Double_t value )
  { fMinPt = value; }

  /// min PyPz
  /* cutting on this variable allows one to select more and more straight tracks
   *  when magnetic field is on
   */
  void SetMinPyPz( Double_t value )
  { fMinPyPz = value; }

  /// min number of stations
  void SetMinStations( Int_t value )
  { fMinStations = value; }

  /// run range
  void SetRunRange( Int_t runNumberMin, Int_t runNumberMax )
  {
    fRunNumberMin = runNumberMin;
    fRunNumberMax = runNumberMax;
  }

  //@}

  /// output data
  virtual void UserCreateOutputObjects();

  /// per-event method
  virtual void UserExec( Option_t* );
  virtual void NotifyRun();

  /// termination cleanup
  virtual void Terminate( const Option_t* )
  {}

  /// end of task execution
  virtual void FinishTaskOutput();

  /// Set default ocdb
  void SetDefaultStorage( TString defaultOCDB )
  { fDefaultStorage = defaultOCDB; }

  /// Set old (misaligned) alignment path for ocdb
  void SetOldAlignStorage( TString oldalignOCDB )
  { fOldAlignStorage = oldalignOCDB; }

  /// Set new (realigned) alignment path for ocdb
  void SetNewAlignStorage( TString newalignOCDB )
  { fNewAlignStorage = newalignOCDB; }

  /// Flag to set OCDB once at first run notify
  void SetLoadOCDBOnce( Bool_t loadOCDBOnce = kTRUE )
  { fLoadOCDBOnce = loadOCDBOnce; }

  protected:

  /// store misalignment matrices from OCDB into geometry transformer
  void SaveMisAlignmentData( AliMUONGeometryTransformer* ) const;

  private:

  /// copy constructor, not implemented
  AliMUONAlignmentTask(const AliMUONAlignmentTask& obj);

  /// asignment operator, not implemented
  AliMUONAlignmentTask& operator=(const AliMUONAlignmentTask& other);

  ///@name flags
  //@{

  /// perform alignment from alignment records is true. Use ESD tracks otherwise
  Bool_t fReadRecords;

  /// write alignment records to AOD if true
  Bool_t fWriteRecords;

  /// perform alignment (from either tracks or records depending on fReadRecords)
  Bool_t fDoAlignment;

  /// merge old and new Align CDBs into the new one.
  Bool_t fMergeAlignmentCDBs;

  /// true if magnetic field was forced to value, instead of reading from GRP
  Bool_t fForceBField;

  /// Flag for Magnetic field On/Off
  Bool_t fBFieldOn;

  /// Flag for refitting straight tracks
  Bool_t fRefitStraightTracks;

  /// Flag for running refit evaluation
  Bool_t fDoEvaluation;

  //@}

  /// track momentum cut
  Double_t fMinPt;

  /// track momentum cut
  Double_t fMinPyPz;

  /// track min number of stations
  Int_t fMinStations;

  /// The MUON alignment object
  AliMUONAlignment *fAlign;

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

  /// set to true if not willing to re-initialize OCDB at every new run
  Bool_t fLoadOCDBOnce;

  /// set to true when OCDB was loaded at least once
  Bool_t fOCDBLoaded;

  //! event number (internal counter)
  Int_t fEvent;

  /// Total number of track read
  Int_t fTrackTot;

  /// Number of tracks used for alignment
  Int_t fTrackOk;

  /// Number of records (should match fTrackOk)
  Int_t fRecordsTot;

  /// run range
  Int_t fRunNumberMin;
  Int_t fRunNumberMax;

  /// Array of alignment parameters
  Double_t fParameters[AliMUONAlignment::fNGlobal];

  /// Array of alignment parameters errors
  Double_t fErrors[AliMUONAlignment::fNGlobal];

  /// Array of alignment parameters pulls
  Double_t fPulls[AliMUONAlignment::fNGlobal];

  /// list of track parameters
  TClonesArray *fTrackParams;

  /// list of track records
  TClonesArray *fRecords;

  /// number of records
  Int_t fRecordCount;

  ClassDef(AliMUONAlignmentTask, 3)

};

#endif
