#ifndef ALIMUONTRACKERQADATAMAKERREC_H
#define ALIMUONTRACKERQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id: AliMUONTrackerQADataMakerRec.h 35760 2009-10-21 21:45:42Z ivana $

/// \ingroup rec
/// \class AliMUONTrackerQADataMakerRec
/// \brief MUON Quality assurance data maker
///

// --- AliRoot header files ---
#include "AliMUONVQADataMakerRec.h"
#include "AliMUONRecoParam.h"
#include "AliQAv1.h"

class AliMUONDigitMaker;
class AliMUONVClusterStore;
class AliMUONVDigitStore;
class AliMUONVStore;
class AliMUONVTrackerData;
class AliMUONVTrackerDataMaker;
class AliMUONCalibrationData;
class AliMUONQAMappingCheck;
class AliMUONLogger;
class AliMUONQADataMakerRec;
class AliRawVEvent;

class AliMUONTrackerQADataMakerRec: public AliMUONVQADataMakerRec {

public:
  AliMUONTrackerQADataMakerRec(AliQADataMakerRec* master);         
  virtual ~AliMUONTrackerQADataMakerRec();
  
  AliMUONVTrackerData* GetTrackerData() const;

  virtual void InitDigits(); 
  virtual void InitESDs(); 
  virtual void InitRaws(); 
  virtual void InitRecPoints(); 
  
  void EndOfDetectorCycleRaws(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleDigits(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleRecPoints(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleESDs(Int_t specie, TObjArray** list);

  virtual void MakeDigits(TTree* dig); 
  virtual void MakeESDs(AliESDEvent* esd) ;
  virtual void MakeRaws(AliRawReader* rawReader); 
  virtual void MakeRecPoints(TTree* recpo); 
  
  void ResetDetectorRaws(TObjArray* list);
  
private:

  AliMUONQADataMakerRec* Master() const;
  
  void BookHistograms(AliQAv1::TASKINDEX_t task);

  void FillReadoutStatus(AliMUONLogger& log, AliMUONVTrackerData* data);
  
  void FillEventSize(const AliRawVEvent* event);
  
  void InitCommon();

  void InsertTrackerData(Int_t specie, TObjArray** list, TObject* object, 
                         Int_t indexNumber, 
                         Bool_t replace=kFALSE);
  
  void ProjectTrackerData(AliMUONVTrackerData* data, 
                          TH1& hbp,
                          TH1& hnevents,
                          TH1& hddl,
                          TH1& hddlevents);

  AliMUONVTrackerDataMaker* TrackerDataMaker(Int_t specie) const;

  AliMUONVTrackerDataMaker* TrackerDataMaker(Int_t specie, Bool_t create);
  
  AliMUONQAMappingCheck* MappingCheckRecPoints(Int_t specie, Bool_t create=kFALSE);
  
  AliMUONVTrackerData* TrackerCalData(Int_t specie, Bool_t create=kFALSE);
  
  AliMUONVTrackerData* TrackerRecData(Int_t specie, Bool_t create=kFALSE);
  
  TObjArray* GetArray(TObjArray*& array, Bool_t create);

private:
  /// Not implemented
  AliMUONTrackerQADataMakerRec(const AliMUONTrackerQADataMakerRec& rhs);
  /// Not implemented
  AliMUONTrackerQADataMakerRec& operator=(const AliMUONTrackerQADataMakerRec& rhs);
  
  AliMUONVDigitStore*   fDigitStore; //!< pointer to digits store
  AliMUONDigitMaker*    fDigitMaker;  //!< pointer to digit maker
  AliMUONVClusterStore* fClusterStore; //!< pointer to cluster store
	
  AliMUONCalibrationData* fCalibrationData; //!< Used to load Local, Regional and Global masks
  
  AliMUONLogger* fLogger; //!< (readout) error logger
  
  TH1* fBusPatchConfig; //!< bus patch configuration
  
  Double_t fBPxmin; //!< min bin value for bus patch
  Double_t fBPxmax; //!< max bin value for bus patch
  Int_t fBPnbins; //!< number of bus patch bins

  TObjArray* fTrackerDataMakerArray; //!< tracker data accumulation (Raws)  
  TObjArray* fTrackerCalDataArray; //!< tracker data accumulation (calibrated digits)  
  TObjArray* fTrackerRecDataArray; //!< tracker data accumulation (only calibrated digits belonging to reconstructed clusters)
  TObjArray* fMappingCheckRecPointsArray; //!< mapping cross-checker (RecPoints)  
  
  ClassDef(AliMUONTrackerQADataMakerRec,5)  // MUON Quality assurance data maker

};
#endif
