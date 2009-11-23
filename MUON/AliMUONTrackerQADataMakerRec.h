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

class AliMUONDigitMaker;
class AliMUONVClusterStore;
class AliMUONVDigitStore;
class AliMUONVStore;
class AliMUONVTrackerData;
class AliMUONVTrackerDataMaker;
class AliMUONCalibrationData;
class AliMUONQAMappingCheck;

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
  void EndOfDetectorCycleRecPoints(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleESDs(Int_t specie, TObjArray** list);

  /// Empty implementation 
  void EndOfDetectorCycleDigits(Int_t, TObjArray**) {}
    
  virtual void MakeDigits(TTree* dig); 
  virtual void MakeESDs(AliESDEvent* esd) ;
  virtual void MakeRaws(AliRawReader* rawReader); 
  virtual void MakeRecPoints(TTree* recpo); 
  
private:
  
  void InsertTrackerData(Int_t specie, TObjArray** list, TObject* object, 
                         Int_t indexNumber, Bool_t replace=kFALSE);

private:
  /// Not implemented
  AliMUONTrackerQADataMakerRec(const AliMUONTrackerQADataMakerRec& rhs);
  /// Not implemented
  AliMUONTrackerQADataMakerRec& operator=(const AliMUONTrackerQADataMakerRec& rhs);
  
  AliMUONVDigitStore*   fDigitStore; //!< pointer to digits store
  AliMUONDigitMaker*    fDigitMaker;  //!< pointer to digit maker
  AliMUONVClusterStore* fClusterStore; //!< pointer to cluster store
	
  AliMUONVTrackerDataMaker* fTrackerDataMaker; //!< tracker data accumulation (Raw)
  
  AliMUONQAMappingCheck* fMappingCheckRecPoints; //!< mapping cross-checker (RecPoints)
  
  AliMUONCalibrationData* fCalibrationData; //!< Used to load Local, Regional and Global masks
  
  ClassDef(AliMUONTrackerQADataMakerRec,1)  // MUON Quality assurance data maker

};
#endif
