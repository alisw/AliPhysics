#ifndef ALIMUONCHAMBERCALIBRATIONTASK_H
#define ALIMUONCHAMBERCALIBRATIONTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calign
/// \class  AliMUONChamberCalibrationTask
/// \brief  Definition of the task to extract cluster information 
///  from MCH tracks after applying the calibration on an aligned ESD
/// \author Andry Rakotozafindrabe CEA/IRFU/SPhN

class AliMUONRecoParam;
class AliMUONClusterInfo;
class AliMUONPadInfo;
class AliMUONCalibrationData;
class AliMUONESDInterface;
class AliMUONVClusterStore;
class AliMUONVDigitStore;
class AliMUONTrack;

class AliESDInputHandler;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"
#include "TTree.h" 
#include "TString.h" 

class AliMUONChamberCalibrationTask : public AliAnalysisTaskSE
{

 public:

  /// enumerate the list of the available modes for the calibration
  enum Calibration_t { 
    kNOGAIN, 
    kGAINCONSTANTCAPA, 
    kGAIN, 
    kINJECTIONGAIN 
  }; 

  // Constructors
  AliMUONChamberCalibrationTask();
  AliMUONChamberCalibrationTask( const char* name, 
				 char* ocdbpath, 
				 const Int_t my_calib_option );

  // Destructor
  virtual ~AliMUONChamberCalibrationTask();

  // Implementation of interface methods
  virtual void CreateOutputObjects(); 
  virtual void LocalInit(); 
  virtual void ConnectInputData( Option_t* option = "" );
  virtual void Exec( Option_t* option  = "" ); 
  virtual void Terminate( Option_t* option = "" ); 

  // Getters
  /// Return TTree filled with the cluster information
  TTree* OutputTree() const { return fClusterInfoTree; }

  UInt_t BuildClusterMap( AliMUONTrack &track );

 private:
  /// Not implemented
  AliMUONChamberCalibrationTask(const AliMUONChamberCalibrationTask& right);
  /// Not implemented
  AliMUONChamberCalibrationTask&  operator = (const AliMUONChamberCalibrationTask& right);

  TString fOCDBPath;                    //!< default path to the condition database
  Calibration_t fCalibChoice;           //!< calibration option
  TTree* fClusterInfoTree;              //!< TTree filled with the cluster information
  AliMUONRecoParam* fMuonRecoParam;     //!< reconstruction parameters for track refitting
  AliMUONClusterInfo* fClusterInfo;     //!< cluster info used to fill the output TTree
  AliMUONCalibrationData* fCalibData;   //!< needed to access to the calibration data for each pad within each cluster
  AliMUONESDInterface* fESDInterface;   //!< interface to easily access to the ESD content
  AliMUONVDigitStore* fDigitStore;      //!< pointer to the digit stored for the current input ESD event 
  AliESDInputHandler* fESDInputHandler; //!< ESD input handler
  AliESDEvent* fESDInputEvent;          //!< pointer to the current input ESD event

  ClassDef( AliMUONChamberCalibrationTask, 1 ) // Task to extract cluster information after applying calibration

};


#endif
