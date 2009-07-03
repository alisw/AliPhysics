#ifndef ALIMUONCHAMBERCALIBRATIONTASK_H
#define ALIMUONCHAMBERCALIBRATIONTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calign
/// \class  AliMUONChamberCalibrationTask
/// \brief  Definition of the task to extract cluster information 
///  from MCH tracks after applying the calibration on aligned ESD
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

class AliMUONChamberCalibrationTask : public AliAnalysisTaskSE
{

 public:

  enum Calibration_t { NOGAIN, GAINCONSTANTCAPA, GAIN, INJECTIONGAIN };
  typedef Calibration_t Calibration_t;

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
  TTree* OutputTree() { return fClusterInfoTree; }

  UInt_t BuildClusterMap( AliMUONTrack &track );

 private:

  char* fOCDBPath; //!< default path to the condition database
  Calibration_t fCalibChoice; //!< calibration option
  TTree* fClusterInfoTree; //!< TTree filled with the cluster information
  AliMUONRecoParam* fMuonRecoParam; //!< reconstruction parameters for track refitting
  AliMUONClusterInfo* fClusterInfo; //!< the cluster info used to fill the output TTree
  AliMUONCalibrationData* fCalibData; 
  AliMUONESDInterface* fESDInterface;
  AliMUONVDigitStore* fDigitStore;
  AliESDInputHandler* fESDInputHandler; //!< ESD input handler
  AliESDEvent* fESDInputEvent; //!< pointer to the current input ESD event

  ClassDef( AliMUONChamberCalibrationTask, 1 ) //Task to extract cluster information after applying calibration

};


#endif
