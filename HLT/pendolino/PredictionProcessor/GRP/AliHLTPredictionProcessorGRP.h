//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPREDICTIONPROCESSORGRP_H
#define ALIHLTPREDICTIONPROCESSORGRP_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTPredictionProcessorGRP.cxx
//  @author Matthias Richter
//  @date   2010-04-12
//  @brief  Prediction processor for the GRP entry
// 

#include "AliHLTPredictionProcessorInterface.h"

/**
 * Predition Processor for the GRP
 *
 */
class AliHLTPredictionProcessorGRP : public AliHLTPredictionProcessorInterface
{
public:
	
  /**
   * Constructor for AliHLTPredictionProcessorGRP
   *
   * @param pendolino pointer to the hosting pendolino
   */
  AliHLTPredictionProcessorGRP(AliHLTPendolino* pendolino);

  /// Destructor for AliHLTPredictionProcessorGRP
  virtual ~AliHLTPredictionProcessorGRP();
		
  /**
   * Virtual function to force the Prediction Processor to implement
   * a function to flag that prediction making is required.
   * This function is called by the Pendolino before the fetched DCS data
   * is handed in for (prediction) processing.
   *
   * @param doPrediction if true, prediction making shall be switched on,
   * 			if false, switched off.
   * 
   * @return 0 on success; a value greater than 0 refers to an error
   */
  virtual UInt_t makePrediction(Bool_t doPrediction = true);

  /**
   * Virtual function, implemented in the detector specific 
   * PredictionProcessors to initialize them.
   *
   * @param run run number
   * @param startTime start time of data
   * @param endTime end time of data
   */
  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

  /**
   * Function called by the Pendolino for each participating subdetector
   * producing the required condition settings. This includes the 
   * encoding of a prediction to the values due to the fact, that only 
   * data up to now can be fetched from DCS (the latest data can also be 
   * up to 2 min old until it is received on HLT side). To write the data
   * to the HCDB, the detector specific implementation of this class has 
   * to call the appropriated storing function provided by the interface.
   *
   * @param dcsAliasMap the map containing aliases and corresponding DCS
   * 			values and timestamps
   *
   * @return 0 on success; a value greater than 0 refers to an error
   */
  virtual UInt_t Process(TMap* dcsAliasMap);

  /**
   * Indicates if DCS data shall be processed.
   * NOTE: should always return true, since it is used as prediction 
   * processor, which will only process DCS data
   *
   * @return true if DCS data can be processed, else false. Note: if false
   * 				the Pendolino would stop, so make sure that it is true.
   */
  virtual Bool_t ProcessDCS();

  /**
   * Function to let the PredictionProcessor produce dummy input data,
   * that can be used in Pendolino tests, where no DCS Archieve DB is
   * contacted. This function is called by the Pendolino, the result is
   * given back to the PredictionProcessor via the Process(...) call.
   * Since the DCSMaps requested from DCS are detector specific, the
   * PredictionProcessor should know, how the maps look like, that it
   * expects.
   * NOTE: The clean - up (delete) of the TMap will be performed by the
   * Pendolino after the usage. The PredictionProcessor has never to
   * call delete on the returned TMap pointer.
   *
   * @param aliasName optional parameter, that can be given in order to
   *          create a DCSMap for dedicated aliases. For a general test
   *          this paramter should be empty. If more than one alias name
   *          shall be specified, they are separated by blanks " ".
   *
   * @return DCSMap containing dummy data for testing the Pendolino.
   */
  virtual TMap* produceTestData(TString aliasName = "");

  /**
   * Create the initial GRP entry.
   * The beam type and run type parameters are propagated form the ECS to
   * the run manager, and provided as arguments to the pendolino.
   * The initial GRP entry is created ignoring the magnetic field, this
   * is updated by the pendolino afterwords.
   */
  static bool CreateInitialGRPEntry(int runno, const char* beamtype,
				    const char* runtype, const char* detectorList);

protected:

private:
  /**
   * Disabled Copy constructor 
   * (parent class is disabled so derived class does the same)
   */
  AliHLTPredictionProcessorGRP(const AliHLTPredictionProcessorGRP& predictPro);

  /**
   * Disabled Assignment operator
   * (parent class is disabled so derived class does the same)
   */
  AliHLTPredictionProcessorGRP& operator=(const AliHLTPredictionProcessorGRP& rhs);
		
  ClassDef(AliHLTPredictionProcessorGRP, 1);
	
};


inline Bool_t AliHLTPredictionProcessorGRP::ProcessDCS() {
	// indicates if DCS values are processed; always true
	return true;
}

#endif


