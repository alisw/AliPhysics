//-*- Mode: C++ -*-
// $Id$

#ifndef ALI_HLT_PREDICT_PROC_TEMP_MONITOR_H
#define ALI_HLT_PREDICT_PROC_TEMP_MONITOR_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTPredicProcTempMonitor.h
    @author Sebastian Bablok
    @date   
    @brief  
*/

#include "AliHLTPredictionProcessorInterface.h"

#include <TObjString.h>


/**
 * Implementation of an AliHLTPredictionProcessor for monitoring DCS 
 * Temperatur values inside the HLT. 
 *
 * @author Sebastian Bablok
 * @author Gaute Oevrebekk (extention for B-Field)
 *
 * @date 2007-12-03
 * @date 2008-03-04
 */
class AliHLTPredicProcTempMonitor : public AliHLTPredictionProcessorInterface {
    public:

		/** Static const defining path part 2 */
		static const TString kPath2; // defining path part 2
		/** Static const defining path part 3 */
		static const TString kPath3; //defining path part 3
		/** Static const defining creator of HCDB object */
		static const TString kCreator; // defining creator of HCDB object
		/** Static const defining comment for HCDB object */
		static const TString kComment; // defining comment for HCDB object
		/** Static const defining aliroot version */
		static const TString kAliRootVersion; // defining aliroot version
		/** Static const defining error ret val for storage error */
		static const UInt_t kUnableToStoreObject; // defining error ret val for storage error
//		static const TString kAmandaTempSensor;
	
		/**
		 * Constructor for AliHLTPredictionProcessorDummy
		 *
		 * @param detector string defining the detector to which the 
		 * 			PredictionProcessor belongs to
		 * @param pendolino pointer to the hosting pendolino (derived from 
		 * 			AliShuttleInterface)
		 */
		AliHLTPredicProcTempMonitor(const char* detector, 
						AliHLTPendolino* pendolino);

		/**
		 * Destructor for AliHLTPredictionProcessorDummy
		 */
		virtual ~AliHLTPredicProcTempMonitor();
		
		/**
		 * Pure virtual function to force the Prediction Processor to implement
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
		 * Pure Virtual function, implemented in the detector specific 
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


	protected:
		
		/**
		 * Function to extract the B-Field from the DCS map.
		 *
		 * @param dcsAliasMap the DCS map, containing the AliDCSValues
		 *
		 * @return 0 in case of success, else an error code
		 */
		UInt_t ExtractBField(TMap* dcsAliasMap);

		/**
		 * Function to extract a float AliDCSValue by its name from the DCS map
		 *
		 * @param dcsAliasMap the DCS value map
		 * @param stringId the name of the desired AliDCSValue - float
		 * @param value pointer to a float, where the extratced value shall be stored
		 * 				(out parameter)
		 * 
		 * @return true in case of success, else false
		 */
		Bool_t GetSensorValue(TMap* dcsAliasMap, const char* stringId, 
						Float_t* value);
		
	private:
		/**
		 * Disabled Copy constructor 
		 * (parent class is disabled so derived class does the same)
		 */
		AliHLTPredicProcTempMonitor(
						const AliHLTPredicProcTempMonitor& predictPro);

		/**
		 * Disabled Assignment operator
		 * (parent class is disabled so derived class does the same)
		 */
		AliHLTPredicProcTempMonitor& operator=(
						const AliHLTPredicProcTempMonitor& rhs);


		/**
		 * Function to extract the actual values from the DCSMap
		 *
		 * @param inputMap pointer to the DCSMap as input
		 * @param outputMap pointer to the Map where to store the extracted 
		 * 			values
		 *
		 * @return true in case of success, else false
		 */
//		Bool_t extractTempvalues(TMap* inputMap, TMap* outputMap);

		/**
		 * Stores if prediction shall be made
		 */
		Bool_t fPredict;  // flag for prediction making

		/**
		 * Stores the run number
		 */
		Int_t fRun;  //  Stores the run number
	
	 	/**
		 * Stores the start time of the to process DCS data
		 */	
		UInt_t fStartTime; // Stores the start time of the to process DCS data
	
	 	/**
		 * Stores the end time of the to process DCS data
		 */	
		UInt_t fEndTime;  // Stores the end time of the to process DCS data

		/**
		 * Stores the extracted B-Field in a TObjString
		 */
		TObjString fBField;  // Stores the extracted B-Field in a TObjString

		/** ClassDef of AliHLTPredicProcTempMonitor for AliRoot */
		ClassDef(AliHLTPredicProcTempMonitor, 4);
	
};


inline Bool_t AliHLTPredicProcTempMonitor::ProcessDCS() {
	// indicates if DCS values are processed; always true
	return true;
}

#endif


