//-*- Mode: C++ -*-
// $Id$

#ifndef ALI_HLT_PREDICTION_PROCESSOR_HLT_H
#define ALI_HLT_PREDICTION_PROCESSOR_HLT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTPredictionProcessorHLT.h
    @author Gaute Ovrebekk
    @date   
    @brief  
*/

#include "AliHLTPredictionProcessorInterface.h"

//new
#include <TObjString.h>


/**
 * Predition Processor for prepering b-field values for general HLT components. 
 *
 * @author Sebastian Bablok, Gaute Ovrebekk
 *
 * @date 2007-10-24
 */
class AliHLTPredictionProcessorHLT : public AliHLTPredictionProcessorInterface {
    public:
	
		/**
		 * Constructor for AliHLTPredictionProcessorHLT
		 *
		 * @param detector string defining the detector to which the 
		 * 			PredictionProcessor belongs to
		 * @param pendolino pointer to the hosting pendolino (derived from 
		 * 			AliShuttleInterface)
		 */
		AliHLTPredictionProcessorHLT(const char* detector, 
						AliHLTPendolino* pendolino);

		/**
		 * Destructor for AliHLTPredictionProcessorHLT
		 */
		virtual ~AliHLTPredictionProcessorHLT();
		
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


        protected:

	private:
		/**
		 * Disabled Copy constructor 
		 * (parent class is disabled so derived class does the same)
		 */
		AliHLTPredictionProcessorHLT(
						const AliHLTPredictionProcessorHLT& predictPro);

		/**
		 * Disabled Assignment operator
		 * (parent class is disabled so derived class does the same)
		 */
		AliHLTPredictionProcessorHLT& operator=(
						const AliHLTPredictionProcessorHLT& rhs);

		/**
		 * Stores if prediction shall be made
		 */
		Bool_t fPredict;  // flag for prediction making

		/**
		 * Stores the run number
		 */
		Int_t fRun;  // Stores the run number
	
	 	/**
		 * Stores the start time of the to process DCS data
		 */	
		UInt_t fStartTime;  // Stores the start time of the to process DCS data
	
	 	/**
		 * Stores the end time of the to process DCS data
		 */	
		UInt_t fEndTime;  // Stores the end time of the to process DCS data

        /**
		 * TObjstring which stores the B-field value
		 */
		TObjString fBField;  // TObjstring which stores the B-field value

		/**
		 * Function to extract the B-field from the DCS value map
		 *
		 * @param dcsAliasMap the retrieved DCS value map
		 *
		 * @return 0 on sucess else an error code
		 */
		UInt_t ExtractBField(TMap* dcsAliasMap);

		/**
		 * Function to rertieve a sensor value from the DCS value map
		 *
		 * @param dcsAliasMap the retrieved DCS value map
		 * @param stringId the alias name of the desired sensor value
		 * @param value [return parameter] - the extracted sensor value
		 *
		 * @return true if sucessful, else false
		 */
		Bool_t GetSensorValue(TMap* dcsAliasMap,const char* stringId, Float_t * value);
		
		ClassDef(AliHLTPredictionProcessorHLT, 1);
	
};


inline Bool_t AliHLTPredictionProcessorHLT::ProcessDCS() {
	// indicates if DCS values are processed; always true
	return true;
}

#endif


