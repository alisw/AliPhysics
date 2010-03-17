//-*- Mode: C++ -*-
// $Id$

#ifndef ALI_HLT_PREDICTION_PROCESSOR_INTERFACE_H
#define ALI_HLT_PREDICTION_PROCESSOR_INTERFACE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTPredictionProcessorInterface.h
    @author Sebastian Bablok
    @date   
    @brief  
*/

#include <AliPreprocessor.h>

class AliHLTPendolino;


/**
 * The Class defines the interface HLT Preprocessor.
 * It inherits from the AliPreprocessorInterface and adds the features for
 * the prediction procession required by the HLT Pendolino.
 *
 * @author Sebastian Bablok
 *
 * @date 2007-10-22
 */
class AliHLTPredictionProcessorInterface : public AliPreprocessor {
    public:
	
		/**
		 * Constructor for AliHLTPredictionProcessorInterface
		 *
		 * @param detector string defining the detector to which the 
		 * 			PredictionProcessor belongs to
		 * @param pendolino pointer to the hosting pendolino (derived from 
		 * 			AliShuttleInterface)
		 */
		AliHLTPredictionProcessorInterface(const char* detector, 
						AliHLTPendolino* pendolino);

		/**
		 * Destructor for AliHLTPredictionProcessorInterface
		 */
		virtual ~AliHLTPredictionProcessorInterface();
		
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
		virtual UInt_t makePrediction(Bool_t doPrediction = true) = 0;

		/**
		 * Function to let the PredictionProcessor produce dummy input data, 
		 * that can be used in Pendolino tests, where no DCS Archieve DB is 
		 * contacted. This function is called by the Pendolino, the result is
		 * given back to the PredictionProcessor via the Process(...) call. 
		 * Since the DCSMaps requested from DCS are detector specific, the 
		 * PredictionProcessor should know, how the maps look like, that it
		 * expects. This function is pure virtual, because it has to be 
		 * implemented by the real PredictionProcessors. 
		 * NOTE: the clean - up (delete) of the TMap will be performed by the
		 * Pendolino after the usage. The PredictionProcessor has never to 
		 * call delete on the returned TMap pointer.
		 *
		 * @param aliasName optional parameter, that can be given in order to 
		 * 			create a DCSMap for dedicated aliases. For a general test
		 * 			this parameter should be empty. If more than one alias name
		 * 			shall be specified, they are separated by blanks " ". 
		 *
		 * @return DCSMap containing dummy data for testing the Pendolino.
		 */
		virtual TMap* produceTestData(TString aliasName = "") = 0;

		/**
		 * Virtual function, implemented in the detector specific 
		 * PredictionProcessors to initialize them.
		 *
		 * @param run run number
		 * @param startTime start time of data
		 * @param endTime end time of data
		 */
//		virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

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
//		virtual UInt_t Process(TMap* dcsAliasMap) = 0;   

		/**
		 * Indicates if DCS data shall be processed.
		 * NOTE: should always return true, since it is used as prediction 
		 * processor, which will only process DCS data
		 *
		 * @return true if DCS data can be processed, else false. Note: if false
		 * 				the Pendolino would stop, so make sure that it is true.
		 */
//		virtual Bool_t ProcessDCS();
		
	protected:
		/**
		 * Helper function to receive the current run number
		 *
		 * @return the current run number
		 */
	   virtual Int_t GetRunNumber();	

       /**
        * Function to let the Penolino add an request for an AliCDBEntry to the
        * according Taxi list. This corresponding object will be request then
        * with the next run of the Taxi. The entry is added to the Pendolino
        * specific list for the Taxi configuration. Note: the name has to
        * include the detector name as well; e.g.: "TPC/Calib/LocalVdrift" .
        *
        * @param entryPath the full path name of the entry, that shall be
        *               included
        *
        * @return  true, when successful included or entry already existing in
        *               list; else false.
        */
       virtual Bool_t includeAliCDBEntryInList(const TString& entryPath);


	private:
		/**
		 * Disabled Standart Constructor
		 */
	   AliHLTPredictionProcessorInterface();
	   
		/**
		 * Disabled Copy constructor 
		 * (parent class is disabled so derived class does the same)
		 */
		AliHLTPredictionProcessorInterface(
						const AliHLTPredictionProcessorInterface& predictPro);

		/**
		 * Disabled Assignment operator
		 * (parent class is disabled so derived class does the same)
		 */
		AliHLTPredictionProcessorInterface& operator=(
						const AliHLTPredictionProcessorInterface& rhs);

		/**
		 * Stores pointer to Pendolino
		 */
		AliHLTPendolino* fpPend;  //  Stores pointer to Pendolino

		
		ClassDef(AliHLTPredictionProcessorInterface, 6);
	
};



#endif


