//-*- Mode: C++ -*-
// $Id$

#ifndef ALI_HLT_PENDOLINO_H
#define ALI_HLT_PENDOLINO_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTPendolino.h
//  @author Sebastian Bablok
//  @date   
//  @brief  
//  @note   maintained by Matthias.Richter@ift.uib.no

//#include <TObject.h>
#include <TString.h>
#include <TMap.h>


#include <AliShuttleInterface.h>

//#include "AliHLTPredictionProcessorInterface.h"

#ifdef SHUTTLE_PRE_REV29388_INTERFACE
#define CONST_PROPERTY const
#else
#define CONST_PROPERTY
#endif

class AliHLTPendolinoLogger;

/**
 * Class that implements the AliShuttleInterface and provides the required
 * features for the contacting the PredictionProcessor and storing the result
 * to the HCDB
 * 
 * @author Sebastian Bablok
 *
 * @date 2007-10-22
 */
class AliHLTPendolino : public AliShuttleInterface {
	public:

		/**
 		 * Static string that defines the local storage indicator in path
 		 */  
		static const TString fgkLocalStorageDefine;  // defines the local storage

		/**
		 * Static string defining the name of this inteface module.
		 */
		static const char* fgkHLTInterfaceModule;  // defines the name of inteface module

		/**
 		 * Static value defining error value for a Pendolino exception.
 		 */ 
		static const Int_t fgkHLTPendolinoException;  // error value for a Pendolino exception

		/**
 		 * Static value defining error value for a bad cast
 		 */ 
		static const Int_t fgkHLTPendolinoBadCast;  // error value for a bad cast

		/**
 		 * Static value defining error value for handed in module is not 
 		 * implementing the PredictionProcessor interface.
 		 */ 
		static const Int_t fgkHLTPendolinoNotPredictProc;  //  error value for "not implementing the PredictionProcessor interface"

		/**
		 * Static value defining error value for module not existing.
		 */ 
		static const Int_t fgkHLTPendolinoModuleNotExisting;  // error value for module not existing

		/**
		 * Static value defining error value for PredictionProc does not
		 * process DCS values.
		 */	
		static const Int_t fgkHLTPendolinoNoDCS; // error value for PredictionProc does not process DCS

		/**
 		 * Static string that defines the base folder for the Taxi list files.
 		 */ 
		static const TString fgkTaxiListBaseFolder;  // defines the base folder for the Taxi list files

		/**
		 * Static string that defines list folder name for taxi list
		 */
		static const TString fgkTaxiListFolderName; // defines list folder name for taxi list

		/**
 		 * Static string that defines the filename for the Taxi list required
 		 * by the Pendolino
 		 */  
		static const TString fgkTaxiListPendolino; // defines the filename for the Taxi list 

		/**
 		 * Static value that defines the max length of a line that can be read
 		 * when browsing through the list file
 		 */ 
		static const Int_t fgkMaxLineLength; // defines the max length of a line
		

		/**
		 * Constructor for AliHLTPendolino.
		 *
		 * @param run the current run number
		 * @param HCDBbase path to the HCDB base folder
		 * @param runType the current run type
		 * @param logger pointer to the desired logger; default is the command
		 * 			line logger.
		 * @param startTime start time of the DCS Archive DB request
		 * @param endTime end time of the DCS Archive DB request
		 */
		AliHLTPendolino(Int_t run, TString HCDBbase, TString runType = "TEST", 
				AliHLTPendolinoLogger* logger = 0, UInt_t startTime = 0, 
				UInt_t endTime = 0);

		/**
		 * Destructor for
		 */
		virtual ~AliHLTPendolino();

		/**
		 * Function to store processed data in the HCDB
		 *
		 * @param path the storage path consisting of "detector/type/objectname"
		 * 			where detector should be the detector name used by offline,
		 * 			type is something like "CALIB" or "ALIGN", while 
		 * 			objectname the name represents, by which it can be requested
		 * 			from the HCDB.
		 * @param object the object to store
		 * @param metaData metaData to store
		 * @param validityStart object is valid from current run number minus 
		 * 			validityStart
		 * @param validityInfinite if true, validity is set to infinity
		 *
		 * @return true on success, else false
		 */
		virtual Bool_t Store(const AliCDBPath& path, TObject* object, 
				AliCDBMetaData* metaData, Int_t validityStart = 0, 
				Bool_t validityInfinite = kFALSE);
		/**
		 * Function is required from interface, but is disbled in its 
		 * implementation since additional refernce data will not be used by 
		 * the DAs inside the HLT. Therefore always returns false
		 *
		 * @return always false, since function is disabled.
		 */
	    virtual Bool_t StoreReferenceData(const AliCDBPath& path, 
				TObject* object, AliCDBMetaData* metaData);
		/**
		 * Function is required from interface, but is disbled in its
         * implementation since additional reference data will not be used by
         * the DAs inside the HLT. Therefore always returns false
         *
         * @return always false, since function is disabled.
         */
	    virtual Bool_t StoreReferenceFile(const char* detector, 
				const char* localFile, const char* gridFileName);

		/**
		 * Function is required from interface, but is disbled in its
		 * implementation since additional run meta data will only be stored
		 * by the the GPR - Preprocessor of the Offline Shuttle.
		 * 
		 * @return always false, since function is disabled.
		 */                  
		virtual Bool_t StoreRunMetadataFile(const char* localFile, 
				const char* gridFileName);

		/**
		 * Function is required from interface, but is disbled in its
         * implementation since additional refernce data will not be used by
         * the DAs inside the HLT. Therefore always returns NULL.
		 * If this feature is lateron required by HLT, inherit from this class
		 * and overwrite this function with the required functionality.
         *
         * @return always NULL, since not used.
         */
//	    virtual const char* GetFile(Int_t system, const char* detector, 
//				const char* id, const char* source);  --> is private now
	    
        /**
         * Function is required from interface, but is disbled in its
         * implementation since additional refernce data will not be used by
         * the DAs inside the HLT. Therefore always returns NULL.
         * If this feature is lateron required by HLT, inherit from this class
         * and overwrite this function with the required functionality.
         *
         * @return always NULL, since not used.
         */
//	    virtual TList* GetFileSources(Int_t system, const char* detector, 
//				const char* id = 0);  -> is priavte now
		
        /**
         * Function is required from interface, but is disbled in its
         * implementation since additional refernce data will not be used by
         * the DAs inside the HLT. Therefore always returns NULL.
         * If this feature is lateron required by HLT, inherit from this class
         * and overwrite this function with the required functionality.
         *
         * @return always NULL, since not used.
         */
//	    virtual TList* GetFileIDs(Int_t system, const char* detector, 
//				const char* source);  -> is private now
	    
		/**
		 * Retrieves current run parameter.
		 * 
		 * @param lbEntry name of the run parameter
		 *
		 * @return value of the run parameter
		 */
	    virtual const char* GetRunParameter(const char* lbEntry);

		/**
		 * Retrieves the current run type.		
		 *
		 * @return the current run type
		 */
		virtual const char* GetRunType();
		
		/**
		 * Returns the HLT status.
		 * Since the Pendolino is part of the HLT, the function will always
		 * return true
		 *
		 * @return always true - see above
		 */
	    virtual Bool_t GetHLTStatus();

		/**
		 * Retrieves a file from the OCDB, but does not uses the OCDB but the 
		 * HCDB (local copy of the OCDB). - Since the HCDB only contains OCDB
		 * objects, that are included in the list the Taxi uses for fetching 
		 * the data from the OCDB, it can happen, that the corresponding 
		 * AliCDBEntry is missing. 
		 * TODO: Think of mechanism to automatically include the object in the
		 * Taxi list.
		 * 
		 * @param detector the detectorname, to which the object belongs
		 * @param path the path of the object
		 *
		 * @return pointer to the fetched HCDB entry
		 */
	    virtual AliCDBEntry* GetFromOCDB(const char* detector, 
				const AliCDBPath& path);
	    
		/**
		 * Function to allow Pendolino and PredictionProcessor to make log 
		 * entries.
		 *
		 * @param detector the detector, that wants to issue this message
		 * @param message the log message
		 */
	    virtual void Log(const char* detector, const char* message);

		/**
		 * Function is required from interface, but is disbled in the Pendolino
		 * because at the moment there is no Trigger configuration available 
		 * for Pendolino. At the moment it just return a NULL pointer and makes
		 * some log output.
		 *
		 * @return NULL pointer, since it should not be used.
		 */
//		virtual const char* GetTriggerConfiguration(); --> is private now

  virtual const char* GetCTPTimeParams() {return "";}
  virtual void Log(const char* name, const char* message, UInt_t level);

		/**
		 * Registers a preprocessor; actually it has to be a PredictionProcessor
		 * since the Pendolino requires a PredictionProcessor. If the registered
		 * preprocessor is not implementing a PredictionProcessor, a log entry
		 * is made and the registration discarded.
		 *
		 * @param preprocessor the PredictionProcessor that shall be registered.
		 * 			For more details please see above !!
		 */ 
		virtual void RegisterPreprocessor(AliPreprocessor* preprocessor);
//					AliHLTPredictionProcessorInterface* preprocessor);

		/**
		 * Function to get the current run number
		 *
		 * @return current run number
		 */
		virtual Int_t GetRunNumber();

		/**
		 * Function to enable prediction making in all registered 
		 * PredictionProcessors, if they implement the required interface
		 *          
		 * @return number of PredictionProcessor, which has been switched to
		 * 				prediction making (NOTE: the internal list may conatin 
		 * 				more PredictionProcessors, but if switching on failed 
		 * 				for one, this one is not counted.)
		 */
		virtual UInt_t SetToPredictMaking(); 

		/**
		 * Function to get the number of registered PredictionProcessors
		 *
		 * @return number of registered PredictionProcessors
		 */
		Int_t GetNumberOfPredictProc();

		/**
		 * Function to check if given PredtionProc allows for processing DCS
		 * and enable prediction making
		 *
		 * @param detector the detector in whose PredictionProcessor the 
		 * 			prediction making shall be enabled.
		 *
		 * @return 0 on success, else an error code is returned
		 */
		virtual Int_t setToPredictMaking(TString detector);

		/**
		 * Function to initlaize a dedicated Prediction Processor
		 * 
		 * @param detector the detector, whose PredictProc shall be initialized
		 * @param run the current run number
		 * @param startTime the start time of the fetched data
		 * @param endTime the end time of the fetched data
		 *
		 * @return 0 on success, else an error code is returned
		 */
		virtual Int_t initPredictProc(TString detector, Int_t run, 
					UInt_t startTime, UInt_t endTime);

		/**
		 * Function to hand in retrieved DCS values. These values are handed to
		 * the corresponding PredictionProcessor of the according detector.
		 * The PredictionProcessor should prepare the data inside and store 
		 * them to the HCDB.
		 *
		 * @param detector the according detector, to whose PredictionProcessor
		 * 			shall prepare the handed-in data
		 * @param DCSValues pointer to the map containing the fetched DCS values
		 *
		 * @return 0 on success, else an error code is returned.
		 */
		virtual Int_t PrepareDCSValues(TString detector, TMap* DCSValues);

		/**
		 * Function to retrieve dummy data for testing the Pendolino from a
		 * given PredictionProcessor. The function called is handed further to
		 * the corresponding PredictionProcessor.
		 * NOTE: The returned TMap can be NULL, if no corresponding 
		 * PredictionProcessor is registered.
		 *
		 * @param detector the according detector, from whom the 
		 * 			PredictionProcessor shall produce the dummy data.
		 * @param aliasName optional parameter to hand in a alias name for 
		 * 			producing a DCSMap for the given alias.
		 *
		 * @return the DCSMap with the dummy data to test (given by the 
		 * 			PredictionProcessor). NOTE: can be NULL, if no corresponding
		 * 			PredictionProcessor is registered.
		 */
		virtual TMap* EmulateDCSMap(TString detector, TString aliasName = "");
   
		/**
 		 * Function to add a entry request to the Taxi lists.
 		 *
 		 * @param entryPath the path entry, that shall be included in the 
 		 * 				list file.
 		 *
 		 * @return true, when successful included or entry already existing in 
 		 * 				list; else false.
 		 */
		virtual Bool_t IncludeAliCDBEntryInList(const TString& entryPath); 

		/**
		 * Function to get the start time of the DCS Archive DB request; in HLT
		 * this is the same like the start time given in the Initialize() call 
		 * to the PredictionProcessors (NOTE: thsi is differnet to the 
		 * implementation in the Offline Shuttle - there the initial start time
		 * is set to the start-of-data for the complete run.)
		 */
#ifdef SHUTTLE_PRE_REV29388_INTERFACE
		virtual const UInt_t GetStartTimeDCSQuery();
#else
		virtual UInt_t GetStartTimeDCSQuery();
#endif

		/**
		 * Function to get the end time of the DCS Archive DB request; in HLT
		 * this is the same like the end time given in the Initialize() call
		 * to the PredictionProcessors (NOTE: thsi is differnet to the
		 * implementation in the Offline Shuttle - there the initial end time
		 * is set to the end-of-data for the complete run.)
		 */
#ifdef SHUTTLE_PRE_REV29388_INTERFACE
		virtual const UInt_t GetEndTimeDCSQuery();
#else
		virtual UInt_t GetEndTimeDCSQuery();
#endif
				
  /**
   * method introduced as pure virtual in r43691
   * needs to be implemented to create pendolino instances
   */
  virtual void SendMLFromDet(const char* /*value*/) {}
  /**
   * method introduced as pure virtual in r45587
   * needs to be implemented to create pendolino instances
   */
  virtual TString* GetLTUConfig(const char* /*det*/) {return NULL;}	
	protected:

		
	private:
        /**
         * Function is required from interface, but is disbled in its
         * implementation since additional refernce data will not be used by
         * the DAs inside the HLT. Therefore always returns NULL.
         * If this feature is lateron required by HLT, inherit from this class
         * and overwrite this function with the required functionality.
         *
         * @return always NULL, since not used.
         */
        virtual TList* GetFileIDs(Int_t system, const char* detector,
                const char* source);

        /**
         * Function is required from interface, but is disbled in its
         * implementation since additional refernce data will not be used by
         * the DAs inside the HLT. Therefore always returns NULL.
         * If this feature is lateron required by HLT, inherit from this class
         * and overwrite this function with the required functionality.
         *
         * @return always NULL, since not used.
         */
        virtual TList* GetFileSources(Int_t system, const char* detector,
                const char* id = 0);

        /**
         * Function is required from interface, but is disbled in its
         * implementation since additional refernce data will not be used by
         * the DAs inside the HLT. Therefore always returns NULL.
         * If this feature is lateron required by HLT, inherit from this class
         * and overwrite this function with the required functionality.
         *
         * @return always NULL, since not used.
         */
        virtual const char* GetFile(Int_t system, const char* detector,
                const char* id, const char* source);

        /**
         * Function is required from interface, but is disbled in the Pendolino
         * because at the moment there is no Trigger configuration available
         * for Pendolino. At the moment it just return a NULL pointer and makes
         * some log output.
         *
         * @return NULL pointer, since it should not be used.
         */
		virtual const char* GetTriggerConfiguration();
                virtual const char* GetTriggerDetectorMask();

		/**
		 * Disable the default constructor.
		 */
		AliHLTPendolino();

		/** copy constructor prohibited */
		AliHLTPendolino(const AliHLTPendolino&);
		/** assignment operator prohibited */
		AliHLTPendolino& operator=(const AliHLTPendolino&);

		/**
		 * Stores the current run type
		 */
		TString fRunType;  // Stores the current run type

		/**
 		 * Stores the current run number
 		 */
		Int_t fRunNumber;  // Stores the current run number

		/**
 		 * Stores the HCDBpath
 		 */
		TString fHCDBPath;	 // Stores the HCDBpath

		/**
		 * Map that stores the all PredictionProcessors with their name
		 * (detector)
		 */
		TMap fPredictionProcessorMap;  // stores the all PredictionProcessors

		/**
		 * Pointer to the used Pendolino logger
		 */
		AliHLTPendolinoLogger* fpLogger; // Pointer to the used Pendolino logger

		/**
		 * Indicates, if Logger is owned by Pendolino
		 */
		Bool_t fOwnLogger;  //  Indicates, if Logger is owned by Pendolino

		/**
		 * Stores the start time of the DCS Archive DB request
		 */
                UInt_t fStartTime; //!

		/**
		 * Stores the end time of the DCS Archive DB request
		 */
                UInt_t fEndTime; //!
		
		ClassDef(AliHLTPendolino, 6);

};


inline const char* AliHLTPendolino::GetRunType() {
	// getter for run type
	return fRunType.Data();
}

inline Int_t AliHLTPendolino::GetRunNumber() {
	// getter for run number
	return fRunNumber;
}

inline Int_t AliHLTPendolino::GetNumberOfPredictProc() {
	// getter for number of registered PredictionProcessors
	return fPredictionProcessorMap.GetSize();
}


#endif

