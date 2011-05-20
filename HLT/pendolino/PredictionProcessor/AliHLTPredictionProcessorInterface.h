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
#include "AliDCSSensor.h"
#include "AliDCSValue.h"
#include "AliSplineFit.h"
#include "TMap.h"
#include <memory>

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

  /**
  *  Create an AliDCSSensor object with specified id and a linear const spline fit
  *  There is no extrapolation in the online reconstruction, only the last value
  *  counts. Some offline code requires the spline fit to be initialized, so it is
  *  added here even though it does not fully make sense.
  */ 
  template <typename T>
  static AliDCSSensor* CreateSensor(const char* id, T value, UInt_t starttime, UInt_t endtime);

  UInt_t GetStartTime() const {return fStartTime;}
  UInt_t GetEndTime() const {return fEndTime;}

	protected:
		/**
		 * Helper function to receive the current run number
		 *
		 * @return the current run number
		 */
	   virtual Int_t GetRunNumber();	
  
  Bool_t IsPredicting() const {return fPredict;}
  UInt_t StartTime() const {return fStartTime;}  
  UInt_t EndTime() const {return fEndTime;}  

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

  /**
   * Function to retrieve a sensor value from the DCS value map
   * The value is extracted from the last entry, because online reconstruction
   * does not need the history.
   *
   * @param dcsAliasMap the retrieved DCS value map
   * @param stringId the alias name of the desired sensor value
   * @param value [return parameter] - the extracted sensor value
   *
   * @return true if sucessful, else false
   */
  template<typename T>
  Bool_t GetSensorValue(TMap* dcsAliasMap,const char* stringId, T * value) const;

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

		/**
		 * Stores if prediction shall be made
		 */
		Bool_t fPredict;  // flag for prediction making

	 	/**
		 * Stores the start time of the to process DCS data
		 */	
		UInt_t fStartTime;  // Stores the start time of the to process DCS data
	
	 	/**
		 * Stores the end time of the to process DCS data
		 */	
		UInt_t fEndTime;  // Stores the end time of the to process DCS data
		
		ClassDef(AliHLTPredictionProcessorInterface, 7);
	
};

template<typename T>
Bool_t AliHLTPredictionProcessorInterface::GetSensorValue(TMap* dcsAliasMap,
							  const char* stringId, T *value) const
{
  // extracts the sensor value
  // return last value read from sensor specified by stringId
  
  TObject* object=dcsAliasMap->FindObject(stringId);
  if (!object) return kFALSE;
  TPair* pair = dynamic_cast<TPair*>(object);
  if (pair && pair->Value()) {
    TObjArray* valueSet = dynamic_cast<TObjArray*>(pair->Value());
    Int_t nentriesDCS = valueSet->GetEntriesFast() - 1;
    if(nentriesDCS>=0 && valueSet->At(nentriesDCS)){
      AliDCSValue *val = dynamic_cast<AliDCSValue *>(valueSet->At(nentriesDCS));
      if (val) {
	*value=*val;
	return kTRUE;
      }
    }
  }
  return kFALSE;
}

template <typename T>
AliDCSSensor* AliHLTPredictionProcessorInterface::CreateSensor(const char* id, T value, UInt_t starttime, UInt_t endtime) 
{
  // create an AliDCSSensor object with specified id and a linear graph
  // There is no extrapolation in the online reconstruction, only the last value
  // counts

  if (!id) return NULL;

  std::auto_ptr<AliDCSSensor> pSensor(new AliDCSSensor);
  if (!pSensor.get()) {
    //HLTFatal("memory allocation failed");
    return NULL;
  }

  // AliDCSSensor allows two types of value representation: a spline fit and
  // a graph. The online system uses a linear graph with const values between
  // start and end time
  // Note: AliDCSSensor::GetValue returns -99 if the requested time is before
  // the start time and the last value if the time is after end time
  // The measurements are stored in fractions of hours (see definition in
  // class AliDCSSensor
  const int points=2;
  const Double_t kSecInHour = 3600.; // seconds in one hour
  T x[points];
  T y[points];
  x[0]=0;
  x[1]=(endtime-starttime)/kSecInHour;
  y[0]=value;
  y[1]=value;
  std::auto_ptr<TGraph> pGraph(new TGraph(2,x,y));
  if (!pGraph.get()) {
    //HLTFatal("can not create graph for id %s", id);
    return NULL;
  }

  AliSplineFit *fit = new AliSplineFit();
  if (!fit) return NULL;

  fit->SetMinPoints(10);
  fit->InitKnots(pGraph.get(),10, 10, 0.0);
  fit->SplineFit(2);

  pSensor->SetStringID(id);
  pSensor->SetStartTime(starttime);
  pSensor->SetEndTime(endtime);
  // note: AliDCSSensor has no correct cleanup in the destructor
  // so the fit object is lost if the sensor is deleted
  pSensor->SetFit(fit);

  return pSensor.release();
}

#endif
