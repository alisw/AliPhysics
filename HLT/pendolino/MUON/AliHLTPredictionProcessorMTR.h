#ifndef ALIHLTPREDICTIONPROCESSORMTR_H
#define ALIHLTPREDICTIONPROCESSORMTR_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTPredictionProcessorMTR.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   2007-11-24
/// @brief  Declaration of the muon trigger prediction processor.
///

#include "AliHLTPredictionProcessorInterface.h"

/**
 * This is the prediction processor for the muon trigger detector.
 * It is required as part of the Pendolino interface so that data fetched from
 * DCS can be transformed and extrapolated into a format suitable for the HCDB.
 *
 * At the moment this class does not really do anything, but will fill the
 * HCDB with data in the future.
 *
 * Refer to:
 * http://wiki.kip.uni-heidelberg.de/ti/HLT/index.php/Pendolino-PredictionProcessor-Interface
 * for more details about the Pendolino and preprocessor interfaces.
 */
class AliHLTPredictionProcessorMTR : public AliHLTPredictionProcessorInterface
{
public:

	/**
	 * The class constructor for the preprocessor.
	 * @param detector  The detector name, which must be "MTR".
	 * @param pendolino A valid pointer to the pendolino instance.
	 */
	AliHLTPredictionProcessorMTR(const char* detector, AliHLTPendolino* pendolino);
	
	virtual ~AliHLTPredictionProcessorMTR();
	
	// The following methods are inherited from AliHLTPredictionProcessorInterface.
	
	/**
	 * This function is called by the Pendolino before the fetched DCS data
	 * is handed in for (prediction) processing.
	 *
	 * @param doPrediction If true then preprocessed data will be extrapolated
	 *           to the current time. Otherwise no extrapolation will be
	 *           performed if false. The default is to perform extrapolations.
	 *
	 * @return Currently always returns 0 for success.
	 */
	virtual UInt_t makePrediction(Bool_t doPrediction = true);

	/**
	 * Performs initialisation of the preprocessor.
	 * At the moment just basic sanity checks are performed.
	 *
	 * @param run The current run number.
	 * @param startTime The start time (earliest timestamp) of the data.
	 * @param endTime The end time (latest timestamp) of the data.
	 */
	virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

	/**
	 * Process the DCS values.
	 * At the moment nothing is done here.
	 *
	 * @param dcsAliasMap The map containing aliases and corresponding DCS
	 *                    values and timestamps
	 *
	 * @return At the moment 0 is always returned for success.
	 */
	virtual UInt_t Process(TMap* dcsAliasMap);

	/**
	 * Indicates that the DCS data shall be processed.
	 * NOTE: should always return true, since it is used as prediction
	 * processor, which will only process DCS data
	 *
	 * @return Always returns the value kTRUE.
	 */
	virtual Bool_t ProcessDCS() { return kTRUE; };

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
	
private:

	// Do not allow copying of this class.
	AliHLTPredictionProcessorMTR(const AliHLTPredictionProcessorMTR& /*obj*/);
	AliHLTPredictionProcessorMTR& operator = (const AliHLTPredictionProcessorMTR& /*obj*/);

	Bool_t fPredict;  //! Flag indicating if the processed values should be extrapolated to current time.
	
	ClassDef(AliHLTPredictionProcessorMTR, 0);  // The prediction processor for the muon trigger detector.
};

#endif // ALIHLTPREDICTIONPROCESSORMTR_H
