/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///
///  @file   AliHLTPredictionProcessorMTR.cxx
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   
///  @brief  Implementation of the AliHLTPredictionProcessorMTR class.
///

#include "AliHLTPredictionProcessorMTR.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include <cstdlib>

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <AliDCSValue.h>

ClassImp(AliHLTPredictionProcessorMTR);


AliHLTPredictionProcessorMTR::AliHLTPredictionProcessorMTR(
		const char* detector, AliHLTPendolino* pendolino
	) :
	AliHLTPredictionProcessorInterface(detector, pendolino),
	fPredict(kTRUE)
{
	/// The class constructor for the preprocessor.
	/// @param detector  The detector name, which must be "MTR".
	/// @param pendolino A valid pointer to the pendolino instance.
	
	if (strcmp(detector, "MTR") != 0)
	{
		Log(Form("Warning: Setting the detector name to an incorrect value = '%s'."
		    " It must be set to 'MTR'.", detector)
		);
	}
}


AliHLTPredictionProcessorMTR::~AliHLTPredictionProcessorMTR()
{
	/// Default destructor.
}


UInt_t AliHLTPredictionProcessorMTR::makePrediction(Bool_t doPrediction)
{
	/// This function is called by the Pendolino before the fetched DCS data
	/// is handed in for (prediction) processing.
	///
	/// @param doPrediction If true then preprocessed data will be extrapolated
	///           to the current time. Otherwise no extrapolation will be
	///           performed if false. The default is to perform extrapolations.
	///
	/// @return Currently always returns 0 for success.

	fPredict = doPrediction;
	return 0;
}


void AliHLTPredictionProcessorMTR::Initialize(
		Int_t run, UInt_t startTime, UInt_t endTime
	)
{
	/// Performs initialisation of the preprocessor.
	/// At the moment just basic sanity checks are performed.
	///
	/// @param run The current run number.
	/// @param startTime The start time (earliest timestamp) of the data.
	/// @param endTime The end time (latest timestamp) of the data.
	
	if (startTime > endTime)
		Log("Error: start time is greater than end time.");
	
	AliPreprocessor::Initialize(run, startTime, endTime);

	if (fPredict)
		Log("Prediction is switched ON.");
	else
		Log("Prediction is switched OFF.");
}


UInt_t AliHLTPredictionProcessorMTR::Process(TMap* /*dcsAliasMap*/)
{
	/// Process the DCS values.
	/// At the moment nothing is done here.
	///
	/// @param dcsAliasMap The map containing aliases and corresponding DCS
	///                    values and timestamps
	///
	/// @return At the moment 0 is always returned for success.
	
	Log("Processing MTR");
	return 0;
}

TMap* AliHLTPredictionProcessorMTR::produceTestData(TString /*aliasName*/) {
    TMap* resultMap = 0;

    // here has to come real dummy data :-)
    resultMap = new TMap();
    TTimeStamp tt;
	Float_t fval = 33.3;
    TObjString* name = new TObjString("DummyData");
    AliDCSValue* val = new AliDCSValue(fval, tt.GetTime());
    TObjArray* arr = new TObjArray();
    arr->Add(val);
    resultMap->Add(name, arr);

    return resultMap;
}


