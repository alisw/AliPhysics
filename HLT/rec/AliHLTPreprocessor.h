//-*- Mode: C++ -*-
// @(#) $Id$

#ifndef ALI_HLT_PREPROCESSOR_H
#define ALI_HLT_PREPROCESSOR_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/**
 * @file   AliHLTPreprocessor.h
 * @author Kenneth Aamodt, Sebastian Bablok
 * @date   2007-12-06
 * @brief  Declaration of the HLT preprocessor (used by the Offline Shuttle) 
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliPreprocessor.h"

/**
 * @class AliHLTPreprocessor
 * Implementation of the HLT version for the Shuttle Preprocessor. 
 *
 * @author Sebastian Bablok, Kenneth Aamodt
 *
 * @date 2007-12-05
 */
class AliHLTPreprocessor : public AliPreprocessor {
    public:
	
		/**
		 * Constructor for AliHLTPreprocessor
		 *
		 * @param shuttle pointer to the hosting shuttle
		 */
		AliHLTPreprocessor(AliShuttleInterface* shuttle);

		/**
		 * Destructor for AliHLTPreprocessor
		 */
		virtual ~AliHLTPreprocessor();
		
		/**
		 * Function to initilaize the Preprocessor.
		 *
		 * @param run run number
		 * @param startTime start time of data
		 * @param endTime end time of data
		 */
		virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

		/**
 		 * Function to process data. Inside the preparation and storing to OCDB
 		 * should be handled.
		 *
		 * @param dcsAliasMap the map containing aliases and corresponding DCS
		 * 			values and timestamps
		 *
		 * @return 0 on success; a value greater than 0 refers to an error
		 */
		virtual UInt_t Process(TMap* dcsAliasMap);

		/**
		 * Indicates if DCS data can be processed.
		 *
		 * @return true if DCS data can be processed, else false. 
		 */
		virtual Bool_t ProcessDCS();

		/** Define for HuffmanTable number */
		static const Int_t fgkHuffmanTablesNum;			// see above

		/** Define for HuffmanFileBase */
		static const char* fgkHuffmanFileBase;			// see above

		/** Define for Detector used for Huffman table */
		static const char* fgkHuffmanFileDetector; 		// see above

		/** Define for Temperature Histogram filename */
		static const char* fgkTempHistoFileName; 		// see above

		/** Define for name of the HLT Preproc */
		static const char* fgkHLTPreproc; 				// see above
		
	protected:

	private:
		/**
		 * Disabled Copy constructor 
		 * (parent class is disabled so derived class does the same)
		 *
		 * @param preproc would be the AliHLTPreproc object to make copy of
		 */
		AliHLTPreprocessor(const AliHLTPreprocessor& preproc);

		/**
		 * Disabled Assignment operator
		 * (parent class is disabled so derived class does the same)
		 *
		 * @param rhs the AliHLTPreproc to assign from
		 *
		 * @return reference to assinged AliHLTPreproc
		 */
		AliHLTPreprocessor& operator=(const AliHLTPreprocessor& rhs);

		/**
		 * Stores the run number
		 */
		Int_t fRun;										// see above
	
	 	/**
		 * Stores the start time 
		 */	
		UInt_t fStartTime;								// see above
	
	 	/**
		 * Stores the end time 
		 */	
		UInt_t fEndTime;								// see above

		
		ClassDef(AliHLTPreprocessor, 2);
	
};


#endif


