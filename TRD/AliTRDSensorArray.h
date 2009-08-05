#ifndef ALITRDSENSORARRAY_H
#define ALITRDSENSORARRAY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TRD sensor array for DCS                                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TString.h>

class TObjArray;
class TClonesArray;
class TMap;

#include "AliDCSSensorArray.h"

//_____________________________________________________________________________
class AliTRDSensorArray : public AliDCSSensorArray
{

	public :

		AliTRDSensorArray();
		AliTRDSensorArray(const char * amanda
				, const char * storeName
                                , Float_t /*diffCut*/
				, TClonesArray * const trdSensor); 
		
		AliTRDSensorArray(const AliTRDSensorArray & source);
		virtual ~AliTRDSensorArray();
		AliTRDSensorArray & operator=(const AliTRDSensorArray & source);
		
		static TObjArray *GetList();
		TString 	  GetStoreName() const           { return fStoreName; }
	  	TMap             *ExtractDCS(TMap * dcsMap);
		void 	     	  SetGraph(TMap * map);
		Int_t 		  GetNGraph( ) const;
		
	private :
		
		TString 	fAmanda;		// amanqda string dcsDatapointAlias
		TString 	fStoreName;		// name used as the name in the storage system
	
		ClassDef(AliTRDSensorArray,2);          // Array of TRD DCS value objects
};

#endif

