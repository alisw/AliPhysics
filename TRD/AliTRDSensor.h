#ifndef ALITRDSENSOR_H
#define ALITRDSENSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TRD sensor for DCS                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include "AliDCSSensor.h"

//_____________________________________________________________________________
class AliTRDSensor : public  AliDCSSensor {

	public :
		AliTRDSensor ();
		~AliTRDSensor();
		AliTRDSensor (const AliTRDSensor & source);
		AliTRDSensor & operator=(const AliTRDSensor &c);
		void Copy(TObject &c) const;
		AliTRDSensor (Int_t dcsId,
					  Double_t x, Double_t y, Double_t z); 
		
				
	

		ClassDef(AliTRDSensor,2);   // Object for TRD DCS values
};

#endif

