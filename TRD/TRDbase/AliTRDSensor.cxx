/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// This class represent a TRD sensor, it is used by AliTRDSensorArray     //
// 
//                                                                        //
// Author:                                                                //
//   W. Monange   (w.monange@gsi.de)                                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDSensor.h"

ClassImp(AliTRDSensor)

//_____________________________________________________________________________
AliTRDSensor::AliTRDSensor ():
		AliDCSSensor ()
{
	//default constructor
	SetIdDCS  (0);
	SetX       (0);
	SetY       (0);
	SetZ       (0);
	//SetStartTime (startTime);
}


//_____________________________________________________________________________
AliTRDSensor::AliTRDSensor (const AliTRDSensor & source) :
		AliDCSSensor (source)
{ 
	//copy constructor

}

//_____________________________________________________________________________
AliTRDSensor::AliTRDSensor (Int_t dcsId,
							Double_t x, Double_t y, Double_t z):
		AliDCSSensor ()
{
	// constructor
	SetIdDCS  (dcsId);
	SetX       (x);
	SetY       (y);
	SetZ       (z);
	//SetStartTime (startTime);
}

//_____________________________________________________________________________
AliTRDSensor::~AliTRDSensor()
{
  //
  // Destructor
  //

}

//_____________________________________________________________________________
AliTRDSensor & AliTRDSensor::operator=(const AliTRDSensor &c)
{
  //
  // Assignment operator
  //

	if (this != &c) ((AliTRDSensor &) c).Copy(*this);
	return *this;
}

//_____________________________________________________________________________
void AliTRDSensor::Copy(TObject &c) const
{
  //
  // Copy function
  //
	TObject::Copy(c);
}
