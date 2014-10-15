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

/* $Id$*/

///
/// \file   AliMUONVRawStreamTracker.cxx
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   28-11-2007
/// \brief  Implementation of the constructors and destructors for AliMUONVRawStreamTracker.
///

//-----------------------------------------------------------------------------
/// \ingroup raw
/// \class AliMUONVRawStreamTracker
/// \brief This class is the base class for raw stream decoders than need to deal with
/// raw data coming from the muon tracking chambers.
///
/// The classes that derive from this abstract class should loops over all MUON
/// digits in the raw data given by the AliRawReader.
/// The Next methods should be overridden so that they step through the all the
/// digits and return kFALSE or zero when done. kTRUE or the number of digits
/// decoded should be returned if any digits were actually found.
///
/// \author Artur Szostak <artursz@iafrica.com>
//-----------------------------------------------------------------------------

#include "AliMUONVRawStreamTracker.h"

/// \cond CLASSIMP
ClassImp(AliMUONVRawStreamTracker)
/// \endcond

const Int_t AliMUONVRawStreamTracker::fgkMaxDDL = 20;


AliMUONVRawStreamTracker::AliMUONVRawStreamTracker() : AliMUONRawStream()
{
	///
	/// Default constructor.
	///
}


AliMUONVRawStreamTracker::AliMUONVRawStreamTracker(AliRawReader* rawReader) :
	AliMUONRawStream(rawReader)
{
	///
	/// Constructor with AliRawReader as argument.
	///
}


AliMUONVRawStreamTracker::~AliMUONVRawStreamTracker()
{
	///
	/// Default destructor.
	///
}
