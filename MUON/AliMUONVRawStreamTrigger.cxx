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

/// \class AliMUONVRawStreamTrigger
///
/// This base class defines what interface all MUON trigger raw data decoders
/// must implement. Thus all trigger decoders inherit from this class.
/// This class is an abstract class.
///
/// \author Artur Szostak <artursz@iafrica.com>

#include "AliMUONVRawStreamTrigger.h"
#include "AliRawReader.h"
#include <cassert>

/// \cond CLASSIMP
ClassImp(AliMUONVRawStreamTrigger)
/// \endcond

//___________________________________________
AliMUONVRawStreamTrigger::AliMUONVRawStreamTrigger()
	: AliMUONRawStream()
{
	///
	/// Create a base object to read MUON raw triggers
	/// Default constructor for monitoring purposes.
	///
}

//_________________________________________________________________
AliMUONVRawStreamTrigger::AliMUONVRawStreamTrigger(AliRawReader* rawReader)
	: AliMUONRawStream(rawReader)
{
	///
	/// Constructor with AliRawReader as argument
	/// for reconstruction purpose.
	///
}

//___________________________________
AliMUONVRawStreamTrigger::~AliMUONVRawStreamTrigger()
{
	///
	/// Default destructor
	///
}

