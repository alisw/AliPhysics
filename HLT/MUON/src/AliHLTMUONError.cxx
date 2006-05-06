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

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
// AliHLTMUONCoreError is the base excpetion class used by the dHLT subsystem.
// All child classes used to throw exception should be derived from this
// class to allow easy catching of classes of errors.
// 
// AliHLTMUONCoreOutOfMemory is also defined to be used when the system runs
// out of memory. Do not throw this object directly but rather use
// AliHLTMUONCoreThrowOutOfMemory which throws a pree allocated static object.
//
////////////////////////////////////////////////////////////////////////////////
 
#include "AliHLTMUONError.h"

namespace
{
	// The one and only pree allocated out of memory error object.
	static AliHLTMUONCoreOutOfMemory gAliOutOfMemObject;

} // end of namespace


const char* AliHLTMUONCoreOutOfMemory::Message() const throw()
{
	return "Out of memory.";
}

Int AliHLTMUONCoreOutOfMemory::ErrorCode() const throw()
{
	return kOutOfMemory;
}

void AliHLTMUONCoreThrowOutOfMemory() throw (AliHLTMUONCoreOutOfMemory)
{
	throw gAliOutOfMemObject;
}

