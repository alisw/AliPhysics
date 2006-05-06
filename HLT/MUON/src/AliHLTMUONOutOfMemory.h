#ifndef ALIHLTMUONOUTOFMEMORY_H
#define ALIHLTMUONOUTOFMEMORY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
// AliHLTMUONOutOfMemory is defined to be used when the system runs
// out of memory. Do not throw this object directly but rather use the routine
// ThrowAliHLTMUONOutOfMemory which throws a pree allocated static object.
//
////////////////////////////////////////////////////////////////////////////////

#include "AliHLTMUONError.h"


// Error code declaration.
enum
{
	kAliHLTMUONOutOfMemory = 0x10000001
};

class AliHLTMUONOutOfMemory;

class AliHLTMUONOutOfMemory : public AliHLTMUONError
{
public:
	virtual const char* Message() const throw ();
	virtual Int ErrorCode() const throw ();
};

/* When one needs to indicate that no more memory is available one should use
   the ThrowAliHLTMUONOutOfMemory routine rather than explicitly using the code
       throw AliHLTMUONOutOfMemory();
   This is because the Throw methods throws a preallocated object so we are
   safe from trying to allocate any more unavailable memory.
 */
void ThrowAliHLTMUONOutOfMemory() throw (AliHLTMUONOutOfMemory);


/* Overload the global new and delete operators to perform better error checking.
   These new operators throw AliHLTMUONOutOfMemory exceptions if the system runs
   out of memory.
 */
void* operator new (size_t size) throw (AliHLTMUONOutOfMemory);
void* operator new [] (size_t size) throw (AliHLTMUONOutOfMemory);
void operator delete (void* memory) throw ();
void operator delete [] (void* memory) throw ();

#endif // ALIHLTMUONOUTOFMEMORY_H
