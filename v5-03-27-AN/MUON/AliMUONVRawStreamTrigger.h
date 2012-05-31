#ifndef ALIMUONVRAWSTREAMTRIGGER_H
#define ALIMUONVRAWSTREAMTRIGGER_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup raw
/// \class AliMUONVRawStreamTrigger
/// \brief Base class for reading MUON raw trigger data.
///
//  MUON trigger decoders must derive from this base class.
//
//  Author Artur Szostak <artursz@iafrica.com>

#include <TObject.h>
#include "AliMUONRawStream.h"

class TArrayS;
class AliRawReader;
class AliMUONDDLTrigger;

class AliMUONVRawStreamTrigger : public AliMUONRawStream
{
public:
	AliMUONVRawStreamTrigger();
	AliMUONVRawStreamTrigger(AliRawReader* rawReader);
	virtual ~AliMUONVRawStreamTrigger();
	
	/// Advance one step in the iteration. Returns false if finished.
	virtual Bool_t Next(UChar_t& id,   UChar_t& dec,     Bool_t& trigY,
	                    UChar_t& yPos, UChar_t& sXDev,   UChar_t& xDev,
	                    UChar_t& xPos, Bool_t& triggerY, Bool_t& triggerX,
	                    TArrayS& xPattern, TArrayS& yPattern) = 0;
	
	/// Return pointer to DDL payload object.
	virtual AliMUONDDLTrigger* GetDDLTrigger() const = 0;
	
	/// Return maximum number of DDLs
	virtual Int_t GetMaxDDL() const = 0;
	/// Return maximum number of regional cards in DATE file
	virtual Int_t GetMaxReg() const = 0;
	/// Return maximum number of local cards in DATE file
	virtual Int_t GetMaxLoc() const = 0;
	
	/// Should set the maximum number of local cards expected in the DDL stream.
	virtual void SetMaxLoc(Int_t loc) = 0;
	
	/// Return number of DDL
	virtual Int_t GetDDL() const = 0;
	
	/// Disable Warnings
	virtual void DisableWarnings() = 0;
	
	/// error numbers
	enum rawStreamTriggerError
	{
		kDarcEoWErr   = 6, ///< end of Darc word error 
		kGlobalEoWErr = 7, ///< end of Global word error
		kRegEoWErr    = 8, ///< end of Regional word error 
		kLocalEoWErr  = 9  ///< end of local word error
	};

private:
	/// Not implemented
	AliMUONVRawStreamTrigger(const AliMUONVRawStreamTrigger& stream);
	/// Not implemented
	AliMUONVRawStreamTrigger& operator = (const AliMUONVRawStreamTrigger& stream);
	
	ClassDef(AliMUONVRawStreamTrigger, 0)  // Base class for MUON trigger rawdata decoders.
};

#endif // ALIMUONVRAWSTREAMTRIGGER_H
