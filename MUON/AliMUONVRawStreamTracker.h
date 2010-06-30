#ifndef ALIMUONVRAWSTREAMTRACKER_H
#define ALIMUONVRAWSTREAMTRACKER_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$*/

///
/// \file   AliMUONVRawStreamTracker.h
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   28-11-2007
/// \brief  Declaration of the abstract base class for muon trigger chamber raw stream decoders.
///

#include "AliMUONRawStream.h"

class AliMUONDDLTracker;

class AliMUONVRawStreamTracker : public AliMUONRawStream
{
public:

	/// Default constructor.
	AliMUONVRawStreamTracker();
	
	/// Constructor setting the raw reader.
	AliMUONVRawStreamTracker(AliRawReader* rawReader);
	
	/// Default destructor.
	virtual ~AliMUONVRawStreamTracker();
	
	/// Advance one step in the iteration. Returns false if finished.
	virtual Bool_t Next(Int_t& busPatchId,
				UShort_t& manuId, UChar_t& manuChannel,
				UShort_t& adc) = 0;
	
	/// Return pointer to DDL payload object.
	virtual AliMUONDDLTracker* GetDDLTracker() const = 0;
	
	/// Return maximum number of DDLs
	static Int_t GetMaxDDL() { return fgkMaxDDL; };
	
	/// Return maximum number of blocks per DDL allowed.
	virtual Int_t GetMaxBlock() const = 0;
	/// Return maximum number of Dsp per block allowed.
	virtual Int_t GetMaxDsp()   const = 0;
	/// Return maximum number of Buspatch per Dsp allowed.
	virtual Int_t GetMaxBus()   const = 0;
	
	/// Set maximum number of blocks per DDL allowed.
	virtual void SetMaxBlock(Int_t blk) = 0;
	/// Set maximum number of Dsp per block allowed.
	virtual void SetMaxDsp(Int_t dsp) = 0;
	/// Set maximum number of Buspatch per Dsp allowed.
	virtual void SetMaxBus(Int_t bus) = 0;
	
	/// Return number of the current DDL.
	virtual Int_t GetDDL() const = 0;
	
	/// check error/Warning presence
	virtual Bool_t IsErrorMessage() const = 0;
	
	/// error numbers
	enum rawStreamTrackerError
	{
		kGlitchErr      = 1, ///< glitch error
		kPaddingWordErr = 2, ///< padding word error
		kParityErr      = 3  ///< parity error
	};

protected:	

	static const Int_t fgkMaxDDL;   //!< maximum number of DDLs

private:

	// Do not allow copying of this class.
        /// Not implemented
	AliMUONVRawStreamTracker(const AliMUONVRawStreamTracker& stream);
        /// Not implemented
	AliMUONVRawStreamTracker& operator = (const AliMUONVRawStreamTracker& stream);

	ClassDef(AliMUONVRawStreamTracker, 0) // Base class for reading MUON raw digits from tracking chambers.
};

#endif  // ALIMUONVRAWSTREAMTRACKER_H
