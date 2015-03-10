#ifndef ALIMUONRAWSTREAMTRIGGERHP_H
#define ALIMUONRAWSTREAMTRIGGERHP_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup raw
/// \class AliMUONRawStreamTriggerHP
/// \brief Higher performance decoder stream class for reading MUON raw trigger data.
///
//  This provides a streamer interface to the high performance decoder which
//  is required for AliRoot.
//
//  Author Artur Szostak <artursz@iafrica.com>

#include "AliMUONVRawStreamTrigger.h"
#include "AliMUONTriggerDDLDecoder.h"
#include "TArrayS.h"

class AliMUONDDLTrigger;


class AliMUONRawStreamTriggerHP : public AliMUONVRawStreamTrigger
{
public:
	class AliLocalStruct;

	AliMUONRawStreamTriggerHP();
	AliMUONRawStreamTriggerHP(AliRawReader* rawReader);
	virtual ~AliMUONRawStreamTriggerHP();
	
	/// Initialize iterator
	virtual void First();
	
	/// DDL iterator 
	virtual Bool_t NextDDL();
	
	/// Whether the iteration is finished or not
	virtual Bool_t IsDone() const;
	
	/// Nothing is actually done in the AddErrorMessage method because we log
	/// the error messages as we find them in AliDecoderEventHandler::OnError().
	virtual void AddErrorMessage() { };
	
	/// Advance one step in the iteration. Returns false if finished.
	virtual Bool_t Next(UChar_t& id,   UChar_t& dec,     Bool_t& trigY,
	                    UChar_t& yPos, UChar_t& sXDev,   UChar_t& xDev,
	                    UChar_t& xPos, Bool_t& triggerY, Bool_t& triggerX,
	                    TArrayS& xPattern, TArrayS& yPattern);
	
	/// Construct and return a pointer to the DDL payload object.
	virtual AliMUONDDLTrigger* GetDDLTrigger() const;
	
	/// Returns the next local trigger structure.
	const AliLocalStruct* Next();
	
	/// Return maximum number of DDLs.
	virtual Int_t GetMaxDDL() const { return fgkMaxDDL; }
	
	/// Return maximum number of regional cards in the DDL.
	virtual Int_t GetMaxReg() const { return (Int_t) fDecoder.MaxRegionals(); }
	
	/// Set the maximum number of regional cards in the DDL.
	virtual void SetMaxReg(Int_t reg);
	
	/// Return maximum number of local cards in the DDL.
	virtual Int_t GetMaxLoc() const { return (Int_t) fDecoder.MaxLocals(); }
	
	/// Sets the maximum number of local cards in the DDL.
	virtual void SetMaxLoc(Int_t loc);
	
	/// Return number of the current DDL being handled in the range [0..1] and -1 if no DDL set.
	virtual Int_t GetDDL() const { return fDDL - 1; }
	
	/// check error/Warning presence
	virtual Bool_t IsErrorMessage() const { return fHadError; }
	
	/// Set warnings flag to disable warnings on data errors.
	virtual void DisableWarnings() { fDecoder.GetHandler().Warnings(kFALSE); }
	/// Set warnings flag to enable warnings on data errors.
	virtual void EnableWarnings() { fDecoder.GetHandler().Warnings(kTRUE); }
	
	/// Get number of end of DARC word errors in the DDL last decoded.
	UInt_t GetDarcEoWErrors() const {return fDecoder.GetHandler().GetDarcEoWErrors();}
	
	/// Get number of end of Global word errors in the DDL last decoded.
	UInt_t GetGlobalEoWErrors() const {return fDecoder.GetHandler().GetGlobalEoWErrors();}
	
	/// Get number of end of regional word errors in the DDL last decoded.
	UInt_t GetRegEoWErrors() const {return fDecoder.GetHandler().GetRegEoWErrors();}
	
	/// Get number of end of local word errors in the DDL last decoded.
	UInt_t GetLocalEoWErrors() const {return fDecoder.GetHandler().GetLocalEoWErrors();}
	
	/// Number of end of DARC word errors since First() was called.
	UInt_t NumberOfDarcEoWErrors() const { return fTotalNumberOfDarcEoWErrors; }
	
	/// Number of end of global word errors since First() was called.
	UInt_t NumberOfGlobalEoWErrors() const { return fTotalNumberOfGlobalEoWErrors; }
	
	/// Number of end of regional word errors since First() was called.
	UInt_t NumberOfRegEoWErrors() const { return fTotalNumberOfRegEoWErrors; }
	
	/// Number of end of local word errors since First() was called.
	UInt_t NumberOfLocalEoWErrors() const { return fTotalNumberOfLocalEoWErrors; }
	
	/// Whether we got any end of DARC word errors or not since calling First().
	Bool_t HasDarcEoWError() const { return NumberOfDarcEoWErrors() > 0; }
	
	/// Whether we got any end of global word errors or not since calling First().
	Bool_t HasGlobalEoWError() const { return NumberOfGlobalEoWErrors() > 0; }
	
	/// Whether we got any end of regional word errors or not since calling First().
	Bool_t HasRegEoWError() const { return NumberOfRegEoWErrors() > 0; }
	
	/// Whether we got any end of local word errors or not since calling First().
	Bool_t HasLocalEoWError() const { return NumberOfLocalEoWErrors() > 0; }
	
	/// Returns the "try to recover from errors" flag.
	Bool_t TryRecover() const { return Bool_t(fDecoder.TryRecover()); }
	
	/// Sets the "try to recover from errors" flag.
	/// i.e. should the decoder try to recover from errors found in the
	/// payload headers.
	void TryRecover(Bool_t value) { fDecoder.TryRecover(bool(value)); }
	
	/// Light weight interface class to the DARC and global header data.
	class AliHeader
	{
	public:
		/// Default constructor.
		AliHeader(
				UInt_t darcHeader = 0,
				const AliMUONDarcScalarsStruct* darcScalars = NULL,
				const AliMUONGlobalHeaderStruct* globalHeader = NULL,
				const AliMUONGlobalScalarsStruct* globalScalars = NULL
			) :
			fDarcHeader(darcHeader), fDarcScalars(darcScalars),
			fGlobalHeader(globalHeader), fGlobalScalars(globalScalars)
		{
		}
		
		/// Implement shallow copying in the copy constructor.
		AliHeader(const AliHeader& o) :
			fDarcHeader(o.fDarcHeader), fDarcScalars(o.fDarcScalars),
			fGlobalHeader(o.fGlobalHeader), fGlobalScalars(o.fGlobalScalars)
		{
		}
		
		/// Implement shallow copying in the assignment operator.
		AliHeader& operator = (const AliHeader& object)
		{
			memcpy(this, &object, sizeof(AliHeader));
			return *this;
		}
	
		/// Default destructor.
		~AliHeader() {};
		
		/// Return first word
		UInt_t   GetWord()               const {return GetDarcHeader();}
		/// Return global input
		UInt_t   GetGlobalInput(Int_t n) const {return fGlobalHeader->fInput[n];}
		/// Return global output
		UChar_t  GetGlobalOutput() const {return AliMUONTriggerDDLDecoderEventHandler::GetGlobalOutput(fGlobalHeader);}
		/// Return global config
		UShort_t GetGlobalConfig() const {return AliMUONTriggerDDLDecoderEventHandler::GetGlobalConfig(fGlobalHeader);}

		/// Return event type
		UChar_t GetEventType()  const {return AliMUONTriggerDDLDecoderEventHandler::GetDarcEventType(GetDarcHeader());}
		/// Returns true if this was a physics event.
		Bool_t IsPhysicsEvent() const {return GetEventType() == 0x1;}
		/// Return Darc type
		UChar_t GetDarcType()   const {return AliMUONTriggerDDLDecoderEventHandler::GetDarcType(GetDarcHeader());}
		/// Return serial number
		UChar_t GetSerialNb()   const {return AliMUONTriggerDDLDecoderEventHandler::GetDarcSerialNb(GetDarcHeader());}
		/// Return version
		UChar_t GetVersion()    const {return AliMUONTriggerDDLDecoderEventHandler::GetDarcVersion(GetDarcHeader());}
		/// Return VME trig
		Bool_t  GetVMETrig()    const {return AliMUONTriggerDDLDecoderEventHandler::GetDarcVMETrig(GetDarcHeader());}
		/// Return global flag
		Bool_t  GetGlobalFlag() const {return AliMUONTriggerDDLDecoderEventHandler::GetDarcGlobalFlag(GetDarcHeader());}
		/// Return CPT trigger
		Bool_t  GetCTPTrig()    const {return AliMUONTriggerDDLDecoderEventHandler::GetDarcCTPTrig(GetDarcHeader());}
		/// Return DAQ flag
		Bool_t  GetDAQFlag()    const {return AliMUONTriggerDDLDecoderEventHandler::GetDarcDAQFlag(GetDarcHeader());}
		/// Return reg pattern
		UChar_t GetRegPattern() const {return AliMUONTriggerDDLDecoderEventHandler::GetDarcRegPattern(GetDarcHeader());}
	
		// DARC get methods
		/// Return DARC L0 received and used
		UInt_t  GetDarcL0R()   const {return (fDarcScalars != NULL) ? fDarcScalars->fL0R : 0;}
		/// Return DARC L1 physics
		UInt_t  GetDarcL1P()   const {return (fDarcScalars != NULL) ? fDarcScalars->fL1P : 0;}
		/// Return DARC L1 software
		UInt_t  GetDarcL1S()   const {return (fDarcScalars != NULL) ? fDarcScalars->fL1S : 0;}
		/// Return DARC L2 accept
		UInt_t  GetDarcL2A()   const {return (fDarcScalars != NULL) ? fDarcScalars->fL2A : 0;}
		/// Return DARC L2 reject
		UInt_t  GetDarcL2R()   const {return (fDarcScalars != NULL) ? fDarcScalars->fL2R : 0;}
		/// Return DARC clock
		UInt_t  GetDarcClock() const {return (fDarcScalars != NULL) ? fDarcScalars->fClk : 0;}
		/// Return DARC hold (dead time)
		UInt_t  GetDarcHold()  const {return (fDarcScalars != NULL) ? fDarcScalars->fHold : 0;}
		
		// global get methods
		/// Return global L0
		UInt_t  GetGlobalL0()    const {return (fGlobalScalars != NULL) ? fGlobalScalars->fL0 : 0;}
		/// Return global clock
		UInt_t  GetGlobalClock() const {return (fGlobalScalars != NULL) ? fGlobalScalars->fClk : 0;}
		/// Return global scalars or NULL if none exist.
		const UInt_t* GetGlobalScaler()  const {return (fGlobalScalars != NULL) ? &fGlobalScalars->fScaler[0] : NULL;}
		/// Return global hold (dead time)
		UInt_t  GetGlobalHold()  const {return (fGlobalScalars != NULL) ? fGlobalScalars->fHold : 0;}
		/// Return global spare
		UInt_t  GetGlobalSpare() const {return (fGlobalScalars != NULL) ? fGlobalScalars->fSpare : 0;}
		
		/// Return true if type for DARC is default.
		Bool_t DarcIsDefaultType() const {return GetDarcType() == AliMUONTriggerDDLDecoder<AliMUONTriggerDDLDecoderEventHandler>::DarcDefaultType();}
		/// Return true if type for DARC is Vadorh.
		Bool_t DarcIsVadohrType()  const {return GetDarcType() == AliMUONTriggerDDLDecoder<AliMUONTriggerDDLDecoderEventHandler>::DarcVadorhType();}
		
		/// Return the DARC header's raw data.
		UInt_t GetDarcHeader() const {return fDarcHeader;}
		
		/// Return the DARC scalars raw data or NULL if none exist.
		const AliMUONDarcScalarsStruct* GetDarcScalars() const {return fDarcScalars;}
		
		/// Return the global header's raw data.
		const AliMUONGlobalHeaderStruct* GetGlobalHeader() const {return fGlobalHeader;}
		
		/// Return the global scalars raw data or NULL if none exist.
		const AliMUONGlobalScalarsStruct* GetGlobalScalars() const {return fGlobalScalars;}
		
		/// Print the contents of the headers to screen.
		void Print() const;
	
	private:
	
		UInt_t fDarcHeader;  ///< Pointer to DARC header in DDL payload.
		const AliMUONDarcScalarsStruct* fDarcScalars;  ///< Pointer to DARC scalars in DDL payload.
		const AliMUONGlobalHeaderStruct* fGlobalHeader;  ///< Pointer to global header in DDL payload.
		const AliMUONGlobalScalarsStruct* fGlobalScalars;  ///< Pointer to global scalars in DDL payload.
	};
	
	/// Light weight interface class to the regional card header data.
	class AliRegionalHeader
	{
	public:
		/// Default constructor.
		AliRegionalHeader(
				const AliLocalStruct* localsArray = NULL,
				const AliMUONRegionalHeaderStruct* header = NULL,
				const AliMUONRegionalScalarsStruct* scalars = NULL
			) :
			fNext(NULL), fLocalsCount(0), fFirstLocal(localsArray),
			fHeader(header), fScalars(scalars)
		{
		}
		
		/// Implement shallow copying in the copy constructor.
		AliRegionalHeader(const AliRegionalHeader& o) :
			fNext(o.fNext), fLocalsCount(o.fLocalsCount),
			fFirstLocal(o.fFirstLocal), fHeader(o.fHeader),
			fScalars(o.fScalars)
		{
		}
		
		/// Implement shallow copying in the assignment operator.
		AliRegionalHeader& operator = (const AliRegionalHeader& object)
		{
			memcpy(this, &object, sizeof(AliRegionalHeader));
			return *this;
		}
		
		/// Default destructor.
		~AliRegionalHeader() {};
		
		/// Return darc word
		UInt_t   GetDarcWord()     const {return fHeader->fDarcWord;}
		/// Return first reg word
		UInt_t   GetWord()         const {return fHeader->fWord;}
		/// Return regional input
		UInt_t   GetInput(Int_t n) const {assert(n < 2); return fHeader->fInput[n];}
		/// Return L0
		UShort_t GetL0()           const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalL0(fHeader);}
		/// Return mask
		UShort_t GetMask()         const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalMask(fHeader);}
		
		//word: phys type:1, reset: 6, serialNb:5, Id:4, version: 8, regional output:8
		//true for phys, false for soft
		/// Return RegPhysFlag
		Bool_t    GetRegPhysFlag() const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalPhysFlag(fHeader);}
		/// Return ResetNb
		UChar_t   GetResetNb()     const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalResetNb(fHeader);}
		/// Return SerialNb
		UChar_t   GetSerialNb()    const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalSerialNb(fHeader);}
		/// Return Id
		UChar_t   GetId()          const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalId(fHeader);}
		/// Return Version
		UChar_t   GetVersion()     const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalVersion(fHeader);}
		/// Return Output
		UChar_t   GetOutput()      const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalOutput(fHeader);}
		
		//Darc Status: error:10, #fpag:3, MBZ:3, phys type:1, present:1, not_full:1
		// not_empty:1, L2Rej:1, L2Acc:1, L1:1, L0:1, #evt:4, busy:4
		/// Return ErrorBits
		UShort_t GetErrorBits()       const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalErrorBits(fHeader);}
		/// Return FPGANumber
		UChar_t  GetFPGANumber()      const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalFPGANumber(fHeader);}
		/// Return DarcPhysFlag
		Bool_t   GetDarcPhysFlag()    const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalDarcPhysFlag(fHeader);}
		/// Return PresentFlag
		Bool_t   GetPresentFlag()     const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalPresentFlag(fHeader);}
		/// Return RamNotFullFlag
		Bool_t   GetRamNotFullFlag()  const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalRamNotFullFlag(fHeader);}
		/// Return RamNotEmptyFlag
		Bool_t   GetRamNotEmptyFlag() const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalRamNotEmptyFlag(fHeader);}
		/// Return L2RejStatus
		Bool_t   GetL2RejStatus()     const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalL2RejStatus(fHeader);}
		/// Return L2AccStatus
		Bool_t   GetL2AccStatus()     const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalL2AccStatus(fHeader);}
		/// Return L1Status
		Bool_t   GetL1Status()        const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalL1Status(fHeader);}
		/// Return L0Status
		Bool_t   GetL0Status()        const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalL0Status(fHeader);}
		/// Return EventInRam
		UChar_t  GetEventInRam()      const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalEventInRam(fHeader);}
		/// Return Busy
		UChar_t  GetBusy()            const {return AliMUONTriggerDDLDecoderEventHandler::GetRegionalBusy(fHeader);}
		
		// scalar methods
		/// Return regional clock
		UInt_t  GetClock()        const {return (fScalars != NULL) ? fScalars->fClk : 0;}
		/// Return regional ouput
		const UInt_t* GetScaler() const {return (fScalars != NULL) ? &fScalars->fScaler[0] : NULL;}
		/// Return regional hold (dead time)
		UInt_t  GetHold()         const {return (fScalars != NULL) ? fScalars->fHold : 0;}
		
		/// Return raw data of regional header.
		const AliMUONRegionalHeaderStruct* GetHeader() const { return fHeader; }
		
		/// Return the raw data of the regional scalars or NULL if none exist.
		const AliMUONRegionalScalarsStruct* GetScalars() const { return fScalars; }
		
		/// Return the next regional structure header.
		const AliRegionalHeader* Next() const { return fNext; }
		
		/// Returns the first AliLocalStruct class in this regional structure.
		const AliLocalStruct* GetFirstLocalStruct() const { return fFirstLocal; }
		
		/// Returns the number of local trigger structures within this regional structure.
		UInt_t GetLocalStructCount() const { return fLocalsCount; }
	
		/// Return the i'th local trigger structure in this regional structure.
		const AliLocalStruct* GetLocalStruct(UInt_t i) const
		{
			return i < fLocalsCount ? GetFirstLocalStruct() + i : NULL;
		}
	
		/// Sets the next regional structure header.
		void SetNext(const AliRegionalHeader* next) { fNext = next; }

		/// Increments the local trigger structure count.
		void IncLocalStructCount() { fLocalsCount++; };
		
		/// Print the contents of the regional header and scalars to screen.
		void Print() const;
	
	private:
	
		const AliRegionalHeader* fNext;  ///< Pointer to next regional header.
		UInt_t fLocalsCount;  ///< The number of AliLocalStruct objects found in the array pointed to by fFirstLocal.
		const AliLocalStruct* fFirstLocal;  ///< The first local trigger structure of this regional structure.
		const AliMUONRegionalHeaderStruct* fHeader;  ///< Pointer to the regional header in the DDL payload.
		const AliMUONRegionalScalarsStruct* fScalars;  ///< Pointer to the regional scalars in the DDL payload.
	};
	
	/// Light weight interface class to the local trigger card data.
	class AliLocalStruct
	{
	public:
		/// Default constructor.
		AliLocalStruct(
				const AliRegionalHeader* regionalHeader = NULL,
				const AliMUONLocalInfoStruct* localStruct = NULL,
				const AliMUONLocalScalarsStruct* scalars = NULL
			) :
			fRegional(regionalHeader), fNext(NULL),
			fLocalStruct(localStruct), fScalars(scalars),
			fCalculatedId(0)
		{
		}
		
		/// Implement shallow copying in the copy constructor.
		AliLocalStruct(const AliLocalStruct& o) :
			fRegional(o.fRegional), fNext(o.fNext),
			fLocalStruct(o.fLocalStruct), fScalars(o.fScalars),
			fCalculatedId(o.fCalculatedId)
		{
		}
		
		/// Implement shallow copying in the assignment operator.
		AliLocalStruct& operator = (const AliLocalStruct& object)
		{
			memcpy(this, &object, sizeof(AliLocalStruct));
			return *this;
		}
	
		/// Default destructor.
		~AliLocalStruct() {};
		
		/// Return local data
		UInt_t GetData(Int_t n) const
		{
			assert(n < 5);
			return (reinterpret_cast<const UInt_t*>(fLocalStruct))[n];
		}
	
		/// Return X2
		UShort_t GetX2() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalX2(fLocalStruct);}
		/// Return X1
		UShort_t GetX1() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalX1(fLocalStruct);}
		/// Return X4
		UShort_t GetX4() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalX4(fLocalStruct);}
		/// Return X3
		UShort_t GetX3() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalX3(fLocalStruct);}
	
		/// Return Y2
		UShort_t GetY2() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalY2(fLocalStruct);}
		/// Return Y1
		UShort_t GetY1() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalY1(fLocalStruct);}
		/// Return Y4
		UShort_t GetY4() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalY4(fLocalStruct);}
		/// Return Y3
		UShort_t GetY3() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalY3(fLocalStruct);}
		
		/// return X pattern array
		void GetXPattern(TArrayS& array) const
		{
			Short_t vec[4] = {static_cast<Short_t>(GetX1()), static_cast<Short_t>(GetX2()), static_cast<Short_t>(GetX3()), static_cast<Short_t>(GetX4())};
			array.Set(4, vec);
		}
	
		/// return Y pattern array
		void GetYPattern(TArrayS& array) const
		{
			Short_t vec[4] = {static_cast<Short_t>(GetY1()), static_cast<Short_t>(GetY2()), static_cast<Short_t>(GetY3()), static_cast<Short_t>(GetY4())};
			array.Set(4, vec);
		}
	
		/// Return Id
		UChar_t GetId() const {return fgOverrideId ? fCalculatedId : AliMUONTriggerDDLDecoderEventHandler::GetLocalId(fLocalStruct);}
		/// Return Dec
		UChar_t GetDec() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalDec(fLocalStruct);}
		/// Return TrigY
		Bool_t GetTrigY() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalTrigY(fLocalStruct);}
		/// Return TriggerY
		Bool_t GetTriggerY() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalTriggerY(fLocalStruct);}
		/// Return Upos
		UChar_t GetYPos() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalYPos(fLocalStruct);}
		/// Get Sign of X deviation 
		Bool_t GetSXDev() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalSXDev(fLocalStruct);}
		/// Get X deviation 
		UChar_t GetXDev() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalXDev(fLocalStruct);}
		/// Return TriggerX
		Bool_t GetTriggerX() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalTriggerX(fLocalStruct);}
		/// Return Xpos
		UChar_t GetXPos() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalXPos(fLocalStruct);}
		/// Return LPT
		UChar_t GetLpt() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalLpt(fLocalStruct);}
		/// Return HPT
		UChar_t GetHpt() const {return AliMUONTriggerDDLDecoderEventHandler::GetLocalHpt(fLocalStruct);}
		
		// Scaler methods
		/// Return local L0
		UInt_t  GetL0()      const {return (fScalars != NULL) ? fScalars->fL0 : 0;}
		/// Return local hold (dead time)
		UInt_t  GetHold()    const {return (fScalars != NULL) ? fScalars->fHold : 0;}
		/// Return local clock
		UInt_t  GetClock()   const {return (fScalars != NULL) ? fScalars->fClk : 0;}
		
		/// Return switch
		UShort_t GetSwitch() const
		{
			return (fScalars != NULL) ?
				AliMUONTriggerDDLDecoderEventHandler::GetLocalSwitch(fScalars) : 0;
		}
		
		/// Return ComptXY
		UChar_t GetComptXY() const
		{
			return (fScalars != NULL) ?
				AliMUONTriggerDDLDecoderEventHandler::GetLocalComptXY(fScalars) : 0;
		}
	
		/// Return XY1
		UShort_t GetXY1(Int_t n) const
		{
			return (fScalars != NULL) ?
				AliMUONTriggerDDLDecoderEventHandler::GetLocalXY1(fScalars, n) : 0;
		}
	
		/// Return XY2
		UShort_t GetXY2(Int_t n) const
		{
			return (fScalars != NULL) ?
				AliMUONTriggerDDLDecoderEventHandler::GetLocalXY2(fScalars, n) : 0;
		}
	
		/// Return XY3
		UShort_t GetXY3(Int_t n) const
		{
			return (fScalars != NULL) ?
				AliMUONTriggerDDLDecoderEventHandler::GetLocalXY3(fScalars, n) : 0;
		}
	
		/// Return XY4
		UShort_t GetXY4(Int_t n) const
		{
			return (fScalars != NULL) ?
				AliMUONTriggerDDLDecoderEventHandler::GetLocalXY4(fScalars, n) : 0;
		}

		/// Return raw data of the local trigger structure.
		const AliMUONLocalInfoStruct* GetData() const {return fLocalStruct;}
		
		/// Return raw data of the local trigger scalars.
		/// Could be NULL if no scalars present.
		const AliMUONLocalScalarsStruct* GetScalars() const {return fScalars;}
		
		/// Return the parent regional header.
		const AliRegionalHeader* GetRegionalHeader() const { return fRegional; }
		
		/// Return the next local trigger structure.
		const AliLocalStruct* Next() const { return fNext; }
		
		/// Sets the next local trigger structure.
		void SetNext(const AliLocalStruct* next) { fNext = next; }
	
		/// Sets the calculated ID value to be returned by GetId if fgOverrideId is true.
		void SetCalculatedId(UChar_t id) { fCalculatedId = id; }
		
		/// Print the contents of the local trigger structure and contents to screen.
		void Print() const;
		
		/// Returns the override flag indicating if the GetId method should return the calculated Id value or not.
		static bool GetOverrideIdFlag() { return fgOverrideId; }
		
		/// Sets the override flag to control what value the GetId method returns.
		static void SetOverrideIdFlag(bool value) { fgOverrideId = value; }
	
	private:
	
		const AliRegionalHeader* fRegional;  ///< The regional structure this local trigger structure belongs to.
		const AliLocalStruct* fNext;  ///< Next local structure object in the regional structure.
		const AliMUONLocalInfoStruct* fLocalStruct;  ///< Pointer to the local trigger structure data in the DDL payload.
		const AliMUONLocalScalarsStruct* fScalars;  ///< Pointer to the local trigger scalars data in the DDL payload.
		UChar_t fCalculatedId;  ///< Calculated ID value returned by GetId() if fgOverrideId == true.
		static bool fgOverrideId; //!<! Flag indicating if we should return a calculated number in the GetId method.
	};
	
	/// Returns the DARC and global headers plus scalars if they exist.
	const AliHeader* GetHeaders() const { return fDecoder.GetHandler().GetHeaders(); }
	
	/// Return the number of regional structures in the DDL payload.
	UInt_t GetRegionalHeaderCount() const
	{
		return fDecoder.GetHandler().RegionalHeaderCount();
	}
	
	/// Return the first regional structure header.
	const AliRegionalHeader* GetFirstRegionalHeader() const
	{
		return fDecoder.GetHandler().RegionalHeader(0);
	}
	
	/// Return the i'th regional header or NULL if not found.
	const AliRegionalHeader* GetRegionalHeader(UInt_t i) const
	{
		return fDecoder.GetHandler().RegionalHeader(i);
	}
	
	/// Returns the number of local trigger structures for the given
	/// regional structure number.
	UInt_t GetLocalStructCount(UInt_t reg) const
	{
		const AliRegionalHeader* r = GetRegionalHeader(reg);
		return r != NULL ? r->GetLocalStructCount() : 0;
	}
	
	/// Returns the i'th local trigger structure for the given regional
	/// structure number or NULL if not found.
	const AliLocalStruct* GetLocalStruct(UInt_t reg, UInt_t i) const
	{
		const AliRegionalHeader* r = GetRegionalHeader(reg);
		return r != NULL ? r->GetLocalStruct(i) : NULL;
	}

	/// Returns the current local struct being decoded or NULL if none found.
	const AliLocalStruct* CurrentLocalStruct() const
	{
		return (fkCurrentLocalStruct != fDecoder.GetHandler().EndOfLocalStructs()) ?
			fkCurrentLocalStruct : NULL;
	}

	/// Returns the current regional structure being decoded
	/// or NULL if none found.
	const AliRegionalHeader* CurrentRegionalHeader() const
	{
		const AliLocalStruct* local = CurrentLocalStruct();
		return (local != NULL) ? local->GetRegionalHeader() : NULL;
	}
	
private:
	/// Not implemented
	AliMUONRawStreamTriggerHP(const AliMUONRawStreamTriggerHP& stream);
	/// Not implemented
	AliMUONRawStreamTriggerHP& operator = (const AliMUONRawStreamTriggerHP& stream);
	
	
	/// This is the custom event handler (callback interface) class which
	/// unpacks raw data words and fills an internal buffer with decoded digits
	/// as they are decoded by the high performance decoder.
	/// Any errors are logged to the parent AliMUONVRawStreamTrigger, so one
	/// must set this pointer appropriately before decoding and DDL payload.
	class AliDecoderEventHandler : public AliMUONTriggerDDLDecoderEventHandler
	{
	public:
	
		/// Default constructor.
		AliDecoderEventHandler();
		/// Default destructor.
		virtual ~AliDecoderEventHandler();
		
		/// Sets the internal arrays based on the maximum number of structures allowed.
		void SetMaxStructs(UInt_t maxRegionals, UInt_t maxLocals);
		
		/// Sets the raw stream object which should be the parent of this class.
		void SetRawStream(AliMUONVRawStreamTrigger* rawStream) { fRawStream = rawStream; }
		
		/// Returns the decoded DARC and global headers (and scalars if they exist).
		const AliHeader* GetHeaders() const { return &fHeaders; }
		
		/// Returns the number of regional headers.
		UInt_t RegionalHeaderCount() const { return fRegionalsCount; }
		
		/// Return the i'th regional structure header.
		const AliRegionalHeader* RegionalHeader(UInt_t i) const
		{
			return i < fRegionalsCount ? &fRegionals[i] : NULL;
		}
		
		/// Return the first local structure decoded.
		const AliLocalStruct* FirstLocalStruct() const { return fLocals; }
		
		/// Returns the marker to the end of local structures.
		/// i.e. one position past the last local structure in fLocals.
		const AliLocalStruct* EndOfLocalStructs() const { return fEndOfLocals; }
		
		/// Get number of end of DARC word errors.
		UInt_t GetDarcEoWErrors() const {return fDarcEoWErrors;}
		
		/// Get number of end of Global word errors.
		UInt_t GetGlobalEoWErrors() const {return fGlobalEoWErrors;}
		
		/// Get number of end of regional word errors.
		UInt_t GetRegEoWErrors() const {return fRegEoWErrors;}
		
		/// Get number of end of local word errors.
		UInt_t GetLocalEoWErrors() const {return fLocalEoWErrors;}

		/// Returns the warnings flag.
		Bool_t Warnings() const { return fWarnings; }
		/// Sets the warnings flag.
		void Warnings(Bool_t value) { fWarnings = value; }
		
		// The following methods are inherited from AliMUONTriggerDDLDecoderEventHandler:
		
		/// New buffer handler.
		void OnNewBuffer(const void* buffer, UInt_t bufferSize);

		/// End of buffer handler marks the end of local trigger structures.
		void OnEndOfBuffer(const void* /*buffer*/, UInt_t /*bufferSize*/)
		{
			fEndOfLocals = fCurrentLocal+1;
		}
		
		/// Handler for the DARC header which just remembers the pointers for later.
		void OnDarcHeader(UInt_t header,
		                  const AliMUONDarcScalarsStruct* scalars,
		                  const void* /*data*/)
		{
			fDarcHeader = header;
			fDarcScalars = scalars;
		}
		
		/// Handler for the global header which stores the pointer to the headers.
		void OnGlobalHeader(const AliMUONGlobalHeaderStruct* header,
		                    const AliMUONGlobalScalarsStruct* scalars,
		                    const void* /*data*/)
		{
			fHeaders = AliHeader(fDarcHeader, fDarcScalars, header, scalars);
		}
		
		/// Handler for new regional card structures.
		void OnNewRegionalStructV2(UInt_t iReg,
		                           const AliMUONRegionalHeaderStruct* header,
		                           const AliMUONRegionalScalarsStruct* scalars,
		                           const void* data);
		
		/// Handler for new local card structures.
		void OnLocalStructV2(UInt_t iLoc,
		                     const AliMUONLocalInfoStruct* localStruct,
		                     const AliMUONLocalScalarsStruct* scalars);
		
		/// Error handler.
		void OnError(ErrorCode error, const void* location);
	
	private:
	
		// Do not allow copying of this class.
                /// Not implemented
		AliDecoderEventHandler(const AliDecoderEventHandler& /*obj*/);
                /// Not implemented
		AliDecoderEventHandler& operator = (const AliDecoderEventHandler& /*obj*/);

		AliMUONVRawStreamTrigger* fRawStream; //!<! Pointer to the parent raw stream object.
		const void* fBufferStart;   //!<! Pointer to the start of the current DDL payload buffer.
		UInt_t fDarcHeader; //!<! Currently decoded DARC header.
		const AliMUONDarcScalarsStruct* fDarcScalars; //!<! Currently decoded DARC scalars.
		AliHeader fHeaders;  //!<! Headers of the DDL payload.
		UInt_t fRegionalsCount; //!<! Number of regional headers filled in fRegionals.
		AliRegionalHeader* fRegionals;  //!<! Array of regional headers. [0..fMaxRegionals-1]
		AliLocalStruct* fLocals; //!<! Array of decoded local structured. [0..fMaxRegionals*fMaxLocals-1]
		AliLocalStruct* fEndOfLocals; //!<! Marker indicating the position just passed the last filled element in fLocals.
		AliRegionalHeader* fCurrentRegional; //!<! Current regional header position.
		AliLocalStruct* fCurrentLocal; //!<! Current local trigger structure.
		UInt_t fDarcEoWErrors;    //!<! Number of end of DARC word errors.
		UInt_t fGlobalEoWErrors;  //!<! Number of end of global word errors.
		UInt_t fRegEoWErrors;     //!<! Number of end of regional word errors.
		UInt_t fLocalEoWErrors;   //!<! Number of end of local word errors.
		Bool_t fWarnings;       //!<! Flag indicating if we should generate a warning for errors.
		
		static const AliMUONRegionalHeaderStruct fgkEmptyHeader;  //!<! Empty header for skipped regional structures.
	};
	
	AliMUONTriggerDDLDecoder<AliDecoderEventHandler> fDecoder;  //!<! The decoder for the DDL payload.
	Int_t fDDL;         //!<! The current DDL number being handled.
	Int_t fBufferSize;  //!<! This is the buffer size in bytes of fBuffer.
	UChar_t* fBuffer;   //!<! This is the buffer in which we store the DDL payload read from AliRawReader.
	const AliLocalStruct* fkCurrentLocalStruct;  //!<! The current local trigger structure being handled by Next().
	Bool_t fHadError;   //!<! Flag indicating if there was a decoding error or not.
	Bool_t fDone;       //!<! Flag indicating if the iteration is done or not.
	mutable AliMUONDDLTrigger* fDDLObject; //!<! Temporary DDL object used by GetDDLTrigger() for caching.
	UInt_t fTotalNumberOfDarcEoWErrors; //!<! The total number of end of DARC word errors since the last call to First().
	UInt_t fTotalNumberOfGlobalEoWErrors; //!<! The total number of end of global word errors since the last call to First().
	UInt_t fTotalNumberOfRegEoWErrors; //!<! The total number of end of regional word errors since the last call to First().
	UInt_t fTotalNumberOfLocalEoWErrors; //!<! The total number of end of local word errors since the last call to First().
	
	static const Int_t  fgkMaxDDL;     //!<! Maximum number of DDLs
	
	ClassDef(AliMUONRawStreamTriggerHP, 0)  // Higher performance decoder class for MUON trigger rawdata.
};

////////////////////////////////////////////////////////////////////////////////

inline const AliMUONRawStreamTriggerHP::AliLocalStruct* AliMUONRawStreamTriggerHP::Next()
{
	/// Iterates through all the local trigger structures for the event.
	/// When no more local triggers are found then NULL is returned.

	do {
		if (fkCurrentLocalStruct != fDecoder.GetHandler().EndOfLocalStructs())
			return fkCurrentLocalStruct++;
	} while (NextDDL());
	return NULL;
}


inline void AliMUONRawStreamTriggerHP::AliDecoderEventHandler::OnNewRegionalStructV2(
		UInt_t iReg,
		const AliMUONRegionalHeaderStruct* header,
		const AliMUONRegionalScalarsStruct* scalars,
		const void* /*data*/
	)
{
	/// New regional structure handler is called by the decoder whenever a
	/// new regional header is found. We just mark the header and increment
	/// the appropriate counters.
	
	assert( header != NULL );
	assert( iReg < fRegionalsCount );
	
	fCurrentRegional = fRegionals+iReg;
	*fCurrentRegional = AliRegionalHeader(fCurrentLocal+1, header, scalars);
	
	// Link to the next regional structure unless this is the last one.
	if (iReg+1 < fRegionalsCount)
	{
		fCurrentRegional->SetNext(fCurrentRegional+1);
	}
}


inline void AliMUONRawStreamTriggerHP::AliDecoderEventHandler::OnLocalStructV2(
		UInt_t iLoc,
		const AliMUONLocalInfoStruct* localStruct,
		const AliMUONLocalScalarsStruct* scalars
	)
{
	/// New local trigger structure handler.
	/// This is called by the high performance decoder when a new local trigger
	/// structure is found within the DDL payload.
	/// We mark the location of the structure in the DDL paytload and increment
	/// appropriate counters.

	assert( localStruct != NULL );
	assert( fCurrentLocal != NULL );
	assert( fCurrentRegional != NULL );
	assert( fCurrentRegional->GetLocalStructCount() < (UInt_t)fRawStream->GetMaxLoc() );
	
	// Link the previous local structure unless it is the first one. 
	if (fCurrentRegional->GetLocalStructCount() > 0)
	{
		fCurrentLocal->SetNext(fCurrentLocal+1);
	}
	
	fCurrentLocal++;
	*fCurrentLocal = AliLocalStruct(fCurrentRegional, localStruct, scalars);
	fCurrentLocal->SetCalculatedId(iLoc);
	fCurrentRegional->IncLocalStructCount();
}

#endif // ALIMUONRAWSTREAMTRIGGERHP_H
