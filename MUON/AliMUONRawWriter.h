#ifndef ALIMUONRAWWRITER_H
#define ALIMUONRAWWRITER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup base
/// \class AliMUONRawWriter
/// \brief Raw data class for trigger and tracker chambers
///
//  Author Christian Finck and Laurent Aphecetche, Subatech

#include <TObject.h>
#include "AliFstream.h"

class AliMUONBlockHeader;
class AliMUONBusStruct;
class AliMUONDarcHeader;
class AliMUONVDigit;
class AliMUONDspHeader;
class AliMUONGlobalTrigger;
class AliMUONLocalStruct;
class AliMUONRegHeader;
class AliMUONVDigitStore;
class AliMUONVTriggerStore;
class AliMpDDLStore;
class AliMpExMap;
class AliRawDataHeaderSim;

class AliMUONRawWriter : public TObject 
{
 public:
  AliMUONRawWriter(); // Constructor
  virtual ~AliMUONRawWriter(); // Destructor
    
  // write raw data
  Int_t Digits2Raw(const AliMUONVDigitStore* digitStore, const AliMUONVTriggerStore* triggerStore);
  
  void SetScalersNumbers();

  /// Set the header of DDL
  void SetHeader(AliRawDataHeaderSim& header) {fHeader = &header;}

private:

  void Digits2BusPatchMap(const AliMUONVDigitStore& digitStore, AliMpExMap& busPatchMap);
  void WriteTrackerDDL(AliMpExMap& busPatchMap, Int_t iDDL);

  //void WriteBusPatch(AliMUONLocalBusStruct* busStruct);
  
  Int_t WriteTriggerDDL(const AliMUONVTriggerStore& triggerStore, AliFstream* file[2]);
  
  Int_t GetBusPatch(const AliMUONVDigit& digit) const;

private:
  /// Not implemented copy constructor
  AliMUONRawWriter (const AliMUONRawWriter& rhs); // copy constructor
  /// Not implemented assignment operator
  AliMUONRawWriter& operator=(const AliMUONRawWriter& rhs);

 static void LocalWordPacking(UInt_t &word, UInt_t locId, UInt_t locDec, 
			      UInt_t trigY, UInt_t posY, UInt_t posX, 
			      UInt_t sdevX, UInt_t devX);

  AliMUONBlockHeader* fBlockHeader;  //!< DDL block header class pointers
  AliMUONDspHeader*   fDspHeader;    //!< DDL Dsp header class pointers
  AliMUONDarcHeader*  fDarcHeader;   //!< DDL darc header class pointers
  AliMUONRegHeader*   fRegHeader;    //!< DDL regional header class pointers
  AliMUONLocalStruct* fLocalStruct;  //!< DDL local structure class pointers

  AliMpDDLStore*            fDDLStore;     //!< DDL store pointer

  Bool_t fScalerEvent;               ///< flag to generates scaler event

  AliRawDataHeaderSim*    fHeader;           ///< header of DDL
  
  Int_t fBufferSize; //!< size of internal data buffer
  Int_t* fBuffer; //!< internal data buffer

  ClassDef(AliMUONRawWriter,5) // MUON cluster reconstructor in ALICE
};
	
#endif
