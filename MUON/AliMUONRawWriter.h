#ifndef ALIMUONRAWWRITER_H
#define ALIMUONRAWWRITER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup sim
/// \class AliMUONRawWriter
/// \brief Raw data class for trigger and tracker chambers
///
/// \author Christian Finck and Laurent Aphecetche, Subatech

#include <TObject.h>
#include "AliRawDataHeader.h"

class AliMUONData;
class AliMUONDigit;
class AliMUONDspHeader;
class AliMUONBlockHeader;
class AliMUONDarcHeader;
class AliMUONRegHeader;
class AliMUONLocalStruct;
class AliMpExMap;
class AliMUONBusStruct;
class AliMUONGlobalTrigger;
class AliMpDDLStore;
class AliMUONTriggerCrateStore;
class TStopwatch;

class AliMUONRawWriter : public TObject 
{
 public:
  AliMUONRawWriter(AliMUONData* data); // Constructor
  virtual ~AliMUONRawWriter(); // Destructor
    
  // write raw data
  Int_t Digits2Raw();

  void SetScalersNumbers();

protected:
  AliMUONRawWriter();                  // Default constructor

  // writing raw data
  Int_t WriteTrackerDDL(Int_t iCh);
  Int_t WriteTriggerDDL();
  
private:

    void AddData(const AliMUONBusStruct& event);

  Int_t GetBusPatch(const AliMUONDigit& digit) const;
  Int_t GetBusPatch(Int_t detElemId, Int_t manuId) const;

  Int_t GetGlobalTriggerPattern(const AliMUONGlobalTrigger* gloTrg) const;

private:

  AliMUONData*  fMUONData;           //!< Data container for MUON subsystem 
 
  FILE*         fFile[4];            //!< DDL binary file pointer one per 1/2 chamber, 4 for one station

  AliMUONBlockHeader* fBlockHeader;  //!< DDL block header class pointers
  AliMUONDspHeader*   fDspHeader;    //!< DDL Dsp header class pointers
  AliMUONDarcHeader*  fDarcHeader;   //!< DDL darc header class pointers
  AliMUONRegHeader*   fRegHeader;    //!< DDL regional header class pointers
  AliMUONLocalStruct* fLocalStruct;  //!< DDL local structure class pointers

  AliMpDDLStore*            fBusPatchManager; //!< buspatch versus DE's & DDL
  AliMUONTriggerCrateStore* fCrateManager;    //!< Crate array

  Bool_t fScalerEvent;               ///< flag to generates scaler event

  AliRawDataHeader    fHeader;           ///< header of DDL

  static Int_t fgManuPerBusSwp1B[12];   //!< array containing the first manuId for each buspatch st1, Bending
  static Int_t fgManuPerBusSwp1NB[12];  //!< array containing the first manuId for each buspatch st1, NBending

  static Int_t fgManuPerBusSwp2B[12];   //!< array containing the first manuId for each buspatch st2, Bending
  static Int_t fgManuPerBusSwp2NB[12];  //!< array containing the first manuId for each buspatch st2, NBending
  
  TStopwatch* fTimers;             //!< time watchers
  
  AliMUONRawWriter (const AliMUONRawWriter& rhs); // copy constructor
  AliMUONRawWriter& operator=(const AliMUONRawWriter& rhs); // assignment operator

  ClassDef(AliMUONRawWriter,2) // MUON cluster reconstructor in ALICE
};
	
#endif
