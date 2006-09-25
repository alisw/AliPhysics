#ifndef ALIMUONRAWWRITER_H
#define ALIMUONRAWWRITER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONRawWriter
/// \brief Raw data class for trigger and tracker chambers
///
/// Writring Raw data class for trigger and tracker chambers

#include <TObject.h>
#include "AliMUONBusStruct.h"
#include "AliRawDataHeader.h"
#include "TStopwatch.h"

class AliMUONData;
class AliMUONDigit;
class AliMUONDspHeader;
class AliMUONBlockHeader;
class AliMUONDarcHeader;
class AliMUONRegHeader;
class AliMUONLocalStruct;
class AliMUONGlobalTrigger;
class AliMpBusPatch;
class AliMUONTriggerCrateStore;
class AliMpSegFactory;
class AliMpExMap;

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

  Int_t GetBusPatch(const AliMUONDigit& digit);

  Int_t GetGlobalTriggerPattern(const AliMUONGlobalTrigger* gloTrg) const;

private:

  AliMUONData*  fMUONData;           //!< Data container for MUON subsystem 
 
  FILE*         fFile[4];            //!< DDL binary file pointer one per 1/2 chamber, 4 for one station
   
  AliMUONBlockHeader* fBlockHeader;  //!< DDL block header class pointers
  AliMUONDspHeader*   fDspHeader;    //!< DDL Dsp header class pointers
  AliMUONBusStruct*   fBusStruct;    //!< DDL bus patch structure class pointers
  AliMUONDarcHeader*  fDarcHeader;   //!< DDL darc header class pointers
  AliMUONRegHeader*   fRegHeader;    //!< DDL regional header class pointers
  AliMUONLocalStruct* fLocalStruct;  //!< DDL local structure class pointers

  AliMpBusPatch*            fBusPatchManager; //!< buspatch versus DE's & DDL
  AliMUONTriggerCrateStore* fCrateManager;    //!< Crate array

  Bool_t fScalerEvent;               ///< flag to generates scaler event

  AliRawDataHeader    fHeader;           ///< header of DDL

  static Int_t fgManuPerBusSwp1B[12];   //!< array containing the first manuId for each buspatch st1, Bending
  static Int_t fgManuPerBusSwp1NB[12];  //!< array containing the first manuId for each buspatch st1, NBending

  static Int_t fgManuPerBusSwp2B[12];   //!< array containing the first manuId for each buspatch st2, Bending
  static Int_t fgManuPerBusSwp2NB[12];  //!< array containing the first manuId for each buspatch st2, NBending
  
  TStopwatch fTrackerTimer;             //!< time watcher for tracker part
  TStopwatch fTriggerTimer;             //!< time watcher for trigger part
  TStopwatch fMappingTimer;             //!< time watcher for mapping-tracker part
  
  AliMpSegFactory* fSegFactory;         //!< mapping segmentation factory
  
  AliMUONRawWriter (const AliMUONRawWriter& rhs); // copy constructor
  AliMUONRawWriter& operator=(const AliMUONRawWriter& rhs); // assignment operator

  ClassDef(AliMUONRawWriter,1) // MUON cluster reconstructor in ALICE
};
	
#endif
