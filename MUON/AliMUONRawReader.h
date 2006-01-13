#ifndef ALIMUONRAWREADER_H
#define ALIMUONRAWREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONRawReader
/// \brief Raw data class for trigger and tracker chambers
///
/// Readding Raw data class for trigger and tracker chambers

#include <TObject.h>
#include "AliMpBusPatch.h"

class TClonesArray;
class TArrayI;
class AliLoader;
class AliMUONData;
class AliMUONDigit;
class AliMUONDDLTracker;
class AliMUONDDLTrigger;
class AliMUONGlobalTrigger;
class AliRawReader;
class AliMpSegFactory;

class AliMUONRawReader : public TObject 
{
 public:
  AliMUONRawReader(AliLoader* loader, AliMUONData* data); // Constructor
  virtual ~AliMUONRawReader(void); // Destructor
    
  // write raw data
  Int_t   Raw2Digits(AliRawReader* rawReader);

  Int_t ReadTrackerDDL(AliRawReader* rawReader);
  Int_t ReadTriggerDDL(AliRawReader* rawReader);

  AliMUONData*   GetMUONData() {return fMUONData;}

  Int_t GetMapping(Int_t buspatchId, UShort_t manuId, 
			  UChar_t channelId, AliMUONDigit* digit );

  AliMUONGlobalTrigger* GetGlobalTriggerPattern(Int_t gloTrg) const;


 protected:
  AliMUONRawReader();                  // Default constructor
  AliMUONRawReader (const AliMUONRawReader& rhs); // copy constructor
  AliMUONRawReader& operator=(const AliMUONRawReader& rhs); // assignment operator

 private:

  AliMUONData*  fMUONData;           //! Data container for MUON subsystem 
 
  AliLoader*    fLoader;             //! alice loader
 
  AliMpSegFactory* fSegFactory;      //! Mapping segmentation factory

   
  AliMUONDDLTracker* fDDLTracker;    //! DDL tracker class pointers
  AliMUONDDLTrigger* fDDLTrigger;    //! DDL trigger class pointers

  AliMpBusPatch* fBusPatchManager;    //! buspatch versus DE's & DDL

  ClassDef(AliMUONRawReader,1) // MUON cluster reconstructor in ALICE
};
	
#endif
