#ifndef ALIMUONTRIGGERELECTRONICS_H
#define ALIMUONTRIGGERELECTRONICS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim 
/// \class AliMUONTriggerElectronics
/// \brief Manager class for muon trigger electronics
///
/// Client of trigger board classes
///
/// \author Rachid Guernane (LPCCFd)

#ifndef ROOT_TTask
#  include "TTask.h"
#endif

#ifndef ROOT_TArrayI
#  include "TArrayI.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONTriggerCrate;
class AliMUONCalibrationData;
class AliMUONData;
class AliMUONGlobalTriggerBoard;
class AliMUONTriggerCrateStore;
class AliMUONLocalTrigger;
class AliMUONGlobalTrigger;

class AliMUONTriggerElectronics : public TTask
{
   public:
      AliMUONTriggerElectronics(AliMUONData* data = 0, 
                                AliMUONCalibrationData* calibData=0);
      virtual ~AliMUONTriggerElectronics();

      virtual void Exec(Option_t*);
      
      /// Set Crate config from ascii file
      virtual void SetDataSource(TString SourceFile = 
                                 "$ALICE_ROOT/MUON/mapping/data/stationTrigger/crate.dat") 
      {fSourceFileName = SourceFile;}

      virtual void Factory(AliMUONCalibrationData* calibData);
      void LoadMasks(AliMUONCalibrationData* calibData);

      virtual void Feed(UShort_t pattern[2][4]);
		  virtual void FeedM();
		  
      virtual void Reset();

      virtual void Scan(Option_t *option);

      virtual void LocalResponse();
      virtual void RegionalResponse();
      virtual void GlobalResponse();

      virtual void DumpOS();

      virtual void Digits2Trigger();
      virtual void Trigger();

   private:
      AliMUONTriggerElectronics(const AliMUONTriggerElectronics& right);
      AliMUONTriggerElectronics&  operator = (const AliMUONTriggerElectronics& right);
     
   private:
      TString                    fSourceFileName;     ///< Source file
      AliMUONTriggerCrateStore  *fCrates;             ///< Crate array
      AliMUONGlobalTriggerBoard *fGlobalTriggerBoard; ///< Global trigger board
      AliMUONData               *fMUONData;           //!< Data container for MUON subsystem
      AliMUONLocalTrigger*       fLocalTrigger;       //!< pointer for local trigger container
      AliMUONGlobalTrigger*      fGlobalTrigger;      //!< pointer for global trigger container

   ClassDef(AliMUONTriggerElectronics,2)
};
#endif
