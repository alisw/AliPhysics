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

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONTriggerCrate;
class AliMUONCalibrationData;
class AliMUONGlobalTriggerBoard;
class AliMUONTriggerCrateStore;
class AliMUONVTriggerStore;
class AliMUONVDigitStore;

class AliMUONTriggerElectronics : public TObject
{
   public:
      AliMUONTriggerElectronics(AliMUONCalibrationData* calibData=0);

      virtual ~AliMUONTriggerElectronics();

      virtual void Factory(AliMUONCalibrationData* calibData);
      void LoadMasks(AliMUONCalibrationData* calibData);

      virtual void Feed(UShort_t pattern[2][4]);
      virtual void Feed(const AliMUONVDigitStore& digitStore);
      virtual void Feed(const AliMUONVTriggerStore& triggerStore);
      virtual void Reset();

      virtual void Scan(const Option_t *option);

      virtual void LocalResponse();
      virtual void RegionalResponse();
      virtual void GlobalResponse();

      virtual void DumpOS();

      virtual void Digits2Trigger(const AliMUONVDigitStore& digitStore,
                                  AliMUONVTriggerStore& triggerStore);

      virtual Bool_t ModifiedLocalResponse(Int_t loCircuit,
					   Bool_t& bendingPlaneResp,
					   Bool_t& nonBendingPlaneResp,
					   Bool_t isCoinc44 = kFALSE,
					   Int_t removeChamber = -1);

      virtual void ResponseRemovingChambers(AliMUONVTriggerStore& triggerStore);

      /// Return Crate array
      AliMUONTriggerCrateStore* GetCrateStore() {return fCrates;}

   private:
      /// Not implemented
      AliMUONTriggerElectronics(const AliMUONTriggerElectronics& right);
      /// Not implemented
      AliMUONTriggerElectronics&  operator = (const AliMUONTriggerElectronics& right);

      void FeedCopyNeighbours();
     
      /// set copy card array
      void SetCopyInput();

   private:
      TList*                     fCopyXInput[2];         ///< list of copy X input from local to local board
      TList*                     fCopyYInput[2];         ///< list of copy Y input from local to local board
      AliMUONTriggerCrateStore  *fCrates;             ///< Crate array
      AliMUONGlobalTriggerBoard *fGlobalTriggerBoard; ///< Global trigger board
      
   ClassDef(AliMUONTriggerElectronics,4) //Manager class for muon trigger electronics
};
#endif
