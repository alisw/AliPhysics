#ifndef ALITENDER_H
#define ALITENDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 31/08/2009

//==============================================================================
//   AliTender - Tender wagon providing access to ESD event and CDB.
//      The tender calls an arbitrary number of user algorithms that add or
//      correct information in ESD based on CDB info that was not available
//      during pass1 reconstruction.
//==============================================================================

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTask.h"
#endif

class AliCDBManager;
class AliESDEvent;
class AliESDInputHandler;
class AliTenderSupply;

class AliTender : public AliAnalysisTask {

private:
  Int_t                     fDebug;          // Debug level
  Int_t                     fRun;            //! Current run
  Bool_t                    fRunChanged;     //! Flag for run change.
  ULong_t                   fCDBkey;         //! Key to unlock CDB manager
  TString                   fDefaultStorage; // Default CDB storage
  AliCDBManager            *fCDB;            //! Pointer to CDB manager
  AliESDInputHandler       *fESDhandler;     //! Pointer to ESD input handler
  AliESDEvent              *fESD;            //! Pointer to current ESD event
  TObjArray                *fSupplies;       // Array of tender supplies
  TObjArray                *fCDBSettings;    // Array with CDB configuration
  
  AliTender(const AliTender &other);
  AliTender& operator=(const AliTender &other);

public:  
  AliTender();
  AliTender(const char *name);
  virtual ~AliTender();

  void                      AddSupply(AliTenderSupply *supply);
  Int_t                     GetRun() const {return fRun;}
  AliCDBManager            *GetCDBManager() const {return fCDB;}
  AliESDInputHandler       *GetESDhandler() const {return fESDhandler;}
  AliESDEvent              *GetEvent() const {return fESD;}
  TObjArray                *GetSupplies() const {return fSupplies;}
  Bool_t                    RunChanged() const {return fRunChanged;}
  // Configuration
  void                      SetDefaultCDBStorage(const char *dbString="local://$ALICE_ROOT/OCDB");
  void                      SetDebugLevel(Int_t debug) {fDebug = debug;}

  // Run control
  virtual void              ConnectInputData(Option_t *option = "");
  virtual void              CreateOutputObjects();
  virtual Bool_t            Notify() {return kTRUE;}
  virtual void              Exec(Option_t *option);
    
  ClassDef(AliTender,2)  // Class describing the tender car for ESD analysis
};
#endif
