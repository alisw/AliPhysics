#ifndef ALIMUONVSUBPROCESSOR_H
#define ALIMUONVSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONVSubprocessor
/// \brief Base class for a shuttle sub-task for MUON (either TRK or TRG)
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif

class TMap;
class TObjectArray;
class AliMUONPreprocessor;

class AliMUONVSubprocessor : public TNamed
{
public:
  AliMUONVSubprocessor(AliMUONPreprocessor* master,
                       const char* name="", const char* title="");
  virtual ~AliMUONVSubprocessor();
  
  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  
  /// Process this sub-task
  virtual UInt_t Process(TMap* dcsAliasMap) = 0;
  
protected:
  /// Return the pointer to our master
  AliMUONPreprocessor* Master() const { return fMaster; }

  Bool_t RemoveValuesOutsideRun ( TObjArray* values );
  
  /// Not implemented
  AliMUONVSubprocessor();
  /// Not implemented
  AliMUONVSubprocessor(const AliMUONVSubprocessor&);
  /// Not implemented
  AliMUONVSubprocessor& operator=(const AliMUONVSubprocessor&);
  
private:
  AliMUONPreprocessor* fMaster; ///< Pointer to our master
  UInt_t fStartTime; ///< Start time of run
  UInt_t fEndTime;   ///< End time of run
  
  ClassDef(AliMUONVSubprocessor,2) // Base class of MUON shuttle sub(pre)processors
};

#endif
