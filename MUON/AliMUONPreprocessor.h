#ifndef ALIMUONPREPROCESSOR_H
#define ALIMUONPREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONPreprocessor
/// \brief Shuttle preprocessor for MUON subsystems (TRK and TRG)
/// 
/// \author Laurent Aphecetche

#ifndef ALI_PREPROCESSOR_H
#  include "AliPreprocessor.h"
#endif

class AliMUONVSubprocessor;
class TObjArray;

class AliMUONPreprocessor : public AliPreprocessor
{
public:
  AliMUONPreprocessor(const TString& detName, AliShuttleInterface* shuttle);
  virtual ~AliMUONPreprocessor();
  
  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t Process(TMap* dcsAliasMap);
  virtual void Print(Option_t* opt="") const;
  
  TList* GetFileSources(Int_t system, const char* id) 
  { return AliPreprocessor::GetFileSources(system,id); }

  UInt_t Store(const char* pathLevel2, const char* pathLevel3, TObject* object,
               AliCDBMetaData* metaData, 
               Int_t validityStart = 0, Bool_t validityInfinite = kFALSE)
  {
    return AliPreprocessor::Store(pathLevel2,pathLevel3,object,metaData,
                                  validityStart,validityInfinite);
  }
  
  const char* GetFile(Int_t system, const char* id, const char* source)
  {
    return AliPreprocessor::GetFile(system,id,source);
  }  
  
private:
  // enum ESubprocessors { kPedestal=0, kGMS=1, kLast };
  
  AliMUONPreprocessor(const AliMUONPreprocessor& rhs);
  AliMUONPreprocessor& operator=(const AliMUONPreprocessor& rhs);
  
  AliMUONVSubprocessor* Subprocessor(Int_t i) const;
  
private:
  static const TString  fgkTrackerDetName;  ///< The tracker detector name
  static const TString  fgkTriggerDetName;  ///< The trigger detector name

  TObjArray* fSubprocessors; ///!< sub processors to execute
  
  ClassDef(AliMUONPreprocessor,1) // MUON Shuttle preprocessor
};

#endif
