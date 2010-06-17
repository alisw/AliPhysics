#ifndef ALIMUONPREPROCESSOR_H
#define ALIMUONPREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONPreprocessor
/// \brief Shuttle preprocessor for MUON subsystems (TRK and TRG)
/// 
//  Author Laurent Aphecetche

#ifndef ALI_PREPROCESSOR_H
#  include "AliPreprocessor.h"
#endif

class AliMUONVSubprocessor;
class TObjArray;

class AliMUONPreprocessor : public AliPreprocessor
{
public:
  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t Process(TMap* dcsAliasMap);
  virtual void Print(Option_t* opt="") const;

  /// Return info whether the current subprocessor(s) needs DCS or not  
  virtual Bool_t ProcessDCS() { return fProcessDCS; }

  /// Publish AliPreprocessor::Log function
  void Log(const char* message) { AliPreprocessor::Log(message); }
  
  /// Publish AliPreprocessor::GetFileSources function
  TList* GetFileSources(Int_t system, const char* id) 
  { return AliPreprocessor::GetFileSources(system,id); }

  /// Publish AliPreprocessor::Store function
  Bool_t Store(const char* pathLevel2, const char* pathLevel3, TObject* object,
               AliCDBMetaData* metaData, 
               Int_t validityStart = 0, Bool_t validityInfinite = kFALSE)
  {
    return AliPreprocessor::Store(pathLevel2,pathLevel3,object,metaData,
                                  validityStart,validityInfinite);
  }
  
  /// Publish AliPreprocessor::GetRunParameter
  const char* GetRunParameter(const char* param)
  {
    return AliPreprocessor::GetRunParameter(param);
  }
  
  /// Publish AliPreprocessor::GetFile function
  const char* GetFile(Int_t system, const char* id, const char* source)
  {
    return AliPreprocessor::GetFile(system,id,source);
  }  

  /// Publish AliPreprocessor::GetFromOCDB function
    AliCDBEntry* GetFromOCDB(const char* pathLevel2, const char* pathLevel3) {
      return AliPreprocessor::GetFromOCDB(pathLevel2,pathLevel3);      
    }

  /// Publish AliPreprocessor::GetFromOCDB function
  AliCDBEntry* GetGeometryFromOCDB()
  {
    return AliPreprocessor::GetGeometryFromOCDB();
  }
  
  /// Whether we can be used (e.g. whether we were properly initialized)
  Bool_t IsValid() const { return fIsValid; }
  
  /// Mark as invalid
  void Invalidate() { fIsValid = kFALSE; }
  
  /// Whether we should do something or not
  Bool_t IsApplicable() { return fIsApplicable; }
  
  /// Return log book parameter
  TString GetLogBookParam(const char* parname)
  { return TString(AliPreprocessor::GetRunParameter(parname)); }
  
protected:
  AliMUONPreprocessor(const char* detName, AliShuttleInterface* shuttle);
  virtual ~AliMUONPreprocessor();
  
  void Add(AliMUONVSubprocessor* subProcessor, Bool_t processDCS=kFALSE); 
  void ClearSubprocessors();
  
  Bool_t fIsValid; //!< whether we were correctly initialized
  Bool_t fIsApplicable; //!< whether we have something to do
  
private:
  /// Not implemented
  AliMUONPreprocessor(const AliMUONPreprocessor& rhs);
  /// Not implemented
  AliMUONPreprocessor& operator=(const AliMUONPreprocessor& rhs);
  
  AliMUONVSubprocessor* Subprocessor(Int_t i) const;
  
private:

  TObjArray* fSubprocessors; //!< sub processors to execute
  Bool_t fProcessDCS; //!< whether the current subprocessor(s) needs DCS or not

  ClassDef(AliMUONPreprocessor,4) // MUON Shuttle preprocessor
};

#endif
