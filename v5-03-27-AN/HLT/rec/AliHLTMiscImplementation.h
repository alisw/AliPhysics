//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTMISCIMPLEMENTATION_H
#define ALIHLTMISCIMPLEMENTATION_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/// @file   AliHLTMiscImplementation.h
/// @author Matthias Richter
/// @date   2009-07-07
/// @brief  Implementation of various glue functions implemented in dynamically
///         loaded libraries

#include "AliHLTMisc.h"

class AliHLTMiscImplementation : public AliHLTMisc
{
 public:
  AliHLTMiscImplementation();
  ~AliHLTMiscImplementation();

  int InitCDB(const char* cdbpath);

  int SetCDBRunNo(int runNo);
  int GetCDBRunNo() const;

  AliCDBEntry* LoadOCDBEntry(const char* path, int runNo=-1) const;

  TObject* ExtractObject(AliCDBEntry* entry) const;
  int CheckOCDBEntries(const TMap* const pMap) const;

  int InitMagneticField() const;

  AliHLTUInt64_t GetTriggerMask(AliRawReader* rawReader) const;

  AliHLTUInt32_t GetTimeStamp(AliRawReader* rawReader) const;
  AliHLTUInt32_t GetEventType(AliRawReader* rawReader) const;
  const char* GetBeamTypeFromGRP() const;

  Double_t GetBz();
  Double_t GetBz(const Double_t *r);
  void GetBxByBz(const Double_t r[3], Double_t b[3]);

  const TClass* IsAliESDHLTDecision() const;
  int Copy(const AliHLTGlobalTriggerDecision* pDecision, TObject* pESDHLTDecision) const;

  int InitStreamerInfos(const char* ocdbEntry) const;
  int InitStreamerInfos(TObjArray* pSchemas) const;

  void SetAliESDtrackOnlineModeFlag(bool mode) const;
  bool GetAliESDtrackOnlineModeFlag() const;

 private:

  ClassDef(AliHLTMiscImplementation, 0)
};

#endif //ALIHLTMISCIMPLEMENTATION_H
