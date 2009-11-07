//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTMISCIMPLEMENTATION_H
#define ALIHLTMISCIMPLEMENTATION_H_H
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

  AliCDBEntry* LoadOCDBEntry(const char* path, int runNo=-1, int version = -1, int subVersion = -1);

  TObject* ExtractObject(AliCDBEntry* entry);

  int InitMagneticField() const;

  AliHLTUInt64_t GetTriggerMask(AliRawReader* rawReader) const;

 private:

  ClassDef(AliHLTMiscImplementation, 0)
};

#endif //ALIHLTMISCIMPLEMENTATION_H
