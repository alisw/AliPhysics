// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTMisc.h
/// @author Matthias Richter
/// @date   2009-07-07
/// @brief  Definition of various glue functions implemented in dynamically
///         loaded libraries

#include "AliHLTMisc.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMisc);

AliHLTMisc::AliHLTMisc()
{
  // see header file for function documentation
}

AliHLTMisc::~AliHLTMisc()
{
  // see header file for function documentation
}

AliHLTMisc* AliHLTMisc::fgInstance=NULL;

AliHLTMisc& AliHLTMisc::Instance()
{
  // see header file for function documentation
  if (!fgInstance) {
    fgInstance=LoadInstance((AliHLTMisc*)NULL, "AliHLTMiscImplementation", ALIHLTMISC_LIBRARY);
  }
  if (!fgInstance) {
    AliHLTLogging log;
    fgInstance=new AliHLTMisc;
    log.Logging(kHLTLogError, "AliHLTMisc::Instance", "HLT Analysis", "falling back to default AliHLTMisc instance");
  }
  return *fgInstance;
}

int AliHLTMisc::InitCDB(const char* /*cdbpath*/)
{
  // default method, functionality is implemented in the child class
  return -EFAULT;
}

int AliHLTMisc::SetCDBRunNo(int /*runNo*/)
{
  // default method, functionality is implemented in the child class
  return -EFAULT;
}

int AliHLTMisc::GetCDBRunNo() const
{
  // default method, functionality is implemented in the child class
  return -1;
}

AliCDBEntry* AliHLTMisc::LoadOCDBEntry(const char* /*path*/, int /*runNo*/) const
{
  // default method, functionality is implemented in the child class
  return NULL;
}

TObject* AliHLTMisc::ExtractObject(AliCDBEntry* /*entry*/) const
{
  // default method, functionality is implemented in the child class
  return NULL;
}

int AliHLTMisc::CheckOCDBEntries(const TMap* const /*pMap*/) const
{
  // default method, functionality is implemented in the child class
  return -ENOENT;
}

int AliHLTMisc::InitMagneticField() const
{
  // default method, functionality is implemented in the child class
  return -EFAULT;
}

AliHLTUInt64_t AliHLTMisc::GetTriggerMask(AliRawReader* /*rawReader*/) const
{
  // default method, functionality is implemented in the child class
  return 0;
}

AliHLTUInt32_t AliHLTMisc::GetTimeStamp(AliRawReader* /*rawReader*/) const
{
  // default method, functionality is implemented in the child class
  return 0;
}

AliHLTUInt32_t AliHLTMisc::GetEventType(AliRawReader* /*rawReader*/) const
{
  // default method, functionality is implemented in the child class
  return 0;
}

Double_t AliHLTMisc::GetBz()
{
  // default method, functionality is implemented in the child class
  return 0.0;
}

Double_t AliHLTMisc::GetBz(const Double_t */*r*/)
{
  // default method, functionality is implemented in the child class
  return 0.0;
}

void AliHLTMisc::GetBxByBz(const Double_t r[3], Double_t b[3])
{
  // default method, functionality is implemented in the child class
  if (!r || !b) return;
  return;
}

const TClass* AliHLTMisc::IsAliESDHLTDecision() const
{
  // default method, functionality is implemented in the child class
  return NULL;
}

int AliHLTMisc::Copy(const AliHLTGlobalTriggerDecision* /*pDecision*/, TObject* /*pESDHLTDecision*/) const
{
  // default method, functionality is implemented in the child class
  return -EFAULT;
}

int AliHLTMisc::InitStreamerInfos(const char* /*ocdbEntry*/) const
{
  // default method, functionality is implemented in the child class
  return -EFAULT;
}

int AliHLTMisc::InitStreamerInfos(TObjArray* /*pSchemas*/) const
{
  // default method, functionality is implemented in the child class
  return -EFAULT;
}

void AliHLTMisc::SetAliESDtrackOnlineModeFlag(bool /*mode*/) const
{
  // default method, functionality is implemented in the child class
}

bool AliHLTMisc::GetAliESDtrackOnlineModeFlag() const
{
  // default method, functionality is implemented in the child class
  return false;
}

AliHLTMisc::AliOnlineGuard::AliOnlineGuard(bool mode)
  : fMode(false)
{
  // store the current value and set the flag
  fMode=AliHLTMisc::Instance().GetAliESDtrackOnlineModeFlag();
  AliHLTMisc::Instance().SetAliESDtrackOnlineModeFlag(mode);
}

AliHLTMisc::AliOnlineGuard::~AliOnlineGuard()
{
  // restore old value of the flag
  AliHLTMisc::Instance().SetAliESDtrackOnlineModeFlag(fMode);
}

ostream  &operator<<(ostream &out, const AliHLTComponentDataType &dt)
{
  // printout of AliHLTComponentDataType struct
  char id[kAliHLTComponentDataTypefIDsize+1];
  strncpy(id, dt.fID, kAliHLTComponentDataTypefIDsize);
  id[kAliHLTComponentDataTypefIDsize]=0;
  char origin[kAliHLTComponentDataTypefOriginSize+1];
  strncpy(origin, dt.fOrigin, kAliHLTComponentDataTypefOriginSize);
  origin[kAliHLTComponentDataTypefOriginSize]=0;
  out << "{" << id << ":" << origin << "}";
  return out;
}
