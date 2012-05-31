// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTDomainEntry.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   20 Nov 2008
/// @brief  Implementation of the AliHLTDomainEntry class.
///
/// The AliHLTDomainEntry class is used to store information identifying a particular
/// HLT internal data block, or set of data blocks using wild card values. This
/// class is used by AliHLTTriggerDomain to store a list of data block classes
/// that should be readout by the HLT. The information identifying a data block is
/// the following:
///  - the data block type
///  - the data block's origin (detector name)
///  - the data block's specification (detector specific bits)
/// Several useful operators and methods are defined to help manipulate this
/// information in the AliHLTTriggerDomain class.

#include "AliHLTDomainEntry.h"
#include "Riostream.h"
#include "TString.h"
#include <cstring>
#include <cerrno>

ClassImp(AliHLTDomainEntry)


AliHLTDomainEntry::AliHLTDomainEntry() :
  TObject(),
  fExclude(kFALSE),
  fUseSpec(kFALSE),
  fType(kAliHLTVoidDataType),
  fSpecification(kAliHLTVoidDataSpec)
{
  // Default constructor.
}

AliHLTDomainEntry::AliHLTDomainEntry(const AliHLTDomainEntry& domain) :
  TObject(domain),
  fExclude(domain.fExclude),
  fUseSpec(domain.fUseSpec),
  fType(domain.fType),
  fSpecification(domain.fSpecification)
{
  // Copy constructor performs a deep copy.
}


AliHLTDomainEntry::AliHLTDomainEntry(const AliHLTComponentDataType& type) :
  TObject(),
  fExclude(kFALSE),
  fUseSpec(kFALSE),
  fType(type),
  fSpecification(kAliHLTVoidDataSpec)
{
  // Constructs a domain entry with a particular data type and any specification.
  // See header file for more information.
}


AliHLTDomainEntry::AliHLTDomainEntry(const char* blocktype, const char* origin) :
  TObject(),
  fExclude(kFALSE),
  fUseSpec(kFALSE),
  fType(),
  fSpecification(kAliHLTVoidDataSpec)
{
  // Constructs a domain entry with a particular data type and any specification.
  // See header file for more information.
  
  char id[kAliHLTComponentDataTypefIDsize];
  memset(&id, 0x0, sizeof(id));
  for (int i = 0; i < kAliHLTComponentDataTypefIDsize && blocktype[i] != '\0'; i++)
  {
    id[i] = blocktype[i];
  }
  fType = AliHLTComponentDataTypeInitializer(id, origin);
}


AliHLTDomainEntry::AliHLTDomainEntry(const AliHLTComponentDataType& type, UInt_t spec) :
  TObject(),
  fExclude(kFALSE),
  fUseSpec(kTRUE),
  fType(type),
  fSpecification(spec)
{
  // Constructs a domain entry with a particular data type and specification.
  // See header file for more information.
}


AliHLTDomainEntry::AliHLTDomainEntry(const char* blocktype, const char* origin, UInt_t spec) :
  TObject(),
  fExclude(kFALSE),
  fUseSpec(kTRUE),
  fType(),
  fSpecification(spec)
{
  // Constructs a domain entry with a particular data type and specification.
  // See header file for more information.
  
  char id[kAliHLTComponentDataTypefIDsize];
  memset(&id, 0x0, sizeof(id));
  for (int i = 0; i < kAliHLTComponentDataTypefIDsize && blocktype[i] != '\0'; i++)
  {
    id[i] = blocktype[i];
  }
  fType = AliHLTComponentDataTypeInitializer(id, origin);
}


AliHLTDomainEntry::AliHLTDomainEntry(Bool_t exclude, const AliHLTDomainEntry& domain) :
  TObject(domain),
  fExclude(exclude),
  fUseSpec(domain.fUseSpec),
  fType(domain.fType),
  fSpecification(domain.fSpecification)
{
  // Constructs a domain entry from an existing one but with the exclude flag set.
  // See header file for more information.
}


AliHLTDomainEntry::AliHLTDomainEntry(Bool_t exclude, const AliHLTComponentDataType& type) :
  TObject(),
  fExclude(exclude),
  fUseSpec(kFALSE),
  fType(type),
  fSpecification(kAliHLTVoidDataSpec)
{
  // Constructs a domain entry with the given data type, any specification
  // and the exclude flag set.
  // See header file for more information.
}


AliHLTDomainEntry::AliHLTDomainEntry(Bool_t exclude, const char* blocktype, const char* origin) :
  TObject(),
  fExclude(exclude),
  fUseSpec(kFALSE),
  fType(),
  fSpecification(kAliHLTVoidDataSpec)
{
  // Constructs a domain entry with a particular data type, any specification
  // and the exclude flag set.
  // See header file for more information.
  
  char id[kAliHLTComponentDataTypefIDsize];
  memset(&id, 0x0, sizeof(id));
  for (int i = 0; i < kAliHLTComponentDataTypefIDsize && blocktype[i] != '\0'; i++)
  {
    id[i] = blocktype[i];
  }
  fType = AliHLTComponentDataTypeInitializer(id, origin);
}


AliHLTDomainEntry::AliHLTDomainEntry(Bool_t exclude, const AliHLTComponentDataType& type, UInt_t spec) :
  TObject(),
  fExclude(exclude),
  fUseSpec(kTRUE),
  fType(type),
  fSpecification(spec)
{
  // Constructs a domain entry with a particular data type and specification,
  // and the exclude flag is set.
  // See header file for more information.
}


AliHLTDomainEntry::AliHLTDomainEntry(Bool_t exclude, const char* blocktype, const char* origin, UInt_t spec) :
  TObject(),
  fExclude(exclude),
  fUseSpec(kTRUE),
  fType(),
  fSpecification(spec)
{
  // Constructs a domain entry with a particular data type and specification,
  // and the exclude flag is set.
  // See header file for more information.
  
  char id[kAliHLTComponentDataTypefIDsize];
  memset(&id, 0x0, sizeof(id));
  for (int i = 0; i < kAliHLTComponentDataTypefIDsize && blocktype[i] != '\0'; i++)
  {
    id[i] = blocktype[i];
  }
  fType = AliHLTComponentDataTypeInitializer(id, origin);
}


AliHLTDomainEntry::~AliHLTDomainEntry()
{
  // Default destructor.
}


AliHLTDomainEntry& AliHLTDomainEntry::operator = (const AliHLTDomainEntry& domain)
{
  // The copy operator performs a deep copy.

  if (this==&domain) return *this;
  TObject::operator = (domain);
  fType = domain.fType;
  fUseSpec = domain.fUseSpec;
  fSpecification = domain.fSpecification;
  return *this;
}


bool AliHLTDomainEntry::IdenticalTo(const AliHLTDomainEntry& rhs) const
{
  // Checks if this domain entry is identical to 'rhs' and do not just have a
  // set intersection.
  // See header file for more information.
  
  if (not MatchExactly(fType, rhs.fType)) return false;
  return (fUseSpec == rhs.fUseSpec) and (fSpecification == rhs.fSpecification);
}


bool AliHLTDomainEntry::SubsetOf(const AliHLTDomainEntry& rhs) const
{
  // Checks if this domain entry is a subset of 'rhs'.
  // See header file for more information.

  if (*this != rhs) return false;
  bool thisTypeIsAny = strncmp(&fType.fID[0], kAliHLTAnyDataTypeID, kAliHLTComponentDataTypefIDsize) == 0;
  bool rhsTypeIsAny = strncmp(&rhs.fType.fID[0], kAliHLTAnyDataTypeID, kAliHLTComponentDataTypefIDsize) == 0;
  if (thisTypeIsAny and not rhsTypeIsAny) return false;
  bool thisOriginIsAny = strncmp(&fType.fOrigin[0], kAliHLTDataOriginAny, kAliHLTComponentDataTypefOriginSize) == 0;
  bool rhsOriginIsAny = strncmp(&rhs.fType.fOrigin[0], kAliHLTDataOriginAny, kAliHLTComponentDataTypefOriginSize) == 0;
  if (thisOriginIsAny and not rhsOriginIsAny) return false;
  bool thisSpecIsAny = not fUseSpec;
  bool rhsSpecIsAny = not rhs.fUseSpec;
  if (thisSpecIsAny and not rhsSpecIsAny) return false;
  return true;
}


bool AliHLTDomainEntry::IntersectWith(const AliHLTDomainEntry& rhs, AliHLTDomainEntry& result) const
{
  // Finds the set intersection between this domain entry and 'rhs'.
  // See header file for more information.

  if (*this != rhs) return false;
  bool thisTypeIsAny = strncmp(&fType.fID[0], kAliHLTAnyDataTypeID, kAliHLTComponentDataTypefIDsize) == 0;
  bool thisOriginIsAny = strncmp(&fType.fOrigin[0], kAliHLTDataOriginAny, kAliHLTComponentDataTypefOriginSize) == 0;
  bool thisSpecIsAny = not fUseSpec;
  const AliHLTComponentDataType& type = (not thisTypeIsAny) ? fType : rhs.fType;
  const AliHLTComponentDataType& origin = (not thisOriginIsAny) ? fType : rhs.fType;
  Bool_t useSpec;
  UInt_t spec;
  if (not thisSpecIsAny)
  {
    useSpec = fUseSpec;
    spec = fSpecification;
  }
  else
  {
    useSpec = rhs.fUseSpec;
    spec = rhs.fSpecification;
  }
  if (useSpec)
  {
    result = AliHLTDomainEntry(type | origin.fOrigin, spec);
  }
  else
  {
    result = AliHLTDomainEntry(type | origin.fOrigin);
  }
  return true;
}


void AliHLTDomainEntry::Print(Option_t* option) const
{
  // Inherited from TObject. Prints the domain entry contents.
  // See header file for more information.
  
  cout << AsString().Data();
  TString opt(option);
  if (opt.Contains("noendl")) return;
  cout << endl;
}


TString AliHLTDomainEntry::AsString() const
{
  // Returns a string representation of the domain entry.
  // See header file for more information.
  
  TString str;
  if (strncmp(&fType.fID[0], kAliHLTAnyDataTypeID, kAliHLTComponentDataTypefIDsize) == 0)
  {
    for (int i = 0; i < kAliHLTComponentDataTypefIDsize; i++) str += "*";
  }
  else
  {
    for (int i = 0; i < kAliHLTComponentDataTypefIDsize; i++)
    {
      if (fType.fID[i] != '\0')
        str += fType.fID[i];
      else
        str += "\\0";
    }
  }
  str += ":";
  if (strncmp(&fType.fOrigin[0], kAliHLTDataOriginAny, kAliHLTComponentDataTypefOriginSize) == 0)
  {
    for (int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++) str += "*";
  }
  else
  {
    for (int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++)
    {
      if (fType.fOrigin[i] != '\0')
        str += fType.fOrigin[i];
      else
        str += "\\0";
    }
  }
  str += ":";
  if (fUseSpec)
  {
    char num[16];
    sprintf(num, "0x%8.8X", fSpecification);
    str += num;
  }
  else
  {
    str += "**********";
  }
  return str;
}

int AliHLTDomainEntry::AsBinary(AliHLTUInt32_t buffer[4]) const
{
  // convert the data type and specification to a 32 byte buffer
  if (!buffer) return -EINVAL;

  AliHLTUInt32_t* tgt=buffer; 
  unsigned ii=0;

  // lower part of the data type id
  *tgt=0;
  for ( ii=0; ii<4; ii++ ) {
    *tgt |= ((AliHLTUInt32_t)(fType.fID[8-1-ii])) << (ii*8);
  }
  tgt++;
	  
  // upper part of the data type id
  *tgt=0;
  for ( ii=0; ii<4; ii++ ) {
    *tgt |= ((AliHLTUInt32_t)(fType.fID[8-5-ii])) << (ii*8);
  }
  tgt++;
  
  // data type origin
  *tgt=0;
  for ( ii=0; ii<4; ii++ ) {
    *tgt |= ((AliHLTUInt32_t)(fType.fOrigin[4-1-ii])) << (ii*8);
  }
  tgt++;
  
  // specification
  if (fUseSpec)
    *tgt = fSpecification;
  else
    *tgt = kAliHLTVoidDataSpec;

  return 0;
}
