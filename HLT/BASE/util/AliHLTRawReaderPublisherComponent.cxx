// $Id$

///**************************************************************************
///* This file is property of and copyright by the ALICE HLT Project        * 
///* ALICE Experiment at CERN, All rights reserved.                         *
///*                                                                        *
///* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
///*                  for The ALICE HLT Project.                            *
///*                                                                        *
///* Permission to use, copy, modify and distribute this software and its   *
///* documentation strictly for non-commercial purposes is hereby granted   *
///* without fee, provided that the above copyright notice appears in all   *
///* copies and that both the copyright notice and this permission notice   *
///* appear in the supporting documentation. The authors make no claims     *
///* about the suitability of this software for any purpose. It is          *
///* provided "as is" without express or implied warranty.                  *
///**************************************************************************/

/// @file   AliHLTRawReaderPublisherComponent.cxx
/// @author Matthias Richter
/// @date   
/// @brief  A general tree publisher component for the AliRawReader.
///

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTRawReaderPublisherComponent.h"
#include "AliRawReader.h"
#include "AliDAQ.h"
#include "AliLog.h"
#include <cerrno>
#include <cassert>
#include <list>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRawReaderPublisherComponent)

AliHLTRawReaderPublisherComponent::AliHLTRawReaderPublisherComponent()
  :
  fMaxSize(5000000),
  fDetector(),
  fMinEquId(-1),
  fMaxEquId(-1),
  fVerbosity(0),
  fDataType(kAliHLTVoidDataType),
  fSpecification(kAliHLTVoidDataSpec),
  fSkipEmpty(kFALSE)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTRawReaderPublisherComponent::~AliHLTRawReaderPublisherComponent()
{
  // see header file for class documentation
}

const char* AliHLTRawReaderPublisherComponent::GetComponentID()
{
  // see header file for class documentation
  return "AliRawReaderPublisher";
}

AliHLTComponentDataType AliHLTRawReaderPublisherComponent::GetOutputDataType()
{
  // see header file for class documentation
  return fDataType;
}

void AliHLTRawReaderPublisherComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=fMaxSize;
  inputMultiplier=1;
}

AliHLTComponent* AliHLTRawReaderPublisherComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTRawReaderPublisherComponent;
}

int AliHLTRawReaderPublisherComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  // scan arguments
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -detector
    if (argument.CompareTo("-detector")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fDetector=argv[i];

      // -equipmentid, -minid
    } else if (argument.CompareTo("-equipmentid")==0 ||
	       argument.CompareTo("-minid")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	fMinEquId=(AliHLTUInt32_t)parameter.Atoi();
      } else {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }

      // -maxid
    } else if (argument.CompareTo("-maxid")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	fMaxEquId=(AliHLTUInt32_t)parameter.Atoi();
      } else {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }

      // -verbose
    } else if (argument.CompareTo("-verbose")==0) {
      fVerbosity++;

      // -silent
    } else if (argument.CompareTo("-silent")==0) {
      fVerbosity=0;

      // -skipempty
    } else if (argument.CompareTo("-skipempty")==0) {
      fSkipEmpty=kTRUE;

      // -datatype
    } else if (argument.CompareTo("-datatype")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      memcpy(&fDataType.fID, argv[i], TMath::Min(kAliHLTComponentDataTypefIDsize, (Int_t)strlen(argv[i])));
      if ((bMissingParam=(++i>=argc))) break;
      memcpy(&fDataType.fOrigin, argv[i], TMath::Min(kAliHLTComponentDataTypefOriginSize, (Int_t)strlen(argv[i])));

      // -dataspec
    } else if (argument.CompareTo("-dataspec")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	fSpecification=(AliHLTUInt32_t)parameter.Atoi();
      } else if (parameter.BeginsWith("0x") &&
		 parameter.Replace(0,2,"",0).IsHex()) {
	sscanf(parameter.Data(),"%x", &fSpecification);
      } else {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
    } else {
      HLTError("unknown argument %s", argument.Data());
      iResult=-EINVAL;
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  if (iResult<0) return iResult;

  if (!fDetector.IsNull()) {
    int ddloffset=-1;
    int ddlcount=-1;
    if ((ddloffset=AliDAQ::DdlIDOffset(fDetector))<0 || 
	(ddlcount=AliDAQ::NumberOfDdls(fDetector))<0) {
      return -EINVAL;
    }
    if (fMinEquId<0) fMinEquId=ddloffset;
    else fMinEquId+=ddloffset;

    if (fMaxEquId<0 || fMaxEquId>ddlcount) fMaxEquId=ddloffset+ddlcount;
    else fMaxEquId+=ddloffset;
  }

  if (fMinEquId>fMaxEquId) fMaxEquId=fMinEquId;

  if (fMinEquId<0) {
    AliErrorStream() << "equipment id required, use \'-equipmentid\' or \'-detector\' option" << endl;
    return -EINVAL;
  }

  AliHLTUInt32_t dummy;
  if (fMinEquId!=fMaxEquId && GetSpecificationFromEquipmentId(0, dummy)==-ENOSYS) {
    AliWarningStream() << "publication of multiple equipment ids needs implementation of a child and function GetSpecificationFromEquipmentId to set correct specifications" << endl;
    //return -EINVAL;
  }

  if (GetRawReader()!=NULL) {
  } else {
    AliErrorStream() << "RawReader instance needed" << endl;
    return -EINVAL;
  }

  return iResult;
}

int AliHLTRawReaderPublisherComponent::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}

int AliHLTRawReaderPublisherComponent::GetEvent(const AliHLTComponentEventData& /*evtData*/, 
						AliHLTComponentTriggerData& /*trigData*/, 
						AliHLTUInt8_t* outputPtr, 
						AliHLTUInt32_t& size, 
						vector<AliHLTComponentBlockData>& outputBlocks)
{
  // see header file for class documentation
  int iResult=0;
  AliHLTUInt32_t capacity=size;
  size=0;

  // process data events only
  if (!IsDataEvent()) return 0;

  unsigned int offset=0;
  assert(outputPtr!=NULL || size==0);
  AliRawReader* pRawReader=GetRawReader();
  if (pRawReader) {
    pRawReader->Reset();
    pRawReader->SelectEquipment(-1, fMinEquId, fMaxEquId);
    if (fVerbosity>1) {
      AliInfo(Form("get event from RawReader %p equipment id range [%d,%d]", pRawReader, fMinEquId, fMaxEquId));
    } else {
      AliDebug(0, Form("get event from RawReader %p equipment id range [%d,%d]", pRawReader, fMinEquId, fMaxEquId));
    }
    list<int> processedIds;
    while (pRawReader->ReadHeader() && (iResult>=0 || iResult==-ENOSPC)) {
      const AliRawDataHeader* pHeader=pRawReader->GetDataHeader();
      if (pHeader==NULL) {
	HLTError("can not get data header from RawReader, skipping data block ...");
	continue;
      }
      unsigned int readSize=pRawReader->GetDataSize()+sizeof(AliRawDataHeader);
      int id=pRawReader->GetEquipmentId();
      if (fVerbosity>0) {
	AliInfo(Form("got header for id %d, size %d", id, readSize));
      } else {
	AliDebug(0, Form("got header for id %d, size %d", id, readSize));
      }
      if (fMinEquId>id || fMaxEquId<id) {
	AliError(Form("id %d returned from RawReader is outside range [%d,%d]", id, fMinEquId, fMaxEquId));
	continue;
      }
      processedIds.push_back(id);
      if (!IsSelected(id)) continue;
      if (readSize+offset<=capacity) {
	memcpy(outputPtr+offset, pHeader, sizeof(AliRawDataHeader));
	if (readSize>sizeof(AliRawDataHeader)) {
	  if (!pRawReader->ReadNext(outputPtr+offset+sizeof(AliRawDataHeader), readSize-sizeof(AliRawDataHeader))) {
	    AliError(Form("error reading %ld bytes from RawReader %p", readSize-sizeof(AliRawDataHeader), pRawReader));
	    iResult=-ENODATA;
	    break;
	  }
	}
	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	bd.fOffset = offset;
	bd.fSize = readSize;
	bd.fDataType = fDataType;
	if (fSpecification == kAliHLTVoidDataSpec) {
	  GetSpecificationFromEquipmentId(id, bd.fSpecification);
	} else {
	  bd.fSpecification=fSpecification;
	}
	outputBlocks.push_back( bd );
      } else {
	// we keep the loop going in order to collect the full size
	fMaxSize=offset+readSize;
	iResult=-ENOSPC;
      }
      offset+=readSize;
    }
    if (!fSkipEmpty && processedIds.size()!=size_t(fMaxEquId-fMinEquId+1)) {
      // add further empty data blocks
      AliRawDataHeader header;
      header.fSize=sizeof(AliRawDataHeader);
      const UInt_t* triggermask=pRawReader->GetTriggerPattern();
      if (triggermask) {
	header.fTriggerClassLow=triggermask[0];
	header.fROILowTriggerClassHigh=triggermask[1];
      }
      processedIds.sort();
      list<int>::iterator curr=processedIds.begin();
      for (int id=fMinEquId; id<=fMaxEquId; id++) {
	if (curr!=processedIds.end() && *curr<=id) {
	  curr++;
	} else {
	  if (sizeof(AliRawDataHeader)<=capacity-offset) {
	    HLTInfo("add empty data block for equipment id %d", id);
	    memcpy(outputPtr+offset, &header, sizeof(AliRawDataHeader));
	    AliHLTComponentBlockData bd;
	    FillBlockData( bd );
	    bd.fOffset = offset;
	    bd.fSize = sizeof(AliRawDataHeader);
	    bd.fDataType = fDataType;
	    if (fSpecification == kAliHLTVoidDataSpec) {
	      GetSpecificationFromEquipmentId(id, bd.fSpecification);
	    } else {
	      bd.fSpecification=fSpecification;
	    }
	    outputBlocks.push_back( bd );
	  } else {
	    // we keep the loop going in order to collect the full size
	    fMaxSize=offset+sizeof(AliRawDataHeader);
	    iResult=-ENOSPC;
	  }
	  offset+=sizeof(AliRawDataHeader);
	}
      }
    }
    if (offset<=capacity) {
      size=offset;
    } else {
      outputBlocks.clear();
    }
  } else {
    AliErrorStream() << "RawReader uninitialized" << endl;
    iResult=-EFAULT;
  }
  return iResult;
}

int AliHLTRawReaderPublisherComponent::GetSpecificationFromEquipmentId(int id, AliHLTUInt32_t& specification) const {
  // see header file for class documentation
  return specification=id;
}
