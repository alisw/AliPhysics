// $Id$

///**************************************************************************
///* This file is property of and copyright by the                          * 
///* ALICE Experiment at CERN, All rights reserved.                         *
///*                                                                        *
///* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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
/// @brief  Publisher component for raw data blocks through the AliRawReader
///         of the offline environment

#include "AliHLTRawReaderPublisherComponent.h"
#include "AliRawReader.h"
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliHLTCDHWrapper.h" 
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
  // constructor
}

AliHLTRawReaderPublisherComponent::~AliHLTRawReaderPublisherComponent()
{
  // destructor
}

const char* AliHLTRawReaderPublisherComponent::GetComponentID()
{
  /// inherited from AliHLTComponent: id of the component
  return "AliRawReaderPublisher";
}

AliHLTComponentDataType AliHLTRawReaderPublisherComponent::GetOutputDataType()
{
  /// inherited from AliHLTComponent: output data type of the component.
  return fDataType;
}

void AliHLTRawReaderPublisherComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  /// inherited from AliHLTComponent: output data size estimator
  constBase=fMaxSize;
  inputMultiplier=1;
}

AliHLTComponent* AliHLTRawReaderPublisherComponent::Spawn()
{
  /// inherited from AliHLTComponent: spawn function.
  return new AliHLTRawReaderPublisherComponent;
}

int AliHLTRawReaderPublisherComponent::DoInit( int argc, const char** argv )
{
  /// inherited from AliHLTComponent: component initialisation and argument scan.
  int iResult=0;

  // scan arguments
  if (argc && (iResult = ConfigureFromArgumentString(argc, argv)) < 0)
    return iResult;

  if (!fDetector.IsNull()) {
    int ddloffset=-1;
    int ddlcount=-1;
    if ((ddloffset=AliDAQ::DdlIDOffset(fDetector))<0 || 
	(ddlcount=AliDAQ::NumberOfDdls(fDetector))<0) {
      return -EINVAL;
    }
    if (fMinEquId<0) fMinEquId=ddloffset;
    else fMinEquId+=ddloffset;

    if (fMaxEquId<0 || fMaxEquId>ddlcount) fMaxEquId=ddloffset+ddlcount-1;
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
  /// inherited from AliHLTComponent: component cleanup
  int iResult=0;
  return iResult;
}

int AliHLTRawReaderPublisherComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  /// inherited from AliHLTComponent: argument scan
  if (argc<1) return 0;
  int bMissingParam=0;
  int i=0;
  TString argument=argv[i];

  do {
    // -detector
    if (argument.CompareTo("-detector")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fDetector=argv[i];
      return 2;

      // -equipmentid, -minid
    } else if (argument.CompareTo("-equipmentid")==0 ||
	       argument.CompareTo("-minid")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	fMinEquId=(AliHLTUInt32_t)parameter.Atoi();
	return 2;
      } else {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	return -EINVAL;
      }

      // -maxid
    } else if (argument.CompareTo("-maxid")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	fMaxEquId=(AliHLTUInt32_t)parameter.Atoi();
	return 2;
      } else {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	return -EINVAL;
      }

      // -verbose
    } else if (argument.CompareTo("-verbose")==0) {
      fVerbosity++;
      return 1;

      // -silent
    } else if (argument.CompareTo("-silent")==0) {
      fVerbosity=0;
      return 1;

      // -skipempty
    } else if (argument.CompareTo("-skipempty")==0) {
      fSkipEmpty=kTRUE;
      return 1;

      // -datatype
    } else if (argument.CompareTo("-datatype")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      memcpy(&fDataType.fID, argv[i], TMath::Min(kAliHLTComponentDataTypefIDsize, (Int_t)strlen(argv[i])));
      if ((bMissingParam=(++i>=argc))) break;
      memcpy(&fDataType.fOrigin, argv[i], TMath::Min(kAliHLTComponentDataTypefOriginSize, (Int_t)strlen(argv[i])));
      return 3;

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
	return -EINVAL;
      }
      return 2;
    }
  } while (0); // using do-while only to have break available
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    return -EINVAL;
  }

  return -EPROTO;
}

int AliHLTRawReaderPublisherComponent::GetEvent(const AliHLTComponentEventData& /*evtData*/, 
						AliHLTComponentTriggerData& /*trigData*/, 
						AliHLTUInt8_t* outputPtr, 
						AliHLTUInt32_t& size, 
						vector<AliHLTComponentBlockData>& outputBlocks)
{
  /// inherited from AliHLTDataSource: get the event
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
    UChar_t headerVersion=0;
    while (pRawReader->ReadHeader() && (iResult>=0 || iResult==-ENOSPC)) {
      const AliRawDataHeader* pHeaderV2=pRawReader->GetDataHeader();
      const AliRawDataHeaderV3* pHeaderV3=pRawReader->GetDataHeaderV3();
      AliHLTCDHWrapper pHeader;
      if(pHeaderV2) {
	pHeader=pHeaderV2;
      } else if (pHeaderV3) {
	pHeader=pHeaderV3;
      } else {
	HLTError("can not get data header from RawReader, skipping data block ...");
	continue;
      }
      // store header version for empty blocks later on
      // any found header will suffice
      headerVersion=pHeader.GetVersion();
      unsigned int headerSize=pHeader.GetHeaderSize();
      unsigned int readSize=pRawReader->GetDataSize()+headerSize;
      int id=pRawReader->GetEquipmentId();
      if (fMinEquId>id || fMaxEquId<id) {
	AliError(Form("id %d returned from RawReader is outside range [%d,%d]", id, fMinEquId, fMaxEquId));
	continue;
      }
      processedIds.push_back(id);
      bool isSelected=IsSelected(id);
      if (fVerbosity>0) {
	AliInfo(Form("got header for id %d, size %d, %s", id, readSize, isSelected?"selected":"discarded"));
      } else {
	AliDebug(0, Form("got header for id %d, size %d", id, readSize));
      }
      if (!isSelected) continue;
      if (readSize+offset<=capacity) {
	memcpy(outputPtr+offset, pHeader.GetHeader(), headerSize);
	if (readSize>headerSize) {
	  if (!pRawReader->ReadNext(outputPtr+offset+headerSize, readSize-headerSize)) {
	    AliError(Form("error reading %d bytes from RawReader %p", readSize-headerSize, pRawReader));
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
      AliHLTCDHWrapper header;
      AliRawDataHeader headerV2;
      AliRawDataHeaderV3 headerV3;
      if(headerVersion==2){
	headerV2.fSize=sizeof(AliRawDataHeader);
	const UInt_t* triggermask=pRawReader->GetTriggerPattern();
	if (triggermask) {
	  headerV2.fTriggerClassLow=triggermask[0];
	  headerV2.fROILowTriggerClassHigh=triggermask[1];
	}
	header=&headerV2;
      } else { //assuming V3 even if no header at all was found above
	if(! headerVersion){
	  AliWarning("No data header found! Creating dummy header for empty blocks assuming CDH v3.");
	}
        headerV3.fSize=sizeof(AliRawDataHeaderV3);
        const UInt_t* triggermask=pRawReader->GetTriggerPattern();
        if (triggermask) {
          headerV3.fTriggerClassLow=triggermask[0];
	  headerV3.fTriggerClassesMiddleLow=triggermask[1] & 0x3FFFF;
	  headerV3.fTriggerClassesMiddleHigh= (triggermask[1]>>18) | (triggermask[2]<<14);
          headerV3.fROILowTriggerClassHigh=(triggermask[2]>>18) | (triggermask[3]<<14);
        }
        header=&headerV3;
      }
      unsigned int headerSize=header.GetHeaderSize();
      processedIds.sort();
      list<int>::iterator curr=processedIds.begin();
      for (int id=fMinEquId; id<=fMaxEquId; id++) {
	if (curr!=processedIds.end() && *curr<=id) {
	  curr++;
	} else {
	  if (offset+headerSize<=capacity) {
	    HLTInfo("add empty data block for equipment id %d", id);
	    memcpy(outputPtr+offset, header.GetHeader(), headerSize);
	    AliHLTComponentBlockData bd;
	    FillBlockData( bd );
	    bd.fOffset = offset;
	    bd.fSize = headerSize;
	    bd.fDataType = fDataType;
	    if (fSpecification == kAliHLTVoidDataSpec) {
	      GetSpecificationFromEquipmentId(id, bd.fSpecification);
	    } else {
	      bd.fSpecification=fSpecification;
	    }
	    outputBlocks.push_back( bd );
	  } else {
	    // we keep the loop going in order to collect the full size
	    fMaxSize=offset+headerSize;
	    iResult=-ENOSPC;
	  }
	  offset+=headerSize;
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
  /// get the data specification from the equipment id
  /// default method just returns the equipment id
  return specification=id;
}
