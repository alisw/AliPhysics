// @(#) $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   AliHLTRawReaderPublisherComponent.cxx
    @author Matthias Richter
    @date   
    @brief  A general tree publisher component for the AliRawReader.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTRawReaderPublisherComponent.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include <cerrno>
#include <cassert>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRawReaderPublisherComponent)

AliHLTRawReaderPublisherComponent::AliHLTRawReaderPublisherComponent()
  :
  fMaxSize(5000000),
  fDetector(),
  fMinEquId(-1),
  fMaxEquId(-1),
  fVerbose(kFALSE),
  fDataType(kAliHLTVoidDataType),
  fSpecification(kAliHLTVoidDataSpec)
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
      fVerbose=kTRUE;

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

  if (fMinEquId>fMaxEquId) fMaxEquId=fMinEquId;

  if (fMinEquId<0) {
    AliErrorStream() << "equipment id required, use \'-equipmentid\' option" << endl;
    return -EINVAL;
  }

  if (!fDetector.IsNull()) {
    AliErrorStream() << "option \'-detector\' not implemented" << endl;
    return -ENOSYS;
  }

  AliHLTUInt32_t dummy;
  if (fMinEquId!=fMaxEquId && GetSpecificationFromEquipmentId(0, dummy)==-ENOSYS) {
    AliWarningStream() << "publication of multiple equipment ids needs implementation of a child and function GetSpecificationFromEquipmentId to set correct specifications" << endl;
    //return -EINVAL;
  }

  AliRawReader* pRawReader=GetRawReader();
  if ((pRawReader=GetRawReader())!=NULL) {
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

int AliHLTRawReaderPublisherComponent::GetEvent(const AliHLTComponentEventData& evtData, 
						AliHLTComponentTriggerData& trigData, 
						AliHLTUInt8_t* outputPtr, 
						AliHLTUInt32_t& size, 
						vector<AliHLTComponentBlockData>& outputBlocks)
{
  // see header file for class documentation
  int iResult=0;
  unsigned int offset=0;
  AliHLTUInt8_t* pTgt=outputPtr;
  assert(outputPtr!=NULL || size==0);
  AliRawReader* pRawReader=GetRawReader();
  if (pRawReader) {
    pRawReader->Reset();
    pRawReader->SelectEquipment(-1, fMinEquId, fMaxEquId);
    AliInfo(Form("get event from RawReader %p equipment id range [%d,%d]", pRawReader, fMinEquId, fMaxEquId));
    while (pRawReader->ReadHeader() && (iResult>=0 || iResult==-ENOSPC)) {
      const AliRawDataHeader* pHeader=pRawReader->GetDataHeader();
      assert(pHeader!=NULL);
      if (pHeader==NULL) continue;
      unsigned int readSize=pRawReader->GetDataSize()+sizeof(AliRawDataHeader);
      int id=pRawReader->GetEquipmentId();
      AliInfo(Form("got header for id %d, size %d", id, readSize));
      if (fMinEquId>id || fMaxEquId<id) {
	AliError(Form("id %d returned from RawReader is outside range [%d,%d]", id, fMinEquId, fMaxEquId));
	continue;
      }
      if (readSize<=size-offset) {
	memcpy(pTgt, pHeader, sizeof(AliRawDataHeader));
	pTgt+=sizeof(AliRawDataHeader);
	if (readSize>0) {
	  if (!pRawReader->ReadNext(pTgt, readSize-sizeof(AliRawDataHeader))) {
	    AliError(Form("error reading %d bytes from RawReader %p", readSize-sizeof(AliRawDataHeader), pRawReader));
	    iResult=-ENODATA;
	    break;
	  }
	  pTgt+=readSize-sizeof(AliRawDataHeader);
	}
      } else {
	// we keep the loop going in order to collect the full size
	fMaxSize=offset+readSize;
	iResult=-ENOSPC;
      }
      if (iResult>=0) {
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
      }
      offset+=readSize;
    }
    if (offset<=size) size=offset;
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
