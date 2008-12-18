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

/** @file   AliHLTDataGenerator.cxx
    @author Matthias Richter
    @date   
    @brief  HLT file publisher component implementation. */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTDataGenerator.h"
#include "TString.h"
#include "TRandom.h"
#include "TDatime.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataGenerator)

AliHLTDataGenerator::AliHLTDataGenerator()
  :
  AliHLTProcessor(),
  fDataType(kAliHLTVoidDataType),
  fSpecification(~(AliHLTUInt32_t)0),
  fSize(0),
  fRange(0),
  fCurrSize(0),
  fDivisor(0),
  fDecrement(0),
  fModulo(0),
  fOffset(0),
  fMultiplier(1.0)
  , fpDice(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  // make the lists owners of their objects in order to automatically
  // de-allocate the objects
}

AliHLTDataGenerator::~AliHLTDataGenerator()
{
  // see header file for class documentation

}

const char* AliHLTDataGenerator::GetComponentID()
{
  // see header file for class documentation
  return "DataGenerator";
}

void AliHLTDataGenerator::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTDataGenerator::GetOutputDataType()
{
  // see header file for class documentation
  return fDataType;
}

int AliHLTDataGenerator::GetOutputDataTypes(vector<AliHLTComponentDataType>& tgtList)
{
  // see header file for class documentation
  int count=0;
  tgtList.clear();
  tgtList.push_back(fDataType);
  return count;
}

void AliHLTDataGenerator::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  if (fSize>0)
    constBase=(unsigned long)(fCurrSize+fRange);
  else
    constBase=(unsigned long)(fOffset+fRange);
  inputMultiplier=fMultiplier;
}

AliHLTComponent* AliHLTDataGenerator::Spawn()
{
  // see header file for class documentation
  return new AliHLTDataGenerator;
}

int AliHLTDataGenerator::DoInit( int argc, const char** argv )
{
  // see header file for class documentation

  //HLTDebug("%d %s", argc, argv[0]);
  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -datatype
    if (argument.CompareTo("-datatype")==0) {
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
      // -size
    } else if (argument.CompareTo("-size")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      if ((iResult=ScanSizeArgument(fSize, argv[i]))==-ERANGE) {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
      // -minsize || -maxsize
    } else if (argument.CompareTo("-minsize")==0 ||
	       argument.CompareTo("-maxsize")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      AliHLTUInt32_t value=0;
      if ((iResult=ScanSizeArgument(value, argv[i]))==-ERANGE) {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      } else {
	if (fSize==0) {
	  fSize=value;
	} else if (fSize<=value) {
	  fRange=value-fSize;
	  fSize=value;
	} else {
	  fRange=fSize-value;
	  fSize=value;
	}
      }
      // -range
    } else if (argument.CompareTo("-range")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      if ((iResult=ScanSizeArgument(fRange, argv[i]))==-ERANGE) {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
      // -divisor
    } else if (argument.CompareTo("-divisor")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      if ((iResult=ScanSizeArgument(fDivisor, argv[i]))==-ERANGE) {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
      // -decrement
    } else if (argument.CompareTo("-decrement")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      if ((iResult=ScanSizeArgument(fDecrement, argv[i]))==-ERANGE) {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
      // -modulo
    } else if (argument.CompareTo("-modulo")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      if ((iResult=ScanSizeArgument(fModulo, argv[i]))==-ERANGE) {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
      // -offset
    } else if (argument.CompareTo("-offset")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      if ((iResult=ScanSizeArgument(fOffset, argv[i]))==-ERANGE) {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
      // -multiplier
    } else if (argument.CompareTo("-multiplier")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      if ((iResult=ScanFloatArgument(fMultiplier, argv[i]))==-ERANGE) {
	HLTError("wrong parameter for argument %s, number expected", argument.Data());
	iResult=-EINVAL;
      }
    } else {
      if ((iResult=ScanArgument(argc-i, &argv[i]))==-EINVAL) {
	HLTError("unknown argument %s", argument.Data());
	break;
      } else if (iResult==-EPROTO) {
	bMissingParam=1;
	break;
      } else if (iResult>=0) {
	i+=iResult;
	iResult=0;
      }
    }
  }

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  fCurrSize=fSize;

  if (iResult>=0) {
    fpDice=new TRandom;
    if (fpDice) {
      TDatime dt;
      // just take the pointer value as seed combined with time 
      unsigned int seed=0;//(int)(this);
      fpDice->SetSeed(seed^dt.Get());
    } else {
       iResult=-ENOMEM;
    }
  }

  return iResult;
}

int AliHLTDataGenerator::ScanSizeArgument(AliHLTUInt32_t &size, const char* arg)
{
  // see header file for class documentation
  int iResult=0;
  if (arg) {
    TString parameter(arg);
    AliHLTUInt32_t base=1;
    parameter.Remove(TString::kLeading, ' '); // remove all blanks
    if (parameter.EndsWith("k")) {
      base=0x400; // one k
      parameter.Remove(TString::kTrailing, 'k');
    } else if (parameter.EndsWith("M")) {
      base=0x100000; // one M
      parameter.Remove(TString::kTrailing, 'M');
    }
    if (parameter.IsDigit()) {
      size=(AliHLTUInt32_t)parameter.Atoi()*base;
    } else {
      iResult=-ERANGE;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataGenerator::ScanFloatArgument(float &value, const char* arg)
{
  // see header file for class documentation
  int iResult=0;
  if (arg) {
    TString parameter(arg);
    parameter.Remove(TString::kLeading, ' '); // remove all blanks
    if (parameter.IsFloat()) {
      value=parameter.Atof();
    } else {
      iResult=-ERANGE;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataGenerator::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation

  // there are no other arguments than the standard ones
  if (argc==0 && argv==NULL) {
    // this is just to get rid of the warning "unused parameter"
  }
  return -EINVAL;
}

int AliHLTDataGenerator::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  if (fpDice) delete fpDice;
  fpDice=NULL;

  return iResult;
}

int AliHLTDataGenerator::DoEvent( const AliHLTComponentEventData& evtData,
				  const AliHLTComponentBlockData* blocks, 
				  AliHLTComponentTriggerData& /*trigData*/,
				  AliHLTUInt8_t* outputPtr, 
				  AliHLTUInt32_t& size,
				  AliHLTComponentBlockDataList& outputBlocks )
{
  // see header file for class documentation
  int iResult=0;

  AliHLTUInt32_t space=size;
  size=0;
  if (!IsDataEvent()) return 0;

  AliHLTUInt32_t generated=0;
  if (fSize>0) {
    // mode 1: fake independent of input data size
    generated=fpDice->Integer(fRange)+fCurrSize;
    if (fModulo>0 && ((GetEventCount()+1)%fModulo)==0) {
      // manipulate the size
      if (fDivisor>0) {
	fCurrSize/=fDivisor;
	if (fCurrSize==0) fCurrSize=fSize; //reset
      }
      if (fDecrement>0) {
	AliHLTUInt32_t backup=fCurrSize;
	if (fCurrSize<fDecrement) {
	  fCurrSize=fSize; // reset
	} else {
	  fCurrSize-=fDecrement;
	}
	HLTDebug("manipulated output size from %d to %d", backup, fCurrSize);
	backup=0; // just to avoid warning 'unused variable' in the release version
      }
    }

  } else {
    for (unsigned int i=0; i<evtData.fBlockCnt; i++) {
      generated+=blocks[i].fSize;
    }
    generated=(AliHLTUInt32_t)(generated*fMultiplier);
    generated+=fOffset;
  }

  if (generated<=space ) {
    HLTDebug("adding block: size %d", generated);
    AliHLTComponentBlockData bd;
    FillBlockData(bd);
    bd.fPtr=outputPtr;
    bd.fOffset=0;
    bd.fSize=generated;
    bd.fDataType=fDataType;
    bd.fSpecification=fSpecification;
    outputBlocks.push_back(bd);
    size=generated;
  } else {
    iResult=-ENOSPC;
  }

  return iResult;
}
