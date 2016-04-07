// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE                    * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

#include "AliHLTFXSWriterComponent.h"
#include "AliHLTReadoutList.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTFXSWriterComponent)

AliHLTFXSWriterComponent::AliHLTFXSWriterComponent()
  : AliHLTCalibrationProcessor()
  , fFXSName("")
  , fFXSDetector("")
  , fDataType(kAliHLTAllDataTypes|kAliHLTDataOriginAny)
  , fRootObject(false)
{
}

AliHLTFXSWriterComponent::~AliHLTFXSWriterComponent()
{
}

void AliHLTFXSWriterComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  list.push_back(kAliHLTAllDataTypes);
}

AliHLTComponentDataType AliHLTFXSWriterComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeFXSCalib;
}

void AliHLTFXSWriterComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  constBase = 0;
  inputMultiplier = 1;
}

void AliHLTFXSWriterComponent::GetOCDBObjectDescription( TMap* const /*targetArray*/)
{
}

int AliHLTFXSWriterComponent::InitCalibration()
{
  int iResult=0;

  return iResult;
}

int AliHLTFXSWriterComponent::DeinitCalibration()
{
  int iResult=0;

  return iResult;
}

Int_t AliHLTFXSWriterComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
    if (fRootObject)
    {
	for (const TObject *obj = GetFirstInputObject(fDataType); obj != NULL; obj = GetNextInputObject())
	{
	    PushToFXS(obj, fFXSDetector, fFXSName, NULL);
	}
    }
    else
    {
	for (const AliHLTComponentBlockData* blk = GetFirstInputBlock(fDataType); blk != NULL; blk = GetNextInputBlock())
	{
	    if (blk->GetDataType() == (kAliHLTAnyDataType | kAliHLTDataOriginPrivate)) continue;
	    PushToFXS(blk->fPtr, blk->fSize, fFXSDetector, fFXSName, NULL);
	}
    }

    return(0);
}

int AliHLTFXSWriterComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/,
						       AliHLTComponentTriggerData& /*trigData*/)
{
  return 0;
}

int AliHLTFXSWriterComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  if (argc<=0) return 0;
  int iRet = 0;
  for( int i=0; i<argc; i++ ){
        TString argument=argv[i];
        if (argument.CompareTo("-FXSName") == 0)
        {
                if (++i >= argc)
                {
                        HLTError("FXSName missing");
                        return(-1);
                }
                fFXSName = argv[i];
                iRet+=2;
        }
        else if (argument.CompareTo("-FXSDetector") == 0)
        {
                if (++i >= argc)
                {
                        HLTError("FXSDetector missing");
                        return(-1);
                }
                fFXSDetector = argv[i];
                iRet+=2;
        }
        else if (argument.CompareTo("-RootObject") == 0)
        {
                fRootObject = true;
                iRet++;
        }
        else if (argument.CompareTo( "-DataType" ) == 0)
        {
                if (i + 2 > argc)
                {
                        HLTError("DataType missing");
                        return(-1);
                }

                fDataType = AliHLTComponentDataTypeInitializerWithPadding(argv[i + 1], argv[i + 2]);
                i += 2;
                iRet += 3;
        }
        else
        {
          iRet = -EINVAL;
          HLTError("Unknown argument %s",argv[i]);
        }
  }
  return iRet;
}
