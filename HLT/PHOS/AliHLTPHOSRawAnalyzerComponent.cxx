/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        * 
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                * 
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSRawAnalyzer.h"
#include "AliHLTPHOSRawAnalyzerComponent.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSMapper.h"
#include "AliHLTPHOSSanityInspector.h"
#include "AliHLTPHOSBaseline.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "AliHLTPHOSDigitMaker.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH2F.h"
#include  "AliAltroDecoder.h"    // decoder for altro payload
#include  "AliAltroData.h"       // container for altro payload
#include  "AliAltroBunch.h"      // container for altro bunches


AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent():AliHLTPHOSRcuProcessor(), 
								 fAnalyzerPtr(0), 
								 fSendChannelData(kFALSE),
								 fOutPtr(0), 
								 fMapperPtr(0), 
								 fSanityInspectorPtr(0),
								 fUseBaselineSubtraction(false), 
								 fDecoderPtr(0),  
								 fAltroDataPtr(0),
								 fAltroBunchPtr(0),
								 fDoPushCellEnergies(false),
								 fDoMakeDigits(false),
								 fDigitMakerPtr(0),
								 fDigitContainerPtr(0),
								 fDoSelectiveReadOut(false),
								 fSelectedChannelsList(0),
								 fDoCheckDataSize(false)
								 //fRawMemoryReader(0), fPHOSRawStream(0) 
{
  //comment
  fMapperPtr = new AliHLTPHOSMapper();
  fAltroDataPtr = new AliAltroData();
  fAltroBunchPtr = new AliAltroBunch();
  fDecoderPtr = new AliAltroDecoder();
  fSanityInspectorPtr = new AliHLTPHOSSanityInspector();
  fSelectedChannelsList = new AliHLTUInt16_t[N_XCOLUMNS_RCU*N_ZROWS_RCU*N_GAINS];
}
AliHLTPHOSRawAnalyzerComponent::~AliHLTPHOSRawAnalyzerComponent()
{
  //comment
  if(fMapperPtr)
    {
      delete  fMapperPtr;
      fMapperPtr = 0;
    }
  if(fAltroDataPtr)
    {
      delete fAltroDataPtr;
      fAltroDataPtr = 0;
    }
  if(fAltroBunchPtr)
    {
      delete fAltroBunchPtr;
      fAltroBunchPtr = 0;
    }
  if(fDecoderPtr)
    {
      delete fDecoderPtr;
      fDecoderPtr = 0;
    }
}


int 
AliHLTPHOSRawAnalyzerComponent::Deinit()
{
  //comment
  if(fMapperPtr)
    {
      delete  fMapperPtr;
      fMapperPtr = 0;
    }
  if(fAltroDataPtr)
    {
      delete fAltroDataPtr;
      fAltroDataPtr = 0;
    }
  if(fAltroBunchPtr)
    {
      delete fAltroBunchPtr;
      fAltroBunchPtr = 0;
    }
  if(fDecoderPtr)
    {
      delete fDecoderPtr;
      fDecoderPtr = 0;
    }
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSRawAnalyzerComponen Deinit");
  return 0;
}


const char* 
AliHLTPHOSRawAnalyzerComponent::GetComponentID()
{
  //comment
  return "AliPhosTestRaw";
}

void
AliHLTPHOSRawAnalyzerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //comment
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSRawAnalyzerComponent::GetOutputDataType()
{
  //comment
  return kAliHLTMultipleDataType;
  //  return AliHLTPHOSDefinitions::fgkDigitDataType;
}

int 
AliHLTPHOSRawAnalyzerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // Added by OD
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(AliHLTPHOSDefinitions::fgkCellEnergyDataType);
  tgtList.push_back(AliHLTPHOSDefinitions::fgkDigitDataType);
  tgtList.push_back(kAliHLTDataTypeHwAddr16);
  return tgtList.size();
}

void
AliHLTPHOSRawAnalyzerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  //comment
  constBase = 30;
  inputMultiplier = 1.2;
}



int 
AliHLTPHOSRawAnalyzerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& /*trigData*/, 
					 AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  //comment
  UInt_t offset            = 0; 
  UInt_t mysize            = 0;
  UInt_t tSize             = 0;
  Float_t baseline         = 0;
  UInt_t digitCount        = 0;

  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;
  Int_t *rawDataBufferPos = (Int_t *)outputPtr; 

  AliHLTPHOSValidCellDataStruct *validCellPtr = 0;

  UInt_t nSamples = 0;

  UInt_t nSelected = 0;

  UInt_t specification = 0;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      Int_t tmpChannelCnt     = 0;
      iter = blocks+ndx;
      mysize = 0;
      offset = tSize;
      Int_t crazyness = 0;
      mysize += sizeof(AliHLTPHOSRcuCellEnergyDataStruct);

      if ( iter->fDataType != AliHLTPHOSDefinitions::fgkDDLPackedRawDataType )
	{
	  continue; 
	}

      specification = specification|iter->fSpecification;
  
      fDecoderPtr->SetMemory(reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize);
      fDecoderPtr->Decode();
      fOutPtr =  (AliHLTPHOSRcuCellEnergyDataStruct*)outBPtr;
      fOutPtr->fRcuX = fRcuX;
      fOutPtr->fRcuZ = fRcuZ;
      fOutPtr->fModuleID =fModuleID;
      Reset(fOutPtr);
      rawDataBufferPos += (mysize)/sizeof(Int_t); 

      while( fDecoderPtr->NextChannel(fAltroDataPtr) == true )
	{
	  
	  nSamples = fAltroDataPtr->GetDataSize() - 2;
	  if(fDoCheckDataSize)
	    {
	      if(nSamples != fNTotalSamples)
		{
		  //cout << "Error, fDataSize = " << nSamples + 2 << endl;
		  //continue;
		}
	      else
		{
		  //cout << "Info, fDataSize = " << fAltroDataPtr->GetDataSize() << endl;
		}
	    }
	 
	  crazyness = fSanityInspectorPtr->CheckInsanity(fAltroDataPtr->GetData(), fAltroDataPtr->GetDataSize() - 2);
	  fAnalyzerPtr->SetData(fAltroDataPtr->GetData());  
	  fAnalyzerPtr->Evaluate(0, fAltroDataPtr->GetDataSize() -2);  

	  Int_t x = fMapperPtr->fHw2geomapPtr[fAltroDataPtr->GetHadd()].fXCol;
	  Int_t z = fMapperPtr->fHw2geomapPtr[fAltroDataPtr->GetHadd()].fZRow;
	  Int_t gain = fMapperPtr->fHw2geomapPtr[fAltroDataPtr->GetHadd()].fGain;	  
	  validCellPtr = &(fOutPtr->fValidData[x][z][gain]);
	  validCellPtr->fX = x;
	  validCellPtr->fZ = z;
	  validCellPtr->fGain = gain;

	  if(fUseBaselineSubtraction)
	    {
	      baseline = fBaselines[validCellPtr->fX][validCellPtr->fZ][ validCellPtr->fGain];
	    }
	  baseline = 0;
	  
	  validCellPtr->fEnergy  = (float)fAnalyzerPtr->GetEnergy() - baseline;
	  validCellPtr->fTime    = (float)fAnalyzerPtr->GetTiming();
	  validCellPtr->fCrazyness = (int)crazyness;
	  validCellPtr->fNSamples = nSamples;
	  validCellPtr->fID = tmpChannelCnt;
	  validCellPtr->fData = rawDataBufferPos;
	  const UInt_t *tmpData =  fAltroDataPtr->GetData();

	  if(fDoPushCellEnergies)
	    {
	      for(UInt_t sample = 0; sample < nSamples; sample++)
		{
		  (validCellPtr->fData)[sample] = tmpData[sample] - (int)baseline;
		}
	    }
	  if(fDoSelectiveReadOut)
	    {
	      if(validCellPtr->fEnergy > fSelectiveReadOutThresholds[x][z][gain])
		{
		  fSelectedChannelsList[nSelected] = (AliHLTUInt16_t)(fAltroDataPtr->GetHadd());
		  nSelected++;
		}
	    }

	  UInt_t tmpSize =  sizeof(Int_t)*(validCellPtr->fNSamples);
	  //	  mysize += sizeof(Int_t)*(validCellPtr->fNSamples);
	  mysize += tmpSize;
	  //	  mysize += tmpSize;
	  rawDataBufferPos += tmpSize/sizeof(Int_t);

	  tmpChannelCnt ++;

	}
      
      fOutPtr->fCnt  = tmpChannelCnt;
      fOutPtr->fSize = mysize;
      
      if(fDoPushCellEnergies)
	{
	  AliHLTComponentBlockData bdCellEnergy;
	  FillBlockData( bdCellEnergy );
	  bdCellEnergy.fOffset = offset;
	  bdCellEnergy.fSize = mysize;
	  bdCellEnergy.fDataType = AliHLTPHOSDefinitions::fgkCellEnergyDataType;
	  bdCellEnergy.fSpecification = specification;
	  outputBlocks.push_back( bdCellEnergy );
	  
	  tSize += mysize;
	  outBPtr += mysize;
	}
      
      //Making Digits
      if(fDoMakeDigits)
	{
	  Int_t digitSize = 0;
	  fDigitMakerPtr->SetDigitContainerStruct((AliHLTPHOSDigitContainerDataStruct*)outBPtr);	  
	  digitCount = fDigitMakerPtr->MakeDigits(fOutPtr);
	  offset = tSize;
	  digitSize += sizeof(AliHLTPHOSDigitContainerDataStruct);
	  
	  AliHLTComponentBlockData bdDigits;
	  FillBlockData(bdDigits);
	  bdDigits.fOffset = offset;
	  bdDigits.fSize = mysize;
	  bdDigits.fDataType = AliHLTPHOSDefinitions::fgkDigitDataType;
	  bdDigits.fSpecification = specification;
	  outputBlocks.push_back( bdDigits );
	  
	  tSize += digitSize;
	  outBPtr += digitSize;
	  fDigitMakerPtr->Reset();
	}      
                  
      //Pushing selected channel addresses
      if(fDoSelectiveReadOut)
	{
	  UInt_t hwAddSize = sizeof(AliHLTUInt16_t);
	  offset = tSize;
	  for(UInt_t n = 0; n < nSelected; n++)
	    {
	      ((AliHLTUInt16_t*)outBPtr)[n] = fSelectedChannelsList[n];
	    }
	  mysize = nSelected*hwAddSize;
	  AliHLTComponentBlockData bdHwAdd;
	  FillBlockData(bdHwAdd);
	  bdHwAdd.fOffset = offset;
	  bdHwAdd.fSize = mysize;
	  bdHwAdd.fDataType = kAliHLTDataTypeHwAddr16;
	  bdHwAdd.fSpecification = specification;
	  outputBlocks.push_back( bdHwAdd );
	  
	  tSize += mysize;
	  outBPtr += mysize;
	}
      
      if( tSize > size )
      	{
	  cout <<"kHLTLogFatal, HLT::AliHLTPHOSRawAnalyzerComponent::DoEvent Too much dataData written over allowed buffer. Amount written:"
	       << tSize << " allowed" << size << endl;
      	  Logging( kHLTLogFatal, "HLT::AliHLTPHOSRawAnalyzerComponent::DoEvent", "Too much data",
      		   "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
      		   , tSize, size );
      	  return EMSGSIZE;
      	}
	
    }

  fPhosEventCount++; 

  if(fPrintInfo == kTRUE)
    {
      if(fPhosEventCount%fPrintInfoFrequncy == 0)
      	{
	  cout << "Analyzing event " <<  fPhosEventCount  << " for Equippment " << fkEquippmentID << endl; 
	  if(fDoSelectiveReadOut) cout << "# of selected channels: " << nSelected << endl;
	  if(fDoMakeDigits) cout << "# of digits: " << digitCount << endl;
	}  
    }

  return 0;
}//end DoEvent



int
AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv )
{
  //See base class for documentation
  
  fSendChannelData = kFALSE;
  fPrintInfo = kFALSE;
  int iResult=0;
  TString argument="";
  iResult = ScanArguments(argc, argv);

  int nSigmas = 3;
  
  fDigitMakerPtr = new AliHLTPHOSDigitMaker();
  
  for(int i = 0; i < argc; i++)
    {
      if(!strcmp("-rmsfilepath", argv[i]))
	{
	  fDigitMakerPtr->SetDigitThresholds(argv[i+1], nSigmas);
	  SetSelectiveReadOutThresholds(argv[i+1], nSigmas);
	}
      if(!strcmp("-baselinefile", argv[i]))
	{
	  SetBaselines(argv[i+1]);
	}
      if(!strcmp("-lowgainfactor", argv[i]))
	{
	  fDigitMakerPtr->SetGlobalLowGainFactor(atof(argv[i+1]));
	}
      if(!strcmp("-highgainfactor", argv[i]))
	{
	  fDigitMakerPtr->SetGlobalHighGainFactor(atof(argv[i+1]));
	}
      if(!strcmp("-selectivereadout", argv[i]))
	{
	  fDoSelectiveReadOut = true;
	}
      if(!strcmp("-makedigits", argv[i]))
	{
	  fDoMakeDigits = true;
	}
      if(!strcmp("-pushcellenergies", argv[i]))
       	{
	  fDoPushCellEnergies = true;
 	}
      if(!strcmp("-checkdatasize", argv[i]))
	{
	  fDoCheckDataSize = true;
	}
    }

  if(fIsSetEquippmentID == kFALSE)
    {
      cout << "The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>" << endl;
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuHistogramProducerComponent::DoInt( int argc, const char** argv )", "Missing argument",
	       "The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>");
      iResult = -3; 
    }
  else
    {
      iResult = 0; 
      //      fRawMemoryReader->SetEquipmentID(fkEquippmentID);
    }

  return iResult;
}


void
AliHLTPHOSRawAnalyzerComponent::Reset(AliHLTPHOSRcuCellEnergyDataStruct* cellDataPtr)
{
  //comment
  //  for(unsigned int mod = 0; mod < N_MODULES; mod ++)
  //{
  for(unsigned int x = 0; x < N_XCOLUMNS_RCU; x ++)
    {
      for(unsigned int z = 0; z < N_ZROWS_RCU; z ++)
	{
	  for(unsigned int gain = 0; gain < N_GAINS; gain ++ )
	    {
	      //fMaxValues[mod][row][col][gain] = 0;
	      cellDataPtr->fValidData[x][z][gain].fEnergy = 0;
	      cellDataPtr->fValidData[x][z][gain].fID = -1;
	    }
	}
    }
  
  ResetDataPtr(0, ALTRO_MAX_SAMPLES);

} // end Reset


void
AliHLTPHOSRawAnalyzerComponent::ResetDataPtr(int startindex, int sampleCnt)
{
  //comment
  for(int i = startindex ; i< sampleCnt; i++)
    {
      fTmpChannelData[i] = 0;
    }
}


void 
AliHLTPHOSRawAnalyzerComponent::SetBaselines(const char* file)
{
  //comment
  fUseBaselineSubtraction = true;
  AliHLTPHOSBaseline *baseline = 0;
  TFile *baselineFile = TFile::Open(file);
  TTree *baselineTree = (TTree*)baselineFile->Get("baselineTree");
  TClonesArray *baselineArray = new TClonesArray("AliHLTPHOSBaseline", 7168);
  baselineTree->SetBranchAddress("Baselines", &baselineArray);
  baselineTree->GetEntry(0);
  for(Int_t i = 0; i < baselineArray->GetEntriesFast(); i++)
    {
      baseline = (AliHLTPHOSBaseline*)baselineArray->At(i);
      if((baseline->GetX() < (Int_t)((fRcuX + 1)*N_XCOLUMNS_RCU)) && (baseline->GetX() >= (Int_t)(fRcuX*N_XCOLUMNS_RCU)))
	{
	  if((baseline->GetZ() < (Int_t)((fRcuZ + 1)*N_ZROWS_RCU)) && (baseline->GetZ() >= (Int_t)(fRcuZ*N_ZROWS_RCU)))
	    {
	      fBaselines[baseline->GetX() - fRcuX*N_XCOLUMNS_RCU][baseline->GetZ() - fRcuZ*N_ZROWS_RCU][baseline->GetGain()] = baseline->GetBaseline();
	    }
	}
    }
  baselineFile->Close();
  delete baselineFile;
  baselineFile = 0;
}

void 
AliHLTPHOSRawAnalyzerComponent::SetSelectiveReadOutThresholds(const char* filepath, Int_t nSigmas)
{
  //See header file for documentation

  TFile *histFile = new TFile(filepath);
  
  TH2F *lgHist = (TH2F*)histFile->Get("RMSLGMapHist");
  TH2F *hgHist = (TH2F*)histFile->Get("RMSHGMapHist");

  for(UInt_t x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(UInt_t z = 0; z < N_ZROWS_MOD; z++)
	{
	  fSelectiveReadOutThresholds[x][z][LOW_GAIN] = lgHist->GetBinContent(x, z) * nSigmas;
	  fSelectiveReadOutThresholds[x][z][HIGH_GAIN] = hgHist->GetBinContent(x, z) * nSigmas;
	}
    }
}
