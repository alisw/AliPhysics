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
//#include "AliHLTPHOSRcuChannelDataStruct.h"
#include "AliHLTPHOSMapper.h"
#include "AliHLTPHOSSanityInspector.h"
#include "AliHLTPHOSBaseline.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include  "AliAltroDecoder.h"    // decoder for altro payload
#include  "AliAltroData.h"       // container for altro payload
#include  "AliAltroBunch.h"      // container for altro bunches


AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent():AliHLTPHOSRcuProcessor(), 
								 fAnalyzerPtr(0), 
								 fSendChannelData(kFALSE),
								 fOutPtr(0), 
								 fMapperPtr(0), 
								 fUseBaselineSubtraction(false), 
								 fDecoderPtr(0),  
								 fAltroDataPtr(0),
								 fAltroBunchPtr(0)   
								 //								 fRawMemoryReader(0), fPHOSRawStream(0) 
{
  //comment
  fMapperPtr = new AliHLTPHOSMapper();
  fAltroDataPtr = new AliAltroData();
  fAltroBunchPtr = new AliAltroBunch();
  fDecoderPtr = new AliAltroDecoder();
  fSanityInspectorPtr = new AliHLTPHOSSanityInspector();
}
 


AliHLTPHOSRawAnalyzerComponent::~AliHLTPHOSRawAnalyzerComponent()
{
  //comment
  delete  fMapperPtr;
}


int 
AliHLTPHOSRawAnalyzerComponent::Deinit()
{
  //comment
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
  return AliHLTPHOSDefinitions::fgkCellEnergyDataType;
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
  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;
  Int_t *rawDataBufferPos = (Int_t *)outputPtr; 

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
      if( fPhosEventCount%100 == 0)
	{
	  cout << "event count = "<< fPhosEventCount <<endl;
	} 
    
      fDecoderPtr->SetMemory(reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize);
      fDecoderPtr->Decode();
      fOutPtr =  (AliHLTPHOSRcuCellEnergyDataStruct*)outBPtr;
      fOutPtr->fRcuX = fRcuX;
      fOutPtr->fRcuZ = fRcuZ;
      fOutPtr->fModuleID =fModuleID;
      rawDataBufferPos += (mysize)/sizeof(Int_t); 

      while( fDecoderPtr->NextChannel(fAltroDataPtr) == true )
	{
	  if(fAltroDataPtr->GetDataSize() != (fNTotalSamples +2))
	    {
	      cout << "Error, fDataSize = " << fAltroDataPtr->GetDataSize() << endl;
	      continue;
	    }
	  else
	    {
	      //	       cout << "Info, fDataSize = " << fAltroDataPtr->GetDataSize() << endl;
	    }

	  crazyness = fSanityInspectorPtr->CheckInsanity(fAltroDataPtr->GetData(), fAltroDataPtr->GetDataSize() - 2);
	  fAnalyzerPtr->SetData(fAltroDataPtr->GetData());  
	  fAnalyzerPtr->Evaluate(0, fAltroDataPtr->GetDataSize() -2);  

	  fOutPtr->fValidData[tmpChannelCnt].fZ  = fMapperPtr->hw2geomapPtr[fAltroDataPtr->GetHadd()].zRow;
	  fOutPtr->fValidData[tmpChannelCnt].fX  = fMapperPtr->hw2geomapPtr[fAltroDataPtr->GetHadd()].xCol; 
	  fOutPtr->fValidData[tmpChannelCnt].fGain  = fMapperPtr->hw2geomapPtr[fAltroDataPtr->GetHadd()].gain; 

	  if(fUseBaselineSubtraction)
	    {
	      baseline = fBaselines[fOutPtr->fValidData[tmpChannelCnt].fX][fOutPtr->fValidData[tmpChannelCnt].fZ][ fOutPtr->fValidData[tmpChannelCnt].fGain];
	    }

	  fOutPtr->fValidData[tmpChannelCnt].fEnergy  = (float)fAnalyzerPtr->GetEnergy() - baseline;
	  fOutPtr->fValidData[tmpChannelCnt].fTime    = (float)fAnalyzerPtr->GetTiming();
	  fOutPtr->fValidData[tmpChannelCnt].fCrazyness = (int)crazyness;
	  fOutPtr->fValidData[tmpChannelCnt].fNSamples = fNTotalSamples;
	  fOutPtr->fValidData[tmpChannelCnt].fData = rawDataBufferPos;
	  const UInt_t *tmpData =  fAltroDataPtr->GetData();

	  for(Int_t sample = 0; sample < fNTotalSamples; sample++)
	    {
	      (fOutPtr->fValidData[tmpChannelCnt].fData)[sample] = tmpData[sample] - (int)baseline;
	    }

	  UInt_t tmpSize =  sizeof(Int_t)*(fOutPtr->fValidData[tmpChannelCnt].fNSamples);
	  mysize += sizeof(Int_t)*(fOutPtr->fValidData[tmpChannelCnt].fNSamples);
	  mysize += tmpSize;
	  rawDataBufferPos += tmpSize/sizeof(Int_t);
	  tmpChannelCnt ++;
	}
      
      fOutPtr->fCnt  =  tmpChannelCnt;
      fOutPtr->fSize = mysize;

      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = offset;
      bd.fSize = mysize;

      bd.fDataType = AliHLTPHOSDefinitions::fgkCellEnergyDataType;
      bd.fSpecification = 0xFFFFFFFF;
      outputBlocks.push_back( bd );

      tSize += mysize;
      outBPtr += mysize;
      
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
	  cout <<"Analyzing event " <<  fPhosEventCount  << "for Equippment " << fkEquippmentID << endl; 
	}  
    }

  return 0;
}//end DoEvent



int
AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv )
{
  //See base class for documentation
  cout <<"AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv ) "<< endl;

  fSendChannelData = kFALSE;
  fPrintInfo = kFALSE;
  Reset();
  int iResult=0;
  TString argument="";
  iResult = ScanArguments(argc, argv);

  for(int i = 0; i < argc; i++)
    {
      if(!strcmp("-baselinefile", argv[i]))
	{
	  cout << "Getting baselines from " << argv[i+1] << endl;
	  SetBaselines(argv[i+1]);
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
AliHLTPHOSRawAnalyzerComponent::Reset()
{
  //comment
  for(unsigned int mod = 0; mod < N_MODULES; mod ++)
    {
      for(unsigned int row = 0; row < N_ZROWS_MOD; row ++)
	{
	  for(unsigned int col = 0; col < N_XCOLUMNS_MOD; col ++)
	    {
	      for(unsigned int gain = 0; gain < N_GAINS; gain ++ )
		{
		  fMaxValues[mod][row][col][gain] = 0;
		}
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
	      //	      cout <<  fBaselines[baseline->GetX() - fRcuX*N_XCOLUMNS_RCU][baseline->GetZ() - fRcuZ*N_ZROWS_RCU][baseline->GetGain()] << endl;
	    }
	}
    }
  baselineFile->Close();
  delete baselineFile;
  baselineFile = 0;
}
