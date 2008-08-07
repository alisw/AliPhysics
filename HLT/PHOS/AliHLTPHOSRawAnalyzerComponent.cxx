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
                                                                 fDoPushRawData(false),
                                                                 fDigitContainerPtr(0),
                                                                 fDoSelectiveReadOut(false),
                                                                 fSelectedChannelsList(0),
                                                                 fDoCheckDataSize(false),
								 fNCorruptedBlocks(0),
								 fNOKBlocks(0)
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
  Deinit();
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


  //  cout << "AliHLTPHOSRawAnalyzerComponent::DoEven. fPhosEventCount = " <<  fPhosEventCount <<endl;

  //  cout << "AliHLTPHOSRawAnalyzerComponent::DoEvent TP0" << endl;

  UInt_t offset            = 0; 
  UInt_t mysize            = 0;
  UInt_t tSize             = 0;
  Float_t baseline         = 0;
  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;
  Int_t *rawDataBufferPos = (Int_t *)outputPtr; 
  AliHLTPHOSValidCellDataStruct *validCellPtr = 0;
  Int_t nSamples = 0;
  UInt_t nSelected = 0;
  UInt_t specification = 0;

  //  cout << "AliHLTPHOSRawAnalyzerComponent::DoEvent TP1" << endl;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      //    cout << "AliHLTPHOSRawAnalyzerComponent::DoEvent TP2" << endl;     

      Int_t tmpChannelCnt     = 0;
      iter = blocks+ndx;
      mysize = 0;
      offset = tSize;
      Int_t crazyness = 0;
      UInt_t tmpSize = 0;
      mysize += sizeof(AliHLTPHOSRcuCellEnergyDataStruct);
      //     cout << "AliHLTPHOSRawAnalyzerComponent::DoEvent TP3" << endl;     
  
      if ( iter->fDataType != AliHLTPHOSDefinitions::fgkDDLPackedRawDataType )
	{
	  //	  cout << "AliHLTPHOSRawAnalyzerComponent::DoEvent TP4" << endl;     
	  continue; 
	}

      //     cout << "AliHLTPHOSRawAnalyzerComponent::DoEvent TP5" << endl;     
      
      specification = specification|iter->fSpecification;

      fDecoderPtr->SetMemory(reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize);
      fDecoderPtr->Decode();
      fOutPtr =  (AliHLTPHOSRcuCellEnergyDataStruct*)outBPtr;

      fOutPtr->fRcuX = fRcuX;
      fOutPtr->fRcuZ = fRcuZ;
      fOutPtr->fModuleID =fModuleID;

      rawDataBufferPos += (mysize)/sizeof(Int_t); 
            
      while( fDecoderPtr->NextChannel(fAltroDataPtr) == true )
	{

	  nSamples = fAltroDataPtr->GetDataSize() - 2;

	  if(fDoCheckDataSize)
	    {
	      if(nSamples != fNTotalSamples)
		{
		  //		  cout <<"processing event " << fPhosEventCount << endl;;  
		  //		  cout << "Wrong number of samples Expected  "<< fNTotalSamples << " samples (assuming non zero supressed data) but recieved  " << nSamples << endl;
		  Logging( kHLTLogError, __FILE__ , "Wrong number of samples", "Expected  %lu samples (assuming non zero supressed data) but recieved %lu", fNTotalSamples,  nSamples); 
		  fNCorruptedBlocks ++;	  
		  continue;  
		}
	      else
		{
		  //	  cout <<"The number of samples is " << fNTotalSamples << endl;
		}

	    }

	  fNOKBlocks ++;
	  
	  if((fPhosEventCount%10 ==0) && fPhosEventCount !=0)
	    {
	      //	      float percent = ((float)(100*fNCorruptedBlocks))/((float)(fNOKBlocks + fNCorruptedBlocks) );
	    }
	  
	  crazyness = fSanityInspectorPtr->CheckInsanity((const UInt_t*)fAltroDataPtr->GetData(), (const Int_t)(fAltroDataPtr->GetDataSize() - 2));

	  fAnalyzerPtr->SetData(fAltroDataPtr->GetData(), fAltroDataPtr->GetDataSize() -2);   
	  fAnalyzerPtr->Evaluate(0, fAltroDataPtr->GetDataSize() -2);  

	  fOutPtr->fValidData[tmpChannelCnt].fZ  = fMapperPtr->fHw2geomapPtr[fAltroDataPtr->GetHadd()].fZRow;
	  fOutPtr->fValidData[tmpChannelCnt].fX  = fMapperPtr->fHw2geomapPtr[fAltroDataPtr->GetHadd()].fXCol; 
	  fOutPtr->fValidData[tmpChannelCnt].fGain  = fMapperPtr->fHw2geomapPtr[fAltroDataPtr->GetHadd()].fGain; 

	  if(fUseBaselineSubtraction)
	    {
	      baseline = fBaselines[fOutPtr->fValidData[tmpChannelCnt].fX][fOutPtr->fValidData[tmpChannelCnt].fZ][ fOutPtr->fValidData[tmpChannelCnt].fGain];
	    }

	  fOutPtr->fValidData[tmpChannelCnt].fEnergy  = (float)fAnalyzerPtr->GetEnergy() - baseline;
	  fOutPtr->fValidData[tmpChannelCnt].fTime    = (float)fAnalyzerPtr->GetTiming();
	  fOutPtr->fValidData[tmpChannelCnt].fCrazyness = (int)crazyness;
	  
	  const UInt_t *tmpData =  fAltroDataPtr->GetData();

	  if(fDoPushRawData)
	    {
	      tmpSize += nSamples + 1;
	      *rawDataBufferPos = nSamples;
	      //	      cout << "# samples: " << *rawDataBufferPos << endl;

	      for(int sample = 0; sample < nSamples; sample++)
		{
		  rawDataBufferPos++;
		  *(rawDataBufferPos) = tmpData[sample] - (Int_t)baseline;
		}

	      rawDataBufferPos++;
	    }
	  if(fDoSelectiveReadOut)
	    {
	      if(validCellPtr->fEnergy > fSelectiveReadOutThresholds[fOutPtr->fValidData[tmpChannelCnt].fX][fOutPtr->fValidData[tmpChannelCnt].fZ][fOutPtr->fValidData[tmpChannelCnt].fGain])
		{
		  fSelectedChannelsList[nSelected] = (AliHLTUInt16_t)(fAltroDataPtr->GetHadd());
		  nSelected++;
		}
	    }

	  //	  UInt_t tmpSize =  sizeof(Int_t)*(fOutPtr->fValidData[tmpChannelCnt].fNSamples);
	  //	  mysize += sizeof(Int_t)*(fOutPtr->fValidData[tmpChannelCnt].fNSamples);
	  //mysize += tmpSize;
	  //rawDataBufferPos += tmpSize/sizeof(Int_t);
	  
	  
	  tmpChannelCnt ++;
	}

      mysize += tmpSize*4;
      //      cout << "mysize: " << mysize << " - tmpSize: " << 

      if(fDoPushRawData)
	{
	  fOutPtr->fHasRawData = true;
	}
      else
	{
	  fOutPtr->fHasRawData = false;
	}
      
      fOutPtr->fCnt  = tmpChannelCnt;
      fOutPtr->fSize = mysize;

      AliHLTComponentBlockData bdCellEnergy;
      FillBlockData( bdCellEnergy );
      bdCellEnergy.fOffset = offset;
      bdCellEnergy.fSize = mysize;
      bdCellEnergy.fDataType = AliHLTPHOSDefinitions::fgkCellEnergyDataType;
      bdCellEnergy.fSpecification = specification;
      //    cout << "Pushing cell energies" << endl;
      outputBlocks.push_back( bdCellEnergy );
      
      tSize += mysize;
      outBPtr += mysize;
                  
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
	  Logging(kHLTLogBenchmark, __FILE__ , IntToChar(  __LINE__ ) , "Analyzing event %lu", fPhosEventCount);
	}  
 
   }

  return 0;
}//end DoEvent



int
AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv )
{ 
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // fAnalyzerPtr->SetCorrectBaselineUsingFirstFiveSamples();
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  //See base class for documentation
  fSendChannelData = kFALSE;
  fPrintInfo = kFALSE;
  int iResult=0;
  TString argument="";
  const int nSigmas = 3;
  fMapperPtr = new AliHLTPHOSMapper();
 
  if(fMapperPtr->GetIsInitializedMapping() == false)
    {
      Logging(kHLTLogFatal, __FILE__ , IntToChar(  __LINE__ ) , "AliHLTPHOSMapper::Could not initial mapping from file %s, aborting", fMapperPtr->GetFilePath());
      return -4;
    }

  char tmpbaselinfile[256];
  if(fUtilitiesPtr->ScanSingleNameArgument(argc, argv, "-baselinefile", tmpbaselinfile) == true  )
    {
      SetBaselines(tmpbaselinfile);
    }
  
  char tmpSelectiveThresholdfile[256];
  if(fUtilitiesPtr->ScanSingleNameArgument(argc, argv, "-baselinefile", tmpSelectiveThresholdfile) == true  )
    {
      fDoSelectiveReadOut = true;
      SetSelectiveReadOutThresholds(tmpSelectiveThresholdfile, nSigmas);
    }


  //fDoSelectiveReadOut = fUtilitiesPtr->ScanSingleNameArgument(argc, argv, "-selectivereadout");

  fDoPushRawData = fUtilitiesPtr->ScanSingleArgument(argc, argv, "-pushrawdata");

  fDoPushRawData = true; //CRAP

  fDoCheckDataSize = fUtilitiesPtr->ScanSingleArgument(argc, argv, "-checkdatasize");
 
  iResult = ScanArguments(argc, argv);
  return iResult;
}


void
AliHLTPHOSRawAnalyzerComponent::Reset(AliHLTPHOSRcuCellEnergyDataStruct* cellDataPtr)
{
  //comment
  //  for(unsigned int mod = 0; mod < N_MODULES; mod ++)
  //{
  /*
  for(int x = 0; x < N_XCOLUMNS_RCU; x ++)
    {

      for(int z = 0; z < N_ZROWS_RCU; z ++)
	{
	  for(int gain = 0; gain < N_GAINS; gain ++ )
	    {
	      //fMaxValues[mod][row][col][gain] = 0;
	      cellDataPtr->fValidData[x][z][gain].fEnergy = 0;
	      cellDataPtr->fValidData[x][z][gain].fID = -1;
	    }
	}
    }
  */
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

  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(int z = 0; z < N_ZROWS_MOD; z++)
	{
	  fSelectiveReadOutThresholds[x][z][LOW_GAIN] = lgHist->GetBinContent(x, z) * nSigmas;
	  fSelectiveReadOutThresholds[x][z][HIGH_GAIN] = hgHist->GetBinContent(x, z) * nSigmas;
	}
    }
}
