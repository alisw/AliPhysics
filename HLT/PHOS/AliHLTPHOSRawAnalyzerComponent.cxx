// $Id$

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
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH2F.h"
#include "AliAltroDecoder.h"    // decoder for altro payload
#include "AliAltroData.h"       // container for altro payload
#include "AliAltroBunch.h"      // container for altro bunches


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
								 fDoMakeDigits(false),
								 fDigitMakerPtr(0),
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
  list.clear();
  list.push_back( AliHLTPHOSDefinitions::fgkDDLPackedRawDataType | kAliHLTDataOriginPHOS);
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
  inputMultiplier = 4;
}




void 
AliHLTPHOSRawAnalyzerComponent::FillDataArray(UInt_t *data, const AliAltroData */*altrodataptr*/, const int /*channel*/)
{
  ResetDataPtr(0, ALTRO_MAX_SAMPLES);
  bool islastbunch = true;

  while( fAltroDataPtr->NextBunch(fAltroBunchPtr) == true)
    {
      const UInt_t *tmpdata  = fAltroBunchPtr->GetData();
   
      if(islastbunch == true)
	{
	  data[0] = fAltroBunchPtr->GetEndTimeBin();
	  islastbunch = false;
	}

      int tmpstartbin =  fAltroBunchPtr->GetStartTimeBin();
      int tmpendbin =  fAltroBunchPtr->GetEndTimeBin();
      int tmplength = tmpendbin -  tmpstartbin;

      for(int i = 0; i < tmplength ; i++)
	{ 
	  data[i+tmpstartbin] = tmpdata[i];
	}
    }


  /*
  cout <<__FILE__ <<" : " <<__LINE__  << "the resulting array is"<<endl;
  
  for(int i=0; i<  data[0]; i++)
    {
      if(i != 0 && i %16 == 0)
	{
	  cout << endl;
	}
      cout <<data[i] << "\t" ;
    }
  cout << endl;
  */
}


void 
AliHLTPHOSRawAnalyzerComponent::GetFirstBunch(AliAltroData */*altrodata*/,  AliAltroBunch */*altrobunch*/)
{
  while( fAltroDataPtr->NextBunch(fAltroBunchPtr) == true)
    {
      
    }
}


int 
AliHLTPHOSRawAnalyzerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& /*trigData*/, 
					 AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  //  cout << "Event" << fPhosEventCount  << endl;

  UInt_t offset            = 0; 
  UInt_t mysize            = 0;
  UInt_t tSize             = 0;
  Float_t baseline         = 0;
  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;
  Int_t *rawDataBufferPos = (Int_t *)outputPtr; 
  Int_t nSamples = 0;
  UInt_t specification = 0;
  bool droppedRaw = true;
  if(fDoPushRawData) {droppedRaw = false;}
  
  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if ( iter->fDataType != AliHLTPHOSDefinitions::fgkDDLPackedRawDataType )
	{
	  continue; 
	}

      Int_t tmpChannelCnt     = 0;
      mysize = 0;
      offset = tSize;
      Int_t crazyness = 0;
      mysize += sizeof(AliHLTPHOSRcuCellEnergyDataStruct);
      tSize += mysize;

      if(tSize > size)
	{
	  HLTError("Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes.", tSize, size);
	  return -ENOBUFS;
	}  

      specification = specification|iter->fSpecification;
      fDecoderPtr->SetMemory(reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize);
      fDecoderPtr->Decode();
      fOutPtr =  (AliHLTPHOSRcuCellEnergyDataStruct*)outBPtr;
      fOutPtr->fRcuX = fRcuX;
      fOutPtr->fRcuZ = fRcuZ;
      fOutPtr->fModuleID =fModuleID;
      
      rawDataBufferPos += (tSize)/sizeof(Int_t); 
           
      while( fDecoderPtr->NextChannel(fAltroDataPtr) == true )
	{          
	  FillDataArray(fTmpChannelData, fAltroDataPtr, tmpChannelCnt); 

	  if(  fAltroDataPtr->GetDataSize() != 0 )
            {
	      GetFirstBunch(fAltroDataPtr, fAltroBunchPtr);
	      nSamples = fAltroBunchPtr->GetBunchSize();
	      //	      cout <<__FILE__ <<" : " <<__LINE__  << ",  the size of the first bunch is " << nSamples <<endl;
	      crazyness = fSanityInspectorPtr->CheckInsanity((const UInt_t*)fAltroBunchPtr->GetData(), (const Int_t)(fAltroBunchPtr->GetBunchSize()));
	      fAnalyzerPtr->SetData(fAltroBunchPtr->GetData(), fAltroBunchPtr->GetBunchSize());   
	      fAnalyzerPtr->Evaluate(0, fAltroBunchPtr->GetBunchSize());  
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

	      if(fDoPushRawData == true && droppedRaw == false)
		{
		  int tmpsize = fTmpChannelData[0];
		  //	  cout << __FILE__ << ":" << __LINE__ << "channel = " << tmpChannelCnt <<  " size  ="<< tmpsize << endl;
		  mysize += (tmpsize + 1)*sizeof(Int_t);
		  tSize += (tmpsize  + 1)*sizeof(Int_t);;

		  if(tSize > size)
		    {
		      HLTError("Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes. Dropping raw data.", tSize, size);
		      droppedRaw = true;
		      tSize -= mysize;
		    }
		  else
		    {
		      *rawDataBufferPos = tmpsize;
		      
		      for(int sample = 0; sample < tmpsize; sample++)
			{
			  rawDataBufferPos++;
			  *(rawDataBufferPos) = fTmpChannelData[sample]; 
			}
		      rawDataBufferPos++;
			      
		    }
		}
	      tmpChannelCnt ++;
	    }
	}

      if(fDoPushRawData && droppedRaw == false)
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
      outputBlocks.push_back( bdCellEnergy );
      outBPtr += mysize;

      if( tSize > size )
	{
	  Logging( kHLTLogFatal, "HLT::AliHLTPHOSRawAnalyzerComponent::DoEvent", "Too much data",
		   "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
		   , tSize, size );
	  return EMSGSIZE;
	}
      
    }
  
  //  *rawDataBufferPos = 0;  
  fPhosEventCount++; 
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
      HLTWarning("baseline file no longer supported");

      //      SetBaselines(tmpbaselinfile);
    }
  
  char tmpSelectiveThresholdfile[256];
  if(fUtilitiesPtr->ScanSingleNameArgument(argc, argv, "-baselinefile", tmpSelectiveThresholdfile) == true  )
    {
      fDoSelectiveReadOut = true;
      SetSelectiveReadOutThresholds(tmpSelectiveThresholdfile, nSigmas);
    }

  //fDoSelectiveReadOut = fUtilitiesPtr->ScanSingleNameArgument(argc, argv, "-selectivereadout");

  fDoPushRawData = fUtilitiesPtr->ScanSingleArgument(argc, argv, "-pushrawdata");

  fDoCheckDataSize = fUtilitiesPtr->ScanSingleArgument(argc, argv, "-checkdatasize");
 
  iResult = ScanArguments(argc, argv);
  return iResult;
}


void
AliHLTPHOSRawAnalyzerComponent::Reset(AliHLTPHOSRcuCellEnergyDataStruct* /*cellDataPtr*/)
{
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
AliHLTPHOSRawAnalyzerComponent::SetBaselines(const char* /*file*/)
{
  //comment
//   fUseBaselineSubtraction = true;
//   AliHLTPHOSBaseline *baseline = 0;
//   TFile *baselineFile = TFile::Open(file);
//   TTree *baselineTree = (TTree*)baselineFile->Get("baselineTree");
//   TClonesArray *baselineArray = new TClonesArray("AliHLTPHOSBaseline", 7168);
//   baselineTree->SetBranchAddress("Baselines", &baselineArray);
//   baselineTree->GetEntry(0);
//   for(Int_t i = 0; i < baselineArray->GetEntriesFast(); i++)
//     {
//       baseline = (AliHLTPHOSBaseline*)baselineArray->At(i);
//       if((baseline->GetX() < (Int_t)((fRcuX + 1)*N_XCOLUMNS_RCU)) && (baseline->GetX() >= (Int_t)(fRcuX*N_XCOLUMNS_RCU)))
// 	{
// 	  if((baseline->GetZ() < (Int_t)((fRcuZ + 1)*N_ZROWS_RCU)) && (baseline->GetZ() >= (Int_t)(fRcuZ*N_ZROWS_RCU)))
// 	    {
// 	      fBaselines[baseline->GetX() - fRcuX*N_XCOLUMNS_RCU][baseline->GetZ() - fRcuZ*N_ZROWS_RCU][baseline->GetGain()] = baseline->GetBaseline();
// 	    }
// 	}
//     }
//   baselineFile->Close();
//   delete baselineFile;
//   baselineFile = 0;
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
