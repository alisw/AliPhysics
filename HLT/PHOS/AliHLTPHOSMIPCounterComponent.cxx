 /**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSMIPCounterComponent.h"
#include "AliHLTPHOSProcessor.h"
#include "AliHLTPHOSMIPCounter.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TFile.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"


const AliHLTComponentDataType AliHLTPHOSMIPCounterComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}};

AliHLTPHOSMIPCounterComponent gAliHLTPHOSMIPCounterComponent;

AliHLTPHOSMIPCounterComponent::AliHLTPHOSMIPCounterComponent()
    : AliHLTPHOSProcessor(),
    fEvtCnt ( 0 ),
    fInterval ( 0 ),
    fMIPCount ( 0 ),
    fTRUThreshold ( 0 ),
    fMIPCountInterval ( 0 ),
    fPath ( 0 ),
    fMIPCounterPtr ( 0 ),
    fHistPtr ( 0 ),
    fIntervalHistPtr ( 0 ),
    fRateHistPtr ( 0 ),
    fChannelHistPtr ( 0 ),
    fRatioHistPtr ( 0 )
{
}


AliHLTPHOSMIPCounterComponent::~AliHLTPHOSMIPCounterComponent()
{
}

int 
AliHLTPHOSMIPCounterComponent::Deinit()
{
  printf("AliHLTPHOSMIPCounterComponent::Deinit()\n");
  char filename[50];
  sprintf(filename, "%s/MIPCount_TRUThreshold%s.root", fPath, fTRUThreshold);
  TFile *outfile = new TFile(filename, "recreate");
  fHistPtr->Write();
  fIntervalHistPtr->Write();
  fRateHistPtr->Write();
  fChannelHistPtr->Write();
  fRatioHistPtr->Write();
  outfile->Close();
  
  printf("Total number of MIPs in %d events: %d\nGives a rate of: %f\n", fEvtCnt, fMIPCounterPtr->GetMIPCountTotal(),
	 ((float)(fMIPCounterPtr->GetMIPCountTotal()))/((float)fEvtCnt));
      
  return 0;
}

const char*
AliHLTPHOSMIPCounterComponent::GetComponentID()
{
  return "PhosMIPCounter";
}

void
    AliHLTPHOSMIPCounterComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{ 
 //Get datatypes for input
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType); 
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSMIPCounterComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkAliHLTMIPDataType;
}


void 
    AliHLTPHOSMIPCounterComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  constBase = 30;
  inputMultiplier = 1;
}

int 
AliHLTPHOSMIPCounterComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
					AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/, //TODO: I think size should be set to zero when returning from this method if not data was written to the output buffer.
					vector<AliHLTComponentBlockData>& /*outputBlocks*/)
{
   //Do event
  int digitCount = 0;
  
  const AliHLTComponentBlockData* iter = 0; 
  unsigned long ndx; 

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
  {
    iter = blocks+ndx;
      
    if(iter->fDataType != AliHLTPHOSDefinitions::fgkAliHLTDigitDataType)
    {
	//  cout << "Warning: data type is not fgkCellEnergyDataType " << endl;
      continue;
    }
    digitCount += (reinterpret_cast<AliHLTPHOSDigitContainerDataStruct*>(iter->fPtr))->fNDigits;
    fMIPCount += fMIPCounterPtr->CountMIPs(reinterpret_cast<AliHLTPHOSDigitContainerDataStruct*>(iter->fPtr));
  }
  fRatioHistPtr->Fill((float)(((float)fMIPCount)/((float)digitCount)));
  fEvtCnt++;
  fMIPCountInterval += fMIPCount;
  fMIPCount = 0;
  
  if(fEvtCnt % fInterval == 0)
  {
    printf("Event #: %d -- Number of MIPs the last %d events: %d -- Which gives a rate of: %f\n",
	   fEvtCnt, fInterval, Int_t(fMIPCountInterval), ((Float_t)fMIPCountInterval/(Float_t)fInterval));   //TODO: check that the proper things are being written to screen.
    fIntervalHistPtr->Fill(fMIPCountInterval);
    fRateHistPtr->Fill((Float_t)fMIPCountInterval/(Float_t)fInterval);
    fMIPCountInterval = 0;
  }
 
  return 0;
}


int
AliHLTPHOSMIPCounterComponent::DoInit(int argc, const char** argv )
{
  //Do initialization
  cout << "Initializing AliHLTPHOSMIPCounterComponent...\n";
  cout << endl;
  Char_t intervalHistName[50];
  Char_t rateHistName[50];
  fPath = new Char_t[50];
  fTRUThreshold = new Char_t[50];
  fMIPCounterPtr = new AliHLTPHOSMIPCounter();
  for(int i = 0; i < argc; i++)
  {
      if(!strcmp("-interval", argv[i]))
      {
	fInterval = atoi(argv[i+1]);
      }
      if(!strcmp("-path", argv[i]))
      {
	strcpy(fPath, argv[i+1]);
      }
      if(!strcmp("-upperbound", argv[i]))
      {
	fMIPCounterPtr->SetUpperBound(atoi(argv[i+1]));
      }
      if(!strcmp("-lowerbound", argv[i]))
      {
	fMIPCounterPtr->SetLowerBound(atoi(argv[i+1]));
      }
      if(!strcmp("-zerothreshold", argv[i]))
      {
	fMIPCounterPtr->SetZeroThreshold(atoi(argv[i+1]));
      }
      if(!strcmp("-truthreshold", argv[i]))
      {
	strcpy(fTRUThreshold, argv[i+1]);
      }
      if(!strcmp("-lowerstarttime", argv[i]))
      {
	fMIPCounterPtr->SetLowerStartTime(atoi(argv[i+1]));
      }
      if(!strcmp("-upperstarttime", argv[i]))
      {
	fMIPCounterPtr->SetUpperStartTime(atoi(argv[i+1]));
      }
  }
  
  sprintf(intervalHistName, "Number of MIPs in %d Events", fInterval);
  sprintf(rateHistName, "MIP Rate for %d Events", fInterval);
  
  fHistPtr = new TH1I("MIPHist", "Number of MIPs in event", 20, 0, 100);
  fIntervalHistPtr = new TH1I("intervalMIPHist", intervalHistName, 100, 0, 500);
  fRateHistPtr = new TH1F("rateHist", rateHistName, 100, 0, 100);
  fChannelHistPtr = new TH2I("channelHist", "MIP Hits in Channels", 64, 0, 63, 56, 0, 55);
  fRatioHistPtr = new TH1F("ratioHist", "Ratio of MIP digits in Event", 100, 0, 0.2);
     
  fMIPCounterPtr->SetChannelHistogram(fChannelHistPtr);
  

  return 0;  
}

AliHLTComponent*
    AliHLTPHOSMIPCounterComponent::Spawn()
{
  return new AliHLTPHOSMIPCounterComponent();
}
