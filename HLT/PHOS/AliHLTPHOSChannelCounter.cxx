
#include "AliHLTPHOSChannelCounter.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSValidCellDataStruct.h"
#include "AliHLTPHOSConstants.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH2D.h"

//ClassImp(AliHLTPHOSChannelCounter);

using namespace PhosHLTConst;

AliHLTPHOSChannelCounter::AliHLTPHOSChannelCounter() :
  AliHLTPHOSBase(),
  // fChannelArrayPtr(0),
  fHistHighGainPtr(0),
  fHistLowGainPtr(0),
  fHistHighRatioPtr(0),
  fHistLowRatioPtr(0)
{
  //  fChannelArrayPtr = new UInt_t[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];
  
  fHistHighGainPtr = new TH2I("highchannelcount", "High gain channel count", 
			      64, 0, 63, 56, 0, 56); 
  fHistLowGainPtr = new TH2I("lowchannelcount", "Low gain channel count", 
			      64, 0, 63, 56, 0, 56); 
  fHistHighRatioPtr = new TH2F("highoutofsync", "High gain channel count divided by number of events", 
			      64, 0, 63, 56, 0, 56); 
  fHistLowRatioPtr = new TH2F("lowoutofsync", "Low gain channel count divided by number of events", 
			     64, 0, 63, 56, 0, 56); 
  
   for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
     {
       for(Int_t z = 0; z < N_ZROWS_MOD; z++)
	 {
	   for(Int_t gain = 0; gain < 2; gain++)
	     {
	       fChannelArrayPtr[x][z][gain] = 0;
	     }
	 }
     }
}

AliHLTPHOSChannelCounter::~AliHLTPHOSChannelCounter()
{
}

void
AliHLTPHOSChannelCounter::CountChannels(AliHLTPHOSRcuCellEnergyDataStruct* channelDataPtr)
{
  Int_t tmp[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];
  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(Int_t z = 0; z < N_ZROWS_MOD; z++)
    {
      for(Int_t gain = 0; gain < 2; gain++)
      {
	tmp[x][z][gain] = 0;
      }
    }
  }
  for(Int_t i = 0; i < channelDataPtr->fCnt; i++)
    {
      AliHLTPHOSValidCellDataStruct *validDataPtr = &(channelDataPtr->fValidData[i]);
      
      tmp[validDataPtr->fX + channelDataPtr->fRcuX*N_XCOLUMNS_RCU]
	  [validDataPtr->fZ + channelDataPtr->fRcuZ*N_ZROWS_RCU]
	  [validDataPtr->fGain]++;
    }
  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
    {
       for(Int_t z = 0; z < N_ZROWS_MOD; z++)
	{
	  for(Int_t gain = 0; gain < N_ZROWS_MOD; gain++)
	  {
	    if(tmp[x][z][gain] > 1)
	      fChannelArrayPtr[x][z][gain] = fChannelArrayPtr[x][z][gain] + 1;
	  }
	}
    }
}

void
AliHLTPHOSChannelCounter::PrintOutOfSyncChannels(Int_t nEvents)
{
  printf("After %d events:\n", nEvents);
  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
  //for(Int_t x = 0; x < 63; x++)
    {
      for(Int_t z = 0; z < N_ZROWS_MOD; z++)
      //for(Int_t z = 0; z < 55; z++)
	{
	  for(Int_t gain = 0; gain < N_GAINS; gain++)
	  //for(Int_t gain = 0; gain < 2; gain++)
	    {
	      printf("x = %d -- z = %d -- gain = %d:   %d\n", x, z, gain, (fChannelArrayPtr[x][z][gain]));
	    }
	}
    }
  printf("\n");
}

void
AliHLTPHOSChannelCounter::FillHistograms(Int_t nEvents)
{
  printf("Filling histograms...");
  //for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
  for(Int_t x = 0; x < 63; x++)
    {
      //  for(Int_t z = 0; z < N_ZROWS_MOD; z++)
      for(Int_t z = 0; z < 55; z++)
	{
	  for(Int_t gain = 0; gain < 2; gain++)
	    {
	      
	      if(gain = 0)
		{
		  fHistHighGainPtr->SetBinContent(x, z, fChannelArrayPtr[x][z][gain]);
		  fHistHighRatioPtr->SetBinContent(x, z, fChannelArrayPtr[x][z][gain]/(float)nEvents);
		  continue;
		}
	      
	      fHistLowGainPtr->SetBinContent(x, z, fChannelArrayPtr[x][z][gain]);
	      fHistLowRatioPtr->SetBinContent(x, z, fChannelArrayPtr[x][z][gain]/(float)nEvents);
	    }
	}
    }
  printf("Done!\n");
}

void 
AliHLTPHOSChannelCounter::WriteHistograms(const char* filename)
{

  printf("Writing histograms to file...");

  TFile *outfile = new TFile(filename,"recreate");

  fHistHighGainPtr->Write();
  fHistLowGainPtr->Write();
  fHistHighRatioPtr->Write();
  fHistLowRatioPtr->Write();

  outfile->Close();

  delete outfile;
  outfile = 0;
  printf("Done!\n");
}
