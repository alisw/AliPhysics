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

#include "AliHLTPHOSBaselineAnalyzer.h"
#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSBaseline.h"
#include "AliHLTPHOSValidCellDataStruct.h"
#include "AliHLTPHOSSanityInspector.h"
#include "AliHLTPHOSConstants.h"

#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace PhosHLTConst;

ClassImp(AliHLTPHOSBaselineAnalyzer); 

AliHLTPHOSBaselineAnalyzer::AliHLTPHOSBaselineAnalyzer() : 
  AliHLTPHOSBase(),
  fSampleStart(5),
  fMaxSignal(0),
  fMaxCrazyDifference(0),	    
  fTreePtr(0),
  fSanityInspector(0)	   
{  
  //comment
  fSanityInspector = new AliHLTPHOSSanityInspector();
  
  char histName[128];
  char histTitle[128];
  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(Int_t z = 0; z < N_ZROWS_MOD; z++)
    {
      for(Int_t gain = 0; gain < N_GAINS; gain++)
      { 
	sprintf(histName, "sample_value_channel_x_col_%d_z_row_%d_gain_%d", x, z, gain);
	sprintf(histTitle, "Distribution of Sample Values for Channel X: %d - Z: %d - Gain: %d", x, z, gain);
	fChannelHistogramsPtr[x][z][gain] = new TH1F(histName, histTitle, 1024, 0, 1023);
	sprintf(histName, "fixed_sample_value_channel_x_col_%d_z_row_%d_gain_%d", x, z, gain);
	sprintf(histTitle, "Distribution of Corrected Sample Values for Channel X: %d - Z: %d - Gain: %d", x, z, gain);
	fFixedChannelHistogramsPtr[x][z][gain] = new TH1F(histName, histTitle, 1024, 0, 1023);
      }
    }
  }
  
  fRMSHistogramPtr = new TH1F("RMSHist", "RMS Values for All Channels", 1023, 0, 1023);
  fRMSMapHighGainHistogramPtr = new TH2F("RMSHGMapHist", "Map of RMS Values for High Gain Channels", 64, 0, 63, 56, 0, 55);
  fRMSMapLowGainHistogramPtr = new TH2F("RMSLGMapHist", "Map of RMS Values for Low Gain Channels", 64, 0, 63, 56, 0, 55);
  
  fFixedRMSHistogramPtr = new TH1F("fixedRMSHist", "Corrected RMS Values for All Channels", 1023, 0, 1023);
  fFixedRMSMapHighGainHistogramPtr = new TH2F("fixedRMSHGMapHist", "Map of Corrected RMS Values for High Gain Channels", 64, 0, 63, 56, 0, 55);
  fFixedRMSMapLowGainHistogramPtr = new TH2F("fixedRMSLGMapHist", "Map of Corrected RMS Values for Low Gain Channels", 64, 0, 63, 56, 0, 55);
  
  ResetBaselines();
  ResetAccumulatedBaselines();
} 

AliHLTPHOSBaselineAnalyzer::~AliHLTPHOSBaselineAnalyzer() 
{ 
  //comment
} 


void
AliHLTPHOSBaselineAnalyzer::CalculateRcuBaselines(AliHLTPHOSRcuCellEnergyDataStruct* rcuData)
{
  //comment
  Int_t xOff = rcuData->fRcuX * N_XCOLUMNS_RCU;
  Int_t zOff = rcuData->fRcuZ * N_ZROWS_RCU;
  Float_t baseline = 0;

  AliHLTPHOSValidCellDataStruct *data = rcuData->fValidData;

  for(Int_t i = 0; i < rcuData->fCnt; i++)
    {
      baseline = CalculateChannelBaseline(&(data[i]), xOff, zOff);
      if(!(baseline < 0))
      {
	CalculateAccumulatedChannelBaseline(data[i].fX + xOff, data[i].fZ + zOff, data[i].fGain, baseline);
      }
    }
}
  


Float_t
AliHLTPHOSBaselineAnalyzer::CalculateChannelBaseline(AliHLTPHOSValidCellDataStruct *cellData, Int_t xOff, Int_t zOff) 
{ 
  //comment
  Int_t crazyness = 0;
  Int_t total = 0;
  Int_t *data = cellData->fData;
  for(Int_t i = fSampleStart; i < fNSamples; i++)
  { 
    /*
    if(cellData->fGain == 0)
    {
      fChannelLowGainHistogramsPtr[cellData->fX + xOff][cellData->fZ + zOff]->Fill(data[i]);
    }
    else
    {
      fChannelHighGainHistogramsPtr[cellData->fX + xOff][cellData->fZ + zOff]->Fill(data[i]);
    }
    */
    fChannelHistogramsPtr[cellData->fX + xOff][cellData->fZ + zOff][cellData->fGain]->Fill(data[i]);
  }
  
  if(cellData->fCrazyness == 0)
  {
       crazyness = fSanityInspector->CheckAndHealInsanity(data, fNSamples);
  }
  if(crazyness < 0)
    return crazyness;
  
  for(Int_t j = fSampleStart; j < fNSamples; j++)
    {
      if(data[j] > fMaxSignal)
	return -data[j];
    
      fFixedChannelHistogramsPtr[cellData->fX + xOff][cellData->fZ + zOff][cellData->fGain]->Fill(data[j]);
      total += data[j];
    }
  fBaselines[cellData->fX + xOff][cellData->fZ + zOff][cellData->fGain] = (float)total / (fNSamples - fSampleStart);
    
  return (float)total/(fNSamples - fSampleStart);
}

void 
AliHLTPHOSBaselineAnalyzer::CalculateAccumulatedChannelBaseline(Int_t x, Int_t z, Int_t gain, Float_t baseline)
{
  //comment
  Float_t oldBaseline = fAccumulatedBaselines[x][z][gain][0];
  Float_t nEntries = fAccumulatedBaselines[x][z][gain][1];
  
  Float_t total = nEntries * oldBaseline + baseline;
  nEntries++;
  fAccumulatedBaselines[x][z][gain][0] = total / nEntries;
  fAccumulatedBaselines[x][z][gain][1] = nEntries;
}

void
AliHLTPHOSBaselineAnalyzer::CalculateChannelsBaselineRMS()
{
  //comment
  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(Int_t z = 0; z < N_ZROWS_MOD; z++)
    {
      for(Int_t gain = 0; gain < N_GAINS; gain++)
      {
	fRMSHistogramPtr->Fill(fChannelHistogramsPtr[x][z][gain]->GetRMS());  
	if(gain == 1)
	  fRMSMapHighGainHistogramPtr->SetBinContent(x, z, fChannelHistogramsPtr[x][z][gain]->GetRMS());
	else
	  fRMSMapLowGainHistogramPtr->SetBinContent(x, z, fChannelHistogramsPtr[x][z][gain]->GetRMS());
		   
	fFixedRMSHistogramPtr->Fill(fFixedChannelHistogramsPtr[x][z][gain]->GetRMS());  
	if(gain == 1)
	  fFixedRMSMapHighGainHistogramPtr->SetBinContent(x, z, fFixedChannelHistogramsPtr[x][z][gain]->GetRMS());
	else
	  fFixedRMSMapLowGainHistogramPtr->SetBinContent(x, z, fFixedChannelHistogramsPtr[x][z][gain]->GetRMS());
      }
    }
  }
}    
 

void 
AliHLTPHOSBaselineAnalyzer::SetRootObjects(TTree *tree, TClonesArray *baselineArray)
{
  //comment
  fTreePtr = tree;
  fBaselineArrayPtr = baselineArray;
  tree->Branch("Baselines", &fBaselineArrayPtr);
}
  

void
AliHLTPHOSBaselineAnalyzer::FillTree()
{
  //comment
  AliHLTPHOSBaseline *baseline;
  Int_t n = 0;

  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(Int_t z = 0; z < N_ZROWS_MOD; z++)
	{
	  for(Int_t gain = 0; gain < N_GAINS; gain++)
	    {
	      baseline = (AliHLTPHOSBaseline*)fBaselineArrayPtr->New(n);
	      baseline->SetX(x);
	      baseline->SetZ(z);
	      baseline->SetGain(gain);
	      baseline->SetBaseline(fAccumulatedBaselines[x][z][gain][0]);
	      baseline->SetEntries((Int_t)fAccumulatedBaselines[x][z][gain][1]);
	      n++;
	    }
	}
    }
  fTreePtr->Fill();
  fBaselineArrayPtr->Clear();
}

void
AliHLTPHOSBaselineAnalyzer::WriteAccumulatedBaselines(const Char_t* filename)
{
  //comment
  TFile *baselineFile = TFile::Open(filename, "recreate");
  fTreePtr->Write();
  baselineFile->Close();
}

void 
AliHLTPHOSBaselineAnalyzer::WriteChannelHistograms(const Char_t* filename)
{
  //comment
  TFile *file = TFile::Open(filename, "recreate");
  
  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(Int_t z = 0; z < N_ZROWS_MOD; z++)
    {
      for(Int_t gain = 0; gain < N_GAINS; gain++)
      {
	fChannelHistogramsPtr[x][z][gain]->Write();
      }
    }
  }
  file->Close();
}


void
AliHLTPHOSBaselineAnalyzer::WriteRMSHistogram(const Char_t* filename)
{
  //comment
  TFile *file = TFile::Open(filename, "recreate");
  fRMSHistogramPtr->Write();
  fRMSMapHighGainHistogramPtr->Write();
  fRMSMapLowGainHistogramPtr->Write();

  fFixedRMSHistogramPtr->Write();
  fFixedRMSMapHighGainHistogramPtr->Write();
  fFixedRMSMapLowGainHistogramPtr->Write();
  
  file->Close();
}

void 
AliHLTPHOSBaselineAnalyzer::ResetBaselines()
{
  //comment
   for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(Int_t z = 0; z < N_ZROWS_MOD; z++)
	{
	  for(Int_t gain = 0; gain < N_GAINS; gain++)
	    {
	      fBaselines[x][z][gain] = 0;    
	    }
	}
    }
}

void 
AliHLTPHOSBaselineAnalyzer::ResetChannelCount()
{
  //comment
  fChannelCount = 0;
}
    
void 
AliHLTPHOSBaselineAnalyzer::ResetAccumulatedBaselines()
{
   for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(Int_t z = 0; z < N_ZROWS_MOD; z++)
	{
	  for(Int_t gain = 0; gain < N_GAINS; gain++)
	    {
	     fAccumulatedBaselines[x][z][gain][0] = 0;
	     fAccumulatedBaselines[x][z][gain][1] = 0;
	    }
	}
    }
}

void 
AliHLTPHOSBaselineAnalyzer::SetMaxCrazyDifference(Int_t diff)
{
  //comment
  fMaxCrazyDifference = diff;  
  fSanityInspector->SetMaxDifference(diff);
}



/*   
void 
AliHLTPHOSBaselineAnalyzer::SetChannelsHistograms(TH1F *channelLowGainHistArray[N_XCOLUMNS_MOD][N_ZROWS_MOD], TH1F *channelHighGainHistArray[N_XCOLUMNS_MOD][N_ZROWS_MOD])
{ 
  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(Int_t z = 0; z < N_ZROWS_MOD; z++)
    {
      for(Int_t gain = 0; gain < N_GAINS; gain++)
      { 
	fChannelLowGainHistogramsPtr[x][z][gain] = channelLowGainHistArray[x][z][gain];
	fChannelHighGainHistogramsPtr[x][z][gain] = channelHighGainHistArray[x][z][gain];
      }
    }
  }
}
*/
