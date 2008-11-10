// $Id$

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


/**
 * Measures the baselines and calculates relevant values
 *
 * @file   AliHLTPHOSBaselineAnalyzer.cxx
 * @author Oystein Djuvsland
 * @date
 * @brief  Baseline analyzer for PHOS HLT
 */
#include <fstream>
#include "AliHLTPHOSBaselineAnalyzer.h"
#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSBaseline.h"
#include "AliHLTPHOSValidCellDataStruct.h"
#include "AliHLTPHOSSanityInspector.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSMapper.h"

#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSSharedMemoryInterface.h"
#include "AliHLTPHOSValidCellDataStruct.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace PhosHLTConst;

ClassImp(AliHLTPHOSBaselineAnalyzer); 

AliHLTPHOSBaselineAnalyzer::AliHLTPHOSBaselineAnalyzer() : 
  AliHLTPHOSBase(),
  fSampleStart(5),
  fMaxCrazyDifference(0),	
  fMaxSignal(0),  
  fChannelCount(0),
  fTreePtr(0),
  fBaselineArrayPtr(0),
  fRMSHistogramPtr(0),
  fRMSMapHighGainHistogramPtr(0),
  fRMSMapLowGainHistogramPtr(0),
  fFixedRMSHistogramPtr(0),
  fFixedRMSMapHighGainHistogramPtr(0),
  fFixedRMSMapLowGainHistogramPtr(0),
  fSanityInspector(0)	   
{  
  //see headerfile for documentation
  fSanityInspector = new AliHLTPHOSSanityInspector();
  
  char histName[128];
  char histTitle[128];
  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(int z = 0; z < N_ZROWS_MOD; z++)
    {
      for(int gain = 0; gain < N_GAINS; gain++)
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
  //see headerfile for documentation
} 


void
AliHLTPHOSBaselineAnalyzer::CalculateRcuBaselines(AliHLTPHOSRcuCellEnergyDataStruct* rcuData)
{

  //see headerfile for documentation
  Int_t xOff = rcuData->fRcuX * N_XCOLUMNS_RCU;
  Int_t zOff = rcuData->fRcuZ * N_ZROWS_RCU;
  Float_t baseline = 0;
  AliHLTPHOSValidCellDataStruct *currentChannel =0;

  AliHLTPHOSSharedMemoryInterface* shmPtr = new AliHLTPHOSSharedMemoryInterface();
  shmPtr->SetMemory(rcuData);
  currentChannel = shmPtr->NextChannel();
  
  while(currentChannel != 0)
    {
      if(rcuData->fHasRawData == true)
	{
	  Int_t* rawPtr = 0;
	  Int_t dummysamples = 0;
	  rawPtr = shmPtr->GetRawData(dummysamples);
	  baseline = CalculateChannelBaseline(currentChannel, rawPtr, xOff, zOff);
	  if(!(baseline < 0))
	    {
	      CalculateAccumulatedChannelBaseline(currentChannel->fX + xOff, currentChannel->fZ + zOff, currentChannel->fGain, baseline);
	    }
	}
      currentChannel = shmPtr->NextChannel();
    }    
      
//   AliHLTPHOSValidCellDataStruct *data = rcuData->fValidData;

//   for(Int_t i = 0; i < rcuData->fCnt; i++)
//     {
//       baseline = CalculateChannelBaseline(&(data[i]), xOff, zOff);
//       if(!(baseline < 0))
//       {
// 	CalculateAccumulatedChannelBaseline(data[i].fX + xOff, data[i].fZ + zOff, data[i].fGain, baseline);
//       }
}
  


Float_t
AliHLTPHOSBaselineAnalyzer::CalculateChannelBaseline(AliHLTPHOSValidCellDataStruct *cellData, Int_t* rawDataPtr, Int_t xOff, Int_t zOff) 
{ 
  //see headerfile for documentation
  Int_t crazyness = 0;
  Int_t total = 0;
  Int_t *data = rawDataPtr;
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
    total += data[i];
  }
  fBaselines[cellData->fX + xOff][cellData->fZ + zOff][cellData->fGain] = (float)total / (fNSamples - fSampleStart);

  //if(cellData->fCrazyness == 0)
  // {
    //       crazyness = fSanityInspector->CheckAndHealInsanity(data, fNSamples);
  // }
  if(crazyness < 0)
    return crazyness;

  total = 0;

  for(Int_t j = fSampleStart; j < fNSamples; j++)
    {
      if(data[j] > fMaxSignal)
	return -data[j];
    
      fFixedChannelHistogramsPtr[cellData->fX + xOff][cellData->fZ + zOff][cellData->fGain]->Fill(data[j]);
      total += data[j];
    }
  //  fBaselines[cellData->fX + xOff][cellData->fZ + zOff][cellData->fGain] = (float)total / (fNSamples - fSampleStart);
    
  return (float)total/(fNSamples - fSampleStart);
}

void 
AliHLTPHOSBaselineAnalyzer::CalculateAccumulatedChannelBaseline(Int_t x, Int_t z, Int_t gain, Float_t baseline)
{
  //see headerfile for documentation
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

  //see headerfile for documentation
  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(int z = 0; z < N_ZROWS_MOD; z++)
    {
      for(int gain = 0; gain < N_GAINS; gain++)
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
  //see headerfile for documentation
  fTreePtr = tree;
  fBaselineArrayPtr = baselineArray;
  tree->Branch("Baselines", &fBaselineArrayPtr);
}
  

void
AliHLTPHOSBaselineAnalyzer::FillTree()
{
  //see headerfile for documentation
  AliHLTPHOSBaseline *baseline;
  Int_t n = 0;

  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(int z = 0; z < N_ZROWS_MOD; z++)
	{
	  for(int gain = 0; gain < N_GAINS; gain++)
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
  //see headerfile for documentation
  TFile *baselineFile = TFile::Open(filename, "recreate");
  fTreePtr->Write();
  baselineFile->Close();
}

void 
AliHLTPHOSBaselineAnalyzer::WriteChannelHistograms(const Char_t* filename)
{
  //see headerfile for documentation
  TH1F* tmpRMSHistPtr = new TH1F("tmp", "tmp", 100, 0, 20);
  char fullfilename[128];
  sprintf(fullfilename, "%s.root", filename);
  TFile *file = TFile::Open(fullfilename, "recreate");
  
  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(int z = 0; z < N_ZROWS_MOD; z++)
    {
      for(int gain = 0; gain < N_GAINS; gain++)
      {
	fChannelHistogramsPtr[x][z][gain]->Write();
      }
    }
  }
  file->Close();

  //  char fullfilename[128];
  char header[128];
  int x = -1;
  int z = -1;
  int gain = -1;
  float maxSigma = 0;
  sprintf(fullfilename, "%s.txt", filename);
  ofstream asciifile(fullfilename);
  AliHLTPHOSMapper* fMapperPtr = new AliHLTPHOSMapper();
  //CRAP need to find module number
  for(int module = 0; module < N_MODULES; module++)
    {
      for(int rcu = 0; rcu < N_RCUS; rcu++)
	{
	  for(int branch = 0; branch < N_BRANCHES; branch++)
	    {
	      for(int card = 0; card < N_FEECS; card++)
		{
		  for(int chip = 0; chip < N_ALTROS+1; chip++)
		    {
		      if(chip==1)continue;
		      //		      sprintf(header, "Module:%d RCU:%d Branch:%d Card:%d Chip:%d", module-1, rcu, branch, card, chip);
		      if(module == 2)
			{
			  sprintf(header, "Module:%d RCU:%d Branch:%d Card:%d Chip:%d", module-1, rcu, branch, card+1, chip);
			  asciifile << header << endl;
			  for(int channel = 0; channel < N_ALTROCHANNELS; channel++)
			    {
			      
			      int hwAddress = ((branch << 11) | (card << 7) | (chip << 4) | channel);
			      int xoff = 0;
			      int zoff = 0;
			      if(rcu == 0 || rcu == 2) xoff = 1;
			      if(rcu == 2 || rcu == 3) zoff = 1;
			      z  = fMapperPtr->fHw2geomapPtr[hwAddress].fZRow + zoff*N_ZROWS_RCU;
			      x  = fMapperPtr->fHw2geomapPtr[hwAddress].fXCol + xoff*N_XCOLUMNS_RCU; 
			      gain = fMapperPtr->fHw2geomapPtr[hwAddress].fGain; 
			      if(x >63 || x < 0 || z > 55 ||  z < 0 || gain < 0 ||  gain >1) 
				{
				  //				  cout << "index out of range!" << " x = " << x << " z = " << z << " gain = " << gain << " xoff = " << xoff << " zoff = " << zoff << endl;
 				  continue;
				}
			      //fChannelHistogramsPtr[x][z][gain]->Fit("gaus", "0");
			      //			      asciifile << (fChannelHistogramsPtr[x][z][gain]->GetFunction("gaus"))->GetParameter(1) << endl;
			      asciifile << (int)fChannelHistogramsPtr[x][z][gain]->GetMean()<< endl;
			      //asciifile << (int)(fAccumulatedBaselines[x][z][gain][0]) << endl;
			      if(fChannelHistogramsPtr[x][z][gain]->GetRMS() > maxSigma) maxSigma = fChannelHistogramsPtr[x][z][gain]->GetRMS();
			      tmpRMSHistPtr->Fill(fChannelHistogramsPtr[x][z][gain]->GetRMS());
			    }
			}
		      else
			{
			  //			      asciifile << 0 << endl;
			  maxSigma = 0;
			}
		    }
		}
	    }
	}
    }
  asciifile << tmpRMSHistPtr->GetMean()*3;
}


void
AliHLTPHOSBaselineAnalyzer::WriteRMSHistogram(const Char_t* filename)
{
  //see headerfile for documentation
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
  //see headerfile for documentation
   for(int x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(int z = 0; z < N_ZROWS_MOD; z++)
	{
	  for(int gain = 0; gain < N_GAINS; gain++)
	    {
	      fBaselines[x][z][gain] = 0;    
	    }
	}
    }
}


void 
AliHLTPHOSBaselineAnalyzer::ResetChannelCount()
{
  //see headerfile for documentation
  fChannelCount = 0;
}
    
void 
AliHLTPHOSBaselineAnalyzer::ResetAccumulatedBaselines()
{
   for(int x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(int z = 0; z < N_ZROWS_MOD; z++)
	{
	  for(int gain = 0; gain < N_GAINS; gain++)
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
  //see headerfile for documentation
  fMaxCrazyDifference = diff;  
  fSanityInspector->SetMaxDifference(diff);
}



/*   
void 
AliHLTPHOSBaselineAnalyzer::SetChannelsHistograms(TH1F *channelLowGainHistArray[N_XCOLUMNS_MOD][N_ZROWS_MOD], TH1F *channelHighGainHistArray[N_XCOLUMNS_MOD][N_ZROWS_MOD])
{ 
  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
  {
    for(int z = 0; z < N_ZROWS_MOD; z++)
    {
      for(int gain = 0; gain < N_GAINS; gain++)
      { 
	fChannelLowGainHistogramsPtr[x][z][gain] = channelLowGainHistArray[x][z][gain];
	fChannelHighGainHistogramsPtr[x][z][gain] = channelHighGainHistArray[x][z][gain];
      }
    }
  }
}
*/
