// $Id$

/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Per Thomas Hille for the ALICE                                *
 * offline/HLT Project. Contributors are mentioned in the code where      *
 * appropriate.                                                           *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSOnlineDisplayFourierTab.h"
#include <iostream>
#include "TGFrame.h"
#include "AliHLTPHOSGetEventButton.h"
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTDataTypes.h"
#include "AliHLTHOMERData.h"
#include "AliHLTHOMERReader.h"
#include "AliHLTHOMERWriter.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
//#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h" 
#include "AliHLTPHOSOnlineDisplay.h"
#include "AliHLTPHOSSharedMemoryInterface.h"
#include "AliHLTPHOSFourier.h"
#include "AliHLTPHOSOnlineDisplayTH2D.h"


#include "AliHLTPHOSRcuFFTDataStruct.h"
#include "TStyle.h"

#define SAMPLINGFREQUENCY 10

using namespace std;

// MT Crap
#include <TMath.h>
//#include <TEveManager.h>
//#include <TEveBoxSet.h>

//TEveBoxSet* gAliEveBoxSet = 0;

AliHLTPHOSOnlineDisplayFourierTab::AliHLTPHOSOnlineDisplayFourierTab()
{
  // See header file for documentation
  cout << "ERROR: You cannot create a onlinedisplay Tab without arguments" << endl;
}


AliHLTPHOSOnlineDisplayFourierTab::AliHLTPHOSOnlineDisplayFourierTab(AliHLTPHOSOnlineDisplay * const onlineDisplayPtr, TGTab  *tabPtr, 
								     const AliHLTHOMERReader * homerSyncPtr, const AliHLTHOMERReader * const homerPtrs[MAXHOSTS], int nHosts) :  AliHLTPHOSOnlineDisplayTab(), fEvtCnt(0)
{     
  // See header file for documentation
  // gStyle->SetOptLogy();
  // gStyle->SetOptStat(false);
  

  fShmPtr = new AliHLTPHOSSharedMemoryInterface();
  fOnlineDisplayPtr =  onlineDisplayPtr;
  fFourierPtr = new AliHLTPHOSFourier();

  for(int gain = 0; gain < NGAINS; gain ++ )
    {
      fFourierHistoNew[gain] = 0;
      fFourierHistoOld[gain] = 0;
      fFourierHistoAccumulated[gain] = 0;
    }

  for(int i=0; i<MAXHOSTS; i++)
    {
       fgHomerReadersPtr[i] = 0;
    }

  fgHomerReaderPtr = const_cast<AliHLTHOMERReader*>(homerSyncPtr);
  
  for(int i=0; i<nHosts; i++)
    {
      fgHomerReadersPtr[i] = const_cast<AliHLTHOMERReader*>(homerPtrs[i]);

    }

  fgNHosts = nHosts;
  InitDisplay(tabPtr);
}


AliHLTPHOSOnlineDisplayFourierTab::~AliHLTPHOSOnlineDisplayFourierTab()

{
  // See header file for documentation
}



int
AliHLTPHOSOnlineDisplayFourierTab::GetNextEvent()
{
  // See header file for documentation
  //  ResetDisplay();
  DoGetNextEvent();
  //  FillHistograms();
  UpdateDisplay();
  fEvtCnt ++;
  // fgEvntCnt ++;

 
}



void 
AliHLTPHOSOnlineDisplayFourierTab::ReadBlockData(AliHLTHOMERReader * const homeReaderPtr)
{ 
  // See header file for documentation
  AliHLTPHOSValidCellDataStruct *currentChannel =0;
  cout << "AliHLTPHOSOnlineDisplayFourierTab::ReadBlockDat, Reading block data, therere are " <<  homeReaderPtr->GetBlockCnt() << " blocks " <<endl;
  //  unsigned long blk = homeReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF );

//   while ( blk != ~(unsigned long)0 ) 
//     {
//       Int_t moduleID;
//       Int_t rcuX = 0;
//       Int_t rcuZ = 0;
//       AliHLTPHOSRcuCellEnergyDataStruct* cellEnergiesPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)homeReaderPtr->GetBlockData( blk ); 
      
//       unsigned int *t = (unsigned int*)cellEnergiesPtr;
      
//       moduleID = cellEnergiesPtr->fModuleID ;
//       rcuX = cellEnergiesPtr->fRcuX;
//       rcuZ = cellEnergiesPtr->fRcuZ;

//       cout << "AliHLTPHOSOnlineDisplayFourierTab::ReadBlockData,  fModuleID =" <<moduleID << endl; 

//       Int_t tmpZ;
//       Int_t tmpX;
//       Int_t tmpGain;
//       int cnt = 0;
//       Int_t* tmpPtr = 0;

//       fShmPtr->SetMemory(cellEnergiesPtr);
//       currentChannel = fShmPtr->NextChannel();

//       while(currentChannel != 0)
// 	{
// 	  cnt ++;
// 	  tmpZ = currentChannel->fZ;
// 	  tmpX = currentChannel->fX;
// 	  tmpGain =  currentChannel->fGain;

// 	  if(cellEnergiesPtr->fHasRawData == true)
// 	    {
// 	      Int_t nSamples = 0;
// 	      Int_t* rawPtr = 0;
// 	      rawPtr = fShmPtr->GetRawData(nSamples);
// 	      fFourierPtr->ProcessFourier(rawPtr, nSamples, tmpZ, tmpX, tmpGain, fEvtCnt);
// 	    }
	  
// 	  currentChannel = fShmPtr->NextChannel();
// 	}
//       blk = homeReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF, blk+1);
//     }

//  FillHistograms(fFourierPtr->GetPSD(), fFourierPtr->GetDataSize());

  unsigned long blk = homeReaderPtr->FindBlockNdx("TREIRUOF","SOHP", 0x0);
  while ( blk != ~(unsigned long)0 ) 
    {  
      AliHLTPHOSRcuFFTDataStruct* fftDataPtr = (AliHLTPHOSRcuFFTDataStruct*)homeReaderPtr->GetBlockData( blk );
      FillHistograms(*fftDataPtr, fftDataPtr->fDataLength);

      blk = homeReaderPtr->FindBlockNdx("TREIRUOF","SOHP", blk+1);
    }
}



void 
AliHLTPHOSOnlineDisplayFourierTab::FillHistograms(const AliHLTPHOSRcuFFTDataStruct psd, const int size)
{
  // See header file for documentation
  //  gStyle->SetOptLogy();
  //  gStyle->SetOptStat(false);

  char tmpname[256];
  char tmptitle[256];

  int linewidth = 0;
  // double  linewidth = 1.2;

  for(int gain = 0; gain < NGAINS; gain ++ )
    {
      if( fFourierHistoNew[gain] == 0)
	{
	  sprintf(tmptitle, "PSD averaged over all %s channels: Most recent event", Gain2Text(gain, ' ')); 
	  sprintf(tmpname,  "PSD_averaged_over_all_%s_channels__most_recent_event", Gain2Text(gain, '_'));  
	  fFourierHistoNew[gain] = new TH1D(tmpname, tmptitle,  (size/2) +1, 0, SAMPLINGFREQUENCY/2);
	  fFourierHistoNew[gain]->GetXaxis()->SetTitle("f/MHz");
	  fFourierHistoNew[gain]->GetYaxis()->SetTitle("Power (arbitrary units)"); 
	  fFourierHistoNew[gain]->SetLineWidth(linewidth);

	}
      if (fFourierHistoOld[gain] == 0)
	{
	  sprintf(tmptitle, "PSD averaged over all %s channels: Previous event", Gain2Text(gain, ' ')); 
	  sprintf(tmpname,  "PSD_averaged_over_all_%s_channels__previous_event", Gain2Text(gain, '_')); 
	  fFourierHistoOld[gain] = new TH1D(tmpname, tmptitle,  (size/2) +1, 0, SAMPLINGFREQUENCY/2);
	  fFourierHistoOld[gain]->GetXaxis()->SetTitle("f/MHz");
	  fFourierHistoOld[gain]->GetYaxis()->SetTitle("Power (arbitrary units)"); 
	  fFourierHistoOld[gain]->SetLineWidth(linewidth);

	}
      if( fFourierHistoAccumulated[gain] == 0 )
	{
	  sprintf(tmptitle, "PSD averaged over all %s channels: All events", Gain2Text(gain, ' ')); 
	  sprintf(tmpname,  "PSD_averaged_over_all_%s_channels__All_events", Gain2Text(gain, '_')); 
	  fFourierHistoAccumulated[gain] = new TH1D(tmpname, tmptitle,  (size/2) +1, 0, SAMPLINGFREQUENCY/2);
	  fFourierHistoAccumulated[gain]->GetXaxis()->SetTitle("f/MHz");
	  fFourierHistoAccumulated[gain]->GetYaxis()->SetTitle("Power (arbitrary units)"); 
	  fFourierHistoAccumulated[gain]->SetLineWidth(linewidth);

	}

//       fFourierHistoNew[gain]->Reset();
//       fFourierHistoOld[gain]->Reset();
//       fFourierHistoAccumulated[gain]->Reset();
      
      for(int i = 0; i <size/2; i++)
	{
	  fFourierHistoOld[gain]->SetBinContent(i+1,  fFourierHistoNew[gain]->GetBinContent(i+1));
	  fFourierHistoNew[gain]->SetBinContent(i+1,  psd.fGlobalLastPSD[gain][i] );
	  fFourierHistoAccumulated[gain]->SetBinContent(i+1,  psd.fGlobalAccumulatedPSD[gain][i] );
	}

    }
}




void
AliHLTPHOSOnlineDisplayFourierTab::InitDisplay(TGTab  *tabPtr)
{
  // See header file for documentatino
  for(int gain=0; gain < NGAINS; gain++)
    {
      char gainLabel[100];
      char label[256];
 
      //     Gain2Text
      fOnlineDisplayPtr->Gain2Text(gain,gainLabel);
      sprintf(label, "PHOS Fourier transform %s", gainLabel);
      fgLegoPlotPtr[gain] = new AliHLTPHOSOnlineDisplayTH2D(fOnlineDisplayPtr, label, label, 
      							    NXCOLUMNSMOD*NMODULES , 0, NXCOLUMNSMOD*NMODULES,  
      							    NZROWSMOD,   0, NZROWSMOD);   
	   //    fgLegoPlotPtr[gain]->SetGain(HIGHGAIN);
      fgLegoPlotPtr[gain]->SetMaximum(1023);
      fgLegoPlotPtr[gain]->Reset();
      fgLegoPlotPtr[gain]->GetXaxis()->SetRange(XRANGESTART, XRANGEEND);
 
   }
  


  TGLayoutHints *fL1 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 15, 1);

  TGCompositeFrame *tf = tabPtr->AddTab("Power spectrum");
  fSubTab1 = new TGTab(tf, 100, 100);
  TGCompositeFrame *tf2 = fSubTab1->AddTab("Most recent event");  
  fSubF1 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
  fEc1 = new TRootEmbeddedCanvas("ecf1", fSubF1, 100, 100);
  fSubF1->AddFrame(fEc1, fL1);
  fEc2 = new TRootEmbeddedCanvas("ecf2", fSubF1, 100, 100);
  fSubF1->AddFrame(fEc2, fL1);
  tf2->AddFrame(fSubF1, fL1);
  
  tf2 = fSubTab1->AddTab("Previous event"); 
  fSubF2 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
  tf2->AddFrame(fSubF2, fL1);
  fEc3 = new TRootEmbeddedCanvas("ecf3", fSubF2, 100, 100);
  fSubF2->AddFrame(fEc3, fL1);
  fEc4 = new TRootEmbeddedCanvas("ecf4", fSubF2, 100, 100);
  fSubF2->AddFrame(fEc4, fL1);
 
  
  tf2 = fSubTab1->AddTab("Accumulated"); 
  fSubF3 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
  tf2->AddFrame(fSubF3, fL1);
  fEc5 = new TRootEmbeddedCanvas("ecf5", fSubF3, 100, 100);
  fSubF3->AddFrame(fEc5, fL1);
  fEc6 = new TRootEmbeddedCanvas("ecf6", fSubF3, 100, 100);
  fSubF3->AddFrame(fEc6, fL1);
  fSubTab1->Resize();
  tf->AddFrame(fSubTab1, fL1);
  

  fgEventButtPtr = new  AliHLTPHOSGetEventButton(fSubF1, "get fourier", 'e');
}



void
AliHLTPHOSOnlineDisplayFourierTab::UpdateDisplay()
{
  // See header file for documentation
  if( fFourierHistoNew[HIGHGAIN])
    {
      fgCanvasPtr[HIGHGAIN] =  fEc1->GetCanvas();
      fgCanvasPtr[HIGHGAIN]->cd();
      gPad->SetLogy();
      //  fgLegoPlotPtr[HIGHGAIN]->Draw("LGZ");
      fFourierHistoNew[HIGHGAIN]->Draw();
      fgCanvasPtr[HIGHGAIN]->Update();
    }

  if( fFourierHistoNew[LOWGAIN])
    {
      fgCanvasPtr[LOWGAIN] = fEc2->GetCanvas();
      fgCanvasPtr[LOWGAIN]->cd();
      gPad->SetLogy();
      //  fgLegoPlotPtr[LOWGAIN]->Draw("HGZ");
      fFourierHistoNew[LOWGAIN]->Draw();
      fgCanvasPtr[LOWGAIN]->Update();
    }
  if( fFourierHistoOld[HIGHGAIN])
    {
      fgCanvasPtr[HIGHGAIN] =  fEc3->GetCanvas();
      fgCanvasPtr[HIGHGAIN]->cd();
      gPad->SetLogy();
      // fgLegoPlotPtr[HIGHGAIN]->Draw("Low gain");
      fFourierHistoOld[HIGHGAIN]->Draw();
      fgCanvasPtr[HIGHGAIN]->Update();
    }

  if( fFourierHistoOld[LOWGAIN])
    {
      fgCanvasPtr[LOWGAIN] = fEc4->GetCanvas();
      fgCanvasPtr[LOWGAIN]->cd();
      //fgLegoPlotPtr[LOWGAIN]->Draw("High gain");
      gPad->SetLogy();
      fFourierHistoOld[LOWGAIN]->Draw();
      fgCanvasPtr[LOWGAIN]->Update();
    }

  if( fFourierHistoAccumulated[HIGHGAIN])
    {
      fgCanvasPtr[HIGHGAIN] =  fEc5->GetCanvas();
      fgCanvasPtr[HIGHGAIN]->cd();
      gPad->SetLogy();
      fFourierHistoAccumulated[HIGHGAIN]->Draw();
      //  fgLegoPlotPtr[HIGHGAIN]->Draw("CONTZ");
      fgCanvasPtr[HIGHGAIN]->Update();
    }

  if( fFourierHistoAccumulated[LOWGAIN])
    {
      fgCanvasPtr[LOWGAIN] = fEc6->GetCanvas();
      fgCanvasPtr[LOWGAIN]->cd();
      gPad->SetLogy();
      //  fgLegoPlotPtr[LOWGAIN]->Draw("CONTZ");
      fFourierHistoAccumulated[LOWGAIN]->Draw();
      fgCanvasPtr[LOWGAIN]->Update();
    }
  
}


const  char* 
AliHLTPHOSOnlineDisplayFourierTab::Gain2Text(const int gain, const char delimeter)
{
  // See header file for documentation
  if(gain ==  LOWGAIN)
    {
      sprintf(fGainText, "low%cgain", delimeter);

    }
  else if(gain ==  HIGHGAIN)
    {
      sprintf(fGainText, "high%cgain", delimeter);
    }
  else
    {
      sprintf(fGainText, "Error, invalid gain");
    }
  return fGainText;
}
