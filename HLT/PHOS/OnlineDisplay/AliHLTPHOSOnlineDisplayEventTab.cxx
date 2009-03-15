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

#include "AliHLTPHOSOnlineDisplayEventTab.h"
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

//#include "TStyle.h"

using namespace std;

// MT Crap
#include <TMath.h>
//#include <TEveManager.h>
//#include <TEveBoxSet.h>

//TEveBoxSet* gAliEveBoxSet = 0;

AliHLTPHOSOnlineDisplayEventTab::AliHLTPHOSOnlineDisplayEventTab()
{
  cout << "ERROR: You cannot create a onlinedisplay Tab without arguments" << endl;
}


AliHLTPHOSOnlineDisplayEventTab::AliHLTPHOSOnlineDisplayEventTab(AliHLTPHOSOnlineDisplay * const onlineDisplayPtr, TGTab  *tabPtr, 
								 const AliHLTHOMERReader * const homerSyncPtr, AliHLTHOMERReader * const homerPtrs[MAXHOSTS], 
								 int nHosts,  const int runnumber) :  AliHLTPHOSOnlineDisplayTab()
{
  //comment
  /*
  if(fIsSetRunNumber == true)
    {
      for(int i=0; i < NGAINS; i++)
	{
	  fgLegoPlotPtr[gain]
	}
   }
  */


  fShmPtr = new AliHLTPHOSSharedMemoryInterface();
  fOnlineDisplayPtr =  onlineDisplayPtr;


  for(int gain = 0; gain < NGAINS; gain ++ )
    {
      fgCanvasPtr[gain] = 0;
      fgLegoPlotPtr[gain] = 0;
 
      for(int mod =0; mod <NMODULES; mod ++)
	{
	  for(int rcuxcoord = 0; rcuxcoord < NZRCUCOORD; rcuxcoord ++)
	    {
	      for(int rcuzcoord = 0; rcuzcoord < NXRCUCOORD; rcuzcoord ++) 
		{
		  for(int z = 0; z < NZROWSRCU; z ++)
		    {
		      for(int x = 0; x < NXCOLUMNSRCU; x ++)
			{
			  fChannelData[mod][rcuzcoord][rcuxcoord][x][z][gain] = new int[ALTROMAXSAMPLES];
			  fNChannelSamples[mod][rcuzcoord][rcuxcoord][x][z][gain] = 0;
			  fChannelEnergy[mod][rcuzcoord][rcuxcoord][x][z][gain] = 0;
			}
		    }	   
		}
	    }
	}
    }

  for(int i=0; i<MAXHOSTS; i++)
    {
       fgHomerReadersPtr[i] = 0;
    }

  fgHomerReaderPtr = homerSyncPtr;
  
  for(int i=0; i<nHosts; i++)
    {
      fgHomerReadersPtr[i] = homerPtrs[i] ;

    }

  fgNHosts = nHosts;
  InitDisplay(tabPtr, runnumber);
}


AliHLTPHOSOnlineDisplayEventTab::~AliHLTPHOSOnlineDisplayEventTab()
{
  //comment
}


Int_t
AliHLTPHOSOnlineDisplayEventTab::GetRawData(TH1D *histPtr, int x, int z, int gain)
{
  int tmpModID = x/64;
  int tmpRcuZ = z/32;
  int tmpRcuX = (x%64)/32;
  int tmpZ = z%28;
  int tmpX = x%32;

  for(int i=0;  i < fNChannelSamples[tmpModID][tmpRcuX][tmpRcuZ][tmpX][tmpZ][gain] ; i++)
    {
      histPtr->SetBinContent(i, fChannelData[tmpModID][tmpRcuX][tmpRcuZ][tmpX][tmpZ][gain][i]);  
    }
  return fNChannelSamples[tmpModID][tmpRcuX][tmpRcuZ][tmpX][tmpZ][gain];
}


int
AliHLTPHOSOnlineDisplayEventTab::GetNextEvent()
{
  //comment
  ResetDisplay();
  // MT crap
  Bool_t is_first = false;
 //  if (gAliEveBoxSet == 0)
//   {
//     is_first = true;
//     gAliEveBoxSet = new TEveBoxSet("PHOS module");
//     // gAliEveBoxSet->SetSecSelectCommand("Draw()");
//     // gAliEveBoxSet->SetSecSelectCommand("phos_histo_draw"); 
//     gEve->AddElement(gAliEveBoxSet);
//   }
//   gAliEveBoxSet->Reset(TEveBoxSet::kBT_AABox, kFALSE, 128);

  DoGetNextEvent();
  UpdateDisplay();

  //  gAliEveBoxSet->ElementChanged();
  // gEve->Redraw3D(is_first);

  fgEvntCnt ++;
}


void
AliHLTPHOSOnlineDisplayEventTab::FindFourierBlocks(AliHLTHOMERReader * const homerReaderPtr) const
{
  //comment
 cout << "AliHLTPHOSOnlineDisplayEventTab::FindFourierBlocks" << endl; 
  // unsigned long blk = homeReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF );
  unsigned long blk = homerReaderPtr->FindBlockNdx(" TREIRUOF","SOHP", 0xFFFFFFFF );

  while ( blk != ~(unsigned long)0 )
    {
      cout << "AliHLTPHOSOnlineDisplayEventTab::FindFourierBlocks(homerReaderPtr) FOUND FOURIER DATA !!!!!!!!!!!!!!" << endl;


      blk = homerReaderPtr->FindBlockNdx("TREIRUOF","SOHP", 0xFFFFFFFF );
    }

}


void 
AliHLTPHOSOnlineDisplayEventTab::ReadBlockData(AliHLTHOMERReader *homeReaderPtr)
{  
  //comment
  AliHLTPHOSValidCellDataStruct *currentChannel =0;
  cout << "AliHLTPHOSOnlineDisplayEventTab::ReadBlockDat, Reading block data, therere are " <<  homeReaderPtr->GetBlockCnt() << " blocks " <<endl;

  FindFourierBlocks(homeReaderPtr);

  unsigned long blk = homeReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF );

  int cnt = 0;

  //CRAP PT
  //  FindFourierBlocks(homeReaderPtr);

  while ( blk != ~(unsigned long)0 ) 
    {
      Int_t moduleID;
      Int_t rcuX = 0;
      Int_t rcuZ = 0;
      AliHLTPHOSRcuCellEnergyDataStruct* cellEnergiesPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)homeReaderPtr->GetBlockData( blk ); 
      
      unsigned int *t = (unsigned int*)cellEnergiesPtr;
      
      moduleID = cellEnergiesPtr->fModuleID ;
      rcuX = cellEnergiesPtr->fRcuX;
      rcuZ = cellEnergiesPtr->fRcuZ;

      cout << "AliHLTPHOSOnlineDisplayEventTab::ReadBlockData,  fModuleID =" <<moduleID << endl; 

      Int_t tmpZ;
      Int_t tmpX;
      Int_t tmpGain;
      int cnt = 0;
      Int_t* tmpPtr = 0;

      fShmPtr->SetMemory(cellEnergiesPtr);
      currentChannel = fShmPtr->NextChannel();

      while(currentChannel != 0)
	{
	  cnt ++;
	  tmpZ = currentChannel->fZ;
	  tmpX = currentChannel->fX;
	  tmpGain =  currentChannel->fGain;
	  fgLegoPlotPtr[tmpGain]->Fill(moduleID*NXCOLUMNSMOD + tmpX +  NXCOLUMNSRCU*cellEnergiesPtr->fRcuX,  
				    tmpZ + NZROWSRCU*cellEnergiesPtr->fRcuZ, currentChannel->fEnergy);
	      
	  // CRAP PTH
	  if(tmpGain == HIGHGAIN)
	    {
	    //   gAliEveBoxSet->AddBox(2.2*(tmpX + N_XCOLUMNS_RCU*cellEnergiesPtr->fRcuX) - 1.1,
// 				    0,
// 				    2.2*(tmpZ + N_ZROWSRCU*cellEnergiesPtr->fRcuZ) - 1.1,
// 				    2.2,
// 				    0.4*140*currentChannel->fEnergy/1024,
// 				    2.2);
// 	      gAliEveBoxSet->DigitValue(TMath::Nint(currentChannel->fEnergy));
	    }
 
	  if(cellEnergiesPtr->fHasRawData == true)
	    //if(0)
	    {
	      Int_t nSamples = 0;
	      Int_t* rawPtr = 0;
	      rawPtr = fShmPtr->GetRawData(nSamples);
	      fNChannelSamples[moduleID][rcuX][rcuZ][tmpX][tmpZ][tmpGain] = nSamples;

	      //	      cout << __FILE__ << ":" <<  __LINE__  <<" gain = " << tmpGain << " z = "<< tmpZ << " x = " << tmpX;
	      //	      cout << " nsamples = " << nSamples;
	      //	      cout << __FILE__ << ":" << __LINE__ << " the address of raw ptr = " << rawPtr  << endl;
	      
	      
	      if(nSamples > ALTROMAXSAMPLES || nSamples < 0 )
		{
		  cout << __FILE__<< ":" <<__LINE__ <<"ERROR, nsamples = "<< nSamples <<" exeeds allowd range, max number of samples is  "<< ALTROMAXSAMPLES << endl;
		}
	      else
		{
		  for(int j= 0; j< nSamples; j++)
		    {
		      //		      cout << __FILE__ << ":" << __LINE__ << " nsamples = " << nSamples << "  j =" << j << endl;
		      fChannelData[moduleID][cellEnergiesPtr->fRcuX][cellEnergiesPtr->fRcuZ][tmpX][tmpZ][tmpGain][j] = rawPtr[j];  
		    }
	      
		}
	    }
	
	  currentChannel = fShmPtr->NextChannel();
	}
      
      blk = homeReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF, blk+1);

    }
}




void
AliHLTPHOSOnlineDisplayEventTab::ResetDisplay()
{
  //comment
  if(fgAccumulate == kFALSE)
    {
      for(int gain=0; gain < NGAINS; gain++)
	{
	  if(fgLegoPlotPtr[gain] !=0)
	    {
	      fgLegoPlotPtr[gain]->Reset(); 
	    }
	}
    } 
}


void
AliHLTPHOSOnlineDisplayEventTab::InitDisplay(const TGTab  * const tabPtr, const int runnumber)
{
  //  gStyle->SetOptLogy();
  ///  gStyle->SetOptStat(true);

  for(int gain=0; gain < NGAINS; gain++)
    {
      char gainLabel[100];
      char label[256];
 
      //     Gain2Text
      fOnlineDisplayPtr->Gain2Text(gain,gainLabel);
      sprintf(label, "PHOS HLT Online Display %s", gainLabel);
      fgLegoPlotPtr[gain] = new AliHLTPHOSOnlineDisplayTH2D(fOnlineDisplayPtr, label, label, 
							    NXCOLUMNSMOD*NMODULES , 0, NXCOLUMNSMOD*NMODULES,  
							    NZROWSMOD,   0, NZROWSMOD);   
      
      //      cout << __FILE__ << ":" << __LINE__ << " Runnumber = "  << runnumber <<endl;
   
      fgLegoPlotPtr[gain]->SetRunNumber(runnumber);
      fgLegoPlotPtr[gain]->SetMaximum(1023);
      fgLegoPlotPtr[gain]->Reset();
      fgLegoPlotPtr[gain]->GetXaxis()->SetRange(XRANGESTART, XRANGEEND);
    }
  

  TGLayoutHints *fL1 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 15, 1);

  TGCompositeFrame *tf = tabPtr->AddTab("Event display");
  fSubTab1 = new TGTab(tf, 100, 100);
  TGCompositeFrame *tf2 = fSubTab1->AddTab("LEGO");  
  fSubF1 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
  fEc1 = new TRootEmbeddedCanvas("ec1", fSubF1, 100, 100);
  fSubF1->AddFrame(fEc1, fL1);
  fEc2 = new TRootEmbeddedCanvas("ec2", fSubF1, 100, 100);
  fSubF1->AddFrame(fEc2, fL1);
  tf2->AddFrame(fSubF1, fL1);
  
  tf2 = fSubTab1->AddTab("SCAT"); 
  fSubF2 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
  tf2->AddFrame(fSubF2, fL1);
  fEc3 = new TRootEmbeddedCanvas("ec3", fSubF2, 100, 100);
  fSubF2->AddFrame(fEc3, fL1);
  fEc4 = new TRootEmbeddedCanvas("ec4", fSubF2, 100, 100);
  fSubF2->AddFrame(fEc4, fL1);
  
  tf2 = fSubTab1->AddTab("SURF"); 
  fSubF3 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
  tf2->AddFrame(fSubF3, fL1);
  fEc5 = new TRootEmbeddedCanvas("ec5", fSubF3, 100, 100);
  fSubF3->AddFrame(fEc5, fL1);
  fEc6 = new TRootEmbeddedCanvas("ec6", fSubF3, 100, 100);
  fSubF3->AddFrame(fEc6, fL1);
  fSubTab1->Resize();
  tf->AddFrame(fSubTab1, fL1);

  fgEventButtPtr = new  AliHLTPHOSGetEventButton(fSubF1, "get event", 'e');
}



void
AliHLTPHOSOnlineDisplayEventTab::UpdateDisplay()
{
  // gStyle->SetOptLogy();
  //  gStyle->SetOptStat(true);

  fgCanvasPtr[HIGHGAIN] =  fEc1->GetCanvas();
  fgCanvasPtr[HIGHGAIN]->cd();
  fgLegoPlotPtr[HIGHGAIN]->Draw("LEGO2Z");
  fgCanvasPtr[HIGHGAIN]->Update();
  fgCanvasPtr[LOWGAIN] = fEc2->GetCanvas();
  fgCanvasPtr[LOWGAIN]->cd();
  fgLegoPlotPtr[LOWGAIN]->Draw("LEGO2Z");
  fgCanvasPtr[LOWGAIN]->Update();

  fgCanvasPtr[HIGHGAIN] =  fEc3->GetCanvas();
  fgCanvasPtr[HIGHGAIN]->cd();
  fgLegoPlotPtr[HIGHGAIN]->Draw("SCAT");
  fgCanvasPtr[HIGHGAIN]->Update();
  fgCanvasPtr[LOWGAIN] = fEc4->GetCanvas();
  fgCanvasPtr[LOWGAIN]->cd();
  fgLegoPlotPtr[LOWGAIN]->Draw("SCAT");
  fgCanvasPtr[LOWGAIN]->Update();

  /* 
 fgCanvasPtr[HIGHGAIN] =  fEc5->GetCanvas();
  fgCanvasPtr[HIGHGAIN]->cd();
  fgLegoPlotPtr[HIGHGAIN]->Draw("CONTZ");
  fgCanvasPtr[HIGHGAIN]->Update();
  fgCanvasPtr[LOWGAIN] = fEc6->GetCanvas();
  fgCanvasPtr[LOWGAIN]->cd();
  fgLegoPlotPtr[LOWGAIN]->Draw("CONTZ");
  fgCanvasPtr[LOWGAIN]->Update();
  */

  fgCanvasPtr[HIGHGAIN] =  fEc5->GetCanvas();
  fgCanvasPtr[HIGHGAIN]->cd();
  fgLegoPlotPtr[HIGHGAIN]->Draw("COLZ");
  fgCanvasPtr[HIGHGAIN]->Update();
  fgCanvasPtr[LOWGAIN] = fEc6->GetCanvas();
  fgCanvasPtr[LOWGAIN]->cd();
  fgLegoPlotPtr[LOWGAIN]->Draw("COLZ");
  fgCanvasPtr[LOWGAIN]->Update();


}
