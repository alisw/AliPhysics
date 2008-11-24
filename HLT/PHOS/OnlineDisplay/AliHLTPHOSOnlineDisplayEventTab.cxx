// $Id$

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


AliHLTPHOSOnlineDisplayEventTab::AliHLTPHOSOnlineDisplayEventTab(AliHLTPHOSOnlineDisplay *onlineDisplayPtr, TGTab  *tabPtr, 
								 AliHLTHOMERReader *homerSyncPtr, AliHLTHOMERReader *homerPtrs[MAX_HOSTS], 
								 int nHosts,  const int runnumber) :  AliHLTPHOSOnlineDisplayTab()
{

  /*
  if(fIsSetRunNumber == true)
    {
      for(int i=0; i < N_GAINS; i++)
	{
	  fgLegoPlotPtr[gain]
	}
   }
  */


  fShmPtr = new AliHLTPHOSSharedMemoryInterface();
  fOnlineDisplayPtr =  onlineDisplayPtr;


  for(int gain = 0; gain < N_GAINS; gain ++ )
    {
      fgCanvasPtr[gain] = 0;
      fgLegoPlotPtr[gain] = 0;
 
      for(int mod =0; mod <N_MODULES; mod ++)
	{
	  for(int rcu_x_coord = 0; rcu_x_coord < N_ZRCU_COORD; rcu_x_coord ++)
	    {
	      for(int rcu_z_coord = 0; rcu_z_coord < N_XRCU_COORD; rcu_z_coord ++) 
		{
		  for(int z = 0; z < N_ZROWS_RCU; z ++)
		    {
		      for(int x = 0; x < N_XCOLUMNS_RCU; x ++)
			{
			  fChannelData[mod][rcu_z_coord][rcu_x_coord][x][z][gain] = new int[ALTRO_MAX_SAMPLES];
			  fNChannelSamples[mod][rcu_z_coord][rcu_x_coord][x][z][gain] = 0;
			  fChannelEnergy[mod][rcu_z_coord][rcu_x_coord][x][z][gain] = 0;
			}
		    }	   
		}
	    }
	}
    }

  for(int i=0; i<MAX_HOSTS; i++)
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
AliHLTPHOSOnlineDisplayEventTab::FindFourierBlocks(AliHLTHOMERReader *homerReaderPtr)
{
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
	  fgLegoPlotPtr[tmpGain]->Fill(moduleID*N_XCOLUMNS_MOD + tmpX +  N_XCOLUMNS_RCU*cellEnergiesPtr->fRcuX,  
				    tmpZ + N_ZROWS_RCU*cellEnergiesPtr->fRcuZ, currentChannel->fEnergy);
	      
	  // CRAP PTH
	  if(tmpGain == HIGH_GAIN)
	    {
	    //   gAliEveBoxSet->AddBox(2.2*(tmpX + N_XCOLUMNS_RCU*cellEnergiesPtr->fRcuX) - 1.1,
// 				    0,
// 				    2.2*(tmpZ + N_ZROWS_RCU*cellEnergiesPtr->fRcuZ) - 1.1,
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
	      
	      
	      if(nSamples > ALTRO_MAX_SAMPLES || nSamples < 0 )
		{
		  cout << __FILE__<< ":" <<__LINE__ <<"ERROR, nsamples = "<< nSamples <<" exeeds allowd range, max number of samples is  "<< ALTRO_MAX_SAMPLES << endl;
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
  if(fgAccumulate == kFALSE)
    {
      for(int gain=0; gain < N_GAINS; gain++)
	{
	  if(fgLegoPlotPtr[gain] !=0)
	    {
	      fgLegoPlotPtr[gain]->Reset(); 
	    }
	}
    } 
}


void
AliHLTPHOSOnlineDisplayEventTab::InitDisplay(TGTab  *tabPtr, const int runnumber)
{
  //  gStyle->SetOptLogy();
  ///  gStyle->SetOptStat(true);

  for(int gain=0; gain < N_GAINS; gain++)
    {
      char gainLabel[100];
      char label[256];
 
      //     Gain2Text
      fOnlineDisplayPtr->Gain2Text(gain,gainLabel);
      sprintf(label, "PHOS HLT Online Display %s", gainLabel);
      fgLegoPlotPtr[gain] = new AliHLTPHOSOnlineDisplayTH2D(fOnlineDisplayPtr, label, label, 
							    N_XCOLUMNS_MOD*N_MODULES , 0, N_XCOLUMNS_MOD*N_MODULES,  
							    N_ZROWS_MOD,   0, N_ZROWS_MOD);   
      
      //      cout << __FILE__ << ":" << __LINE__ << " Runnumber = "  << runnumber <<endl;
   
      fgLegoPlotPtr[gain]->SetRunNumber(runnumber);
      fgLegoPlotPtr[gain]->SetMaximum(1023);
      fgLegoPlotPtr[gain]->Reset();
      fgLegoPlotPtr[gain]->GetXaxis()->SetRange(X_RANGE_START, X_RANGE_END);
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

  fgCanvasPtr[HIGH_GAIN] =  fEc1->GetCanvas();
  fgCanvasPtr[HIGH_GAIN]->cd();
  fgLegoPlotPtr[HIGH_GAIN]->Draw("LEGO2Z");
  fgCanvasPtr[HIGH_GAIN]->Update();
  fgCanvasPtr[LOW_GAIN] = fEc2->GetCanvas();
  fgCanvasPtr[LOW_GAIN]->cd();
  fgLegoPlotPtr[LOW_GAIN]->Draw("LEGO2Z");
  fgCanvasPtr[LOW_GAIN]->Update();

  fgCanvasPtr[HIGH_GAIN] =  fEc3->GetCanvas();
  fgCanvasPtr[HIGH_GAIN]->cd();
  fgLegoPlotPtr[HIGH_GAIN]->Draw("SCAT");
  fgCanvasPtr[HIGH_GAIN]->Update();
  fgCanvasPtr[LOW_GAIN] = fEc4->GetCanvas();
  fgCanvasPtr[LOW_GAIN]->cd();
  fgLegoPlotPtr[LOW_GAIN]->Draw("SCAT");
  fgCanvasPtr[LOW_GAIN]->Update();

  /* 
 fgCanvasPtr[HIGH_GAIN] =  fEc5->GetCanvas();
  fgCanvasPtr[HIGH_GAIN]->cd();
  fgLegoPlotPtr[HIGH_GAIN]->Draw("CONTZ");
  fgCanvasPtr[HIGH_GAIN]->Update();
  fgCanvasPtr[LOW_GAIN] = fEc6->GetCanvas();
  fgCanvasPtr[LOW_GAIN]->cd();
  fgLegoPlotPtr[LOW_GAIN]->Draw("CONTZ");
  fgCanvasPtr[LOW_GAIN]->Update();
  */

  fgCanvasPtr[HIGH_GAIN] =  fEc5->GetCanvas();
  fgCanvasPtr[HIGH_GAIN]->cd();
  fgLegoPlotPtr[HIGH_GAIN]->Draw("COLZ");
  fgCanvasPtr[HIGH_GAIN]->Update();
  fgCanvasPtr[LOW_GAIN] = fEc6->GetCanvas();
  fgCanvasPtr[LOW_GAIN]->cd();
  fgLegoPlotPtr[LOW_GAIN]->Draw("COLZ");
  fgCanvasPtr[LOW_GAIN]->Update();


}
