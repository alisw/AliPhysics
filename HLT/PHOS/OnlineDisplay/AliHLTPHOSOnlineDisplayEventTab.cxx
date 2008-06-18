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

using namespace std;

// MT Crap
#include <TMath.h>
#include <TEveManager.h>
#include <TEveBoxSet.h>

TEveBoxSet* gAliEveBoxSet = 0;

AliHLTPHOSOnlineDisplayEventTab::AliHLTPHOSOnlineDisplayEventTab()
{
  cout << "ERROR: You cannot create a onlinedisplay Tab without arguments" << endl;
}


AliHLTPHOSOnlineDisplayEventTab::AliHLTPHOSOnlineDisplayEventTab(AliHLTPHOSOnlineDisplay *onlineDisplayPtr, TGTab  *tabPtr, 
								 AliHLTHOMERReader *homerSyncPtr, AliHLTHOMERReader *homerPtrs[MAX_HOSTS], int nHosts) :  AliHLTPHOSOnlineDisplayTab()
{
  fShmPtr = new AliHLTPHOSSharedMemoryInterface();
  fOnlineDisplayPtr =  onlineDisplayPtr;

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
		      for(int gain = 0; gain < N_GAINS; gain ++ )
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

  fgCanvasHGPtr = 0;
  fgCanvasLGPtr = 0;
  fgLegoPlotLGPtr = 0;
  fgLegoPlotHGPtr = 0;

  fgNHosts = nHosts;
  InitDisplay(tabPtr);
}


AliHLTPHOSOnlineDisplayEventTab::~AliHLTPHOSOnlineDisplayEventTab()
{

}


Int_t
AliHLTPHOSOnlineDisplayEventTab::GetRawData(TH1D *histPtr, int mod, int rcuX, int rcuZ, int x, int z, int gain)
{
  for(int i=0;  i < fNChannelSamples[mod][rcuX][rcuZ][x][z][gain] ; i++)
    {
       histPtr->SetBinContent(i, fChannelData[mod][rcuX][rcuZ][x][z][gain][i]);
    }
  return fNChannelSamples[mod][rcuX][rcuZ][x][z][gain];
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
  if (gAliEveBoxSet == 0)
  {
    is_first = true;
    gAliEveBoxSet = new TEveBoxSet("PHOS module");
    // gAliEveBoxSet->SetSecSelectCommand("Draw()");
    // gAliEveBoxSet->SetSecSelectCommand("phos_histo_draw"); 
    gEve->AddElement(gAliEveBoxSet);
  }
  gAliEveBoxSet->Reset(TEveBoxSet::kBT_AABox, kFALSE, 128);

  DoGetNextEvent();
  UpdateDisplay();

  gAliEveBoxSet->ElementChanged();
  gEve->Redraw3D(is_first);

  fgEvntCnt ++;
}



void 
AliHLTPHOSOnlineDisplayEventTab::ReadBlockData(AliHLTHOMERReader *homeReaderPtr)
{  
  AliHLTPHOSValidCellDataStruct *currentChannel =0;
  cout << "AliHLTPHOSOnlineDisplayEventTab::ReadBlockDat, Reading block data" << endl;

  unsigned long blk = homeReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF );

  while ( blk != ~(unsigned long)0 ) 
    {
      Int_t moduleID;
      Int_t rcuX = 0;
      Int_t rcuZ = 0;
      AliHLTPHOSRcuCellEnergyDataStruct* cellEnergiesPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)homeReaderPtr->GetBlockData( blk ); 
      
      unsigned int *t = (unsigned int*)cellEnergiesPtr;
      
      for(int i = 0; i < 10000; i++)
	{
	  //	  printf("%f\t", (float)*t);
	  if(i%30 == 0)
	    {
	      //  printf("\ni = %d", i);
	    }
	  t ++;
	}

      for(int gain = 1; gain < N_GAINS; gain ++)
	{
	  for(int x=0; x <N_XCOLUMNS_RCU; x ++ )
	    {
	      //	      printf("\nnewline");
	      for(int z=0; z <N_ZROWS_RCU; z ++ ) 
		{
		  //		  printf("%f\t",cellEnergiesPtr->fValidData[x][z][gain].fEnergy);
	    
		}
	    }
	}

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
	  
	  //	  cout << "AliHLTPHOSOnlineDisplayEventTab::ReadBlockData gain = " << tmpGain << endl;
	  
	  if(cnt < 20)
	    {
	      //cout << "the addresss of fData is " << (void *)currentChannel->fData  << endl;
	    }
	  // We have to unroll the gain for loop because of the hack by MT to display the raw data in AliEve
	  // This is crap !! but it is the dirty one hour solution (and it works) 
	  

	  if(
	     tmpGain == HIGH_GAIN)
	    {
	      fgLegoPlotHGPtr->Fill(moduleID*N_XCOLUMNS_MOD + tmpX +  N_XCOLUMNS_RCU*cellEnergiesPtr->fRcuX,  
	      		    tmpZ + N_ZROWS_RCU*cellEnergiesPtr->fRcuZ, currentChannel->fEnergy);
	      
	      //  cout <<   << endl;

	      // CRAP PTH
	      
	      gAliEveBoxSet->AddBox(2.2*(tmpX + N_XCOLUMNS_RCU*cellEnergiesPtr->fRcuX) - 1.1,
				    0,
				    2.2*(tmpZ + N_ZROWS_RCU*cellEnergiesPtr->fRcuZ) - 1.1,
				    2.2,
				    0.4*140*currentChannel->fEnergy/1024,
				    2.2);
	      gAliEveBoxSet->DigitValue(TMath::Nint(currentChannel->fEnergy));
	      
	      if(cellEnergiesPtr->fHasRawData == true)
		{
		  Int_t nSamples = 0;
		  Int_t* rawPtr = 0;

		  rawPtr = fShmPtr->GetRawData(nSamples);

		  fNChannelSamples[moduleID][rcuX][rcuZ][tmpX][tmpZ][tmpGain] = nSamples;
		  //fChannelEnergy[moduleID][rcuX][rcuZ][tmpX][tmpZ][tmpGain] = currentChannel->fEnergy;

		  for(int j= 0; j< nSamples; j++)
		    {
		     //  if(j == fNTotalSamples) 
// 			{
// 			  break;
// 			}		  
 		      
		      fChannelData[moduleID][cellEnergiesPtr->fRcuX][cellEnergiesPtr->fRcuZ][tmpX][tmpZ][HIGH_GAIN][j] = rawPtr[j];  
		    }
		}
	    }
	  
	  else if(tmpGain == LOW_GAIN)
	    {
	      /*	 
	      gAliEveBoxSet->AddBox(2.2*(tmpX + N_XCOLUMNS_RCU*cellEnergiesPtr->fRcuX) - 1.1,
				    0,
				    2.2*(tmpZ + N_ZROWS_RCU*cellEnergiesPtr->fRcuZ) - 1.1,
				    2.2,
				    0.4*140*currentChannel->fEnergy/1024,
				    2.2);
	      gAliEveBoxSet->DigitValue(TMath::Nint(currentChannel->fEnergy));
	      */

	      fgLegoPlotLGPtr->Fill(moduleID*N_XCOLUMNS_MOD + tmpX +  N_XCOLUMNS_RCU*cellEnergiesPtr->fRcuX,
				    tmpZ + N_ZROWS_RCU*cellEnergiesPtr->fRcuZ,    currentChannel->fEnergy);
	      if(cellEnergiesPtr->fHasRawData == true)
		{
		  Int_t nSamples = 0;
		  Int_t* rawPtr = 0;
		  rawPtr = fShmPtr->GetRawData(nSamples);

		  fNChannelSamples[moduleID][rcuX][rcuZ][tmpX][tmpZ][tmpGain] = nSamples;
		  //fChannelEnergy[moduleID][rcuX][rcuZ][tmpX][tmpZ][tmpGain] = currentChannel->fEnergy;


		  for(int j= 0; j< nSamples; j++)
		    {
// 		      if(j == fNTotalSamples) 
// 			{
// 			  break;
// 			}		  
		      
		      fChannelData[moduleID][cellEnergiesPtr->fRcuX][cellEnergiesPtr->fRcuZ][tmpX][tmpZ][LOW_GAIN][j] = rawPtr[j];   
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
      if(fgLegoPlotHGPtr !=0)
	{
	  fgLegoPlotHGPtr->Reset(); 
	}

      if(fgLegoPlotLGPtr !=0)
	{
	  fgLegoPlotLGPtr->Reset();
	}  
    }
 }


void
AliHLTPHOSOnlineDisplayEventTab::InitDisplay(TGTab  *tabPtr)
{
  fgLegoPlotHGPtr = new AliHLTPHOSOnlineDisplayTH2D(fOnlineDisplayPtr, "Cosmics, High gain", "PHOS HLT: Cosmics", 
						    N_XCOLUMNS_MOD*N_MODULES , 0, N_XCOLUMNS_MOD*N_MODULES,  
						    N_ZROWS_MOD,   0, N_ZROWS_MOD);    
  fgLegoPlotHGPtr->SetGain(HIGH_GAIN);
  fgLegoPlotHGPtr->SetMaximum(1023);
  fgLegoPlotHGPtr->Reset();
  fgLegoPlotHGPtr->GetXaxis()->SetRange(X_RANGE_START, X_RANGE_END);
  fgLegoPlotLGPtr = new AliHLTPHOSOnlineDisplayTH2D(fOnlineDisplayPtr, "Cosmics, Low gain", "PHOS HLT: Cosmics",  
						    N_XCOLUMNS_MOD* N_MODULES , 0, N_XCOLUMNS_MOD* N_MODULES,  
						    N_ZROWS_MOD,          0, N_ZROWS_MOD);
  fgLegoPlotLGPtr->SetGain(LOW_GAIN);
  fgLegoPlotLGPtr->SetMaximum(1023); 
  fgLegoPlotLGPtr->Reset();
  fgLegoPlotLGPtr->GetXaxis()->SetRange(X_RANGE_START, X_RANGE_END);

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
  fgCanvasHGPtr =  fEc1->GetCanvas();
  fgCanvasHGPtr->cd();
  fgLegoPlotHGPtr->Draw("LEGO2Z");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr = fEc2->GetCanvas();
  fgCanvasLGPtr->cd();
  fgLegoPlotLGPtr->Draw("LEGO2Z");
  fgCanvasLGPtr->Update();

  fgCanvasHGPtr =  fEc3->GetCanvas();
  fgCanvasHGPtr->cd();
  fgLegoPlotHGPtr->Draw("SCAT");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr = fEc4->GetCanvas();
  fgCanvasLGPtr->cd();
  fgLegoPlotLGPtr->Draw("SCAT");
  fgCanvasLGPtr->Update();

  fgCanvasHGPtr =  fEc5->GetCanvas();
  fgCanvasHGPtr->cd();
  fgLegoPlotHGPtr->Draw("CONTZ");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr = fEc6->GetCanvas();
  fgCanvasLGPtr->cd();
  fgLegoPlotLGPtr->Draw("CONTZ");
  fgCanvasLGPtr->Update();
}
