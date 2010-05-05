// $Id: AliHLTEMCALOnlineDisplayEventTab.cxx 35108 2009-09-30 01:58:37Z phille $

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

#include "AliHLTEMCALOnlineDisplayEventTab.h"
#include <iostream>
#include "TGFrame.h"
#include "AliHLTEMCALGetEventButton.h"
#include "AliHLTDataTypes.h"
#include "AliHLTHOMERData.h"
#include "AliHLTHOMERReader.h"
#include "AliHLTHOMERWriter.h"
#include "TRootEmbeddedCanvas.h"
#include "AliHLTEMCALOnlineDisplay.h"
#include "AliHLTCaloChannelDataStruct.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
//#include "AliHLTCaloSharedMemoryInterfacev2.h"
#include "AliHLTEMCALSharedMemoryInterface.h"
#include "AliHLTCaloCoordinate.h"
//#include "AliHLTCaloChannelRawDataStruct.h"


using namespace std;

// MT Crap
#include <TMath.h>
#include "AliHLTEMCALOnlineDisplayTH2D.h"

#include <TEveManager.h>
#include <TEveBoxSet.h>

TEveBoxSet* gAliEveBoxSet = 0;

//gEve = new TEveManager(300, 300);

AliHLTEMCALOnlineDisplayEventTab::AliHLTEMCALOnlineDisplayEventTab()
{
  cout << "ERROR: You cannot create a onlinedisplay Tab without arguments" << endl;
}


AliHLTEMCALOnlineDisplayEventTab::AliHLTEMCALOnlineDisplayEventTab(AliHLTEMCALOnlineDisplay * onlineDisplayPtr, TGTab  *tabPtr, 
								 AliHLTHOMERReader * homerSyncPtr, AliHLTHOMERReader * homerPtrs[MAXHOSTS], 
								 int nHosts,  int runnumber) :  AliHLTEMCALOnlineDisplayTab()
{
  
  //  gEve = new TEveManager(300, 300, kFALSE);
  

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


  // fShmPtr = new AliHLTEMCALSharedMemoryInterface();

  // fShmPtr = new AliHLTEMCALSharedMemoryInterfacev2();
 
  fShmPtr = new AliHLTEMCALSharedMemoryInterface();

  

  fOnlineDisplayPtr =  onlineDisplayPtr;


  for(int gain = 0; gain < NGAINS; gain ++ )
    {
      fgCanvasPtr[gain] = 0;
      fgLegoPlotPtr[gain] = 0;
 
      
      for(int mod =0; mod <NMODULES; mod ++)
	{
	  for(int z = 0; z < NZROWSMOD ; z ++)
	      {
		for(int x = 0; x < NXCOLUMNSMOD; x ++)
		  {
		    fChannelData[mod][z][x][gain] = new int[ALTROMAXSAMPLES];
		    fNChannelSamples[mod][z][x][gain] = 0;
		    fChannelEnergy[mod][z][x][gain] = 0;
		  }
	      }
	}
    }

  for(int i=0; i<MAXHOSTS; i++)
    {
       fgHomerReadersPtr[i] = 0;
    }

  fgHomerReaderPtr = const_cast<AliHLTHOMERReader*>(homerSyncPtr);
  
  for(int i=0; i<nHosts; i++)
    {
      fgHomerReadersPtr[i] = homerPtrs[i] ;

    }

  fgNHosts = nHosts;
  InitDisplay(tabPtr, runnumber);
}


AliHLTEMCALOnlineDisplayEventTab::~AliHLTEMCALOnlineDisplayEventTab()
{
  //comment
}



Int_t
AliHLTEMCALOnlineDisplayEventTab::GetRawData(TH1D *histPtr, int x, int z, int gain)
{
  
  //  int tmpModID = x/64;
  
  int tmpModID = (x*z)/(NZROWSMOD*NXCOLUMNSMOD);

  // const int NZROWSMOD      =  48;            /**<Number of rows per module*/       
  //  const int NXCOLUMNSMOD   =  24;     
  //int tmpModID = x/64;

  /* 
     int tmpRcuZ = z/32;
     int tmpRcuX = (x%64)/32;
     int tmpZ = z%28;
     int tmpX = x%32;
  */

  cout << __FILE__ << __LINE__ <<": Getting raw data for mod =" << tmpModID << ", z="<< z << ",x=" << x << endl;

  for(  int i=0;  i <  fNChannelSamples[tmpModID][z][x][gain] ; i++)
    {
      histPtr->SetBinContent(i, fChannelData[tmpModID][z][x][gain][i]);  
    }
  return fNChannelSamples [tmpModID][z][x][gain];
}



int
AliHLTEMCALOnlineDisplayEventTab::GetNextEvent()
{
  ResetDisplay();
  DoGetNextEvent();
  UpdateDisplay();
  fgEvntCnt ++;
}


void
AliHLTEMCALOnlineDisplayEventTab::FindFourierBlocks(AliHLTHOMERReader * const homerReaderPtr) const
{
  //comment
 cout << "AliHLTEMCALOnlineDisplayEventTab::FindFourierBlocks" << endl; 
  // unsigned long blk = homeReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF );
  unsigned long blk = homerReaderPtr->FindBlockNdx(" TREIRUOF","SOHP", 0xFFFFFFFF );

  while ( blk != ~(unsigned long)0 )
    {
      cout << "AliHLTEMCALOnlineDisplayEventTab::FindFourierBlocks(homerReaderPtr) FOUND FOURIER DATA !!!!!!!!!!!!!!" << endl;
      blk = homerReaderPtr->FindBlockNdx("TREIRUOF","SOHP", 0xFFFFFFFF );
    }
}


void 
AliHLTEMCALOnlineDisplayEventTab::ReadBlockData(AliHLTHOMERReader *homeReaderPtr)
{  
  //  AliHLTEMCALChannelDataStruct *currentChannel =0;
  AliHLTCaloChannelDataStruct *currentChannel =0; 

 cout << "AliHLTEMCALOnlineDisplayEventTab::ReadBlockDat, Reading block data, therere are " <<  homeReaderPtr->GetBlockCnt() << " blocks " <<endl;
  FindFourierBlocks(homeReaderPtr);
  // unsigned long blk = homeReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF );
  unsigned long blk = homeReaderPtr->FindBlockNdx("TLENNAHC","SOHP", 0xFFFFFFFF );
  cout << __FILE__ << ":" << __LINE__ << "blk"  << blk  << endl ;
  int cnt = 0;
  //
  //  AliHLTEMCALCoordinate tmpCoord;
  AliHLTCaloCoordinate tmpCoord;
  
  long wcnt = 0;

  //  cout << __FILE__ << ":" << __LINE__ << "TP0"   << endl ;
  
  while ( blk != ~(unsigned long)0 ) 
    {
      //     cout << __FILE__ << ":" << __LINE__ << "TP1"    << endl ;
      //    cout << wcnt << "\t"; 

      wcnt ++;
 
      //   AliHLTEMCALChannelDataHeaderStruct* cellEnergiesPtr = (AliHLTEMCALChannelDataHeaderStruct*)homeReaderPtr->GetBlockData( blk ); 
      AliHLTCaloChannelDataHeaderStruct* cellEnergiesPtr = (AliHLTCaloChannelDataHeaderStruct*)homeReaderPtr->GetBlockData( blk );    

      Int_t* tmpPtr = 0;
      fShmPtr->SetMemory(cellEnergiesPtr);
      currentChannel = fShmPtr->NextChannel();
  
      //    cout << __FILE__ << ":" << __LINE__ << "TP1.2"    << endl ;
    
      while(currentChannel != 0)
	{
	  //	  cout << cnt << "\t"; 

	  cnt ++;
	  //	   cout << __FILE__ << ":" << __LINE__ << "TP1.3"    << endl ;
	  AliHLTEMCALMapper::ChannelId2Coordinate( currentChannel->fChannelID, tmpCoord );
	  //	  cout << __FILE__ << ":" << __LINE__ << "TP1.4"    << endl ; 
	  fgLegoPlotPtr[ tmpCoord.fGain ]->Fill(  tmpCoord.fModuleId*NXCOLUMNSMOD +   tmpCoord.fX,   tmpCoord.fZ, currentChannel->fEnergy );
	  //	  cout << __FILE__ << ":" << __LINE__ << "TP1.5"    << endl ;
	  fChannelEnergy[tmpCoord.fModuleId][tmpCoord.fZ][ tmpCoord.fX][tmpCoord.fGain] =  currentChannel->fEnergy;
	  cout << __FILE__ << ":" << __LINE__ << "TP1.6"    << endl ;
	  
	  /*
	  if(cellEnergiesPtr->fHasRawData == true)
	    {
	      FillRawData(fShmPtr->GetRawData());
	    }
	  */
	  
	  currentChannel = fShmPtr->NextChannel();

	  cout << __FILE__ << ":" << __LINE__ << "TP1.7"    << endl ;
	}
      //      blk = homeReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF, blk+1);
      //   cout << __FILE__ << ":" << __LINE__ << "TP2"    << endl ;
      blk = homeReaderPtr->FindBlockNdx("TLENNAHC","SOHP", 0xFFFFFFFF, blk+1 );
      //   cout << __FILE__ << ":" << __LINE__ << "TP3"    << endl ;
   }
}


void 
///AliHLTEMCALOnlineDisplayEventTab::FillRawData(const AliHLTEMCALChannelRawDataStruct &rawStr)
AliHLTEMCALOnlineDisplayEventTab::FillRawData(const AliHLTCaloChannelRawDataStruct &rawStr)  

{
  fNChannelSamples[ rawStr.fCoordinate.fModuleId ][ rawStr.fCoordinate.fZ ]  [ rawStr.fCoordinate.fX ][ rawStr.fCoordinate.fGain ] = rawStr.nSamplesUsed;
  fChannelEnergy[ rawStr.fCoordinate.fModuleId ][ rawStr.fCoordinate.fZ ]  [ rawStr.fCoordinate.fX ][ rawStr.fCoordinate.fGain ] = rawStr.fEnergy;


  /*
  cout << __FILE__ << __LINE__<< "module ID = " << rawStr.fCoordinate.fModuleId  << endl;
  cout << __FILE__ << __LINE__<< "fZ = " << rawStr.fCoordinate.fZ   << endl;
  cout << __FILE__ << __LINE__<< "fX = " << rawStr.fCoordinate.fX   << endl;
  cout << __FILE__ << __LINE__<< "fGain = " << rawStr.fCoordinate.fGain   << endl; 
  cout << __FILE__ << __LINE__<< "nSamples = " <<    rawStr.nSamplesUsed   << endl; 
  */

  for(int i=0; i <  rawStr.nSamplesUsed; i++ )
    {
      //   cout <<  "i = "  << i << endl;
      fChannelData[ rawStr.fCoordinate.fModuleId ][ rawStr.fCoordinate.fZ ]  [ rawStr.fCoordinate.fX ][ rawStr.fCoordinate.fGain ][i] =  rawStr.fDataPtr[i];  
      //      fChannelData[ rawStr.fCoordinate.fModuleId ][   rawStr.fCoordinate.fX ]  [ rawStr.fCoordinate.fZ ][ rawStr.fCoordinate.fGain ][i] =  rawStr.fDataPtr[i];  
    }


}



void
AliHLTEMCALOnlineDisplayEventTab::ResetDisplay()
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
AliHLTEMCALOnlineDisplayEventTab::InitDisplay(TGTab  * tabPtr, int runnumber)
{
  //  gStyle->SetOptLogy();
  ///  gStyle->SetOptStat(true);

  for(int gain=0; gain < NGAINS; gain++)
    {
      char gainLabel[100];
      char label[256];
 
      //     Gain2Text
      fOnlineDisplayPtr->Gain2Text(gain,gainLabel);
      sprintf(label, "EMCAL HLT Online Display %s", gainLabel);
      fgLegoPlotPtr[gain] = new AliHLTEMCALOnlineDisplayTH2D(fOnlineDisplayPtr, label, label, 
							    NXCOLUMNSMOD*NMODULES , 0, NXCOLUMNSMOD*NMODULES,  
							    NZROWSMOD,   0, NZROWSMOD);   
      
      //      cout << __FILE__ << ":" << __LINE__ << " Runnumber = "  << runnumber <<endl;
   
      fgLegoPlotPtr[gain]->SetRunNumber(runnumber);
      fgLegoPlotPtr[gain]->SetMaximum(1023);
      fgLegoPlotPtr[gain]->Reset();
      //     fgLegoPlotPtr[gain]->GetXaxis()->SetRange(XRANGESTART, XRANGEEND);
    }
  

  TGLayoutHints *fL1 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 15, 1);

  TGCompositeFrame * tf = tabPtr->AddTab("Event display");
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

  fgEventButtPtr = new  AliHLTEMCALGetEventButton(fSubF1, "get event", 'e');
}



void
AliHLTEMCALOnlineDisplayEventTab::UpdateDisplay()
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
