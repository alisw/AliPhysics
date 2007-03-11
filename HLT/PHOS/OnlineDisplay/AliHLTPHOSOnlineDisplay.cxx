#include  "AliHLTPHOSOnlineDisplay.h"
#include  "AliHLTDataTypes.h"
#include  "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include  <vector>
#include  "stdio.h"
#include <string>
#include <sys/ipc.h>
#include <errno.h>
#include "TH2.h"
#include "TCanvas.h"
//#include "AliHLTPHOSModuleCellAccumulatedEnergyDataStruct.h"
#include "AliHLTPHOSRcuCellAccumulatedEnergyDataStruct.h"
#include "AliHLTPHOSCommonDefs.h"

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

#ifndef __CINT__
# include <stdexcept>
# include <TSystem.h>
# include <TApplication.h>
# include "TStyle.h" 
#endif
#include <iostream>


#define MAX_BIN_VALUE 1023


AliHLTPHOSGetEventButton*  AliHLTPHOSOnlineDisplay::fgEventButtPtr     = 0;           /**<Button to get a new event from the HLT online stream*/
AliHLTPHOSOnlineDisplay*   AliHLTPHOSOnlineDisplay::fgInstancePtr      = 0;           /**<The one an only instance of PhosOnlineDisplay*/
HOMERReader*               AliHLTPHOSOnlineDisplay::fgHomerReaderPtr   = 0;           /**<Homer reader that fetches events from the HLT online stream*/
HOMERReader*               AliHLTPHOSOnlineDisplay::fgHomerReadersPtr[MAX_HOSTS];     /**<Homer reader that fetches events from the HLT online stream*/
HOMERReader*               AliHLTPHOSOnlineDisplay::fgCalibReadersPtr[MAX_HOSTS];     /**<Homer reader that fetches histograms from the HLT online stream*/
TH2D*                      AliHLTPHOSOnlineDisplay::fgLegoPlotLGPtr      = 0;         /**<2D histogram for low gain channels*/
TH2D*                      AliHLTPHOSOnlineDisplay::fgLegoPlotHGPtr      = 0;         /**<2D histogram for high gain channels*/
TH2D*                      AliHLTPHOSOnlineDisplay::fgCalibHistPtr[N_GAINS];          /**<2D histogram for low gain channels*/
TH2I*                      AliHLTPHOSOnlineDisplay::fgHitsHistPtr[N_GAINS];           /**<2D histogram for low gain channels*/
TH2D*                      AliHLTPHOSOnlineDisplay::fgAveragePtr[N_GAINS];           /**<Accumuated energy/hits*/
char*                      AliHLTPHOSOnlineDisplay::fgDefaultDet       = "SOHP";      /**<PHOS written backwards*/
char*                      AliHLTPHOSOnlineDisplay::fgDefaultDataType  = "RENELLEC";  /**<CELLENER (Celle energy) written backwards*/  
int                        AliHLTPHOSOnlineDisplay::fgEvntCnt          = 0;           /**<Event Counter*/
TCanvas*                   AliHLTPHOSOnlineDisplay::fgCanvasHGPtr      = 0;           /**<Canvas to plot fgLegoplot for High gain channels*/ 
TCanvas*                   AliHLTPHOSOnlineDisplay::fgCanvasLGPtr      = 0;           /**<Canvas to plot fgLegoplot for Low gain channels*/ 
Bool_t                     AliHLTPHOSOnlineDisplay::fgAccumulate       = kFALSE ;     /**<If set to kFALSE reset fgLegoplot between event, kTRUE adds current energies to previous plot*/
Bool_t                     AliHLTPHOSOnlineDisplay::fgSyncronize       = kFALSE ;
unsigned int               AliHLTPHOSOnlineDisplay::fgNHosts           = 0;
unsigned int               AliHLTPHOSOnlineDisplay::fgNPorts           = 0;
char*                      AliHLTPHOSOnlineDisplay::fgHosts[MAX_HOSTS];
short unsigned int*        AliHLTPHOSOnlineDisplay::fgPorts            =0; 


TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fFrame1            = 0; 
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF1                = 0;         
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF2                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF3                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF4                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF5                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fSubF1             = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fSubF2             = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fSubF3             = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fSubF4             = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fSubF5             = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fSubF6             = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fSubF7             = 0;

TGTab*                     AliHLTPHOSOnlineDisplay::fTab               = 0;
TGTab*                     AliHLTPHOSOnlineDisplay::fSubTab1           = 0;
TGTab*                     AliHLTPHOSOnlineDisplay::fSubTab2           = 0;

TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc1               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc2               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc3               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc4               = 0;
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc5               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc6               = 0;
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc7               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc8               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc9               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc10              = 0;
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc11              = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc12              = 0;
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc13              = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc14              = 0;
using namespace std;


AliHLTPHOSOnlineDisplay*
AliHLTPHOSOnlineDisplay::Instance() 
{
  if (!fgInstancePtr) fgInstancePtr = new AliHLTPHOSOnlineDisplay();
  return fgInstancePtr;
}

AliHLTPHOSOnlineDisplay::AliHLTPHOSOnlineDisplay()
{
  char **tmp;
  cout << "creating new PHOS Onlinedisplay" << endl;
  fgLegoPlotLGPtr  = 0;
  fgLegoPlotHGPtr  = 0;

  fgHomerReaderPtr = new  HOMERReader(fgNHosts, (const char**)fgHosts,  fgPorts);

  for(int i = 0; i <fgNHosts; i++)
    {
      fgHomerReadersPtr[i] =  new  HOMERReader(fgHosts[i], fgPorts[i]); 
      fgCalibReadersPtr[i] =  new  HOMERReader(fgHosts[i], fgPorts[i]);
    }
 
 InitDisplay();

}


AliHLTPHOSOnlineDisplay::~AliHLTPHOSOnlineDisplay()
{

}

void
AliHLTPHOSOnlineDisplay::InitDisplay()
{
  fgLegoPlotHGPtr = new TH2D("Homer","HLT: #pi^{0} 5 - 30Gev, High gain",  
			     N_XCOLUMNS_MOD*N_MODULES , 0, N_XCOLUMNS_MOD*N_MODULES,  
                             N_ZROWS_MOD,               0, N_ZROWS_MOD);
  fgLegoPlotHGPtr->SetMaximum( MAX_BIN_VALUE);
  fgLegoPlotHGPtr->Reset();

  fgLegoPlotLGPtr = new TH2D("Homer","HLT: #pi^{0} 5 - 30Gev, Low gain",  
			     N_XCOLUMNS_MOD* N_MODULES , 0, N_XCOLUMNS_MOD* N_MODULES,  
			     N_ZROWS_MOD,          0, N_ZROWS_MOD);
  fgLegoPlotLGPtr->SetMaximum( MAX_BIN_VALUE); 
  fgLegoPlotLGPtr->Reset();

  for(int gain = 0; gain< N_GAINS; gain ++)
    {
      fgCalibHistPtr[gain] = new TH2D("Homer","HLT:",  
				      N_XCOLUMNS_MOD*N_MODULES , 0, N_XCOLUMNS_MOD*N_MODULES , 
				      N_ZROWS_MOD,         0, N_ZROWS_MOD);
      fgCalibHistPtr[gain]->Reset(); 
     
      fgHitsHistPtr[gain] = new TH2I("Homer","HLT: #pi^{0} 5 - 30Gev",  
				    N_XCOLUMNS_MOD* N_MODULES , 0, N_XCOLUMNS_MOD* N_MODULES,  
				    N_ZROWS_MOD,          0, N_ZROWS_MOD);
      fgHitsHistPtr[gain]->SetMaximum( MAX_BIN_VALUE); 
      fgHitsHistPtr[gain]->Reset();

      fgAveragePtr[gain] = new TH2D("Homer","HLT: #pi^{0} 5 - 30Gev",  
				    N_XCOLUMNS_MOD* N_MODULES , 0, N_XCOLUMNS_MOD* N_MODULES,  
				    N_ZROWS_MOD,          0, N_ZROWS_MOD);
      fgAveragePtr[gain]->SetMaximum( MAX_BIN_VALUE); 
      fgAveragePtr[gain]->Reset();

    }

  gStyle->SetPalette(1);
  fTab = new TGTab(this, 100, 100);
  TGLayoutHints *fL1 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 15, 1);

  TGCompositeFrame *tf = fTab->AddTab("Event display");
           fSubTab1 = new TGTab(tf, 100, 100);
           TGCompositeFrame *tf2 = fSubTab1->AddTab("FGLEGO");  
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

  tf = fTab->AddTab("Calibration data");
  fF1 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
     
           fSubTab2 = new TGTab(tf, 100, 100);

	   tf2 = fSubTab2->AddTab("Accumulated energy");   
	   fSubF4 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
	   fEc7 = new TRootEmbeddedCanvas("ec7", fSubF4, 100, 100);
	   fSubF4->AddFrame(fEc7, fL1);
	   fEc8 = new TRootEmbeddedCanvas("ec8", fSubF4, 100, 100);
	   fSubF4->AddFrame(fEc8, fL1);
	   tf2->AddFrame(fSubF4, fL1);
	   
	   
	   tf2 = fSubTab2->AddTab("SCAT (hits)"); 
	   fSubF5 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
	   tf2->AddFrame(fSubF5, fL1);
	   fEc9 = new TRootEmbeddedCanvas("ec9", fSubF5, 100, 100);
	   fSubF5->AddFrame(fEc9, fL1);
	   fEc10 = new TRootEmbeddedCanvas("ec10", fSubF5, 100, 100);
	   fSubF5->AddFrame(fEc10, fL1);

	   tf2 = fSubTab2->AddTab("SURF"); 
	   fSubF6 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
	   tf2->AddFrame(fSubF6, fL1);
	   fEc11 = new TRootEmbeddedCanvas("ec11", fSubF6, 100, 100);
	   fSubF6->AddFrame(fEc11, fL1);
	   fEc12 = new TRootEmbeddedCanvas("ec12", fSubF6, 100, 100);
	   fSubF6->AddFrame(fEc12, fL1);

	   
	   tf2 = fSubTab2->AddTab("acummulated energy / hits"); 
	   fSubF7 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
	   tf2->AddFrame(fSubF7, fL1);
	   fEc13 = new TRootEmbeddedCanvas("ec13", fSubF7, 100, 100);
	   fSubF7->AddFrame(fEc13, fL1);
	   fEc14 = new TRootEmbeddedCanvas("ec14", fSubF7, 100, 100);
	   fSubF7->AddFrame(fEc14, fL1);
	   

	   fSubTab2->Resize();
	   tf->AddFrame(fSubTab2, fL1);





  tf = fTab->AddTab("Tab 3");
  fF2 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  tf->AddFrame(fF2, fL1);

  tf = fTab->AddTab("Tab 4");
  fF4 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  tf->AddFrame(fF4, fL1);

  AddFrame(fTab, fL1);
  fgEventButtPtr = new  AliHLTPHOSGetEventButton(fSubF1, "get event", 'e');
  fgEventButtPtr = new  AliHLTPHOSGetEventButton(fSubF4, "update histograms", 'h');

  MapSubwindows();
  Resize();
  SetWindowName("online display");
  MapWindow();
  MoveResize(100,100,1200,1000);
}


int
AliHLTPHOSOnlineDisplay::ScanArguments(int argc, char** argv)
{
  for(int i=0; i< MAX_HOSTS; i++)
    {
      fgHosts[i] = new char[256];
    }

  fgPorts = new short unsigned[100];

  Bool_t hostIsSet = kFALSE;
  Bool_t portIsSet = kFALSE;
  int iResult=0;
  TString argument="";

  for (int i=0; i<argc && iResult>=0; i++) 
    {
      argument=argv[i];
      
      if (argument.IsNull()) 
	{
	  continue;
	}
      
      if (argument.CompareTo("-sync")==0)
	{
	  cout << "setting Synchronize to true" << endl;
	  fgSyncronize = kTRUE;
	}  

      
      if (argument.CompareTo("-acc")==0)
	{
	  cout << "setting Accumulate to true" << endl;
	  fgAccumulate = kTRUE;
	}  
      
 
      if (argument.CompareTo("-host")==0) 
	{
	  if(i+1 <= argc)
	    {
	      i++;
	      sprintf(fgHosts[fgNHosts],"%s", argv[i]);
	      fgNHosts ++; 
	      hostIsSet = kTRUE; 
	      if(i+1 <= argc)
		{
		  argument=argv[i+1];
		  if(argument.CompareTo("-port")==0)
		  {
		    i++;
		    if(i+1 <= argc)
		      {
			i++;
			fgPorts[fgNPorts] = atoi(argv[i]);	
			cout << "A setting port to   " << fgPorts[fgNPorts]  <<endl; 
			fgNPorts ++;
			portIsSet = kTRUE;
		      }
		  }
		  
		  else
		    {
		      fgPorts[fgNPorts] =  DEFAULT_EVENT_PORT;	
		      cout << "B setting port to   " << fgPorts[fgNPorts]  <<endl; 
		      fgNPorts ++;
		      portIsSet = kTRUE;
		    }
		}

	    }

	}
    }

  if(hostIsSet != kTRUE ||  portIsSet != kTRUE)
    {
      if(hostIsSet == kFALSE)
	{
	  printf("\nERROR: no hostname is specified\n");
	}
      
      if( portIsSet == kFALSE)
	{
	  printf("ERROR: no port spcified\n");
	}
      printf("\nYou must specify at least one host \n\n");
      printf("*****************************************************************\n");
      printf("\nUsage: ./onlinedisplay  -hostname  <hostname>   -port  <port>");
      printf("\n-port is optional, if not set a default port will be assumed\n");
      printf("*****************************************************************\n\n\n");
      iResult = -1;
    }
  
  else
    {
      iResult = 0;
    }

  return iResult;
}//end ScanArguments


int
AliHLTPHOSOnlineDisplay::GetNextEvent()
{
  HOMERReader* CurrentReaderPtr;

  if(fgAccumulate == kFALSE)
    {
      cout <<"resetting fgLegoplot" << endl;
      if(fgLegoPlotHGPtr !=0)
	{
	  fgLegoPlotHGPtr->Reset(); 
	}

      if(fgLegoPlotLGPtr !=0)
	{
	  fgLegoPlotLGPtr->Reset();
	}  
    }

  int ret = 0;
  unsigned long ndx;
  const AliHLTComponentBlockData* iter = NULL;   
  Bool_t nextSwitch=kTRUE; 
  cout << "homerreader connectionstatus  =" <<fgHomerReaderPtr->GetConnectionStatus() << endl;;

  int nLoops=0;
  if(fgSyncronize == kTRUE)
    {
      nLoops = 1;
    }
  else
    {
      nLoops =  fgNHosts; 
    }
  
  for(int reader = 0; reader <  nLoops; reader ++)
    {
      if(fgSyncronize == kTRUE)
	{
	  CurrentReaderPtr =fgHomerReaderPtr;
	}
      else
	{
	  CurrentReaderPtr =fgHomerReadersPtr[reader];
	}
      ret =CurrentReaderPtr->ReadNextEvent();  
      cout << "Event ID =\t " << reader <<" = "<< CurrentReaderPtr->GetEventID() << endl;;
  
      if( ret ) 
	{
	  int ndx = fgHomerReaderPtr->GetErrorConnectionNdx();
	  printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
	  cout << "HOMER getconncetioNdx  status = " << ndx << endl;
	  return ret; 
	}
  
      unsigned long blockCnt = CurrentReaderPtr->GetBlockCnt();
      
      for ( unsigned long i = 0; i < blockCnt; i++ ) 
	{
	  char tmp1[9], tmp2[5];
	  memset( tmp1, 0, 9 );
	  memset( tmp2, 0, 5);
	  void *tmp11 = tmp1;
	  ULong64_t* tmp12 = (ULong64_t*)tmp11;
	  *tmp12 =fgHomerReadersPtr[reader]->GetBlockDataType( i );
	  void *tmp21 = tmp2;
	  ULong_t* tmp22 = (ULong_t*)tmp21;
	  *tmp22 = CurrentReaderPtr->GetBlockDataOrigin( i );
	  cout << "Dataype for block:  "<< i<<"  is:  "<< tmp1<<tmp2 <<endl;
	}

      unsigned long blk = CurrentReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF );

      while ( blk != ~(unsigned long)0 ) 
	{
	  Int_t moduleID;
	  const AliHLTPHOSRcuCellEnergyDataStruct* cellEnergiesPtr = (const AliHLTPHOSRcuCellEnergyDataStruct*) CurrentReaderPtr->GetBlockData( blk );  
	  moduleID = cellEnergiesPtr->fModuleID ;
	  Int_t tmpCount = cellEnergiesPtr->fCnt;
	  Int_t tmpZ;
	  Int_t tmpX;
	  Int_t tmpGain;
	  
	  //	  for(int i= 0; i<tmpCount; i++)
	  for(int i= 0; i <= tmpCount; i++)
	    {
	      tmpZ = cellEnergiesPtr->fValidData[i].fZ;
	      tmpX = cellEnergiesPtr->fValidData[i].fX;
	      tmpGain =  cellEnergiesPtr->fValidData[i].fGain;
	      
	      if(tmpGain == HIGH_GAIN)
		{
		  fgLegoPlotHGPtr->Fill(moduleID*N_XCOLUMNS_MOD + tmpX +  N_XCOLUMNS_RCU*cellEnergiesPtr->fRcuX,  
				      tmpZ + N_ZROWS_RCU*cellEnergiesPtr->fRcuZ, cellEnergiesPtr->fValidData[i].fEnergy);
		}
		  
	      else if(tmpGain == LOW_GAIN)
		{
		  fgLegoPlotLGPtr->Fill(moduleID*N_XCOLUMNS_MOD + tmpX +  N_XCOLUMNS_RCU*cellEnergiesPtr->fRcuX,
				      tmpZ + N_ZROWS_RCU*cellEnergiesPtr->fRcuZ,    cellEnergiesPtr->fValidData[i].fEnergy);
		}
	    }
		
	  blk = CurrentReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF, blk+1);
	}
    }
  UpdateDisplay();
  fgEvntCnt ++;
}

int
AliHLTPHOSOnlineDisplay::GetHistogram()
{
  fgCalibHistPtr[LOW_GAIN]->Reset(); 
  fgCalibHistPtr[HIGH_GAIN]->Reset();
 
  int ret = 0;
  unsigned long ndx;
  const AliHLTComponentBlockData* iter = NULL;   
  Bool_t nextSwitch=kTRUE; 

  for(int reader = 0; reader <  fgNHosts; reader ++)
    {
      ret =fgCalibReadersPtr[reader]->ReadNextEvent(); ;  
      if( ret ) 
	{
	  int ndx = fgCalibReadersPtr[reader]->GetErrorConnectionNdx();
	  printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
	  cout << "HOMER getconncetioNdx  status = " << ndx << endl;
	  return ret; 
	}
      
      unsigned long blockCnt = fgCalibReadersPtr[reader]->GetBlockCnt();
      cout << "AliHLTPHOSOnlineDisplay::GetHistogram():  blockCnt  = " << blockCnt << endl;

      for ( unsigned long i = 0; i < blockCnt; i++ ) 
	{
	  char tmp1[9], tmp2[5];
	  memset( tmp1, 0, 9 );
	  memset( tmp2, 0, 5);
	  void *tmp11 = tmp1;
	  ULong64_t* tmp12 = (ULong64_t*)tmp11;
	  *tmp12 =fgCalibReadersPtr[reader]->GetBlockDataType( i );
	  void *tmp21 = tmp2;
	  ULong_t* tmp22 = (ULong_t*)tmp21;
	  *tmp22 = fgCalibReadersPtr[reader]->GetBlockDataOrigin( i );
	  cout << "Dataype is: "<< tmp1<<"   "<<tmp2 <<endl;
	}
      
      unsigned long blk = fgCalibReadersPtr[reader]->FindBlockNdx("UCCARENE","SOHP", 0xFFFFFFFF );
      //      int tmpWhileCnt = 0;
  


      while ( blk != ~(unsigned long)0 ) 
	{
	  cout << "GetHistogram: updating block " << endl;
	  AliHLTUInt16_t moduleID;
	  const AliHLTPHOSRcuCellAccumulatedEnergyDataStruct* accCellEnergiesPtr = (const AliHLTPHOSRcuCellAccumulatedEnergyDataStruct*)fgCalibReadersPtr[reader]->GetBlockData( blk ); 
	  moduleID = accCellEnergiesPtr->fModuleID ;
	  //	  int RcuXOffset = (accCellEnergiesPtr->fRcuX)*N_XCOLUMNS_RCU;
	  //	  int RcuZOffset = (accCellEnergiesPtr->fRcuZ)*N_ZROWS_RCU;
	  cout << "(X,Z) =" << "("<< (int)accCellEnergiesPtr->fRcuX <<" , " <<  (int)accCellEnergiesPtr->fRcuZ << ") " << endl;
	  
	  int tmpx;
	  int tmpz;

	  for(int x = 0; x < N_XCOLUMNS_RCU; x ++)
	    for(int z = 0; z <N_ZROWS_RCU; z ++)
	      {
		{
		  for(int gain = 0; gain < N_GAINS; gain ++)
		    {
		      tmpx = moduleID*N_XCOLUMNS_MOD + (accCellEnergiesPtr->fRcuX)*N_XCOLUMNS_RCU + x;
		      tmpz = (accCellEnergiesPtr->fRcuZ)*N_ZROWS_RCU +z;

		      fgCalibHistPtr[gain]->Fill(tmpx, tmpz, accCellEnergiesPtr->fAccumulatedEnergies[x][z][gain] );
		      fgHitsHistPtr[gain]->Fill(tmpx, tmpz, accCellEnergiesPtr->fHits[x][z][gain] );
		      
		      if(fgHitsHistPtr[gain]->GetBinContent(tmpx, tmpz) > 0)
			{
			  fgAveragePtr[gain]->SetBinContent(tmpx, tmpz, fgCalibHistPtr[gain]->GetBinContent(tmpx, tmpz)/fgHitsHistPtr[gain]->GetBinContent(tmpx, tmpz));
			}
		    }
		}
	      }

	  blk = fgCalibReadersPtr[reader]->FindBlockNdx("UCCARENE","SOHP", 0xFFFFFFFF, blk+1);
	  //	  tmpWhileCnt ++;
	} 
    }
  
  UpdateHistograms();
  fgEvntCnt ++;
}

/*
void 
AliHLTPHOSOnlineDisplay::EvaluateAverage()
{
  for(int x = 0; x < N_XCOLUMNS_RCU; x ++)
    for(int z = 0; z <N_ZROWS_RCU; z ++)
      {
	{
	  for(int gain = 0; gain < N_GAINS; gain ++)
	    {
	      tmpx = moduleID*N_XCOLUMNS_MOD + (accCellEnergiesPtr->fRcuX)*N_XCOLUMNS_RCU + x;
	      tmpz = (accCellEnergiesPtr->fRcuZ)*N_ZROWS_RCU +z;
	      
	      fgCalibHistPtr[gain]->Fill(tmpx, tmpz, accCellEnergiesPtr->fAccumulatedEnergies[x][z][gain] );
	      fgHitsHistPtr[gain]->Fill(tmpx, tmpz, accCellEnergiesPtr->fHits[x][z][gain] );
	      
	    }
	}
      } 
}
*/



void
AliHLTPHOSOnlineDisplay::UpdateDisplay()
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

void
AliHLTPHOSOnlineDisplay::UpdateHistograms()
{
 fgCanvasHGPtr =  fEc7->GetCanvas();
  fgCanvasHGPtr->cd();
  fgCalibHistPtr[HIGH_GAIN]->Draw("LEGO2Z");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr = fEc8->GetCanvas();
  fgCanvasLGPtr->cd();
  fgCalibHistPtr[LOW_GAIN]->Draw("LEGO2Z");
  fgCanvasLGPtr->Update();

  fgCanvasHGPtr =  fEc9->GetCanvas();
  fgCanvasHGPtr->cd();
  fgHitsHistPtr[HIGH_GAIN]->Draw("SCAT");
  fgCanvasHGPtr->Update();

  fgCanvasLGPtr = fEc10->GetCanvas();
  fgCanvasLGPtr->cd();
  fgHitsHistPtr[LOW_GAIN]->Draw("SCAT");
  fgCanvasLGPtr->Update();

  fgCanvasHGPtr =  fEc11->GetCanvas();
  fgCanvasHGPtr->cd();
  fgCalibHistPtr[HIGH_GAIN]->Draw("COLZ");
  fgCanvasHGPtr->Update();

  fgCanvasLGPtr = fEc12->GetCanvas();
  fgCanvasLGPtr->cd();
  fgCalibHistPtr[LOW_GAIN]->Draw("COLZ");
  fgCanvasLGPtr->Update();


  fgCanvasLGPtr = fEc13->GetCanvas();
  fgCanvasLGPtr->cd();
  fgAveragePtr[HIGH_GAIN]->Draw("COLZ");
  fgCanvasLGPtr->Update();

  fgCanvasHGPtr = fEc14->GetCanvas();
  fgCanvasHGPtr->cd();
  fgAveragePtr[LOW_GAIN]->Draw("COLZ");
  fgCanvasHGPtr->Update();

}
