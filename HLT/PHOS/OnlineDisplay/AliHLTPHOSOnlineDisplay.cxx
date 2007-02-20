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
HOMERReader*               AliHLTPHOSOnlineDisplay::fgHomerReadersPtr[MAX_HOSTS];
TH2S*                      AliHLTPHOSOnlineDisplay::legoPlotLGPtr      = 0;           /**<2D histogram for low gain channels*/
TH2S*                      AliHLTPHOSOnlineDisplay::legoPlotHGPtr      = 0;           /**<2D histogram for high gain channels*/
char*                      AliHLTPHOSOnlineDisplay::fgDefaultDet       = "SOHP";      /**<PHOS written backwards*/
char*                      AliHLTPHOSOnlineDisplay::fgDefaultDataType  = "RENELLEC";  /**<CELLENER (Celle energy) written backwards*/  
int                        AliHLTPHOSOnlineDisplay::fgEvntCnt          = 0;           /**<Event Counter*/
TCanvas*                   AliHLTPHOSOnlineDisplay::fgCanvasHGPtr      = 0;           /**<Canvas to plot legoplot for High gain channels*/ 
TCanvas*                   AliHLTPHOSOnlineDisplay::fgCanvasLGPtr      = 0;           /**<Canvas to plot legoplot for Low gain channels*/ 
Bool_t                     AliHLTPHOSOnlineDisplay::fgAccumulate       = kTRUE ;     /**<If set to kFALSE reset legoplot between event, kTRUE adds current energies to previous plot*/
char*                      AliHLTPHOSOnlineDisplay::host               = 0;
int                        AliHLTPHOSOnlineDisplay::port               = 0;
unsigned int                        AliHLTPHOSOnlineDisplay::fgNHosts           = 0;
unsigned int                        AliHLTPHOSOnlineDisplay::fgNPorts           = 0;
char*                      AliHLTPHOSOnlineDisplay::fgHosts[MAX_HOSTS];
short unsigned int*                        AliHLTPHOSOnlineDisplay::fgPorts             =0; 

TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fFrame1            = 0; 
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF1                = 0;         
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF2                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF3                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF4                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF5                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fSubF1             = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fSubF2             = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fSubF3             = 0;
TGTab*                     AliHLTPHOSOnlineDisplay::fTab               = 0;
TGTab*                     AliHLTPHOSOnlineDisplay::fSubTab            = 0;
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc1               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc2               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc3               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc4               = 0;
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc5               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc6               = 0;
using namespace std;


AliHLTPHOSOnlineDisplay*
AliHLTPHOSOnlineDisplay::Instance() 
{
  if (!fgInstancePtr) fgInstancePtr = new AliHLTPHOSOnlineDisplay(host, port);
  return fgInstancePtr;
}


AliHLTPHOSOnlineDisplay::AliHLTPHOSOnlineDisplay()
{
  //  if (!fgInstancePtr) fgInstancePtr = new AliHLTPHOSOnlineDisplay(host, port);
  //  return fgInstancePtr;
  //  cout << "ERROR: You canot create Onlinedisplay without parameters" << endl;
  //  cout << "Usage: AliHLTPHOSOnlineDisplay(char *hostname, int port)" << endl;
}

AliHLTPHOSOnlineDisplay::AliHLTPHOSOnlineDisplay(char *hostname, int port)
{
  char **tmp;
  cout << "creating new PHOS Onlinedisplay" << endl;
  legoPlotLGPtr  = 0;
  legoPlotHGPtr  = 0;
  
  //fgHomerReaderPtr = new  HOMERReader(hostname, port);

  fgHomerReaderPtr = new  HOMERReader(fgNHosts, (const char**)fgHosts,  fgPorts);



  //    fgHomerReaderPtr = new  HOMERReader(fgNHosts, fgHosts,  fgPorts);
  for(int i = 0; i <fgNHosts; i++)
    {
      fgHomerReadersPtr[i] =  new  HOMERReader(fgHosts[i], fgPorts[i]); 
    }

 
 InitDisplay();
  int ret = 0;
  Bool_t nextSwitch=kTRUE;
  cout << "The number of ports is" << fgNPorts  << endl;
 cout << "The number of HOSTS is" << fgNHosts  << endl;
 
 for(int i = 0; i< fgNPorts; i++)
   {
     cout << "Port[" << i <<"] = " << fgPorts[i] << endl;  
   }

  for(int i = 0; i< fgNHosts; i++)
   {
     cout << "Host[" << i <<"] = " << fgHosts[i] << endl;  
   }

}


AliHLTPHOSOnlineDisplay::~AliHLTPHOSOnlineDisplay()
{

}

void
AliHLTPHOSOnlineDisplay::InitDisplay()
{
  gStyle->SetPalette(1);
  fTab = new TGTab(this, 100, 100);
  TGLayoutHints *fL1 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 15, 1);


   TGCompositeFrame *tf = fTab->AddTab("Event display");
           fSubTab = new TGTab(tf, 100, 100);
           TGCompositeFrame *tf2 = fSubTab->AddTab("LEGO");  

	   fSubF1 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
	   fEc1 = new TRootEmbeddedCanvas("ec1", fSubF1, 100, 100);
	   fSubF1->AddFrame(fEc1, fL1);

	   fEc2 = new TRootEmbeddedCanvas("ec1", fSubF1, 100, 100);
	   fgCanvasHGPtr = fEc2->GetCanvas();
	   fSubF1->AddFrame(fEc2, fL1);
	   tf2->AddFrame(fSubF1, fL1);

	   tf2 = fSubTab->AddTab("SCAT"); 
	   fSubF2 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
	   tf2->AddFrame(fSubF2, fL1);
	   fEc3 = new TRootEmbeddedCanvas("ec1", fSubF2, 100, 100);
	   fSubF2->AddFrame(fEc3, fL1);
	   fEc4 = new TRootEmbeddedCanvas("ec1", fSubF2, 100, 100);
	   fSubF2->AddFrame(fEc4, fL1);


	   tf2 = fSubTab->AddTab("SURF"); 
	   fSubF3 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
	   tf2->AddFrame(fSubF3, fL1);
	   fEc5 = new TRootEmbeddedCanvas("ec1", fSubF3, 100, 100);
	   fSubF3->AddFrame(fEc5, fL1);
	   fEc6 = new TRootEmbeddedCanvas("ec1", fSubF3, 100, 100);
	   fSubF3->AddFrame(fEc6, fL1);
	   fSubTab->Resize();
	   tf->AddFrame(fSubTab, fL1);

  tf = fTab->AddTab("Tab 2");
  fF1 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);

  tf = fTab->AddTab("Tab 3");
  fF2 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  tf->AddFrame(fF2, fL1);


  tf = fTab->AddTab("Tab 4");
  fF4 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  tf->AddFrame(fF4, fL1);
  //  tf->AddFrame(tf2, fL1);

  AddFrame(fTab, fL1);

  //fgEventButtPtr = new  AliHLTPHOSGetEventButton(fF1, "get event");
  fgEventButtPtr = new  AliHLTPHOSGetEventButton(fSubF1, "get event");

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
  //  int hostCnt = 0; 
 
  int bMissingParam=0;
  printf("Main: the number of argumnets is %d \n", argc);
  for (int i=0; i<argc && iResult>=0; i++) 
    {
      argument=argv[i];
 
      if (argument.IsNull()) 
	{
	  continue;
	}

      if (argument.CompareTo("-host")==0) 
	{
	  if(i+1 <= argc)
	    {
	      i++;
	      host = argv[i];
	      cout << "sprintf " <<host <<endl;
	      sprintf(fgHosts[fgNHosts],"%s", argv[i]);
	      cout << "finnished sprintf " <<host <<endl; 
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

			fgPorts[fgNPorts] =  DEFAULT_PORT;	
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
      printf("\nYou must specify hostname & port as command line arguments\n\n");
      printf("*****************************************************************\n");
      printf("\nUsage: ./onlinedisplay  -hostname  <hostname>   -port  <port>\n\n");
      printf("*****************************************************************\n\n\n");
      iResult = -1;
     
    }
  
  else
    {
      iResult = 0;
    }

  return iResult;
}


int
AliHLTPHOSOnlineDisplay::GetNextEvent()
{
  int whileCnt = 0;

  if(fgEvntCnt == 0)
    {
      legoPlotHGPtr   = new TH2S("Homer","HLT:HOMER: #pi^{0} 5 - 30Gev, High gain",  N_COLUMNS_MOD* N_MODULES , 0, 
				 N_COLUMNS_MOD* N_MODULES ,  N_ROWS_MOD, 0,  N_ROWS_MOD);
      legoPlotHGPtr->SetMaximum( MAX_BIN_VALUE);
      legoPlotLGPtr   = new TH2S("Homer","HLT:HOMER: #pi^{0} 5 - 30Gev, Low gain",  N_COLUMNS_MOD* N_MODULES , 0,  
				 N_COLUMNS_MOD* N_MODULES ,  N_ROWS_MOD, 0,  N_ROWS_MOD);
      legoPlotLGPtr->SetMaximum( MAX_BIN_VALUE);
    }  


  if(fgAccumulate == kFALSE)
    {
      cout <<"restting legoplot" << endl;
      if(legoPlotHGPtr !=0)
	{
	  legoPlotHGPtr->Reset(); 
	}

      if(legoPlotLGPtr !=0)
	{
	  legoPlotLGPtr->Reset();
	}  
    }


  int ret = 0;
  unsigned long ndx;
  const AliHLTComponentBlockData* iter = NULL;   
  Bool_t nextSwitch=kTRUE; 
  
  cout << "homerreader connectionstatus  =" <<fgHomerReaderPtr->GetConnectionStatus() << endl;;

  for(int reader = 0; reader <  fgNHosts; reader ++)
    {
      ret =fgHomerReadersPtr[reader]->ReadNextEvent();  
     
      cout << "Event ID for reader " << reader <<" = "<< fgHomerReadersPtr[reader]->GetEventID() << endl;;
  
      if( ret ) 
	{
	  int ndx = fgHomerReaderPtr->GetErrorConnectionNdx();
	  printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
	  cout << "HOMER getconncetioNdx  status = " << ndx << endl;
	  return ret; 
	}
  
      
      unsigned long blockCnt = fgHomerReadersPtr[reader]->GetBlockCnt();

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
	  *tmp22 =fgHomerReadersPtr[reader]->GetBlockDataOrigin( i );
	}
      unsigned long blk = fgHomerReadersPtr[reader]->FindBlockNdx( fgDefaultDataType, fgDefaultDet, 0xFFFFFFFF );    

  
      while ( blk != ~(unsigned long)0 ) 
	{
	  AliHLTUInt16_t moduleID;
	  const AliHLTPHOSRcuCellEnergyDataStruct* cellEnergiesPtr = (const AliHLTPHOSRcuCellEnergyDataStruct*)fgHomerReadersPtr[reader]->GetBlockData( blk );  
	  moduleID = cellEnergiesPtr->fModuleID ;
	  int tmpCount = cellEnergiesPtr->fCnt;
	  int tmpRow;
	  int tmpCol;
	  int tmpGain;
	  Int_t tmpBin;
	  
	  for(int i= 0; i<tmpCount; i++)
	    {
	      tmpRow = cellEnergiesPtr->fValidData[i].fRow;
	      tmpCol = cellEnergiesPtr->fValidData[i].fCol;
	      tmpGain =  cellEnergiesPtr->fValidData[i].fGain;

	      if(tmpGain == HIGH_GAIN)
		{
		  legoPlotHGPtr->Fill(moduleID*N_COLUMNS_MOD + tmpCol +  N_COLUMNS_RCU*cellEnergiesPtr->fRcuZ,  tmpRow + N_ROWS_RCU*cellEnergiesPtr->fRcuX, cellEnergiesPtr->fValidData[i].fEnergy);
		}
	      
	      else if(tmpGain == LOW_GAIN)
		{
		  legoPlotLGPtr->Fill(moduleID*N_COLUMNS_MOD + tmpCol +  N_COLUMNS_RCU*cellEnergiesPtr->fRcuZ, tmpRow + N_ROWS_RCU*cellEnergiesPtr->fRcuX,    cellEnergiesPtr->fValidData[i].fEnergy);
		}

	    }

	  blk = fgHomerReadersPtr[reader]->FindBlockNdx( fgDefaultDataType, fgDefaultDet, 0xFFFFFFFF, blk+1);  
      
	  whileCnt ++;

	}
 
    }

  UpdateDisplay();

  fgEvntCnt ++;
}

int
AliHLTPHOSOnlineDisplay::GetNextEvent2()
{
  int whileCnt = 0;

  if(fgEvntCnt == 0)
    {
      legoPlotHGPtr   = new TH2S("Homer","HLT:HOMER: #pi^{0} 5 - 30Gev, High gain",  N_COLUMNS_MOD* N_MODULES , 0, 
				 N_COLUMNS_MOD* N_MODULES ,  N_ROWS_MOD, 0,  N_ROWS_MOD);
      legoPlotHGPtr->SetMaximum( MAX_BIN_VALUE);
      legoPlotLGPtr   = new TH2S("Homer","HLT:HOMER: #pi^{0} 5 - 30Gev, Low gain",  N_COLUMNS_MOD* N_MODULES , 0,  
				 N_COLUMNS_MOD* N_MODULES ,  N_ROWS_MOD, 0,  N_ROWS_MOD);
      legoPlotLGPtr->SetMaximum( MAX_BIN_VALUE);
    }  


  if(fgAccumulate == kFALSE)
    {
      cout <<"restting legoplot" << endl;
      if(legoPlotHGPtr !=0)
	{
	  legoPlotHGPtr->Reset(); 
	}

      if(legoPlotLGPtr !=0)
	{
	  legoPlotLGPtr->Reset();
	}  
    }


  int ret = 0;
  unsigned long ndx;
  const AliHLTComponentBlockData* iter = NULL;   
  Bool_t nextSwitch=kTRUE; 
  
  cout << "homerreader connectionstatus sync=" <<fgHomerReaderPtr->GetConnectionStatus() << endl;;

  //  for(int reader = 0; reader <  fgNHosts; reader ++)
  //   {
  ret =fgHomerReaderPtr->ReadNextEvent();  
     
  //  cout << "Event ID for reader " << reader <<" = "<< fgHomerReaderPtr->GetEventID() << endl;;
  
  if( ret ) 
    {
      int ndx = fgHomerReaderPtr->GetErrorConnectionNdx();
      printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
      cout << "HOMER getconncetioNdx status sync= " << ndx << endl;
      return ret; 
    }
  
      
  unsigned long blockCnt = fgHomerReaderPtr->GetBlockCnt();

  for ( unsigned long i = 0; i < blockCnt; i++ ) 
    {
      char tmp1[9], tmp2[5];
      memset( tmp1, 0, 9 );
      memset( tmp2, 0, 5);
      void *tmp11 = tmp1;
      ULong64_t* tmp12 = (ULong64_t*)tmp11;
      *tmp12 =fgHomerReaderPtr->GetBlockDataType( i );
      void *tmp21 = tmp2;
      ULong_t* tmp22 = (ULong_t*)tmp21;
      *tmp22 =fgHomerReaderPtr->GetBlockDataOrigin( i );
    }
  unsigned long blk = fgHomerReaderPtr->FindBlockNdx( fgDefaultDataType, fgDefaultDet, 0xFFFFFFFF );    

  
  while ( blk != ~(unsigned long)0 ) 
    {
      AliHLTUInt16_t moduleID;
      const AliHLTPHOSRcuCellEnergyDataStruct* cellEnergiesPtr = (const AliHLTPHOSRcuCellEnergyDataStruct*)fgHomerReaderPtr->GetBlockData( blk );  
      moduleID = cellEnergiesPtr->fModuleID ;
      int tmpCount = cellEnergiesPtr->fCnt;
      int tmpRow;
      int tmpCol;
      int tmpGain;
      Int_t tmpBin;

      for(int i= 0; i<tmpCount; i++)
	{
	  tmpRow = cellEnergiesPtr->fValidData[i].fRow;
	  tmpCol = cellEnergiesPtr->fValidData[i].fCol;
	  tmpGain =  cellEnergiesPtr->fValidData[i].fGain;

	  if(tmpGain == HIGH_GAIN)
	    {
	      legoPlotHGPtr->Fill(moduleID*N_COLUMNS_MOD + tmpCol +  N_COLUMNS_RCU*cellEnergiesPtr->fRcuZ,  tmpRow + N_ROWS_RCU*cellEnergiesPtr->fRcuX, cellEnergiesPtr->fValidData[i].fEnergy);
	    }

	  else if(tmpGain == LOW_GAIN)
	    {
	      legoPlotLGPtr->Fill(moduleID*N_COLUMNS_MOD + tmpCol +  N_COLUMNS_RCU*cellEnergiesPtr->fRcuZ, tmpRow + N_ROWS_RCU*cellEnergiesPtr->fRcuX,    cellEnergiesPtr->fValidData[i].fEnergy);
	    }

	}

      blk = fgHomerReaderPtr->FindBlockNdx( fgDefaultDataType, fgDefaultDet, 0xFFFFFFFF, blk+1);  
      
      //    whileCnt ++;

      //   }
 
    }

  UpdateDisplay();

  fgEvntCnt ++;
}


void
AliHLTPHOSOnlineDisplay::UpdateDisplay()
{
  fgCanvasHGPtr =  fEc1->GetCanvas();
  fgCanvasHGPtr->cd();
  legoPlotHGPtr->Draw("LEGO2Z");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr = fEc2->GetCanvas();
  fgCanvasLGPtr->cd();
  legoPlotLGPtr->Draw("LEGO2Z");
  fgCanvasLGPtr->Update();

  fgCanvasHGPtr =  fEc3->GetCanvas();
  fgCanvasHGPtr->cd();
  legoPlotHGPtr->Draw("SCAT");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr = fEc4->GetCanvas();
  fgCanvasLGPtr->cd();
  legoPlotLGPtr->Draw("SCAT");
  fgCanvasLGPtr->Update();

  fgCanvasHGPtr =  fEc5->GetCanvas();
  fgCanvasHGPtr->cd();
  legoPlotHGPtr->Draw("CONTZ");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr = fEc6->GetCanvas();
  fgCanvasLGPtr->cd();
  legoPlotLGPtr->Draw("CONTZ");
  fgCanvasLGPtr->Update();
}

