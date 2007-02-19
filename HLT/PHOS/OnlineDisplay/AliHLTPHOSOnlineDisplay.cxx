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
TH2S*                      AliHLTPHOSOnlineDisplay::legoPlotLGPtr      = 0;           /**<2D histogram for low gain channels*/
TH2S*                      AliHLTPHOSOnlineDisplay::legoPlotHGPtr      = 0;           /**<2D histogram for high gain channels*/
char*                      AliHLTPHOSOnlineDisplay::fgDefaultDet       = "SOHP";      /**<PHOS written backwards*/
char*                      AliHLTPHOSOnlineDisplay::fgDefaultDataType  = "RENELLEC";  /**<CELLENER (Celle energy) written backwards*/  
int                        AliHLTPHOSOnlineDisplay::fgEvntCnt          = 0;           /**<Event Counter*/
TCanvas*                   AliHLTPHOSOnlineDisplay::fgCanvasHGPtr      = 0;           /**<Canvas to plot legoplot for High gain channels*/ 
TCanvas*                   AliHLTPHOSOnlineDisplay::fgCanvasLGPtr      = 0;           /**<Canvas to plot legoplot for Low gain channels*/ 
Bool_t                     AliHLTPHOSOnlineDisplay::fgAccumulate       = kTRUE ;     /**<If set to kFALSE reset legoplot between event, kTRUE adds current energies to previous plot*/
Bool_t                     AliHLTPHOSOnlineDisplay::test[17920][2];
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fFrame1            = 0; 
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF1                = 0;         
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF2                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF3                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF4                = 0;
TGCompositeFrame*          AliHLTPHOSOnlineDisplay::fF5                = 0;
TGTab*                     AliHLTPHOSOnlineDisplay::fTab               = 0;
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc1               = 0; 
TRootEmbeddedCanvas*       AliHLTPHOSOnlineDisplay::fEc2               = 0; 
using namespace std;


AliHLTPHOSOnlineDisplay*
AliHLTPHOSOnlineDisplay::Instance(char *hostname, int port) 
{
  if (!fgInstancePtr) fgInstancePtr = new AliHLTPHOSOnlineDisplay(hostname, port);
  return fgInstancePtr;
}


AliHLTPHOSOnlineDisplay::AliHLTPHOSOnlineDisplay()
{
  cout << "ERROR: You canot create Onlinedisplay without parameters" << endl;
  cout << "Usage: AliHLTPHOSOnlineDisplay(char *hostname, int port)" << endl;
}

AliHLTPHOSOnlineDisplay::AliHLTPHOSOnlineDisplay(char *hostname, int port)
{
  cout << "creating new PHOS Onlinedisplay" << endl;
  legoPlotLGPtr  = 0;
  legoPlotHGPtr  = 0;
  fgHomerReaderPtr = new  HOMERReader(hostname, port);
  InitDisplay();
  int ret = 0;
  Bool_t nextSwitch=kTRUE;
}


AliHLTPHOSOnlineDisplay::~AliHLTPHOSOnlineDisplay()
{

}

void
AliHLTPHOSOnlineDisplay::InitDisplay()
{
  gStyle->SetPalette(1);

  fTab = new TGTab(this, 100, 100);
  TGCompositeFrame *tf = fTab->AddTab("Tab 1");
  fF1 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  TGLayoutHints *fL1 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 15, 1);
  fEc1 = new TRootEmbeddedCanvas("ec1", fF1, 100, 100);
  fgCanvasHGPtr = fEc1->GetCanvas();
  fF1->AddFrame(fEc1, fL1);
  fEc2 = new TRootEmbeddedCanvas("ec2", fF1, 100, 100);
  fgCanvasLGPtr = fEc2->GetCanvas();
  fF1->AddFrame(fEc2, fL1);
  tf->AddFrame(fF1, fL1);

  tf = fTab->AddTab("Tab 2");
  fF2 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  tf->AddFrame(fF2, fL1);

  tf = fTab->AddTab("Tab 3");
  fF3 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  tf->AddFrame(fF3, fL1);

  tf = fTab->AddTab("Tab 4");
  fF4 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  tf->AddFrame(fF4, fL1);

  AddFrame(fTab, fL1);
  fgEventButtPtr = new  AliHLTPHOSGetEventButton(fF1, "get event");
  MapSubwindows();
  Resize();
  SetWindowName("online display");
  MapWindow();
  MoveResize(100,100,1200,1000);
}


int
AliHLTPHOSOnlineDisplay::GetNextEvent()
{
  int whileCnt = 0;

  if(fgEvntCnt == 0)
    {
      legoPlotHGPtr   = new TH2S("Homer","HLT:HOMER: #pi^{0} 5 - 30Gev, High gain",  N_COLUMNS_MOD* N_MODULES , 0,  N_COLUMNS_MOD* N_MODULES ,  N_ROWS_MOD, 0,  N_ROWS_MOD);
      legoPlotHGPtr->SetMaximum( MAX_BIN_VALUE);
      legoPlotLGPtr   = new TH2S("Homer","HLT:HOMER: #pi^{0} 5 - 30Gev, Low gain",  N_COLUMNS_MOD* N_MODULES , 0,  N_COLUMNS_MOD* N_MODULES ,  N_ROWS_MOD, 0,  N_ROWS_MOD);
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
  //  cout << "homerreader connectionstatus =" <<fgHomerReaderPtr->GetConnectionStatus() << endl;;

  ret =fgHomerReaderPtr->ReadNextEvent();  
      
  if( ret ) 
    {
      int ndx = fgHomerReaderPtr->GetErrorConnectionNdx();
      printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
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
      
      whileCnt ++;

    }
 
  UpdateDisplay();

  fgEvntCnt ++;
}


void
AliHLTPHOSOnlineDisplay::UpdateDisplay()
{
  fgCanvasHGPtr->cd();
  //  fEc1->cd();
  legoPlotHGPtr->Draw("LEGO2");
  //  legoPlotHGPtr->Draw("COLZ");
  legoPlotHGPtr->SetMarkerColor(4);
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr->cd();
  legoPlotLGPtr->Draw("LEGO2");

  //  legoPlotLGPtr->Draw("COLZ");
  fgCanvasLGPtr->Update();

}

