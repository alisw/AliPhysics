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
#endif
#include <iostream>


AliHLTPHOSOnlineDisplay*   AliHLTPHOSOnlineDisplay::fgInstancePtr      = 0;
HOMERReader*               AliHLTPHOSOnlineDisplay::homerReaderPtr     = 0;
TH2S*                      AliHLTPHOSOnlineDisplay::legoPlotLGPtr      = 0;
TH2S*                      AliHLTPHOSOnlineDisplay::legoPlotHGPtr      = 0;
AliHLTPHOSGetEventButton*  AliHLTPHOSOnlineDisplay::fgEventButtPtr     = 0;
char*                      AliHLTPHOSOnlineDisplay::fgDefaultDet       = "SOHP";      //PHOS written backwards
char*                      AliHLTPHOSOnlineDisplay::fgDefaultDataType  = "RENELLEC";  //CELLENER (Celle energy) written backwards  
int                        AliHLTPHOSOnlineDisplay::fgEvntCnt          = 0;           //CELLENER (Celle energy) written backwards
TCanvas*                   AliHLTPHOSOnlineDisplay::fgCanvasHGPtr[100];
TCanvas*                   AliHLTPHOSOnlineDisplay::fgCanvasLGPtr[100];



using namespace std;


AliHLTPHOSOnlineDisplay*
AliHLTPHOSOnlineDisplay::Instance() 
{
  if (!fgInstancePtr) fgInstancePtr = new AliHLTPHOSOnlineDisplay;
  return fgInstancePtr;
}


AliHLTPHOSOnlineDisplay::AliHLTPHOSOnlineDisplay()
{
  this->MapWindow();
  this->SetWindowName("GetEvent");
  this->MoveResize(100,100,1000,720);
  fgEventButtPtr = new  AliHLTPHOSGetEventButton(this, "online display");
  this->MapSubwindows();  
  homerReaderPtr = new  HOMERReader("mixing", 42001);
  std::vector<unsigned> blockList;
  cout << "creating new PHOS Onlinedisplay" << endl;
  int ret = 0;
  Bool_t nextSwitch=kTRUE;
}


AliHLTPHOSOnlineDisplay::~AliHLTPHOSOnlineDisplay()
{

}

int
AliHLTPHOSOnlineDisplay::GetNextEvent()
{
  int ret = 0;
  unsigned long ndx;
  const AliHLTComponentBlockData* iter = NULL;   
  Bool_t nextSwitch=kTRUE; 
  cout << "homerreader connectionstatus =" << homerReaderPtr->GetConnectionStatus() << endl;;

  ret = homerReaderPtr->ReadNextEvent();  
      
  if( ret ) 
    {
      int ndx =  homerReaderPtr->GetErrorConnectionNdx();
      printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
      return ret; 
  }
      
  unsigned long blockCnt =  homerReaderPtr->GetBlockCnt();

  printf( "Event 0x%016LX (%Lu) with %lu blocks\n", (ULong64_t) homerReaderPtr->GetEventID(), (ULong64_t) homerReaderPtr->GetEventID(), blockCnt );
      
  for ( unsigned long i = 0; i < blockCnt; i++ ) 
    {
      char tmp1[9], tmp2[5];
      memset( tmp1, 0, 9 );
      memset( tmp2, 0, 5);
      void *tmp11 = tmp1;
      ULong64_t* tmp12 = (ULong64_t*)tmp11;
      *tmp12 = homerReaderPtr->GetBlockDataType( i );
      void *tmp21 = tmp2;
      ULong_t* tmp22 = (ULong_t*)tmp21;
      *tmp22 = homerReaderPtr->GetBlockDataOrigin( i );
      printf( "Block %lu length: %lu - type: %s - origin: %s\n",ndx, homerReaderPtr->GetBlockDataLength( i ), tmp1, tmp2 );
	  
    }
  unsigned long blk =  homerReaderPtr->FindBlockNdx( fgDefaultDataType, fgDefaultDet, 0xFFFFFFFF );    


  fgCanvasHGPtr[fgEvntCnt] = new TCanvas();
  fgCanvasHGPtr[fgEvntCnt]->cd();
  legoPlotHGPtr   = new TH2S("Lego plot 1","Phi0 10 - 20Gev, High gain", 56*5, 0, 56*5, 64, 0, 64);
  
  fgCanvasLGPtr[fgEvntCnt] = new TCanvas();
  fgCanvasLGPtr[fgEvntCnt]->cd();
  legoPlotLGPtr   = new TH2S("Lego plot 1","Phi0 10 - 20Gev, Low gain", 56*5, 0, 56*5, 64, 0, 64);



  while ( blk != ~(unsigned long)0 ) 
    {
      AliHLTUInt16_t moduleID;
      //      printf( "Found Cell Energy block %lu\n", blk );
      const AliHLTPHOSRcuCellEnergyDataStruct* cellEnergiesPtr = (const AliHLTPHOSRcuCellEnergyDataStruct*) homerReaderPtr->GetBlockData( blk );  
      moduleID = cellEnergiesPtr->fModuleID ;
      //      cout << endl;
 
      int tmp = cellEnergiesPtr->fCnt;
      int tmpRow;
      int tmpCol;
 
  
      for(int i= 0; i<tmp; i++)
	{
	  tmpRow = cellEnergiesPtr->fValidData[i].fRow;
	  tmpCol = cellEnergiesPtr->fValidData[i].fCol;
	  legoPlotHGPtr->SetBinContent(moduleID*56 + tmpCol + 28*cellEnergiesPtr->fRcuZ,  tmpRow + 32*cellEnergiesPtr->fRcuX , cellEnergiesPtr->fCellEnergies[tmpRow][tmpCol][0]);
	  legoPlotLGPtr->SetBinContent(moduleID*56 + tmpCol + 28*cellEnergiesPtr->fRcuZ,  tmpRow + 32*cellEnergiesPtr->fRcuX , cellEnergiesPtr->fCellEnergies[tmpRow][tmpCol][1]);
	}

      blk =  homerReaderPtr->FindBlockNdx( fgDefaultDataType, fgDefaultDet, 0xFFFFFFFF, blk+1);  
    }
  
  //  cout << "event count =" << fgEvntCnt<<endl;

  fgCanvasLGPtr[fgEvntCnt]->cd();
  legoPlotLGPtr->Draw("LEGO2");

  fgCanvasHGPtr[fgEvntCnt]->cd();
  legoPlotHGPtr->Draw("LEGO2");


  fgEvntCnt ++;
}
