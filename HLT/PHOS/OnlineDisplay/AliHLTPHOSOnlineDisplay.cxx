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
Bool_t                     AliHLTPHOSOnlineDisplay::fgAccumulate       = kFALSE ;     /**<If set to kFALSE reset legoplot between event, kTRUE adds current energies to previous plot*/


using namespace std;


AliHLTPHOSOnlineDisplay*
AliHLTPHOSOnlineDisplay::Instance() 
{
  if (!fgInstancePtr) fgInstancePtr = new AliHLTPHOSOnlineDisplay;
  return fgInstancePtr;
}


AliHLTPHOSOnlineDisplay::AliHLTPHOSOnlineDisplay()
{
  legoPlotLGPtr  = 0;
  legoPlotHGPtr  = 0;
  MapWindow();
  SetWindowName("online display");
  MoveResize(100,100,1000,720);
  fgEventButtPtr = new  AliHLTPHOSGetEventButton(this, "get event");
  MapSubwindows();
  fgHomerReaderPtr = new  HOMERReader("mixing", 42001);
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
  if(fgEvntCnt == 0)
    {
      fgCanvasLGPtr    = new TCanvas();
      fgCanvasHGPtr    = new TCanvas();
      legoPlotHGPtr   = new TH2S("Homer","HLT:HOMER: #pi^{0} 5 - 30Gev, High gain",  N_COLUMNS_MOD* N_MODULES , 0,  N_COLUMNS_MOD* N_MODULES ,  N_ROWS_MOD, 0,  N_ROWS_MOD);
      legoPlotHGPtr->SetMaximum( MAX_BIN_VALUE);
      legoPlotLGPtr   = new TH2S("Homer","HLT:HOMER: #pi^{0} 5 - 30Gev, Low gain",  N_COLUMNS_MOD* N_MODULES , 0,  N_COLUMNS_MOD* N_MODULES ,  N_ROWS_MOD, 0,  N_ROWS_MOD);
      legoPlotLGPtr->SetMaximum( MAX_BIN_VALUE);
    }  


  if(fgAccumulate == kFALSE)
    {
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
  cout << "homerreader connectionstatus =" <<fgHomerReaderPtr->GetConnectionStatus() << endl;;

  ret =fgHomerReaderPtr->ReadNextEvent();  
      
  if( ret ) 
    {
      int ndx = fgHomerReaderPtr->GetErrorConnectionNdx();
      printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
      return ret; 
  }
      
  unsigned long blockCnt = fgHomerReaderPtr->GetBlockCnt();

  printf( "Event 0x%016LX (%Lu) with %lu blocks\n", (ULong64_t)fgHomerReaderPtr->GetEventID(), (ULong64_t)fgHomerReaderPtr->GetEventID(), blockCnt );
      
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
      printf( "Block %lu length: %lu - type: %s - origin: %s\n",ndx,fgHomerReaderPtr->GetBlockDataLength( i ), tmp1, tmp2 );
	  
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
	        tmpBin = legoPlotHGPtr->GetBin(moduleID*N_COLUMNS_MOD + tmpCol +  N_COLUMNS_RCU*cellEnergiesPtr->fRcuZ,  tmpRow + N_ROWS_RCU*cellEnergiesPtr->fRcuX );
		legoPlotHGPtr->AddBinContent(tmpBin, cellEnergiesPtr->fValidData[i].fEnergy);
	
	    }

	  else if(tmpGain == LOW_GAIN)
	    {
	      tmpBin = legoPlotLGPtr->GetBin(moduleID*N_COLUMNS_MOD + tmpCol +  N_COLUMNS_RCU*cellEnergiesPtr->fRcuZ,  tmpRow + N_ROWS_RCU*cellEnergiesPtr->fRcuX );
	      legoPlotLGPtr->AddBinContent(tmpBin ,cellEnergiesPtr->fValidData[i].fEnergy);
	    }

	}

      blk = fgHomerReaderPtr->FindBlockNdx( fgDefaultDataType, fgDefaultDet, 0xFFFFFFFF, blk+1);  
    }
  
  fgCanvasHGPtr->cd();
  legoPlotHGPtr->Draw("LEGO2");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr->cd();
  legoPlotLGPtr->Draw("LEGO2");
  fgCanvasLGPtr->Update();

  fgEvntCnt ++;
}




