#include  "AliHLTPHOSOnlineDisplayRawTab.h"
#include <iostream>
#include "AliHLTPHOSRcuCellAccumulatedEnergyDataStruct.h"
#include "AliHLTPHOSGetEventButton.h"
#include "AliHLTPHOSRcuChannelDataStruct.h"

using namespace std;

AliHLTPHOSOnlineDisplayRawTab::AliHLTPHOSOnlineDisplayRawTab()
{
  cout << "AliHLTPHOSOnlineDisplayRawTab:ERROR: You cannot create a onlinedisplay Tab without arguments" << endl;
}

AliHLTPHOSOnlineDisplayRawTab::AliHLTPHOSOnlineDisplayRawTab(TGTab  *tabPtr, HOMERReader *homerSyncPtr, HOMERReader *homerPtrs[MAX_HOSTS], int nHosts)
{
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
  InitDisplay(tabPtr);
}


AliHLTPHOSOnlineDisplayRawTab::~AliHLTPHOSOnlineDisplayRawTab()
{


}


void
AliHLTPHOSOnlineDisplayRawTab::ReadBlockData(HOMERReader *homerReaderPtr)
{
  unsigned long blk = homerReaderPtr->FindBlockNdx("ATADNAHC","SOHP", 0xeFFFFFFF );
  
  while ( blk != ~(unsigned long)0 ) 
    {
      cout << "GetNextEventRaw(): updating block " << endl;
      AliHLTUInt16_t moduleID;
      const AliHLTPHOSRcuChannelDataStruct* rcuChannelDataPtr = (const AliHLTPHOSRcuChannelDataStruct*)homerReaderPtr->GetBlockData( blk ); 
      moduleID = rcuChannelDataPtr->fModuleID ;
      cout << "GetNextEventRaw() Module ID =" << moduleID << endl; 
      int tmpx;
      int tmpz;
      AliHLTUInt32_t tmpChCnt =0;
      AliHLTUInt16_t tmpSampleCnt =0;
      
      tmpChCnt = rcuChannelDataPtr->fNValidChannels;
      cout << "tmpChCnt = " << tmpChCnt << endl; 
      
      for( AliHLTUInt32_t ch =0; ch < tmpChCnt; ch ++)
	{
	  {
	    tmpz =  rcuChannelDataPtr->fValidData[ch].fZ;
	    tmpx =  rcuChannelDataPtr->fValidData[ch].fX;	
	    tmpSampleCnt =  rcuChannelDataPtr->fValidData[ch].fDataSize -2;

	    //	    fgChannelDataPlotPtr[tmpz][tmpx]->SetFillColor(3);

	    fgChannelDataPlotPtr[tmpz][tmpx]->SetFillColor(1);

	    for(AliHLTUInt16_t sample =0; sample <tmpSampleCnt; sample ++)
	      {
		if(  rcuChannelDataPtr->fValidData[ch].fGain  == 0)
		  {
		    fgChannelDataPlotPtr[tmpz][tmpx]->SetBinContent(sample,  rcuChannelDataPtr->fValidData[ch].fChannelData[sample]);
		  }

	      }
	  }
	}

      blk =  homerReaderPtr->FindBlockNdx("ATADNAHC","SOHP", 0xeFFFFFFF, blk+1);

    }
}


int 
AliHLTPHOSOnlineDisplayRawTab::GetNextEvent()
{
  cout << "AliHLTPHOSOnlineDisplayRawTab::GetNextEvent()" << endl;
  ResetDisplay();
  DoGetNextEvent();
  UpdateDisplay();
  fgEvntCnt ++;
}


void
AliHLTPHOSOnlineDisplayRawTab::ResetDisplay()
{

}


void 
AliHLTPHOSOnlineDisplayRawTab::InitDisplay(TGTab *tabPtr)
{
  char tmpHistoName[256]; 
  //  char tmpHistoName[256];
  char tmpChDtaName[256];

  for(int z = 0; z < N_ZROWS_RCU; z ++)
    {
      for(int x = 0; x < N_XCOLUMNS_RCU; x ++)
	{
	  sprintf(tmpHistoName, "row %d column %d",z, x);
	  fgChannelDataPlotPtr[z][x] = new TH1D(tmpHistoName, tmpHistoName, 300, 0, 299);
	  //	  fgChannelDataPlotPtr[z][x] = new AliHLTPHOSTH1D(tmpHistoName, tmpHistoName, 100, 0, 100);
	  //	  fgChannelDataPlotPtr[z][x]->SetMaximum(MAX_BIN_VALUE); 
	  fgChannelDataPlotPtr[z][x]->SetMaximum(1023); 
	  fgChannelDataPlotPtr[z][x]->Reset();
	}
    }



  TGCompositeFrame  *tf = tabPtr->AddTab("Raw Data Dislay");
 
  fSubTab3 = new TGTab(tf, 100, 100);
  TGLayoutHints *hints = new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0); 
  char tmpTabName[256];
  char tmpCanvasName[256];
  sprintf(tmpTabName, "Raw data zzz");
  TGCompositeFrame	   *tf2 = fSubTab3->AddTab(tmpTabName);   
  fgEventButtPtr = new  AliHLTPHOSGetEventButton(tf, "Get Rawdata2xxxxx", 'r');
  AliHLTPHOSGetEventButton*	   EventButtPtr2 = new  AliHLTPHOSGetEventButton(tf, "Get Rawdata", 'r'); 
  EventButtPtr2->Move(200, 200); 

}


void 
AliHLTPHOSOnlineDisplayRawTab::UpdateDisplay()
{
  fgTestCanvasPtr = new TCanvas("TEST", "testcanvas", 1200, 1000);  

  int startZ = 25;
  int endZ   = 31;
  int startX = 20;
  int endX   = 27;

  int nRows =  endZ - startZ;
  int nCols =  endX - startX;


  //  int start = 0;
  //  fgTestCanvasPtr->Divide(N_XCOLUMNS_RCU, N_ZROWS_RCU, 0, 0);
  cout <<"nRows = "<< nRows << endl;
  cout <<"nCols = "<< nCols << endl;


  /*
  fgTestCanvasPtr->Divide(nCols, nRows, 0, 0);
  //  fgTestCanvasPtr->Divide(nRows,  nCols,  0, 0);
  
   
  for(int z = startZ; z < endZ; z ++)
    {
      for(int x = startX; x < endX; x ++)
	{
	  cout << "z = "<< z << " x= " << x<<endl;
	  //	  fgTestCanvasPtr->cd(x*N_ZROWS_RCU (tmpHistoName, tmpHistoName, 300, 0, 299);+z +1);
	  //	  fgTestCanvasPtr->cd((z - startZ)*nCols +x +1);
	  fgTestCanvasPtr->cd((z - startZ)*nCols +x);
	  //	  fgTestCanvasPtr->cd( z*nCols    +x);

	  cout << "(z - startZ)*nCols +x  ="<< (z - startZ)*nCols +x <<endl;
	  
	  fgChannelDataPlotPtr[x][z]->Draw();
	} 
    }
  */

  
  /*  
  for(int z = 0; z < nRows; z ++)
    {
      for(int x = 0; x < nCols; x ++)
	{
	  //	  fgTestCanvasPtr->cd(x*N_ZROWS_RCU +z + 1);
	  fgTestCanvasPtr->cd(z*nCols +x + 1);
	  fgChannelDataPlotPtr[z][x]->Draw();
	} 
    }
  */



  
  //  fgTestCanvasPtr->Divide(N_XCOLUMNS_RCU, N_ZROWS_RCU, 0, 0);
  
  /*
    for(int z = 0; z < N_ZROWS_RCU; z ++)
    {
    for(int x = 0; x < N_XCOLUMNS_RCU; x ++)
    {
    //	  fgTestCanvasPtr->cd(x*N_ZROWS_RCU +z + 1);
	  
    fgTestCanvasPtr->cd(z*N_XCOLUMNS_RCU +x + 1);
	  
    //   fgTestCanvasPtr->cd(z*N_ZROWS_RCU +x + 1);
	  
    fgChannelDataPlotPtr[z][x]->Draw();
    } 
    }
  */


  
  fgTestCanvasPtr->Divide(N_XCOLUMNS_RCU, N_ZROWS_RCU, 0, 0);

  
  for(int z = 0; z < N_ZROWS_RCU; z ++)
    {
      for(int x = 0; x < N_XCOLUMNS_RCU; x ++)
	{
	  //	  fgTestCanvasPtr->cd(x*N_ZROWS_RCU +z + 1);
	  
	  fgTestCanvasPtr->cd(z*N_XCOLUMNS_RCU +x +1);
	  
	  //   fgTestCanvasPtr->cd(z*N_ZROWS_RCU +x + 1);
	  
	  fgChannelDataPlotPtr[z][x]->Draw();
	} 
    }
  
 

 

  fgTestCanvasPtr->Update();
}


