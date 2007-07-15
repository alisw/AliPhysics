#include "AliHLTPHOSOnlineDisplayEventTab.h"
#include <iostream>
#include "TGFrame.h"
#include "AliHLTPHOSGetEventButton.h"
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTDataTypes.h"
#include "HOMERData.h"
#include "HOMERReader.h"
#include "HOMERWriter.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"


using namespace std;


AliHLTPHOSOnlineDisplayEventTab::AliHLTPHOSOnlineDisplayEventTab()
{
  cout << "ERROR: You cannot create a onlinedisplay Tab without arguments" << endl;
}


AliHLTPHOSOnlineDisplayEventTab::AliHLTPHOSOnlineDisplayEventTab(TGTab  *tabPtr, HOMERReader *homerSyncPtr, HOMERReader *homerPtrs[MAX_HOSTS], int nHosts)
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


int
AliHLTPHOSOnlineDisplayEventTab::GetNextEvent()
{
  ResetDisplay();
  DoGetNextEvent();
  UpdateDisplay();
  fgEvntCnt ++;
}



void 
AliHLTPHOSOnlineDisplayEventTab::ReadBlockData(HOMERReader *homeReaderPtr)
{
  cout << "Reading block data" << endl;

  unsigned long blk = homeReaderPtr->FindBlockNdx("RENELLEC","SOHP", 0xFFFFFFFF );

  while ( blk != ~(unsigned long)0 ) 
    {
      Int_t moduleID;
      const AliHLTPHOSRcuCellEnergyDataStruct* cellEnergiesPtr = (const AliHLTPHOSRcuCellEnergyDataStruct*)homeReaderPtr->GetBlockData( blk );  
      moduleID = cellEnergiesPtr->fModuleID ;

      cout <<"AliHLTPHOSOnlineDisplayEventTab::ReadBlockData,  fModuleID =" <<moduleID << endl; 

      Int_t tmpCount = cellEnergiesPtr->fCnt;
      Int_t tmpZ;
      Int_t tmpX;
      Int_t tmpGain;
	  
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


  fgLegoPlotHGPtr = new TH2D("Homer a eventTAB","xx HLT: #pi^{0} 5 - 30Gev HG, High gain",  
			     N_XCOLUMNS_MOD*N_MODULES , 0, N_XCOLUMNS_MOD*N_MODULES,  
                             N_ZROWS_MOD,               0, N_ZROWS_MOD);
  fgLegoPlotHGPtr->SetMaximum( MAX_BIN_VALUE);
  fgLegoPlotHGPtr->Reset();

  fgLegoPlotLGPtr = new TH2D("Homer b eventTab","x HLT: #pi^{0} 5 - 30Gev LG, Low gain",  
			     N_XCOLUMNS_MOD* N_MODULES , 0, N_XCOLUMNS_MOD* N_MODULES,  
			     N_ZROWS_MOD,          0, N_ZROWS_MOD);
  fgLegoPlotLGPtr->SetMaximum( MAX_BIN_VALUE); 
  fgLegoPlotLGPtr->Reset();

  TGLayoutHints *fL1 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 15, 1);

  TGCompositeFrame *tf = tabPtr->AddTab("Event display TAB");
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
