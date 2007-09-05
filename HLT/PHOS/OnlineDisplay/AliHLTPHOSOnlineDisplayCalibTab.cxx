#include  "AliHLTPHOSOnlineDisplayCalibTab.h"
#include <iostream>
#include "AliHLTPHOSRcuCellAccumulatedEnergyDataStruct.h"
#include "AliHLTPHOSGetEventButton.h"

using namespace std;

AliHLTPHOSOnlineDisplayCalibTab::AliHLTPHOSOnlineDisplayCalibTab()
{
  cout << "AliHLTPHOSOnlineDisplayCalibTab:ERROR: You cannot create a onlinedisplay Tab without arguments" << endl;
}

AliHLTPHOSOnlineDisplayCalibTab::AliHLTPHOSOnlineDisplayCalibTab(TGTab  *tabPtr, HOMERReader *homerSyncPtr, HOMERReader *homerPtrs[MAX_HOSTS], int nHosts)
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


AliHLTPHOSOnlineDisplayCalibTab::~AliHLTPHOSOnlineDisplayCalibTab()
{


}


void
AliHLTPHOSOnlineDisplayCalibTab::ReadBlockData(HOMERReader *homerReaderPtr)
{
  unsigned long blk = homerReaderPtr->FindBlockNdx("UCCARENE","SOHP", 0xFFFFFFFF );

  while ( blk != ~(unsigned long)0 ) 
    {
      cout << "AliHLTPHOSOnlineDisplayCalibTab::ReadBlockDat:GetHistogram: updating block " << endl;
      AliHLTUInt16_t moduleID;
      const AliHLTPHOSRcuCellAccumulatedEnergyDataStruct* accCellEnergiesPtr = (const AliHLTPHOSRcuCellAccumulatedEnergyDataStruct*)homerReaderPtr->GetBlockData( blk ); 
      moduleID = accCellEnergiesPtr->fModuleID ;
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

		  //////////////////////////////////////////////////////
		  //Added for debugging purposes, want a pure dead map//
		  //////////////////////////////////////////////////////
		  
		  if(accCellEnergiesPtr->fHits[x][z][gain] > 0)
		    {
		      fgHitsHistPtr[gain]->SetBinContent(x, z, 10);
		    }
		  
		  ////////////////////////////////////////////////////



		  //fgHitsHistPtr[gain]->Fill(tmpx, tmpz, accCellEnergiesPtr->fHits[x][z][gain] );
		  
		  if(fgHitsHistPtr[gain]->GetBinContent(tmpx, tmpz) > 0)
		    {
		      fgAveragePtr[gain]->SetBinContent(tmpx, tmpz, fgCalibHistPtr[gain]->GetBinContent(tmpx, tmpz)/fgHitsHistPtr[gain]->GetBinContent(tmpx, tmpz));
		    }
		}
	    }
	  }
     
      blk = homerReaderPtr->FindBlockNdx("UCCARENE","SOHP", 0xFFFFFFFF, blk+1);
    }
}


int 
AliHLTPHOSOnlineDisplayCalibTab::GetNextEvent()
{
  ResetDisplay();
  DoGetNextEvent();
  UpdateDisplay();
  fgEvntCnt ++;
}

void
AliHLTPHOSOnlineDisplayCalibTab::ResetDisplay()
{

}

void 
AliHLTPHOSOnlineDisplayCalibTab::InitDisplay(TGTab *tabPtr)
{
  char tmpHistoName[256]; 

  fgLegoPlotHGPtr = new TH2D("Homer","HLT: #pi^{0} 5 - 30Gev HG, High gain",  
			     N_XCOLUMNS_MOD*N_MODULES , 0, N_XCOLUMNS_MOD*N_MODULES,  
                             N_ZROWS_MOD,               0, N_ZROWS_MOD);
  fgLegoPlotHGPtr->SetMaximum( MAX_BIN_VALUE);
  fgLegoPlotHGPtr->Reset();

  fgLegoPlotLGPtr = new TH2D("Homer","HLT: #pi^{0} 5 - 30Gev LG, Low gain",  
			     N_XCOLUMNS_MOD* N_MODULES , 0, N_XCOLUMNS_MOD* N_MODULES,  
			     N_ZROWS_MOD,          0, N_ZROWS_MOD);
  fgLegoPlotLGPtr->SetMaximum( MAX_BIN_VALUE); 
  fgLegoPlotLGPtr->Reset();


 TGLayoutHints *fL1 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 15, 1);


  for(int gain = 0; gain< N_GAINS; gain ++)
    {
      sprintf(tmpHistoName, "TAB a HLT gain %d", gain);
      fgCalibHistPtr[gain] = new TH2D(tmpHistoName, tmpHistoName,  
				      N_XCOLUMNS_MOD*N_MODULES , 0, N_XCOLUMNS_MOD*N_MODULES , 
				      N_ZROWS_MOD,         0, N_ZROWS_MOD);
      fgCalibHistPtr[gain]->Reset(); 
      
      sprintf(tmpHistoName, "TAB b Calibration Data HLT: #pi^{0} 5 - 30GeV gain %d", gain);
      fgHitsHistPtr[gain] = new TH2I(tmpHistoName, tmpHistoName,  
				    N_XCOLUMNS_MOD* N_MODULES , 0, N_XCOLUMNS_MOD*N_MODULES,  
				    N_ZROWS_MOD,          0, N_ZROWS_MOD);
      fgHitsHistPtr[gain]->SetMaximum( MAX_BIN_VALUE); 
      fgHitsHistPtr[gain]->Reset();
      
      sprintf(tmpHistoName, "TAB c Average Data HLT: #pi^{0} 5 - 30GeV gain %d", gain);
      fgAveragePtr[gain] = new TH2D(tmpHistoName,tmpHistoName,  
				    N_XCOLUMNS_MOD* N_MODULES , 0, N_XCOLUMNS_MOD*N_MODULES,  
				    N_ZROWS_MOD,          0, N_ZROWS_MOD);
      fgAveragePtr[gain]->SetMaximum( MAX_BIN_VALUE); 
      fgAveragePtr[gain]->Reset();
    }


   TGCompositeFrame  *tf = tabPtr->AddTab("Calibration data zzz");
     

          fSubTab2 = new TGTab(tf, 100, 100);

	   TGCompositeFrame	   *tf2 = fSubTab2->AddTab("Accumulated energy");   
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

  fgEventButtPtr = new  AliHLTPHOSGetEventButton(fSubF4, "update histograms z", 'h');
}


void 
AliHLTPHOSOnlineDisplayCalibTab::UpdateDisplay()
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
