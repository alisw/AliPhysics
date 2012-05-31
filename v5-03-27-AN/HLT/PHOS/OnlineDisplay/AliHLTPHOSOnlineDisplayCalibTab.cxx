// $Id$

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

#include  "AliHLTPHOSOnlineDisplayCalibTab.h"
#include <iostream>
#include "AliHLTPHOSRcuCellAccumulatedEnergyDataStruct.h"
#include "AliHLTPHOSGetEventButton.h"

#include "TRootEmbeddedCanvas.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TCanvas.h"
#include "AliHLTPHOSGetEventButton.h"

//#include "TStyle.h"

using namespace std;


AliHLTPHOSOnlineDisplayCalibTab::AliHLTPHOSOnlineDisplayCalibTab()
{
  // See header file for documentation
  cout << "AliHLTPHOSOnlineDisplayCalibTab:ERROR: You cannot create a onlinedisplay Tab without arguments" << endl;
}

AliHLTPHOSOnlineDisplayCalibTab::AliHLTPHOSOnlineDisplayCalibTab(TGTab  *tabPtr, HOMERReader *homerSyncPtr, HOMERReader *homerPtrs[MAXHOSTS], int nHosts)
{
  // See header file for documentation
  for(int i=0; i<MAXHOSTS; i++)
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
  // See header file for documentation

}

void
AliHLTPHOSOnlineDisplayCalibTab::ReadBlockData(HOMERReader *homerReaderPtr)
{
  // See header file for documentation
  //  gStyle->SetOptLogy();
  //  gStyle->SetOptStat(true);

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
      int tmpdcsx;
      int tmpdcsz;

      for(int x = 0; x < NXCOLUMNSRCU; x ++)
	for(int z = 0; z <NZROWSRCU; z ++)
	  {
	    {
	      for(int gain = 0; gain < NGAINS; gain ++)
		{
		  tmpx = moduleID*NXCOLUMNSMOD + (accCellEnergiesPtr->fRcuX)*NXCOLUMNSRCU + x;
		  tmpz = (accCellEnergiesPtr->fRcuZ)*NZROWSRCU +z;

		  tmpdcsx =  (accCellEnergiesPtr->fRcuX)*NXCOLUMNSRCU + x  ;
		  tmpdcsz =  (accCellEnergiesPtr->fRcuZ)*NZROWSRCU +z;  

		  //	  if(tmpx < 140 && (tmpz < 6 || tmpz > 50)  )
		  //	  if(tmpx < 132 && (tmpz  > 52)  )
		  //		  if(tmpx < 132 && (tmpz  < 4)  )    
		  //		    {
		  //		      cout << " tmpx = "<< tmpx <<  "   tmpz =" << tmpz << "   x = "<< x << "   z = "<<  z << "  
		  // Energy ="<< accCellEnergiesPtr->fAccumulatedEnergies[x][z][gain] <<endl;
		  //		    }

		  
		  fgHitsHistPtr[gain]->Fill(tmpx , tmpz +1, accCellEnergiesPtr->fHits[x][z][gain] );
		  fgCalibHistPtr[gain]->SetBinContent(tmpx +1, tmpz +1, accCellEnergiesPtr->fAccumulatedEnergies[x][z][gain] );
		  fgHitsHistPtr[gain]->SetBinContent(tmpx +1, tmpz +1, accCellEnergiesPtr->fHits[x][z][gain] );
		  fDeadCannelMapPtr[gain]->SetBinContent(tmpx +1, tmpz +1, accCellEnergiesPtr->fDeadChannelMap[x][z][gain]);
		  
		  
		  //		  fgDCSViewPtr[gain]->SetBinContent(tmpdcsz +1, tmpdcsx +1, accCellEnergiesPtr->fAccumulatedEnergies[x][z][gain]);
		  fgDCSViewPtr[gain]->SetBinContent(tmpdcsz +1, tmpdcsx +1, accCellEnergiesPtr->fAccumulatedEnergies[x][z][gain]);

		  if(fgHitsHistPtr[gain]->GetBinContent(tmpx + 1, tmpz +1) > 0)
		    {
		      fgAveragePtr[gain]->SetBinContent(tmpx + 1, tmpz +1, fgCalibHistPtr[gain]->GetBinContent(tmpx +1, tmpz +1)/fgHitsHistPtr[gain]->GetBinContent(tmpx +1, tmpz +1));
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
  // See header file for documentation
  ResetDisplay();
  DoGetNextEvent();
  UpdateDisplay();
  fgEvntCnt ++;
}

void
AliHLTPHOSOnlineDisplayCalibTab::ResetDisplay() const
{
  // See header file for documentation
}

void 
AliHLTPHOSOnlineDisplayCalibTab::InitDisplay(TGTab *tabPtr)
{
  // See header file for documentation
  //  gStyle->SetOptLogy();
  //  gStyle->SetOptStat(true); 

  char tmpHistoName[256]; 

  // fgLegoPlotHGPtr = new TH2D("a Homer","HLT: #pi^{0} 5 - 30Gev HG, High gain",  
  fgLegoPlotHGPtr = new TH2D("a Homer","HLT: #pi^{0} 5 - 30Gev HG, High gain",  
			     NXCOLUMNSMOD*NMODULES , 0, NXCOLUMNSMOD*NMODULES,  
                             NZROWSMOD,               0, NZROWSMOD);
  fgLegoPlotHGPtr->SetMaximum( MAXBINVALUE);
  fgLegoPlotHGPtr->Reset();
  //  fgLegoPlotHGPtr->GetXaxis()->SetRange(128, 128 + 64);


  fgLegoPlotLGPtr = new TH2D("b Homer","HLT: #pi^{0} 5 - 30Gev LG, Low gain",  
			     NXCOLUMNSMOD* NMODULES , 0, NXCOLUMNSMOD* NMODULES,  
			     NZROWSMOD,          0, NZROWSMOD);
  fgLegoPlotLGPtr->SetMaximum( MAXBINVALUE); 
  fgLegoPlotLGPtr->Reset();
  //  fgLegoPlotLGPtr->GetXaxis()->SetRange(128, 128 + 64);



  TGLayoutHints *fL1 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 15, 1);


  for(int gain = 0; gain< NGAINS; gain ++)
    {
      sprintf(tmpHistoName, "TAB a HLT gain %d", gain);
      fgCalibHistPtr[gain] = new TH2D(tmpHistoName, tmpHistoName,  
				      NXCOLUMNSMOD*NMODULES , 0, NXCOLUMNSMOD*NMODULES , 
				      NZROWSMOD,         0, NZROWSMOD);
      fgCalibHistPtr[gain]->Reset(); 
      fgCalibHistPtr[gain]->GetXaxis()->SetRange(XRANGESTART, XRANGEEND);
      
      
      sprintf(tmpHistoName, "TAB b Calibration Data HLT: #pi^{0} 5 - 30GeV gain %d", gain);
      
      fgHitsHistPtr[gain] = new TH2I(tmpHistoName, tmpHistoName,  
				    NXCOLUMNSMOD* NMODULES , 0, NXCOLUMNSMOD*NMODULES,  
				    NZROWSMOD,          0, NZROWSMOD);
      //      fgHitsHistPtr[gain]->SetMaximum( MAXBINVALUE); 
      fgHitsHistPtr[gain]->Reset();
      fgHitsHistPtr[gain]->GetXaxis()->SetRange(XRANGESTART, XRANGEEND);
     
      sprintf(tmpHistoName, "Average Energy gain %d", gain);
      fgAveragePtr[gain] = new TH2D(tmpHistoName,tmpHistoName,  
				    NXCOLUMNSMOD* NMODULES , 0, NXCOLUMNSMOD*NMODULES,  
				    NZROWSMOD,          0, NZROWSMOD);
      //    fgAveragePtr[gain]->SetMaximum( MAXBINVALUE); 
      fgAveragePtr[gain]->Reset();
      fgAveragePtr[gain]->GetXaxis()->SetRange(XRANGESTART, XRANGEEND);

      //sprintf(tmpHistoName, "Dead Channel Map gain%d", gain);
      sprintf(tmpHistoName,  "Dead Channel Map gain%d", gain);
      fDeadCannelMapPtr[gain] = new TH2D(tmpHistoName,tmpHistoName,  
				    NXCOLUMNSMOD* NMODULES , 0, NXCOLUMNSMOD*NMODULES,  
				    NZROWSMOD,          0, NZROWSMOD);
      //    fDeadCannelMapPtr[gain]->SetMaximum( MAXBINVALUE); 
 
      fDeadCannelMapPtr[gain]->Reset();
      fDeadCannelMapPtr[gain]->GetXaxis()->SetRange(XRANGESTART, XRANGEEND);


      
      sprintf(tmpHistoName, "DCS view gain %d", gain);
      fgDCSViewPtr[gain] = new TH2D(tmpHistoName, tmpHistoName, 
				    NZROWSMOD, 0, NZROWSMOD,  
				    NXCOLUMNSMOD, 0, NXCOLUMNSMOD);
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

	   tf2 = fSubTab2->AddTab("Dead Channel Map"); 
	   fSubF8 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
	   tf2->AddFrame(fSubF8, fL1);
	   fEc15 = new TRootEmbeddedCanvas("ec15", fSubF8, 100, 100);
	   fSubF8->AddFrame(fEc15, fL1);
	   fEc16 = new TRootEmbeddedCanvas("ec16", fSubF8, 100, 100);
	   fSubF8->AddFrame(fEc16, fL1);	   
	  
	   tf2 = fSubTab2->AddTab("SURF (DCS view)"); 
	   fSubF9 = new TGCompositeFrame(tf2, 60, 20, kVerticalFrame);
	   tf2->AddFrame(fSubF9, fL1);
	   
	   fEc17 = new TRootEmbeddedCanvas("ec17", fSubF9, 100, 100);
	   fSubF9->AddFrame(fEc17, fL1);
	   fEc18 = new TRootEmbeddedCanvas("ec18", fSubF9, 100, 100);
	   fSubF9->AddFrame(fEc18, fL1);
	   


	   fSubTab2->Resize();
	   tf->AddFrame(fSubTab2, fL1);

	   fgEventButtPtr = new  AliHLTPHOSGetEventButton(fSubF4, "update histograms z", 'h');
}


void 
AliHLTPHOSOnlineDisplayCalibTab::UpdateDisplay()
{
  // See header file for documentation
  // gStyle->SetOptLogy();
  //  gStyle->SetOptStat(true);
 
  fgCanvasHGPtr =  fEc7->GetCanvas();
  fgCanvasHGPtr->cd();

  fgCalibHistPtr[HIGHGAIN]->Draw("LEGO2Z");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr = fEc8->GetCanvas();
  fgCanvasLGPtr->cd();
  fgCalibHistPtr[LOWGAIN]->Draw("LEGO2Z");
  fgCanvasLGPtr->Update();
  fgCanvasHGPtr =  fEc9->GetCanvas();

  fgCanvasHGPtr->cd();
  fgHitsHistPtr[HIGHGAIN]->Draw("SCAT");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr = fEc10->GetCanvas();
  fgCanvasLGPtr->cd();
  fgHitsHistPtr[LOWGAIN]->Draw("SCAT");
  fgCanvasLGPtr->Update();
  fgCanvasHGPtr =  fEc11->GetCanvas();
  fgCanvasHGPtr->cd();

  fgCalibHistPtr[HIGHGAIN]->Draw("COLZ");
  fgCanvasHGPtr->Update();
  fgCanvasLGPtr = fEc12->GetCanvas();
  fgCanvasLGPtr->cd();
  fgCalibHistPtr[LOWGAIN]->Draw("COLZ");
  fgCanvasLGPtr->Update();


  fgCanvasLGPtr = fEc13->GetCanvas();
  fgCanvasLGPtr->cd();
  fgAveragePtr[HIGHGAIN]->Draw("COLZ");
  fgCanvasLGPtr->Update();

  fgCanvasHGPtr = fEc14->GetCanvas();
  fgCanvasHGPtr->cd();
  fgAveragePtr[LOWGAIN]->Draw("COLZ");
  fgCanvasHGPtr->Update();

  fgCanvasLGPtr = fEc15->GetCanvas();
  fgCanvasLGPtr->cd();
  fDeadCannelMapPtr[HIGHGAIN]->Draw("COL");
  fgCanvasLGPtr->Update();

  fgCanvasHGPtr = fEc16->GetCanvas();
  fgCanvasHGPtr->cd();
  fDeadCannelMapPtr[LOWGAIN]->Draw("COL");
  fgCanvasHGPtr->Update();


  fgCanvasLGPtr = fEc17->GetCanvas();
  fgCanvasLGPtr->cd();
  fgDCSViewPtr[HIGHGAIN]->Draw("COLZ");
  fgCanvasLGPtr->Update();

  fgCanvasHGPtr = fEc18->GetCanvas();
  fgCanvasHGPtr->cd();
  fgDCSViewPtr[LOWGAIN]->Draw("COLZ");
  fgCanvasHGPtr->Update();
  

}
