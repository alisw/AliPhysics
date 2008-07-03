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
#include "AliHLTPHOSOnlineDisplayTH2D.h"
#include "AliHLTPHOSOnlineDisplay.h"
//#include "AliHLTPHOSBase.h"

#include <iostream>

using namespace std;

// #define PX_MIN 85
// #define PX_MAX 688
// #define PZ_MIN 60
// #define PZ_MAX 380

#define PX_MIN 55
#define PX_MAX 420
#define PZ_MIN 34
#define PZ_MAX 251


AliHLTPHOSOnlineDisplayTH2D::AliHLTPHOSOnlineDisplayTH2D()
{

}


AliHLTPHOSOnlineDisplayTH2D::AliHLTPHOSOnlineDisplayTH2D(AliHLTPHOSOnlineDisplay *onlineDisplayPtr, const char* name, 
							 const char* title, Int_t nbinsx, 
							 Double_t xlow, Double_t xup, Int_t nbinsy, 
							 Double_t ylow, Double_t yup) : // AliHLTPHOSBase(),
											 TH2D(name, title, nbinsx, 
											      xlow, xup, nbinsy, 
											      ylow, yup)


{
  fOnlineDisplayPtr = onlineDisplayPtr;

}


AliHLTPHOSOnlineDisplayTH2D::~AliHLTPHOSOnlineDisplayTH2D()
{
  
}

/*
void
AliHLTPHOSOnlineDisplayTH2D::SetGain(int gain)
{
  fGain = gain;
}
*/

void
AliHLTPHOSOnlineDisplayTH2D::ExecuteEvent(Int_t event, Int_t px, Int_t pz)
{
  //  char tmpName[256];
  int tmpZBin = 0;
  int tmpXBin = 0;


  // cout <<"px = "<< px<<" pz =" <<  pz <<endl;


  for(int gain = 0; gain < N_GAINS; gain++)
    {
      if(event == 61)
	{
	  char gainLabel[100];
	  char label[256];
	  fOnlineDisplayPtr->Gain2Text(gain,gainLabel);
	  sprintf(label, "%s_PHOS_HLT_Online_rawdatadisplay",gainLabel);

	  //	  fgRawDataCanvas = new TCanvas("TEST3", "PHOS HLT Raw Data Display", 1200, 1000); ;
	  fgRawDataCanvasPtr[gain] = new TCanvas(label, label, 1200, 1000); ;

	  tmpZBin =  GetZBin(pz);
	  tmpXBin =  GetXBin(px);
	  if(tmpZBin > 54) {tmpZBin = 54;}
	  if(tmpZBin < 1) {tmpZBin = 1;}
	  
	  //sprintf(tmpName, "Z_%d_X_%d", tmpZBin, tmpXBin);
	
	  fgRawDataCanvasPtr[gain]->Divide(Z_ROWS, X_COLS); 
	  int cnt = 0;
    
	  //	  sprintf(label, "Z_%d_X_%d_%s",  tmpZBin, tmpXBin, gainLabel);
  
	  for(int z= 1; z > -2; z--)
	    {
	      for(int x=-1; x < X_COLS -1; x ++)
		{
		  fOnlineDisplayPtr->Gain2Text(gain,gainLabel);
		  sprintf(label, "Z_%d_X_%d_%s",  tmpZBin + z, tmpXBin + x, gainLabel);

		  cnt ++;
		  //  sprintf(tmpName, "Z_%d_X_%d_gain%d", tmpZBin + z, tmpXBin + x, gain);   

		  fgRawDataPlotsPtr[cnt][gain] = new TH1D(label, label, ALTRO_MAX_SAMPLES, 0, ALTRO_MAX_SAMPLES -1);
		  fgRawDataCanvasPtr[gain]->cd(cnt);	 
		  fgRawDataPlotsPtr[cnt][gain]->SetFillColor(1);
		  //fgRawDataPlotsPtr[cnt]->SetMaximum(1023); 
		  fgRawDataPlotsPtr[cnt][gain]->Reset();   
		  Int_t nSamples = fOnlineDisplayPtr->fgEventTabPtr->GetRawData( fgRawDataPlotsPtr[cnt][gain], tmpXBin +x, tmpZBin +z, gain); 
		  
		  //		  cout << "nAliHLTPHOSOnlineDisplayTH2D::ExecuteEvent Samples = " << nSamples  <<endl;
		  
		  if(nSamples == 0)
		    {
		      //      fgRawDataPlotsPtr[cnt][gain]->();
		    }




		  //	  fgRawDataPlotsPtr[cnt][gain]->GetXaxis()->SetRangeUser(0, nSamples - 1);
		  fgRawDataPlotsPtr[cnt][gain]->GetXaxis()->SetRangeUser(0, nSamples);

		  fgRawDataPlotsPtr[cnt][gain]->GetYaxis()->SetRangeUser(fgRawDataPlotsPtr[cnt][gain]->GetMinimum() - fgRawDataPlotsPtr[cnt][gain]->GetMinimum()*0.1, fgRawDataPlotsPtr[cnt][gain]->GetMaximum() + fgRawDataPlotsPtr[cnt][gain]->GetMaximum()*0.1 + 5); 
		  fgRawDataPlotsPtr[cnt][gain]->SetXTitle("Time/samples");
		  fgRawDataPlotsPtr[cnt][gain]->SetYTitle("Amplitude/ADC counts");

		  //	      fgRawDataPlotsPtr[cnt]->GetYaxis()->SetRangeUser(-30, fgRawDataPlotsPtr[cnt]->GetMaximum() + fgRawDataPlotsPtr[cnt]->GetMaximum()*0.1 + 5); 
		  fgRawDataPlotsPtr[cnt][gain]->Draw();
		}
	    }


	  sprintf(label, "%s_PHOS_HLT__Online_rawdatadisplay",gainLabel);
	  fgRawDataCanvasSinglePtr[gain] = new TCanvas(label, label, 1200, 1000); ;  
	  fgRawDataCanvasSinglePtr[gain]->cd();
	  sprintf(label, "Z_%d_X_%d__%s",  tmpZBin, tmpXBin, gainLabel);
	  fgRawDataPlotsSinglePtr[gain] = new TH1D(label, label, ALTRO_MAX_SAMPLES, 0, ALTRO_MAX_SAMPLES -1);
	  fgRawDataPlotsSinglePtr[gain]->SetFillColor(1);
	  fgRawDataPlotsSinglePtr[gain]->Reset();   
	  Int_t nSamples = fOnlineDisplayPtr->fgEventTabPtr->GetRawData(fgRawDataPlotsSinglePtr[gain], tmpXBin, tmpZBin, gain); 
	  fgRawDataPlotsSinglePtr[gain]->GetXaxis()->SetRangeUser(0, nSamples - 1);
	  fgRawDataPlotsSinglePtr[gain]->GetYaxis()->SetRangeUser(fgRawDataPlotsSinglePtr[gain]->GetMinimum() - fgRawDataPlotsSinglePtr[gain]->GetMinimum()*0.1, fgRawDataPlotsSinglePtr[gain]->GetMaximum() + fgRawDataPlotsSinglePtr[gain]->GetMaximum()*0.1 + 5); 
	  //SetXTitle(const char* title)
	  fgRawDataPlotsSinglePtr[gain]->SetXTitle("Time/samples");
	  fgRawDataPlotsSinglePtr[gain]->SetYTitle("Amplitude/ADC counts");

	  fgRawDataPlotsSinglePtr[gain]->Draw();

		  //	      fgRawDataPlotsPtr[cnt]->GetYaxis()->SetRangeUser(-30, fgRawDataPlotsPtr[cnt]->GetMaximum() + fgRawDataPlotsPtr[cnt]->GetMaximum()*0.1 + 5); 

	  fgRawDataCanvasPtr[gain]->Update();
	  
	}
    }
}

int
AliHLTPHOSOnlineDisplayTH2D::GetXBin(Int_t px)
{
  float tmpBinRange = GetXaxis()->GetLast() -  GetXaxis()->GetFirst();
  float tmpPixelRange =   PX_MAX  -   PX_MIN;
  float tmpPixRelative = (px - PX_MIN)/tmpPixelRange ;
  int xBin = ((int)((tmpPixRelative)*tmpBinRange) +  (float)GetXaxis()->GetFirst());
  return xBin;
}

int
AliHLTPHOSOnlineDisplayTH2D::GetZBin(Int_t pz)
{
  float tmpBinRange = GetYaxis()->GetLast() -  GetYaxis()->GetFirst();
  float tmpPixelRange =   PZ_MAX  -   PZ_MIN;
  float tmpPixRelative = (pz - PZ_MIN)/tmpPixelRange;
  int zBin = (int)((1-tmpPixRelative)*tmpBinRange);
  return zBin;
}
