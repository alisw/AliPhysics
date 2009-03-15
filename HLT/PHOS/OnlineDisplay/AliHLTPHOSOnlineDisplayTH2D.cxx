// $Id$

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

#define PXMIN 55
#define PXMAX 420
#define PZMIN 34
#define PZMAX 251


AliHLTPHOSOnlineDisplayTH2D::AliHLTPHOSOnlineDisplayTH2D()
{

}


AliHLTPHOSOnlineDisplayTH2D::AliHLTPHOSOnlineDisplayTH2D(AliHLTPHOSOnlineDisplay *onlineDisplayPtr, const char* name, 
							 const char* title, Int_t nbinsx, 
							 Double_t xlow, Double_t xup, Int_t nbinsy, 
							 Double_t ylow, Double_t yup) : TH2D(name, title, nbinsx, 
											     xlow, xup, nbinsy, 
											     ylow, yup),
											fRunNumber(0) 
											//		fIsSetRunNumber(false)
{
  fOnlineDisplayPtr = onlineDisplayPtr;
}


AliHLTPHOSOnlineDisplayTH2D::~AliHLTPHOSOnlineDisplayTH2D()
{
  
}


/*
  void
  AliHLTPHOS2DHistogram::ExecuteEvent(Int_t event, Int_t px, Int_t pz)
  {
  if(event == 61) // CRAP PTH, use root enumerator instead 
  {
  int zbin;
  int xbin;
  EvaluateBinPosition( GetObjectInfo(px, pz),  &zbin, &xbin);
  printf("\nThe bin that was double clicked was z = %d, x = %d\n", zbin, xbin);
 
  HandleDoubleClick(zbin, xbin);
  }
}
*/

void
AliHLTPHOSOnlineDisplayTH2D::EvaluateBinPosition(const char *info, int *z, int *x)
{
  float zpos;
  float xpos;
  float bincnt;
  sscanf(info,"(x=%f, y=%f, binx=%d, biny=%d, binc=%2.4f)\n", &zpos, &xpos, z, x, &bincnt);
}


void
AliHLTPHOSOnlineDisplayTH2D::ExecuteEvent(Int_t event, Int_t px, Int_t pz)
{
  int tmpZBin = 0;
  int tmpXBin = 0;

  for(int gain = 0; gain < NGAINS; gain++)
    {
      if(event == 61)
	{
	  char gainLabel[100];
	  char label[256];
	  fOnlineDisplayPtr->Gain2Text(gain,gainLabel);
	  sprintf(label, "%s_PHOS_HLT_Online_rawdatadisplay",gainLabel);
	  fgRawDataCanvasPtr[gain] = new TCanvas(label, label, 1200, 1000); ;
	  
	  
	  tmpZBin =  GetZBin(pz);
	  tmpXBin =  GetXBin(px);
	  
	  if(tmpZBin > 54) {tmpZBin = 54;}
	  if(tmpZBin < 1) {tmpZBin = 1;}
	  
	  // int tmpZBin = 0;
	  // int tmpXBin = 0;

	  /*
	  EvaluateBinPosition( GetObjectInfo(px, pz),  &tmpZBin, &tmpXBin);
	  */
	  cout << "px = " << px << "pz =" << pz << endl;
	  cout << __FILE__ << ":" <<__LINE__ << " ObjectInfo IS !!" <<  GetObjectInfo(px, pz)  << endl;

	  fgRawDataCanvasPtr[gain]->Divide(ZROWS, XCOLS); 
	  int cnt = 0;
    
	  for(int z= 1; z > -2; z--)
	    {
	      for(int x=-1; x < XCOLS -1; x ++)
		{
		  fOnlineDisplayPtr->Gain2Text(gain,gainLabel);

		  if( fRunNumber >= 0)
		    {
		      sprintf(label, "(z = %d, x = %d) %s , run %d",  tmpZBin + z, (tmpXBin + x)%64 , gainLabel, fRunNumber);

		    }
		  else
		    {
		      sprintf(label, "(z = %d, x = %d) %s, unknow run number",  tmpZBin + z, (tmpXBin + x)%64,  gainLabel); 
		    }


		  cnt ++;
		  fgRawDataPlotsPtr[cnt][gain] = new TH1D(label, label, ALTROMAXSAMPLES, 0, ALTROMAXSAMPLES -1);
		  fgRawDataCanvasPtr[gain]->cd(cnt);	 
		  fgRawDataPlotsPtr[cnt][gain]->SetFillColor(1);
		  fgRawDataPlotsPtr[cnt][gain]->Reset();   
		  Int_t nSamples = fOnlineDisplayPtr->fgEventTabPtr->GetRawData( fgRawDataPlotsPtr[cnt][gain], tmpXBin +x, tmpZBin +z, gain); 
		  fgRawDataPlotsPtr[cnt][gain]->GetXaxis()->SetRangeUser(0, nSamples);
		  fgRawDataPlotsPtr[cnt][gain]->GetYaxis()->SetRangeUser(fgRawDataPlotsPtr[cnt][gain]->GetMinimum() - fgRawDataPlotsPtr[cnt][gain]->GetMinimum()*0.1, fgRawDataPlotsPtr[cnt][gain]->GetMaximum() + fgRawDataPlotsPtr[cnt][gain]->GetMaximum()*0.1 + 5); 
		  fgRawDataPlotsPtr[cnt][gain]->SetXTitle("Time / samples");
		  fgRawDataPlotsPtr[cnt][gain]->SetYTitle("Amplitude / ADC counts");
		  fgRawDataPlotsPtr[cnt][gain]->Draw();
		}
	    }

	  sprintf(label, "%s_PHOS_HLT__Online_rawdatadisplay",gainLabel);
	  fgRawDataCanvasSinglePtr[gain] = new TCanvas(label, label, 1200, 1000); ;  
	  fgRawDataCanvasSinglePtr[gain]->cd();
	  //	  sprintf(label, "Z_%d_X_%d__%s",  tmpZBin, tmpXBin, gainLabel);

	  
	  if( fRunNumber >= 0)
	    {
	      sprintf(label, "(z = %d, x = %d) %s , run %d",  tmpZBin, tmpXBin%64, gainLabel, fRunNumber);
	    }
	  else
	    {
	      sprintf(label, "(z = %d, x = %d) %s, unknow run number",  tmpZBin, tmpXBin%64, gainLabel); 
	    }



	  fgRawDataPlotsSinglePtr[gain] = new TH1D(label, label, ALTROMAXSAMPLES, 0, ALTROMAXSAMPLES -1);
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
  float tmpPixelRange =   PXMAX  -   PXMIN;
  float tmpPixRelative = (px - PXMIN)/tmpPixelRange ;
  int xBin = ((int)((tmpPixRelative)*tmpBinRange) +  (float)GetXaxis()->GetFirst());
  return xBin;
}

int
AliHLTPHOSOnlineDisplayTH2D::GetZBin(Int_t pz)
{
  float tmpBinRange = GetYaxis()->GetLast() -  GetYaxis()->GetFirst();
  float tmpPixelRange =   PZMAX  -   PZMIN;
  float tmpPixRelative = (pz - PZMIN)/tmpPixelRange;
  int zBin = (int)((1-tmpPixRelative)*tmpBinRange);
  return zBin;
}
