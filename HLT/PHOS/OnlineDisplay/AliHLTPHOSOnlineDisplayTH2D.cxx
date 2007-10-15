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

#define PX_MIN 85
#define PX_MAX 688
#define PZ_MIN 60
#define PZ_MAX 380



AliHLTPHOSOnlineDisplayTH2D::AliHLTPHOSOnlineDisplayTH2D()
{

}


AliHLTPHOSOnlineDisplayTH2D::AliHLTPHOSOnlineDisplayTH2D(AliHLTPHOSOnlineDisplay *onlineDisplayPtr, const char* name, 
							 const char* title, Int_t nbinsx, 
							 Double_t xlow, Double_t xup, Int_t nbinsy, 
							 Double_t ylow, Double_t yup) :  AliHLTPHOSBase(),
											 TH2D(name, title, nbinsx, 
											      xlow, xup, nbinsy, 
											      ylow, yup)
{
  fOnlineDisplayPtr = onlineDisplayPtr;

}


AliHLTPHOSOnlineDisplayTH2D::~AliHLTPHOSOnlineDisplayTH2D()
{
  
}

void
AliHLTPHOSOnlineDisplayTH2D::SetGain(int gain)
{
  fGain = gain;
}

void
AliHLTPHOSOnlineDisplayTH2D::ExecuteEvent(Int_t event, Int_t px, Int_t pz)
{
  char tmpName[256];
  int tmpZBin = 0;
  int tmpXBin = 0;

  //  cout <<"px = "<< px<<" pz ="  <<endl;

  if(event == 61)
    {
      fgRawDataCanvas = new TCanvas("TEST3", "PHOS HLT Raw Data Display", 1200, 1000); ;
      tmpZBin =  GetZBin(pz);
      tmpXBin =  GetXBin(px);
      if(tmpZBin > 54) {tmpZBin = 54;}
      if(tmpZBin < 1) {tmpZBin = 1;}
      sprintf(tmpName, "Z_%d  X_%d", tmpZBin, tmpXBin);
      fgRawDataCanvas->Divide(Z_ROWS, X_COLS); 
      int cnt = 0;
      
      for(int z= 1; z > -2; z--)
	{
	  for(int x=-1; x < X_COLS -1; x ++)
	    {
	      cnt ++;
	      sprintf(tmpName, "Z_%d  X_%d_gain%d", tmpZBin + z, tmpXBin + x, fGain);   
	      fgRawDataPlotsPtr[cnt] = new TH1D(tmpName, tmpName, fNTotalSamples, 0, fNTotalSamples -1);
	      fgRawDataCanvas->cd(cnt);	 
	      fgRawDataPlotsPtr[cnt]->SetFillColor(1);
	      fgRawDataPlotsPtr[cnt]->SetMaximum(1023); 
	      fgRawDataPlotsPtr[cnt]->Reset();   
	      fOnlineDisplayPtr->fgEventTabPtr->GetRawData( fgRawDataPlotsPtr[cnt], tmpXBin +x, tmpZBin +z, fGain); 
	      fgRawDataPlotsPtr[cnt]->Draw();
	    }
	}

      fgRawDataCanvas->Update();

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
