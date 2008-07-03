/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliTRDqaGuiBlackGlobal.cxx 23387 2008-01-17 17:25:16Z cblume $ */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of black (non zero zuppresed) events from TRD. 
// It lets display and browse throu histograms created by the class 
// AliTRDqaBlackEvents.
// The class works in cooperation with AliTRDqaGuiMainBlack.
//
// S. Radomski 
// Uni-Heidelberg
// June 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "AliTRDqaGuiBlackGlobal.h"

#include "TH1D.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TGaxis.h"

#include "TRootEmbeddedCanvas.h"

ClassImp(AliTRDqaGuiBlackGlobal)

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiBlackGlobal::AliTRDqaGuiBlackGlobal() 
{
  //
  // Default constructor
  //

}

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiBlackGlobal::AliTRDqaGuiBlackGlobal(TGWindow *parent) 
  : TGCompositeFrame(parent, 720, 500)
{
  //
  // Main constructor
  //
    
  SetLayoutManager(new TGMatrixLayout(this,2,3,0,0));

  for(Int_t i=0; i<6; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(Form("pos_%d", i), this, 330, 350);
    AddFrame(fCanvasList[i]);
  }
  
  for(Int_t i=0; i<3; i++) {
    fHistList[i] = 0;
    fGraphList[i] = 0;
  }
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiBlackGlobal::SetQAFile(const char *filename) {
  //
  // Set the file with histograms
  //

  TGaxis::SetMaxDigits(3);
  
  const char *names[6] = {
    "smHcPed", "noiseTotal", "peakPeak",
    "trendMCM", "nADCinEvent", ""
  };  

  const char *title[6] = {
    "active half-chambers;super-module;half chamber",
    ";noise (ADC)", 
    ";peak-peak (ADC)",
    ";event number;number of strange MCMs",
    ";event number;number of fired ADC channels",
    ""
  };


  strcpy(fFileName,filename);
 
  for(int i=0; i<3; i++) {

    if (fHistList[i]) delete fHistList[i];
    if (fGraphList[i]) delete fGraphList[i];

    fHistList[i] = 0;
    fGraphList[i] = 0;
  }
  
  TFile *file = new TFile(filename);
  
  for(Int_t i=0; i<3; i++) {

    fHistList[i] = (TH1D*)file->Get(names[i]);
    if (fHistList[i]) {
      fCanvasList[i]->GetCanvas()->cd();
      if (i==0) fHistList[i]->Draw("colz");
      if (i>0) {
	gPad->SetLogy();
	fHistList[i]->Draw();
      }
      fHistList[i]->SetTitle(title[i]);
      
      // draw lines
      if (i==0) {
	for(Int_t stx=1; stx < 5; stx++) {
	  TLine *line = new TLine(-0.5, stx*12-0.5, 17.5, stx*12-0.5);
	  line->Draw(); 
	}
      }
    }
    
    fGraphList[i] = (TGraph*)file->Get(names[i+3]);
    if (fGraphList[i]) {
      fCanvasList[i+3]->GetCanvas()->cd();
      fGraphList[i]->Draw("apl");
      fGraphList[i]->GetHistogram()->SetTitle(title[i+3]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
