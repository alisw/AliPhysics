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

  for (Int_t i = 0; i < 6; i++) {
    fCanvasList[i] = 0x0;
  }
  for (Int_t i = 0; i < 3; i++) {
    fHistList[i]   = 0x0;
  }
  for (Int_t i = 0; i < 2; i++) {
    fGraphList[i]  = 0x0;
  }

  strncpy(fFileName,"",256);

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
    "noiseTotal", "peakPeak", "mcmEvDist_999",
    "trendMCM", "fracPP_0", "nADCinEvent"
  };  

  const char *title[6] = {
    ";noise (ADC)", 
    ";peak-peak (ADC)",
    ";#Delta Event Number from MCM",

    "number of active MCMs;event number",
    "number of ADCs with PP > 10;event number",
    "number of ADC chanels;event number"
  };


  strncpy(fFileName,filename,256);
 
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
      gPad->SetLogy();      
      fHistList[i]->Draw();
      fHistList[i]->SetTitle(title[i]);
    }
    
    
    fGraphList[i] = (TGraph*)file->Get(names[i+3]);
    if (fGraphList[i]) {
      fCanvasList[i+3]->GetCanvas()->cd();
      fGraphList[i]->Draw("apl");
      fGraphList[i]->SetMarkerStyle(7);
      fGraphList[i]->GetHistogram()->SetTitle(title[i+3]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
