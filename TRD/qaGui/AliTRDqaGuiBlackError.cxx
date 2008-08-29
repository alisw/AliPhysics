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

/* $Id: AliTRDqaGuiBlackError.cxx 23387 2008-01-17 17:25:16Z cblume $ */

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

#include "AliTRDqaGuiBlackError.h"

#include "TH1D.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TGaxis.h"

#include "TRootEmbeddedCanvas.h"

ClassImp(AliTRDqaGuiBlackError)

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiBlackError::AliTRDqaGuiBlackError() 
{
  //
  // Default constructor
  //

}

//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiBlackError::AliTRDqaGuiBlackError(TGWindow *parent) 
  : TGCompositeFrame(parent, 720, 500)
{
  //
  // Main constructor
  //
    
  SetLayoutManager(new TGMatrixLayout(this,3,3,0,0));

  for(Int_t i=0; i<9; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(Form("pos_%d", i), this, 333, 245);
    AddFrame(fCanvasList[i]);
  }
  
  for(Int_t i=0; i<3; i++) {
    fHistList[i] = 0;
    fGraphList[i] = 0;
    fHistListSM[i] = 0;
  }

  //AddFrame(fGCanvas);
  /**/
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiBlackError::SetQAFile(const char *filename) {
  //
  // Set the file with histograms
  //

  TGaxis::SetMaxDigits(3);
  
  const char *names[9] = {
    "errorHC", "errorMCM", "errorADC",
    "trendErrorHC", "trendErrorMCM", "trendErrorADC",
    "errorSM_HC", "errorSM_MCM", "errorSM_ADC"
  };  

  const char *title[9] = {
    "half-chamber;HC error ID", "MCM;MCM error ID", "ADC;ADC error ID",
    "fraction with HC error (%);event number",
    "fraction with MCM error (%);event number",
    "fraction with ADC error (%);event number",    
    ";sm id;number of errors",
    ";sm id;number of errors",
    ";sm id;number of errors"
  };


  strcpy(fFileName,filename);
 
  for(int i=0; i<3; i++) {

    if (fHistList[i]) delete fHistList[i];
    if (fGraphList[i]) delete fGraphList[i];
    if (fHistListSM[i]) delete fHistListSM[i];

    fHistList[i] = 0;
    fGraphList[i] = 0;
    fHistListSM[i] = 0;
  }
  
  TFile *file = new TFile(filename);
  
  for(Int_t i=0; i<3; i++) {

    fHistList[i] = (TH1D*)file->Get(names[i]);
    if (fHistList[i]) {
      fCanvasList[i]->GetCanvas()->cd();
      fCanvasList[i]->GetCanvas()->SetLogy(1);
      
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

    fHistListSM[i] = (TH1D*)file->Get(names[i+6]);
    if (fHistListSM[i]) {
      fCanvasList[i+6]->GetCanvas()->cd();
      //fCanvasList[i]->GetCanvas()->SetLogy(1);
      
      fHistListSM[i]->Draw();
      fHistListSM[i]->SetTitle(title[i+6]);
    }
  }

  
  /*


      fHistList[i]->SetMinimum(fRangePed[0]);
      fHistList[i]->SetMaximum(fRangePed[1]);
    }
    
    if ( fHistList[i] && (fIdxType == 1) && fSetRangeNoise) {
      fHistList[i]->SetMinimum(fRangeNoise[0]);
      fHistList[i]->SetMaximum(fRangeNoise[1]);
    }

    fCanvasList[pos]->GetCanvas()->Update();
  }
  */
}

//////////////////////////////////////////////////////////////////////////////////
