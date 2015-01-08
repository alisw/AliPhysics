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

/* $Id: AliTRDqaGuiESDs.cxx 23871 2008-02-12 11:48:20Z hristov $ */

#include "AliTRDqaGuiESDs.h"

#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRootEmbeddedCanvas.h"
#include "TGToolTip.h"

ClassImp(AliTRDqaGuiESDs)

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of ESD (Event Summary Data)
// It displays histograms created by 
// the AliTRDQADataMakerRec run during the reconstruction 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////


const Int_t AliTRDqaGuiESDs::fgkLogList[24] = {
  1,1,0,0,0,0,
  1,1,1,1,1,1,
  1,1,1,0,0,0,
  0,0,0,0,0,0
};

//////////////////////////////////////////////////////////////////////////////////
AliTRDqaGuiESDs::AliTRDqaGuiESDs() 
  :TGCompositeFrame()
  ,fPage(0)
 {
  //
  // Default constructor
  //

   for (Int_t i = 0; i < 6; i++) {
     fCanvasList[i] = 0x0;
     fHistList[i]   = 0x0;
     for (Int_t j = 0; j < 4; j++) {
       fNameList[j*i] = 0x0;
     }
   }

 }

//////////////////////////////////////////////////////////////////////////////////
AliTRDqaGuiESDs::AliTRDqaGuiESDs(TGWindow *parent, Int_t page) 
  :TGCompositeFrame(parent, 720, 400)
  ,fPage(page)
 {
  //
  // main constructor
  //
  
  SetLayoutManager(new TGMatrixLayout(this,2,3,1,1));

  fNameList[0] = "bits";
  fNameList[1] = "ptTRDr";
  fNameList[2] = "ptTRDrTPCo";
  fNameList[3] = "sector";
  fNameList[4] = "trdzTRDr";
  fNameList[5] = "trdzTRDrTPCo";

  fNameList[6] = "clsTRDo";
  fNameList[7] = "clsTRDr";
  fNameList[8] = "clsTRDz";
  fNameList[9] = "quality";
  fNameList[10] = "pidQuality";
  fNameList[11] = "chi2";
  
  fNameList[12] = "pid0";
  fNameList[13] = "pid2";
  fNameList[14] = "pid4";
  fNameList[15] = "tracksStack";
  fNameList[16] = "electronStack";
  fNameList[17] = "elRatioStack";
  
  fNameList[18] = "signalPzone_0";
  fNameList[19] = "signalPzone_1";
  fNameList[20] = "signalPzone_2";
  fNameList[21] = "signalPzone_3";
  fNameList[22] = "";
  fNameList[23] = "";

  for(Int_t i=0; i<6; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(fNameList[i+6*fPage], this, 320, 320);
    AddFrame(fCanvasList[i]);
    fCanvasList[i]->GetCanvas()->SetRightMargin(0.05);
    //TGToolTip *tip = new TGToolTip(this,fCanvasList[i], Form("Wal sie na ryja %d", i),1000);
    //tip->Show((i%3)*320, (i/3)*320);
  }
  
  for(Int_t i=0; i<6; i++) {
    fHistList[i] = 0;
  }
}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiESDs::SetQAFile(const char *filename) {
  //
  // sets a file with histograms
  //


  for(Int_t i=0; i<6; i++) {
    if (fHistList[i]) delete fHistList[i];
  }
  
  //  const char *opt[6] = {"colz", "", "", ""};

  TFile *file = new TFile(filename);
  file->cd("TRD/ESDs");
  
  for(int i=0; i<6; i++) {
    fHistList[i] = (TH1D*)gDirectory->Get(Form("qaTRD_esd_%s", fNameList[i+fPage*6]));
    fCanvasList[i]->GetCanvas()->cd();
 
    if (fPage == 3) {
      if (fHistList[i]) fHistList[i]->Draw("colz");
      gPad->SetLogz(1);
      gPad->SetLogx(1);
    } else {
      gPad->SetLogy(fgkLogList[i+6*fPage]);
      if (fHistList[i]) fHistList[i]->Draw();
    }

    fCanvasList[i]->GetCanvas()->Update();
  }
}

//////////////////////////////////////////////////////////////////////////////////
