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

/* $Id: AliTRDqaGuiClusters.cxx 23871 2008-02-12 11:48:20Z hristov $ */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of clusters on the full detector level.
// It displays histograms created by the AliTRDQADataMakerRec 
// run during the reconstruction 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "AliTRDqaGuiClusters.h"

#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRootEmbeddedCanvas.h"

ClassImp(AliTRDqaGuiClusters)

const Int_t AliTRDqaGuiClusters::fgkLogList[4] = {0,1,0,1};

//////////////////////////////////////////////////////////////////////////////////
AliTRDqaGuiClusters::AliTRDqaGuiClusters(TGWindow *parent) 
  : TGCompositeFrame(parent, 720, 400) {
  //
  // main constructor
  //
  
  SetLayoutManager(new TGMatrixLayout(this,2,2,1,1));

  fNameList[0] = "detMap";
  fNameList[1] = "nCls";
  fNameList[2] = "signal";
  fNameList[3] = "ampMPV";

  for(Int_t i=0; i<4; i++) {
    fCanvasList[i] = new TRootEmbeddedCanvas(fNameList[i], this, 480, 320);
    AddFrame(fCanvasList[i]);
  }
  
  for(Int_t i=0; i<4; i++) {
    fHistList[i] = 0;
  }

}

//////////////////////////////////////////////////////////////////////////////////

void AliTRDqaGuiClusters::SetQAFile(const char *filename) {
  //
  // sets the file with histograms
  //

  for(Int_t i=0; i<4; i++) {
    if (fHistList[i]) delete fHistList[i];
  }
  
  const char *opt[4] = {"colz", "", "", ""};

  TFile *file = new TFile(filename);
  file->cd("TRD/RecPoints");
  
  for(int i=0; i<4; i++) {
    fHistList[i] = (TH1D*)gDirectory->Get(Form("qaTRD_recPoints_%s", fNameList[i]));
    fCanvasList[i]->GetCanvas()->cd();
    gPad->SetLogy(fgkLogList[i]);
    if (fHistList[i]) fHistList[i]->Draw(opt[i]);
    fCanvasList[i]->GetCanvas()->Update();
  }
}

//////////////////////////////////////////////////////////////////////////////////
