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

/* $Id$ */

#include "AliTRDqaGuiMainAnalysis.h"
#include "AliTRDqaGuiJPsi.h"
#include "AliTRDqaGuiEnergyDeposit.h"
#include "TGTab.h"

ClassImp(AliTRDqaGuiMainAnalysis)

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of clusters and ESDs. 
// It displays histograms created by the AliTRDQADataMakerRec 
// run during the reconstruction 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiMainAnalysis::AliTRDqaGuiMainAnalysis() 
  : fGTabPanel(0),
    fGJPsi(0),
    fGED(0)
{
  //
  // Default constructor
  //

}

////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiMainAnalysis::AliTRDqaGuiMainAnalysis(TGWindow *parent) :
  TGMainFrame(parent, 400, 400),
  fGTabPanel(0),
  fGJPsi(0),
  fGED(0)
{
  //
  // main constructor
  //
  
  fGTabPanel = new TGTab(this);
  fGJPsi     = new AliTRDqaGuiJPsi(fGTabPanel);
  fGED       = new AliTRDqaGuiEnergyDeposit(fGTabPanel);
   
  fGTabPanel->AddTab("JPsi", fGJPsi);
  fGTabPanel->AddTab("Energy Deposit", fGED);

  AddFrame(fGTabPanel);
  
  SetWindowName("TRD Analysis QA");
  MapSubwindows();
  MapWindow();
  Resize(GetDefaultSize());
}

////////////////////////////////////////////////////////////////////////////////
 
void AliTRDqaGuiMainAnalysis::SetQAFile() {
  //
  // sets a file with histograms
  // 
  
  fGJPsi->SetQAFile("outJPsi.root");
  fGED->SetQAFile("outEnergyDeposit.root");
}

////////////////////////////////////////////////////////////////////////////////
