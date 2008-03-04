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

/* $Id: AliTRDqaGuiMain.cxx 23871 2008-02-12 11:48:20Z hristov $ */

#include "AliTRDqaGuiMain.h"
#include "qaGui/AliTRDqaGuiClusters.h"
#include "AliTRDqaGuiClustersSM.h"
#include "AliTRDqaGuiClustersStack.h"
#include "AliTRDqaGuiESDs.h"

#include "TGTab.h"

ClassImp(AliTRDqaGuiMain)

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

AliTRDqaGuiMain::AliTRDqaGuiMain() 
  : fGTabPanel(0),
    fGDet(0),
    fGSM(0),
    fGStack(0)
{
  //
  // Default constructor
  //

}

////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiMain::AliTRDqaGuiMain(TGWindow *parent) :
  TGMainFrame(parent, 400, 400),
  fGTabPanel(0),
  fGDet(0),
  fGSM(0),
  fGStack(0)
{
  //
  // main constructor
  //
  
  fGTabPanel = new TGTab(this);
  fGDet      = new AliTRDqaGuiClusters(fGTabPanel);
  fGSM       = new AliTRDqaGuiClustersSM(fGTabPanel);
  fGStack    = new AliTRDqaGuiClustersStack(fGTabPanel);
  fGESDs[0]  = new AliTRDqaGuiESDs(fGTabPanel,0);
  fGESDs[1]  = new AliTRDqaGuiESDs(fGTabPanel,1);
  fGESDs[2]  = new AliTRDqaGuiESDs(fGTabPanel,2);
  
  fGTabPanel->AddTab("Clusters", fGDet);
  fGTabPanel->AddTab("Clusters - Super Module", fGSM);
  fGTabPanel->AddTab("Clusters - Stack", fGStack);
  fGTabPanel->AddTab("ESDs (1)", fGESDs[0]);
  fGTabPanel->AddTab("ESDs (2)", fGESDs[1]);
  fGTabPanel->AddTab("ESDs (3)", fGESDs[2]);

  AddFrame(fGTabPanel);
  
  SetWindowName("TRD QA");
  MapSubwindows();
  MapWindow();
  Resize(GetDefaultSize());
}

////////////////////////////////////////////////////////////////////////////////
 
void AliTRDqaGuiMain::SetQAFile(const char *file) {
  //
  // sets a file with histograms
  // 
  
  fGDet->SetQAFile(file);
  fGSM->SetQAFile(file);
  fGStack->SetQAFile(file);
  fGESDs[0]->SetQAFile(file);
  fGESDs[1]->SetQAFile(file);
  fGESDs[2]->SetQAFile(file);
}

////////////////////////////////////////////////////////////////////////////////
