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

/* $Id: AliTRDqaGuiMainBlack.cxx 23871 2008-02-12 11:48:20Z hristov $ */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of clusters. It lets display and browse throu histograms created by 
// the AliTRDQADataMakerRec run during the reconstruction 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "AliTRDqaGuiMainBlack.h"
#include "AliTRDqaGuiBlackSM.h"
#include "AliTRDqaGuiBlackChamber.h"

#include "TGTab.h"

ClassImp(AliTRDqaGuiMainBlack)

////////////////////////////////////////////////////////////////////////////////

AliTRDqaGuiMainBlack::AliTRDqaGuiMainBlack(TGWindow *parent) :
  TGMainFrame(parent, 400, 400) 
{
  //
  // Main constructor
  //
  
  fGTabPanel = new TGTab(this);
  fGSM       = new AliTRDqaGuiBlackSM(fGTabPanel);
  fGChamber  = new AliTRDqaGuiBlackChamber(fGTabPanel);
  
  fGChamber->SetRangePed(8, 11);
  fGChamber->SetRangeNoise(0.5, 2);

  fGSM->SetRangePed(8, 11);
  fGSM->SetRangeNoise(0.5, 3);

  fGTabPanel->AddTab("Super Module", fGSM);
  fGTabPanel->AddTab("Chamber", fGChamber);
   
  AddFrame(fGTabPanel);
  
  SetWindowName("TRD QA -- Black Events");
  MapSubwindows();
  MapWindow();
  Resize(GetDefaultSize());
}

////////////////////////////////////////////////////////////////////////////////
 
void AliTRDqaGuiMainBlack::SetQAFile(const char *file) 
{
  //
  // Set the QA file
  //
 
  fGChamber->SetQAFile(file);
  fGSM->SetQAFile(file);
}

////////////////////////////////////////////////////////////////////////////////
