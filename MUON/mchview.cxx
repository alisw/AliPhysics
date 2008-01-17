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

// $Id$

/// \ingroup graphics
/// \file mchview.cxx
/// \brief Tracker visualization program
///
/// \author Laurent Aphecetche, Subatech


#include "AliMUONPainterDataSourceFrame.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONPainterMasterFrame.h"
#include "AliMUONPainterRegistry.h"
#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TCanvas.h>
#include <TEnv.h>
#include <TGMenu.h>
#include <TGTab.h>
#include <TROOT.h>
#include <TRint.h>
#include <TString.h>
#include <TStyle.h>

//_____________________________________________________________________________
void CreateMenuBar(TRint* app, TGMainFrame* mainFrame, UInt_t w)
{
/// 

  TGPopupMenu* file = new TGPopupMenu(gClient->GetRoot());
  
  file->AddEntry("&Exit",1);
  
  file->Connect("Activated(Int_t)","TRint",app,"Terminate()");

  TGMenuBar* bar = new TGMenuBar(mainFrame,w);
  
  bar->AddPopup("&File",file,new TGLayoutHints(kLHintsLeft|kLHintsTop));
  
  mainFrame->AddFrame(bar,new TGLayoutHints(kLHintsLeft|kLHintsExpandX));
  
  AliMUONPainterRegistry::Instance()->SetMenuBar(bar);
}


int main(int argc, char** argv)
{
///

  AliWarningGeneral("main","Remove default storage and run number from here...");
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetRun(0);
  
  TRint *theApp = new TRint("mchview", &argc, argv);

  gROOT->SetStyle("Plain");
  
  gStyle->SetPalette(1);
  
  Int_t n = gStyle->GetNumberOfColors();
  
  Int_t* colors = new Int_t[n+2];
  
  for ( Int_t i = 1; i <= n; ++i )
  {
    colors[i] = gStyle->GetColorPalette(i-1);
  }
  
  colors[0] = 0;
  colors[n+1] = 1;
  
  gStyle->SetPalette(n+2,colors);
  
  delete[] colors;
  
  UInt_t dw = gClient->GetDisplayWidth(); 
  UInt_t dh = gClient->GetDisplayHeight(); 
  
  UInt_t w = (UInt_t)(0.7*dw);
  UInt_t h = (UInt_t)(0.90*dh);
    
  TGMainFrame* mainFrame = new TGMainFrame(gClient->GetRoot(),w,h);
  
  const Int_t bs = 2;
  
  CreateMenuBar(theApp,mainFrame,w);

//  h -= 60; // menubar
  
  TGTab* tabs = new TGTab(mainFrame,w,h);
  
  TGCompositeFrame* t = tabs->AddTab("Painter Master Frame");

  AliMUONPainterMasterFrame* pf = 
    new AliMUONPainterMasterFrame(t,t->GetWidth()-bs*2,t->GetHeight()-bs*2);
  
  t->AddFrame(pf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,bs,bs,bs,bs));

  t = tabs->AddTab("Data Sources");
  
  AliMUONPainterDataSourceFrame* dsf = 
    new AliMUONPainterDataSourceFrame(t,t->GetWidth()-bs*2,t->GetHeight()-bs*2);
  
  t->AddFrame(dsf,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,bs,bs,bs,bs));
  
  mainFrame->AddFrame(tabs,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,0,0,0,0));

  mainFrame->SetWindowName("mchview - Visualization of MUON Tracker detector");

  mainFrame->MapSubwindows();
  mainFrame->Resize();
  mainFrame->MapWindow();
  
  mainFrame->Connect("CloseWindow()","TRint",theApp,"Terminate()");
  
  UInt_t x = dw/2 - w/2;
  UInt_t y = 0;
  
  mainFrame->MoveResize(x, y, w, h); 
  mainFrame->SetWMPosition(x, y);
  
  mainFrame->SetWMSizeHints(w,h,w,h,0,0);

  AliCodeTimer::Instance()->Print();

  // --- Start the event loop ---
  theApp->Run(kTRUE);

  AliMUONPainterHelper::Instance()->Save();
}
