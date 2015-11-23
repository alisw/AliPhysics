// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGLayout.h>
#include <TGFrame.h>
#include <TEveManager.h>
#include <TEveBrowser.h>

#include <AliTPCCalibViewerGUI.h>
#endif

void tpc_calib_viewer(const char* file="CalibTree.root")
{
   TEveBrowser* b = gEve->GetBrowser();
   b->StartEmbedding(1);

   TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
   frmMain->SetWindowName("AliTPCCalibViewer GUI");
   frmMain->SetCleanup(kDeepCleanup);

   AliTPCCalibViewerGUI* calibViewer1 = new AliTPCCalibViewerGUI(frmMain, 1000, 600, (char*)file);
   frmMain->AddFrame(calibViewer1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

   frmMain->MapSubwindows();
   frmMain->Resize();
   frmMain->MapWindow();

   b->StopEmbedding();
   b->SetTabTitle("TPC CalibViewer", 1);
}
