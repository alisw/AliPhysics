// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGLViewer.h>
#include <TTree.h>
#include <TStyle.h>
#include <TEveManager.h>
#include <TEveTreeTools.h>

#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#include <AliEveITSModuleStepper.h>
#include <AliEveITSDigitsInfo.h>
#endif

class AliEveITSModuleStepper;

void its_module_stepper(Int_t det = 0)
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadDigits("ITS");
  TTree* dt = rl->GetTreeD("ITS", false);

  AliEveITSDigitsInfo* di = new AliEveITSDigitsInfo();
  di->SetTree(dt);
  di->Dump();

  gEve->DisableRedraw();
  AliEveITSModuleStepper* ms = new AliEveITSModuleStepper(di);
  ms->SetMainColor(8);
  gStyle->SetPalette(1, 0);
  ms->DisplayDet(det, -1);
  gEve->AddElement(ms);
  gEve->Redraw3D(kTRUE); // To enforce camera reset
  gEve->EnableRedraw();

  TGLViewer* v = gEve->GetDefaultGLViewer();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  
  /*
   * Disabling obsolete code
   * 
   */
  //TGLCameraMarkupStyle* mup = v->GetCameraMarkup();
  //if(mup) mup->SetShow(kFALSE);
}
