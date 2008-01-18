// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#include "TGLViewer.h"
namespace Alieve
{
class AliEveITSModuleStepper;
}

void its_module_stepper(Int_t det = 0)
{
  TFile *file = TFile::Open("ITS.Digits.root");
  TDirectory* dir = (TDirectory*) file->Get("Event0");
  TTree* tree =  (TTree*)dir->Get("TreeD");
  AliEveITSDigitsInfo* di = new AliEveITSDigitsInfo();
  di->SetTree(tree);

  gEve->DisableRedraw();
  AliEveITSModuleStepper* ms = new AliEveITSModuleStepper(di);
  ms->SetMainColor(Color_t(8));
  gStyle->SetPalette(1, 0);
  ms->DisplayDet(det, -1);
  gEve->AddElement(ms);
  gEve->Redraw3D(kTRUE); // To enforce camera reset
  gEve->EnableRedraw();

  TGLViewer* v = (TGLViewer *)gEve->GetGLViewer();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  TGLCameraMarkupStyle* mup = v->GetCameraMarkup();
  if(mup) mup->SetShow(kFALSE);
}
