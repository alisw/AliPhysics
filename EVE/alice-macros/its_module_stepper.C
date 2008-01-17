#include "TGLViewer.h"
namespace Alieve
{
class ITSModuleStepper;
}

void its_module_stepper(Int_t det = 0)
{
  TFile *file = TFile::Open("ITS.Digits.root");
  TDirectory* dir = (TDirectory*) file->Get("Event0");
  TTree* tree =  (TTree*)dir->Get("TreeD");
  Alieve::ITSDigitsInfo* di = new Alieve::ITSDigitsInfo();
  di->SetTree(tree);

  gEve->DisableRedraw();
  Alieve::ITSModuleStepper* ms = new Alieve::ITSModuleStepper(di);
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
