#include "TGLViewer.h"

namespace Alieve {
class ITSModuleStepper;
}

Alieve::ITSModuleStepper* stepper = 0;

void its_module_stepper(Int_t col = 4 , Int_t row = 3)
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadDigits("ITS");
  TTree* dt = rl->GetTreeD("ITS", false);
  Alieve::ITSDigitsInfo* di = new Alieve::ITSDigitsInfo();
  di->SetTree(dt);

  gStyle->SetPalette(1, 0);
  gReve->DisableRedraw();

  Alieve::ITSModuleStepper* store = new Alieve::ITSModuleStepper(di);
  store->SetStepper(col, row, 10, 10);
  store->SetMainColor((Color_t)2);
  gReve->AddRenderElement(store);
  stepper = store;
  
  TRandom r(0);
  Int_t module;
  for (Int_t i=0; i<40; ++i) {
    module = r.Integer(51);
    store->AddToList(module);
  }
  store->Start();
  gReve->EnableRedraw();

  TGLViewer * v = (TGLViewer *)gPad->GetViewer3D();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
}

void start()
{
  stepper->Start();
}

void next() 
{
  stepper->Next();
}
