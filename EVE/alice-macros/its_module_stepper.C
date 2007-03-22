#include "TGLViewer.h"

namespace Alieve
{
class ITSModuleStepper;
}

Alieve::ITSModuleStepper* stepper = 0;

void its_module_stepper(Int_t col = 4 , Int_t row = 4)
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadDigits("ITS");
  TTree* dt = rl->GetTreeD("ITS", false);
  Alieve::ITSDigitsInfo* di = new Alieve::ITSDigitsInfo();
  di->SetTree(dt);

  gStyle->SetPalette(1, 0);
  gReve->DisableRedraw();

  Alieve::ITSModuleStepper* store = new Alieve::ITSModuleStepper(di);
  store->SetStepper(col, row, 10, 8);
  store->SetMainColor((Color_t)2);
  gReve->AddRenderElement(store);
  stepper = store;
  
  TRandom r(0);
  Int_t module;
  // SPD
  for (Int_t i=0; i<40; ++i) {
    module = r.Integer(240);
    store->AddToList(module);
  }
  // SDD
  for (Int_t i=0; i<40; ++i) {
    module = 240 + r.Integer(260);
    store->AddToList(module);
  }
  // SSD
  for (Int_t i=0; i<40; ++i) {
    module = 500 + r.Integer(1600);
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
