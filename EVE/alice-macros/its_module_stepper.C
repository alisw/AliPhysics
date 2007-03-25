#include "TGLViewer.h"
namespace Alieve
{
class ITSModuleStepper;
}

void its_module_stepper(Int_t col = 4 , Int_t row = 4)
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadDigits("ITS");
  TTree* dt = rl->GetTreeD("ITS", false);
  Alieve::ITSDigitsInfo* di = new Alieve::ITSDigitsInfo();
  di->SetTree(dt);

  gStyle->SetPalette(1, 0);
  gReve->DisableRedraw();


  Int_t sw = 10;
  Int_t sh = 8;
  Float_t CW = col*sw;
  Float_t CH = col*sh;
   
  Float_t off_x = CW*0.07;
  Float_t off_y = CH*0.1;

  Reve::ZTrans* mx;

  Bool_t rnrFrame = kTRUE;
  Color_t wcol = 41;  
  
  // SPD  
  Alieve::ITSModuleStepper* spd_lay1 = new Alieve::ITSModuleStepper(di);
  spd_lay1->SetStepper(col, row, sw, sh);
  spd_lay1->DisplayDet(0, 1);
  spd_lay1->SetName("SPD 1");
  spd_lay1->SetRnrFrame(rnrFrame);
  spd_lay1->SetWColor(wcol);
  gReve->AddRenderElement(spd_lay1);

  Alieve::ITSModuleStepper* spd_lay2 = new Alieve::ITSModuleStepper(di);
  spd_lay2->SetStepper(col, row, sw, sh);
  mx = spd_lay2->PtrMainHMTrans();
  mx->SetPos(CW+off_x, 0, 0);
  spd_lay2->DisplayDet(0, 2);
  spd_lay2->SetName("SPD 2");
  spd_lay2->SetRnrFrame(rnrFrame);
  spd_lay2->SetWColor(wcol);
  gReve->AddRenderElement(spd_lay2);
  
  // SDD  
  Alieve::ITSModuleStepper* sdd_lay1 = new Alieve::ITSModuleStepper(di);
  sdd_lay1->SetStepper(col, row, 10, 8);
  mx = sdd_lay1->PtrMainHMTrans();
  mx->SetPos(0, CH+off_y, 0);
  sdd_lay1->DisplayDet(1, 3);
  sdd_lay1->SetName("SDD 1");
  sdd_lay1->SetRnrFrame(rnrFrame);
  sdd_lay1->SetWColor(wcol);
  gReve->AddRenderElement(sdd_lay1);

  Alieve::ITSModuleStepper* sdd_lay2 = new Alieve::ITSModuleStepper(di);
  sdd_lay2->SetStepper(col, row, 10, 8); 
  mx = sdd_lay2->PtrMainHMTrans();
  mx->SetPos(CW+off_x, CH+off_y, 0);
  sdd_lay2->DisplayDet(1, 4);
  sdd_lay2->SetName("SDD 2");
  sdd_lay2->SetRnrFrame(rnrFrame);
  sdd_lay2->SetWColor(wcol);
  gReve->AddRenderElement(sdd_lay2);
  

  // SSD  
  Alieve::ITSModuleStepper* ssd_lay1 = new Alieve::ITSModuleStepper(di);
  ssd_lay1->SetStepper(col, row, 10, 8);
  mx = ssd_lay1->PtrMainHMTrans();
  mx->SetPos(0, 2*(CH+off_y), 0);
  ssd_lay1->DisplayDet(2, 5);
  ssd_lay1->SetName("SSD 1");
  ssd_lay1->SetRnrFrame(rnrFrame);
  ssd_lay1->SetWColor(wcol);
  gReve->AddRenderElement(ssd_lay1);

  Alieve::ITSModuleStepper* ssd_lay2 = new Alieve::ITSModuleStepper(di);
  ssd_lay2->SetStepper(col, row, 10, 8); 
  mx = ssd_lay2->PtrMainHMTrans();
  mx->SetPos(CW+off_x, 2*(CH+off_y), 0);
  ssd_lay2->DisplayDet(2, 6);
  ssd_lay2->SetName("SSD 2");
  ssd_lay2->SetRnrFrame(rnrFrame);
  ssd_lay2->SetWColor(wcol);
  gReve->AddRenderElement(ssd_lay2);

  Alieve::DigitScaleInfo* si =  Alieve::ITSScaledModule::fgDigitScaleInfo;
  si->ScaleChanged(2);

  gReve->EnableRedraw();
  TGLViewer * v = (TGLViewer *)gPad->GetViewer3D();
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
}

