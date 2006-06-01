void muon_cpads()
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();

  rl->LoadTracks("MUON");

  TTree* tt = rl->GetTreeT("MUON", false);

  Alieve::MUONDigitsInfo* di = new Alieve::MUONDigitsInfo();
  di->SetTTree(tt);

  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  AliMUON *pMUON = (AliMUON*)gAlice->GetModule("MUON");
  const AliMUONGeometryTransformer* kGeomTransformer = pMUON->GetGeometryTransformer();

  gStyle->SetPalette(1, 0);

  gReve->DisableRedraw();

  char name[128];
  char title[128];

  /* CHAMBERS */

  sprintf(name,"M-Chamber pads");
  Reve::RenderElementList* lc = new Reve::RenderElementList(name);
  lc->SetTitle(title);
  lc->SetMainColor((Color_t)2);
  TGListTreeItem *tic = gReve->AddRenderElement(lc);
  /*
  for (Int_t i = 1; i <= 14; i++) {

    Alieve::MUONModule* m = new Alieve::MUONModule(i, -2, di, 0, 0, (Color_t)2);
    lc->AddElement(m);
    gReve->AddRenderElement(tic,m);      
       
  }
  */
  // draw pads in the trigger chambers
  for (Int_t i = 11; i <= 14; i++) {

    Alieve::MUONModule* m1 = new Alieve::MUONModule(i, -3, di, 0, 0, (Color_t)2);
    lc->AddElement(m1);
    gReve->AddRenderElement(tic,m1);      
       
    Alieve::MUONModule* m2 = new Alieve::MUONModule(i, -4, di, 0, 0, (Color_t)2);
    lc->AddElement(m2);
    gReve->AddRenderElement(tic,m2);      
       
  }
  
  gReve->DrawRenderElement(lc);

  gReve->EnableRedraw();

}
