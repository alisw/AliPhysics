void muon_digits()
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();

  rl->LoadDigits("MUON");

  TTree* dt = rl->GetTreeD("MUON", false);

  Alieve::MUONDigitsInfo* di = new Alieve::MUONDigitsInfo();
  di->SetDTree(dt);

  gStyle->SetPalette(1, 0);

  gReve->DisableRedraw();

  rl->LoadgAlice();

  char name[128];
  char title[128];

  /* DIGITS */

  for (Int_t iSta = 1; iSta <= 7; iSta++) {

    for (Int_t iCha = 1; iCha <= 2; iCha++) {

      sprintf(name,"M-ST%1dCH%1d/Digits",iSta,iCha);

      Reve::RenderElementList* l = new Reve::RenderElementList(name);
      
      if (iSta <= 5) {
	sprintf(title,"Station %1d chamber %1d (tracking)",iSta,iCha);
      } else {
	sprintf(title,"Station %1d chamber %1d (trigger)",iSta,iCha);
      }
      
      l->SetTitle(title);
      l->SetMainColor((Color_t)4);
      TGListTreeItem *ti = gReve->AddRenderElement(l);
    
      Int_t iChamber = (iSta-1) * 2 + iCha; 

      AliMpDEIterator it;
      for ( it.First(iChamber-1); ! it.IsDone(); it.Next() ) {
	
	Int_t detElemId = it.CurrentDE();
	
	for (Int_t iCat = 1; iCat <=2; iCat++) {

	  Alieve::MUONModule* m = new Alieve::MUONModule(detElemId, iCat, di, 1, 0, (Color_t)2);
	  l->AddElement(m);
	  gReve->AddRenderElement(ti,m);      
	
	}

      }

      gReve->DrawRenderElement(l);
      
    }
    
  }
  
  gReve->EnableRedraw();

}
