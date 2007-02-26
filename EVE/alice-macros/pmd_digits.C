// ---------------------------------------------------------------------- //
void pmd_digits()
{
  gStyle->SetPalette(1, 0);

  gReve->DisableRedraw();


  TString spl;

  TString sddl;
  TString bsddl="DDL";


  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadDigits("PMD");
  TTree* pmdt = rl->GetTreeD("PMD", false);

  //  cout << pmdt->GetEntries() << endl;

  Reve::RenderElementList* l = new Reve::RenderElementList("PMD");
  // l->SetTitle("tooltip");
  // l->SetMainColor((Color_t)3);
  gReve->AddRenderElement(l);
  
  Int_t istartDDL = 0;
  Int_t iendDDL    = 0;
  Int_t modnumber = 0;
  Int_t NSM       = 0;

  for (Int_t ipl = 0; ipl < 1; ipl++)
    {

      if (ipl == 0)
	{
	  spl = "PRE";
	  istartDDL = 0;
	  iendDDL   = 4;
	}
      if (ipl == 1)
	{
	  spl = "CPV";
	  istartDDL = 4;
	  iendDDL   = 6;
	}
      
      Reve::RenderElementList* lplane = new Reve::RenderElementList(spl.Data());
      //  l->SetMainColor((Color_t)3);
      gReve->AddRenderElement(l,lplane);
      
      for (Int_t iddl = istartDDL; iddl < iendDDL; iddl++)
	{
	  sddl = bsddl;
	  sddl += iddl;
	  Reve::RenderElementList* lddl = new Reve::RenderElementList(sddl.Data());
	  //  l->SetMainColor((Color_t)3);
	  gReve->AddRenderElement(lplane,lddl);

	  modnumber = iddl*6;

	  if (iddl < 4)
	    {
	      NSM = 6;
	    }
	  else if (iddl >=4 && iddl < 6)
	    {
	      NSM = 12;
	    }

	  for (Int_t ism = 0; ism < NSM; ism++)
	    {
	      
	      Alieve::PMDModule *lmodule = new Alieve::PMDModule();
	      lmodule->SetPosition(0.,0.,360.);
	      lmodule->DisplayDigitsData(modnumber, pmdt);
	      gReve->AddRenderElement(lddl, lmodule);
	      modnumber++;
	      if (iddl == 4 && modnumber == 30) modnumber = 42;
	    }

	}

    }








  gReve->EnableRedraw();
  
  gReve->Redraw3D();
}

// ---------------------------------------------------------------------- //
