// ---------------------------------------------------------------------- //
void pmd_digits(Int_t mode = 0)
{
  Bool_t drawBackPolygons = kFALSE;

  gStyle->SetPalette(1, 0);


  TString spl;

  TString sddl;
  TString bsddl="DDL";


  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadDigits("PMD");
  TTree* pmdt = rl->GetTreeD("PMD", false);

  //  cout << pmdt->GetEntries() << endl;

  gReve->DisableRedraw();

  Reve::RenderElementList* l = new Reve::RenderElementList("PMD");
  // l->SetTitle("tooltip");
  // l->SetMainColor((Color_t)3);
  gReve->AddRenderElement(l);

  Reve::RGBAPalette* pal = new Reve::RGBAPalette(20, 1000);
  pal->SetLimits(0, 1024);

  Int_t NSM         = 0;
  Int_t istartDDL   = 0;
  Int_t iendDDL     = 0;
  Int_t modnumber   = 0;
  Int_t istartPlane = 0;
  Int_t iendPlane   = 0;
  Float_t zpos      = 0;

  switch(mode)
    {
    case 0:
      istartPlane = 0;
      iendPlane   = 1;
      printf("--- Visualization is set for PREshower Plane ---\n");
      break;

    case 1:
      istartPlane = 1;
      iendPlane   = 2;
      printf("--- Visualization is set for CPV Plane ---\n");
      break;

    case 2:
      istartPlane = 0;
      iendPlane   = 2;
      printf("--- Visualization is set for both Planes ---\n");
      break;

    default:
      printf("--- Not set for any Plane ---\n");
    }

  for (Int_t ipl = istartPlane; ipl < iendPlane; ipl++)
    {

      if (ipl == 0)
	{
	  spl       = "PRE";
	  istartDDL = 0;
	  iendDDL   = 4;
	  zpos      = 365.;
	}
      if (ipl == 1)
	{
	  spl       = "CPV";
	  istartDDL = 4;
	  iendDDL   = 6;
	  zpos      = 360.;
	}
      
      Reve::RenderElementList* lplane = new Reve::RenderElementList(spl.Data());
      //  l->SetMainColor((Color_t)3);
      gReve->AddRenderElement(lplane, l);
      
      for (Int_t iddl = istartDDL; iddl < iendDDL; iddl++)
	{
	  sddl = bsddl;
	  sddl += iddl;
	  Reve::RenderElementList* lddl = new Reve::RenderElementList(sddl.Data());
	  //  l->SetMainColor((Color_t)3);
	  gReve->AddRenderElement(lddl, lplane);

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
	      lmodule->SetPosition(0.,0.,zpos);
	      lmodule->DisplayDigitsData(modnumber, pmdt);
              lmodule->SetPalette(pal);
              if (drawBackPolygons)
              {
                Reve::FrameBox* b = lmodule->GetFrame();
                b->SetFrameWidth(1.5);
                b->SetFrameColor((Color_t) 1);
                b->SetBackColor ((Color_t) (kTeal - 9));
                b->SetFrameFill (kFALSE);
                b->SetDrawBack  (kTRUE);
              }
	      gReve->AddRenderElement(lmodule, lddl);
	      modnumber++;
	      if (iddl == 4 && modnumber == 30) modnumber = 42;
	    }

	}

    }


  gReve->EnableRedraw();
  
  gReve->Redraw3D();
}

// ---------------------------------------------------------------------- //
