void ViewFMD()
{
  gMC->Gsatt("FMD1","seen",0);
  gMC->Gsatt("FMD2","seen",0);
  gMC->Gsatt("FMD3","seen",0);

  TString name;
  // Rings
  for (Int_t i = 0; i < 2; i++) {
    char c;
    switch (i) {
    case 0: c = 'I'; break;
    case 1: c = 'O'; break;
    }
     
    name = Form("FMD%c", c);
    gMC->Gsatt(name.Data(),"seen",0); // Ring volume	     

    name = Form("F%cFV", c);
    gMC->Gsatt(name.Data(),"seen",0); // Virtual volume front

    name = Form("F%cBV", c);
    gMC->Gsatt(name.Data(),"seen",0); // Virtual volume back

    name = Form("F%cAC", c);
    gMC->Gsatt(name.Data(),"seen",-2); // Active volume

    name =  Form("F%cSL", c);
    gMC->Gsatt(name.Data() ,"seen",1);

    name =  Form("F%cLL", c);
    gMC->Gsatt(name.Data() ,"seen",1);

    // name = Form("F%cAP", c);
    // gMC->Gsatt(name.Data(),"seen",-1); // Phi segmentation of active

    // name = Form("F%cAR", c);
    // gMC->Gsatt(name.Data(),"seen",-1); // R segmentation of active

    name = Form("F%cPT", c);
    gMC->Gsatt(name.Data(),"seen",1); // Top of print-board

    name = Form("F%cPB", c);
    gMC->Gsatt(name.Data(),"seen",1); // Bottom of print board
  }
   
  for (Int_t i = 1; i <= 3; i++) {
    for (Int_t j = 0;  j < 2; j++) {
      if (i == 1 && j == 1) break;
      char c;
      switch (j) {
      case 0: c = 'I'; break;
      case 1: c = 'O'; break;
      }
       
      name = Form("F%d%cI", i, c);
      gMC->Gsatt(name.Data(),"seen",-2); // Honeycomp top 

      name = Form("F%d%cJ", i, c);
      gMC->Gsatt(name.Data(),"seen",-2); // Honeycomp bottom

      name = Form("F%d%cK", i, c);
      gMC->Gsatt(name.Data(),"seen",0); // Honeycomp inner top 

      name = Form("F%d%cL", i, c);
      gMC->Gsatt(name.Data(),"seen",0); // Honeycomp inner bottom 
    }
  }

  gMC->Gsatt("F3SN", "seen", 1); // Nose of FMD3 Cone
  gMC->Gsatt("F3SB", "seen", 1); // Back of FMD3 Cone
  gMC->Gsatt("F3SL", "seen", 1); // Beams of FMD3 Cone
  gMC->Gsatt("F3SF", "seen", 1); // Flanges on FMD3 Cone
}
