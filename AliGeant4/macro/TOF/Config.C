void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTOF");

  AliTOF *TOF = 0;
  switch (version) {
    case 0: TOF  = new AliTOFv0("TOF", "TOFv0 detector"); break;
    case 1: TOF  = new AliTOFv1("TOF", "TOFv1 detector"); break;
    case 2: TOF  = new AliTOFv2("TOF", "TOFv2 detector"); break;
    case 3: TOF  = new AliTOFv3("TOF", "TOFv3 detector"); break;
    case 4: TOF  = new AliTOFv4("TOF", "TOFv4 detector"); break;
    case 5: if (AliRunConfiguration::Holes())
              TOF = new AliTOFv2FHoles("TOF", "TOF with Holes");
	    else  
	      TOF = new AliTOFv4T0("TOF", "normal TOF");
	    break;
  }   

//=================== TOF parameters ============================
}  
