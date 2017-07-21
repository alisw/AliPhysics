AliGenerator* AddMCGenPythia8(Float_t e_cms) {
    
   gSystem->Load("libpythia6");
   gSystem->Load("liblhapdf");
   // loaded automatically -> gSystem->Load("libpythia8");

   gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
   gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
   gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));


  AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
  
  //Centre of mass energy 
  gener->SetEnergyCMS(e_cms); // in GeV

  // Event list
  gener->SetEventListRange(-1, -1);

  //tune of PYTHIA
  //(AliPythia8::Instance())->ReadString("Tune:ee = 0");

  //random seed based on time
  AliPythia8::Instance()->ReadString("Random:setSeed = on");
  AliPythia8::Instance()->ReadString("Random:seed = 0");

  //color reconnection
  //if(kCR)             
  //(AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = on");
  //else
  //(AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = off");
  
 
  //AliPythia8::Instance()->ReadString(Form("MultipartonInteractions:kFactor = %i", kF));
  
  return gener;
}
