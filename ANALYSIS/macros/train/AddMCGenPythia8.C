AliGenerator* AddMCGenPythia8(Float_t e_cms = 2760., Bool_t kCR = kTRUE, Int_t kF = 1) 
{
  // Add Pythia 8 generator: 
  //    - Color reconnection = ON/OFF
  //    - Set k factor, default = 1; range of possible values in xmldoc/CouplingsAndScales.xml

  gSystem->Load("liblhapdf");
 
  AliGenerator *genP = NULL;
  genP = CreatePythia8Gen(e_cms, kCR, kF);
  
  return genP;
}

AliGenerator* CreatePythia8Gen(Float_t e_cms, Bool_t kCR, Int_t kF) {
    
   gSystem->Load("libpythia6");
   gSystem->Load("libEGPythia6");
   gSystem->Load("libAliPythia6");
   gSystem->Load("libpythia8");
   gSystem->Load("libAliPythia8");
   gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8175/xmldoc"));
   gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
   gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));


  AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());

  // set process (MB)
  gener->SetProcess(kPyMbDefault);
  
  //   Centre of mass energy 
  gener->SetEnergyCMS(e_cms); // in GeV

  // Event list
  gener->SetEventListRange(-1, 2);

  // color reconnection
  (AliPythia8::Instance())->ReadString("Tune:pp = 5");//CR

  //random seed based on time
  AliPythia8::Instance()->ReadString("Random:setSeed = on");
  AliPythia8::Instance()->ReadString("Random:seed = 0");

  if(kCR)             
    (AliPythia8::Instance())->ReadString("BeamRemnants:reconnectColours = on");
  else
    (AliPythia8::Instance())->ReadString("BeamRemnants:reconnectColours = off");
  
 
  AliPythia8::Instance()->ReadString(Form("MultipartonInteractions:kFactor = %i", kF));
  
  return gener;
}
