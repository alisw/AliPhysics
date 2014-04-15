AliGenerator* AddMCGenPythia8(Float_t e_cms = 2760., Bool_t kCR = kTRUE) 
{
  // Add Pythia 8 generator: 
  //    - Color reconnection = ON/OFF

  gSystem->Load("liblhapdf.so");
 
  AliGenerator *genP = NULL;
  genP = CreatePythia8Gen(e_cms, kCR);
  
  return genP;
}

AliGenerator* CreatePythia8Gen(Float_t e_cms, Bool_t kCR) {
    
   gSystem->Load("libpythia6.so");
   gSystem->Load("libAliPythia6.so");
   gSystem->Load("libpythia8.so");
   gSystem->Load("libAliPythia8.so");
   gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8175/xmldoc"));
   gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
   gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));


  AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());

  // set process (MB)
  gener->SetProcess(kPyMbDefault);
  // //   charm, beauty, charm_unforced, beauty_unforced, jpsi, jpsi_chi, mb
  // if(ptHardMin>0.) {
  //   gener->SetProcess(kPyJets);
  //   gener->SetPtHard((float)ptHardMin,(float)ptHardMax);
  // } else
  //   gener->SetProcess(kPyMb); // Minimum Bias  
  
  //   Centre of mass energy 
  gener->SetEnergyCMS(e_cms); // in GeV

  // Event list
  gener->SetEventListRange(-1, 2);

  // color reconnection
  (AliPythia8::Instance())->ReadString("Tune:pp = 5");//CR
  //gener->UseNewMultipleInteractionsScenario(); // for all Pythia versions >= 6.3

  if(kCR)             
    (AliPythia8::Instance())->ReadString("BeamRemnants:reconnectColours = on");
  else
    (AliPythia8::Instance())->ReadString("BeamRemnants:reconnectColours = off");


  //   vertex position and smearing 
  //gener->SetVertexSmear(kPerEvent);

  // structure function
  // use kCTEQ5l for Perugia tunes 
  // except for tunes: Perugia * (325, MRSTLO*), Perugia 6 (326, CTEQ6L),  
  // Perugia 11 M (355, MRST LO**), Perugia 11 C (356, CTEQ6L1) 
  //gener->SetStrucFunc(kCTEQ5L); 
   
  return gener;
}
