void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libTRD");

  AliTRD *TRD = 0;
  switch (version) {
    case 0: TRD  = new AliTRDv0("TRD", "TRDv0 detector"); break;
    case 1: TRD  = new AliTRDv1("TRD","TRD version 0");   break;
  }

//=================== TRD parameters ============================
  
  //TRD->SetHits();
  
  //AliTRD *TRD  = new AliTRDv1("TRD","TRD slow simulator");
  //TRD->SetSensPlane(0);
  //TRD->SetSensChamber(2);
  //TRD->SetSensSector(17);
  
  // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
  TRD->SetGasMix(1);
  
  // With hole in front of PHOS
  TRD->SetPHOShole();
  // With hole in front of RICH
  TRD->SetRICHhole();
}
