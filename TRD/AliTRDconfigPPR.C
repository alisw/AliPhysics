if (iTRD) {

  //=================== TRD parameters ============================

  AliTRD *TRD  = new AliTRDv1("TRD","TRD slow simulator");

  // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 90% Xe + 10% CO2)
  TRD->SetGasMix(1);

  // Switch on TR
  AliTRDsim *TRDsim = TRD->CreateTR();

}
