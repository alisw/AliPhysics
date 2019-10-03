AliGenerator* AddMCGenDPMJET(Float_t e_cms = 8.)
{
  // Add Pythia generator: pt-hard bin or min bias

  AliGenDPMjet* dpmjet = new AliGenDPMjet(-1);
  dpmjet->SetMomentumRange(0, 999999.);
  dpmjet->SetThetaRange(0., 180.);
  dpmjet->SetYRange(-12., 12.);
  dpmjet->SetPtRange(0, 1000.);
  dpmjet->SetProcess(kDpmMb);
  dpmjet->SetEnergyCMS(e_cms);
  dpmjet->SetCrossingAngle(0, 0.);

  dpmjet->SetProjectile("P", 1, 1);
  dpmjet->SetTarget    ("A", 208, 82);
  dpmjet->SetProjectileBeamEnergy(0.5 * e_cms * TMath::Sqrt(208./82.));
  dpmjet->SetImpactParameterRange(0., 20.);

  return dpmjet;
}
