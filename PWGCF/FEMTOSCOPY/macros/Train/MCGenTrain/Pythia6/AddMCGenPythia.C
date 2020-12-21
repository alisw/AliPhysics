//PYTHIA tunes
// tune 320 - Perugia-0
// tune 350 - Perugia 2011
// tune 370 - Perugia 2012

AliGenerator* AddMCGenPythia(Double_t energy = 14000, TString projectile = "p", TString target = "p", Int_t tune = 0, Int_t ntrig = 0, Int_t *trig = NULL) 
{
  //Add Pythia generator: pt-hard bin or min bias

  gSystem->Load("liblhapdf");
  
  if(tune==370)  gSystem->Load("libpythia6_4_28"); // for Perugia 2012
  if(tune<370)  gSystem->Load("libpythia6_4_25");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6");

  
  // Pythia
  AliGenPythia* pythia = new AliGenPythia(-1); 
  pythia->SetMomentumRange(0, 999999.);
  pythia->SetThetaRange(0., 180.);
  pythia->SetYRange(-12.,12.);
  pythia->SetPtRange(0,1000.);
  //pythia->SetProcess(kPyOldPopcorn); //kPyMb
  pythia->SetEnergyCMS(energy);
  pythia->SetCrossingAngle(0,0);
  pythia->SetCollisionSystem(projectile,target);
  //
  // Tune
  if (tune > 0) {
    pythia->SetTune(tune); 
    pythia->UseNewMultipleInteractionsScenario();
  }
  //
  // Trigger particles
  if (ntrig > 0) {
    Int_t pdg = trig[gRandom->Integer(ntrig)];
    comment = comment.Append(Form(" | %s enhanced",
				  TDatabasePDG::Instance()->GetParticle(pdg)->GetName()));
    pythia->SetTriggerParticle(pdg, 1.2);
  }
  //
  return pythia;
}

