enum ECalorimeterAcceptance_t {
  kCalorimeterAcceptance_FullDetector = 0,
  kCalorimeterAcceptance_EMCRun1,
  kCalorimeterAcceptance_PHSRun1,
  kCalorimeterAcceptance_EMCRun2,
  kCalorimeterAcceptance_PHSRun2,
  kCalorimeterAcceptance_PHSDMC
};

void 
GetCalorimeterAcceptance(Int_t acceptance, Float_t &etaMax, Float_t &phiMin, Float_t &phiMax)
{
  switch (acceptance) 
    {
    case kCalorimeterAcceptance_FullDetector:
      etaMax = 1.5 ; phiMin =   0.; phiMax = 360.;
      break;   
    case kCalorimeterAcceptance_EMCRun1:
      etaMax = 0.7 ; phiMin =  80.; phiMax = 180.;
      break;
    case kCalorimeterAcceptance_EMCRun2:
      etaMax = 0.7 ; phiMin =  80.; phiMax = 187.;
      break;    
    case kCalorimeterAcceptance_PHSRun1:
      etaMax = 0.13; phiMin = 260.; phiMax = 320.;
      break;
    case kCalorimeterAcceptance_PHSRun2:
      etaMax = 0.13; phiMin = 250.; phiMax = 320.;
      break;    
    case kCalorimeterAcceptance_PHSDMC:
      etaMax = 0.7 ; phiMin = 250.; phiMax = 327.;
      break;
    default:
      printf("Error in acceptance switch\n");
    }
  
  printf("\t Calorimeter acceptance for %d: |eta|<%2.2f - %2.2f<phi<%2.2f\n",acceptance,etaMax,phiMin,phiMax);
}

AliGenerator *
GeneratorPythia8(Int_t tune, Float_t energyConfig)
{
  //
  // Libraries
  gSystem->Load("libpythia8.so");
  gSystem->Load("libAliPythia8.so");
  //
  // Environment settings
  gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
  gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
  gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
  //
  // Pythia
  AliGenPythiaPlus *pythia = new AliGenPythiaPlus(AliPythia8::Instance()); 
  pythia->SetMomentumRange(0, 999999.);
  pythia->SetThetaRange(0., 180.);
  pythia->SetYRange(-12.,12.);
  pythia->SetPtRange(0,1000.);
  pythia->SetProcess(kPyMbDefault); // pythia->SetProcess(kPyMb);
  pythia->SetEnergyCMS(energyConfig);
  //
  // Initialize
  //pythia->SetEventListRange(-1, 2); 
  (AliPythia8::Instance())->ReadString("Random:setSeed = on");
  (AliPythia8::Instance())->ReadString("Random:seed = 0");
  (AliPythia8::Instance())->ReadString("111:mayDecay = on");
  //
  // Tune
  if (tune > 0) {
    pythia->SetTune(tune);
  }
  //
  return pythia;
}

AliGenerator* Add_MCGenPythia8_JetJet(  Float_t e_cms       = 2760., 
                                        Int_t tune          = 5,
                                        Int_t acceptance    = 0,
                                        Double_t pthardminConfig = 0., 
                                        Double_t pthardmaxConfig  = 1.
                                    ) {
  // Pythia
  AliGenPythiaPlus *pythia = (AliGenPythiaPlus*) GeneratorPythia8(tune, e_cms);
  //
  // jets settings
  pythia->SetProcess(kPyJets);
  Float_t etaMax, phiMin, phiMax;
  GetCalorimeterAcceptance(acceptance, etaMax, phiMin, phiMax);
  pythia->SetJetEtaRange(-etaMax, etaMax); // Final state kinematic cuts
  pythia->SetJetPhiRange(phiMin, phiMax);
  pythia->SetJetEtRange(0., 1000.);
  pythia->SetPtHard(pthardminConfig, pthardmaxConfig); // Pt transfer of the hard scattering
  pythia->SetStrucFunc(kCTEQ5L);
  //
  return pythia;
}
