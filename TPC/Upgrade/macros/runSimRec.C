void runSimRec(Int_t simtype, Int_t SCtype, Int_t nevents, Int_t ntracks, Int_t rate=50)
{
  //rate is in kHz

  Int_t recoType=simtype/10;
  recoType%=10;
  Int_t subRecoType=simtype/100;
  simtype%=10;
  
  //simulation part
  AliToyMCEventGeneratorSimple s;

  TString outputFile="toyMC";

  //for simtype also below
  switch (simtype) {
    case 0:
      outputFile.Append(Form("_fixed_%dkHz",rate));
      break;
    case 1:
      outputFile.Append(Form("_train_%dkHz",rate));
      break;
  }
  
  switch (SCtype) {
    case 0:
      s.SetSpaceCharge(AliToyMCEventGeneratorSimple::kEps5);
      outputFile.Append("_eps05");
      break;
    case 1:
      s.SetSpaceCharge(AliToyMCEventGeneratorSimple::kEps10);
      outputFile.Append("_eps10");
      break;
    case 2:
      s.SetSpaceCharge(AliToyMCEventGeneratorSimple::kEps20);
      outputFile.Append("_eps20");
      break;
  }

  outputFile.Append(Form("_%04dev_%04dtr",nevents,ntracks));
  outputFile.Append(".root");
  s.SetOutputFileName(outputFile.Data());

  //TODO: Add other types
  switch (simtype) {
    case 0:
      s.RunSimulation(nevents,ntracks,rate);
      break;
    case 1:
      s.RunSimulationBunchTrain(nevents,ntracks);
      break;
  }

  //reconstruction part
  AliToyMCReconstruction rec;
  // rec.SetUseMaterialBudget(kTRUE)

  if (recoType==0){
    rec.SetRecoSettings(1,0,AliToyMCReconstruction::kNoCorrection);
    if (!subRecoType||subRecoType==1) rec.RunReco(outputFile.Data());

    rec.SetRecoSettings(1,1,AliToyMCReconstruction::kIdeal);
    if (!subRecoType||subRecoType==2) rec.RunReco(outputFile.Data());

    rec.SetRecoSettings(0,1,AliToyMCReconstruction::kIdeal);
    if (!subRecoType||subRecoType==3) rec.RunReco(outputFile.Data());

    rec.SetRecoSettings(0,1,AliToyMCReconstruction::kAverageEta);
    if (!subRecoType||subRecoType==4) rec.RunReco(outputFile.Data());

    rec.SetRecoSettings(0,1,AliToyMCReconstruction::kNoCorrection);
    if (!subRecoType||subRecoType==5) rec.RunReco(outputFile.Data());

    rec.SetRecoSettings(0,0,AliToyMCReconstruction::kNoCorrection);
    if (!subRecoType||subRecoType==6) rec.RunReco(outputFile.Data());
  }

  if (recoType==1) {
    rec.SetRecoSettings(0,1,AliToyMCReconstruction::kNoCorrection);
    if (!subRecoType||subRecoType==1) rec.RunFullTracking(outputFile.Data());
    
    rec.SetRecoSettings(0,0,AliToyMCReconstruction::kNoCorrection);
    if (!subRecoType||subRecoType==2) rec.RunFullTracking(outputFile.Data());
  }
  
}
