void runSimRec(Int_t simtype, Int_t SCtype, Int_t nevents, Int_t ntracks, Int_t rate=50)
{
  //rate is in kHz

  //simulation part
  AliToyMCEventGeneratorSimple s;

  TString outputFile="toyMC";

  //for simtype also below
  switch (simtype) {
    case 0:
      outputFile.Append(Form("_fixed_%dkHz",rate));
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
  }

  //reconstruction part
  AliToyMCReconstruction rec;
  // rec.SetUseMaterialBudget(kTRUE)
  
  rec.SetRecoSettings(1,0,AliToyMCReconstruction::kNoCorrection);
  rec.RunFullTracking(outputFile.Data());
  
  rec.SetRecoSettings(1,1,AliToyMCReconstruction::kIdeal);
  rec.RunFullTracking(outputFile.Data());
  
  rec.SetRecoSettings(0,1,AliToyMCReconstruction::kIdeal);
  rec.RunFullTracking(outputFile.Data());
  
  rec.SetRecoSettings(0,1,AliToyMCReconstruction::kAverageEta);
  rec.RunFullTracking(outputFile.Data());
  
  rec.SetRecoSettings(0,1,AliToyMCReconstruction::kNoCorrection);
  rec.RunFullTracking(outputFile.Data());
  
  rec.SetRecoSettings(0,0,AliToyMCReconstruction::kNoCorrection);
  rec.RunFullTracking(outputFile.Data());
  
}
