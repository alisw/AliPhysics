void runSimRec(Int_t simtype, Int_t SCtype, Int_t nevents, Int_t ntracks)
{

  //simulation part
  AliToyMCEventGeneratorSimple s;

  TString outputFile="toyMC";
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

  outputFile.Append(Form("_%04dev_%04dtr"));
  outputFile.Append(".root");
  s.SetOutputFileName(outputFile.Data());

  //TODO: Add other types
  switch (simtype) {
    case 0:
      s.RunSimulation(nevents,ntracks);
      break;
  }

  //reconstruction part
  AliToyMCReconstruction rec;
  // rec.SetUseMaterialBudget(kTRUE)
  
  rec.SetRecoSettings(1,0,AliToyMCReconstruction::kNoCorrection);
  rec.RunReco(outputFile.Data(),-1);
  
  rec.SetRecoSettings(1,1,AliToyMCReconstruction::kIdeal);
  rec.RunReco(outputFile.Data(),-1);
  
  rec.SetRecoSettings(0,1,AliToyMCReconstruction::kIdeal);
  rec.RunReco(outputFile.Data(),-1);
  
  rec.SetRecoSettings(0,1,AliToyMCReconstruction::kAverageEta);
  rec.RunReco(outputFile.Data(),-1);
  
  rec.SetRecoSettings(0,1,AliToyMCReconstruction::kNoCorrection);
  rec.RunReco(outputFile.Data(),-1);
  
  rec.SetRecoSettings(0,0,AliToyMCReconstruction::kNoCorrection);
  rec.RunReco(outputFile.Data(),-1);
  
}
