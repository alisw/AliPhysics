void MakeQAThresholdsEntry(const char* storageUri="local://$ALICE_ROOT/../AliRoot/OCDB", Int_t firstRun=0, Int_t lastRun=999999999)
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(storageUri);
  // QAThresholds
  TObjArray* qaThrArray = new TObjArray();
  for (Int_t idet = 0; idet < AliDAQ::kNDetectors; idet++){
    TString detName = AliDAQ::OnlineName(idet);
    if (detName == "TRI" || detName == "HLT" || detName == "TST") continue;   // skipping TRI, HLT, TST since they do not produce QAThresholds
    Printf("Processing QAThreshold for detector %s",detName.Data()); 
    TString inFile(gSystem->ExpandPathName("$ALICE_ROOT/../AliRoot/GRP/ShuttleInput/"));
    inFile += "run000168322_";
    inFile += detName;
    inFile += "_DQM_QAThresholds";
    Printf("Opening QAThreshold file %s", inFile.Data());
    TFile dqmFile(inFile.Data(),"READ");
    if (dqmFile.IsOpen()) {
      AliQAThresholds* qaThr = dynamic_cast<AliQAThresholds*>(dqmFile.Get(detName.Data()));
      if (qaThr){
        Int_t qaThrId = qaThr->GetDetectorId();
        if (qaThrId != idet){
          Printf("ERROR: Expecting QA threshold for detector %s, but found that for detector %s, skipping",detName.Data(), AliDAQ::OnlineName(qaThrId));
          continue;
        }
        else{
          qaThrArray->AddAtAndExpand(qaThr, qaThrId);
        }
      }
      else {
        Printf("ERROR: No QAThresholds object found in the file for detector %s, skipping",detName.Data());
        continue;
      }
    }
    else {
      Printf("ERROR: Can't open QAThreshold file for detector %s, skipping",detName.Data());
      continue;					
    }
  }
  if (qaThrArray->GetEntries() > 0){
    AliCDBMetaData md;
    md.SetResponsible("Barthélémy von Haller");
    md.SetComment("QA Threshold TObjArray");
    AliCDBId id("GRP/Calib/QAThresholds", firstRun, lastRun);
    cdb->Put(qaThrArray, id, &md); 
  }
  else{
    Printf("No valid QAThresholds entries found, storing nothing in the OCDB");
  }

}
