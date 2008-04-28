void simqa()
{
  const char * kYear = "08" ; 
  gEnv->SetValue("Root.Stacktrace","no");
  gEnv->SetValue("Root.Stacktrace","no");
  AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), kYear)) ;
  AliQA::SetQARefDataDirName("Sim") ; //Data, Pedestals, BlackEvent, .....
  
  TString detectors("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD PMD ZDC T0 VZERO"); 

  AliQADataMakerSteer qas ; 
  qas.Run(detectors.Data(), AliQA::kHITS);
  qas.Run(detectors.Data(), AliQA::kSDIGITS);
  qas.Run(detectors.Data(), AliQA::kDIGITS);
}
