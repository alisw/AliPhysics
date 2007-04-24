void monrec() {
  // MonaLisa monitoring
  gSystem->Load("libNet.so");
  gSystem->Load("libMonaLisa.so");

  SysInfo_t info;
  gSystem->GetSysInfo(&info);

  TString platform(info.fOS);
  platform += ".";
  platform += info.fCpuType;
  new TMonaLisaWriter(platform.Data(),"Reconstruction pp","aliendb3.cern.ch");


  gROOT->LoadMacro("rec.C");
  rec();
  gMonitoringWriter->SendProcessingProgress(1,1,kTRUE);  

  // Send the size of the AliESDs.root file

  FileStat_t buf;
  gSystem->GetPathInfo("./AliESDs.root",buf);

  TList *valuelist = new TList();
  valuelist->SetOwner(kTRUE);

  TMonaLisaValue* valdouble = new TMonaLisaValue("AliESDs.root size",buf.fSize);
  valuelist->Add(valdouble);

  gMonitoringWriter->SendParameters(valuelist);
  delete valuelist;

}
