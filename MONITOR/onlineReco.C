UInt_t onlineReco(const char* param = "listen",const char *recMacroPath = "$ALICE_ROOT/test/cosmic/rec.C") {

  TString paramStr(param);
	
  UInt_t run = 0;

  TStopwatch stopwatch;
  stopwatch.Start();

  if (paramStr.IsDigit()) {
    run = paramStr.Atoi();
  } else if (paramStr == "listen") {
    AliOnlineRecoTrigger trigger;
    run = trigger.Run();
  } else {
    cout<<"Bad parameter: "<<param<<endl;
    cout<<"Parameter options: "<<endl;
    cout<<"<run> - run online reconstruction for the given run"<<endl;
    cout<<"listen - start listening for ECS SOR notification"<<endl;
    cout<<"<empty parameter> - the same as 'listen'"<<endl;
  }

  if (run > 0) {
    TString gdcList;
    if (grp(run) > 0) {

      // "t" stores the token on this disk, otherwise the alien connection is not thread/process-safe
      TGrid::Connect("alien://", 0, 0, "t");

      TObjArray *gdcs = gdcList.Tokenize(" ");
      Int_t ngdcs = tokens->GetEntriesFast();
      if (ngdcs > 0) {

	TString dataSource = ((TObjString*)gdcs->At(0))->String();
	dataSource.Prepend("mem://@");
	datasource.Append(":");

	// Setting CDB
	AliCDBManager * man = AliCDBManager::Instance();
	man->SetDefaultStorage("raw://");
	man->SetSpecificStorage("GRP/GRP/Data",
			      Form("local://%s",gSystem->pwd()));
	man->SetSpecificStorage("GRP/CTP/Config",
			      Form("local://%s",gSystem->pwd()));
	man->Lock();

	gSystem->mkdir(Form("run%d",run));
	gSystem->cd(Form("run%d",run));

	gROOT->LoadMacro(gSystem->ExpandPathName(recMacroPath));
	rec(dataSource.Data());

	AliCDBManager::Destroy();
      }
      else {
	cout << "No GDCs defined in the logbook entry for run " << run << endl;
      }
      delete gdcs;
    }
  }

  printf("Execution time: R:%.2fs C:%.2fs \n",
	 stopwatch.RealTime(),stopwatch.CpuTime());

  return run;

}

Int_t grp(UInt_t run, TString &gdcList) {

  Int_t ret=AliGRPPreprocessor::ReceivePromptRecoParameters(run, "aldaqdb", 0, "LOGBOOK", "logbook", "alice",
							    Form("local://%s",gSystem->pwd()),
							    gdcList);
  if(ret>0) cout << "Last run of the same type is: " << ret << endl;
  else if(ret==0) cout << "No previous run of the same type found" << endl;
  else if(ret<0) cout << "Error code while retrieving GRP parameters returned: " << ret << endl;
  return(ret);
}
