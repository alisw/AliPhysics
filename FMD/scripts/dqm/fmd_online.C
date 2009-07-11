void fmd_online(const char *input = "^FMD")
{
  gROOT->LoadMacro("alieve_init.C");
  gROOT->LoadMacro("fmd_raw.C");
  //          path, event, esdFile, aodFile, rawFile, cdbUri,  
  alieve_init("",   0,     0,       0,       input,   0, 
	      //rl,   esd,    aod,    raw
 	      kFALSE, kFALSE, kFALSE, kTRUE);

  gAliEveEvent->SetAutoLoad(kTRUE);
  // gAliEveEvent->SetAutoLoad(kFALSE);
  gAliEveEvent->AddNewEventCommand("fmd_raw()");

  //  while(true) 
  //gSystem->Sleep(100);


}
