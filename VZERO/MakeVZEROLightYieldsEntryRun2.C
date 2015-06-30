void MakeVZEROLightYieldsEntryRun2()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://./OCDB");

  // Creation of the light yields OCDB object
  Double_t lightYieldCorr[66] = {0.0,
				 0.013680 , 0.012620 , 0.013706 , 0.013465 , 0.017078 , 0.013017 , 0.011998 , 0.010337 ,
				 0.017869 , 0.009923 , 0.014161 , 0.014068 , 0.009942 , 0.012284 , 0.013803 , 0.014526 ,
				 0.014905 , 0.025752 , 0.013444 , 0.011385 , 0.017576 , 0.017524 , 0.014404 , 0.017145 ,
				 0.015326 , 0.019219 , 0.021629 , 0.023111 , 0.019571 , 0.030764 , 0.016716 , 0.019850 ,
				 0.001608 , 0.002109 , 0.001913 , 0.001575 , 0.000520 , 0.000701 , 0.000849 , 0.002509 ,
				 0.002153 , 0.002154 , 0.002430 , 0.001578 , 0.001347 , 0.001711 , 0.001976 , 0.001385 ,
				 0.002110 , 0.001277 , 0.000676 , 0.002059 , 0.001763 , 0.002143 , 0.002693 , 0.002029 ,
				 0.005020 , 0.003182 , 0.002699 , 0.004431 , 0.002670 , 0.002217 , 0.004805 , 0.004317 ,
				 0.0};

  TH1F *yields = new TH1F("VZEROLightYields","VZERO Light Yields",64,-0.5,63.5);
  yields->SetContent(lightYieldCorr);
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Light Yields channel by channel for Run2");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("VZERO/Calib/LightYields",0,AliCDBRunRange::Infinity());

  storLoc->Put(yields, id, md);

  storLoc->Delete();
  delete md;

}
