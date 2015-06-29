void MakeVZEROLightYieldsEntryRun2()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://./OCDB");

  // Creation of the light yields OCDB object
  Double_t lightYieldCorr[66] = {0.0,
				 0.020519 , 0.018930 , 0.020560 , 0.020197 , 0.025617 , 0.019526 , 0.017997 , 0.015505 ,
				 0.026803 , 0.014885 , 0.021242 , 0.021103 , 0.014913 , 0.018426 , 0.020704 , 0.021788 ,
				 0.022357 , 0.038628 , 0.020165 , 0.017078 , 0.026363 , 0.026287 , 0.021607 , 0.025717 ,
				 0.022989 , 0.028829 , 0.032443 , 0.034666 , 0.029357 , 0.046145 , 0.025073 , 0.029776 ,
				 0.002412 , 0.003164 , 0.002870 , 0.002362 , 0.000780 , 0.001051 , 0.001273 , 0.003763 ,
				 0.003230 , 0.003231 , 0.003645 , 0.002367 , 0.002021 , 0.002566 , 0.001976 , 0.002077 ,
				 0.003166 , 0.001915 , 0.001013 , 0.003088 , 0.002644 , 0.003215 , 0.004039 , 0.003043 ,
				 0.007530 , 0.004773 , 0.004048 , 0.006647 , 0.004005 , 0.003325 , 0.007207 , 0.006475 ,
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
