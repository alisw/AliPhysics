void MakeADPMGainsEntry(const char *outputCDB = "local://$ALICE_ROOT/OCDB")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  // Creation of the PM gains OCDB object
  //2015 before PM termination resistors - up to TS1
  /*/
  Double_t a[16] = {1085.6, 972.9, 1039.2, 986.8, 1085.2, 997.3, 951.9, 1012.7,
  		    999.9, 1055.2, 1080.6, 987.7, 1001.6, 1061.4, 986.1, 938.0};
		    
  Double_t b[16] = {7.052, 5.835, 6.913, 6.755, 7.508, 7.113, 6.547, 6.342, 
  		    6.747, 6.704, 6.413, 7.340, 6.146, 5.957, 7.160, 6.250};
  /*/	    

  //2015 after PM termination resistors 
  Double_t a[16] = {1202.8 ,1083.9 ,1140.7 ,1194.4 ,1135.3 ,1110.0 ,1065.3 ,1102.6,
	   	    1104.9 ,1170.4 ,1210.2 ,1061.4 ,1117.9 ,1204.6 ,1068.4 ,1045.7}; 
		    
  Double_t b[16] = {7.147 ,5.681 ,6.794 ,8.848 ,6.503 ,7.394 ,6.568 ,5.928,
  		    6.739 ,6.815 ,6.510 ,6.866 ,6.136 ,6.118 ,6.799 ,6.230}; 
		    
  TH2F *gains = new TH2F("ADPMGains", "AD PM gain factors", 16, -0.5, 15.5, 2, -0.5, 1.5);
  for(Int_t channel = 0; channel < 16; ++channel) {
    gains->SetBinContent(channel+1,1,a[channel]);
    gains->SetBinContent(channel+1,2,b[channel]);
  }
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("PM gain factors channel by channel");
  md->AddDateToComment();
  md->PrintMetaData();

  AliCDBId id("AD/Calib/PMGains", 0, AliCDBRunRange::Infinity());
  man->Put(gains, id, md);

  delete md;

}
