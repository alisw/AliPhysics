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
  /*/
  Double_t a[16] = {1202.8 ,1083.9 ,1140.7 ,1194.4 ,1135.3 ,1110.0 ,1065.3 ,1102.6,
	   	    1104.9 ,1170.4 ,1210.2 ,1061.4 ,1117.9 ,1204.6 ,1068.4 ,1045.7}; 
		    
  Double_t b[16] = {7.147 ,5.681 ,6.794 ,8.848 ,6.503 ,7.394 ,6.568 ,5.928,
  		    6.739 ,6.815 ,6.510 ,6.866 ,6.136 ,6.118 ,6.799 ,6.230};
   
		    
  //2015 after run  233912 (LHC15h period).
  Double_t a[16] = {1518.0, 1522.1, 1565.3, 1630.1, 1517.8, 1465.8, 1164.0, 1512.2, 
  	            1482.2, 1606.7, 1824.3, 1420.4, 1508.8, 1678.5, 1453.7, 1422.9};
		    
		    
  Double_t b[16] = {7.664, 5.681, 6.794, 8.848, 6.503, 7.394, 5.609, 5.928, 
                    6.739, 6.815, 6.510, 6.866, 6.136, 6.118, 6.799, 6.230};
		    
  //2015 after run 238785 		    
  Double_t a[16] = {1272.8, 1283.7, 1325.5, 1511.5, 1498.4, 1449.6, 1349.2, 1337.8, 
  	            1348.3, 1476.9, 1243.3, 1403.6, 1337.3, 1347.0, 1376.4, 1487.0};
  
  Double_t b[16] = { 6.230, 6.866, 6.799, 6.510, 6.118, 6.815, 6.136, 6.739,
  		     5.928, 8.848, 6.568, 6.794, 7.394, 5.681, 6.503, 7.147};

   		    
  //2015 after run 245683
  Double_t a[16] = {1272.8, 1283.7, 1325.5, 1511.5, 1498.4, 1449.6, 1104.4, 1337.8, 
  	            1348.3, 1476.9, 1243.3, 1403.6, 1337.3, 1347.0, 1376.4, 1487.0};
  
  Double_t b[16] = { 6.230, 6.866, 6.799, 6.510, 6.118, 6.815, 9.027, 6.739,
  		     5.928, 8.848, 6.568, 6.794, 7.394, 5.681, 6.503, 7.147};
  /*/		     
  //2016 after run 252234		     
  Double_t a[16] = {1584.3, 1472.4, 1422.5, 1376.6, 1417.6, 1316.6, 1251.3, 1437.6,
		    1419.6, 1488.5, 1599.5, 1320.5, 1438.0, 1557.8, 1396.2, 1352.5};
		    
  Double_t b[16] = { 8.252, 7.496, 7.229, 7.232, 7.488, 7.284, 7.150, 7.613,
		     7.296, 8.003, 9.106, 7.477, 7.225, 7.254, 8.036, 7.600};  
		    
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

  AliCDBId id("AD/Calib/PMGains", 252234, AliCDBRunRange::Infinity());
  man->Put(gains, id, md);

  delete md;

}
