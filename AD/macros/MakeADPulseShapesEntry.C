void MakeADPulseShapesEntry(const char *outputCDB = "local://$ALICE_ROOT/OCDB")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  // Creation of the PM gains OCDB object
  Double_t offset[16] = {6,6,6,6,6,6,6,6,
  		         6,6,6,6,6,6,6,6};
		    
  Double_t tau[16] = {21.3, 24.5, 24.2, 24.3, 21.6, 26.7, 23.8, 24.1,
		      23.5, 24.1, 22.2, 25.8, 24.9, 23.0, 25.0, 23.9};
		    
  Double_t sigma[16] = {0.57, 0.44, 0.44, 0.43, 0.54, 0.41, 0.44, 0.43, 
		        0.49, 0.43, 0.53, 0.42, 0.45, 0.49, 0.43, 0.45};
		  

  TH2F *pulseParams = new TH2F("ADPulseShapes", "AD Pulse Shape parameters", 16, -0.5, 15.5, 3, -0.5, 2.5);
  for(Int_t channel = 0; channel < 16; ++channel) {
    pulseParams->SetBinContent(channel+1,1,offset[channel]);
    pulseParams->SetBinContent(channel+1,2,tau[channel]);
    pulseParams->SetBinContent(channel+1,3,sigma[channel]);

  }
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("AD Pulse Shape parameters");
  md->AddDateToComment();
  md->PrintMetaData();

  AliCDBId id("AD/Calib/PulseShapes", 0, AliCDBRunRange::Infinity());
  man->Put(pulseParams, id, md);

  delete md;

}
