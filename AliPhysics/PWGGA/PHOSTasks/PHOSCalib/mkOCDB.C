void mkOCDB(const Bool_t toGRID=kFALSE, const Bool_t debug=kFALSE)
{
  // This script reads PHOS-EMC calibration object from OCDB raw://,
  // then reads corrections to gains O(1) from local://OCDB_corr
  // and create the new OCDB object with the new gains as a product of old gains
  // and gain corrections.
  // The new calibration object is written either to local://OCDB_new
  // --
  // Gain corrections were found by equalizing the mean cell energy
  // from the reconstructed data LHC10h by B.Polishchuk
  // --
  // Yuri Kharlov. 7-Feb-2011
  gROOT->cd() ;

  AliCDBManager::Instance()->SetDefaultStorage("raw://");

  AliCDBManager::Instance()->SetRun(182724); //LHC12c
  AliPHOSCalibData db1(182724);

  TFile fgains("Calibration_pass1.root") ;
  char key[55] ;
  Int_t module=3;
  sprintf(key,"Mass_mod%d",module) ;
  TH2D * h=(TH2D*)fgains.Get(key) ;
 
  
  char key[55] ;
  for(Int_t module=3; module<4; module++) {

    for(Int_t ix=48; ix<64; ix++) {
      for(Int_t iz=28; iz<56; iz++) {

	Float_t cc_i = db1.GetADCchannelEmc(module+1,iz+1,ix+1);
	Float_t corr = h->GetBinContent(ix+1,iz+1) ;
	printf("ix=%d, iz=%d, corr=%f \n",ix,iz,corr) ;
	  printf("    old=%f \n",db1.GetADCchannelEmc(module+1,iz+1,ix+1));
        if(corr>0){
   	  db1.SetADCchannelEmc(module+1,iz+1,ix+1,cc_i*corr);
	  printf("    new=%f \n",db1.GetADCchannelEmc(module+1,iz+1,ix+1));
	}
      }
    }
  }
  
  //Writing new calibration coefficients (CCs) to OCDB

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://OCDB");

  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Dmitri Peressounko");
  md->SetComment("PHOS gains calculated from pi0 peak position (LHC12c). Corrected for HV adjustment. Widht 3.8-4.8 MeV");
  AliCDBId id("PHOS/Calib/EmcGainPedestals",177295,AliCDBRunRange::Infinity());
//   cdb->Put(&db1,id, md);
  db1.WriteEmc(177295,AliCDBRunRange::Infinity(),md);


  
}
