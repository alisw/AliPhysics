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

  TH1F *hGainOrigM1 = new TH1F("hGainOrigM1","hGainM1OrigM1",1000,0.,0.01);
  TH1F *hGainOrigM2 = new TH1F("hGainOrigM2","hGainM1OrigM2",1000,0.,0.01);
  TH1F *hGainOrigM3 = new TH1F("hGainOrigM3","hGainM1OrigM3",1000,0.,0.01);
  TH1F *hGainOrigM4 = new TH1F("hGainOrigM4","hGainM1OrigM4",1000,0.,0.01);
  TH1F *hGainOrigM5 = new TH1F("hGainOrigM5","hGainM1OrigM5",1000,0.,0.01);

  TH2F *hOrigNewM1 = new TH2F("hOrigNewM1","hOrigNewM1",100,0.,0.01,100,0.,0.01);
  TH2F *hOrigNewM2 = new TH2F("hOrigNewM2","hOrigNewM2",100,0.,0.01,100,0.,0.01);
  TH2F *hOrigNewM3 = new TH2F("hOrigNewM3","hOrigNewM3",100,0.,0.01,100,0.,0.01);
  TH2F *hOrigNewM4 = new TH2F("hOrigNewM4","hOrigNewM4",100,0.,0.01,100,0.,0.01);
  TH2F *hOrigNewM5 = new TH2F("hOrigNewM5","hOrigNewM5",100,0.,0.01,100,0.,0.01);

  TH1F *hGainNewM1 = new TH1F("hGainNewM1","hGainM1NewM1",1000,0.,0.01);
  TH1F *hGainNewM2 = new TH1F("hGainNewM2","hGainM1NewM2",1000,0.,0.01);
  TH1F *hGainNewM3 = new TH1F("hGainNewM3","hGainM1NewM3",1000,0.,0.01);
  TH1F *hGainNewM4 = new TH1F("hGainNewM4","hGainM1NewM4",1000,0.,0.01);
  TH1F *hGainNewM5 = new TH1F("hGainNewM5","hGainM1NewM5",1000,0.,0.01);

  AliCDBManager::Instance()->SetDefaultStorage("raw://");

  AliCDBManager::Instance()->SetRun(114783);
  AliPHOSCalibData db1(114783);

  AliCDBManager::Instance()->SetDefaultStorage("local://OCDB_corr");
  AliCDBManager::Instance()->SetSpecificStorage("PHOS/Calib/EmcGainPedestals",
						"local://OCDB_corr");
  AliCDBManager::Instance()->SetRun(114783);
  AliPHOSCalibData db2(114783);  

  for(Int_t module=0; module<5; module++) {
    for(Int_t ix=0; ix<64; ix++) {
      for(Int_t iz=0; iz<56; iz++) {

	Float_t cc_i = db1.GetADCchannelEmc(module+1,iz+1,ix+1);
	Float_t corr = db2.GetADCchannelEmc(module+1,iz+1,ix+1);
        if (module <= 2)
          corr *= 0.135/0.1223;

	if      (module == 0) {
	  hGainOrigM1->Fill(cc_i);
	  hGainNewM1 ->Fill(cc_i*corr);
	  hOrigNewM1 ->Fill(cc_i,cc_i*corr);
	}
	else if (module == 1) {
	  hGainOrigM2->Fill(cc_i);
	  hGainNewM2 ->Fill(cc_i*corr);
	  hOrigNewM2 ->Fill(cc_i,cc_i*corr);
	}
	else if (module == 2) {
	  hGainOrigM3->Fill(cc_i);
	  hGainNewM3 ->Fill(cc_i*corr);
	  hOrigNewM3 ->Fill(cc_i,cc_i*corr);
	}
	else if (module == 3) {
	  hGainOrigM4->Fill(cc_i);
	  hGainNewM4 ->Fill(cc_i*corr);
	  hOrigNewM4 ->Fill(cc_i,cc_i*corr);
	}
	else if (module == 4) {
	  hGainOrigM5->Fill(cc_i);
	  hGainNewM5 ->Fill(cc_i*corr);
	  hOrigNewM5 ->Fill(cc_i,cc_i*corr);
	}

	db1.SetADCchannelEmc(module+1,iz+1,ix+1,cc_i*corr);
	if (debug) {
	  if (TMath::Abs(corr-1.)>0.01)
	    printf("CC for mod%d x%d z%d set to %f (orig %f, corr %f)\n",
		   module+1,ix+1,iz+1,
		   db1.GetADCchannelEmc(module+1,iz+1,ix+1),
		   cc_i,corr);
	}

      }
    }
  }
  
  //Writing new calibration coefficients (CCs) to OCDB

  AliCDBManager* cdb = AliCDBManager::Instance();
  if (toGRID)
    cdb->SetDefaultStorage("alien://Folder=/alice/cern.ch/user/k/kharlov/RAW/2010/OCDB");
  else
    cdb->SetDefaultStorage("local://OCDB_new");

  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Yuri Kharlov / Boris Polishchuk");
  md->SetComment("PHOS gains calculated from mean cell energy equalization. All scaled to shift pi0 to 135 MeV. HL/LG ratio is calculated from LED runs offline. The pi0 width is 5-7 MeV");
  AliCDBId id("PHOS/Calib/EmcGainPedestals",114783,AliCDBRunRange::Infinity());
//   cdb->Put(&db1,id, md);
  db1.WriteEmc(136833,AliCDBRunRange::Infinity(),md);

  TFile *gainHisto = new TFile("gainHisto.root","recreate");
  hGainOrigM1->Write();
  hGainOrigM2->Write();
  hGainOrigM3->Write();
  hGainOrigM4->Write();
  hGainOrigM5->Write();
  hGainNewM1->Write();
  hGainNewM2->Write();
  hGainNewM3->Write();
  hGainNewM4->Write();
  hGainNewM5->Write();
  hOrigNewM1->Write();
  hOrigNewM2->Write();
  hOrigNewM3->Write();
  hOrigNewM4->Write();
  hOrigNewM5->Write();
  gainHisto->Close();
}
