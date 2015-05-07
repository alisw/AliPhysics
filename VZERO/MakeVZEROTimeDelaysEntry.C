void MakeVZEROTimeDelaysEntry()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // Creation of the time delays OCDB object
  const Double_t timeShift[66] = {0.0 , 0.477957 , 0.0889999 , 0.757669 , 0.205439 , 0.239666 , -0.183705 , 0.442873 , -0.281366 , 0.260976 , 0.788995 , 0.974758 , 0.548532 , 0.495023 , 0.868472 , 0.661167 , 0.358307 , 0.221243 , 0.530179 , 1.26696 , 1.33082 , 1.27086 , 1.77133 , 1.10253 , 0.634806+0.885 , 2.14838 , 1.50212 , 1.59253 , 1.66122+0.740 , 1.16957 , 1.52056 , 1.47791 , 1.81905 , -1.94123 , -1.29124-0.350 , -2.16045 , -1.78939 , -3.11111 , -1.87178 , -1.57671-0.560 , -1.70311 , -1.81208 , -1.94475 , -2.53058+0.930 , -1.7042 , -2.08109 , -1.84416 , -0.61073 , -1.77145 , 0.16999 , -0.0585339 , 0.00401133 , 0.397726 , 0.851111 , 0.264187 , 0.59573 , -0.158263 , 0.584362 , 1.20835 , 0.927573 , 1.13895 , 0.64648 , 2.18747 , 1.68909 , 0.451194 , 0.0};
  TH1F *delays = new TH1F("VZEROTimeDelays","VZERO Time delays",64,-0.5,63.5);
  delays->SetContent(timeShift);
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Time delays channel by channel for year >= 2012");
  md->PrintMetaData();

  AliCDBId id("VZERO/Calib/TimeDelays",0,AliCDBRunRange::Infinity());

  man->Put(delays, id, md);

  delete md;

}
