void MakeOCDBTempTree()
{
  AliCDBEntry *fCDBEntry = new AliCDBEntry();
  TEnv *fConfEnv = new TEnv();

  // variables for tree
  Int_t Sensor=0; // 0 to 169
  Int_t Sec=0; // 0 to 5
  char  Side[10]; // A or C
  Int_t Num=0; // 0 to 8, per SM
  Int_t ECha=0; // DCS id; same as Sensor (I hope)..
  const char *sideStr[] = {"A", "C"};

  TTree *treeT = new TTree("treeT","treeT");
  treeT->Branch("Sensor",&Sensor,"Sensor/I");
  treeT->Branch("Side",Side,"Side/C");
  treeT->Branch("Sec",&Sec,"Sec/I");
  treeT->Branch("Num",&Num,"Num/I");
  treeT->Branch("ECha",&ECha,"ECha/I");

  for (int iSM=0; iSM<20; iSM++) { // SuperModule
    Sec = iSM/2;
    int iside = iSM%2; 
    sprintf(Side,"%s",sideStr[iside]);

    for (int ip=0; ip<8; ip++) { // points/sensors per SM
      Num = ip;
      Sensor = iSM*8 + ip; // Id
      ECha = Sensor; // IdDCS
      if (Sec != 5) { 
	treeT->Fill();
      }
      else { // SMA5 and SMC5 only has 4 sensors each
	if (ip<4) {
	  treeT->Fill();
	}
      }
    }
  }

  fCDBEntry->SetObject(treeT);

  // done; now save it..; add some metadata business etc.
  Int_t firstRun   =  217000; // a runno from between Run1 and Run2 (LS1)
  Int_t lastRun    =  999999999;
  Int_t version = 0;
  Int_t subversion = 0;
  Int_t beamPeriod =  1;

  char filename[200];
  sprintf(filename, "Temperature/Run%d_%d_v%d_s%d.root",
	  firstRun, lastRun, version, subversion);

  AliCDBMetaData md;
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("David Silvermyr");
  
  AliCDBId id("EMCAL/Config/Temperature", firstRun, lastRun, version, subversion);

  fCDBEntry->SetId(id);
  fCDBEntry->SetMetaData(&md);

  // ok, write the file
  TFile f(filename, "recreate");
  if (!f.IsZombie()) {
    f.cd();
    fCDBEntry->Write("AliCDBEntry");
    f.Close();
  }

}
