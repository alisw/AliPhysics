void timeStamp(const char *listfile = "temperatureValuesPerSMPerRun.txt")
  
{

  Int_t nruns = 908;
  const Int_t useRuns=nruns;
  Int_t currentRunNumber;
  Int_t row=0;
  Int_t count=-1;
  const Int_t kNSM=20;
  char id[100];
  char title[100];

  ifstream fin(listfile);
  Int_t runno = 0;
  Int_t iSM;
  Float_t minTempAvg;

  TH1F *hTemperature[kNSM];
  for (int iSM=0; iSM<kNSM; iSM++) {
    sprintf(id, "hTemperatureSM_%d", iSM);
    sprintf(title, "avg min Temp: SM %d ", iSM);
    hTemperature[iSM]= new TH1F(id, title, nruns, -0.5, useRuns-0.5);
  }

  TH1F *hRunNumber=new TH1F("hRunNumber", " ", nruns, -0.5, useRuns-0.5);

  while ( fin.good() ) {
    fin >> runno >> iSM >> minTempAvg;
//    cout<<runno<< iSM  <<minTempAvg<<endl;
    //row=row+1;
    if (row%20==0) {
      count=count+1;
      hRunNumber->Fill(count, runno);
    }
    hTemperature[iSM]->Fill(count, minTempAvg);
    
    row=row+1;
  }

  TFile destFile("tempDependenceAllRuns.root", "recreate");
  for (int iSM=0; iSM<kNSM; iSM++) {
    hTemperature[iSM]->Write();
  }
  hRunNumber->Write();
}
