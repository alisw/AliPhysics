//

Int_t FindKrClusterCheck(const char *fileName="data.root");


Int_t FindKrClustersRaw(const char *fileName="data.root"){

  

  //
  // remove Altro warnings
  //
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);
  //
  // Get calibration
  //
  char *ocdbpath = gSystem->Getenv("OCDB_PATH");
  //char *ocdbpath ="local:///afs/cern.ch/alice/tpctest/OCDB";
  if (ocdbpath==0){
    ocdbpath="alien://folder=/alice/data/2007/LHC07w/OCDB/";
  }
  printf("OCDB PATH = %s\n",ocdbpath); 
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbpath);
  man->SetRun(100000000);

  AliTPCCalPad * noiseTPC = AliTPCcalibDB::Instance()->GetPadNoise();
  AliTPCAltroMapping** mapping =AliTPCcalibDB::Instance()->GetMapping();
  //
  //


  //define tree
  TFile *hfile=new TFile("adc.root","RECREATE","ADC file");
  // Create a ROOT Tree
  TTree *mytree = new TTree("Kr","Krypton cluster tree");


  Int_t debugLevel=1;
  if(debugLevel>0){
    TH1F *histoRow   =new TH1F("histoRow","rows",100,0.,100.);
    TH1F *histoPad   =new TH1F("histoPad","pads",150,0.,150.);
    TH1F *histoTime  =new TH1F("histoTime","timebins",100,0.,1000.);
    TH2F *histoRowPad=new TH2F("histoRowPad","pads-vs-rows" ,150,0.,150.,100,0.,100.);
  }


  //one general output
  AliTPCclustererKr *clusters = new AliTPCclustererKr();
  clusters->SetOutput(mytree);
  clusters->SetRecoParam(0);

  if(debugLevel>0){
    clusters->SetDebugLevel(debugLevel);
    clusters->SetHistoRow(histoRow   );
    clusters->SetHistoPad(histoPad   );
    clusters->SetHistoTime(histoTime  );
    clusters->SetHistoRowPad(histoRowPad);
  }


  AliTPCParamSR *param=new AliTPCParamSR();
  //only for geometry parameters loading - temporarly
//  AliRunLoader* rl = AliRunLoader::Open("galice.root");
//  AliTPCParam *param=(AliTPCParamSR *)gDirectory->Get("75x40_100x60_150x60");
  //if (!param) {cerr<<"TPC parameters have not been found !\n"; return 4;}

  clusters->SetParam(param);

  //set cluster finder parameters (from data)
  clusters->SetZeroSup(param->GetZeroSup());//zero suppression parameter
  clusters->SetFirstBin(60);//first bin
  clusters->SetLastBin(950);//last bin
  clusters->SetMaxNoiseAbs(2);//maximal noise
  clusters->SetMaxNoiseSigma(3);//maximal amount of sigma of noise

  //set cluster finder parameters (from MC)
  clusters->SetMinAdc(3);//signal threshold (everything below is treated as 0)
  clusters->SetMinTimeBins(2);//number of neighbouring timebins
  clusters->SetMaxPadRangeCm(5.);//distance of the cluster center to the center of a pad (in cm)
  clusters->SetMaxRowRangeCm(5.);//distance of the cluster center to the center of a padrow (in cm)
  clusters->SetMaxTimeRange(5.);//distance of the cluster center to the max time bin on a pad (in tackts)
  //ie. fabs(centerT - time)<7
  clusters->SetValueToSize(7.);//cut reduce peak at 0
  clusters->SetIsolCut(3);//set isolation cut threshold

  AliRawReader *reader = new AliRawReaderRoot(fileName);
  reader->Reset();

  TStopwatch timer;
  timer.Start();

  AliAltroRawStreamFast* stream = new AliAltroRawStreamFast(reader);
  stream->SelectRawData("TPC");

  Int_t evtnr=0;
  while (reader->NextEvent()) {
    //output for each event

    //if(evtnr>4) break;
    cout<<"Evt = "<<evtnr<<endl;
    clusters->FinderIO(reader);
    evtnr++;
    AliSysInfo::AddStamp(Form("Event%d",evtnr),evtnr);
  }


  mytree->Print();//print rootuple summary 
  // Save all objects in this file
  hfile->Write();
  // Close the file
  hfile->Close();

  timer.Stop();
  timer.Print();
  printf("Deleting clusterer\n");
  delete clusters;
  printf("Deleting stream\n");
  delete stream;
  printf("Deleting raw reader\n");
  delete reader;

//   TCanvas *c2=new TCanvas("c2","title",800,800);
//   c2->SetHighLightColor(2);
//   c2->Range(-458.9552,-2948.238,3296.642,26856.6);
//   c2->SetBorderSize(2);
//   c2->SetLeftMargin(0.15);
//   c2->SetRightMargin(0.06);
//   c2->SetFrameFillColor(0);

//   gStyle->SetOptStat(111111);
//   histoRow->Draw();
//   c2->Print("rows.ps");
//   histoPad->Draw();
//   c2->Print("pads.ps");
//   histoTime->Draw();
//   c2->Print("timebins.ps");
//   histoRowPad->Draw();
//   c2->Print("row-pad.ps");

  return 0;
}


Int_t FindKrClusterCheck(const char *fileName){
  //
  //
  gSystem->Load("$ROOTSYS/lib/libGui.so");
  gSystem->Load("$ROOTSYS/lib/libTree.so");
  gSystem->Load("$MEMSTAT/libMemStat.so");
  {
    TMemStat memstat(1000000,100000,kTRUE);
    AliSysInfo::AddCallBack(TMemStatManager::GetInstance()->fStampCallBack);
    FindKrClustersRaw(fileName); 
  }
  // the output memstat.root file
  TMemStat draw("memstat.root");
  // Print some information
  // code info
  draw.MakeReport(0,0);

}
