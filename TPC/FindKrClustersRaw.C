Int_t FindKrClustersRaw(const char *fileName="data.root"){



  //
  // remove Altro warnings
  //
  AliLog::SetClassDebugLevel("AliTPCRawStream",-5);
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);
  //
  // Get calibration
  //
  //  char *ocdbpath = gSystem->Getenv("OCDB_PATH");
  char *ocdbpath ="local:///afs/cern.ch/alice/tpctest/OCDB";
  if (ocdbpath==0){
    ocdbpath="alien://folder=/alice/data/2007/LHC07w/OCDB/";
  }
  printf("OCDB PATH = %s\n",ocdbpath); 
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbpath);
  man->SetRun(0);

  AliTPCCalPad * noiseTPC = AliTPCcalibDB::Instance()->GetPadNoise();
  AliTPCAltroMapping** mapping =AliTPCcalibDB::Instance()->GetMapping();
  //
  //


  //define tree
  TFile *hfile=new TFile("adc.root","RECREATE","ADC file");
  // Create a ROOT Tree
  TTree *mytree = new TTree("Kr","Krypton cluster tree");


  AliRawReader *reader = new AliRawReaderRoot(fileName);
  //AliRawReader *reader = new AliRawReaderDate(fileName);
  reader->Reset();

  TStopwatch timer;
  timer.Start();

  AliAltroRawStreamFast* stream = new AliAltroRawStreamFast(reader);
  stream->SelectRawData("TPC");

  //one general output
  AliTPCclustererKr *clusters = new AliTPCclustererKr();
  clusters->SetOutput(mytree);
  clusters->SetRecoParam(0);


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
  clusters->SetMaxPadRangeCm(2.5);//distance of the cluster center to the center of a pad (in cm)
  clusters->SetMaxRowRangeCm(3.5);//distance of the cluster center to the center of a padrow (in cm)
  clusters->SetMaxTimeRange(7);//distance of the cluster center to the max time bin on a pad (in tackts)
  //ie. fabs(centerT - time)<7
  clusters->SetValueToSize(3.5);//cut reduce peak at 0




  Int_t evtnr=0;
  while (reader->NextEvent()) {
    //output for each event
  //  AliTPCclustererKr *clusters = new AliTPCclustererKr();
  //  clusters->SetOutput(mytree);
  //  clusters->SetRecoParam(0);
  //  clusters->SetParam(param);

    // if(evtnr++>5) break;
    cout<<"Evt = "<<evtnr<<endl;
    clusters->FinderIO(reader);
    evtnr++;

  }

  mytree->Print();//print rootuple summary 
  // Save all objects in this file
  hfile->Write();
  // Close the file
  hfile->Close();

  timer.Stop();
  timer.Print();

  delete stream;

  return 0;
}
