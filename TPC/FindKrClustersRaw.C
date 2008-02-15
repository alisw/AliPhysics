Int_t FindKrClustersRaw(const char *fileName="data.root"){

  //define tree
  TFile *hfile=new TFile("KryptonCl.root","RECREATE","ADC file");
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

  //only for geometry parameters loading - temporarly
  //AliRunLoader* rl = AliRunLoader::Open("galice.root");
  //AliTPCParam *param=(AliTPCParamSR *)gDirectory->Get("75x40_100x60_150x60");
  AliTPCParam *param=new AliTPCParam;
  param->SetTSample(1.00000002337219485e-07);
  param->Update();
  //if (!param) {cerr<<"TPC parameters have not been found !\n"; return 4;}
  clusters->SetParam(param);
  

  Int_t evtnr=0;
  while (reader->NextEvent()) {
    //output for each event
    cout<<"Evt = "<<evtnr<<endl;
    clusters->finderIO(reader);
    evtnr++;
  }

  //mytree->Print();//print rootuple summary 
  // Save all objects in this file
  hfile->Write();
  // Close the file
  hfile->Close();

  timer.Stop();
  timer.Print();

  delete stream;

  return 0;
}
