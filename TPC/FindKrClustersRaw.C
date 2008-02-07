Int_t FindKrClustersRaw(const char *fileName="data.root"){

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

  //only for geometry parameters loading - temporarly
  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  AliTPCParam *param=(AliTPCParamSR *)gDirectory->Get("75x40_100x60_150x60");
  //if (!param) {cerr<<"TPC parameters have not been found !\n"; return 4;}
  clusters->SetParam(param);


  Int_t evtnr=0;
  while (reader->NextEvent()) {
    //output for each event
    //  AliTPCclustererKr *clusters = new AliTPCclustererKr();
    //  clusters->SetOutput(mytree);

    //if(evtnr++<35)continue;
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
