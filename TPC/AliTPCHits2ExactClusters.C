/****************************************************************************
 *           Origin: M.Ivanov, CERN, Marian.Ivanov@cern.ch                 *
 ****************************************************************************/

Int_t AliTPCHits2ExactClusters(Int_t nevent=1)
{
  //
  // loop over hits and find itersection of the trak with pad rows
  // so called exact clusters stored in the structure similar to hits
  // it's used for comparison purpose - to have benchmark of cluster
  // finder and tracker


  // Connect the Root Galice file containing Geometry, Kine and Hits

  const char * inFile_old = "galice.root"; 
  const char * inFile_new = "galice.root";
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile_old);
  if (file) {file->Close(); delete file;}
  file =  TFile::Open(inFile_new,"UPDATE");
  if (!file->IsOpen()) {
    cerr<<"Can't open "<<inFile_new<<" !\n";
    return 1;
  }

  // Get AliRun object from file or create it if not on file
  if (gAlice) delete gAlice;
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
    cerr<<"AliTPCHits2Digits.C : AliRun object not found on file\n";
    return 2;
  }



  // gAlice->GetEvent(0);
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");      
  AliTPCParam * param = TPC->GetParam();
  TStopwatch timer;
  timer.Start();

  for(Int_t eventn =0;eventn<nevent;eventn++){
    printf("Processing event %d\n",eventn);
    gAlice->GetEvent(eventn);
    //setup AliTPCClustersArray
    AliTPCClustersArray * arr=new AliTPCClustersArray;
    arr->SetClusterType("AliComplexCluster");
    arr->Setup(param);
    arr->MakeTree();    
    TPC->SetClustersArray(arr);    
    for (Int_t i=0;i<72;i++)
      TPC->Hits2ExactClustersSector(i);
    char treeName[100];
    sprintf(treeName,"TreeCExact_%d",eventn);
    TPC->GetClustersArray()->GetTree()->Write(treeName);
  }

  delete gAlice; gAlice=0;
  file->Close(); delete file;
  timer.Stop();
  timer.Print();

  return 0;
};




