//#include "alles.h"

void TPCHits2ExactClusters()
{
  if (gTPCParam==0) {
    cout<<"TPC Param Not Initialised\n";
    return;
  }
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }  
  // Connect the Root Galice file containing Geometry, Kine and Hits
  const char * inFile = "galice.root";  
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
  if (file) file->Close();
  file = new TFile(inFile,"UPDATE");
  // Get AliRun object from file or create it if not on file
  //  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    //  };  
  gAlice->GetEvent(0);
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");      

  //setup AliTPCClustersArray
  AliTPCClustersArray * arr=new AliTPCClustersArray;
  arr->SetClusterType("AliCluster");
  arr->Setup(gTPCParam);
  TPC->SetParam(gTPCParam);
  arr->MakeTree();

  TPC->SetClustersArray(arr); 
  TPC->Hits2ExactClustersSector(0);  
  //write results
  char treeName[100];
  sprintf(treeName,"TreeCExact_%s",gTPCParam->GetTitle());
  TPC->GetClustersArray()->GetTree()->Write(treeName);
}


AliTPCClustersArray *  GetExactClustersArray(Bool_t newtree=kFALSE, const char* name=0 ) 
{
  AliTPCClustersArray * arr;
  if ( (gAlice!=0) && (gAlice->GetDetector("TPC")!=0) ) {    
    arr = ((AliTPC*)gAlice->GetDetector("TPC"))->GetClustersArray();
    if (arr!=0) arr->Update();
  }
  if (arr==0) {    
    arr = new AliTPCClustersArray;
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    if (file==0) file = new TFile("galice.root","update");   
    arr->SetClusterType("AliCluster");         
    arr->Setup(gTPCParam);
    cout<<"Update status : "<<arr->Update()<<"\n";
    char treeName[100];
    if (name==0)
      sprintf(treeName,"TreeCExact_%s",gTPCParam->GetTitle());
    else sprintf(treeName,"TreeCExact_%s",name); 
    if (newtree!=kTRUE){
      cout<<"Connect tree status : "<<arr->ConnectTree(treeName)<<"\n";
    }
    else 
      arr->MakeTree();
  }
  return arr;
}



     
