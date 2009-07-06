// Macro to merge TPC performance train output.
// One can merge all performance objects stored in the 
// output files or one of them:
// 
// RES 
// EFF
// DEDX
// DCA
// TPC
//
// It is recommeneded to merge particular components 
// for a big statistcs due to too big size of the output. For small
// dedicated productions (e.g. flat-p) all components can be merged at once.
//

void MergeOutput(const char* filename = "list.txt", const char* comp = "ALL", Int_t max=2)
{
  // Performance object 
  TString str(comp);

  // Create new collections
  TList *collTPC = new TList;
  TList *collRes = new TList;
  TList *collResTPCInner = new TList;
  TList *collEff = new TList;
  TList *colldEdxTPCInner = new TList;
  TList *collDCA = new TList;

  // Open the input stream
  ifstream in;
  in.open(filename);

  // Read the input list of files and add them to the chain
  TString perfile;
  Int_t counter=0;
  while(!in.eof()) {
    in >> perfile;
    if (!perfile.Contains("root")) continue; // protection
    counter++;
    if(counter>max) break;
    printf("%s \n",perfile.Data());
    TFile::Open(perfile.Data());
    gFile->cd();
    TList *coutput = gFile->Get("coutput"); 

    if(str.CompareTo("ALL")==0) {
      collTPC->Add(coutput->FindObject("AliPerformanceTPC"));
      collRes->Add(coutput->FindObject("AliPerformanceRes"));
      collResTPCInner->Add(coutput->FindObject("AliPerformanceResTPCInner"));
      collEff->Add(coutput->FindObject("AliPerformanceEff"));
      colldEdxTPCInner->Add(coutput->FindObject("AliPerformanceDEdxTPCInner"));
      collDCA->Add(coutput->FindObject("AliPerformanceDCA"));
    } 
    else if(str.CompareTo("TPC")==0) {
      collTPC->Add(coutput->FindObject("AliPerformanceTPC"));
    }
    else if(str.CompareTo("RES")==0) {
      collRes->Add(coutput->FindObject("AliPerformanceRes"));
      collResTPCInner->Add(coutput->FindObject("AliPerformanceResTPCInner"));
    }
    else if(str.CompareTo("EFF")==0) {
      collEff->Add(coutput->FindObject("AliPerformanceEff"));
    }
    else if(str.CompareTo("DEDX")==0) {
      colldEdxTPCInner->Add(coutput->FindObject("AliPerformanceDEdxTPCInner"));
    }
    else if(str.CompareTo("DCA")==0) {
      collDCA->Add(coutput->FindObject("AliPerformanceDCA"));
    }
    else {
     continue;
    }
  }
  in.close();

  //
  AliPerformanceTPC *tpc = new AliPerformanceTPC("AliPerformanceTPC","AliPerformanceTPC",0,kFALSE);
  AliPerformanceRes *res = new AliPerformanceRes("AliPerformanceRes","res",0,kFALSE);
  AliPerformanceRes *resTPCInner = new AliPerformanceRes("AliPerformanceResTPCInner","AliPerformanceResTPCInner",3,kFALSE);
  AliPerformanceEff *eff = new AliPerformanceEff("AliPerformanceEff","AliPerformanceEff",0,kFALSE);
  AliPerformanceDEdx *dedx = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner","AliPerformanceDEdxTPCInner",3,kFALSE);
  AliPerformanceDCA *dca = new AliPerformanceDCA("AliPerformanceDCA","AliPerformanceDCA",0,kFALSE);

  // Create output list
  TFolder *outList = new TFolder();
  outList->SetOwner();
  outList->SetName("coutput");

  // merge objects and add to the list
  if(!str.CompareTo("ALL")) {
    tpc->Merge(collTPC); 
    outList->Add(tpc);
    res->Merge(collRes);
    outList->Add(res);
    resTPCInner->Merge(collResTPCInner);
    outList->Add(resTPCInner);
    eff->Merge(collEff);
    outList->Add(eff);
    dedx->Merge(colldEdxTPCInner);
    outList->Add(dedx);
    dca->Merge(collDCA);
    outList->Add(dca);
  }
  else if(!str.CompareTo("TPC")) {
    tpc->Merge(collTPC); 
    outList->Add(tpc);
  }
  else if(!str.CompareTo("RES")) {
    res->Merge(collRes);
    outList->Add(res);
    resTPCInner->Merge(collResTPCInner);
    outList->Add(resTPCInner);
  }
  else if(!str.CompareTo("EFF")) {
    eff->Merge(collEff);
    outList->Add(eff);
  }
  else if(!str.CompareTo("DEDX")) {
    dedx->Merge(colldEdxTPCInner);
    outList->Add(dedx);
  }
  else if(!str.CompareTo("DCA")) {
    dca->Merge(collDCA);
    outList->Add(dca);
  }

  //
  TFile *outFile = new TFile("TPC.Performance.Merged.root","recreate");
  outFile->cd();
  outList->Write();
  outFile->Close();
}
