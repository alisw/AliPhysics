// Macro to make final QA spectra from information stored
// in QA performance objects. It is done in 2 steps:
// 1. Merging TPC QA performance train output
// 2. Analyse generic histograms to make final QA spectra 
//
// 
// One can merge all performance objects stored in the 
// output files or one of them:
// 
// ALL - all components only for simulaltion
// NO_MC - all components only for realdata
// RES (At DCA, at InnerTPC, at outerTPC)
// EFF
// DEDX
// DCA
// TPC
// Match (TPC+ITS, TPC+TRD, TPC Eff wrt ITS)
//
// It is recommeneded to merge particular components 
// for a big statistcs due to too big size of the output. For small
// dedicated productions (e.g. flat-p) all components can be merged at once.
//
//void MergeOutput(const char* filename = "list.txt", const char* comp = "ALL", Int_t max=2);
//void MakeFinalSpectra(const char* filename = "merged.root", const char* comp = "ALL");
void MergeOutput(const char* filename, const char* comp, Int_t max);
//void MakeFinalSpectra(const char* filename, const char* comp);
void MakeFinalSpectra(const char* comp);

void MakeResults(const char* filename = "list.txt", const char* comp = "NO_MC", Int_t max=-1)
{
  // Merge performance objects
  MergeOutput(filename, comp, max);

  // Make final spectra
  MakeFinalSpectra(comp);
}


void MergeOutput(const char* filename, const char* comp, Int_t max)
{
  //
  if(!filename) 
  { 
    printf("ERROR file %s does not exist \n",filename);
    return;
  }

  // Performance object 
  TString str(comp);

  // Create new collections
  TList *collTPC = new TList;
  TList *collRes = new TList;
  TList *collResTPCInner = new TList;
  TList *collResTPCOuter = new TList;
  TList *collEff = new TList;
  TList *colldEdxTPCInner = new TList;
  TList *collDCA = new TList;
  TList *collMatchTPCITS = new TList;
  TList *collMatchTPCTRD = new TList;
  TList *collMatchTPCEFF = new TList;

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
    if(max>0 && counter>max) break;
    printf("%s \n",perfile.Data());
    TFile::Open(perfile.Data());
    gFile->cd();
    TList *coutput = gFile->Get("TPC"); 

    if(str.CompareTo("NO_MC")==0) {
      collTPC->Add(coutput->FindObject("AliPerformanceTPC"));
      colldEdxTPCInner->Add(coutput->FindObject("AliPerformanceDEdxTPCInner"));
      collDCA->Add(coutput->FindObject("AliPerformanceDCA"));
      collMatchTPCITS->Add(coutput->FindObject("AliPerformanceMatchTPCITS"));
      collMatchTPCTRD->Add(coutput->FindObject("AliPerformanceMatchTPCTRD"));
      collMatchTPCEFF->Add(coutput->FindObject("AliPerformanceMatchTPCEFF"));
    } 

    if(str.CompareTo("ALL")==0) {
      collTPC->Add(coutput->FindObject("AliPerformanceTPC"));
      collRes->Add(coutput->FindObject("AliPerformanceRes"));
      collResTPCInner->Add(coutput->FindObject("AliPerformanceResTPCInner"));
      collResTPCOuter->Add(coutput->FindObject("AliPerformanceResTPCOuter"));
      collEff->Add(coutput->FindObject("AliPerformanceEff"));
      colldEdxTPCInner->Add(coutput->FindObject("AliPerformanceDEdxTPCInner"));
      collDCA->Add(coutput->FindObject("AliPerformanceDCA"));
      collMatchTPCITS->Add(coutput->FindObject("AliPerformanceMatchTPCITS"));
      collMatchTPCTRD->Add(coutput->FindObject("AliPerformanceMatchTPCTRD"));
      collMatchTPCEFF->Add(coutput->FindObject("AliPerformanceMatchTPCEFF"));
    } 
    else if(str.CompareTo("TPC")==0) {
      collTPC->Add(coutput->FindObject("AliPerformanceTPC"));
    }
    else if(str.CompareTo("RES")==0) {
      collRes->Add(coutput->FindObject("AliPerformanceRes"));
      collResTPCInner->Add(coutput->FindObject("AliPerformanceResTPCInner"));
      collResTPCOuter->Add(coutput->FindObject("AliPerformanceResTPCOuter"));
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
    else if(str.CompareTo("MATCH")==0) {
      collMatchTPCITS->Add(coutput->FindObject("AliPerformanceMatchTPCITS"));
      collMatchTPCTRD->Add(coutput->FindObject("AliPerformanceMatchTPCTRD"));
      collMatchTPCEFF->Add(coutput->FindObject("AliPerformanceMatchTPCEFF"));
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
  AliPerformanceRes *resTPCOuter = new AliPerformanceRes("AliPerformanceResTPCOuter","AliPerformanceResTPCOuter",4,kFALSE);
  AliPerformanceEff *eff = new AliPerformanceEff("AliPerformanceEff","AliPerformanceEff",0,kFALSE);
  AliPerformanceDEdx *dedx = new AliPerformanceDEdx("AliPerformanceDEdxTPCInner","AliPerformanceDEdxTPCInner",3,kFALSE);
  AliPerformanceDCA *dca = new AliPerformanceDCA("AliPerformanceDCA","AliPerformanceDCA",0,kFALSE);
  AliPerformanceMatch *matchTPCITS = new AliPerformanceMatch("AliPerformanceMatchTPCITS","AliPerformanceTPCITS",0,kFALSE);
  AliPerformanceMatch *matchTPCTRD = new AliPerformanceMatch("AliPerformanceMatchTPCTRD","AliPerformanceTPCTRD",1,kFALSE);
  AliPerformanceMatch *matchTPCEFF = new AliPerformanceMatch("AliPerformanceMatchTPCEFF","AliPerformanceTPCEFF",2,kFALSE);

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
    resTPCOuter->Merge(collResTPCOuter);
    outList->Add(resTPCOuter);
    eff->Merge(collEff);
    outList->Add(eff);
    dedx->Merge(colldEdxTPCInner);
    outList->Add(dedx);
    dca->Merge(collDCA);
    outList->Add(dca);
    matchTPCITS->Merge(collMatchTPCITS);
    outList->Add(matchTPCITS);
    matchTPCTRD->Merge(collMatchTPCTRD);
    outList->Add(matchTPCTRD);
    matchTPCEFF->Merge(collMatchTPCEFF);
    outList->Add(matchTPCEFF);
  }
  else if(!str.CompareTo("NO_MC")) {
    tpc->Merge(collTPC); 
    outList->Add(tpc);
    dedx->Merge(colldEdxTPCInner);
    outList->Add(dedx);
    dca->Merge(collDCA);
    outList->Add(dca);
    matchTPCITS->Merge(collMatchTPCITS);
    outList->Add(matchTPCITS);
    matchTPCTRD->Merge(collMatchTPCTRD);
    outList->Add(matchTPCTRD);
    matchTPCEFF->Merge(collMatchTPCEFF);
    outList->Add(matchTPCEFF);
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
    resTPCOuter->Merge(collResTPCOuter);
    outList->Add(resTPCOuter);
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
  else if(!str.CompareTo("MATCH")) {
    matchTPCITS->Merge(collMatchTPCITS);
    outList->Add(matchTPCITS);
    matchTPCTRD->Merge(collMatchTPCTRD);
    outList->Add(matchTPCTRD);
    matchTPCEFF->Merge(collMatchTPCEFF);
    outList->Add(matchTPCEFF);
  }

  //
  
  TFile *outFile = new TFile( Form("TPC.Performance.Merged.%s.root",comp),"recreate");
  outFile->cd();
  outList->Write();
  outFile->Close();
}

void MakeFinalSpectra(const char* comp) 
{
  // open proper input file
  TFile *inFile = TFile::Open(Form("TPC.Performance.Merged.%s.root",comp));
  inFile->cd();

  // Performance object 
  TString str(comp);

  if(str.CompareTo("ALL") == 0) 
  {
    AliPerformanceTPC * compObjTPC = (AliPerformanceTPC*)coutput->FindObject("AliPerformanceTPC");
    compObjTPC->Analyse();
    compObjTPC->GetAnalysisFolder()->ls("*");
    compObjTPC->PrintHisto(kTRUE,"PerformanceTPCQA.ps");
    TFile fout("PerformanceTPCQA.root","recreate");
    compObjTPC->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceEff * compObjEff = (AliPerformanceEff*)coutput->FindObject("AliPerformanceEff");
    compObjEff->Analyse();
    compObjEff->GetAnalysisFolder()->ls("*");
    compObjEff->PrintHisto(kTRUE,"PerformanceEffQA.ps");
    TFile fout("PerformanceEffQA.root","recreate");
    compObjEff->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceDCA * compObjDCA = (AliPerformanceDCA*)coutput->FindObject("AliPerformanceDCA");
    compObjDCA->Analyse();
    compObjDCA->GetAnalysisFolder()->ls("*");
    compObjDCA->PrintHisto(kTRUE,"PerformanceDCAQA.ps");
    TFile fout("PerformanceDCAQA.root","recreate");
    compObjDCA->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceRes * compObjRes = (AliPerformanceRes*)coutput->FindObject("AliPerformanceRes");
    compObjRes->Analyse();
    compObjRes->GetAnalysisFolder()->ls("*");
    compObjRes->PrintHisto(kTRUE,"PerformanceResQA.ps");
    TFile fout("PerformanceResQA.root","recreate");
    compObjRes->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceRes * compObjRes = (AliPerformanceRes*)coutput->FindObject("AliPerformanceResTPCInner");
    compObjRes->Analyse();
    compObjRes->GetAnalysisFolder()->ls("*");
    compObjRes->PrintHisto(kTRUE,"PerformanceResTPCInnerQA.ps");
    TFile fout("PerformanceResTPCInnerQA.root","recreate");
    compObjRes->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceRes * compObjRes = (AliPerformanceRes*)coutput->FindObject("AliPerformanceResTPCOuter");
    compObjRes->Analyse();
    compObjRes->GetAnalysisFolder()->ls("*");
    compObjRes->PrintHisto(kTRUE,"PerformanceResTPCOuterQA.ps");
    TFile fout("PerformanceResTPCOuterQA.root","recreate");
    compObjRes->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceDEdx* compObjDEdx = (AliPerformanceDEdx*)coutput->FindObject("AliPerformanceDEdxTPCInner");
    compObjDEdx->Analyse();
    compObjDEdx->GetAnalysisFolder()->ls("*");
    compObjDEdx->PrintHisto(kTRUE,"PerformanceDEdxTPCInnerQA.ps");
    TFile fout("PerformanceDEdxTPCInnerQA.root","recreate");
    compObjDEdx->GetAnalysisFolder()->Write();
    fout.Close();

 inFile->cd();
    AliPerformanceMatch * compObjMatchTPCITS = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCITS");
    compObjMatchTPCITS->Analyse();
    compObjMatchTPCITS->GetAnalysisFolder()->ls("*");
    compObjMatchTPCITS->PrintHisto(kTRUE,"PerformanceMatchTPCITSQA.ps");
    TFile fout("PerformanceMatchTPCITSQA.root","recreate");
    compObjMatchTPCITS->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceMatch * compObjMatchTPCTRD = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCTRD");
    compObjMatchTPCTRD->Analyse();
    compObjMatchTPCTRD->GetAnalysisFolder()->ls("*");
    compObjMatchTPCTRD->PrintHisto(kTRUE,"PerformanceMatchTPCTRDQA.ps");
    TFile fout("PerformanceMatchTPCTRDQA.root","recreate");
    compObjMatchTPCTRD->GetAnalysisFolder()->Write();
    fout.Close();

inFile->cd();
    AliPerformanceMatch * compObjMatchTPCEFF = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCEFF");
    compObjMatchTPCEFF->Analyse();
    compObjMatchTPCEFF->GetAnalysisFolder()->ls("*");
    compObjMatchTPCEFF->PrintHisto(kTRUE,"PerformanceMatchTPCEFFQA.ps");
    TFile fout("PerformanceMatchTPCEFFQA.root","recreate");
    compObjMatchTPCEFF->GetAnalysisFolder()->Write();
    fout.Close();
  }
  else if(str.CompareTo("NO_MC") == 0) 
  {
    AliPerformanceTPC * compObjTPC = (AliPerformanceTPC*)coutput->FindObject("AliPerformanceTPC");
    compObjTPC->Analyse();
    compObjTPC->GetAnalysisFolder()->ls("*");
    compObjTPC->PrintHisto(kTRUE,"PerformanceTPCQA.ps");
    TFile fout("PerformanceTPCQA.root","recreate");
    compObjTPC->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceDCA * compObjDCA = (AliPerformanceDCA*)coutput->FindObject("AliPerformanceDCA");
    compObjDCA->Analyse();
    compObjDCA->GetAnalysisFolder()->ls("*");
    compObjDCA->PrintHisto(kTRUE,"PerformanceDCAQA.ps");
    TFile fout("PerformanceDCAQA.root","recreate");
    compObjDCA->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceDEdx* compObjDEdx = (AliPerformanceDEdx*)coutput->FindObject("AliPerformanceDEdxTPCInner");
    compObjDEdx->Analyse();
    compObjDEdx->GetAnalysisFolder()->ls("*");
    compObjDEdx->PrintHisto(kTRUE,"PerformanceDEdxTPCInnerQA.ps");
    TFile fout("PerformanceDEdxTPCInnerQA.root","recreate");
    compObjDEdx->GetAnalysisFolder()->Write();
    fout.Close();
    
  inFile->cd();
    AliPerformanceMatch * compObjMatchTPCITS = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCITS");
    compObjMatchTPCITS->Analyse();
    compObjMatchTPCITS->GetAnalysisFolder()->ls("*");
    compObjMatchTPCITS->PrintHisto(kTRUE,"PerformanceMatchTPCITSQA.ps");
    TFile fout("PerformanceMatchTPCITSQA.root","recreate");
    compObjMatchTPCITS->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceMatch * compObjMatchTPCTRD = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCTRD");
    compObjMatchTPCTRD->Analyse();
    compObjMatchTPCTRD->GetAnalysisFolder()->ls("*");
    compObjMatchTPCTRD->PrintHisto(kTRUE,"PerformanceMatchTPCTRDQA.ps");
    TFile fout("PerformanceMatchTPCTRDQA.root","recreate");
    compObjMatchTPCTRD->GetAnalysisFolder()->Write();
    fout.Close();
   
inFile->cd();
    AliPerformanceMatch * compObjMatchTPCEFF = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCEFF");
    compObjMatchTPCEFF->Analyse();
    compObjMatchTPCEFF->GetAnalysisFolder()->ls("*");
    compObjMatchTPCEFF->PrintHisto(kTRUE,"PerformanceMatchTPCEFFQA.ps");
    TFile fout("PerformanceMatchTPCEFFQA.root","recreate");
    compObjMatchTPCEFF->GetAnalysisFolder()->Write();
    fout.Close();
  }
  else if(str.CompareTo("TPC") == 0) 
  {
  inFile->cd();
    AliPerformanceTPC * compObjTPC = (AliPerformanceTPC*)coutput->FindObject("AliPerformanceTPC");
    compObjTPC->Analyse();
    compObjTPC->GetAnalysisFolder()->ls("*");
    compObjTPC->PrintHisto(kTRUE,"PerformanceTPCQA.ps");
    TFile fout("PerformanceTPCQA.root","recreate");
    compObjTPC->GetAnalysisFolder()->Write();
    fout.Close();
  }
  else if(str.CompareTo("EFF") == 0) 
  {
  inFile->cd();
    AliPerformanceEff * compObjEff = (AliPerformanceEff*)coutput->FindObject("AliPerformanceEff");
    compObjEff->Analyse();
    compObjEff->GetAnalysisFolder()->ls("*");
    compObjEff->PrintHisto(kTRUE,"PerformanceEffQA.ps");
    TFile fout("PerformanceEffQA.root","recreate");
    compObjEff->GetAnalysisFolder()->Write();
    fout.Close();
  }
  else if(str.CompareTo("RES") == 0) 
  {
  inFile->cd();
    AliPerformanceRes * compObjRes = (AliPerformanceRes*)coutput->FindObject("AliPerformanceRes");
    compObjRes->Analyse();
    compObjRes->GetAnalysisFolder()->ls("*");
    compObjRes->PrintHisto(kTRUE,"PerformanceResQA.ps");
    TFile fout("PerformanceResQA.root","recreate");
    compObjRes->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceRes * compObjRes = (AliPerformanceRes*)coutput->FindObject("AliPerformanceResTPCInner");
    compObjRes->Analyse();
    compObjRes->GetAnalysisFolder()->ls("*");
    compObjRes->PrintHisto(kTRUE,"PerformanceResTPCInnerQA.ps");
    TFile fout("PerformanceResTPCInnerQA.root","recreate");
    compObjRes->GetAnalysisFolder()->Write();
    fout.Close();
  }
  else if(str.CompareTo("DEDX") == 0) 
  {
  inFile->cd();
    AliPerformanceDEdx* compObjDEdx = (AliPerformanceDEdx*)coutput->FindObject("AliPerformanceDEdxTPCInner");
    compObjDEdx->Analyse();
    compObjDEdx->GetAnalysisFolder()->ls("*");
    compObjDEdx->PrintHisto(kTRUE,"PerformanceDEdxTPCInnerQA.ps");
    TFile fout("PerformanceDEdxTPCInnerQA.root","recreate");
    compObjDEdx->GetAnalysisFolder()->Write();
    fout.Close();
  }
  else if(str.CompareTo("DCA") == 0) 
  {
  inFile->cd();
    AliPerformanceDCA * compObjDCA = (AliPerformanceDCA*)coutput->FindObject("AliPerformanceDCA");
    compObjDCA->Analyse();
    compObjDCA->GetAnalysisFolder()->ls("*");
    compObjDCA->PrintHisto(kTRUE,"PerformanceDCAQA.ps");
    TFile fout("PerformanceDCAQA.root","recreate");
    compObjDCA->GetAnalysisFolder()->Write();
    fout.Close();
  }
  else if(str.CompareTo("MATCH") == 0) 
  {
  inFile->cd();
    AliPerformanceMatch * compObjMatchTPCITS = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCITS");
    compObjMatchTPCITS->Analyse();
    compObjMatchTPCITS->GetAnalysisFolder()->ls("*");
    compObjMatchTPCITS->PrintHisto(kTRUE,"PerformanceMatchTPCITSQA.ps");
    TFile fout("PerformanceMatchTPCITSQA.root","recreate");
    compObjMatchTPCITS->GetAnalysisFolder()->Write();
    fout.Close();

  inFile->cd();
    AliPerformanceMatch * compObjMatchTPCTRD = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCTRD");
    compObjMatchTPCTRD->Analyse();
    compObjMatchTPCTRD->GetAnalysisFolder()->ls("*");
    compObjMatchTPCTRD->PrintHisto(kTRUE,"PerformanceMatchTPCTRDQA.ps");
    TFile fout("PerformanceMatchTPCTRDQA.root","recreate");
    compObjMatchTPCTRD->GetAnalysisFolder()->Write();
    fout.Close();

inFile->cd();
    AliPerformanceMatch * compObjMatchTPCEFF = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCEFF");
    compObjMatchTPCEFF->Analyse();
    compObjMatchTPCEFF->GetAnalysisFolder()->ls("*");
    compObjMatchTPCEFF->PrintHisto(kTRUE,"PerformanceMatchTPCEFFQA.ps");
    TFile fout("PerformanceMatchTPCEFFQA.root","recreate");
    compObjMatchTPCEFF->GetAnalysisFolder()->Write();
    fout.Close();
  }
}
