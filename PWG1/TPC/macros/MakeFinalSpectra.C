void MakeFinalSpectra(const char *file,  const char* comp) 
{
  gStyle->SetPalette(1);
  gStyle->SetNumberContours(50);

  // open proper input file
  TFile *inFile = TFile::Open(file);
  inFile->cd();

  TList *coutput = gFile->Get("TPC");

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
