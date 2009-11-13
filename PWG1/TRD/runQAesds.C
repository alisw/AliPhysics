Int_t ExtractRunNumber(TString filename);
void Plot(Int_t run);
void runQAesds(const char *esdlist)
{
  ifstream input(esdlist);
  string str;
  getline(input, str);
  TString tstr(str.c_str());
  if(tstr.Contains("alien://")){
    printf("Data from alien\n");
    gSystem->Load("libRAliEn");
    TGrid::Connect("alien://");
  }

  Int_t runnr = ExtractRunNumber(tstr);
  printf("runnr %d\n", runnr);

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWG1");
  AliLog::SetGlobalLogLevel(AliLog::kError);
  gROOT->LoadMacro(Form("%s/PWG0/CreateESDChain.C", gSystem->ExpandPathName("$ALICE_ROOT")));
  TChain *chain = CreateESDChain(esdlist, -1);
  chain->SetBranchStatus("*FMD*",0);
  chain->SetBranchStatus("*Calo*",0);
  chain->SetBranchStatus("Tracks", 1);
  chain->SetBranchStatus("ESDfriend*",1);
  chain->Lookup();
  chain->GetListOfFiles()->Print();
  printf("\n ----> CHAIN HAS %d ENTRIES <----\n\n", (Int_t)chain->GetEntries());
    
  AliAnalysisManager *TRDqa = new AliAnalysisManager("TRD Cosmics QA");
  TRDqa->SetInputEventHandler(new AliESDInputHandler);
  AliTRDinfoGen *info = new AliTRDinfoGen();
  AliTRDcheckDET *detector = new AliTRDcheckDET();
  detector->SetDebugLevel(10);
  detector->SetMCdata(kFALSE);
  detector->UseClustersOutsideChamber();
  TRDqa->AddTask(info);
  TRDqa->AddTask(detector);
  
  AliAnalysisDataContainer *co[3];
  TRDqa->ConnectInput( info, 0, TRDqa->GetCommonInputContainer());
  co[0] = TRDqa->CreateContainer("trackInfo", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[1] = TRDqa->CreateContainer("eventInfo", AliTRDeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  co[2] = TRDqa->CreateContainer("v0Info", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  TRDqa->ConnectOutput(info, 0, co[0]);
  TRDqa->ConnectOutput(info, 1, co[1]);
  TRDqa->ConnectOutput(info, 2, co[2]);
  TRDqa->ConnectInput(detector, 0, co[0]);
  TRDqa->ConnectInput(detector, 1, co[1]);
  TRDqa->ConnectOutput(detector, 0, TRDqa->CreateContainer(detector->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Perf%d.root", runnr)));
  
  if(TRDqa->InitAnalysis()){
	  TRDqa->Print();
	  TRDqa->StartAnalysis("local", chain);
  }
  gSystem->Exec(Form("mv TRD.DebugPerformance.root TRD.DbgPerf%d.root", runnr));
  Plot(runnr);
}

void Plot(Int_t run){
  TFile *f=new TFile(Form("TRD.Perf%d.root",run));
  TH1F *hNtrSec = (TH1F*)checkDET->FindObject("hNtrksSector");
  TH1F *hNtrks = (TH1F*)checkDET->FindObject("hNtrks"); // tracks/event
  TH1F *hNtlsSTA = (TH1F*)checkDET->FindObject("hNtlsSTA");
  TH1F *hNtlsBAR = (TH1F*)checkDET->FindObject("hNtlsBAR");
  TH1F *hNclTrk = (TH1F*)checkDET->FindObject("hNclTls");
  TProfile *hPHp = (TProfile*) checkDET->FindObject("<PH>")->FindObject("hPHt");
  f->Close();

  TH2F *hPH = new TH2F("hPH", "PH", 30, 0, 30, 100, 0, 1000);
  TString filename = Form("TRD.DbgPerf%d.root",run);
  TChain *debchain = new TChain("PHt");
  if(filename.EndsWith(".root")) {  
    debchain->AddFile(filename.Data());
    debchain->Draw("Charge:Timebin>>hPH", "","colz");
  } else{
    cout << "No such file... " << filename <<endl;
  }

  const Float_t Ls1=0.04;
  TCanvas *tc = new TCanvas ("Monitor", "Monitor",50,50,600,900);
  tc->Divide(2,3);

  tc->cd(1); 
  pad = tc->cd(1); pad->SetLogy();
  hNtrks->GetXaxis()->SetTitle("Tracks per event");
  hNtrks->GetXaxis()->SetTitleOffset(.9);
  hNtrks->GetXaxis()->SetTitleSize(0.05);
  hNtrks->GetYaxis()->SetTitle("");
  //hNtrks->SetTitle(Form("%d",run));
  hNtrks->SetLineWidth(2);
  hNtrks->SetLabelSize(Ls1,"xy");
  hNtrks->GetXaxis()->SetRangeUser(0,20);    
  hNtrks->Draw();
  int Ntracks=hNtrks->Integral();
  TLatex *text=new TLatex(8,10,Form("%d tracks",Ntracks));
  text->SetTextSize(0.045); text->SetTextFont(42);
  text->Draw("same");

  tc->cd(2);
  hNtrSec->GetXaxis()->SetTitle("Sector");
  hNtrSec->GetXaxis()->SetTitleOffset(.9);
  hNtrSec->GetXaxis()->SetTitleSize(0.05);
  hNtrSec->GetYaxis()->SetTitle("Tracks");
  //hNtrSec->GetYaxis()->SetTitleSize(0.05);
  hNtrSec->GetYaxis()->SetTitleOffset(1.2);
  hNtrSec->SetLineWidth(2);
  hNtrSec->SetLabelSize(Ls1,"xy");
  hNtrSec->Draw();

  tc->cd(3);
  hNtlsBAR->GetXaxis()->SetTitle("Tracklets per track");
  hNtlsBAR->GetXaxis()->SetTitleOffset(.9);
  hNtlsBAR->GetXaxis()->SetTitleSize(0.05);
  hNtlsBAR->GetYaxis()->SetTitle("");
  //hNtlsBAR->SetStats(1111);
  hNtlsBAR->SetLabelSize(Ls1,"xy");
  hNtlsBAR->SetLineWidth(2);
  //hNtlsBAR->GetYaxis()->SetRangeUser(0,500);    
  hNtlsSTA->GetXaxis()->SetTitle("Tracklets per track");
  hNtlsSTA->GetYaxis()->SetTitle("");
  hNtlsBAR->SetLineColor(2); hNtlsBAR->Draw("");
  hNtlsSTA->SetLineColor(4); hNtlsSTA->Draw("same");
  leg = new TLegend(0.14,0.35,0.65,0.5);
  leg->AddEntry(hNtlsBAR,"Barrel","l");
  leg->AddEntry(hNtlsSTA,"Standalone","l");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetMargin(0.15); //separation symbol-text
  leg->SetEntrySeparation(0.1); //this seems not to do much...
  leg->Draw();

  tc->cd(4);
  hNclTrk->GetXaxis()->SetTitle("Clusters per tracklet");
  hNclTrk->GetXaxis()->SetTitleOffset(.9);
  hNclTrk->GetXaxis()->SetTitleSize(0.05);
  hNclTrk->GetYaxis()->SetTitle("");
  hNclTrk->SetLabelSize(Ls1,"xy");
  hNclTrk->SetLineWidth(2);
  hNclTrk->Draw();

  tc->cd(5);
  hPHp->GetXaxis()->SetTitle("Time bin (100 ns)");
  hPHp->GetXaxis()->SetTitleOffset(.9);
  hPHp->GetXaxis()->SetTitleSize(0.05);
  hPHp->GetYaxis()->SetTitle("<PH> (ADC ch.)");
  hPHp->GetYaxis()->SetTitleOffset(1.3);
  hPHp->SetLabelSize(Ls1,"xy");
  hPHp->SetLineWidth(2);
  hPHp->Draw();

  tc->cd(6);
  hPH->GetXaxis()->SetTitle("Time bin (100 ns)");
  hPH->GetXaxis()->SetTitleOffset(.9);
  hPH->GetXaxis()->SetTitleSize(0.05);
  hPH->GetYaxis()->SetTitle("Pulse Height (ADC ch.)");
  hPH->GetYaxis()->SetTitleOffset(1.3);
  hPH->SetLabelSize(Ls1,"xy");
  hPH->SetLineWidth(2);
  hPH->Draw("colz");

  //tc->SaveAs(Form("Mon_run%dp%d.eps",run,pass));
  tc->SaveAs(Form("Mon_run%d.gif",run));
}

Int_t ExtractRunNumber(TString filename){
  TObjArray *toks = filename.Tokenize("/");
  TObjString *file = dynamic_cast<TObjString *>(toks->At(toks->GetEntries() - 2));
  TString &filestr = file->String();
  TString runstr = filestr(2,9);
  return runstr.Atoi();
}