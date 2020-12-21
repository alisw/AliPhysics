void CompareTime3Stages_v1(){
  //
  // QA macro after step1 of Time Calibration
  // thre is a memory leak somewhere - to be checked where
  //
  char tmpin[256];
  Int_t runno[]={
    259867, 259703, 259091, 260187, 259954, 259471, 259470, 259379, 259095

  };
  Int_t nruns=9;//160;//172;

  for (Int_t i=0;i<nruns;i++){
    sprintf(tmpin,"%d/AnalysisResults.root",runno[i]);
    CompareTime3Stages_1file(tmpin,Form("%d",runno[i]));
 }
}

void CompareTime3Stages_1file(TString filename="timeResults.root",TString period="LHC15i"){
  TFile *file=new TFile(filename.Data());
  TList *list=file->Get("chistolist");
  //  list->ls();
  TH2F *hRawTime[4];
  TH2F *hRawTimeCorr[4];
  TH2F *hTime[4];

  char text[256];
  for(Int_t i=0;i<4;i++){//BC loop
    sprintf(text,"RawTimeVsIdBC%d",i);
    hRawTime[i]=(TH2F *)list->FindObject(text);
    sprintf(text,"RawCorrTimeVsIdBC%d",i);
    hRawTimeCorr[i]=(TH2F *)list->FindObject(text);
    sprintf(text,"TimeVsIdBC%d",i);
    hTime[i]=(TH2F *)list->FindObject(text);

    if(hRawTime[i]!=0) hRawTime[i]->SetMarkerColor(i+1);
    if(hRawTimeCorr[i]!=0) hRawTimeCorr[i]->SetMarkerColor(i+1);
    if(hTime[i]!=0) hTime[i]->SetMarkerColor(i+1);

  }

  TCanvas *c1=new TCanvas("c1","c1",1200,400); 
  c1->Divide(3,1);
  c1->cd(1);
  if(hRawTime[0]!=0) hRawTime[0]->Draw();
  c1->cd(2);
  if(hRawTimeCorr[0]!=0) hRawTimeCorr[0]->Draw();
  c1->cd(3);
  if(hTime[0]!=0) hTime[0]->Draw();
  for(Int_t i=1;i<4;i++){//BC loop
    c1->cd(1);
    if(hRawTime[i]!=0) hRawTime[i]->Draw("same");
    c1->cd(2);
    if(hRawTimeCorr[i]!=0) hRawTimeCorr[i]->Draw("same");
    c1->cd(3);
    if(hTime[0]!=0) hTime[i]->Draw("same");
  }

//  sprintf(text,"plots/RawTimeCompare_%s.pdf",period.Data());
//  c1->Print(text);
  sprintf(text,"plots/RawTimeCompare_%s.jpg",period.Data());
  c1->Print(text);

  TCanvas *c2=new TCanvas("c2","c2",800,800); 
  c2->Divide(2,2);
  for(Int_t i=0;i<4;i++){//BC loop
    c2->cd(i+1)->SetLogz();
    if(hRawTime[i]!=0) hRawTime[i]->Draw("colz");
  }
//  sprintf(text,"plots/RawTimeCompareBC_%s.pdf",period.Data());
//  c2->Print(text);
  sprintf(text,"plots/RawTimeCompareBC_%s.jpg",period.Data());
  c2->Print(text);

  delete c1;
  delete c2;

  delete list;

  file->Close();
  delete file;

}
