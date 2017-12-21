const int kNSM = 20;

void plotProf(const float minRunDiv1k=216,const float maxRunDiv1k=250)
{ // select range of runnumbers to plot
  gStyle->SetOptStat(0);
  TFile *_file0 = TFile::Open("all.root");

  TString cutStr = "nOk>=3&&nPTot>=3&&aveMax>aveMin&&(aveMax-aveMin)<5";

  tree->Draw("aveMax-aveMin", cutStr.Data());
  c1->SetLogz();
  c1->SaveAs("png/diff.png");

  char plotFileName[200];
  char title[100];
  char iSMCut[200];

  for (int iSM=0; iSM<kNSM; iSM++) {
    sprintf(iSMCut,"%s&&iSM==%d", cutStr.Data(), iSM);
    sprintf(title, "aveMin vs runno: SM==%d", iSM);
    
    TProfile *hp = new TProfile("hp",title,100, minRunDiv1k, maxRunDiv1k);
    tree->Draw("aveMin:runno/1000.>>hp", iSMCut, "colz");
    hp->SetMinimum(20); 
    sprintf(plotFileName,"png/SM%dProf.png", iSM); 
    c1->SaveAs(plotFileName);
    delete hp;
  }

}