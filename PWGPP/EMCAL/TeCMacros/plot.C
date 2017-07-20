const int kNSM = 20;

void plot()
{
  TFile *_file0 = TFile::Open("all.root");

  TString cutStr = "nOk>=3&&nPTot>=3&&aveMax>aveMin&&(aveMax-aveMin)<5";

  tree->Draw("aveMin:runno", cutStr.Data(),"colz");
  c1->SetLogz();
  c1->SaveAs("png/all.png");

  char plotFileName[200];
  char iSMCut[200];

  for (int iSM=0; iSM<kNSM; iSM++) {
    sprintf(iSMCut,"%s&&iSM==%d", cutStr.Data(), iSM);

    tree->Draw("aveMin:runno", iSMCut,"colz");
    sprintf(plotFileName,"png/SM%d.png", iSM); 
    c1->SaveAs(plotFileName);

    tree->Draw("aveMin" , iSMCut );
    sprintf(plotFileName,"png/SM%d_1D.png", iSM); 
    c1->SaveAs(plotFileName);

    tree->Draw("aveMax" , iSMCut );
    sprintf(plotFileName,"png/SM%d_1D_aveMax.png", iSM); 
    c1->SaveAs(plotFileName);
  }

}
