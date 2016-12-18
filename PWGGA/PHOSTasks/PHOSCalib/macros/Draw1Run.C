void Draw1Run(const Int_t runnumber=244628,TString trigger="kINT7", const Int_t gain=1)
{
	TCanvas *c1 = new TCanvas("c1",Form("run%d_%s_g%d",runnumber,trigger.Data(),gain),0,0,1500,900);
	c1->Divide(5,3);

  gStyle->SetOptStat(0);
	TString infile = Form("CellTime_%d_%s.root",runnumber,trigger.Data());
  TFile *rootfile = TFile::Open(infile,"READ");



	for(Int_t iddl=6;iddl<20;iddl++){
    TString histname = Form("h2BC4vsCorrTime_DDL%d_g%d",iddl,gain);
    TH2F *h2 = (TH2F*)rootfile->Get(histname);
    h2->SetTitle(Form("%d %s DDL%d g%d",runnumber,trigger.Data(),iddl,gain));
    c1->cd(iddl-5);
    h2->GetXaxis()->SetRangeUser(-200,200);
    h2->Draw("colz");
	}

}
