void DrawBC4vsTime(const TString trigger="kINT7", const Int_t ddl=11, const Int_t gain=1)
{
  gStyle->SetOptStat(0);
	TCanvas *c1 = new TCanvas("c1",Form("DDL%d_%s_g%d",ddl,trigger.Data(),gain),0,0,1200,800);
	c1->Divide(6,4);

	TString filename = "runlist_LHC15n.txt";
	TString infile="";
	//TFile *rootfile=0x0;

	Int_t runnumber=0;
	Int_t counter=0;
	ifstream fin;
  fin.open(filename);
	while(1){
    counter++;
    fin >> runnumber;
    if(fin.eof()) break;

		//cout << "runnumber = " << runnumber << endl;
		infile = Form("CellTime_%d_%s.root",runnumber,trigger.Data());
		cout << "infile = " << infile << endl;
	
		TFile *rootfile = TFile::Open(infile,"READ");
		TH2F *h2 = (TH2F*)rootfile->Get(Form("hBC4vsRawTimeDDL%d_G%d",ddl,gain));

 //   hBC4vsRawTimeDDL16_G0


		h2->GetXaxis()->SetRangeUser(-50,150);
    h2->SetTitle(Form("%d %s DDL%d g%d",runnumber,trigger.Data(),ddl,gain));
		c1->cd(counter);
		h2->Draw("colz");
	
    //rootfile->Close();
	
	}
  fin.close();

  //c1->SaveAs(Form("WithoutL1corr_DDL%d_%s_g%d.eps",ddl,trigger.Data(),gain));

}

