THnSparseD *s;


void showHistos(const char * fn="histosBenchmark.root" ){
  
  
  TFile * file = TFile::Open(fn);
  s = (THnSparseD*)file->Get("benchmarkInformation");
  
  
  if(!s){
	cout<<"Cannot open THnSparse!"<<endl;
	return;
  }
  
    s->GetAxis(1)->SetRange(1,160);
    s->GetAxis(2)->SetRange(1,600);
    s->GetAxis(3)->SetRange(1,350);
    s->GetAxis(5)->SetRange(1,400);
  
  createHisto(3,6,"realTimeVsTime");
  createHisto(3,1,"realTimeVsSize","col");
  createHisto(3,5,"realTimeVsnTracks","col");
  createHisto(3,2,"realTimeVsOutsize","col");
  createHisto(3,7,"realTimeVsV0s","col");
  createHisto(2,1,"outVsIn","col");
  createHisto(1,5,"sizeVsnTracks","col");
  

    s->GetAxis(7)->SetRange(2,10);
  
    TH2D * realTimeVsSizeV0s = s->Projection(3,1);
    TCanvas *c = new TCanvas("benchmarks","benchmarks",500,400);
	realTimeVsSizeV0s->SetStats(0);
	realTimeVsSizeV0s->Draw("col");
	c->SaveAs("realTimeVsSizeV0s.pdf");
	c->SaveAs("realTimeVsSizeV0s.root");

  
}

void createHisto( Int_t a1, Int_t a2, const char* name,const char* opt=""){
  
	TH2D*h=s->Projection(a1,a2);
  
    TCanvas *c = new TCanvas("benchmarks","benchmarks",500,400);
	h->SetStats(0);
	h->Draw(opt);
	c->SaveAs( Form("%s.pdf", name ) );
	c->SaveAs( Form("%s.root", name ));
	delete c;
}