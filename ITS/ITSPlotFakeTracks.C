#include <iostream.h>
#include <fstream.h>

void ITSPlotFakeTracks(){  
  
  ifstream in ("AliITSTra.out");
  
  TVector DataOut(10);
  
  ///////////////////////////////// Histograms definition ///////////////////////////////////////////

  
  TH1F *hp=new TH1F("hp","PHI resolution",25,-80.,80.); hp->SetFillColor(4);     
  TH1F *hl=new TH1F("hl","LAMBDA resolution",25,-80.,80.); hl->SetFillColor(4);  
  
  // TH1F *hp=new TH1F("hp","PHI resolution",50,-5.,5.); hp->SetFillColor(4);     
  //TH1F *hl=new TH1F("hl","LAMBDA resolution",50,-5.,5.); hl->SetFillColor(4); 
      
 // TH1F *hpt=new TH1F("hpt","Relative Pt resolution",50,-5.,5.);
 TH1F *hpt=new TH1F("hpt","Relative Pt resolution",40,-20.,20.); 
  hpt->SetFillColor(2); 
  TH1F *hd=new TH1F("hd","Impact parameter distribution ",100,0,2000); 
  hd->SetFillColor(6);
  
  TH1F *hdr=new TH1F("hdr","Dr ",50,-1000,1000);  
  hdr->SetFillColor(kGreen);
  TH1F *hdz=new TH1F("hdz","Dz ",50,-1000,1000);  
  hdz->SetFillColor(kBlue);
  
  /*
  TH1F *hdr=new TH1F("hdr","Dr ",50,-100,100);  
  hdr->SetFillColor(kGreen);
  TH1F *hdz=new TH1F("hdz","Dz ",50,-1000,1000);  
  hdz->SetFillColor(kBlue);  
  */ 

  TH1F *hgood=new TH1F("hgood","Good tracks",10,0,2);
  hgood->Sumw2();  // aggiunto il 15-01-01
  TH1F *hfound=new TH1F("hfound","Found tracks",10,0,2);
  hfound->Sumw2();  // aggiunto il 15-01-01
  TH1F *hfake=new TH1F("hfake","Fake tracks",10,0,2);
  hfake->Sumw2();    
  TH1F *hg=new TH1F("hg","",10,0,2); //efficiency for good tracks  
  hg->SetLineColor(4); hg->SetLineWidth(2);
  TH1F *hf=new TH1F("hf","Efficiency for fake tracks",10,0,2);
  /*hf->SetFillColor(1); hf->SetFillStyle(3013);*/ hf->SetLineColor(4); hf->SetLineWidth(2);

  /////////////////////////////////////////////////////////////////////////////////////////////////// 

  ifstream in1 ("AliITSTrag.out");
  Double_t ptg; 
  for(;;) {
    in1 >> ptg; 
    if( in1.eof() ) break;
    hgood->Fill(ptg);
  }
  in1.close();
   
  Int_t neglabs=0;	 // provvisoria
  for (;;){    
    for (int r=0; r<9; r++) in>>DataOut(r);
    if( in.eof() ) break;

    Double_t ptg=DataOut(0); Double_t labITS=DataOut(1); Double_t labTPC=DataOut(2); 
	 Double_t ptperc=DataOut(3);	  
    Double_t deltalam=DataOut(4); Double_t deltaphi=DataOut(5);
    Double_t Dtot=DataOut(6); Double_t Dr=DataOut(7); Double_t Dz=DataOut(8);	

     if(labITS<0) neglabs++;    // provvisoria
	  if(labITS>=0) hfound->Fill(ptg); else 	{   hfake->Fill(ptg);}

	  if(labITS<0 ) {	      // >=
      hpt->Fill(ptperc);
      hl->Fill(deltalam);
      hp->Fill(deltaphi);
      hd->Fill(Dtot);
      hdr->Fill(Dr);
      hdz->Fill(Dz);
   }    
  }
  

  in.close();
  Stat_t ngood=hgood->GetEntries(); cerr<<"Good tracks "<<ngood<<endl;  
  Stat_t nfound=hfound->GetEntries(); cerr<<"Found tracks "<<nfound<<endl;
  Stat_t nfake=hfake->GetEntries(); cerr<<"Fake tracks "<<nfake<<endl;       
  gStyle->SetOptStat(111110);
  gStyle->SetOptFit(1);	   
  TCanvas *c1=new TCanvas("c1","",0,0,700,700);
  TPad *p1=new TPad("p1","",0,0.5,0.5,1); p1->Draw(); hp->SetXTitle("(mrad)");
  p1->cd(); hp->Draw();  hp->Fit("gaus"); c1->cd();
  TPad *p2=new TPad("p2","",0.5,0.5,1,1); p2->Draw(); hl->SetXTitle("(mrad)");
  p2->cd(); hl->Draw(); hl->Fit("gaus"); c1->cd();
  TPad *p3=new TPad("p3","",0,0,0.5,0.5); p3->Draw(); hpt->SetXTitle("(%)");
  p3->cd(); hpt->Draw(); hpt->Fit("gaus"); c1->cd();
  TPad *p4=new TPad("p4","",0.5,0,1,0.5); p4->Draw(); hd->SetXTitle("(micron)");
  p4->cd(); hd->Draw(); c1->cd();
   
  TCanvas *c3=new TCanvas("c3","",200,200,800,500);
  hfound->Print("all");  // aggiunto il 16-01-01
  hgood->Print("all");  // aggiunto il 16-01-01
  TPad *p7=new TPad("p7","",0,0,0.333,1); p7->Draw(); p7->cd(); hfound->Draw(); c3->cd(); 
  TPad *p8=new TPad("p8","",0.333,0,0.666,1); p8->Draw(); p8->cd(); hfake->Draw(); c3->cd();
  TPad *p9=new TPad("p9","",0.666,0,1,1); p9->Draw(); p9->cd(); hgood->Draw(); c3->cd();
     
  TCanvas *c4=new TCanvas("c4","",300,300,800,500);
  hg->Divide(hfound,hgood,1.,1.); //,"b");
  hf->Divide(hfake,hgood,1.,1.); //,"b");
  hg->SetMaximum(1.4);
  hg->SetYTitle("Tracking efficiency");
  hg->SetXTitle("Pt (GeV/c)");
  hg->Print("all");  // aggiunto il 16-01-01
  hg->Draw();  // to not plot the erro bars    hg->Draw("histo");
  TLine *line1 = new TLine(0,1.0,2,1.0); line1->SetLineStyle(4);
  line1->Draw("same");
  TLine *line2 = new TLine(0,0.9,2,0.9); line2->SetLineStyle(4);
  line2->Draw("histosame");
  
   
  hf->SetFillColor(1);
  hf->SetFillStyle(3013);
  hf->SetLineColor(2);
  hf->SetLineWidth(2);
  hf->Draw("same");  // to not plot the error bars  hf->Draw("histosame");
  
  TText *text = new TText(0.461176,0.248448,"Fake tracks");
  text->SetTextSize(0.05);
  text->Draw();
  text = new TText(0.453919,1.11408,"Good tracks");
  text->SetTextSize(0.05);
  text->Draw();
         
  TCanvas *c2=new TCanvas("c2","",100,100,700,400);
  TPad *p5=new TPad("p5","",0,0,0.5,1); p5->Draw(); hdr->SetXTitle("(micron)");
  p5->cd(); hdr->Draw();  hdr->Fit("gaus"); c2->cd();
  TPad *p6=new TPad("p6","",0.5,0,1,1); p6->Draw(); hdz->SetXTitle("(micron)");
  p6->cd(); hdz->Draw(); hdz->Fit("gaus"); c2->cd();
  
  cout<<"neglabs = "<<neglabs<<"\n";  // provvisoria
}
