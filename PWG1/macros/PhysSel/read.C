#include "QALHC11c_QA69.h"
#include "DefinePlots.h"
#include "TGraphErrors.h"

TH1F *fAccOverAll[3];
TH1F *fBGOverAll[3];
TH1F *fV0BGOverAll[3];
TH1F *fV0ABGOverAll[3];
TH1F *fV0CBGOverAll[3];
TH1F *fAfterOverBefore[3];
TH1F *fV0AOverV0C[3];
TH1F *fF02OverAll[3];
TH1F *fF01OverAll[3];

TFile fou;
TLine *l1=new TLine(153776, 0, 153776, 1);
TH1F * GetEmpty(const char * name, Int_t nfile);
TGraphErrors * GetGraphRej(TGraphErrors * gr, TList * rejRunList, const char * reason, Float_t &mean, Bool_t doDraw) ;
double GetMedian(double* arr, int n);
double meanMed(double* vec, int np, double nsigmaCut, int &nrej, int *rejList);
TGraphErrors * graph[kNGraphs];


void read(){
  
  for(Int_t i=0;i<3;++i){
    fAccOverAll[i]=new TH1F(Form("fAccOverAll_%d",i),Form("fAccOverAll_%d",i),100,0,1);
    fBGOverAll[i]=new TH1F(Form("fBGOverAll_%d",i),Form("fBGOverAll_%d",i),50,0,0.5);
    fV0BGOverAll[i]=new TH1F(Form("fV0BGOverAll_%d",i),Form("fV0BGOverAll_%d",i),50,0,0.1);
    fV0ABGOverAll[i]=new TH1F(Form("fV0ABGOverAll_%d",i),Form("fV0ABGOverAll_%d",i),50,0,0.1);
    fV0CBGOverAll[i]=new TH1F(Form("fV0CBGOverAll_%d",i),Form("fV0CBGOverAll_%d",i),50,0,0.1);
    fAfterOverBefore[i]=new TH1F(Form("fAfterOverBefore_%d",i),Form("fAfterOverBefore_%d",i),100,0.5,1);
    fV0AOverV0C[i]=new TH1F(Form("fV0AOverV0C_%d",i),Form("fV0AOverV0C_%d",i),100,0.7,1.2);
    fF01OverAll[i]=new TH1F(Form("fF01OverAll_%d",i),Form("fF01OverAll_%d",i),200,0.6,1.1);
    fF02OverAll[i]=new TH1F(Form("fF02OverAll_%d",i),Form("fF02OverAll_%d",i),200,0.6,1.1);
  }




  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("read.C","");
  dir.ReplaceAll("/./","/");
  
  // read file and add to fit object
  Double_t *run = new Double_t[1330];
  Double_t *ratio1 = new Double_t[1330];
  Double_t *nofill = new Double_t[1330];
  Double_t vX, vY, vZ;
  Int_t vNData = 0;
  ifstream vInput;
  vInput.open(Form("%srct.dat",dir.Data()));
  while (1) {
    vInput >> vX >> vY >> vZ;
    if (!vInput.good()) break;
    run[vNData] = vX;
    ratio1[vNData] = vY;
    nofill[vNData] = vZ;
    vNData++;
  }//while
  vInput.close();

  grin = new TGraph(vNData,run,ratio1);
  cout<<"n points="<<grin->GetN()<<endl;


  Int_t NoFillings=0;
  for(Int_t i=1;i<grin->GetN();++i){
    if(nofill[i]!=nofill[i-1])
      NoFillings++;
  }
  cout<<"No. Dif. de Fill="<<NoFillings<<endl;

  Int_t *boundary_run=new Double_t[NoFillings];
  Int_t *boundary_fill=new Double_t[NoFillings];
  Int_t counter=0;
  for(Int_t i=1;i<grin->GetN();++i){
    if(nofill[i]!=nofill[i-1]){
      boundary_run[counter]=run[i-1];
      boundary_fill[counter]=nofill[i-1];
      counter++;
    }
    
  }


  const Int_t NoOfFill =NoFillings;
  TLine *l2[NoOfFill];

  for(Int_t i=0;i<NoFillings;++i){
    cout<<"interfaz run="<<boundary_run[i]<<";  Fill #"<<boundary_fill[i]<<endl;
    l2[i]=new TLine(boundary_run[i], 0, boundary_run[i], 20);
  }




  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  





  loadlibs();
  TList * listRejectedRuns = new TList(); // keep track of rejected runes

  
  
  
  
  
  for(Int_t i=0;i<3;++i)
    cout<<gnames[i]<<endl;
  
  fou = new TFile("qaLHC10apass2.root");
  for(Int_t igraph = 0; igraph < kNGraphs; igraph++){
    graph[igraph]=(TGraphErrors *)  fou.Get(gnames[igraph]);
  }
  
  GetHisto(graph[kGraphBGOverAllAC],fBGOverAll[0]);
  GetHisto(graph[kGraphBGOverAllAC2],fBGOverAll[1]);
  GetHisto(graph[kGraphBGOverAllSMH],fBGOverAll[2]);
  
  GetHisto(graph[kGraphV0BGOverAllAC],fV0BGOverAll[0]);
  GetHisto(graph[kGraphV0BGOverAllAC2],fV0BGOverAll[1]);
  GetHisto(graph[kGraphV0BGOverAllSMH],fV0BGOverAll[2]);
  
  
  GetHisto(graph[kGraphV0ABGOverAllAC],fV0ABGOverAll[0]);
  GetHisto(graph[kGraphV0ABGOverAllAC2],fV0ABGOverAll[1]);
  GetHisto(graph[kGraphV0ABGOverAllSMH],fV0ABGOverAll[2]);
  
  GetHisto(graph[kGraphV0CBGOverAllAC],fV0CBGOverAll[0]);
  GetHisto(graph[kGraphV0CBGOverAllAC2],fV0CBGOverAll[1]);
  GetHisto(graph[kGraphV0CBGOverAllSMH],fV0CBGOverAll[2]);
  
  GetHisto(graph[kGraphNevACratioAfter],fAfterOverBefore[0]);
  GetHisto(graph[kGraphNevAC2ratioAfter],fAfterOverBefore[1]);
  GetHisto(graph[kGraphNevSMHratioAfter],fAfterOverBefore[2]);

  GetHisto(graph[kGraphV0AOverV0CAC],fV0AOverV0C[0]);
  GetHisto(graph[kGraphV0AOverV0CAC2],fV0AOverV0C[1]);
  GetHisto(graph[kGraphV0AOverV0CSMH],fV0AOverV0C[2]);

  GetHisto(graph[kGraphFO1OverAllAC],fF01OverAll[0]);
  GetHisto(graph[kGraphFO1OverAllAC2],fF01OverAll[1]);
  GetHisto(graph[kGraphFO1OverAllSMH],fF01OverAll[2]);

  GetHisto(graph[kGraphFO2OverAllAC],fF02OverAll[0]);
  GetHisto(graph[kGraphFO2OverAllAC2],fF02OverAll[1]);
  GetHisto(graph[kGraphFO2OverAllSMH],fF02OverAll[2]);



  Float_t meanDummy;
  



  /*
    THREE TRIGGERS IN THE SAME CANVAS
    
  */
  c4a = new TCanvas("GraphACCOverAll","GraphACCOverAll",1000,700);
  

  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();


  c4a->cd();
  pad1->Draw();


  pad1->cd();


  graph[kGraphACCOverAllAC]->Draw("PA");
  graph[kGraphACCOverAllACRej] = GetGraphRej(graph[kGraphACCOverAllAC] , listRejectedRuns, "Acc/All AC" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }

  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  GetHisto(graph[kGraphACCOverAllAC],fAccOverAll[0]);
  fAccOverAll[0]->Draw();





  pad2->cd();


  graph[kGraphACCOverAllACS2]->Draw("PA");
  graph[kGraphACCOverAllACS2Rej] = GetGraphRej(graph[kGraphACCOverAllACS2] , listRejectedRuns, "Acc/All [ACS2]" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }


  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  GetHisto(graph[kGraphACCOverAllACS2],fAccOverAll[1]);
  fAccOverAll[1]->Draw();



  pad3->cd();


  graph[kGraphACCOverAllSMH]->Draw("PA");
  graph[kGraphACCOverAllSMHRej] = GetGraphRej(graph[kGraphACCOverAllSMH] , listRejectedRuns, "Acc/All SMH" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }


  //fAccOverAll[i]

  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  GetHisto(graph[kGraphACCOverAllSMH],fAccOverAll[2]);
  fAccOverAll[2]->Draw();
  fAccOverAll[2]->GetXaxis()->SetTitle("Accepted / All");

  /*
  for(Int_t ipoint = 0; ipoint < graph[kGraphACCOverAllAC]->GetN(); ipoint++){
    fhistTest->Fill(graph[kGraphACCOverAllAC]->GetY()[ipoint]);
  }

  pad1->cd(4);
  fhistTest->Sumw2();
  fhistTest->SetMarkerStyle(25);
  fhistTest->GetXaxis()->SetTitle(ylabels[kGraphACCOverAllAC]);
  fhistTest->GetYaxis()->SetTitle("Entries");
  fhistTest->Draw();
  */
  //pad1->cd(2);





  c4a->Update();
  gSystem->ProcessEvents();
  c4a->SaveAs(Form("picturesLHC11hAOD50/c4a_%s.png",c4a->GetName()));



  c5a = new TCanvas("GraphNev","GraphNev",1000,700);
  
  TPad *pad1 = new TPad("pad1",
			"The pad with the function",0.01,0.01,0.99,0.94,0);
  pad1->Draw();


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();




  pad1->Divide(1,3);
  pad1->cd(1);
  //pad1->cd(1)->SetLogy(1);
  graph[kGraphNevAC]->Draw("PA");

  pad1->cd(2);
  //pad1->cd(2)->SetLogy(1);
  graph[kGraphNevAC2]->Draw("PA");


  pad1->cd(3);
  //pad1->cd(3)->SetLogy(1);
  graph[kGraphNevSMH]->Draw("PA");
 


  c5a->Update();
  gSystem->ProcessEvents();
  c5a->SaveAs(Form("picturesLHC11hAOD50/c5a_%s.png",c5a->GetName()));





  c7a = new TCanvas("GraphBGOverAll","GraphACCOverAll",1000,700);
  

  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();


  c7a->cd();
  pad1->Draw();


  pad1->cd();

  graph[kGraphBGOverAllAC]->Draw("PA");
  graph[kGraphBGOverAllAC] = GetGraphRej(graph[kGraphBGOverAllAC] , listRejectedRuns, "BG/All AC" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }


  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fBGOverAll[0]->Draw();
 






  pad2->cd();

  graph[kGraphBGOverAllAC2]->Draw("PA");
  graph[kGraphBGOverAllAC2] = GetGraphRej(graph[kGraphBGOverAllAC2] , listRejectedRuns, "BG/All AC2" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }


  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fBGOverAll[1]->Draw();



  pad3->cd();

 
  graph[kGraphBGOverAllSMH]->Draw("PA");
  graph[kGraphBGOverAllSMHRej] = GetGraphRej(graph[kGraphBGOverAllSMH] , listRejectedRuns, "BG/All SMH" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }


  pad3a->cd();
  pad3a->cd()->SetGridx(1);

  fBGOverAll[2]->GetXaxis()->SetTitle("BG / All");
  fBGOverAll[2]->Draw();


  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }

  c7a->Update();
  gSystem->ProcessEvents();
  c7a->SaveAs(Form("picturesLHC11hAOD50/c7a_%s.png",c7a->GetName()));



  c7b = new TCanvas("GraphV0BGOverAll","GraphV0BGOverAll",1000,700);
    TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();


  c7b->cd();



  pad1->cd();

  graph[kGraphV0BGOverAllAC]->Draw("PA");
  graph[kGraphV0BGOverAllAC] = GetGraphRej(graph[kGraphV0BGOverAllAC] , listRejectedRuns, "V0BG/All AC" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }



  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fV0BGOverAll[0]->Draw();


  pad2->cd();
  graph[kGraphV0BGOverAllAC2]->Draw("PA");
  graph[kGraphV0BGOverAllAC2] = GetGraphRej(graph[kGraphV0BGOverAllAC2] , listRejectedRuns, "V0BG/All AC2" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }


  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fV0BGOverAll[1]->Draw();

  pad3->cd();
  graph[kGraphV0BGOverAllSMH]->Draw("PA");
  graph[kGraphV0BGOverAllSMHRej] = GetGraphRej(graph[kGraphV0BGOverAllSMH] , listRejectedRuns, "V0BG/All SMH" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }


  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fV0BGOverAll[2]->GetXaxis()->SetTitle("V0BG / All");
  fV0BGOverAll[2]->Draw();



  c7b->Update();
  gSystem->ProcessEvents();
  c7b->SaveAs(Form("picturesLHC11hAOD50/c7b_%s.png",c7b->GetName()));










  c7c = new TCanvas("GraphV0A_BGOverAll","GraphV0A_BGOverAll",1000,700);
   TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();


  c7c->cd();

  pad1->cd(); 


  graph[kGraphV0ABGOverAllAC]->Draw("PA");
  graph[kGraphV0ABGOverAllAC] = GetGraphRej(graph[kGraphV0ABGOverAllAC] , listRejectedRuns, "V0A_BG/All AC" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();


  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fV0ABGOverAll[0]->Draw();



  pad2->cd(); 
  graph[kGraphV0ABGOverAllAC2]->Draw("PA");
  graph[kGraphV0ABGOverAllAC2] = GetGraphRej(graph[kGraphV0ABGOverAllAC2] , listRejectedRuns, "V0A_BG/All AC2" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();


  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fV0ABGOverAll[1]->Draw();

  pad3->cd(); 
  graph[kGraphV0ABGOverAllSMH]->Draw("PA");
  graph[kGraphV0ABGOverAllSMHRej] = GetGraphRej(graph[kGraphV0ABGOverAllSMH] , listRejectedRuns, "V0A_BG/All SMH" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();


  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fV0ABGOverAll[2]->Draw();
  fV0ABGOverAll[2]->GetXaxis()->SetTitle("V0ABG / All");


  c7c->Update();
  gSystem->ProcessEvents();
  c7c->SaveAs(Form("picturesLHC11hAOD50/c7c_%s.png",c7c->GetName()));





  c7d = new TCanvas("GraphV0C_BGOverAll","GraphV0C_BGOverAll",1000,700);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();


  c7d->cd();

  pad1->cd();
  graph[kGraphV0CBGOverAllAC]->Draw("PA");
  graph[kGraphV0CBGOverAllAC] = GetGraphRej(graph[kGraphV0CBGOverAllAC] , listRejectedRuns, "V0C_BG/All AC" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();


  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fV0CBGOverAll[0]->Draw();


  pad2->cd();
  graph[kGraphV0CBGOverAllAC2]->Draw("PA");
  graph[kGraphV0CBGOverAllAC2] = GetGraphRej(graph[kGraphV0CBGOverAllAC2] , listRejectedRuns, "V0C_BG/All AC2" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fV0CBGOverAll[1]->Draw();


	   pad3->cd();
  graph[kGraphV0CBGOverAllSMH]->Draw("PA");
  graph[kGraphV0CBGOverAllSMHRej] = GetGraphRej(graph[kGraphV0CBGOverAllSMH] , listRejectedRuns, "V0A_BG/All SMH" ,  meanDummy, 1);


  l1->SetLineStyle(3);
  l1->Draw();

  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fV0CBGOverAll[2]->Draw();
  fV0CBGOverAll[2]->GetXaxis()->SetTitle("V0CBG / All");



  c7d->Update();
  gSystem->ProcessEvents();
  c7d->SaveAs(Form("picturesLHC11hAOD50/c7d_%s.png",c7d->GetName()));





  c8a = new TCanvas("GraphAfterOverBefore","GraphAfterOverBefore",1000,700);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();
  
  pad1->cd();
  
  graph[kGraphNevACratioAfter]->Draw("PA");
  graph[kGraphNevACratioAfterRej] = GetGraphRej(graph[kGraphNevACratioAfter] , listRejectedRuns, "V0BG/All AC" ,  meanDummy, 1);

  l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }




  
  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fAfterOverBefore[0]->Draw();


  pad2->cd();
  graph[kGraphNevAC2ratioAfter]->Draw("PA");
  graph[kGraphNevAC2ratioAfterRej] = GetGraphRej(graph[kGraphNevAC2ratioAfter] , listRejectedRuns, "V0BG/All AC2" ,  meanDummy, 1);

 l1->SetLineStyle(3);
  l1->Draw();


  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }


  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fAfterOverBefore[1]->Draw();

  pad3->cd();
  graph[kGraphNevSMHratioAfter]->Draw("PA");
  graph[kGraphNevSMHratioAfterRej] = GetGraphRej(graph[kGraphNevSMHratioAfter] , listRejectedRuns, "V0BG/All SMH" ,  meanDummy, 1);

 l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }



  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fAfterOverBefore[2]->Draw();
  fAfterOverBefore[2]->GetXaxis()->SetTitle("After_PS / Before_PS");


  c8a->Update();
  gSystem->ProcessEvents();
  c8a->SaveAs(Form("picturesLHC11hAOD50/c8a_%s.png",c8a->GetName()));





  c8b = new TCanvas("GraphNevBefore","GraphNevBefore",1000,700);
  
  TPad *pad1 = new TPad("pad1",
			"The pad with the function",0.01,0.01,0.99,0.94,0);
  pad1->Draw();


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();




  pad1->Divide(1,3);

  pad1->cd(1);
  graph[kGraphNevACBefore]->Draw("PA");

  pad1->cd(2);
  graph[kGraphNevAC2Before]->Draw("PA");

  pad1->cd(3);
  graph[kGraphNevSMHBefore]->Draw("PA");



  c8b->Update();
  gSystem->ProcessEvents();
  c8b->SaveAs(Form("picturesLHC11hAOD50/c8b_%s.png",c8b->GetName()));





  c8c = new TCanvas("GraphNevAfter","GraphNevAfter",1000,700);
  
  TPad *pad1 = new TPad("pad1",
			"The pad with the function",0.01,0.01,0.99,0.94,0);
  pad1->Draw();


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();




  pad1->Divide(1,3);

  pad1->cd(1);
  graph[kGraphNevACAfter]->Draw("PA");


  pad1->cd(2);
  graph[kGraphNevAC2After]->Draw("PA");


  pad1->cd(3);

  graph[kGraphNevSMHAfter]->Draw("PA");



  c8c->Update();
  gSystem->ProcessEvents();
  c8c->SaveAs(Form("picturesLHC11hAOD50/c8c_%s.png",c8c->GetName()));








  c9a = new TCanvas("GraphV0AOverV0C","GraphV0AOverV0C",1000,700);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();

  pad1->cd();
  graph[kGraphV0AOverV0CAC]->Draw("PA");
  graph[kGraphV0AOverV0CACRej] = GetGraphRej(graph[kGraphV0AOverV0CAC] , listRejectedRuns, "Nev_V0A/Nev_V0C AC" ,  meanDummy, 1);


 l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }


  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fV0AOverV0C[0]->Draw();


  pad2->cd();
  graph[kGraphV0AOverV0CAC2]->Draw("PA");
  graph[kGraphV0AOverV0CAC2Rej] = GetGraphRej(graph[kGraphV0AOverV0CAC2] , listRejectedRuns, "Nev_V0A/Nev_V0C AC2" ,  meanDummy, 1);

 l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }



  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fV0AOverV0C[1]->Draw();


  pad3->cd();
  graph[kGraphV0AOverV0CSMH]->Draw("PA");
  graph[kGraphV0AOverV0CSMHRej] = GetGraphRej(graph[kGraphV0AOverV0CSMH] , listRejectedRuns, "Nev_V0A/Nev_V0C SMH" ,  meanDummy, 1);

 l1->SetLineStyle(3);
  l1->Draw();


  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }


  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fV0AOverV0C[2]->Draw();
  fV0AOverV0C[2]->GetXaxis()->SetTitle("V0A / VOC");



  c9a->Update();
  gSystem->ProcessEvents();
  c9a->SaveAs(Form("picturesLHC11hAOD50/c9a_%s.png",c9a->GetName()));



  c10a = new TCanvas("GraphFO1OverAll","GraphFO1OverAll",1000,700);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();  


  pad1->cd();
  graph[kGraphFO1OverAllAC]->Draw("PA");
  graph[kGraphFO1OverAllACRej] = GetGraphRej(graph[kGraphFO1OverAllAC] , listRejectedRuns, "Nev_FO1/All AC" ,  meanDummy, 1);

 l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }



  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fF01OverAll[0]->Draw();


  pad2->cd();
  graph[kGraphFO1OverAllAC2]->Draw("PA");
  graph[kGraphFO1OverAllAC2Rej] = GetGraphRej(graph[kGraphFO1OverAllAC2] , listRejectedRuns, "Nev_FO1/All AC2" ,  meanDummy, 1);

 l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }



  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fF01OverAll[1]->Draw();



  pad3->cd();
  graph[kGraphFO1OverAllSMH]->Draw("PA");
  graph[kGraphFO1OverAllSMHRej] = GetGraphRej(graph[kGraphFO1OverAllSMH] , listRejectedRuns, "Nev_FO1/All SMH" ,  meanDummy, 1);

 l1->SetLineStyle(3);
  l1->Draw();

  for(Int_t i=25;i<NoFillings;++i){
    l2[i]->SetLineStyle(3);
    l2[i]->SetLineColor(4);
    l2[i]->Draw();
  }



  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fF01OverAll[2]->Draw();
  fF01OverAll[2]->GetXaxis()->SetTitle("F01 / All");




  c10a->Update();
  gSystem->ProcessEvents();
  c10a->SaveAs(Form("picturesLHC11hAOD50/c10a_%s.png",c10a->GetName()));




  c10b = new TCanvas("GraphFO2OverAll","GraphFO2OverAll",1000,700);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.635,0.7,0.94,0);
  TPad * pad2=new TPad("pad2","pad2",0.01,0.33,0.7,0.635,0);
  TPad * pad3=new TPad("pad3","pad3",0.01,0.01,0.7,0.33,0);

  TPad * pad1a=new TPad("pad1a","pad1a",0.7,0.635,0.99,0.94,0);
  TPad * pad2a=new TPad("pad2a","pad2a",0.7,0.33,0.99,0.635,0);
  TPad * pad3a=new TPad("pad3a","pad3a",0.7,0.01,0.99,0.33,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1a->Draw();
  pad2a->Draw();
  pad3a->Draw();


  pad1->SetBottomMargin(0);
  pad1->SetBorderSize(0);
  pad1->SetRightMargin(0.01);

  pad2->SetBottomMargin(0.0);
  pad2->SetTopMargin(0);
  pad2->SetRightMargin(0.01);
  pad2->SetBorderSize(0);

  pad3->SetBottomMargin(0.2);
  pad3->SetTopMargin(0);
  pad3->SetRightMargin(0.01);
  pad3->SetBorderSize(0);

  pad1a->SetBottomMargin(0);
  pad1a->SetBorderSize(0);
  pad1a->SetRightMargin(0.01);

  pad2a->SetBottomMargin(0.0);
  pad2a->SetTopMargin(0);
  pad2a->SetRightMargin(0.01);
  pad2a->SetBorderSize(0);

  pad3a->SetBottomMargin(0.2);
  pad3a->SetTopMargin(0);
  pad3a->SetRightMargin(0.01);
  pad3a->SetBorderSize(0);


  // Draw a global picture title
  TPaveLabel *title = new TPaveLabel(0.01,0.95,0.99,0.99,
				     Form("%s",period));

  title->SetFillColor(0);
  title->SetTextFont(52);
  title->SetTextColor(4);
  title->SetTextSize(0.7);
  title->Draw();  


  pad1->cd();
  graph[kGraphFO2OverAllAC]->Draw("PA");
  graph[kGraphFO2OverAllACRej] = GetGraphRej(graph[kGraphFO2OverAllAC] , listRejectedRuns, "Nev_FO2/All AC" ,  meanDummy, 1);

  pad1a->cd();
  pad1a->cd()->SetGridx(1);
  fF02OverAll[0]->Draw();
  fF02OverAll[0]->GetXaxis()->SetTitle("F02 / All");


  pad2->cd();
  graph[kGraphFO2OverAllAC2]->Draw("PA");
  graph[kGraphFO2OverAllAC2Rej] = GetGraphRej(graph[kGraphFO2OverAllAC2] , listRejectedRuns, "Nev_FO2/All AC2" ,  meanDummy, 1);

  pad2a->cd();
  pad2a->cd()->SetGridx(1);
  fF02OverAll[1]->Draw();
  fF02OverAll[1]->GetXaxis()->SetTitle("F02 / All");

  pad3->cd();
  graph[kGraphFO2OverAllSMH]->Draw("PA");
  graph[kGraphFO2OverAllSMHRej] = GetGraphRej(graph[kGraphFO2OverAllSMH] , listRejectedRuns, "Nev_FO2/All SMH" ,  meanDummy, 1);

  pad3a->cd();
  pad3a->cd()->SetGridx(1);
  fF02OverAll[2]->Draw();
  fF02OverAll[2]->GetXaxis()->SetTitle("F02 / All");





  c10b->Update();
  gSystem->ProcessEvents();
  c10b->SaveAs(Form("picturesLHC11hAOD50/c10b_%s.png",c10b->GetName()));








}
double meanMed(double* vec, int np, double nsigmaCut, int &nrej, int *rejList)
{
  // compute the mean of the array "vec" rejecting the outliers
  // if rejlist array is provided, fill indices of rejected values
  // This method works assuming that the "good" points are nearly normally 
  // distrubted around the mean (i.e. there is no significant trend)
  //
  // create copy of the vector
  double *vec1 = new Double_t[np];
  memcpy(vec1, vec, np*sizeof(double));
  //
  // get median of the vector as a robust initial estimate of the mean
  cout << "Points " << np << endl;
  double md = GetMedian(vec1,np);
  //
  // compute squared differences to median, 
  for (int i=0;i<np;i++) {
    vec1[i] = TMath::Abs(vec1[i] - md);
  }
  //
  // compute median squared difference for robust estimate of the sigma
  double sigmd = GetMedian(vec1,np);
  //
  printf("Median Value: %+e | Median abs residual: %e\n",md,sigmd);
  //
  sigmd /= 0.6745; // go to sigma assuming normal distribution
  printf("Estimate of sigma: %+e\n",sigmd);
  // if the estimate of the sigma is null, do not look for outliers

  cout<<"md="<<md<<endl;

  if(!sigmd) return md;

 


  // compute mean rejecting more than nsigmaCut deviations
  double mean = 0;
  int npok = 0;
  for (int i=0;i<np;i++) {
    double dev =  TMath::Abs(vec[i]-md)/sigmd;
    if ( dev < nsigmaCut ) {
      mean += vec[i];
      npok++;
    }
    else {
      printf("Reject#%d (Y=%+e) : deviation: %e\n",i,vec[i],dev);
      if (rejList) rejList[nrej] = i; 
      nrej++;
    }
  }
  //
  delete vec1;
  return npok ? mean/npok : 0;
  //


  cout<<"ending meanMed"<<endl;

}

double GetMedian(double* arr, int n)
{
  // get median by Wirths algroithm
  int i=0,j=0,l=0,m=0;
  double x;
  l=0 ; m=n-1;
  int k = (n&1) ? n/2 : n/2-1;
  while (l<m) {
    x=arr[k] ;
    i=l ;
    j=m ;
    do {
      //      cout << arr[i] << " " << arr[j] << i<<","<<j << " " << arr[k]  << " " << k << " " << x<< endl;
      
      while (arr[i]<x) i++ ;
      while (x<arr[j]) j-- ;
      if (i<=j) { // swap elements	  
	//	cout << "Swapping" << endl;	  
	double t = arr[i];
	arr[i] = arr[j];
	arr[j] = t;
	i++ ; j-- ;
      }      
    } while (i<=j) ;
    if (j<k) l=i ;
    if (k<i) m=j ;
  }
  return arr[k] ;
}

// double GetMedian(double* arr, int n)
// {
//   // get median by Wirths algroithm
//   int i,j,l,m;
//   double x;
//   l=0 ; m=n-1;
//   int k = (n&1) ? n/2 : n/2-1;
//   while (l<m) {
//     x=arr[k] ;
//     i=l ;
//     j=m ;
//     do {
//       //      cout << i << " " << j << endl;
      
//       //skip runs which were not set (value is -1)
//       while (arr[i]<x) i++ ;
//       while (x<arr[j]) j-- ;
//       if(i<=j) { // swap elements
// 	// if((TMath::Abs(arr[i]+1)< 0.0001) && ((TMath::Abs(arr[j]+1)<0.0001))){
// 	//   i++ ; j-- ;
// 	// } else {
// 	  double t = arr[i];
// 	  arr[i] = arr[j];
// 	  arr[j] = t;
// 	  i++ ; j-- ;
// 	  //	}
//       }      
//     } while (i<=j) ;
//     if (j<k) l=i ;
//     if (k<i) m=j ;
    
//   }
//   return arr[k] ;
// }

TGraphErrors * GetGraphRej(TGraphErrors * gr, TList * rejRunList, const char * reason, Float_t &mean, Bool_t doDraw) {

  //Returns a new graph containing only rejected points

  const Double_t nsigmaCut = 3;

  int *rejList = new Int_t[gr->GetN()];
  int nrej = 0;

  Double_t * array = new Double_t[gr->GetN()];
  Int_t * correspondenceFullArray = new Int_t[gr->GetN()];
  Int_t npoint = 0; 
  //exclude from the array all the -1
  for(Int_t ipoint = 0; ipoint < gr->GetN(); ipoint++){
    if (TMath::Abs(gr->GetY()[ipoint]+1)>0.001) { // skip points for which there is no file (==-1)
      array[npoint] = gr->GetY()[ipoint];
      correspondenceFullArray[npoint] = ipoint;
      npoint++;
    } 
  }

  cout<<"calling meanMed"<<endl;
  mean = meanMed(array,npoint,nsigmaCut, nrej, rejList);

  cout<<"nrej="<<nrej<<"  mean="<<mean<<endl;
  
 
    
  TGraphErrors *grrej = new TGraphErrors(nrej);
  for (int i=0;i<nrej;i++) {
    grrej->SetPoint(i, gr->GetX()[correspondenceFullArray[rejList[i]]], gr->GetY()[correspondenceFullArray[rejList[i]]]);
    grrej->SetPointError(i, gr->GetEX()[correspondenceFullArray[rejList[i]]], gr->GetEY()[correspondenceFullArray[rejList[i]]]);
    if(!knownProblems.Contains(Form("%d",(int)gr->GetX()[correspondenceFullArray[rejList[i]]])))
      rejRunList->Add(new TObjString(Form("[%d], %s",  (int)gr->GetX()[correspondenceFullArray[rejList[i]]], reason)));
  }
  grrej->SetMarkerColor(kRed);
  //grrej->SetMarkerStyle(gr->GetMarkerStyle());
  grrej->SetMarkerStyle(29);
  grrej->SetMarkerSize(1.5);


  delete rejList;

  if(doDraw) {
    Float_t minXDraw = gr->GetXaxis()->GetBinLowEdge(1);
    Float_t maxXDraw = gr->GetXaxis()->GetBinLowEdge(gr->GetXaxis()->GetNbins());

    grrej->Draw("P");
    TLine * l = new TLine (minXDraw, mean,maxXDraw, mean);
    l->SetLineStyle(kDashed);
    l->Draw();


  }

  return grrej;



} 

void loadlibs()
{
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libMinuit");
  gSystem->Load("libPWG2spectra");
  gSystem->Load("libPWG0base"); 
}
void GetHisto(TGraphErrors * gr, TH1F *ftemp)
{

  for(Int_t ipoint = 0; ipoint < gr->GetN(); ipoint++){
  //for(Int_t ipoint = 0; ipoint < gr->GetN(); ipoint++){
    ftemp->Fill(gr->GetY()[ipoint]);
  }
  

  ftemp->Sumw2();
  ftemp->SetMarkerColor(kRed+2);
  ftemp->SetLineColor(kRed+2);
  ftemp->SetMarkerStyle(30);
  ftemp->GetYaxis()->SetTitle("Entries");
  ftemp->GetXaxis()->SetLabelSize(0.07);
  ftemp->GetXaxis()->SetTitleSize(0.09);
  ftemp->GetYaxis()->SetLabelSize(0.04);
  ftemp->GetYaxis()->SetTitleSize(0.05);
  ftemp->GetYaxis()->SetTitleOffset(0.7);
  ftemp->GetYaxis()->CenterTitle(0);
  ftemp->GetYaxis()->SetLabelFont(42);
  ftemp->GetYaxis()->SetTitleFont(42);
  ftemp->GetYaxis()->SetNoExponent(kTRUE);
  ftemp->GetXaxis()->SetLabelFont(42);
  ftemp->GetXaxis()->SetTitleFont(42);
  ftemp->GetXaxis()->SetNdivisions(5);


 
  return ftemp;
 
}
