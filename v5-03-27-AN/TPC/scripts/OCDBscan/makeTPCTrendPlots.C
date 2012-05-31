/*
  Make a TPC basic calibration trend plots:
  Input  - calibTime.root  tree with summary info per run expected to be in the local directory
         - optional parameters run range can be specified - startRun -endRun
  Output - default plots are saved in the pwd/pic/
  macro to define the picture style (NimStyle.C)  expected to be in the current directory 

  
  Example usage:
  
  aliroot -b -q /u/miranov/AliRoot/trunk/TPC/scripts/OCDBscan/makeTPCTrendPlots.C

*/

TTree * tree=0;
const Double_t kmin=0.01;
const Double_t kmax=0.99;
const Double_t kEpsilon=0.0000001;
Int_t run0=0;  // value to set from outside
Int_t run1=10000000;  // vaue to set  from outside

void makeTPCTrendPlots(Int_t startRun=0, Int_t endRun=1000000){
  //
  // make trend plots of the basic TPC calibration parameters
  //
  run0=startRun;
  run1=endRun;
  gROOT->Macro("NimStyle.C");
  TFile f("calibTime.root");
  tree = (TTree*)f.Get("dcs");
  tree->SetMarkerStyle(25);
  tree->SetMarkerSize(0.4);
  gStyle->SetMarkerSize(0.4);
  tree->SetAlias("isValidITS","abs(dits)<3600");
  tree->SetAlias("isValidCE","min(abs(dcea),abs(dcec))<3600&&max(tdriftCE.fElements[72],tdriftCE.fElements[73])>100");
  tree->SetAlias("isValidCEB","max(abs(dcea),abs(dcec))<3600&&min(tdriftCE.fElements[72],tdriftCE.fElements[73])>100");
  printf("makeTPCTrendPlots.C\n\n");
  printf("Processing DrawDriftTime();\n\n");
  DrawDriftTime();
  printf("DrawDriftRun();\n\n");
  DrawDriftRun();
  printf("DrawDriftCorel();\n\n");
  DrawDriftCorel();
}



void DrawDriftTime(){
  //
  // Draw drift velocity trend grapsh - as function of time
  //
  TCut cutRun=Form("run>%d&&run<%d",run0,run1);
  TCut cutCE="(tdriftCE.fElements[72]>100||tdriftCE.fElements[72]>100)&&min(abs(dcea),abs(dcec))<3600";
  TCut cutITS="abs(dits)<3600";

  Int_t entries=0;
  Double_t max=6;
  Double_t min= -6;
  Double_t maxT=0,minT=0;
  Double_t dmaxmin= 200000;
  TCanvas * canvasDrift = new TCanvas("canvasDriftTime","canvasDriftTime",2000,900);
  canvasDrift->Divide(1,3);
  // P/T part
  entries=tree->Draw("100*ptrel0:time",cutRun,"goff");
  TGraph * graphPTA = new TGraph(entries,tree->GetV2(),tree->GetV1());
  graphPTA->SetTitle("P/T correction A side");
  graphPTA->SetName("PTcorrectionAside");
  //
  entries=tree->Draw("100*ptrel1:time",cutRun,"goff");
  TGraph * graphPTC = new TGraph(entries,tree->GetV2(),tree->GetV1());
  graphPTC->SetTitle("P/T correction C side");
  graphPTC->SetName("PTcorrectionCside");
  //
  entries=tree->Draw("1000*(ptrel0-ptrel1):time",cutRun,"goff");
  TGraph * graphPTAMC = new TGraph(entries,tree->GetV2(),tree->GetV1());
  graphPTAMC->SetTitle("P/T correction A-C side");
  graphPTAMC->SetName("PTcorrectionAMCside");
  //
  //
  entries=tree->Draw("isValidCE*(-tdriftCE.fElements[72]+990)/10:time",cutRun,"goff");
  TGraph * graphCEA = new TGraph(entries,tree->GetV2(),tree->GetV1());
  graphCEA->SetName("CEAside");
  //
  entries=tree->Draw("isValidCE*(-tdriftCE.fElements[73]+990)/10:time",cutRun,"goff");
  TGraph * graphCEC = new TGraph(entries,tree->GetV2(),tree->GetV1());
  graphCEC->SetName("CECside");
  //
  entries=tree->Draw("isValidCE*(tdriftCE.fElements[73]-tdriftCE.fElements[72]):time",cutRun,"goff");
  TGraph * graphCEAMC = new TGraph(entries,tree->GetV2(),tree->GetV1());
  graphCEAMC->SetName("CEAMCside");
  //
  entries=tree->Draw("isValidITS*vdriftITS*100:time",cutRun,"goff");
  TGraph * graphITSTPCA = new TGraph(entries,tree->GetV2(),tree->GetV1());
  graphITSTPCA->SetName("ITSTPCAside");
  //
  entries=tree->Draw("isValidITS*vdriftITS*100:time",cutRun,"goff");
  TGraph * graphITSTPCC = new TGraph(entries,tree->GetV2(),tree->GetV1());
  graphITSTPCC->SetName("ITSTPCCside");
  //
  // Drawing part
  //
  min=-6;
  max=6;
  graphCEA->SetMinimum(min);
  graphCEA->SetMaximum(max);
  graphPTA->SetMinimum(min);
  graphPTA->SetMaximum(max);
  graphITSTPCA->SetMinimum(min);
  graphITSTPCA->SetMaximum(max);
  graphPTA->GetXaxis()->SetRangeUser(minT,maxT);
  graphCEA->GetXaxis()->SetRangeUser(minT,maxT);
  graphITSTPCA->GetXaxis()->SetRangeUser(minT,maxT);
  //
  //
  canvasDrift->cd(1); 
  graphPTA->GetXaxis()->SetTimeDisplay(kTRUE);
  graphPTA->GetYaxis()->SetTitle("#Delta P/T (%)");
  graphPTA->SetMarkerColor(2);graphPTA->SetMarkerStyle(25);
  graphPTC->SetMarkerColor(4);graphPTC->SetMarkerStyle(27);
  graphPTA->Draw("alp");
  graphPTC->Draw("lp");
  graphPTAMC->SetMarkerColor(3);graphPTAMC->SetMarkerStyle(25); 
  graphPTAMC->Draw("lp");
  TLegend *legend = new TLegend(0.11,0.11,0.4,0.35,"P/T correction");
  legend->AddEntry(graphPTA,"A side (%)");
  legend->AddEntry(graphPTC,"C side (%)");
  legend->AddEntry(graphPTAMC,"A-C side (0.1%)");
  legend->Draw();
  //
  canvasDrift->cd(2); 
  graphITSTPCA->GetXaxis()->SetTimeDisplay(kTRUE);
  graphITSTPCA->GetYaxis()->SetTitle("v_{dcorr} (%)");
  graphITSTPCA->SetMarkerColor(2);graphITSTPCA->SetMarkerStyle(25);
  graphITSTPCC->SetMarkerColor(4);graphITSTPCC->SetMarkerStyle(27); 
  graphITSTPCA->Draw("ap");
  graphITSTPCC->Draw("p");
  TLegend *legend = new TLegend(0.11,0.11,0.4,0.35,"Drift correction (TPC-ITS)");
  legend->AddEntry(graphITSTPCA,"A side (%)");
  legend->AddEntry(graphITSTPCC,"C side (%)");
  legend->Draw();
  //
  canvasDrift->cd(3); 
  graphCEA->GetXaxis()->SetTimeDisplay(kTRUE);
  graphCEA->GetYaxis()->SetTitle("(T_{CE0}-T_{CE})/T_{CE0}");
  graphCEA->SetMarkerColor(2);graphCEA->SetMarkerStyle(25);
  graphCEC->SetMarkerColor(4);graphCEC->SetMarkerStyle(27);
  graphCEA->Draw("ap");
  graphCEC->Draw("p");
  graphCEAMC->SetMarkerColor(3);graphCEAMC->SetMarkerStyle(25); 
  graphCEAMC->Draw("p");
  TLegend *legend = new TLegend(0.11,0.11,0.4,0.35,"CE laser time (T_{CE0}=990)");
  legend->AddEntry(graphCEA,"A side (%)");
  legend->AddEntry(graphCEC,"C side (%)");
  legend->AddEntry(graphCEAMC,"A-C side (0.1%)");
  legend->Draw();
  //
  canvasDrift->SaveAs("pic/canvasDriftTime.gif");
}


void DrawDriftRun(){
  //
  //
  TCut cutRun=Form("run>%d&&run<%d",run0,run1);
  TCut cutCE="(tdriftCE.fElements[72]>100||tdriftCE.fElements[72]>100)&&min(abs(dcea),abs(dcec))<3600";
  TCut cutITS="abs(dits)<3600";

  Int_t entries=0;
  Double_t max=-100000;
  Double_t min= 100000;
  Double_t maxT=0,minT=0;
  Double_t dmaxmin= 200000;
  TCanvas * canvasDrift = new TCanvas("canvasDriftRun","canvasDriftRun",2000,900);
  canvasDrift->Divide(1,3);
  //
  // P/T part
  TGraph * graphPTA = TStatToolkit::MakeGraphSparse(tree,"100*ptrel0:run",cutRun);
  graphPTA->SetTitle("P/T correction A side");
  graphPTA->SetName("PTcorrectionAside");
  TGraph * graphPTC = TStatToolkit::MakeGraphSparse(tree,"100*ptrel1:run",cutRun);
  graphPTC->SetTitle("P/T correction C side");
  graphPTC->SetName("PTcorrectionCside");
  TGraph * graphPTAMC = TStatToolkit::MakeGraphSparse(tree,"1000*(ptrel0-ptrel1):run",cutRun);
  graphPTAMC->SetTitle("P/T correction A-C side");
  graphPTAMC->SetName("PTcorrectionAMCside");
  //
  TGraph * graphCEA = TStatToolkit::MakeGraphSparse(tree,"isValidCE*(-tdriftCE.fElements[72]+990)/10:run",cutRun+cutITS);
  graphCEA->SetTitle("P/T correction A side");
  graphCEA->SetName("CEcorrectionAside");
  TGraph * graphCEC = TStatToolkit::MakeGraphSparse(tree,"isValidCE*(-tdriftCE.fElements[73]+990)/10:run",cutRun+cutITS);
  graphCEC->SetTitle("P/T correction C side");
  graphCEC->SetName("CEcorrectionCside");
  TGraph * graphCEAMC = TStatToolkit::MakeGraphSparse(tree,"(tdriftCE.fElements[73]-tdriftCE.fElements[72]):run",cutRun+cutITS);
  graphCEAMC->SetTitle("P/T correction A-C side");
  graphCEAMC->SetName("CEcorrectionAMCside");
  //
  TGraph * graphITSTPCA = TStatToolkit::MakeGraphSparse(tree,"isValidITS*vdriftITS*100:run",cutRun+cutITS);
  graphITSTPCA->SetTitle("P/T correction A side");
  graphITSTPCA->SetName("ITSTPCcorrectionAside");
  TGraph * graphITSTPCC = TStatToolkit::MakeGraphSparse(tree,"isValidITS*vdriftITS*100:run",cutRun+cutITS);
  graphITSTPCC->SetTitle("P/T correction C side");
  graphITSTPCC->SetName("ITSTPCcorrectionCside");
  //
  // Drawing part
  //
  min=-6;
  max=6;
  graphCEA->SetMinimum(min);
  graphCEA->SetMaximum(max);
  graphPTA->SetMinimum(min);
  graphPTA->SetMaximum(max);
  graphITSTPCA->SetMinimum(min);
  graphITSTPCA->SetMaximum(max);
  //
  //
  canvasDrift->cd(1); 
  graphPTA->GetXaxis()->SetTimeDisplay(kTRUE);
  graphPTA->GetYaxis()->SetTitle("#Delta P/T (%)");
  graphPTA->SetMarkerColor(2);graphPTA->SetMarkerStyle(25);
  graphPTC->SetMarkerColor(4);graphPTC->SetMarkerStyle(27);
  graphPTA->Draw("alp");
  graphPTC->Draw("lp");
  graphPTAMC->SetMarkerColor(3);graphPTAMC->SetMarkerStyle(25); 
  graphPTAMC->Draw("lp");
  TLegend *legend = new TLegend(0.11,0.11,0.4,0.35,"P/T correction");
  legend->AddEntry(graphPTA,"A side (%)");
  legend->AddEntry(graphPTC,"C side (%)");
  legend->AddEntry(graphPTAMC,"A-C side (0.1%)");
  legend->Draw();
  //
  canvasDrift->cd(2); 
  graphITSTPCA->GetXaxis()->SetTimeDisplay(kTRUE);
  graphITSTPCA->GetYaxis()->SetTitle("v_{dcorr} (%)");
  graphITSTPCA->SetMarkerColor(2);graphITSTPCA->SetMarkerStyle(25);
  graphITSTPCC->SetMarkerColor(4);graphITSTPCC->SetMarkerStyle(27); 
  graphITSTPCA->Draw("ap");
  graphITSTPCC->Draw("p");
  TLegend *legend = new TLegend(0.11,0.11,0.4,0.35,"Drift correction (TPC-ITS)");
  legend->AddEntry(graphITSTPCA,"A side (%)");
  legend->AddEntry(graphITSTPCC,"C side (%)");
  legend->Draw();
  //
  canvasDrift->cd(3); 
  graphCEA->GetXaxis()->SetTimeDisplay(kTRUE);
  graphCEA->GetYaxis()->SetTitle("(T_{CE0}-T_{CE})/T_{CE0}");
  graphCEA->SetMarkerColor(2);graphCEA->SetMarkerStyle(25);
  graphCEC->SetMarkerColor(4);graphCEC->SetMarkerStyle(27);
  graphCEA->Draw("ap");
  graphCEC->Draw("p");
  graphCEAMC->SetMarkerColor(3);graphCEAMC->SetMarkerStyle(25); 
  graphCEAMC->Draw("p");
  TLegend *legend = new TLegend(0.11,0.11,0.4,0.35,"CE laser time (T_{CE0}=990)");
  legend->AddEntry(graphCEA,"A side (%)");
  legend->AddEntry(graphCEC,"C side (%)");
  legend->AddEntry(graphCEAMC,"A-C side (0.1%)");
  legend->Draw();
  //
  canvasDrift->SaveAs("pic/canvasDriftRun.gif");
}

void DrawDriftCorel(){
  //
  //
  //
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;
  TString *strDelta=0;
  Int_t npointsMax=100000;
  //
  tree->SetAlias("tCEB","(tdriftCE.fElements[72]+tdriftCE.fElements[73]-2000)/2000");
  tree->SetAlias("tCEA","(tdriftCE.fElements[72]-1000)/1000");
  tree->SetAlias("tCEC","(tdriftCE.fElements[73]-1000)/1000");
  
  tree->SetAlias("tCEE","1-((1+(ptrel0+ptrel1)*0.5)*(1+vdriftITS))");
  strDelta= TStatToolkit::FitPlaneConstrain(tree,"tCEB", "tCEE","isValidCEB&&isValidITS", chi2,npoints,param,covar,-1,0, npointsMax, 20);
  //
  strDelta->Tokenize("++")->Print();
  tree->SetAlias("tCEF",strDelta->Data());

  TGraph * graphLT = TStatToolkit::MakeGraphSparse(tree,"2.64*1000*(tCEB-tCEF):run","isValidCEB&&isValidITS"); 
  TGraph * graphLTA = TStatToolkit::MakeGraphSparse(tree,"2.64*1000*(tCEA-tCEF):run","isValidCEB&&isValidITS"); 
  TGraph * graphLTC = TStatToolkit::MakeGraphSparse(tree,"2.64*1000*(tCEC-tCEF):run","isValidCEB&&isValidITS"); 
  graphLT->GetYaxis()->SetTitle("#Delta (mm)");
  graphLT->SetMarkerStyle(25);
  graphLTA->SetMarkerStyle(26); graphLTA->SetMarkerColor(2);
  graphLTC->SetMarkerStyle(27); graphLTC->SetMarkerColor(4);
  graphLT->SetMaximum(10);
  graphLT->SetMinimum(-10);

  TCanvas * canvasDrift = new TCanvas("canvasDriftDiff","canvasDriftDiff",2000,500);
  graphLT->Draw("alp");
  graphLTA->Draw("lp");
  graphLTC->Draw("lp");
  TLegend *legend = new TLegend(0.11,0.11,0.4,0.35,"CE plane  - corrected with tracks");
  legend->AddEntry(graphLT,"AC side mean (mm)");
  legend->AddEntry(graphLTA,"A side (mm)");
  legend->AddEntry(graphLTC,"C side (mm)");
  legend->Draw();
  //
  canvasDrift->SaveAs("pic/canvasDriftDiffRun.gif");
}
