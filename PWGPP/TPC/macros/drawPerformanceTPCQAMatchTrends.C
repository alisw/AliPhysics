
// 05.2013 new functionality: Status Bar
// to set the variables, do the following. more info & examples in the TStatToolkit.
/*{
  TStatToolkit::SetStatusAlias(tree, "meanTPCncl",    "", "varname_Out:(abs(varname-MeanEF)>6.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "tpcItsMatchA",  "", "varname_Out:(abs(varname-MeanEF)>6.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "meanTPCncl",    "", "varname_Warning:(abs(varname-MeanEF)>3.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "tpcItsMatchA",  "", "varname_Warning:(abs(varname-MeanEF)>3.*RMSEF):0.8");
  TObjArray* oaMultGr = new TObjArray(); int igr=0;
  oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, "tpcItsMatchA:run",  "", "(1):(varname_Out==0):(varname_Out):(varname_Warning)", igr) ); igr++;
  oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, "meanTPCncl:run",    "", "(1):(varname_Out==0):(varname_Out):(varname_Warning)", igr) ); igr++;
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  }*/


drawPerformanceTPCQAMatchTrends(const char* inFile = "tree_list", const char* runType="pp") {
  //
  if (gSystem->Exec("ls qaConfig.C")==0)
    gROOT->LoadMacro("qaConfig.C"); //if you run only this macro
  else {
    printf("now loading source/qaConfig.C\n");
    gROOT->LoadMacro("source/qaConfig.C");
  }
  //
  // colors & markers:
  Int_t colPosA=kGreen+2; //kRed;
  Int_t colPosC=kMagenta+1; //kAzure-4;
  Int_t colNegA=kGreen+2; //kRed;//kOrange;
  Int_t colNegC=kMagenta+1; //kAzure-4;//kGreen;
  Int_t colSum=kBlack;
  Int_t marA1=20;
  Int_t marA2=24;
  Int_t marC1=20;
  Int_t marC2=24;
  Int_t marSum=34;//full cross
  Int_t marCorr=31;//snowflake
  // shifting of graphs for one run:
	Float_t sh_gr0=-0.3;
	Float_t sh_gr1=-0.1;
	Float_t sh_gr2=+0.1;
	Float_t sh_gr3=+0.3;

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetLabelSize(0.04,"x");
  gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

  float  ncl_min = 80, ncl_max = 140;
  float  ratio_min = 0.6, ratio_max = 1.2;
  float  mip_min = 30, mip_max = 60;
  float  mipr_min = 0, mipr_max = 0.3;
  float  vx_min = -0.3, vx_max = 0.3;
  float  vy_min = -0.45, vy_max = 0.45;
  float  vz_min = -3, vz_max = 3;
  float  dca_min = -1.2, dca_max = 1.2;
  float  mult_min = 1, mult_max = 35;
  float  pt_min = 0, pt_max = 1.6;
  
  if( strcmp(runType,"PbPb") == 0){
    ncl_min = 90;  ncl_max = 140;
    ratio_min = 0.5;  ratio_max = 1.1;
    mip_min = 35;  mip_max = 70;
    mipr_min = 0;  mipr_max = 0.15;
    vx_min = -0.3;  vx_max = 0.3;
    vy_min = -0.45;  vy_max = 0.45;
    vz_min = -3;  vz_max = 3;
    dca_min = -1;  dca_max = 1;
    mult_min = 5; mult_max = 700;
  }
  if ( strcmp(runType,"PbPbCentr0") == 0){
    ncl_min = 95;  ncl_max = 125;
    ratio_min = 0.5;  ratio_max = 1.1;
    mip_min = 45;  mip_max = 60;
    mipr_min = 0;  mipr_max = 0.15;
    vx_min = -0.3;  vx_max = 0.3;
    vy_min = -0.45;  vy_max = 0.45;
    vz_min = -3;  vz_max = 3;
    dca_min = -0.6;  dca_max = 0.6;
    mult_min = 5; mult_max = 2000;
  }
  if ( strcmp(runType,"PbPbCentr30") == 0){
    ncl_min = 95;  ncl_max = 125;
    ratio_min = 0.5;  ratio_max = 1.1;
    mip_min = 45;  mip_max = 60;
    mipr_min = 0;  mipr_max = 0.15;
    vx_min = -0.3;  vx_max = 0.3;
    vy_min = -0.45;  vy_max = 0.45;
    vz_min = -3;  vz_max = 3;
    dca_min = -0.6;  dca_max = 0.6;
    mult_min = 5; mult_max = 500;
  }
  if ( strcmp(runType,"PbPbCentr70") == 0){
    ncl_min = 95;  ncl_max = 125;
    ratio_min = 0.5;  ratio_max = 1.1;
    mip_min = 45;  mip_max = 60;
    mipr_min = 0;  mipr_max = 0.15;
    vx_min = -0.3;  vx_max = 0.3;
    vy_min = -0.45;  vy_max = 0.45;
    vz_min = -3;  vz_max = 3;
    dca_min = -0.6;  dca_max = 0.6;
    mult_min = 5; mult_max = 50;
  }
  
  // open input file
  //
  TFile *file = TFile::Open(inFile);
  if(!file) return;
  file->cd();
  
  //
  TTree *tree = (TTree*)file->Get("tpcQA");
  //TTree *tree = ch->GetTree();
  if(!tree) return;
  int const entries_tree = tree->GetEntries();
  cout<<"number of tree entries: "<<entries_tree<<endl; 

  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanTPCncl:run","");
  int const entries = gr->GetN();
  cout<<"number of graph entries: "<<entries<<endl;
  cout<<"(multiple occurences of runs removed)"<<endl;

  if(entries<3)
    const float norm_runs = 8.0; 
  else if(entries<20)
    const float norm_runs = 20.0;
  else if(entries<35)
    const float norm_runs = 35.0;
  else 
    const float norm_runs = 50.0;
  // 50 is the max number of runs that shall be viewed on a 1700-wide canvas
  // only when there are more runs, the canvas shall become wider than that.

  int const canvas_width  = int(((entries*1.0)/norm_runs)*1700.0);
  int const canvas_height = 600;
  gStyle->SetPadLeftMargin(0.1*900/canvas_width);
  gStyle->SetPadRightMargin(0.01);
  
  if(entries>50){
    gStyle->SetTickLength(0.03*norm_runs/(entries*1.0),"Y");
    gStyle->SetTitleYOffset((norm_runs/(entries*1.0))*0.8);
    gStyle->SetPadLeftMargin(0.1*norm_runs/(entries*1.0));
    gStyle->SetPadRightMargin(0.01*norm_runs/(entries*1.0));
  }

  TCanvas *c1 = new TCanvas("can","can",canvas_width,canvas_height);
  c1->SetGrid(3);
  c1->cd();

  //
  // compute TPC Status graphs
  TStatToolkit::SetStatusAlias(tree, "meanTPCncl",    "", "varname_Out:(abs(varname-MeanEF)>6.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "meanMIP",       "", "varname_Out:(abs(varname-MeanEF)>6.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "meanVertZ",     "", "varname_Out:(abs(varname-MeanEF)>6.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "tpcItsMatchA",  "", "varname_Out:(abs(varname-MeanEF)>6.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "tpcItsMatchC",  "", "varname_Out:(abs(varname-MeanEF)>6.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "meanTPCncl",    "", "varname_Warning:(abs(varname-MeanEF)>3.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "meanMIP",       "", "varname_Warning:(abs(varname-MeanEF)>3.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "meanVertZ",     "", "varname_Warning:(abs(varname-MeanEF)>3.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "tpcItsMatchA",  "", "varname_Warning:(abs(varname-MeanEF)>3.*RMSEF):0.8");
  TStatToolkit::SetStatusAlias(tree, "tpcItsMatchC",  "", "varname_Warning:(abs(varname-MeanEF)>3.*RMSEF):0.8");
  TObjArray* oaMultGr = new TObjArray();
  int igr=0;
  // the order in the status bar (top->down) will be the opposite of this order:
  oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, "tpcItsMatchC:run",  "", "(1):(varname_Out==0):(varname_Out):(varname_Warning)", igr) ); igr++;
  oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, "tpcItsMatchA:run",  "", "(1):(varname_Out==0):(varname_Out):(varname_Warning)", igr) ); igr++;
  oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, "meanVertZ:run",     "", "(1):(varname_Out==0):(varname_Out):(varname_Warning)", igr) ); igr++;
  oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, "meanMIP:run",       "", "(1):(varname_Out==0):(varname_Out):(varname_Warning)", igr) ); igr++;
  oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, "meanTPCncl:run",    "", "(1):(varname_Out==0):(varname_Out):(varname_Warning)", igr) ); igr++;
  // configure the pad in which they are plotted
  Float_t statPadHeight=0.25; //fraction of canvas height
  Float_t statPadBotMar=0.40; //bottom margin of pad for run numbers
  //


/****** 1/Pt  ******/
// THIS CAN BE USED FOR TESTING!!!
 /*
  TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"deltaPt:run:deltaPt_Err","",1,1,1);
  DrawPlot(gr0, "deltaPt:run", marSum, 1.2, colSum, "AP");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"deltaPtA:run:deltaPtA_Err","",1,1,1, +sh_gr1);
  DrawPlot(gr2, "deltaPtA:run", marA1, 1.0, colPosA, "P");
  TGraphErrors *gr4 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"deltaPtC:run:deltaPtC_Err","",1,1,1, -sh_gr1);
  DrawPlot(gr4, "deltaPtC:run", marC1, 1.0, colPosC, "P");

  gr0->GetHistogram()->SetYTitle("delta (1/pt) ");
  gr0->GetHistogram()->SetTitle(". delta (1/pt)");
  gr0->GetHistogram()->SetMinimum(-0.008);
  gr0->GetHistogram()->SetMaximum(0.008);
  gr0->GetXaxis()->LabelsOption("v");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.9,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("deltaPt:run","deltaPt: both sides","ap");
  leg->AddEntry("deltaPtA:run","deltaPtA: A side only","p");
  leg->AddEntry("deltaPtC:run","deltaPtC: C side only","p");
  leg->Draw();

  PlotThresholds(deltaPt:run); 
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("1overPt_vs_run.png");
  //c1->Clear();

break;
*/



  /****** Number of TPC Clusters vs run number ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanTPCncl:run","",marA1,colPosA,1.0);
  gr->SetName("meanTPCncl:run");
  gr->GetHistogram()->SetYTitle("Number of TPC Clusters");
  gr->GetHistogram()->SetTitle("p_{T} > 0.25GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0");
  gr->GetHistogram()->SetMinimum(ncl_min);
  gr->GetHistogram()->SetMaximum(ncl_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
 
	PlotThresholds(meanTPCncl:run);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanTPCncl_vs_run.png");
  c1->Clear();

  /****** Ratio of findable TPC clusters vs run number ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanTPCnclF:run","",marA1,colPosA,1.0);
  gr->SetName("meanTPCnclF:run");
  gr->GetHistogram()->SetYTitle("# of Found Clusters/ # of Findable Clusters");
  gr->GetHistogram()->SetTitle("p_{T} > 0.25GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0");
  gr->GetHistogram()->SetMinimum(ratio_min);
  gr->GetHistogram()->SetMaximum(ratio_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  
	PlotThresholds(meanTPCnclF:run);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanTPCnclF_vs_run.png");
  c1->Clear();


  /****** Mean MIPs ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanMIP:run","",marA1,colPosA,1.0);
  gr->SetName("meanMIP:run");
  gr->GetHistogram()->SetYTitle("Mean of MIPs");
  gr->GetHistogram()->SetTitle("0,4<p<0.55GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0, 80<#Cluster<160, 35<dE/dx<60");
  gr->GetHistogram()->SetMinimum(mip_min);
  gr->GetHistogram()->SetMaximum(mip_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  
	PlotThresholds(meanMIP:run);  
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanMIP_vs_run.png");
  c1->Clear();

  /****** Mean MIP Resolution ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"resolutionMIP:run","",marA1,colPosA,1.0);
  gr->SetName("resolutionMIP:run");
  gr->GetHistogram()->SetYTitle("Resolution of MIPs");
  gr->GetHistogram()->SetTitle("0,4<p<0.55GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0, 80<#Cluster<160, 35<dE/dx<60");
  gr->GetHistogram()->SetMinimum(mipr_min);
  gr->GetHistogram()->SetMaximum(mipr_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  
  PlotThresholds(resolutionMIP:run); 
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("resolutionMIP_vs_run.png");
  c1->Clear();
  
  /****** Mean energy loss for electrons ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanMIPele:run","",marA1,colPosA,1.0);
  gr->SetName("meanMIPele:run");
  gr->GetHistogram()->SetYTitle("Mean of electron dEdx");
  gr->GetHistogram()->SetTitle("0,32<p<0.38GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0, 80<#Cluster<160, 70<dE/dx<100");
  gr->GetHistogram()->SetMinimum(40);
  gr->GetHistogram()->SetMaximum(110);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  //PlotThresholds(meanMIP:run);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meandEdxele_vs_run.png");
  c1->Clear();

  /****** Mean Energy loss electron Resolutio ******/
 
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"resolutionMIPele:run","",marA1,colPosA,1.0);
  gr->SetName("resolutionMIPele:run");
  gr->GetHistogram()->SetYTitle("Resolution of electrons dEdx");
  //gr->GetHistogram()->SetTitle("0,4<p<0.55GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0, 80<#Cluster<160, 35<dE/dx<60");
  gr->GetHistogram()->SetTitle("0,32<p<0.38GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0, 80<#Cluster<160, 70<dE/dx<100");
  gr->GetHistogram()->SetMinimum(mipr_min);
  gr->GetHistogram()->SetMaximum(mipr_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  //PlotThresholds(resolutionMIPele:run);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("resolutionMeandEdxEle_vs_run.png");
  c1->Clear();
  
////////////////////////////////////////////////////////////////////////////////////////////////
  
  

  /****** Mean VertX ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanVertX:run","",marA1,colPosA,1.0);
  gr->SetName("meanVertX:run");
  gr->GetHistogram()->SetYTitle("Mean of Vert_{X} / [cm]");
  gr->GetHistogram()->SetTitle("");
  gr->GetHistogram()->SetMinimum(vx_min);
  gr->GetHistogram()->SetMaximum(vx_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  
  PlotThresholds(meanVertX:run);  
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanVertX_vs_run.png");
  c1->Clear();

  /****** Mean VertY ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanVertY:run","",marA1,colPosA,1.0);
  gr->SetName("meanVertY:run");
  gr->GetHistogram()->SetYTitle("Mean of Vert_{Y} / [cm]");
  gr->GetHistogram()->SetTitle("");
  gr->GetHistogram()->SetMinimum(vy_min);
  gr->GetHistogram()->SetMaximum(vy_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  
  PlotThresholds(meanVertY:run); 
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanVertY_vs_run.png");
  c1->Clear();

  /****** Mean VertZ ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanVertZ:run","",marA1,colPosA,1.0);
  gr->SetName("meanVertZ:run");
  gr->GetHistogram()->SetYTitle("Mean of Vert_{Z} / [cm]");
  gr->GetHistogram()->SetTitle("");
  gr->GetHistogram()->SetMinimum(vz_min);
  gr->GetHistogram()->SetMaximum(vz_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  
  PlotThresholds(meanVertZ:run);  
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanVertZ_vs_run.png");
  c1->Clear();
  

  /****** Offset DCA  ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"offsetdRA:run:offsetdRAErr","",marA1,colPosA,1.0,sh_gr0);
	gr->SetName("offsetdRA:run");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"offsetdZA:run:offsetdZAErr","",marA2,colNegA,1.0,sh_gr1);
  gr1->SetName("gr1");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"offsetdRC:run:offsetdRCErr","",marC1,colPosC,1.0,sh_gr2);
  gr2->SetName("gr2");
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"offsetdZC:run:offsetdZCErr","",marC2,colNegC,1.0,sh_gr3);
  gr3->SetName("gr3");
  
  gr->GetHistogram()->SetYTitle("DCAs / [cm]");
  gr->GetHistogram()->SetTitle("p_{T} > 0.25GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 0.8");
  gr->GetHistogram()->SetMinimum(dca_min);
  gr->GetHistogram()->SetMaximum(dca_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  gr2->Draw("P");
  gr3->Draw("P");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("offsetdRA:run","DCA_{R}, A Side","p");
  leg->AddEntry("gr1","DCA_{Z}, A Side","p");
  leg->AddEntry("gr2","DCA_{R}, C Side","p");
  leg->AddEntry("gr3","DCA_{Z}, C Side","p");
  leg->Draw();

  PlotThresholds(offsetdRA:run);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
	c1->SaveAs("DCAOffset_vs_run.png");
  c1->Clear();


  /****** Mean Mult  ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanMultPos:run:rmsMultPos","",marA1,colPosA,1.0,sh_gr1);
  gr->SetName("meanMultPos:run");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanMultNeg:run:rmsMultPos","",marA2,colNegA,1.0,sh_gr2);
  gr1->SetName("gr1");

  gr->GetHistogram()->SetYTitle("Multiplicites of Primary Tracks");
  gr->GetHistogram()->SetTitle("|DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, #Cluster > 70");
  //gr->GetHistogram()->SetMinimum(mult_min);
  //gr->GetHistogram()->SetMaximum(mult_max);
  gr->GetHistogram()->SetMinimum(0);
  gr->GetHistogram()->SetMaximum(100);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  
  TLegend *leg = new TLegend(0.6,0.85,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("meanMultPos:run","Positive Charged Tracks","p");
  leg->AddEntry("gr1","Negative Charged Tracks","p");
  leg->Draw();

  PlotThresholds(meanMultPos:run);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanMult_vs_run.png");
  c1->Clear();
  
  /****** TPC-ITS matching Efficiency  ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"tpcItsMatchA:run","",marA1,colPosA,1.0,sh_gr0);
  gr->SetName("tpcItsMatchA:run");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"tpcItsMatchHighPtA:run","",marA2,colNegA,1.0,sh_gr1);
  gr1->SetName("gr1");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"tpcItsMatchC:run","",marC1,colPosC,1.0,sh_gr2);
  gr2->SetName("gr2");
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"tpcItsMatchHighPtC:run","",marC2,colNegC,1.0,sh_gr3);
  gr3->SetName("gr3");

  gr->GetHistogram()->SetYTitle("Matching Efficiencies");
  gr->GetHistogram()->SetTitle("");
  gr->GetHistogram()->SetMinimum(0.0);
  gr->GetHistogram()->SetMaximum(1.2);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  gr2->Draw("P");
  gr3->Draw("P");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("tpcItsMatchA:run","TPC-ITS matching, A Side","p");
  leg->AddEntry("gr1","TPC-ITS matching ( Pt>4GeV/c ), A Side","p");
  leg->AddEntry("gr2","TPC-ITS matching, C Side","p");
  leg->AddEntry("gr3","TPC-ITS matching ( Pt>4GeV/c ), C Side","p");
  leg->Draw();
  
  PlotThresholds(tpcItsMatchA:run);   
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("TPC-ITS_vs_run.png");
  c1->Clear();

  /****** ITS-TPC matching Efficiency  ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"lambdaPull:run","",marA1,colPosA,1.0,sh_gr0);
  gr->SetName("lambdaPull:run");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ptPull:run","",marA2,colNegA,1.0,sh_gr1);
  gr1->SetName("gr1");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"yPull:run","",marC1,colPosC,1.0,sh_gr2);
  gr2->SetName("gr2");
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"zPull:run","",marC2,colNegC,1.0,sh_gr3);
  gr3->SetName("gr3");

  gr->GetHistogram()->SetYTitle("Pulls");
  gr->GetHistogram()->SetTitle("");
  gr->GetHistogram()->SetMinimum(-3);
  gr->GetHistogram()->SetMaximum(3);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  gr2->Draw("P");
  gr3->Draw("P");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("lambdaPull:run","Tan#lambda pull","p");
  leg->AddEntry("gr1","Pt pull","p");
  leg->AddEntry("gr2","y pull","p");
  leg->AddEntry("gr3","z pull","p");
  leg->Draw();

  PlotThresholds(lambdaPull:run);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
	c1->SaveAs("ITS-TPC_vs_run.png");
  c1->Clear();


  /****** pullPhi for TPC Constrain  ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"tpcConstrainPhiA:run","",marA1,colPosA,1.0,sh_gr1);
  gr->SetName("tpcConstrainPhiA:run");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"tpcConstrainPhiC:run","",marC1,colPosC,1.0,sh_gr2);
  gr1->SetName("gr1");

  gr->GetHistogram()->SetYTitle("(sin#phi_{TPC} - sin#phi_{Global})/#sigma");
  gr->GetHistogram()->SetTitle("");
  gr->GetHistogram()->SetMinimum(-1);
  gr->GetHistogram()->SetMaximum(1);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.9,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("tpcConstrainPhiA:run","A Side","p");
  leg->AddEntry("gr1","C Side","p");
  leg->Draw();
  
  PlotThresholds(tpcConstrainPhiA:run); 
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("pullPhiConstrain_vs_run.png");
  c1->Clear();

  
	/****** 1/Pt  ******/
  //  TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"deltaPt:run","");
  TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"deltaPt:run:deltaPt_Err","",1,1,1);
  DrawPlot(gr0, "deltaPt:run", marSum, 1.2, colSum, "AP");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"deltaPtA:run:deltaPtA_Err","",1,1,1, +sh_gr1);
  DrawPlot(gr2, "deltaPtA:run", marA1, 1.0, colPosA, "P");
  TGraphErrors *gr4 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"deltaPtC:run:deltaPtC_Err","",1,1,1, -sh_gr1);
  DrawPlot(gr4, "deltaPtC:run", marC1, 1.0, colPosC, "P");

  gr0->GetHistogram()->SetYTitle("delta (q/pt) ");
  gr0->GetHistogram()->SetTitle("delta (q/pt)");
  gr0->GetHistogram()->SetMinimum(-0.008);
  gr0->GetHistogram()->SetMaximum(0.008);
  gr0->GetXaxis()->LabelsOption("v");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.9,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("deltaPt:run","deltaPt: both sides","ap");
  leg->AddEntry("deltaPtA:run","deltaPtA: A side only","p");
  leg->AddEntry("deltaPtC:run","deltaPtC: C side only","p");
  leg->Draw();

  PlotThresholds(deltaPt:run); 
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("1overPt_vs_run.png");
  c1->Clear();


  /****** DCAr fitting parameters  ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcarAP0:run","",marA1,colPosA,1.0,sh_gr0);
  gr->SetName("dcarAP0:run");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcarAP1:run","",marA2,colNegA,1.0,sh_gr1);
  gr1->SetName("gr1");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcarCP0:run","",marC1,colPosC,1.0,sh_gr2);
  gr2->SetName("gr2");
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcarCP1:run","",marC2,colNegC,1.0,sh_gr3);
  gr3->SetName("gr3");

  gr->GetHistogram()->SetYTitle("DCAR Fitting Parameters");
  gr->GetHistogram()->SetTitle("sqrt(P0^{2} + P1^{2}/(pT^{2}))");
  gr->GetHistogram()->SetMinimum(-1);
  gr->GetHistogram()->SetMaximum(1);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  gr2->Draw("P");
  gr3->Draw("P");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcarAP0:run","P0, A Side","p");
  leg->AddEntry("gr1","P1, A Side","p");
  leg->AddEntry("gr2","P0, C Side","p");
  leg->AddEntry("gr3","P1, C Side","p");
  leg->Draw();

  PlotThresholds(dcarAP0:run);  
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("dcar_fitting_run.png");
  c1->Clear();
 
  ////////////////////////////////////////////////////////////////////////
  //test DCAR plots
  //DCAr first parameter _0

  TCanvas *c2  = new  TCanvas("can2","can2",canvas_width,canvas_height); 
  c2->cd();
  c2->Update();
  c2->SetGrid(3);
   
  //TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_posA_0:run","",1,1,1,sh_gr0);
  TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_posA_0:run:dcar_posA_0_Err","",1,1,1,sh_gr0);
  DrawPlot(gr0, "dcar_posA_0:run", marA1, 1.0, colPosA, "AP");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_negA_0:run:dcar_negA_0_Err","",1,1,1,sh_gr1);
  DrawPlot(gr1, "dcar_negA_0:run", marA2, 1.0, colNegA, "P");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_posC_0:run:dcar_posC_0_Err","",1,1,1,sh_gr2);
  DrawPlot(gr2, "dcar_posC_0:run", marC1, 1.0, colPosC, "P");    
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_negC_0:run:dcar_negC_0_Err","",1,1,1,sh_gr3);
  DrawPlot(gr3, "dcar_negC_0:run", marC2, 1.0, colNegC, "P");    

  gr0->GetHistogram()->SetYTitle("DCARs");
  gr0->GetHistogram()->SetMinimum(-0.2);
  gr0->GetHistogram()->SetMaximum(0.2);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcar_0:run");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.9,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcar_0:run","dcar_posA_0","p");
  leg->AddEntry("dcar_negA_0:run","dcar_negA_0","p");
  leg->AddEntry("dcar_posC_0:run","dcar_posC_0","p");
  leg->AddEntry("dcar_negC_0:run","dcar_negC_0","p");
  leg->Draw();
  
  PlotThresholds(dcar_0:run);
  TStatToolkit::AddStatusPad(c2, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c2->SaveAs("dcar_0_vs_run.png");//for C,A side and pos/neg particle
  c2->Update();

  ////////////////////////////////////////////////////////////////////////
  //DCAr second parameter _1
  
  TCanvas *c3  = new  TCanvas("can3","can3",canvas_width,canvas_height); 
  c3->cd();
  c3->Update();
  c3->SetGrid(3);
  
  TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_posA_1:run:dcar_posA_1_Err","",1,1,1,sh_gr0);
  DrawPlot(gr0, "dcar_posA_1:run", marA1, 1.0, colPosA, "AP");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_negA_1:run:dcar_negA_1_Err","",1,1,1,sh_gr1);
  DrawPlot(gr1, "dcar_negA_1:run", marA2, 1.0, colNegA, "P");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_posC_1:run:dcar_posC_1_Err","",1,1,1,sh_gr2);
  DrawPlot(gr2, "dcar_posC_1:run", marC1, 1.0, colPosC, "P");    
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_negC_1:run:dcar_negC_1_Err","",1,1,1,sh_gr3);
  DrawPlot(gr3, "dcar_negC_1:run", marC2, 1.0, colNegC, "P");    
  
  gr0->GetHistogram()->SetYTitle("DCARs");
  gr0->GetHistogram()->SetMinimum(-0.1);
  gr0->GetHistogram()->SetMaximum(0.1);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcar_1:run");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.9,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcar_1:run","dcar_posA_1","p");
  leg->AddEntry("dcar_negA_1:run","dcar_negA_1","p");
  leg->AddEntry("dcar_posC_1:run","dcar_posC_1","p");
  leg->AddEntry("dcar_negC_1:run","dcar_negC_1","p");
  leg->Draw();
  
  PlotThresholds(dcar_1:run);
  TStatToolkit::AddStatusPad(c3, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c3->SaveAs("dcar_1_vs_run.png");//for C,A side and pos/neg particle
  c3->Update();
  
  ////////////////////////////////////////////////////////////////////////
  //DCAr third parameter _2
  
  TCanvas *c5  = new  TCanvas("can5","can5",canvas_width,canvas_height); 
  c5->cd();
  c5->Update();
  c5->SetGrid(3);
   
  TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_posA_2:run:dcar_posA_2_Err","",1,1,1,sh_gr0);
  DrawPlot(gr0, "dcar_posA_2:run", marA1, 1.0, colPosA, "AP");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_negA_2:run:dcar_negA_2_Err","",1,1,1,sh_gr1);
  DrawPlot(gr1, "dcar_negA_2:run", marA2, 1.0, colNegA, "P");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_posC_2:run:dcar_posC_2_Err","",1,1,1,sh_gr2);
  DrawPlot(gr2, "dcar_posC_2:run", marC1, 1.0, colPosC, "P");    
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar_negC_2:run:dcar_negC_2_Err","",1,1,1,sh_gr3);
  DrawPlot(gr3, "dcar_negC_2:run", marC2, 1.0, colNegC, "P");    
  
  gr0->GetHistogram()->SetYTitle("DCARs");
  gr0->GetHistogram()->SetMinimum(-0.1);
  gr0->GetHistogram()->SetMaximum(0.1);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcar_2:run");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.9,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcar_2:run","dcar_posA_2","p");
  leg->AddEntry("dcar_negA_2:run","dcar_negA_2","p");
  leg->AddEntry("dcar_posC_2:run","dcar_posC_2","p");
  leg->AddEntry("dcar_negC_2:run","dcar_negC_2","p");
  leg->Draw();
  
  PlotThresholds(dcar_2:run);
  TStatToolkit::AddStatusPad(c5, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c5->SaveAs("dcar_2_vs_run.png");//for C,A side and pos/neg particle
  c5->Update();

  ////////////////////////////////////////////////////////////////////////
  //DCAz parameters
  //Dcaz first parameter _0

  TCanvas *c6  = new  TCanvas("can6","can6",canvas_width,canvas_height); 
  c6->cd();
  c6->Update();
  c6->SetGrid(3);
   
  TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_posA_0:run:dcaz_posA_0_Err","",1,1,1,sh_gr0);
  DrawPlot(gr0, "dcaz_posA_0:run", marA1, 1.0, colPosA, "AP");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_negA_0:run:dcaz_negA_0_Err","",1,1,1,sh_gr1);
  DrawPlot(gr1, "dcaz_negA_0:run", marA2, 1.0, colNegA, "P");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_posC_0:run:dcaz_posC_0_Err","",1,1,1,sh_gr2);
  DrawPlot(gr2, "dcaz_posC_0:run", marC1, 1.0, colPosC, "P");    
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_negC_0:run:dcaz_negC_0_Err","",1,1,1,sh_gr3);
  DrawPlot(gr3, "dcaz_negC_0:run", marC2, 1.0, colNegC, "P");    
 
  gr0->GetHistogram()->SetYTitle("DCAZs");
  gr0->GetHistogram()->SetMinimum(-2.);
  gr0->GetHistogram()->SetMaximum(2.);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcaz_0:run");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.9,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcaz_0:run","dcaz_posA_0","p");
  leg->AddEntry("dcaz_negA_0:run","dcaz_negA_0","p");
  leg->AddEntry("dcaz_posC_0:run","dcaz_posC_0","p");
  leg->AddEntry("dcaz_negC_0:run","dcaz_negC_0","p");
  leg->Draw();
  
  PlotThresholds(dcaz_0:run);
  TStatToolkit::AddStatusPad(c6, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c6->SaveAs("dcaz_0_vs_run.png");//for C,A side and pos/neg particle
  c6->Update();

  ////////////////////////////////////////////////////////////////////////
  //Dcaz second parameter _1
  
  TCanvas *c7  = new  TCanvas("can7","can7",canvas_width,canvas_height); 
  c7->cd();
  c7->Update();
  c7->SetGrid(3);

  TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_posA_1:run:dcaz_posA_1_Err","",1,1,1,sh_gr0);
  DrawPlot(gr0, "dcaz_posA_1:run", marA1, 1.0, colPosA, "AP");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_negA_1:run:dcaz_negA_1_Err","",1,1,1,sh_gr1);
  DrawPlot(gr1, "dcaz_negA_1:run", marA2, 1.0, colNegA, "P");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_posC_1:run:dcaz_posC_1_Err","",1,1,1,sh_gr2);
  DrawPlot(gr2, "dcaz_posC_1:run", marC1, 1.0, colPosC, "P");    
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_negC_1:run:dcaz_negC_1_Err","",1,1,1,sh_gr3);
  DrawPlot(gr3, "dcaz_negC_1:run", marC2, 1.0, colNegC, "P");    
 
  gr0->GetHistogram()->SetYTitle("DCAZs");
  gr0->GetHistogram()->SetMinimum(-0.2);
  gr0->GetHistogram()->SetMaximum(0.2);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcaz_1:run");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.9,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcaz_1:run","dcaz_posA_1","p");
  leg->AddEntry("dcaz_negA_1:run","dcaz_negA_1","p");
  leg->AddEntry("dcaz_posC_1:run","dcaz_posC_1","p");
  leg->AddEntry("dcaz_negC_1:run","dcaz_negC_1","p");
  leg->Draw();
  
  PlotThresholds(dcaz_1:run);
  TStatToolkit::AddStatusPad(c7, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c7->SaveAs("dcaz_1_vs_run.png");//for C,A side and pos/neg particle
  c7->Update();
  
  ////////////////////////////////////////////////////////////////////////
  //Dcaz third parameter _2
  
  TCanvas *c8  = new  TCanvas("can8","can8",canvas_width,canvas_height); 
  c8->cd();
  c8->Update();
  c8->SetGrid(3);

  TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_posA_2:run:dcaz_posA_2_Err","",1,1,1,sh_gr0);
  DrawPlot(gr0, "dcaz_posA_2:run", marA1, 1.0, colPosA, "AP");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_negA_2:run:dcaz_negA_2_Err","",1,1,1,sh_gr1);
  DrawPlot(gr1, "dcaz_negA_2:run", marA2, 1.0, colNegA, "P");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_posC_2:run:dcaz_posC_2_Err","",1,1,1,sh_gr2);
  DrawPlot(gr2, "dcaz_posC_2:run", marC1, 1.0, colPosC, "P"); 
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz_negC_2:run:dcaz_negC_2_Err","",1,1,1,sh_gr3);
  DrawPlot(gr3, "dcaz_negC_2:run", marC2, 1.0, colNegC, "P");

  gr0->GetHistogram()->SetYTitle("DCAZs");
  gr0->GetHistogram()->SetMinimum(-0.1);
  gr0->GetHistogram()->SetMaximum(0.1);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcaz_2:run");
  
  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.9,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcaz_2:run","dcaz_posA_2","p");
  leg->AddEntry("dcaz_negA_2:run","dcaz_negA_2","p");
  leg->AddEntry("dcaz_posC_2:run","dcaz_posC_2","p");
  leg->AddEntry("dcaz_negC_2:run","dcaz_negC_2","p");
  leg->Draw();
  
  PlotThresholds(dcaz_2:run);
  TStatToolkit::AddStatusPad(c8, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c8->SaveAs("dcaz_2_vs_run.png");//for C,A side and pos/neg particle
  c8->Update();

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // Plot Occupancy IROC, OROC, A,C side
  
  TCanvas *c9  = new  TCanvas("can9","can9",canvas_width,canvas_height);
  c9->cd();
  c9->Update();
  c9->SetGrid(3);

  TGraphErrors *gr0 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"iroc_A_side:run","",marA1,colPosA,1.0,sh_gr0);
  gr0->SetName("iroc_A_side:run");

  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"oroc_A_side:run","",marA2,colNegA,1.0,sh_gr1);
  gr1->SetName("gr1");

  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"iroc_C_side:run","",marC1,colPosC,1.0,sh_gr2);
  gr2->SetName("gr2");

  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"oroc_C_side:run","",marC2,colNegC,1.0,sh_gr3);
  gr3->SetName("gr3");

  gr0->GetHistogram()->SetYTitle("(nr_Chamber) - (nr_Chamber_lowOcc)");
  gr0->GetHistogram()->SetMinimum(14.0);
  gr0->GetHistogram()->SetMaximum(20.0);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("occ_AC_Side_IROC_OROC:run");
  
  gr0->Draw("AP");
  gr1->Draw("P");
  gr2->Draw("P");
  gr3->Draw("P");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.9,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("occ_AC_Side_IROC_OROC:run","iroc_A_side","p");
  leg->AddEntry("gr1","oroc_A_side","p");
  leg->AddEntry("gr2","iroc_C_side","p");
  leg->AddEntry("gr3","oroc_C_side","p");
  leg->Draw();

  TStatToolkit::AddStatusPad(c9, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c9->SaveAs("occ_AC_Side_IROC_OROC_vs_run.png");//for C,A side and IROC,OROC                                                                                                 
  c9->Update();
                                                                                                                                                 
  /****** attachemnt parameters for A and C side ******/
  
  TCanvas *c10  = new  TCanvas("can10","can10",canvas_width,canvas_height);
  c10->cd();
  c10->Update();
  c10->SetGrid(3);
  
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"MIPattachSlopeA:run","",marA1,colPosA,1.0);
  gr->SetName("MIPattachSlopeA:run");
  gr->GetHistogram()->SetYTitle("Attachment parameter p1");
  gr->GetHistogram()->SetMinimum(-20);
  gr->GetHistogram()->SetMaximum(20);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  TStatToolkit::AddStatusPad(c10, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c10->SaveAs("MIPattachSlopeA_vs_run.png");
  c10->Clear();

  //C side
  TCanvas *c11  = new  TCanvas("can11","can11",canvas_width,canvas_height);
  c11->cd();
  c11->Update();
  c11->SetGrid(3);

  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"MIPattachSlopeC:run","",marA1,colPosA,1.0);
  gr->SetName("MIPattachSlopeC:run");
  gr->GetHistogram()->SetYTitle("Attachment parameter p1");
  gr->GetHistogram()->SetMinimum(-20);
  gr->GetHistogram()->SetMaximum(20);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  TStatToolkit::AddStatusPad(c11, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c11->SaveAs("MIPattachSlopeC_vs_run.png");
  c11->Clear();

  /****** electron and MIPs separation ******/

  TCanvas *c12  = new  TCanvas("can12","can12",canvas_width,canvas_height);
  c12->cd();
  c12->Update();
  c12->SetGrid(3);

  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"electroMIPSeparation:run","",marA1,colPosA,1.0);
  gr->SetName("electroMIPSeparation:run");
  gr->GetHistogram()->SetYTitle("Electron - MIP");
  gr->GetHistogram()->SetMinimum(0);
  gr->GetHistogram()->SetMaximum(120);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  TStatToolkit::AddStatusPad(c12, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c12->SaveAs("ElectroMIPSeparation_vs_run.png");
  c12->Clear();
  
}

//the function plots thresholds
Int_t PlotThresholds(TGraphErrors* h)
{
  TMap* configMap = qaConfig();
 
  Float_t min = ConfigEntryMin(configMap,h);
  Float_t max = ConfigEntryMax(configMap,h);
  
  cout<<"min, max "<< min << " "<< max << endl;
  
  TString desc = ConfigEntryDescription(configMap,h);
  TAxis *xaxis = h->GetXaxis();
  Float_t x1 = xaxis->GetBinLowEdge(1);
  Float_t x2 = xaxis->GetBinUpEdge(xaxis->GetLast());
  TLine *lineMin = new TLine(x1,min,x2,min); lineMin->SetLineColor(kRed);
  lineMin->Draw();
  TLine *lineMax = new TLine(x1,max,x2,max); lineMax->SetLineColor(kRed);
  lineMax->Draw();

  return 1;
}

//the function draws the plots
Int_t DrawPlot(TGraphErrors* gr, TString nameHisto, Int_t markerStyle, Int_t markerSize, Int_t markerColor, TString drawMode)
{
  gr->SetName(nameHisto);
  gr->SetMarkerStyle(markerStyle);
  gr->SetMarkerSize(markerSize);
  gr->SetMarkerColor(markerColor);
  gr->Draw(drawMode);
  return 1;
}


