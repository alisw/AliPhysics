
// 05.2013 new functionality: Status Bar
// 07.2014 updated to steer most of this from the qaConfig.C
// 09.2014 new functions to produce status lines from aliases and to write status infos to a tree for later use (see TStatToolkit).
// To create the Status Bar, the following is done in principle. more info & examples in the TStatToolkit and qaConfig.C.
/*{
 TStatToolkit::SetStatusAlias(tree, "meanTPCncl",    "", "varname_Out:(abs(varname-MeanEF)>6.*RMSEF):0.8");
 TStatToolkit::SetStatusAlias(tree, "tpcItsMatchA",  "", "varname_Out:(abs(varname-MeanEF)>6.*RMSEF):0.8");
 TStatToolkit::SetStatusAlias(tree, "meanTPCncl",    "", "varname_Warning:(abs(varname-MeanEF)>3.*RMSEF):0.8");
 TStatToolkit::SetStatusAlias(tree, "tpcItsMatchA",  "", "varname_Warning:(abs(varname-MeanEF)>3.*RMSEF):0.8");
 TObjArray* oaMultGr = new TObjArray(); int igr=0;
 oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, "tpcItsMatchA:run",  "", "(1):(meanTPCncl>0):(varname_Warning):(varname_Outlier):", igr) ); igr++;
 oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, "meanTPCncl:run",    "", "(1):(meanTPCncl>0):(varname_Warning):(varname_Outlier):", igr) ); igr++;
 TCanvas *c1 = new TCanvas("c1","c1");
 TStatToolkit::AddStatusPad(c1, 0.30, 0.40);
 TStatToolkit::DrawStatusGraphs(oaMultGr);
 }*/


TTree *tree;
TTree *statusTree;

drawPerformanceTPCQAMatchTrends(const char* inFile = "trending.root", const char* runType="pp") {
  //

  if (gSystem->Exec("ls qaConfig.C")==0)
    gROOT->LoadMacro(  "qaConfig.C");
  else {
    printf("now loading $ALICE_PHYSICS/PWGPP/TPC/macros/qaConfig.C\n");
    gROOT->LoadMacro(  "$ALICE_PHYSICS/PWGPP/TPC/macros/qaConfig.C");
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
  // shifting of graphs within one run for better visibility:
  Float_t sh_gr0=-0.3;
  Float_t sh_gr1=-0.1;
  Float_t sh_gr2=+0.1;
  Float_t sh_gr3=+0.3;
  // properties of status lines:
  // currently set in 'MakeStatusLines()'

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

  // make backup of rootfile
  //
  TString sBackupfile(inFile);
  sBackupfile.ReplaceAll(".root",".backup.root");
  gSystem->Exec(Form("cp %s %s", inFile, sBackupfile.Data()));

  // open input file
  //
  TFile *_file0 = TFile::Open(inFile, "UPDATE");
  if(!_file0) return;
  _file0->cd();

  //
  tree = (TTree*)_file0->Get("tpcQA");
  if (!tree) tree = (TTree*)_file0->Get("trending");
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

  double SpaceForLegend = 0.;
  int const canvas_width  = SpaceForLegend + int(((entries*1.0)/norm_runs)*1700.0);
  int const canvas_height = 600;
  gStyle->SetPadLeftMargin(0.12*900/canvas_width);
  // gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadRightMargin(SpaceForLegend/canvas_width); // SpaceForLegend pixels of whitespace right next to the histograms

  if(entries>50){
    gStyle->SetTickLength(0.03*norm_runs/(entries*1.0),"Y");
    gStyle->SetTitleYOffset((norm_runs/(entries*1.0))*0.8);
    gStyle->SetPadLeftMargin(0.12*norm_runs/(entries*1.0));
    gStyle->SetPadRightMargin(0.01*norm_runs/(entries*1.0));
  }

  TCanvas *c1 = new TCanvas("can","can",canvas_width,canvas_height);
  c1->SetGrid(3);
  c1->cd();
  gPad->SetTicks(1,2);
  // c1->SetRightMargin(SpaceForLegend/canvas_width);
  double leftlegend  = 1 - 180./c1->GetWw();
  double rightlegend = 1 - 10./c1->GetWw();

  //
  // process config file qaConfig.C to initialize status aliases (outliers etc.), status bar criteria, status lines, ...
  //
  TString returnStrings[3];
  qaConfig(tree, returnStrings);
  // configures outlier criteria and descriptions for the needed TPC variables, as specified in the qaConfig.C.
  // defines aliases according to these criteria.

  TString sStatusbarVars  = returnStrings[0];
  TString sStatusbarNames = returnStrings[1];
  TString sCriteria       = returnStrings[2];
  cout << "sStatusbarVars = " << sStatusbarVars.Data() << endl;
  cout << "sCriteria      = " << sCriteria.Data() << endl;

  //
  // compute TPC status graphs
  //
  TObjArray* oaStatusbarVars = sStatusbarVars.Tokenize(";");
  TObjArray* oaStatusbarNames = sStatusbarNames.Tokenize(";");
  TObjArray* oaMultGr = new TObjArray();
  int igr=0;

  for (Int_t vari=oaStatusbarVars->GetEntriesFast()-1; vari>=0; vari--) // invert the order of the status graphs
  {
    TString sVar = Form("%s:run", oaStatusbarVars->At(vari)->GetName()); //e.g. -> dcar:run
    oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, sVar.Data(),  "", sCriteria.Data(), igr) );
    TString sYtitle = oaStatusbarNames->At(vari)->GetName(); // set better name for y axis of statuspad
    ((TMultiGraph*) oaMultGr->At(igr))->SetTitle(sYtitle.Data());
    igr++;
  }

  //
  // save status into Tree and write to rootfile
  // we update the original rootfile trending.root, as it is complicated to 'copy-paste' a TTree...
  //
  statusTree = TStatToolkit::WriteStatusToTree(oaMultGr);
  statusTree->BuildIndex("run");
  tree->AddFriend(statusTree,"Tstatus");
//  tree->Write("", TObject::kOverwrite);
//  statusTree->Write();
  // if we save statusTree to file here, then the run number in the plots will be always the same run. no idea why.
  // so we do it at the end...
  //
  // alternative: write statusTree to different rootfile: (same problem)
//  TFile* file_out = new TFile("trendingStatusTree.root","RECREATE");
//  file_out->cd();
//  statusTree->Write();
//  file_out->Close();
//  Printf("Status tree written to file '%s'", file_out->GetName());

  //afterwards one can open the rootfile and correlate the trees:
  /*
   // read the trees and draw tests:
   // [terminal]$ aliroot -l trending.root
   TTree* tree = (TTree*)_file0->Get("tpcQA");
   tree->Draw("meanMIP:run","run>0","*");
   TTree* statusTree = (TTree*)_file0->Get("statusTree");
   statusTree->Draw("MIPquality_Warning:run","run>0","*");

   // correlate:
   tree->Draw("meanMIP:Tstatus.MIPquality_Warning","run>0","*");

   TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"Tstatus.MIPquality_Warning:run","",20,kRed,1.0);
   gr->Draw("AP");
  */


  cout << "Start plotting of trending graphs... " << endl;
  // configure the pad in which the status graphs are plotted ('status bar')
  Float_t statPadHeight=0.30; //fraction of canvas height (was 0.25 until Aug 2014)
  Float_t statPadBotMar=0.40; //bottom margin of pad for run numbers
  //
  // automatic plot ranges based on outlier bands, computed for each variable later.
  Float_t plotmean;
  Float_t plotoutlier;
  //
  c1->cd();

  /****** Number of TPC Clusters vs run number ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanTPCncl:run","",marA1,colPosA,1.0);
  gr->SetName("meanTPCncl:run");
  gr->GetHistogram()->SetYTitle("Number of TPC Clusters");
  gr->GetHistogram()->SetTitle("p_{T} > 0.25GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0");
  ComputeRange(tree, "meanTPCncl", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
//  gr->GetHistogram()->SetMinimum(ncl_min);
//  gr->GetHistogram()->SetMaximum(ncl_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"meanTPCncl:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanTPCncl_vs_run.png");
  c1->Clear();

  /****** Ratio of findable TPC clusters vs run number ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanTPCnclF:run","",marA1,colPosA,1.0);
  gr->SetName("meanTPCnclF:run");
  gr->GetHistogram()->SetYTitle("# of Found Clusters/ # of Findable Clusters");
  gr->GetHistogram()->SetTitle("p_{T} > 0.25GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0");
  ComputeRange(tree, "meanTPCnclF", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
//  gr->GetHistogram()->SetMinimum(ratio_min);
//  gr->GetHistogram()->SetMaximum(ratio_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"meanTPCnclF:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanTPCnclF_vs_run.png");
  c1->Clear();

  /****** Mean MIPs ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanMIP:run","",marA1,colPosA,1.0);
  gr->SetName("meanMIP:run");
  gr->GetHistogram()->SetYTitle("Mean of MIPs");
  gr->GetHistogram()->SetTitle("0,4<p<0.55GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0, 80<#Cluster<160, 35<dE/dx<60");
  ComputeRange(tree, "meanMIP", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
//  gr->GetHistogram()->SetMinimum(mip_min);
//  gr->GetHistogram()->SetMaximum(mip_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"meanMIP:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanMIP_vs_run.png");
  c1->Clear();


  /****** Mean MIP Resolution ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"resolutionMIP:run","",marA1,colPosA,1.0);
  gr->SetName("resolutionMIP:run");
  gr->GetHistogram()->SetYTitle("Resolution of MIPs");
  gr->GetHistogram()->SetTitle("0,4<p<0.55GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0, 80<#Cluster<160, 35<dE/dx<60");
  ComputeRange(tree, "resolutionMIP", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
//  gr->GetHistogram()->SetMinimum(mipr_min);
//  gr->GetHistogram()->SetMaximum(mipr_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"resolutionMIP:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("resolutionMIP_vs_run.png");
  c1->Clear();

  /****** Mean energy loss for electrons ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanMIPele:run","",marA1,colPosA,1.0);
  gr->SetName("meanMIPele:run");
  gr->GetHistogram()->SetYTitle("Mean of electron dEdx");
  gr->GetHistogram()->SetTitle("0,32<p<0.38GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0, 80<#Cluster<160, 70<dE/dx<100");
  ComputeRange(tree, "meanMIPele", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
//  gr->GetHistogram()->SetMinimum(40);
//  gr->GetHistogram()->SetMaximum(110);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"meanMIPele:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meandEdxele_vs_run.png");
  c1->Clear();

  /****** Mean Energy loss electron Resolution ******/

  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"resolutionMIPele:run","",marA1,colPosA,1.0);
  gr->SetName("resolutionMIPele:run");
  gr->GetHistogram()->SetYTitle("Resolution of electrons dEdx");
  //gr->GetHistogram()->SetTitle("0,4<p<0.55GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0, 80<#Cluster<160, 35<dE/dx<60");
  gr->GetHistogram()->SetTitle("0,32<p<0.38GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 1.0, 80<#Cluster<160, 70<dE/dx<100");
  ComputeRange(tree, "resolutionMIPele", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
//  gr->GetHistogram()->SetMinimum(mipr_min);
//  gr->GetHistogram()->SetMaximum(mipr_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"resolutionMIPele:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("resolutionMeandEdxEle_vs_run.png");
  c1->Clear();

  /****** Separation Power ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"PIDSepPow_comb2:run","",marA1,colPosA,1.0);
  // TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"(meanMIPele-meanMIP)/(0.5*(resolutionMIP*meanMIP+resolutionMIPele*meanMIPele)):run","",marA1,colPosA,1.0);
  gr->SetName("(meanMIPele-meanMIP)/(0.5*(resolutionMIP*meanMIP+resolutionMIPele*meanMIPele)):run");
  gr->GetHistogram()->SetYTitle("Separation Power");
  gr->GetHistogram()->SetTitle("(MIP_ele-meanMIP)/(0.5*(sigma(MIP)+sigma(MIP_ele))");
  // ComputeRange(tree, "meanMIP", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(0);
  gr->GetHistogram()->SetMaximum(20);
//  gr->GetHistogram()->SetMinimum(mip_min);
//  gr->GetHistogram()->SetMaximum(mip_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"PIDSepPow_comb2:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("SeparationPower_vs_run.png");
  c1->Clear();

  ////////////////////////////////////////////////////////////////////////////////////////////////



  /****** Mean VertX ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanVertX:run","",marA1,colPosA,1.0);
  gr->SetName("meanVertX:run");
  gr->GetHistogram()->SetYTitle("Mean of Vert_{X} / [cm]");
  gr->GetHistogram()->SetTitle("");
//  ComputeRange(tree, "meanVertX", plotmean, plotoutlier);
//  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
//  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
  gr->GetHistogram()->SetMinimum(vx_min);
  gr->GetHistogram()->SetMaximum(vx_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"meanVertX:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanVertX_vs_run.png");
  c1->Clear();

  /****** Mean VertY ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanVertY:run","",marA1,colPosA,1.0);
  gr->SetName("meanVertY:run");
  gr->GetHistogram()->SetYTitle("Mean of Vert_{Y} / [cm]");
  gr->GetHistogram()->SetTitle("");
//  ComputeRange(tree, "meanVertY", plotmean, plotoutlier);
//  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
//  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
  gr->GetHistogram()->SetMinimum(vy_min);
  gr->GetHistogram()->SetMaximum(vy_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"meanVertY:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("meanVertY_vs_run.png");
  c1->Clear();

  /****** Mean VertZ ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanVertZ:run","",marA1,colPosA,1.0);
  gr->SetName("meanVertZ:run");
  gr->GetHistogram()->SetYTitle("Mean of Vert_{Z} / [cm]");
  gr->GetHistogram()->SetTitle("");
//  ComputeRange(tree, "meanVertZ", plotmean, plotoutlier);
//  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
//  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
  gr->GetHistogram()->SetMinimum(vz_min);
  gr->GetHistogram()->SetMaximum(vz_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"meanVertZ:run","");
  PlotTimestamp(entries,entries_tree,c1);
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
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"offsetd_comb4:run","",marSum,colSum,1.2);
  grComb->SetName("grComb");

  gr->GetHistogram()->SetYTitle("DCAs / [cm]");
  gr->GetHistogram()->SetTitle("p_{T} > 0.25GeV/c, |DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, |#eta| < 0.8");
  ComputeRange(tree, "offsetd_comb4", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(-plotmean-2*plotoutlier);
  gr->GetHistogram()->SetMaximum(+plotmean+2*plotoutlier);
//  gr->GetHistogram()->SetMinimum(dca_min);
//  gr->GetHistogram()->SetMaximum(dca_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  gr2->Draw("P");
  gr3->Draw("P");
  grComb->Draw("P");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("offsetdRA:run","DCA_{R}, A Side","p");
  leg->AddEntry("gr1","DCA_{Z}, A Side","p");
  leg->AddEntry("gr2","DCA_{R}, C Side","p");
  leg->AddEntry("gr3","DCA_{Z}, C Side","p");
  leg->AddEntry("grComb","combined = #sqrt{#Sigma x_{i}^{2}}","p");
  leg->Draw();

  PlotStatusLines(tree,"offsetd_comb4:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("DCAOffset_vs_run.png");
  c1->Clear();


  /****** Mean Mult  ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanMultPos:run:errorMultPos","",marA1,colPosA,1.0,sh_gr1);
  gr->SetName("meanMultPos:run");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanMultNeg:run:errorMultPos","",marA2,colNegA,1.0,sh_gr2);
  gr1->SetName("gr1");
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"meanMult_comb2:run","",marSum,colSum,1.2);
  grComb->SetName("grComb");

  gr->GetHistogram()->SetYTitle("Multiplicites of Primary Tracks");
  gr->GetHistogram()->SetTitle("|DCA_{R}| < 3cm, |DCA_{Z}| < 3cm, #Cluster > 70");
  ComputeRange(tree, "meanMult_comb2", plotmean, plotoutlier);
  // gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMinimum(0);
  gr->GetHistogram()->SetMaximum(plotmean+2*plotoutlier);
 // gr->GetHistogram()->SetMinimum(10);  //gr->GetHistogram()->SetMinimum(mult_min);
 // gr->GetHistogram()->SetMaximum(20);  //gr->GetHistogram()->SetMaximum(mult_max);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  grComb->Draw("P");

  TLegend *leg = new TLegend(0.6,0.80,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("meanMultPos:run","Positive Charged Tracks","p");
  leg->AddEntry("gr1","Negative Charged Tracks","p");
  leg->AddEntry("grComb","combined = (#Sigma x_{i})/N","p");
  leg->Draw();

  PlotStatusLines(tree,"meanMult_comb2:run","");
  PlotTimestamp(entries,entries_tree,c1);
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
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"tpcItsMatch_comb4:run","",marSum,colSum,1.2);
  grComb->SetName("grComb");

  gr->GetHistogram()->SetYTitle("Matching Efficiencies");
  gr->GetHistogram()->SetTitle("TPC-ITS Matching Efficiency");
  ComputeRange(tree, "tpcItsMatch_comb4", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMaximum(1.2);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  gr2->Draw("P");
  gr3->Draw("P");
  grComb->Draw("P");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("tpcItsMatchA:run","TPC-ITS matching, A Side","p");
  leg->AddEntry("gr1","TPC-ITS matching ( p_{T}>4GeV/c ), A Side","p");
  leg->AddEntry("gr2","TPC-ITS matching, C Side","p");
  leg->AddEntry("gr3","TPC-ITS matching ( p_{T}>4GeV/c ), C Side","p");
  leg->AddEntry("grComb","combined = (#Sigma x_{i})/N","p");
  leg->Draw();

  PlotStatusLines(tree,"tpcItsMatch_comb4:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("TPC-ITS-matching-efficiency_vs_run.png");
  c1->Clear();

  /****** ITS-TPC matching quality  ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"lambdaPull:run","",marA1,colPosA,1.0,sh_gr0);
  gr->SetName("lambdaPull:run");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ptPull:run","",marA2,colNegA,1.0,sh_gr1);
  gr1->SetName("gr1");
  TGraphErrors *gr2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"yPull:run","",marC1,colPosC,1.0,sh_gr2);
  gr2->SetName("gr2");
  TGraphErrors *gr3 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"zPull:run","",marC2,colNegC,1.0,sh_gr3);
  gr3->SetName("gr3");
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"itsTpcPulls_comb4:run","",marSum,colSum,1.2);
  grComb->SetName("grComb");

  gr->GetHistogram()->SetYTitle("Pulls");
  gr->GetHistogram()->SetTitle("ITS-TPC Matching Quality");
//  ComputeRange(tree, "itsTpcPulls_comb4", plotmean, plotoutlier);
//  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
//  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
  gr->GetHistogram()->SetMinimum(-3);
  gr->GetHistogram()->SetMaximum(3);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  gr2->Draw("P");
  gr3->Draw("P");
  grComb->Draw("P");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("lambdaPull:run","Tan#lambda pull bias","p");
  leg->AddEntry("gr1","q/p_{T} pull bias","p");
  leg->AddEntry("gr2","y pull bias","p");
  leg->AddEntry("gr3","z pull bias","p");
  leg->AddEntry("grComb","combined = #sqrt{#Sigma x_{i}^{2}}","p");
  leg->Draw();

  PlotStatusLines(tree,"itsTpcPulls_comb4:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c1, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c1->SaveAs("ITS-TPC-matching-quality_vs_run.png");
  c1->Clear();

  /****** pullPhi for TPC Constrain  ******/
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"tpcConstrainPhiA:run","",marA1,colPosA,1.0,sh_gr1);
  gr->SetName("tpcConstrainPhiA:run");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"tpcConstrainPhiC:run","",marC1,colPosC,1.0,sh_gr2);
  gr1->SetName("gr1");
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"tpcConstrainPhi_comb2:run","",marSum,colSum,1.2);
  grComb->SetName("grComb");

  gr->GetHistogram()->SetYTitle("(sin#phi_{TPC} - sin#phi_{Global})/#sigma");
  gr->GetHistogram()->SetTitle("");
  ComputeRange(tree, "tpcConstrainPhi_comb2", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(-plotmean-2*plotoutlier);
  gr->GetHistogram()->SetMaximum(+plotmean+2*plotoutlier);
//  gr->GetHistogram()->SetMinimum(-1);
//  gr->GetHistogram()->SetMaximum(1);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  grComb->Draw("P");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("tpcConstrainPhiA:run","A Side","p");
  leg->AddEntry("gr1","C Side","p");
  leg->AddEntry("grComb","combined = #sqrt{#Sigma x_{i}^{2}}","p");
  leg->Draw();

  PlotStatusLines(tree,"tpcConstrainPhi_comb2:run","");
  PlotTimestamp(entries,entries_tree,c1);
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
  ComputeRange(tree, "deltaPt", plotmean, plotoutlier);
  gr0->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr0->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
//  gr0->GetHistogram()->SetMinimum(-0.008);
//  gr0->GetHistogram()->SetMaximum(0.008);
  gr0->GetXaxis()->LabelsOption("v");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("deltaPt:run","deltaPt: both sides","ap");
  leg->AddEntry("deltaPtA:run","deltaPtA: A side only","p");
  leg->AddEntry("deltaPtC:run","deltaPtC: C side only","p");
  leg->Draw();

  PlotStatusLines(tree,"deltaPt:run","");
  PlotTimestamp(entries,entries_tree,c1);
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
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcarFitpar_comb4:run","",marSum,colSum,1.2);
  grComb->SetName("grComb");

  gr->GetHistogram()->SetYTitle("DCAR Fitting Parameters");
  gr->GetHistogram()->SetTitle("sqrt(P0^{2} + P1^{2}/(pT^{2}))");
  ComputeRange(tree, "dcarFitpar_comb4", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(-plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMaximum(+plotmean+3*plotoutlier);
  gr->GetHistogram()->SetMinimum(-1);
  gr->GetHistogram()->SetMaximum(1);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  gr2->Draw("P");
  gr3->Draw("P");
  grComb->Draw("P");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcarAP0:run","P0, A Side","p");
  leg->AddEntry("gr1","P1, A Side","p");
  leg->AddEntry("gr2","P0, C Side","p");
  leg->AddEntry("gr3","P1, C Side","p");
  leg->AddEntry("grComb","combined = (#Sigma x_{i})/N","p");
  leg->Draw();

  PlotStatusLines(tree,"dcarFitpar_comb4:run","");
  PlotTimestamp(entries,entries_tree,c1);
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
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar0_comb4:run","",1,1,1);
  DrawPlot(grComb, "grComb", marSum, 1.4, colSum, "P");

  gr0->GetHistogram()->SetYTitle("DCARs");
  ComputeRange(tree, "dcar0_comb4", plotmean, plotoutlier);
  gr0->GetHistogram()->SetMinimum(-plotmean-2*plotoutlier);
  gr0->GetHistogram()->SetMaximum(+plotmean+2*plotoutlier);
//  gr0->GetHistogram()->SetMinimum(-0.2);
//  gr0->GetHistogram()->SetMaximum(0.2);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcar_posA_0:run");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcar_posA_0:run","dcar_posA_0","p");
  leg->AddEntry("dcar_negA_0:run","dcar_negA_0","p");
  leg->AddEntry("dcar_posC_0:run","dcar_posC_0","p");
  leg->AddEntry("dcar_negC_0:run","dcar_negC_0","p");
  leg->AddEntry("grComb","combined = #sqrt{#Sigma x_{i}^{2}}","p");
  leg->Draw();

  PlotStatusLines(tree,"dcar0_comb4:run","");
  PlotTimestamp(entries,entries_tree,c1);
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
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar1_comb4:run","",1,1,1);
  DrawPlot(grComb, "grComb", marSum, 1.4, colSum, "P");

  gr0->GetHistogram()->SetYTitle("DCARs");
  ComputeRange(tree, "dcar1_comb4", plotmean, plotoutlier);
  gr0->GetHistogram()->SetMinimum(-plotmean-2*plotoutlier);
  gr0->GetHistogram()->SetMaximum(+plotmean+2*plotoutlier);
//  gr0->GetHistogram()->SetMinimum(-0.1);
//  gr0->GetHistogram()->SetMaximum(0.1);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcar_posA_1:run");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcar_posA_1:run","dcar_posA_1","p");
  leg->AddEntry("dcar_negA_1:run","dcar_negA_1","p");
  leg->AddEntry("dcar_posC_1:run","dcar_posC_1","p");
  leg->AddEntry("dcar_negC_1:run","dcar_negC_1","p");
  leg->AddEntry("grComb","combined = #sqrt{#Sigma x_{i}^{2}}","p");
  leg->Draw();

  PlotStatusLines(tree,"dcar1_comb4:run","");
  PlotTimestamp(entries,entries_tree,c1);
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
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcar2_comb4:run","",1,1,1);
  DrawPlot(grComb, "grComb", marSum, 1.4, colSum, "P");

  gr0->GetHistogram()->SetYTitle("DCARs");
  ComputeRange(tree, "dcar2_comb4", plotmean, plotoutlier);
  gr0->GetHistogram()->SetMinimum(-plotmean-2*plotoutlier);
  gr0->GetHistogram()->SetMaximum(+plotmean+2*plotoutlier);
//  gr0->GetHistogram()->SetMinimum(-0.1);
//  gr0->GetHistogram()->SetMaximum(0.1);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcar_posA_2:run");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcar_posA_2:run","dcar_posA_2","p");
  leg->AddEntry("dcar_negA_2:run","dcar_negA_2","p");
  leg->AddEntry("dcar_posC_2:run","dcar_posC_2","p");
  leg->AddEntry("dcar_negC_2:run","dcar_negC_2","p");
  leg->AddEntry("grComb","combined = #sqrt{#Sigma x_{i}^{2}}","p");
  leg->Draw();

  PlotStatusLines(tree,"dcar2_comb4:run","");
  PlotTimestamp(entries,entries_tree,c1);
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
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz0_comb4:run","",1,1,1);
  DrawPlot(grComb, "grComb", marSum, 1.4, colSum, "P");

  gr0->GetHistogram()->SetYTitle("DCAZs");
  ComputeRange(tree, "dcaz0_comb4", plotmean, plotoutlier);
  gr0->GetHistogram()->SetMinimum(-plotmean-2*plotoutlier);
  gr0->GetHistogram()->SetMaximum(+plotmean+2*plotoutlier);
//  gr0->GetHistogram()->SetMinimum(-2.);
//  gr0->GetHistogram()->SetMaximum(2.);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcaz_posA_0:run");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcaz_posA_0:run","dcaz_posA_0","p");
  leg->AddEntry("dcaz_negA_0:run","dcaz_negA_0","p");
  leg->AddEntry("dcaz_posC_0:run","dcaz_posC_0","p");
  leg->AddEntry("dcaz_negC_0:run","dcaz_negC_0","p");
  leg->AddEntry("grComb","combined = #sqrt{#Sigma x_{i}^{2}}","p");
  leg->Draw();

  PlotStatusLines(tree,"dcaz0_comb4:run","");
  PlotTimestamp(entries,entries_tree,c1);
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
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz1_comb4:run","",1,1,1);
  DrawPlot(grComb, "grComb", marSum, 1.4, colSum, "P");

  gr0->GetHistogram()->SetYTitle("DCAZs");
  ComputeRange(tree, "dcaz1_comb4", plotmean, plotoutlier);
  gr0->GetHistogram()->SetMinimum(-plotmean-2*plotoutlier);
  gr0->GetHistogram()->SetMaximum(+plotmean+2*plotoutlier);
//  gr0->GetHistogram()->SetMinimum(-0.2);
//  gr0->GetHistogram()->SetMaximum(0.2);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcaz_posA_1:run");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcaz_posA_1:run","dcaz_posA_1","p");
  leg->AddEntry("dcaz_negA_1:run","dcaz_negA_1","p");
  leg->AddEntry("dcaz_posC_1:run","dcaz_posC_1","p");
  leg->AddEntry("dcaz_negC_1:run","dcaz_negC_1","p");
  leg->AddEntry("grComb","combined = #sqrt{#Sigma x_{i}^{2}}","p");
  leg->Draw();

  PlotStatusLines(tree,"dcaz1_comb4:run","");
  PlotTimestamp(entries,entries_tree,c1);
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
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"dcaz2_comb4:run","",1,1,1);
  DrawPlot(grComb, "grComb", marSum, 1.4, colSum, "P");

  gr0->GetHistogram()->SetYTitle("DCAZs");
  ComputeRange(tree, "dcaz2_comb4", plotmean, plotoutlier);
  gr0->GetHistogram()->SetMinimum(-plotmean-2*plotoutlier);
  gr0->GetHistogram()->SetMaximum(+plotmean+2*plotoutlier);
//  gr0->GetHistogram()->SetMinimum(-0.1);
//  gr0->GetHistogram()->SetMaximum(0.1);
  gr0->GetHistogram()->SetTitleOffset(10);
  gr0->GetXaxis()->LabelsOption("v");
  gr0->SetName("dcaz_posA_2:run");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("dcaz_posA_2:run","dcaz_posA_2","p");
  leg->AddEntry("dcaz_negA_2:run","dcaz_negA_2","p");
  leg->AddEntry("dcaz_posC_2:run","dcaz_posC_2","p");
  leg->AddEntry("dcaz_negC_2:run","dcaz_negC_2","p");
  leg->AddEntry("grComb","combined = #sqrt{#Sigma x_{i}^{2}}","p");
  leg->Draw();

  PlotStatusLines(tree,"dcaz2_comb4:run","");
  PlotTimestamp(entries,entries_tree,c1);
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

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("occ_AC_Side_IROC_OROC:run","iroc_A_side","p");
  leg->AddEntry("gr1","oroc_A_side","p");
  leg->AddEntry("gr2","iroc_C_side","p");
  leg->AddEntry("gr3","oroc_C_side","p");
  leg->Draw();

  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c9, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c9->SaveAs("occ_AC_Side_IROC_OROC_vs_run.png");//for C,A side and IROC,OROC
  c9->Update();

  /****** attachment parameters for A and C side ******/

  TCanvas *c10  = new  TCanvas("can10","can10",canvas_width,canvas_height);
  c10->cd();
  c10->Update();
  c10->SetGrid(3);

  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"MIPattachSlopeA:run","",marA1,colPosA,1.0,sh_gr1);
  gr->SetName("MIPattachSlopeA:run");
  TGraphErrors *gr1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"MIPattachSlopeC*(-1):run","",marC1,colPosC,1.0,sh_gr2);
  gr1->SetName("gr1");
  TGraphErrors *grComb = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"MIPattachSlope_comb2:run","",marSum,colSum,1.2);
  grComb->SetName("grComb");

  gr->GetHistogram()->SetYTitle("Attachment parameter p1");
  gr->GetHistogram()->SetTitle("showing p1 of fit: p0 + p1 * tan(#theta)"); // info from Marian, 19.11.2014. to be checked in code that produces the tree.
  ComputeRange(tree, "MIPattachSlope_comb2", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");
  gr1->Draw("P");
  grComb->Draw("P");

  TLegend *leg = new TLegend(0.6,0.75,0.62*sqrt(norm_runs/entries),0.95,"","brNDC");
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry("MIPattachSlopeA:run","A Side","p");
  leg->AddEntry("gr1","C Side *(-1)","p");
  leg->AddEntry("grComb","combined = (#Sigma x_{i})/N","p");
  leg->Draw();

  PlotStatusLines(tree,"MIPattachSlope_comb2:run","");
  PlotTimestamp(entries,entries_tree,c1);
  TStatToolkit::AddStatusPad(c10, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c10->SaveAs("MIPattachSlopes_vs_run.png");
  c10->Clear();

//  //C side
//  TCanvas *c11  = new  TCanvas("can11","can11",canvas_width,canvas_height);
//  c11->cd();
//  c11->Update();
//  c11->SetGrid(3);
//
//  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"MIPattachSlopeC*(-1):run","",marA1,colPosA,1.0);
//  gr->SetName("MIPattachSlopeC:run");
//  gr->GetHistogram()->SetYTitle("Attachment parameter p1");
//  gr->GetHistogram()->SetMinimum(-10);
//  gr->GetHistogram()->SetMaximum(+10);
//  gr->GetXaxis()->LabelsOption("v");
//  gr->Draw("AP");
//
//  PlotStatusLines(tree,"MIPattachSlopeC:run","");
//  PlotTimestamp(entries,entries_tree,c1);
//  TStatToolkit::AddStatusPad(c11, statPadHeight, statPadBotMar);
//  TStatToolkit::DrawStatusGraphs(oaMultGr);
//  c11->SaveAs("MIPattachSlopeC_vs_run.png");
//  c11->Clear();

  /****** electron and MIPs separation ******/

  TCanvas *c12  = new  TCanvas("can12","can12",canvas_width,canvas_height);
  c12->cd();
  c12->Update();
  c12->SetGrid(3);

  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"electroMIPSeparation:run","",marA1,colPosA,1.0);
  gr->SetName("electroMIPSeparation:run");
  gr->GetHistogram()->SetYTitle("Electron - MIP");
  ComputeRange(tree, "electroMIPSeparation", plotmean, plotoutlier);
  gr->GetHistogram()->SetMinimum(plotmean-3*plotoutlier);
  gr->GetHistogram()->SetMaximum(plotmean+3*plotoutlier);
  gr->GetXaxis()->LabelsOption("v");
  gr->Draw("AP");

  PlotStatusLines(tree,"electroMIPSeparation:run","");
  PlotTimestamp(entries,entries_tree,c12);
  TStatToolkit::AddStatusPad(c12, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  c12->SaveAs("ElectroMIPSeparation_vs_run.png");
  c12->Clear();

  //
  // save status into Tree and write to rootfile
  //
  if (statusTree) {
    cout << "updating trending rootfile with status tree... ";
    //tree->AddFriend(statusTree,"Tstatus");
    tree->Write("", TObject::kOverwrite);
    statusTree->Write();
    cout << " successful." << endl;
  }

  cout << "...done with trending." << endl;
  return 1;
}


Int_t PlotStatusLines(TTree * tree, const char * expr, const char * cut)
{
  //the function plots status lines
  char* alias = "varname_OutlierMin:varname_OutlierMax:varname_WarningMin:varname_WarningMax:varname_PhysAccMin:varname_PhysAccMax:varname_RobustMean";
  TMultiGraph* mgStatusLines = TStatToolkit::MakeStatusLines(tree,expr,cut,alias);

  if (mgStatusLines) mgStatusLines->Draw("l");
  else { cout << " no mgStatusLines available!" << endl; return 0; }

  return 1;
}

Int_t PlotTimestamp(const int nruns=0, const int nentries=0, TCanvas* c1)
{
  double rightlegend = 1 - 10./c1->GetWw();
  //the function plots a timestamp, the used Aliroot version, and the number of runs
  TString sTimestamp  = gSystem->GetFromPipe("date");
  TString sAlirootVer;
  if (gSystem->GetFromPipe("echo $ALICE_VER") == "master"){
    sAlirootVer = "AliRoot: " + gSystem->GetFromPipe("wdir=`pwd`; cd $ALICE_ROOT/../src; git describe; cd $wdir;");
  }
  else {
    sAlirootVer = "AliRoot: " + gSystem->GetFromPipe("echo $ALIROOT_VERSION");
  }

  TString sAliphysicsVer;
  if (gSystem->GetFromPipe("echo $ALIPHYSICS_VER") == "master"){
    sAliphysicsVer = "AliPhysics: " + gSystem->GetFromPipe("wdir=`pwd`; cd $ALICE_PHYSICS/../src; git describe; cd $wdir;");
  }
  else {
    sAliphysicsVer = "AliPhysics: " + gSystem->GetFromPipe("echo $ALIPHYSICS_VERSION");
  }
  TLatex* latTime = new TLatex(rightlegend,0.13,sTimestamp.Data());
  latTime->SetTextSize(0.03);
  latTime->SetTextAlign(31);
  latTime->SetNDC();
  latTime->Draw("same");
  TLatex* latAliroot = new TLatex(rightlegend,0.09,sAlirootVer.Data());
  latAliroot->SetTextSize(0.03);
  latAliroot->SetTextAlign(31);
  latAliroot->SetNDC();
  latAliroot->Draw("same");
  TLatex* latAliphysics = new TLatex(rightlegend,0.05,sAliphysicsVer.Data());
  latAliphysics->SetTextSize(0.03);
  latAliphysics->SetTextAlign(31);
  latAliphysics->SetNDC();
  latAliphysics->Draw("same");
  TLatex* latNruns = new TLatex(rightlegend,0.01,Form("N shown runs: %i (tree entries: %i)",nruns,nentries));
  latNruns->SetTextSize(0.03);
  latNruns->SetTextAlign(31);
  latNruns->SetNDC();
  if (nruns>0) latNruns->Draw("same");
return 1;
}

Int_t ComputeRange(TTree* tree, const char* varname, Float_t &plotmean, Float_t &plotoutlier)
{
  //the function computes useful numbers for plot ranges from the outlier criteria
  plotmean    = (Float_t) TFormula("fcn", tree->GetAlias(Form("%s_RobustMean",varname))).Eval(0);
  plotoutlier = (Float_t) TFormula("fcn", tree->GetAlias(Form("%s_OutlierMax",varname))).Eval(0) - plotmean;
  return 1;
}

Int_t DrawPlot(TGraphErrors* gr, TString nameHisto, Int_t markerStyle, Int_t markerSize, Int_t markerColor, TString drawMode)
{
  //the function draws the plots
  gr->SetName(nameHisto);
  gr->SetMarkerStyle(markerStyle);
  gr->SetMarkerSize(markerSize);
  gr->SetMarkerColor(markerColor);
  gr->Draw(drawMode);
  return 1;
}
