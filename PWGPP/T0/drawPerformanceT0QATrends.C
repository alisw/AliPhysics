
TGraphErrors * MakeGraphSparse(TTree * tree, const Char_t * expr="mean-fdelta:run", const Char_t * cut="isTPC&&ptype==0&&theta>0"){
  // Format of expr: Var:Error:Run
  //
  // Make a sparse draw of the variables
  //  Format of expr : Var:Error:Run
  //

  const Int_t entries =  tree->Draw(expr,cut,"goff");
  if (entries<=0) return 0x0;
  Double_t *graphX, *graphY, *graphError;

  //check whether error argument exists in expr
  if(!tree->GetV3()){
      graphX = tree->GetV2();
      graphY = tree->GetV1();
  }
  else{
      graphX = tree->GetV3();
      graphY = tree->GetV1();
      graphError = tree->GetV2();
  }

  // sort according to run number
  Int_t *index = new Int_t[entries];
  TMath::Sort(entries,graphX,index,false);

  // define arrays for the new graph
  Double_t *tempArray = new Double_t[entries];
  Double_t *xError = new Double_t[entries];
  Double_t *yError = new Double_t[entries];
  Int_t *vrun = new Int_t[entries]; 
  Double_t count = 0.5;

  // evaluate arrays for the new graph accroding to the run-number
  Int_t icount=0;
  tempArray[index[0]] = count;
  xError[0] = 0;
  yError[0] = 0;
  if(tree->GetV3())
      yError[index[0]] = graphError[index[0]];
  vrun[0] = graphX[index[0]];

  // loop the rest entries
  for(Int_t i=1;i<entries;i++){
      xError[i] = 0;
      yError[i] = 0;

      if(tree->GetV3())
	  yError[i] = graphError[index[i]];


      if(graphX[index[i]]==graphX[index[i-1]])
	  tempArray[index[i]] = count; 
      else if((graphX[index[i]]!=graphX[index[i-1]])){
	      count++;
	      icount++;
	      tempArray[index[i]] = count;
	      vrun[icount]=graphX[index[i]];
      }
  }

  // count the number of xbins (run-wise) for the new graph
  const Int_t newNbins = int(count+0.5);
  Double_t *newBins = new Double_t[newNbins+1];
  for(Int_t i=0; i<=count+1;i++){
    newBins[i] = i;
  }
  
  // define and fill the new graph
  TGraphErrors *graphNew = new TGraphErrors(entries,tempArray,graphY,xError,yError);
  graphNew->GetXaxis()->Set(newNbins,newBins);
  
  // set the bins for the x-axis
  Char_t xName[50];
  for(Int_t i=0;i<count;i++){
    sprintf(xName,"%d",Int_t(vrun[i]));
    graphNew->GetXaxis()->SetBinLabel(i+1,xName);
  }
  graphNew->GetHistogram()->SetTitle("");
  
  // memory clearing
  delete [] xError; 
  delete [] yError; 
  delete [] tempArray;
  delete [] index;
  delete [] newBins;
  delete [] vrun;
  return graphNew;

}
//------------------------------------------------------------------------------------------------
void drawPerformanceT0QATrends(const char* inFile = "trending.root", const char* runType="pp") {
  //
  //
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetLabelSize(0.04,"x");

 
  // open input file
  //
  TFile *file = TFile::Open(inFile);
  if(!file) return;
  file->cd();
  
  TTree *tree = (TTree*)file->Get("trending");
  if(!tree) return;
  int const entries = tree->GetEntries();
  cout<<"number of entries   "<<entries<<endl; 
  TH1F *frame = new TH1F();

  float norm_runs = 10.0; 
  if(entries<8)
    norm_runs =10.0;
  else if(entries>=8&&entries<16)
    norm_runs =20.0;
  else
    norm_runs =50.0;
  int const canvas_width = int(((entries*1.0)/norm_runs)*2000.0);
  if(entries>50){
    gStyle->SetTickLength(0.03*norm_runs/(entries*1.0),"Y");
    gStyle->SetTitleYOffset((norm_runs/(entries*1.0))*0.8);
    gStyle->SetPadLeftMargin(0.1*norm_runs/(entries*1.0));
    gStyle->SetPadRightMargin(0.1*norm_runs/(entries*1.0));
  }

  //Define ranges of you trending plots
  double resolutionMin = 10,   resolutionMax = 100;    // OR A - OR C
  double oraplusorcMin = -100, oraplusorcMax = 100;   // OR A + OR CA
  double oraMin        = -100, oraMax        = 100;   // OR A  
  double orcMin        = -100, orcMax        = 100;   // OR C
  //  double amplMin       =0,    amplMax    =3   ;     // amplitude in each PMT
  //  double timeMin     =  2900,   timeMax =3300       ;     // amplitude in each PMT
  //-----> add ranges of your new trending plot


  TCanvas *c1  = new  TCanvas("can","can",canvas_width,500);
  c1->SetGridy(1);
  c1->SetGridx(1);
  c1->SetBottomMargin(0.17);
  
  /****** T0 ORA+ORC ******/
  double orMinCalib    = -20,  orMaxCalib    = 20;    // OR limits for GOOD run with calibration
  double orMinNoCalib  = -40,  orMaxNoCalib  = 40;    // OR limits for GOOD run without calibration
  
  TGraphErrors *grSum = MakeGraphSparse(tree,"tzeroOrAPlusOrC:run","");
  grSum->SetMarkerStyle(20);
  grSum->SetMarkerSize(1.0);
  grSum->SetMarkerColor(2);
  TGraphErrors *grORA = MakeGraphSparse(tree,"tzeroOrA:run","");
  grORA->SetMarkerStyle(28);
  grORA->SetMarkerSize(1.0);
  grORA->SetMarkerColor(4);
  TGraphErrors *grORC = MakeGraphSparse(tree,"tzeroOrC:run","");
  grORC->SetMarkerStyle(25);
  grORC->SetMarkerSize(1.0);
  grORC->SetMarkerColor(1);
 
  grSum->GetHistogram()->SetYTitle("mean [ps]");
  grSum->GetHistogram()->SetTitle("T0 ORA, ORC and (ORA+ORC)/2");
  grSum->GetHistogram()->SetMinimum(oraplusorcMin);
  grSum->GetHistogram()->SetMaximum(oraplusorcMax);
  
  
  int nRuns = grSum->GetN();
  double *x =  grSum->GetX();
  double min = x[0];
  double max = x[0];
  for(int irun =1; irun<nRuns;irun++){
      if(min > x[irun] && x[irun]>0) min = x[irun];
      if(max < x[irun]) max = x[irun];
  }
  min-=0.5; max+=0.5;
  TBox* outOfLimits = new TBox(min,oraplusorcMin,max,oraplusorcMax);
  outOfLimits->SetFillColor(kOrange);
  TBox* limitsNoCalib = new TBox(min,orMinNoCalib,max,orMaxNoCalib);
  limitsNoCalib->SetFillColor(kYellow);
  TBox* limitsCalib = new TBox(min,orMinCalib,max,orMaxCalib);
  limitsCalib->SetFillColor(kTeal);
  TLine* limitsOR[4];
  limitsOR[0] = new TLine(min,orMinCalib,max,orMinCalib);    limitsOR[0]->SetLineColor(kGreen);limitsOR[0]->SetLineWidth(3);
  limitsOR[1] = new TLine(min,orMaxCalib,max,orMaxCalib);    limitsOR[1]->SetLineColor(kGreen);limitsOR[1]->SetLineWidth(3);
  limitsOR[2] = new TLine(min,orMinNoCalib,max,orMinNoCalib);limitsOR[2]->SetLineColor(kRed);  limitsOR[2]->SetLineWidth(3);
  limitsOR[3] = new TLine(min,orMaxNoCalib,max,orMaxNoCalib);limitsOR[3]->SetLineColor(kRed);  limitsOR[3]->SetLineWidth(3);
  
  grSum->Draw("AP");
  outOfLimits->Draw("same");
  limitsNoCalib->Draw("same");
  limitsCalib->Draw("same");
//   for(Int_t i=0;i<nRuns-1;i++)
//   {
//       ((TLine*)(TLine::DrawLine(min,oraplusorcMin,min,oraplusorcMax)));//x[i]+0.5,oraplusorcMin,x[i]+0.5,oraplusorcMax)));//->SetLineStyle(2);
//   }
  grSum->Draw("psame");
  grORA->Draw("psame");
  grORC->Draw("psame");
  for(Int_t i=0;i<4;i++)limitsOR[i]->Draw("same");
  TLegend *leg = new TLegend(0.1,0.85,0.3,0.95," ","brNDC");
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);leg->SetNColumns(3);leg->SetColumnSeparation(1);
  leg->AddEntry(grORA,"ORA","p"); 
  leg->AddEntry(grORC,"ORC","p"); 
  leg->AddEntry(grSum,"(ORA+ORC)/2","p"); 
  leg->Draw(); 
 
  grSum->GetXaxis()->LabelsOption("v");
  c1->SaveAs("meanT0OrAPlusOrC_vs_run.png");
 
  /****** T0 Resolution ******/
  TGraphErrors *gr = MakeGraphSparse(tree,"resolution:run","");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.0);
  gr->SetMarkerColor(2);
  gr->GetHistogram()->SetMinimum(resolutionMin);
  gr->GetHistogram()->SetMaximum(resolutionMax);
  gr->Draw("AP");
  gr->GetXaxis()->LabelsOption("v");
  gr->GetHistogram()->SetYTitle("sigma [ps]");
  gr->GetHistogram()->SetTitle("T0 resolution (ORA -ORC)/2");
  c1->SaveAs("sigmaResolutionT0_vs_run.png");
 
  /****** Mean Amplitude in PMT ******/
  const int kNPMTs = 24;
  char name[200];

  for(int ipmt=1;ipmt<=kNPMTs; ipmt++){
    sprintf(name,"amplPMT%d:run",ipmt);
    TString cutamp = Form("amplPMT%i>0",ipmt);
    TGraphErrors *gramp = MakeGraphSparse(tree, name, cutamp.Data());
    if (!gramp)
      continue;
    // gr = MakeGraphSparse(tree,name,"");
    gramp->SetMarkerStyle(20);
    gramp->SetMarkerSize(1.0);
    gramp->SetMarkerColor(6);
    gramp->GetHistogram()->SetYTitle("mean");
    gramp->GetHistogram()->SetTitle(Form("Amplitude PMT%d", ipmt));

    int nRuns = gramp->GetN();
    double *y = gramp->GetY();
    double min = y[0];
    double max = y[0];
    for (int irun = 1; irun < nRuns; irun++) {
      if (min > y[irun] & y[irun] > 0)
        min = y[irun];
      if (max < y[irun])
        max = y[irun];
    }
    //  amplMin = min - 2;
    //  amplMax = max + 2;

    //    gr->GetHistogram()->SetMinimum(amplMin);
    //   gr->GetHistogram()->SetMaximum(amplMax);
    gramp->Draw("AP");
    gramp->GetXaxis()->LabelsOption("v");
    c1->SaveAs(Form("meanAmplPMT%d_vs_run.png", ipmt));
    }
    /****** Mean Time in PMT ******/
    for (int ipmt = 1; ipmt <= kNPMTs; ipmt++) {
      sprintf(name, "timePMT%d:run", ipmt);
      TString cut = Form("timePMT%i>0", ipmt);
      TGraphErrors *grtime = MakeGraphSparse(tree, name, cut.Data());
      if (!grtime)
        continue;
      // regular run
      int nRuns = grtime->GetN();
      double *y = grtime->GetY();
      double min = y[0];
      double max = y[0];
      for (int irun = 1; irun < nRuns; irun++) {
        if (min > y[irun] && y[irun] > 0)
          min = y[irun];
        if (max < y[irun])
          max = y[irun];
      }
      grtime->SetMarkerStyle(20);
      grtime->SetMarkerSize(1.0);
      grtime->SetMarkerColor(2);

      grtime->GetHistogram()->SetYTitle("mean [channels]");
      grtime->GetHistogram()->SetTitle(Form("Time PMT%d", ipmt));
      //    gr->GetHistogram()->SetMinimum(timeMin);
      // gr->GetHistogram()->SetMaximum(timeMax);
      grtime->GetXaxis()->LabelsOption("v");
      grtime->Draw("AP");

      TLegend *leg = new TLegend(0.1, 0.85, 0.3, 0.95, " ", "brNDC");
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.05);
      leg->SetNColumns(3);
      leg->AddEntry(grtime, "mean time", "p");
      leg->Draw();

      c1->SaveAs(Form("meanTimePMT%d_vs_run.png", ipmt));
    }
    //   TGraphErrors *grEfficiencySPD =
    //   MakeGraphSparse(tree,"efficiencySPD:run","");
    //   grEfficiencySPD->SetMarkerStyle(25);
    //   grEfficiencySPD->SetMarkerSize(1.0);
    //   grEfficiencySPD->SetMarkerColor(1);
    //   grEfficiencySPD->GetXaxis()->LabelsOption("v");
    //   grEfficiencySPD->SetTitle("T0 to SPD
    //   ratio;run;Nevts_{T0}/Nevts_{SPD}"); grEfficiencySPD->Draw("AP");
    // c1->SaveAs("efficiencyT0toSPD_vs_run.gif");
    if (tree->GetBranch("efficiency0TVX_CINT7")) {
      TGraphErrors *grEfficiency0TVX_CINT7 =
          MakeGraphSparse(tree, "efficiency0TVX_CINT7:run", "");
      if (grEfficiency0TVX_CINT7) {
        grEfficiency0TVX_CINT7->GetYaxis()->SetRangeUser(0., 1.);
        grEfficiency0TVX_CINT7->SetMarkerStyle(25);
        grEfficiency0TVX_CINT7->SetMarkerSize(1.0);
        grEfficiency0TVX_CINT7->SetMarkerColor(1);
        grEfficiency0TVX_CINT7->GetXaxis()->LabelsOption("v");
        grEfficiency0TVX_CINT7->SetTitle("T0 to V0 ratio;run;0TVX/CINT7");
        grEfficiency0TVX_CINT7->Draw("AP");
        c1->SaveAs("efficiencyT0toV0_vs_run.png");
      }
  }
  if ( tree->GetBranch("efficiency0TVX_CADAND")) {
    TGraphErrors *grEfficiency0TVX_CADAND =
        MakeGraphSparse(tree, "efficiency0TVX_CADAND:run", "");
    if (grEfficiency0TVX_CADAND) {
      grEfficiency0TVX_CADAND->GetYaxis()->SetRangeUser(0., 1.);
      grEfficiency0TVX_CADAND->SetMarkerStyle(22);
      grEfficiency0TVX_CADAND->SetMarkerSize(1.0);
      grEfficiency0TVX_CADAND->SetMarkerColor(2);
      grEfficiency0TVX_CADAND->GetXaxis()->LabelsOption("v");
      grEfficiency0TVX_CADAND->SetTitle("T0 to AD ratio;run;0TVX/CADAND");
      grEfficiency0TVX_CADAND->Draw("AP");
      c1->SaveAs("efficiencyT0toAD_vs_run.png");
    }
  }
  // leg->Clear();
  //   TLegend *leg = new TLegend(0.1,0.8,0.3,0.95," ","efficiency");
  //   leg->SetFillStyle(0); leg->SetBorderSize(0);
  //   leg->SetTextSize(0.05);leg->SetNColumns(1);//leg->SetColumnSeparation(1);
  //   leg->AddEntry(grEfficiency0TVX_CINT7,"T0 to V0 ratio","p");
  //   leg->AddEntry(grEfficiency0TVX_CADAND,"T0 to AD ratio","p");
  // leg->AddEntry(grSum,"(ORA+ORC)/2","p");
  // leg->Draw();
  // c1->SaveAs("efficiency_vs_run.gif");

  //-----> draw your new trending plot here
  c1->Close();
}

