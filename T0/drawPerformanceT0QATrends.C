
TGraphErrors * MakeGraphSparse(TTree * tree, const Char_t * expr="mean-fdelta:run", const Char_t * cut="isTPC&&ptype==0&&theta>0"){
  // Format of expr: Var:Error:Run
  //
  // Make a sparse draw of the variables
  //  Format of expr : Var:Error:Run
  //

  const Int_t entries =  tree->Draw(expr,cut,"goff");
  
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
drawPerformanceT0QATrends(const char* inFile = "trending.root", const char* runType="pp") {
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
  
  TTree *tree = (TTree*)file->Get("t0QA");
  if(!tree) return;
  int const entries = tree->GetEntries();
  cout<<"number of entries   "<<entries<<endl; 
  TH1F *frame = new TH1F();

  //const float norm_runs = 8.0; 
  if(entries<8)
    const float norm_runs =10.0;
  else if(entries>=8&&entries<16)
    const float norm_runs =20.0;
 else
   const float norm_runs =50.0;
  int const canvas_width = int(((entries*1.0)/norm_runs)*2000.0);
  if(entries>50){
    gStyle->SetTickLength(0.03*norm_runs/(entries*1.0),"Y");
    gStyle->SetTitleYOffset((norm_runs/(entries*1.0))*0.8);
    gStyle->SetPadLeftMargin(0.1*norm_runs/(entries*1.0));
    gStyle->SetPadRightMargin(0.1*norm_runs/(entries*1.0));
  }

  //Define ranges of you trending plots
  double resolutionMin = 10,   resolutionMax = 50;    // OR A - OR C
  double oraplusorcMin = -100, oraplusorcMax = 100;   // OR A + OR CA
  double oraMin        = -100, oraMax        = 100;   // OR A  
  double orcMin        = -100, orcMax        = 100;   // OR C
  double amplMin       =0,    amplMax    =3   ;     // amplitude in each PMT
  double timeMin       ,   timeMax       ;     // amplitude in each PMT
  //-----> add ranges of your new trending plot


  TCanvas *c1  = new  TCanvas("can","can",canvas_width,500);
  c1->SetGridy(1);
  c1->SetGridx(1);
  c1->SetBottomMargin(0.17);
  
  /****** T0 ORA+ORC ******/
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
  grSum->Draw("AP");
  grORA->Draw("psame");
  grORC->Draw("psame");
  TLegend *leg = new TLegend(0.1,0.85,0.3,0.95," ","brNDC");
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);leg->SetNColumns(3);leg->SetColumnSeparation(1);
  leg->AddEntry(grORA,"ORA","p"); 
  leg->AddEntry(grORC,"ORC","p"); 
  leg->AddEntry(grSum,"(ORA+ORC)/2","p"); 
  leg->Draw(); 
 
  grSum->GetXaxis()->LabelsOption("v");
  c1->SaveAs("meanT0OrAPlusOrC_vs_run.gif");
 
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
  c1->SaveAs("sigmaResolutionT0_vs_run.gif");
 

   /****** Mean T0 OR A ******/
 /* TGraphErrors *gr = MakeGraphSparse(tree,"tzeroOrA:run","");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.0);
  gr->SetMarkerColor(2);
  gr->GetHistogram()->SetYTitle("mean [ps]");
  gr->GetHistogram()->SetTitle("T0 OR A");
  gr->GetHistogram()->SetMinimum(oraMin);
  gr->GetHistogram()->SetMaximum(oraMax);
  gr->Draw("AP");
  gr->GetXaxis()->LabelsOption("v");
  c1->SaveAs("meanT0OrA_vs_run.gif");
*/

  /****** Mean T0 OR C ******/
 /* TGraphErrors *gr = MakeGraphSparse(tree,"tzeroOrC:run","");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.0);
  gr->SetMarkerColor(2);
  gr->GetHistogram()->SetYTitle("mean [ps]");
  gr->GetHistogram()->SetTitle("T0 OR C");
  gr->GetHistogram()->SetMinimum(orcMin);
  gr->GetHistogram()->SetMaximum(orcMax);
  gr->Draw("AP");
  gr->GetXaxis()->LabelsOption("v");
  c1->SaveAs("meanT0OrC_vs_run.gif");
*/

  /****** Mean Amplitude in PMT ******/
  const int kNPMTs = 24;
  char name[200];

  for(int ipmt=1;ipmt<=kNPMTs; ipmt++){
    sprintf(name,"amplPMT%d:run",ipmt);
    gr = MakeGraphSparse(tree,name,"");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);
    gr->SetMarkerColor(6);
    gr->GetHistogram()->SetYTitle("mean");
    gr->GetHistogram()->SetTitle(Form("Amplitude PMT%d",ipmt));

    int nRuns = gr->GetN();
    double *y =  gr->GetY();
	double min = y[0];
	double max = y[0];
    for(int irun =1; irun<nRuns;irun++){
      if(min > y[irun] & y[irun]>0) min = y[irun];
      if(max < y[irun]) max = y[irun];
    }
    //  amplMin = min - 2; 
    //  amplMax = max + 2; 

    gr->GetHistogram()->SetMinimum(amplMin);
    gr->GetHistogram()->SetMaximum(amplMax);
    gr->Draw("AP");
    gr->GetXaxis()->LabelsOption("v");
    c1->SaveAs(Form("meanAmplPMT%d_vs_run.gif",ipmt));
  }
  /****** Mean Time in PMT ******/
  for(int ipmt=1;ipmt<=kNPMTs; ipmt++){
    sprintf(name,"timePMT%d:run",ipmt);
    gr = MakeGraphSparse(tree,name,"");


    sprintf(name,"timeDelayPMT%d:run",ipmt);
    TGraphErrors *grDelay = MakeGraphSparse(tree,name,"");
    //regular run
    int nRuns = gr->GetN();
    double *y =  gr->GetY();
	double min = y[0];
	double max = y[0];
    for(int irun =1; irun<nRuns;irun++){
      if(min > y[irun] && y[irun]>0) min = y[irun];
      if(max < y[irun]) max = y[irun];
    }
    //Delay
    //   double *yDelay =  grDelay->GetY();
    //  nRuns = grDelay->GetN();
    //   for(int irun =0; irun<nRuns;irun++){
    //      if(min > yDelay[irun] && yDelay[irun]>0) min = yDelay[irun];
    //    if(max < yDelay[irun]) max = yDelay[irun];
    //  }
    
    timeMin = min - 2; 
    timeMax = max + 2; 

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);
    gr->SetMarkerColor(2);
    grDelay->SetMarkerStyle(24);
    grDelay->SetMarkerSize(1.0);
    grDelay->SetMarkerColor(1);

    gr->GetHistogram()->SetYTitle("mean [channels]");
    gr->GetHistogram()->SetTitle(Form("Time PMT%d",ipmt));
    gr->GetHistogram()->SetMinimum(timeMin);
    gr->GetHistogram()->SetMaximum(timeMax);
    gr->GetXaxis()->LabelsOption("v");
    gr->Draw("AP");
    grDelay->Draw("Psame");

    TLegend *leg = new TLegend(0.1,0.85,0.3,0.95," ","brNDC");
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.05);leg->SetNColumns(3);
    leg->AddEntry(gr,"mean time","p"); 
    //leg->AddEntry(grDelay,"Time Delay OCDB","p"); 
    leg->Draw(); 
 
    c1->SaveAs(Form("meanTimePMT%d_vs_run.gif",ipmt));
  }


  //-----> draw your new trending plot here
}

