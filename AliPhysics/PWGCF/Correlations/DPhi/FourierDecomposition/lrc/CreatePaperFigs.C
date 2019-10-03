// --- Switches and settings ---
bool saveOutput = 1;
bool resave = 0;  // Sometimes pdf gets corrupted...retry then quit.
const char* tag = "-Dec2"; // identifier appended to filenames
const char* outputDir = "../analysis-output";
const char* inputHistFile = "$MYDATA/gsi07_etamin08-withmeanpt.root";
TString grFile   = Form("%s/root-objects/graphs%s.root",
			outputDir, tag);
TString canvFile = Form("%s/root-objects/canvases%s.root",
			outputDir, tag);
TString pdfFile  = Form("%s/plots/all-figs%s",
			outputDir, tag);

// --- Globals ---
TFile* pwg2File   = new TFile("comparisons/PWG2points.root");
TFile* alexv2File = new TFile("comparisons/v2_Eta04_pass2.root");
TFile* alexv3File = new TFile("comparisons/v3_EP_eta04.root");
TLatex ltx;
TObjArray* cList = new TObjArray();

// --- Functions ---
void PlotCollection(int k, int i, int j, TString opt);
void SaveCanvases(TObjArray* canvases, const char* fileName);
TObjArray* GetObjectsFromFile(TFile& file, TString clname, TString dir = "");
void SaveCanvasesFromFile(const char* rootFile, const char* targetDir, const char* fileType);
TGraphErrors* PWG2v(int n, double cen1, double cen2, double etaGap=1.);

void CreatePaperFigs()
{
  ltx.SetNDC();
  int NCENT = 8;
  bool isBatch = gROOT->IsBatch();
  gROOT->LoadMacro("./common/Utils.C+");
  PlotUtils utils;
  gROOT->LoadMacro("FourierPlus.C+g");
  initialize(inputHistFile);

  if (resave) {
    SaveCanvasesFromFile(canvFile.Data(), Form("%s/plots", outputDir), "pdf");
    return;
  }

  int ncb_set = 6;
  int centbins_set[] = { CentBin(0,2), 
			 CentBin(2,10), 
			 CentBin(10,20), 
			 CentBin(20,30), 
			 CentBin(30,40), 
			 CentBin(40,50) };
  
  // SingleDraw2D(CentBin(0,10), PtBin(2, 2.5), PtBin(1.5, 2));
  // SingleDraw2D(CentBin(0,10), PtBin(3, 4), PtBin(2, 2.5));
  // //  cList->Add(SingleDraw2D(CentBin(0,2), PtBin(2, 2.5), PtBin(1.5, 2)));
  // SingleDraw(CentBin(0,2), PtBin(6, 8), PtBin(1., 1.5), "color");
  // SingleDraw(CentBin(0,10), PtBin(6, 8), PtBin(1., 1.5), "color");

    // cList->Add(Drawv2to5(ncb_set, centbins_set, 0, 999, "") );
  //  DrawAgreement(CentBin(0, 10), "contour");
  //DrawChi2(CentBin(0, 10), "");

  int t1[] = {PtBin(2., 2.5)};
  int a1[] = {PtBin(1.5, 2.)};
  int t2[] = {PtBin(8,15)};
  int a2[] = {PtBin(6, 8)};
  int c2[] = {CentBin(40,50), CentBin(0,20)};
  int c3[] = {CentBin(0,2)};
  
  int ncb = 5;
  int centbins[] = { CentBin(40,50), 
		     CentBin(20,30), 
		     CentBin(10,20), 
		     CentBin(2,10), 
		     CentBin(0,2) };

  // int cbin05[] = { CentBin(0,5) };
  // int cbin010[] = { CentBin(0,10) };
  
  // cList->Add(DrawVnDeltaVsN(1,  t1, 1, a1, 1,cbin05, "05_nobars_ext"));
  // cList->Add(DrawVnDeltaVsN(1,  t1, 1, a1, 1,cbin010, "010_nobars_ext"));
  
  // cList->Add(DrawVnDeltaVsN(1,  t1, 1, a1, ncb, centbins, "nobars_ext")); // all centralities
  // cList->Add(DrawVnDeltaVsN(1,  t1, 1, a1, 1,   c3, "cent02_nobars_ext")); // just 0-2%

  // cList->Add(DrawVnDeltaVsN(1,  t1, 1, a1, ncb, centbins,
  // 			    "sine_ext")); // for sys errors.
  // cList->Add(DrawVnDeltaVsN(1,  t2, 1, a2, 2, c2, "sine_ext"));
  // cList->Add(DrawVnDeltaVsN(1, t2, 1, a2, 2, c2, "nobars_ext"));
  
  // DrawQ(CentBin(0,10), 2);
  // DrawQ(CentBin(0,10), 3);
  //  return;

  // For CERN courier article
  cList->Add(SingleDrawPlain(CentBin(0,2), PtBin(2, 2.5), PtBin(1.5, 2), 
			     "harmonics_sum") );

  int aa1[1] = {PtBin(0.5, 0.75)};
  int aa2[3] = {PtBin(0.5, 0.75), PtBin(1.5, 2), PtBin(2.5, 3)};
  int cc1[1] = {CentBin(0,20)};
  int cc2[5] = {CentBin(0,2), CentBin(2,10), CentBin(10,20),
  CentBin(20,30), CentBin(30,40) };

  //  cList->Add( DrawChi2vsPtTrig(1, aa2, 1, cc1, "") );

  if (1) { // Draw Correlation functions

    if (0) {// Draw 3 CFs overlaid together
      int ptt = PtBin(2., 2.5);
      int pta = PtBin(1.5, 2.);
      TLegend* ml = new TLegend(0.5, 0.58, 0.98, 0.8, "Pb-Pb 2.76 TeV");
      ml->SetFillColor(kNone);
      ml->SetBorderSize(0);
      int col[] = {18, kAzure-9, kBlack, kYellow, kOrange, kRed };
      cList->Add(SingleDrawPlain(CentBin(3,5), ptt, pta, "ALL"));
      TH1* h0 = (TH1*)(gPad->GetListOfPrimitives())->FindObject("h");
      utils.set_hist_props(h0, kNone, kNone, kNone, kOpenCircle, 1.4);
      TH1* hUltCnt[3];
      int cntBins[] = {CentBin(3,5), CentBin(2,3), CentBin(0,1)};
      for (int n=0; n<3; n++) {
	hUltCnt[n] = Hist("ALL", "cA", ptt, pta, cntBins[n], Form("_%d", n));
	utils.set_hist_props(hUltCnt[n], kBlack, col[n], col[n], kFullCircle, 1.4);
	hUltCnt[n]->SetLineWidth(2); 
	hUltCnt[n]->Draw("same");
	TH1* hc = (TH1*)hUltCnt[n]->Clone();
	utils.set_hist_props(hc, kBlack, kBlack, kBlack, kOpenCircle, 1.4);
	hc->Draw("same");
      }
      for (int n=2; n>=0; n--) {
	ml->AddEntry(hUltCnt[n], centLabel(cntBins[n]), "epl");
      }
      ml->Draw();
    }
    if (0) { // "plain" CFs (no lower ratio panel)
      int centEvol[] = { CentBin(0,5), 
			 CentBin(0,1), 
			 CentBin(1,2), 
			 CentBin(2,3), 
			 CentBin(3,4), 
			 CentBin(4,5)};
      for (int k=0;k<6;k++) {
	int i = PtBin(2., 2.5), j = PtBin(1.5, 2);
	int kk = centEvol[k]; 
	
	cList->Add(SingleDrawPlain(kk,i,j, "ALL"));
	TH1* h0 = (TH1*)(gPad->GetListOfPrimitives())->FindObject("h");
	h0->GetYaxis()->SetRangeUser(0.98, 1.03);
	
	cList->Add(SingleDrawPlain(kk,i,j, "harmonics_ALL"));
	TH1* h0 = (TH1*)(gPad->GetListOfPrimitives())->FindObject("h");
	h0->GetYaxis()->SetRangeUser(0.98, 1.03);
      }
    }

    // New: High ptt, low pta...
    PlotCollection(CentBin(0,2), PtBin(6, 8), PtBin(1., 1.5), "2D");
    PlotCollection(CentBin(0,10), PtBin(6, 8), PtBin(1., 1.5), "2D");
  
    PlotCollection(CentBin(0,1),  PtBin(2., 2.5), PtBin(1.5,2.), "");
    PlotCollection(CentBin(0,2),  PtBin(2., 2.5), PtBin(1.5,2.), "");
    PlotCollection(CentBin(0,10), PtBin(2, 2.5),  PtBin(1.5, 2), "2D");
    PlotCollection(CentBin(0,10), PtBin(3.0, 4.), PtBin(2.,2.5), "2D");
    PlotCollection(CentBin(0,20), PtBin(8., 15.), PtBin(6.,8), "2D");
    PlotCollection(CentBin(30,40),PtBin(6., 8.),  PtBin(1,1.5), "");
    cList->Add(SingleDraw2D(CentBin(0,20),PtBin(8, 15), PtBin(6, 8), "zoom"));
    ltx.DrawLatex(0.4, 0.7, "#splitline{zoomed to }{0 < C(#Delta#phi) < 5}");
  }
  //  return;

  if (1) { // v_n{GF} for v1 to v5
    cList->Add(Drawv1to5(ncb_set, centbins_set, 0, PtBin(2.0, 2.5), "fitbelow2.5") );
    cList->Add(Drawv1to5(ncb_set, centbins_set, 0, 999, "") );
    cList->Add(Drawv2to5(ncb_set, centbins_set, 0, 999, "") );
    cout << "Drawv1to5() ok" << endl;
  }

  if (1) { // Global fit plots
    DrawQ(CentBin(0,10), 6); // bug somewhere...w/o this, the next canvas is messed up
    for (int n=1; n<=5; n++) {
      cList->Add( DrawQ(CentBin(0,10), n) );
    }
    for (int cb=0; cb<ncb_set; cb++) {
      for (int n=1; n<=5; n++) {
	int k = centbins_set[cb];
	cout << Form("DrawQ() %s, n = %d / 5", centLabel(k), n) << endl;
	cList->Add( DrawQ(k, n) );
      }
    }
  }
  
  if (1) { // "Power spectrum" plots
    int t1[] = {PtBin(2., 2.5)};
    int t2[] = {PtBin(8,15)};
    int t3[] = {PtBin(6,8)};
    int a1[] = {PtBin(1.5, 2.)};
    int a2[] = {PtBin(6, 8)};
    int a3[] = {PtBin(1., 1.5)};

    int c2[] = {CentBin(40,50), CentBin(0,20)};
    int c3[] = {CentBin(0,2)};

    int ncb = 5;
    int centbins[] = { CentBin(40,50), 
		       CentBin(20,30), 
		       CentBin(10,20), 
		       CentBin(2,10), 
		       CentBin(0,2) };

    int ncb2 = 4;
    int centbins2[] = { CentBin(40,50), 
		       CentBin(20,30), 
		       CentBin(10,20), 
		       CentBin(0,10)};

    
    // VnDelta vs n...
    if (0) {
    // Fig 2,3 with color bars
    cList->Add(DrawVnDeltaVsN(1,t1, 1, a1, ncb, centbins, "")); // all
    cList->Add(DrawVnDeltaVsN(1,t1, 1, a1, 1,   c3, "cent02")); // just 0-2%
    cList->Add(DrawVnDeltaVsN(1,t2, 1, a2, 2, c2, "ext"));
    }
    // And without
    cList->Add(DrawVnDeltaVsN(1,t1, 1, a1, ncb, centbins, "nobars"));
    cList->Add(DrawVnDeltaVsN(1,t1, 1, a1, 1,   c3, "cent02_nobars"));
    cList->Add(DrawVnDeltaVsN(1,t2, 1, a2, 2, c2, "ext_nobars"));

    // up to n=12
    cList->Add(DrawVnDeltaVsN(1,t1, 1, a1, ncb, centbins,"sine_ext")); // sys err.
    cList->Add(DrawVnDeltaVsN(1,t2, 1, a2, 2, c2, "sine_ext"));


    // High ptt x low pta
    cList->Add(DrawVnDeltaVsN(1,a2, 1, a1, ncb2, centbins2, ""));
    cList->Add(DrawVnDeltaVsN(1,t3, 1, a3, ncb2, centbins2, "nobars"));
    cList->Add(DrawVnDeltaVsN(1,t2, 1, a3, ncb2, centbins2, "nobars"));

    // Additional
    cList->Add(DrawVnDeltaVsN(1,t1, 1, a1, ncb2, centbins2, "sine"));


    // vn{GF} vs n...
    cList->Add(DrawGlobalvnVsN(1, t1, ncb, centbins, ""));
    cList->Add(DrawGlobalvnVsN(1, a1, ncb, centbins, ""));
    //    cList->Add(DrawVnDeltaVsN(1, t1, 1, a1, ncb_set, centbins_set, "zoom"));

  }

  if (1) { // VnDelta vs. trigger pt -- 1 canvas/centrality bin
    int cbins[] = { CentBin(0,10) , CentBin(10,20), CentBin(20,30), 
		    CentBin(30,40), CentBin(40,50) };
    int ptabins02[] = {PtBin(0.25, 0.5), PtBin(1.0, 1.5)};
    int ptabins10[] = {PtBin(0.25, 0.5), PtBin(1.0, 1.5), PtBin(2.0, 2.5)};

    // cList->Add(DrawVnVsPtTrig(CentBin(0,2), PtBin(0.25, 0.5),
    // PtBin(1.0, 1.5)));
    cList->Add(DrawVnVsPtTrig(CentBin(0,2), 2, ptabins02, ""));

    for (int i=0; i<5; i++) {

      // cList->Add(DrawVnVsPtTrig(cbins[i], PtBin(0.25, 0.5), PtBin(2.5,
      // 3.0)));
      cList->Add(DrawVnVsPtTrig(cbins[i], 3, ptabins10, ""));
    }
  }

 if (1) { // individual v_n{GF} canvases
    // v1 comparison with Luzum
    int ncbv1 = 2;
    int centbinsv1[] = {CentBin(0,10), CentBin(40,50)};

    //    Just our points
    for (int n=1; n<=5; n++) {
      cList->Add(DrawVnFromGlobalFit(n,0, 999, ncb_set,centbins_set,""));
      cList->Add(DrawVnFromGlobalFit(n,0, PtBin(2.0, 2.5), ncb_set,centbins_set,""));
    }
    int ncb = 2;
    int centbins[] = {CentBin(0,2), CentBin(30, 40)};
    for (int n=2; n<=5; n++) {
      TCanvas* c = DrawVnFromGlobalFit(n,0,999,ncb,centbins,"_pwg2");
      TGraphErrors* g0002 = PWG2v(n, 0,  2 );
      TGraphErrors* g3040 = PWG2v(n, 30, 40);
      utils.set_tgraph_props(g0002, kBlack, kBlack, kOpenSquare, 1.5);
      utils.set_tgraph_props(g3040, kRed,   kRed,   kFullSquare, 0.8);
      g0002->Draw("epsame");
      g3040->Draw("epsame");
      if (0)
	if(n==2 || n==3)
	  Alexv(n, 30, 40)->Draw("epsame"); // Great agmt, but I need <pt>
      // ltx.SetTextSize(0.04);
      // ltx.DrawLatex(0.17, 0.95, Form("Global fit"));
      // ltx.DrawLatex(0.45, 0.95, Form("S.P. (CERN-PH-EP-2011-073)"));
      // ltx.SetTextSize(0.05);
      // ltx.DrawLatex(0.45, 0.88, Form("1.0 < |#Delta#eta| < 1.6", n));
      TObject* obj = c->GetListOfPrimitives()->FindObject("lv");
      TLegend* leg = obj;
      leg->AddEntry(g0002, Form("SP 0-2%%", n), "epl");
      leg->AddEntry(g3040, Form("SP 30-40%%", n), "epl");
      cList->Add(c);
    }
  }

  if (saveOutput) { 
    SaveGraphs(grFile.Data());
    SaveCanvases(cList, canvFile.Data()); 

    if (1) // individual pdfs + multipage pdf
      SaveCanvasesFromFile(canvFile.Data(), Form("%s/plots", outputDir), "pdf");
    if (0) // individual plot macros
      SaveCanvasesFromFile(canvFile.Data(), Form("%s/plots", outputDir), "C");

  }
  
  return; // end CreatePaperFigs()
}

TGraphErrors* Alexv(int n, double cen1, double cen2)
{
  // 0-5% -> 0, 5-10% ->1, 10-20%->2
  int bin = -1;
  if (cen1==0  && cen2==5)   bin = 0;
  else if (cen1==5  && cen2==10)  bin = 1;
  else if (cen1==10  && cen2==20) bin = 2;
  else if (cen1==20  && cen2==30) bin = 3;
  else if (cen1==30  && cen2==40) bin = 4;
  else if (cen1==40  && cen2==50) bin = 5;
  else if (cen1==50  && cen2==60) bin = 6;
  else if (cen1==60  && cen2==70) bin = 7;
  else if (cen1==70  && cen2==80) bin = 8;

  if (bin<0)
    return 0;

  TFile* f = 0;
  if (n==2)
    f = alexv2File;
  else if (n==3)
    f = alexv3File;
  else
    return 0;

  const char* xtra = cen1==30 ? "ALICE" : "";
  TGraphErrors* g = 0;
  const char* name = "";

  if (n==2)
    name = Form("v2_all_Eta04_hist_%d", bin);
  if (n==3)
    name = Form("v3_eta04_%d", bin);

   if (0) cout << name << endl;

  g = (TGraphErrors*)f->Get(name);
  if (!g)
    Error("Alexv()", "%s not found");

  return g;
}

TGraphErrors* PWG2v(int n, double cen1, double cen2, double etaGap)
{
  const char* xtra = cen1==30 ? "ALICE" : "";
  TGraphErrors* g = 0;
   const char* name = Form("GrSP_%02d%02d%s_v%d_etaGap%.0f", 
			   (int)cen1, (int)cen2, xtra, n, 10*etaGap);

   if (0) cout << name << endl;

  g = (TGraphErrors*)pwg2File->Get(name);
  if (!g)
    Error("PWG2v()", "%s not found");
  
  return g;
}

void SaveCanvases(TObjArray* canvases, const char* fileName)
{
  TFile* f = new TFile(fileName, "recreate");

  for (int n=0; n<canvases->GetEntries(); n++) {
    TCanvas* c = (TCanvas*)canvases->At(n);
    c->Write(c->GetTitle());
  }
  f->Close();
  return;
}

void PlotCollection(int k, int i, int j, TString opt)
{
  // cList->Add(SingleDrawPlain(k, i, j, "pl"));
  cList->Add(SingleDrawPlain(k, i, j, "harmonics"));
  cList->Add(SingleDraw(k, i, j, "color"));
  //  cList->Add(SingleDraw(k, i, j, "gray"));
  //cList->Add(SingleDraw(k, i, j, "gray_upto10"));
  cList->Add(SingleDraw(k, i, j, "global"));
  if (opt.Contains("2D")) {
    cList->Add(SingleDraw2D(k, i, j, "pl"));
  }
  return;
}

void SaveCanvasesFromFile(const char* rootFile, const char* targetDir, const char* fileType)
{
  // Get a list of canvases from rootFile into array, then save each
  // to its own file in targetDir/. fileType = "eps", "pdf", "C",
  // "png", etc. Not all formats have been tested.
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TString name = "";
  TString base(targetDir);
  TFile *cFile = new TFile(rootFile, "read");
  cout << cFile->GetName() << endl;
  TObjArray* cList = GetObjectsFromFile(*cFile, "TCanvas");

  for (int n=0; n<cList->GetEntries(); n++) {
    TCanvas* c = (TCanvas*)cList->At(n);
    if (c) {
      name = "";
      name = base;
      name += TString("/");
      name += TString(fileType);
      name += TString("/");
      name += TString(c->GetTitle());
      name += TString(".");
      name += TString(fileType);
      cout<<name.Data()<<endl;

      c->Draw();
      c->Modified();
      c->Update();
      c->SaveAs(name.Data());
    }
    else
      Error("SaveCanvasesFromFile()", "!c");
  }

  if (1) {
    utils.print_pdf(cList, Form("%s/all-figs%s", base.Data(), tag), "pdf");
  }
  
  return;
}

TObjArray* GetObjectsFromFile(TFile& file, TString clname, TString dir)
{
  file.cd(dir.Data());

  TObjArray* objList = new TObjArray();
  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  
  while ((key=(TKey*)next())) {
    TString className(key->GetClassName());
    TString keyName(key->GetName());
    if (1) 
      printf("%10s %20s\n", className.Data(), keyName.Data());
    
    if (className.Contains(clname)) {
      objList->Add(gDirectory->Get(keyName.Data()));
    }
  }

  cout << objList->GetEntries() << " objects retrieved from "
       << file.GetName()  << "/" << gDirectory->GetName() 
       << endl;

  return objList;
}
