// --- Global Switches and settings ---
bool saveOutput = 1;
bool resave = 1;  // Sometimes pdf gets corrupted...retry then quit.
const char* tag = "-v3.4"; // identifier appended to filenames
const char* outputDir = "../test";
const char* inputHistFile = "$MYDATA/gsi07_etamin08-withmeanpt.root";
TString grFile   = Form("%s/root-objects/graphs%s.root",
			outputDir, tag);
TString canvFile = Form("%s/root-objects/canvases%s.root",
			outputDir, tag);
TString pdfFile  = Form("%s/plots/all-figs%s",
			outputDir, tag);
int NCENT = 8;
TLatex ltx;
TObjArray* cList = new TObjArray();
bool isBatch = gROOT->IsBatch();

void PlotCollection(int k, int i, int j, TString opt);

int Setup()
{
  ltx.SetNDC();
  gROOT->LoadMacro("~/Dropbox/ALICE/common/Utils.C+");
  gROOT->LoadMacro("~/Dropbox/ALICE/common/IOUtilFns.C");
  gROOT->LoadMacro("FourierPlus.C+g");
  initialize(inputHistFile);
  return 0;
}

void RunFourierPlus()
{
  Setup();

  if (resave) {
    SaveCanvasesFromFile(canvFile.Data(), Form("%s/plots", outputDir), tag, "pdf");
    return;
  }

  int t1[] = {PtBin(2., 2.5)};
  int t2[] = {PtBin(8,15)};
  int t3[] = {PtBin(6,8)};
  int t4[] = {PtBin(4,5)};
  int a1[] = {PtBin(1.5, 2.)};
  int a2[] = {PtBin(6, 8)};
  int a3[] = {PtBin(1., 1.5)};
  int a4[] = {PtBin(3, 4)};
  int a5[] = {PtBin(0.5., 0.75)};
  int c2[] = {CentBin(40,50),CentBin(0,20)};
  int c3[] = {CentBin(0,2)};
  int c4[] = {CentBin(0,10)};
  int c5[] = {CentBin(0,5)};

  int ncb_set = 5;
  int centbins_set[] = { CentBin(0,2), 
			 CentBin(2,10), 
			 CentBin(10,20), 
			 CentBin(20,30), 
			 CentBin(40,50) };
  int ncb_1 = 6;
  int centbins_1[] = { CentBin(0,2), 
		       CentBin(2,10), 
		       CentBin(10,20), 
		       CentBin(20,30), 
		       CentBin(30,40), 
		       CentBin(40,50) };
  int ncb_2 = 2;
  int centbins_2[] = { CentBin(0,20), 
		       CentBin(40,50) };
  int ncb_3 = 2;
  int centbins_3[] = { CentBin(0,10), 
		       CentBin(40,50) };
  int ncb_4 = 5;
  int centbins_4[] = { CentBin(40,50), 
		     CentBin(20,30), 
		     CentBin(10,20), 
		     CentBin(2,10), 
		     CentBin(0,2) };
  
  // For CERN courier article
  cList->Add(SingleDrawPlain(CentBin(0,2), PtBin(2, 2.5), PtBin(1.5, 2), 
			     "harmonics_sum") );
  
  // Correlation functions
  PlotCollection(CentBin(0,1),  PtBin(2., 2.5), PtBin(1.5,2.), "");
  PlotCollection(CentBin(0,2),  PtBin(2., 2.5), PtBin(1.5,2.), "");
  PlotCollection(CentBin(0,10), PtBin(2, 2.5),  PtBin(1.5, 2), "2D");
  PlotCollection(CentBin(0,10), PtBin(3.0, 4.), PtBin(2.,2.5), "2D");
  PlotCollection(CentBin(0,20), PtBin(8., 15.), PtBin(6.,8), "2D");
  PlotCollection(CentBin(30,40),PtBin(6., 8.),  PtBin(1,1.5), "");
  cList->Add(SingleDraw2D(CentBin(0,20),PtBin(8, 15), PtBin(6, 8), "zoom"));
  ltx.DrawLatex(0.4, 0.7, "#splitline{zoomed to }{0 < C(#Delta#phi) < 5}");
  
  if (1) { // v_n{GF} for v1 to v5
    cList->Add(Drawv1to5(ncb_1, centbins_1, 0, PtBin(2.0, 2.5), "fitbelow2.5") );
    cList->Add(Drawv1to5(ncb_1, centbins_1, 0, 999, "") );
    cList->Add(Drawv2to5(ncb_1, centbins_1, 0, 999, "") );
    cout << "Drawv1to5() ok" << endl;
  }
  //  cList->Add(Drawv2to5(ncb_set, centbins_set, 0, 999, "") );  

  if (1) { // v_n{GF} for v1 to v5 - superimposed with two fit ranges
    cList->Add(Drawv1to5(ncb_2, centbins_2, PtBin(5., 6.), 999));
    cList->Add(Drawv1to5(ncb_2, centbins_2, 0, PtBin(3., 4.), "split") );
    cList->Add(Drawv1to5(ncb_2, centbins_2, 0, PtBin(4., 5.), "split") );
    cout << "Drawv1to5() ok" << endl;
  }

  if (1) { // Global fit plots
    DrawQ(CentBin(0,10), 6); // bug somewhere...w/o this, the next canvas is messed up
    for (int n=1; n<=5; n++) {
      cList->Add( DrawQ(CentBin(0,10), n) );
    }
    for (int cb=0; cb<ncb_set; cb++) {
      for (int n=1; n<=5; n++) {
	int k = centbins_1[cb];
	cout << Form("DrawQ() %s, n = %d / 5", centLabel(k), n) << endl;
	cList->Add( DrawQ(k, n) );
      }
    }
  }
  if (1) { // Global fit plots at high pt only
    for (int n=1; n<=5; n++) {
      cout << Form("DrawQ() %s, n = %d / 5", "highptfit", n) << endl;
      cList->Add( DrawQ(CentBin(0,10), n, "highptfit") );   // pta > 5
      cList->Add( DrawQ(CentBin(0,20), n, "highptfit") );   // pta > 5
    }
    cout << "DrawQ() ok" << endl;
  }

  cList->Add( DrawQ(CentBin(0,20), 1, 
   		    Form("ptabin%dto%d_%s",PtBin(0.25,0.5),PtBin(0.75,1.0), "RIDGE")));
  cList->Add( DrawQ(CentBin(0,20), 1, 
		    Form("ptabin%dto%d_%s",PtBin(2.,2.5),PtBin(3,4), "RIDGE")));
  
  // individual v_n{GF} canvases
  cList->Add(DrawVnFromGlobalFit(1,0,999,ncb_3,centbins_3,"multi"));
  cList->Add(DrawVnFromGlobalFit(2,0,999,ncb_3,centbins_3,"multi"));

  // VnDelta vs n...
  cList->Add(DrawVnDeltaVsN(1,t1, 1, a1, ncb_4, centbins_4, "nobars"));
  cList->Add(DrawVnDeltaVsN(1,t1, 1, a1, 1,   c3, "cent02_nobars"));
  cList->Add(DrawVnDeltaVsN(1,t2, 1, a2, 2, c2, "ext_nobars"));
  cout << "DrawVnDeltaVsN() ok" << endl;
  
  // // vn{GF} vs n...
  // cList->Add(DrawGlobalvnVsN(1, t1, ncb_4, centbins_4, ""));
  // cList->Add(DrawGlobalvnVsN(1, a1, ncb_4, centbins_4, ""));

  if (1) { // VnDelta vs. trigger pt -- 1 canvas/centrality bin
    int cbins[] = { CentBin(0,10) , CentBin(10,20), CentBin(20,30), 
		    CentBin(30,40), CentBin(40,50) };
    int ptabins02[] = {PtBin(0.25, 0.5), PtBin(1.0, 1.5)};
    int ptabins10[] = {PtBin(0.25, 0.5), PtBin(1.0, 1.5), PtBin(2.0, 2.5)};

    cList->Add(DrawVnVsPtTrig(CentBin(0,2), 2, ptabins02, ""));

    for (int i=0; i<5; i++) {
      cList->Add(DrawVnVsPtTrig(cbins[i], 3, ptabins10, ""));
    }
      cout << Form("DrawVnVsPtTrig() ok") << endl;
  }

  if (saveOutput) {
    SaveGraphs(grFile.Data());
    cout << Form("SaveGraphs() ok") << endl;
    SaveCanvases(cList, canvFile.Data()); 
    cout << Form("SaveCanvases() ok") << endl;

    if (1) // individual pdfs + multipage pdf
      SaveCanvasesFromFile(canvFile.Data(), Form("%s/plots", outputDir), tag, "pdf");
    if (0) // individual plot macros
      SaveCanvasesFromFile(canvFile.Data(), Form("%s/plots", outputDir), tag, "C");
  }

  return;
}


void PlotCollection(int k, int i, int j, TString opt)
{
  cList->Add(SingleDrawPlain(k, i, j, "harmonics"));
  cList->Add(SingleDraw(k, i, j, "color"));
  cList->Add(SingleDraw(k, i, j, "global"));
  if (opt.Contains("2D")) {
    cList->Add(SingleDraw2D(k, i, j, "pl"));
  }
  return;
}
