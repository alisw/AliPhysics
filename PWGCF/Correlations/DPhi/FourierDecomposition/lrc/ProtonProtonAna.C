// --- Global Switches and settings ---
bool saveOutput = 1;
bool resave = 0;  // Sometimes pdf gets corrupted...retry then quit.
const char* tag = "-pp-nov28-etamin12"; // identifier appended to filenames
const char* outputDir = "../test";
const char* inputHistFile = "$MYDATA/lhc11a_etamin12.root";
TString grFile   = Form("%s/root-objects/graphs%s.root",
			outputDir, tag);
TString canvFile = Form("%s/root-objects/canvases%s.root",
			outputDir, tag);
TString pdfFile  = Form("%s/plots/all-figs%s",
			outputDir, tag);
int NCENT = 1;
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
  SetProtonProtonFlag(1);
  SetMinRidgeDEta(1.2);
  return 0;
}

void ProtonProtonAna()
{
  Setup();
  cout << IsProtonProton() << endl;
  if (resave) {
    SaveCanvasesFromFile(canvFile.Data(), Form("%s/plots", outputDir), tag, "pdf");
    return;
  }

  int nt = 5, na = 6;
  int t1[] = {PtBin(1., 1.5), PtBin(2., 2.5), PtBin(4,5), PtBin(6,8), PtBin(8,15)};
  int a1[] = {PtBin(0.25., 0.5), PtBin(0.5., 0.75), PtBin(1., 1.5), PtBin(2., 2.5), 
	      PtBin(3, 4), PtBin(6, 8)};

  //    cList->Add( DrawQ(0, 1, Form("ptabin%dto%d_%s",PtBin(0.25,0.5),PtBin(4,5), "ptcons")));

    if (1) {
      for (int n=1; n<=5; n++) {
	cList->Add( DrawQ(0, n, Form("ptabin%dto%d_%s",PtBin(0.25,0.5),PtBin(8,15), "RIDGE")));
      }
      for (int n=1; n<=5; n++) {
	cList->Add( DrawQ(0, n, Form("ptabin%dto%d_%s",PtBin(0.25, 0.5),PtBin(1,1.5), "RIDGE")));
      }
    }
    
  int ppcb[] = {CentBin(0,0)};
  cList->Add(Drawv1to5(1, ppcb, PtBin(0.25, 0.5), PtBin(1., 1.5), "ptcons" ));

  if (1) {
    for (int i=0; i<nt; i++) {
      for (int j=0; j<na; j++) {
	int ti = t1[i], aj = a1[j];
	if (ti >= aj)
	  cList->Add(PlotCollection(ti, aj, "2D"));
      }
    }
  }
  

  if (saveOutput) {
    SaveGraphs(grFile.Data());
    SaveCanvases(cList, canvFile.Data()); 

    if (1) // individual pdfs + multipage pdf
      SaveCanvasesFromFile(canvFile.Data(), Form("%s/plots", outputDir), tag, "pdf");
    if (0) // individual plot macros
      SaveCanvasesFromFile(canvFile.Data(), Form("%s/plots", outputDir), tag, "C");
  }

  return;
}

void PlotCollection(int i, int j, TString opt)
{
  //  cList->Add(SingleDraw(0, i, j, ""));
  if (opt.Contains("2D")) {
    cList->Add(SingleDraw2D(0, i, j, "pl"));
  }
  cList->Add(SingleDrawPlain(0, i, j, "harmonics_sum"));
  return;
}
