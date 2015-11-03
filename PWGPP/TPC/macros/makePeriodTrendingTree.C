makePeriodTrendingTree(const char* inFile = "trending.root", const char* runType="PbPb")
{

 TTreeSRedirector pcstream("periodTrending.root");
 TCut selection;
 TObjString runTypeStr(runType);
 TFile* f = TFile::Open(inFile);
 gROOT->cd();
 TTree *t = (TTree*)f->Get("trending");
 TStatToolkit::MakeSummaryTree(t, &pcstream, runTypeStr, selection);
 delete t;
 delete f;
}
