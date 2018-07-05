{
  TFile::SetCacheFileDir(".");
  TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
  TTree *tree = (TTree *) finput->Get("hisPtAll");
  hisArray = new TObjArray();
  TList *keys = finput->GetListOfKeys();
  for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
  TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
  hisArray->AddLast(o);
  }
  TCanvas *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);
  AliPainter::DivideTPad(canvasQA, "<horizontal>[1,1]", "Canvas1");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E,ylim=[<min>+0.5*<mean>,<max>-0.5*<mean>], div=1)");
}
