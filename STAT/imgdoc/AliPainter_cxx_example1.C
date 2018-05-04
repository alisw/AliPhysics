{
TFile::SetCacheFileDir(".");
TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
TTree *tree = (TTree *) finput.Get("hisPtAll");
hisArray = new TObjArray();
TList *keys = finput->GetListOfKeys();
for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
  TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
  hisArray->AddLast(o);
}

TCanvas *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);
AliPainter::DivideTPad("<horizontal>[1,1,1,1]", "Canvas41", "", canvasQA);
canvasQA->cd(1);
AliPainter::DrawHistogram("hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
AliPainter::DrawHistogram("hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)", hisArray);
AliPainter::DrawHistogram("hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)", hisArray);
AliPainter::DrawHistogram("hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)", hisArray);
AliDrawStyle::RegisterCssStyle("figTemplateHex", AliDrawStyle::ReadCSSFile("../test/figTemplateHex.css"));
AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
}