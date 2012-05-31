
void testQAGuiAnalysis() {
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
 
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleY(0.95);
  gStyle->SetTitleW(0.45);
  gStyle->SetTitleH(0.05);


  gSystem->Load("libGui.so");
  gSystem->Load("libTRDgui.so");
  
  AliTRDqaGuiMainAnalysis *fMain = new AliTRDqaGuiMainAnalysis(gClient->GetRoot());
  fMain->SetQAFile();
}
