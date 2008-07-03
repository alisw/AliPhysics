
void testQAGui(const char *filename) {
  
  // gROOT->ForceStyle();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
 
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleY(0.95);
  gStyle->SetTitleW(0.45);
  gStyle->SetTitleH(0.05);

  gStyle->SetLabelFont(52, "XY");


  gSystem->Load("libGui.so");
  gSystem->Load("libTRDgui.so");
  
  AliTRDqaGuiMain *fMain = new AliTRDqaGuiMain(gClient->GetRoot());
  fMain->SetQAFile(filename);
}
