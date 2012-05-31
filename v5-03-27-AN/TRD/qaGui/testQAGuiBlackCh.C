
void testQAGuiBlackCh(const char *filename) {
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  

  gSystem->Load("libGui.so");
  gSystem->Load("libTRDGui.so");
  
  TGMainFrame *fMain = new TGMainFrame(gClient->GetRoot(), 400, 400);
  //TGuiDetector *det = new TGuiDetector((TGWindow*)fMain);
  //fMain->AddFrame(det);
  
  //AliTRDqaGuiBlackSM *sm = new AliTRDqaGuiBlackSM((TGWindow*)fMain);
  AliTRDqaGuiBlackChamber *chamber = new AliTRDqaGuiBlackChamber((TGWindow*)fMain);
  //fMain->AddFrame(sm);
  fMain->AddFrame(chamber);
  
  chamber->SetQAFile(filename);
  //sm->SetQAFile(filename);
  
  /**/
  chamber->SetRangePed(8, 11);
  chamber->SetRangeNoise(0.5, 2);
  chamber->SetQAFile(filename);
  /**/

  fMain->SetWindowName("TRD QA");
  fMain->MapSubwindows();
  fMain->MapWindow();
  fMain->Resize(fMain->GetDefaultSize());
}
