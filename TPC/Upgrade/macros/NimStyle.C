Int_t kmicolors[10]={1,2,3,4,6,7,8,9,10,11};
Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};


void NimStyle(){
  TStyle *stNimTPC = new TStyle("tpcNimStyle","tpcNimStyle");
  gROOT->GetStyle("Plain")->Copy((*stNimTPC));
  Int_t nimTPCFont=22; //or 62 for sans serif font
  //default definitions
  stNimTPC->SetTextFont(nimTPCFont);
  stNimTPC->SetTitleFont(nimTPCFont, "T");
  stNimTPC->SetTitleFont(nimTPCFont, "XYZ");
  stNimTPC->SetLabelFont(nimTPCFont,"XYZ");
  stNimTPC->SetStatFont(nimTPCFont);
  stNimTPC->SetOptTitle(0);
  stNimTPC->SetNumberContours(50);
  stNimTPC->SetPalette(1,0);
  stNimTPC->SetStatBorderSize(1);
  new TColor(101,1,1,1);
  stNimTPC->SetFillColor(101);
  if (nimTPCFont==22){
    stNimTPC->SetLabelSize(0.045,"XYZ");
    stNimTPC->SetTitleSize(0.05,"XYZ");
    stNimTPC->SetTitleOffset(0.92,"XYZ");
  }
  //depending on taste
  stNimTPC->SetTitleX(0.1);
  stNimTPC->SetTitleW(0.8);
  stNimTPC->SetTitleH(0.08);
  stNimTPC->SetStatX(.9);
  stNimTPC->SetStatY(.9);
  stNimTPC->SetOptStat(1100);
  stNimTPC->cd();
  //
  //
}


