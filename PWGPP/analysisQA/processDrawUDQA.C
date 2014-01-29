void processDrawUDQA(const char *infile="AnalysisResults.root", 
		     TString suffix = "eps", 
		     const char *outfile="DrawUDQA_outfile.root") {

TFile *f = TFile::Open(infile, "read");

if(!f) {
 Printf("FATAL: File doesn't exist");
 return;
}

TList *flistQA = (TList*)f->Get("UpcTree/ListHist");

if(!flistQA) {
 Printf("FATAL: TList - UpcTree/ListHist -  doesn't exist");
 return;
}



TH1D *fHistNeventsJPsi = dynamic_cast<TH1D*> (flistQA->At(0));
fHistNeventsJPsi->GetYaxis()->SetTitle("Counts (Events)");

TH2D *fHistTPCsignalJPsi = dynamic_cast<TH2D*> (flistQA->At(1));
fHistTPCsignalJPsi->GetXaxis()->SetTitle("d#it{E}/d#it{x}^{TPC} (#it{l}^{+}) (a.u.)");
fHistTPCsignalJPsi->GetYaxis()->SetTitle("d#it{E}/d#it{x}^{TPC} (#it{l}^{-}) (a.u.)");

TH2D *fHistDiLeptonPtJPsi = dynamic_cast<TH2D*> (flistQA->At(2));
fHistDiLeptonPtJPsi->GetXaxis()->SetTitle("#it{p}_{T} (#it{l}^{+}) (GeV/#it{c})");
fHistDiLeptonPtJPsi->GetYaxis()->SetTitle("#it{p}_{T} (#it{l}^{-}) (GeV/#it{c})");

TH1D *fHistDiElectronMass = dynamic_cast<TH1D*> (flistQA->At(3));
fHistDiElectronMass->GetYaxis()->SetTitle("dN/dM");

TH1D *fHistDiMuonMass = dynamic_cast<TH1D*> (flistQA->At(4));
fHistDiMuonMass->GetYaxis()->SetTitle("dN/dM");

//___________________________________
  // Added by sjena
TFile *fout = TFile::Open(outfile,"UPDATE");
  fout->ls();
  
  TDirectoryFile *cdd = NULL;
  cdd = (TDirectoryFile*)fout->Get("UD");
  if(!cdd) {
    Printf("Warning: UD <dir> doesn't exist, creating a new one");
    cdd = (TDirectoryFile*)fout->mkdir("UD");
  }
  cdd->cd();
  cdd->ls();

  fHistNeventsJPsi->Write();
  fHistTPCsignalJPsi->Write();
  fHistDiLeptonPtJPsi->Write();
  fHistDiElectronMass->Write();
  fHistDiMuonMass->Write();
  fout->Close();
//___________________________________




  myOptions();
  gROOT->ForceStyle();

  TDatime now;
  int iDate = now.GetDate();
  int iYear=iDate/10000;
  int iMonth=(iDate%10000)/100;
  int iDay=iDate%100;
  char* cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
		    "Jul","Aug","Sep","Oct","Nov","Dec"};
  char cStamp1[25],cStamp2[25];
  sprintf(cStamp1,"%i %s %i",iDay, cMonth[iMonth-1], iYear);
  sprintf(cStamp2,"%i/%.2d/%i",iDay, iMonth, iYear);

myHistSetUp(fHistNeventsJPsi);
fHistNeventsJPsi->SetLineWidth(3.0);
myHistSetUp(fHistTPCsignalJPsi);
myHistSetUp(fHistDiLeptonPtJPsi);
myHistSetUp(fHistDiElectronMass);
fHistDiElectronMass->SetLineWidth(2.0);
myHistSetUp(fHistDiMuonMass);
fHistDiMuonMass->SetLineWidth(2.0);

TCanvas *c1 = new TCanvas("fig_ud_NeventsJPsi"," ",800,400);
c1->Draw();
c1->cd();
TPad *myPad1 = new TPad("myPad1", "The pad",0,0,1,1);
myPadSetUp(myPad1,0.15,0.1,0.04,0.15);
myPad1->SetLogy();
myPad1->Draw();
myPad1->cd();

fHistNeventsJPsi->Draw();
 c1->SaveAs(Form("fig_ud_NeventsJPsi.%s",suffix.Data()));


TCanvas *c21 = new TCanvas("fig_ud_TPCsignalJPsi"," ",500,500);
c21->Draw();
c21->cd();
TPad *myPad21 = new TPad("myPad21", "The pad",0,0,1,1);
myPadSetUp(myPad21,0.15,0.15,0.15,0.15);
myPad21->Draw();
myPad21->cd();
fHistTPCsignalJPsi->Draw("COLZ");
 c21->SaveAs(Form("fig_ud_TPCsignalJPsi.%s",suffix.Data()));


TCanvas *c22 = new TCanvas("fig_ud_DiLeptonPtJPsi"," ",500,500);
c22->Draw();
c22->cd();
TPad *myPad22 = new TPad("myPad22", "The pad",0,0,1,1);
myPadSetUp(myPad22,0.15,0.15,0.15,0.15);
myPad22->Draw();
myPad22->cd();
fHistDiLeptonPtJPsi->Draw("COLZ");
 c22->SaveAs(Form("fig_ud_DiLeptonPtJPsi.%s",suffix.Data()));
 

TCanvas *c31 = new TCanvas("fig_ud_DiElectronMassJPsi"," ",800,500);
c31->Draw();
c31->cd();
TPad *myPad31 = new TPad("myPad31", "The pad",0,0,1,1);
myPadSetUp(myPad31,0.15,0.1,0.04,0.15);
myPad31->Draw();
myPad31->cd();
fHistDiElectronMass->Draw();
 c31->SaveAs(Form("fig_ud_DiElectronMassJPsi.%s",suffix.Data()));
 

TCanvas *c32 = new TCanvas("fig_ud_DiMuonMassJPsi"," ",800,500);
c32->Draw();
c32->cd();
TPad *myPad32 = new TPad("myPad32", "The pad",0,0,1,1);
myPadSetUp(myPad32,0.15,0.1,0.04,0.15);
myPad32->Draw();
myPad32->cd();
fHistDiMuonMass->Draw();
 c32->SaveAs(Form("fig_ud_DiMuonMassJPsi.%s",suffix.Data()));
 

/*
   TFile* outfile = new TFile("fig_ud.root","recreate");
   c1->Write();
   c21->Write();
   c31->Write();
   c22->Write();
   c32->Write();
 */

}

void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}


void myHistSetUp(TH1 *currentGraph=0){
 
  currentGraph->SetLabelSize(0.05,"xyz");
  currentGraph->SetLabelFont(42,"xyz"); 
  currentGraph->SetLabelOffset(0.01,"xyz");
  currentGraph->SetTitleFont(42,"xyz");   
  currentGraph->GetXaxis()->SetTitleOffset(1.1);
  currentGraph->GetYaxis()->SetTitleOffset(1.1);  
  currentGraph->SetTitleSize(0.06,"xyz");
  currentGraph->SetStats(kFALSE); 
  currentGraph->SetTitle("");
  return;
}
void myHistSetUp(TH2 *currentGraph=0){
 
  currentGraph->SetLabelSize(0.05,"xyz");
  currentGraph->SetLabelFont(42,"xyz"); 
  currentGraph->SetLabelOffset(0.01,"xyz");
  currentGraph->SetTitleFont(42,"xyz"); 
  currentGraph->GetXaxis()->SetTitleOffset(1.3);
  currentGraph->GetYaxis()->SetTitleOffset(1.3);  
  currentGraph->SetTitleSize(0.05,"xyz");
  currentGraph->SetStats(kFALSE);
  currentGraph->SetTitle(""); 
  return;
}

void myOptions(Int_t lStat=0){
  // Set gStyle
  int font = 42;
  // From plain
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.0,"xyz");  
  gStyle->SetTitleSize(0.06,"xyz");  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0); 
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }
}

