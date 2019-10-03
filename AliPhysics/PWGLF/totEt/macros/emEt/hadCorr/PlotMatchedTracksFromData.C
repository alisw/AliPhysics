void PlotMatchedTracksFromData(TString datafilename="rootFiles/LHC11a4_bis/Et.ESD.simPbPb.EMCAL.LHC11a4_bis.root",TString simfilename="rootFiles/LHC11a4_bis/Et.ESD.simPbPb.EMCAL.LHC11a4_bis.root", int bin = 10, int binLast = 10, TString det = "EMCal"){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString outname = "";
  TString outnamebin = Form("%iTo%i",bin,binLast);

  TFile *fsim = TFile::Open(simfilename, "READ");
  TList *lsim = (TList*)fsim->Get("out1");
  TH3F  *fHistMatchedTracksEvspTBkgdvsMult = lsim->FindObject("fHistMatchedTracksEvspTBkgdMult");
  fHistMatchedTracksEvspTBkgdvsMult->GetZaxis()->SetRange(bin,binLast);
  TH2D *hBkgd2D = (TProfile2D*) fHistMatchedTracksEvspTBkgdvsMult->Project3D("yx");
  TProfile * profBkgd2D = hBkgd2D->ProfileX();
  profBkgd2D->SetLineColor(2);
  profBkgd2D->SetLineWidth(2);

  TFile *f = TFile::Open(datafilename, "READ");

  f->cd();  

  TList *l = (TList*)f->Get("out1");
  TH3F  *fHistMatchedTracksEvspTvsMult = l->FindObject("fHistMatchedTracksEvspTvsMult");
  fHistMatchedTracksEvspTvsMult->GetZaxis()->SetRange(bin,binLast);
  TH2D *hMatchedTracks2D = (TProfile2D*) fHistMatchedTracksEvspTvsMult->Project3D("yx");
  TProfile * profMatchedTracks2D = hMatchedTracks2D->ProfileX();
  hMatchedTracks2D->GetXaxis()->SetTitle("p");
  hMatchedTracks2D->GetYaxis()->SetTitle("E^{cluster}");
  profMatchedTracks2D->SetLineWidth(2);
  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  hMatchedTracks2D->Draw("colz");
  profMatchedTracks2D->Draw("same");
  profBkgd2D->Draw("same");
  outname = "/tmp/TrackMatchingData2D"+det+outnamebin+".png";
  c1->SaveAs(outname.Data());
}
