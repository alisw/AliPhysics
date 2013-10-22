void MakeSensitivityPlots() {
  //
  //  -> THIS MACRO SHOULD BE COMPILABLE.
  //  -> ALL PLOTS SHOULD BE LABELED (ESPECIALLY THE AXES).
  //  -> DATA RED AND MC BLUE.
  //
  TH1D * nclAcceptedData80 = GetAcceptedFractionNclCut(80, "output/Data_LHC13b.root");
  nclAcceptedData80->SetNameTitle("nclAcceptedData80","nclAcceptedData80");
  nclAcceptedData80->SetLineColor(kRed -3);
  //
  TH1D * nclAcceptedMc80   = GetAcceptedFractionNclCut(80, "output/MC_LHC13b.root");
  nclAcceptedMc80->SetNameTitle("nclAcceptedMc80","nclAcceptedMc80");
  nclAcceptedMc80->SetLineColor(kBlue -3);
  //
  TCanvas * canvNclCut = new TCanvas("canvNclCut","sensitivity to ncl cut",600,800);
  canvNclCut->Divide(1,2);
  canvNclCut->cd(1)->SetLogx();
  //
  nclAcceptedData80->DrawCopy();
  nclAcceptedMc80->DrawCopy("SAME");
  //
  //
  //
  canvNclCut->cd(2)->SetLogx();;
  nclAcceptedData80->Divide(nclAcceptedMc80);
  nclAcceptedData80->DrawCopy();

}


TH1D * GetAcceptedFractionNclCut(Int_t nclCut = 80, const Char_t * inFileName = "output/Data_LHC13b.root") {
  //
  // accepted fraction of tracks for ncl cut vs. pT
  //
  TFile * inFileData = TFile::Open(inFileName);
  TList * l = (TList * ) inFileData->Get("akalweit_TrackingUncert");
  THnF * histNcl = (THnF *) l->FindObject("histNcl");
  //  histNcl->GetListOfAxes()->Print();
  //
  // determine sensitivities
  //
  TH1D * hAll = histNcl->Projection(1);
  hAll->SetNameTitle("hAll","hAll");  
  //
  const Int_t kVeryBig = 10000;
  histNcl->GetAxis(0)->SetRangeUser(nclCut, kVeryBig);
  TH1D * hAccepted = histNcl->Projection(1);
  hAccepted->SetNameTitle("hAccepted","hAccepted");
  //
  histNcl->GetAxis(0)->SetRangeUser(0,nclCut);
  TH1D * hRejected = histNcl->Projection(1);
  hRejected->SetNameTitle("hRejected","hRejected");
  //
  //
  hAccepted->Divide(hAll);
  hRejected->Divide(hAll);
  //
  // some cosmetics
  //
  hAccepted->SetLineWidth(2);
  //
  return hAccepted;

  
}
