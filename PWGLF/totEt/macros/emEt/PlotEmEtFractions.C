void SetStyles(TH1 *histo,int marker, int color){
  histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  //histo->GetXaxis()->SetTitle(xtitle);
  //histo->GetYaxis()->SetTitle(ytitle);
}
void PlotEmEtFractions(Bool_t isPhos = kFALSE){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  Float_t min = 0;
  float max = 1;
  TString filename, detname;
  if(isPhos){
    min = 0.655;
    max = 0.785;
    detname = "PHOS";
    filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.Run139465.root";
  }
  else{
    min = 0.58;
    max = 0.725;
    filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root";
    detname = "EMCal";
  }
  
  TFile *f = TFile::Open(filename, "READ");
  TList *l = dynamic_cast<TList*>(f->Get("out1"));
  TString part[] = {"Signal","Hadrons","Neutrons","Kaons","Secondaries"};
  Int_t nparts = 5;
  TString var[] = {"NClusters","NMultiplicity","NMatchedTracks","NTotalTracks"};
  Int_t nvar = 4;
  Int_t partnum = 0;
  Int_t varnum = 0;

  for(Int_t partnum = 0; partnum<nparts;partnum++){
    for(Int_t varnum = 0;varnum<nvar; varnum++){
      TString histoname = "fHistFrac"+part[partnum]+"Vs"+var[varnum];
      
      TH2F *histo = l->FindObject(histoname.Data());
      histo->GetXaxis()->SetTitle(var[varnum].Data());
      TString ytitle = part[partnum]+" fraction";
      histo->GetYaxis()->SetTitle(ytitle.Data());
      //this line here keeps the zeros from plotting
      histo->GetYaxis()->SetRange(2,histo->GetYaxis()->GetNbins());
      
      TCanvas *c1 = new TCanvas("c1","Simulation",600,400);
      c1->SetTopMargin(0.02);
      c1->SetRightMargin(0.03);
      c1->SetLeftMargin(0.11745);
      c1->SetBottomMargin(0.11745);
      c1->SetBorderSize(0);
      c1->SetFillColor(0);
      c1->SetFillColor(0);
      c1->SetBorderMode(0);
      c1->SetFrameFillColor(0);
      c1->SetFrameBorderMode(0);
      histo->Draw("colz");
      TH1 * prof = histo->ProfileX();
      prof->Draw("same");
      TString canvasname = "/tmp/"+histoname+detname+".png";
      c1->SaveAs(canvasname.Data());
      //return;
      delete c1;
      delete prof;
    }
  }
}
