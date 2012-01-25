//Christine Nattrass, University of Tennessee at Knoxville
//This macro is for investigating the number of particles misidentified by the PID algorithm used for the transverse energy analysis
//Misidentified particles are plotted as a percentage of the *total* number of particles.
//Uses the output of AliAnalysisTaskHadEt

void Misidentified(char *prodname = "Enter Production Name", char *shortprodname = "LHC10d4", bool TPC = true, char *filename="Et.ESD.new.sim.LHC10d4.pp.merged.root"){
  char *myname = "TPCITS";
  if(TPC) myname = "TPC";
  TFile *file =  new TFile(filename);
  if(!file){
    cerr<<"Error, no file found"<<endl;
    return;
  }
  TH2F *all = out2->FindObject(Form("MisidentifiedPIDs%s",myname));
  Float_t totaln = ((TH2F*)out2->FindObject(Form("dEdxAll%s",myname)))->GetEntries();
  all->Scale(100.0/totaln);
  gStyle->SetPalette(1);

  all->GetXaxis()->SetRange(all->GetXaxis()->FindBin(0.5),all->GetXaxis()->FindBin(4.));
  all->GetYaxis()->SetRange(all->GetYaxis()->FindBin(0.0));
  all->SetMarkerSize(2);
  all->GetXaxis()->SetLabelSize(0.0);
  all->GetYaxis()->SetLabelSize(0.0);
  all->GetXaxis()->SetTitleSize(.06);
  all->GetYaxis()->SetTitleSize(.06);
  all->GetYaxis()->SetTitle("PID real");
  all->GetXaxis()->SetTitle("PID identified");


  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",600,400);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.124161);
  c->SetBottomMargin(0.147849);
  c->SetLeftMargin(0.129195);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  all->Draw("colz");
  all->Draw("textsame");

  TLatex *pi = new TLatex(0.937284,-0.85435,"#pi^{#pm}");
  pi->Draw();
  TLatex *p = new TLatex(1.937284,-0.85435,"p,#bar{p}");
  p->Draw();
  TLatex *K = new TLatex(2.937284,-0.85435,"K^{#pm}");
  K->Draw();
  TLatex *E = new TLatex(3.937284,-0.85435,"e^{#pm}");
  E->Draw();
  TLatex *pi2 = new TLatex(0.297751,0.937284,"#pi^{#pm}");
  pi2->Draw();
  TLatex *p2 = new TLatex(0.297751,1.937284,"p,#bar{p}");
  p2->Draw();
  TLatex *K2 = new TLatex(0.297751,2.937284,"K^{#pm}");
  K2->Draw();
  TLatex *E2 = new TLatex(0.297751,3.937284,"e^{#pm}");
  E2->Draw();
  TLatex *other = new TLatex(0.297751,0.937284-1.0,"#mu^{#pm}");
  other->Draw();


  TLatex *tex = new TLatex(-0.0932606,-1.19453,prodname);
 tex->SetTextSize(0.0537634);
 tex->Draw();


  if(TPC){
    c->SaveAs(Form("pics/%s/Misidentified.eps",shortprodname));
    c->SaveAs(Form("pics/%s/Misidentified.png",shortprodname));
  }
  else{
    c->SaveAs(Form("pics/%s/MisidentifiedITS.eps",shortprodname));
    c->SaveAs(Form("pics/%s/MisidentifiedITS.png",shortprodname));
  }
}
