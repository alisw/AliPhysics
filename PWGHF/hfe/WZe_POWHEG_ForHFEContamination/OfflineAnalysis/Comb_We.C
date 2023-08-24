void Comb_We()
{
 
 TFile *file0 = TFile::Open("Wp_e.root");
 TH1D *Wp = Wpe->Clone("Wp");
 TFile *file1 = TFile::Open("Wm_e.root");
 TH1D *Wm = Wme->Clone("Wm");

 Wp->Add(Wm);
 Wp->Scale(0.5);
 Wp->Scale(1.0,"width");

 Wp->Draw(); 

 
 TGraphErrors *Gwe = new TGraphErrors(50);
 for(int i=0; i<50; i++)
    {
     Double_t pT = Wp->GetBinCenter(i+1);
     Double_t yWe = Wp->GetBinContent(i+1)/pT;
     //Double_t yWe_err = Wm->GetBinError(i+1)/pT;

     Gwe->SetPoint(i,pT,yWe);   
     Gwe->SetPointError(i,0.0,0.0);   
    }

 Gwe->SetTitle("electrons from W;p_{T}(GeV/c);d#sigma^{2}/(2#pi p_{T} #Delta p_{T} #Delta #eta)");
 Gwe->SetMarkerStyle(20);
 Gwe->SetMarkerColor(2);
 Gwe->Draw("P");

 TFile *fout = new TFile("W_decay_e.root","recreate");
 Gwe->Write("W_decay_e");

}
