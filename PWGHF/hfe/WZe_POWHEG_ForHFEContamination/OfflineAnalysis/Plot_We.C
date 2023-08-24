void Plot_We()
{

 
 TCanvas *c1 = new TCanvas();
 gPad->SetLogy();

  // Ze
  //TFile *file0 = TFile::Open("We_pp_all.root");
  TFile *file0 = TFile::Open("We_pp.root");
 // TH1F *ze=(TH1F*) file0->Get("ze");
  Double_t norm = 2.0*acos(-1.0)*1.2*15337.0/1.17e-05;

  ze->Scale(1.0/19.0); //scale # of files
  ze->Scale(1.0/norm);
  /*
  Double_t xbins[14] = {1.5,2,3,4,5,6,8,10,12,14,16,20,26,35};
  ze->Rebin(13,"ze_rebin",xbins);
  ze_rebin->Scale(1,"width");

  ze_rebin->Draw();
  */
 
  TH1D *we = ze->Clone("we"); 

  ze->Draw();
  /*
  TFile *fout = new TFile("Wp_e.root","recreate");  
  //ze_rebin->Write("Wpe");
  we->Write("Wpe");
  */
  TFile *fout = new TFile("Wm_e.root","recreate");  

  //ze_rebin->Write("Wpe");

  we->Write("Wme");
}
