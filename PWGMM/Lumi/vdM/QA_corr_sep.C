#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// accessing the rates

void GetRate(Double_t *rate, const char *rate_name, const char *rate_type, Int_t scan, Int_t scan_type, Int_t bc)
{
  char *file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/%sRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_type,rate_name,scan);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/%sRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_type,rate_name,scan);
  TFile *rate_file = new TFile(file_name);
  TTree *rate_tree = (TTree *) rate_file->Get("Rate");
  rate_tree->ResetBranchAddresses();
  rate_tree->SetBranchAddress("rate",rate);
  rate_tree->GetEntry(bc);
  delete [] file_name;
}

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the corrections to the rate
// Note that it uses fBCT corrected rates.
//-------------------------------------------------------

void QA_corr_vs_sep(Int_t Fill, const char *rate_name, Int_t scan, Int_t scan_type, Int_t bc)
// scan_type: 1 => x-scan; 2 => y-scan

{
  
  // initialize
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // reserve space for rated
  Int_t n_sep = FindNumberSeparations(scan_type, scan);
  Double_t *raw = new Double_t[n_sep];
  Double_t *bkgd = new Double_t[n_sep];
  Double_t *pu = new Double_t[n_sep];
  Double_t *all = new Double_t[n_sep];      

  // get the rates
  GetRate(raw,rate_name,"Raw",scan,scan_type,bc);
  GetRate(bkgd,rate_name,"BkgdCorr",scan,scan_type,bc);
  GetRate(pu,rate_name,"PileupCorr",scan,scan_type,bc);
  GetRate(all,rate_name,"IntensityCorrFBCT",scan,scan_type,bc);  
 
  // get the separations
  char *file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/NomSep_x_Scan_%d.root",g_vdm_Fill,scan);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/NomSep_y_Scan_%d.root",g_vdm_Fill,scan);
  TFile *sep_file = new TFile(file_name);
  TTree *sep_tree = (TTree *) sep_file->Get("Separations");
  Double_t *sep = new Double_t[n_sep];
  sep_tree->ResetBranchAddresses();
  sep_tree->SetBranchAddress("separation",sep);
  sep_tree->GetEntry(bc);

  // print if needed
  // for(Int_t i=0;i<n_sep;i++) cout << i << " " << sep[i] << " " << raw[i] << " " << bkgd[i]<< " " << pu[i]<< " " << all [i] << endl;

  // make ratios
  Double_t *raw_bkgd = new Double_t[n_sep];
  Double_t *bkgd_pu = new Double_t[n_sep];
  Double_t *pu_all = new Double_t[n_sep];  
  for(Int_t i=0;i<n_sep;i++) {
    if (raw[i]>0 && bkgd[i]>0) raw_bkgd[i] = bkgd[i]/raw[i]; else raw_bkgd[i] = -1;
    if (bkgd[i]>0 && pu[i]>0) bkgd_pu[i] = pu[i]/bkgd[i]; else bkgd_pu[i] = -1;
    if (pu[i]>0 && all[i]>0) pu_all[i] = all[i]/pu[i]; else pu_all[i] = -1;    
  }

  // define the limits for the plot
  // --> separation
  Double_t sep_min = 0;
  Double_t sep_max = 0;
  for(Int_t i=0;i<n_sep;i++) {
    if(sep[i]<sep_min) sep_min=sep[i];
    if(sep[i]>sep_max) sep_max=sep[i];
  }
  sep_min*=1.2;
  sep_max*=1.2;  
  // --> ratio
  Double_t ratio_max = 0;
  Double_t ratio_min = 2;  
  for(Int_t i=0;i<n_sep;i++) {
    if(raw_bkgd[i] > -1 && raw_bkgd[i] > ratio_max) ratio_max=raw_bkgd[i];
    if(bkgd_pu[i] > -1 && bkgd_pu[i] > ratio_max) ratio_max=bkgd_pu[i];
    if(pu_all[i] > -1 && pu_all[i] > ratio_max) ratio_max=pu_all[i];    
    if(raw_bkgd[i] > -1 && raw_bkgd[i] < ratio_min) ratio_min=raw_bkgd[i];
    if(bkgd_pu[i] > -1 && bkgd_pu[i] < ratio_min) ratio_min=bkgd_pu[i];
    if(pu_all[i] > -1 && pu_all[i] < ratio_min) ratio_min=pu_all[i];    
  }
  ratio_max *= 1.5;
  ratio_min *= 0.8; 

  // make graphs
  TGraph *gr_raw_bkgd = new TGraph(n_sep,sep,raw_bkgd);
  gr_raw_bkgd->SetMarkerStyle(20);gr_raw_bkgd->SetMarkerColor(1);
  TGraph *gr_bkgd_pu = new TGraph(n_sep,sep,bkgd_pu);
  gr_bkgd_pu->SetMarkerStyle(20);gr_bkgd_pu->SetMarkerColor(2);  
  TGraph *gr_pu_all = new TGraph(n_sep,sep,pu_all);
  gr_pu_all->SetMarkerStyle(20);gr_pu_all->SetMarkerColor(4);  
  
  // plot graphs
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *corr_C = new TCanvas("corr_C","rate versus separation",600,400);
  TH1F* frame = gPad->DrawFrame(sep_min,ratio_min,sep_max,ratio_max);
  frame->SetTitle(";separation (mm); correction factor");
  gr_raw_bkgd->Draw("p,e1,same");
  gr_bkgd_pu->Draw("p,e1,same");
  gr_pu_all->Draw("p,e1,same");  
  TLegend *legend = new TLegend(0.3,0.7,0.7,0.9);
  legend->AddEntry(gr_raw_bkgd,"Raw/Bkgd","p");
  legend->AddEntry(gr_bkgd_pu,"Bkgd/Pileup","p");
  legend->AddEntry(gr_pu_all,"Pileup/Intensity","p");  
  legend->Draw();
  
  // clean up
  delete [] raw;
  delete [] bkgd;
  delete [] pu;
  delete [] all;
  delete [] raw_bkgd;
  delete [] bkgd_pu;
  delete [] pu_all;  

}
