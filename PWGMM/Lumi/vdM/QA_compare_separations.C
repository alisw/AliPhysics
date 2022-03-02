
//
// QA macro to compare the different separations as a
// function of the nominal separation
//

//-------------------------------------------------------
// headers

#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// accessing the separations

void GetSeparation(Double_t *sep, const char *separation_name, Int_t scan, Int_t scan_type, Int_t bc)
{
  char *file_name = new char[kg_string_size];
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/%sSep_x_Scan_%d.root",g_vdm_Fill,separation_name,scan);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/%sSep_y_Scan_%d.root",g_vdm_Fill,separation_name,scan);
  TFile *separation_file = new TFile(file_name);
  TTree *separation_tree = (TTree *) separation_file->Get("Separations");
  separation_tree->ResetBranchAddresses();
  separation_tree->SetBranchAddress("separation",sep);
  separation_tree->GetEntry(bc);
  delete [] file_name;
}

//-------------------------------------------------------
// entry point

void QA_compare_separations(Int_t Fill, Int_t scan, Int_t scan_type, Int_t bc)
// scan_type: 1 => x-scan; 2 => y-scan
{
  // initialize
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();
  
  // reserve space for separations
  Int_t n_sep = FindNumberSeparations(scan_type, scan);
  Double_t *nom = new Double_t[n_sep];
  Double_t *odc = new Double_t[n_sep];
  Double_t *bbd = new Double_t[n_sep];
  Double_t *hys = new Double_t[n_sep];

  // get the separations
  GetSeparation(nom,"Nom",scan,scan_type,bc);
  GetSeparation(odc,"ODC",scan,scan_type,bc);
  GetSeparation(bbd,"ODCBBD",scan,scan_type,bc);
  GetSeparation(hys,"ODCBBDHyst",scan,scan_type,bc);  
  
  // make differences
  Double_t *odc_nom = new Double_t[n_sep];
  Double_t *bbd_odc = new Double_t[n_sep];
  Double_t *hys_bbd = new Double_t[n_sep];
  for(Int_t i=0;i<n_sep;i++) {
    // factor of 1000 to go to mum
    odc_nom[i] = (odc[i]-nom[i])*1000;
    bbd_odc[i] = (bbd[i]-odc[i])*1000;
    hys_bbd[i] = (hys[i]-bbd[i])*1000;    
  }

  // define the limits for the plot
  // --> separation
  Double_t sep_min = 0;
  Double_t sep_max = 0;
  // --> separation diffs
  Double_t dif_min = 0;
  Double_t dif_max = 0;
  for(Int_t i=0;i<n_sep;i++) {
    if(nom[i]<sep_min) sep_min=nom[i];
    if(nom[i]>sep_max) sep_max=nom[i];
    if(odc_nom[i]<dif_min) dif_min=odc_nom[i];
    if(odc_nom[i]>dif_max) dif_max=odc_nom[i];
    if(bbd_odc[i]<dif_min) dif_min=bbd_odc[i];
    if(bbd_odc[i]>dif_max) dif_max=bbd_odc[i];
    if(hys_bbd[i]<dif_min) dif_min=hys_bbd[i];
    if(hys_bbd[i]>dif_max) dif_max=hys_bbd[i];
  }
  sep_min*=1.2;
  sep_max*=1.2;  
  dif_min*=1.2;
  dif_max*=1.2;  
  
  // make graphs
  TGraph *gr_odc_nom = new TGraph(n_sep,nom,odc_nom);
  gr_odc_nom->SetMarkerStyle(20);gr_odc_nom->SetMarkerColor(1);
  TGraph *gr_bbd_odc = new TGraph(n_sep,nom,bbd_odc);
  gr_bbd_odc->SetMarkerStyle(20);gr_bbd_odc->SetMarkerColor(2);
  TGraph *gr_hys_bbd = new TGraph(n_sep,nom,hys_bbd);
  gr_hys_bbd->SetMarkerStyle(20);gr_hys_bbd->SetMarkerColor(4);

  // plot graphs
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *sep_C = new TCanvas("sep_C","separation differences versus separation",600,400);
  TH1F* frame = gPad->DrawFrame(sep_min,dif_min,sep_max,dif_max);
  frame->SetTitle(";separation (mm); separation difference (#mum)");
  gr_odc_nom->Draw("p,e1,same");
  gr_bbd_odc->Draw("p,e1,same");
  gr_hys_bbd->Draw("p,e1,same");  
  TLegend *legend = new TLegend(0.15,0.7,0.35,0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(gr_odc_nom,"odc-nom","p");
  legend->AddEntry(gr_bbd_odc,"bbd-odc","p");
  legend->AddEntry(gr_hys_bbd,"hys-bbd","p");  
  legend->Draw();

  // clean up
  delete [] nom;
  delete [] bbd;
  delete [] odc;
  delete [] hys;
  delete [] odc_nom;
  delete [] bbd_odc;
  delete [] hys_bbd;
  
}
