#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// ghost and satellite charges
//-------------------------------------------------------
//Global scan time variables
UInt_t gScan0StartX = 0;
UInt_t gScan0EndX   = 0;
UInt_t gScan0StartY = 0;
UInt_t gScan0EndY   = 0;
UInt_t gScan1StartX = 0;
UInt_t gScan1EndX   = 0;
UInt_t gScan1StartY = 0;
UInt_t gScan1EndY   = 0;


// Set the start and end times of the scans
void SetScanTimes(){
    //set set branches
    Double_t time;
    g_vdm_Tree->ResetBranchAddresses();
    g_vdm_Tree->SetBranchAddress("time",&time);
    
    g_vdm_Tree->GetEntry(g_Idx_Start_Scan_x[0]);
    gScan0StartX = (UInt_t) time;
    cout << "Start X0: " << gScan0StartX << endl;
    
    g_vdm_Tree->GetEntry(g_Idx_End_Scan_x[0]);
    gScan0EndX = (UInt_t) time;
    cout << "End X0: " << gScan0EndX << endl;
    
    g_vdm_Tree->GetEntry(g_Idx_Start_Scan_y[0]);
    gScan0StartY = (UInt_t) time;
    cout << "Start Y0: " << gScan0StartY << endl;
    
    g_vdm_Tree->GetEntry(g_Idx_End_Scan_y[0]);
    gScan0EndY = (UInt_t) time;
    cout << "End Y0: " << gScan0EndY << endl;
    
    g_vdm_Tree->GetEntry(g_Idx_Start_Scan_x[1]);
    gScan1StartX = (UInt_t) time;
    cout << "Start X1: " << gScan1StartX << endl;
    
    g_vdm_Tree->GetEntry(g_Idx_End_Scan_x[1]);
    gScan1EndX = (UInt_t) time;
    cout << "End X1: " << gScan1EndX << endl;
    
    g_vdm_Tree->GetEntry(g_Idx_Start_Scan_y[1]);
    gScan1StartY = (UInt_t) time;
    cout << "Start Y1: " << gScan1StartY << endl;
    
    g_vdm_Tree->GetEntry(g_Idx_End_Scan_y[1]);
    gScan1EndY = (UInt_t) time;
    cout << "End Y1: " << gScan1EndY << endl;
}

void Get_satelliteCharge(Int_t Fill, Int_t opt, Bool_t stat, Bool_t save)
// opt: 0 => beam 1 device 1
// opt: 1 => beam 1 device 2
// opt: 2 => beam 2

{
  // create names for files
  char *ghost_file_name = new char[kg_string_size];
  if (opt == 0)
    sprintf(ghost_file_name,"../Fill-%d/ghostCharge/ghostCharge_B1_A.root",
	    Fill);
  else if (opt == 1)
    sprintf(ghost_file_name,"../Fill-%d/ghostCharge/ghostCharge_B1_B.root",
	    Fill);
  else
    sprintf(ghost_file_name,"../Fill-%d/ghostCharge/ghostCharge_B2_A.root",
	    Fill);

  //Open files and get trees
  TFile *ghost_file = new TFile(ghost_file_name);
  TTree *satellite_tree = (TTree *) ghost_file->Get("SatelliteCharge");

  // reserve space
  Int_t ng =  satellite_tree->GetEntries();
  Double_t sat;
  Double_t satErr;
  UInt_t time;
  satellite_tree->SetBranchAddress("Q_sat_value",&sat);
  if (stat)
    satellite_tree->SetBranchAddress("sigma_sat_value_stat",&satErr);
  else
    satellite_tree->SetBranchAddress("sigma_sat_value_syst",&satErr);
  satellite_tree->SetBranchAddress("Q_sat_time",&time);

  Double_t *satArray = new Double_t[ng];
  Double_t *satErrArray = new Double_t[ng];
  Double_t *timeArray = new Double_t[ng];
  for (Int_t i=0; i<ng; i++){
    satellite_tree->GetEntry(i);
    timeArray[i] = time;
    satArray[i] = sat;
    satErrArray[i] = satErr;
  }

  // fill graph
  TGraphErrors *sat_gr = new TGraphErrors(ng,timeArray,satArray,NULL,satErrArray);
  sat_gr->SetMarkerStyle(20);

  // define the limits for the plot
  // --> time
  Double_t time_min = 100e18;
  Double_t time_max = 0;
  for(Int_t i=0;i<ng;i++) {
    if(timeArray[i]<time_min) time_min=timeArray[i];
    if(timeArray[i]>time_max) time_max=timeArray[i];
  }
  time_min-=1e3;
  time_max+=1e3;

  Double_t sat_max = 0.5;
  Double_t sat_min = 0;
  if (Fill == 4937) sat_max = 0.3e-3;
  if (Fill == 6864) sat_max = 0.5e-3;
  if (Fill == 6012  ) {
    sat_max = 0.25e-3;
    sat_min = -0.1e-3;
  }
  if (Fill == 7483) {
    sat_max = 0.035;
    sat_min = 0.015;
  }
  // plot TGraph
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *sat_C = new TCanvas("sat_C","satellite charge",600,400);
  TH1F* frame = gPad->DrawFrame(time_min,sat_min,time_max,sat_max);
  TString device = "";
  if (opt == 0) device = "B1_A";
  if (opt == 1) device = "B1_B";
  if (opt == 2) device = "B2";
  frame->SetTitle(Form("Fraction of beam in Satellite Bucket (%s);timestamp ; satellite charge",device.Data()));
  sat_gr->Draw("p,e1,same");

  if (Fill == 4937 || Fill == 6012 || Fill == 6864){
    TF1 *fit1 = new TF1("fit1","pol0");
    sat_gr->Fit(fit1,"","",gScan0StartX,gScan0EndY);
    fit1->Draw("same");
     
    TF1 *fit2 = new TF1("fit2","pol0");
    fit2->SetLineColor(kBlue);
    sat_gr->Fit(fit2,"+","",gScan1StartX,gScan1EndY);
    fit2->Draw("same");

    TLatex* lat = new TLatex();
      lat->SetNDC();
      lat->SetTextSize(0.04);
      lat->SetTextFont(42);
      lat->SetTextColor(1);
    lat->DrawLatex(0.6,0.85,"Scan 0:");
    lat->DrawLatex(0.6,0.80,Form("#chi^{2}/ndf = %.2f/%d",fit1->GetChisquare(),fit1->GetNDF()));
    lat->DrawLatex(0.6,0.75,Form("p_{0} = %f #pm %f",fit1->GetParameter(0),fit1->GetParError(0)));
    lat->DrawLatex(0.6,0.70,"Scan 1:");
    lat->DrawLatex(0.6,0.65,Form("#chi^{2}/ndf = %.2f/%d",fit2->GetChisquare(),fit2->GetNDF()));
    lat->DrawLatex(0.6,0.60,Form("p_{0} = %f #pm %f",fit2->GetParameter(0),fit2->GetParError(0)));
    
  }

  //save plot
  if (save){
    TString unc = "syst";
    if (stat)
        unc = "stat";
	TString plotName = "";
	if (opt == 0) plotName = Form("../Fill-%d/ghostCharge/satelliteCharge_B1_A_%s.pdf",Fill,unc.Data());
    if (opt == 1) plotName = Form("../Fill-%d/ghostCharge/satelliteCharge_B1_B_%s.pdf",Fill,unc.Data());
    if (opt == 2) plotName = Form("../Fill-%d/ghostCharge/satelliteCharge_B2_A_%s.pdf",Fill,unc.Data());
	sat_C->SaveAs(plotName.Data());
  }


  // clean up
  delete [] satArray;
  delete [] satErrArray;
  delete [] timeArray;
}

void QA_satelliteCharge(Int_t Fill, Bool_t save = kFALSE) {
    
    // initialize
    Set_input_file_names(Fill);
    Set_pointers_to_input_files_and_trees();
    
    // find indices for start and end of scans
    cout << "Determining the start and end times of each scan, may take some time" << endl;
    Find_start_and_end_of_scans();
    SetScanTimes();


    //Statistical uncertainty
    //-----------------------
    //beam B1_A
    Get_satelliteCharge(Fill,0,kTRUE,save);
    //beam B1_B
    Get_satelliteCharge(Fill,1,kTRUE,save);
    //beam B2_A
    Get_satelliteCharge(Fill,2,kTRUE,save);
    
    //Systematic uncertainty
    //----------------------
    //beam B1_A
    Get_satelliteCharge(Fill,0,kFALSE,save);
    //beam B1_B
    Get_satelliteCharge(Fill,1,kFALSE,save);
    //beam B2_A
    Get_satelliteCharge(Fill,2,kFALSE,save);
}
