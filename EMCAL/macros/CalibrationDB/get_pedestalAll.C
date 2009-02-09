typedef struct
{        
  int row;
  int col;
  int gain;
  int csp;
  int fee;
  int chip;
  int chan;
} mapData;

#include "constants.h"

void get_pedestalAll(int runNumber = 33, 
		     const int save = 0, const int saveDb=0)
{
  //Extracts pedestals from the file RUN_0{runNumber}.root.
  //If saveDb != 0, saves pedestals to the local OCDB, separately for low and high gains.
  //Used constant.h file with the constants and map.txt file with mapping related 
  //to the Sept-Oct 2007 beam test at CERN.

  TString DBFolder  ="local://";
  TString DBPathLow = "Calib/Data/LowGain";
  TString DBPathHigh = "Calib/Data/HighGain";
    
  Int_t firstRun = runNumber;
  Int_t lastRun = 999999999;

  // setup
  gROOT->Time();
  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  
  gStyle->SetFillColor(10);
  gStyle->SetCanvasColor(10);
  
  gStyle->SetPadBorderSize(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
 
  gStyle->SetTitleOffset(1.2,"X");
  gStyle->SetTitleOffset(0.9,"Y");                  
  gStyle->SetTitleSize(0.045,"X");
  gStyle->SetTitleSize(0.045,"Y");

  char FECstr0[20];
  char FECstr1[20];
  sprintf(FECstr0, "FEC %d", FEC[0]);
  sprintf(FECstr1, "FEC %d", FEC[1]);

  // get root file
  char rootFile[80];
  sprintf(rootFile,"RUN_%04d.root",runNumber);
  cout<<"\n root file: "<<rootFile<<endl;
  TFile *f = new TFile(rootFile);
  if (f==NULL) 
    {
      cout<<"\n ---> No root file.\n"<<endl;
      return 0;
    }


  // get map info
  ifstream fin;
  fin.open("map.txt");
    
  mapData chmap[CHANNELS];

  char dummy[20];
  fin >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
 
  int ichan = 0;
  int idum = 0;
  int nchan = 0;
  while ( fin.good() ) {
    fin >> idum >> ichan;
    if (! (fin.good() && ichan<CHANNELS) ) break;
    fin >> chmap[ichan].row >> chmap[ichan].col >> chmap[ichan].gain >> chmap[ichan].csp
	>> chmap[ichan].fee >> chmap[ichan].chip >> chmap[ichan].chan;
    nchan++;
  }

  cout << " info for " << nchan << " channels read from map " << endl;

  // get hSig and hPed profile histograms; the indices are channel id's
  TProfile *hPed = (TProfile *) gROOT->FindObject("hPed");
  TProfile *hMean = (TProfile *) gROOT->FindObject("hMean");
  TProfile *hSig = (TProfile *) gROOT->FindObject("hSig");
  TProfile *hAll = (TProfile *) gROOT->FindObject("hAll");
  hPed->Approximate(false); // prevent ROOT from smoothing. profiles.
  hMean->Approximate(false);
  hSig->Approximate(false);
  hAll->Approximate(false);

  hPed->SetErrorOption("s"); // don't divide error with sqrt(entries)
  hMean->SetErrorOption("s");
  hSig->SetErrorOption("s");
  hAll->SetErrorOption("s");

  TH1F *hRMS = new TH1F("hRMS", "hRMS", CHANNELS, -0.5, CHANNELS - 0.5);
  TH1F *hRMS1D = new TH1F("hRMS1D", "hRMS1D", 100, 0, 1);

  // define canvas : we'll have plots with info vs channel and info vs CSP (high and low)
  bool extra = false;
  if (extra) {
    TCanvas *cSignal = new TCanvas("cSignal","Signal",500,0,500,620);  
    cSignal->Divide(1,3);
  }
  TCanvas *cPedestal = new TCanvas("cPedestal","Pedestal",400,50,500,620);  
  cPedestal->Divide(1,3); 
  TCanvas *cMean = new TCanvas("cMean","Mean",400,50,500,620);  
  cMean->Divide(1,3); 
  TCanvas *cRMS = new TCanvas("cRMS","Pedestal RMS",400,50,500,620);  
  cRMS->Divide(1,3); 

  TCanvas *cRMS1D = new TCanvas("cRMS1D","Pedestal RMS - 1D",10,10,400,300);  

  const int NGAIN = 2;
  const char * GAINSTR[] = {"Low", "High"};

  // plot results
  TH1F* hRMSCSP[NGAIN];
  TH1F* hPedCSP[NGAIN]; 
  TH1F* hSigCSP[NGAIN]; 
  TH1F* hMeanCSP[NGAIN]; 

  char info[80];
  char id[80];

  for (int igain=0; igain<NGAIN; igain++) { 
    sprintf(id, "hRMSCSP%d",igain);
    sprintf(info, "RMS vs CSP # - %s gain",GAINSTR[igain]);
    hRMSCSP[igain] = new TH1F(id, info, NCSP, -0.5, NCSP-0.5);

    sprintf(id, "hPedCSP%d",igain);
    sprintf(info, "Pedestal vs CSP # - %s gain",GAINSTR[igain]);
    hPedCSP[igain] = new TH1F(id, info, NCSP, -0.5, NCSP-0.5);

    sprintf(id, "hSigCSP%d",igain);
    sprintf(info, "Signal vs CSP # - %s gain",GAINSTR[igain]);
    hSigCSP[igain] = new TH1F(id, info, NCSP, -0.5, NCSP-0.5);

    sprintf(id, "hMeanCSP%d",igain);
    sprintf(info, "Mean vs CSP # - %s gain",GAINSTR[igain]);
    hMeanCSP[igain] = new TH1F(id, info, NCSP, -0.5, NCSP-0.5);

    // setup the way the plots should look
    hRMSCSP[igain]->SetXTitle("CSP channels");
    hRMSCSP[igain]->SetYTitle("RMS (ADC counts)");
    hRMSCSP[igain]->SetMarkerStyle(20);
    hRMSCSP[igain]->SetMarkerColor(2);
    hRMSCSP[igain]->SetMarkerSize(1);
    hRMSCSP[igain]->GetYaxis()->SetNdivisions(505);

    hPedCSP[igain]->SetXTitle("CSP channels");
    hPedCSP[igain]->SetYTitle("Pedestal (ADC counts)");
    hPedCSP[igain]->SetMarkerStyle(20);
    hPedCSP[igain]->SetMarkerColor(2);
    hPedCSP[igain]->SetMarkerSize(1);
    hPedCSP[igain]->GetYaxis()->SetNdivisions(505);

    hSigCSP[igain]->SetXTitle("CSP channels");
    hSigCSP[igain]->SetYTitle("Signal (ADC counts)");
    hSigCSP[igain]->SetMarkerStyle(20);
    hSigCSP[igain]->SetMarkerColor(2);
    hSigCSP[igain]->SetMarkerSize(1);
    hSigCSP[igain]->GetYaxis()->SetNdivisions(505);

    hMeanCSP[igain]->SetXTitle("CSP channels");
    hMeanCSP[igain]->SetYTitle("Signal (ADC counts)");
    hMeanCSP[igain]->SetMarkerStyle(20);
    hMeanCSP[igain]->SetMarkerColor(2);
    hMeanCSP[igain]->SetMarkerSize(1);
    hMeanCSP[igain]->GetYaxis()->SetNdivisions(505);
  }

  //Creates calibration objects
  AliEMCALCalibData* cdb[2];
  cdb[0] = new AliEMCALCalibData(); //for low gain
  cdb[1] = new AliEMCALCalibData(); //for high gain

  // fill CSP histograms too
  for (Int_t k=0; k<CHANNELS; k++) {
    int igain = chmap[k].gain;
    int icsp = chmap[k].csp;
    int col = chmap[k].col;
    int row = chmap[k].row;
    if (chmap[k].fee != FEC[0]) { icsp += NCSP/2; } // offset for 2nd card

    // report suspicious channels
    if ( hPed->GetBinContent(k+1) == 0 ) {
      printf("Zero pedestal channel : %d, FEE %d Chip %d Chan %d : CSP %d gain %d\n", 
	     k, chmap[k].fee, chmap[k].chip, chmap[k].chan,  icsp, igain); 
    }

    hRMSCSP[igain]->SetBinContent(icsp+1, hPed->GetBinError(k+1));
    hRMS->SetBinContent(k+1, hPed->GetBinError(k+1));
    if (igain == 1) {
      hRMS1D->Fill(hPed->GetBinError(k+1));
    }

    Int_t supermod=0;
    cdb[igain]->SetADCpedestal(supermod,col,row,hPed->GetBinContent(k+1));
    printf("--- col%d row%d gain%d ped %.3f\n",col,row,igain,
	   cdb[igain]->GetADCpedestal(supermod,col,row));

    hPedCSP[igain]->SetBinContent(icsp+1, hPed->GetBinContent(k+1));
    hPedCSP[igain]->SetBinError(icsp+1, hPed->GetBinError(k+1));

    hSigCSP[igain]->SetBinContent(icsp+1, hSig->GetBinContent(k+1));
    hSigCSP[igain]->SetBinError(icsp+1, hSig->GetBinError(k+1));

    hMeanCSP[igain]->SetBinContent(icsp+1, hMean->GetBinContent(k+1));
    hMeanCSP[igain]->SetBinError(icsp+1, hMean->GetBinError(k+1));
  }

  //-------------------------------------------------------------------------------
  
  //Save calibration objects to OCDB.
  if(saveDb) {

    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    AliCDBManager::Instance()->SetSpecificStorage("EMCAL/*",DBFolder.Data());

    AliCDBMetaData md;
    md.SetComment("EMCAL beam test at CERN (Oct 2007)");
   
    AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("EMCAL/*");
    if(storage) {
 
      AliCDBId idLow(DBPathLow.Data(),firstRun,lastRun);
      storage->Put(cdb[0],idLow, &md);

      AliCDBId idHigh(DBPathHigh.Data(),firstRun,lastRun);
      storage->Put(cdb[1],idHigh, &md);

    }
    
  }
  
  
  // draw RMS
  cRMS->cd(1);
  hRMS->Draw("hist");
  cRMS->cd(2);
  hRMSCSP[0]->Draw("P");
  cRMS->cd(3);
  hRMSCSP[1]->Draw("P");

  // lines
  double min = hRMSCSP[1]->GetMinimum();
  double max = hRMSCSP[1]->GetMaximum();
  int ntcards = NCSP/8; // 8 CSPs per T-card
  for (int tcard = 0; tcard<ntcards; tcard++) {
    double xpos = tcard*8 - 0.5;
    TLine *line = new TLine(xpos, min, xpos, max); 
    line->SetLineColor(3);
    line->SetLineWidth(1);
    line->Draw("same");

    //text
    double xtex = tcard*8 + 1;
    double ytex = max * 0.1;
    char texstr[20];
    sprintf(texstr,"T-card %d", tcard%4);        
    TLatex *tex = new TLatex(xtex, ytex, texstr);
    tex->SetTextSize(0.05);
    tex->SetTextColor(4);
    tex->SetLineWidth(2);
    tex->Draw();
  }

  char Run[40];
  sprintf(Run,"Run: %d", runNumber);
  TLatex *tex = new TLatex(8, max*1.1, Run);    // 8 is just an arb. position    
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();
   
  TLatex *tex = new TLatex(16, max*1.1, ", Gain: high");        
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  cRMS->cd(1);
  min = hRMS->GetMinimum();
  max = hRMS->GetMaximum();
  double fecpos = CHANNELS/2 - 0.5;
  TLine *line = new TLine(fecpos, min, fecpos, max); 
  line->SetLineColor(4);
  line->SetLineWidth(1);
  line->Draw("same");

  TLatex *tex = new TLatex(fecpos*0.5, max*1.1, FECstr0);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex *tex = new TLatex(fecpos*1.5, max*1.1, FECstr1);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  cRMS->cd(2);
  min = hRMSCSP[0]->GetMinimum();
  max = hRMSCSP[0]->GetMaximum();
  fecpos = NCSP/2 - 0.5;
  TLine *line = new TLine(fecpos, min, fecpos, max); 
  line->SetLineColor(4);
  line->SetLineWidth(1);
  line->Draw("same");

  TLatex *tex = new TLatex(fecpos*0.5, max*1.1, FECstr0);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex *tex = new TLatex(fecpos*1.5, max*1.1, FECstr1);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  // make Pedestal plot
  //-------------------------------------------------------------------------------
  cPedestal->cd(1);
  hPed->Draw("hist");
  cPedestal->cd(2);
  hPedCSP[0]->Draw("P");
  cPedestal->cd(3);
  hPedCSP[1]->Draw("P");

  // lines
  min = hPedCSP[1]->GetMinimum();
  max = hPedCSP[1]->GetMaximum();
  for (int tcard = 0; tcard<ntcards; tcard++) {
    double xpos = tcard*8 - 0.5;
    TLine *line = new TLine(xpos, min, xpos, max); 
    line->SetLineColor(3);
    line->SetLineWidth(1);
    line->Draw("same");

    //text
    double xtex = tcard*8 + 1;
    double ytex = max * 0.1;
    char texstr[20];
    sprintf(texstr,"T-card %d", tcard%4);        
    TLatex *tex = new TLatex(xtex, ytex, texstr);
    tex->SetTextSize(0.05);
    tex->SetTextColor(4);
    tex->SetLineWidth(2);
    tex->Draw();
  }

  sprintf(Run,"Run: %d", runNumber);
  TLatex *tex = new TLatex(8, max*1.1, Run);    // 8 is just an arb. position    
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();
   
  TLatex *tex = new TLatex(16, max*1.1, ", Gain: high");        
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  cPedestal->cd(1);
  min = hPed->GetMinimum();
  max = hPed->GetMaximum();
  fecpos = CHANNELS/2 - 0.5;
  TLine *line = new TLine(fecpos, min, fecpos, max); 
  line->SetLineColor(4);
  line->SetLineWidth(1);
  line->Draw("same");

  TLatex *tex = new TLatex(fecpos*0.5, max*1.1, FECstr0);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex *tex = new TLatex(fecpos*1.5, max*1.1, FECstr1);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  cPedestal->cd(2);
  min = hPedCSP[0]->GetMinimum();
  max = hPedCSP[0]->GetMaximum();
  fecpos = NCSP/2 - 0.5;
  TLine *line = new TLine(fecpos, min, fecpos, max); 
  line->SetLineColor(4);
  line->SetLineWidth(1);
  line->Draw("same");

  TLatex *tex = new TLatex(fecpos*0.5, max*1.1, FECstr0);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex *tex = new TLatex(fecpos*1.5, max*1.1, FECstr1);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  // make Signal plot 
  //-------------------------------------------------------------------------------
  if (extra) {
  cSignal->cd(1);
  hSig->Draw("hist");
  cSignal->cd(2);
  hSigCSP[0]->Draw("P");
  cSignal->cd(3);
  hSigCSP[1]->Draw("P");

  // lines
  min = hSigCSP[1]->GetMinimum();
  max = hSigCSP[1]->GetMaximum();
  for (int tcard = 0; tcard<ntcards; tcard++) {
    double xpos = tcard*8 - 0.5;
    TLine *line = new TLine(xpos, min, xpos, max); 
    line->SetLineColor(3);
    line->SetLineWidth(1);
    line->Draw("same");

    //text
    double xtex = tcard*8 + 1;
    double ytex = max * 0.1;
    char texstr[20];
    sprintf(texstr,"T-card %d", tcard%4);        
    TLatex *tex = new TLatex(xtex, ytex, texstr);
    tex->SetTextSize(0.05);
    tex->SetTextColor(4);
    tex->SetLineWidth(2);
    tex->Draw();
  }

  sprintf(Run,"Run: %d", runNumber);
  TLatex *tex = new TLatex(8, max*1.1, Run);    // 8 is just an arb. position    
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();
   
  TLatex *tex = new TLatex(16, max*1.1, ", Gain: high");        
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  cSignal->cd(1);
  min = hSig->GetMinimum();
  max = hSig->GetMaximum();
  fecpos = CHANNELS/2 - 0.5;
  TLine *line = new TLine(fecpos, min, fecpos, max); 
  line->SetLineColor(4);
  line->SetLineWidth(1);
  line->Draw("same");

  TLatex *tex = new TLatex(fecpos*0.5, max*1.1, FECstr0);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex *tex = new TLatex(fecpos*1.5, max*1.1, FECstr1);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  cSignal->cd(2);
  min = hSigCSP[0]->GetMinimum();
  max = hSigCSP[0]->GetMaximum();
  fecpos = NCSP/2 - 0.5;
  TLine *line = new TLine(fecpos, min, fecpos, max); 
  line->SetLineColor(4);
  line->SetLineWidth(1);
  line->Draw("same");

  TLatex *tex = new TLatex(fecpos*0.5, max*1.1, FECstr0);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex *tex = new TLatex(fecpos*1.5, max*1.1, FECstr1);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();
  }

  // make Mean plot 
  //-------------------------------------------------------------------------------
  cMean->cd(1);
  hMean->Draw("hist");
  cMean->cd(2);
  hMeanCSP[0]->Draw("P");
  cMean->cd(3);
  hMeanCSP[1]->Draw("P");

  // lines
  min = hMeanCSP[1]->GetMinimum();
  max = hMeanCSP[1]->GetMaximum();
  for (int tcard = 0; tcard<ntcards; tcard++) {
    double xpos = tcard*8 - 0.5;
    TLine *line = new TLine(xpos, min, xpos, max); 
    line->SetLineColor(3);
    line->SetLineWidth(1);
    line->Draw("same");

    //text
    double xtex = tcard*8 + 1;
    double ytex = max * 0.1;
    char texstr[20];
    sprintf(texstr,"T-card %d", tcard%4);        
    TLatex *tex = new TLatex(xtex, ytex, texstr);
    tex->SetTextSize(0.05);
    tex->SetTextColor(4);
    tex->SetLineWidth(2);
    tex->Draw();
  }

  sprintf(Run,"Run: %d", runNumber);
  TLatex *tex = new TLatex(8, max*1.1, Run);    // 8 is just an arb. position    
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();
   
  TLatex *tex = new TLatex(16, max*1.1, ", Gain: high");        
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  cMean->cd(1);
  min = hMean->GetMinimum();
  max = hMean->GetMaximum();
  fecpos = CHANNELS/2 - 0.5;
  TLine *line = new TLine(fecpos, min, fecpos, max); 
  line->SetLineColor(4);
  line->SetLineWidth(1);
  line->Draw("same");

  TLatex *tex = new TLatex(fecpos*0.5, max*1.1, FECstr0);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex *tex = new TLatex(fecpos*1.5, max*1.1, FECstr1);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  cMean->cd(2);
  min = hMeanCSP[0]->GetMinimum();
  max = hMeanCSP[0]->GetMaximum();
  fecpos = NCSP/2 - 0.5;
  TLine *line = new TLine(fecpos, min, fecpos, max); 
  line->SetLineColor(4);
  line->SetLineWidth(1);
  line->Draw("same");

  TLatex *tex = new TLatex(fecpos*0.5, max*1.1, FECstr0);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex *tex = new TLatex(fecpos*1.5, max*1.1, FECstr1);
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();

  // make unique files 
  if (save) {
  char gifFile[80];
  if (extra) {
    sprintf(gifFile,"gif/Run%d_Signal.gif",runNumber);
    cSignal->SaveAs(gifFile);
  }
  sprintf(gifFile,"gif/Run%d_Pedestal.gif",runNumber);
  cPedestal->SaveAs(gifFile);
  sprintf(gifFile,"gif/Run%d_RMS.gif",runNumber);
  cRMS->SaveAs(gifFile);
  sprintf(gifFile,"gif/Run%d_Mean.gif",runNumber);
  cMean->SaveAs(gifFile);
  }

  cRMS1D->cd();
  hRMS1D->Fit("gaus");
  hRMS1D->Draw();

}

