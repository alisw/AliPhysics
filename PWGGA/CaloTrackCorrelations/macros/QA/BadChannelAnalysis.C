// This macro has been developed to find badcells candidates in EMCal based on cells amplitude distributions
//
// Input needed :
// can be either outputs QA from AliAnaCalorimeterQA task (BadChannelAnalysis() function)
// Or from merged output of AliAnalysisTaskCaloCellsQA (use BCAnalysis() function)
//
// Ouput:
// in Results.txt file and cell control plots for bad channel candidates will be created in BadCandidate.pdf
// 
// To classify between bad/warm/good (this  last case can also occur in bad candidates) 
// the user will have to check manually (by eye)control plots. In those in red the bad cell candidate energy distrib in black a reference choosen
// just to guide the eye
// 
//
// Author : Alexis Mas (SUBATECH), based on getCellsRunQA.C from Olga Driga (SUBATECH)

#if !defined(__CINT__) || defined(__MAKECINT__) 
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3D.h>
#include <TLine.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TStyle.h>
#include <TFileMerger.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TString.h>
#include <TGridCollection.h>
#include <TGridResult.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#endif
using namespace std;


void Draw(Int_t cell[], Int_t iBC, Int_t nBC, const Int_t cellref=151){
  //Allow to produce a pdf file with badcells candidates (red) compared to a refence cell (black)
  
  gROOT ->SetStyle("Plain");
  gStyle->SetOptStat(0); 
  gStyle->SetFillColor(kWhite);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetPalette(1);
  
  char out[120]; char title[100]; char name[100];
  TString slide = "Cells ";
  slide+=cell[iBC];
  slide+="-";
  slide+=cell[iBC+nBC-1];
  sprintf(out,"%d-%d.gif",cell[iBC],cell[iBC+nBC-1]); 
  
  TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude"); 
  TH1 *hCellref       = hCellAmplitude->ProjectionY("badcells",cellref+1,cellref+1);
  
  Int_t i;
  TCanvas *c1 = new TCanvas("badcells","badcells",1000,750) ;
  c1->SetLogy();
  if      (nBC > 6) c1->Divide(3,3) ;
  else if (nBC > 3) c1->Divide(3,2) ;
  else              c1->Divide(3,1);  
  // hCellref->Rebin(3);
  
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  for(i=0; i<nBC ; i++){
    sprintf(name,"Cell %d",cell[iBC+i]) ;
    TH1 *hCell = hCellAmplitude->ProjectionY(name,cell[iBC+i]+1,cell[iBC+i]+1);
    c1->cd(i%9 + 1) ;
    c1->cd(i%9 + 1)->SetLogy(); 
    sprintf(title,"Cell %d      Entries : %d",cell[iBC+i], (Int_t)hCell->GetEntries()) ;
    hCell->SetLineColor(2)  ; 
    // cout<<title<<endl ; 
    hCell->SetMaximum(1e5);
    // hCell->Rebin(3);
    hCell->SetAxisRange(0.,4.);
    hCell->GetXaxis()->SetTitle("E (GeV)"); 
    hCell->GetYaxis()->SetTitle("N Entries");
    hCellref->SetAxisRange(0.,4.);
    hCell->SetLineWidth(1) ;
    hCellref->SetLineWidth(1) ;
    hCell->SetTitle(title);
    hCellref->SetLineColor(1)  ;  
    if(i==0){
      leg->AddEntry(hCellref,"reference","l"); 
      leg->AddEntry(hCell,"current","l");; 
    }
    hCell->Draw() ;
    hCellref->Draw("same") ; 
    leg->Draw();
  }
  
  //CREATE A PDF FILE 
  TString PdfName = "BadCellsCandidate.pdf";
  if(nBC<9) {
    PdfName +=")"; 
    c1->Print(PdfName.Data());}
  else if(iBC==0) {  
    PdfName +="("; 
    c1->Print(PdfName.Data());}
  else  c1->Print(PdfName.Data());
  
  
  //  c1->SaveAs(out);  
  delete hCellref ; 
  delete c1 ; 
  delete leg ;
  
}

//_________________________________________________________________________
//_________________________________________________________________________

void Convert(TString fCalorimeter = "EMCAL", TString period = "LHC11h", TString pass = "pass1_HLT", TString trigger= "default"){ 
  
  // Create one file for the analysis from several outputs QA files listed in runlist.txt
  // You need :
  //   runlist.txt with runs listed 
  //   outputsQA  e.g  period/pass/123456.root 
  
  TH2F * hCellAmplitude = new TH2F("hCellAmplitude","Cell Amplitude",11520,0,11520,200,0,10);
  TH1D * hNEventsProcessedPerRun = new TH1D("hNEventsProcessedPerRun","Number of processed events vs run number",200000,100000,300000);
  FILE * pFile;
  TString file = "/scratch/alicehp2/mas/analyse/QA/"+period+"/"+ pass + "/runlistMB.txt" ;
  cout<<file<<endl;
  pFile = fopen(file.Data(), "r"); //open the text file where include the run list and correct run index
  
  cout<<file<<endl;
  cout << " fcalo: " << fCalorimeter << "; period: " << period << "; pass: " << pass << "  trigger "<<trigger<<  endl; 
  
  Int_t ix,iy;
  Int_t Nentr;
  Int_t p;
  Int_t q;
  Int_t ncols;
  Int_t nlines = 0 ;
  Int_t RunId[500] ;
  
  Double_t x[500] ;
  Double_t xrun[500] ;
  
  while (1){
    ncols = fscanf(pFile,"%d  %d ",&p,&q);
    if (ncols< 0) break;
    x[nlines]=p;
    RunId[nlines]=q;
    xrun[nlines]=1.*q;
    nlines++;
  }
  
  fclose(pFile);
  
  const Int_t nRun = nlines ;
  Double_t content;
  TString base ;
  TString BCfile ;
  TString direct = "CaloQA_"; 
  direct += trigger;
  
  for(Int_t i = 0 ; i < nRun ; i++) { 
    base = "/scratch/alicehp2/germain/QA/";
    BCfile = base + period ;
    BCfile += trigger ;
    BCfile += ".root";
    base += period ;
    base += "/";
    base += pass ;
    base += "/";
    base += RunId[i] ;
    TString infile ;
    infile = base + ".root" ;
    TFile *f = TFile::Open(infile);
    base += "/" ;
    base += trigger ; 
    TDirectoryFile *dir = (TDirectoryFile *)f->Get(direct);
    TList *outputList = (TList*)dir->Get(direct);
    TH2F *hAmpId;
    TH2F *hNEvents;
    hAmpId =(TH2F *)outputList->FindObject("EMCAL_hAmpId");
    hNEvents =(TH2F *)outputList->FindObject("hNEvents");
    Nentr =  (Int_t)hNEvents->GetEntries();
    if (Nentr<100) continue ;
    hNEventsProcessedPerRun->SetBinContent(RunId[i]-100000,(Double_t)Nentr);
    
    for(ix=1;ix<=200;ix++){ 
      for(iy=1; iy<=11520; iy++){ 
        content = 0.0 ;
        content +=  hAmpId->GetBinContent(ix,iy);
        content +=  hCellAmplitude->GetBinContent(iy,ix);
        if(content > 0.5) 
          hCellAmplitude->SetBinContent(iy,ix,content);
      }
    }  
    
    //cout<<i<<endl;
    if(i==0){ cout<<"Merging/Converting procedure ..." ; cout.flush();}
    else { cout<<"..." ; cout.flush();}
    outputList->Delete();
    dir->Delete();
    f->Close();
    delete f;
  }
  
  
  TFile *BCF = TFile::Open(BCfile,"recreate");
  hNEventsProcessedPerRun->Write();
  hCellAmplitude->Write();
  BCF->Close();
  cout<<"DONE !"<<endl;
  
}

//_________________________________________________________________________
//_________________________________________________________________________

void Process(Int_t *pflag[11520][7], TH1* inhisto, Double_t Nsigma = 4., Int_t dnbins = 200, Double_t dmaxval = -1.)
{  
  //  1) create a distribution for the input histogram;
  //  2) fit the distribution with a gaussian
  //  3) define good area within +-Nsigma to identfy badcells
  // 
  // inhisto -- input histogram;
  // dnbins -- number of bins in distribution;
  // dmaxval -- maximum value on distribution histogram.
  
  Int_t crit = *pflag[0][0] ; //identify the criterum processed
  Double_t goodmax= 0. ; 
  Double_t goodmin= 0. ;
  *pflag[0][0] =1;
  if (dmaxval < 0.) {
    dmaxval = inhisto->GetMaximum()*1.01;  // 1.01 - to see the last bin
    if(crit==2 && dmaxval > 1) dmaxval =1. ;
  }    
  
  TH1 *distrib = new TH1F(Form("%sDistr",inhisto->GetName()), "", dnbins, inhisto->GetMinimum(), dmaxval);
  distrib->SetXTitle(inhisto->GetYaxis()->GetTitle());
  distrib->SetYTitle("Entries");
  
  // fill distribution
  for (Int_t c = 1; c <= inhisto->GetNbinsX(); c++)
    distrib->Fill(inhisto->GetBinContent(c));
  
  // draw histogram + distribution
  TCanvas *c1 = new TCanvas(inhisto->GetName(),inhisto->GetName(), 800,400);
  c1->Divide(2,1);
  
  c1->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetLogy();
  inhisto->SetTitleOffset(1.7,"Y");
  inhisto->Draw();
  
  c1->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.10);
  gPad->SetLogy();
  distrib->Draw();
  
  Int_t higherbin=0,i;
  for (i = 2; i <= distrib->GetNbinsX(); i++){
    if(distrib->GetBinContent(higherbin) < distrib->GetBinContent(i))
      higherbin = i ;
  }
  
  for(i=higherbin ; i<=dnbins ; i++)
    if(distrib->GetBinContent(i)<2) break ;
  goodmax = distrib->GetBinCenter(i);
  
  for(i=higherbin ; i>0 ; i--)
    if(distrib->GetBinContent(i)<2) break ;
  goodmin = distrib->GetBinLowEdge(i);
  
  TF1 *fit2 = new TF1("fit2", "gaus");
  
  distrib->Fit(fit2, "0LQEM", "", goodmin, goodmax);
  Double_t sig, mean, chi2ndf; 
  mean = fit2->GetParameter(1);
  sig = fit2->GetParameter(2);
  chi2ndf = fit2->GetChisquare()/fit2->GetNDF();
  goodmin = mean - Nsigma*sig ;
  goodmax = mean + Nsigma*sig ;
  
  // lines
  TLine *lline = new TLine(goodmin, 0, goodmin, distrib->GetMaximum());
  lline->SetLineColor(kOrange);
  lline->SetLineStyle(7);
  lline->Draw();
  
  TLine *rline = new TLine(goodmax, 0, goodmax, distrib->GetMaximum());
  rline->SetLineColor(kOrange);
  rline->SetLineStyle(7);
  rline->Draw();
  
  // legend
  TLegend *leg = new TLegend(0.60,0.82,0.9,0.88);
  leg->AddEntry(lline, "Good region boundary","l");
  leg->Draw("same");
  fit2->Draw("same");
  
  c1->Update();
  TString name = "criteria-" ;
  name+= crit;
  name+= ".gif";
  
  c1->SaveAs(name); 
  
  Int_t ntot = 0, cel;
  
  for (Int_t c = 1; c <= inhisto->GetNbinsX(); c++) {
    cel=(Int_t)(inhisto->GetBinLowEdge(c)+0.1);
    if (inhisto->GetBinContent(c) < goodmin) { 
      ntot++;
      *pflag[cel][crit]=0; 
    }
    else if (inhisto->GetBinContent(c) > goodmax) { 
      ntot++; 
      *pflag[cel][crit]=2;
    }
  }
}

//_________________________________________________________________________
//_________________________________________________________________________

void TestCellEandN(Int_t *pflag[11520][7], Double_t Emin = 0.1, Double_t Nsigma = 4., char* hname = "hCellAmplitude", Int_t dnbins = 200)
{
  
  
  // Three more tests for bad cells:
  //  1) total deposited energy;
  //  2) total number of entries;
  //  3) average energy = [total deposited energy]/[total number of entries].
  //
  
  
  // input; X axis -- absId numbers
  TH2 *hCellAmplitude = (TH2*) gFile->Get(hname);
  
  // binning parameters
  Int_t ncells = hCellAmplitude->GetNbinsX();
  Double_t amin = hCellAmplitude->GetXaxis()->GetXmin();
  Double_t amax = hCellAmplitude->GetXaxis()->GetXmax();
  
  TH1* hCellEtotal = new TH1F(Form("%s_hCellEtotal_E%.2f",hname,Emin),
                              Form("Total deposited energy, E > %.2f GeV",Emin), ncells,amin,amax);
  hCellEtotal->SetXTitle("AbsId");
  hCellEtotal->SetYTitle("Energy, GeV");
  
  TH1F *hCellNtotal = new TH1F(Form("%s_hCellNtotal_E%.2f",hname,Emin),
                               Form("Number of entries per events, E > %.2f GeV",Emin), ncells,amin,amax);
  hCellNtotal->SetXTitle("AbsId");
  hCellNtotal->SetYTitle("Entries");
  
  TH1F *hCellEtoNtotal = new TH1F(Form("%s_hCellEtoNtotal_E%.2f",hname,Emin),
                                  Form("Average energy per hit, E > %.2f GeV",Emin), ncells,amin,amax);
  hCellEtoNtotal->SetXTitle("AbsId");
  hCellEtoNtotal->SetYTitle("Energy, GeV");
  
  TH1* hNEventsProcessedPerRun = (TH1*) gFile->Get("hNEventsProcessedPerRun");
  Double_t totalevents = hNEventsProcessedPerRun->Integral(1, hNEventsProcessedPerRun->GetNbinsX());
  
  // fill cells
  for (Int_t c = 1; c <= ncells; c++) {
    Double_t Esum = 0;
    Double_t Nsum = 0;
    
    
    for (Int_t j = 1; j <= hCellAmplitude->GetNbinsY(); j++) {
      Double_t E = hCellAmplitude->GetYaxis()->GetBinCenter(j);
      Double_t N = hCellAmplitude->GetBinContent(c, j);
      if (E < Emin) continue;
      Esum += E*N;
      Nsum += N;
    }
    
    hCellEtotal->SetBinContent(c, Esum);
    hCellNtotal->SetBinContent(c, Nsum/totalevents);
    
    if (Nsum > 0.5)  // number of entries >= 1
      hCellEtoNtotal->SetBinContent(c, Esum/Nsum);
    
  }
  
  delete hCellAmplitude;
  
  // Process(hCellEtotal,   dnbins );
  if(*pflag[0][0]==1)
    Process(pflag, hCellEtoNtotal, Nsigma, dnbins );
  if(*pflag[0][0]==2)
    Process(pflag, hCellNtotal, Nsigma,  dnbins );
}

//_________________________________________________________________________
//_________________________________________________________________________

void TestCellShapes(Int_t *pflag[11520][7], Double_t fitEmin, Double_t fitEmax, Double_t Nsigma =4., char* hname = "hCellAmplitude", Int_t dnbins = 1000)
{
  // Test cells shape using fit function f(x)=A*exp(-B*x)/x^2.
  // Produce values per cell + distributions for A,B and chi2/ndf parameters.
  
  TH2 *hCellAmplitude = (TH2*) gFile->Get(Form("%s",hname));
  
  // binning parameters
  Int_t  ncells = hCellAmplitude->GetNbinsX();
  Double_t amin = hCellAmplitude->GetXaxis()->GetXmin();
  Double_t amax = hCellAmplitude->GetXaxis()->GetXmax();
  
  // initialize histograms
  TH1 *hFitA = new TH1F(Form("hFitA_%s",hname),"Fit A value", ncells,amin,amax);
  hFitA->SetXTitle("AbsId");
  hFitA->SetYTitle("A");
  
  TH1 *hFitB = new TH1F(Form("hFitB_%s",hname),"Fit B value", ncells,amin,amax);
  hFitB->SetXTitle("AbsId");
  hFitB->SetYTitle("B");
  
  TH1 *hFitChi2Ndf = new TH1F(Form("hFitChi2Ndf_%s",hname),"Fit #chi^{2}/ndf value", ncells,amin,amax);
  hFitChi2Ndf->SetXTitle("AbsId");
  hFitChi2Ndf->SetYTitle("#chi^{2}/ndf");
  
  Double_t maxval1=0., maxval2=0., maxval3=0.;
  Double_t prev=0., MSA=0., AvA = 0. ; //those param are used to automaticaly determined a reasonable maxval1
  Double_t prev2=0., MSB=0., AvB = 0.  ; //those param are used to automaticaly determined a reasonable maxval2
  Double_t prev3=0., MSki2=0., Avki2 = 0. ; //those param are used to automaticaly determined a reasonable maxval3
  Double_t ki2=0.0 ;
  
  for (Int_t k = 1; k <= ncells; k++) { 
    TF1 *fit = new TF1("fit", "[0]*exp(-[1]*x)/x^2");
    TH1 *hCell = hCellAmplitude->ProjectionY("",k,k);
    if (hCell->GetEntries() == 0) continue;
    // hCell->Rebin(3);
    hCell->Fit(fit, "0QEM", "", fitEmin, fitEmax);
    delete hCell; 
    
    if(fit->GetParameter(0) < 5000.){ 
      hFitA->SetBinContent(k, fit->GetParameter(0));
      if(k<3000) {
        AvA +=  fit->GetParameter(0);
        if(k==2999)  maxval1  = AvA/3000. ;
        if (prev < fit->GetParameter(0)) MSA += fit->GetParameter(0) - prev;
        else MSA -= (fit->GetParameter(0) - prev) ;
        prev = fit->GetParameter(0);
      }
      
      else 
      {
        
        if((fit->GetParameter(0) - maxval1) > 0. && (fit->GetParameter(0) - maxval1) < (MSA/1000.))
        {
          maxval1 = fit->GetParameter(0); 
        }
      }
    }
    else hFitA->SetBinContent(k, 5000.);
    
    
    
    if(fit->GetParameter(1) < 5000.){ 
      hFitB->SetBinContent(k, fit->GetParameter(1));
      if(k<3000) {
        AvB +=  fit->GetParameter(1);
        if(k==2999)  maxval2  = AvB/3000. ;
        if (prev2 < fit->GetParameter(1)) MSB += fit->GetParameter(1) - prev2;
        else MSB -= (fit->GetParameter(1) - prev2) ;
        prev2 = fit->GetParameter(1);
      }
      
      else 
      {
        
        if((fit->GetParameter(1) - maxval2) > 0. && (fit->GetParameter(1) - maxval2) < (MSB/1000.))
        {
          maxval2 = fit->GetParameter(1); 
        }
      }
    }
    else hFitB->SetBinContent(k, 5000.);
    
    
    if (fit->GetNDF() != 0 ) ki2 =  fit->GetChisquare()/fit->GetNDF();
    else ki2 = 1000.;
    
    if(ki2 < 1000.){ 
      hFitChi2Ndf->SetBinContent(k, ki2);
      if(k<3000) {
        Avki2 +=  ki2;
        if(k==2999)  maxval3  = Avki2/3000. ;
        if (prev3 < ki2) MSki2 += ki2 - prev3;
        else MSki2 -= (ki2 - prev3) ;
        prev3 = ki2;
      }
      
      else 
      {
        
        if((ki2 - maxval3) > 0. && (ki2 - maxval3) < (MSki2/1000.))
        {
          maxval3 = ki2; 
        }
      }
    }
    else hFitChi2Ndf->SetBinContent(k, 1000.);
    
    
    delete fit ;
  }
  
  delete hCellAmplitude;
  
  // if you have problem with automatic parameter :
  //  maxval1 = 
  //  maxval2 =
  //  maxval3 =
  
  
  if(*pflag[0][0]==3)
    Process(pflag, hFitChi2Ndf, Nsigma, dnbins, maxval3); 
  
  
  if(*pflag[0][0]==4)
    Process(pflag, hFitA, Nsigma, dnbins,  maxval1); 
  
  
  if(*pflag[0][0]==5)
    Process(pflag, hFitB, Nsigma, dnbins, maxval2);
  
}

//_________________________________________________________________________
//_________________________________________________________________________

void ExcludeCells(Int_t *pexclu[11520]) {
  //find the cell with 0 entrie for excluding
  TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude"); 
  
  
  for (Int_t c = 1; c <= 11520; c++) {
    Double_t Nsum = 0;
    
    for (Int_t l = 1; l <= hCellAmplitude->GetNbinsY(); l++) {
      Double_t N = hCellAmplitude->GetBinContent(c, l);
      Nsum += N;
    }
    if(Nsum < 0.5 && *pexclu[c-1]!=1) *pexclu[c-1]=1; //trick for criterum 7
    //if(Nsum < 0.5 ) *pexclu[c-1]=1; 
    else *pexclu[c-1]=0;
  }
  delete hCellAmplitude;
}

//_________________________________________________________________________
//_________________________________________________________________________

void KillCells(Int_t filter[], Int_t nbc) {
  // kill a cell : put it to 0 entrie 
  TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");
  
  for(Int_t i =0; i<nbc; i++){
    for(Int_t j=0; j<= hCellAmplitude->GetNbinsY() ;j++){
      hCellAmplitude->SetBinContent(filter[i]+1,j,0) ; }}
  
  TH1* hNEventsProcessedPerRun = (TH1*) gFile->Get("hNEventsProcessedPerRun");
  
  TFile *tf = new TFile("filter.root","recreate");
  hCellAmplitude->Write(); 
  hNEventsProcessedPerRun->Write(); 
  tf->Write();
  tf->Close();
  delete hCellAmplitude; delete hNEventsProcessedPerRun;
}

//_________________________________________________________________________
//_________________________________________________________________________

void PeriodAnalysis(Int_t criterum=7, Double_t Nsigma = 4.0, Double_t Emin=0.1, Double_t Emax=1.0, TString file ="none") { 
  
  // what it does in function of criterum value 
  
  // 1 : average E for E>Emin
  // 2 : entries for E>Emin
  // 3 : ki²/ndf  (from fit of each cell Amplitude between Emin and Emax) 
  // 4 : A parameter (from fit of each cell Amplitude between Emin and Emax) 
  // 5 : B parameter (from fit of each cell Amplitude between Emin and Emax) 
  // 6 : 
  // 7 : give bad + dead list
  
  Int_t newBC[11520]; // newBC[0] donne l'id de la premiere BC trouvée
  Int_t *pexclu[11520] ;
  Int_t exclu[11520];
  Int_t *pflag[11520][7] ;
  Int_t flag[11520][7];
  Int_t bad[1000] ;
  Int_t i, j, nb=0; 
  
  //INIT
  TString output, bilan;
  if(criterum == 7) bilan = "Results.txt" ;
  output.Form("Criterum-%d_Emin-%.2f.txt",criterum,Emin); 
  for(i=0;i<11520;i++) { exclu[i]=0; pexclu[i] =&exclu[i];
    for(j=0;j<7;j++) { flag[i][j] =1 ; pflag[i][j] = &flag[i][j];}}
  flag[0][0]=criterum ; //to identify the criterum tested
  
  
  //CELLS EXCLUDED
  ExcludeCells(pexclu); //exclude cells from analysis (will not appear in results)
  if(criterum < 7){
    cout<<"Excluded/dead cells : "<<endl;
    for(i=0;i<11520;i++) {if(exclu[i]!=0) {cout<<i<<", " ; nb++;}}
    cout<<"("<<nb<<")"<<endl; nb=0;}
  
  
  //CRITERUM 7 : FINAL RESULT
  if(criterum ==7) { 
    cout<<"FINAL RESULTS"<<endl;
    ofstream fichier(bilan, ios::out | ios::trunc);  
    if(fichier){
      fichier<<"Dead cells : "<<endl;  
      cout<<"Dead cells : "<<endl;
      for(i=0;i<11520;i++) {
        if(exclu[i]!=0) {fichier<<i<<", " ; cout<<i<<", " ; nb++;}}
      fichier<<"("<<nb<<")"<<endl; cout<<"("<<nb<<")"<<endl; nb=0;
      
      TFile::Open("filter.root");
      ExcludeCells(pexclu); 
      fichier<<"Bad cells candidates : "<<endl; cout<<"Bad cells candidates : "<<endl;
      for(i=0;i<11520;i++) {
        if(exclu[i]!=0) {bad[nb]=i; fichier<<i<<", " ; cout<<i<<", " ;
          nb++; if(nb==999){ cout<<"TO MUCH BAD CELLS"<<endl ; break;}}}
      fichier<<"("<<nb<<")"<<endl; cout<<"("<<nb<<")"<<endl;}
    fichier.close();
    
    if(file!="none"){
      TFile::Open(file);
      Int_t w=0 ;
      Int_t c;   
      for(w=0; (w*9)<=nb; w++) {
        if(9<=(nb-w*9)) c = 9 ; 
        else c = nb-9*w ;
        Draw(bad, w*9, c) ;
      }}
    
    nb=0;
  }
  
  
  //ANALYSIS
  if (criterum < 3)  TestCellEandN(pflag, Emin, Nsigma);
  else if (criterum < 6)
    TestCellShapes(pflag, Emin, Emax, Nsigma);
  
  
  //RESULTS
  if(criterum < 6) { nb=0;
    cout<<"bad by lower value : "<<endl;
    for(i=0;i<11520;i++) {
      if(flag[i][criterum]==0 && exclu[i]==0){nb++;
        cout<<i<<", " ;}} cout<<"("<<nb<<")"<<endl; nb=0;
    
    cout<<"bad by higher value : "<<endl;
    for(i=0;i<11520;i++) {
      if(flag[i][criterum]==2 && exclu[i]==0) {nb++;
        cout<<i<<", " ;}} cout<<"("<<nb<<")"<<endl; nb=0;
    
    cout<<"total bad "<<endl;
    for(i=0;i<11520;i++) {
      if(flag[i][criterum]!=1 && exclu[i]==0) {
        newBC[nb]=i;
        nb++;
        cout<<i<<", " ; }} cout<<"("<<nb<<")"<<endl;
    
    
    //create a filtered file
    KillCells(newBC,nb) ; nb=0;
    
    //write in a file the results
    ofstream fichier(output, ios::out | ios::trunc);  
    if(fichier)
    {
      fichier <<"criterum : "<<criterum<<", Emin = "<<Emin<<" GeV"<<", Emax = "<<Emax<<" GeV"<<endl;
      fichier<<"bad by lower value : "<<endl;
      for(i=0;i<11520;i++) {
        if(flag[i][criterum]==0 && exclu[i]==0){nb++;
          fichier<<i<<", " ;}} fichier<<"("<<nb<<")"<<endl; nb=0;
      
      fichier<<"bad by higher value : "<<endl;
      for(i=0;i<11520;i++) {
        if(flag[i][criterum]==2 && exclu[i]==0) {nb++;
          fichier<<i<<", " ;}} fichier<<"("<<nb<<")"<<endl; nb=0;
      
      fichier<<"total bad "<<endl;
      for(i=0;i<11520;i++) {
        if(flag[i][criterum]!=1 && exclu[i]==0) {
          newBC[nb]=i;
          nb++;
          fichier<<i<<", " ; }} fichier<<"("<<nb<<")"<<endl;   
      fichier.close();  
    }
    else  
      cerr << "opening error" << endl; 
    
  }
  
}


//_________________________________________________________________________
//_________________________________________________________________________

void BCAnalysis(TString file, TString trigger = "default"){
  
  //Configure a complete analysis with different criteria, it provides bad+dead cells lists
  //You can manage criteria used and their order, the first criteria will use the original 
  //output file from AliAnalysisTaskCaloCellsQA task, then after each criteria it will use a filtered file without the badchannel previously identified
  
  if(trigger=="default"){
    
    TFile::Open(file);
    PeriodAnalysis(3, 8., 0.1, 2.5);
    TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 0.1, 2.5); 
    TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 0.5, 2.5);
    TFile::Open("filter.root");
    PeriodAnalysis(1, 4., 0.1, 2.5);
    TFile::Open("filter.root");
    PeriodAnalysis(4, 4., 0.1,2.5); 
    
  }
  
  else { //you have the possibility to change analysis configuration  in function of trigger type
    
    TFile::Open(file);
    PeriodAnalysis(3, 8., 0.1, 2.5);
    TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 0.1, 2.5); 
    TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 0.5, 2.5);
    TFile::Open("filter.root");
    PeriodAnalysis(1, 4., 0.1, 2.5);
    TFile::Open("filter.root");
    PeriodAnalysis(4, 4., 0.1, 2.5); 
    
  }
  
  TFile::Open(file);
  PeriodAnalysis(7,0.,0.,0.,file); //provide dead cells list from original file and draw bad cells candidate from indicated file
  
}

//_________________________________________________________________________
//_________________________________________________________________________

void BadChannelAnalysis(TString fCalorimeter = "EMCAL", TString period = "LHC11h", TString pass = "pass1_HLT", TString trigger= "default"){
  Convert(fCalorimeter, period, pass, trigger);
  TString inputfile =  "/scratch/alicehp2/germain/QA/" + period;
  inputfile += trigger;
  inputfile += ".root";
  BCAnalysis(inputfile, trigger);
}
