/// \file BadChannelAnalysis.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Do EMCal bad channel analysis
///
/// This macro has been developed to find badcells candidates 
/// in EMCal based on cells amplitude distributions
/// Input needed can be either outputs QA from AliAnaCalorimeterQA task 
/// (BadChannelAnalysis() function)
/// Or from merged output of AliAnalysisTaskCaloCellsQA (use BCAnalysis() function)
///
///-----------------
/// Main method:
///---------------
/// BadChannelAnalysis:
///
///    step 1 : Convert()
///       read list of mergeable runs  in your working directory (in example below the $workdir is "/scratch/alicehp2/germain/QANew2/"
///       The QAresults.root files should be aleady copied from alien and be in $workdir/period/pass/runnb.root
///       read/merge the histos"EMCAL_hAmpId" and "EMCAL_hTimeId"  from QAresults.root file and write them in
///       $workdir/period/pass/period_pass_Runlist0New.root"  !!! this is hardcoded !!!!!
///       step 1 has to be called only the first time runing on a new list
///
///    step 2 BCanalysis() main method to analyse previously created file (hardcoded)
///
///       call of different Periodanalysis(criterium,..) functions according to the wanted tests (critreria)
///           1 : average E for E>Emin
///           2 : entries for E>Emin
///           3 : ki²/ndf  (from fit of each cell Amplitude between Emin and Emax)
///           4 : A parameter (from fit of each cell Amplitude between Emin and Emax)
///           5 : B parameter (from fit of each cell Amplitude between Emin and Emax)
///           6 :
///           7 : give bad + dead list
///
///       Mainly used: 1 and 2 (with different settings (chi2, intervals of energy : this is quite dependent of the stat you may have )
///       Further improvement: implement tests 1 and 2 on time distribustion histogram
///
/// ---------------------
///  Running the macro
/// ---------------------
///  root [2] .L BadChannelAnalysis.C++
///  root [2] BadChannelAnalysis("EMCAL","LHC15o","muon_caloLego","AnyINTnoBC",trial=0)
///
///  !!! pay attention the trigger name depends on the caloQA_triggername you want to analyse check first in QAresults.root what is abvailable
///  !!! it is generally not good to run it on triggered data for the following reasons:
///
/// ------------------
/// outputs
///-------------------
/// the output of this analysis povides you:
/// intermediate steps files: (those will be recreated each time you rerun so pa attention to save the different files when changing period/listof runs....
///   - Criterum-xx_Emin-xx_Emax-xx.txt : list of identified bad for the different test/Emin/Emax
///  - period_pass.txt file with list of dead/bad cells identified
///  - a pdf file with all energy distributions plots of all bad cells candidates (compared to a reference one (hard coded see Draw function to change)
///
///    Further improvement: implement tests 1 and 2 on time distribution histogram
///
///
/// \author Alexis Mas, SUBATECH
/// \author Marie Germain,  SUBATECH
/// based on getCellsRunQA.C from 
/// \author Olga Driga (SUBATECH)
///

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
#include <alloca.h>
#include <string>
#include <cstring>
#endif
using namespace std;

void Draw2(Int_t cell, Int_t cellref=400) {
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0); 
  gStyle->SetFillColor(kWhite);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetPalette(1);
  char out[120]; char title[100]; char name[100];char name2[100];
  TString slide(Form("Cells %d-%d",cell,cell));

  sprintf(out,"%d.gif",cell);
  TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");
  TH1 *hCellref = hCellAmplitude->ProjectionX("badcells",cellref+1,cellref+1);

  TCanvas *c1 = new TCanvas("badcells","badcells",600,600) ;
  c1->SetLogy();

  // hCellref->Rebin(3);
   TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

    sprintf(name,"Cell %d",cell) ;
    TH1 *hCell = hCellAmplitude->ProjectionX(name,cell+1,cell+1);

    sprintf(title,"Cell %d      Entries : %d Enties ref: %d",cell, (Int_t)hCell->GetEntries(),(Int_t)hCellref->GetEntries()) ;
    hCell->SetLineColor(2)  ;
    // cout<<title<<endl ;
    hCell->SetMaximum(1e5);
    // hCell->Rebin(3);
    hCell->SetAxisRange(0.,10.);
    hCell->GetXaxis()->SetTitle("E (GeV)");
    hCell->GetYaxis()->SetTitle("N Entries");
    hCellref->SetAxisRange(0.,10.);
    hCell->SetLineWidth(1) ;
    hCellref->SetLineWidth(1) ;
    hCell->SetTitle(title);
    hCellref->SetLineColor(1)  ;
    leg->AddEntry(hCellref,"reference","l");
   leg->AddEntry(hCell,"current","l");
    hCell->Draw() ;
    hCellref->Draw("same") ;
    leg->Draw();
    sprintf(name2,"Cell%dLHC13MB.gif",cell) ;
    c1->SaveAs(name2);
  
}

void Draw(Int_t cell[], Int_t iBC, Int_t nBC,TString datapath="/scratch/alicehp2/germain/QANew2/", TString period="LHC15f", TString pass="pass2", Int_t trial=0,const Int_t cellref=2377){
  //Allow to produce a pdf file with badcells candidates (red) compared to a refence cell (black)

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetFillColor(kWhite);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetPalette(1);
  char out[120]; char title[100]; char name[100];
  TString slide(Form("Cells %d-%d",cell[iBC],cell[iBC+nBC-1]));

  TString reflegend =  Form("reference Cell %i",cellref);
  sprintf(out,"%d-%d.gif",cell[iBC],cell[iBC+nBC-1]); 
  TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude"); 
  TH1 *hCellref = hCellAmplitude->ProjectionX("badcells",cellref+1,cellref+1);
  Int_t i;
  TCanvas *c1 = new TCanvas("badcells","badcells",1000,750) ;
  c1->SetLogy();
  if(nBC > 6) c1->Divide(3,3) ;
  else if (nBC > 3)  c1->Divide(3,2) ;
  else  c1->Divide(3,1);
  // hCellref->Rebin(3);
   TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  for(i=0; i<nBC ; i++){
    sprintf(name,"Cell %d",cell[iBC+i]) ;
   TH1 *hCell = hCellAmplitude->ProjectionX(name,cell[iBC+i]+1,cell[iBC+i]+1);

   c1->cd(i%9 + 1) ;
    c1->cd(i%9 + 1)->SetLogy(); 
    sprintf(title,"Cell %d      Entries : %d  Ref : %d",cell[iBC+i], (Int_t)hCell->GetEntries(), (Int_t)hCellref->GetEntries() ) ;
    hCell->SetLineColor(2)  ; 
    // cout<<title<<endl ; 
    hCell->SetMaximum(1e6);
    // hCell->Rebin(3);
    hCell->SetAxisRange(0.,10.);
    hCell->GetXaxis()->SetTitle("E (GeV)"); 
    hCell->GetYaxis()->SetTitle("N Entries");
    hCellref->SetAxisRange(0.,8.);
    hCell->SetLineWidth(1) ;
    hCellref->SetLineWidth(1) ;
    hCell->SetTitle(title);
    hCellref->SetLineColor(1)  ;  



    if(i==0){

      leg->AddEntry(hCellref,reflegend,"l");
      leg->AddEntry(hCell,"current","l");
    }
    hCell->Draw() ;
    hCellref->Draw("same") ; 
    leg->Draw();
  }

  //CREATE A PDF FILE 
  TString PdfName =  Form("%s/BadChannelNew/2016/%s%sList1Test%i.pdf",datapath.Data(),period.Data(),pass.Data(),trial);

  //  TString PdfName =  Form("/scratch/alicehp2/germain/QANew2/BadChannelNew/2016/%s%sList0Test%i.pdf",period.Data(),pass.Data(),trial);
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

void Convert(TString datapath= "/scratch/alicehp2/germain/QANew2",TString fCalorimeter = "EMCAL", TString period = "LHC11h", TString pass = "pass1_HLT", TString trigger= "default"){

  //Create one file for the analysis from several outputs QA files listed in runlist.txt
  //You need :
  // runlist.txt with runs listed
  // outputsQA  e.g  period/pass/123456.root

  TH2F *hCellAmplitude = new TH2F("hCellAmplitude","Cell Amplitude",200,0,10,23040,0,23040);
  TH2F *hCellTime = new TH2F("hCellTime","Cell Time",250,-275,975,23040,0,23040);

  TH1D *hNEventsProcessedPerRun = new TH1D("hNEventsProcessedPerRun","Number of processed events vs run number",200000,100000,300000);
  FILE * pFile;

  TString file = Form("%s/%s%sBC1.txt",datapath.Data(),period.Data(),pass.Data());

  //  TString file = Form("/scratch/alicehp2/germain/QANew2/%s%sBC0.txt",period.Data(),pass.Data());
cout << " file = " << file << endl;;
  pFile = fopen(file.Data(), "r"); //open the text file where include the run list and correct run index


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
    ncols = fscanf(pFile,"  %d ",&q);
    if (ncols< 0) break;
    //    x[nlines]=p;
    RunId[nlines]=q;
    //xrun[nlines]=1.*q;
    nlines++;
  }
  fclose(pFile);
  const Int_t nRun = nlines ;
  Double_t content;
  TString base ;
  TString base2 ;
  TString BCfile ;


 TString direct(Form("CaloQA_%s",trigger.Data()));
 // TString direct  ="EMCAL_TrigAnyINT_Cl" ;

  for(Int_t i = 0 ; i < nRun ; i++) { 
    // base2 = Form("/scratch/alicehp2/germain/QANew2/%s/%s/",period.Data(),pass.Data());
    // base = Form("/scratch/alicehp2/germain/QANew2/%s/%s/%d",period.Data(),pass.Data(),RunId[i]);
    base2 = Form("%s/%s/%s/",datapath.Data(),period.Data(),pass.Data());
    base = Form("%s/%s/%s/%d",datapath.Data(),period.Data(),pass.Data(),RunId[i]);
    BCfile=Form("%s%s%sRunlist1New.root",base2.Data(),period.Data(),pass.Data());

    TString infile ;

    // ICI on met le nom period/pass/runblabla.root
    if ((pass=="cpass1_pass2")||(pass=="cpass1-2")){
      if (trigger=="default"){
	infile = Form("%s_barrel.root",base.Data());}
      else {infile = Form("%s_outer.root",base.Data());}
    }
    else
      infile = Form("%s.root",base.Data()) ;
    cout<<"file : "<<infile<<endl;
    TFile *f = TFile::Open(infile);
    //   base += "/" ;
    // base += trigger ;
    base=Form("%s/%s",base.Data(),trigger.Data());
    if (!f) continue;
    cout << " jusqu'ici ca va "<< endl;
    TDirectoryFile *dir = (TDirectoryFile *)f->Get(direct);
      if (!dir) continue;
     cout << " jusqu'ici ca va dir OK "<< endl;
   TList *outputList = (TList*)dir->Get(direct);
    if (!outputList) continue;
    //  TList *outputList = (TList *)f->Get(direct.Data());
      cout << " jusqu'ici ca va List OK "<< endl;

    TH2F *hAmpId;
    TH2F *hNEvents;
    TH2F * hTimeId;

    hAmpId =(TH2F *)outputList->FindObject("EMCAL_hAmpId");
    hAmpId->Draw();

   hTimeId =(TH2F *)outputList->FindObject("EMCAL_hTimeId");
   hTimeId->Draw();

   cout<<"file : "<<infile<<endl;
    hNEvents =(TH2F *)outputList->FindObject("hNEvents");
    if (!hNEvents)continue;
    Nentr =  (Int_t)hNEvents->GetEntries();
    cout << " N entries " << Nentr << endl;
    if (Nentr<100) continue ;
    hNEventsProcessedPerRun->SetBinContent(RunId[i]-100000,(Double_t)Nentr);



      hCellAmplitude->Add(hAmpId);
      hCellTime->Add(hTimeId);


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
  hCellTime->Write();
  BCF->Close();
  cout<<"DONE !"<<endl;

}

//_________________________________________________________________________
//_________________________________________________________________________

void Process(Int_t *pflag[23040][7], TH1* inhisto, Double_t Nsigma = 4., Int_t dnbins = 200, Double_t dmaxval = -1., Int_t compteur = 1)
{
  //  1) create a distribution for the input histogram;
  //  2) fit the distribution with a gaussian
  //  3) define good area within +-Nsigma to identfy badcells
  //
  // inhisto -- input histogram;
  // dnbins -- number of bins in distribution;
  // dmaxval -- maximum value on distribution histogram.

       gStyle->SetOptStat(1); // MG modif
       gStyle->SetOptFit(1); // MG modif
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
  // Marie midif to take into account very non gaussian distrig
  mean = fit2->GetParameter(1);
  sig = fit2->GetParameter(2);
  chi2ndf = fit2->GetChisquare()/fit2->GetNDF();


//   mean = distrib->GetMean(); // MG
//   sig = distrib->GetRMS(); // MG

 // cout << "----------------------------------------------"<< endl;
//   cout <<" pass " << compteur <<  " mean " << mean << " rms" << sig << endl;
//   cout << "----------------------------------------------"<< endl;
  if (mean <0.) mean=0.;

    goodmin = mean - Nsigma*sig ;
    goodmax = mean + Nsigma*sig ;

    if (goodmin<0) goodmin=0.;
 //    if (compteur==0){
//       goodmin = 0.;
//       goodmax = 2.*mean;
//     }

  cout << "-----------------------------------------------------------------"<< endl;
  cout  << " pass " << compteur <<  " mean " << mean << " rms" << sig << " goodmin " << goodmin<< " goodmax" << goodmax <<  endl;
  cout << "-----------------------------------------------------------------"<< endl;

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
  //  fit2->Draw("same");

  c1->Update();
  TString name = "criteria-" ;
  name+= crit;
  name+= ".gif";
  
  c1->SaveAs(name); 

  Int_t ntot = 0, cel;

  for (Int_t c = 1; c <= inhisto->GetNbinsX(); c++) {
    cel=(Int_t)(inhisto->GetBinLowEdge(c)+0.1);
    if (inhisto->GetBinContent(c) <= goodmin) {
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

void TestCellEandN(Int_t *pflag[23040][7], Double_t Emin = 0.1, Double_t Emax=2., Double_t Nsigma = 4.,  Int_t compteur = 1, char const * hname = "hCellAmplitude", Int_t dnbins = 200)
{


  // Three more tests for bad  cells:
  //  1) total deposited energy;
  //  2) total number of entries;
  //  3) average energy = [total deposited energy]/[total number of entries].
  //

 //  Int_t count;
//   count = compteur;

  // input; X axis -- absId numbers

   TH2 *hCellAmplitude = (TH2*) gFile->Get(hname);

  // binning parameters
  Int_t ncells = hCellAmplitude->GetNbinsY();
  Double_t amin = hCellAmplitude->GetYaxis()->GetXmin();
  Double_t amax = hCellAmplitude->GetYaxis()->GetXmax();

  cout << "ncells " << ncells << " amin = " << amin << "amax = " << amax<< endl;


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
    

    for (Int_t j = 1; j <= hCellAmplitude->GetNbinsX(); j++) {
    //   for (Int_t j =  hCellAmplitude->GetXaxis()->FindBin(Emin); j <= hCellAmplitude->GetXaxis()->FindBin(Emax); j++) {
      Double_t E = hCellAmplitude->GetXaxis()->GetBinCenter(j);
      Double_t N = hCellAmplitude->GetBinContent(j, c);
      if (E < Emin || E>Emax) continue;
      // if (E > Emin && E< Emax) {
      Esum += E*N;
      Nsum += N;
      //}
    }

    hCellEtotal->SetBinContent(c, Esum);
    hCellNtotal->SetBinContent(c, Nsum/totalevents);
    //  hCellNtotal->SetBinContent(c, Nsum);

    if (Nsum > 0.)  // number of entries >= 1
      hCellEtoNtotal->SetBinContent(c, Esum/Nsum);
    
  }

  delete hCellAmplitude;

  // Process(hCellEtotal,   dnbins );
  if(*pflag[0][0]==1)
    Process(pflag, hCellEtoNtotal, Nsigma, dnbins, -1,compteur);
  if(*pflag[0][0]==2)
    Process(pflag, hCellNtotal, Nsigma,  dnbins,-1, compteur);
}

//_________________________________________________________________________
//_________________________________________________________________________

void TestCellShapes(Int_t *pflag[23040][7], Double_t fitEmin, Double_t fitEmax, Double_t Nsigma =4., Int_t compteur= 1, char const * hname = "hCellAmplitude", Int_t dnbins = 1000)
{
  // Test cells shape using fit function f(x)=A*exp(-B*x)/x^2.
  // Produce values per cell + distributions for A,B and chi2/ndf parameters.
//   Int_t count;
//   count = compteur;

  TH2 *hCellAmplitude = (TH2*) gFile->Get(Form("%s",hname));

  // binning parameters
  Int_t  ncells = hCellAmplitude->GetNbinsY();
  Double_t amin = hCellAmplitude->GetYaxis()->GetXmin();
  Double_t amax = hCellAmplitude->GetYaxis()->GetXmax();
  cout << "ncells " << ncells << " amin = " << amin << "amax = " << amax<< endl;

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
    TH1 *hCell = hCellAmplitude->ProjectionX("",k,k);
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
    Process(pflag, hFitChi2Ndf, Nsigma, dnbins, maxval3,compteur);

  
  if(*pflag[0][0]==4)
    Process(pflag, hFitA, Nsigma, dnbins,  maxval1,compteur);


  if(*pflag[0][0]==5)
    Process(pflag, hFitB, Nsigma, dnbins, maxval2,compteur);

}

//_________________________________________________________________________
//_________________________________________________________________________

void ExcludeCells(Int_t *pexclu[23040]) {
  //find the cell with 0 entrie for excluding
  TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");


  for (Int_t c = 1; c <= 23040; c++) {
    Double_t Nsum = 0;
    //   cout << "exclude cells ca va "<<endl
    for (Int_t l = 1; l <= hCellAmplitude->GetNbinsX(); l++) {
      Double_t N = hCellAmplitude->GetBinContent(l, c);
      Nsum += N;
    }
    //  cout << "2emem exclude cells ca va "<< c-1 << endl;
    if(Nsum < 0.5 && *pexclu[c-1]!=1)
      { *pexclu[c-1]=1;    }//trick for criterum 7
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
    for(Int_t j=0; j<= hCellAmplitude->GetNbinsX() ;j++){
      hCellAmplitude->SetBinContent(j,filter[i]+1,0) ; }}

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

void PeriodAnalysis(Int_t criterum=7, Double_t Nsigma = 4.0, Double_t Emin=0.1, Double_t Emax=2.0, Int_t compteur = 1, TString datapath="/scratch/alicehp2/germain/QANew2",TString period = "LHC15f", TString pass = "pass2", Int_t trial=0, TString file ="none"){
  
  // what it does in function of criterum value

  // 1 : average E for E>Emin
  // 2 : entries for E>Emin
  // 3 : ki²/ndf  (from fit of each cell Amplitude between Emin and Emax)
  // 4 : A parameter (from fit of each cell Amplitude between Emin and Emax)
  // 5 : B parameter (from fit of each cell Amplitude between Emin and Emax)
  // 6 :
  // 7 : give bad + dead list

  Int_t newBC[23040]; // newBC[0] donne l'id de la premiere BC trouvée
  Int_t *pexclu[23040] ;
  Int_t exclu[23040];
  Int_t *pflag[23040][7] ;
  Int_t flag[23040][7];
  Int_t bad[10000] ;
  Int_t i, j, nb=0;

  //INIT
  TString output, bilan;
  //bilan = Form("%s%sBC0Test%s.txt",period.Data());
  if(criterum == 7) bilan = Form("%s%sBC0Test%i.txt",period.Data(),pass.Data(),trial);
  // if(criterum == 7) bilan = "ResultsLHC15oBCRunList0Test1.txt" ;
  output.Form("Criterum-%d_Emin-%.2f_Emax-%.2f.txt",criterum,Emin,Emax);
  for(i=0;i<23040;i++) { exclu[i]=0; pexclu[i] =&exclu[i];
    for(j=0;j<7;j++) { flag[i][j] =1 ; pflag[i][j] = &flag[i][j];}}
  flag[0][0]=criterum ; //to identify the criterum tested
  

  //CELLS EXCLUDED
   ExcludeCells(pexclu); //exclude cells from analysis (will not appear in results)
  // if(criterum < 7){
  //   cout<<"Excluded/dead cells : "<<endl;
  //   for(i=0;i<23040;i++) {if(exclu[i]!=0) {cout<<i<<", " ; nb++;}}
  //   cout<<"("<<nb<<")"<<endl; nb=0;}
  

  //CRITERUM 7 : FINAL RESULT
  if(criterum ==7) { 
    cout<<"FINAL RESULTS"<<endl;
    ofstream fichier(bilan, ios::out | ios::trunc);  
    if(fichier){
      fichier<<"Dead cells : "<<endl;  
      cout<<"Dead cells : "<<endl;
      // MG modif
      //for(i=0;i<23040;i++) {
      for(i=0;i<17665;i++) {
	//if(exclu[i]!=0) {fichier<<i<<", " ; cout<<i<<", " ; nb++;}}
	    if(exclu[i]!=0) {fichier<<i<<"\n" ; cout<<i<<", " ; nb++;}}
  // fichier<<i<<"\n" ; cout<<i<<", " ; nb++;}
      fichier<<"("<<nb<<")"<<endl; cout<<"("<<nb<<")"<<endl; nb=0;

      TFile::Open("filter.root");
      ExcludeCells(pexclu); 
      fichier<<"Bad cells : "<<endl; cout<<"Bad cells : "<<endl;
      for(i=0;i<17665;i++) {
	//	if(exclu[i]!=0) {bad[nb]=i; fichier<<i<<", " ; cout<<i<<", " ;
	if(exclu[i]!=0) {bad[nb]=i; fichier<<i<<"\n" ; cout<<i<<", " ;
	  nb++;
	  //if(nb==999){ cout<<"TO MUCH BAD CELLS"<<endl ; break;}
	}
      }
      fichier<<"("<<nb<<")"<<endl; cout<<"("<<nb<<")"<<endl;}
    fichier.close();

    if(file!="none"){
      TFile::Open(file);
      Int_t w=0 ;
      Int_t c;   
      for(w=0; (w*9)<=nb; w++) {
	if(9<=(nb-w*9)) c = 9 ;
	else c = nb-9*w ;
	Draw(bad, w*9, c,datapath,period,pass,trial) ;
      }}

    nb=0;
  }

  
  //ANALYSIS
  if (criterum < 3)  TestCellEandN(pflag, Emin, Emax,Nsigma,compteur);
  else if (criterum < 6)
    TestCellShapes(pflag, Emin, Emax, Nsigma,compteur);
  
  
  //RESULTS
  if(criterum < 6) { nb=0;
    cout<<"bad by lower value  Emin : "<< Emin << "  Emax = " << Emax << endl;
    for(i=0;i<17665;i++) {
      if(flag[i][criterum]==0 && exclu[i]==0){nb++;
	cout<<i<<", " ;}} cout<<"("<<nb<<")"<<endl; nb=0;

    cout<<"bad by higher value  Emin : "<< Emin << "  Emax = " << Emax  <<endl;
    for(i=0;i<17665;i++) {
      if(flag[i][criterum]==2 && exclu[i]==0) {nb++;
	cout<<i<<", " ;}} cout<<"("<<nb<<")"<<endl; nb=0;

    cout<<"total bad "<<endl;
    for(i=0;i<17665;i++) {
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
	for(i=0;i<17665;i++) {
	  if(flag[i][criterum]==0 && exclu[i]==0){nb++;
	    fichier<<i<<", " ;}} fichier<<"("<<nb<<")"<<endl; nb=0;

	fichier<<"bad by higher value : "<<endl;
	for(i=0;i<17665;i++) {
	  if(flag[i][criterum]==2 && exclu[i]==0) {nb++;
	    fichier<<i<<", " ;}} fichier<<"("<<nb<<")"<<endl; nb=0;
	fichier<<"total bad "<<endl;
	for(i=0;i<17665;i++) {
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

  void BCAnalysis(TString file, TString datapath="scratch/alicehp2/germain/QANew2",TString trigger = "default",TString period = "LHC15f", TString pass = "pass2",Int_t trial = 0){

  //Configure a complete analysis with different criteria, it provides bad+dead cells lists
  //You can manage criteria used and their order, the first criteria will use the original output file from AliAnalysisTaskCaloCellsQA task, then after each criteria it will use a filtered file without the badchannel previously identified



    Int_t criter;
    Double_t  Emini, Emaxi, Nsig;




// Default Configuration:

    //if(trigger=="default"){
 if(trigger=="default"||trigger=="INT7"||trigger=="DMC7"||trigger=="AnyINTnoBC"){


    TFile::Open(file);
    PeriodAnalysis(2, 4., 0.2, 0.5,1,datapath,period,pass,trial); // nb ent emin emax
    TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 0.5, 1.,1,datapath,period,pass,trial); // nb ent emin emax
    TFile::Open("filter.root");
    PeriodAnalysis(1, 6., 0.5, 1.,1,datapath,period,pass,trial); // energy mea emin emax
    TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 1., 2.,1,datapath,period,pass,trial); // nb ent emin emax
    TFile::Open("filter.root");
    PeriodAnalysis(1, 6., 1., 2.,1,datapath,period,pass,trial); // energy mea emin emax
    TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 1., 10.,1,datapath, period,pass,trial); //nb ent emin emax
    //    TFile::Open("filter.root");
    // PeriodAnalysis(1, 6., 1., 10.,1,datapath,period,pass,trial); //energy mea emin emax

   TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 2., 3.,1,datapath,period,pass,trial); //entriesmea emin emax
   TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 3., 4.,1,datapath,period,pass,trial); //entries  mea emin emax
  TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 4., 5.,1,datapath,period,pass,trial); //entriesmea emin emax
   TFile::Open("filter.root");
    PeriodAnalysis(2, 4., 5., 10.,1,datapath,period,pass,trial); //entries  mea emin emax

//     TFile::Open("filter.root");
//     PeriodAnalysis(4, 4., 0.1, 2.,1,period,pass,trial);  // param A fit E = A exp(-Bx)/x^2

  }


  else { //you have the possibility to change analysis configuration  in function of trigger type

  TFile::Open(file);
  PeriodAnalysis(2, 6., 0.5, 2.,1,datapath,period,pass,trial); // nb ent emin emax
    TFile::Open("filter.root");
    PeriodAnalysis(1, 6., 0.5, 2.,1,datapath,period,pass,trial); // energy mea emin emax
    TFile::Open("filter.root");
    PeriodAnalysis(2, 6., 2., 5.,1,datapath,period,pass,trial); // nb ent emin emax
    TFile::Open("filter.root");
    PeriodAnalysis(1, 6., 2., 5.,1,datapath,period,pass,trial); // energy mea emin emax
   TFile::Open("filter.root");
   PeriodAnalysis(2, 6., 5., 10.,1,datapath,period,pass,trial); // nb ent emin emax
    TFile::Open("filter.root");
    PeriodAnalysis(1, 6., 5., 10.,1,datapath,period,pass,trial); // energy mea emin emax
  TFile::Open("filter.root");
    PeriodAnalysis(1, 6., 5., 10.,1,datapath,period,pass,trial); // energy mea emin emax

    // this was old settings for EMC riggers checks
//     TFile::Open(file);
// //    PeriodAnalysis(3, 8., 1., 3.,1,period,pass,trial);


  }

  TFile::Open(file);
  PeriodAnalysis(7,0.,0.,0.,1,datapath,period,pass,trial,file); //provide dead cells list from original file and draw bad cells candidate from indicated file

}

//_________________________________________________________________________
//________________________________________________________________________

void BadChannelAnalysis(TString datapath= "/scratch/alicehp2,germain/QANew2",TString fCalorimeter = "EMCAL", TString period = "LHC15f", TString pass = "pass2", TString trigger= "default",Int_t trial=0){
  //Convert(datapath,fCalorimeter, period, pass, trigger);
  //   TString inputfile(Form( "/scratch/alicehp2/germain/QANew2/%s/%s/%s%sRunlist0New.root",period.Data(),pass.Data(),period.Data(),pass.Data(),trigger.Data()));

     TString inputfile(Form( "%s/%s/%s/%s%sRunlist1New.root",datapath.Data(),period.Data(),pass.Data(),period.Data(),pass.Data(),trigger.Data()));

     BCAnalysis(inputfile,datapath,trigger,period,pass,trial);
}
