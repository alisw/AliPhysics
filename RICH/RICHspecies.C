#include "iostream.h"

void RICHspecies (Int_t pion=1,Int_t kaon=0, Int_t proton=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

//Create canvases and plot style

  TCanvas *c1 = 0;
  TCanvas *c2 = 0;
  /*TCanvas *c3 = 0;
  TCanvas *c4 = 0;
   TCanvas *c5 = 0;
   TCanvas *c6 = 0;
   TCanvas *c7 = 0;
   TCanvas *c8 = 0;
   TCanvas *c9 = 0;
   TCanvas *c10 = 0;
   TCanvas *c11 = 0;
   TCanvas *c12 = 0;
   TCanvas *c13 = 0;*/

   //TF1* expo = 0;
   //TF1* gaus = 0;
   
  TStyle *mystyle=new TStyle("Plain","mystyle");
  mystyle->SetPalette(1,0);
  //mystyle->SetTitleYSize(0.2);
   //mystyle->SetStatW(0.19);
  //mystyle->SetStatH(0.1);
  //mystyle->SetStatFontSize(0.01);
   //mystyle->SetTitleYSize(0.3);
  mystyle->SetFuncColor(2);
  mystyle->SetOptStat(0000);
  mystyle->SetDrawBorder(0);
  mystyle->SetTitleBorderSize(0);
  mystyle->SetOptFit(0000);
  mystyle->cd();


  //defining the histos

  TH1F *piondata = new TH1F("ckovangle1","Produced Cerenkov angle per photon",100,.35,1);
  TH1F *kaondata = new TH1F("ckovangle2","Produced Cerenkov angle per photon",100,.35,1);
  TH1F *protondata = new TH1F("ckovangle3","Produced Cerenkov angle per photon",100,.35,1);

  TH1F *pionrec = new TH1F("omega3D1","Reconstructed Cerenkov angle per photon",100,.35,1);
  TH1F *kaonrec = new TH1F("omega3D2","Reconstructed Cerenkov angle per photon",100,.35,1);
  TH1F *protonrec = new TH1F("omega3D3","Reconstructed Cerenkov angle per photon",100,.35,1);

  TH2F *pionid = new TH2F("identification1","Particle Identification",100,1,5,100,0,.8);
  TH2F *kaonid = new TH2F("identification2","Particle Identification",100,1,5,100,0,.8);
  TH2F *protonid = new TH2F("identification3","Particle Identification",100,1,5,100,0,.8);

// Connect the Root Galice files containing Geometry, Kine and Hits

  //pion data file

  if(pion)
    {
      TFile *filepion = (TFile*)gROOT->GetListOfFiles()->FindObject("pion.root");
      if (filepion) filepion->Close(); 
      filepion = new TFile("pion.root","UPDATE");

      
      piondata = (TH1F*) filepion->Get("ckovangle");
      if (piondata) printf("Pion data object found on filepion\n");
      if (!piondata) printf("Pion data object not found on filepion\n");
      //printf(" %f %f %f %d\n", data->GetMean(), data->GetMaximum(), data->GetRMS(),data->GetEntries() );
      pionrec = (TH1F*) filepion->Get("omega");
 
      pionid = (TH2F*) filepion->Get("identification");
    }


  //kaon data file

  if(kaon)
    {
      filekaon = (TFile*)gROOT->GetListOfFiles()->FindObject("kaon.root");
      if (filekaon) filekaon->Close(); 
      filekaon = new TFile("kaon.root","UPDATE");
      

      kaondata = (TH1F*) filekaon->Get("ckovangle");
      if (kaondata) printf("Kaon data object found on filekaon\n");
      if (!kaondata) printf("Kaon data object not found on filekaon\n");
      //printf(" %f %f %f %d\n", data->GetMean(), data->GetMaximum(), data->GetRMS(),data->GetEntries() );
      kaonrec = (TH1F*) filekaon->Get("omega");
      
      kaonid = (TH2F*) filekaon->Get("identification");
    }



  //proton data file
  
  if(proton)
    {
      TFile *fileproton = (TFile*)gROOT->GetListOfFiles()->FindObject("proton.root");
      if (fileproton) fileproton->Close(); 
      fileproton = new TFile("proton.root","UPDATE");
      
      
      TH1F *protondata = (TH1F*) fileproton->Get("ckovangle");
      if (protondata) printf("Proton data object found on fileproton\n");
      if (!protondata) printf("Proton data object not found on fileproton\n");
      //printf(" %f %f %f %d\n", data->GetMean(), data->GetMaximum(), data->GetRMS(),data->GetEntries() );
      TH1F *protonrec = (TH1F*) fileproton->Get("omega");
      
      protonid = (TH2F*) fileproton->Get("identification");
    }






  c1 = new TCanvas("c1","Cerenkov angle",50,50,300,700);
  c1->Divide(1,2);

  c1->cd(1);
  piondata->SetFillColor(5);
  piondata->SetXTitle("(rad)");
  piondata->Draw();

  kaondata->SetFillColor(4);
  kaondata->SetXTitle("(rad)");
  kaondata->Draw("same");

      
  protondata->SetFillColor(3);
  protondata->SetXTitle("(rad)");
  protondata->Draw("same");


  c1->cd(2);
  pionrec->SetFillColor(5);
  pionrec->SetXTitle("(rad)");
  pionrec->Draw();

  kaonrec->SetFillColor(4);
  kaonrec->SetXTitle("(rad)");
  kaonrec->Draw("same");

      
  protonrec->SetFillColor(3);
  protonrec->SetXTitle("(rad)");
  protonrec->Draw("same");

  c1 = new TCanvas("c12","Cerenkov angle vs. Momentum",150,150,550,350);

  
  TF1 *pionplot = new TF1("pion","acos(sqrt((.139*.139+x*x)/(x*x*1.285*1.285)))",1,5);
  TF1 *kaonplot = new TF1("kaon","acos(sqrt((.439*.439+x*x)/(x*x*1.285*1.285)))",1,5);
  TF1 *protonplot = new TF1("proton","acos(sqrt((.938*.938+x*x)/(x*x*1.285*1.285)))",1,5);
	

  pionplot->SetLineColor(5);
  pionplot->Draw("same");

  kaonplot->SetLineColor(4);
  kaonplot->Draw("same");
  
  protonplot->SetLineColor(3);
  protonplot->Draw("same");

  pionid->SetXTitle("Momentum (GeV/c)");
  pionid->SetYTitle("Cherenkov angle (radians)");
  
  pionid->Draw("cont0");
  
  kaonid->Draw("cont0 same");
  
  protonid->Draw("cont0 same");


  pionplot->SetLineColor(5);
  pionplot->Draw("same");

  kaonplot->SetLineColor(4);
  kaonplot->Draw("same");
  
  protonplot->SetLineColor(3);
  protonplot->Draw("same");


  pionplot->SetLineColor(5);
  pionplot->Draw("same");

  kaonplot->SetLineColor(4);
  kaonplot->Draw("same");
  
  protonplot->SetLineColor(3);
  protonplot->Draw("same");
	   

  //filepion->Close();
  
  //delete gAlice;
  printf("\nEnd of Macro  *************************************\n");
}



