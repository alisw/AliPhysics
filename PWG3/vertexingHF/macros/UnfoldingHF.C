//-----------------------------------------------------//
//                                                     //
//      Base Macro to evaluate pt-eta unfolding        //
//      for charm mesons                               //
//                                                     //
// Usage:                                              //
//                                                     //
// 1) By default need output.root file from AliCF...   //
//    task                                             //
// 2) ./output.root - used to evaluate the eff         //
// 3) To correct, the data sample should be spitted in //
//    2, in this case in GetMeasuredSpectrum() connect //
//    the second container ./output_data.root          //
// 4) Set the number of iterations for Baysian unf     //
// 5) Set min. Chi2                                    //
// 6) Select if unfold at Acceptance level or          //
//    after PPR cuts                                   //
// 7) Gaussian error propagation assumed, may be       //
//    errors overestimated - to be studied             //
//                                                     //
//-----------------------------------------------------//
//                                                     //
//  a.grelli@uu.nl- Utrecht for R.Vernaut example      //
//                                                     //
//-----------------------------------------------------//
//                                                     //
//   NOTE: by default the macro is computing the       //
//   unfolding after acceptance step. To evaluate      //
//   after PPR switch on the option in the CF task!    //
//   and change in this macro the selection steps      //
//                                                     //
//-----------------------------------------------------//

#include "TString.h"
#include <iostream>
#include <fstream>

void UnfoldingHF() {

  // not jet impemented the "systematics" mode!

  //TString smode(analysis_mode);
  //TString hist(method);
  
  Load();

  printf("==================================================================\n");
  printf("===========        RUNNING D MESONS UNFOLDING           ==========\n");
  printf("==================================================================\n");

  // --------------------- get variables from container ----------------------//

  AliCFDataGrid  *measured         = (AliCFDataGrid*) GetMeasuredSpectrum();
  AliCFDataGrid  *reconstructed    = (AliCFDataGrid*) GetReconstructedSpectrum();
  AliCFDataGrid  *generated        = (AliCFDataGrid*) GetGeneratedSpectrum();
  AliCFEffGrid   *efficiency       = (AliCFEffGrid*)  GetEfficiency();
  
  // the response matrix

  THnSparse      *response  = (THnSparse*) GetResponseMatrix();
    
  // -------------------------------------------------------------------------//
  //                                                                          //
  //    The HF correction framework has 12 dimensions. This is bad for        //
  //    pt-eta unfolding so reduce the dimesions from 12 to 2.                //
  //                                                                          //
  //--------------------------------------------------------------------------//
  

  THnSparse* guessed = CreateGuessed(((AliCFGridSparse*)generated->GetData())->GetGrid()) ;
  THnSparse* data    = CreateGuessed(((AliCFGridSparse*)measured->GetData())->GetGrid());
  THnSparse* Reco    = CreateGuessed(((AliCFGridSparse*)reconstructed->GetData())->GetGrid()) ;
  
 
  // best guess, traslate the 12 dimentions container over 2 (eta and pt)
  
  Int_t dim[2];

  dim[0] = 0;
  dim[1] = 1;
 
  THnSparse* BestGuess      = guessed ->Projection(2,dim);
  THnSparse* DataSample     = data    ->Projection(2,dim);
  THnSparse* Reconstructed  = Reco    ->Projection(2,dim);
  THnSparse* SimpleEff      = Reco    ->Projection(2,dim);

  SimpleEff->Divide(Reconstructed,BestGuess,1.,1.);
 

  // here I do the unfolding and I set the number of iterations and the chi2

  AliCFUnfolding unfolding("unfolding","",2,response,SimpleEff,DataSample,BestGuess);
  unfolding.SetMaxNumberOfIterations(100);
  unfolding.SetMaxChi2PerDOF(0.00000000005);
  unfolding.Unfold();

  THnSparse* result = unfolding.GetUnfolded();
  

  TH2D* h_gen     =  generated   ->Project(1,0);
  TH2D* h_meas    =  DataSample  ->Projection(1,0);
  TH2D* h_reco    =  Reconstructed ->Projection(1,0);
  TH2D* h_guessed =  BestGuess     ->Projection(1,0);
  TH2D* h_eff     =  SimpleEff     ->Projection(1,0);
  TH2D* h_unf     =  result        ->Projection(1,0);

  // Apply simple efficiency correction 


  TH2D* simpleC = (TH2D*)h_meas->Clone();
  simpleC->Divide(h_eff);
  
  TCanvas * canvas3 = new TCanvas("Unfolded efficiency map","canvas 3",1000,700);

  canvas3->cd();

  TH2D* ratio = (TH2D*)h_unf->Clone();
  ratio->SetTitle("corrected using unfolding");
  ratio->Divide(h_unf,h_guessed,1,1);
 
  ratio->Draw("cont4z");

  TCanvas * canvas4 = new TCanvas("Simple Efficiency map","",1000,700);

  canvas4->cd();
  h_meas->SetTitle("simple correction");
  h_meas->Divide(h_eff);
  h_meas->Divide(h_guessed);
  h_meas->Draw("cont4z");  

  TCanvas * distribution2 = new TCanvas("dist2","",1000,700);

  distribution2->cd();

  TH1D* gen2 = generated->Project(1); // generated pt
  gen2->SetLineColor(1);
  gen2->SetTitle("D0 #eta distribution");
  gen2->GetXaxis()->SetTitle("eta");
  gen2->SetMarkerStyle(21);
  gen2->SetMarkerSize(1.);
  gen2->SetMarkerColor(1);
  gen2->Draw();

  TH1D* meas2 = measured->Project(1); // generated pt
  meas2->SetTitle("measurements");
  meas2->SetLineColor(2);
  meas2->SetMarkerStyle(22);
  meas2->SetMarkerSize(1.);
  meas2->SetMarkerColor(2);
  meas2->DrawCopy("same");

  TH1D* unf2 = result->Projection(1); // generated pt
  unf2->SetTitle("measurements");
  unf2->SetLineColor(3);
  unf2->SetMarkerStyle(23);
  unf2->SetMarkerSize(1.);
  unf2->SetMarkerColor(3);
  unf2->DrawCopy("SAME");

  corr2 = new TH1D();
  corr2->Sumw2();
  corr2 = simpleC->ProjectionY();
  corr2->SetTitle("measurements");
  corr2->SetLineColor(4);
  corr2->SetMarkerStyle(24);
  corr2->SetMarkerSize(1.);
  corr2->SetMarkerColor(4);
  corr2->Sumw2();
  corr2->DrawCopy("SAME");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->SetHeader("Distributions");
  leg->AddEntry(gen2,"Generated PYTHIA","l");
  leg->AddEntry(meas2,"Reconstructed","l");
  leg->AddEntry(corr2,"Corrected for eff","l");
  leg->AddEntry(unf2,"Unfolded","l");

  leg->Draw();

  TCanvas * distribution3 = new TCanvas("dist3","",1000,700);

  distribution3->cd();


  TH1D* gen3 = generated->Project(0); // generated pt
  Int_t Entries = gen3->GetEntries();
  gen3->SetLineColor(1);
  gen3->SetTitle("D0 pt distribution");
  gen3->GetXaxis()->SetTitle("pt (GeV/c)");
  gen3->SetMarkerStyle(21);
  gen3->SetMarkerSize(1.);
  gen3->SetMarkerColor(1);
  gen3->Draw();


  distribution3->cd();
  TH1D* meas3 = measured->Project(0); // generated pt
  meas3->SetTitle("measurements");
  meas3->SetLineColor(2);
  meas3->SetMarkerStyle(22);
  meas3->SetMarkerSize(1.);
  meas3->SetMarkerColor(2);
  meas3->DrawCopy("SAME");

  TH1D* unf3 = result->Projection(0); // 
  unf3->SetTitle("unfolded");
  unf3->SetLineColor(3);
  unf3->SetMarkerStyle(23);
  unf3->SetMarkerSize(1.);
  unf3->SetMarkerColor(3);
  unf3->Sumw2();
  unf3->DrawCopy("SAME");

  TH1D* corr3 =  simpleC->ProjectionX();
  corr3->SetTitle("corrected");
  corr3->SetLineColor(4);
  corr3->SetMarkerStyle(24);
  corr3->SetMarkerSize(1.);
  corr3->SetMarkerColor(4);
  corr3->Sumw2();
  corr3->DrawCopy("SAME");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg -> SetHeader("Distributions");
  leg -> AddEntry(gen3,"Generated PYTHIA","l");
  leg -> AddEntry(meas3,"Reconstructed","l");
  leg -> AddEntry(corr3,"Corrected for eff","l");
  leg -> AddEntry(unf3,"Unfolded","l");

  leg->Draw();

  return;
}

// ====================================================================


void Load(){
  gSystem->Load("libANALYSIS");
  gSystem->Load("libCORRFW");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
}

//___________________________________________________________-

TObject* GetMeasuredSpectrum() {
  TFile * f = TFile::Open("output.root","read");
  AliCFContainer* c = (AliCFContainer*)f->Get("container");
  AliCFDataGrid* data1 = new AliCFDataGrid("data1","",*c);
  data1->SetMeasured(3);
  return data1;
}
TObject* GetReconstructedSpectrum() {
  TFile * f = TFile::Open("output.root","read");
  AliCFContainer* c = (AliCFContainer*)f->Get("container");
  AliCFDataGrid* data2 = new AliCFDataGrid("data2","",*c);
  data2->SetMeasured(3);
 
  return data2;
}
TObject* GetGeneratedSpectrum() {
  TFile * f = TFile::Open("output.root","read");
  AliCFContainer* c = (AliCFContainer*)f->Get("container");
  AliCFDataGrid* data3 = new AliCFDataGrid("data3","",*c);
  data3->SetMeasured(1);
  return data3;
}
TObject* GetEfficiency() {
  TFile * f = TFile::Open("output.root","read");
  AliCFContainer* c = (AliCFContainer*)f->Get("container");
  AliCFEffGrid* eff = new AliCFEffGrid("eff","",*c);
  eff->CalculateEfficiency(3,1);
  return eff;
}
THnSparse* GetResponseMatrix() {
  TFile * f = TFile::Open("output.root","read");
  THnSparse* h = (THnSparse*)f->Get("correlation");
  return h;
}

THnSparse* GetResponseMatrixPPR() {
  TFile * f = TFile::Open("output.root","read");
  THnSparse* h = (THnSparse*)f->Get("correlationPPR");
  return h;
}

THnSparse* CreateGuessed(const THnSparse* h) {
  THnSparse* guessed = (THnSparse*) h->Clone();
  return guessed ;
}


