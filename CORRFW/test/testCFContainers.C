#include <Riostream.h>

extern TRandom *gRandom;
extern TBenchmark *gBenchmark;
extern TSystem *gSystem;

void testCFContainers(){

  // simple example macros for usage of a N-dim container (AliCFContainer)
  // handling a set of grids to accumulate data at different 
  // selection steps, & derive efficiency (this is stored in AliCFEffGrid) 
  // book, fill and draw some histos
  // The efficiency is then used to correct the data (trivially self-correct, 
  // in this example)

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(111110);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);

  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include  -I$ROOTSYS/include");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("$ALICE_ROOT/CORRFW/libCORRFW.so") ;
 
  //Setting up the container grid... 

  const Int_t nstep=2; //number of selection steps  (just 2 in this ex)

  const Int_t nvar=4; //number of variables on the grid:pt,vtx

  const Int_t nbin1=6; //bins in pt
  const Int_t nbin2=10; //bins in eta 
  const Int_t nbin3=18; //bins in phi
  const Int_t nbin4=10; //bins in vertex


  //Flag the sel steps. In this example, we have two, may be any nstep
  Int_t stepGen=0;
  Int_t stepRec=1;

  //the sensitive variables, their indeces
  Int_t ipt =0;
  Int_t ieta=1;
  Int_t iphi=2;
  Int_t ivtx=3;

  //arrays for the number of bins in each dimension
  const Int_t iBin[nvar] ={nbin1,nbin2,nbin3,nbin4};

  //arrays for bin limits
  Double_t binLim1[nbin1+1];
  Double_t binLim2[nbin2+1];
  Double_t binLim3[nbin3+1];
  Double_t binLim4[nbin4+1];
  
  for(Int_t i=0;i<=nbin1;i++){
    // pt [0,3] GeV/c
    binLim1[i]=i*0.5; 
  }

  for(Int_t i=0;i<=nbin2;i++){
    // eta [-1,1] 
    binLim2[i]=i*0.2-1.;
  }
  for(Int_t i=0;i<=nbin3;i++){
    //phi [0,360]
    binLim3[i]=i*20.;
  }
  for(Int_t i=0;i<=nbin4;i++){
    //vertex [-20,20] cm
    binLim4[i]=i*4.-20.;
  }
  
  
  //the nstep grids "container" 
  //   AliCFContainer *cont = new AliCFContainer("cont","example of  container",nstep,nvar,iBin);
  AliCFContainer *cont = new AliCFContainer("cont","example of  container",nstep,nvar,iBin, 0x0,kTRUE);

  //setting the bin limits
  cont->SetBinLimits(ipt,binLim1);
  cont->SetBinLimits(ieta,binLim2);
  cont->SetBinLimits(iphi,binLim3);
  cont->SetBinLimits(ivtx,binLim4);

  //Start filling the mc and the data

  //data sample (1M tracks)
  Int_t nev=1000000;
  Int_t seed =1234;
  gRandom->SetSeed(seed);
  Double_t Value[nvar];
  for(Int_t iev =0;iev<nev;iev++){
    Float_t y=gRandom->Rndm();
    Float_t pt=-TMath::Log(y)/0.5; //pt, exponential
    Double_t eta=2.*gRandom->Rndm()-1.;//flat in eta
    Double_t phi=360.*gRandom->Rndm(); //flat in phi
    Float_t vtx=gRandom->Gaus( 0,5.);//gaussian in vertex
    Value[ipt]=pt;
    Value[ieta]=eta;
    Value[iphi]=phi;
    Value[ivtx]=vtx;    
    cont->Fill(Value, stepGen); //fill the efficiency denominator, sel step=0
    Float_t rndm=gRandom->Rndm();
    //simulate 80% constant efficiency everywhere
    if(rndm<0.8){
      cont->Fill(Value,stepRec); //fill the efficiency denominator, sel step =1
    }		
  }   

// Save it to a file
   cont->Save("container.root");
  //delete it
   delete cont;

// Read the  container from file
   TFile *file = new TFile("container.root");
   AliCFContainer *data = (AliCFContainer*) (file->Get("cont"));

  // Make some 1 & 2-D projections..
  // pt and vertex, generator and reconstructed level
  TCanvas *cmc =new TCanvas("cmc","The  distributions",0,300,900,900);
  cmc->Divide(2,2);
  cmc->cd(1);
  TH1D *hpt1a = data->ShowProjection(ipt, stepGen);
  hpt1a->SetMinimum(0.01);
  hpt1a->Draw();
  cmc->cd(2);
  TH1D *hpt1b = data->ShowProjection(ipt, stepRec);
  hpt1b->SetMinimum(0.01);
  hpt1b->Draw();
  cmc->cd(3);
  TH2D *hptvtx1a = data->ShowProjection(ipt,ivtx, stepGen);
  hptvtx1a->SetMinimum(0.01);
  hptvtx1a->Draw("lego");
  cmc->cd(4);
  TH2D *hptvtx1b = data->ShowProjection(ipt,ivtx, stepRec);
  hptvtx1b->SetMinimum(0.01);
  hptvtx1b->Draw("lego");
  cmc->Print("data.gif");

 
  //construct the efficiency grid from the data container 
  AliCFEffGrid *eff = new AliCFEffGrid("eff"," The efficiency",*data);
  eff->CalculateEfficiency(stepRec,stepGen); //eff= step1/step0

  //The efficiency along pt and vertex, and 2-D projection
  TCanvas *ceff =new TCanvas("ceff"," Efficiency",0,300,900,300);
  ceff->Divide(3,1);
  ceff->cd(1);
  TH1D *hpt2a = eff->Project(ipt); //the efficiency vs pt
  hpt2a->SetMinimum(0.01);
  hpt2a->Draw();
  ceff->cd(2);
  TH1D *hvtx2a = eff->Project(ivtx); //the efficiency vs vtx
  hvtx2a->SetMinimum(0.01);
  hvtx2a->Draw();
  ceff->cd(3);
  TH2D *hptvtx2a = eff->Project(ipt,ivtx); //look at the numerator
  hptvtx2a->SetMinimum(0.01);
  hptvtx2a->SetMinimum(0.01);
  hptvtx2a->Draw("lego");
  
  ceff->Print("eff.gif");

  //get the corrected data grid  
  AliCFDataGrid *corrdata = new AliCFDataGrid("corrdata","corrected data",*data);
  //correct selection step "reconstructed"
  corrdata->SetMeasured(stepRec); //set data to be corrected
  corrdata->ApplyEffCorrection(*eff);//apply the correction for efficiency

  //The observed data, the corrected ones and the "MC truth" distributions 
  //vs pt and vtx
  TCanvas *ccorrdata =new TCanvas("ccorrdata"," corrected data",0,300,900,900);
  ccorrdata->Divide(2,2);
  ccorrdata->cd(1);
  TH1D *hpt3a = corrdata->GetData()->Project(ipt); //uncorrected data
  hpt3a->SetMinimum(0.01);
  hpt3a->Draw();
  ccorrdata->cd(2);
  TH1D *hpt3b = corrdata->Project(ipt); //corrected data
  hpt3b->SetMinimum(0.01);
  hpt3b->Draw();
  ccorrdata->cd(3);
  TH1D *hvtx3a = corrdata->GetData()->Project(ivtx); //uncorrected data
  hvtx3a->SetMinimum(0.01);
  hvtx3a->Draw();
  ccorrdata->cd(4);
  TH1D *hvtx3b = corrdata->Project(ivtx); //corrected data
  hvtx3b->SetMinimum(0.01);
  hvtx3b->Draw();
  ccorrdata->Print("corrdata.gif");

}
