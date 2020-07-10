#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPythia6.h"
#include "AliTPythia8.h"
#include "AliDecayerPythia8.h"
#include "TPythia6Decayer.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TParticle.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "Riostream.h"
#include "TObjArray.h"
#include "TMCParticle.h"
#include "TStyle.h"
#include "TTree.h"
#include "TRandom3.h"
#endif


/// Macro to compute the pt-differential cross section of D mesons from B decays
///   using FONLL cross sections for B hadrons and Pythia decayer
///
/// Parameters:
///   nGener = number of B mesons to generate (in |y|<1)
///   pythiaver = decayer used
///                  6 --> TPythia6Decayer
///                  8 --> AliDecayerPythia8
///   fonllCase = FONLL prediction to be used
///                  0 --> central value
///                  1 --> max
///                  -1 --> min
///   fileNameFONLL* = files from FONLL with dsigma/dpt of B hadrons, D0, D+ and D* mesons
///                    (as taken from http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html)
///   opt4ff = set of fragmentation fractions f(b->B) to be used
///                  0 --> ppbar fractions from PDG
///                  1 --> Z decay fractions from PDG
///                  2 --> pT-dependent fractions from LHCb - cent hypothesis
///                  3 --> pT-dependent fractions from LHCb - min hypothesis
///                  4 --> pT-dependent fractions from LHCb - max hypothesis
///                  5 --> pT-dependent fractions from LHCb - min hypothesis + Lb energy scaling to 5 TeV
///   optForNorm = treatment of rapidity cut and normalization of xsec
///                  0 --> no cut on y(D) and normalisation ot xsec of B in |y|<0.5
///                  1 --> generate B in |yB|<1, cut on |yD|<0.5, count B in |yB|<0.5 for normalisation to xsec
///   writeTree = flag to control writing of Tree of decay kinematics


TH1D* ReadFONLL(TString filename, Int_t fonllCase, Int_t nPtBins, Double_t histominpt, Double_t histomaxpt);

void ComputeBtoDdecay(Int_t nGener=10000000,
		      Int_t pythiaver=8,
		      Int_t fonllCase=0,
		      TString fileNameFONLLb="FONLL-Bhadron-dsdpt-sqrts5020-100GeV-50MeVbins.txt",
		      TString fileNameFONLLd0="FONLL-D0-dsdpt-sqrts5020-100GeV-50MeVbins.txt",
		      TString fileNameFONLLdplus="FONLL-Dplus-dsdpt-sqrts5020-100GeV-50MeVbins.txt",
		      TString fileNameFONLLdave="FONLL-D0DplusAv-sqrts5020-100GeV-50MeVbins.txt",
		      TString fileNameFONLLdstar="FONLL-Dstar-dsdpt-sqrts5020-100GeV-50MeVbins.txt",
		      Int_t opt4ff=0,
		      Int_t optForNorm=1,
		      Bool_t writeTree=kFALSE){

  const Int_t nBeautyHadSpecies=4;
  Int_t pdgArrB[nBeautyHadSpecies]={511,521,531,5122};
  TString bhadrname[nBeautyHadSpecies]={"B0","Bplus","Bs","Lb"};
  Double_t fracB[4]={0.401,0.401,0.105,0.093};
  TF1 *fracU[15];
  TF1 *fracBs[15];
  TF1 *fracLb[15];
  TF1 *enScal;

  if(opt4ff==0){
    // ppbar fractions
    fracB[0]=0.340;
    fracB[1]=0.340;
    fracB[2]=0.101;
    fracB[3]=0.219;
  }else if(opt4ff==1){
    // e+e- fractions
    fracB[0]=0.412;
    fracB[1]=0.412;
    fracB[2]=0.088;
    fracB[3]=0.088;    
  }
  else if(opt4ff>=2){
    // pt-dependent fractions - evaluate when b hadron pt is calculated in the gen. loop
    // parameters for pT-dependence from LHCb beauty fraction measurement from https://arxiv.org/pdf/1902.06794.pdf
    fracB[0]=0.;
    fracB[1]=0.;
    fracB[2]=0.;
    fracB[3]=0.;    

    // FF uncertainty defined as the envelope of the variations of the fit parameters in LHCb paper
    // ipar=0 : central
    // 1 <= ipar <= 7 : upper uncertainties 
    // 8 <= ipar <= 14 : lower uncertainties 
    for(Int_t ipar=0;ipar<15;ipar++){

      fracU[ipar] = new TF1("fracU","1 /  (2 * (([0] * ([1] + [2] * (x - [3])))  + ([4] * ([5] + exp([6] + [7] * x))) + 1) ) ",0,50 );
      fracLb[ipar] = new TF1("fracLb","([4] * ([5] + exp([6] + [7] * x))) /  (([0] * ([1] + [2] * (x - [3])))  + ([4] * ([5] + exp([6] + [7] * x))) + 1)  ",0,50 );
      fracBs[ipar] = new TF1("fracBs","([0] * ([1] + [2] * (x - [3]))) /  (([0] * [1] + [2] * (x - [3]))  + ([4] * ([5] + exp([6] + [7] * x))) + 1)  ",0,50 );

      Double_t parLbA = 1;
      Double_t parLbp1 = 0.0793;
      Double_t parLbp2 = -1.022;
      Double_t parLbp3 = -0.107;
      Double_t parBsA = 1;
      Double_t parBsp1 = 0.119;
      Double_t parBsp2 = -0.00091;
      Double_t parBsAvePt = 10.1;

      // positive
      if(ipar==1) parLbA += 0.061;
      if(ipar==2) parLbp1 += 0.0141;
      if(ipar==3) parLbp2 += 0.0047;
      if(ipar==4) parLbp3 += 0.002;
      if(ipar==5) parBsA += 0.043;
      if(ipar==6) parBsp1 += 0.001;
      if(ipar==7) parBsp2 += 0.00025;
      // negative errors
      if(ipar==8) parLbA -= 0.061;
      if(ipar==9) parLbp1 -= 0.0141;
      if(ipar==10) parLbp2 -= 0.0047;
      if(ipar==11) parLbp3 -= 0.002;
      if(ipar==12) parBsA -= 0.043;
      if(ipar==13) parBsp1 -= 0.001;
      if(ipar==14) parBsp2 -= 0.00025;

      fracU[ipar]->SetParameter(0,parBsA);
      fracU[ipar]->SetParameter(1,parBsp1);
      fracU[ipar]->SetParameter(2,parBsp2);
      fracU[ipar]->SetParameter(3,parBsAvePt);
      fracU[ipar]->SetParameter(4,parLbA);
      fracU[ipar]->SetParameter(5,parLbp1);
      fracU[ipar]->SetParameter(6,parLbp2);
      fracU[ipar]->SetParameter(7,parLbp3);
      fracBs[ipar]->SetParameter(0,parBsA);
      fracBs[ipar]->SetParameter(1,parBsp1);
      fracBs[ipar]->SetParameter(2,parBsp2);
      fracBs[ipar]->SetParameter(3,parBsAvePt);
      fracBs[ipar]->SetParameter(4,parLbA);
      fracBs[ipar]->SetParameter(5,parLbp1);
      fracBs[ipar]->SetParameter(6,parLbp2);
      fracBs[ipar]->SetParameter(7,parLbp3);
      fracLb[ipar]->SetParameter(0,parBsA);
      fracLb[ipar]->SetParameter(1,parBsp1);
      fracLb[ipar]->SetParameter(2,parBsp2);
      fracLb[ipar]->SetParameter(3,parBsAvePt);
      fracLb[ipar]->SetParameter(4,parLbA);
      fracLb[ipar]->SetParameter(5,parLbp1);
      fracLb[ipar]->SetParameter(6,parLbp2);
      fracLb[ipar]->SetParameter(7,parLbp3);

    }
    if(opt4ff==5) {
      enScal = new TF1("enScal","[0]",0,50);
      enScal->SetParameter(0,0.796357); // average difference between 13 TeV/7 TeV Lb/B ratio, + 50%
    }
  }

  const Int_t nCharmHadSpecies=5;
  Int_t pdgArrC[nCharmHadSpecies]={421,411,431,4122,413};
  TString chadrname[nCharmHadSpecies]={"D0","Dplus","Ds","Lc","Dstar"};
  Double_t fracC[nCharmHadSpecies]={0.542,0.225,0.092,0.057,0.236}; // Values from e+e- ARXIV:1404.3888 (D0, D+, Ds, Lc, D*+)
  Int_t cols[nCharmHadSpecies]={2,4,kGreen+1,kMagenta+1,kYellow+1};
  
  Int_t nPtBins=2001;
  Double_t ptmin=0.;
  Double_t ptmax=100.05;

  TVirtualMCDecayer* pdec=0x0;
  
  if(pythiaver==6){
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    //    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
    pdec=new TPythia6Decayer();
  }else{
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    pdec=new AliDecayerPythia8();
  }

  TDatabasePDG* db=TDatabasePDG::Instance();
  pdec->Init();

  TH1D *hBptDistr = ReadFONLL(fileNameFONLLb.Data(),fonllCase,nPtBins,ptmin,ptmax);
  hBptDistr->SetName("hfonllB");
  hBptDistr->Scale(1.e-6); // convert to ub
  TH1D *hbFragmFrac = new  TH1D("hbFragmFrac"," ; ; Fragmentation Fraction",nBeautyHadSpecies,-0.5,nBeautyHadSpecies-0.5);
  hbFragmFrac->SetStats(0);
  for(Int_t ib=0; ib<nBeautyHadSpecies; ib++){
    hbFragmFrac->GetXaxis()->SetBinLabel(ib+1,bhadrname[ib].Data());
    hbFragmFrac->SetBinContent(ib+1,fracB[ib]);
  }

  TH1D** hpromptDpt=new TH1D*[nCharmHadSpecies];
  TH1D *hcFragmFrac = new  TH1D("hcFragmFrac"," ; ; Fragmentation Fraction",nCharmHadSpecies,-0.5,nCharmHadSpecies-0.5);
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    if(pdgArrC[ic]==411) hpromptDpt[ic] = ReadFONLL(fileNameFONLLdplus.Data(),fonllCase,nPtBins,ptmin,ptmax);
    else if(pdgArrC[ic]==421) hpromptDpt[ic] = ReadFONLL(fileNameFONLLd0.Data(),fonllCase,nPtBins,ptmin,ptmax);
    else if(pdgArrC[ic]==413) hpromptDpt[ic] = ReadFONLL(fileNameFONLLdstar.Data(),fonllCase,nPtBins,ptmin,ptmax);
    else hpromptDpt[ic] = ReadFONLL(fileNameFONLLdave.Data(),fonllCase,nPtBins,ptmin,ptmax);
    hpromptDpt[ic]->Scale(fracC[ic]);
    hpromptDpt[ic]->SetName(Form("hfonllPrompt%s",chadrname[ic].Data()));
    hpromptDpt[ic]->Scale(1.e-6); // convert to ub
    hcFragmFrac->GetXaxis()->SetBinLabel(ic+1,chadrname[ic].Data());
    hcFragmFrac->SetBinContent(ic+1,fracC[ic]);
  }
  
  Double_t xsecb=0;
  for(Int_t i=1; i<=hBptDistr->GetNbinsX(); i++){
    xsecb+=(hBptDistr->GetBinContent(i)*hBptDistr->GetBinWidth(i));
  }
  Double_t xsecpromptD[nCharmHadSpecies];
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    xsecpromptD[ic]=0;
    for(Int_t i=1; i<=hpromptDpt[ic]->GetNbinsX(); i++){
      xsecpromptD[ic]+=(hpromptDpt[ic]->GetBinContent(i)*hpromptDpt[ic]->GetBinWidth(i));
    }
  }

  
  TH1D** hnonpromptDorigin=new TH1D*[nCharmHadSpecies];
  TH1D** hnonpromptDpt=new TH1D*[nCharmHadSpecies];
  TH2D** hnonpromptDptByOrigin=new TH2D*[nCharmHadSpecies];
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    hnonpromptDorigin[ic] = new TH1D(Form("hnonprompt%sOrigin",chadrname[ic].Data()),Form("%s mother ; ; Entries",chadrname[ic].Data()),nBeautyHadSpecies,-0.5,nBeautyHadSpecies-0.5);
    for(Int_t ib=0; ib<nBeautyHadSpecies; ib++)  hnonpromptDorigin[ic]->GetXaxis()->SetBinLabel(ib+1,bhadrname[ib].Data());
    hnonpromptDpt[ic]  = new TH1D(Form("hnonprompt%spt",chadrname[ic].Data())," ; p_{T} (GeV) ; d#sigma/dp_{T} (#mub/GeV)",nPtBins,ptmin,ptmax);
    hnonpromptDptByOrigin[ic]  = new TH2D(Form("hnonprompt%sptByOrigin",chadrname[ic].Data()),"",nBeautyHadSpecies,-0.5,nBeautyHadSpecies-0.5,nPtBins,ptmin,ptmax);
    for(Int_t ib=0; ib<nBeautyHadSpecies; ib++) hnonpromptDptByOrigin[ic]->GetXaxis()->SetBinLabel(ib+1,bhadrname[ib].Data());
    hnonpromptDptByOrigin[ic]->GetYaxis()->SetTitle("p_{T} (GeV)");
  }

  TH1D** hBhadDau = new TH1D*[nBeautyHadSpecies];
  for(Int_t ib=0; ib<nBeautyHadSpecies; ib++){
    hBhadDau[ib] = new TH1D(Form("h%sdau",bhadrname[ib].Data())," ; ; Entries",nCharmHadSpecies+1,-1.5,nCharmHadSpecies-0.5);
    hBhadDau[ib]->GetXaxis()->SetBinLabel(1,Form("All %s",bhadrname[ib].Data()));
    for(Int_t ic=0; ic<nCharmHadSpecies; ic++) hBhadDau[ib]->GetXaxis()->SetBinLabel(2+ic,chadrname[ic].Data());
  }
  
  TH2F** hDptVsBpt = new TH2F*[nCharmHadSpecies*nBeautyHadSpecies];
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    for(Int_t ib=0; ib<nBeautyHadSpecies; ib++){
      hDptVsBpt[ib*nCharmHadSpecies+ic] = new TH2F(Form("h%sptVs%spt",chadrname[ic].Data(),bhadrname[ib].Data()),Form(" ; p_{T}(%s) ; p_{T}(%s)",bhadrname[ib].Data(),chadrname[ic].Data()),200,0.,100.,200.,0.,100.);
    }
  }
  TH1D *hUfrac = new TH1D("hUfrac","hUfrac",250,0,50);
  TH1D *hBsfrac = new TH1D("hBsfrac","hBsfrac",250,0,50);
  TH1D *hLbfrac = new TH1D("hLbfrac","hLbfrac",250,0,50);

  TTree* fTreeDecays = 0x0;
  Int_t pdgB = -9999;
  Double_t ptB = -1.;
  Double_t pB = -1.;
  Double_t yB = -1.;
  vector<float> arrptD;
  vector<float> arrpD;
  vector<float> arryD;
  vector<int> arrpdgD;
  Double_t norm = xsecb;

  if(writeTree){
    fTreeDecays = new TTree("fTreeDecays", "fTreeDecays");
    fTreeDecays->Branch("pdgB", &pdgB);
    fTreeDecays->Branch("ptB", &ptB);
    fTreeDecays->Branch("pB", &pB);
    fTreeDecays->Branch("yB", &yB);
    fTreeDecays->Branch("ptD", &arrptD);
    fTreeDecays->Branch("pD", &arrpD);
    fTreeDecays->Branch("yD", &arryD);
    fTreeDecays->Branch("pdgD", &arrpdgD);
    fTreeDecays->Branch("norm", &norm);
  }
  
  TRandom3* gener=new TRandom3(0);
  TClonesArray *array = new TClonesArray("TParticle",100);
  TLorentzVector* vec=new TLorentzVector();

  Double_t countB=0;

  for(Int_t itry=0; itry<nGener; itry++){
    if(itry%10000==0) printf("Particle %d\n",itry);

    Int_t iBhad=0;
    Double_t value=gener->Rndm();
    TH2F* hptD0tofill=0x0;
    TH2F* hptDptofill=0x0;
    TH2F* hptDstofill=0x0;
    TH2F* hptLctofill=0x0;
    ptB=hBptDistr->GetRandom();
    if(opt4ff>=2){
      if(opt4ff==2){
        fracB[0] = fracU[0]->Eval(ptB>5?ptB:5);  
        fracB[1] = fracU[0]->Eval(ptB>5?ptB:5);  
        fracB[2] = fracBs[0]->Eval(ptB>5?ptB:5); 
        fracB[3] = fracLb[0]->Eval(ptB>5?ptB:5); 
      }
      else{ 
        fracB[0] = opt4ff==4?0:1;
        fracB[1] = opt4ff==4?0:1;
        fracB[2] = opt4ff==4?0:1;
        fracB[3] = opt4ff==4?0:1;
        for(Int_t ipar=1;ipar<15;ipar++){
          if(opt4ff==3){ // minimum
            Double_t fracLbtry = fracLb[ipar]->Eval(ptB>5?ptB:5);
            if(fracLbtry<fracB[3]){
              fracB[0] = fracU[ipar]->Eval(ptB>5?ptB:5);
              fracB[1] = fracU[ipar]->Eval(ptB>5?ptB:5);
              fracB[2] = fracBs[ipar]->Eval(ptB>5?ptB:5);
              fracB[3] = fracLb[ipar]->Eval(ptB>5?ptB:5);
            }
          }
          else if(opt4ff==4){ // maximum
            Double_t fracLbtry = fracLb[ipar]->Eval(ptB);
            if(fracLbtry>fracB[3]){
              fracB[0] = fracU[ipar]->Eval(ptB);
              fracB[1] = fracU[ipar]->Eval(ptB);
              fracB[2] = fracBs[ipar]->Eval(ptB);
              fracB[3] = fracLb[ipar]->Eval(ptB);
            }
          }
          else if(opt4ff==5){ // minimum, with additional energy uncertainty
            Double_t fracLbtry = fracLb[ipar]->Eval(ptB>5?ptB:5);
            Double_t scaleLb; 
            if(opt4ff==5) {
              scaleLb = enScal->Eval(ptB);
            }
            if(fracLbtry<fracB[3]){
              Double_t diffLb = fracLb[ipar]->Eval(ptB>5?ptB:5) * (1. -  scaleLb);
              fracB[0] = fracU[ipar]->Eval(ptB>5?ptB:5) ;
              fracB[1] = fracU[ipar]->Eval(ptB>5?ptB:5) ;
              fracB[2] = fracBs[ipar]->Eval(ptB>5?ptB:5) + diffLb;
              fracB[3] = fracLb[ipar]->Eval(ptB>5?ptB:5) * scaleLb;
            }
          }
        }
      }
      hUfrac->SetBinContent(hUfrac->FindBin(ptB),fracB[0]);
      hBsfrac->SetBinContent(hBsfrac->FindBin(ptB),fracB[2]);
      hLbfrac->SetBinContent(hLbfrac->FindBin(ptB),fracB[3]);
    }
    if(value<fracB[0]){ 
      pdgB=pdgArrB[0];
      iBhad=0;
   }else if(value<(fracB[0]+fracB[1])){ 
      pdgB=pdgArrB[1];
      iBhad=1;
    }else if(value<(fracB[0]+fracB[1]+fracB[2])){ 
      pdgB=pdgArrB[2];
      iBhad=2;
    }else{
      pdgB=pdgArrB[3];
      iBhad=3;
    }
    hBhadDau[iBhad]->Fill(-1);
    
    Double_t mass=db->GetParticle(pdgB)->Mass();
    Double_t phiB=gener->Rndm()*2*TMath::Pi();
    yB=gener->Rndm()*2.-1.; // flat in -1<y<1
    Double_t px=ptB*TMath::Cos(phiB);
    Double_t py=ptB*TMath::Sin(phiB);
    Double_t mt=TMath::Sqrt(mass*mass+ptB*ptB);
    Double_t pz=mt*TMath::SinH(yB);
    pB=TMath::Sqrt(ptB*ptB+pz*pz);
    Double_t E=TMath::Sqrt(mass*mass+pB*pB);
    vec->SetPxPyPzE(px,py,pz,E);
    pdec->Decay(pdgB,vec);
    if(optForNorm==0) countB+=1.;
    else if(optForNorm==1 && TMath::Abs(yB)<0.5) countB+=1.;
    Int_t nentries = pdec->ImportParticles(array);
    //    TParticle* bmes=(TParticle*)array->At(0);

    for(int j=0; j<nentries; j++){
      TParticle * part = (TParticle*)array->At(j);
      Int_t pdgdau=TMath::Abs(part->GetPdgCode());
      Double_t ptD=-999;
      Double_t yD=-999;
      Int_t iChad=-999;
      for(Int_t ic=0; ic<nCharmHadSpecies; ic++) if(pdgdau==pdgArrC[ic]) iChad=ic;
      if(iChad>=0){
	Double_t ptD=part->Pt();
	Double_t yD=part->Y();
	arrptD.push_back(ptD);
	arrpD.push_back(part->P());
	arryD.push_back(yD);
	arrpdgD.push_back(pdgdau);
	hBhadDau[iBhad]->Fill(iChad); // filled outside the cut on yD to recover correctly the BRs
	if(optForNorm==0 || (optForNorm==1 && TMath::Abs(yD)<0.5)){
	  hnonpromptDorigin[iChad]->Fill(iBhad);
	  hnonpromptDpt[iChad]->Fill(ptD);
	  hnonpromptDptByOrigin[iChad]->Fill(iBhad,ptD);
	  hDptVsBpt[iBhad*nCharmHadSpecies+iChad]->Fill(ptB,ptD);
	}
      }
    }
    if(arrptD.size() == 0){
      arrptD.push_back(-1);
      arrpD.push_back(-1);
      arryD.push_back(-999);
      arrpdgD.push_back(-1);
    }
    if(fTreeDecays) fTreeDecays->Fill();
    arrptD.clear();
    arrpD.clear();
    arryD.clear();
    arrpdgD.clear();
    array->Clear();
  }
  
  delete vec;
  delete array;
  delete gener;

  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    hnonpromptDpt[ic]->Scale(xsecb/countB/hnonpromptDpt[ic]->GetBinWidth(1));
    hnonpromptDpt[ic]->SetStats(0);
    hnonpromptDptByOrigin[ic]->Scale(xsecb/countB/hnonpromptDptByOrigin[ic]->GetYaxis()->GetBinWidth(1));
  }

  printf("Cross sections for prompt charm:");
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    printf(" %s = %f ub   ",chadrname[ic].Data(),xsecpromptD[ic]);
  }
  printf("\n");

  printf("Cross section for B = %f ub \n",xsecb);
  Double_t xsecD[nCharmHadSpecies];
  printf("Cross sections for feeddown charm:");
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    xsecD[ic]=0.;
    for(Int_t i=1; i<=hnonpromptDpt[ic]->GetNbinsX(); i++){
      xsecD[ic]+=(hnonpromptDpt[ic]->GetBinContent(i)*hnonpromptDpt[ic]->GetBinWidth(i));
    }
    printf(" %s = %f ub   ",chadrname[ic].Data(),xsecD[ic]);
  }
  printf("\n");

  TH1D** hfpromptD = new TH1D*[nCharmHadSpecies];
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    hfpromptD[ic] = new TH1D(Form("hfprompt%s",chadrname[ic].Data())," ; p_{T} (GeV) ; f_{prompt}",hnonpromptDpt[ic]->GetNbinsX(),hnonpromptDpt[ic]->GetXaxis()->GetXmin(),hnonpromptDpt[ic]->GetXaxis()->GetXmax());
    hfpromptD[ic]->SetStats(0);
    for(Int_t i=1; i<=hnonpromptDpt[ic]->GetNbinsX(); i++){
      Double_t fp=hpromptDpt[ic]->GetBinContent(i)/(hpromptDpt[ic]->GetBinContent(i)+hnonpromptDpt[ic]->GetBinContent(i));
      Double_t efp=fp/hpromptDpt[ic]->GetBinContent(i)*hnonpromptDpt[ic]->GetBinError(i);
      hfpromptD[ic]->SetBinContent(i,fp);
      hfpromptD[ic]->SetBinError(i,efp);
    }
  }
  
  hBptDistr->SetStats(0);
  hBptDistr->GetYaxis()->SetTitle("d#sigma/dp_{T} (#mub/GeV)");
  hBptDistr->GetXaxis()->SetTitle("p_{T} (GeV)");

  TH1D* hnonpromptDsKKpipt=(TH1D*)hnonpromptDpt[2]->Clone("hnonpromptDsKKpipt");
  hnonpromptDsKKpipt->Scale(0.0227);
  hnonpromptDsKKpipt->GetYaxis()->SetTitle("d#sigma/dp_{T}xBR (#mub/GeV)");
    
  TCanvas* c1=new TCanvas("c1","B mother",1500,900);
  c1->Divide(2,2);
  c1->cd(1);
  hnonpromptDorigin[0]->Draw();
  c1->cd(2);
  hnonpromptDorigin[1]->Draw();
  c1->cd(3);
  hnonpromptDorigin[2]->Draw();
  c1->cd(4);
  hnonpromptDorigin[3]->Draw();

  
  TCanvas* c2=new TCanvas("c2","PtD vs PtB",1500,1000);
  c2->Divide(nCharmHadSpecies,nBeautyHadSpecies);
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    for(Int_t ib=0; ib<nBeautyHadSpecies; ib++){      
      c2->cd(1+ib*nCharmHadSpecies+ic);
      hDptVsBpt[ib*nCharmHadSpecies+ic]->Draw("colz");
    }
  }
  c2->SaveAs(Form("DecayKine_Pythia%d.png",pythiaver));
  
  TCanvas* c3=new TCanvas("c3","pt-diff xsec",900,800);
  gPad->SetLogy();
  hBptDistr->Draw();
  hpromptDpt[0]->SetLineColor(kRed-9);
  hpromptDpt[0]->Draw("lsame");
  hpromptDpt[1]->SetLineColor(kBlue-9);
  hpromptDpt[1]->Draw("lsame");
  hpromptDpt[4]->SetLineColor(kYellow);
  hpromptDpt[4]->Draw("lsame");
  TLegend* leg=new TLegend(0.5,0.6,0.89,0.89);
  leg->AddEntry(hBptDistr,"B hadron, FONLL","L")->SetTextColor(hBptDistr->GetLineColor());
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    hnonpromptDpt[ic]->SetLineColor(cols[ic]);
    hnonpromptDpt[ic]->Draw("same");
    leg->AddEntry(hnonpromptDpt[ic],Form("%s #leftarrowB (FONLL+PYTHIA%d)",chadrname[ic].Data(),pythiaver),"L")->SetTextColor(cols[ic]);
  }
  leg->AddEntry(hpromptDpt[0],"Prompt D^{0}, FONLL","L")->SetTextColor(hpromptDpt[0]->GetLineColor());
  leg->AddEntry(hpromptDpt[1],"Prompt D^{+}, FONLL","L")->SetTextColor(hpromptDpt[1]->GetLineColor());
  leg->AddEntry(hpromptDpt[4],"Prompt D^{*+}, FONLL","L")->SetTextColor(hpromptDpt[4]->GetLineColor());
  leg->Draw();
  c3->SaveAs(Form("XsecBandDfromB_FONLLPythia%d.png",pythiaver));

  TCanvas* c4=new TCanvas("c4","fprompt",900,800);
  TH2F* hframe=new TH2F("hframe"," ; p_{T} (GeV) ; f_{prompt}",nPtBins,ptmin,ptmax,500,0.,1.);
  hframe->SetStats(0);
  hframe->Draw();
  TLegend* legf=new TLegend(0.15,0.15,0.35,0.4);
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    hfpromptD[ic]->SetLineColor(cols[ic]);
    hfpromptD[ic]->Draw("same");
    legf->AddEntry(hfpromptD[ic],chadrname[ic].Data(),"L")->SetTextColor(cols[ic]);
  }
  legf->Draw();


  
  TString outfilnam="DfromB_FONLL";
  if(fonllCase==1) outfilnam.Append("max");
  else if(fonllCase==-1) outfilnam.Append("min");
  else outfilnam.Append("cent");
  outfilnam.Append(Form("Pythia%d",pythiaver));
  if(opt4ff==0) outfilnam.Append("_FFppbar");
  else if(opt4ff==1) outfilnam.Append("_FFee");
  else if(opt4ff==2) outfilnam.Append("_FFptDepcent");
  else if(opt4ff==3) outfilnam.Append("_FFptDepmin");
  else if(opt4ff==4) outfilnam.Append("_FFptDepmax");
  else if(opt4ff==5) outfilnam.Append("_FFptDepminEnScaleConst");
  else outfilnam.Append("_FFold");
  if(optForNorm==1) outfilnam.Append("_yDcut");
  outfilnam.Append(".root");
  TFile* outfil=new TFile(outfilnam.Data(),"recreate");
  hcFragmFrac->Write();
  hbFragmFrac->Write();
  hBptDistr->Write();
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    hpromptDpt[ic]->Write();
    hnonpromptDpt[ic]->Write();
    hnonpromptDorigin[ic]->Write();
    hnonpromptDptByOrigin[ic]->Write();
  }
  hnonpromptDsKKpipt->Write();
  for(Int_t ib=0; ib<nBeautyHadSpecies; ib++){
    hBhadDau[ib]->Write();
  }
  for(Int_t ic=0; ic<nCharmHadSpecies; ic++){
    for(Int_t ib=0; ib<nBeautyHadSpecies; ib++){      
      hDptVsBpt[ib*nCharmHadSpecies+ic]->Write();
    }
  }
  if(fTreeDecays) fTreeDecays->Write();
  if(opt4ff>=2) {
    hUfrac->Write();
    hBsfrac->Write();
    hLbfrac->Write();
  }
  outfil->Close();
}


//----------------------------------------------
TH1D* ReadFONLL(TString filename, Int_t fonllCase, Int_t nPtBins, Double_t histominpt, Double_t histomaxpt){
  FILE* infil=fopen(filename.Data(),"r");
  Char_t line[200];
  Char_t* rc;
  for(Int_t il=0; il<16; il++){
    rc=fgets(line,200,infil);
    if(strstr(line,"central")) break;
  }
  Float_t pt,csc,csmin,csmax,dum;
  Double_t ptmin=999,ptmax=0.;
  Int_t iPt=0;
  Double_t x[2000],y[2000];
  Bool_t ok;
  while(!feof(infil)){
    ok=fscanf(infil,"%f %f %f %f",&pt,&csc,&csmin,&csmax);
    for(Int_t i=0; i<12;i++) ok=fscanf(infil,"%f",&dum);
    if(feof(infil)) break;
    if(pt==0.) continue;
    if(pt<ptmin) ptmin=pt;
    if(pt>ptmax) ptmax=pt;
    x[iPt]=pt;
    if(fonllCase==-1) y[iPt]=csmin;
    else if(fonllCase==1) y[iPt]=csmax;
    else  y[iPt]=csc;
    iPt++;
  }
  fclose(infil);
  Double_t binw=(ptmax-ptmin)/(iPt-1);
  if(iPt!=nPtBins){
    printf("ERROR: different number of pt bins in FONLL (%d) and histos (%d)\n",iPt,nPtBins);
  }
  Double_t lowlim=ptmin-0.5*binw;
  if(TMath::Abs(lowlim-histominpt)>0.0001){
    printf("ERROR: different lower limit for pt in FONLL (%f) and histos (%f)\n",lowlim,histominpt);
  }else{
    lowlim=histominpt; // to avoid numerical precision problems with axis limits
  }
  Double_t higlim=ptmax+0.5*binw;
  if(TMath::Abs(higlim-histomaxpt)>0.0001){
    printf("ERROR: different upper limit for pt in FONLL (%f) and histos (%f)\n",higlim,histomaxpt);
  }else{
    higlim=histomaxpt; // to avoid numerical precision problems with axis limits
  }
  TH1D* hfonll=new TH1D("hfonll","",iPt,lowlim,higlim);
  for(Int_t iBin=0; iBin<iPt; iBin++){
    hfonll->SetBinContent(iBin+1,y[iBin]);
    hfonll->SetBinError(iBin+1,0.);
  }
  return hfonll;

}
