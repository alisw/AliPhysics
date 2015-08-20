#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TPythia6Decayer.h>
#include <TPaveStats.h>
#endif

enum EDDecay{kD0Kpi,kDplusKpipi,kDstarD0pi,kDsKKpi};
enum EFidY{kFixedY,kPtDepY};
enum EPtShape{kFlat,kFONLL7TeV,kPythia7TeV};

// Configuration
Int_t fDDecay=kD0Kpi;
Double_t fPtMinDau=0.1;
Double_t fEtaMaxDau=0.9;
Int_t fOptionYFiducial=kFixedY;
Double_t fYMaxFidAccCut=0.8;
Int_t fPtShape=kFONLL7TeV;
TString fDecayTableFileName="$ALICE_PHYSICS/../src/PWGHF/vertexingHF/macros/decaytable_acc.dat"; 
Int_t fDebugLevel=0;
Int_t totTrials=1000000;


Bool_t CountKpi(TClonesArray *array, Int_t nentries, Int_t &nPions, Int_t &nKaons, Int_t &nPionsInAcc, Int_t &nKaonsInAcc);
Bool_t IsInFiducialAcceptance(Double_t pt, Double_t y);

void ComputeAcceptance(){
  // main function
  
  gSystem->Load("liblhapdf.so");      // Parton density functions
  gSystem->Load("libEGPythia6.so");   // TGenerator interface
  gSystem->Load("libpythia6.so");     // Pythia
  
  TPythia6Decayer* pdec=TPythia6Decayer::Instance();
  if(fDecayTableFileName.CompareTo("")!=0){
    if(fDecayTableFileName.Contains("ALICE_PHYSICS")){
      gSystem->Exec(Form("cp %s .",fDecayTableFileName.Data()));
      fDecayTableFileName.ReplaceAll("$ALICE_PHYSICS/../src/PWGHF/vertexingHF/macros/","./");
    }
    pdec->SetDecayTableFile(fDecayTableFileName.Data());
    pdec->ReadDecayTable();
  }
  pdec->Init();

  Int_t pdgCode=0;
  Int_t nPionDau=-1;
  Int_t nKaonDau=-1;
  TString outFileName="Acceptance_Toy_";
  if(fDDecay==kD0Kpi){
    pdgCode=421;
    nPionDau=1;
    nKaonDau=1;
    outFileName.Append("D0Kpi_");
  }else if(fDDecay==kDplusKpipi){
    pdgCode=411;
    nPionDau=2;
    nKaonDau=1;
    outFileName.Append("DplusKpipi_");
  }else if(fDDecay==kDstarD0pi){
    pdgCode=413;
    nPionDau=2;
    nKaonDau=1;
    outFileName.Append("DStarD0pi_");
  }else if(fDDecay==kDsKKpi){
    pdgCode=431;
    nPionDau=1;
    nKaonDau=2;
    outFileName.Append("DsKKpi_");
  }else{
    printf("ERROR: Wrong decay selected\n");
    return;
  }
  if(fOptionYFiducial==kFixedY) outFileName.Append(Form("yfid%02d_",(Int_t)(fYMaxFidAccCut*10)));
  else outFileName.Append("yfidPtDep_");
  outFileName.Append(Form("etaDau%02d_",(Int_t)(fEtaMaxDau*10)));
  outFileName.Append(Form("ptDau%d_",(Int_t)(fPtMinDau*1000)));
  TDatabasePDG* db=TDatabasePDG::Instance();
  Float_t massD=db->GetParticle(pdgCode)->Mass();
  TClonesArray *array = new TClonesArray("TParticle",100);

  TH2D* hPtVsYGen=new TH2D("hPtVsYGen","",400,0.,40.,20.,-1.,1.);
  hPtVsYGen->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtVsYGen->GetYaxis()->SetTitle("y");
  TH2D* hPtVsYGenLimAcc=new TH2D("hPtVsYGenLimAcc","",400,0.,40.,20.,-1.,1.);
  hPtVsYGenLimAcc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtVsYGenLimAcc->GetYaxis()->SetTitle("y");
  TH2D* hPtVsYGenAcc=new TH2D("hPtVsYGenAcc","",400,0.,40.,20.,-1.,1.);
  hPtVsYGenAcc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtVsYGenAcc->GetYaxis()->SetTitle("y");

  TF1* funcPt=0x0;
  if(fPtShape==kFONLL7TeV){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.322643,2.96275,2.30301,2.5);
    outFileName.Append("FONLL7ptshape.root");
  }else if(fPtShape==kPythia7TeV){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.322643,1.94635,1.40463,2.5);
    outFileName.Append("PYTHIA7ptshape.root");  
  }else{
    funcPt=new TF1("fFlat","pol0",0.,40.);
    funcPt->SetParameter(0,1.);
    outFileName.Append("flatpt.root");
  }
  TRandom3* gener=new TRandom3(0);

  for(Int_t itry=0; itry<totTrials; itry++){
    if(itry%10000==0) printf("Event %d\n",itry);
    Float_t ptD=funcPt->GetRandom();
    Float_t phiD=gener->Rndm()*2*TMath::Pi();
    Float_t yD=gener->Rndm()*2.-1.; // flat in -1<y<1
    Float_t px=ptD*TMath::Cos(phiD);
    Float_t py=ptD*TMath::Sin(phiD);
    Float_t mt=TMath::Sqrt(massD*massD+ptD*ptD);
    Float_t pz=mt*TMath::SinH(yD);
    Float_t E=TMath::Sqrt(massD*massD+px*px+py*py+pz*pz);

    TLorentzVector* vec=new TLorentzVector(px,py,pz,E);
    pdec->Decay(pdgCode,vec);
    array->Clear();
    Int_t nentries = pdec->ImportParticles(array);
    TParticle* dmes=(TParticle*)array->At(0);
    Int_t nDaughters=dmes->GetNDaughters();
    if(fDDecay==kD0Kpi && nDaughters!=2) return;
    Int_t nPionsInAcc=0;
    Int_t nKaonsInAcc=0;
    Int_t nPions=0;
    Int_t nKaons=0;
    Bool_t isOk=CountKpi(array,nentries,nPions,nKaons,nPionsInAcc,nKaonsInAcc);

    if(isOk){
      if(nPions==nPionDau && nKaons==nKaonDau){
	hPtVsYGen->Fill(ptD,yD);
	if(TMath::Abs(yD)<0.5){
	  hPtVsYGenLimAcc->Fill(ptD,yD);
	}
	if(IsInFiducialAcceptance(ptD,yD)){	  
	  if(nPionsInAcc==nPionDau && nKaonsInAcc==nKaonDau){ 
	    hPtVsYGenAcc->Fill(ptD,yD);
	  }
	}
      }
      delete vec;
    }
  }

  TH1D* hPtGenAcc=(TH1D*)hPtVsYGenAcc->ProjectionX("hPtGenAcc"); 
  hPtGenAcc->GetYaxis()->SetTitle("Entries");
  TH1D* hPtGenLimAcc=(TH1D*)hPtVsYGenLimAcc->ProjectionX("hPtGenLimAcc"); 
  hPtGenLimAcc->GetYaxis()->SetTitle("Entries");
  hPtGenAcc->Sumw2();
  hPtGenLimAcc->Sumw2();
  TH1D* hAccVsPt=(TH1D*)hPtGenAcc->Clone("hAccVsPt");
  hAccVsPt->Divide(hPtGenAcc,hPtGenLimAcc,1,1,"B");
  hAccVsPt->GetYaxis()->SetTitle("Acceptance");
  hAccVsPt->SetStats(0);

  TCanvas* c2d=new TCanvas("c2d","Pt vs y",1200,600);
  c2d->Divide(3,1);
  c2d->cd(1);
  hPtVsYGen->Draw("colz");
  c2d->cd(2);
  hPtVsYGenLimAcc->Draw("colz");
  c2d->cd(3);
  hPtVsYGenAcc->Draw("colz");

  TCanvas* c1d=new TCanvas("c1d","Acceptance",1200,600);
  c1d->Divide(2,1);
  c1d->cd(1);
  hPtGenLimAcc->Draw("");
  Double_t ymax=1.2*TMath::Max(hPtGenLimAcc->GetMaximum(),hPtGenAcc->GetMaximum());
  hPtGenLimAcc->SetMaximum(ymax);
  gPad->Update();
  TPaveStats *st1=(TPaveStats*)hPtGenLimAcc->GetListOfFunctions()->FindObject("stats");
  st1->SetY1NDC(0.71);
  st1->SetY2NDC(0.9);
  hPtGenAcc->SetLineColor(kRed+1);
  hPtGenAcc->Draw("sames");
  gPad->Update();
  TPaveStats *st2=(TPaveStats*)hPtGenAcc->GetListOfFunctions()->FindObject("stats");
  st2->SetY1NDC(0.51);
  st2->SetY2NDC(0.7);
  st2->SetTextColor(hPtGenAcc->GetLineColor());
  gPad->Modified();
  c1d->cd(2);
  hAccVsPt->Draw();


  TFile* outfil=new TFile(outFileName.Data(),"recreate");
  hPtVsYGen->Write();
  hPtVsYGenLimAcc->Write();
  hPtVsYGenAcc->Write();
  hPtGenLimAcc->Write();
  hPtGenAcc->Write();
  hAccVsPt->Write();
  outfil->Close();

}

//___________________________________________________
Bool_t IsInFiducialAcceptance(Double_t pt, Double_t y){
  // check fiducial acceptance

  if(fOptionYFiducial==kFixedY){
    if(TMath::Abs(y) > fYMaxFidAccCut) return kFALSE;
    else return kTRUE;
  }

  if(pt > 5.) {
    if (TMath::Abs(y) > 0.8) return kFALSE;
  } else {
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5; 
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;		
    if (y < minFiducialY || y > maxFiducialY) return kFALSE;
  }

  return kTRUE;
}
//___________________________________________________
Bool_t CountKpi(TClonesArray *array, Int_t nentries, Int_t &nPions, Int_t &nKaons, Int_t &nPionsInAcc, Int_t &nKaonsInAcc){
  // count K and pi in Acc

  TParticle* dmes=(TParticle*)array->At(0);
  Double_t sumPx=0;
  Double_t sumPy=0;
  Double_t sumPz=0;
  for(int j=0; j<nentries; j++){
    TParticle * o = (TParticle*)array->At(j);
    Int_t pdgdau=TMath::Abs(o->GetPdgCode());
    if(fDebugLevel>0) printf("%d ",pdgdau);
    Float_t ptdau=TMath::Sqrt(o->Px()*o->Px()+o->Py()*o->Py());      
    Float_t etadau=o->Eta();
    if(pdgdau==211){ 
      nPions++;
      sumPx+=o->Px();
      sumPy+=o->Py();
      sumPz+=o->Pz();
    }
    if(pdgdau==321){ 
      nKaons++;
      sumPx+=o->Px();
      sumPy+=o->Py();
      sumPz+=o->Pz();
    }
    if(TMath::Abs(etadau)<fEtaMaxDau && ptdau>fPtMinDau){
      if(pdgdau==211) nPionsInAcc++;
      if(pdgdau==321) nKaonsInAcc++;
    }
  }
  if(fDebugLevel>0) printf("\n");
  if(TMath::Abs(sumPx-dmes->Px())>0.001 ||
     TMath::Abs(sumPy-dmes->Py())>0.001 ||
     TMath::Abs(sumPz-dmes->Pz())>0.001){
    printf("Momentum conservation violation\n");
    return kFALSE;
  }
  return kTRUE;
}


