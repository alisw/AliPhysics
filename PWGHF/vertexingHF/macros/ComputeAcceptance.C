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

enum EDDecay{kD0Kpi,kDplusKpipi,kDstarD0pi,kDsKKpi,kLcpKpi,kLcK0Sp};
enum EFidY{kFixedY,kPtDepY};
enum EPtShape{kFlat,kFONLL8TeV,kFONLL8TeVfeeddown,kFONLL7TeV,kPythia7TeV,kFONLL5TeV,kFONLL13TeVprompt,kPythia13TeVprompt};

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
Bool_t CountPKpi(TClonesArray *array, Int_t nentries, Int_t &nPions, Int_t &nKaons, Int_t &nProtons, Int_t &nPionsInAcc, Int_t &nKaonsInAcc, Int_t &nProtonsInAcc);


// Pt-shape histograms
TH1D* LoadFONLL13TeV_promptD0();
TH1D* LoadFONLL13TeV_promptDplus();
TH1D* LoadFONLL13TeV_promptDstar();
TH1D* LoadPYTHIA13TeV_promptD0();
TH1D* LoadPYTHIA13TeV_promptDplus();
TH1D* LoadPYTHIA13TeV_promptDstar();
TH1D* LoadPYTHIA13TeV_promptDs();



void ComputeAcceptance(){
  // main function
  
  gSystem->Load("liblhapdf.so");      // Parton density functions
  gSystem->Load("libEGPythia6.so");   // TGenerator interface
  gSystem->Load("libpythia6.so");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations

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
  Int_t nProtonDau=-1;
  Int_t nKaonDau=-1;
  TString outFileName="Acceptance_Toy_";
  if(fDDecay==kD0Kpi){
    pdgCode=421;
    nPionDau=1;
    nKaonDau=1;
    nProtonDau=0;
    outFileName.Append("D0Kpi_");
  }else if(fDDecay==kDplusKpipi){
    pdgCode=411;
    nPionDau=2;
    nKaonDau=1;
    nProtonDau=0;
    outFileName.Append("DplusKpipi_");
  }else if(fDDecay==kDstarD0pi){
    pdgCode=413;
    nPionDau=2;
    nKaonDau=1;
    nProtonDau=0;
    outFileName.Append("DStarD0pi_");
  }else if(fDDecay==kDsKKpi){
    pdgCode=431;
    nPionDau=1;
    nKaonDau=2;
    nProtonDau=0;
    outFileName.Append("DsKKpi_");
  }else if(fDDecay==kLcpKpi){
    pdgCode=4122;
    nPionDau=1;
    nKaonDau=1;
    nProtonDau=1;
    outFileName.Append("LcpKpi_");
  }else if(fDDecay==kLcK0Sp){
    pdgCode=4122;
    nPionDau=2;
    nKaonDau=0;
    nProtonDau=1;
    outFileName.Append("LcK0Sp_");
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

  TF1*  funcPt=0x0;
  TH1D* histPt=0x0;
  if(fPtShape==kFONLL8TeV){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.518046,3.01138,3.38914,1.75899); // Prompt
    outFileName.Append("FONLL8ptshape.root");
  }else if(fPtShape==kFONLL8TeVfeeddown){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.398252, 3.9603, 3.915, 1.51853); // FeedDown
    outFileName.Append("FONLL8ptshapeFeedDown.root");
  }else if(fPtShape==kFONLL7TeV){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.322643,2.96275,2.30301,2.5);
    outFileName.Append("FONLL7ptshape.root");
  }else if(fPtShape==kFONLL5TeV){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.302879,2.9750,3.68139,1.68855);
    outFileName.Append("FONLL5ptshape.root");
  }else if(fPtShape==kPythia7TeV){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.322643,1.94635,1.40463,2.5);
    outFileName.Append("PYTHIA7ptshape.root");  
  }else if(fPtShape==kFONLL13TeVprompt){
    if(fDDecay==kDplusKpipi){
      histPt = LoadFONLL13TeV_promptDplus();
      outFileName.Append("promptDplus");
    }else if (fDDecay==kDstarD0pi){
      histPt = LoadFONLL13TeV_promptDstar();
      outFileName.Append("promptDstar");
    }else{
      histPt = LoadFONLL13TeV_promptD0();
      outFileName.Append("promptD0");
    }
    outFileName.Append("FONLL13ptshape.root");
  }else if(fPtShape==kPythia13TeVprompt){
    if(fDDecay==kDplusKpipi){
      histPt = LoadPYTHIA13TeV_promptDplus();
      outFileName.Append("promptDplus");
    }else if (fDDecay==kDstarD0pi){
      histPt = LoadPYTHIA13TeV_promptDstar();
      outFileName.Append("promptDstar");
    }else if (fDDecay==kDsKKpi){
      histPt = LoadPYTHIA13TeV_promptDs();
      outFileName.Append("promptDs");
    }else{
      histPt = LoadPYTHIA13TeV_promptD0();
      outFileName.Append("promptD0");
    }
    outFileName.Append("PYTHIA13ptshape.root");
  }else{
    funcPt=new TF1("fFlat","pol0",0.,40.);
    funcPt->SetParameter(0,1.);
    outFileName.Append("flatpt.root");
  }

  if (funcPt) funcPt->SetNpx(10000);

  TRandom3* gener=new TRandom3(0);
  TLorentzVector* vec=new TLorentzVector();


  for(Int_t itry=0; itry<totTrials; itry++){
    if(itry%10000==0) printf("Event %d\n",itry);
    Float_t ptD = funcPt ? funcPt->GetRandom() : histPt->GetRandom();
    Float_t phiD=gener->Rndm()*2*TMath::Pi();
    Float_t yD=gener->Rndm()*2.-1.; // flat in -1<y<1
    Float_t px=ptD*TMath::Cos(phiD);
    Float_t py=ptD*TMath::Sin(phiD);
    Float_t mt=TMath::Sqrt(massD*massD+ptD*ptD);
    Float_t pz=mt*TMath::SinH(yD);
    Float_t E=TMath::Sqrt(massD*massD+px*px+py*py+pz*pz);

    // TLorentzVector* vec=new TLorentzVector(px,py,pz,E);
    vec->SetPxPyPzE(px,py,pz,E);
    pdec->Decay(pdgCode,vec);
    array->Clear();
    Int_t nentries = pdec->ImportParticles(array);
    TParticle* dmes=(TParticle*)array->At(0);
    Int_t nDaughters=dmes->GetNDaughters();
    if(fDDecay==kD0Kpi && nDaughters!=2) continue;
    if(fDDecay==kLcK0Sp && nentries>6) continue;
    Int_t nPionsInAcc=0;
    Int_t nProtonsInAcc=0;
    Int_t nKaonsInAcc=0;
    Int_t nPions=0;
    Int_t nProtons=0;
    Int_t nKaons=0;
    Bool_t isOk=CountPKpi(array,nentries,nPions,nKaons,nProtons,nPionsInAcc,nKaonsInAcc,nProtonsInAcc);

    if(isOk){
      if(nPions==nPionDau && nKaons==nKaonDau && nProtons==nProtonDau){
	hPtVsYGen->Fill(ptD,yD);
	if(TMath::Abs(yD)<0.5){
	  hPtVsYGenLimAcc->Fill(ptD,yD);
	}
	if(IsInFiducialAcceptance(ptD,yD)){	  
	  if(nPionsInAcc==nPionDau && nKaonsInAcc==nKaonDau && nProtonsInAcc==nProtonDau){ 
	    hPtVsYGenAcc->Fill(ptD,yD);
	  }
	}
      }
    }
    // delete vec;
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
    if(pdgdau==130) {
      if(fDebugLevel>0) printf("K0 dacaying into K0L\n");
      return kFALSE;
    }
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


//___________________________________________________
Bool_t CountPKpi(TClonesArray *array, Int_t nentries, Int_t &nPions, Int_t &nKaons, Int_t &nProtons, Int_t &nPionsInAcc, Int_t &nKaonsInAcc, Int_t &nProtonsInAcc){
  // count K and pi in Acc

  TParticle* dmes=(TParticle*)array->At(0);
  Double_t sumPx=0;
  Double_t sumPy=0;
  Double_t sumPz=0;
  
  for(int j=0; j<nentries; j++){
    TParticle * o = (TParticle*)array->At(j);
    Int_t pdgdau=TMath::Abs(o->GetPdgCode());
    if(fDebugLevel>0) printf("%d ",pdgdau);
    if(pdgdau==130) {
      if(fDebugLevel>0) printf("K0 dacaying into K0L\n");
      return kFALSE;
    }
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
    if(pdgdau==2212){ 
      nProtons++;
      sumPx+=o->Px();
      sumPy+=o->Py();
      sumPz+=o->Pz();
    }
    if(TMath::Abs(etadau)<fEtaMaxDau && ptdau>fPtMinDau){
      if(pdgdau==211) nPionsInAcc++;
      if(pdgdau==321) nKaonsInAcc++;
      if(pdgdau==2212) nProtonsInAcc++;
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



//___________________________________________________
TH1D* LoadFONLL13TeV_promptD0()
{
  TH1D *hFONLL13 = new TH1D("hFONLL13TeV_D0", "", 80, 0., 40.);
  Float_t val[80] = {
    1.4686e+08, 3.9542e+08, 5.4901e+08, 5.2166e+08, 4.1083e+08, 2.9968e+08, 2.1299e+08, 1.5057e+08, 1.0701e+08, 7.6919e+07,
    5.6121e+07, 4.1546e+07, 3.1184e+07, 2.3715e+07, 1.8253e+07, 1.4206e+07, 1.1175e+07, 8.8774e+06, 7.1169e+06, 5.7544e+06,
    4.6899e+06, 3.8509e+06, 3.1841e+06, 2.6499e+06, 2.2189e+06, 1.8687e+06, 1.5823e+06, 1.3467e+06, 1.1516e+06, 9.8933e+05,
    8.5356e+05, 7.3942e+05, 6.4302e+05, 5.6122e+05, 4.9154e+05, 4.3193e+05, 3.8074e+05, 3.3662e+05, 2.9846e+05, 2.6535e+05,
    2.3653e+05, 2.1136e+05, 1.8932e+05, 1.6997e+05, 1.5293e+05, 1.3789e+05, 1.2458e+05, 1.1277e+05, 1.0227e+05, 9.2913e+04,
    8.4560e+04, 7.7087e+04, 7.0388e+04, 6.4371e+04, 5.8956e+04, 5.4074e+04, 4.9667e+04, 4.5680e+04, 4.2068e+04, 3.8791e+04,
    3.5813e+04, 3.3103e+04, 3.0633e+04, 2.8380e+04, 2.6321e+04, 2.4436e+04, 2.2710e+04, 2.1127e+04, 1.9673e+04, 1.8336e+04,
    1.7106e+04, 1.5972e+04, 1.4926e+04, 1.3961e+04, 1.3069e+04, 1.2243e+04, 1.1479e+04, 1.0771e+04, 1.0113e+04, 9.5031e+03
  };
  for (Int_t ibin=0; ibin<80; ++ibin) hFONLL13->SetBinContent(ibin+1, val[ibin]);

  return hFONLL13;
}



//___________________________________________________
TH1D* LoadFONLL13TeV_promptDplus()
{
  TH1D *hFONLL13 = new TH1D("hFONLL13TeV_Dplus", "", 80, 0., 40.);
  Float_t val[80] = {
    1.5242e+08, 3.9396e+08, 5.3767e+08, 5.1263e+08, 4.0706e+08, 2.9926e+08, 2.1414e+08, 1.5230e+08, 1.0879e+08, 7.8507e+07,
    5.7460e+07, 4.2651e+07, 3.2093e+07, 2.4462e+07, 1.8865e+07, 1.4708e+07, 1.1588e+07, 9.2185e+06, 7.3998e+06, 5.9900e+06,
    4.8871e+06, 4.0167e+06, 3.3240e+06, 2.7686e+06, 2.3200e+06, 1.9552e+06, 1.6566e+06, 1.4107e+06, 1.2070e+06, 1.0374e+06,
    8.9550e+05, 7.7610e+05, 6.7519e+05, 5.8953e+05, 5.1652e+05, 4.5404e+05, 4.0037e+05, 3.5409e+05, 3.1405e+05, 2.7929e+05,
    2.4902e+05, 2.2258e+05, 1.9943e+05, 1.7909e+05, 1.6117e+05, 1.4535e+05, 1.3135e+05, 1.1892e+05, 1.0787e+05, 9.8024e+04,
    8.9230e+04, 8.1359e+04, 7.4302e+04, 6.7963e+04, 6.2256e+04, 5.7112e+04, 5.2465e+04, 4.8261e+04, 4.4452e+04, 4.0996e+04,
    3.7854e+04, 3.4995e+04, 3.2390e+04, 3.0012e+04, 2.7838e+04, 2.5849e+04, 2.4027e+04, 2.2355e+04, 2.0820e+04, 1.9408e+04,
    1.8108e+04, 1.6910e+04, 1.5805e+04, 1.4785e+04, 1.3842e+04, 1.2969e+04, 1.2161e+04, 1.1411e+04, 1.0716e+04, 1.0070e+04
  };
  for (Int_t ibin=0; ibin<80; ++ibin) hFONLL13->SetBinContent(ibin+1, val[ibin]);

  return hFONLL13;
}



//___________________________________________________
TH1D* LoadFONLL13TeV_promptDstar()
{
  TH1D *hFONLL13 = new TH1D("hFONLL13TeV_Dstar", "", 80, 0., 40.);
  Float_t val[80] = {
    1.2433e+08, 3.4512e+08, 5.0662e+08, 5.1020e+08, 4.2016e+08, 3.1661e+08, 2.3064e+08, 1.6632e+08, 1.2007e+08, 8.7340e+07,
    6.4329e+07, 4.8000e+07, 3.6287e+07, 2.7772e+07, 2.1492e+07, 1.6808e+07, 1.3278e+07, 1.0589e+07, 8.5178e+06, 6.9084e+06,
    5.6462e+06, 4.6479e+06, 3.8519e+06, 3.2124e+06, 2.6951e+06, 2.2738e+06, 1.9285e+06, 1.6437e+06, 1.4076e+06, 1.2108e+06,
    1.0459e+06, 9.0711e+05, 7.8968e+05, 6.8993e+05, 6.0484e+05, 5.3197e+05, 4.6933e+05, 4.1529e+05, 3.6850e+05, 3.2786e+05,
    2.9246e+05, 2.6152e+05, 2.3441e+05, 2.1058e+05, 1.8958e+05, 1.7104e+05, 1.5461e+05, 1.4003e+05, 1.2706e+05, 1.1550e+05,
    1.0517e+05, 9.5920e+04, 8.7626e+04, 8.0172e+04, 7.3461e+04, 6.7409e+04, 6.1940e+04, 5.6992e+04, 5.2508e+04, 4.8437e+04,
    4.4736e+04, 4.1367e+04, 3.8296e+04, 3.5492e+04, 3.2929e+04, 3.0583e+04, 2.8433e+04, 2.6460e+04, 2.4648e+04, 2.2981e+04,
    2.1446e+04, 2.0032e+04, 1.8727e+04, 1.7522e+04, 1.6407e+04, 1.5376e+04, 1.4420e+04, 1.3535e+04, 1.2713e+04, 1.1950e+04
  };
  for (Int_t ibin=0; ibin<80; ++ibin) hFONLL13->SetBinContent(ibin+1, val[ibin]);

  return hFONLL13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_promptD0()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13TeV_D0", "", 40, 0., 40.);
  Float_t val[40] = {
    2617743, 3763836, 2235903, 1140807, 587707, 317881, 180756, 107762, 67299,42882,
    28550, 19690, 13817, 9754, 7205, 5299, 4100, 3177, 2386,1822,
    1484, 1248, 933, 798, 639, 534, 435, 405, 304,260,
    210, 196, 161, 143, 139, 99, 89, 80, 66,63
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_promptDplus()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13TeV_Dplus", "", 40, 0., 40.);
  Float_t val[40] = {
    1192379, 1736283, 1047485, 540499, 281523, 152153, 87150, 52200, 32709, 20991,
    13986, 9611, 6832, 4705, 3484, 2650, 2001, 1518, 1197, 924,
    735, 576, 497, 402, 292, 261, 214, 173, 169, 150,
    106, 102, 95, 70, 67, 47, 53, 46, 35, 35
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_promptDstar()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13TeV_Dstar", "", 40, 0., 40.);
  Float_t val[40] = {
    922985, 1419749, 904587, 484407, 257862, 142386, 82177, 49839, 31103, 20422,
    13493, 9272, 6549, 4746, 3344, 2575, 1959, 1493, 1161, 935,
    672, 585, 456, 415, 328, 260, 199, 180, 161, 134,
    109, 97, 103, 72, 63, 56, 40, 46, 38, 27
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_promptDs()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13TeV_Ds", "", 40, 0., 40.);
  Float_t val[40] = {
    346381, 519143, 320488, 167944, 87484, 47325, 26932, 16376, 10058, 6527,
    4347, 3041, 2112, 1521, 1069, 849, 609, 468, 359, 275,
    226, 197, 152, 124, 111, 90, 63, 60, 52, 46,
    39, 35, 26, 29, 16, 20, 10, 13, 18, 9
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}

