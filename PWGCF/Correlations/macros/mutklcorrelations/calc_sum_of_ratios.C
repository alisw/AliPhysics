#include "TFile.h"
#include "TCanvas.h"
#include "THnBase.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TStyle.h"
#include "env.h"


Int_t debug=0; //plots for debugging
Int_t SaveMixed=0; //save mixed event distribution
Char_t drawOption[]="surf1";//"colz";
const Int_t nc = 6; 
//Double_t binscmin[] = { 0, 20, 40,  60,  70, 60};
//Double_t binscmax[] = {20, 40, 60, 100, 100, 80};
Double_t binscmin[] = { 100, 100,  100,  5,   0,   0};
Double_t binscmax[] = { 150, 200, 1000, 50,  50, 100};

void calc_sum_of_ratios(
    TString inFileName = "corr/corr_mu_tkl_U0S_bcde.root",
    const Bool_t isFB =  1, //central-central correlations
    const Bool_t isTkl = 1, //for analysis with tracklets
    const Int_t gStudySystematic = 0,//check env.h for the meaning
    const Bool_t isAMPTgen = 0 //fine binning for AMPT all
){
  gStyle->SetOptStat(0);

  Int_t na = 0; Double_t* binsa;
  Int_t nt = 0; Double_t* binst;

  if (isFB){
    if(isTkl){
      na = nbins_tkl;      binsa = bins_tkl;
      nt = nbins_muon_tkl; binst = bins_muon_tkl;
    }else{
      na = nbins_trk;      binsa = bins_trk;
      nt = nbins_muon_trk; binst = bins_muon_trk;
      if(isAMPTgen){
        nt = nbins_muon_ampt_gen;
        binst = bins_muon_ampt_gen;}
    }
  } else {
    if(isTkl){
      na = nbins_tkl; binsa = bins_tkl;
      nt = nbins_tkl; binst = bins_tkl;
    }else{
      na = nbins_trk; binsa = bins_trk;
      nt = nbins_trk; binst = bins_trk;
    }
  }

  const char* strPtTrg = (isTkl && !isFB) ? "#Delta#varphi^{trig} (mrad)"  : "p_{T}^{trig} (GeV/c)";
  const char* strPtAss = (isTkl)          ? "#Delta#varphi^{assoc} (mrad)" : "p_{T}^{assoc} (GeV/c)";

  TString outFileName = TString("dphi/dphi_")+TString(inFileName).ReplaceAll("corr/","");
  if(gStudySystematic == k_vertex05)outFileName.ReplaceAll(".root","_vertex05.root");
  if(gStudySystematic == k_vertex01)outFileName.ReplaceAll(".root","_vertex01.root");
  if(gStudySystematic == k_vertex02)outFileName.ReplaceAll(".root","_vertex02.root");
  if(gStudySystematic == k_vertex03)outFileName.ReplaceAll(".root","_vertex03.root");
  if(gStudySystematic == k_vertex04)outFileName.ReplaceAll(".root","_vertex04.root");
  if(gStudySystematic == k_mixed_norm)outFileName.ReplaceAll(".root","_mixed_norm.root");
  if(gStudySystematic == k_mixed_integrated)outFileName.ReplaceAll(".root","_mixed_integrated.root");
  if(gStudySystematic ==k_no_mixed)outFileName.ReplaceAll(".root","_no_mixed.root");

  THnBase* hTrackHistS; // same:  deta,ptassoc,pttrig,cent,dphi,zvtx
  THnBase* hTrackHistM; // mixed: deta,ptassoc,pttrig,cent,dphi,zvtx
  THnBase* hEventHistS; // same:  pttrig,cent,zvtx
  THnBase* hEventHistM; // mixed: pttrig,cent,zvtx
  TH3D*    hEventCount; // event count step, cent, run

  TFile* inputFile = new TFile(inFileName);
  if (0){
    //    TList* list = (TList*) inputFile->Get("histosPhiCorrelations");
    //    AliUEHistograms* hHistS = (AliUEHistograms*) list->FindObject("AliUEHistogramsSame");
    //    AliUEHistograms* hHistM = (AliUEHistograms*) list->FindObject("AliUEHistogramsMixed");
    //    hTrackHistS = hHistS->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetGrid(6)->GetGrid();
    //    hTrackHistM = hHistM->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetGrid(6)->GetGrid();
    //    hEventHistS = hHistS->GetUEHist(2)->GetEventHist()->GetGrid(6)->GetGrid();
    //    hEventHistM = hHistM->GetUEHist(2)->GetEventHist()->GetGrid(6)->GetGrid();
    //    hEventCount = (TH3D*) list->FindObject("hEventCount");
  } else {
    hTrackHistS = (THnBase*) inputFile->Get("hTrackHistS");
    hTrackHistM = (THnBase*) inputFile->Get("hTrackHistM");
    hEventHistS = (THnBase*) inputFile->Get("hEventHistS");
    hEventHistM = (THnBase*) inputFile->Get("hEventHistM");
    hEventCount = (TH3D*)    inputFile->Get("hEventCount");
  }
  Int_t nz = hTrackHistS->GetAxis(kTrackZvtx)->GetNbins();
  printf("zbins = %i\n",nz);

  Double_t binWidthEta = hTrackHistS->GetAxis(kTrackDeta)->GetBinWidth(1);
  Double_t finiteBinCorrection = 1;//1 - 1./detaMax*binWidthEta/2;

  TFile* fout = new TFile(outFileName,"recreate");

  TCanvas* cTriggerStat = debug ? new TCanvas(Form("cTriggerStat"),Form("cTriggerStat"),1400,500) : 0;
  if(cTriggerStat)cTriggerStat->Divide(nc,1);

  for (Int_t ic=0;ic<nc;ic++){
    //centrality selection
    hTrackHistS->GetAxis(kTrackCent)->SetRangeUser(binscmin[ic]+0.001,binscmax[ic]-0.001);
    hEventHistS->GetAxis(kEventCent)->SetRangeUser(binscmin[ic]+0.001,binscmax[ic]-0.001);
    //only set the range for mixed event if not integrated
    if(gStudySystematic!=k_mixed_integrated)hTrackHistM->GetAxis(kTrackCent)->SetRangeUser(binscmin[ic]+0.001,binscmax[ic]-0.001);
    if(gStudySystematic!=k_mixed_integrated)hEventHistM->GetAxis(kEventCent)->SetRangeUser(binscmin[ic]+0.001,binscmax[ic]-0.001);

    Double_t cmin = hTrackHistS->GetAxis(kTrackCent)->GetBinLowEdge(hTrackHistS->GetAxis(kTrackCent)->GetFirst());
    Double_t cmax = hTrackHistS->GetAxis(kTrackCent)->GetBinUpEdge (hTrackHistS->GetAxis(kTrackCent)->GetLast());
    hEventHistS->GetAxis(kEventZvtx)->SetRange(1,hEventHistS->GetAxis(kEventZvtx)->GetNbins());//reset z bin
    //count the number of trigger particles
    TH1D* hTriggerStat = hEventHistS->Projection(kEventPtTr); hTriggerStat->SetName(Form("TriggerStat_%i",ic));
    TH1D* hNeV = hEventCount->ProjectionY(Form("hNeV_%i",ic),hEventCount->GetXaxis()->FindBin(cmin+0.001),hEventCount->GetXaxis()->FindBin(cmax-0.001));
    hEventCount->SetName(Form("hNeV_%i",ic));
    Printf("N muon in [%.1f-%.1f]: %.0f",cmin,cmax,hTriggerStat->Integral());
    // 7 = after selection on the muon
    Printf("N events after mu cut in [%.1f-%.1f]: %.0f",cmin,cmax,hNeV->GetBinContent(7));
    hTriggerStat->Scale(1./hNeV->GetBinContent(7));//after cut on muon
    if(cTriggerStat){
      cTriggerStat->cd(ic+1);
      hTriggerStat->DrawCopy();
    }
    hTriggerStat->Write();

    Printf("\033[1;31m -- Multiplicity class %d = %.0f-%.0f\033[m",ic,cmin,cmax);
    TCanvas* cCent = debug ? new TCanvas(Form("cCent_%i",ic),Form("cCent_%i",ic),1900,1000) : 0;
    if (cCent) cCent->Divide(nt,na);
    for (Int_t it=0;it<nt;it++) {
      //pt or dphi trigger loop
      Int_t ibintMin = (!isFB && isTkl) ? 0 : it;
      hTrackHistS->GetAxis(kTrackPtTr)->SetRangeUser(binst[ibintMin]+0.001,binst[it+1]-0.001);
      hTrackHistM->GetAxis(kTrackPtTr)->SetRangeUser(binst[ibintMin]+0.001,binst[it+1]-0.001);
      hEventHistS->GetAxis(kEventPtTr)->SetRangeUser(binst[ibintMin]+0.001,binst[it+1]-0.001);
      hEventHistM->GetAxis(kEventPtTr)->SetRangeUser(binst[ibintMin]+0.001,binst[it+1]-0.001);
      Double_t tmin = hTrackHistS->GetAxis(kTrackPtTr)->GetBinLowEdge(hTrackHistS->GetAxis(kTrackPtTr)->GetFirst());
      Double_t tmax = hTrackHistS->GetAxis(kTrackPtTr)->GetBinUpEdge (hTrackHistS->GetAxis(kTrackPtTr)->GetLast());
      TCanvas* cS  =  debug ? new TCanvas(Form("cS_%i_%i",ic,it),Form("cS_%i_%i",ic,it),1900,1000) : 0;
      TCanvas* cM  =  debug ? new TCanvas(Form("cM_%i_%i",ic,it),Form("cM_%i_%i",ic,it),1900,1000) : 0;
      TCanvas* cSM =  debug ? new TCanvas(Form("cSM_%i_%i",ic,it),Form("cSM_%i_%i",ic,it),1900,1000) : 0;
      TH1D* hTriggersS = hEventHistS->Projection(kEventZvtx); hTriggersS->SetName(Form("TriggersS_%i_%i",ic,it));
      TH1D* hTriggersM = hEventHistM->Projection(kEventZvtx); hTriggersM->SetName(Form("TriggersM_%i_%i",ic,it));
      if (cS)  cS->Divide(nz,na,0.001,0.001);
      if (cM)  cM->Divide(nz,na,0.001,0.001);
      if (cSM) cSM->Divide(nz,na,0.001,0.001);
      for (Int_t ia=0;ia<na;ia++) {
        //pt or dphi assoc loop
        Int_t ibinaMin = isTkl ? 0 : ia;
        hTrackHistS->GetAxis(kTrackPtAs)->SetRangeUser(binsa[ibinaMin]+0.001,binsa[ia+1]-0.001);
        hTrackHistM->GetAxis(kTrackPtAs)->SetRangeUser(binsa[ibinaMin]+0.001,binsa[ia+1]-0.001);
        Double_t amin = hTrackHistS->GetAxis(kTrackPtAs)->GetBinLowEdge(hTrackHistS->GetAxis(kTrackPtAs)->GetFirst());
        Double_t amax = hTrackHistS->GetAxis(kTrackPtAs)->GetBinUpEdge (hTrackHistS->GetAxis(kTrackPtAs)->GetLast());
        TH2D* hPhiEtaSMsum=0;
        Double_t nTriggersS =0.;
        Double_t nTriggersM =0.;
        for (Int_t iz=1;iz<=nz;iz++){
          //z vtx loop
          if(gStudySystematic== k_vertex05 && TMath::Abs(hTrackHistS->GetAxis(kTrackZvtx)->GetBinCenter(iz))> 5. ){
            Printf("-----gStudySystematic=%d skipping zvtx bin %d",gStudySystematic,iz);
            continue;
          }
          if(gStudySystematic== k_vertex01 && TMath::Abs(hTrackHistS->GetAxis(kTrackZvtx)->GetBinCenter(iz))> 1. ){
            Printf("-----gStudySystematic=%d skipping zvtx bin %d",gStudySystematic,iz);
            continue;
          }

          if(gStudySystematic== k_vertex02 && TMath::Abs(hTrackHistS->GetAxis(kTrackZvtx)->GetBinCenter(iz))> 2. ){
            Printf("-----gStudySystematic=%d skipping zvtx bin %d",gStudySystematic,iz);
            continue;
          }

          if(gStudySystematic== k_vertex03 && TMath::Abs(hTrackHistS->GetAxis(kTrackZvtx)->GetBinCenter(iz))> 3. ){
            Printf("-----gStudySystematic=%d skipping zvtx bin %d",gStudySystematic,iz);
            continue;
          }

          if(gStudySystematic== k_vertex04 && TMath::Abs(hTrackHistS->GetAxis(kTrackZvtx)->GetBinCenter(iz))> 4. ){
            Printf("-----gStudySystematic=%d skipping zvtx bin %d",gStudySystematic,iz);
            continue;
          }

          nTriggersS+=hTriggersS->Integral(iz,iz);
          nTriggersM+=hTriggersM->Integral(iz,iz);
          hTrackHistS->GetAxis(kTrackZvtx)->SetRange(iz,iz);
          hTrackHistM->GetAxis(kTrackZvtx)->SetRange(iz,iz);
          Double_t zmin = hTrackHistS->GetAxis(kTrackZvtx)->GetBinLowEdge(iz);
          Double_t zmax = hTrackHistS->GetAxis(kTrackZvtx)->GetBinUpEdge(iz);
          TH2D* hPhiEtaS  = hTrackHistS->Projection(kTrackDeta,kTrackDphi,"E");  hPhiEtaS->SetName(Form("hPhiEtaS_%i_%i_%i_%i",ic,it,ia,iz));
          TH2D* hPhiEtaM  = hTrackHistM->Projection(kTrackDeta,kTrackDphi,"E");  hPhiEtaM->SetName(Form("hPhiEtaM_%i_%i_%i_%i",ic,it,ia,iz));
          hPhiEtaS->SetTitle(Form("%.1f < %s < %.1f - %.1f< z(cm) <%.1f",amin,strPtAss,amax,zmin,zmax));
          hPhiEtaM->SetTitle(Form("%.1f < %s < %.1f - %.1f< z(cm) <%.1f",amin,strPtAss,amax,zmin,zmax));
          // draw same and mixed in each z bin
          if (cS)  { cS->cd(ia*nz+iz); hPhiEtaS->DrawCopy(drawOption); }
          if (cM)  { cM->cd(ia*nz+iz); hPhiEtaM->DrawCopy(drawOption); }
          // scale mixed to 1 at (0,0)
          //          Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-0.01);
          //          Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(+0.01);
          //          Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(-0.01);
          //          Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(+0.01);
          //          Int_t nNormBins = (binEta2-binEta1+1)*(binPhi2-binPhi1+1);
          //          Double_t norm = hPhiEtaM->Integral(binPhi1,binPhi2,binEta1,binEta2)/nNormBins;
          Double_t norm = 1.;
          if(gStudySystematic == k_mixed_norm){//scale mixed to the integral
            norm=hPhiEtaM->Integral();
          }else{//scale mixed to the region of full acceptance
            Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin( -pi/2+0.0001);
            Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*pi/2-0.0001);
            Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin( isFB ? -3.5+0.0001 : -0.0001);
            Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin( isFB ? -3.0-0.0001 : +0.0001);
            Int_t nNormBins = (binEta2-binEta1+1)*(binPhi2-binPhi1+1);
            norm = hPhiEtaM->Integral(binPhi1,binPhi2,binEta1,binEta2)/nNormBins;
          }
          hPhiEtaM->Scale(1./norm/finiteBinCorrection);
          if (debug) printf("nTriggersS(%i)=%.0f nTriggersM(%i)=%.0f norm=%.0f\n",iz,hTriggersS->GetBinContent(iz),iz,hTriggersM->GetBinContent(iz),norm);
          if (SaveMixed) hPhiEtaM->Write();
          // same/mixed
          TH2D* hPhiEtaSM = (TH2D*) hPhiEtaS->Clone(Form("hPhiEtaSM_%i_%i_%i_%i",ic,it,ia,iz));
          if(gStudySystematic!=k_no_mixed)hPhiEtaSM->Divide(hPhiEtaM);
          if (!hPhiEtaSMsum) hPhiEtaSMsum = (TH2D*) hPhiEtaSM->Clone(Form("dphi_%i_%i_%i",it,ia,ic)); else hPhiEtaSMsum->Add(hPhiEtaSM);
          // draw same/mixed
          if (cSM) {
            // scale same by the number of triggers in iz bin (for drawing purposes only)
            hPhiEtaSM->SetName("");
            hPhiEtaSM->Scale(1./hTriggersS->GetBinContent(iz));
            cSM->cd(ia*nz+iz); hPhiEtaSM->DrawCopy(drawOption);
          }
        } // zvtx
        printf("nTriggersS=%.0f nTriggersM=%.0f\n",nTriggersS,nTriggersM);
        hPhiEtaSMsum->Scale(1./nTriggersS);
        hPhiEtaSMsum->Scale(1./hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum->SetTitle(Form("%.1f < %s < %.1f - %.1f < %s < %.1f - %.0f-%.0f%%",tmin,strPtTrg,tmax,amin,strPtAss,amax,cmin,cmax));
        printf("tmax = %f\n",tmax);
        hPhiEtaSMsum->Write();
        if (cCent) { cCent->cd(ia*nt+it+1); hPhiEtaSMsum->Draw(drawOption); }
      } // pt-a
    } // pt-t
  } // centrality
  fout->Close();
  delete fout;
}
