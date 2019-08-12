#include "TFile.h"
#include "TCanvas.h"
#include "THnBase.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TStyle.h"
#include "AliUEHistograms.h"
#include "AliCFContainer.h"
#include "env.h"
#include "utils.h"

#include <iostream>
using namespace std;


void calc_sum_of_ratios(
			TString inFileName = "corr/corr_mu_tkl_U0S_bcde.root",
			Int_t anaType =  kTrkTrk,
			const Int_t gStudySystematic = 0,//check env.h for the meaning
			const Int_t StepPhiCorr = 0 //fine binning for AMPT all
			){
  gStyle->SetOptStat(0);

  Int_t na = 0; Double_t* binsa;
  Int_t nt = 0; Double_t* binst;

  switch (anaType){
  case kTrkTrk: 
    na = ntbins_trktrk[kTrackPtAs]; binsa = bins_trk;
    nt = ntbins_trktrk[kTrackPtTr]; binst = bins_trk;
    break;
  case kTrkTrkITS:
    na = ntbins_trktrkits[kTrackPtAs]; binsa = bins_trkits;
    nt = ntbins_trktrkits[kTrackPtTr]; binst = bins_trkits;
    break;
  case kTrkTrkGen: 
    na = ntbins_trktrk[kTrackPtAs]; binsa = bins_trk;
    nt = ntbins_trktrk[kTrackPtTr]; binst = bins_trk;
    break;
  case kTklTkl: 
    na = ntbins_tkltkl[kTrackPtAs]; binsa = bins_tkl;
    nt = ntbins_tkltkl[kTrackPtTr]; binst = bins_tkl;
    break;
  case kTklTklMC:
    na = ntbins_tkltkl[kTrackPtAs]; binsa = bins_tkl;
    nt = ntbins_tkltkl[kTrackPtTr]; binst = bins_tkl;
    break;
  case kTklTklGen:
    na = ntbins_tkltklGen[kTrackPtAs]; binsa = bins_tklGen;
    nt = ntbins_tkltklGen[kTrackPtTr]; binst = bins_tklGen;
    break;
  default: cout << "Invalid analysis type choice!" << endl; 
  }

  const char* strPtTrg = (anaType==kTklTkl || anaType==kTklTklMC) ? "#Delta#varphi^{trig} (mrad)"  : "p_{T}^{trig} (GeV/c)";
  const char* strPtAss = (anaType==kTklTkl  || anaType==kTklTklMC)          ? "#Delta#varphi^{assoc} (mrad)" : "p_{T}^{assoc} (GeV/c)";
  TString outFileName = TString(inFileName).ReplaceAll("corr_","dphi_");
  if( (anaType==kTklTklGen || anaType==kTrkTrkGen) && outFileName.Contains("AnalysisResults")){
    outFileName.ReplaceAll("AnalysisResults","dphi");//match other names
    outFileName.ReplaceAll("_zvtx","");//match other names  
    if(anaType==kTklTklGen){
      outFileName.ReplaceAll("CL1.root","tkl_tkl_CL1.root");
      outFileName.ReplaceAll("V0M.root","tkl_tkl_V0M.root");
    }
    if(anaType==kTrkTrkGen){
      outFileName.ReplaceAll("CL1.root","trk_trk_CL1.root");
      outFileName.ReplaceAll("V0M.root","trk_trk_V0M.root");
    }
  }
  if(gStudySystematic == k_vertex)outFileName.ReplaceAll(".root","_vertex02.root");
  if(gStudySystematic == k_mixed_norm)outFileName.ReplaceAll(".root","_mixed_norm.root");
  if(gStudySystematic == k_mixed_integrated)outFileName.ReplaceAll(".root","_mixed_integrated.root");
  if(gStudySystematic ==k_no_mixed)outFileName.ReplaceAll(".root","_no_mixed.root");
  
  THnBase* hTrackHistS; // same:  deta,ptassoc,pttrig,cent,dphi,zvtx
  THnBase* hTrackHistM; // mixed: deta,ptassoc,pttrig,cent,dphi,zvtx
  THnBase* hEventHistS; // same:  pttrig,cent,zvtx
  THnBase* hEventHistM; // mixed: pttrig,cent,zvtx
  TH3D*    hEventCount; // event count step, cent, run

  TFile* inputFile = new TFile(inFileName);
  TDirectoryFile *df=(TDirectoryFile*)inputFile->Get("PWG4_PhiCorrelations");
  if (df){
    TList* list = (TList*) df->Get("histosPhiCorrelations");
    AliUEHistograms* hHistS = (AliUEHistograms*) list->FindObject("AliUEHistogramsSame");
    AliUEHistograms* hHistM = (AliUEHistograms*) list->FindObject("AliUEHistogramsMixed");
    hTrackHistS = hHistS->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetGrid(StepPhiCorr)->GetGrid();
    hTrackHistM = hHistM->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetGrid(StepPhiCorr)->GetGrid();
    hEventHistS = hHistS->GetUEHist(2)->GetEventHist()->GetGrid(StepPhiCorr)->GetGrid();
    hEventHistM = hHistM->GetUEHist(2)->GetEventHist()->GetGrid(StepPhiCorr)->GetGrid();
    //    hEventCount = (TH3D*) list->FindObject("hEventCount");
  } 
  else {
    hTrackHistS = (THnBase*) inputFile->Get("hTrackHistS");
    hTrackHistM = (THnBase*) inputFile->Get("hTrackHistM");
    hEventHistS = (THnBase*) inputFile->Get("hEventHistS");
    hEventHistM = (THnBase*) inputFile->Get("hEventHistM");
    //    hEventCount = (TH3D*)    inputFile->Get("hEventCount");
  }
  Int_t nz = hTrackHistS->GetAxis(kTrackZvtx)->GetNbins();
  printf("zbins = %i\n",nz);

  Double_t binWidthEta = hTrackHistS->GetAxis(kTrackDeta)->GetBinWidth(1);
  Double_t detaMax = 1.;
  switch (anaType){
  case kTrkTrk:    detaMax = xtmax_trktrk[kTrackDeta]; break;
  case kTrkTrkITS: detaMax = xtmax_trktrkits[kTrackDeta]; break;
  case kTrkTrkGen: detaMax = xtmax_trktrk[kTrackDeta]; break;
  case kTklTkl:    detaMax = xtmax_tkltkl[kTrackDeta]; break;
  case kTklTklMC:  detaMax = xtmax_tkltkl[kTrackDeta]; break;
  case kTklTklGen: detaMax = xtmax_tkltklGen[kTrackDeta]; break;
  default: cout << "Invalid analysis type choice!" << endl; 
  }
  if(inFileName.Contains("fb89"))detaMax = 1.8;

  Double_t finiteBinCorrection = 1. - 1./detaMax*binWidthEta/2.;
  if(anaType==kTrkTrk) finiteBinCorrection = 1.;

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

    /////////////////////////////////////////
    //count the number of trigger particles
    //    TH1D* hTriggerStat = hEventHistS->Projection(kEventPtTr); hTriggerStat->SetName(Form("TriggerStat_%i",ic)); // is this right?

    //    TH1D* hNeV = hEventCount->ProjectionY(Form("hNeV_%i",ic),hEventCount->GetXaxis()->FindBin(cmin+0.001),hEventCount->GetXaxis()->FindBin(cmax-0.001));
    //    hEventCount->SetName(Form("hNeV_%i",ic));
    //    Printf("N muon in [%.1f-%.1f]: %.0f",cmin,cmax,hTriggerStat->Integral());
    //    // 7 = after selection on the muon
    //    Printf("N events after mu cut in [%.1f-%.1f]: %.0f",cmin,cmax,hNeV->GetBinContent(7));
    //    hTriggerStat->Scale(1./hNeV->GetBinContent(7));//after cut on muon
    //    if(cTriggerStat){
    //      cTriggerStat->cd(ic+1);
    //      hTriggerStat->DrawCopy();
    //    }
    //    hTriggerStat->Write();
    Printf("\033[1;31m -- Multiplicity class %d = %.2f-%.2f\033[m",ic,cmin,cmax);
    ///////////////////////////////////////////
    TCanvas* cCent = debug ? new TCanvas(Form("cCent_%i",ic),Form("cCent_%i",ic),1900,1000) : 0;
    if (cCent) cCent->Divide(nt,na);
    for (Int_t it=0;it<nt;it++) {
      //pt or dphi trigger loop
      Int_t ibintMin = (anaType==kTklTkl  || anaType==kTklTklMC) ? 0 : it; // Leonardo: check this!
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
        Int_t ibinaMin = (anaType==kTklTkl || anaType==kTklTklMC) ? 0 : ia; // Leonardo: check this!
        hTrackHistS->GetAxis(kTrackPtAs)->SetRangeUser(binsa[ibinaMin]+0.001,binsa[ia+1]-0.001);
        hTrackHistM->GetAxis(kTrackPtAs)->SetRangeUser(binsa[ibinaMin]+0.001,binsa[ia+1]-0.001);
        Double_t amin = hTrackHistS->GetAxis(kTrackPtAs)->GetBinLowEdge(hTrackHistS->GetAxis(kTrackPtAs)->GetFirst());
        Double_t amax = hTrackHistS->GetAxis(kTrackPtAs)->GetBinUpEdge (hTrackHistS->GetAxis(kTrackPtAs)->GetLast());
        TH2D* hPhiEtaSMsum=0;
        Double_t nTriggersS =0.;
        Double_t nTriggersM =0.;
        for (Int_t iz=1;iz<=nz;iz++){
          //          if(anaType==kTklTklGen && TMath::Abs(hTrackHistS->GetAxis(kTrackZvtx)->GetBinCenter(iz))> .1)continue;
          //z vtx loop
          if(gStudySystematic== k_vertex && TMath::Abs(hTrackHistS->GetAxis(kTrackZvtx)->GetBinCenter(iz))> 2. ){
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
          hPhiEtaS->SetTitle(Form("%.2f < %s < %.2f - %.2f< z(cm) <%.2f",amin,strPtAss,amax,zmin,zmax));
          hPhiEtaM->SetTitle(Form("%.2f < %s < %.2f - %.2f< z(cm) <%.2f",amin,strPtAss,amax,zmin,zmax));
          // draw same and mixed in each z bin
          if (cS)  { cS->cd(ia*nz+iz); if(hPhiEtaS->Integral() > 0) hPhiEtaS->DrawCopy(drawOption); }
          if (cM)  { cM->cd(ia*nz+iz); if(hPhiEtaM->Integral() > 0) hPhiEtaM->DrawCopy(drawOption); }
          Double_t norm = 1.;
          if(gStudySystematic == k_mixed_norm){//scale mixed to the integral
            norm=hPhiEtaM->Integral();
          }
          else{
	    if(anaType==kTrkTrk || anaType==kTrkTrkGen || anaType==kTrkTrkITS) {
              // scale mixed to 1 at eta=0 (integrated over phi)
              Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(-0.01);
              Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(+0.01);
              Int_t nNormBins = (binEta2-binEta1+1)*(hPhiEtaM->GetXaxis()->GetNbins());
              norm = hPhiEtaM->Integral(1,hPhiEtaM->GetXaxis()->GetNbins(),binEta1,binEta2)/nNormBins;
	    }
            else if(anaType==kTklTkl || anaType==kTklTklMC || anaType==kTklTklGen) {
              // scale mixed to 1 at (0,0)
              Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-0.01);
              Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(+0.01);
              Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(-0.01);
              Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(+0.01);
              Int_t nNormBins = (binEta2-binEta1+1)*(binPhi2-binPhi1+1);
              norm = hPhiEtaM->Integral(binPhi1,binPhi2,binEta1,binEta2)/nNormBins;
            }
            else {
              //scale mixed to the region of full acceptance
              Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin( -pi/2+0.0001);
              Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*pi/2-0.0001);
              Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin( -0.0001);
              Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin( +0.0001);
              Int_t nNormBins = (binEta2-binEta1+1)*(binPhi2-binPhi1+1);
              norm = hPhiEtaM->Integral(binPhi1,binPhi2,binEta1,binEta2)/nNormBins;
            }
          }
          if(norm!=0.)hPhiEtaM->Scale(1./norm*finiteBinCorrection);//enough stat in mixed events
          else hPhiEtaS->Reset("M"); //reset also the same evt. distr.
          if (debug) printf("nTriggersS(%i)=%.0f nTriggersM(%i)=%.0f norm=%.0f\n",iz,hTriggersS->GetBinContent(iz),iz,hTriggersM->GetBinContent(iz),norm);
          if (SaveMixed) hPhiEtaM->Write();
	  if (SaveSame)  hPhiEtaS->Write();
          // same/mixed
          TH2D* hPhiEtaSM = (TH2D*) hPhiEtaS->Clone(Form("hPhiEtaSM_%i_%i_%i_%i",ic,it,ia,iz));
          if (!hPhiEtaSMsum){
            hPhiEtaSMsum = (TH2D*) hPhiEtaS->Clone(Form("dphi_%i_%i_%i",it,ia,ic));
            hPhiEtaSMsum->Reset("M");
          }
          if(gStudySystematic!=k_no_mixed && norm!=0.)hPhiEtaSM->Divide(hPhiEtaM);

	  // some code to judge the relative contribution of the different correlation functions to the overall uncertainty
          Double_t errors[] = { 0, 0, 0 };
	  Double_t sums[] = { 0, 0, 0 }; //S M S/M
          // From AliUEHist
          //          for (Int_t x=1; x<=hPhiEtaS->GetNbinsX(); x++)
          //            for (Int_t y=1; y<=hPhiEtaS->GetNbinsY(); y++)
          //            {
          //              sums[0] += hPhiEtaS->GetBinContent(x, y);
          //              errors[0] += hPhiEtaS->GetBinError(x, y);
          //              sums[1] += hPhiEtaM->GetBinContent(x, y);
          //              errors[1] += hPhiEtaM->GetBinError(x, y);
          //            }
          //          for (Int_t x=1; x<=hPhiEtaSM->GetNbinsX(); x++)
          //            for (Int_t y=1; y<=hPhiEtaSM->GetNbinsY(); y++)
          //            {
          //              sums[2] += hPhiEtaSM->GetBinContent(x, y);
          //              errors[2] += hPhiEtaSM->GetBinError(x, y);
          //            }
          //          for (Int_t x=0; x<3; x++)
          //            if (sums[x] > 0)
          //              errors[x] /= sums[x];

          for (Int_t x=1; x<=hPhiEtaS->GetNbinsX(); x++)
            for (Int_t y=1; y<=hPhiEtaS->GetNbinsY(); y++)
	      {
		sums[0] += hPhiEtaS->GetBinContent(x, y);
		sums[1] += hPhiEtaM->GetBinContent(x, y);
		sums[2] += hPhiEtaSM->GetBinContent(x, y);
		if( hPhiEtaS->GetBinContent(x, y)!=0)errors[0] += TMath::Abs(hPhiEtaS->GetBinError(x, y)/hPhiEtaS->GetBinContent(x, y));
		if( hPhiEtaM->GetBinContent(x, y)!=0)errors[1] += TMath::Abs(hPhiEtaM->GetBinError(x, y)/hPhiEtaM->GetBinContent(x, y));
		if(hPhiEtaSM->GetBinContent(x, y)!=0)errors[2] += TMath::Abs(hPhiEtaSM->GetBinError(x, y)/hPhiEtaSM->GetBinContent(x, y));
		//              Printf("%f   %f   %f",errors[0],errors[1],errors[2]);
	      }
          for (Int_t x=0; x<3; x++)errors[x] /= (hPhiEtaS->GetNbinsX()*hPhiEtaS->GetNbinsY());
	  if(sums[0]<1.e-9 || sums[1]<1.e-9){
	    Printf("\033[1;31m skipping (ic=%d iz=%d) because histograms are empty (sums = %f  %f)\033[m", ic, iz,sums[0],sums[1]);
	    continue;
          }
	  
          if(debug)
            Printf("The correlation function (ic=%d iz=%d) has uncertainties S=%f M=%f S/M= %f (Ratio S/M %f)", ic, iz, errors[0], errors[1], errors[2], errors[0] / errors[1]);
          if(errors[0]/errors[1] < gCutStatErrorSM){
            Printf("\033[1;31m skipping (ic=%d iz=%d) because Ratio S/M=%f\033[m", ic, iz,errors[0]/errors[1]);
            continue;
          }

	  

          hPhiEtaSMsum->Add(hPhiEtaSM);
          // draw same/mixed
          if (cSM) {
            // scale same by the number of triggers in iz bin (for drawing purposes only)
            hPhiEtaSM->SetName("");
            if(hTriggersS->GetBinContent(iz) > 0) hPhiEtaSM->Scale(1./hTriggersS->GetBinContent(iz));
            cSM->cd(ia*nz+iz); if(hPhiEtaSM->Integral() > 0) hPhiEtaSM->DrawCopy(drawOption);
          }
        } // zvtx
        printf("nTriggersS=%.0f nTriggersM=%.0f\n",nTriggersS,nTriggersM);
        if(hPhiEtaSMsum){
          if(nTriggersS!=0)hPhiEtaSMsum->Scale(1./nTriggersS);
          hPhiEtaSMsum->Scale(1./hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
          hPhiEtaSMsum->Scale(1./hPhiEtaSMsum->GetYaxis()->GetBinWidth(1));
          hPhiEtaSMsum->SetTitle(Form("%.2f < %s < %.2f - %.2f < %s < %.2f - %.2f-%.2f",tmin,strPtTrg,tmax,amin,strPtAss,amax,cmin,cmax));
          printf("trigger = [%f,%f] associated = [%f,%f]\n",tmin,tmax,amin,amax);
          if(anaType==kTrkTrkGen) rotateDphi(hPhiEtaSMsum);
          hPhiEtaSMsum->Write();
          if (cCent) { cCent->cd(ia*nt+it+1); hPhiEtaSMsum->Draw(drawOption); }
        }
      } // pt-a
    } // pt-t
  } // centrality
  fout->Close();
  delete fout;
  if(inFileName.Contains("fb89")) Printf("Automatically changed detaMax to 1.8 for hybrid tracks!");
}
