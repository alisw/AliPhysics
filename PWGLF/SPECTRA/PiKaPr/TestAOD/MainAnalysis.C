#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <iostream>
#include "TFile.h"
#include "TH1I.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"

TString Charge[]={"Pos","Neg"};
TString Sign[]={"Plus","Minus"};
TString Particle[]={"Pion","Kaon","Proton"};
Int_t Color[3]={1,2,4};
Int_t Marker[6]={20,21,22,24,25,26};
Double_t projl[3]={0,20,80};
Double_t proju[3]={10,40,90};


void MainAnalysis() {

  //TString fold="5SigmaPID_AOD086";
  TString fold="5SigmaPID_AOD046";
  Int_t ibinToCompare=-1;
  
  TString sname="Cent0to5_QVec0.0to100.0";ibinToCompare=0;
  
  //TString sname="Cent5to10_QVec0.0to100.0";ibinToCompare=1;
  //TString sname="Cent10to20_QVec0.0to100.0";ibinToCompare=2;
  //TString sname="Cent20to30_QVec0.0to100.0";ibinToCompare=3;
  //TString sname="Cent30to40_QVec0.0to100.0";ibinToCompare=4;
  
  TString dataFile = Form("output/%s/Pt.AOD.1._data_ptcut_%s.root",fold.Data(),sname.Data());
  TString mcFile =Form("output/%s/Pt.AOD.1._MC_%s.root",fold.Data(),sname.Data());
  
  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libProof");
  gSystem->Load("libMatrix");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  gSystem->Load("libCORRFW");
  //gSystem->Load("libPWG0base");
  gSystem->Load("libMinuit");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGLFSPECTRA");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gStyle->SetPalette(1);
  
  // Open root MC file and get classes
  cout << "Analysis Macro" << endl;
  cout << "  > Reading MC data" << endl;
  TFile *_mc = TFile::Open(mcFile.Data());
  AliSpectraAODHistoManager* hman_mc = (AliSpectraAODHistoManager*) _mc->Get("SpectraHistos");
  AliSpectraAODEventCuts* ecuts_mc = (AliSpectraAODEventCuts*) _mc->Get("Event Cuts");
  AliSpectraAODTrackCuts* tcuts_mc = (AliSpectraAODTrackCuts*) _mc->Get("Track Cuts");
  // print info about mc track and Event cuts
  cout << " -- Info about MC -- "<< endl;
  ecuts_mc->PrintCuts();
  // tcuts_mc->PrintCuts();
  // proceed likewise for data
  TFile *_data = TFile::Open(dataFile.Data());
  AliSpectraAODHistoManager* hman_data = (AliSpectraAODHistoManager*) _data->Get("SpectraHistos");
  AliSpectraAODEventCuts* ecuts_data = (AliSpectraAODEventCuts*) _data->Get("Event Cuts");
  AliSpectraAODTrackCuts* tcuts_data = (AliSpectraAODTrackCuts*) _data->Get("Track Cuts");
  // print info about track and Event cuts
  cout << " -- Info about data -- " << endl;
  ecuts_data->PrintCuts();
  tcuts_data->PrintCuts();
  
  //dedx in data and MC
  TCanvas *cPIDSig=new TCanvas("cPIDSig","cPIDSig",700,500);
  cPIDSig->Divide(2,2);
  cPIDSig->cd(1);
  TH2F *PIDSig_data = (TH2F*)((TH2F*)hman_data->GetPIDHistogram("histPIDTPC"))->Clone();
  gPad->SetLogz();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_data->DrawClone("colz");
  cPIDSig->cd(2);
  TH2F *PIDSig_mc = (TH2F*)((TH2F*)hman_mc->GetPIDHistogram("histPIDTPC"))->Clone();
  gPad->SetLogz();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_mc->DrawClone("colz");
  cPIDSig->cd(3);
  TH2F *PIDSig_data = (TH2F*)((TH2F*)hman_data->GetPIDHistogram("histPIDTOF"))->Clone();
  gPad->SetLogz();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_data->DrawClone("colz");
  cPIDSig->cd(4);
  TH2F *PIDSig_mc = (TH2F*)((TH2F*)hman_mc->GetPIDHistogram("histPIDTOF"))->Clone();
  gPad->SetLogz();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_mc->DrawClone("colz");

  //nsig in data and MC
  for(Int_t ipart=0;ipart<3;ipart++){
    TCanvas *cnsig=new TCanvas(Form("cnsig%s",Particle[ipart].Data()),Form("cnsig%s",Particle[ipart].Data()),700,500);
    cnsig->Divide(2,3);
    cnsig->cd(1);
    TH2F *nsig_data = (TH2F*)((TH2F*)hman_data->GetNSigHistogram(Form("histNSig%sPtTPC",Particle[ipart].Data())))->Clone();
    nsig_data->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_data->DrawClone("colz");
    cnsig->cd(2);
    TH2F *nsig_mc = (TH2F*)((TH2F*)hman_mc->GetNSigHistogram(Form("histNSig%sPtTPC",Particle[ipart].Data())))->Clone();
    nsig_mc->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_mc->DrawClone("colz");
    cnsig->cd(3);
    TH2F *nsig_data = (TH2F*)((TH2F*)hman_data->GetNSigHistogram(Form("histNSig%sPtTOF",Particle[ipart].Data())))->Clone();
    nsig_data->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_data->DrawClone("colz");
    cnsig->cd(4);
    TH2F *nsig_mc = (TH2F*)((TH2F*)hman_mc->GetNSigHistogram(Form("histNSig%sPtTOF",Particle[ipart].Data())))->Clone();
    nsig_mc->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_mc->DrawClone("colz");
    cnsig->cd(5);
    TH2F *nsig_data = (TH2F*)((TH2F*)hman_data->GetNSigHistogram(Form("histNSig%sPtTPCTOF",Particle[ipart].Data())))->Clone();
    nsig_data->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_data->DrawClone("colz");
    cnsig->cd(6);
    TH2F *nsig_mc = (TH2F*)((TH2F*)hman_mc->GetNSigHistogram(Form("histNSig%sPtTPCTOF",Particle[ipart].Data())))->Clone();
    nsig_mc->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_mc->DrawClone("colz");
  }
  //qvec in data and MC
  for (Int_t icharge=0;icharge<2;icharge++){
    TCanvas *cqvec=new TCanvas(Form("cqvec%s",Charge[icharge].Data()),Form("cqvec%s",Charge[icharge].Data()),700,500);
    cqvec->Divide(2,1);
    cqvec->cd(1);
    TString hname=Form("histq%s",Charge[icharge].Data());
    TH2F *qvec_data = (TH2F*)((TH2F*)hman_data->GetqVecHistogram(hname.Data()))->Clone();
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    qvec_data->DrawClone("colz");
    cqvec->cd(2);
    TH2F *qvec_mc = (TH2F*)((TH2F*)hman_mc->GetqVecHistogram(hname.Data()))->Clone();
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    qvec_mc->DrawClone("colz");
    
    TCanvas *cProjqvec=new TCanvas(Form("cProjqvec%s",Charge[icharge].Data()),Form("cProjqvec%s",Charge[icharge].Data()),700,500);
    cProjqvec->Divide(3,2);
    for(Int_t iproj=0;iproj<3;iproj++){
      TH1F *proj_data=(TH1F*)qvec_data->ProjectionX("data",qvec_data->GetYaxis()->FindBin(projl[iproj]),qvec_data->GetYaxis()->FindBin(proju[iproj]));
      if(proj_data->GetEntries()==0)continue;
      proj_data->Scale(1/proj_data->GetEntries());
      proj_data->SetTitle(Form("data q%s [%.0f,%.0f]",Charge[icharge].Data(),projl[iproj],proju[iproj]));
      TH1F *proj_mc=(TH1F*)qvec_mc->ProjectionX("mc",qvec_mc->GetYaxis()->FindBin(projl[iproj]),qvec_mc->GetYaxis()->FindBin(proju[iproj]));
      proj_mc->Scale(1/proj_mc->GetEntries());
      proj_mc->SetTitle(Form("mc q%s [%.0f,%.0f]",Charge[icharge].Data(),projl[iproj],proju[iproj]));
      proj_mc->SetLineColor(2);
      cProjqvec->cd(iproj+1);
      gPad->SetGridy();
      gPad->SetGridx();
      proj_data->DrawClone();
      proj_mc->DrawClone("same");
      gPad->BuildLegend();
      proj_data->Divide(proj_mc);
      cProjqvec->cd(iproj+4);
      proj_data->DrawClone();
    }
  }
  
  //efficiencies
  Printf("\n\n-> Calculating MC Correction Factors");
  TH1F *CorrFact[6];
  TCanvas *ceff=new TCanvas("ceff","ceff",700,500);
  Bool_t UseMCDCACorrection=kTRUE;
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      //BE CAREFUL! depending on the efficiency you choose, you must change the DCA correction (data or data/mc)
      
      TString hname=Form("histPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()); //MC correction for Prim+Cont+Eff
      //TString hname=Form("histPtRecTrue%s%s",Particle[ipart].Data(),Sign[icharge].Data());//MC correction for Prim+Eff
      //TString hname=Form("histPtRecSigmaPrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());UseMCDCACorrection=kFALSE; //MC correction for Cont+Eff
      //TString hname=Form("histPtRecTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());UseMCDCACorrection=kFALSE; // Pure MC efficiency for Prim. BE CAREFUL WITH MUONS!!!
      Printf("Getting %s",hname.Data());
      CorrFact[index]=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
      CorrFact[index]->SetName(Form("CorrFact_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      CorrFact[index]->SetTitle(Form("CorrFact %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      CorrFact[index]->SetMarkerStyle(Marker[index]);
      CorrFact[index]->SetMarkerColor(Color[ipart]);
      CorrFact[index]->SetLineColor(Color[ipart]);
      hname=Form("histPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
      Printf("... and divide it by %s",hname.Data());
      CorrFact[index]->Divide(CorrFact[index],(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),1,1))->Clone(),1,1,"B");//binomial error
      if(index==0)CorrFact[index]->DrawClone();
      else CorrFact[index]->DrawClone("same");
    }
  } 
  TFile *fESD=new TFile("EffAlex/pionEffPbPb.root");
  TH1F *hEffESD=(TH1F*)fESD->Get("effMapPionTpcOnlyNeg0");
  hEffESD->DrawClone("same");
  gPad->BuildLegend();
  
    //Normalization
  printf("\n\n-> Spectra Normalization");
  AliSpectraAODEventCuts* ecuts_data = (AliSpectraAODEventCuts*) _data->Get("Event Cuts");
  Double_t events_data =  ecuts_data->NumberOfEvents();
  Printf(": accepted events: %.0f",events_data);
  
  //divide RAW for Correction Factor
  Printf("\n\n-> Using MC correction factor to correct RAW spectra");
  TH1F *Spectra[6];
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      TString hname=Form("histPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
      printf("Getting %s",hname.Data());
      Spectra[index] =(TH1F*)((TH1F*) hman_data->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
      Spectra[index]->SetName(Form("Spectra_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      Spectra[index]->SetTitle(Form("Spectra %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      Spectra[index]->SetMarkerStyle(Marker[index]);
      Spectra[index]->SetMarkerColor(Color[ipart]);
      Spectra[index]->SetLineColor(Color[ipart]);
      Printf("... and divide it by %s",hname.Data());
      Spectra[index]->Divide(CorrFact[index]);//////////////////////////////////////////////////////////////////////////////////////////FIXME
      // if(index!=3)Spectra[index]->Divide(CorrFact[index]);
      // else{
      // 	Spectra[index]=AliPWGHistoTools::MyDivideHistosDifferentBins(Spectra[index],hEffESD);
      // }
      Spectra[index]->Scale(1./events_data,"width");//NORMALIZATION
    }
  } 
  
  
  //Geant/Fluka Correction
  Printf("\nGF correction for Kaons");
  TString fnameGeanFlukaK="GFCorrection/correctionForCrossSection.321.root";
  TFile *fGeanFlukaK= new TFile(fnameGeanFlukaK.Data());
  TH1F *hGeantFlukaKPos=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionParticles");
  TH1F *hGeantFlukaKNeg=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionAntiParticles");
  for(Int_t binK=0;binK<=Spectra[1]->GetNbinsX();binK++){
    Float_t FlukaCorrKPos=hGeantFlukaKPos->GetBinContent(hGeantFlukaKPos->FindBin(Spectra[1]->GetBinCenter(binK)));
    Float_t FlukaCorrKNeg=hGeantFlukaKNeg->GetBinContent(hGeantFlukaKNeg->FindBin(Spectra[4]->GetBinCenter(binK)));
    Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPos);
    Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNeg);
    Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPos);
    Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNeg);
  }
  
  //Geant Fluka for P
  Printf("\nGF correction for Protons");
  const Int_t kNCharge=2;
  Int_t kPos=0;
  Int_t kNeg=1;
  TFile* fGFProtons = new TFile ("GFCorrection/correctionForCrossSection.root");
  TH2D * hCorrFluka[kNCharge];
  TH2D * hCorrFluka[2];
  hCorrFluka[kPos] = (TH2D*)fGFProtons->Get("gHistCorrectionForCrossSectionProtons");
  hCorrFluka[kNeg] = (TH2D*)fGFProtons->Get("gHistCorrectionForCrossSectionAntiProtons");
  for(Int_t icharge = 0; icharge < kNCharge; icharge++){
    Int_t nbins = Spectra[2]->GetNbinsX();
    Int_t nbinsy=hCorrFluka[icharge]->GetNbinsY();
    for(Int_t ibin = 0; ibin < nbins; ibin++){
      Float_t pt = Spectra[2]->GetBinCenter(ibin);
      Float_t minPtCorrection = hCorrFluka[icharge]->GetYaxis()->GetBinLowEdge(1);
      Float_t maxPtCorrection = hCorrFluka[icharge]->GetYaxis()->GetBinLowEdge(nbinsy+1);
      if (pt < minPtCorrection) pt = minPtCorrection+0.0001;
      if (pt > maxPtCorrection) pt = maxPtCorrection;
      Float_t correction = hCorrFluka[icharge]->GetBinContent(1,hCorrFluka[icharge]->GetYaxis()->FindBin(pt));
      //cout<<correction<<"     charge "<<icharge<<endl;
      if(icharge==0){
  	if (correction != 0) {// If the bin is empty this is a  0
  	  Spectra[2]->SetBinContent(ibin,Spectra[2]->GetBinContent(ibin)*correction);
  	  Spectra[2]->SetBinError(ibin,Spectra[2]->GetBinError  (ibin)*correction);
  	}else if (Spectra[2]->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
  	  cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for protons secondaries," << endl;
  	  cout << " Bin content: " << Spectra[2]->GetBinContent(ibin)  << endl;
  	}
      }
      if(icharge==1){
  	if (correction != 0) {// If the bin is empty this is a  0
  	  Spectra[5]->SetBinContent(ibin,Spectra[5]->GetBinContent(ibin)*correction);
  	  Spectra[5]->SetBinError(ibin,Spectra[5]->GetBinError  (ibin)*correction);
  	}else if (Spectra[5]->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
  	  cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for Antiprotons secondaries," << endl;
  	  cout << " Bin content: " << Spectra[5]->GetBinContent(ibin)  << endl;
  	}
      }
    }	
  }
  
  //DCA Correction
  printf("\n\n-> DCA Correction");
  DCACorrection(Spectra,hman_data,hman_mc,UseMCDCACorrection);
  
  
  TH1F *MCTruth[6];
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      hname=Form("histPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
      MCTruth[index]=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),0.9,1.1))->Clone();
      MCTruth[index]->Scale(1./events_data,"width");//NORMALIZATION
    }
  }
  
  //Contamination
  Printf("\n\n-> Contamination from MC");
  TH1F *Cont[6];
  TCanvas *ccont=new TCanvas("ccont","ccont",700,500);
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      TString hname=Form("histPtRecTrue%s%s",Particle[ipart].Data(),Sign[icharge].Data());
      Printf("Getting %s",hname.Data());
      Cont[index] =(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
      Cont[index]->SetName(Form("Cont_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      Cont[index]->SetTitle(Form("RecTrue/RecSigma %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      Cont[index]->SetMarkerStyle(Marker[index]);
      Cont[index]->SetMarkerColor(Color[ipart]);
      Cont[index]->SetLineColor(Color[ipart]);
      hname=Form("histPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
      Printf("... and divide it by %s",hname.Data());
      Cont[index]->Divide((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1));
      if(index==0)Cont[index]->DrawClone();
      else Cont[index]->DrawClone("same");
      //Spectra[index]->Multiply(Cont[index]);
    }
  } 
  gPad->BuildLegend();
  
  
  //Drawing Final Spectra
  TCanvas *cspectra=new TCanvas("cspectra","cspectra",700,500);
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      if(index==0)Spectra[index]->DrawClone();
      else Spectra[index]->DrawClone("same");
      //MCTruth[index]->DrawClone("same");
    }
  } 
  gPad->BuildLegend();
  
  //saving spectra
  Printf("\n\n-> Saving spectra in Out%s_%s.root",sname.Data(),fold.Data());
  TFile * fout=new TFile(Form("results/Out_%s_%s.root",sname.Data(),fold.Data()),"RECREATE");
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      Spectra[index]->Write();
      CorrFact[index]->Write();
      MCTruth[index]->Write();
    }
  }
  
  
  //if Bin 0-5% with no cut ratio with combined analysis
  if(ibinToCompare!=-1){
    TCanvas *CratioComb=new TCanvas("CratioComb","CratioComb",700,500);
    CratioComb->Divide(3,2);
    TString nameComb[6]={Form("cent%d_pion_plus",ibinToCompare),Form("cent%d_kaon_plus",ibinToCompare),Form("cent%d_proton_plus",ibinToCompare),
  			 Form("cent%d_pion_minus",ibinToCompare),Form("cent%d_kaon_minus",ibinToCompare),Form("cent%d_proton_minus",ibinToCompare)};
    TFile *fComb=new TFile("Combined05/SPECTRA_COMB_20120323.root");
    for(Int_t icharge=0;icharge<2;icharge++){
      for(Int_t ipart=0;ipart<3;ipart++){
  	Int_t index=ipart+3*icharge;
  	Spectra[index]->GetXaxis()->SetRangeUser(0,2.5);
  	TH1F *hcomb=fComb->Get(nameComb[index].Data())->Clone();
  	CratioComb->cd(ipart+1);
  	gPad->SetGridy();
  	gPad->SetGridx();
  	if(icharge==0)Spectra[index]->DrawClone();
  	else Spectra[index]->DrawClone("same");
  	hcomb->DrawClone("same");
  	Spectra[index]->Divide(hcomb);
  	Spectra[index]->SetMaximum(1.3);
  	Spectra[index]->SetMinimum(0.7);
  	CratioComb->cd(ipart+4);
  	gPad->SetGridy();
  	gPad->SetGridx();
  	if(icharge==0)Spectra[index]->DrawClone();
  	else Spectra[index]->DrawClone("same");
      }
    }
  }	
  

  
  //comparison with charged hadron
  Printf("\n\n-> ChargedHadron comparison");
  TCanvas *cAllCh=new TCanvas("cAllCh","cAllCh",700,500);
  cAllCh->Divide(1,4);
  TH1F *hChHad_data=(TH1F*)((TH1F*)hman_data->GetPtHistogram1D("histPtRec",-1,-1))->Clone();
  //fraction of sec in MC
  TH1F *hPrimRec_mc=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("histPtRecPrimary",-1,-1))->Clone();
  TH1F *hAllRec_mc=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("histPtRec",-1,-1))->Clone();
  for(Int_t ibin=1;ibin<=hChHad_data->GetNbinsX();ibin++){
    Double_t en_data=hChHad_data->GetBinContent(ibin);
    Double_t en_mc=hAllRec_mc->GetBinContent(ibin);
    Double_t prim_mc=hPrimRec_mc->GetBinContent(ibin);
    if(en_mc!=0)hChHad_data->SetBinContent(ibin,en_data-(en_data*(en_mc-prim_mc)*1.2/en_mc));
    Printf("Before: %.0f After: %.0f  fraction: %.1f",en_data,hChHad_data->GetBinContent(ibin),hChHad_data->GetBinContent(ibin)/en_data);
  }
  cAllCh->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hPrimRec_mc->Divide(hAllRec_mc);
  hPrimRec_mc->DrawClone();
  //efficiency for primaries
  TH1F *hEff_mc=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("histPtRecPrimary",-1,-1))->Clone();
  hEff_mc->Divide((TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("histPtGen",1,1))->Clone());
  Printf("+-------------------------------------%f ",((TH1F*)hman_mc->GetPtHistogram1D("histPtGen",1,1))->GetEntries()/1.6/ecuts_mc->NumberOfEvents());
  hChHad_data->Scale(1./events_data,"width");//NORMALIZATION
  hChHad_data->Divide(hEff_mc);//Efficiency
  hChHad_data->Scale(1./(2*tcuts_data->GetEta()));
  cAllCh->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hEff_mc->DrawClone();
  cAllCh->cd(3);
  gPad->SetGridy();
  gPad->SetGridx();
  hChHad_data->DrawClone();
  TFile *fCh=new TFile("ChargedHadron/SPECTRA_UNID_110906.root");
  TH1F *hCh=fCh->Get("hpt_c0_5");
  //invariant yield
  for(Int_t ibin=0;ibin<hCh->GetNbinsX();ibin++){
    hCh->SetBinContent(ibin,hCh->GetBinContent(ibin)*(2*TMath::Pi()*hCh->GetBinCenter(ibin)));
    hCh->SetBinError(ibin,hCh->GetBinError(ibin)*(2*TMath::Pi()*hCh->GetBinCenter(ibin)));
  }
  hCh->DrawClone("same");
  TH1F *gRatio=AliPWGHistoTools::MyDivideHistosDifferentBins(hChHad_data,hCh);
  gRatio->SetMaximum(1.5);
  gRatio->SetMinimum(.5);
  cAllCh->cd(4);
  gPad->SetGridy();
  gPad->SetGridx();
  gRatio->DrawClone("");
  
  // //Comparison of efficiency with TPCTOF ESD analysis
  Printf("\n\n-> Calculating Efficiency to be compared with ESD analysis");
  TH1F *EffTRUEPions;
  TH1F *EffSIGMAPions;
  TCanvas *cEffESD=new TCanvas("cEffESD","cEffESD",700,500);
  cEffESD->Divide(1,2);
  Int_t icharge=1;
  Int_t ipart=0;
  //using MC truth
  TString hname=Form("histPtRecTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
  Printf("Getting %s",hname.Data());
  EffTRUEPions=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
  EffTRUEPions->SetName(Form("Eff TRUE_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
  EffTRUEPions->SetTitle(Form("Eff TRUE %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
  EffTRUEPions->SetMarkerStyle(Marker[icharge]);
  EffTRUEPions->SetMarkerColor(Color[ipart]);
  EffTRUEPions->SetLineColor(Color[ipart]);
  hname=Form("histPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
  Printf("... and divide it by %s",hname.Data());
  EffTRUEPions->Divide(EffTRUEPions,(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),1,1))->Clone(),1,1,"B");//binomial error
  //using NSigma
  hname=Form("histPtRecSigmaPrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
  Printf("Getting %s",hname.Data());
  EffSIGMAPions=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
  EffSIGMAPions->SetName(Form("Eff SIGMA_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
  EffSIGMAPions->SetTitle(Form("Eff SIGMA %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
  EffSIGMAPions->SetMarkerStyle(Marker[icharge]);
  EffSIGMAPions->SetMarkerColor(Color[ipart+1]);
  EffSIGMAPions->SetLineColor(Color[ipart+1]);
  hname=Form("histPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
  Printf("... and divide it by %s",hname.Data());
  EffSIGMAPions->Divide(EffSIGMAPions,(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),1,1))->Clone(),1,1,"B");//binomial error
  cEffESD->cd(1);
  EffTRUEPions->DrawClone("lhist");
  EffSIGMAPions->DrawClone("lhistsame");
  hEffESD->DrawClone("lhistsame");
  gPad->BuildLegend();
  cEffESD->cd(2);
  TH1F *hRatioTRUE=AliPWGHistoTools::MyDivideHistosDifferentBins(EffTRUEPions,hEffESD);
  hRatioTRUE->DrawClone("lhist");
  TH1F *hRatioSIGMA=AliPWGHistoTools::MyDivideHistosDifferentBins(EffSIGMAPions,hEffESD);
  hRatioSIGMA->DrawClone("lhistsame");
  
  
  // //Muon over Pion Ratio
  // TCanvas *cMu=new TCanvas("cMu","cMu");
  // TH1F *hMuOverPi[2];
  // for(Int_t icharge=0;icharge<2;icharge++){
  //   TString hname=Form("histPtRecTruePrimaryMuon%s",Sign[icharge].Data());
  //   hMuOverPi[icharge]=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
  //   hname=Form("histPtRecTruePrimaryPion%s",Sign[icharge].Data());
  //   hMuOverPi[icharge]->Divide((TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone());
  //   if(icharge==0)hMuOverPi[icharge]->DrawClone();
  //   else hMuOverPi[icharge]->DrawClone("same");
  // }
}



void DCACorrection(TH1F **Spectra, AliSpectraAODHistoManager* hman_data, AliSpectraAODHistoManager* hman_mc,Bool_t UseMCDCACorrection){
  
  Double_t FitRange[2]={-3,3};
  Int_t nrebin=20;
  Printf("\DCACorr");
  TH1F *hcorrection[2];
  TCanvas *ccorrection=new TCanvas("DCAcorrection","DCAcorrection",700,500);
  ccorrection->Divide(2,1);
  TString sample[2]={"data","mc"};
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      for(Int_t isample=0;isample<2;isample++){
	TCanvas *cDCA=new TCanvas(Form("cDCA%d%s",index,sample[isample].Data()),Form("cDCA%d%s",index,sample[isample].Data()),700,500);
	hcorrection[isample]=(TH1F*)Spectra[index]->Clone();
	hcorrection[isample]->Reset("all");
	cDCA->Divide(8,4);
	for(Int_t ibin_data=6;ibin_data<38;ibin_data++){
	  Double_t lowedge=Spectra[index]->GetBinLowEdge(ibin_data);
	  Double_t binwidth=Spectra[index]->GetBinWidth(ibin_data);
	  if(isample==0)TH1F *hToFit =(TH1F*) ((TH1F*)hman_data->GetDCAHistogram1D(Form("histPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge+binwidth))->Clone();
	  if(isample==1)TH1F *hToFit =(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("histPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge+binwidth))->Clone();
	  TH1F *hmc1=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("histPtRecSigmaPrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge+binwidth))->Clone();
	  TH1F *hmc2=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("histPtRecSigmaSecondaryWeakDecay%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge+binwidth))->Clone();
	  TH1F *hmc3=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("histPtRecSigmaSecondaryMaterial%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge+binwidth))->Clone();
	  Double_t minentries=5;
	  if(hToFit->GetEntries()<=minentries || hmc1->GetEntries()<=minentries || hmc2->GetEntries()<=minentries || hmc3->GetEntries()<=minentries)continue;
	  hmc1->Rebin(nrebin);
	  hmc2->Rebin(nrebin);
	  hmc3->Rebin(nrebin);
	  hToFit->Rebin(nrebin);
	  //Data and MC can have different stat
	  hToFit->Sumw2();
	  hmc1->Sumw2();
	  hmc2->Sumw2();
	  hmc3->Sumw2();
	  hToFit->Scale(1./hToFit->GetEntries());
	  Double_t normMC=hmc1->GetEntries()+hmc2->GetEntries()+hmc3->GetEntries();
	  hmc1->Scale(1./normMC);
	  hmc2->Scale(1./normMC);
	  hmc3->Scale(1./normMC);
	  cDCA->cd(ibin_data-5);
	  gPad->SetGridy();
	  gPad->SetGridx();
	  gPad->SetLogy();
	  hToFit->DrawClone("lhist");
	  hmc3->DrawClone("lhist");
	  TObjArray *mc = new TObjArray(3);        // MC histograms are put in this array
	  mc->Add(hmc1);
	  mc->Add(hmc2);
	  mc->Add(hmc3);
	  TFractionFitter* fit = new TFractionFitter(hToFit,mc); // initialise
	  fit->Constrain(0,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	  fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	  fit->Constrain(2,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	  fit->SetRangeX(hToFit->GetXaxis()->FindBin(FitRange[0]),hToFit->GetXaxis()->FindBin(FitRange[1]));
	  hToFit->GetXaxis()->SetRange(hToFit->GetXaxis()->FindBin(FitRange[0]),hToFit->GetXaxis()->FindBin(FitRange[1]));
	  hToFit->SetTitle(Form("DCA distr - %s",sample[isample].Data()));
	  Int_t status = fit->Fit();               // perform the fit
	  cout << "fit status: " << status << endl;
	  if (status == 0) {                       // check on fit status
	    TH1F* result = (TH1F*) fit->GetPlot();
	    TH1F* PrimMCPred=(TH1F*)fit->GetMCPrediction(0);
	    TH1F* secStMCPred=(TH1F*)fit->GetMCPrediction(1);
	    TH1F* secMCPred=(TH1F*)fit->GetMCPrediction(2);
	    //Drawing section
	    PrimMCPred->SetLineColor(2);
	    secStMCPred->SetLineColor(6);
	    secMCPred->SetLineColor(4);
	    hToFit->DrawClone("lhist");
	    result->SetTitle("Fit result");
	    result->DrawClone("lhistsame");
	    PrimMCPred->DrawClone("lhistsame");
	    secStMCPred->DrawClone("lhistsame");
	    secMCPred->DrawClone("lhistsame");
	    Double_t v1=0,v2=0,v3=0;
	    Double_t ev1=0,ev2=0,ev3=0;
	    fit->GetResult(0,v1,ev1);
	    fit->GetResult(1,v2,ev2);
	    fit->GetResult(2,v3,ev3);
	    Printf("\n\n\n\n\nv1:%f  v2:%f  v3:%f   ev1:%f  ev2:%f  ev3:%f   ",v1,v2,v3,ev1,ev2,ev3);
	    hcorrection[isample]->SetBinContent(ibin_data,v1/(v1+v2+v3));
	    //hcorrection[isample]->SetBinError(ibin_data,TMath::Sqrt(ev1*ev1+ev2*ev2+ev3*ev3));
	    hcorrection[isample]->SetBinError(ibin_data,0);
	  }
	  else{
	    hcorrection[isample]->SetBinContent(ibin_data,1);
	    hcorrection[isample]->SetBinError(ibin_data,0);
	  }
	  Printf("deleting");
	  delete hToFit;
	}
	
	ccorrection->cd(icharge+1);
	gPad->SetGridy();
	gPad->SetGridx();
	hcorrection[isample]->SetTitle(Form("DCA corr %s %s %s",Particle[ipart].Data(),Charge[icharge].Data(),sample[isample].Data()));
	hcorrection[isample]->SetLineWidth(2);
	hcorrection[isample]->SetLineColor(Color[ipart]);
	hcorrection[isample]->SetLineStyle(isample+1);
	hcorrection[isample]->SetMarkerColor(Color[ipart]);
	hcorrection[isample]->GetXaxis()->SetRangeUser(0.2,2.5);
	if(ipart==0 && isample==0)hcorrection[isample]->DrawClone("lhist");
	else hcorrection[isample]->DrawClone("lhistsame");
	//else hcorrection[isample]->DrawClone("psame");
	// smooth the DCA correction
	// TF1 *fitfun = new TF1("fitfun","[0]+[1]/x^[2]",0.2,0.8);
	// hcorrection[isample]->Fit(fitfun,"WRVMN","",0.35,0.8);
	// fitfun->SetLineWidth(1.5);
	// fitfun->SetLineColor(ipart+1);
	// fitfun->SetLineStyle(ipart+1);
	// fitfun->DrawClone("same");
	// for(Int_t ibin=1;ibin<30;ibin++){
	// 	hcorrection[isample]->SetBinContent(ibin,fitfun->Eval(hcorrection[isample]->GetBinCenter(ibin)));
	// }
      }
      Spectra[index]->Multiply(hcorrection[0]);//multiply for data
      Printf("DCACorrection for DATA used: Spectra[index]->Multiply(hcorrection[0])")
	if(UseMCDCACorrection){
	  Spectra[index]->Divide(hcorrection[1]); //divide by Monte Carlo
	  Printf("DCACorrection for MC used: Spectra[index]->Divide(hcorrection[1]")
	    }
    }
  }
  ccorrection->cd(1);
  gPad->BuildLegend();
  ccorrection->cd(2);
  gPad->BuildLegend();
}


