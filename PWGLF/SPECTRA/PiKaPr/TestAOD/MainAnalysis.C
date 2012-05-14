class AliPID;
#include "QAPlots.C"

TString Charge[]={"Pos","Neg"};
TString Sign[]={"Plus","Minus"};
TString Particle[]={"Pion","Kaon","Proton"};
Int_t Color[3]={1,2,4};
Int_t Marker[6]={20,21,22,24,25,26};
Double_t projl[3]={0,20,80};
Double_t proju[3]={10,40,90};
Double_t Range[3]={0.3,0.3,0.5}; // LowPt range for pi k p
enum ECharge_t {
  kPositive,
  kNegative,
  kNCharges
};


void MainAnalysis()  {
  
  gSystem->Load("libCore.so");  
  //gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  //gSystem->Load("libVMC");
  gSystem->Load("libTree");
  //gSystem->Load("libProof");
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
  //gSystem->Load("libMinuit");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGLFSPECTRA");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
 
  // Set Masses
  Double_t mass[3];
  mass[0]   = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  mass[1]   = TDatabasePDG::Instance()->GetParticle("K+")->Mass();
  mass[2] = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
  
  //TString fold="3SigmaPID_AOD086-090_FilterBit10";
  TString fold="5SigmaPID_AOD046_FilterBit6";
  Int_t ibinToCompare=-1;
  
  TString sname="Cent0to5_QVec0.0to100.0";ibinToCompare=0;
  
  TString dataFile = Form("output/%s/Pt.AOD.1._data_ptcut_%s.root",fold.Data(),sname.Data());
  TString mcFile =Form("output/%s/Pt.AOD.1._MC_%s.root",fold.Data(),sname.Data());
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
  
  //QAPlots(hman_data,hman_mc);

  //efficiencies
  Printf("\n\n-> Calculating MC Correction Factors");
  TH1F *CorrFact[6];
  TCanvas *ceff=new TCanvas("ceff","ceff",700,500);
  Bool_t UseMCDCACorrection=kTRUE;
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      //BE CAREFUL! depending on the efficiency you choose, you must change the DCA correction (data or data/mc)
      
      TString hname=Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()); //MC correction for Prim+Cont+Eff
      //TString hname=Form("hHistPtRecTrue%s%s",Particle[ipart].Data(),Sign[icharge].Data());//MC correction for Prim+Eff
      //TString hname=Form("hHistPtRecSigmaPrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());UseMCDCACorrection=kFALSE; //MC correction for Cont+Eff
      //TString hname=Form("hHistPtRecTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());UseMCDCACorrection=kFALSE; // Pure MC efficiency for Prim. BE CAREFUL WITH MUONS!!!
      Printf("Getting %s",hname.Data());
      CorrFact[index]=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
      CorrFact[index]->SetName(Form("CorrFact_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      CorrFact[index]->SetTitle(Form("CorrFact %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      CorrFact[index]->SetMarkerStyle(Marker[index]);
      CorrFact[index]->SetMarkerColor(Color[ipart]);
      CorrFact[index]->SetLineColor(Color[ipart]);
      hname=Form("hHistPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
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
  Double_t events_mc =  ecuts_mc->NumberOfEvents();
  Printf(": accepted events (MC): %.0f",events_mc);
  
  //divide RAW for Correction Factor
  Printf("\n\n-> Using MC correction factor to correct RAW spectra");
  TH1F *Spectra[6];
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      TString hname=Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
      printf("Getting %s",hname.Data());
      Spectra[index] =(TH1F*)((TH1F*) hman_data->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
      Spectra[index]->SetName(Form("Spectra_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      Spectra[index]->SetTitle(Form("Spectra %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      Spectra[index]->SetMarkerStyle(Marker[index]);
      Spectra[index]->SetMarkerColor(Color[ipart]);
      Spectra[index]->SetLineColor(Color[ipart]);
      Printf("... and divide it by %s",hname.Data());
      Spectra[index]->Divide(CorrFact[index]);
      Spectra[index]->Scale(1./events_data,"width");//NORMALIZATION
    }
  } 
  
  //Put Bin Content = 0
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      for(Int_t ibin=0;ibin<Spectra[index]->GetNbinsX();ibin++){
	if(Spectra[index]->GetBinCenter(ibin)<Range[ipart]){
	  Spectra[index]->SetBinContent(ibin,0);
	  Spectra[index]->SetBinError(ibin,0);
	}
      }
    }
  }
  //DCA Correction with the "right" DCA sample
  DCACorrection(Spectra,hman_data,hman_mc,UseMCDCACorrection);
  
  //DCA Correction forcing loose DCA
  // TString fold_LooseDCA="5SigmaPID_AOD046_FilterBit5";
  // TString dataFile_LooseDCA = Form("output/%s/Pt.AOD.1._data_ptcut_%s.root",fold_LooseDCA.Data(),sname.Data());
  // TString mcFile_LooseDCA =Form("output/%s/Pt.AOD.1._MC_%s.root",fold_LooseDCA.Data(),sname.Data());
  // TFile *_mc_LooseDCA = TFile::Open(mcFile_LooseDCA.Data());
  // AliSpectraAODHistoManager* hman_mc_LooseDCA = (AliSpectraAODHistoManager*) _mc_LooseDCA->Get("SpectraHistos");
  // TFile *_data_LooseDCA = TFile::Open(dataFile_LooseDCA.Data());
  // AliSpectraAODHistoManager* hman_data_LooseDCA = (AliSpectraAODHistoManager*) _data_LooseDCA->Get("SpectraHistos");
  // DCACorrection(Spectra,hman_data_LooseDCA,hman_mc_LooseDCA,UseMCDCACorrection);
  
  
  
  
  
  //GFCorrection
  GFCorrection(Spectra,tcuts_data);
  
  TH1F *MCTruth[6];
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      hname=Form("hHistPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
      MCTruth[index]=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),1,1))->Clone();
      MCTruth[index]->SetName(Form("MCTruth_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      MCTruth[index]->SetTitle(Form("MCTruth_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      MCTruth[index]->Scale(1./events_mc,"width");//NORMALIZATION
    }
  }
  
  //Drawing Final Spectra
  TCanvas *cspectra=new TCanvas("cspectra","cspectra",700,500);
  gPad->SetGridy();
  gPad->SetGridy();
  gPad->SetLogy();
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      if(index==0)Spectra[index]->DrawClone();
      else Spectra[index]->DrawClone("same");
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
    TFile *fComb=new TFile("Combined05/SPECTRA_COMB_20120412.root");
    TH1F *Spectra_copy[6]=0x0;
    for(Int_t icharge=0;icharge<2;icharge++){
      for(Int_t ipart=0;ipart<3;ipart++){
	Int_t index=ipart+3*icharge;
	TH1F *htmp=(TH1F*)((TH1F*)Spectra[index])->Clone("");
	htmp->GetXaxis()->SetRangeUser(0,2.5);
  	TH1F *hcomb=fComb->Get(nameComb[index].Data())->Clone();
  	CratioComb->cd(ipart+1);
  	gPad->SetGridy();
  	gPad->SetGridx();
  	if(icharge==0)htmp->DrawClone();
  	else htmp->DrawClone("same");
  	hcomb->DrawClone("same");
  	htmp->Divide(hcomb);
  	htmp->SetMaximum(1.3);
  	htmp->SetMinimum(0.7);
  	CratioComb->cd(ipart+4);
  	gPad->SetGridy();
  	gPad->SetGridx();
  	if(icharge==0)htmp->DrawClone();
  	else htmp->DrawClone("same");
      }
    }
  }	
  
  //comparison with charged hadron
  Printf("\n\n-> ChargedHadron comparison");
  TH1F *hChHad_data=(TH1F*)((TH1F*)hman_data->GetPtHistogram1D("hHistPtRec",-1,-1))->Clone();
  //fraction of sec in MC
  TH1F *hPrimRec_mc=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("hHistPtRecPrimary",-1,-1))->Clone();
  TH1F *hAllRec_mc=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("hHistPtRec",-1,-1))->Clone();
  for(Int_t ibin=1;ibin<=hChHad_data->GetNbinsX();ibin++){
    Double_t en_data=hChHad_data->GetBinContent(ibin);
    Double_t en_mc=hAllRec_mc->GetBinContent(ibin);
    Double_t prim_mc=hPrimRec_mc->GetBinContent(ibin);
    if(en_mc!=0)hChHad_data->SetBinContent(ibin,en_data-(en_data*(en_mc-prim_mc)*1.1/en_mc));
    //Printf("Before: %.0f After: %.0f  fraction: %.1f",en_data,hChHad_data->GetBinContent(ibin),hChHad_data->GetBinContent(ibin)/en_data);
  }
  hPrimRec_mc->Divide(hAllRec_mc);
  //efficiency for primaries
  TH1F *hEff_mc=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("hHistPtRecPrimary",-1,-1))->Clone();
  hEff_mc->Divide((TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("hHistPtGen",1,1))->Clone());
  TCanvas *cAllChFactors=new TCanvas("cAllChFactors","cAllChFactors",700,500);
  cAllChFactors->Divide(1,2);
  cAllChFactors->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hPrimRec_mc->SetTitle("Prim/All, charged hadron pure MC");
  hPrimRec_mc->DrawClone("lhist");
  gPad->BuildLegend();
  cAllChFactors->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hEff_mc->SetTitle("Efficiency for Primaries, charged hadron pure MC");
  hEff_mc->DrawClone("lhist");
  gPad->BuildLegend();
  //Printf("--------%f ",((TH1F*)hman_mc->GetPtHistogram1D("hHistPtGen",1,1))->GetEntries()/1.6/ecuts_mc->NumberOfEvents());
  hChHad_data->Scale(1./events_data,"width");//NORMALIZATION
  hChHad_data->Divide(hEff_mc);//Efficiency
  hChHad_data->Scale(1./(2*tcuts_data->GetEta()));
  hChHad_data->SetTitle("All Ch from AOD");
  TCanvas *cAllCh=new TCanvas("cAllCh","cAllCh",700,500);
  cAllCh->Divide(1,2);
  cAllCh->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hChHad_data->DrawClone();
  TFile *fCh=new TFile("ChargedHadron/SPECTRA_UNID_110906.root");
  TH1F *hCh=fCh->Get("hpt_c0_5");
  hCh->SetTitle("All Ch from Jacek");
  hCh->SetMarkerColor(2);
  hCh->SetLineColor(2);
  //invariant yield
  for(Int_t ibin=0;ibin<hCh->GetNbinsX();ibin++){
    hCh->SetBinContent(ibin,hCh->GetBinContent(ibin)*(2*TMath::Pi()*hCh->GetBinCenter(ibin)));
    hCh->SetBinError(ibin,hCh->GetBinError(ibin)*(2*TMath::Pi()*hCh->GetBinCenter(ibin)));
  }
  hCh->DrawClone("same");
  gPad->BuildLegend();
  TH1F *gRatio=AliPWGHistoTools::MyDivideHistosDifferentBins(hChHad_data,hCh);
  gRatio->SetMaximum(1.3);
  gRatio->SetMinimum(.7);
  cAllCh->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  gRatio->DrawClone("");
  //fitting
  TCanvas *cFitChargHad=new TCanvas("cFitChargHad","cFitChargHad",700,500);
  gPad->SetGridy();
  gPad->SetGridx();
  hChHad_data->DrawClone();
  //Fitting the sum of all particles
  AliPWGFunc * fm = new AliPWGFunc;
  fm->SetVarType(AliPWGFunc::kdNdpt);
  Float_t fitmin = 0.3;
  Float_t fitmax = 3;
  TF1 * func = 0;
  Int_t normPar = 0;
  func = fm->GetBGBW(0.13,0.6,0.3, 1, 1e5);
  func->SetParLimits(1, 0.1, 0.99);
  func->SetParLimits(2, 0.01, 1);
  func->SetParLimits(3, 0.01, 2);
  TH1F * hToFit = hChHad_data;
  hToFit->Fit(func,"N","VMRN",fitmin,fitmax);
  // // if(!AliPWGHistoTools::Fit(hToFit,func,fitmin,fitmax)) {
  // //   cout << " FIT ERROR " << endl;
  // //   return;      
  // // }
  Double_t yieldTools, yieldETools;
  Double_t partialYields[3],partialYieldsErrors[3]; 
  AliPWGHistoTools::GetYield(hToFit, func, yieldTools, yieldETools,0, 100, partialYields,partialYieldsErrors);
  func->DrawClone("same");   
  Printf("TOTAL YIELD (AOD Charged Hadron) : %f +- %f",yieldTools,yieldETools);
  //Fit All Charged
  hToFit = hCh;
  hToFit->Fit(func,"N","VMRN",fitmin,fitmax);
  AliPWGHistoTools::GetYield(hToFit, func, yieldTools, yieldETools,0, 100, partialYields,partialYieldsErrors);
  func->SetLineColor(2);
  hCh->DrawClone("same");
  func->DrawClone("same");   
  gPad->BuildLegend();
  Printf("TOTAL YIELD (JACEK): %f +- %f",yieldTools,yieldETools);
  //sumID vs AllCh
  //Convert spectra to dNdeta and sum
  TH1F * hsum = (TH1F*) Spectra[0]->Clone();
  hsum->Reset("all");
  Double_t epsilon = 0.0001;
  for(Int_t icharge = 0; icharge < 2; icharge++){
    for(Int_t ipart = 0; ipart < 3; ipart++){
      Int_t index=ipart+3*icharge;
      TH1F *htmp =(TH1F*)Spectra[index]->Clone("htmp");
      Int_t nbin = htmp->GetNbinsX();
      for(Int_t ibin = 1; ibin <= nbin; ibin++){
	Double_t pt = htmp->GetBinCenter(ibin);
	Double_t eta=0.8;//////////////////eta range///////////////////////////////////////
	Double_t jacobian =eta2y(pt,mass[ipart],eta)/eta;
	//Printf("jacobian: %f pt:%f   BinContent:%f  mass:%f",jacobian,pt,htmp->GetBinContent(ibin),mass[ipart]);
	htmp->SetBinContent(ibin,htmp->GetBinContent(ibin)*jacobian);
	htmp->SetBinError(ibin,htmp->GetBinError(ibin)*jacobian);
	Int_t ibinSum = hsum->FindBin(pt);
	if ( htmp->GetBinContent(ibin) > 0 && 
	     (TMath::Abs(htmp->GetBinLowEdge(ibin)   - hsum->GetBinLowEdge(ibinSum)) > epsilon || 
	      TMath::Abs(htmp->GetBinLowEdge(ibin+1) - hsum->GetBinLowEdge(ibinSum+1)) )
	     ) {
	  cout << "DISCREPANCY IN BIN RANGES" << endl;
	  cout << pt << " " << ibinSum << " " << ibin  << "; " << h->GetBinContent(ibin) << endl
	       << h->GetBinLowEdge(ibin) << "-"  << h->GetBinLowEdge(ibin+1) << endl
	       << hsum->GetBinLowEdge(ibinSum) << "-"  << hsum->GetBinLowEdge(ibinSum+1) << endl;
	  cout << "" << endl;	    
	}
      }
      hsum->Add(htmp);
    }
  }
  hsum->SetTitle("Sum ID from AOD");
  TCanvas *cChargHadComp=new TCanvas("cChargHadComp","cChargHadComp",700,500);
  cChargHadComp->Divide(1,2);
  cChargHadComp->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hsum->DrawClone();
  hToFit = hsum;
  hToFit->Fit(func,"N","VMRN",fitmin,fitmax);
  AliPWGHistoTools::GetYield(hToFit, func, yieldTools, yieldETools,0, 100, partialYields,partialYieldsErrors);
  func->SetLineColor(2);
  Printf("TOTAL YIELD (Pi+K+p): %f +- %f",yieldTools,yieldETools);
  hChHad_data->SetMarkerColor(2);
  hChHad_data->SetLineColor(2);
  hChHad_data->DrawClone("same");
  gPad->BuildLegend();
  cChargHadComp->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hsum->Divide(hChHad_data);
  hsum->SetMaximum(1.2);
  hsum->SetMinimum(.8);
  hsum->DrawClone("");
  
  return;
  // //Comparison of efficiency with TPCTOF ESD analysis
  // Printf("\n\n-> Calculating Efficiency to be compared with ESD analysis");
  // TH1F *EffTRUEPions;
  // TH1F *EffSIGMAPions;
  // TCanvas *cEffESD=new TCanvas("cEffESD","cEffESD",700,500);
  // cEffESD->Divide(1,2);
  // Int_t icharge=1;
  // Int_t ipart=0;
  // //using MC truth
  // TString hname=Form("hHistPtRecTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
  // Printf("Getting %s",hname.Data());
  // EffTRUEPions=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
  // EffTRUEPions->SetName(Form("Eff TRUE_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
  // EffTRUEPions->SetTitle(Form("Eff TRUE %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
  // EffTRUEPions->SetMarkerStyle(Marker[icharge]);
  // EffTRUEPions->SetMarkerColor(Color[ipart]);
  // EffTRUEPions->SetLineColor(Color[ipart]);
  // hname=Form("hHistPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
  // Printf("... and divide it by %s",hname.Data());
  // EffTRUEPions->Divide(EffTRUEPions,(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),1,1))->Clone(),1,1,"B");//binomial error
  // //using NSigma
  // hname=Form("hHistPtRecSigmaPrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
  // Printf("Getting %s",hname.Data());
  // EffSIGMAPions=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
  // EffSIGMAPions->SetName(Form("Eff SIGMA_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
  // EffSIGMAPions->SetTitle(Form("Eff SIGMA %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
  // EffSIGMAPions->SetMarkerStyle(Marker[icharge]);
  // EffSIGMAPions->SetMarkerColor(Color[ipart+1]);
  // EffSIGMAPions->SetLineColor(Color[ipart+1]);
  // hname=Form("hHistPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
  // Printf("... and divide it by %s",hname.Data());
  // EffSIGMAPions->Divide(EffSIGMAPions,(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),1,1))->Clone(),1,1,"B");//binomial error
  // cEffESD->cd(1);
  // EffTRUEPions->DrawClone("lhist");
  // EffSIGMAPions->DrawClone("lhistsame");
  // hEffESD->DrawClone("lhistsame");
  // gPad->BuildLegend();
  // cEffESD->cd(2);
  // TH1F *hRatioTRUE=AliPWGHistoTools::MyDivideHistosDifferentBins(EffTRUEPions,hEffESD);
  // hRatioTRUE->DrawClone("lhist");
  // TH1F *hRatioSIGMA=AliPWGHistoTools::MyDivideHistosDifferentBins(EffSIGMAPions,hEffESD);
  // hRatioSIGMA->DrawClone("lhistsame");
  // return;
  
}




void DCACorrection(TH1F **Spectra, AliSpectraAODHistoManager* hman_data, AliSpectraAODHistoManager* hman_mc,Bool_t UseMCDCACorrection){
  printf("\n\n-> DCA Correction");
  
  Double_t FitRange[2]={-1.,1.};
  Int_t nrebin=2;
  Printf("\DCACorr");
  TH1F *hcorrection[2];
  TCanvas *ccorrection=new TCanvas("DCAcorrection","DCAcorrection",700,500);
  TCanvas *cRatiocorrection=new TCanvas("DCARatiocorrection","DCARatiocorrection",700,500);
  cRatiocorrection->Divide(2,1);
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
	  if(isample==0)TH1F *hToFit =(TH1F*) ((TH1F*)hman_data->GetDCAHistogram1D(Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge+binwidth))->Clone();
	  if(isample==1)TH1F *hToFit =(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge+binwidth))->Clone();
	  TH1F *hmc1=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigmaPrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge+binwidth))->Clone();
	  TH1F *hmc2=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigmaSecondaryWeakDecay%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge+binwidth))->Clone();
	  TH1F *hmc3=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigmaSecondaryMaterial%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge+binwidth))->Clone();
	  //Double_t minentries=200;
	  //if(hToFit->GetEntries()<=minentries || hmc1->GetEntries()<=minentries || hmc2->GetEntries()<=minentries || hmc3->GetEntries()<=minentries)continue;
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
	    
	    Double_t v1=0,v2=0,v3=0;
	    Double_t ev1=0,ev2=0,ev3=0;
	    //first method, use directly the fit result
	    fit->GetResult(0,v1,ev1);
	    fit->GetResult(1,v2,ev2);
	    fit->GetResult(2,v3,ev3);
	    result->Scale(hToFit->GetSumOfWeights()/result->GetSumOfWeights());
	    PrimMCPred->Scale(hToFit->GetSumOfWeights()*v1/PrimMCPred->GetSumOfWeights());
	    secStMCPred->Scale(hToFit->GetSumOfWeights()*v2/secStMCPred->GetSumOfWeights());
	    secMCPred->Scale(hToFit->GetSumOfWeights()*v3/secMCPred->GetSumOfWeights());
	    //second method, integrated the MC predisction, it should give the same as the first method
	    // v1=PrimMCPred->Integral();
	    // v2=secStMCPred->Integral();
	    // v3=secMCPred->Integral();
	    //Printf("\n\n\n\n\nv1:%f  v2:%f  v3:%f   ev1:%f  ev2:%f  ev3:%f   ",v1,v2,v3,ev1,ev2,ev3);
	    hcorrection[isample]->SetBinContent(ibin_data,v1/(v1+v2+v3));
	    hcorrection[isample]->SetBinError(ibin_data,0);
	    //Drawing section
	    PrimMCPred->SetLineColor(2);
	    secStMCPred->SetLineColor(6);
	    secMCPred->SetLineColor(4);
	    hToFit->SetMinimum(0.0001);
	    hToFit->DrawClone("lhist");
	    result->SetTitle("Fit result");
	    result->DrawClone("lhistsame");
	    PrimMCPred->DrawClone("lhistsame");
	    secStMCPred->DrawClone("lhistsame");
	    secMCPred->DrawClone("lhistsame");
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
	// smooth the DCA correction
	// TF1 *fitfun = new TF1("fitfun","[0]+[1]*x^[2]+[3]*x^[4]",0.2,2.5);
	// hcorrection[isample]->Fit(fitfun,"WRN","N",0.35,2);
	// fitfun->SetLineWidth(1.5);
	// fitfun->SetLineColor(Color[ipart]);
	// fitfun->SetLineStyle(isample);
	// fitfun->DrawClone("same");
	// for(Int_t ibin=1;ibin<30;ibin++){
	//   hcorrection[isample]->SetBinContent(ibin,fitfun->Eval(hcorrection[isample]->GetBinCenter(ibin)));
	// }
      }
      Spectra[index]->Multiply(hcorrection[0]);//multiply for data
      Printf("DCACorrection for DATA used: Spectra[index]->Multiply(hcorrection[0])")
	if(UseMCDCACorrection){
	  Spectra[index]->Divide(hcorrection[1]); //divide by Monte Carlo
	  Printf("DCACorrection for MC used: Spectra[index]->Divide(hcorrection[1]")
	    }
      cRatiocorrection->cd(icharge+1);
      gPad->SetGridy();
      gPad->SetGridx();
      hcorrection[0]->Divide(hcorrection[1]);
      if(ipart==0)hcorrection[0]->DrawClone("lhist");
      else hcorrection[0]->DrawClone("lhistsame");
    }
  }
  ccorrection->cd(1);
  gPad->BuildLegend();
  ccorrection->cd(2);
  gPad->BuildLegend();
}


//////////////

Double_t eta2y(Double_t pt, Double_t mass, Double_t eta){
  Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
  return TMath::ASinH(pt / mt * TMath::SinH(eta));
}
///////////////////////
void GFCorrection(TH1F **Spectra,AliSpectraAODTrackCuts* tcuts_data){
  //Geant/Fluka Correction
  Printf("\nGF correction for Kaons");
  //Getting GF For Kaons in TPC
  TGraph *gGFCorrectionKaonPlus=new TGraph();
  gGFCorrectionKaonPlus->SetName("gGFCorrectionKaonPlus");
  gGFCorrectionKaonPlus->SetTitle("gGFCorrectionKaonPlus");
  TGraph *gGFCorrectionKaonMinus=new TGraph();
  gGFCorrectionKaonMinus->SetName("gGFCorrectionKaonMinus");
  gGFCorrectionKaonMinus->SetTitle("gGFCorrectionKaonMinus");
  TString fnameGeanFlukaK="GFCorrection/correctionForCrossSection.321.root";
  TFile *fGeanFlukaK= new TFile(fnameGeanFlukaK.Data());
  TH1F *hGeantFlukaKPos=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionParticles");
  TH1F *hGeantFlukaKNeg=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionAntiParticles");
  //getting GF func for Kaons with TOF
  TF1 *fGFKPosTracking;
  fGFKPosTracking = TrackingEff_geantflukaCorrection(3,kPositive);
  TF1 *fGFKNegTracking;
  fGFKNegTracking = TrackingEff_geantflukaCorrection(3,kNegative);
  TF1 *fGFKPosMatching;
  fGFKPosMatching = TOFmatchMC_geantflukaCorrection(3,kPositive);
  TF1 *fGFKNegMatching;
  fGFKNegMatching = TOFmatchMC_geantflukaCorrection(3,kNegative);
  for(Int_t binK=0;binK<=Spectra[1]->GetNbinsX();binK++){
    if(Spectra[1]->GetBinCenter(binK)<tcuts_data->GetPtTOFMatching()){//use TPC GeantFlukaCorrection
      Float_t FlukaCorrKPos=hGeantFlukaKPos->GetBinContent(hGeantFlukaKPos->FindBin(Spectra[1]->GetBinCenter(binK)));
      Float_t FlukaCorrKNeg=hGeantFlukaKNeg->GetBinContent(hGeantFlukaKNeg->FindBin(Spectra[4]->GetBinCenter(binK)));
      Printf("TPC Geant/Fluka: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPos,FlukaCorrKNeg);
      Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPos);
      Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNeg);
      Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPos);
      Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNeg);
      gGFCorrectionKaonPlus->SetPoint(binK,Spectra[1]->GetBinCenter(binK),FlukaCorrKPos);
      gGFCorrectionKaonMinus->SetPoint(binK,Spectra[4]->GetBinCenter(binK),FlukaCorrKNeg);
    }else{
      gGFCorrectionKaonPlus->SetPoint(binK,Spectra[1]->GetBinCenter(binK),0);
      gGFCorrectionKaonMinus->SetPoint(binK,Spectra[4]->GetBinCenter(binK),0);
      Float_t FlukaCorrKPosTracking=fGFKPosTracking->Eval(Spectra[1]->GetBinCenter(binK));
      Float_t FlukaCorrKNegTracking=fGFKNegTracking->Eval(Spectra[1]->GetBinCenter(binK));
      Printf("TPC/TOF Geant/Fluka Tracking: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPosTracking,FlukaCorrKNegTracking);
      Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPosTracking);
      Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNegTracking);
      Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPosTracking);
      Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNegTracking);
      Float_t FlukaCorrKPosMatching=fGFKPosMatching->Eval(Spectra[1]->GetBinCenter(binK));
      Float_t FlukaCorrKNegMatching=fGFKNegMatching->Eval(Spectra[1]->GetBinCenter(binK));
      Printf("TPC/TOF Geant/Fluka Matching: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPosMatching,FlukaCorrKNegMatching);
      Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPosMatching);
      Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNegMatching);
      Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPosMatching);
      Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNegMatching);
    }
  }
  
  //Geant Fluka for P in TPC
  Printf("\nGF correction for Protons");
  const Int_t kNCharge=2;
  Int_t kPos=0;
  Int_t kNeg=1;
  TFile* fGFProtons = new TFile ("GFCorrection/correctionForCrossSection.root");
  TH2D * hCorrFluka[kNCharge];
  TH2D * hCorrFluka[2];
  hCorrFluka[kPos] = (TH2D*)fGFProtons->Get("gHistCorrectionForCrossSectionProtons");
  hCorrFluka[kNeg] = (TH2D*)fGFProtons->Get("gHistCorrectionForCrossSectionAntiProtons");
  TGraph *gGFCorrectionProtonPlus=new TGraph();
  gGFCorrectionProtonPlus->SetName("gGFCorrectionProtonPlus");
  gGFCorrectionProtonPlus->SetTitle("gGFCorrectionProtonPlus");
  TGraph *gGFCorrectionProtonMinus=new TGraph();
  gGFCorrectionProtonMinus->SetName("gGFCorrectionProtonMinus");
  gGFCorrectionProtonMinus->SetTitle("gGFCorrectionProtonMinus");
   //getting GF func for Kaons with TPCTOF
  TF1 *fGFpPosTracking;
  fGFpPosTracking = TrackingEff_geantflukaCorrection(4,kPositive);
  TF1 *fGFpNegTracking;
  fGFpNegTracking = TrackingEff_geantflukaCorrection(4,kNegative);
  TF1 *fGFpPosMatching;
  fGFpPosMatching = TOFmatchMC_geantflukaCorrection(4,kPositive);
  TF1 *fGFpNegMatching;
  fGFpNegMatching = TOFmatchMC_geantflukaCorrection(4,kNegative);
  
  for(Int_t icharge = 0; icharge < kNCharge; icharge++){
    Int_t nbins = Spectra[2]->GetNbinsX();
    Int_t nbinsy=hCorrFluka[icharge]->GetNbinsY();
    for(Int_t ibin = 0; ibin < nbins; ibin++){
      if(Spectra[2]->GetBinCenter(ibin)<tcuts_data->GetPtTOFMatching()){//use TPC GeantFlukaCorrection
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
	  gGFCorrectionProtonPlus->SetPoint(ibin,pt,correction);
	}else if (Spectra[2]->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
  	  cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for protons secondaries, Pt:"<< pt<< endl;
  	  cout << " Bin content: " << Spectra[2]->GetBinContent(ibin)  << endl;
  	}
      }
      if(icharge==1){
  	if (correction != 0) {// If the bin is empty this is a  0
  	  Spectra[5]->SetBinContent(ibin,Spectra[5]->GetBinContent(ibin)*correction);
  	  Spectra[5]->SetBinError(ibin,Spectra[5]->GetBinError  (ibin)*correction);
	  gGFCorrectionProtonMinus->SetPoint(ibin,pt,correction);
  	}else if (Spectra[5]->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
  	  cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for Antiprotons secondaries, Pt:"<< pt<< endl;
  	  cout << " Bin content: " << Spectra[5]->GetBinContent(ibin)  << endl;
  	}
      }
      }else{
	gGFCorrectionProtonPlus->SetPoint(binK,Spectra[1]->GetBinCenter(binK),0);
	gGFCorrectionProtonMinus->SetPoint(binK,Spectra[4]->GetBinCenter(binK),0);
	Float_t FlukaCorrpPosTracking=fGFpPosTracking->Eval(Spectra[1]->GetBinCenter(binK));
	Float_t FlukaCorrpNegTracking=fGFpNegTracking->Eval(Spectra[1]->GetBinCenter(binK));
	Printf("TPC/TOF Geant/Fluka Tracking: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrpPosTracking,FlukaCorrpNegTracking);
	Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrpPosTracking);
	Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrpNegTracking);
	Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrpPosTracking);
	Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrpNegTracking);
	Float_t FlukaCorrpPosMatching=fGFpPosMatching->Eval(Spectra[1]->GetBinCenter(binK));
	Float_t FlukaCorrpNegMatching=fGFpNegMatching->Eval(Spectra[1]->GetBinCenter(binK));
	Printf("TPC/TOF Geant/Fluka Matching: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrpPosMatching,FlukaCorrpNegMatching);
	Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrpPosMatching);
	Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrpNegMatching);
	Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrpPosMatching);
	Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrpNegMatching);
      }
    }//end loop on bins	
  }
  gGFCorrectionKaonPlus->SetLineColor(kRed);
  gGFCorrectionKaonMinus->SetLineColor(kRed+2);
  gGFCorrectionProtonPlus->SetLineColor(kGreen);
  gGFCorrectionProtonMinus->SetLineColor(kGreen+2);
  fGFKPosTracking->SetLineColor(kRed);
  fGFKNegTracking->SetLineColor(kRed+2);
  fGFKPosMatching->SetLineColor(kRed);
  fGFKNegMatching->SetLineColor(kRed+2);
  fGFpPosTracking->SetLineColor(kGreen);
  fGFpNegTracking->SetLineColor(kGreen+2);
  fGFpPosMatching->SetLineColor(kGreen);
  fGFpNegMatching->SetLineColor(kGreen+2);
  fGFKPosTracking->SetLineStyle(2);
  fGFKNegTracking->SetLineStyle(2);
  fGFKPosMatching->SetLineStyle(3);
  fGFKNegMatching->SetLineStyle(3);
  fGFpPosTracking->SetLineStyle(2);
  fGFpNegTracking->SetLineStyle(2);
  fGFpPosMatching->SetLineStyle(3);
  fGFpNegMatching->SetLineStyle(3);
  fGFKPosTracking->SetRange(.6,5);
  fGFKNegTracking->SetRange(.6,5);
  fGFKPosMatching->SetRange(.6,5);
  fGFKNegMatching->SetRange(.6,5);
  fGFpPosTracking->SetRange(.6,5);
  fGFpNegTracking->SetRange(.6,5);
  fGFpPosMatching->SetRange(.6,5);
  fGFpNegMatching->SetRange(.6,5);
 
  TCanvas * GFCorrection = new TCanvas ("GFCorrection","GFCorrection",700,500);
  gPad->SetGridx();
  gPad->SetGridy();
  gGFCorrectionKaonPlus->DrawClone("al");
  gGFCorrectionKaonMinus->DrawClone("lsame");
  gGFCorrectionProtonPlus->DrawClone("lsame");
  gGFCorrectionProtonMinus->DrawClone("lsame");
  fGFKPosTracking->DrawClone("lsame");
  fGFKNegTracking->DrawClone("lsame");
  fGFKPosMatching->DrawClone("lsame");
  fGFKNegMatching->DrawClone("lsame");
  fGFpPosTracking->DrawClone("lsame");
  fGFpNegTracking->DrawClone("lsame");
  fGFpPosMatching->DrawClone("lsame");
  fGFpNegMatching->DrawClone("lsame");
  gPad->BuildLegend();
}



///////////
TF1 *
TrackingEff_geantflukaCorrection(Int_t ipart, Int_t icharge)
{

  if (ipart == 3 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "TrackingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    return f;
  }
  else if (ipart == 4 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "TrackingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
  }
  else
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "TrackingPtGeantFlukaCorrectionNull(x)", 0., 5.);

  return f;
}

Double_t
TrackingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
  return 1.;
}

Double_t
TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  return (1 - 0.129758 *TMath::Exp(-pTmc*0.679612));
}

Double_t
TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  return TMath::Min((0.972865 + 0.0117093*pTmc), 1.);
}
///////////////////////////////////////////
TF1 *
TOFmatchMC_geantflukaCorrection(Int_t ipart, Int_t icharge)
{

  if (ipart == 3 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "MatchingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    return f;
  }
  else if (ipart == 4 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "MatchingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
  }
  else
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "MatchingPtGeantFlukaCorrectionNull(x)", 0., 5.);

  return f;
}


Double_t
MatchingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
  return 1.;
}

Double_t 
MatchingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  Float_t ptTPCoutP =pTmc*(1-6.81059e-01*TMath::Exp(-pTmc*4.20094));
  return (TMath::Power(1 - 0.129758*TMath::Exp(-ptTPCoutP*0.679612),0.07162/0.03471));
}

Double_t
MatchingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  Float_t ptTPCoutK=pTmc*(1- 3.37297e-03/pTmc/pTmc - 3.26544e-03/pTmc);
  return TMath::Min((TMath::Power(0.972865 + 0.0117093*ptTPCoutK,0.07162/0.03471)), 1.);
}
