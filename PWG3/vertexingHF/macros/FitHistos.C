void FitHistos(){
  Bool_t useParFiles=kFALSE;
  Int_t load=gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/macros/LoadLibraries.C");
  LoadLibraries(useParFiles);
  Int_t rebin=4;//4;
  
    
  Int_t typeb=1;
  Int_t types=0;
  Int_t factor4refl=0;
  TFile *f=new TFile("AnalysisResults.root");
  TDirectory *dir = (TDirectory*)f->Get("PWG3_D2H_InvMassDplus");
  TList *list = (TList*)dir->Get("coutputDplus");
  TList *listc = (TList*)dir->Get("coutputDplusCuts");
  AliRDHFCutsDplustoKpipi* cutsDp=(AliRDHFCutsDplustoKpipi*)listc->At(0);
  Int_t nPtBins=cutsDp->GetNPtBins();
  printf("Number of pt bins = %d\n",nPtBins);
  Float_t *ptlims=cutsDp->GetPtBinLimits();
  ptlims[nPtBins]=ptlims[nPtBins-1]+4.;
  TH1F* hSignal=new TH1F("hSignal","hSignal",nPtBins,ptlims);
  TH1F* hBackground=new TH1F("hBackground","hBackground",nPtBins,ptlims);
  TH1F* hSignificance=new TH1F("hSignificance","hSignificance",nPtBins,ptlims);
  TH1F** hmass=new TH1F*[nPtBins];
  for(Int_t i=0;i<nPtBins;i++){
    hmass[i]=(TH1F*)list->FindObject(Form("hMassPt%dTC",i));
  }
  Int_t nMassBins=hmass[1]->GetNbinsX();
  Double_t hmin=hmass[1]->GetBinLowEdge(3);
  Double_t hmax=hmass[1]->GetBinLowEdge(nMassBins-2)+hmass[1]->GetBinWidth(nMassBins-2);
  TCanvas** c1= new TCanvas*[nPtBins];
  Int_t iPad=1;
  //  Int_t iBin=2;

  AliHFMassFitter** fitter=new AliHFMassFitter*[nPtBins];
  Double_t sig,errsig,s,errs,b,errb;
  for(Int_t iBin=0; iBin<nPtBins; iBin++){
    c1[iBin]=new TCanvas(Form("cbin%d",iBin),Form("cbin%d",iBin));    
    if(hmass[iBin]->GetEntries()>1000.){
      fitter[iBin]=new AliHFMassFitter(hmass[iBin],hmin, hmax,rebin,typeb,types);
      Bool_t out=fitter[iBin]->MassFitter(0);
      hmass[iBin]->Rebin(rebin);
      hmass[iBin]->Draw();
      TF1* fB1=fitter[iBin]->GetBackgroundFullRangeFunc();
      TF1* fB2=fitter[iBin]->GetBackgroundRecalcFunc();
      TF1* fM=fitter[iBin]->GetMassFunc();
      if(fB1) fB1->DrawClone("same");
      if(fB2) fB2->DrawClone("same");
      if(fM) fM->DrawClone("same");
      //      fitter[iBin]->DrawHere(c1[iBin]);    
      fitter[iBin]->Signal(3,s,errs);
      fitter[iBin]->Background(3,b,errb);
      fitter[iBin]->Significance(3,sig,errsig);
      hSignal->SetBinContent(iBin+1,s);
      hSignal->SetBinError(iBin+1,errs);
      hBackground->SetBinContent(iBin+1,b);
      hBackground->SetBinError(iBin+1,errb);
      hSignificance->SetBinContent(iBin+1,sig);
      hSignificance->SetBinError(iBin+1,errsig);
    }
  }
  for(Int_t iBin=0; iBin<nPtBins; iBin++){
    printf("Bin %d  Signal %.2f+-%.2f\n",iBin,sig[iBin],errsig[iBin]);
  }

  TCanvas* csig=new TCanvas("csig","Results",1200,600);
  csig->Divide(3,1);
  csig->cd(1);
  hSignal->SetMarkerStyle(20);
  hSignal->Draw("P");
  csig->cd(2);
  hBackground->SetMarkerStyle(20);
  hBackground->Draw("P");
  csig->cd(3);
  hSignificance->SetMarkerStyle(20);
  hSignificance->Draw("P");

  TFile* outf=new TFile("RawYield.root","recreate");
  outf->cd();
  hSignal->Write();
  hBackground->Write();
  hSignificance->Write();
}
