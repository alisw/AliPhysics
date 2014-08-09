#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFractionFitter.h>
#endif

enum species {kPion, kKaon, kProton};
enum chargesign {kPositive,kNegative};

Float_t PerformFit(TH1F* fHistDCA, TObjArray *mc,TString partname,Float_t xyMax);

void FitDCADistributions(TString period="LHC10d",
			 TString MCname="LHC10f6a",  
			 Int_t iSpecies=2, 
			 Int_t cSign=1
			 )
{
  gROOT->SetStyle("Plain");
  gStyle->SetFillColor(0);
  gStyle->SetOptStat(0000);
  //binning
  const Int_t nptbins = 22;
  Double_t ptbins[nptbins+1]={0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0};
  //Final correction
  TH1F *hPrimAllDATA=new TH1F("hPrimAllDATA","hPrimAllDATA",nptbins,ptbins);
  TH1F *hPrimAllMC=new TH1F("hPrimAllMC","hPrimAllMC",nptbins,ptbins);
  TH1F *hPrimAllDATAMC=new TH1F("hPrimAllDATAMC","hPrimAllDATAMC",nptbins,ptbins);
  
  
  TString fnameMC=Form("%s/AnalysisResults.root",MCname.Data());
  TString lnameMC="clistITSsaMult-1to-1";
  TString lnameCutMC="DCACutMult-1to-1";
  TString fnameDATA=Form("%s/AnalysisResults.root",period.Data());
  TString lnameDATA="clistITSsaMult-1to-1";
  TString lnameCutDATA="DCACutMult-1to-1";

  const Int_t firstbin=1;
  Int_t rebin=10;
  Bool_t optSmooth=kTRUE;


  TString pCharge="Pos";
  TString pcharge="pos";
  if(cSign==kNegative){
    pCharge="Neg";
    pcharge="neg";
  }
  TString species[3]={"Pi","K","P"};
  TString speciesRoot[6]={"#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#bar{p}"};
  Int_t idPart=iSpecies*2+cSign;
  TString partname=Form("%s%s",pCharge.Data(),species[iSpecies].Data());
  printf("%s\n",partname.Data());

  // MC histograms 
  TFile *fMC=new TFile(fnameMC,"READ");  
  TH1F *fHistDCAMC[nptbins]; 
  TH1F *fHistDCAPrim[nptbins]; 
  TH1F *fHistDCAsec[nptbins]; 
  TH1F *fHistDCAsecSt[nptbins]; 
  TDirectory *dirFileMC=(TDirectory*)fMC->Get("PWG2SpectraITSsa");
  TList *liMC = (TList*)dirFileMC->Get(lnameMC.Data());
  for(Int_t m=0;m<nptbins;m++){
    fHistDCAMC[m] = (TH1F*)liMC->FindObject(Form("fHistDCA%s%i",partname.Data(),m));
    fHistDCAPrim[m] = (TH1F*)liMC->FindObject(Form("fHistMCPrimDCA%s%i",partname.Data(),m));
    fHistDCAsec[m] = (TH1F*)liMC->FindObject(Form("fHistMCSecMatDCA%s%i",partname.Data(),m));
    fHistDCAsecSt[m] = (TH1F*)liMC->FindObject(Form("fHistMCSecStDCA%s%i",partname.Data(),m));
    fHistDCAMC[m]->Rebin(rebin);
    fHistDCAPrim[m]->Rebin(rebin);
    fHistDCAsec[m]->Rebin(rebin);
    fHistDCAsecSt[m]->Rebin(rebin);
  }
  TList *liCutMC = (TList*)dirFileMC->Get(lnameCutMC.Data());
  Bool_t setDCA=kFALSE;

  Float_t DCAcut=7.;
  Double_t xyPMC[3];   //parameters of DCA cut
  xyPMC[0]=36.; //MC LHC10d1
  xyPMC[1]=43.9;
  xyPMC[2]=1.3;
  TF1* fDCAxyCutMC=0x0;
  TF1* fDCAzCutMC=0x0;
  Float_t nSigmaDCAxyMC=0.;
  Float_t nSigmaDCAzMC=0.;
  if(liCutMC){
    fDCAxyCutMC=(TF1*)liCutMC->FindObject("fDCAxyCutFunc");
    fDCAzCutMC=(TF1*)liCutMC->FindObject("fDCAzCutFunc");
    nSigmaDCAxyMC=fDCAxyCutMC->GetParameter(3);
    nSigmaDCAzMC=fDCAzCutMC->GetParameter(3);
    printf("DCA cut parameters from TF1: %f %f %f --- DCA cut = %f\n",fDCAxyCutMC->GetParameter(0),fDCAxyCutMC->GetParameter(1),fDCAxyCutMC->GetParameter(2),nSigmaDCAxyMC);
    setDCA=kTRUE;
  }
  if(!setDCA){
    printf("DCA cut values for MC not found in the output file, use dafault values. DCAcut=%f with default resolution\n",DCAcut);
    fDCAxyCutMC = new TF1("fDCAxyCutFunc","[3]*([0]+[1]/TMath::Power(TMath::Abs(x),[2]))",0.05,10.);
    for(Int_t ipar=0; ipar<3; ipar++) fDCAxyCutMC->SetParameter(ipar,xyPMC[ipar]);
    fDCAxyCutMC->SetParameter(3,DCAcut);
    nSigmaDCAxyMC=DCAcut;
  }

  //DATA histograms
  TFile *fDATA=new TFile(fnameDATA,"READ");
  TDirectory *dirFileDATA=(TDirectory*)fDATA->Get("PWG2SpectraITSsa");
  TList *liDATA = (TList*)dirFileDATA->Get(lnameDATA.Data());
  TH1F *fHistDCADATA[nptbins]; 
  for(Int_t m=0;m<nptbins;m++){
    fHistDCADATA[m] = (TH1F*)liDATA->FindObject(Form("fHistDCA%s%i",partname.Data(),m));
    fHistDCADATA[m]->Rebin(rebin);
  }
  TList *liCutDATA = (TList*)dirFileDATA->Get(lnameCutDATA.Data());
  setDCA=kFALSE;
  Float_t nSigmaDCAxyDATA=0.;
  Float_t nSigmaDCAzDATA=0.;
  Double_t xyP[3];   // parameters of DCA cut
  xyP[0]=32.7;  //DATA 7 TeV pass2
  xyP[1]=44.8;
  xyP[2]=1.3;
  TF1* fDCAxyCutDATA=0x0;
  TF1* fDCAzCutDATA=0x0;
  if(liCutDATA){
    fDCAxyCutDATA=(TF1*)liCutDATA->FindObject("fDCAxyCutFunc");
    fDCAzCutDATA=(TF1*)liCutDATA->FindObject("fDCAzCutFunc");
    nSigmaDCAxyDATA=fDCAxyCutDATA->GetParameter(3);
    nSigmaDCAzDATA=fDCAzCutDATA->GetParameter(3);
    printf("DCA cut parameters from TF1: %f %f %f --- DCA cut = %f\n",fDCAxyCutDATA->GetParameter(0),fDCAxyCutDATA->GetParameter(1),fDCAxyCutDATA->GetParameter(2),nSigmaDCAxyDATA);
    setDCA=kTRUE;
  }
  if(!setDCA){
    printf("DCA cut values for DATA not found in the output file, use dafault values. DCAcut=%f with default resolution\n",DCAcut);
    fDCAxyCutDATA = new TF1("fDCAxyCutFunc","[3]*([0]+[1]/TMath::Power(TMath::Abs(x),[2]))",0.05,10.);
    for(Int_t ipar=0; ipar<3; ipar++) fDCAxyCutDATA->SetParameter(ipar,xyP[ipar]);
    fDCAxyCutDATA->SetParameter(3,DCAcut);
    nSigmaDCAxyDATA=DCAcut;
  }
  
  if(TMath::Abs(nSigmaDCAxyDATA-nSigmaDCAxyMC)<0.01) DCAcut=nSigmaDCAxyDATA;
  else{
    printf("ERROR: DCAxy cuts do not match between data and MC\n");
    return;
  }
  if(TMath::Abs(nSigmaDCAzDATA-nSigmaDCAzMC)>0.01){
    printf("ERROR: DCAz cuts do not match between data (%f) and MC (%f) \n",nSigmaDCAzDATA,nSigmaDCAzMC);
    return;
  }

  TCanvas *cfitDATA=new TCanvas("cfitDATA","cfitDATA");
  cfitDATA->Divide(6,4);
  TCanvas *cfitMC=new TCanvas("cfitMC","cfitMC");
  cfitMC->Divide(6,4);

  for(Int_t ibin=firstbin;ibin<nptbins;ibin++){    
    
    Double_t xyMax =fDCAxyCutDATA->Eval(TMath::Abs((ptbins[ibin+1]+ptbins[ibin])/2))/10000.;
    Double_t xyMaxMC =fDCAxyCutMC->Eval(TMath::Abs((ptbins[ibin+1]+ptbins[ibin])/2))/10000.;

    TObjArray *mcTemplates = 0x0;
    if(partname=="PosP"){
      mcTemplates = new TObjArray(3);        // MC histograms are put in this array
      mcTemplates->Add(fHistDCAPrim[ibin]);
      mcTemplates->Add(fHistDCAsecSt[ibin]);
      mcTemplates->Add(fHistDCAsec[ibin]);
    }else{
      mcTemplates = new TObjArray(2);        // MC histograms are put in this array
      mcTemplates->Add(fHistDCAPrim[ibin]);
      mcTemplates->Add(fHistDCAsecSt[ibin]);
    }
    //    if(fHistDCADATA[ibin]->Integral()==0) continue;

    ////////////////// Fit on DATA
    cfitDATA->cd(ibin);
    gPad->SetLogy();
    Double_t primAllDATA=PerformFit(fHistDCADATA[ibin],mcTemplates,partname,xyMax);
    
    ////////////////// Fit on MC
    cfitMC->cd(ibin);
    gPad->SetLogy();
    Double_t primAllMC=PerformFit(fHistDCAMC[ibin],mcTemplates,partname,xyMaxMC);
    
    if(hPrimAllMC!=0 && hPrimAllDATA!=0){
      hPrimAllDATA->Fill((ptbins[ibin]+ptbins[ibin+1])/2,primAllDATA);
      hPrimAllMC->Fill((ptbins[ibin]+ptbins[ibin+1])/2,primAllMC);
      hPrimAllDATAMC->Fill((ptbins[ibin]+ptbins[ibin+1])/2,primAllDATA/primAllMC);
    }
  }
  
  //Adding MC truth
  TH1F* hAll=(TH1F*)liMC->FindObject(Form("hHist%sNSigmaMean%d",pCharge.Data(),iSpecies));
  TH1F* hPrim=(TH1F*)liMC->FindObject(Form("hHist%sNSigmaPrimMean%d",pCharge.Data(),iSpecies));
  TH1F* hPrimRecMC=(TH1F*)liMC->FindObject(Form("fHistPrimMC%sReco%d",pcharge.Data(),iSpecies));
  TH1F* hSecMatRecMC=(TH1F*)liMC->FindObject(Form("fHistSecMatMC%sReco%d",pcharge.Data(),iSpecies));
  TH1F* hSecStrRecMC=(TH1F*)liMC->FindObject(Form("fHistSecStrMC%sReco%d",pcharge.Data(),iSpecies));
  hSecStrRecMC->Add(hSecMatRecMC);
  hSecStrRecMC->Add(hPrimRecMC);
  hPrim->SetLineColor(8);
  hPrim->SetLineWidth(2);
  hPrim->SetTitle("Prim/All from MC Truth");
  hPrim->Divide(hAll);
  hPrimRecMC->SetLineColor(11);
  hPrimRecMC->SetLineWidth(2);
  hPrimRecMC->SetTitle("Prim/All from MC Truth Reco");
  hPrimRecMC->Divide(hSecStrRecMC);

  hPrimAllDATAMC->GetXaxis()->SetRangeUser(0.09,0.79);
  hPrimAllMC->GetXaxis()->SetRangeUser(0.09,0.79);
  hPrimAllDATA->GetXaxis()->SetRangeUser(0.09,0.79);
  if(optSmooth) hPrimAllDATAMC->Smooth(1,"R");
  
  hPrimAllDATA->SetTitle("Prim/all data");
  hPrimAllDATA->SetLineColor(2);
  hPrimAllDATA->SetLineWidth(2);
  hPrimAllMC->SetTitle("Prim/all mc");
  hPrimAllMC->SetLineColor(4);
  hPrimAllMC->SetLineWidth(2);
  hPrimAllDATAMC->SetTitle("Prim/all data/mc");
  hPrimAllDATAMC->SetLineColor(1);
  hPrimAllDATAMC->SetLineWidth(2);

  TCanvas *cPrimAll=new TCanvas("cPrimAll","cPrimAll");
  hPrimAllDATA->SetMinimum(0.7);
  hPrimAllDATA->SetMaximum(1.1);
  hPrimAllDATA->Draw("l");
  hPrimAllMC->Draw("lsame");
  hPrimAllDATAMC->Draw("lsame");
  hPrim->Draw("lsame");
  hPrimRecMC->Draw("lsame");
  gPad->BuildLegend();
  TLatex* tsp=new TLatex(0.4,0.75,speciesRoot[idPart].Data());
  tsp->SetTextFont(43);
  tsp->SetTextSize(28);
  tsp->SetNDC();
  tsp->Draw();
  cPrimAll->Update();

  TString fout;
  fout=Form("DCACorr%s_%s_%.0fDCA_%s_TFraction.root",period.Data(),MCname.Data(),DCAcut,partname.Data());
  TFile *out=new TFile(fout.Data(),"recreate");
  hPrimAllDATA->Write();
  hPrimAllMC->Write();
  hPrimAllDATAMC->Write();
  out->Close();
  delete out;
  
}

Float_t PerformFit(TH1F* fHistDCA, TObjArray *mc,TString partname,Float_t xyMax){
  
  Double_t prim=0,secSt=0,sec=0;
  TFractionFitter* fit = new TFractionFitter(fHistDCA, mc); // initialise
  fit->Constrain(0,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
  fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
  if(partname=="PosP")fit->Constrain(2,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
  
  //fit->SetRangeX(200,1000);                    // use only the first 15 bins in the fit
  Int_t status = fit->Fit();               // perform the fit
  cout << "fit status: " << status << endl;
  
  if (status == 0) {                       // check on fit status
    TH1F* result = (TH1F*) fit->GetPlot();
    TH1F* PrimMCPred=(TH1F*)fit->GetMCPrediction(0);
    TH1F* secStMCPred=(TH1F*)fit->GetMCPrediction(1);
    TH1F* secMCPred=0x0;
    if(partname=="PosP"){
      secMCPred=(TH1F*)fit->GetMCPrediction(2);
      secMCPred->SetLineColor(4);
    }
    //Drawing section
    PrimMCPred->SetLineColor(2);
    secStMCPred->SetLineColor(6);
    fHistDCA->SetTitle("DCA distr - data");
    fHistDCA->DrawCopy("Ep");
    result->SetTitle("Fit result");
    result->DrawCopy("same");
    Double_t value,error;
    
    fit->GetResult(0,value,error);
    PrimMCPred->Scale(fHistDCA->GetSumOfWeights()*value/PrimMCPred->GetSumOfWeights());
    PrimMCPred->SetTitle("Primaries");
    PrimMCPred->DrawCopy("same");
    fit->GetResult(1,value,error);
    secStMCPred->Scale(fHistDCA->GetSumOfWeights()*value/secStMCPred->GetSumOfWeights());
    secStMCPred->SetTitle("Sec from strangeness");
    secStMCPred->DrawCopy("same");
    if(partname=="PosP" && secMCPred){
      fit->GetResult(2,value,error);
      secMCPred->Scale(fHistDCA->GetSumOfWeights()*value/secMCPred->GetSumOfWeights());
      secMCPred->SetTitle("Sec from material");
      secMCPred->DrawCopy("same");
    }
    prim=PrimMCPred->Integral(PrimMCPred->FindBin(-xyMax),PrimMCPred->FindBin(xyMax));
    secSt=secStMCPred->Integral(secStMCPred->FindBin(-xyMax),secStMCPred->FindBin(xyMax));
    if(partname=="PosP")sec=secMCPred->Integral(secMCPred->FindBin(-xyMax),secMCPred->FindBin(xyMax));
  }
  else{
    prim=1;
    secSt=1;
    if(partname=="PosP")sec=1;
  }
  Printf("Yields:  primary=%f material=%f strange=%f\n",prim,sec,secSt);
  return prim/(prim+secSt+sec);
}
