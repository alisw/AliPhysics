#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPaveStats.h>
#include <TMath.h>
#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TLegendEntry.h>
#endif

TString configFileName="configfile4lowptanalysis.txt";
TString fileNameMC="";
TString suffix="";
TString fileNameToy="";

const Int_t maxPtBins=30;
Int_t nPtBins=8;
Double_t binLims[maxPtBins+1]={0.,1.,2.,3.,4.,5.,6.,8.,12.};
Int_t ptcol[maxPtBins]={1,kRed+1,kRed,kGreen+2,kCyan,4,kOrange+2,kMagenta,kMagenta+2,kBlue+1,kGray,kGray+2,kGreen,kYellow+7};

TString ptDWeight="";
TString ptBWeight="";
Double_t maxMult=-1;//200;

TH1F* hAccToyFine=0x0;
TH1F* hAccToy=0x0;
TH2D* hPtVsYGenAccToy=0x0;
TH2D* hPtVsYGenLimAccToy=0x0;

TF1* funcPtWeight=0x0;
TF1* funcPtBWeight=0x0;

TString fileWeightName="";
TString histoWeightName="";
TH1F* hMultWeight=0x0;

Bool_t useToyMC=kFALSE;

void ComputeAndWriteEff(TList* l, TString dCase);
void ComputeAndWriteEffOld(TList* l, TString dCase, TString var3="Mult");
Bool_t ReadConfig(TString configName);

void ComputeEfficiencyFromCombinHF(TString configInput=""){

  if(configInput!="") configFileName=configInput.Data();

  if(configFileName.Length()>0){
    if(gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",configFileName.Data()))==0){
      printf("Read configuration from file %s\n",configFileName.Data());
      Bool_t readOK=ReadConfig(configFileName);
      if(!readOK){
	printf("Error in reading configuration file\n");
	return;
      }
    }
  }

  // multiplicity weights
  if(fileWeightName.Length()>0 && gSystem->AccessPathName(fileWeightName.Data())==0){
    printf("Open file with multiplicity weights\n");
    TFile* filw = new TFile(fileWeightName.Data());
    hMultWeight=(TH1F*)filw->Get(histoWeightName.Data());
    if(hMultWeight){
      hMultWeight->SetStats(0);
      hMultWeight->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
      hMultWeight->GetYaxis()->SetTitle("Weight");
      hMultWeight->GetYaxis()->SetTitleOffset(1.4);
      hMultWeight->SetMaximum(3.);
      hMultWeight->SetMinimum(0.);
    }
  }


  // pt weights  
  if(ptDWeight=="FONLL5overLHC13d3"){
    funcPtWeight=new TF1("funcPtWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,30.);
    funcPtWeight->SetParameters(2.94999e+00,3.47032e+00,2.81278e+00,2.5,1.93370e-02,3.86865e+00,-1.54113e-01,8.86944e-02,2.56267e-02);
  }else if (ptDWeight=="FONLL7overLHC10f6a"){
    funcPtWeight=new TF1("funcPtWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
    funcPtWeight->SetParameters(2.41522e+01,4.92146e+00,6.72495e+00,2.5,6.15361e-03,4.78995e+00,-4.29135e-01,3.99421e-01,-1.57220e-02);
  }else if (ptDWeight=="FONLL7overLHC10f7a"){
    funcPtWeight=new TF1("funcPtWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
    funcPtWeight->SetParameters(3.59525,2.67871,2.70402,1.72578,4.78167e-03,4.90992,-1.26424e-01,8.21269e-02,-1.26425e-01);
  }else if(ptDWeight=="FLAToverLHC10f7a"){
    funcPtWeight=new TF1("funcPtWeight","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.15,40.);
    funcPtWeight->SetParameters(1.99498e-01,-9.90532e-02,3.03645e-02,7.42483e-04);
  }else if(ptDWeight=="FONLL5overLHC20g2"){
    funcPtWeight=new TF1("funcPtWeight","[0]+[1]*x+[2]*TMath::Exp(-(x-[3])*(x-[3])/2/[4]/[4])+[5]*TMath::Exp(-(x-[6])*(x-[6])/2/[7]/[7])+[8]*TMath::Exp(-x/[9])",0.,40.);
    funcPtWeight->SetParameters(-1.604e+00,7.346e-02,5.901e-01,1.279e+01,4.785e+00,2.508e+00,2.051e+00,5.081e+00,7.675e-01,2.037e-01);
  }else if(ptDWeight=="DATA010overLHC20g2a"){
    funcPtWeight=new TF1("funcPtWeight","[0]+[1]*x+[2]*TMath::Exp(-(x-[3])*(x-[3])/2/[4]/[4])+[5]*TMath::Exp(-(x-[6])*(x-[6])/2/[7]/[7])+[8]*TMath::Exp(-x/[9])",0.,25.);
    funcPtWeight->SetParameters(5.890e-02,-2.140e-03,4.851e-01,-1.162e+00,4.878e+00,7.749e-01,2.108e+00,1.314e+00,6.717e-01,2.573e-01,2.573e-01);
  }else if(ptDWeight=="DATA3050overLHC20g2b"){
    funcPtWeight=new TF1("funcPtWeight","[0]+[1]*x+[2]*TMath::Exp(-(x-[3])*(x-[3])/2/[4]/[4])+[5]*TMath::Exp(-(x-[6])*(x-[6])/2/[7]/[7])+[8]*TMath::Exp(-x/[9])",0.,25.);
    funcPtWeight->SetParameters(9.774e-02,-3.386e-03,6.433e-01,7.425e-01,4.889e+00,4.575e-01,2.223e+00,1.254e+00,5.560e-01,2.083e-01,2.083e-01);
  }else{
    funcPtWeight=new TF1("funcPtWeight","[0]");
    funcPtWeight->SetParameter(0,1.);
  }
  if(ptBWeight=="FONLL5overLHC19c3"){
    funcPtBWeight=new TF1("ff","[0]+[1]*x+[2]*TMath::Exp(-(x-[3])*(x-[3])/2/[4]/[4])+[5]*TMath::Exp(-(x-[6])*(x-[6])/2/[7]/[7])+[8]*x*x",0.,50.);
    funcPtBWeight->SetParameters(5.359e-01,-1.921e-02,7.247e-01,6.899e+00,5.310e+00,2.273e-01,3.474e+00,1.444e+00,1.800e-04);
  }else if(ptBWeight=="FONLL5overLHC20g2"){
    funcPtBWeight=new TF1("ff","[0]+[1]*x+[2]*TMath::Exp(-(x-[3])*(x-[3])/2/[4]/[4])+[5]*TMath::Exp(-(x-[6])*(x-[6])/2/[7]/[7])+[8]*TMath::Exp(-x/[9])",0.,50.);
    funcPtBWeight->SetParameters(4.866e-01,-9.217e-03,6.428e+00,-2.662e+00,1.144e+01,1.100e+01,-2.319e+00,2.725e+00,-1.358e+01,4.622e+00);
  }else{
    funcPtBWeight=new TF1("funcPtWeight","[0]");
    funcPtBWeight->SetParameter(0,1.);
  }

  TString dirName=Form("PWG3_D2H_InvMassDzeroLowPt%s",suffix.Data());
  TString lstName=Form("coutputDzero%s",suffix.Data());
  
  if(gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",fileNameMC.Data())) !=0){
    printf("File %s with raw data results does not exist -> exiting\n",fileNameMC.Data());
    return;
  }
  TFile* fil=new TFile(fileNameMC.Data());
  TDirectoryFile* df=(TDirectoryFile*)fil->Get(dirName.Data());
  if(!df){
    printf("TDirectoryFile %s not found in TFile\n",dirName.Data());
    fil->ls();
    return;
  }
  TList* l=(TList*)df->Get(lstName.Data());


  // aceptance from toy MC
  TFile* fileAccToy=new TFile(fileNameToy.Data());
  hPtVsYGenAccToy=(TH2D*)fileAccToy->Get("hPtVsYGenAcc");
  hPtVsYGenLimAccToy=(TH2D*)fileAccToy->Get("hPtVsYGenLimAcc");

  TH1F* hPtGenAccToy=(TH1F*)fileAccToy->Get("hPtGenAcc");
  TH1F* hPtGenLimAccToy=(TH1F*)fileAccToy->Get("hPtGenLimAcc");
  hAccToyFine=(TH1F*)fileAccToy->Get("hAccVsPt");
  TH1F* hPtGenAccToyR=(TH1F*)hPtGenAccToy->Rebin(nPtBins,"hPtGenAccToyReb",binLims);
  TH1F* hPtGenLimAccToyR=(TH1F*)hPtGenLimAccToy->Rebin(nPtBins,"hPtGenLimAccToyReb",binLims);
  hAccToy=(TH1F*)hPtGenAccToyR->Clone("hAccToy");
  hPtGenAccToyR->Sumw2();
  hPtGenLimAccToyR->Sumw2();
  hAccToy->Divide(hPtGenAccToyR,hPtGenLimAccToyR,1,1,"B");
  hAccToy->SetLineColor(kGreen+2);
  hAccToy->SetLineWidth(3);
  hAccToy->SetMarkerStyle(27);
  hAccToy->SetMarkerSize(1.8);
  hAccToy->SetMarkerColor(kGreen+2);
  hAccToy->SetStats(0);
  hAccToyFine->SetLineColor(kGreen+2);

  TFile* out=new TFile(Form("outputEff%s.root",suffix.Data()),"recreate");
  hAccToy->Write();
  out->Close();



  TH2F* hEventMultZv=(TH2F*)l->FindObject("hEventMultZv");
  if(hEventMultZv){
    TH2F* hEventMultZvEvSel=(TH2F*)l->FindObject("hEventMultZvEvSel");
    hEventMultZv->GetXaxis()->SetTitle("Z_{vertex} (cm)");
    hEventMultZv->GetYaxis()->SetTitle("N_{tracklets}");
    hEventMultZvEvSel->GetXaxis()->SetTitle("Z_{vertex} (cm)");
    hEventMultZvEvSel->GetYaxis()->SetTitle("N_{tracklets}");
    Int_t binzm10=hEventMultZvEvSel->GetXaxis()->FindBin(-9.999);
    Int_t binzp10=hEventMultZvEvSel->GetXaxis()->FindBin(9.999);
    printf("%d %f    --- %d %f\n",binzm10,hEventMultZvEvSel->GetXaxis()->GetBinLowEdge(binzm10),
	   binzp10,hEventMultZvEvSel->GetXaxis()->GetBinUpEdge(binzp10));
    TH1D* hMultAllEv=hEventMultZv->ProjectionY("hMultAllEv");
    TH1D* hMultEvZ10=hEventMultZv->ProjectionY("hMultEvZ10",binzm10,binzp10);
    TH1D* hMultEvSel=hEventMultZvEvSel->ProjectionY("hMultEvSel");
    hMultAllEv->SetLineColor(1);
    hMultEvZ10->SetLineColor(kGreen+2);
    hMultEvSel->SetLineColor(6);
    hMultAllEv->GetYaxis()->SetTitle("Entries");
    
    TH1D* hZvAllEv=hEventMultZv->ProjectionX("hZvAllEv");
    TH1D* hZvEvSel=hEventMultZvEvSel->ProjectionX("hZvEvSel");
    hZvAllEv->SetLineColor(1);
    hZvEvSel->SetLineColor(6);
    hZvAllEv->GetYaxis()->SetTitle("Entries");
    
    TH1D* hRatioMultEv=(TH1D*)hMultEvSel->Clone("hRatioMultEv");
    hRatioMultEv->Divide(hMultEvSel,hMultEvZ10,1.,1.,"B");
    hRatioMultEv->SetStats(0);
    hRatioMultEv->GetYaxis()->SetTitle("Ratio");
    hRatioMultEv->SetLineColor(4);
    TH1D* hRatioMultEvAll=(TH1D*)hMultEvSel->Clone("hRatioMultEvAll");
    hRatioMultEvAll->Divide(hMultEvSel,hMultAllEv,1.,1.,"B");
    hRatioMultEvAll->SetStats(0);
    hRatioMultEvAll->GetYaxis()->SetTitle("Ratio");
    hRatioMultEvAll->SetLineColor(6);
    TH1D* hRatioMultEvZ10=(TH1D*)hMultEvZ10->Clone("hRatioMultEvZ10");
    hRatioMultEvZ10->Divide(hMultEvZ10,hMultAllEv,1.,1.,"B");
    hRatioMultEvZ10->SetStats(0);
    hRatioMultEvZ10->GetYaxis()->SetTitle("Ratio");
    TH1D* hRatioZvEv=(TH1D*)hZvEvSel->Clone("hRatioZvEv");
    hRatioZvEv->Divide(hZvEvSel,hZvAllEv,1.,1.,"B");
    hRatioZvEv->SetStats(0);
    hRatioZvEv->GetYaxis()->SetTitle("Ratio");
    
    TCanvas* cev=new TCanvas("cev","Event selection",1300,950);
    cev->Divide(2,2);
    cev->cd(1);
    gPad->SetLogy();
    gPad->SetTickx();
    gPad->SetTicky();
    if(maxMult>0) hMultAllEv->GetXaxis()->SetRangeUser(0.,maxMult);
    hMultAllEv->Draw();
    gPad->Update();
    TPaveStats *sta=(TPaveStats*)hMultAllEv->GetListOfFunctions()->FindObject("stats");
    sta->SetY1NDC(0.7);
    sta->SetY2NDC(0.89);
    gPad->Modified();
    hMultEvZ10->Draw("sames");
    gPad->Update();
    TPaveStats *stz=(TPaveStats*)hMultEvZ10->GetListOfFunctions()->FindObject("stats");
    stz->SetY1NDC(0.5);
    stz->SetY2NDC(0.69);
    stz->SetTextColor(hMultEvZ10->GetLineColor());
    gPad->Modified();
    hMultEvSel->Draw("sames");
    gPad->Update();
    TPaveStats *stb=(TPaveStats*)hMultEvSel->GetListOfFunctions()->FindObject("stats");
    stb->SetY1NDC(0.3);
    stb->SetY2NDC(0.49);
    stb->SetTextColor(hMultEvSel->GetLineColor());
    gPad->Modified();
    cev->cd(3);
    gPad->SetTickx();
    gPad->SetTicky();
    hRatioMultEvAll->SetMinimum(0.);
    if(maxMult>0) hRatioMultEvAll->GetXaxis()->SetRangeUser(0.,maxMult);
    hRatioMultEvAll->Draw();
    hRatioMultEvZ10->Draw("same");
    hRatioMultEv->Draw("same");
    TLegend* leg=new TLegend(0.45,0.3,0.89,0.7);
    leg->AddEntry(hRatioMultEvAll,"EvSel / AllEv","L")->SetTextColor(hRatioMultEvAll->GetLineColor());
    leg->AddEntry(hRatioMultEvZ10,"Ev(|z|<10) / AllEv","L")->SetTextColor(hRatioMultEvZ10->GetLineColor());
    leg->AddEntry(hRatioMultEv,"EvSel / Ev(|z|<10)","L")->SetTextColor(hRatioMultEv->GetLineColor());
    leg->Draw();
    cev->cd(2);
    gPad->SetTickx();
    gPad->SetTicky();
    hZvAllEv->Draw();
    gPad->Update();
    TPaveStats *stc=(TPaveStats*)hZvAllEv->GetListOfFunctions()->FindObject("stats");
    stc->SetY1NDC(0.7);
    stc->SetY2NDC(0.89);
    gPad->Modified();
    hZvEvSel->Draw("sames");
    gPad->Update();
    TPaveStats *std=(TPaveStats*)hZvEvSel->GetListOfFunctions()->FindObject("stats");
    std->SetY1NDC(0.5);
    std->SetY2NDC(0.69);
    std->SetTextColor(hZvEvSel->GetLineColor());
    gPad->Modified();
    cev->cd(4);
    gPad->SetTickx();
    gPad->SetTicky();
    hRatioZvEv->SetMinimum(0.85);
    hRatioZvEv->SetMaximum(1.01);
    hRatioZvEv->Draw();
    cev->SaveAs(Form("figures/FracEvSel%s.eps",suffix.Data()));
  }
  
  ComputeAndWriteEff(l,"Prompt");
  ComputeAndWriteEff(l,"Feeddw");

  TFile* outup=new TFile(Form("outputEff%s.root",suffix.Data()),"update");
  TH1D* hEffPr=(TH1D*)outup->Get("hEffPromptVsPtMultAndPtWeight");
  TH1D* hEffFd=(TH1D*)outup->Get("hEffFeeddwVsPtMultWeight");
  if(ptBWeight!="") hEffFd=(TH1D*)outup->Get("hEffFeeddwVsPtMultAndPtBWeight");
  TH1D* hAccPr=(TH1D*)outup->Get("hAccPromptVsPtMultAndPtWeight");
  TH1D* hAccFd=(TH1D*)outup->Get("hAccFeeddwVsPtMultWeight");
  if(ptBWeight!="") hAccFd=(TH1D*)outup->Get("hAccFeeddwVsPtMultAndPtBWeight");
  TH1D* hAxePr=(TH1D*)outup->Get("hAxePromptVsPtMultAndPtWeight");
  TH1D* hAxeFd=(TH1D*)outup->Get("hAxeFeeddwVsPtMultWeight");
  if(ptBWeight!="") hAxeFd=(TH1D*)outup->Get("hAxeFeeddwVsPtMultAndPtBWeight");

  hEffPr->SetLineColor(kRed+1);
  hEffFd->SetLineColor(kBlue+1);
  hEffPr->SetMarkerColor(kRed+1);
  hEffFd->SetMarkerColor(kBlue+1);
  hEffPr->SetMarkerStyle(20);
  hEffFd->SetMarkerStyle(25);
  hAccPr->SetLineColor(kRed+1);
  hAccFd->SetLineColor(kBlue+1);
  hAccPr->SetMarkerColor(kRed+1);
  hAccFd->SetMarkerColor(kBlue+1);
  hAccPr->SetMarkerStyle(20);
  hAccFd->SetMarkerStyle(25);
  
  TCanvas* capf=new TCanvas("capf","Acc Prompt vs Feeddown Toy vs MC",1200,600);
  capf->Divide(2,1);
  capf->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hAccPr->GetYaxis()->SetTitle("GenAcc / GenLimAcc");
  hAccPr->GetYaxis()->SetTitleOffset(1.3);
  hAccPr->SetStats(0);
  hAccFd->SetStats(0);
  hAccPr->DrawCopy();
  hAccFd->DrawCopy("same");
  hAccToy->DrawCopy("same");
  TLegend* legapf=new TLegend(0.5,0.16,0.89,0.36);
  legapf->AddEntry(Form("%s_copy",hAccPr->GetName()),"Prompt (full MC)","P");
  legapf->AddEntry(Form("%s_copy",hAccFd->GetName()),"Feeddown (full MC)","P");
  legapf->AddEntry(hAccToy,"Toy MC","P");
  legapf->Draw();  
  capf->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  TH1D* hRatioPrToy=(TH1D*)hAccPr->Clone("hRatioPrToy");
  hRatioPrToy->Divide(hAccToy);
  TH1D* hRatioFdToy=(TH1D*)hAccFd->Clone("hRatioFdToy");
  hRatioFdToy->Divide(hAccToy);
  hRatioPrToy->SetStats(0);
  hRatioFdToy->SetStats(0);
  hRatioPrToy->SetMinimum(0.98);
  hRatioPrToy->SetMaximum(1.02);
  hRatioPrToy->GetYaxis()->SetTitle("Ratio to Toy MC");
  hRatioPrToy->GetYaxis()->SetTitleOffset(1.2);
  hRatioPrToy->DrawCopy();
  hRatioFdToy->DrawCopy("same");

  
  TCanvas* cepf=new TCanvas("cepf","Eff, Prompt vs Feeddown",1200,600);
  cepf->Divide(2,1);
  cepf->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hEffPr->GetYaxis()->SetTitle("Reco / GenAcc");
  hEffPr->GetYaxis()->SetTitleOffset(1.3);
  hEffPr->SetStats(0);
  hEffFd->SetStats(0);
  hEffPr->SetMinimum(0.2);
  hEffPr->SetMaximum(0.6);
  hEffPr->DrawCopy();
  hEffFd->DrawCopy("same");
  TLegend* legepf=new TLegend(0.5,0.16,0.89,0.36);
  legepf->AddEntry(Form("%s_copy",hEffPr->GetName()),"Prompt","P");
  legepf->AddEntry(Form("%s_copy",hEffFd->GetName()),"Feeddown","P");
  legepf->Draw(); 
  cepf->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  TH1D* hRatioEff=(TH1D*)hEffFd->Clone("hRatioEff");
  hRatioEff->Divide(hEffPr);
  hRatioEff->SetStats(0);
  hRatioEff->SetMarkerColor(kGray+2);
  hRatioEff->SetLineColor(kGray+2);
  hRatioEff->SetMarkerStyle(20);
  hRatioEff->SetMinimum(0.9);
  hRatioEff->SetMaximum(1.04);
  hRatioEff->GetYaxis()->SetTitle("Ratio efficiency feeddown/prompt");
  hRatioEff->GetYaxis()->SetTitleOffset(1.3);
  hRatioEff->DrawCopy();
  
  TH1D* hEffD=0x0;
  TH1D* hEffB=0x0;
  TH1D* hRatioAxe=0x0;
  if(useToyMC){
    hEffD=(TH1D*)hEffPr->Clone("hEffD");
    hEffD->Reset("imen");
    hEffB=(TH1D*)hEffFd->Clone("hEffB");
    hEffB->Reset("imen");
    hRatioAxe=(TH1D*)hEffFd->Clone("hRatioAxe");
    hRatioAxe->Reset("imen");
    
    for(Int_t iBin=1; iBin<=hEffPr->GetNbinsX(); iBin++){
      Double_t lowPt=hEffPr->GetXaxis()->GetBinLowEdge(iBin);
      Double_t highPt=hEffPr->GetXaxis()->GetBinUpEdge(iBin);
      Double_t pt=hEffPr->GetBinCenter(iBin);
      
      Int_t theBinB=hEffFd->FindBin(pt);
      Double_t lowPtB=hEffFd->GetXaxis()->GetBinLowEdge(theBinB);
      Double_t highPtB=hEffFd->GetXaxis()->GetBinUpEdge(theBinB);
      
      Int_t theBinA=hAccToy->FindBin(pt);
      Double_t lowPtA=hAccToy->GetXaxis()->GetBinLowEdge(theBinA);
      Double_t highPtA=hAccToy->GetXaxis()->GetBinUpEdge(theBinA);
      
      if(TMath::Abs(lowPtA-lowPt)>0.001 ||
	 TMath::Abs(lowPtB-lowPt)>0.001 ||
	 TMath::Abs(highPtA-highPt)>0.001 ||
	 TMath::Abs(highPtB-highPt)>0.001){
	printf("ERROR IN BINNING: %f %f    %f %f   %f %f\n",lowPt,highPt,lowPtA,highPtA,lowPtB,highPtB);
	return;
      }
      
      Double_t effD=hEffPr->GetBinContent(iBin);
      Double_t acc=hAccToy->GetBinContent(theBinA);
      Double_t effB=hEffFd->GetBinContent(theBinB);
      
      Double_t effDerr=hEffPr->GetBinError(iBin);
      Double_t accerr=hAccToy->GetBinError(theBinA);
      Double_t effBerr=hEffFd->GetBinError(theBinB);
      
      Double_t acceffD=effD*acc;
      Double_t acceffB=effB*acc;
      Double_t acceffDerr=TMath::Sqrt(acc*acc*effDerr*effDerr+effD*effD*accerr*accerr);
      Double_t acceffBerr=TMath::Sqrt(acc*acc*effBerr*effBerr+effB*effB*accerr*accerr);
      
      Double_t r=-999.;
      Double_t er=0.;
      if(effD>0 && effB>0){
	r=effB/effD;
	er=r*TMath::Sqrt(effDerr/effD*effDerr/effD+effBerr/effB*effBerr/effB);
      }
      printf("Bin %d   acc=%f+-%f\n",iBin,acc,accerr);
      printf("         effD=%f+-%f   acc*effD=%f+-%f\n",effD,effDerr,acceffD,acceffDerr);
      printf("         effB=%f+-%f   acc*effB=%f+-%f\n",effB,effBerr,acceffB,acceffBerr);
      printf("         r=%f+-%f\n",r,er);
      hEffD->SetBinContent(iBin,acceffD);
      hEffD->SetBinError(iBin,acceffDerr);
      hEffB->SetBinContent(iBin,acceffB);
      hEffB->SetBinError(iBin,acceffBerr);
      hRatioAxe->SetBinContent(iBin,r);
      hRatioAxe->SetBinError(iBin,er);
    }
  }else{
    hEffD=(TH1D*)hAxePr->Clone("hEffD");
    hEffB=(TH1D*)hAxeFd->Clone("hEffB");
    hRatioAxe=(TH1D*)hAxeFd->Clone("hRatioAxe");
    hRatioAxe->Divide(hAxePr);
    for(Int_t iBin=1; iBin<=hAxePr->GetNbinsX(); iBin++){
      Double_t acceffD=hAxePr->GetBinContent(iBin);
      Double_t acceffB=hAxeFd->GetBinContent(iBin);
      Double_t acceffDerr=hAxePr->GetBinError(iBin);
      Double_t acceffBerr=hAxeFd->GetBinError(iBin);
      Double_t r=-999.;
      Double_t er=0.;
      if(acceffD>0 && acceffB>0){
	r=acceffB/acceffD;
	er=r*TMath::Sqrt(acceffDerr/acceffD*acceffDerr/acceffD+acceffBerr/acceffB*acceffBerr/acceffB);
      }
      printf("Bin %d   \n",iBin);
      printf("         acc*effD=%f+-%f\n",acceffD,acceffDerr);
      printf("         acc*effB=%f+-%f\n",acceffB,acceffBerr);
      printf("         r=%f+-%f\n",r,er);      
    }
  }
  hEffD->SetStats(0);
  hEffB->SetStats(0);

  TCanvas* cpf=new TCanvas("cpf","Axe Prompt vs Feeddown",1200,600);
  cpf->Divide(2,1);
  cpf->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hEffD->GetYaxis()->SetTitle("Acceptance x Efficiency");
  hEffD->GetYaxis()->SetTitleOffset(1.3);
  hEffD->SetMarkerStyle(20);
  hEffD->SetMarkerColor(kRed+1);
  hEffD->SetLineColor(kRed+1);
  hEffD->DrawCopy();
  hEffB->SetMarkerStyle(25);
  hEffB->SetMarkerColor(kBlue+1);
  hEffB->SetLineColor(kBlue+1);
  hEffB->DrawCopy("same");
  TLegend* legpf=new TLegend(0.6,0.16,0.89,0.36);
  legpf->AddEntry(Form("%s_copy",hEffD->GetName()),"Prompt","P");
  legpf->AddEntry(Form("%s_copy",hEffB->GetName()),"Feeddown","P");
  legpf->Draw();
  cpf->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hRatioAxe->SetStats(0);
  hRatioAxe->SetMarkerColor(kGray+2);
  hRatioAxe->SetLineColor(kGray+2);
  hRatioAxe->SetMarkerStyle(20);
  hRatioAxe->SetMinimum(0.9);
  hRatioAxe->SetMaximum(1.04);
  hRatioAxe->GetYaxis()->SetTitle("Ratio efficiency feeddown/prompt");
  hRatioAxe->GetYaxis()->SetTitleOffset(1.3);
  hRatioAxe->DrawCopy();
  cpf->SaveAs(Form("figures/EfficVsPt_PromptFd_%s.eps",suffix.Data()));

  outup->cd();  
  hEffD->Write();
  hEffB->Write();
  outup->Close();



}

void ComputeAndWriteEff(TList* l, TString dCase){

  TString var1="Pt";
  TString var2="Y";
  TString var3="Mult";
  TString yTitle="y";
  TString zTitle="N_{tracklets} in |#eta|<1";
  if(dCase=="Feeddw"){
    var2="PtB";
    yTitle="B-hadron p_{T} (GeV/c)";
    TH3F* h4test=(TH3F*)l->FindObject(Form("h%sVs%sVs%sReco%s",var1.Data(),var2.Data(),var3.Data(),dCase.Data()));
    if(!h4test){
      printf("WARNING: pt(B) and mult weights cannot be both applied because histos are not there -> re-run with more recent version of the task\n");
      printf(" ---> Try pt(B) only weights\n");
      var2="Y";
      var3="PtB";
      yTitle="y";
      zTitle="B-hadron p_{T} (GeV/c)";
      TH3F* h4test2=(TH3F*)l->FindObject(Form("h%sVs%sVs%sReco%s",var1.Data(),var2.Data(),var3.Data(),dCase.Data()));
      if(!h4test2){
	printf("WARNING: pt(B) weights cannot be applied because histos are not there -> re-run with more recent version of the task\n");
	printf(" ---> use only mult weights\n");
	var3="Mult";
	zTitle="N_{tracklets} in |#eta|<1";
      }
    }
  }
  printf("Case: %s  Variables used for projections: %s %s %s\n",dCase.Data(),var1.Data(),var2.Data(),var3.Data());
  TH3F* h3dForAcc=(TH3F*)l->FindObject(Form("hPtVsYVsMultGenAcc%s",dCase.Data()));
  TH3F* h3dForLimAcc=(TH3F*)l->FindObject(Form("hPtVsYVsMultGenLimAcc%s",dCase.Data()));
  TH3F* h3dReco=(TH3F*)l->FindObject(Form("h%sVs%sVs%sReco%s",var1.Data(),var2.Data(),var3.Data(),dCase.Data()));
  TH3F* h3dGenAccEvSel=(TH3F*)l->FindObject(Form("h%sVs%sVs%sGenAccEvSel%s",var1.Data(),var2.Data(),var3.Data(),dCase.Data()));
  TH3F* h3dGenAcc=(TH3F*)l->FindObject(Form("h%sVs%sVs%sGenAcc%s",var1.Data(),var2.Data(),var3.Data(),dCase.Data()));
  TH3F* h3dGenLimAcc=(TH3F*)l->FindObject(Form("h%sVs%sVs%sGenLimAcc%s",var1.Data(),var2.Data(),var3.Data(),dCase.Data()));
  
  TH2D* hypt=(TH2D*)h3dForAcc->Project3D("yx");
  hypt->SetTitle(Form("Generated in acceptance, all multiplicities, %s",dCase.Data()));
  TH2D* hptmult=(TH2D*)h3dForLimAcc->Project3D("xz");
  hptmult->SetTitle(Form("Generated in |y|<0.5, %s",dCase.Data()));
  hptmult->SetStats(0);
  hypt->SetStats(0);
  hypt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hypt->GetYaxis()->SetTitle("y");
  hptmult->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  hptmult->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
  if(maxMult>0) hptmult->GetXaxis()->SetRangeUser(0.,maxMult);
  TProfile* hMeanPtMult=hptmult->ProfileX(Form("hMeanPtMult%s",dCase.Data()));
  TH2D* hptbmult=0x0;
  TProfile* hMeanPtbMult=0x0;
  Int_t cwid=1200;
  Int_t npads=2;
  if(dCase=="Feeddw" && var2=="PtB"){
    hptbmult=(TH2D*)h3dGenLimAcc->Project3D("yz");
    hptbmult->SetTitle(Form("Generated in |y|<0.5, %s",dCase.Data()));
    hptbmult->SetStats(0);
    hptbmult->GetYaxis()->SetTitle("B-hadron p_{T} (GeV/c)");
    hptbmult->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
    hMeanPtbMult=hptbmult->ProfileX(Form("hMeanPtbMult%s",dCase.Data()));
    cwid=1500;
    npads=3;
  }
				      
     
  TCanvas* c0=new TCanvas(Form("c0%s",dCase.Data()),Form("%s - 2D plots",dCase.Data()),cwid,600);
  c0->Divide(npads,1);
  c0->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  hypt->Draw("colz");
  c0->cd(2);
  gPad->SetLogz();
  gPad->SetRightMargin(0.13);
  hptmult->Draw("colz");
  hMeanPtMult->Draw("same");
  if(npads==3){
    c0->cd(3);
    gPad->SetLogz();
    hptbmult->Draw("colz");
    hMeanPtbMult->Draw("same");
  }
  
  TH1D* hPtReco=h3dReco->ProjectionX(Form("hPtReco%s",dCase.Data()));
  TH1D* hPtGenAcc=h3dGenAcc->ProjectionX(Form("hPtGenAcc%s",dCase.Data()));
  TH1D* hPtGenAccEvSel=0x0;
  if(h3dGenAccEvSel) hPtGenAccEvSel=h3dGenAccEvSel->ProjectionX(Form("hPtGenAccEvSel%s",dCase.Data()));
  TH1D* hPtGenLimAcc=h3dGenLimAcc->ProjectionX(Form("hPtGenLimAcc%s",dCase.Data()));
  TH1D* hPtRecoR=(TH1D*)hPtReco->Rebin(nPtBins,Form("hPtReco%sReb",dCase.Data()),binLims);
  TH1D* hPtGenAccR=(TH1D*)hPtGenAcc->Rebin(nPtBins,Form("hPtGenAcc%sReb",dCase.Data()),binLims);
  TH1D* hPtGenLimAccR=(TH1D*)hPtGenLimAcc->Rebin(nPtBins,Form("hPtGenLimAcc%sReb",dCase.Data()),binLims);
  TH1D* hEffVsPt=(TH1D*)hPtReco->Clone(Form("hEff%s",dCase.Data()));
  hEffVsPt->Divide(hPtReco,hPtGenAcc,1,1,"B");
  TH1D* hEffVsPtR=(TH1D*)hPtRecoR->Clone(Form("hEff%sRebin",dCase.Data()));
  hEffVsPtR->Divide(hPtRecoR,hPtGenAccR,1,1,"B");
  hEffVsPt->SetStats(0);
  hEffVsPtR->SetStats(0);
  TH1D* hAccVsPt=(TH1D*)hPtGenAcc->Clone(Form("hAcc%s",dCase.Data()));
  hAccVsPt->Divide(hPtGenAcc,hPtGenLimAcc,1,1,"B");
  TH1D* hAccVsPtR=(TH1D*)hPtGenAccR->Clone(Form("hAcc%sReb",dCase.Data()));
  hAccVsPtR->Divide(hPtGenAccR,hPtGenLimAccR,1,1,"B");
  hAccVsPt->SetStats(0);
  hAccVsPtR->SetStats(0);
  TH1D* hRatioAccToyMC=(TH1D*)hAccVsPt->Clone("hRatioAccToyMC");
  hRatioAccToyMC->Reset("IMEN");
  for(Int_t j=1; j<hRatioAccToyMC->GetNbinsX(); j++){
    Double_t pt=hAccVsPt->GetBinCenter(j);
    Double_t pt1=hAccVsPt->GetBinLowEdge(j);
    Double_t pt2=hAccVsPt->GetBinLowEdge(j+1);
    Int_t jbinF=hAccToyFine->FindBin(pt);    
    Double_t ptToy1=hAccToyFine->GetBinLowEdge(jbinF);
    Double_t ptToy2=hAccToyFine->GetBinLowEdge(jbinF+1);
    if(TMath::Abs(ptToy1-pt1)<0.001 &&
       TMath::Abs(ptToy2-pt2)<0.001){
      Double_t r=hAccVsPt->GetBinContent(j)/hAccToyFine->GetBinContent(jbinF);
      Double_t er=r*TMath::Sqrt(hAccVsPt->GetBinError(j)/hAccVsPt->GetBinContent(j)*hAccVsPt->GetBinError(j)/hAccVsPt->GetBinContent(j)+hAccToyFine->GetBinError(jbinF)/hAccToyFine->GetBinContent(jbinF)*hAccToyFine->GetBinError(jbinF)/hAccToyFine->GetBinContent(jbinF));
      hRatioAccToyMC->SetBinContent(j,r);
      hRatioAccToyMC->SetBinError(j,er);     
    }
 }

  
  TH1D* hRatioAccToyMCR=(TH1D*)hAccVsPtR->Clone("hRatioAccToyMCR");
  hRatioAccToyMCR->Divide(hAccToy);
  
  TH1D* hEvSelEffVsPt=0x0;
  if(hPtGenAccEvSel){
    hEvSelEffVsPt=(TH1D*)hPtGenAccEvSel->Clone(Form("hEvSelEff%s",dCase.Data()));
    hEvSelEffVsPt->Divide(hPtGenAccEvSel,hPtGenAcc,1,1,"B");
    hEvSelEffVsPt->SetStats(0);
  }


  TCanvas* c1a=new TCanvas(Form("c1a%s",dCase.Data()),Form("%s - AccVsPt",dCase.Data()),1200,800);
  c1a->Divide(2,2);
  c1a->cd(1);
  gPad->SetLogy();
  hPtGenLimAcc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtGenLimAcc->GetYaxis()->SetTitle("Entries");
  hPtGenLimAcc->SetLineColor(1);
  hPtGenLimAcc->Draw();
  gPad->Update();
  TPaveStats *st1=(TPaveStats*)hPtGenLimAcc->GetListOfFunctions()->FindObject("stats");
  st1->SetY1NDC(0.7);
  st1->SetY2NDC(0.89);
  hPtGenAcc->SetLineColor(2);
  hPtGenAcc->Draw("sames");
  gPad->Update();
  TPaveStats *st2=(TPaveStats*)hPtGenAcc->GetListOfFunctions()->FindObject("stats");
  st2->SetY1NDC(0.5);
  st2->SetY2NDC(0.69);
  st2->SetTextColor(2);
  if(hPtGenAccEvSel){
    hPtGenAccEvSel->SetLineColor(6);
    hPtGenAccEvSel->Draw("sames");
    gPad->Update();
    TPaveStats *st2s=(TPaveStats*)hPtGenAccEvSel->GetListOfFunctions()->FindObject("stats");
    st2s->SetY1NDC(0.3);
    st2s->SetY2NDC(0.49);
    st2s->SetTextColor(6);
  }
  hPtReco->SetLineColor(4);
  hPtReco->Draw("sames");
  gPad->Update();
  TPaveStats *st13=(TPaveStats*)hPtReco->GetListOfFunctions()->FindObject("stats");
  st13->SetY1NDC(0.1);
  st13->SetY2NDC(0.29);
  st13->SetTextColor(4);
  gPad->Modified();
  c1a->cd(2);
  hEffVsPt->SetLineColor(4);
  hEffVsPt->SetMinimum(0);
  hEffVsPt->SetMaximum(1.95);
  hEffVsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hEffVsPt->GetYaxis()->SetTitle("Ratio");
  hEffVsPt->Draw();
  hAccVsPt->SetLineColor(2);
  hAccVsPt->Draw("same");
  hAccVsPtR->SetMarkerStyle(20);
  hAccVsPtR->SetMarkerColor(kRed+1);
  hAccVsPtR->SetLineColor(kRed+1);
  hAccVsPtR->Draw("same");
  hAccToyFine->Draw("same");
  hAccToy->Draw("same");
  TLatex* tacc=new TLatex(0.16,0.83,"Acceptance (CombinHF)");
  tacc->SetNDC();
  tacc->SetTextColor(hAccVsPt->GetLineColor());
  tacc->Draw();
  TLatex* tacct=new TLatex(0.16,0.76,"Acceptance (ToyMC)");
  tacct->SetNDC();
  tacct->SetTextColor(hAccToyFine->GetLineColor());
  tacct->Draw();
  TLatex* te=new TLatex(0.16,0.69,"Efficiency");
  te->SetNDC();
  te->SetTextColor(hEffVsPt->GetLineColor());
  te->Draw();
  c1a->cd(3);
  TLatex* t1=new TLatex(6.,0.25,"GenAccInEvSel / GenAcc");
  t1->SetTextColor(6);
  if(hPtGenAccEvSel){
    hEvSelEffVsPt->SetLineColor(6);
    hEvSelEffVsPt->SetMinimum(0);
    hEvSelEffVsPt->SetMaximum(1.05);
    hEvSelEffVsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hEvSelEffVsPt->GetYaxis()->SetTitle("Ratio");
    hEvSelEffVsPt->Draw();
    t1->Draw();  
  }
  hEffVsPt->Draw("same");
  TLatex* t3=new TLatex(6.,0.15,"Reco / GenAcc");
  t3->SetTextColor(4);
  t3->Draw();
  c1a->cd(4);
  hRatioAccToyMC->SetLineColor(kGray+1);
  hRatioAccToyMCR->SetLineColor(1);
  hRatioAccToyMCR->SetMarkerColor(1);
  hRatioAccToyMCR->SetMarkerStyle(21);
  hRatioAccToyMC->SetMinimum(0.95);
  hRatioAccToyMC->SetMaximum(1.05);
  hRatioAccToyMC->GetYaxis()->SetTitle("Acceptance(CombinHF) / Acceptance(ToyMC)");
  hRatioAccToyMC->Draw();
  hRatioAccToyMCR->Draw("same");
  c1a->SaveAs(Form("figures/EfficVsPt_%s_%s.eps",suffix.Data(),dCase.Data()));

  TH1D* hVar2RecoAllPt=h3dReco->ProjectionY(Form("hVar2Reco%s",dCase.Data()));
  TH1D* hVar2GenAccAllPt=h3dGenAcc->ProjectionY(Form("hVar2GenAcc%s",dCase.Data()));
  TH1D* hVar2GenAccEvSelAllPt=0x0;
  if(h3dGenAccEvSel) hVar2GenAccEvSelAllPt=h3dGenAccEvSel->ProjectionY(Form("hVar2GenAccEvSel%s",dCase.Data()));
  TH1D* hVar2GenLimAccAllPt=h3dGenLimAcc->ProjectionY(Form("hVar2GenLimAcc%s",dCase.Data()));
  TH1D* hEffVsVar2AllPt=(TH1D*)hVar2RecoAllPt->Clone(Form("hEff%s",dCase.Data()));
  hEffVsVar2AllPt->Divide(hVar2RecoAllPt,hVar2GenAccAllPt,1,1,"B");
  TH1D* hAccEffVsVar2AllPt=(TH1D*)hVar2RecoAllPt->Clone(Form("hAccEff%s",dCase.Data()));
  hAccEffVsVar2AllPt->Divide(hVar2RecoAllPt,hVar2GenLimAccAllPt,1,1,"B");
  TH1D* hAccVsVar2AllPt=(TH1D*)hVar2GenAccAllPt->Clone(Form("hAcc%s",dCase.Data()));
  hAccVsVar2AllPt->Divide(hVar2GenAccAllPt,hVar2GenLimAccAllPt,1,1,"B");
  hEffVsVar2AllPt->SetStats(0);
  hAccVsVar2AllPt->SetStats(0);
  hAccEffVsVar2AllPt->SetStats(0);
  TH1D* hEvSelEffVsVar2AllPt=0x0;
  if(hVar2GenAccEvSelAllPt){
    hEvSelEffVsVar2AllPt=(TH1D*)hVar2GenAccEvSelAllPt->Clone(Form("hEvSelEff%s%s",var2.Data(),dCase.Data()));
    hEvSelEffVsVar2AllPt->Divide(hVar2GenAccEvSelAllPt,hVar2GenAccAllPt,1,1,"B");
    hEvSelEffVsVar2AllPt->SetStats(0);
  }
  
  TCanvas* c2a=new TCanvas(Form("c2a%s",dCase.Data()),Form("%s - AccVs%s",dCase.Data(),var2.Data()),1500,600);
  c2a->Divide(3,1);
  c2a->cd(1);
  gPad->SetLogy();
  hVar2GenLimAccAllPt->SetLineColor(1);
  hVar2GenLimAccAllPt->GetXaxis()->SetTitle(yTitle.Data());
  hVar2GenLimAccAllPt->GetYaxis()->SetTitle("Entries");
  if(var2=="Y") hVar2GenLimAccAllPt->SetMinimum(0.9);
  hVar2GenLimAccAllPt->Draw();
  gPad->Update();
  TPaveStats *st21=(TPaveStats*)hVar2GenLimAccAllPt->GetListOfFunctions()->FindObject("stats");
  st21->SetY1NDC(0.7);
  st21->SetY2NDC(0.89);
  hVar2GenAccAllPt->SetLineColor(2);
  hVar2GenAccAllPt->Draw("sames");
  gPad->Update();
  TPaveStats *st22=(TPaveStats*)hVar2GenAccAllPt->GetListOfFunctions()->FindObject("stats");
  st22->SetY1NDC(0.5);
  st22->SetY2NDC(0.69);
  st22->SetTextColor(2);
  if(hVar2GenAccEvSelAllPt){
    hVar2GenAccEvSelAllPt->SetLineColor(6);
    hVar2GenAccEvSelAllPt->Draw("sames");
    gPad->Update();
    TPaveStats *st22s=(TPaveStats*)hVar2GenAccEvSelAllPt->GetListOfFunctions()->FindObject("stats");
    st22s->SetY1NDC(0.3);
    st22s->SetY2NDC(0.49);
    st22s->SetTextColor(6);
  }
  hVar2RecoAllPt->SetLineColor(4);
  hVar2RecoAllPt->Draw("sames");
  gPad->Update();
  TPaveStats *st23=(TPaveStats*)hVar2RecoAllPt->GetListOfFunctions()->FindObject("stats");
  st23->SetY1NDC(0.1);
  st23->SetY2NDC(0.29);
  st23->SetTextColor(4);
  gPad->Modified();
  c2a->cd(2);
  hEffVsVar2AllPt->SetLineColor(4);
  hEffVsVar2AllPt->SetMinimum(0);
  hEffVsVar2AllPt->SetMaximum(1.6);
  hEffVsVar2AllPt->GetXaxis()->SetTitle(yTitle.Data());
  hEffVsVar2AllPt->GetYaxis()->SetTitle("Ratio");
  hEffVsVar2AllPt->Draw();
  hAccVsVar2AllPt->SetLineColor(2);
  hAccVsVar2AllPt->Draw("same");
  //  hAccEffVsVar2AllPt->SetLineColor(6);
  // hAccEffVsVar2AllPt->Draw("same");
  TLatex* tacc2=new TLatex(0.16,0.8,"Acceptance (CombinHF)");
  tacc2->SetNDC();
  tacc2->SetTextColor(hAccVsVar2AllPt->GetLineColor());
  tacc2->Draw();
  TLatex* te2=new TLatex(0.16,0.72,"Efficiency");
  te2->SetNDC();
  te2->SetTextColor(hEffVsVar2AllPt->GetLineColor());
  te2->Draw();
  c2a->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  if(hVar2GenAccEvSelAllPt){
    hEvSelEffVsVar2AllPt->SetLineColor(6);
    hEvSelEffVsVar2AllPt->SetMinimum(0);
    hEvSelEffVsVar2AllPt->SetMaximum(1.05);
    hEvSelEffVsVar2AllPt->GetXaxis()->SetTitle(yTitle.Data());
    hEvSelEffVsVar2AllPt->GetYaxis()->SetTitle("Ratio");
    hEvSelEffVsVar2AllPt->Draw();
    t1->Draw();
  }
  hEffVsVar2AllPt->SetLineColor(4);
  hEffVsVar2AllPt->Draw("same");
  t3->Draw();
  c2a->SaveAs(Form("figures/EfficVs%s_%s_%s.eps",var2.Data(),suffix.Data(),dCase.Data()));

  
  TH1D* hVar3RecoAllPt=h3dReco->ProjectionZ(Form("hVar3Reco%s",dCase.Data()));
  TH1D* hVar3GenAccAllPt=h3dGenAcc->ProjectionZ(Form("hVar3GenAcc%s",dCase.Data()));
  TH1D* hVar3GenAccEvSelAllPt=0x0;
  if(h3dGenAccEvSel) hVar3GenAccEvSelAllPt=h3dGenAccEvSel->ProjectionZ(Form("hVar3GenAccEvSel%s",dCase.Data()));
  TH1D* hVar3GenLimAccAllPt=h3dGenLimAcc->ProjectionZ(Form("hVar3GenLimAcc%s",dCase.Data()));
  TH1D* hEffVsVar3AllPt=(TH1D*)hVar3RecoAllPt->Clone(Form("hEff%s",dCase.Data()));
  hEffVsVar3AllPt->Divide(hVar3RecoAllPt,hVar3GenAccAllPt,1,1,"B");
  TH1D* hAccEffVsVar3AllPt=(TH1D*)hVar3RecoAllPt->Clone(Form("hAccEff%s",dCase.Data()));
  hAccEffVsVar3AllPt->Divide(hVar3RecoAllPt,hVar3GenLimAccAllPt,1,1,"B");
  TH1D* hAccVsVar3AllPt=(TH1D*)hVar3GenAccAllPt->Clone(Form("hAcc%s",dCase.Data()));
  hAccVsVar3AllPt->Divide(hVar3GenAccAllPt,hVar3GenLimAccAllPt,1,1,"B");
  hEffVsVar3AllPt->SetStats(0);
  hAccVsVar3AllPt->SetStats(0);
  hAccEffVsVar3AllPt->SetStats(0);
  TH1D* hEvSelEffVsVar3AllPt=0x0;
  if(hVar3GenAccEvSelAllPt){
    hEvSelEffVsVar3AllPt=(TH1D*)hVar3GenAccEvSelAllPt->Clone(Form("hEvSelEff%s%s",var3.Data(),dCase.Data()));
    hEvSelEffVsVar3AllPt->Divide(hVar3GenAccEvSelAllPt,hVar3GenAccAllPt,1,1,"B");
    hEvSelEffVsVar3AllPt->SetStats(0);
  }

  TCanvas* c3a=new TCanvas(Form("c3a%s",dCase.Data()),Form("%s - AccVs%s",dCase.Data(),var3.Data()),1500,600);
  c3a->Divide(3,1);
  c3a->cd(1);
  gPad->SetLogy();
  hVar3GenLimAccAllPt->SetLineColor(1);
  hVar3GenLimAccAllPt->GetXaxis()->SetTitle(zTitle.Data());
  hVar3GenLimAccAllPt->GetYaxis()->SetTitle("Entries");
  if(var3=="Mult" && maxMult>0) hVar3GenLimAccAllPt->GetXaxis()->SetRangeUser(0.,maxMult);
  hVar3GenLimAccAllPt->Draw();
  gPad->Update();
  TPaveStats *st31=(TPaveStats*)hVar3GenLimAccAllPt->GetListOfFunctions()->FindObject("stats");
  st31->SetY1NDC(0.7);
  st31->SetY2NDC(0.89);
  hVar3GenAccAllPt->SetLineColor(2);
  hVar3GenAccAllPt->Draw("sames");
  gPad->Update();
  TPaveStats *st32=(TPaveStats*)hVar3GenAccAllPt->GetListOfFunctions()->FindObject("stats");
  st32->SetY1NDC(0.5);
  st32->SetY2NDC(0.69);
  st32->SetTextColor(2);
  if(hVar3GenAccEvSelAllPt){
    hVar3GenAccEvSelAllPt->SetLineColor(6);
    hVar3GenAccEvSelAllPt->Draw("sames");
    gPad->Update();
    TPaveStats *st32s=(TPaveStats*)hVar3GenAccEvSelAllPt->GetListOfFunctions()->FindObject("stats");
    st32s->SetY1NDC(0.3);
    st32s->SetY2NDC(0.49);
    st32s->SetTextColor(6);
  }
  hVar3RecoAllPt->SetLineColor(4);
  hVar3RecoAllPt->Draw("sames");
  gPad->Update();
  TPaveStats *st33=(TPaveStats*)hVar3RecoAllPt->GetListOfFunctions()->FindObject("stats");
  st33->SetY1NDC(0.1);
  st33->SetY2NDC(0.29);
  st33->SetTextColor(4);
  gPad->Modified();
  c3a->cd(2);
  hEffVsVar3AllPt->SetLineColor(4);
  hEffVsVar3AllPt->SetMinimum(0);
  hEffVsVar3AllPt->SetMaximum(1.6);
  hEffVsVar3AllPt->GetXaxis()->SetTitle(zTitle.Data());
  hEffVsVar3AllPt->GetYaxis()->SetTitle("Ratio");
  if(var3=="Mult" && maxMult>0) hEffVsVar3AllPt->GetXaxis()->SetRangeUser(0.,maxMult);
  hEffVsVar3AllPt->Draw();
  hAccVsVar3AllPt->SetLineColor(2);
  hAccVsVar3AllPt->Draw("same");
  //  hAccEffVsVar3AllPt->SetLineColor(6);
  // hAccEffVsVar3AllPt->Draw("same");
  TLatex* tacc3=new TLatex(0.16,0.8,"Acceptance (CombinHF)");
  tacc3->SetNDC();
  tacc3->SetTextColor(hAccVsVar3AllPt->GetLineColor());
  tacc3->Draw();
  TLatex* te3=new TLatex(0.16,0.72,"Efficiency");
  te3->SetNDC();
  te3->SetTextColor(hEffVsVar3AllPt->GetLineColor());
  te3->Draw();
  c3a->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  if(hVar3GenAccEvSelAllPt){
    hEvSelEffVsVar3AllPt->SetLineColor(6);
    hEvSelEffVsVar3AllPt->SetMinimum(0);
    hEvSelEffVsVar3AllPt->SetMaximum(1.05);
    hEvSelEffVsVar3AllPt->GetXaxis()->SetTitle(zTitle.Data());
    hEvSelEffVsVar3AllPt->GetYaxis()->SetTitle("Ratio");
    if(var3=="Mult" && maxMult>0) hEvSelEffVsVar3AllPt->GetXaxis()->SetRangeUser(0.,maxMult);
    hEvSelEffVsVar3AllPt->Draw();
    t1->Draw();
  }
  hEffVsVar3AllPt->SetLineColor(4);
  hEffVsVar3AllPt->Draw("same");
  t3->Draw();
  c3a->SaveAs(Form("figures/EfficVs%s_%s_%s.eps",var3.Data(),suffix.Data(),dCase.Data()));

  TH1D* hEffVsMult[nPtBins];

  if(var3=="Mult"){
    // double differential acceptance plot pt/mult
    Int_t colMult[20]={kMagenta+1,kMagenta,kBlue+1,kBlue,kBlue-9,
		       kGreen+2,kGreen+1,kGreen,kYellow+1,kYellow,
		       kOrange+2,kOrange+1,kRed-9,kRed,kRed+1};
    TH1D* hAccVsPtMultBin[300];
    TH1D* hPtGenLimAccMultBin[300];
    TH1D* hPtGenAccMultBin[300];
    TGraphErrors* gMeanPtGenLimAccVsMult=new TGraphErrors(0);
    TGraphErrors* gMeanPtGenAccVsMult=new TGraphErrors(0);
    TCanvas* c2dch=new TCanvas(Form("c2dch%s",dCase.Data()),Form("%s - AccVsPt and %s",dCase.Data(),var3.Data()),1400,800);
    c2dch->Divide(3,2);
    Int_t nhp=0;
    TLegend* leg=new TLegend(0.1,0.1,0.6,0.9);
    for(Int_t iBinm=0; iBinm<h3dGenAcc->GetNbinsZ(); iBinm++){
      hPtGenAccMultBin[iBinm]=(TH1D*)h3dGenAcc->ProjectionX(Form("hPtGenAccMultBin%d",iBinm),0,-1,iBinm+1,iBinm+1);
      hPtGenLimAccMultBin[iBinm]=(TH1D*)h3dGenLimAcc->ProjectionX(Form("hPtGenLimAccMultBin%d",iBinm),0,-1,iBinm+1,iBinm+1);
      hAccVsPtMultBin[iBinm]=(TH1D*)hPtGenAccMultBin[iBinm]->Clone(Form("hAccVsPtMultBin%d",iBinm));
      hAccVsPtMultBin[iBinm]->Divide(hPtGenAccMultBin[iBinm],hPtGenLimAccMultBin[iBinm],1,1,"B");
      Double_t minmul=h3dGenAcc->GetZaxis()->GetBinLowEdge(iBinm+1);
      Double_t maxmul=h3dGenAcc->GetZaxis()->GetBinUpEdge(iBinm+1);
      Double_t centmul=0.5*(minmul+maxmul);
      Double_t meanptlimacc=hPtGenLimAccMultBin[iBinm]->GetMean();
      Double_t meanptacc=hPtGenAccMultBin[iBinm]->GetMean();
      Double_t emeanptlimacc=hPtGenLimAccMultBin[iBinm]->GetMeanError();
      Double_t emeanptacc=hPtGenAccMultBin[iBinm]->GetMeanError();
      hPtGenAccMultBin[iBinm]->SetTitle("GenAcc");
      hPtGenLimAccMultBin[iBinm]->SetTitle("GenLimAcc");
      hAccVsPtMultBin[iBinm]->SetTitle(" ");
      hPtGenAccMultBin[iBinm]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hPtGenLimAccMultBin[iBinm]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hAccVsPtMultBin[iBinm]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hAccVsPtMultBin[iBinm]->GetYaxis()->SetTitle("Acceptance");
      hPtGenLimAccMultBin[iBinm]->GetYaxis()->SetTitle("Counts (normalized)");
      hPtGenAccMultBin[iBinm]->GetYaxis()->SetTitle("Counts (normalized)");
    
      hPtGenAccMultBin[iBinm]->SetStats(0);
      hPtGenLimAccMultBin[iBinm]->SetStats(0);
      hAccVsPtMultBin[iBinm]->SetStats(0);
      if(hPtGenLimAccMultBin[iBinm]->GetEntries()>1000 && nhp<20){
	gMeanPtGenLimAccVsMult->SetPoint(nhp,centmul,meanptlimacc);
	gMeanPtGenLimAccVsMult->SetPointError(nhp,centmul-minmul,emeanptlimacc);
	gMeanPtGenAccVsMult->SetPoint(nhp,centmul,meanptacc);
	gMeanPtGenAccVsMult->SetPointError(nhp,centmul-minmul,emeanptacc);
	hPtGenLimAccMultBin[iBinm]->SetLineColor(colMult[nhp]);
	hPtGenAccMultBin[iBinm]->SetLineColor(colMult[nhp]);
	hAccVsPtMultBin[iBinm]->SetLineColor(colMult[nhp]);
	leg->AddEntry(hPtGenLimAccMultBin[iBinm],Form("%.0f<mult<%.0f",minmul,maxmul),"L")->SetTextColor(colMult[nhp]);
	++nhp;
	c2dch->cd(2);
	gPad->SetLogy();
	if(iBinm==0) hPtGenLimAccMultBin[iBinm]->Draw();
	else hPtGenLimAccMultBin[iBinm]->DrawNormalized("same");
	c2dch->cd(3);
	gPad->SetLogy();
	if(iBinm==0) hPtGenAccMultBin[iBinm]->Draw();
	else hPtGenAccMultBin[iBinm]->DrawNormalized("same");
	c2dch->cd(5);
	if(iBinm==0){
	  hAccVsPtMultBin[iBinm]->SetMinimum(0.2);
	  hAccVsPtMultBin[iBinm]->SetMaximum(2.2);
	  hAccVsPtMultBin[iBinm]->Draw();
	}
	else hAccVsPtMultBin[iBinm]->Draw("same");
      }
    }
    c2dch->cd(6);
    leg->Draw();
    c2dch->cd(4);
    gMeanPtGenLimAccVsMult->SetMarkerStyle(20);
    gMeanPtGenLimAccVsMult->SetTitle(" ");
    gMeanPtGenLimAccVsMult->GetXaxis()->SetTitle("N_{tracklets}");
    gMeanPtGenLimAccVsMult->GetYaxis()->SetTitle("<p_{T}> (Gev/c)");
    gMeanPtGenLimAccVsMult->SetMinimum(1.5);
    gMeanPtGenLimAccVsMult->SetMaximum(5.);
    gMeanPtGenLimAccVsMult->Draw("AP");
    gMeanPtGenAccVsMult->SetMarkerStyle(25);
    gMeanPtGenAccVsMult->Draw("PSAME");
    TLegend* legmp=new TLegend(0.15,0.7,0.4,0.89);
    legmp->AddEntry(gMeanPtGenLimAccVsMult,"GenLimAcc","P");
    legmp->AddEntry(gMeanPtGenAccVsMult,"GenAcc","P");
    legmp->Draw();
    c2dch->cd(1);
    TH1D* hCopy=(TH1D*)hAccVsVar3AllPt->Clone("hAccVsVar3AllPtCopy");
    hCopy->SetLineColor(1);
    hCopy->SetMinimum(0.55);
    hCopy->SetMaximum(1.05);
    hCopy->SetMarkerStyle(22);
    hCopy->SetTitle(" ");
    hCopy->GetXaxis()->SetTitle("N_{tracklets}");
    hCopy->GetYaxis()->SetTitle("Acceptance (p_{T} integrated)");
    hCopy->Draw();

    TH1D* hAccVsMultPtBin[10];
    for(Int_t iBinp=0; iBinp<10; iBinp++){
      Double_t minptBin=0.5+iBinp;
      Double_t maxptBin=minptBin+0.1;
      Int_t binMinPt=h3dGenAcc->GetXaxis()->FindBin(minptBin+0.00001);
      Int_t binMaxPt=h3dGenAcc->GetXaxis()->FindBin(maxptBin-0.00001);
      TH1D* htmpm1=(TH1D*)h3dGenAcc->ProjectionZ(Form("hMultGenAccPtBin%d",iBinp),binMinPt,binMaxPt,0,-1);
      TH1D* htmpm2=(TH1D*)h3dGenLimAcc->ProjectionZ(Form("hMultGenLimAccPtBin%d",iBinp),binMinPt,binMaxPt,0,-1);
      hAccVsMultPtBin[iBinp]=(TH1D*)htmpm1->Clone(Form("hAccVsMultPtBin%d",iBinp));
      hAccVsMultPtBin[iBinp]->Divide(htmpm1,htmpm2,1,1,"B");
      minptBin=h3dGenAcc->GetXaxis()->GetBinLowEdge(binMinPt);
      maxptBin=h3dGenAcc->GetXaxis()->GetBinUpEdge(binMaxPt);
      hAccVsMultPtBin[iBinp]->SetTitle(Form("%.2f<pt<%.2f",minptBin,maxptBin));
      hAccVsMultPtBin[iBinp]->GetXaxis()->SetTitle("N_{tracklets}");
      hAccVsMultPtBin[iBinp]->GetYaxis()->SetTitle("Acceptance");
      hAccVsMultPtBin[iBinp]->SetStats(0);
      hAccVsMultPtBin[iBinp]->SetLineWidth(2);
      hAccVsMultPtBin[iBinp]->GetXaxis()->SetRangeUser(1000.,5000.);
    }
    TCanvas* c2dch2=new TCanvas(Form("c2dch2%s",dCase.Data()),Form("%s - AccVs%s and Pt",dCase.Data(),var3.Data()),1400,800);
    c2dch2->Divide(5,2);
    for(Int_t iBinp=0; iBinp<10; iBinp++){
      c2dch2->cd(iBinp+1);
      hAccVsMultPtBin[iBinp]->Draw();
    }
  
    const Int_t nMultBins=6;
    Double_t mulLims[nMultBins+1]={0.,5.,12.,20.,40.,80.,200.};
    if(h3dReco->GetZaxis()->GetXmax()>300. && maxMult<0){
      for(Int_t jb=0; jb<=nMultBins; jb++) mulLims[jb]=(Double_t)jb/nMultBins*h3dReco->GetZaxis()->GetXmax();
    }
    Int_t mulcol[nMultBins]={1,kRed+1,kGreen+2,4,kOrange+2,kMagenta};

    TH1D* hPtRecoM[nMultBins];
    TH1D* hPtGenAccM[nMultBins];
    //  TH1D* hPtGenLimAccM[nMultBins];
    TH1D* hEffVsPtM[nMultBins];
    for(Int_t j=0; j<nMultBins; j++){
      Int_t lowBin=h3dReco->GetZaxis()->FindBin(mulLims[j]);
      Int_t hiBin=h3dReco->GetZaxis()->FindBin(mulLims[j+1]-0.001);
      //    printf("%d (%f)  %d(%f)\n",lowBin,hPtVsYVsMultReco->GetZaxis()->GetBinLowEdge(lowBin),hiBin,hPtVsYVsMultReco->GetZaxis()->GetBinUpEdge(hiBin));
    
      hPtRecoM[j]=h3dReco->ProjectionX(Form("hPtRecoM%d",j),0,-1,lowBin,hiBin);
      hPtGenAccM[j]=h3dGenAcc->ProjectionX(Form("hPtGenAccM%d",j),0,-1,lowBin,hiBin);
      //    hPtGenLimAccM[j]=h3dGenLimAcc->ProjectionX(Form("hPtGenLimAccM%d",j),0,-1,lowBin,hiBin);
      hEffVsPtM[j]=(TH1D*)hPtRecoM[j]->Clone(Form("hEffM%d",j));
      hEffVsPtM[j]->Divide(hPtRecoM[j],hPtGenAccM[j],1,1,"B");
    }

    hEffVsPtM[0]->SetStats(0);
    TCanvas* cwp=new TCanvas(Form("cwp%s",dCase.Data()),Form("%s - Eff and Wei vs. Pt",dCase.Data()),1200,600);
    cwp->Divide(2,1);
    cwp->cd(1);
    gPad->SetLeftMargin(0.12);
    hEffVsPtM[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hEffVsPtM[0]->GetYaxis()->SetTitle("Efficiency");
    hEffVsPtM[0]->GetYaxis()->SetTitleOffset(1.4);
    hEffVsPtM[0]->SetLineColor(mulcol[0]);
    hEffVsPtM[0]->SetMinimum(0);
    hEffVsPtM[0]->SetMaximum(1.8);
    hEffVsPtM[0]->Draw();
    TLegend* legp=new TLegend(0.16,0.5,0.55,0.89);
    legp->SetFillStyle(0);
    legp->SetBorderSize(0);
    legp->AddEntry(hEffVsPtM[0],Form("%.0f<N_{tracklets}<%.0f",mulLims[0],mulLims[1]),"L")->SetTextColor(mulcol[0]);
    for(Int_t j=1; j<nMultBins; j++){
      hEffVsPtM[j]->SetLineColor(mulcol[j]);
      hEffVsPtM[j]->Draw("same");
      legp->AddEntry(hEffVsPtM[j],Form("%.0f<N_{tracklets}<%.0f",mulLims[j],mulLims[j+1]),"L")->SetTextColor(mulcol[j]);
    }
    legp->Draw();
    cwp->cd(2);
    gPad->SetLeftMargin(0.12);
    funcPtWeight->SetTitle("");
    funcPtWeight->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    funcPtWeight->GetYaxis()->SetTitle("Weight");
    funcPtWeight->Draw();
    cwp->SaveAs(Form("figures/EfficVsPtMultBins_%s_%s.eps",suffix.Data(),dCase.Data()));
  
    TH1D* hMultReco[nPtBins];
    TH1D* hMultGenAcc[nPtBins];
    //  TH1D* hMultGenLimAcc[nPtBins];
    for(Int_t j=0; j<nPtBins; j++){
      Int_t lowBin=h3dReco->GetXaxis()->FindBin(binLims[j]);
      Int_t hiBin=h3dReco->GetXaxis()->FindBin(binLims[j+1]-0.001);
      //printf("%d (%f)  %d(%f)\n",lowBin,h3dReco->GetXaxis()->GetBinLowEdge(lowBin),hiBin,h3dReco->GetXaxis()->GetBinUpEdge(hiBin));
    
      hMultReco[j]=h3dReco->ProjectionZ(Form("hMultReco%s%d",dCase.Data(),j),lowBin,hiBin);
      hMultGenAcc[j]=h3dGenAcc->ProjectionZ(Form("hMultGenAcc%s%d",dCase.Data(),j),lowBin,hiBin);
      //    hMultGenLimAcc[j]=h3dGenLimAcc->ProjectionZ(Form("hMultGenLimAcc%d",j),lowBin,hiBin);
      hEffVsMult[j]=(TH1D*)hMultReco[j]->Clone(Form("hEff%sVsMultPtBin%d",dCase.Data(),j));
      hEffVsMult[j]->Divide(hMultReco[j],hMultGenAcc[j],1,1,"B");
      hEffVsMult[j]->SetStats(0);
      hEffVsMult[j]->SetTitle(Form("%.1f<p_{T}<%.1f GeV/c",binLims[j],binLims[j+1]));
      hEffVsMult[j]->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
      hEffVsMult[j]->GetYaxis()->SetTitle("Efficiency");
    }

    // TCanvas* c2e=new TCanvas("c2e");
    // hEffVsMult[0]->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
    // hEffVsMult[0]->GetYaxis()->SetTitle("Efficiency");
    // hEffVsMult[0]->GetYaxis()->SetTitleOffset(1.4);
    // hEffVsMult[0]->SetLineColor(ptcol[0]);
    // hEffVsMult[0]->SetMinimum(0);
    // hEffVsMult[0]->SetMaximum(1.6);
    // hEffVsMult[0]->Draw();
    // TLegend* legm=new TLegend(0.16,0.5,0.4,0.89);
    // legm->SetFillStyle(0);
    // legm->SetBorderSize(0);
    // legm->AddEntry(hEffVsMult[0],Form("%.0f<p_{T}<%.0f GeV/c",binLims[0],binLims[1]),"L")->SetTextColor(ptcol[0]);
    // for(Int_t j=1; j<nPtBins; j++){
    //   hEffVsMult[j]->SetLineColor(ptcol[j]);
    //   hEffVsMult[j]->Draw("same");
    //   legm->AddEntry(hEffVsMult[j],Form("%.0f<p_{T}<%.0f GeV/c",binLims[j],binLims[j+1]),"L")->SetTextColor(ptcol[j]);
    // }
    // legm->Draw();



    hEffVsMult[0]->SetStats(0);
    TCanvas* cw=new TCanvas(Form("cw%s",dCase.Data()),Form("%s - Eff and wei vs. mult",dCase.Data()),1200,600);
    cw->Divide(2,1);
    cw->cd(1);
    gPad->SetLeftMargin(0.12);
    if(maxMult>0) hEffVsMult[0]->GetXaxis()->SetRangeUser(0.,maxMult);
    hEffVsMult[0]->GetYaxis()->SetTitleOffset(1.4);
    hEffVsMult[0]->SetLineColor(ptcol[0]);
    hEffVsMult[0]->SetMinimum(0);
    hEffVsMult[0]->SetMaximum(1.6);
    hEffVsMult[0]->Draw();
    TLegend* legm=new TLegend(0.16,0.5,0.4,0.89);
    legm->SetFillStyle(0);
    legm->SetBorderSize(0);
    legm->AddEntry(hEffVsMult[0],Form("%.0f<p_{T}<%.0f GeV/c",binLims[0],binLims[1]),"L")->SetTextColor(ptcol[0]);
    for(Int_t j=1; j<nPtBins; j++){
      hEffVsMult[j]->SetLineColor(ptcol[j]);
      hEffVsMult[j]->Draw("same");
      legm->AddEntry(hEffVsMult[j],Form("%.0f<p_{T}<%.0f GeV/c",binLims[j],binLims[j+1]),"L")->SetTextColor(ptcol[j]);
    }
    legm->Draw();
    cw->cd(2);
    if(hMultWeight){
      gPad->SetLeftMargin(0.12);
      if(maxMult>0) hMultWeight->GetXaxis()->SetRangeUser(0.,maxMult);
      hMultWeight->Draw();
      TLegend* legw=new TLegend(0.2,0.65,0.55,0.89);
      legw->SetFillStyle(0);
      legw->AddEntry(hMultWeight,hMultWeight->GetName(),"L");
      legw->Draw();
    }
    cw->SaveAs(Form("figures/Effic%sVsMultPtBins_%s.eps",dCase.Data(),suffix.Data()));
  }

  if(var2=="PtB" || var3=="PtB"){
    TCanvas* cwb=new TCanvas("cwb",Form("%s - ptB weight",dCase.Data()),800,800);
    funcPtBWeight->SetTitle("");
    funcPtBWeight->GetXaxis()->SetTitle("B-hadron p_{T} (GeV/c)");
    funcPtBWeight->GetYaxis()->SetTitle("Weight (FONLL/MC)");
    funcPtBWeight->Draw();
  }
  
  TH1D *hpteffNoWeight =new TH1D(Form("hEff%sVsPtNoWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hpteffMultWeight = new TH1D(Form("hEff%sVsPtMultWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hpteffMultPtWeight =new TH1D(Form("hEff%sVsPtMultAndPtWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hpteffMultPtBWeight =new TH1D(Form("hEff%sVsPtMultAndPtBWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hpteffPtBWeight =new TH1D(Form("hEff%sVsPtPtBWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hptaccNoWeight =new TH1D(Form("hAcc%sVsPtNoWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hptaccMultWeight = new TH1D(Form("hAcc%sVsPtMultWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hptaccMultPtWeight =new TH1D(Form("hAcc%sVsPtMultAndPtWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hptaccMultPtBWeight =new TH1D(Form("hAcc%sVsPtMultAndPtBWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hptaccPtBWeight =new TH1D(Form("hAcc%sVsPtPtBWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hptaxeNoWeight =new TH1D(Form("hAxe%sVsPtNoWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hptaxeMultWeight = new TH1D(Form("hAxe%sVsPtMultWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hptaxeMultPtWeight =new TH1D(Form("hAxe%sVsPtMultAndPtWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hptaxeMultPtBWeight =new TH1D(Form("hAxe%sVsPtMultAndPtBWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hptaxePtBWeight =new TH1D(Form("hAxe%sVsPtPtBWeight",dCase.Data()),"",nPtBins,binLims);

  TH1D* hMultRecoAllPtW=h3dReco->ProjectionZ(Form("hMultRecoW%s",dCase.Data()));
  TH1D* hMultGenAccAllPtW=h3dGenAcc->ProjectionZ(Form("hMultGenAccW%s",dCase.Data()));
  TH1D* hMultGenLimAccAllPtW=h3dGenLimAcc->ProjectionZ(Form("hMultGenLimAccW%s",dCase.Data()));
  hMultRecoAllPtW->Reset("ines");
  hMultGenAccAllPtW->Reset("ines");
  hMultGenLimAccAllPtW->Reset("ines");

  Double_t countNumer[nPtBins];
  Double_t countDenom[nPtBins];
  Double_t countDenomLimAcc[nPtBins];
  Double_t countNumerPtBWei[nPtBins];
  Double_t countDenomPtBWei[nPtBins];
  Double_t countDenomLimAccPtBWei[nPtBins];

  Double_t countNumerMultWei[nPtBins];
  Double_t countDenomMultWei[nPtBins];
  Double_t countDenomLimAccMultWei[nPtBins];
  Double_t countNumerMultPtWei[nPtBins];
  Double_t countDenomMultPtWei[nPtBins];
  Double_t countDenomLimAccMultPtWei[nPtBins];
  Double_t countNumerMultPtBWei[nPtBins];
  Double_t countDenomMultPtBWei[nPtBins];
  Double_t countDenomLimAccMultPtBWei[nPtBins];

  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){
    countNumer[iPtBin]=0.0;
    countDenom[iPtBin]=0.0;
    countDenomLimAcc[iPtBin]=0.0;
    countNumerPtBWei[iPtBin]=0.0;
    countDenomPtBWei[iPtBin]=0.0;
    countDenomLimAccPtBWei[iPtBin]=0.0;
    countNumerMultWei[iPtBin]=0.0;
    countDenomMultWei[iPtBin]=0.0;
    countDenomLimAccMultWei[iPtBin]=0.0;
    countNumerMultPtWei[iPtBin]=0.0;
    countDenomMultPtWei[iPtBin]=0.0;
    countDenomLimAccMultPtWei[iPtBin]=0.0;
    countNumerMultPtBWei[iPtBin]=0.0;
    countDenomMultPtBWei[iPtBin]=0.0;
    countDenomLimAccMultPtBWei[iPtBin]=0.0;
  }

  for(Int_t ibx=0; ibx<=h3dGenAcc->GetNbinsX()+1; ibx++){
    Double_t pt=h3dGenAcc->GetXaxis()->GetBinCenter(ibx);
    Double_t wpt=funcPtWeight->Eval(pt);
    Int_t jPtBin=TMath::BinarySearch(nPtBins+1,binLims,pt);
    if(jPtBin>=0 && jPtBin<nPtBins){
      for(Int_t iby=0; iby<=h3dGenAcc->GetNbinsY()+1; iby++){
	Double_t y=h3dGenAcc->GetYaxis()->GetBinCenter(iby);
	Double_t ptB=-999.;
	Double_t wptB=1.;
	if(var2=="PtB"){
	  ptB=h3dGenAcc->GetYaxis()->GetBinCenter(iby);
	  wptB=funcPtBWeight->Eval(ptB);
	}
	for(Int_t ibz=0; ibz<=h3dGenAcc->GetNbinsZ()+1; ibz++){
	  Double_t crec=h3dReco->GetBinContent(ibx,iby,ibz);
	  Double_t cgen=h3dGenAcc->GetBinContent(ibx,iby,ibz);
	  Double_t cgenla=h3dGenLimAcc->GetBinContent(ibx,iby,ibz);
	  Double_t wmult=1.;
	  if(var3=="PtB"){
	    ptB=h3dGenAcc->GetZaxis()->GetBinCenter(ibz);
	    wptB=funcPtBWeight->Eval(ptB);
	  }else if(var3=="Mult"){
	    Double_t mult=h3dGenAcc->GetZaxis()->GetBinCenter(ibz);
	    if(hMultWeight){
	      Int_t binw=hMultWeight->FindBin(mult);
	      if(binw>=1 && binw<hMultWeight->GetNbinsX()+1){
		wmult=hMultWeight->GetBinContent(binw);
		//		printf("mult %.0f   bin %d   wei %f\n",mult,binw,w);
	      }else{
		if(cgen>0){
		  printf("mult %.0f   bin %d   wei %f\n",mult,binw,wmult);
		  getchar();
		}
	      }
	    }
	    hMultRecoAllPtW->Fill(mult,crec*wmult);
	    hMultGenAccAllPtW->Fill(mult,cgen*wmult);
	    hMultGenLimAccAllPtW->Fill(mult,h3dGenLimAcc->GetBinContent(ibx,iby,ibz)*wmult);
	  }
	  countNumer[jPtBin]+=crec;
	  countDenom[jPtBin]+=cgen;
	  countDenomLimAcc[jPtBin]+=cgenla;
	  countNumerPtBWei[jPtBin]+=crec*wptB;
	  countDenomPtBWei[jPtBin]+=cgen*wptB;
	  countDenomLimAccPtBWei[jPtBin]+=cgenla*wptB;
	  countNumerMultWei[jPtBin]+=crec*wmult;
	  countDenomMultWei[jPtBin]+=cgen*wmult;
	  countDenomLimAccMultWei[jPtBin]+=cgenla*wmult;
	  countNumerMultPtWei[jPtBin]+=crec*wmult*wpt;
	  countDenomMultPtWei[jPtBin]+=cgen*wmult*wpt;
	  countDenomLimAccMultPtWei[jPtBin]+=cgenla*wmult*wpt;
	  countNumerMultPtBWei[jPtBin]+=crec*wmult*wptB;
	  countDenomMultPtBWei[jPtBin]+=cgen*wmult*wptB;
	  countDenomLimAccMultPtBWei[jPtBin]+=cgenla*wmult*wptB;
	}
      }
    }
  }

  TH1F* hErrEff1=new TH1F("hErrEff1","Reco/GenAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrEff1b=new TH1F("hErrEff1b","Reco/GenAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAcc1=new TH1F("hErrAcc1","GenAcc/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAcc1b=new TH1F("hErrAcc1b","GenAcc/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAxe1=new TH1F("hErrAxe1","Reco/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAxe1b=new TH1F("hErrAxe1b","Reco/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAxe1c=new TH1F("hErrAxe1c","Reco/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrEff2=new TH1F("hErrEff2","Reco/GenAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAcc2=new TH1F("hErrAcc2","GenAcc/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAxe2=new TH1F("hErrAxe2","Reco/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrEff3=new TH1F("hErrEff3","Reco/GenAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAcc3=new TH1F("hErrAcc3","GenAcc/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAxe3=new TH1F("hErrAxe3","Reco/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrEff4=new TH1F("hErrEff4","Reco/GenAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAcc4=new TH1F("hErrAcc4","GenAcc/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAxe4=new TH1F("hErrAxe4","Reco/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrEff5=new TH1F("hErrEff5","Reco/GenAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAcc5=new TH1F("hErrAcc5","GenAcc/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  TH1F* hErrAxe5=new TH1F("hErrAxe5","Reco/GenLimAcc ; p_{T} (GeV/c) ; Relative Uncertainty (%)",nPtBins,binLims);
  hErrEff1->SetStats(0);
  hErrEff1b->SetStats(0);
  hErrEff2->SetStats(0);
  hErrEff3->SetStats(0);
  hErrEff4->SetStats(0);
  hErrEff5->SetStats(0);
  hErrAcc1->SetStats(0);
  hErrAcc1b->SetStats(0);
  hErrAcc2->SetStats(0);
  hErrAcc3->SetStats(0);
  hErrAcc4->SetStats(0);
  hErrAcc5->SetStats(0);
  hErrAxe1->SetStats(0);
  hErrAxe1b->SetStats(0);
  hErrAxe1c->SetStats(0);
  hErrAxe2->SetStats(0);
  hErrAxe3->SetStats(0);
  hErrAxe4->SetStats(0);
  hErrAxe5->SetStats(0);

  Int_t iybin1=hPtVsYGenAccToy->GetYaxis()->FindBin(-0.499999);
  Int_t iybin2=hPtVsYGenAccToy->GetYaxis()->FindBin(0.499999);
  TH1D* hptga=(TH1D*)hPtVsYGenAccToy->ProjectionX("hptga");
  TH1D* hptgay05=(TH1D*)hPtVsYGenAccToy->ProjectionX("hptga05",iybin1,iybin2);
  TH1D* hptgla=(TH1D*)hPtVsYGenLimAccToy->ProjectionX("hptgla");
  TH1D* hptgaR=(TH1D*)hptga->Rebin(nPtBins,"hptgaR",binLims);
  TH1D* hptglaR=(TH1D*)hptgla->Rebin(nPtBins,"hptglaR",binLims);
  TH1D* hptgay05R=(TH1D*)hptgay05->Rebin(nPtBins,"hptgay05R",binLims);
  TH1D* hr1=(TH1D*)hptgay05R->Clone("hr1");
  hr1->Divide(hptgay05R,hptgaR,1,1,"B");
  TH1D* hr2=(TH1D*)hptgay05R->Clone("hr2");
  hr2->Divide(hptgay05R,hptglaR,1,1,"B");
  TH1D* hr3=(TH1D*)hr1->Clone("hr3");
  hr3->Multiply(hr1,hr2);
  TGraph* g1=new TGraph(0);
  TGraph* g2=new TGraph(0);
  TGraph* g3=new TGraph(0);
  for(Int_t ip=0; ip<hptgay05R->GetNbinsX(); ip++){
    double ptcent=hptgay05R->GetBinCenter(ip+1);
    Double_t yfid=0.8;
    if(ptcent<5) yfid=-0.2/15*ptcent*ptcent+1.9/15*ptcent+0.5;
    g1->SetPoint(ip,ptcent,0.5/yfid);
    double acc=hAccToy->GetBinContent(hAccToy->FindBin(ptcent));
    g2->SetPoint(ip,ptcent,acc/1.6);
    g3->SetPoint(ip,ptcent,0.5/yfid*acc/1.6);
  }

  hr1->SetStats(0);
  hr2->SetStats(0);
  hr3->SetStats(0);
  hr1->SetMinimum(0.3);
  hr1->SetMaximum(1.03);
  hr2->SetMinimum(0.3);
  hr2->SetMaximum(1.03);
  hr3->SetMinimum(0.3);
  hr3->SetMaximum(1.03);
  hr1->SetLineWidth(2);
  hr2->SetLineWidth(2);
  hr3->SetLineWidth(2);
  hr1->GetYaxis()->SetTitle("N(GenAcc && |y|<0.5) / N(GenAcc)");
  hr2->GetYaxis()->SetTitle("N(GenAcc && |y|<0.5) / N(|y|<0.5)");
  hr3->GetYaxis()->SetTitle("N(GenAcc && |y|<0.5)^{2} /  (N(GenAcc) * N(|y|<0.5))");
  hr1->GetYaxis()->SetTitleOffset(1.3);
  hr2->GetYaxis()->SetTitleOffset(1.3);
  hr3->GetYaxis()->SetTitleOffset(1.3);

  if(dCase=="Prompt"){
    TCanvas* crho=new TCanvas("crho","coveriance",1600,600);
    crho->Divide(3,1);
    crho->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    hr1->Draw();
    g1->SetLineColor(2);
    g1->SetLineWidth(2);
    g1->Draw("csame");
    TLegend* l1=new TLegend(0.5,0.7,0.9,0.89);
    l1->SetTextFont(43);
    l1->SetTextSize(22);
    l1->AddEntry(hr1,"Toy MC, Fonll y","L");
    l1->AddEntry(g1,"0.5/y_{fid}","L")->SetTextColor(2);
    l1->Draw();
    crho->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    hr2->Draw();
    g2->SetLineColor(2);
    g2->SetLineWidth(2);
    g2->Draw("csame");
    TLegend* l2=new TLegend(0.3,0.2,0.9,0.39);
    l2->SetTextFont(43);
    l2->SetTextSize(22);
    l2->AddEntry(hr2,"Toy MC, Fonll y","L");
    l2->AddEntry(g2,"(GenAcc/GenLimAcc)/1.6","L")->SetTextColor(2);
    l2->Draw();
    crho->cd(3);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    hr3->Draw();
    g3->SetLineColor(2);
    g3->SetLineWidth(2);
    g3->Draw("csame");
    TLegend* l3=new TLegend(0.2,0.7,0.9,0.89);
    l3->SetMargin(0.15);
    l3->SetTextFont(43);
    l3->SetTextSize(22);
    l3->AddEntry(hr3,"Toy MC, Fonll y","L");
    l3->AddEntry(g3,"0.5/y_{fid}*(GenAcc/GenLimAcc)/1.6","L")->SetTextColor(2);
    l3->Draw();
  }
  
  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){
    printf("---- Pt range %.1f - %.1f ----\n",binLims[iPtBin],binLims[iPtBin+1]);
    Double_t eff1=countNumer[iPtBin]/countDenom[iPtBin];
    Double_t erreff1=TMath::Sqrt(eff1*(1-eff1)/countDenom[iPtBin]);
    // sqrt(n/d*(1-n/d)/d) = sqrt(n/d * (d-n)/d / d) = sqrt( n*(d-n)/d/d/d) = 1/d * sqrt(n*(d-n)/d) = eff * sqrt((d-n)/dn)
    Double_t rho=TMath::Sqrt(eff1);
    Double_t erreff1b=eff1*TMath::Sqrt(1./countNumer[iPtBin]+1./countDenom[iPtBin]-2.*rho*1./TMath::Sqrt(countNumer[iPtBin]*countDenom[iPtBin]));
    // eff * sqrt(1/n + 1/d -2 sqrt(eff)/sqrt(nd)) = eff * sqrt( (n+d)/nd -2sqrt(1/d^2)) = eff*sqrt( (n+d)/nd -2/d) = eff*sqrt( (n-d)/2d 
    Double_t acc1=countDenom[iPtBin]/countDenomLimAcc[iPtBin];
    Double_t yfid=0.8;
    Double_t ptcent=0.5*(binLims[iPtBin]+binLims[iPtBin+1]);
    if(ptcent<5) yfid=-0.2/15*ptcent*ptcent+1.9/15*ptcent+0.5;
    rho=TMath::Sqrt(0.5/yfid*acc1/1.6);
    Double_t rhoToy=TMath::Sqrt(hr3->GetBinContent(hr3->FindBin(ptcent)));
    Double_t erracc1=acc1*TMath::Sqrt(1./countDenom[iPtBin]+1./countDenomLimAcc[iPtBin]-2*rho*1/TMath::Sqrt(countDenom[iPtBin]*countDenomLimAcc[iPtBin]));
    Double_t erracc1b=acc1*TMath::Sqrt(1./countDenom[iPtBin]+1./countDenomLimAcc[iPtBin]-2*rhoToy*1/TMath::Sqrt(countDenom[iPtBin]*countDenomLimAcc[iPtBin]));
    Double_t axe1=countNumer[iPtBin]/countDenomLimAcc[iPtBin];
    Double_t erraxe1=TMath::Sqrt(axe1*(1-axe1)/countDenomLimAcc[iPtBin]);
    rho=TMath::Sqrt(0.5/yfid*acc1/1.6*eff1);
    rhoToy=TMath::Sqrt(hr3->GetBinContent(hr3->FindBin(ptcent))*eff1);
    Double_t erraxe1b=axe1*TMath::Sqrt(1./countNumer[iPtBin]+1./countDenomLimAcc[iPtBin]-2*rho*1/TMath::Sqrt(countNumer[iPtBin]*countDenomLimAcc[iPtBin]));
    Double_t erraxe1c=axe1*TMath::Sqrt(1./countNumer[iPtBin]+1./countDenomLimAcc[iPtBin]-2*rhoToy*1/TMath::Sqrt(countNumer[iPtBin]*countDenomLimAcc[iPtBin]));
    
    printf("Eff from Projection = %f/%f = %f+-%f\n",hPtRecoR->GetBinContent(iPtBin+1),hPtGenAccR->GetBinContent(iPtBin+1),hEffVsPtR->GetBinContent(iPtBin+1),hEffVsPtR->GetBinError(iPtBin+1));
    printf("Eff No weights      = %f/%f = %f+-%f\n",countNumer[iPtBin],countDenom[iPtBin],eff1,erreff1);
    hpteffNoWeight->SetBinContent(iPtBin+1,eff1);
    hpteffNoWeight->SetBinError(iPtBin+1,erreff1);
    hptaccNoWeight->SetBinContent(iPtBin+1,acc1);
    hptaccNoWeight->SetBinError(iPtBin+1,erracc1);
    hptaxeNoWeight->SetBinContent(iPtBin+1,axe1);
    hptaxeNoWeight->SetBinError(iPtBin+1,erraxe1);
    hErrEff1->SetBinContent(iPtBin+1,100*erreff1/eff1);
    hErrEff1b->SetBinContent(iPtBin+1,100*erreff1b/eff1);
    hErrAcc1->SetBinContent(iPtBin+1,100*erracc1/acc1);
    hErrAcc1b->SetBinContent(iPtBin+1,100*erracc1b/acc1);
    hErrAxe1->SetBinContent(iPtBin+1,100*erraxe1/axe1);
    hErrAxe1b->SetBinContent(iPtBin+1,100*erraxe1b/axe1);
    hErrAxe1c->SetBinContent(iPtBin+1,100*erraxe1c/axe1);
    if(dCase=="Feeddw"){
      Double_t eff4=countNumerPtBWei[iPtBin]/countDenomPtBWei[iPtBin];
      Double_t erreff4=TMath::Sqrt(eff4*(1-eff4)/countDenom[iPtBin]);// countDenom is NOT a typo, it has to be like this to get proper statistical errors from the no-weight case
      Double_t acc4=countDenomPtBWei[iPtBin]/countDenomLimAccPtBWei[iPtBin];
      Double_t erracc4=erracc1b*acc4/acc1;
      Double_t axe4=countNumerPtBWei[iPtBin]/countDenomLimAccPtBWei[iPtBin];
      Double_t erraxe4=erraxe1c*axe4/axe1;
      printf("Eff With pt(B) weights = %f/%f = %f+-%f\n",countNumerPtBWei[iPtBin],countDenomPtBWei[iPtBin],eff4,erreff4);
      hpteffPtBWeight->SetBinContent(iPtBin+1,eff4);
      hpteffPtBWeight->SetBinError(iPtBin+1,erreff4);
      hptaccPtBWeight->SetBinContent(iPtBin+1,acc4);
      hptaccPtBWeight->SetBinError(iPtBin+1,erracc4);
      hptaxePtBWeight->SetBinContent(iPtBin+1,axe4);
      hptaxePtBWeight->SetBinError(iPtBin+1,erraxe4);
      hErrEff4->SetBinContent(iPtBin+1,100*erreff4/eff4);
      hErrAcc4->SetBinContent(iPtBin+1,100*erracc4/acc4);
      hErrAxe4->SetBinContent(iPtBin+1,100*erraxe4/axe4);
    }
    if(var3=="Mult"){
      Double_t eff2=countNumerMultWei[iPtBin]/countDenomMultWei[iPtBin];
      Double_t erreff2=TMath::Sqrt(eff2*(1-eff2)/countDenom[iPtBin]);// countDenom is NOT a typo, it has to be like this to get proper statistical errors from the no-weight case
      Double_t acc2=countDenomMultWei[iPtBin]/countDenomLimAccMultWei[iPtBin];
      Double_t erracc2=erracc1b*acc2/acc1;
      Double_t axe2=countNumerMultWei[iPtBin]/countDenomLimAccMultWei[iPtBin];
      Double_t erraxe2=erraxe1c*axe2/axe1;
      printf("Eff With mult weights  = %f/%f = %f+-%f\n",countNumerMultWei[iPtBin],countDenomMultWei[iPtBin],eff2,erreff2);
      hpteffMultWeight->SetBinContent(iPtBin+1,eff2);
      hpteffMultWeight->SetBinError(iPtBin+1,erreff2);
      hptaccMultWeight->SetBinContent(iPtBin+1,acc2);
      hptaccMultWeight->SetBinError(iPtBin+1,erracc2);
      hptaxeMultWeight->SetBinContent(iPtBin+1,axe2);
      hptaxeMultWeight->SetBinError(iPtBin+1,erraxe2);
      hErrEff2->SetBinContent(iPtBin+1,100*erreff2/eff2);
      hErrAcc2->SetBinContent(iPtBin+1,100*erracc2/acc2);
      hErrAxe2->SetBinContent(iPtBin+1,100*erraxe2/axe2);
      Double_t eff3=countNumerMultPtWei[iPtBin]/countDenomMultPtWei[iPtBin];
      Double_t erreff3=TMath::Sqrt(eff3*(1-eff3)/countDenom[iPtBin]);// countDenom is NOT a typo, it has to be like this to get proper statistical errors from the no-weight case
      Double_t acc3=countDenomMultPtWei[iPtBin]/countDenomLimAccMultPtWei[iPtBin];
      Double_t erracc3=erracc1b*acc3/acc1;
      Double_t axe3=countNumerMultPtWei[iPtBin]/countDenomLimAccMultPtWei[iPtBin];
      Double_t erraxe3=erraxe1c*axe3/axe1;
      printf("Eff With mult+pt weights   = %f/%f = %f+-%f\n",countNumerMultPtWei[iPtBin],countDenomMultPtWei[iPtBin],eff3,erreff3);
      hpteffMultPtWeight->SetBinContent(iPtBin+1,eff3);
      hpteffMultPtWeight->SetBinError(iPtBin+1,erreff3);
      hptaccMultPtWeight->SetBinContent(iPtBin+1,acc3);
      hptaccMultPtWeight->SetBinError(iPtBin+1,erracc3);
      hptaxeMultPtWeight->SetBinContent(iPtBin+1,axe3);
      hptaxeMultPtWeight->SetBinError(iPtBin+1,erraxe3);
      hErrEff3->SetBinContent(iPtBin+1,100*erreff3/eff3);
      hErrAcc3->SetBinContent(iPtBin+1,100*erracc3/acc3);
      hErrAxe3->SetBinContent(iPtBin+1,100*erraxe3/axe3);
      if(dCase=="Feeddw"){
	Double_t eff5=countNumerMultPtBWei[iPtBin]/countDenomMultPtBWei[iPtBin];
	Double_t erreff5=TMath::Sqrt(eff5*(1-eff5)/countDenom[iPtBin]);// countDenom is NOT a typo, it has to be like this to get proper statistical errors from the no-weight case
	Double_t acc5=countDenomMultPtBWei[iPtBin]/countDenomLimAccMultPtBWei[iPtBin];
	Double_t erracc5=erracc1b*acc5/acc1;
	Double_t axe5=countNumerMultPtBWei[iPtBin]/countDenomLimAccMultPtBWei[iPtBin];
	Double_t erraxe5=erraxe1c*axe5/axe1;
 	printf("Eff With mult+pt(B) weights   = %f/%f = %f+-%f\n",countNumerMultPtBWei[iPtBin],countDenomMultPtBWei[iPtBin],eff5,erreff5);
	hpteffMultPtBWeight->SetBinContent(iPtBin+1,eff5);
	hpteffMultPtBWeight->SetBinError(iPtBin+1,erreff5);
	hptaccMultPtBWeight->SetBinContent(iPtBin+1,acc5);
	hptaccMultPtBWeight->SetBinError(iPtBin+1,erracc5);
	hptaxeMultPtBWeight->SetBinContent(iPtBin+1,axe5);
	hptaxeMultPtBWeight->SetBinError(iPtBin+1,erraxe5);
	hErrEff5->SetBinContent(iPtBin+1,100*erreff5/eff5);
	hErrAcc5->SetBinContent(iPtBin+1,100*erracc5/acc5);
	hErrAxe5->SetBinContent(iPtBin+1,100*erraxe5/axe5);
      }
    }
  }

  TCanvas* ceer=new TCanvas(Form("ceer%s",dCase.Data()),Form("%s - Unc on Eff",dCase.Data()),1600,600);
  ceer->Divide(3,1);
  ceer->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetTickx();
  gPad->SetTicky();
  hErrEff1b->SetLineColor(kCyan);
  hErrEff1b->SetLineWidth(4);
  hErrEff1b->GetYaxis()->SetTitleOffset(1.8);
  hErrEff1b->Draw();
  hErrEff1->SetLineColor(4);
  hErrEff1->Draw("same");
  hErrEff2->SetLineColor(2);
  hErrEff2->SetLineStyle(2);
  hErrEff2->SetLineWidth(2);
  hErrEff2->Draw("same");
  hErrEff3->SetLineColor(kGreen+1);
  hErrEff3->SetLineStyle(3);
  hErrEff3->SetLineWidth(3);
  hErrEff3->Draw("same");
  TLegend* legErr=new TLegend(0.4,0.15,0.92,0.4);
  legErr->AddEntry(hErrEff1b,"Covariance, no weights","L");
  legErr->AddEntry(hErrEff1,"Binomial, no weights","L");
  legErr->AddEntry(hErrEff2,"Binomial, mult weight","L");
  legErr->AddEntry(hErrEff3,"Binomial, mult+pt weights","L");
  if(dCase=="Feeddw"){
    hErrEff4->SetLineColor(kOrange+1);
    hErrEff4->SetLineStyle(4);
    hErrEff4->SetLineWidth(2);
    hErrEff4->Draw("same");
    hErrEff5->SetLineColor(6);
    hErrEff5->SetLineStyle(5);
    hErrEff5->SetLineWidth(2);
    hErrEff5->Draw("same");
    legErr->AddEntry(hErrEff4,"Binomial, pt(B) weights","L");
    legErr->AddEntry(hErrEff5,"Binomial, mult+pt(B) weights","L");
  }
  legErr->Draw();
  ceer->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetTickx();
  gPad->SetTicky();
  hErrAcc1->SetLineColor(2);
  hErrAcc1->SetLineWidth(2);
  hErrAcc1->GetYaxis()->SetTitleOffset(1.8);
  hErrAcc1->SetMinimum(0.9*TMath::Min(hErrAcc1->GetMinimum(),hErrAcc1b->GetMinimum()));
  hErrAcc1->Draw();
  hErrAcc1b->SetLineWidth(2);
  hErrAcc1b->SetLineColor(1);
  hErrAcc1b->Draw("same");
  TLegend* legErrA=new TLegend(0.45,0.15,0.92,0.25);
  legErrA->SetMargin(0.15);
  legErrA->AddEntry(hErrAcc1,"Covariance, simple formula","L");
  legErrA->AddEntry(hErrAcc1b,"Covariance, ToyMC","L");
  legErrA->Draw();
  // hErrAcc2->SetLineColor(2);
  // hErrAcc2->SetLineStyle(2);
  // hErrAcc2->SetLineWidth(2);
  // hErrAcc2->Draw("same");
  // hErrAcc3->SetLineColor(kGreen+1);
  // hErrAcc3->SetLineStyle(3);
  // hErrAcc3->SetLineWidth(3);
  // hErrAcc3->Draw("same");
  // if(dCase=="Feeddw"){
  //   hErrAcc4->SetLineColor(4);
  //   hErrAcc4->SetLineStyle(4);
  //   hErrAcc4->SetLineWidth(2);
  //   hErrAcc4->Draw("same");
  //   hErrAcc5->SetLineColor(6);
  //   hErrAcc5->SetLineStyle(5);
  //   hErrAcc5->SetLineWidth(2);
  //   hErrAcc5->Draw("same");
  // }
  ceer->cd(3);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetTickx();
  gPad->SetTicky();
  hErrAxe1b->SetLineColor(2);
  hErrAxe1b->SetLineWidth(2);
  hErrAxe1b->GetYaxis()->SetTitleOffset(1.8);
  hErrAxe1b->Draw();
  hErrAxe1c->SetLineColor(1);
  hErrAxe1c->SetLineWidth(2);
  hErrAxe1c->Draw("same");
  hErrAxe1->SetLineColor(4);
  hErrAxe1->Draw("same");
  // hErrAxe2->SetLineColor(2);
  // hErrAxe2->SetLineStyle(2);
  // hErrAxe2->SetLineWidth(2);
  // hErrAxe2->Draw("same");
  // hErrAxe3->SetLineColor(kGreen+1);
  // hErrAxe3->SetLineStyle(3);
  // hErrAxe3->SetLineWidth(3);
  // hErrAxe3->Draw("same");
  // if(dCase=="Feeddw"){
  //   hErrAxe4->SetLineColor(kOrange+1);
  //   hErrAxe4->SetLineStyle(4);
  //   hErrAxe4->SetLineWidth(2);
  //   hErrAxe4->Draw("same");
  //   hErrAxe5->SetLineColor(6);
  //   hErrAxe5->SetLineStyle(5);
  //   hErrAxe5->SetLineWidth(2);
  //   hErrAxe5->Draw("same");
  // }
  
  hEffVsPtR->SetMarkerStyle(0);
  hEffVsPtR->SetMarkerColor(0);
  hEffVsPtR->SetMarkerSize(1.2);
  hEffVsPtR->SetLineColor(kGray);
  hEffVsPtR->SetLineWidth(4);
  hEffVsPtR->SetStats(0);
  hpteffNoWeight->SetMarkerStyle(20);
  hptaccNoWeight->SetMarkerStyle(20);
  
  if(var3=="Mult"){
    TH1D* hEffVsMultAllPtW=(TH1D*)hMultRecoAllPtW->Clone(Form("hEff%s",dCase.Data()));
    TH1D* hAccVsMultAllPtW=(TH1D*)hMultGenAccAllPtW->Clone(Form("hAcc%s",dCase.Data()));
    hEffVsMultAllPtW->Divide(hMultRecoAllPtW,hMultGenAccAllPtW,1,1,"B");
    hAccVsMultAllPtW->Divide(hMultGenAccAllPtW,hMultGenLimAccAllPtW,1,1,"B");
    hEffVsMultAllPtW->SetStats(0);
    hAccVsMultAllPtW->SetStats(0);

    TCanvas* c2w=new TCanvas(Form("c2w%s",dCase.Data()),Form("%s - EffVsMultW",dCase.Data()),1200,600);
    c2w->Divide(2,1);
    c2w->cd(1);
    gPad->SetLogy();
    hMultGenLimAccAllPtW->SetLineColor(1);
    if(maxMult>0) hMultGenLimAccAllPtW->GetXaxis()->SetRangeUser(0.,maxMult);
    hMultGenLimAccAllPtW->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
    hMultGenLimAccAllPtW->GetYaxis()->SetTitle("Entries");
    hMultGenLimAccAllPtW->Draw();
    hMultGenAccAllPtW->SetLineColor(2);
    hMultGenAccAllPtW->Draw("same");
    hMultRecoAllPtW->SetLineColor(4);
    hMultRecoAllPtW->Draw("same");
    c2w->cd(2);
    hEffVsMultAllPtW->SetLineColor(4);
    hEffVsMultAllPtW->SetMinimum(0);
    hEffVsMultAllPtW->SetMaximum(1.6);
    if(maxMult>0) hEffVsMultAllPtW->GetXaxis()->SetRangeUser(0.,maxMult);
    hEffVsMultAllPtW->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
    hEffVsMultAllPtW->GetYaxis()->SetTitle("Ratio");
    hEffVsMultAllPtW->Draw();
    hAccVsMultAllPtW->SetLineColor(2);
    hAccVsMultAllPtW->Draw("same");
    tacc2->Draw();
    te2->Draw();
    

    TH1F* hRatioMultWei=(TH1F*)hpteffMultWeight->Clone(Form("hRatioMultWei%s",dCase.Data()));
    hRatioMultWei->Divide(hpteffMultWeight,hpteffNoWeight);
    hRatioMultWei->GetYaxis()->SetTitle("With Mult weight / Without weight");
    hRatioMultWei->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hRatioMultWei->GetYaxis()->SetTitleOffset(1.4);
    hRatioMultWei->SetStats(0);    
    TH1F* hRatioPtWei=(TH1F*)hpteffMultPtWeight->Clone(Form("hRatioPtWei%s",dCase.Data()));
    hRatioPtWei->Divide(hpteffMultPtWeight,hpteffMultWeight);
    hRatioPtWei->GetYaxis()->SetTitle("With Mult+p_{T} weight / With only Mult weight");
    hRatioPtWei->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hRatioPtWei->GetYaxis()->SetTitleOffset(1.4);
    hRatioPtWei->SetStats(0);
    TH1F* hRatioPtBWei=(TH1F*)hpteffMultPtBWeight->Clone(Form("hRatioPtBWei%s",dCase.Data()));
    hRatioPtBWei->Divide(hpteffMultPtBWeight,hpteffMultWeight);
    hRatioPtBWei->GetYaxis()->SetTitle("With Mult+p_{T}(B) weight / With only Mult weight");
    hRatioPtBWei->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hRatioPtBWei->GetYaxis()->SetTitleOffset(1.4);
    hRatioPtBWei->SetStats(0);

    TCanvas* ceff=new TCanvas(Form("ceff%s",dCase.Data()),Form("%s - Eff Mult Wei",dCase.Data()),1200,600);
    ceff->Divide(2,1);
    ceff->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetTickx();
    gPad->SetTicky();
    hEffVsPtR->SetMinimum(0);
    hEffVsPtR->SetMaximum(1);
    hEffVsPtR->Draw();
    hEffVsPtR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hEffVsPtR->GetYaxis()->SetTitle("Efficiency");
    hEffVsPtR->GetYaxis()->SetTitleOffset(1.4);
    hpteffNoWeight->Draw("same");
    hpteffMultWeight->SetLineColor(kGreen+1);
    hpteffMultWeight->SetMarkerColor(kGreen+1);
    hpteffMultWeight->SetMarkerStyle(25);
    hpteffMultWeight->Draw("same");
    TLegend* leg=new TLegend(0.55,0.15,0.89,0.45);
    leg->SetFillStyle(0);
    leg->AddEntry(hEffVsPtR,"TH3F::Project","L");
    leg->AddEntry(hpteffNoWeight,"Multiplcity slices - No Weight","PL");
    leg->AddEntry(hpteffMultWeight,"Multiplcity slices - MultWei","PL");
    leg->Draw();
    ceff->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetTickx();
    gPad->SetTicky();
    hRatioMultWei->SetMarkerStyle(25);
    hRatioMultWei->SetLineWidth(2);
    hRatioMultWei->Draw();
    hRatioMultWei->SetMinimum(0.95);
    hRatioMultWei->SetMaximum(1.05);
    ceff->SaveAs(Form("figures/Effic%sWithMultWeights%s.eps",dCase.Data(),suffix.Data()));

    if(dCase=="Prompt"){
      TCanvas* ceff2=new TCanvas(Form("ceff2%s",dCase.Data()),Form("%s - Eff Mult+Pt wei",dCase.Data()),1200,600);
      ceff2->Divide(2,1);
      ceff2->cd(1);
      gPad->SetLeftMargin(0.12);
      gPad->SetTickx();
      gPad->SetTicky();
      hEffVsPtR->SetMinimum(0);
      hEffVsPtR->SetMaximum(1);
      hEffVsPtR->Draw();
      hEffVsPtR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hEffVsPtR->GetYaxis()->SetTitle("Efficiency");
      hEffVsPtR->GetYaxis()->SetTitleOffset(1.4);
      hpteffNoWeight->Draw("same");
      hpteffMultWeight->Draw("same");
      hpteffMultPtWeight->SetLineColor(2);
      hpteffMultPtWeight->SetMarkerColor(2);
      hpteffMultPtWeight->SetMarkerStyle(33);
      hpteffMultPtWeight->Draw("same");
      TLegend* leg2=new TLegend(0.55,0.15,0.89,0.45);
      leg2->SetFillStyle(0);
      leg2->AddEntry(hEffVsPtR,"TH3F::Project","L");
      leg2->AddEntry(hpteffNoWeight,"Multiplcity slices - No Weight","PL");
      leg2->AddEntry(hpteffMultWeight,"Multiplcity slices - Mult Weight","PL");
      leg2->AddEntry(hpteffMultPtWeight,"Multiplicity slices - Pt+Mult Weight","PL");
      leg2->Draw();
      ceff2->cd(2);
      gPad->SetLeftMargin(0.12);
      gPad->SetTickx();
      gPad->SetTicky();
      hRatioPtWei->SetMarkerStyle(33);
      hRatioPtWei->SetLineWidth(2);
      hRatioPtWei->Draw();
      hRatioPtWei->SetMinimum(0.95);
      hRatioPtWei->SetMaximum(1.05);
      ceff2->SaveAs(Form("figures/Effic%sWithMultAndPtWeights%s.eps",dCase.Data(),suffix.Data()));
    }
    if(var2=="PtB"){
      TCanvas* ceff2B=new TCanvas(Form("ceff2B%s",dCase.Data()),Form("%s - Eff Mult+Pt(B) wei",dCase.Data()),1200,600);
      ceff2B->Divide(2,1);
      ceff2B->cd(1);
      gPad->SetLeftMargin(0.12);
      gPad->SetTickx();
      gPad->SetTicky();
      hEffVsPtR->SetMinimum(0);
      hEffVsPtR->SetMaximum(1);
      hEffVsPtR->Draw();
      hEffVsPtR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hEffVsPtR->GetYaxis()->SetTitle("Efficiency");
      hEffVsPtR->GetYaxis()->SetTitleOffset(1.4);
      hpteffNoWeight->Draw("same");
      hpteffMultWeight->Draw("same");
      hpteffMultPtBWeight->SetLineColor(2);
      hpteffMultPtBWeight->SetMarkerColor(2);
      hpteffMultPtBWeight->SetMarkerStyle(33);
      hpteffMultPtBWeight->Draw("same");
      TLegend* leg2=new TLegend(0.55,0.15,0.89,0.45);
      leg2->SetFillStyle(0);
      leg2->AddEntry(hEffVsPtR,"TH3F::Project","L");
      leg2->AddEntry(hpteffNoWeight,"Multiplcity slices - No Weight","PL");
      leg2->AddEntry(hpteffMultWeight,"Multiplcity slices - Mult Weight","PL");
      leg2->AddEntry(hpteffMultPtBWeight,"Multiplicity slices - Pt(B)+Mult Weight","PL");
      leg2->Draw();
      ceff2B->cd(2);
      gPad->SetLeftMargin(0.12);
      gPad->SetTickx();
      gPad->SetTicky();
      hRatioPtBWei->SetMarkerStyle(29);
      hRatioPtBWei->SetLineWidth(2);
      hRatioPtBWei->Draw();
      hRatioPtBWei->SetMinimum(0.95);
      hRatioPtBWei->SetMaximum(1.05);
      ceff2B->SaveAs(Form("figures/Effic%sWithMultAndPtBWeights%s.eps",dCase.Data(),suffix.Data()));
    }
  }else{
    hpteffPtBWeight->SetMarkerStyle(21);
    hpteffPtBWeight->SetMarkerColor(kRed+1);
    hpteffPtBWeight->SetLineColor(kRed+1);

    TCanvas* ceff3=new TCanvas(Form("ceff3%s",dCase.Data()),Form("%s - Eff PtB wei",dCase.Data()),1200,600);
    ceff3->Divide(2,1);
    ceff3->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetTickx();
    gPad->SetTicky();
    hEffVsPtR->SetMinimum(0);
    hEffVsPtR->SetMaximum(1);
    hEffVsPtR->Draw();
    hEffVsPtR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hEffVsPtR->GetYaxis()->SetTitle("Efficiency");
    hEffVsPtR->GetYaxis()->SetTitleOffset(1.4);
    hpteffNoWeight->Draw("same");
    hpteffPtBWeight->Draw("same");
    TLegend* leg3=new TLegend(0.55,0.15,0.89,0.45);
    leg3->SetFillStyle(0);
    leg3->AddEntry(hEffVsPtR,"TH3F::Project","L");
    leg3->AddEntry(hpteffNoWeight,"Counting - No Weight","PL");
    leg3->AddEntry(hpteffPtBWeight,"Counting - p_{T}(B) Weight","PL");
    leg3->Draw();
    ceff3->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetTickx();
    gPad->SetTicky();
    TH1F* hRatioPtB=(TH1F*)hpteffPtBWeight->Clone(Form("hRatioPtB%s",dCase.Data()));
    hRatioPtB->Divide(hpteffPtBWeight,hpteffNoWeight);
    hRatioPtB->GetYaxis()->SetTitle("With/without p_{T}(B) weight");
    hRatioPtB->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hRatioPtB->GetYaxis()->SetTitleOffset(1.4);
    hRatioPtB->SetStats(0);
    hRatioPtB->Draw();
  }
  
  TFile* out=new TFile(Form("outputEff%s.root",suffix.Data()),"update");
  if(var3=="Mult") for(Int_t j=0; j<nPtBins; j++) if(hEffVsMult[j]) hEffVsMult[j]->Write();
  hEffVsPtR->Write();
  hpteffNoWeight->Write();
  hpteffMultWeight->Write();
  hpteffMultPtWeight->Write();
  hptaccNoWeight->Write();
  hptaccMultWeight->Write();
  hptaccMultPtWeight->Write();
  hptaxeNoWeight->Write();
  hptaxeMultWeight->Write();
  hptaxeMultPtWeight->Write();
  if(dCase=="Feeddw"){
    hpteffPtBWeight->Write();
    hpteffMultPtBWeight->Write();
    hptaccPtBWeight->Write();
    hptaccMultPtBWeight->Write();
    hptaxePtBWeight->Write();
    hptaxeMultPtBWeight->Write();
  }
  hEvSelEffVsPt->Write();
  out->Close();

}


Bool_t ReadConfig(TString configName){
  FILE* confFil=fopen(configName.Data(),"r");
  char line[100];
  char name[200];
  int n;
  float x;
  bool readok;
  while(!feof(confFil)){
    readok=fscanf(confFil,"%s:",line);
    if(strstr(line,"MCFile")){
      readok=fscanf(confFil,"%s",name);
      fileNameMC=name;
    }
    else if(strstr(line,"SuffixMC")){
      readok=fscanf(confFil,"%s",name);
      suffix=name;
    }
    else if(strstr(line,"AcceptanceFile")){
      readok=fscanf(confFil,"%s",name);
      fileNameToy=name;
    }
    else if(strstr(line,"MultWeiFile")){
      readok=fscanf(confFil,"%s",name);
      fileWeightName=name;
    }
    else if(strstr(line,"HistoWeiName")){
      readok=fscanf(confFil,"%s",name);
      histoWeightName=name;
    }
    else if(strstr(line,"PtDWei")){
      readok=fscanf(confFil,"%s",name);
      ptDWeight=name;
    }
    else if(strstr(line,"PtBWei")){
      readok=fscanf(confFil,"%s",name);
      ptBWeight=name;
    }
    else if(strstr(line,"UseAccFromToyMC")){
      readok=fscanf(confFil,"%d",&n);
      if(n>0) useToyMC=kTRUE;
      else useToyMC=kFALSE;
    }
    else if(strstr(line,"NumOfPtBins")){
      readok=fscanf(confFil,"%d",&n);
      nPtBins=n;
    }
    else if(strstr(line,"BinLimits")){
      readok=fscanf(confFil," [ ");
      for(int j=0; j<nPtBins; j++){
	readok=fscanf(confFil,"%f,",&x);
	binLims[j]=x;
	if(j>0 && binLims[j]<=binLims[j-1]){
	  printf("ERROR in array of pt bin limits\n");
	  return kFALSE;
	}
      }
      readok=fscanf(confFil,"%f",&x);
      binLims[nPtBins]=x;
      readok=fscanf(confFil," ] ");
    }
  }
  return kTRUE;
}

