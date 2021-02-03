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

enum EPtWei{kFONLL5overLHC13d3,kFONLL7overLHC10f7a,kFONLL7overLHC10f6a,kFLAToverLHC10f7a,kNoWei};
enum EPtBWei{kFONLL5overLHC19c3,kFONLL5overLHC20g2,kNoPtBWei};

TString configFileName="configfile4lowptanalysis.txt";
TString fileNameMC="";
TString suffix="";
TString fileNameToy="";

const Int_t maxPtBins=30;
Int_t nPtBins=8;
Double_t binLims[maxPtBins+1]={0.,1.,2.,3.,4.,5.,6.,8.,12.};
Int_t ptcol[maxPtBins]={1,kRed+1,kRed,kGreen+2,kCyan,4,kOrange+2,kMagenta,kMagenta+2,kBlue+1,kGray,kGray+2,kGreen,kYellow+7};

Int_t ptWeight=kNoWei;
Int_t ptBWeight=kFONLL5overLHC20g2;
Bool_t useMultWeight=kTRUE;
Double_t maxMult=-1;//200;

TH1F* hAccToyFine=0x0;
TH1F* hAccToy=0x0;
TF1* funcPtWeight=0x0;
TF1* funcPtBWeight=0x0;

TH1F** hWeight=new TH1F*[3];
Int_t wcol[3]={kRed+1,kGreen+1,4};
Int_t wmark[3]={22,23,26};

void ComputeAndWriteEff(TList* l, TString dCase, TString var3="Mult");
Bool_t ReadConfig(TString configName);

void ComputeEfficiencyFromCombinHF(){


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
  TString fileWeightName="trackletsWeightsMultInt_LHC13d3_08092014.root";
  TString histoWeightName[3];
  histoWeightName[0]="hNtrUnCorrEvWithCand";
  histoWeightName[1]="hNtrUnCorrEvWithD";
  histoWeightName[2]="hNtrUnCorrEvSel";
  for(Int_t iw=0;iw<3; iw++) hWeight[iw]=0x0;

  if(gSystem->AccessPathName(fileWeightName.Data())==0){
    printf("Open file with multiplicity weights\n");
    TFile* filw = new TFile(fileWeightName.Data());
    for(Int_t iw=0;iw<3; iw++){
      hWeight[iw]=(TH1F*)filw->Get(histoWeightName[iw].Data());
      hWeight[iw]->SetLineColor(wcol[iw]);
      hWeight[iw]->SetStats(0);
      hWeight[iw]->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
      hWeight[iw]->GetYaxis()->SetTitle("Weight");
      hWeight[iw]->GetYaxis()->SetTitleOffset(1.4);
      hWeight[iw]->SetMaximum(3.);
      hWeight[iw]->SetMinimum(0.);
    }
  }


  // pt weights  
  if(ptWeight==kFONLL5overLHC13d3){
    funcPtWeight=new TF1("funcPtWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,30.);
    funcPtWeight->SetParameters(2.94999e+00,3.47032e+00,2.81278e+00,2.5,1.93370e-02,3.86865e+00,-1.54113e-01,8.86944e-02,2.56267e-02);
  }else if (ptWeight==kFONLL7overLHC10f6a){
    funcPtWeight=new TF1("funcPtWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
    funcPtWeight->SetParameters(2.41522e+01,4.92146e+00,6.72495e+00,2.5,6.15361e-03,4.78995e+00,-4.29135e-01,3.99421e-01,-1.57220e-02);
  }else if (ptWeight==kFONLL7overLHC10f7a){
    funcPtWeight=new TF1("funcPtWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
    funcPtWeight->SetParameters(3.59525,2.67871,2.70402,1.72578,4.78167e-03,4.90992,-1.26424e-01,8.21269e-02,-1.26425e-01);
  }else if(ptWeight==kFLAToverLHC10f7a){
    funcPtWeight=new TF1("funcPtWeight","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.15,40.);
    funcPtWeight->SetParameters(1.99498e-01,-9.90532e-02,3.03645e-02,7.42483e-04);
  }else{
    funcPtWeight=new TF1("funcPtWeight","[0]");
    funcPtWeight->SetParameter(0,1.);
  }
  if(ptBWeight==kFONLL5overLHC19c3){
    funcPtBWeight=new TF1("ff","[0]+[1]*x+[2]*TMath::Exp(-(x-[3])*(x-[3])/2/[4]/[4])+[5]*TMath::Exp(-(x-[6])*(x-[6])/2/[7]/[7])+[8]*x*x",0.,50.);
    funcPtBWeight->SetParameters(5.359e-01,-1.921e-02,7.247e-01,6.899e+00,5.310e+00,2.273e-01,3.474e+00,1.444e+00,1.800e-04);
  }else if(ptBWeight==kFONLL5overLHC20g2){
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
  hAccToy->SetMarkerStyle(25);
  hAccToy->SetMarkerColor(kGreen+2);
  hAccToy->SetStats(0);
  hAccToyFine->SetLineColor(kGreen+2);

  TFile* out=new TFile(Form("outputEff%s.root",suffix.Data()),"recreate");
  hAccToy->Write();
  out->Close();



  TH2F* hEventMultZv=(TH2F*)l->FindObject("hEventMultZv");
  if(hEventMultZv){
    TH2F* hEventMultZvEvSel=(TH2F*)l->FindObject("hEventMultZvEvSel");
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
  
  ComputeAndWriteEff(l,"Prompt","Mult");
  if(ptBWeight==kNoPtBWei){
    ComputeAndWriteEff(l,"Feeddw","Mult");
  }else{
    ComputeAndWriteEff(l,"Feeddw","PtB");
  }
  
  TFile* outup=new TFile(Form("outputEff%s.root",suffix.Data()),"update");
  outup->ls();
  TH1D* hEffPr=(TH1D*)outup->Get("hEffPromptVsPtNoWeight");
  TH1D* hEffFd=(TH1D*)outup->Get("hEffFeeddwVsPtNoWeight");
  if(ptBWeight!=kNoPtBWei) hEffFd=(TH1D*)outup->Get("hEffFeeddwVsPtPtBWeight");
  hEffFd->SetLineColor(kGray+1);
  hEffFd->SetMarkerColor(kGray+1);
  hEffFd->SetMarkerStyle(24);

  TCanvas* cpf=new TCanvas("cpf","Prompt vs Feeddown",1200,600);
  cpf->Divide(2,1);
  cpf->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hAccToy->GetYaxis()->SetTitle("Efficiency, acceptance");
  hAccToy->GetYaxis()->SetTitleOffset(1.3);
  hAccToy->Draw();
  hAccToy->SetMinimum(0);
  hEffPr->DrawCopy("same");
  hEffFd->DrawCopy("same");
  TH1D* hEffD=(TH1D*)hEffPr->Clone("hEffD");
  hEffD->Reset("imen");
  TH1D* hEffB=(TH1D*)hEffFd->Clone("hEffB");
  hEffB->Reset("imen");
  TH1D* hRatioEff=(TH1D*)hEffFd->Clone("hRatioEff");
  hRatioEff->Reset("imen");

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
    hRatioEff->SetBinContent(iBin,r);
    hRatioEff->SetBinError(iBin,er);
    
  }

  hEffD->SetMarkerStyle(21);
  hEffD->SetMarkerColor(2);
  hEffD->SetLineColor(2);
  hEffD->DrawCopy("same");
  hEffB->SetMarkerStyle(25);
  hEffB->SetMarkerColor(kRed+1);
  hEffB->SetLineColor(kRed+1);
  hEffB->DrawCopy("same");

  TLegend* legpf=new TLegend(0.6,0.16,0.89,0.36);
  legpf->AddEntry(Form("%s_copy",hEffPr->GetName()),"Effic. prompt","P");
  legpf->AddEntry(Form("%s_copy",hEffFd->GetName()),"Effic. feeddown","P");
  legpf->AddEntry(hAccToy,"Acceptance","P");
  legpf->AddEntry(Form("%s_copy",hEffD->GetName()),"Acc x eff prompt","P");
  legpf->AddEntry(Form("%s_copy",hEffB->GetName()),"Acc x eff feeddown","P");
  legpf->Draw();
  cpf->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hRatioEff->SetStats(0);
  hRatioEff->SetMarkerColor(kBlue+1);
  hRatioEff->SetLineColor(kBlue+1);
  hRatioEff->SetMarkerStyle(20);
  hRatioEff->SetMinimum(0.9);
  hRatioEff->SetMaximum(1.1);
  hRatioEff->GetYaxis()->SetTitle("Ratio efficiency feeddown/prompt");
  hRatioEff->GetYaxis()->SetTitleOffset(1.3);
  hRatioEff->DrawCopy();
  cpf->SaveAs(Form("figures/EfficVsPt_PromptFd_%s.eps",suffix.Data()));

  outup->cd();  
  hEffD->Write();
  hEffB->Write();
  outup->Close();



}

void ComputeAndWriteEff(TList* l, TString dCase, TString var3){

  TH3F* hPtVsYVsVar3Reco=(TH3F*)l->FindObject(Form("hPtVsYVs%sReco%s",var3.Data(),dCase.Data()));
  if(var3=="PtB" && !hPtVsYVsVar3Reco){
    printf("WARNING: pt(B) weight cannot be applied because histos are not there -> re-run with more recent version of the task\n");
    printf(" ---> Resort to mult weights\n");
    var3="Mult";
    hPtVsYVsVar3Reco=(TH3F*)l->FindObject(Form("hPtVsYVs%sReco%s",var3.Data(),dCase.Data()));
  }
  TH3F* hPtVsYVsVar3GenAccEvSel=(TH3F*)l->FindObject(Form("hPtVsYVs%sGenAccEvSel%s",var3.Data(),dCase.Data()));
  TH3F* hPtVsYVsVar3GenAcc=(TH3F*)l->FindObject(Form("hPtVsYVs%sGenAcc%s",var3.Data(),dCase.Data()));
  TH3F* hPtVsYVsVar3GenLimAcc=(TH3F*)l->FindObject(Form("hPtVsYVs%sGenLimAcc%s",var3.Data(),dCase.Data()));
  TString zTitle=hPtVsYVsVar3GenLimAcc->GetZaxis()->GetTitle();
  if(var3=="PtB"){
    zTitle="B-hadron p_{T} (GeV/c)";
    maxMult=-1;
  }
  
  TH2D* hypt=(TH2D*)hPtVsYVsVar3GenAcc->Project3D("yx");
  hypt->SetTitle(Form("Generated in acceptance, all multiplcities, %s",dCase.Data()));
  TH2D* hptmult=(TH2D*)hPtVsYVsVar3GenLimAcc->Project3D("xz");
  hptmult->SetTitle(Form("Generated in |y|<0.5, %s",dCase.Data()));
  hptmult->SetStats(0);
  hypt->SetStats(0);
  hypt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hypt->GetYaxis()->SetTitle("y");
  hptmult->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  hptmult->GetXaxis()->SetTitle(zTitle.Data());
  if(maxMult>0) hptmult->GetXaxis()->SetRangeUser(0.,maxMult);
  TProfile* hMeanPtVar3=hptmult->ProfileX(Form("hMeanPt%s%s",var3.Data(),dCase.Data()));

  TCanvas* c0=new TCanvas(Form("c0%s",dCase.Data()),Form("%s - 2D plots",dCase.Data()),1200,600);
  c0->Divide(2,1);
  c0->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  hypt->Draw("colz");
  c0->cd(2);
  gPad->SetLogz();
  gPad->SetRightMargin(0.13);
  hptmult->Draw("colz");
  hMeanPtVar3->Draw("same");

  TH1D* hPtReco=hPtVsYVsVar3Reco->ProjectionX(Form("hPtReco%s",dCase.Data()));
  TH1D* hPtGenAcc=hPtVsYVsVar3GenAcc->ProjectionX(Form("hPtGenAcc%s",dCase.Data()));
  TH1D* hPtGenAccEvSel=0x0;
  if(hPtVsYVsVar3GenAccEvSel) hPtGenAccEvSel=hPtVsYVsVar3GenAccEvSel->ProjectionX(Form("hPtGenAccEvSel%s",dCase.Data()));
  TH1D* hPtGenLimAcc=hPtVsYVsVar3GenLimAcc->ProjectionX(Form("hPtGenLimAcc%s",dCase.Data()));
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
  TH1D* hEvSelEffVsPt=0x0;
  if(hPtGenAccEvSel){
    hEvSelEffVsPt=(TH1D*)hPtGenAccEvSel->Clone(Form("hEvSelEff%s",dCase.Data()));
    hEvSelEffVsPt->Divide(hPtGenAccEvSel,hPtGenAcc,1,1,"B");
    hEvSelEffVsPt->SetStats(0);
  }


  TCanvas* c1a=new TCanvas(Form("c1a%s",dCase.Data()),Form("%s - AccVsPt",dCase.Data()),1200,600);
  c1a->Divide(2,1);
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

  
  TCanvas* c1e=new TCanvas(Form("c1e%s",dCase.Data()),Form("%s - EffVsPt",dCase.Data()),1200,600);
  c1e->Divide(2,1);
  c1e->cd(1);
  gPad->SetLogy();
  hPtGenLimAcc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtGenLimAcc->GetYaxis()->SetTitle("Entries");
  hPtGenLimAcc->SetLineColor(1);
  hPtGenLimAcc->Draw();
  gPad->Update();
  TPaveStats *st11=(TPaveStats*)hPtGenLimAcc->GetListOfFunctions()->FindObject("stats");
  st11->SetY1NDC(0.7);
  st11->SetY2NDC(0.89);
  hPtGenAcc->SetLineColor(2);
  hPtGenAcc->Draw("sames");
  gPad->Update();
  TPaveStats *st12=(TPaveStats*)hPtGenAcc->GetListOfFunctions()->FindObject("stats");
  st12->SetY1NDC(0.5);
  st12->SetY2NDC(0.69);
  st12->SetTextColor(2);
  if(hPtGenAccEvSel){
    hPtGenAccEvSel->SetLineColor(6);
    hPtGenAccEvSel->Draw("sames");
    gPad->Update();
    TPaveStats *st12s=(TPaveStats*)hPtGenAccEvSel->GetListOfFunctions()->FindObject("stats");
    st12s->SetY1NDC(0.3);
    st12s->SetY2NDC(0.49);
    st12s->SetTextColor(6);
  }
  hPtReco->SetLineColor(4);
  hPtReco->Draw("sames");
  gPad->Update();
  TPaveStats *st3=(TPaveStats*)hPtReco->GetListOfFunctions()->FindObject("stats");
  st3->SetY1NDC(0.1);
  st3->SetY2NDC(0.29);
  st3->SetTextColor(4);
  gPad->Modified();
  c1e->cd(2);
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
  c1e->SaveAs(Form("figures/EfficVsPt_%s_%s.eps",suffix.Data(),dCase.Data()));

  TH1D* hVar3RecoAllPt=hPtVsYVsVar3Reco->ProjectionZ(Form("hVar3Reco%s",dCase.Data()));
  TH1D* hVar3GenAccAllPt=hPtVsYVsVar3GenAcc->ProjectionZ(Form("hVar3GenAcc%s",dCase.Data()));
  TH1D* hVar3GenAccEvSelAllPt=0x0;
  if(hPtVsYVsVar3GenAccEvSel) hVar3GenAccEvSelAllPt=hPtVsYVsVar3GenAccEvSel->ProjectionZ(Form("hVar3GenAccEvSel%s",dCase.Data()));
  TH1D* hVar3GenLimAccAllPt=hPtVsYVsVar3GenLimAcc->ProjectionZ(Form("hVar3GenLimAcc%s",dCase.Data()));
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

  if(var3=="Mult"){
    // double differential acceptance plot pt/mult
    Int_t colMult[20]={kMagenta+1,kMagenta,kBlue+1,kBlue,kBlue-9,
		       kGreen+2,kGreen+1,kGreen,kYellow+1,kYellow,
		       kOrange+2,kOrange+1,kRed-9,kRed,kRed+1};
    TH1D* hAccVsPtMultBin[100];
    TH1D* hPtGenLimAccMultBin[100];
    TH1D* hPtGenAccMultBin[100];
    TGraphErrors* gMeanPtGenLimAccVsMult=new TGraphErrors(0);
    TGraphErrors* gMeanPtGenAccVsMult=new TGraphErrors(0);
    TCanvas* c2dch=new TCanvas(Form("c2dch%s",dCase.Data()),Form("%s - AccVsPt and %s",dCase.Data(),var3.Data()),1400,800);
    c2dch->Divide(3,2);
    Int_t nhp=0;
    TLegend* leg=new TLegend(0.1,0.1,0.6,0.9);
    for(Int_t iBinm=0; iBinm<hPtVsYVsVar3GenAcc->GetNbinsZ(); iBinm++){
      hPtGenAccMultBin[iBinm]=(TH1D*)hPtVsYVsVar3GenAcc->ProjectionX(Form("hPtGenAccMultBin%d",iBinm),0,-1,iBinm+1,iBinm+1);
      hPtGenLimAccMultBin[iBinm]=(TH1D*)hPtVsYVsVar3GenLimAcc->ProjectionX(Form("hPtGenLimAccMultBin%d",iBinm),0,-1,iBinm+1,iBinm+1);
      hAccVsPtMultBin[iBinm]=(TH1D*)hPtGenAccMultBin[iBinm]->Clone(Form("hAccVsPtMultBin%d",iBinm));
      hAccVsPtMultBin[iBinm]->Divide(hPtGenAccMultBin[iBinm],hPtGenLimAccMultBin[iBinm],1,1,"B");
      Double_t minmul=hPtVsYVsVar3GenAcc->GetZaxis()->GetBinLowEdge(iBinm+1);
      Double_t maxmul=hPtVsYVsVar3GenAcc->GetZaxis()->GetBinUpEdge(iBinm+1);
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
    gMeanPtGenLimAccVsMult->GetXaxis()->SetTitle(hPtVsYVsVar3GenAcc->GetZaxis()->GetTitle());
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
    hCopy->GetXaxis()->SetTitle(hPtVsYVsVar3GenAcc->GetZaxis()->GetTitle());
    hCopy->GetYaxis()->SetTitle("Acceptance (p_{T} integrated)");
    hCopy->Draw();

    // TH1D* hAccVsMultPtBin[10];
    // for(Int_t iBinp=0; iBinp<10; iBinp++){
    //   Double_t minptBin=0.5+iBinp;
    //   Double_t maxptBin=minptBin+0.1;
    //   Int_t binMinPt=hPtVsYVsVar3GenAcc->GetXaxis()->FindBin(minptBin+0.00001);
    //   Int_t binMaxPt=hPtVsYVsVar3GenAcc->GetXaxis()->FindBin(maxptBin-0.00001);
    //   TH1D* htmpm1=(TH1D*)hPtVsYVsVar3GenAcc->ProjectionZ(Form("hMultGenAccPtBin%d",iBinp),binMinPt,binMaxPt,0,-1);
    //   TH1D* htmpm2=(TH1D*)hPtVsYVsVar3GenLimAcc->ProjectionZ(Form("hMultGenLimAccPtBin%d",iBinp),binMinPt,binMaxPt,0,-1);
    //   hAccVsMultPtBin[iBinp]=(TH1D*)htmpm1->Clone(Form("hAccVsMultPtBin%d",iBinp));
    //   hAccVsMultPtBin[iBinp]->Divide(htmpm1,htmpm2,1,1,"B");
    //   minptBin=hPtVsYVsVar3GenAcc->GetXaxis()->GetBinLowEdge(binMinPt);
    //   maxptBin=hPtVsYVsVar3GenAcc->GetXaxis()->GetBinUpEdge(binMaxPt);
    //   hAccVsMultPtBin[iBinp]->SetTitle(Form("%.2f<pt<%.2f",minptBin,maxptBin));
    //   hAccVsMultPtBin[iBinp]->GetXaxis()->SetTitle(hPtVsYVsVar3GenAcc->GetZaxis()->GetTitle());
    //   hAccVsMultPtBin[iBinp]->GetYaxis()->SetTitle("Acceptance");
    //   hAccVsMultPtBin[iBinp]->SetStats(0);
    //   hAccVsMultPtBin[iBinp]->SetLineWidth(2);
    //   hAccVsMultPtBin[iBinp]->GetXaxis()->SetRangeUser(1000.,5000.);
    // }
    // TCanvas* c2dch2=new TCanvas(Form("c2dch2%s",dCase.Data()),Form("%s - AccVs%s and Pt",dCase.Data(),var3.Data()),1400,800);
    // c2dch2->Divide(5,2);
    // for(Int_t iBinp=0; iBinp<10; iBinp++){
    //   c2dch2->cd(iBinp+1);
    //   hAccVsMultPtBin[iBinp]->Draw();
    // }
  }

  TCanvas* c2a=new TCanvas(Form("c2a%s",dCase.Data()),Form("%s - AccVs%s",dCase.Data(),var3.Data()),1200,600);
  c2a->Divide(2,1);
  c2a->cd(1);
  gPad->SetLogy();
  hVar3GenLimAccAllPt->SetLineColor(1);
  hVar3GenLimAccAllPt->GetXaxis()->SetTitle(zTitle.Data());
  hVar3GenLimAccAllPt->GetYaxis()->SetTitle("Entries");
  if(maxMult>0) hVar3GenLimAccAllPt->GetXaxis()->SetRangeUser(0.,maxMult);
  hVar3GenLimAccAllPt->Draw();
  gPad->Update();
  TPaveStats *st21=(TPaveStats*)hVar3GenLimAccAllPt->GetListOfFunctions()->FindObject("stats");
  st21->SetY1NDC(0.7);
  st21->SetY2NDC(0.89);
  hVar3GenAccAllPt->SetLineColor(2);
  hVar3GenAccAllPt->Draw("sames");
  gPad->Update();
  TPaveStats *st22=(TPaveStats*)hVar3GenAccAllPt->GetListOfFunctions()->FindObject("stats");
  st22->SetY1NDC(0.5);
  st22->SetY2NDC(0.69);
  st22->SetTextColor(2);
  if(hVar3GenAccEvSelAllPt){
    hVar3GenAccEvSelAllPt->SetLineColor(6);
    hVar3GenAccEvSelAllPt->Draw("sames");
    gPad->Update();
    TPaveStats *st22s=(TPaveStats*)hVar3GenAccEvSelAllPt->GetListOfFunctions()->FindObject("stats");
    st22s->SetY1NDC(0.3);
    st22s->SetY2NDC(0.49);
    st22s->SetTextColor(6);
  }
  hVar3RecoAllPt->SetLineColor(4);
  hVar3RecoAllPt->Draw("sames");
  gPad->Update();
  TPaveStats *st23=(TPaveStats*)hVar3RecoAllPt->GetListOfFunctions()->FindObject("stats");
  st23->SetY1NDC(0.1);
  st23->SetY2NDC(0.29);
  st23->SetTextColor(4);
  gPad->Modified();
  c2a->cd(2);
  hEffVsVar3AllPt->SetLineColor(4);
  hEffVsVar3AllPt->SetMinimum(0);
  hEffVsVar3AllPt->SetMaximum(1.6);
  hEffVsVar3AllPt->GetXaxis()->SetTitle(zTitle.Data());
  hEffVsVar3AllPt->GetYaxis()->SetTitle("Ratio");
  if(maxMult>0) hEffVsVar3AllPt->GetXaxis()->SetRangeUser(0.,maxMult);
  hEffVsVar3AllPt->Draw();
  hAccVsVar3AllPt->SetLineColor(2);
  hAccVsVar3AllPt->Draw("same");
  //  hAccEffVsVar3AllPt->SetLineColor(6);
  // hAccEffVsVar3AllPt->Draw("same");
  TLatex* tacc2=new TLatex(0.16,0.8,"Acceptance (CombinHF)");
  tacc2->SetNDC();
  tacc2->SetTextColor(hAccVsVar3AllPt->GetLineColor());
  tacc2->Draw();
  TLatex* te2=new TLatex(0.16,0.72,"Efficiency");
  te2->SetNDC();
  te2->SetTextColor(hEffVsVar3AllPt->GetLineColor());
  te2->Draw();

  TCanvas* c2e=new TCanvas(Form("c2e%s",dCase.Data()),Form("%s - EffVs%s",dCase.Data(),var3.Data()),1200,600);
  c2e->Divide(2,1);
  c2e->cd(1);
  gPad->SetLogy();
  hVar3GenLimAccAllPt->GetXaxis()->SetTitle(zTitle.Data());
  hVar3GenLimAccAllPt->GetYaxis()->SetTitle("Entries");
  if(maxMult>0) hVar3GenLimAccAllPt->GetXaxis()->SetRangeUser(0.,maxMult);
  hVar3GenLimAccAllPt->SetLineColor(1);
  hVar3GenLimAccAllPt->Draw();
  gPad->Update();
  TPaveStats *st1e=(TPaveStats*)hVar3GenLimAccAllPt->GetListOfFunctions()->FindObject("stats");
  st1e->SetY1NDC(0.7);
  st1e->SetY2NDC(0.89);
  hVar3GenAccAllPt->SetLineColor(2);
  hVar3GenAccAllPt->Draw("sames");
  gPad->Update();
  TPaveStats *st2e=(TPaveStats*)hVar3GenAccAllPt->GetListOfFunctions()->FindObject("stats");
  st2e->SetY1NDC(0.5);
  st2e->SetY2NDC(0.69);
  st2e->SetTextColor(2);
  if(hVar3GenAccEvSelAllPt){
    hVar3GenAccEvSelAllPt->SetLineColor(6);
    hVar3GenAccEvSelAllPt->Draw("sames");
    gPad->Update();
    TPaveStats *st2se=(TPaveStats*)hVar3GenAccEvSelAllPt->GetListOfFunctions()->FindObject("stats");
    st2se->SetY1NDC(0.3);
    st2se->SetY2NDC(0.49);
    st2se->SetTextColor(6);
  }
  hVar3RecoAllPt->SetLineColor(4);
  hVar3RecoAllPt->Draw("sames");
  gPad->Update();
  TPaveStats *st3s=(TPaveStats*)hVar3RecoAllPt->GetListOfFunctions()->FindObject("stats");
  st3s->SetY1NDC(0.1);
  st3s->SetY2NDC(0.29);
  st3s->SetTextColor(4);
  gPad->Modified();
  c2e->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  if(hVar3GenAccEvSelAllPt){
    hEvSelEffVsVar3AllPt->SetLineColor(6);
    hEvSelEffVsVar3AllPt->SetMinimum(0);
    hEvSelEffVsVar3AllPt->SetMaximum(1.05);
    hEvSelEffVsVar3AllPt->GetXaxis()->SetTitle(zTitle.Data());
    hEvSelEffVsVar3AllPt->GetYaxis()->SetTitle("Ratio");
    if(maxMult>0) hEvSelEffVsVar3AllPt->GetXaxis()->SetRangeUser(0.,maxMult);
    hEvSelEffVsVar3AllPt->Draw();
    t1->Draw();
  }
  hEffVsVar3AllPt->SetLineColor(4);
  hEffVsVar3AllPt->Draw("same");
  t3->Draw();
  c2e->SaveAs(Form("figures/EfficVs%s_%s_%s.eps",var3.Data(),suffix.Data(),dCase.Data()));

  if(var3=="Mult"){
  
    const Int_t nMultBins=6;
    Double_t mulLims[nMultBins+1]={0.,5.,12.,20.,40.,80.,200.};
    if(hPtVsYVsVar3Reco->GetZaxis()->GetXmax()>300. && maxMult<0){
      for(Int_t jb=0; jb<=nMultBins; jb++) mulLims[jb]=(Double_t)jb/nMultBins*hPtVsYVsVar3Reco->GetZaxis()->GetXmax();
    }
    Int_t mulcol[nMultBins]={1,kRed+1,kGreen+2,4,kOrange+2,kMagenta};

    TH1D* hPtRecoM[nMultBins];
    TH1D* hPtGenAccM[nMultBins];
    //  TH1D* hPtGenLimAccM[nMultBins];
    TH1D* hEffVsPtM[nMultBins];
    for(Int_t j=0; j<nMultBins; j++){
      Int_t lowBin=hPtVsYVsVar3Reco->GetZaxis()->FindBin(mulLims[j]);
      Int_t hiBin=hPtVsYVsVar3Reco->GetZaxis()->FindBin(mulLims[j+1]-0.001);
      //    printf("%d (%f)  %d(%f)\n",lowBin,hPtVsYVsMultReco->GetZaxis()->GetBinLowEdge(lowBin),hiBin,hPtVsYVsMultReco->GetZaxis()->GetBinUpEdge(hiBin));
    
      hPtRecoM[j]=hPtVsYVsVar3Reco->ProjectionX(Form("hPtRecoM%d",j),0,-1,lowBin,hiBin);
      hPtGenAccM[j]=hPtVsYVsVar3GenAcc->ProjectionX(Form("hPtGenAccM%d",j),0,-1,lowBin,hiBin);
      //    hPtGenLimAccM[j]=hPtVsYVsVar3GenLimAcc->ProjectionX(Form("hPtGenLimAccM%d",j),0,-1,lowBin,hiBin);
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
    TH1D* hEffVsMult[nPtBins];
    for(Int_t j=0; j<nPtBins; j++){
      Int_t lowBin=hPtVsYVsVar3Reco->GetXaxis()->FindBin(binLims[j]);
      Int_t hiBin=hPtVsYVsVar3Reco->GetXaxis()->FindBin(binLims[j+1]-0.001);
      //printf("%d (%f)  %d(%f)\n",lowBin,hPtVsYVsVar3Reco->GetXaxis()->GetBinLowEdge(lowBin),hiBin,hPtVsYVsVar3Reco->GetXaxis()->GetBinUpEdge(hiBin));
    
      hMultReco[j]=hPtVsYVsVar3Reco->ProjectionZ(Form("hMultReco%s%d",dCase.Data(),j),lowBin,hiBin);
      hMultGenAcc[j]=hPtVsYVsVar3GenAcc->ProjectionZ(Form("hMultGenAcc%s%d",dCase.Data(),j),lowBin,hiBin);
      //    hMultGenLimAcc[j]=hPtVsYVsVar3GenLimAcc->ProjectionZ(Form("hMultGenLimAcc%d",j),lowBin,hiBin);
      hEffVsMult[j]=(TH1D*)hMultReco[j]->Clone(Form("hEff%s%d",dCase.Data(),j));
      hEffVsMult[j]->Divide(hMultReco[j],hMultGenAcc[j],1,1,"B");
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
    hEffVsMult[0]->GetYaxis()->SetTitle("Efficiency");
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
    if(hWeight[0]){
      gPad->SetLeftMargin(0.12);
      if(maxMult>0) hWeight[0]->GetXaxis()->SetRangeUser(0.,maxMult);
      hWeight[0]->Draw();
      TLegend* legw=new TLegend(0.2,0.65,0.55,0.89);
      legw->SetFillStyle(0);
      legw->AddEntry(hWeight[0],"Candidate Weight","L");
      if(hWeight[1]){
	hWeight[1]->Draw("same");
	legw->AddEntry(hWeight[1],"D Weight","L");
      }
      if(hWeight[2]){ 
	hWeight[2]->Draw("same");
	legw->AddEntry(hWeight[2],"EvSel Weight","L");
      }
      legw->Draw();
    }
    cw->SaveAs(Form("figures/Effic%sVsMultPtBins_%s.eps",dCase.Data(),suffix.Data()));
  }else{
    TCanvas* cwb=new TCanvas("cwb",Form("%s - ptB weight",dCase.Data()),800,800);
    funcPtBWeight->SetTitle("");
    funcPtBWeight->GetXaxis()->SetTitle("B-hadron p_{T} (GeV/c)");
    funcPtBWeight->GetYaxis()->SetTitle("Weight (FONLL/MC)");
    funcPtBWeight->Draw();
  }
  
  TH1D *hpteffNoWeight =new TH1D(Form("hEff%sVsPtNoWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D **hpteffMultWeight =new TH1D*[3];
  hpteffMultWeight[0]=new TH1D(Form("hEff%sVsPtCandWeight",dCase.Data()),"",nPtBins,binLims);
  hpteffMultWeight[1]=new TH1D(Form("hEff%sVsPtDWeight",dCase.Data()),"",nPtBins,binLims);
  hpteffMultWeight[2]=new TH1D(Form("hEff%sVsPtEvSelWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D **hpteffMultPtWeight =new TH1D*[3];
  hpteffMultPtWeight[0]=new TH1D(Form("hEff%sVsPtCandFONLLWeight",dCase.Data()),"",nPtBins,binLims);
  hpteffMultPtWeight[1]=new TH1D(Form("hEff%sVsPtDFONLLWeight",dCase.Data()),"",nPtBins,binLims);
  hpteffMultPtWeight[2]=new TH1D(Form("hEff%sVsPtEvSelFONLLWeight",dCase.Data()),"",nPtBins,binLims);
  TH1D *hpteffPtBWeight =new TH1D(Form("hEff%sVsPtPtBWeight",dCase.Data()),"",nPtBins,binLims);

  TH1D* hMultRecoAllPtW=hPtVsYVsVar3Reco->ProjectionZ(Form("hMultRecoW%s",dCase.Data()));
  TH1D* hMultGenAccAllPtW=hPtVsYVsVar3GenAcc->ProjectionZ(Form("hMultGenAccW%s",dCase.Data()));
  TH1D* hMultGenLimAccAllPtW=hPtVsYVsVar3GenLimAcc->ProjectionZ(Form("hMultGenLimAccW%s",dCase.Data()));
  hMultRecoAllPtW->Reset("ines");
  hMultGenAccAllPtW->Reset("ines");
  hMultGenLimAccAllPtW->Reset("ines");

  Double_t countNumer[nPtBins];
  Double_t countDenom[nPtBins];
  Double_t countNumerPtBWei[nPtBins];
  Double_t countDenomPtBWei[nPtBins];

  Double_t countNumerMultWei[nPtBins][3];
  Double_t countDenomMultWei[nPtBins][3];
  Double_t countNumerMultPtWei[nPtBins][3];
  Double_t countDenomMultPtWei[nPtBins][3];

  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){
    countNumer[iPtBin]=0.0;
    countDenom[iPtBin]=0.0;
    countNumerPtBWei[iPtBin]=0.0;
    countDenomPtBWei[iPtBin]=0.0;
    for(Int_t iw=0; iw<3; iw++){
      countNumerMultWei[iPtBin][iw]=0.0;
      countDenomMultWei[iPtBin][iw]=0.0;
      countNumerMultPtWei[iPtBin][iw]=0.0;
      countDenomMultPtWei[iPtBin][iw]=0.0;
    }
  }

  for(Int_t ibx=0; ibx<=hPtVsYVsVar3GenAcc->GetNbinsX()+1; ibx++){
    for(Int_t iby=0; iby<=hPtVsYVsVar3GenAcc->GetNbinsY()+1; iby++){
      for(Int_t ibz=0; ibz<=hPtVsYVsVar3GenAcc->GetNbinsZ()+1; ibz++){
	Double_t pt=hPtVsYVsVar3GenAcc->GetXaxis()->GetBinCenter(ibx);
	//	Double_t y=hPtVsYVsVar3GenAcc->GetYaxis()->GetBinCenter(iby);
	Double_t v3=hPtVsYVsVar3GenAcc->GetZaxis()->GetBinCenter(ibz);
	//	  printf("%f %f %f\n",pt,y,mult);
	Int_t jPtBin=TMath::BinarySearch(nPtBins+1,binLims,pt);
	if(jPtBin>=0 && jPtBin<nPtBins){
	  Double_t crec=hPtVsYVsVar3Reco->GetBinContent(ibx,iby,ibz);
	  Double_t cgen=hPtVsYVsVar3GenAcc->GetBinContent(ibx,iby,ibz);
	  Double_t wpt=funcPtWeight->Eval(pt);
	  countNumer[jPtBin]+=crec;
	  countDenom[jPtBin]+=cgen;
	  if(var3=="PtB"){
	    Double_t w=funcPtBWeight->Eval(v3);
	    countNumerPtBWei[jPtBin]+=crec*w;
	    countDenomPtBWei[jPtBin]+=cgen*w;
	    
	  }else{
	    for(Int_t iw=0; iw<3; iw++){
	      Double_t w=0;
	      if(hWeight[iw]){
		Int_t binw=hWeight[iw]->FindBin(v3);
		if(binw>=1 && binw<hWeight[iw]->GetNbinsX()+1){
		  w=hWeight[iw]->GetBinContent(binw);
		  //		printf("mult %.0f   bin %d   wei %f\n",mult,binw,w);
		}else{
		  if(cgen>0){
		    printf("mult %.0f   bin %d   wei %f\n",v3,binw,w);
		    getchar();
		  }
		}
	      }else{
		w=1;
	      }
	      if(!useMultWeight) w=1;
	      countNumerMultWei[jPtBin][iw]+=crec*w;
	      countDenomMultWei[jPtBin][iw]+=cgen*w;
	      countNumerMultPtWei[jPtBin][iw]+=crec*w*wpt;
	      countDenomMultPtWei[jPtBin][iw]+=cgen*w*wpt;
	      hMultRecoAllPtW->Fill(v3,crec*w);
	      hMultGenAccAllPtW->Fill(v3,cgen*w);
	      hMultGenLimAccAllPtW->Fill(v3,hPtVsYVsVar3GenLimAcc->GetBinContent(ibx,iby,ibz)*w);
	    }
	  }
	}
      }
    }
  }

  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){
    Double_t eff1=countNumer[iPtBin]/countDenom[iPtBin];
    Double_t erreff1=TMath::Sqrt(eff1*(1-eff1)/countDenom[iPtBin]);
    printf("---- Pt range %.0f - %.0f ----\n",binLims[iPtBin],binLims[iPtBin+1]);
    printf("Eff from Projection = %f/%f = %f+-%f\n",hPtRecoR->GetBinContent(iPtBin+1),hPtGenAccR->GetBinContent(iPtBin+1),hEffVsPtR->GetBinContent(iPtBin+1),hEffVsPtR->GetBinError(iPtBin+1));
    printf("Eff No weights      = %f/%f = %f+-%f\n",countNumer[iPtBin],countDenom[iPtBin],eff1,erreff1);
    hpteffNoWeight->SetBinContent(iPtBin+1,eff1);
    hpteffNoWeight->SetBinError(iPtBin+1,erreff1);
    if(var3=="PtB"){
      Double_t eff4=countNumerPtBWei[iPtBin]/countDenomPtBWei[iPtBin];
      Double_t erreff4=TMath::Sqrt(eff4*(1-eff4)/countDenom[iPtBin]);// countDenom is NOT a typo, it has to be like this to get proper statistical errors from the no-weight case
      printf("Eff With pt(B) weights = %f/%f = %f+-%f\n",countNumerPtBWei[iPtBin],countDenomPtBWei[iPtBin],eff4,erreff4);
      hpteffPtBWeight->SetBinContent(iPtBin+1,eff4);
      hpteffPtBWeight->SetBinError(iPtBin+1,erreff4);
    }else{
      for(Int_t iw=0;iw<3; iw++){
	Double_t eff2=countNumerMultWei[iPtBin][iw]/countDenomMultWei[iPtBin][iw];
	Double_t erreff2=TMath::Sqrt(eff2*(1-eff2)/countDenom[iPtBin]);// countDenom is NOT a typo, it has to be like this to get proper statistical errors from the no-weight case
	printf("Eff With mult weights %d   = %f/%f = %f+-%f\n",iw,countNumerMultWei[iPtBin][iw],countDenomMultWei[iPtBin][iw],eff2,erreff2);
	hpteffMultWeight[iw]->SetBinContent(iPtBin+1,eff2);
	hpteffMultWeight[iw]->SetBinError(iPtBin+1,erreff2);
	Double_t eff3=countNumerMultPtWei[iPtBin][iw]/countDenomMultPtWei[iPtBin][iw];
	Double_t erreff3=TMath::Sqrt(eff3*(1-eff3)/countDenom[iPtBin]);// countDenom is NOT a typo, it has to be like this to get proper statistical errors from the no-weight case
	printf("Eff With mult+pt weights %d   = %f/%f = %f+-%f\n",iw,countNumerMultPtWei[iPtBin][iw],countDenomMultPtWei[iPtBin][iw],eff3,erreff3);
	hpteffMultPtWeight[iw]->SetBinContent(iPtBin+1,eff3);
	hpteffMultPtWeight[iw]->SetBinError(iPtBin+1,erreff3);
      }
    }
  }

  hEffVsPtR->SetMarkerStyle(0);
  hEffVsPtR->SetMarkerColor(0);
  hEffVsPtR->SetMarkerSize(1.2);
  hEffVsPtR->SetLineColor(kGray);
  hEffVsPtR->SetLineWidth(4);
  hEffVsPtR->SetStats(0);
  hpteffNoWeight->SetMarkerStyle(20);
  
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
    hEffVsMultAllPtW->GetYaxis()->SetTitle("Ratio");
    hEffVsMultAllPtW->Draw();
    hAccVsMultAllPtW->SetLineColor(2);
    hAccVsMultAllPtW->Draw("same");
    tacc2->Draw();
    te2->Draw();
    


    TH1F** hRatio=new TH1F*[3];
    for(Int_t iw=0;iw<3; iw++){
      hRatio[iw]=(TH1F*)hpteffMultWeight[iw]->Clone(Form("hRatio%s",dCase.Data()));
      hRatio[iw]->Divide(hpteffMultWeight[iw],hpteffNoWeight);
      hRatio[iw]->GetYaxis()->SetTitle("With Mult weight / Without weight");
      hRatio[iw]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hRatio[iw]->GetYaxis()->SetTitleOffset(1.4);
      hRatio[iw]->SetLineColor(wcol[iw]);
      hRatio[iw]->SetMarkerColor(wcol[iw]);
      hRatio[iw]->SetMarkerStyle(wmark[iw]);
      hRatio[iw]->SetStats(0);
      hpteffMultWeight[iw]->SetLineColor(wcol[iw]);
      hpteffMultWeight[iw]->SetMarkerColor(wcol[iw]);
      hpteffMultWeight[iw]->SetMarkerStyle(wmark[iw]);
    }
    TH1F** hRatioPtW=new TH1F*[3];
    for(Int_t iw=0;iw<3; iw++){
      hRatioPtW[iw]=(TH1F*)hpteffMultPtWeight[iw]->Clone(Form("hRatioWei%s",dCase.Data()));
      hRatioPtW[iw]->Divide(hpteffMultPtWeight[iw],hpteffMultWeight[iw]);
      hRatioPtW[iw]->GetYaxis()->SetTitle("With Mult+p_{T} weight / With only Mult weight");
      hRatioPtW[iw]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hRatioPtW[iw]->GetYaxis()->SetTitleOffset(1.4);
      hRatioPtW[iw]->SetLineColor(wcol[iw]);
      hRatioPtW[iw]->SetMarkerColor(wcol[iw]);
      hRatioPtW[iw]->SetMarkerStyle(wmark[iw]);
      hRatioPtW[iw]->SetStats(0);
      hpteffMultPtWeight[iw]->SetLineColor(wcol[iw]);
      hpteffMultPtWeight[iw]->SetMarkerColor(wcol[iw]);
      hpteffMultPtWeight[iw]->SetMarkerStyle(wmark[iw]);
    }

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
    for(Int_t iw=0;iw<3; iw++) hpteffMultWeight[iw]->Draw("same");
    TLegend* leg=new TLegend(0.55,0.15,0.89,0.45);
    leg->SetFillStyle(0);
    leg->AddEntry(hEffVsPtR,"TH3F::Project","L");
    leg->AddEntry(hpteffNoWeight,"Multiplcity slices - No Weight","PL");
    leg->AddEntry(hpteffMultWeight[0],"Multiplcity slices - Cand Weight","PL");
    leg->AddEntry(hpteffMultWeight[1],"Multiplcity slices - D Weight","PL");
    leg->AddEntry(hpteffMultWeight[2],"Multiplcity slices - EvSel Weight","PL");
    leg->Draw();
    ceff->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetTickx();
    gPad->SetTicky();
    hRatio[0]->Draw();
    hRatio[0]->SetMinimum(0.95);
    hRatio[0]->SetMaximum(1.05);
    for(Int_t iw=1;iw<3; iw++) hRatio[iw]->Draw("same");
    ceff->SaveAs(Form("figures/Effic%sWithMultWeights%s.eps",dCase.Data(),suffix.Data()));

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
    for(Int_t iw=0;iw<3; iw++) hpteffMultPtWeight[iw]->Draw("same");
    TLegend* leg2=new TLegend(0.55,0.15,0.89,0.45);
    leg2->SetFillStyle(0);
    leg2->AddEntry(hEffVsPtR,"TH3F::Project","L");
    leg2->AddEntry(hpteffNoWeight,"Multiplcity slices - No Weight","PL");
    leg2->AddEntry(hpteffMultPtWeight[0],"Multiplicity slices - FONLL+Cand Weight","PL");
    leg2->AddEntry(hpteffMultPtWeight[1],"Multiplicity slices - FONLL+D Weight","PL");
    leg2->AddEntry(hpteffMultPtWeight[2],"Multiplicity slices - FONLL+EvSel Weight","PL");
    leg2->Draw();
    ceff2->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetTickx();
    gPad->SetTicky();
    hRatioPtW[0]->Draw();
    hRatioPtW[0]->SetMinimum(0.95);
    hRatioPtW[0]->SetMaximum(1.05);
    for(Int_t iw=1;iw<3; iw++) hRatioPtW[iw]->Draw("same");
    ceff2->SaveAs(Form("figures/Effic%sWithMultAndPtWeights%s.eps",dCase.Data(),suffix.Data()));
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
  hEffVsPtR->Write();
  hpteffNoWeight->Write();
  hpteffMultWeight[0]->Write();
  hpteffMultWeight[1]->Write();
  hpteffMultWeight[2]->Write();
  hpteffMultPtWeight[0]->Write();
  hpteffMultPtWeight[1]->Write();
  hpteffMultPtWeight[2]->Write();
  hpteffPtBWeight->Write();
  hEvSelEffVsPt->Write();
  out->Close();

}


Bool_t ReadConfig(TString configName){
  FILE* confFil=fopen(configName.Data(),"r");
  char line[50];
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

