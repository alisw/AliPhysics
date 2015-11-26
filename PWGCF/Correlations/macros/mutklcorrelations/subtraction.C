#include <TPad.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include "env.h"

void subtraction(
    //        TString inFileName  = "dphi/dphi_corr_mu_tkl_V0S_f.root",
    //    TString inFileName  = "dphi/dphi_corr_mu_tkl_V0L_f.root",
    //    TString inFileName  = "dphi/dphi_corr_tkl_tkl_V0S_f_vertex05.root",
    //    TString inFileName  = "dphi/dphi_corr_tkl_tkl_V0S_bcde_vertex01.root",
    //    TString inFileName  = "dphi/dphi_corr_tkl_tkl_V0S_bcde.root",
    //    TString inFileName  = "dphi/dphi_corr_mu_tkl_V0S_bcde.root",
    //    TString inFileName  = "dphi/dphi_corr_mu_tkl_AMA_f.root",
    //    TString inFileName  = "dphi/dphi_corr_mu_tkl_LEE_f.root",
    //    Int_t i=0,    // trigger
    //    Int_t j=4,    // associate
//        TString inFileName  = "dphi/dphi_corr_mu_tkl_AMA_bcde.root",
    TString inFileName  = "dphi/dphi_corr_tkl_tkl_DPM_bcde.root",
    Int_t i=0,    // trigger
    Int_t j=0,    // associate
    Int_t c1 = 0,
    Int_t c2 = 3,
//    Float_t etaMin=-5.,
//    Float_t etaMax=-1.5,
//    Float_t exclusionMin = -10,
//    Float_t exclusionMax = +10,
        Float_t etaMin=-2,
        Float_t etaMax=2,
        Float_t exclusionMin = -1.2,
        Float_t exclusionMax = 1.2,
    Bool_t debug = 1,
    Double_t* res = NULL
){
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(2);
  Bool_t isForward = !(exclusionMin>=etaMin && exclusionMin<=etaMax && exclusionMax>=etaMin && exclusionMax<=etaMax);
  Int_t rebinDphi = 6;
  Int_t rebinDeta;

  if(inFileName.Contains("_ET"))rebinDeta = isForward ? 5 : 1;
  else rebinDeta = isForward ? 7 : 1;

  TFile *infile=new TFile(inFileName);
  TH2F* hist1 = (TH2F*) infile->Get(Form("dphi_%d_%d_%d", i, j, c1));
  TH2F* hist2 = (TH2F*) infile->Get(Form("dphi_%d_%d_%d", i, j, c2));

  if (!isForward) {
    for (Int_t h=1;h<=hist1->GetNbinsX();h++){
      for (Int_t k=1;k<=hist1->GetNbinsY();k++){
        hist1->SetBinError(h,k,hist1->GetBinError(h,k)*sqrt(2));
        hist2->SetBinError(h,k,hist2->GetBinError(h,k)*sqrt(2));
      }
    }
  }

  hist1->GetYaxis()->SetRangeUser(etaMin+0.001,etaMax-0.001);
  hist2->GetYaxis()->SetRangeUser(etaMin+0.001,etaMax-0.001);
  hist1->Scale(1./hist1->GetYaxis()->GetBinWidth(1));//only Dphi normalization done in calc_sum_of_ratio
  hist2->Scale(1./hist2->GetYaxis()->GetBinWidth(1));
  hist1->GetXaxis()->SetTitleOffset(1.6);
  hist1->GetYaxis()->SetTitleOffset(1.6);
  hist1->GetZaxis()->SetTitleOffset(2.2);
  hist2->GetXaxis()->SetTitleOffset(1.6);
  hist2->GetYaxis()->SetTitleOffset(1.6);
  hist2->GetZaxis()->SetTitleOffset(2.2);

  TString label1(hist1->GetTitle());
  TString label2(hist2->GetTitle());
  TObjArray* objArray1 = label1.Tokenize("-");
  TObjArray* objArray2 = label2.Tokenize("-");
  printf("label1=%s\n",label1.Data());
  printf("label2=%s\n",label2.Data());

  TString textTrg = objArray1->At(0)->GetName();
  TString textAss = objArray1->At(1)->GetName();
  TString cent11 = objArray1->At(2)->GetName();
  TString cent12 = objArray1->At(3)->GetName();
  TString cent21 = objArray2->At(2)->GetName();
  TString cent22 = objArray2->At(3)->GetName();
  cent11.ReplaceAll(" ","");
  cent22.ReplaceAll(" ","");
  TString textCent1   = Form("%s-%s", cent11.Data(),cent12.Data());
  TString textCent2   = Form("%s-%s", cent21.Data(),cent22.Data());
  TString textCentDif = Form("(%s) - (%s)", textCent1.Data(),textCent2.Data());
  hist1->SetTitle(Form(";%s;%s;%s",kTitlePhi,kTitleEta,kCorrFuncTitle));
  hist2->SetTitle(Form(";%s;%s;%s",kTitlePhi,kTitleEta,kCorrFuncTitle));

  TLatex* latex = 0;
  if (debug) {
    latex = new TLatex();
    latex->SetTextSize(fontSize);
    latex->SetTextFont(42);
    latex->SetTextAlign(21);
    latex->SetNDC();
  }
  hist1->Rebin2D(rebinDphi,rebinDeta);
  hist2->Rebin2D(rebinDphi,rebinDeta);
  hist1->Scale(1./rebinDphi/rebinDeta);
  hist2->Scale(1./rebinDphi/rebinDeta);

  if (debug) {
    TCanvas* cCent1 = new TCanvas("cCent1","cCent1",1000,1000);
    PadFor2DCorr();
    hist1->GetYaxis()->SetRangeUser(etaMin+0.001,etaMax-0.001);
    //    hist1->SetMaximum(1.5*hist1->GetBinContent(hist1->FindBin(3.14159,0))-0.5*hist1->GetMinimum());
    hist1->DrawCopy("surf1");
    latex->DrawLatex(0.25,0.93,kSystemEnergy);
    latex->DrawLatex(0.25,0.86,textCent1.Data());
    latex->DrawLatex(0.80,0.93,textTrg.Data());
    latex->DrawLatex(0.80,0.86,textAss.Data());
    gPad->Print(Form("cent1%i%s",i,format));

    TCanvas* cCent2 = new TCanvas("cCent2","cCent2",1000,1000);
    PadFor2DCorr();
    hist2->GetYaxis()->SetRangeUser(etaMin+0.001,etaMax-0.001);
    //    hist2->SetMaximum(1.5*hist2->GetBinContent(hist2->FindBin(3.14159,0))-0.5*hist2->GetMinimum());
    hist2->DrawCopy("surf1");
    latex->DrawLatex(0.25,0.93,kSystemEnergy);
    latex->DrawLatex(0.25,0.86,textCent2.Data());
    latex->DrawLatex(0.80,0.93,textTrg.Data());
    latex->DrawLatex(0.80,0.86,textAss.Data());
    gPad->Print(Form("cent2%i%s",i,format));
    TFile* f = new TFile(inFileName.ReplaceAll("dphi/dphi_corr","debug"),"recreate");
    hist1->Write();
    hist2->Write();
    f->Close();
    delete f;
  }
  Float_t phiBaselineMin = pi/2 - 0.2 + 0.0001;
  Float_t phiBaselineMax = pi/2 + 0.2 - 0.0001;
  Int_t binEtaMin  = hist1->GetYaxis()->FindBin(etaMin + 0.001);
  Int_t binEtaMax  = hist1->GetYaxis()->FindBin(etaMax - 0.001);
  Int_t binExclusionMin  = hist1->GetYaxis()->FindBin(exclusionMin - 0.001);
  Int_t binExclusionMax  = hist1->GetYaxis()->FindBin(exclusionMax + 0.001);
  Int_t binPhiBaselineMin = hist1->GetXaxis()->FindBin(phiBaselineMin);
  Int_t binPhiBaselineMax = hist1->GetXaxis()->FindBin(phiBaselineMax);
  Float_t actualExclusionMin = hist1->GetYaxis()->GetBinUpEdge(binExclusionMin);
  Float_t actualExclusionMax = hist1->GetYaxis()->GetBinLowEdge(binExclusionMax);

  Double_t baseLineE;
  Float_t baseLine = hist1->IntegralAndError(binPhiBaselineMin, binPhiBaselineMax, binEtaMin, binEtaMax, baseLineE);
  baseLine  /= (binPhiBaselineMax - binPhiBaselineMin + 1) * (binEtaMax - binEtaMin + 1);
  baseLineE /= (binPhiBaselineMax - binPhiBaselineMin + 1) * (binEtaMax - binEtaMin + 1);

  Double_t baseLineE2;
  Float_t baseLine2 = hist2->IntegralAndError(binPhiBaselineMin, binPhiBaselineMax, binEtaMin, binEtaMax, baseLineE2);
  baseLine2  /= (binPhiBaselineMax - binPhiBaselineMin + 1) * (binEtaMax - binEtaMin + 1);
  baseLineE2 /= (binPhiBaselineMax - binPhiBaselineMin + 1) * (binEtaMax - binEtaMin + 1);

  TH1D* proj_c1  = 0;
  TH1D* proj2_c1 = 0;
  TH1D* proj_c2  = 0;
  TH1D* proj2_c2 = 0;

  TF1* v2gaus1 = 0;
  TF1* v2gaus2 = 0;
  if (gStudySystematic==k_LowMultScale || gStudySystematic==k_baselineGausFit) {
    if (!isForward){
      //c1
      proj_c1  = hist1->ProjectionX(Form("%s_proj1x", hist1->GetName()), binEtaMin, binExclusionMin);
      proj2_c1 = hist1->ProjectionX(Form("%s_proj2x", hist1->GetName()), binExclusionMax, binEtaMax);
      proj_c1->Add(proj2_c1, 1);
      proj_c1->Scale(1.0 / (binEtaMax - binExclusionMax + 1 + binExclusionMin - binEtaMin + 1));  
      proj_c2  = hist2->ProjectionX(Form("%s_proj1x", hist2->GetName()), binEtaMin, binExclusionMin);
      proj2_c2 = hist2->ProjectionX(Form("%s_proj2x", hist2->GetName()), binExclusionMax, binEtaMax);
      proj_c2->Add(proj2_c2, 1);
      proj_c2->Scale(1.0 / (binEtaMax - binExclusionMax + 1 + binExclusionMin - binEtaMin + 1));
    } else {
      //      proj_c1  = hist1->ProjectionX(Form("%s_proj1x", hist1->GetName()), binEtaMin, binEtaMax);
      //      proj_c2  = hist2->ProjectionX(Form("%s_proj1x", hist2->GetName()), binEtaMin, binEtaMax);
      Int_t binEtaMaxProj  = hist1->GetYaxis()->FindBin(-2.5 - 0.001);//reduced range to get rid of the NS jet
      proj_c1  = hist1->ProjectionX(Form("%s_proj1x", hist1->GetName()), binEtaMin, binEtaMaxProj);
      proj_c2  = hist2->ProjectionX(Form("%s_proj1x", hist2->GetName()), binEtaMin, binEtaMaxProj);
      proj_c1->Scale(1.0 / (binEtaMaxProj - binEtaMin + 1));
      proj_c2->Scale(1.0 / (binEtaMaxProj - binEtaMin + 1));
    }
    SymmetrizeTo0toPi(proj_c1);
    SymmetrizeTo0toPi(proj_c2);
    v2gaus2 = new TF1("v2gaus2","[0]+2*[1]*cos(2*x)+[2]*exp(-0.5*((x-TMath::Pi())/[3])**2)",-5,5);
    v2gaus2->SetParLimits(3,0.,1000);
    //    v2gaus2->SetParameter(0,baseLine2);//initial value of the baseline
    v2gaus2->SetParameter(0,proj_c2->GetMinimum());//initial value of the baseline
    v2gaus2->FixParameter(1,0);//no v2
    v2gaus2->SetParameter(2,0.1);//initial value of the normalization
    v2gaus2->SetParameter(3,proj_c2->GetRMS());//initial value of sigma
    proj_c2->Fit(v2gaus2,"MN","0",0,pi);
    v2gaus1 = new TF1("v2gaus1","[0]+2*[1]*cos(2*x)+[2]*exp(-0.5*((x-TMath::Pi())/[3])**2)",-5,5);
    v2gaus1->SetParLimits(3,0.,1000);
    v2gaus1->SetParameter(0,baseLine);//initial value of the baseline
    v2gaus1->SetParameter(1,0.1);//initial value of v2
    v2gaus1->SetParameter(2,0.1);//initial value of normalization
    v2gaus1->FixParameter(3,v2gaus2->GetParameter(3));//sigma fixed from low mult
    proj_c1->Fit(v2gaus1,"MN","0",0.00001,TMath::Pi()-0.00001);
    //now fix v2 and fit again
    if(0){//release the sigma and fix the v2 to see the evolution with multiplicity
      v2gaus1->FixParameter(1,v2gaus1->GetParameter(1));//fixed v2
      v2gaus1->ReleaseParameter(3);//release sigma
      proj_c1->Fit(v2gaus1,"MN","0",0.00001,TMath::Pi()-0.00001);
    }
    if(debug){
      TCanvas *cGausFit=new TCanvas("cGausFit","cGausFit");
      cGausFit->Divide(2,1);
      cGausFit->cd(1);
      proj_c1->GetXaxis()->SetRangeUser(0.001,pi-0.001);
      proj_c1->DrawCopy();
      v2gaus1->DrawCopy("same");
      TF1 *c_v2gaus1=(TF1*)v2gaus1->Clone();
      c_v2gaus1->SetLineColor(kBlue);
      c_v2gaus1->SetParameter(2,0);
      c_v2gaus1->DrawCopy("same");
      c_v2gaus1->SetLineStyle(2);
      c_v2gaus1->SetParameter(1,0);
      c_v2gaus1->DrawCopy("same");
      latex->DrawLatex(0.35,0.8,kSystemEnergy);
      latex->DrawLatex(0.35,0.73,textCentDif.Data());
      latex->DrawLatex(0.35,0.66,"High Multiplciity");
      cGausFit->cd(2);
      proj_c2->GetXaxis()->SetRangeUser(0.001,pi-0.001);
      proj_c2->DrawCopy();
      v2gaus2->DrawCopy("same");
      TF1 *c_v2gaus2=(TF1*)v2gaus2->Clone();
      c_v2gaus2->SetLineColor(kBlue);
      c_v2gaus2->SetParameter(2,0);
      c_v2gaus2->DrawCopy("same");
      latex->DrawLatex(0.35,0.8,textTrg.Data());
      latex->DrawLatex(0.35,0.73,textAss.Data());
      latex->DrawLatex(0.35,0.66,"Low Multiplciity");
      cGausFit->SaveAs("cGausFit.eps");
      cGausFit->SaveAs("cGausFit.pdf");
    }

  }

  if(gStudySystematic == k_baselineGausFit){
    //get projection with gap if in central barrel
    Printf("\033[1;31m Setting new baseLine2 from gaus fit: %f -> %f\033[m",baseLine2,v2gaus2->GetParameter(0));
    baseLine2 =v2gaus2->GetParameter(0);
    baseLineE2=v2gaus2->GetParError(0);
  }

  Float_t factor = 1.;
  Float_t Sigma1 = 0.,Sigma2 = 0.;
  Float_t factorErr = 0.;
  if(gStudySystematic == k_LowMultScale){
    //calculating the ratio of the yields in high/low multiplciity
    //areas of gaussian as a*|c|*sart(pi)
    Float_t p2_1=v2gaus1->GetParameter(2);
    Float_t p2_2=v2gaus2->GetParameter(2);
    Float_t ep2_1=v2gaus1->GetParError(2);
    Float_t ep2_2=v2gaus2->GetParError(2);
    Float_t p3_1=v2gaus1->GetParameter(3);
    Float_t p3_2=v2gaus2->GetParameter(3);
    Float_t ep3_1=v2gaus1->GetParError(3);
    Float_t ep3_2=v2gaus2->GetParError(3);

    Float_t yield1val = p2_1*p3_1;
    Float_t yield2val = p2_2*p3_2;
    Float_t yield1err = yield1val*TMath::Sqrt(ep2_1*ep2_1/p2_1/p2_1+ep3_1*ep3_1/p3_1/p3_1);
    Float_t yield2err = yield2val*TMath::Sqrt(ep2_2*ep2_2/p2_2/p2_2+ep3_2*ep3_2/p3_2/p3_2);

    factor = yield1val/yield2val;
    factorErr = factor*TMath::Sqrt(yield1err/yield1val*yield1err/yield1val + yield2err/yield2val*yield2err/yield2val);
    Sigma1=v2gaus1->GetParameter(3);
    Sigma2=v2gaus2->GetParameter(3);

    //    if (debug){
    //      new TCanvas;
    //      proj_c1->DrawCopy();
    //      v2gaus1->DrawCopy("same");
    //      v2gaus1->SetParameter(2,0);
    //      v2gaus1->SetLineColor(kBlue);
    //      v2gaus1->DrawCopy("same");
    //
    //      new TCanvas;
    //      proj_c2->DrawCopy();
    //      v2gaus2->DrawCopy("same");
    //      v2gaus2->SetParameter(2,0);
    //      v2gaus2->SetLineColor(kBlue);
    //      v2gaus2->DrawCopy("same");
    //      v2gaus2->SetParameter(1,0);
    //      v2gaus2->SetLineStyle(2);
    //      v2gaus2->DrawCopy("same");
    //    }
    Printf("\033[1;31m Scaling factors calculated on the fly\033[m");
    Printf("\033[1;31m Scaling Low Mult by %f  +- %f\033[m",factor,factorErr);
    hist2->Scale(factor);
    baseLine2*=factor;
    baseLineE2*=factor;
  }

  if(gStudySystematic!=k_no_subtraction)hist1->Add(hist2, -1);

  if (debug) {
    TCanvas* cSubt = new TCanvas("cSubt", "cSubt",1000,1000);
    PadFor2DCorr();
    TH2D* histSub = (TH2D*) hist1->Clone("histSub");
    histSub->GetYaxis()->SetRangeUser(etaMin+0.001,etaMax-0.001);
    histSub->DrawCopy("surf1");
    latex->DrawLatex(0.25,0.93,kSystemEnergy);
    latex->DrawLatex(0.25,0.86,textCentDif.Data());
    latex->DrawLatex(0.80,0.93,textTrg.Data());
    latex->DrawLatex(0.80,0.86,textAss.Data());
    gPad->Print(Form("subt%i%s",i,format));
    TFile* f = new TFile(inFileName.ReplaceAll("dphi/dphi_corr","debug"),"update");
    histSub->Write();
    f->Close();
    delete f;
  }

  TH1D* proj=0x0;
  if (!isForward){
    //standard projection

    //    TH1D* proj2=0x0;
    //    proj  = LinearFitProjectionX(hist1,Form("%s_proj1x", hist1->GetName()), binEtaMin, binExclusionMin,1);
    //    proj2 = LinearFitProjectionX(hist1,Form("%s_proj2x", hist1->GetName()), binExclusionMax, binEtaMax,1);
    //    proj->Add(proj2, 1);
    //    proj->Scale(1./(binEtaMax - binExclusionMax + 1 + binExclusionMin - binEtaMin + 1));

    //        projection symmetrized in eta

    proj  = LinearFitProjectionX(hist1,Form("%s_proj1x", hist1->GetName()), binEtaMin, binExclusionMin,0,1);
    proj->Scale(1./(binExclusionMin - binEtaMin + 1));

  } else {
    proj = LinearFitProjectionX(hist1,Form("%s_proj1x", hist1->GetName()), binEtaMin, binEtaMax);
    proj->Scale(1./(binEtaMax - binEtaMin + 1));
  }
  proj->SetStats(0);
  proj->GetYaxis()->SetNdivisions(505);
  proj->GetXaxis()->SetTitle("#Delta#varphi (rad)");
  proj->GetYaxis()->SetTitle(kProjYieldTitlePhi);
  proj->GetXaxis()->SetTitleOffset(1.1);
  proj->GetYaxis()->SetTitleOffset(1.1);
  proj->GetYaxis()->SetLabelSize(fontSize);
  proj->GetXaxis()->SetLabelSize(fontSize);
  proj->GetXaxis()->SetTitleSize(fontSize);
  proj->GetYaxis()->SetTitleSize(fontSize);
  proj->SetTitle("");
  proj->SetMarkerStyle(21);
  proj->SetMarkerSize(0.7);
  Float_t range = proj->GetMaximum()-proj->GetMinimum();
  proj->SetMaximum(proj->GetMinimum()+range*2);
  proj->SetMinimum(proj->GetMinimum()-range*0.1);

  TF1* vn = new TF1("vn","[0]+2*[1]*cos(1*x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)",-5,5);
  TF1* v1 = new TF1("v1","[0]+2*[1]*cos(1*x)",-5,5);
  TF1* v2 = new TF1("v2","[0]+2*[1]*cos(2*x)",-5,5);
  TF1* v3 = new TF1("v3","[0]+2*[1]*cos(3*x)",-5,5);
  TF1* bl = new TF1("bl","[0]",-5,5);
  vn->SetLineStyle(kSolid);       vn->SetLineColor(kBlack);
  bl->SetLineStyle(kSolid);       bl->SetLineColor(kBlue);
  v1->SetLineStyle(kDotted);      v1->SetLineColor(kBlue);
  v2->SetLineStyle(kDashed);      v2->SetLineColor(kRed);
  v3->SetLineStyle(kDashDotted);  v3->SetLineColor(kMagenta);
  if (gStudySystematic==k_no_v3) vn->FixParameter(3,0);

  proj->Fit(vn,gStudySystematic==k_fit_opt ? "0Q" : fitOption,"");
  Double_t* pval = vn->GetParameters();
  Double_t* perr = vn->GetParErrors();
  Float_t min = vn->GetMinimum();

  Float_t v0val = baseLine2+pval[0];
  Float_t v0err = TMath::Sqrt(baseLineE2*baseLineE2+perr[0]*perr[0]);

  TF1 *parabola=0x0;
  if (gStudySystematic == k_baselineHighMult || gStudySystematic == k_baselineHighMultParab) {
    // alternative way for baseline
    Printf("\033[1;31m Using the minum to calculate the v2: %f\033[m",min);
    if(gStudySystematic == k_baselineHighMultParab){
      parabola = new TF1("parabola", "[0] + [1]*(x - [2])**2",0,pi);
      parabola->SetLineColor(3);
      parabola->SetParameters(min, -0.1,pi/2);
      proj->Fit(parabola, fitOption, "", pi/2-1,pi/2 + 1);
      min = parabola->GetParameter(0);
      Printf("\033[1;31m Setting new minimum from parabolic fit: %f\033[m",min);
      delete parabola;
    }
    Float_t diffMinParam0 = pval[0] - min;
    v0val = baseLine + diffMinParam0;
    v0err = baseLineE;
  }

  Float_t chi2 = vn->GetChisquare();
  Int_t ndf  = vn->GetNDF();
  Printf("Chi2: %f ndf: %d chi2/ndf: %f Min: %f",chi2,ndf,chi2/ndf, min);

  bl->SetParameter(0,min);
  v1->SetParameters(pval[0],pval[1]);
  v2->SetParameters(pval[0],pval[2]);
  v3->SetParameters(pval[0],pval[3]);

  Float_t v1val = pval[1]/v0val;
  Float_t v2val = pval[2]/v0val;
  Float_t v3val = pval[3]/v0val;
  Float_t v1err = v1val*sqrt(pow(v0err/v0val,2.)+pow(perr[1]/pval[1],2.));
  Float_t v2err = v2val*sqrt(pow(v0err/v0val,2.)+pow(perr[2]/pval[2],2.));
  Float_t v3err = v3val*sqrt(pow(v0err/v0val,2.)+pow(perr[3]/pval[3],2.));

  if(gStudySystematic==k_v2_analytical){
    Printf("\033[1;31m v2 calculated analytically: %f\033[m",min);
    v2val = CalculateVnAnalytically(proj,2,baseLine2);//baseline2 is the peripheral one
    v2err  = CalculateVnErrAnalytically(proj,2,baseLine2,baseLineE2);
  }

  Printf("Baseline: %f +- %f; v1 = %f +- %f; v2 = %f +- %f; v3 = %f +- %f", baseLine, baseLineE, v1val, v1err, v2val, v2err, v3val, v3err);

  if (debug) {
    TCanvas* c = new TCanvas("c2", "c2", 1200, 800);
    gPad->SetMargin(0.12,0.01,0.12,0.01);

    proj->DrawCopy("E0 X0");
    vn->DrawCopy("same");
    bl->DrawCopy("same");
    v1->DrawCopy("same");
    v2->DrawCopy("same");
    v3->DrawCopy("same");

    TLegend* legend = new TLegend(0.48, 0.68, 0.98, 0.98);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(fontSize);

    legend->AddEntry(proj, "Data", "P");
    legend->AddEntry(vn, "a_{0} + #Sigma_{n=1}^{3} a_{n} cos(n#Delta#varphi)", "L");
    legend->AddEntry(v1, Form("a_{0} + a_{1} cos(#Delta#varphi), V_{1}#times10^{3}=%.2f#pm%.2f",v1val*1000,TMath::Abs(v1err*1000)), "L");
    legend->AddEntry(v2, Form("a_{0} + a_{2} cos(2#Delta#varphi), V_{2}#times10^{3}=%.2f#pm%.2f",v2val*1000,TMath::Abs(v2err*1000)), "L");
    legend->AddEntry(v3, Form("a_{0} + a_{3} cos(3#Delta#varphi), V_{3}#times10^{3}=%.2f#pm%.2f",v3val*1000,TMath::Abs(v3err*1000)), "L");
    legend->AddEntry(bl, "Baseline", "L");
    legend->Draw();

    latex->DrawLatex(0.3,0.93,kSystemEnergy);
    latex->DrawLatex(0.3,0.86,textCentDif.Data());
    latex->DrawLatex(0.3,0.79,textTrg.Data());
    latex->DrawLatex(0.3,0.72,textAss.Data());
    latex->DrawLatex(0.3,0.65,Form("%.1f < #Delta#eta < %.1f",etaMin,etaMax));
    if (exclusionMin>etaMin) latex->DrawLatex(0.3,0.58,Form("Ex: %.1f < #Delta#eta < %.1f",actualExclusionMin,actualExclusionMax));

    gPad->Print(Form("proj%i%s",i,format));
    TFile* f = new TFile(inFileName.ReplaceAll("dphi/dphi_corr","debug"),"update");
    proj->Write();
    vn->Write();
    bl->Write();
    v1->Write();
    v2->Write();
    v3->Write();
    f->Close();
    delete f;

  }

  if (res) {
    res[0] = baseLine;
    res[1] = v1val;
    res[2] = v2val;
    res[3] = v3val;
    res[4] = baseLineE;
    res[5] = v1err;
    res[6] = v2err;
    res[7] = v3err;
    res[8] = chi2;
    res[9] = ndf;
    res[10] = factor;
    res[11] = factorErr;
    res[12] = Sigma1;
    res[13] = Sigma2;
  }
  infile->Close();
  delete infile;
  delete vn;
  delete v1;
  delete v2;
  delete v3;
  delete bl;
}



Float_t CalculateVnAnalyticallyRaw(TH1* proj,Int_t n, Float_t BaseLineShift){
  //proj->Print("all");
  Float_t num=0;
  Float_t den=0;
  Printf("BaseLineShift %f",BaseLineShift);
  for(Int_t ibin=1;ibin<=proj->GetNbinsX();ibin++){
    if(proj->GetBinContent(ibin)==0)continue;
    num+=TMath::Cos(n*proj->GetBinCenter(ibin))*(proj->GetBinContent(ibin)+BaseLineShift);
    den+=proj->GetBinContent(ibin)+BaseLineShift;
    //Printf("num %f   den %f BaselineShift %f",num,den,BaseLineShift);
  }
  Float_t vn=num/den;
  return vn;
}

Float_t CalculateVnAnalytically(TH1* proj,Int_t n, Float_t BaseLineShift){
  Float_t vn=CalculateVnAnalyticallyRaw(proj, n, BaseLineShift);
  if (vn < 0)
    return 0.;
  Printf("v2 Value: %f",vn);
  return vn;
}

Float_t CalculateVnErrAnalytically(TH1* proj,Int_t n, Float_t BaseLineShift, Float_t BaseLineShiftError){
  //TODO cross check error calculation
  Float_t I=0.;
  Float_t v2=CalculateVnAnalytically(proj,n,BaseLineShift);
  if (v2 <= 0.)
    return 0.;

  for(Int_t ibin=1;ibin<=proj->GetNbinsX();ibin++){
    if(proj->GetBinContent(ibin)==0)continue;//skip empty bins (if any)
    I+=proj->GetBinContent(ibin)+BaseLineShift;//calculate the integral
  }
  Float_t sigma2=0.;
  for(Int_t ibin=1;ibin<=proj->GetNbinsX();ibin++){
    if(proj->GetBinContent(ibin)==0)continue;
    //error calculation assuming v2=sum_i wi*pi / sum_j pj
    Float_t wi  =TMath::Cos(n*proj->GetBinCenter(ibin));
    Float_t pi  = proj->GetBinContent(ibin)+BaseLineShift;
    Float_t epi  =proj->GetBinError(ibin)+BaseLineShiftError;
    sigma2+=TMath::Power( (wi - v2)/I * epi ,2);
  }
  Float_t vnerr=TMath::Sqrt(sigma2);
  return vnerr;
}

void SymmetrizeTo0toPi(TH1D* h){
  for(Int_t k=1; k<=h->GetNbinsX(); k++){//symmetrize in 0-Pi
    Int_t ibin = -999.;
    Double_t phi = h->GetBinCenter(k); 
    if      (phi<0)  ibin = h->FindBin(TMath::Abs(phi));
    else if (phi>pi) ibin = h->FindBin(2*pi-phi);
    if (ibin<0) continue;//here we are within 0 an Pi
    h->SetBinContent(ibin,(h->GetBinContent(ibin)+h->GetBinContent(k))/2);
    h->SetBinError(ibin,TMath::Sqrt(h->GetBinError(ibin)*h->GetBinError(ibin)+h->GetBinError(k)*h->GetBinError(k))/2);
    h->SetBinContent(k,0.);
    h->SetBinError(k,0.);
  }
}

void SymmetrizeToEta0(TH1D* h){
  for(Int_t k=1; k<=h->GetNbinsX(); k++){//symmetrize in eta
    Int_t ibin = -999.;
    Double_t eta = h->GetBinCenter(k);
    if(eta<0.){
      ibin = h->FindBin(TMath::Abs(eta));
      h->SetBinContent(k,(h->GetBinContent(ibin)+h->GetBinContent(k))/2);
      h->SetBinError(k,TMath::Sqrt(h->GetBinError(ibin)*h->GetBinError(ibin)+h->GetBinError(k)*h->GetBinError(k))/2);
    }else{
      h->SetBinContent(k,0.);
      h->SetBinError(k,0.);
    }
  }
}

TH1D* LinearFitProjectionX(TH2F* h2, TString name, const Int_t bin1, const Int_t bin2, Bool_t UseStandardProj=kFALSE ,Bool_t symmetrizeEta=kFALSE ,Int_t debug=0){
  Float_t ymin = h2->GetYaxis()->GetBinLowEdge(bin1);
  Float_t ymax = h2->GetYaxis()->GetBinUpEdge(bin2);
  printf("Fit in %f %f\n",ymin,ymax);
  TH1D* h = h2->ProjectionX(name,bin1,bin2);
  TCanvas* cdebug = 0;
  TLatex* latex = 0;
  if (debug>1) {
    cdebug = new TCanvas("linproj","linproj",1900,1000);
    cdebug->Divide(h2->GetNbinsX()/3,3,0.001,0.001);
    latex = new TLatex();
    latex->SetTextSize(fontSize);
    latex->SetTextFont(42);
    latex->SetTextAlign(11);
    latex->SetNDC();
  }
  for (Int_t ix=1;ix<=h2->GetNbinsX();ix++) {
    TH1D* hbin = h2->ProjectionY(Form("hbin%i",ix),ix,ix);
    if(symmetrizeEta)SymmetrizeToEta0(hbin);
    hbin->SetLineColor(kBlue);
    hbin->SetLineWidth(2);
    //    TF1* pol = new TF1("pol",gStudySystematic == k_parabolic_fit ? "pol2":"pol1" ,ymin,ymax);
    TF1* pol=0x0;
    if(gStudySystematic == k_parabolic_fit)pol = new TF1("pol","pol2",ymin,ymax);
    else if(gStudySystematic == k_constant_fit)pol = new TF1("pol","pol0",ymin,ymax);
    else pol = new TF1("pol","pol1",ymin,ymax);
    hbin->Fit(pol,"Q0","",ymin,ymax);
    Double_t value = pol->Integral(ymin,ymax)/hbin->GetXaxis()->GetBinWidth(1);
    Double_t error = pol->IntegralError(ymin,ymax)/hbin->GetXaxis()->GetBinWidth(1);
    if (gStudySystematic == k_hist_integral || UseStandardProj) value = hbin->IntegralAndError(hbin->FindFixBin(ymin+0.001),hbin->FindFixBin(ymax-0.001),error);
    h->SetBinContent(ix,value);
    h->SetBinError(ix,error);
    if (cdebug) {
      cdebug->cd(ix);
      gPad->SetRightMargin(0.03);
      gPad->SetTopMargin(0.02);
      gPad->SetBottomMargin(0.11);
      hbin->SetMaximum(h2->GetMaximum()*1.1);
      hbin->SetMinimum(h2->GetMinimum()*0.98);
      hbin->DrawCopy();
      latex->DrawLatex(0.2,0.90,Form("#varphi bin = %i",ix));
      latex->DrawLatex(0.2,0.84,Form("#chi^{2}/ndf = %.2f/%i",pol->GetChisquare(),pol->GetNDF()));
      pol->SetLineWidth(1);
      pol->DrawCopy("same");
      value/=(bin2-bin1+1);
      error/=(bin2-bin1+1);
      latex->DrawLatex(0.2,0.78,Form("fit integral = %.4f #pm %.4f",value,error));
      value = hbin->IntegralAndError(hbin->FindFixBin(ymin+0.001),hbin->FindFixBin(ymax-0.001),error);
      value/=(bin2-bin1+1);
      error/=(bin2-bin1+1);
      latex->DrawLatex(0.2,0.72,Form("hist integral = %.4f #pm %.4f",value,error));
    }
    delete pol;
    delete hbin;
  }
  if (cdebug) cdebug->Print("proj_debug.png");
  if (cdebug) cdebug->Print(Form("proj_debug_%s.eps",name.Data()));
  return h;
}
