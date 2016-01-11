//////////////////////////////////////////////////////////////
//                                                          //
// This macro produces all final results and QA plots for   //
// the DEUTERON nuclei analysis.                            //
//                                                          //
//                                                          //
//                                                          //
//////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFractionFitter.h"
#include "TColor.h"



using namespace std;

//
// global variables
//

Char_t * fileNameMcEnhancedNew = "MC/McEnhancedNew.root";
Char_t * fileNameMcEnhancedOld = "MC/McEnhancedOld2.root";
Char_t * fileNameMcUnEnhanced  = "MC/McUnEnhanced.root";
Char_t * fileNameData          = "data/dataFinal.root";

//
// definition of functions
//
void       MakeEfficiency();
TH1D    *  ExtractEfficiency(TList * listMC, Int_t sign = -1, Int_t centrality = 0, Bool_t isTOF = kTRUE);
TF1     *  MakeEtaRapidityGenCorrection(Bool_t drawQA = kFALSE);
//
void       MakeFinalSpectra();
TH1D    *  MakeRawSpectraTPC(TList * list, Int_t sign = -1, Int_t centralityBin = 0);
TH1D    *  MakeRawSpectra(TList * list, Int_t sign = -1, Int_t centralityBin = 0);
void       NormalizeSpectrum(TH1D * spectrum, Float_t dy, Float_t numberOfEvents);
void       PlotSpectraComparison();
void       CombinationAndSystematics();
TH1D    *  GetSpectrumWithSyst(TH1D * histStatError);
//
void       MakeQAplots();
TCanvas *  PlotTpcQA(Char_t * fileName = "", Bool_t isMC = kTRUE,  Bool_t isMCold = kFALSE);
void       PlotSpectraQA(TList * list, Int_t particleType = 0, Int_t sign = -1);
//
void       MakeMaterialCorrection();
Float_t    GetMaterialCorrection(Float_t ptBin, Int_t centralityBin, Bool_t pureMC = kFALSE);



//
// implementation
//
//_______________________________________________________________________
void CombinationAndSystematics() {
  //
  // produce the final spectra for the paper
  //
  Int_t kMaxCentrality = 5;
  //
  TFile fileOut("output/spectraDeuteronCombined.root", "RECREATE");
  fileOut.Close();
  //
  TFile * fileIn = TFile::Open("output/spectraDeuteron.root");
  //
  TCanvas * canvDeuteron = new TCanvas("canvDeuteron","deuteron spectra");
  //
  for(Int_t iCentr = 0; iCentr < kMaxCentrality; iCentr++) {
    //
    TH1D * tpc = (TH1D *) fileIn->Get(Form("cent%i_sign%s_TPC",iCentr,"Pos"));
    TH1D * tof = (TH1D *) fileIn->Get(Form("cent%i_sign%s",iCentr,"Pos"));
    //
    TH1D * deuteron = (TH1D *) tof->Clone();
    deuteron->SetNameTitle(Form("deuteron_centrality%i",iCentr), Form("deuteron_centrality%i",iCentr));
    deuteron->GetYaxis()->SetTitle("d#it{N}/(d#it{y}d#it{p_{T}})");
    deuteron->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
    //
    // merge the two
    //
    deuteron->SetBinContent(4, tpc->GetBinContent(4));  deuteron->SetBinError(4, tpc->GetBinError(4));
    deuteron->SetBinContent(5, tpc->GetBinContent(5));  deuteron->SetBinError(5, tpc->GetBinError(4));
    //
    // get systematics
    //
    TH1D * deuteronSyst = GetSpectrumWithSyst(deuteron);
    //
    // drawing and saving
    //
    if (iCentr == 0) {
      deuteronSyst->DrawCopy("E2");
      deuteron->DrawCopy("SAME");
    } else {
      deuteronSyst->DrawCopy("E2,SAME");
      deuteron->DrawCopy("SAME");
    }
    //

    //
    TFile fileIn2("output/fitsBW.root");
    TF1 * fit = (TF1 * ) fileIn2.Get(Form("fit%i_signNeg", iCentr));
    if (fit) {
      fit->Draw("SAME");
    }
    fileIn2.Close();
    //
    TFile fileOut("output/spectraDeuteronCombined.root", "UPDATE");
    deuteron->Write();
    deuteronSyst->Write();
    if (fit) fit->Write();
    fileOut.Close();

  }



}

//_______________________________________________________________________
TH1D * GetSpectrumWithSyst(TH1D * histStatError) {
  //
  // add systematic error to the points
  //
  TH1D * histSyst = (TH1D *) histStatError->Clone();
  histSyst->SetNameTitle(Form("%s_SYST", histStatError->GetName()), 
			 Form("%s_SYST", histStatError->GetName()));
  histSyst->SetFillColor(histStatError->GetMarkerColor());
  histSyst->SetFillStyle(0);
  histSyst->SetDrawOption("E2");
  //
  for(Int_t i=0; i <  histSyst->GetXaxis()->GetNbins(); i++) { // begin loop over pt-bins
      Float_t syst = 0.;
      Float_t pt = histSyst->GetXaxis()->GetBinCenter(i);
      //
      // (1.) tracking and matching normally 4%, but we add 2% for the pt-correction
      //
      syst += 0.06*0.06;
      //
      // (2.) PID -- difference between bin counting and fit in TOF
      //
      if (pt > 1.4) syst += 0.05*0.05;
      //
      // (3.) material knock-out -- 20% of correction
      //
      Float_t corr = 1.5*TMath::Exp(1.27259 - 2.527*pt);
      syst += (0.2*corr)*(0.2*corr);
      //
      // (4.) absorption in material (shifted proton exponential, difference between Eulogio's and geant production)
      // ==> for anti-deuterons correspondingly more. ==> see plot by Natasha (backup slide APW)
      //
      syst += 0.03*0.03;
      //
      // (5.) matching efficiency
      //
      //if (pt > 1.) syst += 0.05*0.05;
      //
      // put it to the spectrum
      //
      Float_t  sum = TMath::Sqrt(syst)*histSyst->GetBinContent(i);
      histSyst->SetBinError(i,sum);
    } // end loop over pt-bins

  //
  return histSyst;


}


//_______________________________________________________________________
void PlotSpectraComparison() {
  //
  // plot TPC & TOF spectra for comparison
  //
  Int_t kMaxCentrality = 5;
  //
  //
  TFile * fileIn = TFile::Open("output/spectraDeuteron.root");
  //
  TCanvas * canvRatio = new TCanvas("canvRatio","canvRatio");
  //
  for(Int_t iCentr = 0; iCentr < kMaxCentrality; iCentr++) {
    //
    TCanvas * canvComparison = new TCanvas(Form("canvComparison_%i",iCentr),Form("spectra in centrality %i", iCentr));
    //
    TH1D * tpc = (TH1D *) fileIn->Get(Form("cent%i_sign%s_TPC",iCentr,"Pos"));
    TH1D * tof = (TH1D *) fileIn->Get(Form("cent%i_sign%s",iCentr,"Pos"));
    //
    TH1D * tpcNeg = (TH1D *) fileIn->Get(Form("cent%i_sign%s_TPC",iCentr,"Neg"));
    TH1D * tofNeg = (TH1D *) fileIn->Get(Form("cent%i_sign%s",iCentr,"Neg"));
    tpcNeg->SetMarkerStyle(21);
    tofNeg->SetMarkerStyle(25);
    //
    tof->DrawCopy();
    tpc->DrawCopy("SAME");
    tpcNeg->DrawCopy("SAME");
    tofNeg->DrawCopy("SAME");
    //
    TFile fileIn2("output/fitsBW.root");
    TF1 * fit = (TF1 * ) fileIn2.Get(Form("fit%i_signNeg", iCentr));
    fit->Draw("SAME");
    fileIn2.Close();
    //
    //TCanvas * canvRatio = new TCanvas(Form("canvRatio_%i",iCentr),Form("Ratio in centrality %i", iCentr));
    canvRatio->cd();
    tofNeg->Divide(tof);
    if (iCentr == 0) {
      tofNeg->DrawCopy();
    } else {
      tofNeg->DrawCopy("SAME");
    }
  }


}


//_______________________________________________________________________
void MakeFinalSpectra() {
  //
  // make final spectra for all centrality bins
  //
  Int_t kMaxCentrality = 5;
  TFile fileOut("output/spectraDeuteron.root", "RECREATE");
  fileOut.Close();
  //
  // get efficiencies and data
  //
  TFile * inFileEff = TFile::Open("output/efficiencies.root");
  TH1D  * efficiencyTrackingNegNew = (TH1D*) inFileEff->Get("efficiencyTrackingNegNew");
  TH1D  * efficiencyTrackingPosNew = (TH1D*) inFileEff->Get("efficiencyTrackingPosNew");
  TH1D  * efficiencyTofPosNew      = (TH1D*) inFileEff->Get("efficiencyTofPosNew");
  TH1D  * efficiencyTofNegNew      = (TH1D*) inFileEff->Get("efficiencyTofNegNew");
  //
  TFile inFileData(fileNameData);
  TList * list = (TList * ) inFileData.Get("akalweit_Nuclei");
  //
  // positive particles
  //
  TCanvas * canvAllPos = new TCanvas("canvAllPos","deuteron spectra");
  for(Int_t iCentr = 0; iCentr < kMaxCentrality; iCentr++) {
    //
    TH1D * rawSpecTpc = MakeRawSpectraTPC(list, +1, iCentr);
    TH1D * rawSpecTof = MakeRawSpectra(list, +1, iCentr);
    rawSpecTpc->Divide(efficiencyTrackingPosNew);
    rawSpecTof->Divide(efficiencyTofPosNew);
    //
    canvAllPos->cd();
    if (iCentr == 0) {
      rawSpecTpc->DrawCopy("EP");
      rawSpecTof->DrawCopy("EPSAME");
    } else {
      rawSpecTpc->DrawCopy("EPSAME");
      rawSpecTof->DrawCopy("EPSAME");
    }
    //
    // apply material correction
    //
    TFile matFile("output/materialCorrection.root");
    TH1D * materialCorrection = (TH1D *) matFile.Get(Form("matCorr_%i",iCentr));
    for(Int_t iBin = 0; iBin < rawSpecTpc->GetXaxis()->GetNbins(); iBin++) {
      Double_t pT = rawSpecTpc->GetXaxis()->GetBinCenter(iBin);
      if (pT > 0.4 && pT < 2.) {
	Float_t primFactor = 1. - materialCorrection->GetBinContent(iBin);
	rawSpecTof->SetBinContent(iBin, rawSpecTof->GetBinContent(iBin)*primFactor);
	rawSpecTof->SetBinError(iBin, rawSpecTof->GetBinError(iBin)*primFactor);
	rawSpecTpc->SetBinContent(iBin, rawSpecTpc->GetBinContent(iBin)*primFactor);
	rawSpecTpc->SetBinError(iBin, rawSpecTpc->GetBinError(iBin)*primFactor);
      }
    }
    //
    matFile.Close();
    //
    TFile fileOutput("output/spectraDeuteron.root", "UPDATE");
    rawSpecTpc->Write();
    rawSpecTof->Write();
    fileOutput.Close();
  }
  canvAllPos->Print("QAplots/posSpectra.pdf");
  //
  // negative particles
  //
  TCanvas * canvAllNeg = new TCanvas("canvAllNeg","anti-deuteron spectra");
  for(Int_t iCentr = 0; iCentr < kMaxCentrality; iCentr++) {
    //
    TH1D * rawSpecTpc = MakeRawSpectraTPC(list, -1, iCentr);
    TH1D * rawSpecTof = MakeRawSpectra(list, -1, iCentr);
    rawSpecTpc->Divide(efficiencyTrackingNegNew);
    rawSpecTof->Divide(efficiencyTofNegNew);
    //
    canvAllPos->cd();
    if (iCentr == 0) {
      rawSpecTpc->DrawCopy("EP");
      rawSpecTof->DrawCopy("EPSAME");
    } else {
      rawSpecTpc->DrawCopy("EPSAME");
      rawSpecTof->DrawCopy("EPSAME");
    }
    TFile fileOutput("output/spectraDeuteron.root", "UPDATE");
    rawSpecTpc->Write();
    rawSpecTof->Write();
    fileOutput.Close();
  }
  canvAllNeg->Print("QAplots/negSpectra.pdf");
  


}


//_______________________________________________________________________
TH1D * MakeRawSpectraTPC(TList * list, Int_t sign, Int_t centralityBin) {
  //
  // Make raw spectra for given sign and centrality bin
  //
  // (0.) assumed particle: 0. deuteron, 1. triton, 2. He-3
  // (1.) multiplicity or centrality -- number of accepted ESD tracks per events (deprecated), but now classes from 1 to 10, 0: Min. Bias
  // (2.) pT
  // (3.) sign
  // (4.) rapidity --> filled 4xW
  // (5.)  pull TPC dEx --> filled 4x
  // (6.) has valid TOF pid signal
  // (7.) nsigma TOF --> filled 4x XXXXXXXXX no mass*mass
  // (8..) dca_xy
  // (9.) CODE -- only MC 0-generated, 1-true rec. primaries, 2-misident prim, 3-second weak, 4-second material, 5-misident sec
  //
  TH1D * histCentr = (TH1D *)  list->FindObject("fHistCentrality");
  //
  THnSparseF * hist = (THnSparseF *) list->FindObject("fHistRealTracks");
  hist->GetAxis(0)->SetRangeUser(0,0);                     // select deuterons
  //
  // centrality 0: 0-5% / 1: 5-10% / 2: 10-20% / 3: 20-30% / 4: 30-40% / 5: 40-50% / 6: 50-60% / 7: 60-70% / 8: 70-80%
  //
  Int_t centrality = 0;
  if (centralityBin == 0) centrality = 0; //  0-10%
  if (centralityBin == 1) centrality = 2; // 10-20%
  if (centralityBin == 2) centrality = 3; // 20-40%
  if (centralityBin == 3) centrality = 5; // 40-60%
  if (centralityBin == 4) centrality = 7; // 60-80%
  //
  hist->GetAxis(1)->SetRangeUser(centrality,centrality);   // select centrality bin
  Float_t norm = histCentr->GetBinContent(centrality+2);
  if (centralityBin != 1) {
    hist->GetAxis(1)->SetRangeUser(centrality,centrality+1);   // select centrality bin
    norm += histCentr->GetBinContent(centrality+3);
  }
  //
  // apply cuts -- Has to be consistent with efficiency
  //
  hist->GetAxis(4)->SetRangeUser(-0.49, 0.49);             // select rapidity range
  hist->GetAxis(3)->SetRangeUser(sign,sign);               // select anti-deuterons
  hist->GetAxis(8)->SetRangeUser(-0.5,0.5);  // DCA-cut / TODO: do proper unfolding
  hist->GetAxis(5)->SetRangeUser(-2.5,2.5);  // TPC-PID cut  
  //
  TH1D * spec    = hist->Projection(2);
  spec->SetNameTitle("spec","SPEC");
  TH1D * rawSpec = hist->Projection(2);
  rawSpec->Reset();
  for(Int_t iBin = 0; iBin < spec->GetXaxis()->GetNbins(); iBin++) {
    Double_t pT = spec->GetXaxis()->GetBinCenter(iBin);
    if (pT > 0.6 && pT < 1.2) {
      rawSpec->SetBinContent(iBin, spec->GetBinContent(iBin));
      rawSpec->SetBinError(iBin, spec->GetBinError(iBin));
    }
  }
  //
  //
  NormalizeSpectrum(rawSpec,1,norm);
  //
  if (sign < 0) rawSpec->SetNameTitle(Form("cent%i_sign%s_TPC",centralityBin,"Neg"),
				      Form("cent%i_sign%s_TPC",centralityBin,"Neg"));
  if (sign > 0) rawSpec->SetNameTitle(Form("cent%i_sign%s_TPC",centralityBin,"Pos"),
				      Form("cent%i_sign%s_TPC",centralityBin,"Pos"));
  //
  Int_t colorList[9] = {600+4, 600+3, 600+2, 600+1, 600, 600-4, 600-7, 600-9, 600-10};
  rawSpec->SetMarkerStyle(20);
  rawSpec->SetMarkerSize(1.4);
  rawSpec->SetMarkerColor(colorList[centrality]);
  rawSpec->SetLineColor(colorList[centrality]);
  //
  hist->GetAxis(5)->SetRangeUser(-1000.,1000);  // TPC-PID cut  
  hist->GetAxis(6)->SetRangeUser(-1000.,1000);  // TOF-PID cut  
  hist->GetAxis(7)->SetRangeUser(-1000.,1000);  // TOF-PID cut    
  hist->GetAxis(8)->SetRangeUser(-1000.,1000);  // TOF-PID cut    

  //
  return rawSpec;

}



//_______________________________________________________________________
TH1D * MakeRawSpectra(TList * list, Int_t sign, Int_t centralityBin) {
  //
  // Make raw spectra for given sign and centrality bin
  //
  //
  // (0.) assumed particle: 0. deuteron, 1. triton, 2. He-3
  // (1.) multiplicity or centrality -- number of accepted ESD tracks per events (deprecated), but now classes from 1 to 10, 0: Min. Bias
  // (2.) pT
  // (3.) sign
  // (4.) rapidity --> filled 4xW
  // (5.)  pull TPC dEx --> filled 4x
  // (6.) has valid TOF pid signal
  // (7.) nsigma TOF --> filled 4x XXXXXXXXX no mass*mass
  // (8..) dca_xy
  // (9.) CODE -- only MC 0-generated, 1-true rec. primaries, 2-misident prim, 3-second weak, 4-second material, 5-misident sec
  //
  TH1D * histCentr = (TH1D *)  list->FindObject("fHistCentrality");
  //
  THnSparseF * hist = (THnSparseF *) list->FindObject("fHistRealTracks");
  hist->GetAxis(0)->SetRangeUser(0,0);                     // select deuterons
  //
  // centrality 0: 0-5% / 1: 5-10% / 2: 10-20% / 3: 20-30% / 4: 30-40% / 5: 40-50% / 6: 50-60% / 7: 60-70% / 8: 70-80%
  //
  Int_t centrality = 0;
  if (centralityBin == 0) centrality = 0; //  0-10%
  if (centralityBin == 1) centrality = 2; // 10-20%
  if (centralityBin == 2) centrality = 3; // 20-40%
  if (centralityBin == 3) centrality = 5; // 40-60%
  if (centralityBin == 4) centrality = 7; // 60-80%
  //
  hist->GetAxis(1)->SetRangeUser(centrality,centrality);   // select centrality bin
  Float_t norm = histCentr->GetBinContent(centrality+2);
  if (centralityBin != 1) {
    hist->GetAxis(1)->SetRangeUser(centrality,centrality+1);   // select centrality bin
    norm += histCentr->GetBinContent(centrality+3);
  }
  //
  //
  //
  hist->GetAxis(4)->SetRangeUser(-0.49, 0.49);             // select rapidity range
  hist->GetAxis(3)->SetRangeUser(sign,sign);               // select anti-deuterons
  //
  hist->GetAxis(8)->SetRangeUser(-0.5,0.5);  // DCA-cut / TODO: do proper unfolding
  //
  hist->GetAxis(5)->SetRangeUser(-2.5,2.5);  // TPC-PID cut  
  //hist->GetAxis(7)->SetRangeUser(-2.,2.5);
  //
  // PID in the TOF
  //
  hist->GetAxis(6)->SetRangeUser(1,1);        // TOF-PID cut: require hasTOF
  //
  TH1D * rawSpec = hist->Projection(2);
  rawSpec->Reset();
  if (sign < 0) rawSpec->SetNameTitle(Form("cent%i_sign%s",centralityBin,"Neg"),
				      Form("cent%i_sign%s",centralityBin,"Neg"));
  if (sign > 0) rawSpec->SetNameTitle(Form("cent%i_sign%s",centralityBin,"Pos"),
				      Form("cent%i_sign%s",centralityBin,"Pos"));
  //
  // --> simple PID: hist->GetAxis(7)->SetRangeUser(-0.8.,0.8);  // TOF-PID cut
  //TF1 fitFunc("fitFunc","gaus(0) + [3] + [4]*x",-1.5,1.5);
  TF1 fTOFsignal("fTOFsignal", "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * TMath::Exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))    +    [4] + [5]*x   +   [6]*TMath::Exp(-[7]*x) ", -2.2,2.2); //0:norm, 1:mean, 2:sigma, 3:tail
  //
  //
  TH2D * m2fit = (TH2D *) hist->Projection(7,2);
  m2fit->SetNameTitle(Form("m2fits_%i_%i",sign,centralityBin),
		      Form("m**2 fits for sign %i and centrality %i",sign,centralityBin));
  hist->GetAxis(5)->SetRangeUser(-1000.,1000);  // TPC-PID cut  
  hist->GetAxis(6)->SetRangeUser(-1000.,1000);  // TOF-PID cut  
  hist->GetAxis(7)->SetRangeUser(-1000.,1000);  // TOF-PID cut    
  hist->GetAxis(8)->SetRangeUser(-1000.,1000);  // TOF-PID cut    
  //
  // rebin because of statistics
  //
  m2fit->RebinY(2);
  //
  TCanvas * canvM2fits = new TCanvas(Form("canvM2fits_%i_%i",sign,centralityBin),
				     Form("m**2 fits for sign %i and centrality %i",sign,centralityBin));
  canvM2fits->Print(Form("QA/M2fits_%i_%i.pdf(",sign,centralityBin));
  //
  Double_t oldParStartValue[8] = {664.78 , 0.0178744 , 0.212441 , 0.229438 , 171.889, 0. , 0.0203292};
  //
  for(Int_t iBin = 1; iBin < m2fit->GetXaxis()->GetNbins(); iBin++) {
    //
    Float_t ptBin = m2fit->GetXaxis()->GetBinCenter(iBin);
    if (ptBin < 0.6) continue;
    //
    TH1D * m2distr = m2fit->ProjectionY(Form("m2_%f",ptBin),iBin,iBin);
    //    m2distr->RebinX(2);
    //
    m2distr->SetMarkerColor(kBlue);
    m2distr->SetMarkerStyle(20);
    m2distr->SetMarkerSize(1.4);
    //
    // initialize signal with background parameters
    //
    fTOFsignal.SetParameters(oldParStartValue);
    //
    // improve fit convergence
    //
    fTOFsignal.SetParLimits(1,-0.2,0.2);
    fTOFsignal.SetParLimits(2,0.05,0.8);
    fTOFsignal.SetParLimits(3,0.,1.5);
    //
    fTOFsignal.SetParameter(4, m2distr->GetBinContent(m2distr->GetXaxis()->FindBin(-1.)));
    fTOFsignal.SetParLimits(4, 0.,  m2distr->GetBinContent(m2distr->GetXaxis()->FindBin(0.)));
    //
    fTOFsignal.SetParLimits(7,0.,10.);
    fTOFsignal.SetParLimits(6,0., 5.);
    //
    if (centralityBin < 3 && ptBin < 2.4) {
      fTOFsignal.SetRange(-1.3,2.2);
    } else{
      fTOFsignal.SetRange(-2.2,2.2);
    }
    //
    Int_t fitStatus = m2distr->Fit(&fTOFsignal, "QNR");
    cout << fTOFsignal.GetParameter(3) << endl;
    //
    for(Int_t ival =0; ival < 7; ival++) oldParStartValue[ival] = fTOFsignal.GetParameter(ival);
    //
    canvM2fits->cd();
    m2distr->GetYaxis()->SetRangeUser(0.1, m2distr->GetBinContent(m2distr->GetXaxis()->FindBin(0.))*1.5);
    m2distr->DrawCopy("EP");
    fTOFsignal.DrawCopy("same");
    //
    //
    // bin counting at low pt and subtract the background at high pt
    //
    Float_t yield =1;//  m2distr->Integral(m2distr->GetXaxis()->FindBin(-0.5), m2distr->GetXaxis()->FindBin(+0.5)) - 
      (m2distr->Integral(m2distr->GetXaxis()->FindBin(-1.), m2distr->GetXaxis()->FindBin(-0.5)) + m2distr->Integral(m2distr->GetXaxis()->FindBin(0.5), m2distr->GetXaxis()->FindBin(1.)));
    //
    //
    //
    cout << "pT: " << ptBin 
	 << " , " << fTOFsignal.GetParameter(0)
	 << " , " << fTOFsignal.GetParameter(1)
	 << " , " << fTOFsignal.GetParameter(2)
	 << " , " << fTOFsignal.GetParameter(3)
	 << " , " << fTOFsignal.GetParameter(4)
	 << " , " << fTOFsignal.GetParameter(5)
	 << " , " << fTOFsignal.GetParameter(6)
	 << " , " << fTOFsignal.GetParameter(7) << endl;
    //
    //
    fTOFsignal.SetParameter(4,0.);
    fTOFsignal.SetParameter(5,0.);
    fTOFsignal.SetParameter(6,0.);
    //
    //if (ptBin > 2.3) 
    yield = fTOFsignal.Integral(-1.,1.5)/m2distr->GetBinWidth(10);
    //     yield = m2distr->Integral(m2distr->GetXaxis()->FindBin(-1.), m2distr->GetXaxis()->FindBin(+1.));
    //
    Float_t yieldErr = TMath::Sqrt(yield);
    //    if (fTOFsignal.GetParameter(0) != 0 ) yieldErr = yield*fTOFsignal.GetParError(0)/fTOFsignal.GetParameter(0);
    //
    Float_t maxPt = 4.2;
    if (centralityBin > 2) maxPt = 3.2;
    if (ptBin > 0.5 && ptBin < maxPt && yield > 0) {
      rawSpec->SetBinContent(iBin, yield);
      rawSpec->SetBinError(iBin, yieldErr);
    }
    //
    //
    canvM2fits->Print(Form("QA/M2fits_%i_%i.pdf",sign,centralityBin));
  }
  canvM2fits->Print(Form("QA/M2fits_%i_%i.pdf)",sign,centralityBin));
  //
  rawSpec->Sumw2();
  //
  NormalizeSpectrum(rawSpec,1,norm);
  //
  Int_t colorList[9] = {600+4, 600+3, 600+2, 600+1, 600, 600-4, 600-7, 600-9, 600-10};
  rawSpec->SetMarkerStyle(24);
  rawSpec->SetMarkerSize(1.4);
  rawSpec->SetMarkerColor(colorList[centrality]);
  rawSpec->SetLineColor(colorList[centrality]);
  //
  return rawSpec;

}

//_______________________________________________________________________
void NormalizeSpectrum(TH1D * spectrum, Float_t dy, Float_t numberOfEvents) {
  //
  // make (1 / Nev)*  dN/(dy*dpt)
  //
  spectrum->Scale(1./dy);
  spectrum->Scale(1./numberOfEvents);
  //    
  Int_t nBins = spectrum->GetNbinsX(); 
  for(Int_t i=0; i < nBins+1; i++) {
    Double_t Content =spectrum->GetBinContent(i);
    Double_t error =spectrum->GetBinError(i);
    if (spectrum->GetBinWidth(i) == 0) continue;
    spectrum->SetBinContent(i, Content/spectrum->GetBinWidth(i));
    spectrum->SetBinError(i, error/spectrum->GetBinWidth(i));
  }

}



//_______________________________________________________________________
void MakeEfficiency() {
  //
  // make efficiencies for TPC and TOF part.
  //
  TFile inFileMCold(fileNameMcEnhancedOld);
  TList * listMCold = (TList * ) inFileMCold.Get("akalweit_Nuclei");
  //
  //
  // TPC and tracking as well as matching efficiencies
  //
  //
  TH1D * efficiencyTrackingNegOld = ExtractEfficiency(listMCold, -1, 0, kFALSE);
  TH1D * efficiencyTrackingPosOld = ExtractEfficiency(listMCold, +1, 0, kFALSE);
  efficiencyTrackingNegOld->SetNameTitle("efficiencyTrackingNegOld","efficiencyTrackingNegOld");
  efficiencyTrackingPosOld->SetNameTitle("efficiencyTrackingPosOld","efficiencyTrackingPosOld");
  //
  TH1D * efficiencyTofNegOld = ExtractEfficiency(listMCold, -1, 0, kTRUE);
  TH1D * efficiencyTofPosOld = ExtractEfficiency(listMCold, +1, 0, kTRUE);
  efficiencyTofNegOld->SetNameTitle("efficiencyTofNegOld","efficiencyTofNegOld");
  efficiencyTofPosOld->SetNameTitle("efficiencyTofPosOld","efficiencyTofPosOld");
  //
  TF1 * etaGenCorrection = MakeEtaRapidityGenCorrection(kFALSE);
  efficiencyTrackingPosOld->Divide(etaGenCorrection);
  efficiencyTrackingNegOld->Divide(etaGenCorrection);
  efficiencyTofPosOld->Divide(etaGenCorrection);
  efficiencyTofNegOld->Divide(etaGenCorrection);
  //
  // new simulation
  //
  TFile inFileMCnew(fileNameMcEnhancedNew);
  TList * listMCnew = (TList * ) inFileMCnew.Get("akalweit_Nuclei");
  //
  TH1D * efficiencyTrackingNegNew = ExtractEfficiency(listMCnew, -1, 0, kFALSE);
  TH1D * efficiencyTrackingPosNew = ExtractEfficiency(listMCnew, +1, 0, kFALSE);
  efficiencyTrackingNegNew->SetNameTitle("efficiencyTrackingNegNew","efficiencyTrackingNegNew");
  efficiencyTrackingPosNew->SetNameTitle("efficiencyTrackingPosNew","efficiencyTrackingPosNew");
  efficiencyTrackingNegNew->SetMarkerStyle(21);
  efficiencyTrackingPosNew->SetMarkerStyle(21);
  //
  TH1D * efficiencyTofNegNew = ExtractEfficiency(listMCnew, -1, 0, kTRUE);
  TH1D * efficiencyTofPosNew = ExtractEfficiency(listMCnew, +1, 0, kTRUE);
  efficiencyTofNegNew->SetNameTitle("efficiencyTofNegNew","efficiencyTofNegNew");
  efficiencyTofPosNew->SetNameTitle("efficiencyTofPosNew","efficiencyTofPosNew");
  efficiencyTofNegNew->SetMarkerStyle(24);
  efficiencyTofPosNew->SetMarkerStyle(24);  
  //
  TCanvas * canvEfficiencyTracking = new TCanvas("canvEfficiencyTracking","Efficiencies Tracking");
  efficiencyTrackingNegOld->GetYaxis()->SetRangeUser(0.,1.3);
  efficiencyTrackingNegOld->DrawCopy("E");
  efficiencyTrackingPosOld->DrawCopy("ESame");
  //
  efficiencyTrackingNegNew->DrawCopy("ESame");
  efficiencyTrackingPosNew->DrawCopy("ESame");
  //
  TCanvas * canvEfficiencyTof = new TCanvas("canvEfficiencyTof","Efficiencies Tracking");
  efficiencyTofNegNew->GetYaxis()->SetRangeUser(0.,1.3);
  efficiencyTofNegNew->DrawCopy("E");
  efficiencyTofPosNew->DrawCopy("ESame");
  efficiencyTofNegOld->DrawCopy("ESAME");
  efficiencyTofPosOld->DrawCopy("ESame");
  //
  canvEfficiencyTracking->Print("QAplots/efficiencyTracking.pdf");
  canvEfficiencyTof->Print("QAplots/efficiencyTof.pdf");
  //
  // store output
  //
  TFile * outFile = new TFile("output/efficiencies.root","RECREATE");
  efficiencyTrackingNegNew->Write();
  efficiencyTrackingPosNew->Write();
  efficiencyTofNegNew->Write();
  efficiencyTofPosNew->Write();
  outFile->Close();
  //
  //
  TCanvas * canvEfficiencyMatching = new TCanvas("canvEfficiencyMatching","matching efficiency");
  efficiencyTofNegNew->Divide(efficiencyTrackingNegNew);
  efficiencyTofPosNew->Divide(efficiencyTrackingPosNew);
  efficiencyTofPosNew->DrawCopy();
  efficiencyTofNegNew->DrawCopy("SAME");
  //
  TCanvas * canvEfficiencyRatio = new TCanvas("canvEfficiencyRatio","Efficiencies ratio");
  efficiencyTrackingPosNew->Divide(efficiencyTrackingPosOld);
  efficiencyTrackingPosNew->DrawCopy();
  efficiencyTrackingNegNew->Divide(efficiencyTrackingNegOld);
  efficiencyTrackingNegNew->DrawCopy("SAME");

  

}



//_______________________________________________________________________
TH1D * ExtractEfficiency(TList * listMC, Int_t sign, Int_t centrality, Bool_t isTOF) {
  //
  // Extract the efficiency
  //
  //
  // (0.) assumed particle: 0. deuteron, 1. triton, 2. He-3
  // (1.) multiplicity or centrality -- number of accepted ESD tracks per events (deprecated), but now classes from 1 to 10, 0: Min. Bias
  // (2.) pT
  // (3.) sign
  // (4.) rapidity --> filled 4x
  // (5.)  pull TPC dEx --> filled 4x
  // (6.) has valid TOF pid signal
  // (7.) nsigma TOF --> filled 4x XXXXXXXXX no mass*mass
  // (8..) dca_xy
  // (9.) CODE -- only MC 0-generated, 1-true rec. primaries, 2-misident prim, 3-second weak, 4-second material, 5-misident sec
  //
  THnSparseF * hist = (THnSparseF *) listMC->FindObject("fHistMCparticles");
  //
  // common selections for generated and reconstructed particles
  // 
  hist->GetAxis(0)->SetRangeUser(0,0);                     // select deuterons
  hist->GetAxis(4)->SetRangeUser(-0.49, 0.49);             // select rapidity range
  hist->GetAxis(3)->SetRangeUser(sign,sign);               // select anti-deuterons
  //assume eff. indepedent of centrality  ---> hist->GetAxis(1)->SetRangeUser(centrality,centrality);   // select centrality bin
  //
  hist->GetAxis(9)->SetRangeUser(0,0);
  TH1D * generated = hist->Projection(2);
  generated->SetNameTitle(Form("generated_cent%i_sign%i",centrality,sign),
			  Form("generated_cent%i_sign%i",centrality,sign));
  //
  // further cuts on reconstructed particles
  //
  hist->GetAxis(8)->SetRangeUser(-0.5,0.5);  // DCA-cut / TODO: do proper unfolding
  //
  hist->GetAxis(5)->SetRangeUser(-4.5,4.5);  // TPC-PID cut  
  //
  if (isTOF) {
    //    hist->GetAxis(5)->SetRangeUser(-2.,2.);  // TPC-PID cut  
    hist->GetAxis(6)->SetRangeUser(1.,1.);  // TOF-PID cut  
    hist->GetAxis(7)->SetRangeUser(-1.2,1.2);  // TOF-PID cut  
  }
  //
  hist->GetAxis(9)->SetRangeUser(1,1);
  TH1D * rawSpec = hist->Projection(2);
  rawSpec->Sumw2();
  rawSpec->SetNameTitle(Form("eff_cent%i_sign%i",centrality,sign),
			Form("eff_cent%i_sign%i",centrality,sign));
  //
  rawSpec->Divide(generated);
  //
  if (sign > 0) {
    rawSpec->SetMarkerColor(kRed);
    rawSpec->SetLineColor(kRed);
  } else{
    rawSpec->SetMarkerColor(kBlue);
    rawSpec->SetLineColor(kBlue);
  }
  //
  // reset ranges
  //
  hist->GetAxis(5)->SetRangeUser(-1000.,1000);  // TPC-PID cut  
  hist->GetAxis(6)->SetRangeUser(-1000.,1000);  // TOF-PID cut  
  hist->GetAxis(7)->SetRangeUser(-1000.,1000);  // TOF-PID cut    
  hist->GetAxis(8)->SetRangeUser(-1000.,1000);  // TOF-PID cut    
  //
  //
  //
  return rawSpec;

}

//__________________________________________________________________
TF1 * MakeEtaRapidityGenCorrection(Bool_t drawQA) {
  //
  // the enhanced sample LHC11b9_1 has a stupid cut on the generated level
  //
  TF1 * funcDeut = new TF1("funcDeut","TMath::ASinH(TMath::Sqrt(1 + (1.876*1.876)/(x*x))*TMath::SinH(0.5))",0.5,5);
  TF1 * etaCut = new TF1("etaCut","0.9",0.5,5);
  //
  funcDeut->SetLineColor(kBlue);
  //
  if (drawQA) {
    TCanvas * canv1 = new TCanvas("canv1","original functions");
    funcDeut->DrawCopy();
    funcDeut->GetHistogram()->GetYaxis()->SetRange(0,1.8);
    etaCut->DrawCopy("same");
  }
  //
  TF1 * funcCorr = new TF1("funcCorr","funcDeut/etaCut*(x - 1.1 < 0) + 1*(x - 1.1 >= 0)",0.5,2);
  if (drawQA) {
    TCanvas * canv2 = new TCanvas("canv2","correction function");
    funcCorr->Draw();
  }
  //
  return funcCorr;

}


//__________________________________________________________________
void MakeQAplots() {
  //
  // make the QA plots
  //
  TCanvas * canvData          =  PlotTpcQA(fileNameData,          kFALSE);
  TCanvas * canvMcEnhancedNew =  PlotTpcQA(fileNameMcEnhancedNew, kTRUE);
  TCanvas * canvMcUnEnhanced  =  PlotTpcQA(fileNameMcUnEnhanced,  kTRUE, kTRUE);
  TCanvas * canvMcEnhancedOld =  PlotTpcQA(fileNameMcEnhancedOld, kTRUE, kTRUE);
  //
  TCanvas * canvQA = new TCanvas("canvQA","canvQA");
  canvQA->Print("QAplots/QAplots.pdf(");
  //
  // add all QA plots to file
  //
  canvData->Print("QAplots/QAplots.pdf");
  canvMcEnhancedNew->Print("QAplots/QAplots.pdf");
  canvMcUnEnhanced->Print("QAplots/QAplots.pdf");
  canvMcEnhancedOld->Print("QAplots/QAplots.pdf");
  //
  canvQA->Print("QAplots/QAplots.pdf)");
  //
  // spectra QA
  //
  TFile inFile(fileNameData);
  TList * list = (TList * ) inFile.Get("akalweit_Nuclei");
  PlotSpectraQA(list);


}


//__________________________________________________________________
TCanvas * PlotTpcQA(Char_t * fileName, Bool_t isMC, Bool_t isMCold) {
  //
  // plot QA of TPC-PID
  //
  TFile inFileMC(fileName);
  TList * list = (TList * ) inFileMC.Get("akalweit_Nuclei");
  //
  TF1 * fBBlines[6];
  const Double_t masses[6] = {0.000511, 0.138, 0.4936, 0.938, 2*0.938, 3*0.938};
  //
  // set BB parameters
  //
  Double_t parData[5] = {1.45802, 27.4992, 4.00313e-15, 2.48485, 8.31768};
  Double_t parMC[5] = {1.17329, 27.4992, 4.00313e-15, 2.1204316, 4.1373729}; 
  Double_t parMCold[5] = {1.17329, 27.4992, 4.00313e-15, 2.35563, 9.47569}; // OLD FOR LHC11b9_1 !!
  //
  for(Int_t ii = 0; ii < 6; ii++) {
    fBBlines[ii] = new TF1(Form("fBB_%i",ii),"AliExternalTrackParam::BetheBlochAleph(x/([0]),[1],[2],[3], [4], [5])",0.3,10);
    
    fBBlines[ii]->SetParameters(masses[ii], parData[0], parData[1], parData[2], parData[3], parData[4]);
    if (isMC) {
      fBBlines[ii]->SetParameters(masses[ii], parMC[0], parMC[1], parMC[2], parMC[3], parMC[4]);
      if (isMCold) fBBlines[ii]->SetParameters(masses[ii], parMCold[0], parMCold[1], parMCold[2], parMCold[3], parMCold[4]);
    }
    fBBlines[ii]->SetLineColor(kRed);
    fBBlines[ii]->SetLineWidth(2);
  }
  //
  TCanvas * canvDeDx = new TCanvas(Form("canvDeDxTpcQA_%s",fileName),"control histogram for dE/dx");
  canvDeDx->Divide(1,2);
  canvDeDx->cd(1);
  TH3D * histTPC = (TH3D *) list->FindObject("fHistPidQA");
  histTPC->GetZaxis()->SetRangeUser(-1,-1);
  TH2D * histNegTPC = (TH2D *) histTPC->Project3D("yx");
  histNegTPC->SetNameTitle("histNegTPC", Form("histNegTPC_%s", fileName));
  gPad->SetLogx();
  histNegTPC->DrawCopy("colZ");
  for(Int_t ii = 0; ii < 6; ii++) fBBlines[ii]->DrawCopy("same");
  //
  canvDeDx->cd(2);
  histTPC->GetZaxis()->SetRangeUser(+1,+1);
  TH2D * histPosTPC = (TH2D*) histTPC->Project3D("yx");
  histPosTPC->SetNameTitle("histPosTPC", Form("histPosTPC_%s", fileName));
  gPad->SetLogx();
  histPosTPC->DrawCopy("colZ");
  //  for(Int_t ii = 0; ii < 6; ii++) fBBlines[ii]->DrawCopy("same");
  fBBlines[4]->DrawCopy("same"); // Draw deuteron line
  //
  // He-3 goes extra
  //
  TF1 * funcHel3 = new TF1("funcHel3","4*AliExternalTrackParam::BetheBlochAleph(2*x/(0.938*3),1.74962,27.4992,4.00313e-15,2.42485,8.31768)",0.2,6.);
  if (isMC) funcHel3 = new TF1("funcHel3MC","4*AliExternalTrackParam::BetheBlochAleph(2*x/(0.938*3), 1.17329, 27.4992, 4.00313e-15, 2.1204316, 4.1373729)",0.2,6.);
  //
  fBBlines[5]->DrawCopy("same"); // Draw deuteron line
  funcHel3->DrawCopy("same");

  return canvDeDx;

}



//_______________________________________________________________________
void PlotSpectraQA(TList * list, Int_t particleType, Int_t sign) {
  //
  // Make some basic QA plots for TOF related PID
  //
  THnSparseF * hist = (THnSparseF *) list->FindObject("fHistRealTracks");
  hist->GetAxis(4)->SetRangeUser(-0.49, 0.49);             // select Rapidity range
  hist->GetAxis(0)->SetRangeUser(particleType,particleType); // select deuterons
  hist->GetAxis(3)->SetRangeUser(sign,sign); // select deuterons
  //
  TCanvas * canvSpectraQA1 = new TCanvas("canvSpectraQA1","QA for nsigma TPC");
  TH2F * histTpcNsigma = (TH2F *) hist->Projection(5,2);
  histTpcNsigma->SetNameTitle("histTpcNsigma","histTpcNsigma");
  histTpcNsigma->DrawCopy("colZ");
  //
  TCanvas * canvSpectraQA2 = new TCanvas("canvSpectraQA2","QA for nsigma TOF");
  hist->GetAxis(6)->SetRangeUser(1,1);   // require proper match to TOF
  hist->GetAxis(5)->SetRangeUser(-4.,4); // pre-select tracks with TPC-PID
  TH2F * histTofNsigma = (TH2F *) hist->Projection(7,2);
  histTofNsigma->RebinY(4);
  histTofNsigma->SetNameTitle("histTofNsigma","histTofNsigma");
  histTofNsigma->DrawCopy("colZ");

}



//_______________________________________________________________________
void MakeMaterialCorrection(){
  //
  // compare fitted values with MC pt-dependence
  //
  TFile outFile("output/materialCorrection.root","RECREATE");
  outFile.Close();
  //
  Int_t kMaxCentrality = 5;
  TFile * fileIn = TFile::Open("output/efficiencies.root");
  //
  //
  for(Int_t iCentr = 0; iCentr < kMaxCentrality; iCentr++) {
    //
    TH1D * matCorrMC = (TH1D *) fileIn->Get("efficiencyTrackingNegNew");
    matCorrMC->Reset();
    matCorrMC->SetNameTitle(Form("matCorrMC_%i",iCentr), Form("matCorrMC_i",iCentr));
    //
    TH1D * matCorr = (TH1D *) fileIn->Get("efficiencyTrackingNegNew");
    matCorr->SetLineColor(kRed);
    matCorr->Reset();
    matCorr->SetNameTitle(Form("matCorr_%i",iCentr), Form("matCorr_i",iCentr));    
    //
    for(Int_t iBin = 0; iBin < matCorrMC->GetXaxis()->GetNbins(); iBin++) {
      Double_t pT = matCorrMC->GetXaxis()->GetBinCenter(iBin);
      if (pT > 0.6 && pT < 2.8) {
	Float_t corr = GetMaterialCorrection(pT,iCentr,kTRUE);
	matCorrMC->SetBinContent(iBin, corr);
	if (pT < 2.) {
	  Float_t corrReal = GetMaterialCorrection(pT,iCentr,kFALSE);
	  matCorr->SetBinContent(iBin, corrReal);
	}
      }
    }
    //
    //
    TCanvas * canvMaterial = new TCanvas("canvMaterial","canvMaterial");
    canvMaterial->cd();
    if (iCentr == 0) {
      matCorr->DrawCopy();
      matCorrMC->DrawCopy("SAME");
    } else {
      matCorrMC->DrawCopy("SAME");
    }
    TFile * outputFile = new TFile("output/materialCorrection.root","UPDATE");
    matCorr->Write();
    outputFile->Close();
    delete outputFile;
  }


}



//_______________________________________________________________________
Float_t GetMaterialCorrection(Float_t ptBin, Int_t centralityBin, Bool_t pureMC) {
  //
  // subtract the material contamination from positive particles
  //
  //
  Float_t dcaCutInAnalysis = 0.5;
  //
  TFile * inFileData = TFile::Open("data/dataFinal.root");
  TList * listData = (TList *) inFileData->Get("akalweit_Nuclei");
  THnSparse * histData = (THnSparse *) listData->FindObject("fHistRealTracks");
  //
  TFile * inFileMC = TFile::Open("MC/McCombined.root");
  TList * listMC = (TList *) inFileMC->Get("akalweit_Nuclei");
  THnSparse * histMC = (THnSparse *) listMC->FindObject("fHistMCparticles");
  //
  //
  // (1.) select deuterons
  //
  histData->GetAxis(0)->SetRangeUser(0,0);
  histMC->GetAxis(0)->SetRangeUser(0,0);
  histData->GetAxis(3)->SetRangeUser(+1,+1); // sign
  histMC->GetAxis(3)->SetRangeUser(+1,+1);   // sign
  //
  histData->GetAxis(4)->SetRangeUser(-0.49, 0.49);             // select rapidity range
  histMC->GetAxis(4)->SetRangeUser(-0.49, 0.49);             // select rapidity range
  //
  //  hist->GetAxis(8)->SetRangeUser(-1.,1.);  // DCA-cut / TODO: do proper unfolding
  histData->GetAxis(5)->SetRangeUser(-2.5,2.5);  // TPC-PID cut  
  histData->GetAxis(6)->SetRangeUser(1.,1.);  // TOF-PID cut  
  histData->GetAxis(7)->SetRangeUser(-0.5,0.5);  // TOF-PID cut  
  //
  // (2.) select pt-range
  //
  histData->GetAxis(2)->SetRangeUser(ptBin,ptBin);
  histMC->GetAxis(2)->SetRangeUser(ptBin,ptBin);
  //
  // (4.a) get different MC templates
  //
  histMC->GetAxis(9)->SetRangeUser(1,1);
  TH1D * prim = (TH1D*) histMC->Projection(8);
  prim->SetNameTitle("prim","prim");
  prim->SetLineColor(kRed);
  //
  histMC->GetAxis(9)->SetRangeUser(4,4);
  //
  // centrality can oonly be selected after the primary yield
  //
  Int_t centrality = 0;
  if (centralityBin == 0) centrality = 0; //  0-10%
  if (centralityBin == 1) centrality = 2; // 10-20%
  if (centralityBin == 2) centrality = 3; // 20-40%
  if (centralityBin == 3) centrality = 5; // 40-60%
  if (centralityBin == 4) centrality = 7; // 60-80%
  //
  histData->GetAxis(1)->SetRangeUser(centrality,centrality);       // select centrality bin
  histMC->GetAxis(1)->SetRangeUser(centrality,centrality);       // select centrality bin
  if (centralityBin != 1) {
    histData->GetAxis(1)->SetRangeUser(centrality,centrality+1);   // select centrality bin
    histMC->GetAxis(1)->SetRangeUser(centrality,centrality+1);     // select centrality bin
  }
  //
  TH1D * material = (TH1D*) histMC->Projection(8);
  histMC->GetAxis(2)->SetRangeUser(ptBin,ptBin);
  material->SetNameTitle("material","material");  
  material->SetLineColor(kGreen);
  //
  TH1D * dcaDistr = (TH1D*) histData->Projection(8);
  //
  // if only MC, we just return the values
  //
  if (pureMC) {
    Int_t fractionBinLow = dcaDistr->FindBin(-dcaCutInAnalysis);
    Int_t fractionBinUp  = dcaDistr->FindBin(+dcaCutInAnalysis);
    Float_t materialFraction = 0.;
    prim->Add(material); // has to be scaled with number of events somehow
    materialFraction = material->Integral(fractionBinLow,fractionBinUp)/prim->Integral(fractionBinLow,fractionBinUp);
    //
    delete listData;
    delete listMC;
    inFileData->Close();
    inFileMC->Close();
    delete inFileData;
    delete inFileMC;
    //
    return materialFraction;
  }
  //
  // (6.) do the fit
  //
  Float_t dcaUp = 0.5;//0.0182 + 0.035/TMath::Power(ptBin,1.01);
  Float_t dcaLow  = -1.*dcaUp;
  Int_t binLow = dcaDistr->GetXaxis()->FindBin(dcaLow);
  Int_t binUp  = dcaDistr->GetXaxis()->FindBin(dcaUp);
  //
  TObjArray *mcShapes = new TObjArray(3);        // MC histograms are put in this array
  mcShapes->Add(prim);
  mcShapes->Add(material);
  //
  TFractionFitter* fit = new TFractionFitter(dcaDistr, mcShapes); // initialise
  fit->SetRangeX(binLow,binUp);
  fit->Constrain(1,0.,1.);
  fit->Constrain(2,0.,1.);
  //
  //
  //
  Int_t status = fit->Fit();
  if (status==0) fit->GetPlot()->Draw("same");
  //
  // (7.) Rescaling of histograms and final plots
  //
  TH1D * result = 0x0;
  Double_t yieldPrim, yieldSec, yieldMat, error;
  if (status == 0) {
    fit->GetResult(0,yieldPrim,error);
    fit->GetResult(1,yieldMat,error);
    result = (TH1D*) fit->GetPlot();
    result->SetLineWidth(3);
  }
  //
  prim->Scale(1./prim->Integral(binLow,binUp));
  if (status == 0) prim->Scale(yieldPrim*result->Integral(binLow,binUp));
  //
  material->Scale(1./material->Integral(binLow,binUp));
  if (status == 0) material->Scale(yieldMat*result->Integral(binLow,binUp));
  //
  NormalizeSpectrum(dcaDistr,1,1);
  NormalizeSpectrum(prim,1,1);
  NormalizeSpectrum(material,1,1);
  if (status == 0) NormalizeSpectrum(result,1,1);
  //
  TH1D * sumHist = (TH1D*) prim->Clone(); // result is not equal to weighted sum ---> see FractionFitter doku and paper (?)
  sumHist->SetNameTitle("sumHist","sumHist");
  sumHist->Add(material);
  sumHist->SetLineColor(kOrange);
  //
  TCanvas * canvFeed = new TCanvas("canvFeed", "feed down and material correction");
  dcaDistr->SetMarkerStyle(24);  dcaDistr->SetMarkerSize(1.4);
  dcaDistr->DrawCopy("Ep");
  prim->DrawCopy("same");
  material->DrawCopy("same");
  sumHist->DrawCopy("same");
  //
  // Final-Results
  //
  Int_t fractionBinLow = dcaDistr->FindBin(-dcaCutInAnalysis);
  Int_t fractionBinUp  = dcaDistr->FindBin(+dcaCutInAnalysis);
  //
  Float_t materialFraction = 0.;
  if (status == 0) {
    result->DrawCopy("same");
    materialFraction = material->Integral(fractionBinLow,fractionBinUp)/dcaDistr->Integral(fractionBinLow,fractionBinUp);
  }
  Printf("RESULTING MATERIAL FRACTION: %f", materialFraction);
  Printf("TEMPLATE  MATERIAL FRACTION: %f", yieldMat);
  //
  delete sumHist;
  delete result;
  delete mcShapes;
  delete fit;
  //  delete histData;
  //delete histMC;
  delete listData;
  delete listMC;
  inFileData->Close();
  inFileMC->Close();
  delete inFileData;
  delete inFileMC;
  //
  return materialFraction;

}
