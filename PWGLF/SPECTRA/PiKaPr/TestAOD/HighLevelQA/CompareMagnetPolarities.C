//////////////////////////////////////////////////////////////////////////////////////
// CompareMagnetPolarities.C                                                        //
// computes and draws the ratio of pions from runs with different magnet polarities //
// written by John Groh                                                             //
//////////////////////////////////////////////////////////////////////////////////////

const Int_t nCharge = 2;
TString Sign[nCharge] = {"Plus","Minus"};
const Int_t nRuns = 2;
const Int_t Runs[nRuns] = {138275, 139465};
TString Names[nCharge] = {"#pi^{+}","#pi^{-}"};


void CompareMagnetPolarities()
{
  // get pion spectra histos for run 138275
  TH1F * hSpectra138275[nCharge];
  TFile * fin138275 = TFile::Open("results/SeeIfPionsVaryWithMagnetPolarity/Res_OutputAODSpectraTask_Data_Cent30to40_QVec0.0to100.0_Eta-0.8to0.8_3.0SigmaPID_TrBit1024_SeeIfPionsVaryWithMagnetPolarity_Run138275.root");
  TCanvas * cSpectra138275 = (TCanvas*)fin138275->Get("cSpectra");
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      hSpectra138275[icharge] = (TH1F*)cSpectra138275->FindObject(Form("Spectra_Pion%s",Sign[icharge].Data()));
    }
  fin138275->Close();

  // repeat for run 139465
  TH1F * hSpectra139465[nCharge];
  TFile * fin139465 = TFile::Open("results/SeeIfPionsVaryWithMagnetPolarity/Res_OutputAODSpectraTask_Data_Cent30to40_QVec0.0to100.0_Eta-0.8to0.8_3.0SigmaPID_TrBit1024_SeeIfPionsVaryWithMagnetPolarity_Run139465.root");
  TCanvas * cSpectra139465 = (TCanvas*)fin139465->Get("cSpectra");
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      hSpectra139465[icharge] = (TH1F*)cSpectra139465->FindObject(Form("Spectra_Pion%s",Sign[icharge].Data()));
    }
  fin139465->Close();
  
  // repeat for run 138364
  TH1F * hSpectra138364[nCharge];
  TFile * fin138364 = TFile::Open("results/SeeIfPionsVaryWithMagnetPolarity/Res_OutputAODSpectraTask_Data_Cent30to40_QVec0.0to100.0_Eta-0.8to0.8_3.0SigmaPID_TrBit1024_SeeIfPionsVaryWithMagnetPolarity_Run138364.root");
  TCanvas * cSpectra138364 = (TCanvas*)fin138364->Get("cSpectra");
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      hSpectra138364[icharge] = (TH1F*)cSpectra138364->FindObject(Form("Spectra_Pion%s",Sign[icharge].Data()));
    }
  fin138364->Close();

  // ratio 1 - (late negative) / (late positive)
  TCanvas * cRatioPions_1 = new TCanvas("cRatioPions_1","cRatioPions_1");
  TLegend * lRatioPions_1 = new TLegend(.79,.79,.99,.99);
  TH1F * hRatio[nCharge];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      hRatio[icharge] = (TH1F*)hSpectra138275[icharge]->Clone();
      hRatio[icharge]->Divide(hSpectra139465[icharge]);
      hRatio[icharge]->SetTitle(Form("Ratio of Spectra: #frac{Run %i}{Run %i};p_{T} (GeV/c);",Runs[0],Runs[1]));
      hRatio[icharge]->SetLineStyle(icharge+1);
      hRatio[icharge]->GetYaxis()->SetRangeUser(.7,1.3);
      if (icharge == 0) hRatio[icharge]->DrawCopy("hist][");
      else hRatio[icharge]->DrawCopy("hist][same");
      lRatioPions_1->AddEntry(hRatio[icharge],Names[icharge].Data(),"l");
    }
  lRatioPions_1->DrawClone();

  // another ratio - (early positive) / (late positive)
  TCanvas * cRatioPions_2 = new TCanvas("cRatioPions_2","cRatioPions_2");
  TLegend * lRatioPions_2 = new TLegend(.79,.79,.99,.99);
  TH1F * hRatio[nCharge];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      hRatio[icharge] = (TH1F*)hSpectra138364[icharge]->Clone();
      hRatio[icharge]->Divide(hSpectra139465[icharge]);
      hRatio[icharge]->SetTitle(Form("Ratio of Spectra: #frac{Run %i}{Run %i};p_{T} (GeV/c);",Runs[0],Runs[1]));
      hRatio[icharge]->SetLineStyle(icharge+1);
      hRatio[icharge]->GetYaxis()->SetRangeUser(.7,1.3);
      if (icharge == 0) hRatio[icharge]->DrawCopy("hist][");
      else hRatio[icharge]->DrawCopy("hist][same");
      lRatioPions_2->AddEntry(hRatio[icharge],Names[icharge].Data(),"l");
    }
  lRatioPions_2->DrawClone();


}
