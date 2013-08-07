

void epsilonScanFineBinningForPerformanceNote(Int_t radiusScale=1, Int_t epsScale=1,Int_t iOmegaTau = 0, Double_t x = 85., Double_t y = 0., Double_t z = 10.,Double_t integrateStep = 0.,Bool_t bLowE = kFALSE){
  //
  // do for one gas (with omegaTau values accordingly) the correction
  // use the integrate along drift line option if integrateStep > 0
  // else old method
  //
  // 1. Initialzation form space charge maps
  //
  AliTPCSpaceCharge3D *spaceCharge = new AliTPCSpaceCharge3D;

  // omega tau parameters and TF1 
  const Int_t nOmegaTau = 1;
  Double_t omegaTau[nOmegaTau] = {0.32};
  TString tGas[nOmegaTau] = {"NeCO2_2"}; 
  TString tGasOut[nOmegaTau] = {"NeCO2_2"}; 
  TString sGas[nOmegaTau] = {"Ne-CO_{2}-N_{2} (90-10-5)"};
  TString sLowE = "";
  //if(bLowE) sLowE = "_lowE";
  if(bLowE) sLowE = "_lowE2"; //second iteration

  const Int_t nEps = 51;
  Double_t epsilon[nEps];
  for(Int_t iEps = 0; iEps < nEps; iEps++){
    epsilon[iEps] = (Double_t)iEps;
  }
  TH1F *hR    = new TH1F("hR",Form("distortion at x = %.0f, y = %.0f, z = %.0f",x,y,z),nEps,epsilon[0],epsilon[nEps-1]+1);
  TH1F *hPhiR = new TH1F("hPhiR",Form("distortion at x = %.0f, y = %.0f, z = %.0f",x,y,z),nEps,epsilon[0],epsilon[nEps-1]+1);

  //use always the integrate option here
  Double_t rate = 50.;

  for(Int_t iEps = 0; iEps < nEps; iEps++){

    // select gas 
    // 0 = NeCO2N2 
    cout<<"Open file = "<<Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%.0f_50kHz_radiusScaling%d_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),epsilon[iEps],radiusScale,epsScale)<<endl;
    spaceCharge->SetSCDataFileName(Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%.0f_50kHz_radiusScaling%d_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),epsilon[iEps],radiusScale,epsScale));
    
    // select omegaTau value
    if(iOmegaTau ==  0){
      spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],1,0.99); // Ne CO2 N2
    }
 
    //
    // init and add to corrections
    spaceCharge->InitSpaceCharge3DDistortion();
    spaceCharge->AddVisualCorrection(spaceCharge,1);
    
    //
    // 2. get TH1F with differences 
    //
    // get corrections (at x, y and z) for visual correction 1 (last argument)
    // 0 = dR
    // 1 = dPhiR
    // 2 = dZ?    
    if(integrateStep < 0.0001){// == 0. --> old method
      hR->Fill(epsilon[iEps],AliTPCCorrection::GetDistXYZ(x,y,z,0,1));
      hPhiR->Fill(epsilon[iEps],AliTPCCorrection::GetDistXYZ(x,y,z,1,1));
    }
    else{// > 0. --> following driftlines
      hR->Fill(epsilon[iEps],AliTPCCorrection::GetDistXYZIntegrateZ(x,y,z,0,1,integrateStep));
      hPhiR->Fill(epsilon[iEps],AliTPCCorrection::GetDistXYZIntegrateZ(x,y,z,1,1,integrateStep));
    }
  }


  //
  // 3. Plot and store TH1Fs
  //

  TCanvas *cEpsilon = new TCanvas("cEpsilon","cEpsilon",1200,500);
  cEpsilon->Divide(2,1);

  TLegend *legend = new TLegend(0.6,0.3,0.85,0.75,"","brNDC");
  setupLegend(legend,0.04);

  TString outfilename = Form("epsilonScanFineBinning_performanceNote_radiusScaling%d_epsScaling%d_%s%s_x%0.f_y%.0f_z%.0f_int%.0f.root",radiusScale,epsScale,tGasOut[iOmegaTau].Data(),sLowE.Data(),x,y,z,integrateStep);
  TFile *fOut = TFile::Open(Form("%s",outfilename.Data()),"recreate");
 
  cEpsilon->cd(1);
  hR->Draw();

  cEpsilon->cd(2);
  hPhiR->Draw();


  cEpsilon->SaveAs(Form("epsilonScanFineBinning_performanceNote_radiusScaling%d_epsScaling%d_%s%s_x%0.f_y%.0f_z%.0f_int%.0f.eps",radiusScale,epsScale,tGasOut[iOmegaTau].Data(),sLowE.Data(),x,y,z,integrateStep));
  
  hR->Write();
  hPhiR->Write();
  fOut->Close();
}

void epsilonScanFineBinningPlotOnly(Int_t radiusScale=1, Int_t epsScale=1,Int_t iOmegaTau = 0, Double_t x = 85., Double_t y = 0., Double_t z = 10.){
  //
  // do for one gas (with omegaTau values accordingly) the correction
  // use the integrate along drift line option

  gStyle->SetOptStat(0);

  //
  // 1. Plot and store TH1Fs
  //

  // omega tau parameters and TF1 
  const Int_t nOmegaTau = 1;
  Double_t omegaTau[nOmegaTau] = {0.32};
  TString tGas[nOmegaTau] = {"NeCO2_2"}; 
  TString sGas[nOmegaTau] = {"Ne/CO_{2}/N_{2} (90-10-5)"};

 
  TCanvas *cEpsilon = new TCanvas("cEpsilon","cEpsilon",1200,500);
  cEpsilon->Divide(2,1);

  TLegend *legend = new TLegend(0.15,0.6,0.6,0.85,Form("%s: 50 kHz",sGas[iOmegaTau].Data(),"brNDC"));
  setupLegend(legend,0.05);

  TString infilename0 = Form("epsilonScanFineBinning_performanceNote_radiusScaling%d_epsScaling%d_%s_x%0.f_y%.0f_z%.0f_int0.root",radiusScale,epsScale,tGas[iOmegaTau].Data(),x,y,z);
  TString infilename1 = Form("epsilonScanFineBinning_performanceNote_radiusScaling%d_epsScaling%d_%s_x%0.f_y%.0f_z%.0f_int1.root",radiusScale,epsScale,tGas[iOmegaTau].Data(),x,y,z);
  TFile *fIn0 = TFile::Open(Form("%s",infilename0.Data()),"read");
  TFile *fIn1 = TFile::Open(Form("%s",infilename1.Data()),"read");

  TH1F *hR0    = (TH1F*)fIn0->Get("hR");
  TH1F *hPhiR0 = (TH1F*)fIn0->Get("hPhiR");
 
  TH1F *hR1    = (TH1F*)fIn1->Get("hR");
  TH1F *hPhiR1 = (TH1F*)fIn1->Get("hPhiR");

  hR0->SetTitle(Form("x = %.0f, y = %.0f, z = %.0f",x,y,z)); 
  hPhiR0->SetTitle(Form("x = %.0f, y = %.0f, z = %.0f",x,y,z)); 

  hR1->SetTitle(Form("x = %.0f, y = %.0f, z = %.0f",x,y,z)); 
  hPhiR1->SetTitle(Form("x = %.0f, y = %.0f, z = %.0f",x,y,z)); 

  cEpsilon->cd(1);
  hR0->GetXaxis()->SetTitle("#varepsilon");
  hR0->GetYaxis()->SetTitle("dr (cm)");
  hR0->SetLineWidth(2);
  hR0->DrawCopy("lp");
  hR1->SetLineWidth(2);
  hR1->SetLineColor(2);
  hR1->DrawCopy("lp,same");

  legend->AddEntry(hR0,"linear","l");
  legend->AddEntry(hR1,"integrated in z","l");

  legend->Draw();

  cEpsilon->cd(2);
  hPhiR0->GetXaxis()->SetTitle("#varepsilon");
  hPhiR0->GetYaxis()->SetTitle("d(r#varphi) (cm)");
  hPhiR0->SetLineWidth(2);
  hPhiR0->DrawCopy("lp");
  hPhiR1->SetLineWidth(2);
  hPhiR1->SetLineColor(2);
  hPhiR1->DrawCopy("lp,same");


  cEpsilon->SaveAs(Form("epsilonScanFineBinning_performanceNote_radiusScaling%d_epsScaling%d_%s_x%0.f_y%.0f_z%.0f.eps",radiusScale,epsScale,tGas[iOmegaTau].Data(),x,y,z));
  cEpsilon->SaveAs(Form("epsilonScanFineBinning_performanceNote_radiusScaling%d_epsScaling%d_%s_x%0.f_y%.0f_z%.0f.pdf",radiusScale,epsScale,tGas[iOmegaTau].Data(),x,y,z));
  cEpsilon->SaveAs(Form("epsilonScanFineBinning_performanceNote_radiusScaling%d_epsScaling%d_%s_x%0.f_y%.0f_z%.0f.png",radiusScale,epsScale,tGas[iOmegaTau].Data(),x,y,z));
  
  fIn0->Close();
  fIn1->Close();
}




  
//____________________________________________________________//
void setupLegend(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}
