void SCdistortionsForPerformanceNoteAll(){

  SCdistortionsForPerformanceNote(1,1);
  SCdistortionsForPerformanceNote(2,1);
  SCdistortionsForPerformanceNote(1.5,1);

}

void SCdistortionsForPerformanceNote(Double_t radiusScale=1.5, Int_t epsScale=1,Int_t iOmegaTau = 0){
  //
  // do for one file (given by directory) the correction (specify the gas = iOmegaTau)
  // use the integrate along drift line option
  //
  // 1. Initialzation form space charge maps
  //
  AliTPCSpaceCharge3D *spaceCharge = new AliTPCSpaceCharge3D;
  const Double_t fgke0 = 8.854187817e-12; // vacuum permittivity [A·s/(V·m)]

  // omega tau parameters and TF1
  const Int_t nEps = 2;
  Int_t eps[nEps] = {20,10};
  Int_t col[nEps] = {kBlack,kRed};
  Int_t sty[nEps] = {1,2};

  const Int_t nOmegaTau = 1;
  Double_t omegaTau[nOmegaTau] = {0.32};
  TString tGas[nOmegaTau] = {"NeCO2_2"}; 
  TString sGas[nOmegaTau] = {"Ne-CO_{2}-N_{2} (90-10-5)"};
  TF1 * fdiffR[nEps];
  TF1 * fdiffPhiR[nEps];
  TH1F * hdiffR[nEps];
  TH1F * hdiffPhiR[nEps];
  TH2F * hMap[nEps];

  //use always the integrate option here
  Double_t integrateStep = 1.;

  TCanvas *cMap[nEps];

  TCanvas *cOmegaTau = new TCanvas("cOmegaTau","cOmegaTau",1200,500);
  cOmegaTau->Divide(2,1);

  TString outfilename = Form("SCdistortions_SC_NeCO2_50kHz_radiusScaling%.0f_epsScaling%d",radiusScale,epsScale);
  if(radiusScale>1.1 && radiusScale < 1.9) outfilename = Form("SCdistortions_SC_NeCO2_50kHz_radiusScaling%.1f_epsScaling%d",radiusScale,epsScale);

  TLegend *legend = new TLegend(0.35,0.6,0.8,0.85,Form("%s: 50 kHz",sGas[iOmegaTau].Data(),"brNDC"));
  setupLegend(legend,0.05);

  //loop over epsilons
  for(Int_t iEps = 0; iEps < nEps; ++iEps){

    cMap[iEps] = new TCanvas(Form("cMap%d",iEps),Form("cMap%d",iEps),600,500);
    cMap[iEps]->SetPhi(150);
    
    // select gas 
    // 0 = NeCO2N2 

    if(radiusScale>1.1 && radiusScale < 1.9){
      cout<<"Open file = "<<Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.1f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iEps],radiusScale,epsScale)<<endl;
      spaceCharge->SetSCDataFileName(Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.1f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iEps],radiusScale,epsScale));
    }
    else{
      cout<<"Open file = "<<Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.0f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iEps],radiusScale,epsScale)<<endl;
      spaceCharge->SetSCDataFileName(Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.0f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iEps],radiusScale,epsScale));
    }

    // select omegaTau value
    if(iOmegaTau ==  0){
      spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],1,0.99); // Ne CO2 N2
    }
  
    //
    //
    // init and add to corrections
    spaceCharge->InitSpaceCharge3DDistortion();
    spaceCharge->AddVisualCorrection(spaceCharge,1);
    
    // draw the map
    
    hMap[iEps] = (TH2F*)spaceCharge->CreateHistoSCinZR(0.);
    hMap[iEps]->Scale(fgke0/1e6*1e15); // C/m^3/e0 --> fC/cm^3
    hMap[iEps]->GetXaxis()->SetTitle("z (cm)");
    hMap[iEps]->GetYaxis()->SetTitle("r (cm)");
    hMap[iEps]->GetZaxis()->SetTitle("#rho_{SC} (fC/cm^{3})");
    hMap[iEps]->SetTitleSize(0.05,"XYZ");
    hMap[iEps]->SetTitleOffset(1.5,"XY");
    hMap[iEps]->SetTitleOffset(0.9,"Z");
    hMap[iEps]->SetTitle(Form("%s: 50 kHz, #varepsilon = %d",sGas[iOmegaTau].Data(),eps[iEps]));
    hMap[iEps]->DrawCopy("surf2fb");
    
    
    //
    // 2. get TF1 with differences 
    //
    // get corrections (at y = 0 and z = 10) for visual correction 1 (last argument)
    // 0 = dR
    // 1 = dPhiR
    // 2 = dZ?    
    fdiffR[iEps]       = new TF1(Form("fdiffR%d",iEps), Form("AliTPCCorrection::GetDistXYZIntegrateZ(x,0,10,0,1,%f)",integrateStep),85,245);
    fdiffPhiR[iEps]    = new TF1(Form("fdiffPhiR%d",iEps), Form("AliTPCCorrection::GetDistXYZIntegrateZ(x,0,10,1,1,%f)",integrateStep),85,245);
    
    hdiffR[iEps] = (TH1F*)fdiffR[iEps]->GetHistogram();
    hdiffPhiR[iEps] = (TH1F*)fdiffPhiR[iEps]->GetHistogram();
    
    
    hdiffR[iEps]->SetName(fdiffR[iEps]->GetName());
    hdiffPhiR[iEps]->SetName(fdiffPhiR[iEps]->GetName());
    
    
    
    //
    // 3. Plot and store TH1Fs
    //
    
    hdiffR[iEps]->GetXaxis()->SetTitle("r (cm)");
    hdiffPhiR[iEps]->GetXaxis()->SetTitle("r (cm)");
    
    hdiffR[iEps]->GetYaxis()->SetTitle("dr (cm)");
    hdiffPhiR[iEps]->GetYaxis()->SetTitle("d(r#varphi) (cm)");
    
    hdiffR[iEps]->SetTitle("dr (cm)");
    hdiffPhiR[iEps]->SetTitle("d(r#varphi) (cm)");
  
    hdiffR[iEps]->SetMarkerColor(col[iEps]);
    hdiffPhiR[iEps]->SetMarkerColor(col[iEps]);
    
    hdiffR[iEps]->SetLineColor(col[iEps]);
    hdiffPhiR[iEps]->SetLineColor(col[iEps]);

    hdiffR[iEps]->SetLineStyle(sty[iEps]);
    hdiffPhiR[iEps]->SetLineStyle(sty[iEps]);
    
    hdiffR[iEps]->SetFillColor(col[iEps]);
    hdiffPhiR[iEps]->SetFillColor(col[iEps]);
    
    hdiffR[iEps]->SetTitleSize(0.05,"XYZ");
    hdiffPhiR[iEps]->SetTitleSize(0.05,"XYZ");
    
    legend->AddEntry(hdiffR[iEps],Form("#varepsilon = %d",eps[iEps]),"lp");
    
    cOmegaTau->cd(1);
    if(iEps==0)
      hdiffR[iEps]->DrawCopy("lf");
    else{
      hdiffR[iEps]->DrawCopy("lf,same");
      legend->Draw();
    }

    cOmegaTau->cd(2);
    if(iEps==0)
      hdiffPhiR[iEps]->DrawCopy("fl");
    else{
      hdiffPhiR[iEps]->DrawCopy("fl,same");
    }
  }
  
  for(Int_t iEps = 0; iEps < nEps; iEps++){
    cMap[iEps]->SaveAs(Form("%s_epsilon%d_SC_performanceNote.eps",outfilename.Data(),eps[iEps]));
    cMap[iEps]->SaveAs(Form("%s_epsilon%d_SC_performanceNote.pdf",outfilename.Data(),eps[iEps]));
  }
  cOmegaTau->SaveAs(Form("%s_Distortions_performanceNote.eps",outfilename.Data()));
  cOmegaTau->SaveAs(Form("%s_Distortions_performanceNote.pdf",outfilename.Data()));
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
