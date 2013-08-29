void SCdistortionsForTDRfinalAllGases(Double_t radiusScale=1.5, Int_t epsScale=1){
  //
  // do for one file (given by directory) the correction (specify the gas = iOmegaTau)
  // use the integrate along drift line option
  //
  // 1. Initialzation form space charge maps
  //
  AliTPCSpaceCharge3D *spaceCharge = new AliTPCSpaceCharge3D;
  const Double_t fgke0 = 8.854187817e-12; // vacuum permittivity [A·s/(V·m)]

  // omega tau parameters and TF1
  const Int_t nOmegaTau = 4;
  Double_t omegaTau[nOmegaTau] = {0.32};
  Double_t omegaTau[nOmegaTau] = {0.32,0.43,1.77,1.84};
  TString tGas[nOmegaTau] = {"NeCO2_2","ArCO2","NeCF4","NeCF4_2"}; // CF4 is the same as CO2 here, but different omegaTau
  TString sGas[nOmegaTau] = {"Ne-CO_{2}-N_{2} (90-10-5)","Ar-CO_{2} (90-10)","Ne-CF_{4} (90-10)","Ne-CF_{4} (80-20)"};
  Int_t eps[nOmegaTau] = {20,10,20,20};
  Int_t col[nOmegaTau] = {kBlack,kRed,kOrange-3,kGreen+2};
  TF1 * fdiffR[nOmegaTau];
  TF1 * fdiffPhiR[nOmegaTau];
  TH1F * hdiffR[nOmegaTau];
  TH1F * hdiffPhiR[nOmegaTau];
  TH2F * hMap[nOmegaTau];
  TLegend *legend = new TLegend(0.25,0.6,0.85,0.85,Form("#rho_{SC} ~ r^{-%.1f} for 50 kHz",radiusScale),"brNDC");
  setupLegend(legend,0.05);

  //use always the integrate option here
  Double_t integrateStep = 1.;

  TCanvas *cMap = new TCanvas("cMap","cMap",1200,900);
  cMap->Divide(2,2);

  TCanvas *cOmegaTau = new TCanvas("cOmegaTau","cOmegaTau",1200,500);
  cOmegaTau->Divide(2,1);

  TString outfilename = Form("SCdistortions_SC_NeCO2_50kHz_radiusScaling%.0f_epsScaling%d",radiusScale,epsScale);
  if(radiusScale>1.1 && radiusScale < 1.9) outfilename = Form("SCdistortions_SC_NeCO2_50kHz_radiusScaling%.1f_epsScaling%d",radiusScale,epsScale);


  //loop over gases
  for(Int_t iOmegaTau = 0; iOmegaTau < nOmegaTau; ++iOmegaTau){

    cMap->cd(iOmegaTau+1)->SetPhi(150);
    
    // select gas 
    // 0 = NeCO2N2 

    if(radiusScale>1.1 && radiusScale < 1.9){
      cout<<"Open file = "<<Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.1f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iOmegaTau],radiusScale,epsScale)<<endl;
      spaceCharge->SetSCDataFileName(Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.1f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iOmegaTau],radiusScale,epsScale));
    }
    else{
      cout<<"Open file = "<<Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.0f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iOmegaTau],radiusScale,epsScale)<<endl;
      spaceCharge->SetSCDataFileName(Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.0f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iOmegaTau],radiusScale,epsScale));
    }

  // select omegaTau value
  if(iOmegaTau ==  0){
    spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],1.00,0.99); // Ne CO2 N2 (90-10-5)
  }
  else if(iOmegaTau ==  1){
    spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],0.99,1.03); // Ar CO2 (90-10)
  }
  else if(iOmegaTau == 2){
    spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],0.41,0.70); // Ne CF4 (90-10)
  }
  else if(iOmegaTau == 3){
    spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],0.41,0.70); // Ne CF4 (80-20) (not in table use same as for other)
  }
  else{
    cout<<"wrong gas (iOmegaTau = "<<iOmegaTau<<")"<<endl;
    return;
  }
    //
    //
    // init and add to corrections
    spaceCharge->InitSpaceCharge3DDistortion();
    spaceCharge->AddVisualCorrection(spaceCharge,1);
    
    // draw the map
    
    hMap[iOmegaTau] = (TH2F*)spaceCharge->CreateHistoSCinZR(0.);
    hMap[iOmegaTau]->Scale(fgke0/1e6*1e15); // C/m^3/e0 --> fC/cm^3
    hMap[iOmegaTau]->GetXaxis()->SetTitle("z (cm)");
    hMap[iOmegaTau]->GetYaxis()->SetTitle("r (cm)");
    hMap[iOmegaTau]->GetZaxis()->SetTitle("#rho_{SC} (fC/cm^{3})");
    hMap[iOmegaTau]->SetTitleSize(0.05,"XYZ");
    hMap[iOmegaTau]->SetTitleOffset(1.5,"XY");
    hMap[iOmegaTau]->SetTitleOffset(0.9,"Z");
    hMap[iOmegaTau]->SetTitle(Form("%s: 50 kHz, #varepsilon = %d",sGas[iOmegaTau].Data(),eps[iOmegaTau]));
    hMap[iOmegaTau]->DrawCopy("surf2fb");
    
    
    //
    // 2. get TF1 with differences 
    //
    // get corrections (at y = 0 and z = 10) for visual correction 1 (last argument)
    // 0 = dR
    // 1 = dPhiR
    // 2 = dZ?    
    fdiffR[iOmegaTau]       = new TF1(Form("fdiffR%d",iOmegaTau), Form("AliTPCCorrection::GetDistXYZIntegrateZ(x,0,10,0,1,%f)",integrateStep),85,245);
    fdiffPhiR[iOmegaTau]    = new TF1(Form("fdiffPhiR%d",iOmegaTau), Form("AliTPCCorrection::GetDistXYZIntegrateZ(x,0,10,1,1,%f)",integrateStep),85,245);
    
    hdiffR[iOmegaTau] = (TH1F*)fdiffR[iOmegaTau]->GetHistogram();
    hdiffPhiR[iOmegaTau] = (TH1F*)fdiffPhiR[iOmegaTau]->GetHistogram();
    
    
    hdiffR[iOmegaTau]->SetName(fdiffR[iOmegaTau]->GetName());
    hdiffPhiR[iOmegaTau]->SetName(fdiffPhiR[iOmegaTau]->GetName());
    
    
    
    //
    // 3. Plot and store TH1Fs
    //
    
    hdiffR[iOmegaTau]->GetXaxis()->SetTitle("r (cm)");
    hdiffPhiR[iOmegaTau]->GetXaxis()->SetTitle("r (cm)");
    
    hdiffR[iOmegaTau]->GetYaxis()->SetTitle("dr (cm)");
    hdiffPhiR[iOmegaTau]->GetYaxis()->SetTitle("d(r#varphi) (cm)");
    
    hdiffR[iOmegaTau]->SetTitle("dr (cm)");
    hdiffPhiR[iOmegaTau]->SetTitle("d(r#varphi) (cm)");
  
    hdiffR[iOmegaTau]->SetMarkerColor(col[iOmegaTau]);
    hdiffPhiR[iOmegaTau]->SetMarkerColor(col[iOmegaTau]);
    
    hdiffR[iOmegaTau]->SetLineColor(col[iOmegaTau]);
    hdiffPhiR[iOmegaTau]->SetLineColor(col[iOmegaTau]);

    hdiffR[iOmegaTau]->SetLineWidth(2);
    hdiffPhiR[iOmegaTau]->SetLineWidth(2);
    
    hdiffR[iOmegaTau]->SetFillColor(col[iOmegaTau]);
    hdiffPhiR[iOmegaTau]->SetFillColor(col[iOmegaTau]);

    hdiffR[iOmegaTau]->SetMaximum(35);
    hdiffR[iOmegaTau]->SetMinimum(-25);
    hdiffPhiR[iOmegaTau]->SetMaximum(10);
    hdiffPhiR[iOmegaTau]->SetMinimum(-22);

    hdiffR[iOmegaTau]->SetTitleSize(0.05,"XYZ");
    hdiffPhiR[iOmegaTau]->SetTitleSize(0.05,"XYZ");
    
    legend->AddEntry(hdiffR[iOmegaTau],Form("%s (#varepsilon = %d)",sGas[iOmegaTau].Data(),eps[iOmegaTau]),"lp");
    
    cOmegaTau->cd(1);
    if(iOmegaTau==0)
      hdiffR[iOmegaTau]->DrawCopy("lf");
    else{
      hdiffR[iOmegaTau]->DrawCopy("lf,same");
    }
    if(iOmegaTau == nOmegaTau-1)
      legend->Draw();
    

    cOmegaTau->cd(2);
    if(iOmegaTau==0)
      hdiffPhiR[iOmegaTau]->DrawCopy("fl");
    else{
      hdiffPhiR[iOmegaTau]->DrawCopy("fl,same");
    }
  }
  
  
  cMap->SaveAs(Form("%s_SC_TDRallGases.eps",outfilename.Data()));
  cMap->SaveAs(Form("%s_SC_TDRallGases.pdf",outfilename.Data()));
  cOmegaTau->SaveAs(Form("%s_Distortions_TDRallGases.eps",outfilename.Data()));
  cOmegaTau->SaveAs(Form("%s_Distortions_TDRallGases.pdf",outfilename.Data()));
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
