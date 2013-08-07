void SCdistortionsForTDRfinalAll(){

  SCdistortionsForTDRfinalsingle(1,1,0);
  SCdistortionsForTDRfinalsingle(1,1,1);
  SCdistortionsForTDRfinalsingle(1,1,2);
  SCdistortionsForTDRfinalsingle(1,1,3);
  SCdistortionsForTDRfinalsingle(1,1,4);

}

void SCdistortionsForTDRfinal(Double_t radiusScale=1.5, Int_t epsScale=1,Int_t iOmegaTau = 1){
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
  Int_t col[nEps] = {kBlack,kBlack};

  const Int_t nOmegaTau = 5;
  Double_t omegaTau[nOmegaTau] = {0.34,0.32,0.43,1.77,1.84};
  TString tGas[nOmegaTau] = {"NeCO2","NeCO2_2","ArCO2","NeCF4","NeCF4_2"}; // CF4 is the same as CO2 here, but different omegaTau
  TString sGas[nOmegaTau] = {"Ne-CO_{2} (90-10)","Ne-CO_{2}-N_{2} (90-10-5)","Ar-CO_{2} (90-10)","Ne-CF_{4} (90-10)","Ne-CF_{4} (80-20)"};
  TF2 * fdiffR[nEps];
  TF2 * fdiffPhiR[nEps];
  TH2F * hdiffR[nEps];
  TH2F * hdiffPhiR[nEps];
  TF2 * fdiffIntR[nEps];
  TF2 * fdiffIntPhiR[nEps];
  TF2 * fdiffIntZ[nEps];
  TH2F * hdiffIntR[nEps];
  TH2F * hdiffIntPhiR[nEps];
  TH2F * hdiffIntZ[nEps];
  TH2F * hMap[nEps];
  TH2F * hDistRMap[nEps];
  TH2F * hDistRPMap[nEps];

  //use always the integrate option here
  Double_t integrateStep = 1.;

  TCanvas *cMap = new TCanvas("cMap","cMap",1200,500);
  cMap->Divide(2,1);

  TCanvas *cDistRMap = new TCanvas("cDistRMap","cDistRMap",1200,500);
  cDistRMap->Divide(2,1);

  TCanvas *cDistRPMap = new TCanvas("cDistRPMap","cDistRPMap",1200,500);
  cDistRPMap->Divide(2,1);

  TCanvas *cDistRNonIntMap = new TCanvas("cDistRNonIntMap","cDistRNonIntMap",1200,500);
  cDistRNonIntMap->Divide(2,1);

  TCanvas *cDistRPNonIntMap = new TCanvas("cDistRPNonIntMap","cDistRPNonIntMap",1200,500);
  cDistRPNonIntMap->Divide(2,1);

  TCanvas *cDistRIntMap = new TCanvas("cDistRIntMap","cDistRIntMap",1200,500);
  cDistRIntMap->Divide(2,1);

  TCanvas *cDistRPIntMap = new TCanvas("cDistRPIntMap","cDistRPIntMap",1200,500);
  cDistRPIntMap->Divide(2,1);

  TCanvas *cDistZIntMap = new TCanvas("cDistZIntMap","cDistZIntMap",1200,500);
  cDistZIntMap->Divide(2,1);



  TString outfilename = Form("SCdistortions_%s_50kHz_radiusScaling%.0f_epsScaling%d",tGas[iOmegaTau].Data(),radiusScale,epsScale);
  if(radiusScale>1.1 && radiusScale < 1.9) outfilename = Form("SCdistortions_%s_50kHz_radiusScaling%.1f_epsScaling%d",tGas[iOmegaTau].Data(),radiusScale,epsScale);

  //loop over epsilons
  for(Int_t iEps = 0; iEps < nEps; ++iEps){

    cMap->cd(iEps+1);

  // select gas 
  // 0 = NeCO2 
  // 1 = NeCO2N2 
  // 2 = ArCO2 
  // 3 = NeCF4 
  // 4 = ArCF4 
     if(radiusScale>1.1 && radiusScale < 1.9){
      cout<<"Open file = "<<Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.1f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iEps],radiusScale,epsScale)<<endl;
      spaceCharge->SetSCDataFileName(Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.1f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iEps],radiusScale,epsScale));
    }
    else{
      cout<<"Open file = "<<Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.0f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iEps],radiusScale,epsScale)<<endl;
      spaceCharge->SetSCDataFileName(Form("/Users/physics/ALICE/TPCupgrade/SpaceCharge/Maps/SC_%s_eps%d_50kHz_radiusScaling%.0f_epsScaling%d/SpaceChargeMap.root",tGas[iOmegaTau].Data(),eps[iEps],radiusScale,epsScale));
    }

  // select omegaTau value
  if(iOmegaTau ==  1){
    spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],1.00,0.99); // Ne CO2 N2 (90-10-5)
  }
  if(iOmegaTau ==  2){
    spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],0.99,1.03); // Ar CO2 (90-10)
  }
  else if(iOmegaTau == 3){
    spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],0.41,0.70); // Ne CF4 (90-10)
  }
  else if(iOmegaTau == 4){
    spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],0.41,0.70); // Ne CF4 (80-20) (not in table use same as for other)
  }
  else{
    spaceCharge->SetOmegaTauT1T2(omegaTau[iOmegaTau],1,1.01); // Ne CO2 (90-10)
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
  hMap[iEps]->GetZaxis()->SetTitle("");
  //hMap[iEps]->GetZaxis()->SetTitle("#rho_{SC} (fC/cm^{3})");
  hMap[iEps]->SetTitleSize(0.05,"XYZ");
  //hMap[iEps]->SetTitleOffset(1.5,"XY");
  //hMap[iEps]->SetTitleOffset(0.9,"Z");
  hMap[iEps]->SetTitle(Form("#rho_{SC} (fC/cm^{3}) for %s, 50 kHz, #varepsilon = %d",sGas[iOmegaTau].Data(),eps[iEps]));
  hMap[iEps]->DrawCopy("colz");

  // draw the distortion maps
  cDistRMap->cd(iEps+1);
  hDistRMap[iEps] = (TH2F*)spaceCharge->CreateHistoDRinZR(0.);
  hDistRMap[iEps]->SetMaximum(25);
  hDistRMap[iEps]->SetMinimum(-12);
  hDistRMap[iEps]->GetXaxis()->SetTitle("z (cm)");
  hDistRMap[iEps]->GetYaxis()->SetTitle("r (cm)");
  hDistRMap[iEps]->GetZaxis()->SetTitle("");
  hDistRMap[iEps]->SetTitleSize(0.05,"XYZ");
  //hDistRMap[iEps]->SetTitleOffset(1.5,"XY");
  //hDistRMap[iEps]->SetTitleOffset(0.9,"Z");
  hDistRMap[iEps]->SetTitle(Form("dr (cm) for %s, 50 kHz, #varepsilon = %d",sGas[iOmegaTau].Data(),eps[iEps]));
  hDistRMap[iEps]->DrawCopy("colz");

  cDistRPMap->cd(iEps+1);
  hDistRPMap[iEps] = (TH2F*)spaceCharge->CreateHistoDRinZR(0.);
  hDistRPMap[iEps]->SetMaximum(25);
  hDistRPMap[iEps]->SetMinimum(-12);
  hDistRPMap[iEps]->GetXaxis()->SetTitle("z (cm)");
  hDistRPMap[iEps]->GetYaxis()->SetTitle("r (cm)");
  hDistRPMap[iEps]->GetZaxis()->SetTitle("");
  hDistRPMap[iEps]->SetTitleSize(0.05,"XYZ");
  //hDistRPMap[iEps]->SetTitleOffset(1.5,"XY");
  //hDistRPMap[iEps]->SetTitleOffset(0.9,"Z");
  hDistRPMap[iEps]->SetTitle(Form("d(r#varphi) (cm) for %s, 50 kHz, #varepsilon = %d",sGas[iOmegaTau].Data(),eps[iEps]));
  hDistRPMap[iEps]->DrawCopy("colz");

  
  //
  // 2. get TF2 with differences 
  //
  // get corrections (at y = 0) for visual correction 1 (last argument)
  // 0 = dR
  // 1 = dPhiR
  // 2 = dZ?    
  fdiffR[iEps]       = new TF2(Form("fdiffR%d",iEps), Form("AliTPCCorrection::GetDistXYZ(y,0,x,0,1)"),-250,250,85,250);
  fdiffPhiR[iEps]    = new TF2(Form("fdiffPhiR%d",iEps), Form("AliTPCCorrection::GetDistXYZ(y,0,x,1,1)"),-250,250,85,250);
  hdiffR[iEps] = (TH2F*)fdiffR[iEps]->GetHistogram();
  hdiffPhiR[iEps] = (TH2F*)fdiffPhiR[iEps]->GetHistogram();

  hdiffR[iEps]->SetName(fdiffR[iEps]->GetName());
  hdiffPhiR[iEps]->SetName(fdiffPhiR[iEps]->GetName());

  fdiffIntR[iEps]       = new TF2(Form("fdiffIntR%d",iEps), Form("AliTPCCorrection::GetDistXYZIntegrateZ(y,0,x,0,1,%f)",integrateStep),-250,250,85,250);
  fdiffIntPhiR[iEps]    = new TF2(Form("fdiffIntPhiR%d",iEps), Form("AliTPCCorrection::GetDistXYZIntegrateZ(y,0,x,1,1,%f)",integrateStep),-250,250,85,250);
  fdiffIntZ[iEps]    = new TF2(Form("fdiffIntZ%d",iEps), Form("AliTPCCorrection::GetDistXYZIntegrateZ(y,0,x,2,1,%f)",integrateStep),-250,250,85,250);
  hdiffIntR[iEps] = (TH2F*)fdiffIntR[iEps]->GetHistogram();
  hdiffIntPhiR[iEps] = (TH2F*)fdiffIntPhiR[iEps]->GetHistogram();
  hdiffIntZ[iEps] = (TH2F*)fdiffIntZ[iEps]->GetHistogram();

  hdiffIntR[iEps]->SetName(fdiffIntR[iEps]->GetName());
  hdiffIntPhiR[iEps]->SetName(fdiffIntPhiR[iEps]->GetName());
  hdiffIntZ[iEps]->SetName(fdiffIntZ[iEps]->GetName());

  

  //
  // 3. Plot and store TH1Fs
  //
  
  // draw the distortion maps
  cDistRNonIntMap->cd(iEps+1);
  hdiffR[iEps]->SetMaximum(25);
  hdiffR[iEps]->SetMinimum(-12);
  hdiffR[iEps]->GetXaxis()->SetTitle("z (cm)");
  hdiffR[iEps]->GetYaxis()->SetTitle("r (cm)");
  //hdiffR[iEps]->GetZaxis()->SetTitle("dr (cm)");
  hdiffR[iEps]->SetTitleSize(0.05,"XYZ");
  //hdiffR[iEps]->SetTitleOffset(1.5,"XY");
  //hdiffR[iEps]->SetTitleOffset(0.9,"Z");
  hdiffR[iEps]->SetTitle(Form("dr (cm) for %s, 50 kHz, #varepsilon = %d",sGas[iOmegaTau].Data(),eps[iEps]));
  hdiffR[iEps]->DrawCopy("colz");


  cDistRPNonIntMap->cd(iEps+1);
  hdiffPhiR[iEps]->SetMaximum(25);
  hdiffPhiR[iEps]->SetMinimum(-12);
  hdiffPhiR[iEps]->GetXaxis()->SetTitle("z (cm)");
  hdiffPhiR[iEps]->GetYaxis()->SetTitle("r (cm)");
  hdiffPhiR[iEps]->SetTitleSize(0.05,"XYZ");
  //hdiffPhiR[iEps]->SetTitleOffset(1.5,"XY");
  //hdiffPhiR[iEps]->SetTitleOffset(0.9,"Z");
  hdiffPhiR[iEps]->SetTitle(Form("d(r#varphi) (cm) for %s, 50 kHz, #varepsilon = %d",sGas[iOmegaTau].Data(),eps[iEps]));
  hdiffPhiR[iEps]->DrawCopy("colz");

  cDistRIntMap->cd(iEps+1);
  hdiffIntR[iEps]->SetMaximum(25);
  hdiffIntR[iEps]->SetMinimum(-12);
  hdiffIntR[iEps]->GetXaxis()->SetTitle("z (cm)");
  hdiffIntR[iEps]->GetYaxis()->SetTitle("r (cm)");
  //hdiffIntR[iEps]->GetZaxis()->SetTitle("dr (cm)");
  hdiffIntR[iEps]->SetTitleSize(0.05,"XYZ");
  //hdiffIntR[iEps]->SetTitleOffset(1.5,"XY");
  //hdiffIntR[iEps]->SetTitleOffset(0.9,"Z");
  hdiffIntR[iEps]->SetTitle(Form("dr (cm) for %s, 50 kHz, #varepsilon = %d",sGas[iOmegaTau].Data(),eps[iEps]));
  hdiffIntR[iEps]->DrawCopy("colz");


  cDistRPIntMap->cd(iEps+1);
  hdiffIntPhiR[iEps]->SetMaximum(5);
  hdiffIntPhiR[iEps]->SetMinimum(-9);
  hdiffIntPhiR[iEps]->GetXaxis()->SetTitle("z (cm)");
  hdiffIntPhiR[iEps]->GetYaxis()->SetTitle("r (cm)");
  hdiffIntPhiR[iEps]->SetTitleSize(0.05,"XYZ");
  //hdiffIntPhiR[iEps]->SetTitleOffset(1.5,"XY");
  //hdiffIntPhiR[iEps]->SetTitleOffset(0.9,"Z");
  hdiffIntPhiR[iEps]->SetTitle(Form("d(r#varphi) (cm) for %s, 50 kHz, #varepsilon = %d",sGas[iOmegaTau].Data(),eps[iEps]));
  hdiffIntPhiR[iEps]->DrawCopy("colz");

  cDistZIntMap->cd(iEps+1);
  hdiffIntZ[iEps]->SetMaximum(5);
  hdiffIntZ[iEps]->SetMinimum(-5);
  hdiffIntZ[iEps]->GetXaxis()->SetTitle("z (cm)");
  hdiffIntZ[iEps]->GetYaxis()->SetTitle("r (cm)");
  hdiffIntZ[iEps]->SetTitleSize(0.05,"XYZ");
  //hdiffIntZ[iEps]->SetTitleOffset(1.5,"XY");
  //hdiffIntZ[iEps]->SetTitleOffset(0.9,"Z");
  hdiffIntZ[iEps]->SetTitle(Form("dz (cm) for %s, 50 kHz, #varepsilon = %d",sGas[iOmegaTau].Data(),eps[iEps]));
  hdiffIntZ[iEps]->DrawCopy("colz");
  }
  
  // just control histograms (not saved)
  // cDistRNonIntMap->SaveAs(Form("%s_TDR_DistortR.eps",outfilename.Data()));
  // cDistRNonIntMap->SaveAs(Form("%s_TDR_DistortR.png",outfilename.Data()));
  // cDistRNonIntMap->SaveAs(Form("%s_TDR_DistortR.pdf",outfilename.Data()));
  // cDistRPNonIntMap->SaveAs(Form("%s_TDR_DistortRPhi.eps",outfilename.Data()));
  // cDistRPNonIntMap->SaveAs(Form("%s_TDR_DistortRPhi.png",outfilename.Data()));
  // cDistRPNonIntMap->SaveAs(Form("%s_TDR_DistortRPhi.pdf",outfilename.Data()));

  // histograms for TDR (saved)
  cMap->SaveAs(Form("%s_TDR_SpaceCharge.eps",outfilename.Data()));
  cMap->SaveAs(Form("%s_TDR_SpaceCharge.png",outfilename.Data()));
  cMap->SaveAs(Form("%s_TDR_SpaceCharge.pdf",outfilename.Data()));
  cDistRIntMap->SaveAs(Form("%s_TDR_DistortIntR.eps",outfilename.Data()));
  cDistRIntMap->SaveAs(Form("%s_TDR_DistortIntR.png",outfilename.Data()));
  cDistRIntMap->SaveAs(Form("%s_TDR_DistortIntR.pdf",outfilename.Data()));
  cDistRPIntMap->SaveAs(Form("%s_TDR_DistortIntRPhi.eps",outfilename.Data()));
  cDistRPIntMap->SaveAs(Form("%s_TDR_DistortIntRPhi.png",outfilename.Data()));
  cDistRPIntMap->SaveAs(Form("%s_TDR_DistortIntRPhi.pdf",outfilename.Data()));
  cDistZIntMap->SaveAs(Form("%s_TDR_DistortIntZ.eps",outfilename.Data()));
  cDistZIntMap->SaveAs(Form("%s_TDR_DistortIntZ.png",outfilename.Data()));
  cDistZIntMap->SaveAs(Form("%s_TDR_DistortIntZ.pdf",outfilename.Data()));
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
