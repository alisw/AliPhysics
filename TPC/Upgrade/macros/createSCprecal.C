void createSCprecal(TString input, Int_t gas=1)
{
  AliTPCSpaceCharge3D *spaceCharge = new AliTPCSpaceCharge3D;
  spaceCharge->SetSCDataFileName(input.Data());
  
  const Int_t nOmegaTau = 5;
  Double_t omegaTau[nOmegaTau] = {0.34, 0.32, 0.43, 1.77, 1.84};
  Double_t T1[nOmegaTau]       = {1.00, 1.00, 0.99, 0.41, 0.41};
  Double_t T2[nOmegaTau]       = {1.01, 0.99, 1.03, 0.70, 0.70};
  
  TString tGas[nOmegaTau] = {"NeCO2","NeCO2_2","ArCO2","NeCF4","NeCF4_2"}; // CF4 is the same as CO2 here, but different omegaTau
  TString sGas[nOmegaTau] = {"Ne-CO_{2} (90-10)","Ne-CO_{2}-N_{2} (90-10-5)","Ar-CO_{2} (90-10)","Ne-CF_{4} (90-10)","Ne-CF_{4} (80-20)"};

  spaceCharge->SetOmegaTauT1T2(omegaTau[gas], T1[gas] , T2[gas]);
  spaceCharge->InitSpaceCharge3DDistortion();

  TString outName=input;
  outName.ReplaceAll(".root","_precal.root");

  TFile fout(outName,"recreate");
  spaceCharge->Write("map");
  fout.Write();
  fout.Close();

  delete spaceCharge;
}
