AliTPCParam * SetTPCParam()
{
  TDirectory *dirsave=gDirectory;

  AliTPCParamSR  * param = new AliTPCParamSR;
  AliTPCPRF2D    * prfinner   = new AliTPCPRF2D;
  AliTPCPRF2D    * prfouter   = new AliTPCPRF2D;
  AliTPCRF1D     * rf    = new AliTPCRF1D(kTRUE);
  param->SetTitle("75x40_100x60");


  param->SetSectorAngles(20.,0.,20.,0.);
  param->SetInnerRadiusLow(87.35);
  param->SetInnerRadiusUp(130.6);
  param->SetOuterRadiusLow(131.2);
  param->SetOuterRadiusUp(252.2);
  param->SetInnerPadPitchLength(0.75);
  param->SetInnerPadPitchWidth(0.40);
  param->SetOuterPadPitchLength(1.00);
  param->SetOuterPadPitchWidth(0.60);
  param->SetInnerNWires(3);
  param->SetOuterNWires(4);

  param->SetZeroSup(2); //3 is included !
  param->SetDriftV(2.83e6);
  param->SetDiffT(0.022);
  param->SetDiffL(0.022);
  param->SetNoise(1000);
  param->SetGasGain(2.e4);
  param->SetTFWHM(1.9e-7);
    param->SetTSample(2.0e-7);
    param->SetMaxTBin(445);
  param->SetChipGain(12);      
  param->SetChipNorm(0.4);
  param->SetNCrossRows(1);
  param->SetFacSigmaPadRow(3.);
  param->SetFacSigmaPad(3.);
  param->SetFacSigmaTime(3.);
  param->Update();
  //Set z (time) response function
  rf->SetGauss(param.GetZSigma(),param.GetZWidth(),1.);
  rf->SetOffset(3*param.GetZSigma());
  rf->Update();
  //Set two dimensional pad response function
  TFile f("$ALICE_ROOT/TPC/AliTPCprf2d.root");
  prfinner->Read("prf_07504_Gati_056068_d02");
  prfouter->Read("prf_10006_Gati_047051_d03");
  param->SetInnerPRF(prfinner); //param object is responsible for destroying objects
  param->SetOuterPRF(prfouter); 
  param->SetTimeRF(rf);
  f.Close();
  
  //gTPCParam =param;
  dirsave->cd();
  return param;
};
