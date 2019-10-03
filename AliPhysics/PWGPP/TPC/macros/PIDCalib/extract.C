void extract(const char* input, const char* splines, const char* DataOrMc, const char* period, const char* pass, const char* collsys){
  //gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGPP/TPC/AliTPCcalibResidualPID.cxx+g");
  TFile f(input);
  hsp2=(THnSparseF*)TPCresPID->FindObject("fHistPidQA");
  //hsp2->GetAxis(8)->SetRangeUser(-0.2,0.2);
  gROOT->cd();
  Double_t parametersBB[6] = { 0. }; 
  parametersBB[0] = 59;
  parametersBB[1] = 19;
  parametersBB[2] = 15;
  parametersBB[3] = 0.8;
  parametersBB[4] = 73;
  parametersBB[5] = -3;
  
  //2.6,14.3,-15,2.2,2.7, 0.06};
  AliTPCcalibResidualPID::ExtractResidualPID(hsp2,kTRUE,splines,DataOrMc,period,pass,collsys,0x0/*or parametersBB to set initial values manually*/,"", AliTPCcalibResidualPID::kAleph);
}

