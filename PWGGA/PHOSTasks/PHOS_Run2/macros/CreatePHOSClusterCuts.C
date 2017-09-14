AliPHOSClusterCuts *CreatePHOSClusterCuts(Bool_t useCoreDisp, Double_t NsigmaCPV, Double_t NsigmaDisp)
{
  //floating point can not be used as a list name. index goes to listname.
  TString index="";
  if(NsigmaCPV > 0) index += Form("_CPV%d",(Int_t)(NsigmaCPV*10));
  if(NsigmaDisp > 0){
    if(useCoreDisp) index += Form("_CoreDisp%d",(Int_t)(NsigmaDisp*10));
    else            index += Form("_FullDisp%d",(Int_t)(NsigmaDisp*10));
  }

  AliPHOSClusterCuts *cuts = new AliPHOSClusterCuts(index.Data());
  cuts->SetUseCoreDispersion(useCoreDisp);
  cuts->SetNsigmaCPV(NsigmaCPV);
  cuts->SetNsigmaDisp(NsigmaDisp);

  return cuts;
};

