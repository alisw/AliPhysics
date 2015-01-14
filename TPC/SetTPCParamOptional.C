/// \file SetTPCParamOptional.C
/// \brief Optional set of parameters

AliTPCParam *SetTPCParamOptional(){

  AliTPCParam *par=new AliTPCParamSR;
  par->SetTitle("optional"); //This is a dummy code ! Just to change something
  return par;
}
