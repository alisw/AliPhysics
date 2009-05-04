// hinvMass -> invariant mass histogram to fit
// typeb -> type of fit function for background. Can be 0, 1, 2, 3 (see AliHFMassFitter.cxx for details)
// types -> type of fit function for signal. Can be 0, 1 (see AliHFMassFitter.cxx for details)
// factor4refl -> sigmaRefl=factor4refl*sigmaSgn (have a look to AliHFMassFitter.cxx for details). Set it if types=1

void AliHFMassFitterTest(TH1F *hinvMass, Int_t typeb, Int_t types, Int_t factor4refl=1){

  Bool_t useParFiles=kFALSE;
  Int_t load=gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  LoadLibraries(useParFiles);

  Int_t nbin=hinvMass->GetNbinsX();
  Double_t min=hinvMass->GetBinLowEdge(1),max=hinvMass->GetBinLowEdge(nbin)+hinvMass->GetBinWidth(nbin);
  //TH1F *hMass=new TH1F("hMass","Invariant Mass",nbin,min,max);
  
  AliHFMassFitter *fitter=new AliHFMassFitter(hinvMass,min, max,1,typeb,types);
  
  fitter->SetReflectionSigmaFactor(factor4refl);
  //fitter->SetRangeFit(min,max);
  fitter->MassFitter();
  
  TNtuple *ntu=fitter->NtuParamOneShot("ntuInvMass");
  Double_t err1,err3,err6;
  Double_t significance1sigma,significance3sigma,significance6sigma;
  fitter->Significance(1,significance1sigma,err1);
  fitter->Significance(3,significance3sigma,err3);
  fitter->Significance(6,significance6sigma,err6);

  cout<<"N sigma \t significance\n";
  cout<<"1\t"<<significance1sigma<<endl;
  cout<<"3\t"<<significance3sigma<<endl;
  cout<<"6\t"<<significance6sigma<<endl;

}
