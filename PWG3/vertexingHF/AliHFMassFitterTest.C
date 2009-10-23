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
  fitter->MassFitter(kTRUE); //kFALSE do not draw the histogram
  
  TNtuple *ntu=fitter->NtuParamOneShot("ntuInvMass");

  Double_t mD0pdg=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t sigmaval= 0.027;
  Double_t min=mD0pdg-sigmaval, max=mD0pdg+sigmaval;
  Double_t significance3sigma,significanceRange,signal3sigma,signalRange,bkg3sigma,bkgRange;
  Double_t err3,errRange,errS3,errSRange,errB3,errBRange;

  fitter->Signal(3,signal3sigma,errS3);
  fitter->Signal(min, max,signalRange,errSRange);
  fitter->Background(3,bkg3sigma,errB3);
  fitter->Background(min, max,bkgRange,errBRange);
  fitter->Significance(3,significance3sigma,err3);
  fitter->Significance(min, max,significanceRange,errRange);

  cout<<"# 3 sigma: \n signal = "<<signal3sigma<<" +/- "<<errS3<<"\t backgdround = "<<bkg3sigma<<" +/- "<<errB3<<"\tsignificance = "<<significance3sigma<<" +/- "<<err3<<"\n";

  cout<<"# range ("<<min<<", "<<max<<"): \n signal = "<<signalRange<<" +/- "<<errSRange<<"\t backgdround = "<<bkgRange<<" +/- "<<errBRange<<"\tsignificance = "<<significanceRange<<" +/- "<<errRange<<"\n";


}
