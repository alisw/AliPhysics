/**********************************************
 *
 *   Class meant to do Glauber+NBD fits
 *
 *   This class makes full use of analytical
 *   properties of the Negative Binomial Function
 *
 *   Only the Glauber component is MC, while
 *   the NBD is evaluated probabilistically
 *
 *  - bugs, comments, suggestions, complaints?
 *  - Feel free to write to:
 *     david.dobrigkeit.chinellato@cern.ch
 *
 **********************************************/

#include "AliMultGlauberNBDFitter.h"
#include "TList.h"
#include "TFile.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TVirtualFitter.h"

ClassImp(AliMultGlauberNBDFitter);

AliMultGlauberNBDFitter::AliMultGlauberNBDFitter() : TNamed(), 
fNBD(0x0),
fhNanc(0x0),
fhNpNc(0x0),
ffChanged(kTRUE),
fCurrentf(-1),
fNpart(0x0),
fNcoll(0x0),
fContent(0x0),
fNNpNcPairs(-1),
fMaxNpNcPairs(1000000),
fMu(45),
fk(1.5),
ff(0.8),
fnorm(100),
fFitOptions("R0")
{
  // Constructor
  fNpart = new Double_t[fMaxNpNcPairs];
  fNcoll = new Double_t[fMaxNpNcPairs];
  fContent = new Long_t[fMaxNpNcPairs];
  
  //Ancestor histo
  fhNanc = new TH1D("fhNanc", "", 1000, -0.5, 999.5);
  
  //NBD
  fNBD = new TF1("fNBD","ROOT::Math::negative_binomial_pdf(x,[0],[1])",0,800);
  
  //master function
  fGlauberNBD = new TF1("fGlauberNBD", this, &AliMultGlauberNBDFitter::ProbDistrib,
                        0, 50000, 4 , "AliMultGlauberNBDFitter", "ProbDistrib");
  fGlauberNBD->SetParameter(0,fMu);
  fGlauberNBD->SetParameter(1,fk);
  fGlauberNBD->SetParameter(2,ff);
  fGlauberNBD->SetParameter(3,fnorm);
}

AliMultGlauberNBDFitter::AliMultGlauberNBDFitter(const char * name, const char * title): TNamed(name,title),
fNBD(0x0),
fhNanc(0x0),
fhNpNc(0x0),
ffChanged(kTRUE),
fCurrentf(-1),
fNpart(0x0),
fNcoll(0x0),
fContent(0x0),
fNNpNcPairs(-1),
fMaxNpNcPairs(1000000),
fMu(45),
fk(1.5),
ff(0.8),
fnorm(100),
fFitOptions("R0")
{
  //Named constructor
  fNpart = new Double_t[fMaxNpNcPairs];
  fNcoll = new Double_t[fMaxNpNcPairs];
  fContent = new Long_t[fMaxNpNcPairs];
  
  //Ancestor histo
  fhNanc = new TH1D("fhNanc", "", 1000, -0.5, 999.5);
  
  //NBD
  fNBD = new TF1("fNBD","ROOT::Math::negative_binomial_pdf(x,[0],[1])",0,800);
  
  //master function
  fGlauberNBD = new TF1("fGlauberNBD", this, &AliMultGlauberNBDFitter::ProbDistrib,
                        0, 50000, 4 , "AliMultGlauberNBDFitter", "ProbDistrib");
  fGlauberNBD->SetParameter(0,fMu);
  fGlauberNBD->SetParameter(1,fk);
  fGlauberNBD->SetParameter(2,ff);
  fGlauberNBD->SetParameter(3,fnorm);
}
//________________________________________________________________
AliMultGlauberNBDFitter::~AliMultGlauberNBDFitter() {
  // Destructor
  if (fNBD) {
    delete fNBD;
    fNBD = 0x0;
  }
  if (fhNanc) {
    delete fhNanc;
    fhNanc = 0x0;
  }
  if (fhNpNc) {
    delete fhNpNc;
    fhNpNc = 0x0;
  }
  if (fNpart) delete [] fNpart;
  if (fNcoll) delete [] fNcoll;
  if (fContent) delete [] fContent;
}

//______________________________________________________
Double_t AliMultGlauberNBDFitter::ProbDistrib(Double_t *x, Double_t *par)
//Master fitter function
{
  Double_t lMultValue = TMath::Floor(x[0]+0.5);
  Double_t lProbability = 0.0;
  ffChanged = kTRUE;

  //Comment this line in order to make the code evaluate Nancestor all the time
  if ( TMath::Abs( fCurrentf - par[2] ) < kAlmost0 ) ffChanged = kFALSE ;
  
  //______________________________________________________
  //Recalculate the ancestor distribution in case f changed
  if( ffChanged ){
    fCurrentf = par[2];
    fhNanc->Reset();
    
    for(int ibin=0;ibin<fNNpNcPairs;ibin++){
      //Atentar-se à normalização de Nanc
      fhNanc->Fill(TMath::Floor(fNpart[ibin]*par[2] + fNcoll[ibin]*(1-par[2]) + 0.5),fContent[ibin]);
    }
    fhNanc->Scale(1./fhNanc->Integral());
  }
  //______________________________________________________
  //Actually ealuate function
  for(Long_t iNanc = 1; iNanc<900; iNanc++){
    Double_t lThisMu = ((Double_t)iNanc)*par[0];
    Double_t lThisk = ((Double_t)iNanc)*par[1];
    Double_t lpval = TMath::Power(1+lThisMu/lThisk,-1);
    fNBD->SetParameter(1,lThisk);
    fNBD->SetParameter(0,lpval);
    Double_t lMult = fNBD->Eval(lMultValue);
    lProbability += fhNanc->GetBinContent(fhNanc->FindBin(iNanc))*lMult;
  }
  //______________________________________________________
  return par[3]*lProbability;
}

//________________________________________________________________
Bool_t AliMultGlauberNBDFitter::SetNpartNcollCorrelation(TH2 *hNpNc){
  Bool_t lReturnValue = kTRUE;
  if( hNpNc ){
    fhNpNc = (TH2*) hNpNc;
  }else{
    lReturnValue = kFALSE;
  }
  return lReturnValue;
}

//________________________________________________________________
Bool_t AliMultGlauberNBDFitter::SetInputV0M(TH1 *hV0M){
  Bool_t lReturnValue = kTRUE;
  if( hV0M ){
    fhV0M = (TH1*) hV0M;
  }else{
    lReturnValue = kFALSE;
  }
  return lReturnValue;
}

//________________________________________________________________
TF1 *AliMultGlauberNBDFitter::GetNBD(){
  return fNBD;
}

//________________________________________________________________
TF1 *AliMultGlauberNBDFitter::GetGlauberNBD(){
  return fGlauberNBD;
}

//________________________________________________________________
void AliMultGlauberNBDFitter::SetFitRange( Double_t lMin, Double_t lMax){
  fGlauberNBD -> SetRange(lMin, lMax);
}

//________________________________________________________________
void AliMultGlauberNBDFitter::SetFitOptions(TString lOpt){
  fFitOptions = lOpt;
}

//________________________________________________________________
Bool_t AliMultGlauberNBDFitter::DoFit(){
  //Try very hard, please
  TVirtualFitter::SetMaxIterations(5000000);
  if( !InitializeNpNc() ) return kFALSE ;
  
  TStopwatch* timer = new TStopwatch();
  timer->Start ( kTRUE );
  cout<<"---> Now fitting, please wait..."<<endl;
  
  fGlauberNBD->SetNpx(100);
  fhV0M->Fit("fGlauberNBD",fFitOptions.Data());
  
  timer->Stop();
  Double_t lTotalTime = timer->RealTime();
  cout<<"---> Fitting took "<<lTotalTime<<" seconds"<<endl;
  
  fMu   = fGlauberNBD -> GetParameter(0);
  fk    = fGlauberNBD -> GetParameter(1);
  ff    = fGlauberNBD -> GetParameter(2);
  fnorm = fGlauberNBD -> GetParameter(3);
  
  return kTRUE;
}

//________________________________________________________________
Bool_t AliMultGlauberNBDFitter::InitializeNpNc(){
  //This function initializes fhNpNc
  //Warning: X == Npart, Y == Ncoll
  Bool_t lReturnValue = kFALSE;
  if(fhNpNc){
    fNNpNcPairs = 0;
    //Sweep all allowed values of Npart, Ncoll; find counters
    for(int xbin=1;xbin<500;xbin++){
      for(int ybin=1;ybin<3000;ybin++){
        if(fhNpNc->GetBinContent(fhNpNc->FindBin(xbin,ybin)) != 0){
          fNpart[fNNpNcPairs] = xbin;
          fNcoll[fNNpNcPairs] = ybin;
          fContent[fNNpNcPairs] = fhNpNc->GetBinContent(fhNpNc->FindBin(xbin,ybin));
          fNNpNcPairs++;
        }
      }
    }
    cout<<"Initialized with number of (Npart, Ncoll) pairs: "<<fNNpNcPairs<<endl;
    lReturnValue = kTRUE;
  }else{
    cout<<"Failed to initialize! Please provide input histogram with (Npart, Ncoll) info!"<<endl;
    cout<<"Please remember to call SetNpartNcollCorrelation before doing fit!"<<endl;
  }
  return lReturnValue;
}
