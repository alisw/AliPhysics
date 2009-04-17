/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////
//
// AliHFMassFitter for the fit of invariant mass distribution
// of charmed mesons
//
// Author: C.Bianchin, chiara.bianchin@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TCanvas.h>

#include "AliHFMassFitter.h"


ClassImp(AliHFMassFitter)

  //constructors
AliHFMassFitter::AliHFMassFitter() : TNamed(),
  fhistoInvMass(0),
  fminMass(0),
  fmaxMass(0),
  fNbin(1),
  fWithBkg(0),
  ftypeOfFit4Bkg(0),
  ftypeOfFit4Sgn(0),
  ffactor(1),
  fntuParam(0),
  fMass(1.85),
  fSigmaSgn(0.012),
  fSideBands(0)
{
  // default constructor

  cout<<"Default constructor"<<endl;
}

//___________________________________________________________________________

AliHFMassFitter::AliHFMassFitter (TH1F *histoToFit, Double_t minvalue, Double_t maxvalue, Int_t rebin,Int_t fittypeb,Int_t fittypes): 
 TNamed(),
 fhistoInvMass(0),
 fminMass(0),
 fmaxMass(0),
 fNbin(1),
 fWithBkg(0),
 ftypeOfFit4Bkg(0),
 ftypeOfFit4Sgn(0),
 ffactor(1),
 fntuParam(0),
 fMass(1.85),
 fSigmaSgn(0.012),
 fSideBands(0)
{
  // standard constructor

  fhistoInvMass=histoToFit;
  fminMass=minvalue; 
  fmaxMass=maxvalue;
  if(rebin!=1) RebinMass(rebin); 
  else fNbin=(Int_t)fhistoInvMass->GetNbinsX();
  ftypeOfFit4Bkg=fittypeb;
  ftypeOfFit4Sgn=fittypes;
  if(ftypeOfFit4Bkg!=0 && ftypeOfFit4Bkg!=1 && ftypeOfFit4Bkg!=2) fWithBkg=kFALSE;
  else fWithBkg=kTRUE;
  if (!fWithBkg) cout<<"Fit Histogram of Signal only"<<endl;
  else  cout<<"Type of fit For Background = "<<ftypeOfFit4Bkg<<endl;
}

//___________________________________________________________________________

Double_t AliHFMassFitter::FitFunction4MassDistr (Double_t *x, Double_t *par){
  // Fit function for signal+background


  //exponential or linear fit
  //
  // par[0] = tot integral
  // par[1] = slope
  // par[2] = gaussian integral
  // par[3] = gaussian mean
  // par[4] = gaussian sigma
  
  Double_t total,bkg=0,sgn=0;
  
  if (ftypeOfFit4Bkg==0 || ftypeOfFit4Bkg==1) {
    if(ftypeOfFit4Sgn == 0) {

      Double_t parbkg[2] = {par[0]-par[2], par[1]};
      bkg = FitFunction4Bkg(x,parbkg);
    }
    if(ftypeOfFit4Sgn == 1) {
      Double_t parbkg[5] = {par[2],par[3],ffactor*par[4],par[0]-2*par[2], par[1]};
      bkg = FitFunction4Bkg(x,parbkg);
    }

    sgn = FitFunction4Sgn(x,&par[2]);  

  }

  //polynomial fit

    // par[0] = tot integral
    // par[1] = coef1
    // par[2] = coef2
    // par[3] = gaussian integral
    // par[4] = gaussian mean
    // par[5] = gaussian sigma

  if (ftypeOfFit4Bkg==2) {
    
    if(ftypeOfFit4Sgn == 0) {
      //parbkg = new Double_t[2];
      Double_t parbkg[3] = {par[0]-par[3], par[1], par[2]};
      bkg = FitFunction4Bkg(x,parbkg);
    }
    if(ftypeOfFit4Sgn == 1) {
      //parbkg = new Double_t[6];
      Double_t parbkg[6] = {par[3],par[4],ffactor*par[5],par[0]-2*par[3], par[1], par[2]};
      bkg = FitFunction4Bkg(x,parbkg);
    }
    
    sgn = FitFunction4Sgn(x,&par[3]);
  }

  if (ftypeOfFit4Bkg==3) {
   
    if(ftypeOfFit4Sgn == 0) {
	bkg=FitFunction4Bkg(x,par);
	sgn=FitFunction4Sgn(x,&par[1]);
    }
    if(ftypeOfFit4Sgn == 1) {
      Double_t parbkg[4]={par[1],par[2],ffactor*par[3],par[0]};
      bkg=FitFunction4Bkg(x,parbkg);
      sgn=FitFunction4Sgn(x,&par[1]);
    }
  }

  total = bkg + sgn;
  
  return  total;
}

//_________________________________________________________________________
Double_t AliHFMassFitter::FitFunction4Sgn (Double_t *x, Double_t *par){
  // Fit function for the signal

  //gaussian = A/(sigma*sqrt(2*pi))*exp(-(x-mean)^2/2/sigma^2)
  //Par:
  // * [0] = integralSgn
  // * [1] = mean
  // * [2] = sigma
  //gaussian = [0]/TMath::Sqrt(2.*TMath::Pi())/[2]*exp[-(x-[1])*(x-[1])/(2*[2]*[2])]

  return par[0]/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]);

}

//__________________________________________________________________________

Double_t AliHFMassFitter::FitFunction4Bkg (Double_t *x, Double_t *par){
  // Fit function for the background

  Double_t maxDeltaM = 4.*fSigmaSgn;
  if(fSideBands && TMath::Abs(x[0]-fMass) < maxDeltaM) {
    TF1::RejectPoint();
    return 0;
  }
  Int_t firstPar=0;
  Double_t gaus2=0,total=-1;
  if(ftypeOfFit4Sgn == 1){
    firstPar=3;
    //gaussian = A/(sigma*sqrt(2*pi))*exp(-(x-mean)^2/2/sigma^2)
    //Par:
    // * [0] = integralSgn
    // * [1] = mean
    // * [2] = sigma
    //gaussian = [0]/TMath::Sqrt(2.*TMath::Pi())/[2]*exp[-(x-[1])*(x-[1])/(2*[2]*[2])]
    gaus2 = FitFunction4Sgn(x,par);
  }

  switch (ftypeOfFit4Bkg){
  case 0:
    //exponential
    //exponential = A*exp(B*x) -> integral(exponential)=A/B*exp(B*x)](min,max)
    //-> A = B*integral/(exp(B*max)-exp(B*min)) where integral can be written
    //as integralTot- integralGaus (=par [2])
    //Par:
    // * [0] = integralBkg;
    // * [1] = B;
    //exponential = [1]*[0]/(exp([1]*max)-exp([1]*min))*exp([1]*x)
    total = par[0+firstPar]*par[1+firstPar]/(TMath::Exp(par[1+firstPar]*fmaxMass)-TMath::Exp(par[1+firstPar]*fminMass))*TMath::Exp(par[1+firstPar]*x[0]);
    break;
  case 1:
    //linear
    //y=a+b*x -> integral = a(max-min)+1/2*b*(max^2-min^2) -> a = (integral-1/2*b*(max^2-min^2))/(max-min)=integral/(max-min)-1/2*b*(max+min)
    // * [0] = integralBkg;
    // * [1] = b;
    total= par[0+firstPar]/(fmaxMass-fminMass)+par[1+firstPar]*(x[0]-0.5*(fmaxMass+fminMass));
    break;
  case 2:
    //polynomial
    //y=a+b*x+c*x**2 -> integral = a(max-min) + 1/2*b*(max^2-min^2) +
    //+ 1/3*c*(max^3-min^3) -> 
    //a = (integral-1/2*b*(max^2-min^2)-1/3*c*(max^3-min^3))/(max-min)
    // * [0] = integralBkg;
    // * [1] = b;
    // * [2] = c;
    total = par[0+firstPar]/(fmaxMass-fminMass)+par[1]*(x[0]-0.5*(fmaxMass+fminMass))+par[2+firstPar]*(x[0]*x[0]-1/3.*(fmaxMass*fmaxMass*fmaxMass-fminMass*fminMass*fminMass)/(fmaxMass-fminMass));
    break;
  case 3:
    total=par[0+firstPar];
    break;
//   default:
//     Types of Fit Functions for Background:
//     * 0 = exponential;
//     * 1 = linear;
//     * 2 = polynomial 2nd order
//     * 3 = no background"<<endl;

  }
  return total+gaus2;
}

//__________________________________________________________________________

void AliHFMassFitter::MassFitter(Bool_t draw){  
  // Main method of the class: performs the fit of the histogram
      
  Int_t nFitPars=0; //total function's number of parameters
  switch (ftypeOfFit4Bkg){
  case 0:
    nFitPars=5; //3+2
    break;
  case 1:
    nFitPars=5; //3+2
    break;
  case 2:
    nFitPars=6; //3+3
    break;
  case 3:
    nFitPars=4; //3+1
    break;
  }

  Int_t bkgPar = nFitPars-3; //background function's number of parameters

  cout<<"nFitPars = "<<nFitPars<<"\tbkgPar = "<<bkgPar<<endl;

  //Total integral
  Double_t totInt = fhistoInvMass->Integral("width");

  fSideBands = kTRUE;
  Double_t width;
  Int_t binleft,binright;
  fNbin=fhistoInvMass->GetNbinsX();
  width=fhistoInvMass->GetBinWidth(8);
  //width=(fmaxMass-fminMass)/(Double_t)fNbin;
  binleft=(Int_t)((fMass-4.*fSigmaSgn-fminMass)/width);
  binright=(Int_t)((fMass+4.*fSigmaSgn-fminMass)/width);

  //sidebands integral - first approx (from histo)
  Double_t sideBandsInt=(Double_t)fhistoInvMass->Integral(1,binleft,"width") + (Double_t)fhistoInvMass->Integral(binright,fNbin,"width");
  cout<<"------nbin = "<<fNbin<<"\twidth = "<<width<<"\tbinleft = "<<binleft<<"\tbinright = "<<binright<<endl;
  cout<<"------sideBandsInt - first approx = "<<sideBandsInt<<endl;

  /*Fit Bkg*/

  TF1 *funcbkg = new TF1("funcbkg",this,&AliHFMassFitter::FitFunction4Bkg,fminMass,fmaxMass,bkgPar,"AliHFMassFitter","FitFunction4Bkg");

  funcbkg->SetLineColor(2); //red

  //first fit for bkg: approx bkgint
 
  switch (ftypeOfFit4Bkg) {
  case 0: //gaus+expo
    funcbkg->SetParNames("BkgInt","Slope"); 
    funcbkg->SetParameters(sideBandsInt,-2.); 
    break;
  case 1:
    funcbkg->SetParNames("BkgInt","Slope");
    funcbkg->SetParameters(sideBandsInt,-100.); 
    break;
  case 2:
    funcbkg->SetParNames("BkgInt","Coef1","Coef2");
    funcbkg->SetParameters(sideBandsInt,-10.,5);
    break;
  case 3:
    if(ftypeOfFit4Sgn==0){
      //in principle it doesn't have effects
      funcbkg->SetParNames("Const");
      funcbkg->SetParameter(0,0.);
      funcbkg->FixParameter(0,0.);
    }
    break;
  default:
    cout<<"Wrong choise of ftypeOfFit4Bkg ("<<ftypeOfFit4Bkg<<")"<<endl;
    return;
    break;
  }
  cout<<"\nBACKGROUND FIT - only combinatorial"<<endl;
  Int_t ftypeOfFit4SgnBkp=ftypeOfFit4Sgn;
  
  Double_t intbkg1=0,slope1=0,conc1=0;
  //if only signal and reflection skip
  if (!(ftypeOfFit4Bkg==3 && ftypeOfFit4Sgn==1)) {
    ftypeOfFit4Sgn=0;
    fhistoInvMass->Fit("funcbkg","R,L,E,0");
   
    for(Int_t i=0;i<bkgPar;i++){
      fFitPars[i]=funcbkg->GetParameter(i);
      //cout<<i<<"\t"<<funcbkg->GetParameter(i)<<"\t";
      fFitPars[nFitPars+2*bkgPar+3+i]= funcbkg->GetParError(i);
      //cout<<nFitPars+2*bkgPar+3+i<<"\t"<< funcbkg->GetParError(i)<<endl;
    }
    
    intbkg1 = funcbkg->GetParameter(0);
    if(ftypeOfFit4Bkg!=3) slope1 = funcbkg->GetParameter(1);
    if(ftypeOfFit4Bkg==2) conc1 = funcbkg->GetParameter(2);
    cout<<"Primo fit: \nintbkg1 = "<<intbkg1<<"\nslope1= "<<slope1<<"\nconc1 = "<<conc1<<endl;
  } else cout<<"\t\t//"<<endl;
  
  ftypeOfFit4Sgn=ftypeOfFit4SgnBkp;
  TF1 *funcbkg1=0;
  if (ftypeOfFit4Sgn == 1) {
    cout<<"\nBACKGROUND FIT WITH REFLECTION"<<endl;
    bkgPar+=3;
    
    cout<<"nFitPars = "<<nFitPars<<"\tbkgPar = "<<bkgPar<<endl;   
    funcbkg1 = new TF1("funcbkg1",this,&AliHFMassFitter::FitFunction4Bkg,fminMass,fmaxMass,bkgPar,"AliHFMassFitter","FitFunction4Bkg");
    funcbkg1->SetLineColor(2); //red

    if(ftypeOfFit4Bkg==2){
      cout<<"*** Polynomial Fit ***"<<endl;
      funcbkg1->SetParNames("IntGB","MeanGB","SigmaGB","BkgInt","Coef1","Coef2");
      funcbkg1->SetParameters(0.5*(totInt-intbkg1),fMass,ffactor*fSigmaSgn,intbkg1,slope1,conc1);
//     cout<<"Parameters set to: "<<0.5*(totInt-intbkg1)<<"\t"<<fMass<<"\t"<<ffactor*fSigmaSgn<<"\t"<<intbkg1<<"\t"<<slope1<<"\t"<<conc1<<"\t"<<endl;
//     cout<<"Limits: ("<<fminMass<<","<<fmaxMass<<")\tnPar = "<<bkgPar<<"\tgsidebands = "<<fSideBands<<endl;
    } else{
      if(ftypeOfFit4Bkg==3) //no background: gaus sign+ gaus broadened
	{
	  cout<<"*** No background Fit ***"<<endl;
	  funcbkg1->SetParNames("IntGB","MeanGB","SigmaGB","Const");
	  funcbkg1->SetParameters(0.5*totInt,fMass,ffactor*fSigmaSgn,0.); 
	  funcbkg1->FixParameter(3,0.);
	} else{ //expo or linear
	  if(ftypeOfFit4Bkg==0) cout<<"*** Exponential Fit ***"<<endl;
	  if(ftypeOfFit4Bkg==1) cout<<"*** Linear Fit ***"<<endl;
	  funcbkg1->SetParNames("IntGB","MeanGB","SigmaGB","BkgInt","Slope");
	  funcbkg1->SetParameters(0.5*(totInt-intbkg1),fMass,ffactor*fSigmaSgn,intbkg1,slope1);
	}
    }
    fhistoInvMass->Fit("funcbkg1","R,L,E,+,0");
  
    for(Int_t i=0;i<bkgPar;i++){
      fFitPars[bkgPar-3+i]=funcbkg1->GetParameter(i);
      //cout<<bkgPar-3+i<<"\t"<<funcbkg1->GetParameter(i);
      fFitPars[nFitPars+3*bkgPar-6+i]= funcbkg1->GetParError(i);
      //cout<<"\t"<<nFitPars+3*bkgPar-6+i<<"\t"<<funcbkg1->GetParError(i)<<endl;
      
    }

    intbkg1=funcbkg1->GetParameter(3);
    if(ftypeOfFit4Bkg!=3) slope1 = funcbkg1->GetParameter(4);
    if(ftypeOfFit4Bkg==2) conc1 = funcbkg1->GetParameter(5);

  } else {
    bkgPar+=3;

    for(Int_t i=0;i<3;i++){
      fFitPars[bkgPar-3+i]=0.;
      //cout<<bkgPar-3+i<<"\t"<<0.<<"\t";
      fFitPars[nFitPars+3*bkgPar-6+i]= 0.;
      //cout<<nFitPars+3*bkgPar-6+i<<"\t"<<0.<<endl;
    }
  
    for(Int_t i=0;i<bkgPar-3;i++){
      fFitPars[bkgPar+i]=funcbkg->GetParameter(i);
      //cout<<bkgPar+i<<"\t"<<funcbkg->GetParameter(i)<<"\t";
      fFitPars[nFitPars+3*bkgPar-3+i]= funcbkg->GetParError(i);
      //cout<<nFitPars+3*bkgPar-3+i<<"\t"<< funcbkg->GetParError(i)<<endl;
    }

   
  }

  //sidebands integral - second approx (from fit)
  fSideBands = kFALSE;
  Double_t bkgInt;
  
  if(ftypeOfFit4Sgn == 1) bkgInt=funcbkg1->Integral(fminMass,fmaxMass);
  else bkgInt=funcbkg->Integral(fminMass,fmaxMass);
  cout<<"------BkgInt(Fit) = "<<bkgInt<<endl;

  //Signal integral - first approx
  Double_t sgnInt;
  sgnInt = totInt-bkgInt;
  cout<<"------TotInt = "<<totInt<<"\tsgnInt = "<<sgnInt<<endl;

  /*Fit All Mass distribution with exponential + gaussian (+gaussiam braodened) */
  TF1 *funcmass = new TF1("funcmass",this,&AliHFMassFitter::FitFunction4MassDistr,fminMass,fmaxMass,nFitPars,"AliHFMassFitter","FitFunction4MassDistr");

  funcmass->SetLineColor(4); //blue

  //Set parameters
  cout<<"\nTOTAL FIT"<<endl;

  if(nFitPars==5){
    funcmass->SetParNames("TotInt","Slope","SgnInt","Mean","Sigma");
    funcmass->SetParameters(totInt,slope1,sgnInt,fMass,fSigmaSgn);
    // cout<<"Parameters set to: "<<totInt<<"\t"<<slope1<<"\t"<<sgnInt<<"\t"<<fMass<<"\t"<<fSigmaSgn<<"\t"<<endl;
//     cout<<"Limits: ("<<fminMass<<","<<fmaxMass<<")\tnPar = "<<nFitPars<<"\tgsidebands = "<<fSideBands<<endl;
    funcmass->FixParameter(0,totInt);
  }
  if (nFitPars==6){
    funcmass->SetParNames("TotInt","Coef1","Coef2","SgnInt","Mean","Sigma");
    funcmass->SetParameters(totInt,slope1,conc1,sgnInt,fMass,fSigmaSgn);
//     cout<<"Parameters set to: "<<totInt<<"\t"<<slope1<<"\t"<<conc1<<"\t"<<sgnInt<<"\t"<<fMass<<"\t"<<fSigmaSgn<<"\t"<<endl;
//     cout<<"Limits: ("<<fminMass<<","<<fmaxMass<<")\tnPar = "<<nFitPars<<"\tgsidebands = "<<fSideBands<<endl;
    funcmass->FixParameter(0,totInt);
  }
  if(nFitPars==4){
    funcmass->SetParNames("Const","SgnInt","Mean","Sigma");
    funcmass->SetParameters(0.,0.5*totInt,fMass,fSigmaSgn);
    funcmass->FixParameter(0,0.);
    //cout<<"Parameters set to: "<<0.5*totInt<<"\t"<<fMass<<"\t"<<fSigmaSgn<<"\t"<<endl;
    //cout<<"Limits: ("<<fminMass<<","<<fmaxMass<<")\tnPar = "<<nFitPars<<"\tgsidebands = "<<fSideBands<<endl;

  }

  fhistoInvMass->Fit("funcmass","R,L,E,+,0");
  cout<<"fit done"<<endl;
  
  for(Int_t i=0;i<nFitPars;i++){
    fFitPars[i+2*bkgPar-3]=funcmass->GetParameter(i);
    fFitPars[nFitPars+4*bkgPar-6+i]= funcmass->GetParError(i);
    cout<<i+2*bkgPar-3<<"\t"<<funcmass->GetParameter(i)<<"\t\t"<<nFitPars+4*bkgPar-6+i<<"\t"<<funcmass->GetParError(i)<<endl;
  }
  /*
  //check: cout parameters  
  for(Int_t i=0;i<2*(nFitPars+2*bkgPar-3);i++){
    cout<<i<<"\t"<<fFitPars[i]<<endl;
    }
  */
  /*
  if(draw){
    TCanvas *canvas=new TCanvas("canvas",fhistoInvMass->GetName());
    TH1F *fhistocopy=new TH1F(*fhistoInvMass);
    canvas->cd();
    fhistocopy->Draw();
    if(ftypeOfFit4Sgn == 1) {
      cout<<"funcbkg1 "<<funcbkg1<<endl;
      funcbkg1->Draw("sames");
    }
    else {
      cout<<"funcbkg "<<funcbkg<<endl;
      funcbkg->Draw("sames");
    }
    cout<<"funcmass "<<funcmass<<endl;
    funcmass->Draw("sames");
    
  }
  */
}

//_________________________________________________________________________
void  AliHFMassFitter::GetFitPars(Float_t *vector25) const {
  // Return fit parameters

  for(Int_t i=0;i<25;i++){
    vector25[i]=fFitPars[i];
  }
}

//_________________________________________________________________________

void AliHFMassFitter::InitNtuParam(char *ntuname) {
  // Create ntuple to keep fit parameters

  fntuParam=new TNtuple(ntuname,"Contains fit parameters","intbkg1:slope1:conc1:intGB:meanGB:sigmaGB:intbkg2:slope2:conc2:inttot:slope3:conc3:intsgn:meansgn:sigmasgn:intbkg1Err:slope1Err:conc1Err:intGBErr:meanGBErr:sigmaGBErr:intbkg2Err:slope2Err:conc2Err:inttotErr:slope3Err:conc3Err:intsgnErr:meansgnErr:sigmasgnErr");
  
}

//_________________________________________________________________________

void AliHFMassFitter::FillNtuParam() {
  // Fill ntuple with fit parameters

  Float_t nothing=0.;

  if (ftypeOfFit4Bkg==2) {
      fntuParam->SetBranchAddress("intbkg1",&fFitPars[0]);
      fntuParam->SetBranchAddress("slope1",&fFitPars[1]);
      fntuParam->SetBranchAddress("conc1",&fFitPars[2]);
      fntuParam->SetBranchAddress("intGB",&fFitPars[3]);
      fntuParam->SetBranchAddress("meanGB",&fFitPars[4]);
      fntuParam->SetBranchAddress("sigmaGB",&fFitPars[5]);
      fntuParam->SetBranchAddress("intbkg2",&fFitPars[6]);
      fntuParam->SetBranchAddress("slope2",&fFitPars[7]);
      fntuParam->SetBranchAddress("conc2",&fFitPars[8]);
      fntuParam->SetBranchAddress("inttot",&fFitPars[9]);
      fntuParam->SetBranchAddress("slope3",&fFitPars[10]);
      fntuParam->SetBranchAddress("conc3",&fFitPars[11]);
      fntuParam->SetBranchAddress("intsgn",&fFitPars[12]);
      fntuParam->SetBranchAddress("meansgn",&fFitPars[13]);
      fntuParam->SetBranchAddress("sigmasgn",&fFitPars[14]);

      fntuParam->SetBranchAddress("intbkg1Err",&fFitPars[15]);
      fntuParam->SetBranchAddress("slope1Err",&fFitPars[16]);
      fntuParam->SetBranchAddress("conc1Err",&fFitPars[17]);
      fntuParam->SetBranchAddress("intGBErr",&fFitPars[18]);
      fntuParam->SetBranchAddress("meanGBErr",&fFitPars[19]);
      fntuParam->SetBranchAddress("sigmaGBErr",&fFitPars[20]);
      fntuParam->SetBranchAddress("intbkg2Err",&fFitPars[21]);
      fntuParam->SetBranchAddress("slope2Err",&fFitPars[22]);
      fntuParam->SetBranchAddress("conc2Err",&fFitPars[23]);
      fntuParam->SetBranchAddress("inttotErr",&fFitPars[24]);
      fntuParam->SetBranchAddress("slope3Err",&fFitPars[25]);
      fntuParam->SetBranchAddress("conc3Err",&fFitPars[26]);
      fntuParam->SetBranchAddress("intsgnErr",&fFitPars[27]);
      fntuParam->SetBranchAddress("meansgnErr",&fFitPars[28]);
      fntuParam->SetBranchAddress("sigmasgnErr",&fFitPars[29]);
    
  } else {
    
    if(ftypeOfFit4Bkg==3){
      fntuParam->SetBranchAddress("intbkg1",&fFitPars[0]);
      fntuParam->SetBranchAddress("slope1",&nothing);
      fntuParam->SetBranchAddress("conc1",&nothing);
      fntuParam->SetBranchAddress("intGB",&fFitPars[1]);
      fntuParam->SetBranchAddress("meanGB",&fFitPars[2]);
      fntuParam->SetBranchAddress("sigmaGB",&fFitPars[3]);
      fntuParam->SetBranchAddress("intbkg2",&fFitPars[4]);
      fntuParam->SetBranchAddress("slope2",&nothing);
      fntuParam->SetBranchAddress("conc2",&nothing);
      fntuParam->SetBranchAddress("inttot",&fFitPars[6]);
      fntuParam->SetBranchAddress("slope3",&nothing);
      fntuParam->SetBranchAddress("conc3",&nothing);
      fntuParam->SetBranchAddress("intsgn",&fFitPars[6]);
      fntuParam->SetBranchAddress("meansgn",&fFitPars[7]);
      fntuParam->SetBranchAddress("sigmasgn",&fFitPars[8]);

      fntuParam->SetBranchAddress("intbkg1Err",&fFitPars[9]);
      fntuParam->SetBranchAddress("slope1Err",&nothing);
      fntuParam->SetBranchAddress("conc1Err",&nothing);
      fntuParam->SetBranchAddress("intGBErr",&fFitPars[10]);
      fntuParam->SetBranchAddress("meanGBErr",&fFitPars[11]);
      fntuParam->SetBranchAddress("sigmaGBErr",&fFitPars[12]);
      fntuParam->SetBranchAddress("intbkg2Err",&fFitPars[13]);
      fntuParam->SetBranchAddress("slope2Err",&nothing);
      fntuParam->SetBranchAddress("conc2Err",&nothing);
      fntuParam->SetBranchAddress("inttotErr",&fFitPars[15]);
      fntuParam->SetBranchAddress("slope3Err",&nothing);
      fntuParam->SetBranchAddress("conc3Err",&nothing);
      fntuParam->SetBranchAddress("intsgnErr",&fFitPars[15]);
      fntuParam->SetBranchAddress("meansgnErr",&fFitPars[16]);
      fntuParam->SetBranchAddress("sigmasgnErr",&fFitPars[17]);

    }
    else{
      fntuParam->SetBranchAddress("intbkg1",&fFitPars[0]);
      fntuParam->SetBranchAddress("slope1",&fFitPars[1]);
      fntuParam->SetBranchAddress("conc1",&nothing);
      fntuParam->SetBranchAddress("intGB",&fFitPars[2]);
      fntuParam->SetBranchAddress("meanGB",&fFitPars[3]);
      fntuParam->SetBranchAddress("sigmaGB",&fFitPars[4]);
      fntuParam->SetBranchAddress("intbkg2",&fFitPars[5]);
      fntuParam->SetBranchAddress("slope2",&fFitPars[6]);
      fntuParam->SetBranchAddress("conc2",&nothing);
      fntuParam->SetBranchAddress("inttot",&fFitPars[7]);
      fntuParam->SetBranchAddress("slope3",&fFitPars[8]);
      fntuParam->SetBranchAddress("conc3",&nothing);
      fntuParam->SetBranchAddress("intsgn",&fFitPars[9]);
      fntuParam->SetBranchAddress("meansgn",&fFitPars[10]);
      fntuParam->SetBranchAddress("sigmasgn",&fFitPars[11]);

      fntuParam->SetBranchAddress("intbkg1Err",&fFitPars[12]);
      fntuParam->SetBranchAddress("slope1Err",&fFitPars[13]);
      fntuParam->SetBranchAddress("conc1Err",&nothing);
      fntuParam->SetBranchAddress("intGBErr",&fFitPars[14]);
      fntuParam->SetBranchAddress("meanGBErr",&fFitPars[15]);
      fntuParam->SetBranchAddress("sigmaGBErr",&fFitPars[16]);
      fntuParam->SetBranchAddress("intbkg2Err",&fFitPars[17]);
      fntuParam->SetBranchAddress("slope2Err",&fFitPars[18]);
      fntuParam->SetBranchAddress("conc2Err",&nothing);
      fntuParam->SetBranchAddress("inttotErr",&fFitPars[19]);
      fntuParam->SetBranchAddress("slope3Err",&fFitPars[20]);
      fntuParam->SetBranchAddress("conc3Err",&nothing);
      fntuParam->SetBranchAddress("intsgnErr",&fFitPars[21]);
      fntuParam->SetBranchAddress("meansgnErr",&fFitPars[22]);
      fntuParam->SetBranchAddress("sigmasgnErr",&fFitPars[23]);
    }
     
  }
  fntuParam->TTree::Fill();
}

//_________________________________________________________________________

TNtuple* AliHFMassFitter::NtuParamOneShot(char *ntuname){
  // Create, fill and return ntuple with fit parameters

  InitNtuParam(ntuname);
  FillNtuParam();
  return fntuParam;
}
//_________________________________________________________________________

void AliHFMassFitter::RebinMass(Int_t bingroup){
  // Rebin invariant mass histogram

  if(bingroup<1){
    cout<<"Error! Cannot group "<<bingroup<<" bins\n";
    fNbin=fhistoInvMass->GetNbinsX();
    cout<<"Kept original number of bins: "<<fNbin<<endl;
  } else{
    fhistoInvMass->Rebin(bingroup);
    fNbin = fhistoInvMass->GetNbinsX();
    cout<<"New number of bins: "<<fNbin<<endl;
  } 
  
       
}

//_________________________________________________________________________

void AliHFMassFitter::Signal(Double_t nOfSigma,Double_t &signal,Double_t &errsignal) const {
  // Return signal integral in mean+- n sigma

  TF1 *funcmass=fhistoInvMass->GetFunction("funcmass");
  Double_t mean=0,sigma=1;
  Double_t intS,intSerr;

  if(ftypeOfFit4Bkg == 2) { //polynomial
    mean=funcmass->GetParameter(4); //mean
    sigma=funcmass->GetParameter(5); //sigma
    intS=fFitPars[12];
    intSerr=fFitPars[27];
  } else if(ftypeOfFit4Bkg == 3){ //no background
    mean=funcmass->GetParameter(2); //mean
    sigma=funcmass->GetParameter(3); //sigma
    intS=fFitPars[6];
    intSerr=fFitPars[15];
  } else { //expo or linear
    mean=funcmass->GetParameter(3); //mean
    sigma=funcmass->GetParameter(4); //sigma
    intS=fFitPars[9];
    intSerr=fFitPars[21];
  }

  TF1 *funcbkg=0;
  if(ftypeOfFit4Sgn == 0) funcbkg=fhistoInvMass->GetFunction("funcbkg");
  else funcbkg=fhistoInvMass->GetFunction("funcbkg1");
  Double_t min=mean-nOfSigma*sigma;
  Double_t max=mean+nOfSigma*sigma;
  signal=funcmass->Integral(min,max)-funcbkg->Integral(min,max);
  errsignal=intSerr/intS*signal; // assume relative error is the same as for total integral
  
  return;
}
//_________________________________________________________________________

void AliHFMassFitter::Background(Double_t nOfSigma,Double_t &background,Double_t &errbackground) const {
  // Return background integral in mean+- n sigma

  TF1 *funcmass=fhistoInvMass->GetFunction("funcmass");
  Double_t mean=0,sigma=1;
  Double_t intB,intBerr,err;

  if(ftypeOfFit4Bkg == 2) { //polynomial
    mean=funcmass->GetParameter(4); //mean
    sigma=funcmass->GetParameter(5); //sigma
    intB=fFitPars[9]-fFitPars[12];
    intBerr=fFitPars[27];
  } else if(ftypeOfFit4Bkg == 3){ //no background
    mean=funcmass->GetParameter(2); //mean
    sigma=funcmass->GetParameter(3); //sigma
    if (ftypeOfFit4Sgn == 1){ //reflection
      intB=fFitPars[1]; 
      intBerr=fFitPars[9];
    } else {
      intB=-1;intBerr=-1;
      err=0;
    }
  } else { //expo or linear
    mean=funcmass->GetParameter(3); //mean
    sigma=funcmass->GetParameter(4); //sigma
    intB=fFitPars[7]-fFitPars[9];
    intBerr=fFitPars[21];
  }

  TF1 *funcbkg=0;
  if(ftypeOfFit4Sgn == 0) funcbkg=fhistoInvMass->GetFunction("funcbkg");
  else funcbkg=fhistoInvMass->GetFunction("funcbkg1");
  Double_t min=mean-nOfSigma*sigma;
  Double_t max=mean+nOfSigma*sigma;
  background=funcbkg->Integral(min,max);
  errbackground=intBerr/intB*background;

  return;

}

//__________________________________________________________________________

void AliHFMassFitter::Significance(Double_t nOfSigma,Double_t &significance,Double_t &errsignificance)  const {
  // Return significance in mean+- n sigma

  Double_t signal,errsignal,background,errbackground;
  Signal(nOfSigma,signal,errsignal);
  Background(nOfSigma,background,errbackground);

  significance =  signal/TMath::Sqrt(signal+background);
  errsignificance = TMath::Sqrt(significance*significance/(signal+background)/(signal+background)*(1/4.*errsignal*errsignal+errbackground*errbackground)+significance*significance/signal/signal*errsignal*errsignal);

  return;
}

//__________________________________________________________________________
