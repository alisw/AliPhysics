#include "AliHBTCorrectQInvCorrelFctn.h"
//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTCorrectQInvCorrelFctn                       //
//                                                   //
// Class for calculating Q Invariant correlation     //
// taking to the account resolution of the           //
// detector and coulomb effects.                     //
// Implemented on basis of STAR procedure            //
// implemented and described by Michael A. Lisa      //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

//Parameters fit from pi+ pi+ resolution analysis for pair with qinv<50MeV
// chi2/NFD = 97/157
// FCN=98.0971 FROM MIGRAD    STATUS=CONVERGED     332 CALLS         333 TOTAL
//   EDM=2.13364e-12    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.4 per cent
//  EXT PARAMETER                                   STEP         FIRST
//  NO.   NAME         VALUE          ERROR         SIZE      DERIVATIVE
//   1   fThetaA     2.72546e-03   1.43905e-04  -0.00000e+00   5.54375e-02
//   2   fThetaB     1.87116e-04   5.11862e-05   0.00000e+00   3.66500e-01
//   3  fThetaAlpha -2.36868e+00   1.83230e-01   0.00000e+00  -7.01301e-05

// FCN=120.603 FROM MIGRAD    STATUS=CONVERGED     117 CALLS         118 TOTAL
//  EDM=2.5153e-12    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.4 per cent
//  EXT PARAMETER                                   STEP         FIRST
//  NO.   NAME         VALUE          ERROR         SIZE      DERIVATIVE
//   1   fPhiA       1.93913e-03   1.31059e-04  -0.00000e+00  -5.87280e-02
//   2   fPhiA       2.48687e-04   5.41251e-05   0.00000e+00  -1.87361e-01
//   3  fPhiAlpha   -2.22649e+00   1.44503e-01   0.00000e+00   8.97538e-06

#include <TH1.h>
#include <TH3.h>
#include <TF1.h>
#include <TRandom.h>
#include <AliAODParticle.h>


ClassImp(AliHBTCorrectQInvCorrelFctn)

AliHBTCorrectQInvCorrelFctn::AliHBTCorrectQInvCorrelFctn(const char* name,const char* title):
  AliHBTOnePairFctn1D(name,title),
  fMeasCorrelFctn(0x0),
  fMeasNumer(0x0),
  fMeasDenom(0x0),
  fSmearedNumer(0x0),
  fSmearedDenom(0x0),
  fDPtOverPtRMS(0.004),
  fThetaA(2.72e-03),
  fThetaB(1.87e-04),
  fThetaAlpha(-2.4),
  fPhiA(1.94e-03),
  fPhiB(2.5e-04),
  fPhiAlpha(-2.2),
  fR2(0.0),
  fLambda(0.0),
  fRConvergenceTreshold(0.3),
  fLambdaConvergenceTreshold(0.05)
{

}
/******************************************************************/

AliHBTCorrectQInvCorrelFctn::AliHBTCorrectQInvCorrelFctn(TH1D* measqinv,const char* name,const char* title):
  AliHBTOnePairFctn1D(name,title,
                      measqinv->GetNbinsX(),
	  measqinv->GetXaxis()->GetXmax(),
	  measqinv->GetXaxis()->GetXmin()),
  fMeasCorrelFctn(measqinv),
  fMeasNumer(0x0),
  fMeasDenom(0x0),
  fSmearedNumer(0x0),
  fSmearedDenom(0x0),
  fDPtOverPtRMS(0.004),
  fThetaA(2.72e-03),
  fThetaB(1.87e-04),
  fThetaAlpha(-2.4),
  fPhiA(1.94e-03),
  fPhiB(2.5e-04),
  fPhiAlpha(-2.2),
  fR2(0.0),
  fLambda(0.0),
  fRConvergenceTreshold(0.3),
  fLambdaConvergenceTreshold(0.05)
{
  
}
/******************************************************************/

AliHBTCorrectQInvCorrelFctn::AliHBTCorrectQInvCorrelFctn(const char* name, const char* title,
	            Int_t nbins, Float_t maxXval, Float_t minXval):
  AliHBTOnePairFctn1D(name,title,nbins,maxXval,minXval),
  fMeasCorrelFctn(0x0),
  fMeasNumer(0x0),
  fMeasDenom(0x0),
  fSmearedNumer(0x0),
  fSmearedDenom(0x0),
  fDPtOverPtRMS(0.004),
  fThetaA(2.72e-03),
  fThetaB(1.87e-04),
  fThetaAlpha(-2.4),
  fPhiA(1.94e-03),
  fPhiB(2.5e-04),
  fPhiAlpha(-2.2),
  fR2(0.0),
  fLambda(0.0),
  fRConvergenceTreshold(0.3),
  fLambdaConvergenceTreshold(0.05)
{
//ctor
}	            
/******************************************************************/
AliHBTCorrectQInvCorrelFctn::AliHBTCorrectQInvCorrelFctn(const AliHBTCorrectQInvCorrelFctn& in):
  AliHBTOnePairFctn1D(in),
  fMeasCorrelFctn(0x0),
  fMeasNumer(0x0),
  fMeasDenom(0x0),
  fSmearedNumer(0x0),
  fSmearedDenom(0x0),
  fDPtOverPtRMS(0),
  fThetaA(0),
  fThetaB(0),
  fThetaAlpha(0),
  fPhiA(0),
  fPhiB(0),
  fPhiAlpha(0),
  fR2(0.0),
  fLambda(0.0),
  fRConvergenceTreshold(0),
  fLambdaConvergenceTreshold(0)
{
//cpy ;ctor
 in.Copy(*this);
}

AliHBTCorrectQInvCorrelFctn::~AliHBTCorrectQInvCorrelFctn()
{
//dtor
  delete fMeasCorrelFctn;
  delete fSmearedNumer;
  delete fSmearedDenom;
  delete fMeasNumer;
  delete fMeasDenom;
}
/******************************************************************/

void AliHBTCorrectQInvCorrelFctn::BuildHistos(Int_t nbins, Float_t max, Float_t min)
{
  AliHBTFunction1D::BuildHistos(nbins,max,min);
  
  TString numstr = fName + " Smeared Numerator";  //title and name of the numerator histogram
  TString denstr = fName + " Smeared Denominator";//title and name of the denominator histogram
  
  fSmearedNumer = new TH1D(numstr.Data(),numstr.Data(),nbins,min,max);
  fSmearedDenom = new TH1D(denstr.Data(),denstr.Data(),nbins,min,max);
  fSmearedNumer->Sumw2();
  fSmearedDenom->Sumw2();
  
  if (fMeasCorrelFctn == 0x0)
   { 
     numstr = fName + " Measured Numerator";  //title and name of the numerator histogram
     denstr = fName + " Measured Denominator";//title and name of the denominator histogram
    
     fMeasNumer = new TH1D(numstr.Data(),numstr.Data(),nbins,min,max);
     fMeasDenom = new TH1D(denstr.Data(),denstr.Data(),nbins,min,max);
     fMeasNumer->Sumw2();
     fMeasDenom->Sumw2();
   }
}
/******************************************************************/

void AliHBTCorrectQInvCorrelFctn::Init()
{
//Init
  AliHBTOnePairFctn1D::Init();
  Info("Init","");
  fSmearedNumer->Reset();
  fSmearedDenom->Reset();
  if (fMeasNumer) fMeasNumer->Reset();
  if (fMeasDenom) fMeasDenom->Reset();
  fFittedR = -1.0;
  fFittedLambda = 0.0;
}
/******************************************************************/

void AliHBTCorrectQInvCorrelFctn::ProcessSameEventParticles(AliHBTPair* pair)
{
 //Processes particles that originates from the same event
  if (fMeasNumer == 0x0) return;
  pair = CheckPair(pair);
  if( pair == 0x0) return;
  fMeasNumer->Fill(pair->GetQInv());
}
/******************************************************************/

void AliHBTCorrectQInvCorrelFctn::ProcessDiffEventParticles(AliHBTPair* pair)
{
//Process different events 
  static AliAODParticle part1, part2;
  static AliHBTPair smearedpair(&part1,&part2);
  
  pair = CheckPair(pair);
  if( pair == 0x0) return;
  Double_t cc = GetCoulombCorrection(pair);
  Double_t qinv = pair->GetQInv();

  //measured histogram -> if we are interested 
  //only if fMeasCorrelFctn is not specified by user
  if (fMeasDenom) fMeasDenom->Fill(qinv,cc);

  Smear(pair,smearedpair);
  Double_t modelqinv = GetModelValue(qinv);  
  //Ideal histogram
  fNumerator->Fill(qinv,modelqinv*cc);
  fDenominator->Fill(qinv,cc);
  
  //Smeared histogram
  Double_t smearedqinv = smearedpair.GetQInv();
  fSmearedNumer->Fill(smearedqinv,modelqinv);
  Double_t smearedcc = GetCoulombCorrection(&smearedpair);
  fSmearedDenom->Fill(smearedqinv,smearedcc);
  
}
/******************************************************************/
void AliHBTCorrectQInvCorrelFctn::Smear(AliHBTPair* pair,AliHBTPair& smeared)
{
//Smears pair
  Smear(pair->Particle1(),smeared.Particle1());
  Smear(pair->Particle2(),smeared.Particle2());
  smeared.Changed();
}
/******************************************************************/

void AliHBTCorrectQInvCorrelFctn::Smear(AliVAODParticle* part, AliVAODParticle* smeared)
{
 //Smears momenta
  Double_t sin2theta = TMath::Sin(part->Theta());
  sin2theta = sin2theta*sin2theta;
  Double_t pt = part->Pt();

  double dPtDivPt = gRandom->Gaus(0.0,fDPtOverPtRMS);
  double dphi = gRandom->Gaus(0.0,fPhiA+fPhiB*TMath::Power(pt,fPhiAlpha));
  double dtheta = gRandom->Gaus(0.0,fPhiA+fPhiB*TMath::Power(pt,fThetaAlpha));
  
  Double_t smearedPx = part->Px()*(1.0+dPtDivPt) - part->Py()*dphi;
//  fourmom.setX(px*(1.0+dPtDivPt) - py*dphi);
  Double_t smearedPy = part->Py()*(1.0+dPtDivPt) - part->Px()*dphi;
//  fourmom.setY(py*(1.0+dPtDivPt) + px*dphi);
  Double_t smearedPz = part->Pz()*(1.0+dPtDivPt) - pt*dtheta/sin2theta;
//  fourmom.setZ(pz*(1.0+dPtDivPt) - pT*dtheta/sin2theta);
  
  Double_t mass2 = part->Mass()*part->Mass();
  Double_t e = mass2 + smearedPx*smearedPx + 
                       smearedPy*smearedPy + 
                       smearedPz*smearedPz;
	   
  smeared->SetMomentum(smearedPx,smearedPy,smearedPz,TMath::Sqrt(e));
}
/******************************************************************/

void AliHBTCorrectQInvCorrelFctn::SetInitialValues(Double_t lambda, Double_t r)
{
 //Sets Initial Values
  fLambda = lambda;
  fR2 = r*r;
}
/******************************************************************/
void AliHBTCorrectQInvCorrelFctn::MakeMeasCF()
{
  //makes measured correlation function
  delete fMeasCorrelFctn;
  fMeasCorrelFctn = 0x0;
  
  if (fMeasNumer&&fMeasDenom)
   {
     Double_t measscale = Scale(fMeasNumer,fMeasDenom);
     if (measscale == 0.0)
      {
        Error("GetResult","GetRatio for measured CF returned 0.0");
        return;
      }
     TString str = fName + "measured ratio";
     fMeasCorrelFctn = (TH1D*)fMeasNumer->Clone(str.Data());
     fMeasCorrelFctn->SetTitle(str.Data());
     fMeasCorrelFctn->Divide(fMeasNumer,fMeasDenom,measscale);
   }
}

TH1* AliHBTCorrectQInvCorrelFctn::GetResult()
{
  //In case we don't have yet Measured Correlation Function
  //Try to get it 
  //result is
  //        N[meas]   N[ideal]/D[ideal]
  //C(Q) =  ------- * -----------------
  //        D[meas]   N[smear]/D[smear]
 
  TString str;
  if (fMeasCorrelFctn == 0x0) MakeMeasCF();
  
  if (fMeasCorrelFctn == 0x0)
   {
     Error("GetResult",  
           "Measured Correlation Function is not defined and measured numeraor and/or denominator are/is null");
     return 0x0;
   }
  
  TH1D* ideal = (TH1D*)GetRatio(Scale());
  if (ideal == 0x0)
   {
     Error("GetResult","Ratio of ideal histograms is null");
     return 0x0;
   }
  str = fName + " smeared ratio";
  TH1D* smearedCF = (TH1D*)fSmearedNumer->Clone(str.Data());
  smearedCF->SetTitle(str.Data());
  Double_t smearedscale = Scale(fSmearedNumer,fSmearedDenom);
  smearedCF->Divide(fSmearedNumer,fSmearedDenom,smearedscale);
  
  str = fName + " product meas ideal CF";
  TH1D* measideal = (TH1D*)ideal->Clone(str.Data());
  measideal->Multiply(ideal,fMeasCorrelFctn);
  
  str = fName + " Corrected Result";
  TH1D* result = (TH1D*)fSmearedNumer->Clone(str.Data());
  result->SetTitle(str.Data());
  
  Double_t resultscale = Scale(measideal,smearedCF);
  result->Divide(measideal,smearedCF,resultscale);
  return result;
  
}
/******************************************************************/

void AliHBTCorrectQInvCorrelFctn::Fit()
{
//fits resuting histogram with function 1.0 + [0]*exp([1]*[1]*x*x/(-0.038936366329))
//where [0] is lambda
//      [1] is radius
//   0.038936366329 - constant needed for units transformation eV(c=1,etc.) -> SI

  Info("Fit","Before fFittedLambda = %f",fFittedLambda);
  Info("Fit","Before fFittedR = %f",fFittedR);
  TH1D* result = (TH1D*)GetResult();
  if (result == 0x0)
   {
     Error("Fit","Can not get result");
     return;
   }
  TF1* fitfctn = new TF1("fitfctn","1.0 + [0]*exp([1]*[1]*x*x/(-0.038936366329))");
  fitfctn->SetParameter(0,1.0);
  fitfctn->SetParameter(1,6.0);
  Float_t max = result->GetXaxis()->GetXmax();
  Info("Fit","Max is %f",max);
  result->Fit(fitfctn,"0","",0.008,max);
  fFittedLambda = fitfctn->GetParameter(0);
  fFittedR = fitfctn->GetParameter(1);
  Info("Fit","After fFittedLambda = %f",fFittedLambda);
  Info("Fit","After fFittedR = %f",fFittedR);
  delete fitfctn;
  delete result;
}
/******************************************************************/

Bool_t AliHBTCorrectQInvCorrelFctn::IsConverged()
{
  //check if fitting was performed
  if (fFittedR <= 0.0)
   {
     Fit();//if not do fit first
     if (fFittedR <= 0.0)
      {
        Error("IsConverged","Fitting failed");
        return kFALSE;
      }
   }

  Double_t guessedR = TMath::Sqrt(fR2);
  
  Info("IsConverged","Fitted   lambda            : %8f Fitted  Radius             : %8f",fFittedLambda,fFittedR);
  Info("IsConverged","Guessed  lambda            : %8f Guessed Radius             : %8f",fLambda,guessedR);
  Info("IsConverged","Demanded lambda convergence: %8f Demanded Radius convergence: %8f",
       fLambdaConvergenceTreshold,fRConvergenceTreshold);
  
  if ( (TMath::Abs(fLambda-fFittedLambda)<fLambdaConvergenceTreshold) && 
       (TMath::Abs(fFittedR-guessedR)<fRConvergenceTreshold) )
    {
      Info("IsConverged","Cnvergence reached");
      return kTRUE;
    }
  else
   {
      Info("IsConverged","Cnvergence NOT reached");
      return kFALSE;
   }
}
/******************************************************************/

Double_t AliHBTCorrectQInvCorrelFctn::GetFittedRadius()
{
 //Returns Fitted radius
  if (fFittedR <= 0.0) Fit();
  return fFittedR;
}
/******************************************************************/

Double_t AliHBTCorrectQInvCorrelFctn::GetFittedLambda()
{
 //Returns Fitted Intercept paramter
  if (fFittedR <= 0.0) Fit();
  return fFittedLambda;
}
/******************************************************************/

void AliHBTCorrectQInvCorrelFctn::WriteAll()
{
//Writes function and all additional information
  Write();
  if (fMeasCorrelFctn) fMeasCorrelFctn->Write();
  if (fMeasNumer ) fMeasNumer->Write();
  if (fMeasDenom ) fMeasDenom->Write();
  if (fSmearedNumer) fSmearedNumer->Write();
  if (fSmearedDenom) fSmearedDenom->Write();
  if (fSmearedNumer && fSmearedDenom)
   {
    TString str = fName + " smeared ratio";
    TH1D* smearedCF = (TH1D*)fSmearedNumer->Clone(str.Data());
    smearedCF->SetTitle(str.Data());
    Double_t smearedscale = Scale(fSmearedNumer,fSmearedDenom);
    smearedCF->Divide(fSmearedNumer,fSmearedDenom,smearedscale);
   }
}

