#include "AliHBTCorrectOSLCorrelFctn.h"
//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTCorrectOSLCorrelFctn                        //
//                                                   //
// Class for calculating Q Invariant correlation     //
// taking to the account resolution of the           //
// detector and coulomb effects.                     //
//                                                   //
///////////////////////////////////////////////////////

#include <TH3.h>

#include <AliAODParticle.h>
#include <AliHBTPair.h>

AliHBTCorrectOSLCorrelFctn::AliHBTCorrectOSLCorrelFctn(const char* name, const char* title):
 AliHBTOnePairFctn3D(name,title),
 AliHBTCorrectedCorrelFctn(),
 fMeasCorrelFctn(0x0),
 fSmearedNumer(0x0),
 fSmearedDenom(0x0),
 fMeasNumer(0x0),
 fMeasDenom(0x0),
 fLambda(0.0),
 fROutSq(0.0),
 fRSideSq(0.0),
 fRLongSq(0.0)
{
//ctor
  fWriteNumAndDen = kTRUE;//change default behaviour
}
/******************************************************************/

AliHBTCorrectOSLCorrelFctn::AliHBTCorrectOSLCorrelFctn(const Char_t *name, const Char_t *title,
	         Int_t nXbins, Double_t maxXval, Double_t minXval, 
	         Int_t nYbins, Double_t maxYval, Double_t minYval, 
	         Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTOnePairFctn3D(name,title,nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval),
 AliHBTCorrectedCorrelFctn(),
 fMeasCorrelFctn(0x0),
 fSmearedNumer(0x0),
 fSmearedDenom(0x0),
 fMeasNumer(0x0),
 fMeasDenom(0x0),
 fLambda(0.0),
 fROutSq(0.0),
 fRSideSq(0.0),
 fRLongSq(0.0)
{
//ctor
  fWriteNumAndDen = kTRUE;//change default behaviour}
}

/******************************************************************/
AliHBTCorrectOSLCorrelFctn::AliHBTCorrectOSLCorrelFctn(const AliHBTCorrectOSLCorrelFctn& in):
 AliHBTOnePairFctn3D(in),
 AliHBTCorrectedCorrelFctn(),
 fMeasCorrelFctn(0x0),
 fSmearedNumer(0x0),
 fSmearedDenom(0x0),
 fMeasNumer(0x0),
 fMeasDenom(0x0),
 fLambda(0.0),
 fROutSq(0.0),
 fRSideSq(0.0),
 fRLongSq(0.0)
{
//cpy constructor
 in.Copy(*this);
}
/******************************************************************/

AliHBTCorrectOSLCorrelFctn::~AliHBTCorrectOSLCorrelFctn()
{
 //dtor
 delete fMeasCorrelFctn;
 delete fSmearedNumer;
 delete fSmearedDenom;
 delete fMeasNumer;
 delete fMeasDenom;
}
/******************************************************************/

void AliHBTCorrectOSLCorrelFctn::BuildHistos(Int_t nxbins, Float_t xmax, Float_t xmin,
                    	     Int_t nybins, Float_t ymax, Float_t ymin,
		     Int_t nzbins, Float_t zmax, Float_t zmin)
{
//build histograms
  if (AliVAODParticle::GetDebug()>1) Info("BuildHistos","Enetered BuildHistos(...)");
  
  AliHBTFunction3D::BuildHistos(nxbins,xmax,xmin,nybins,ymax,ymin,nzbins,zmax,zmin);
  
  TString numstr = fName + " Smeared Numerator";  //title and name of the numerator histogram
  TString denstr = fName + " Smeared Denominator";//title and name of the denominator histogram
  
  fSmearedNumer = new TH3F(numstr.Data(),numstr.Data(),nxbins,xmin,xmax,nybins,ymin,ymax,nzbins,zmin,zmax);
  fSmearedDenom = new TH3F(denstr.Data(),denstr.Data(),nxbins,xmin,xmax,nybins,ymin,ymax,nzbins,zmin,zmax);
  fSmearedNumer->Sumw2();
  fSmearedDenom->Sumw2();
  
  if (fMeasCorrelFctn == 0x0)
   { 
     numstr = fName + " Measured Numerator";  //title and name of the numerator histogram
     denstr = fName + " Measured Denominator";//title and name of the denominator histogram
    
     fMeasNumer = new TH3F(numstr.Data(),numstr.Data(),nxbins,xmin,xmax,nybins,ymin,ymax,nzbins,zmin,zmax);
     fMeasDenom = new TH3F(denstr.Data(),denstr.Data(),nxbins,xmin,xmax,nybins,ymin,ymax,nzbins,zmin,zmax);
     fMeasNumer->Sumw2();
     fMeasDenom->Sumw2();
   }

}
/******************************************************************/

void AliHBTCorrectOSLCorrelFctn::Init()
{
//Init
  AliHBTOnePairFctn3D::Init();
  Info("Init","");
  if ( (fSmearedNumer == 0x0) || (fSmearedDenom == 0x0) )
   {
     if (fNumerator == 0x0) Fatal("Init","Sth. goes wrong");
     Int_t nxbins = fNumerator->GetNbinsX();
     Float_t xmax = fNumerator->GetXaxis()->GetXmax();
     Float_t xmin = fNumerator->GetXaxis()->GetXmin();
     Int_t nybins = fNumerator->GetNbinsY();
     Float_t ymax = fNumerator->GetYaxis()->GetXmax();
     Float_t ymin = fNumerator->GetYaxis()->GetXmin();
     Int_t nzbins = fNumerator->GetNbinsZ();
     Float_t zmax = fNumerator->GetZaxis()->GetXmax();
     Float_t zmin = fNumerator->GetZaxis()->GetXmin();
     BuildHistos(nxbins,xmax,xmin, nybins,ymax,ymin, nzbins,zmax,zmin);
   }
   
  fSmearedNumer->Reset();
  fSmearedDenom->Reset();
  if (fMeasNumer) fMeasNumer->Reset();
  if (fMeasDenom) fMeasDenom->Reset();
}
/******************************************************************/

void AliHBTCorrectOSLCorrelFctn::SetInitialValues(Double_t lambda, Double_t rout, Double_t rside, Double_t rlong)
{
  //Sets assumed parameters
 fLambda = lambda;
 fROutSq = rout*rout;
 fRSideSq = rside*rside;
 fRLongSq = rlong*rlong;
}

/******************************************************************/

void AliHBTCorrectOSLCorrelFctn::ProcessSameEventParticles(AliHBTPair* pair)
{
 //Processes particles that originates from the same event
  
  return; //we already heave the measured in hand
 
 
  if (fMeasNumer == 0x0) return;
  pair = CheckPair(pair);
  if( pair == 0x0) return;
  fMeasNumer->Fill(pair->GetQInv());
}
/******************************************************************/

void AliHBTCorrectOSLCorrelFctn::ProcessDiffEventParticles(AliHBTPair* pair)
{
//Process different events 
  static AliAODParticle part1, part2;
  static AliHBTPair smearedpair(&part1,&part2);
  
  pair = CheckPair(pair);
  if( pair == 0x0) return;
  
  Double_t cc = GetCoulombCorrection(pair);

  Double_t qout,qside,qlong;
  GetValues(pair,qout,qside,qlong);

  //measured histogram -> if we are interested 
  //only if fMeasCorrelFctn is not specified by user

//  if (fMeasDenom) fMeasDenom->Fill(qout,qside,qlong,cc);

  Smear(pair,smearedpair);
  Double_t modelval = GetModelValue(qout,qside,qlong);  
  //Ideal histogram
  fNumerator->Fill(qout,qside,qlong,modelval*cc);
  fDenominator->Fill(qout,qside,qlong,cc);
  
  //Smeared histogram

  Double_t smearedqout,smearedqside,smearedqlong;
  
  GetValues(&smearedpair,smearedqout,smearedqside,smearedqlong);
  
  fSmearedNumer->Fill(smearedqout,smearedqside,smearedqlong,modelval);
  
  Double_t smearedcc = GetCoulombCorrection(&smearedpair);
  fSmearedDenom->Fill(smearedqout,smearedqside,smearedqlong,smearedcc);
  
}
/******************************************************************/

void AliHBTCorrectOSLCorrelFctn::WriteFunction()
{
  AliHBTFunction::WriteFunction();
  if (fSmearedNumer) fSmearedNumer->Write();
  if (fSmearedDenom) fSmearedDenom->Write();

}
/******************************************************************/

TH1* AliHBTCorrectOSLCorrelFctn::GetResult()
{
  //reuturns result histogram
   delete fRatio;
   fRatio = GetRatio(Scale());
   return fRatio;
}

void AliHBTCorrectOSLCorrelFctn::GetValues(AliHBTPair* pair, Double_t& x, Double_t& y, Double_t& z) const 
{
    //calculates values of that function
  //qout qside and qlong
  
  x=pair->GetQOutLCMS(); 
  y=pair->GetQSideLCMS(); 
  z=pair->GetQLongLCMS();
  if (fAbs)
   {
     x = TMath::Abs(x);
     y = TMath::Abs(y);
     z = TMath::Abs(z);
   }

}
