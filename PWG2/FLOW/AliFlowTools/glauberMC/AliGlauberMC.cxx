/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////////////
//
//  Glauber MC implementation
//
//  origin: PHOBOS experiment
//  alification: Mikolaj Krzewicki, Nikhef, mikolaj.krzewicki@cern.ch
//
////////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>
#include <TEllipse.h>
#include <TRandom.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TArray.h>

#include "AliGlauberNucleon.h"
#include "AliGlauberNucleus.h"
#include "AliGlauberMC.h"

ClassImp(AliGlauberMC)

//______________________________________________________________________________
AliGlauberMC::AliGlauberMC(Option_t* NA, Option_t* NB, Double_t xsect) :
  TNamed(),
  fANucleus(NA),
  fBNucleus(NB),
  fXSect(xsect),
  fNucleonsA(0),
  fNucleonsB(0),
  fAN(0),
  fBN(0),
  fnt(0),
  fMeanX2(0),
  fMeanY2(0),
  fMeanXY(0),
  fMeanXParts(0),
  fMeanYParts(0),
  fMeanXColl(0),
  fMeanYColl(0),
  fMeanX2Coll(0),
  fMeanY2Coll(0),
  fMeanXYColl(0),
  fMeanXSystem(0),
  fMeanYSystem(0),
  fMeanX_A(0),
  fMeanY_A(0),
  fMeanX_B(0),
  fMeanY_B(0),
  fB_MC(0),
  fEvents(0),
  fTotalEvents(0),
  fBMin(0.),
  fBMax(20.),
  fMaxNpartFound(0),
  fNpart(0),
  fNcoll(0),
  fSx2(0.),
  fSy2(0.),
  fSxy(0.),
  fSx2Coll(0.),
  fSy2Coll(0.),
  fSxyColl(0.),
  fX(0.13),
  fNpp(8.)
{
  fdNdEtaParam[0] = 8.0;
  fdNdEtaParam[1] = 0.13;
  fdNdEtaGBWParam[0] = 0.79;
  fdNdEtaGBWParam[1] = 0.288;
  fdNdEtaGBWParam[2] = 30.25;
  fdNdEtaNBDParam[0] = 3.0;
  fdNdEtaNBDParam[1] = 4.0;
  fdNdEtaNBDParam[2] = 0.13;
  fdNdEtaTwoNBDParam[0] = 3.0;
  fdNdEtaTwoNBDParam[1] = 4.0;
  fdNdEtaTwoNBDParam[2] = 2.0;
  fdNdEtaTwoNBDParam[3] = 11.0;
  fdNdEtaTwoNBDParam[4] = 0.4;
  fdNdEtaTwoNBDParam[5] = 0.13;

  SetName(Form("Glauber_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
  SetTitle(Form("Glauber %s+%s Version",fANucleus.GetName(),fBNucleus.GetName()));
}

//______________________________________________________________________________
AliGlauberMC::~AliGlauberMC()
{
  //dtor
  delete fnt;
}

//______________________________________________________________________________
AliGlauberMC::AliGlauberMC(const AliGlauberMC& in):
  TNamed(in),
  fANucleus(in.fANucleus),
  fBNucleus(in.fBNucleus),
  fXSect(in.fXSect),
  fNucleonsA(in.fNucleonsA),
  fNucleonsB(in.fNucleonsB),
  fAN(in.fAN),
  fBN(in.fBN),
  fnt(in.fnt),
  fMeanX2(in.fMeanX2),
  fMeanY2(in.fMeanY2),
  fMeanXY(in.fMeanXY),
  fMeanXParts(in.fMeanXParts),
  fMeanYParts(in.fMeanYParts),
  fMeanXColl(in.fMeanXColl),
  fMeanYColl(in.fMeanYColl),
  fMeanX2Coll(in.fMeanX2Coll),
  fMeanY2Coll(in.fMeanY2Coll),
  fMeanXYColl(in.fMeanXYColl),
  fMeanXSystem(in.fMeanXSystem),
  fMeanYSystem(in.fMeanYSystem),
  fMeanX_A(in.fMeanX_A),
  fMeanY_A(in.fMeanY_A),
  fMeanX_B(in.fMeanX_B),
  fMeanY_B(in.fMeanY_B),
  fB_MC(in.fB_MC),
  fEvents(in.fEvents),
  fTotalEvents(in.fTotalEvents),
  fBMin(in.fBMin),
  fBMax(in.fBMax),
  fMaxNpartFound(in.fMaxNpartFound),
  fNpart(in.fNpart),
  fNcoll(in.fNcoll),
  fSx2(in.fSx2),
  fSy2(in.fSy2),
  fSxy(in.fSxy),
  fSx2Coll(in.fSx2Coll),
  fSy2Coll(in.fSy2Coll),
  fSxyColl(in.fSxyColl),
  fX(in.fX),
  fNpp(in.fNpp)
{
  //copy ctor
  memcpy(fdNdEtaParam,in.fdNdEtaParam,2*sizeof(Double_t));
  memcpy(fdNdEtaGBWParam,in.fdNdEtaGBWParam,3*sizeof(Double_t));
  memcpy(fdNdEtaNBDParam,in.fdNdEtaNBDParam,3*sizeof(Double_t));
  memcpy(fdNdEtaTwoNBDParam,in.fdNdEtaTwoNBDParam,6*sizeof(Double_t));
}

//______________________________________________________________________________
AliGlauberMC& AliGlauberMC::operator=(const AliGlauberMC& in)
{
  //assignment
  fANucleus=in.fANucleus; 
  fBNucleus=in.fBNucleus; 
  fXSect=in.fXSect;    
  fNucleonsA=in.fNucleonsA;
  fNucleonsB=in.fNucleonsB;
  fAN=in.fAN;       
  fBN=in.fBN;       
  fnt=in.fnt;       
  fMeanX2=in.fMeanX2;   
  fMeanY2=in.fMeanY2;   
  fMeanXY=in.fMeanXY;   
  fMeanXParts=in.fMeanXParts;
  fMeanYParts=in.fMeanYParts;
  fMeanXColl=in.fMeanXColl;
  fMeanYColl=in.fMeanYColl;
  fMeanX2Coll=in.fMeanX2Coll;
  fMeanY2Coll=in.fMeanY2Coll;
  fMeanXYColl=in.fMeanXYColl;
  fMeanXSystem=in.fMeanXSystem;
  fMeanYSystem=in.fMeanYSystem;
  fMeanX_A=in.fMeanX_A;  
  fMeanY_A=in.fMeanY_A;  
  fMeanX_B=in.fMeanX_B;  
  fMeanY_B=in.fMeanY_B;  
  fB_MC=in.fB_MC;     
  fEvents=in.fEvents;   
  fTotalEvents=in.fTotalEvents;
  fBMin=in.fBMin;     
  fBMax=in.fBMax;     
  memcpy(fdNdEtaParam,in.fdNdEtaParam,2*sizeof(Double_t));
  memcpy(fdNdEtaGBWParam,in.fdNdEtaGBWParam,3*sizeof(Double_t));
  memcpy(fdNdEtaNBDParam,in.fdNdEtaNBDParam,3*sizeof(Double_t));
  memcpy(fdNdEtaTwoNBDParam,in.fdNdEtaTwoNBDParam,6*sizeof(Double_t));
  fMaxNpartFound=in.fMaxNpartFound;
  fNpart=in.fNpart;   
  fNcoll=in.fNcoll;    
  fSx2=in.fSx2;      
  fSy2=in.fSy2;      
  fSxy=in.fSxy;      
  fSx2Coll=in.fSx2Coll;  
  fSy2Coll=in.fSy2Coll;  
  fSxyColl=in.fSxyColl;  
  fX=in.fX;      
  fNpp=in.fNpp;       
  return *this;
}

//______________________________________________________________________________
Bool_t AliGlauberMC::CalcEvent(Double_t bgen)
{
  // prepare event
  fANucleus.ThrowNucleons(-bgen/2.);
  fNucleonsA = fANucleus.GetNucleons();
  fAN = fANucleus.GetN();
  for (Int_t i = 0; i<fAN; i++)
  {
    AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->At(i));
    nucleonA->SetInNucleusA();
  }
  fBNucleus.ThrowNucleons(bgen/2.);
  fNucleonsB = fBNucleus.GetNucleons();
  fBN = fBNucleus.GetN();
  for (Int_t i = 0; i<fBN; i++)
  {
    AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->At(i));
    nucleonB->SetInNucleusB();
  }

  // "ball" diameter = distance at which two balls interact
  Double_t d2 = (Double_t)fXSect/(TMath::Pi()*10); // in fm^2

  // for each of the A nucleons in nucleus B
  for (Int_t i = 0; i<fBN; i++)
  {
    AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->At(i));
    for (Int_t j = 0 ; j < fAN ; j++)
    {
      AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->At(j));
      Double_t dx = nucleonB->GetX()-nucleonA->GetX();
      Double_t dy = nucleonB->GetY()-nucleonA->GetY();
      Double_t dij = dx*dx+dy*dy;
      if (dij < d2)
      {
        nucleonB->Collide();
        nucleonA->Collide();
      }
    }
  }

  return CalcResults(bgen);
}

//______________________________________________________________________________
Bool_t AliGlauberMC::CalcResults(Double_t bgen)
{
  // calc results for the given event
  //return true if we have participants

  fNpart=0;
  fNcoll=0;
  fMeanX2=0;
  fMeanY2=0;
  fMeanXY=0;
  fMeanXParts=0;
  fMeanYParts=0;
  fMeanXColl=0;
  fMeanYColl=0;
  fMeanX2Coll=0;
  fMeanY2Coll=0;
  fMeanXYColl=0;
  fMeanXSystem=0;
  fMeanYSystem=0;
  fMeanX_A=0;
  fMeanY_A=0;
  fMeanX_B=0;
  fMeanY_B=0;

  for (Int_t i = 0; i<fAN; i++)
  {
    AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->At(i));
    Double_t xA=nucleonA->GetX();
    Double_t yA=nucleonA->GetY();
    fMeanXSystem  += xA;
    fMeanYSystem  += yA;
    fMeanX_A  += xA;
    fMeanY_A  += yA;

    if(nucleonA->IsWounded())
    {
      fNpart++;
      fMeanXParts  += xA;
      fMeanYParts  += yA;
      fMeanX2 += xA * xA;
      fMeanY2 += yA * yA;
      fMeanXY += xA * yA;
    }
  }

  for (Int_t i = 0; i<fBN; i++)
  {
    AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->At(i));
    Double_t xB=nucleonB->GetX();
    Double_t yB=nucleonB->GetY();
    fMeanXSystem  += xB;
    fMeanYSystem  += yB;
    fMeanX_B  += xB;
    fMeanY_B  += yB;

    if(nucleonB->IsWounded())
    {
      Int_t ncoll = nucleonB->GetNColl();
      fNpart++;
      fMeanXParts  += xB;
      fMeanXColl  += xB*ncoll;
      fMeanYParts  += yB;
      fMeanYColl += yB*ncoll;
      fMeanX2 += xB * xB;
      fMeanX2Coll += xB*xB*ncoll*ncoll;
      fMeanY2 += yB * yB;
      fMeanY2Coll += yB*yB*ncoll*ncoll;
      fMeanXY += xB * yB;
      fMeanXYColl += xB*yB*ncoll*ncoll;
      fNcoll += nucleonB->GetNColl();
    }
  }

  if (fNpart>0)
  {
    fMeanXParts /= fNpart;
    fMeanYParts /= fNpart;
    fMeanX2 /= fNpart;
    fMeanY2 /= fNpart;
    fMeanXY /= fNpart;
  }
  else
  {
    fMeanXParts = 0;
    fMeanYParts = 0;
    fMeanX2 = 0;
    fMeanY2 = 0;
    fMeanXY = 0;
  }

  if (fNcoll>0)
  {
    fMeanXColl /= fNcoll;
    fMeanYColl /= fNcoll;
    fMeanX2Coll /= fNcoll;
    fMeanY2Coll /= fNcoll;
    fMeanXYColl /= fNcoll;
  }
  else
  {
    fMeanXColl = 0;
    fMeanYColl = 0;
    fMeanX2Coll = 0;
    fMeanY2Coll = 0;
    fMeanXYColl = 0;
  }

  if(fAN+fBN>0)
  {
    fMeanXSystem /= (fAN + fBN);
    fMeanYSystem /= (fAN + fBN);
  }
  else
  {
    fMeanXSystem = 0;
    fMeanYSystem = 0;
  }
  if(fAN>0)
  {
    fMeanX_A /= fAN;
    fMeanY_A /= fAN;
  }
  else
  {
    fMeanX_A = 0;
    fMeanY_A = 0;
  }

  if(fBN>0)
  {
    fMeanX_B /= fBN;
    fMeanY_B /= fBN;
  }
  else
  {
    fMeanX_B = 0;
    fMeanY_B = 0;
  }

  fSx2=fMeanX2-(fMeanXParts*fMeanXParts);
  fSy2=fMeanY2-(fMeanYParts*fMeanYParts);
  fSxy=fMeanXY-fMeanXParts*fMeanYParts;

  fSx2Coll=fMeanX2Coll-(fMeanXColl*fMeanXColl);
  fSy2Coll=fMeanY2Coll-(fMeanYColl*fMeanYColl);
  fSxyColl=fMeanXYColl-fMeanXColl*fMeanYColl;

  fB_MC = bgen;

  fTotalEvents++;
  if (fNpart>0) fEvents++;

  if (fNpart==0) return kFALSE;
  if (fNpart > fMaxNpartFound) fMaxNpartFound = fNpart;

  return kTRUE;
}

//______________________________________________________________________________
void AliGlauberMC::Draw(Option_t* /*option*/)
{
  fANucleus.Draw(fXSect, 2);
  fBNucleus.Draw(fXSect, 4);

  TEllipse e;
  e.SetFillColor(0);
  e.SetLineColor(1);
  e.SetLineStyle(2);
  e.SetLineWidth(1);
  e.DrawEllipse(GetB()/2,0,fBNucleus.GetR(),fBNucleus.GetR(),0,360,0);
  e.DrawEllipse(-GetB()/2,0,fANucleus.GetR(),fANucleus.GetR(),0,360,0);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetTotXSect() const
{
  return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBMax*fBMax/100;
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetTotXSectErr() const
{
  return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) *
         TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

//______________________________________________________________________________
TObjArray *AliGlauberMC::GetNucleons()
{
  if(!fNucleonsA || !fNucleonsB) return 0;
  fNucleonsA->SetOwner(0);
  fNucleonsB->SetOwner(0);
  TObjArray *allnucleons=new TObjArray(fAN+fBN);
  allnucleons->SetOwner();
  for (Int_t i = 0; i<fAN; i++)
  {
    allnucleons->Add(fNucleonsA->At(i));
  }
  for (Int_t i = 0; i<fBN; i++)
  {
    allnucleons->Add(fNucleonsB->At(i));
  }
  return allnucleons;
}
//______________________________________________________________________________
Double_t AliGlauberMC::NegativeBinomialDistribution(Int_t x, Int_t k, Double_t nmean)
{
  if(k<=0)
  {
    cout << "Error, zero or negative K" << endl;
    return 0;
  }
  return (TMath::Binomial(x+k-1,x))
         *TMath::Power(((nmean/Double_t(k))/(1+nmean/Double_t(k))),Double_t(x))
         *(1/(TMath::Power((1+nmean/Double_t(k)),Double_t(k))));
}
//______________________________________________________________________________
Int_t AliGlauberMC::NegativeBinomialRandom(Int_t k, Double_t nmean)
//return random integer from a Negative Binomial Distribution
{
  static const Int_t fMaxPlot = 1000;
  Double_t array[fMaxPlot];
  array[0] = NegativeBinomialDistribution(0,k,nmean);
  for (Int_t i=1; i<fMaxPlot; i++)
  {
    array[i] = NegativeBinomialDistribution(i,k,nmean) + array[i-1];
  }
  Double_t r = gRandom->Uniform(0,1);
  return TMath::BinarySearch(fMaxPlot,array,r)+2;
}
//______________________________________________________________________________
Int_t AliGlauberMC::DoubleNegativeBinomialRandom( Int_t k,
                                                  Double_t nmean,
                                                  Int_t k2,
                                                  Double_t nmean2,
                                                  Double_t alpha )
{
  //return random integer from a Double Negative Binomial Distribution
  static const Int_t fMaxPlot = 1000;
  Double_t array[fMaxPlot];
  array[0] = alpha*NegativeBinomialDistribution(0,k,nmean)+(1-alpha)*NegativeBinomialDistribution(0,k2,nmean2);
  for (Int_t i=1; i<fMaxPlot; i++)
  {
    array[i] = alpha*NegativeBinomialDistribution(i,k,nmean)+(1-alpha)*NegativeBinomialDistribution(i,k2,nmean2) + array[i-1];
  }
  Double_t r = gRandom->Uniform(0,1);
  return TMath::BinarySearch(fMaxPlot,array,r)+2;
}


//______________________________________________________________________________
void AliGlauberMC::SetdNdEtaParam(Double_t nnp, Double_t x)
{
  // set parameters used for calculating multiplicity, see GetdNdEta() for commments
  fdNdEtaParam[0]=nnp;
  fdNdEtaParam[1]=x;
}

//______________________________________________________________________________
void AliGlauberMC::SetdNdEtaGBWParam(Double_t delta, Double_t lambda, Double_t snn)
{
  // set parameters used for calculating multiplicity see GetdNdEtaGBW() for commments
  fdNdEtaGBWParam[0]=delta;
  fdNdEtaGBWParam[1]=lambda;
  fdNdEtaGBWParam[2]=snn;
}
//______________________________________________________________________________
void AliGlauberMC::SetdNdEtaNBDParam(Double_t k, Double_t nmean, Double_t beta)
{
  // set parameters used for calculating multiplicity see GetdNdEtaNBD() for commments
  fdNdEtaNBDParam[0]=k;
  fdNdEtaNBDParam[1]=nmean;
  fdNdEtaNBDParam[2]=beta;
}
//______________________________________________________________________________
void AliGlauberMC::SetdNdEtaTwoNBDParam( Double_t k1, 
                                         Double_t nmean1, 
                                         Double_t k2, 
                                         Double_t nmean2, 
                                         Double_t alpha,
                                         Double_t beta)
{
  // set parameters used for calculating multiplicity see GetdNdEtaTwoNBD() for commments
  fdNdEtaTwoNBDParam[0]=k1;
  fdNdEtaTwoNBDParam[1]=nmean1;
  fdNdEtaTwoNBDParam[2]=k2;
  fdNdEtaTwoNBDParam[3]=nmean2;
  fdNdEtaTwoNBDParam[4]=alpha;
  fdNdEtaTwoNBDParam[5]=beta;
}
//______________________________________________________________________________
Double_t AliGlauberMC::GetdNdEta(Double_t nnp, Double_t x)
{
  //Get particle density per unit of rapidity
  //using two component model
  //Parameters: npp, x
  return nnp*((1.-x)*fNpart/2.+x*fNcoll);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetdNdEtaGBW( Double_t delta,
                                     Double_t lambda,
                                     Double_t  snn)
{
  //Get particle density per unit of rapidity
  //using the GBW model
  //Parameters: delta, lambda, snn
  return fNpart*0.47*TMath::Sqrt(TMath::Power(snn,lambda))
         * TMath::Power(fNpart,(1.-delta)/3./delta);
}
//_______________________________________________________________________________

Double_t AliGlauberMC::GetdNdEtaNBD ( Int_t k, Double_t nmean, Double_t beta)
{
  //Get particle density per unit of rapidity
  //using a aandomized number from a negative binomial distrubution
  //Parameters:   k  = related to distribition width
  //              nmean = mean of distribution
  //              beta = set contribution of participants / binary collisions to multiplicity
  Double_t mulnp=0.0;
  for(Int_t i = 0; i<fNpart; i++)
  {
    mulnp+=NegativeBinomialRandom(k,nmean);
  }
  Double_t mulnb=0.0;
  for(Int_t i = 0; i<fNcoll; i++)
  {
    mulnb+=NegativeBinomialRandom(k,nmean);
  }
  return (1-beta)*mulnp/2+beta*mulnb;
}
//______________________________________________________________________________

Double_t AliGlauberMC::GetdNdEtaTwoNBD ( Int_t k1,
                                         Double_t nmean1,
                                         Int_t k2,
                                         Double_t nmean2,
                                         Double_t alpha,
                                         Double_t beta )
{
  //Get particle density per unit of rapidity
  //using random numbers from two negative binomial distributions
  //Parameters:   k1 = related to distribition width of distribution 1
  //              nmean1 = mean of distribution 1
  //              k2 = related to distribition width of distribution 2
  //              nmean2 = mean of distribution 2
  //              alpha = set contributions of distrubitin 1 / distribution 2
  //              beta = set contribution of participants / binary collisions to multiplicity
  Double_t mulnp=0.0;
  for(Int_t i = 0; i<fNpart; i++)
  {
    mulnp+=DoubleNegativeBinomialRandom(k1,nmean1,k2,nmean2,alpha);
  }
  Double_t mulnb=0.0;
  for(Int_t i = 0; i<fNcoll; i++)
  {
    mulnb+=DoubleNegativeBinomialRandom(k1,nmean1,k2,nmean2,alpha);
  }
  Double_t mul=(1-beta)*mulnp/2+beta*mulnb;
  return mul;

}
//______________________________________________________________________________
Double_t AliGlauberMC::GetEccentricityPart() const
{
  //get participant eccentricity of participants
  if (fNpart<2) return 0.0;
  return (TMath::Sqrt((fSy2-fSx2)*(fSy2-fSx2)+4*fSxy*fSxy)/(fSy2+fSx2));
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEccentricityPartColl() const
{
  //get participant eccentricity of binary collisions
  if (fNcoll<2) return 0.0;
  return (TMath::Sqrt((fSy2Coll-fSx2Coll)*(fSy2Coll-fSx2Coll)+4*fSxyColl*fSxyColl)/(fSy2Coll+fSx2Coll));
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEccentricity() const
{
  //get standard eccentricity of participants
  return ((fSy2-fSx2)/(fSy2+fSx2));
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEccentricityColl() const
{
  //get standard eccentricity of binary collisions
  return ((fSy2Coll-fSx2Coll)/(fSy2Coll+fSx2Coll));
}
//______________________________________________________________________________
Bool_t AliGlauberMC::NextEvent(Double_t bgen)
{
  //Make a new event
  Int_t nAttempts = 10; // set indices, max attempts and max comparisons
  Bool_t succes = kFALSE;
  for(Int_t j=0; j<nAttempts; j++)
  {
    if(bgen<0||!succes) //get impactparameter
    {
      bgen = TMath::Sqrt((fBMax*fBMax-fBMin*fBMin)*gRandom->Rndm()+fBMin*fBMin);
    }
    if ( succes=CalcEvent(bgen) ) break; //ends if we have particparts
  }
  return succes;
}

//______________________________________________________________________________
void AliGlauberMC::Run(Int_t nevents)
{
  cout << "Generating " << nevents << " events..." << endl;
  TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
  TString title(Form("%s + %s (x-sect = %d mb)",fANucleus.GetName(),fBNucleus.GetName(),(Int_t) fXSect));
  if (fnt == 0)
  {
    fnt = new TNtuple(name,title,
                      "Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB:VarE:VarEColl:VarEPart:VarEPartColl:dNdEta:dNdEtaGBW:dNdEtaNBD:dNdEtaTwoNBD");
    fnt->SetDirectory(0);
  }
  Int_t q = 0;
  Int_t u = 0;
  for (Int_t i = 0; i<nevents; i++)
  {

    if(!NextEvent())
    {
      u++;
      continue;
    }

    q++;
    Float_t v[25];
    v[0]  = GetNpart();
    v[1]  = GetNcoll();
    v[2]  = fB_MC;
    v[3]  = fMeanXParts;
    v[4]  = fMeanYParts;
    v[5]  = fMeanX2;
    v[6]  = fMeanY2;
    v[7]  = fMeanXY;
    v[8]  = fSx2;
    v[9]  = fSy2;
    v[10] = fSxy;
    v[11] = fMeanXSystem;
    v[12] = fMeanYSystem;
    v[13] = fMeanX_A;
    v[14] = fMeanY_A;
    v[15] = fMeanX_B;
    v[16] = fMeanY_B;
    v[17] = GetEccentricity();
    v[18] = GetEccentricityColl();
    v[19] = GetEccentricityPart();
    v[20] = GetEccentricityPartColl();
    v[21] = GetdNdEta( fdNdEtaParam[0],fdNdEtaParam[1] );
    v[22] = GetdNdEtaGBW( fdNdEtaGBWParam[0],fdNdEtaGBWParam[1],fdNdEtaGBWParam[2] );
    v[23] = GetdNdEtaNBD( TMath::Nint(fdNdEtaNBDParam[0]),
                          fdNdEtaNBDParam[1],
                          fdNdEtaNBDParam[2] );
    v[24] = GetdNdEtaTwoNBD( TMath::Nint(fdNdEtaTwoNBDParam[0]),
                             fdNdEtaTwoNBDParam[1],
                             TMath::Nint(fdNdEtaTwoNBDParam[2]),
                             fdNdEtaTwoNBDParam[3],
                             fdNdEtaTwoNBDParam[4],
                             fdNdEtaTwoNBDParam[5] );
    fnt->Fill(v);

    if ((i%50)==0) std::cout << "Generating Event # " << i << "... \r" << flush;
  }
  std::cout << "Generating Event # " << nevents << "... \r" << endl << "Done! Succesfull events:  " << q << "  discarded events:  " << u <<"."<< endl;
}

//---------------------------------------------------------------------------------
void AliGlauberMC::runAndSaveNtuple( Int_t n,
                                     Option_t *sysA,
                                     Option_t *sysB,
                                     Double_t signn,
                                     Double_t mind,
                                     const char *fname)
{
  AliGlauberMC mcg(sysA,sysB,signn);
  mcg.SetMinDistance(mind);
  mcg.Run(n);
  TNtuple  *nt=mcg.GetNtuple();
  TFile out(fname,"recreate",fname,9);
  if(nt) nt->Write();
  out.Close();
}

//---------------------------------------------------------------------------------
void AliGlauberMC::runAndSaveNucleons( Int_t n,
                                       Option_t *sysA,
                                       Option_t *sysB,
                                       Double_t signn,
                                       Double_t mind,
                                       Bool_t verbose,
                                       const char *fname)
{
  AliGlauberMC mcg(sysA,sysB,signn);
  mcg.SetMinDistance(mind);
  TFile *out=0;
  if(fname)
    out=new TFile(fname,"recreate",fname,9);

  for(Int_t ievent=0; ievent<n; ievent++)
  {

    //get an event with at least one collision
    while(!mcg.NextEvent()) {}

    //access, save and (if wanted) print out nucleons
    TObjArray* nucleons=mcg.GetNucleons();
    if(!nucleons) continue;
    if(out)
      nucleons->Write(Form("nucleonarray%d",ievent),TObject::kSingleKey);

    if(verbose)
    {
      cout<<endl<<endl<<"EVENT NO: "<<ievent<<endl;
      cout<<"B = "<<mcg.GetB()<<"  Npart = "<<mcg.GetNpart()<<endl<<endl;
      printf("Nucleus\t X\t Y\t Z\tNcoll\n");
      Int_t nNucls=nucleons->GetEntries();
      for(Int_t iNucl=0; iNucl<nNucls; iNucl++)
      {
        AliGlauberNucleon *nucl=(AliGlauberNucleon *)nucleons->At(iNucl);
        Char_t nucleus='A';
        if(nucl->IsInNucleusB()) nucleus='B';
        Double_t x=nucl->GetX();
        Double_t y=nucl->GetY();
        Double_t z=nucl->GetZ();
        Int_t ncoll=nucl->GetNColl();
        printf("   %c\t%2.2f\t%2.2f\t%2.2f\t%3d\n",nucleus,x,y,z,ncoll);
      }
    }
  }
  if(out) delete out;
}

//---------------------------------------------------------------------------------
void AliGlauberMC::Reset()
{
  //delete the ntuple
  delete fnt; fnt=NULL;
}
