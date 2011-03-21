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
//  update:      You Zhou, Nikhef, yzhou@nikhef.nl
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
  fMeanXA(0),
  fMeanYA(0),
  fMeanXB(0),
  fMeanYB(0),
  fBMC(0),
  fEvents(0),
  fTotalEvents(0),
  fBMin(0.),
  fBMax(20.),
  fMultType(kNBDSV),
  fMaxNpartFound(0),
  fNpart(0),
  fNcoll(0),
  fMeanr2(0),
  fMeanr3(0),
  fMeanr4(0),
  fMeanr5(0),
  fMeanr2Cos2Phi(0),
  fMeanr2Sin2Phi(0),
  fMeanr2Cos3Phi(0),
  fMeanr2Sin3Phi(0),
  fMeanr2Cos4Phi(0),
  fMeanr2Sin4Phi(0),
  fMeanr2Cos5Phi(0),
  fMeanr2Sin5Phi(0),
  fMeanr3Cos3Phi(0),
  fMeanr3Sin3Phi(0),
  fMeanr4Cos4Phi(0),
  fMeanr4Sin4Phi(0),
  fMeanr5Cos5Phi(0),
  fMeanr5Sin5Phi(0),
  fMeanr2Coll(0),
  fMeanr3Coll(0),
  fMeanr4Coll(0),
  fMeanr5Coll(0),
  fMeanr2Cos2PhiColl(0),
  fMeanr2Sin2PhiColl(0),
  fMeanr2Cos3PhiColl(0),
  fMeanr2Sin3PhiColl(0),
  fMeanr2Cos4PhiColl(0),
  fMeanr2Sin4PhiColl(0),
  fMeanr2Cos5PhiColl(0),
  fMeanr2Sin5PhiColl(0),
  fMeanr3Cos3PhiColl(0),
  fMeanr3Sin3PhiColl(0),
  fMeanr4Cos4PhiColl(0),
  fMeanr4Sin4PhiColl(0),
  fMeanr5Cos5PhiColl(0),
  fMeanr5Sin5PhiColl(0),
  fSx2(0.),
  fSy2(0.),
  fSxy(0.),
  fSx2Coll(0.),
  fSy2Coll(0.),
  fSxyColl(0.),
  fX(0.13),
  fNpp(8.),
  fDoPartProd(kFALSE)
{
  //ctor
  for (UInt_t i=0; i<(sizeof(fdNdEtaParam)/sizeof(fdNdEtaParam[0])); i++)
  {
    fdNdEtaParam[i]=0.0;
  }

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
  fMeanXA(in.fMeanXA),
  fMeanYA(in.fMeanYA),
  fMeanXB(in.fMeanXB),
  fMeanYB(in.fMeanYB),
  fBMC(in.fBMC),
  fEvents(in.fEvents),
  fTotalEvents(in.fTotalEvents),
  fBMin(in.fBMin),
  fBMax(in.fBMax),
  fMultType(in.fMultType),
  fMaxNpartFound(in.fMaxNpartFound),
  fNpart(in.fNpart),
  fNcoll(in.fNcoll),
  fMeanr2(in.fMeanr2),
  fMeanr3(in.fMeanr3),
  fMeanr4(in.fMeanr4),
  fMeanr5(in.fMeanr5),
  fMeanr2Cos2Phi(in.fMeanr2Cos2Phi),
  fMeanr2Sin2Phi(in.fMeanr2Sin2Phi),
  fMeanr2Cos3Phi(in.fMeanr2Cos3Phi),
  fMeanr2Sin3Phi(in.fMeanr2Sin3Phi),
  fMeanr2Cos4Phi(in.fMeanr2Cos4Phi),
  fMeanr2Sin4Phi(in.fMeanr2Sin4Phi),
  fMeanr2Cos5Phi(in.fMeanr2Cos5Phi),
  fMeanr2Sin5Phi(in.fMeanr2Sin5Phi),
  fMeanr3Cos3Phi(in.fMeanr3Cos3Phi),
  fMeanr3Sin3Phi(in.fMeanr3Sin3Phi),
  fMeanr4Cos4Phi(in.fMeanr4Cos4Phi),
  fMeanr4Sin4Phi(in.fMeanr4Sin4Phi),
  fMeanr5Cos5Phi(in.fMeanr5Cos5Phi),
  fMeanr5Sin5Phi(in.fMeanr5Sin5Phi),
  fMeanr2Coll(in.fMeanr2Coll),
  fMeanr3Coll(in.fMeanr3Coll),
  fMeanr4Coll(in.fMeanr4Coll),
  fMeanr5Coll(in.fMeanr5Coll),
  fMeanr2Cos2PhiColl(in.fMeanr2Cos2PhiColl),
  fMeanr2Sin2PhiColl(in.fMeanr2Sin2PhiColl),
  fMeanr2Cos3PhiColl(in.fMeanr2Cos3PhiColl),
  fMeanr2Sin3PhiColl(in.fMeanr2Sin3PhiColl),
  fMeanr2Cos4PhiColl(in.fMeanr2Cos4PhiColl),
  fMeanr2Sin4PhiColl(in.fMeanr2Sin4PhiColl),
  fMeanr2Cos5PhiColl(in.fMeanr2Cos5PhiColl),
  fMeanr2Sin5PhiColl(in.fMeanr2Sin5PhiColl),
  fMeanr3Cos3PhiColl(in.fMeanr3Cos3PhiColl),
  fMeanr3Sin3PhiColl(in.fMeanr3Sin3PhiColl),
  fMeanr4Cos4PhiColl(in.fMeanr4Cos4PhiColl),
  fMeanr4Sin4PhiColl(in.fMeanr4Sin4PhiColl),
  fMeanr5Cos5PhiColl(in.fMeanr5Cos5PhiColl),
  fMeanr5Sin5PhiColl(in.fMeanr5Sin5PhiColl),
  fSx2(in.fSx2),
  fSy2(in.fSy2),
  fSxy(in.fSxy),
  fSx2Coll(in.fSx2Coll),
  fSy2Coll(in.fSy2Coll),
  fSxyColl(in.fSxyColl),
  fX(in.fX),
  fNpp(in.fNpp),
  fDoPartProd(kFALSE)
{
  //copy ctor
  memcpy(fdNdEtaParam,in.fdNdEtaParam,sizeof(fdNdEtaParam));
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
  fMeanr2=in.fMeanr2;
  fMeanr3=in.fMeanr3;
  fMeanr4=in.fMeanr4;
  fMeanr5=in.fMeanr5;
  fMeanXParts=in.fMeanXParts;
  fMeanYParts=in.fMeanYParts;
  fMeanXColl=in.fMeanXColl;
  fMeanYColl=in.fMeanYColl;
  fMeanX2Coll=in.fMeanX2Coll;
  fMeanY2Coll=in.fMeanY2Coll;
  fMeanXYColl=in.fMeanXYColl;
  fMeanr2Coll=in.fMeanr2Coll;
  fMeanr3Coll=in.fMeanr3Coll;
  fMeanr4Coll=in.fMeanr4Coll;
  fMeanr5Coll=in.fMeanr5Coll;
  fMeanXSystem=in.fMeanXSystem;
  fMeanYSystem=in.fMeanYSystem;
  fMeanXA=in.fMeanXA;
  fMeanYA=in.fMeanYA;
  fMeanXB=in.fMeanXB;
  fMeanYB=in.fMeanYB;
  fBMC=in.fBMC;
  fEvents=in.fEvents;
  fTotalEvents=in.fTotalEvents;
  fBMin=in.fBMin;
  fBMax=in.fBMax;
  fMultType=in.fMultType,
  memcpy(fdNdEtaParam,in.fdNdEtaParam,sizeof(fdNdEtaParam));
  fMaxNpartFound=in.fMaxNpartFound;
  fNpart=in.fNpart;
  fNcoll=in.fNcoll;
  fMeanr2Cos2Phi=in.fMeanr2Cos2Phi;
  fMeanr2Sin2Phi=in.fMeanr2Sin2Phi;
  fMeanr2Cos3Phi=in.fMeanr2Cos3Phi;
  fMeanr2Sin3Phi=in.fMeanr2Sin3Phi;
  fMeanr2Cos4Phi=in.fMeanr2Cos4Phi;
  fMeanr2Sin4Phi=in.fMeanr2Sin4Phi;
  fMeanr2Cos5Phi=in.fMeanr2Cos5Phi;
  fMeanr2Sin5Phi=in.fMeanr2Sin5Phi;
  fMeanr3Cos3Phi=in.fMeanr3Cos3Phi;
  fMeanr3Sin3Phi=in.fMeanr3Sin3Phi;
  fMeanr4Cos4Phi=in.fMeanr4Cos4Phi;
  fMeanr4Sin4Phi=in.fMeanr4Sin4Phi;
  fMeanr5Cos5Phi=in.fMeanr5Cos5Phi;
  fMeanr5Sin5Phi=in.fMeanr5Sin5Phi;
  fMeanr2Cos2PhiColl=in.fMeanr2Cos2PhiColl;
  fMeanr2Sin2PhiColl=in.fMeanr2Sin2PhiColl;
  fMeanr2Cos3PhiColl=in.fMeanr2Cos3PhiColl;
  fMeanr2Sin3PhiColl=in.fMeanr2Sin3PhiColl;
  fMeanr2Cos4PhiColl=in.fMeanr2Cos4PhiColl;
  fMeanr2Sin4PhiColl=in.fMeanr2Sin4PhiColl;
  fMeanr2Cos5PhiColl=in.fMeanr2Cos5PhiColl;
  fMeanr2Sin5PhiColl=in.fMeanr2Sin5PhiColl;
  fMeanr3Cos3PhiColl=in.fMeanr3Cos3PhiColl;
  fMeanr3Sin3PhiColl=in.fMeanr3Sin3PhiColl;
  fMeanr4Cos4PhiColl=in.fMeanr4Cos4PhiColl;
  fMeanr4Sin4PhiColl=in.fMeanr4Sin4PhiColl;
  fMeanr5Cos5PhiColl=in.fMeanr5Cos5PhiColl;
  fMeanr5Sin5PhiColl=in.fMeanr5Sin5PhiColl;
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
    AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->UncheckedAt(i));
    nucleonA->SetInNucleusA();
  }
  fBNucleus.ThrowNucleons(bgen/2.);
  fNucleonsB = fBNucleus.GetNucleons();
  fBN = fBNucleus.GetN();
  for (Int_t i = 0; i<fBN; i++)
  {
    AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->UncheckedAt(i));
    nucleonB->SetInNucleusB();
  }

  // "ball" diameter = distance at which two balls interact
  Double_t d2 = (Double_t)fXSect/(TMath::Pi()*10); // in fm^2

  // for each of the A nucleons in nucleus B
  for (Int_t i = 0; i<fBN; i++)
  {
    AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->UncheckedAt(i));
    for (Int_t j = 0 ; j < fAN ; j++)
    {
      AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->UncheckedAt(j));
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
  fMeanX2=0.;
  fMeanY2=0.;
  fMeanXY=0.;
  fMeanXParts=0.;
  fMeanYParts=0.;
  fMeanXColl=0.;
  fMeanYColl=0.;
  fMeanX2Coll=0.;
  fMeanY2Coll=0.;
  fMeanXYColl=0.;
  fMeanXSystem=0.;
  fMeanYSystem=0.;
  fMeanXA=0.;
  fMeanYA=0.;
  fMeanXB=0.;
  fMeanYB=0.;
  fMeanr2=0.;
  fMeanr3=0.;
  fMeanr4=0.;
  fMeanr5=0.;
  fMeanr2Cos2Phi=0.;
  fMeanr2Sin2Phi=0.;
  fMeanr2Cos3Phi=0.;
  fMeanr2Sin3Phi=0.;
  fMeanr2Cos4Phi=0.;
  fMeanr2Sin4Phi=0.;
  fMeanr2Cos5Phi=0.;
  fMeanr2Sin5Phi=0.;
  fMeanr3Cos3Phi=0.;
  fMeanr3Sin3Phi=0.;
  fMeanr4Cos4Phi=0.;
  fMeanr4Sin4Phi=0.;
  fMeanr5Cos5Phi=0.;
  fMeanr5Sin5Phi=0.;
  fMeanr2Coll=0.;
  fMeanr3Coll=0.;
  fMeanr4Coll=0.;
  fMeanr5Coll=0.;
  fMeanr2Cos2PhiColl=0.;
  fMeanr2Sin2PhiColl=0.;
  fMeanr2Cos3PhiColl=0.;
  fMeanr2Sin3PhiColl=0.;
  fMeanr2Cos4PhiColl=0.;
  fMeanr2Sin4PhiColl=0.;
  fMeanr2Cos5PhiColl=0.;
  fMeanr2Sin5PhiColl=0.;
  fMeanr3Cos3PhiColl=0.;
  fMeanr3Sin3PhiColl=0.;
  fMeanr4Cos4PhiColl=0.;
  fMeanr4Sin4PhiColl=0.;
  fMeanr5Cos5PhiColl=0.;
  fMeanr5Sin5PhiColl=0.;

  for (Int_t i = 0; i<fAN; i++)
  {
    AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->UncheckedAt(i));
    Double_t xA = nucleonA->GetX();
    Double_t yA = nucleonA->GetY();
    Double_t r2A = xA*xA+yA*yA;
    Double_t rA = TMath::Sqrt(r2A);
    Double_t r3A = r2A*rA;
    Double_t r4A = r3A*rA;
    Double_t r5A = r4A*rA;
    Double_t phiA = TMath::ATan2(yA,xA);
    Double_t sin2PhiA = TMath::Sin(2.*phiA);
    Double_t cos2PhiA = TMath::Cos(2.*phiA);
    Double_t sin3PhiA = TMath::Sin(3.*phiA);
    Double_t cos3PhiA = TMath::Cos(3.*phiA);
    Double_t sin4PhiA = TMath::Sin(4.*phiA);
    Double_t cos4PhiA = TMath::Cos(4.*phiA);
    Double_t sin5PhiA = TMath::Sin(5.*phiA);
    Double_t cos5PhiA = TMath::Cos(5.*phiA);
    //printf("x: %.3f, y:%.3f, r:%.3f, phi:%.3f, sp:%.6f,cp:%.6f\n",xA,yA,TMath::Sqrt(r2A),PhiA,Sin2PhiA, Cos2PhiA);
    fMeanXSystem  += xA;
    fMeanYSystem  += yA;
    fMeanXA  += xA;
    fMeanYA  += yA;

    if(nucleonA->IsWounded())
    {
      fNpart++;
      fMeanXParts  += xA;
      fMeanYParts  += yA;
      fMeanX2 += xA * xA;
      fMeanY2 += yA * yA;
      fMeanXY += xA * yA;
      fMeanr2 += r2A;
      fMeanr3 += r3A;
      fMeanr4 += r4A;
      fMeanr5 += r5A;
      fMeanr2Cos2Phi += r2A*cos2PhiA;
      fMeanr2Sin2Phi += r2A*sin2PhiA;
      fMeanr2Cos3Phi += r2A*cos3PhiA;
      fMeanr2Sin3Phi += r2A*sin3PhiA;
      fMeanr2Cos4Phi += r2A*cos4PhiA;
      fMeanr2Sin4Phi += r2A*sin4PhiA;
      fMeanr2Cos5Phi += r2A*cos5PhiA;
      fMeanr2Sin5Phi += r2A*sin5PhiA;
      fMeanr3Cos3Phi += r3A*cos3PhiA;
      fMeanr3Sin3Phi += r3A*sin3PhiA;
      fMeanr4Cos4Phi += r4A*cos4PhiA;
      fMeanr4Sin4Phi += r4A*sin4PhiA;
      fMeanr5Cos5Phi += r5A*cos5PhiA;
      fMeanr5Sin5Phi += r5A*sin5PhiA;
    }
  }

  for (Int_t i = 0; i<fBN; i++)
  {
    AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->UncheckedAt(i));
    Double_t xB=nucleonB->GetX();
    Double_t yB=nucleonB->GetY();
    Double_t r2B=xB*xB+yB*yB;
    Double_t rB = TMath::Sqrt(r2B);
    Double_t r3B = r2B*rB;
    Double_t r4B = r3B*rB;
    Double_t r5B = r4B*rB;
    Double_t phiB = TMath::ATan2(yB,xB);
    Double_t sin2PhiB = TMath::Sin(2*phiB);
    Double_t cos2PhiB = TMath::Cos(2*phiB);
    Double_t sin3PhiB = TMath::Sin(3*phiB);
    Double_t cos3PhiB = TMath::Cos(3*phiB);
    Double_t sin4PhiB = TMath::Sin(4*phiB);
    Double_t cos4PhiB = TMath::Cos(4*phiB);
    Double_t sin5PhiB = TMath::Sin(5*phiB);
    Double_t cos5PhiB = TMath::Cos(5*phiB);
    fMeanXSystem  += xB;
    fMeanYSystem  += yB;
    fMeanXB  += xB;
    fMeanYB  += yB;

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
      fMeanr2 += r2B;
      fMeanr3 += r3B;
      fMeanr4 += r4B;
      fMeanr5 += r5B;
      fMeanr2Cos2Phi += r2B*cos2PhiB;
      fMeanr2Sin2Phi += r2B*sin2PhiB;
      fMeanr2Cos3Phi += r2B*cos3PhiB;
      fMeanr2Sin3Phi += r2B*sin3PhiB;
      fMeanr2Cos4Phi += r2B*cos4PhiB;
      fMeanr2Sin4Phi += r2B*sin4PhiB;
      fMeanr2Cos5Phi += r2B*cos5PhiB;
      fMeanr2Sin5Phi += r2B*sin5PhiB;
      fMeanr3Cos3Phi += r3B*cos3PhiB;
      fMeanr3Sin3Phi += r3B*sin3PhiB;
      fMeanr4Cos4Phi += r4B*cos4PhiB;
      fMeanr4Sin4Phi += r4B*sin4PhiB;
      fMeanr5Cos5Phi += r5B*cos5PhiB;
      fMeanr5Sin5Phi += r5B*sin5PhiB;
      fMeanr2Coll += r2B*ncoll*ncoll;
      fMeanr3Coll += r3B*ncoll*ncoll;
      fMeanr4Coll += r4B*ncoll*ncoll;
      fMeanr5Coll += r5B*ncoll*ncoll;
      fMeanr2Cos2PhiColl += r2B*cos2PhiB*ncoll*ncoll;
      fMeanr2Sin2PhiColl += r2B*sin2PhiB*ncoll*ncoll;
      fMeanr2Cos3PhiColl += r2B*cos3PhiB*ncoll*ncoll;
      fMeanr2Sin3PhiColl += r2B*sin3PhiB*ncoll*ncoll;
      fMeanr2Cos4PhiColl += r2B*cos4PhiB*ncoll*ncoll;
      fMeanr2Sin4PhiColl += r2B*sin4PhiB*ncoll*ncoll;
      fMeanr2Cos5PhiColl += r2B*cos5PhiB*ncoll*ncoll;
      fMeanr2Sin5PhiColl += r2B*sin5PhiB*ncoll*ncoll;
      fMeanr3Cos3PhiColl += r3B*cos3PhiB*ncoll*ncoll;
      fMeanr3Sin3PhiColl += r3B*sin3PhiB*ncoll*ncoll;
      fMeanr4Cos4PhiColl += r4B*cos4PhiB*ncoll*ncoll;
      fMeanr4Sin4PhiColl += r4B*sin4PhiB*ncoll*ncoll;
      fMeanr5Cos5PhiColl += r5B*cos5PhiB*ncoll*ncoll;
      fMeanr5Sin5PhiColl += r5B*sin5PhiB*ncoll*ncoll;
    }
  }

  if (fNpart>0)
  {
    fMeanXParts /= fNpart;
    fMeanYParts /= fNpart;
    fMeanX2 /= fNpart;
    fMeanY2 /= fNpart;
    fMeanXY /= fNpart;
    fMeanr2 /= fNpart;
    fMeanr3 /= fNpart;
    fMeanr4 /= fNpart;
    fMeanr5 /= fNpart;
    fMeanr2Cos2Phi /= fNpart;
    fMeanr2Sin2Phi /= fNpart;
    fMeanr2Cos3Phi /= fNpart;
    fMeanr2Sin3Phi /= fNpart;
    fMeanr2Cos4Phi /= fNpart;
    fMeanr2Sin4Phi /= fNpart;
    fMeanr2Cos5Phi /= fNpart;
    fMeanr2Sin5Phi /= fNpart;
    fMeanr3Cos3Phi /= fNpart;
    fMeanr3Sin3Phi /= fNpart;
    fMeanr4Cos4Phi /= fNpart;
    fMeanr4Sin4Phi /= fNpart;
    fMeanr5Cos5Phi /= fNpart;
    fMeanr5Sin5Phi /= fNpart;
  }
  else
  {
    fMeanXParts = 0;
    fMeanYParts = 0;
    fMeanX2 = 0;
    fMeanY2 = 0;
    fMeanXY = 0;
    fMeanr2 = 0;
    fMeanr3 = 0;
    fMeanr4 = 0;
    fMeanr5 = 0;
    fMeanr2Cos2Phi = 0;
    fMeanr2Sin2Phi = 0;
    fMeanr2Cos3Phi = 0;
    fMeanr2Sin3Phi = 0;
    fMeanr2Cos4Phi = 0;
    fMeanr2Sin4Phi = 0;
    fMeanr2Cos5Phi = 0;
    fMeanr2Sin5Phi = 0;
    fMeanr3Cos3Phi = 0;
    fMeanr3Sin3Phi = 0;
    fMeanr4Cos4Phi = 0;
    fMeanr4Sin4Phi = 0;
    fMeanr5Cos5Phi = 0;
    fMeanr5Sin5Phi = 0;
  }

  if (fNcoll>0)
  {
    fMeanXColl /= fNcoll;
    fMeanYColl /= fNcoll;
    fMeanX2Coll /= fNcoll;
    fMeanY2Coll /= fNcoll;
    fMeanXYColl /= fNcoll;
    fMeanr2Coll /= fNcoll;
    fMeanr3Coll /= fNcoll;
    fMeanr4Coll /= fNcoll;
    fMeanr5Coll /= fNcoll;
    fMeanr2Cos2PhiColl /= fNcoll;
    fMeanr2Sin2PhiColl /= fNcoll;
    fMeanr2Cos3PhiColl /= fNcoll;
    fMeanr2Sin3PhiColl /= fNcoll;
    fMeanr2Cos4PhiColl /= fNcoll;
    fMeanr2Sin4PhiColl /= fNcoll;
    fMeanr2Cos5PhiColl /= fNcoll;
    fMeanr2Sin5PhiColl /= fNcoll;
    fMeanr3Cos3PhiColl /= fNcoll;
    fMeanr3Sin3PhiColl /= fNcoll;
    fMeanr4Cos4PhiColl /= fNcoll;
    fMeanr4Sin4PhiColl /= fNcoll;
    fMeanr5Cos5PhiColl /= fNcoll;
    fMeanr5Sin5PhiColl /= fNcoll;
  }
  else
  {
    fMeanXColl = 0;
    fMeanYColl = 0;
    fMeanX2Coll = 0;
    fMeanY2Coll = 0;
    fMeanXYColl = 0;
    fMeanr2Cos2PhiColl =0;
    fMeanr2Sin2PhiColl =0;
    fMeanr2Cos3PhiColl =0;
    fMeanr2Sin3PhiColl =0;
    fMeanr2Cos4PhiColl =0;
    fMeanr2Sin4PhiColl =0;
    fMeanr2Cos5PhiColl =0;
    fMeanr2Sin5PhiColl =0;
    fMeanr3Cos3PhiColl =0;
    fMeanr3Sin3PhiColl =0;
    fMeanr4Cos4PhiColl =0;
    fMeanr4Sin4PhiColl =0;
    fMeanr5Cos5PhiColl =0;
    fMeanr5Sin5PhiColl =0;
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
    fMeanXA /= fAN;
    fMeanYA /= fAN;
  }
  else
  {
    fMeanXA = 0;
    fMeanYA = 0;
  }

  if(fBN>0)
  {
    fMeanXB /= fBN;
    fMeanYB /= fBN;
  }
  else
  {
    fMeanXB = 0;
    fMeanYB = 0;
  }

  fSx2=fMeanX2-(fMeanXParts*fMeanXParts);
  fSy2=fMeanY2-(fMeanYParts*fMeanYParts);
  fSxy=fMeanXY-fMeanXParts*fMeanYParts;
  fSx2Coll=fMeanX2Coll-(fMeanXColl*fMeanXColl);
  fSy2Coll=fMeanY2Coll-(fMeanYColl*fMeanYColl);
  fSxyColl=fMeanXYColl-fMeanXColl*fMeanYColl;
  fBMC = bgen;
  fTotalEvents++;
  if (fNpart>0) fEvents++;
  if (fNpart==0) return kFALSE;
  if (fNpart > fMaxNpartFound) fMaxNpartFound = fNpart;

  return kTRUE;
}

//______________________________________________________________________________
void AliGlauberMC::Draw(Option_t* /*option*/)
{
  //draw method
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
  //total xsection
  return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBMax*fBMax/100;
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetTotXSectErr() const
{
  //total xsection error
  return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) *
         TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

//______________________________________________________________________________
TObjArray *AliGlauberMC::GetNucleons()
{
  //get array of nucleons
  if(!fNucleonsA || !fNucleonsB) return 0;
  fNucleonsA->SetOwner(0);
  fNucleonsB->SetOwner(0);
  TObjArray *allnucleons=new TObjArray(fAN+fBN);
  allnucleons->SetOwner();
  for (Int_t i = 0; i<fAN; i++)
  {
    allnucleons->Add(fNucleonsA->UncheckedAt(i));
  }
  for (Int_t i = 0; i<fBN; i++)
  {
    allnucleons->Add(fNucleonsB->UncheckedAt(i));
  }
  return allnucleons;
}

//______________________________________________________________________________
Double_t AliGlauberMC::NegativeBinomialDistribution(Int_t x, Int_t k, Double_t nmean)
{
  //produce a number distributed acc. neg.bin.dist
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
Int_t AliGlauberMC::NegativeBinomialRandom(Int_t k, Double_t nmean) const
{
  //return random integer from a Negative Binomial Distribution
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
Int_t AliGlauberMC::NegativeBinomialRandomSV(Double_t k, Double_t nbar) const
{
  // negative binomial distribution generator, S. Voloshin, 09-May-2007
  Double_t sum=0.;
  Int_t i=0;
  Double_t ran=gRandom->Rndm();
  Double_t trm=1./pow(1.+nbar/k,k);
  if (trm==0.)
  {
    cout<<"NBD overflow"<<"  nbar="<<nbar<<"   k="<<k<<endl;
    return -1;
  }
  for(i=0; i<2000 && sum<ran ; i++)
  {
    sum += trm;
    trm *= (k+i)/(i+1.)*(nbar/(nbar+k));
  }
  return i-1;
}

//______________________________________________________________________________
Int_t AliGlauberMC::DoubleNegativeBinomialRandom( Int_t k,
    Double_t nmean,
    Int_t k2,
    Double_t nmean2,
    Double_t alpha ) const
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
Double_t AliGlauberMC::GetdNdEta() const
{
  switch (fMultType)
  {
    case kSimple:
      return GetdNdEtaSimple(fdNdEtaParam);
    case kNBD:
      return GetdNdEtaNBD(fdNdEtaParam);
    case kNBDSV:
      return GetdNdEtaNBDSV(fdNdEtaParam);
    case kTwoNBD:
      return GetdNdEtaTwoNBD(fdNdEtaParam);
    case kGBW:
      return GetdNdEtaGBW(fdNdEtaParam);
    default:
      return 0.0;
  }
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetdNdEtaSimple(const Double_t* p) const 
{
  //Get particle density per unit of rapidity
  //using two component model
  //Parameters: npp, x
  Double_t nnp = p[0]; //=8.0
  Double_t x = p[1]; //=0.13

  return nnp*((1.-x)*fNpart/2.+x*fNcoll);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetdNdEtaGBW( const Double_t* p ) const
{
  //Get particle density per unit of rapidity
  //using the GBW model
  //Parameters: delta, lambda, snn
  Double_t delta = p[0]; //=0.79
  Double_t lambda = p[1]; //=0.288
  Double_t  snn = p[2]; //=30.25

  return fNpart*0.47*TMath::Sqrt(TMath::Power(snn,lambda))
         * TMath::Power(fNpart,(1.-delta)/3./delta);
}

//_______________________________________________________________________________
Double_t AliGlauberMC::GetdNdEtaNBDSV ( const Double_t* p ) const 
{
  //Get particle density per unit of rapidity (from Sergei)
  Double_t npp = p[0];          //=2.49
  Double_t ratioSgm2Mu = p[1];  //=1.7
  Double_t xhard = p[2];        //=0.13

  double knb=npp/(ratioSgm2Mu-1.);  // sgm^2/mu = 1+ mu/k
  double scale = (1.-xhard)*fNpart/2.+xhard*fNcoll; // effectively get number of pp collisions
  float nch1=-99.;
  if (knb*scale <1000.)
  {
    nch1=(float) NegativeBinomialRandomSV( npp*scale, knb*scale );
  }
  else
  {
    nch1=(float) NegativeBinomialRandomSV( npp*scale/2., knb*scale/2. ) +
         (float) NegativeBinomialRandomSV( npp*scale/2., knb*scale/2. );
  }
  return nch1;
}

//_______________________________________________________________________________
Double_t AliGlauberMC::GetdNdEtaNBD ( const Double_t* p ) const
{
  //Get particle density per unit of rapidity
  //using a aandomized number from a negative binomial distrubution
  //Parameters:   k  = related to distribition width=3
  //              nmean = mean of distribution=4
  //              beta = set contribution of participants / binary collisions to multiplicity=0.13
  Int_t k = TMath::Nint(p[0]);
  Double_t nmean = p[1];
  Double_t beta = p[2];

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
Double_t AliGlauberMC::GetdNdEtaTwoNBD ( const Double_t* p ) const
{
  //Get particle density per unit of rapidity
  //using random numbers from two negative binomial distributions
  //Parameters:   k1 = related to distribition width of distribution 1=3
  //              nmean1 = mean of distribution 1=4
  //              k2 = related to distribition width of distribution 2=2
  //              nmean2 = mean of distribution 2=11
  //              alpha = set contributions of distrubitin 1 / distribution 2=0.4
  //              beta = set contribution of participants / binary collisions to multiplicity =0.13
  Int_t k1 = TMath::Nint(p[0]);
  Double_t nmean1 = p[1];
  Int_t k2 = TMath::Nint(p[2]);
  Double_t nmean2 = p[3];
  Double_t alpha = p[4];
  Double_t beta = p[6];

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

//_____________________________________________________________________________


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
  if (fNpart<2) return 0.0;
  return ((fSy2-fSx2)/(fSy2+fSx2));
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEccentricityColl() const
{
  //get standard eccentricity of binary collisions
  if (fNcoll<2) return 0.0;
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
    if ( (succes=CalcEvent(bgen)) ) break; //ends if we have particparts
  }
  return succes;
}
//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon2Part() const
{
  //get participant eccentricity of participants
  if (fNpart<2) return 0.0;
  return (TMath::Sqrt(fMeanr2Cos2Phi*fMeanr2Cos2Phi+fMeanr2Sin2Phi*fMeanr2Sin2Phi)/fMeanr2);
  //return (TMath::Sqrt((fSy2-fSx2)*(fSy2-fSx2)+4*fSxy*fSxy)/(fSy2+fSx2));
}
//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon3Part() const
{
  //get participant eccentricity of participants
  if (fNpart<2) return 0.0;
  return (TMath::Sqrt(fMeanr2Cos3Phi*fMeanr2Cos3Phi+fMeanr2Sin3Phi*fMeanr2Sin3Phi)/fMeanr2);
  //return (TMath::Sqrt(fMeanr3Cos3Phi*fMeanr3Cos3Phi+fMeanr3Sin3Phi*fMeanr3Sin3Phi)/fMeanr3);
}
//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon4Part() const
{
  //get participant eccentricity of participants
  if (fNpart<2) return 0.0;
  return (TMath::Sqrt(fMeanr2Cos4Phi*fMeanr2Cos4Phi+fMeanr2Sin4Phi*fMeanr2Sin4Phi)/fMeanr2);
  //return (TMath::Sqrt(fMeanr4Cos4Phi*fMeanr4Cos4Phi+fMeanr4Sin4Phi*fMeanr4Sin4Phi)/fMeanr4);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon5Part() const
{
  //get participant eccentricity of participants
  if (fNpart<2) return 0.0;
  return (TMath::Sqrt(fMeanr2Cos5Phi*fMeanr2Cos5Phi+fMeanr2Sin5Phi*fMeanr2Sin5Phi)/fMeanr2);
  //return (TMath::Sqrt(fMeanr5Cos5Phi*fMeanr5Cos5Phi+fMeanr5Sin5Phi*fMeanr5Sin5Phi)/fMeanr5);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon2Coll() const
{
  //get epsilon2 of binary collisions
  if (fNcoll<2) return 0.0;
  return (TMath::Sqrt(fMeanr2Cos2PhiColl*fMeanr2Cos2PhiColl+fMeanr2Sin2PhiColl*fMeanr2Sin2PhiColl)/fMeanr2Coll);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon3Coll() const
{
  //get epsilon3 of binary collisions
  if (fNcoll<2) return 0.0;
  return (TMath::Sqrt(fMeanr2Cos3PhiColl*fMeanr2Cos3PhiColl+fMeanr2Sin3PhiColl*fMeanr2Sin3PhiColl)/fMeanr2Coll);
  //return (TMath::Sqrt(fMeanr3Cos3PhiColl*fMeanr3Cos3PhiColl+fMeanr3Sin3PhiColl*fMeanr3Sin3PhiColl)/fMeanr3Coll);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon4Coll() const
{
  //get epsilon4 of binary collisions
  if (fNcoll<2) return 0.0;
  return (TMath::Sqrt(fMeanr2Cos4PhiColl*fMeanr2Cos4PhiColl+fMeanr2Sin4PhiColl*fMeanr2Sin4PhiColl)/fMeanr2Coll);
  //return (TMath::Sqrt(fMeanr4Cos4PhiColl*fMeanr4Cos4PhiColl+fMeanr4Sin4PhiColl*fMeanr4Sin4PhiColl)/fMeanr4Coll);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon5Coll() const
{
  //get epsilon5 of binary collisions
  if (fNcoll<2) return 0.0;
  return (TMath::Sqrt(fMeanr2Cos5PhiColl*fMeanr2Cos5PhiColl+fMeanr2Sin5PhiColl*fMeanr2Sin5PhiColl)/fMeanr2Coll);
  //return (TMath::Sqrt(fMeanr5Cos5PhiColl*fMeanr5Cos5PhiColl+fMeanr5Sin5PhiColl*fMeanr5Sin5PhiColl)/fMeanr5Coll);
}

//______________________________________________________________________________
void AliGlauberMC::Run(Int_t nevents)
{
  //example run
  cout << "Generating " << nevents << " events..." << endl;
  TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
  TString title(Form("%s + %s (x-sect = %d mb)",fANucleus.GetName(),fBNucleus.GetName(),(Int_t) fXSect));
  if (fnt == 0)
  {
    fnt = new TNtuple(name,title,
                      "Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:"
                      "VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB:VarE:VarEColl:VarEPart:"
                      "VarEPartColl:dNdEta1:dNdEta2:dNdEta:xsect:tAA:Epsl2:Epsl3:Epsl4:Epsl5:"
                      "E2Coll:E3Coll:E4Coll:E5Coll" );
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
    //Float_t v[27];
    Float_t v[34];
    v[0]  = GetNpart();
    v[1]  = GetNcoll();
    v[2]  = fBMC;
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
    v[13] = fMeanXA;
    v[14] = fMeanYA;
    v[15] = fMeanXB;
    v[16] = fMeanYB;
    v[17] = GetEccentricity();
    v[18] = GetEccentricityColl();
    v[19] = GetEccentricityPart();
    v[20] = GetEccentricityPartColl();
    if (fDoPartProd)
    {
      v[21] = GetdNdEta();
      v[22] = GetdNdEta();
      v[23] = v[21]+v[22];
    }
    else
    {
      v[21] = 0;
      v[22] = 0;
      v[23] = 0;
    }
    v[24]=fXSect;

    Float_t mytAA=-999;
    if (GetNcoll()>0) mytAA=GetNcoll()/fXSect;
    v[25]=mytAA;
    //_____________epsilon2,3,4,4_______
    v[26] = GetEpsilon2Part();
    v[27] = GetEpsilon3Part();
    v[28] = GetEpsilon4Part();
    v[29] = GetEpsilon5Part();
    v[30] = GetEpsilon2Coll();
    v[31] = GetEpsilon3Coll();
    v[32] = GetEpsilon4Coll();
    v[33] = GetEpsilon5Coll();
    //always at the end
    fnt->Fill(v);

    if ((i%100)==0) std::cout << "Generating Event # " << i << "... \r" << flush;
  }
  std::cout << "Generating Event # " << nevents << "... \r" << endl << "Done! Succesfull events:  " << q << "  discarded events:  " << u <<"."<< endl;
}

//---------------------------------------------------------------------------------
void AliGlauberMC::RunAndSaveNtuple( Int_t n,
                                     const Option_t *sysA,
                                     const Option_t *sysB,
                                     Double_t signn,
                                     Double_t mind,
                                     Double_t r,
                                     Double_t a,
                                     const char *fname)
{
  //example run
  AliGlauberMC mcg(sysA,sysB,signn);
  mcg.SetMinDistance(mind);
  mcg.Setr(r);
  mcg.Seta(a);
  mcg.Run(n);
  TNtuple  *nt=mcg.GetNtuple();
  TFile out(fname,"recreate",fname,9);
  if(nt) nt->Write();
  printf("total cross section with a nucleon-nucleon cross section \t%f is \t%f",signn,mcg.GetTotXSect());
  out.Close();
}

//---------------------------------------------------------------------------------
void AliGlauberMC::RunAndSaveNucleons( Int_t n,
                                       const Option_t *sysA,
                                       Option_t *sysB,
                                       const Double_t signn,
                                       Double_t mind,
                                       Bool_t verbose,
                                       const char *fname)
{
  //example run
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
        AliGlauberNucleon *nucl=(AliGlauberNucleon *)nucleons->UncheckedAt(iNucl);
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
  delete fnt;
  fnt=NULL;
}
