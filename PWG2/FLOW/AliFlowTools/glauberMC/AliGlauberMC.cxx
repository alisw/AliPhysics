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
  fMeanX2Parts(0),
  fMeanY2Parts(0),
  fMeanXYParts(0),
  fMeanXParts(0),
  fMeanYParts(0),
  fMeanOXParts(0),
  fMeanOYParts(0),
  fMeanXColl(0),
  fMeanYColl(0),
  fMeanOXColl(0),
  fMeanOYColl(0),
  fMeanX2Coll(0),
  fMeanY2Coll(0), 
  fMeanXYColl(0),
  fMeanXCom(0),
  fMeanYCom(0),
  fMeanOXCom(0),
  fMeanOYCom(0),
  fMeanX2Com(0),
  fMeanY2Com(0), 
  fMeanXYCom(0),
  fMeanXSystem(0),
  fMeanYSystem(0),
  fMeanXA(0),
  fMeanYA(0),
  fMeanXB(0),
  fMeanYB(0),
  fMeanOXA(0),
  fMeanOYA(0),
  fMeanOXB(0),
  fMeanOYB(0),
  fBMC(0),
  fEvents(0),
  fTotalEvents(0),
  fBMin(0.),
  fBMax(20.),
  fMultType(kNBDSV),
  fMaxNpartFound(0),
  fONpart(0),
  fONcoll(0),
  fONcom(0),
  fNpart(0),
  fNcoll(0),
  fNcom(0),
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
  fMeanr2Com(0),
  fMeanr3Com(0),
  fMeanr4Com(0),
  fMeanr5Com(0),
  fMeanr2Cos2PhiCom(0),
  fMeanr2Sin2PhiCom(0),
  fMeanr2Cos3PhiCom(0),
  fMeanr2Sin3PhiCom(0),
  fMeanr2Cos4PhiCom(0),
  fMeanr2Sin4PhiCom(0),
  fMeanr2Cos5PhiCom(0),
  fMeanr2Sin5PhiCom(0),
  fMeanr3Cos3PhiCom(0),
  fMeanr3Sin3PhiCom(0),
  fMeanr4Cos4PhiCom(0),
  fMeanr4Sin4PhiCom(0),
  fMeanr5Cos5PhiCom(0),
  fMeanr5Sin5PhiCom(0),
  fSx2Parts(0.),
  fSy2Parts(0.),
  fSxyParts(0.),
  fSx2Coll(0.),
  fSy2Coll(0.),
  fSxyColl(0.),
  fSx2Com(0.),
  fSy2Com(0.),
  fSxyCom(0.),
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
  fMeanX2Parts(in.fMeanX2Parts),
  fMeanY2Parts(in.fMeanY2Parts),
  fMeanXYParts(in.fMeanXYParts),
  fMeanXParts(in.fMeanXParts),
  fMeanYParts(in.fMeanYParts),
  fMeanOXParts(in.fMeanOXParts),
  fMeanOYParts(in.fMeanOYParts),
  fMeanXColl(in.fMeanXColl),
  fMeanYColl(in.fMeanYColl),
  fMeanOXColl(in.fMeanOXColl),
  fMeanOYColl(in.fMeanOYColl),
  fMeanX2Coll(in.fMeanX2Coll),
  fMeanY2Coll(in.fMeanY2Coll),
  fMeanXYColl(in.fMeanXYColl),
  fMeanXCom(in.fMeanXCom),
  fMeanYCom(in.fMeanYCom),
  fMeanOXCom(in.fMeanOXCom),
  fMeanOYCom(in.fMeanOYCom),
  fMeanX2Com(in.fMeanX2Com),
  fMeanY2Com(in.fMeanY2Com),
  fMeanXYCom(in.fMeanXYCom),
  fMeanXSystem(in.fMeanXSystem),
  fMeanYSystem(in.fMeanYSystem),
  fMeanXA(in.fMeanXA),
  fMeanYA(in.fMeanYA),
  fMeanXB(in.fMeanXB),
  fMeanYB(in.fMeanYB),
  fMeanOXA(in.fMeanOXA),
  fMeanOYA(in.fMeanOYA),
  fMeanOXB(in.fMeanOXB),
  fMeanOYB(in.fMeanOYB),
  fBMC(in.fBMC),
  fEvents(in.fEvents),
  fTotalEvents(in.fTotalEvents),
  fBMin(in.fBMin),
  fBMax(in.fBMax),
  fMultType(in.fMultType),
  fMaxNpartFound(in.fMaxNpartFound),
  fONpart(in.fONpart),
  fONcoll(in.fONcoll),
  fONcom(in.fONcom),
  fNpart(in.fNpart),
  fNcoll(in.fNcoll),
  fNcom(in.fNcom),
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
  fMeanr2Com(in.fMeanr2Com),
  fMeanr3Com(in.fMeanr3Com),
  fMeanr4Com(in.fMeanr4Com),
  fMeanr5Com(in.fMeanr5Com),
  fMeanr2Cos2PhiCom(in.fMeanr2Cos2PhiCom),
  fMeanr2Sin2PhiCom(in.fMeanr2Sin2PhiCom),
  fMeanr2Cos3PhiCom(in.fMeanr2Cos3PhiCom),
  fMeanr2Sin3PhiCom(in.fMeanr2Sin3PhiCom),
  fMeanr2Cos4PhiCom(in.fMeanr2Cos4PhiCom),
  fMeanr2Sin4PhiCom(in.fMeanr2Sin4PhiCom),
  fMeanr2Cos5PhiCom(in.fMeanr2Cos5PhiCom),
  fMeanr2Sin5PhiCom(in.fMeanr2Sin5PhiCom),
  fMeanr3Cos3PhiCom(in.fMeanr3Cos3PhiCom),
  fMeanr3Sin3PhiCom(in.fMeanr3Sin3PhiCom),
  fMeanr4Cos4PhiCom(in.fMeanr4Cos4PhiCom),
  fMeanr4Sin4PhiCom(in.fMeanr4Sin4PhiCom),
  fMeanr5Cos5PhiCom(in.fMeanr5Cos5PhiCom),
  fMeanr5Sin5PhiCom(in.fMeanr5Sin5PhiCom),
  fSx2Parts(in.fSx2Parts),
  fSy2Parts(in.fSy2Parts),
  fSxyParts(in.fSxyParts),
  fSx2Coll(in.fSx2Coll),
  fSy2Coll(in.fSy2Coll),
  fSxyColl(in.fSxyColl),
  fSx2Com(in.fSx2Com),
  fSy2Com(in.fSy2Com),
  fSxyCom(in.fSxyCom),
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
  fMeanOXParts=in.fMeanOXParts;
  fMeanOYParts=in.fMeanOYParts;
  fMeanXColl=in.fMeanXColl;
  fMeanYColl=in.fMeanYColl;
  fMeanOXColl=in.fMeanOXColl;
  fMeanOYColl=in.fMeanOYColl;
  fMeanX2Coll=in.fMeanX2Coll;
  fMeanY2Coll=in.fMeanY2Coll;
  fMeanXYColl=in.fMeanXYColl;
  fMeanr2Coll=in.fMeanr2Coll;
  fMeanr3Coll=in.fMeanr3Coll;
  fMeanr4Coll=in.fMeanr4Coll;
  fMeanr5Coll=in.fMeanr5Coll;
  fMeanXCom=in.fMeanXColl;
  fMeanYCom=in.fMeanYColl;
  fMeanOXCom=in.fMeanOXCom;
  fMeanOYCom=in.fMeanOYCom;
  fMeanX2Com=in.fMeanX2Com;
  fMeanY2Com=in.fMeanY2Com;
  fMeanXYCom=in.fMeanXYCom;
  fMeanr2Com=in.fMeanr2Com;
  fMeanr3Com=in.fMeanr3Com;
  fMeanr4Com=in.fMeanr4Com;
  fMeanr5Com=in.fMeanr5Com;
  fMeanXSystem=in.fMeanXSystem;
  fMeanYSystem=in.fMeanYSystem; 
  fMeanXA=in.fMeanXA;
  fMeanYA=in.fMeanYA;
  fMeanXB=in.fMeanXB;
  fMeanYB=in.fMeanYB;
  fMeanOXA=in.fMeanOXA;
  fMeanOYA=in.fMeanOYA;
  fMeanOXB=in.fMeanOXB;
  fMeanOYB=in.fMeanOYB;
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
  fNcom=in.fNcom;
  fONpart=in.fONpart;
  fONcoll=in.fONcoll;
  fONcom=in.fONcom;
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
  fMeanr2Cos2PhiCom=in.fMeanr2Cos2PhiCom;
  fMeanr2Sin2PhiCom=in.fMeanr2Sin2PhiCom;
  fMeanr2Cos3PhiCom=in.fMeanr2Cos3PhiCom;
  fMeanr2Sin3PhiCom=in.fMeanr2Sin3PhiCom;
  fMeanr2Cos4PhiCom=in.fMeanr2Cos4PhiCom;
  fMeanr2Sin4PhiCom=in.fMeanr2Sin4PhiCom;
  fMeanr2Cos5PhiCom=in.fMeanr2Cos5PhiCom;
  fMeanr2Sin5PhiCom=in.fMeanr2Sin5PhiCom;
  fMeanr3Cos3PhiCom=in.fMeanr3Cos3PhiCom;
  fMeanr3Sin3PhiCom=in.fMeanr3Sin3PhiCom;
  fMeanr4Cos4PhiCom=in.fMeanr4Cos4PhiCom;
  fMeanr4Sin4PhiCom=in.fMeanr4Sin4PhiCom;
  fMeanr5Cos5PhiCom=in.fMeanr5Cos5PhiCom;
  fMeanr5Sin5PhiCom=in.fMeanr5Sin5PhiCom;
  fSx2Parts=in.fSx2Parts;
  fSy2Parts=in.fSy2Parts;
  fSxyParts=in.fSxyParts;
  fSx2Coll=in.fSx2Coll;
  fSy2Coll=in.fSy2Coll;
  fSxyColl=in.fSxyColl;
  fSx2Com=in.fSx2Com;
  fSy2Com=in.fSy2Com;
  fSxyCom=in.fSxyCom;
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
  fNcom=0;
  fONpart=0;
  fONcoll=0;
  fONcom=0;
  fMeanX2=0.;
  fMeanY2=0.;
  fMeanXY=0.;
  fMeanXParts=0.;
  fMeanYParts=0.;
  fMeanOXParts=0.;
  fMeanOYParts=0.;
  fMeanXColl=0.;
  fMeanYColl=0.;
  fMeanOXColl=0.;
  fMeanOYColl=0.;
  fMeanX2Coll=0.;
  fMeanY2Coll=0.;
  fMeanXYColl=0.;
  fMeanXCom=0.;
  fMeanYCom=0.;
  fMeanOXCom=0.;
  fMeanOYCom=0.;
  fMeanX2Com=0.;
  fMeanY2Com=0.;
  fMeanXYCom=0.;
  fMeanXSystem=0.;
  fMeanYSystem=0.;
  fMeanXA=0.;
  fMeanYA=0.;
  fMeanXB=0.;
  fMeanYB=0.;
  fMeanOXA=0.;
  fMeanOYA=0.;
  fMeanOXB=0.;
  fMeanOYB=0.;
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
  fMeanr2Com=0.;
  fMeanr3Com=0.;
  fMeanr4Com=0.;
  fMeanr5Com=0.;
  fMeanr2Cos2PhiCom=0.;
  fMeanr2Sin2PhiCom=0.;
  fMeanr2Cos3PhiCom=0.;
  fMeanr2Sin3PhiCom=0.;
  fMeanr2Cos4PhiCom=0.;
  fMeanr2Sin4PhiCom=0.;
  fMeanr2Cos5PhiCom=0.;
  fMeanr2Sin5PhiCom=0.;
  fMeanr3Cos3PhiCom=0.;
  fMeanr3Sin3PhiCom=0.;
  fMeanr4Cos4PhiCom=0.;
  fMeanr4Sin4PhiCom=0.;
  fMeanr5Cos5PhiCom=0.;
  fMeanr5Sin5PhiCom=0.;

  for (Int_t i = 0; i<fAN; i++)
  {
    AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->UncheckedAt(i));
    Double_t oXA = nucleonA->GetX();
    Double_t oYA = nucleonA->GetY();
    //fMeanOXSystem  += oXA;
    //fMeanOYSystem  += oYA;
    fMeanOXA  += oXA;
    fMeanOYA  += oYA;

    if(nucleonA->IsWounded())
    {
      fONpart++;
      fMeanOXParts  += oXA;
      fMeanOYParts  += oYA;
    }
  }

  for (Int_t i = 0; i<fBN; i++)
  {
    AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->UncheckedAt(i));
    Double_t oXB=nucleonB->GetX();
    Double_t oYB=nucleonB->GetY();
    
    if(nucleonB->IsWounded())
    {
      Int_t oNcoll = nucleonB->GetNColl();
      fONpart++;
      fMeanOXParts  += oXB;
      fMeanOXColl  += oXB*oNcoll;
      fMeanOXCom  += oXB*((1-0.15)+0.15*oNcoll);
      fMeanOYParts  += oYB;
      fMeanOYColl += oYB*oNcoll;
      fMeanOYColl += oYB*((1-0.15)+0.15*oNcoll);
      fONcoll += oNcoll;
      fONcom += ((1-0.15)+0.15*oNcoll);
    }
  }

  if (fONpart>0)
  {
    fMeanOXParts /= fONpart;
    fMeanOYParts /= fONpart;
  }
  else
  {
    fMeanOXParts = 0;
    fMeanOYParts = 0;
  }

  if (fONcoll>0)
  {
    fMeanOXColl /= fONcoll;
    fMeanOYColl /= fONcoll;
  }
  else
  {
    fMeanOXColl = 0;
    fMeanOYColl = 0;
  }

 if (fONcom>0)
  {
    fMeanOXCom /= fONcom;
    fMeanOYCom /= fONcom;
  }
  else
  {
    fMeanOXCom = 0;
    fMeanOYCom = 0;
  }
  
  //////////////////////////////////////////////////////////////////
  for (Int_t i = 0; i<fAN; i++)
  {
    AliGlauberNucleon *nucleonA=(AliGlauberNucleon*)(fNucleonsA->UncheckedAt(i));
    Double_t xAA = nucleonA->GetX(); // X
    Double_t yAA = nucleonA->GetY(); // Y
    Double_t xAPart = xAA - fMeanOXParts; // X'
    Double_t yAPart = yAA - fMeanOYParts; // Y'
    Double_t r2APart = xAPart *xAPart+yAPart*yAPart;     // r'^2
    Double_t rAPart = TMath::Sqrt(r2APart);  // r'
    Double_t r3APart = r2APart*rAPart;
    Double_t r4APart = r3APart*rAPart;
    Double_t r5APart = r4APart*rAPart;
    Double_t phiAPart = TMath::ATan2(yAPart,xAPart);
    Double_t sin2PhiAPart = TMath::Sin(2.*phiAPart);
    Double_t cos2PhiAPart = TMath::Cos(2.*phiAPart);
    Double_t sin3PhiAPart = TMath::Sin(3.*phiAPart);
    Double_t cos3PhiAPart = TMath::Cos(3.*phiAPart);
    Double_t sin4PhiAPart = TMath::Sin(4.*phiAPart);
    Double_t cos4PhiAPart = TMath::Cos(4.*phiAPart);
    Double_t sin5PhiAPart = TMath::Sin(5.*phiAPart);
    Double_t cos5PhiAPart = TMath::Cos(5.*phiAPart);
   
    fMeanXSystem  += xAA;
    fMeanYSystem  += yAA;
    fMeanXA  += xAA;
    fMeanYA  += yAA;
    fMeanX2 += xAA * xAA;
    fMeanY2 += yAA * yAA;
    fMeanXY += xAA * yAA;
    
    if(nucleonA->IsWounded())
     {
      fNpart++;
      fMeanXParts  += xAPart;
      fMeanYParts  += yAPart;
      fMeanX2Parts += xAPart * xAPart;
      fMeanY2Parts += yAPart * yAPart;
      fMeanXYParts += xAPart * yAPart;
      fMeanr2 += r2APart;
      fMeanr3 += r3APart;
      fMeanr4 += r4APart;
      fMeanr5 += r5APart;
      fMeanr2Cos2Phi += r2APart*cos2PhiAPart;
      fMeanr2Sin2Phi += r2APart*sin2PhiAPart;
      fMeanr2Cos3Phi += r2APart*cos3PhiAPart;
      fMeanr2Sin3Phi += r2APart*sin3PhiAPart;
      fMeanr2Cos4Phi += r2APart*cos4PhiAPart;
      fMeanr2Sin4Phi += r2APart*sin4PhiAPart;
      fMeanr2Cos5Phi += r2APart*cos5PhiAPart;
      fMeanr2Sin5Phi += r2APart*sin5PhiAPart;
      fMeanr3Cos3Phi += r3APart*cos3PhiAPart;
      fMeanr3Sin3Phi += r3APart*sin3PhiAPart;
      fMeanr4Cos4Phi += r4APart*cos4PhiAPart;
      fMeanr4Sin4Phi += r4APart*sin4PhiAPart;
      fMeanr5Cos5Phi += r5APart*cos5PhiAPart;
      fMeanr5Sin5Phi += r5APart*sin5PhiAPart;
    }
  }
  
  for (Int_t i = 0; i<fBN; i++)
    {
      AliGlauberNucleon *nucleonB=(AliGlauberNucleon*)(fNucleonsB->UncheckedAt(i));
      Double_t xBB = nucleonB->GetX();
      Double_t yBB = nucleonB->GetY();
      // for Wounded
      Double_t xBPart = xBB - fMeanOXParts; // X'
      Double_t yBPart = yBB - fMeanOYParts; // Y'
      Double_t r2BPart = xBPart*xBPart+yBPart*yBPart;     // r'^2
      Double_t rBPart = TMath::Sqrt(r2BPart);  // r'
      Double_t r3BPart = r2BPart*rBPart;
      Double_t r4BPart = r3BPart*rBPart;
      Double_t r5BPart = r4BPart*rBPart;
      Double_t phiBPart = TMath::ATan2(yBPart,xBPart);
      Double_t sin2PhiBPart = TMath::Sin(2.*phiBPart);
      Double_t cos2PhiBPart = TMath::Cos(2.*phiBPart);
      Double_t sin3PhiBPart = TMath::Sin(3.*phiBPart);
      Double_t cos3PhiBPart = TMath::Cos(3.*phiBPart);
      Double_t sin4PhiBPart = TMath::Sin(4.*phiBPart);
      Double_t cos4PhiBPart = TMath::Cos(4.*phiBPart);
      Double_t sin5PhiBPart = TMath::Sin(5.*phiBPart);
      Double_t cos5PhiBPart = TMath::Cos(5.*phiBPart);
      // for Binary
      Double_t xBColl = xBB - fMeanOXColl; // X'-Binary
      Double_t yBColl = yBB - fMeanOYColl; // Y'-B
      Double_t r2BColl = xBColl*xBColl+yBColl*yBColl;     // r'^2-B
      Double_t rBColl = TMath::Sqrt(r2BColl);  // r'-B
      Double_t r3BColl = r2BColl*rBColl;
      Double_t r4BColl = r3BColl*rBColl;
      Double_t r5BColl = r4BColl*rBColl;
      Double_t phiBColl = TMath::ATan2(yBColl,xBColl);
      Double_t sin2PhiBColl = TMath::Sin(2.*phiBColl);
      Double_t cos2PhiBColl = TMath::Cos(2.*phiBColl);
      Double_t sin3PhiBColl = TMath::Sin(3.*phiBColl);
      Double_t cos3PhiBColl = TMath::Cos(3.*phiBColl);
      Double_t sin4PhiBColl = TMath::Sin(4.*phiBColl);
      Double_t cos4PhiBColl = TMath::Cos(4.*phiBColl);
      Double_t sin5PhiBColl = TMath::Sin(5.*phiBColl);
      Double_t cos5PhiBColl = TMath::Cos(5.*phiBColl);
      // for combine
      Double_t xBCom = xBB - fMeanOXCom; // X'-Combine
      Double_t yBCom = yBB - fMeanOYCom; // Y'-C
      Double_t r2BCom = xBCom *xBCom+yBCom*yBCom;     // r'^2-C
      Double_t rBCom = TMath::Sqrt(r2BCom);  // r'-C
      Double_t r3BCom = r2BCom*rBCom;
      Double_t r4BCom = r3BCom*rBCom;
      Double_t r5BCom = r4BCom*rBCom;
      Double_t phiBCom = TMath::ATan2(yBCom,xBCom);
      Double_t sin2PhiBCom = TMath::Sin(2.*phiBCom);
      Double_t cos2PhiBCom = TMath::Cos(2.*phiBCom);
      Double_t sin3PhiBCom = TMath::Sin(3.*phiBCom);
      Double_t cos3PhiBCom = TMath::Cos(3.*phiBCom);
      Double_t sin4PhiBCom = TMath::Sin(4.*phiBCom);
      Double_t cos4PhiBCom = TMath::Cos(4.*phiBCom);
      Double_t sin5PhiBCom = TMath::Sin(5.*phiBCom);
      Double_t cos5PhiBCom = TMath::Cos(5.*phiBCom);   
      
      fMeanXSystem  += xBB;
      fMeanYSystem  += yBB;
      fMeanXB  += xBB;
      fMeanYB  += yBB;
      fMeanX2 += xBB*xBB;
      fMeanY2 += yBB*yBB;
      fMeanXY += xBB*yBB;
      
      if(nucleonB->IsWounded())
	{
	  Int_t ncoll = nucleonB->GetNColl();
	  fNpart++;
	  fMeanXParts  += xBPart;
	  fMeanXColl  += xBColl*ncoll;
	  fMeanXCom  += xBCom*((1-0.15)+0.15*ncoll);
	  fMeanYParts  += yBPart;
	  fMeanYColl += yBColl*ncoll;
	  fMeanYCom += yBCom*((1-0.15)+0.15*ncoll);
	  fMeanX2Parts += xBPart * xBPart;
	  fMeanX2Coll += xBColl*xBColl*ncoll;
	  fMeanX2Com += xBCom*xBCom*((1-0.15)+0.15*ncoll);
	  fMeanY2Parts += yBPart * yBPart;
	  fMeanY2Coll += yBColl*yBColl*ncoll;
	  fMeanY2Com += yBCom*yBCom*((1-0.15)+0.15*ncoll);
	  fMeanXYParts += xBPart * yBPart;
	  fMeanXYColl += xBColl*yBColl*ncoll;
	  fMeanXYCom += xBCom*yBCom*((1-0.15)+0.15*ncoll);
	  fNcoll += ncoll;
	  fNcom += ((1-0.15)+0.15*ncoll);
	  fMeanr2 += r2BPart;
	  fMeanr3 += r3BPart;
	  fMeanr4 += r4BPart;
	  fMeanr5 += r5BPart;
	  fMeanr2Cos2Phi += r2BPart*cos2PhiBPart;
	  fMeanr2Sin2Phi += r2BPart*sin2PhiBPart;
	  fMeanr2Cos3Phi += r2BPart*cos3PhiBPart;
	  fMeanr2Sin3Phi += r2BPart*sin3PhiBPart;
	  fMeanr2Cos4Phi += r2BPart*cos4PhiBPart;
	  fMeanr2Sin4Phi += r2BPart*sin4PhiBPart;
	  fMeanr2Cos5Phi += r2BPart*cos5PhiBPart;
	  fMeanr2Sin5Phi += r2BPart*sin5PhiBPart;
	  fMeanr3Cos3Phi += r3BPart*cos3PhiBPart;
	  fMeanr3Sin3Phi += r3BPart*sin3PhiBPart;
	  fMeanr4Cos4Phi += r4BPart*cos4PhiBPart;
	  fMeanr4Sin4Phi += r4BPart*sin4PhiBPart;
	  fMeanr5Cos5Phi += r5BPart*cos5PhiBPart;
	  fMeanr5Sin5Phi += r5BPart*sin5PhiBPart;
	  fMeanr2Coll += r2BColl*ncoll;
	  fMeanr3Coll += r3BColl*ncoll;
	  fMeanr4Coll += r4BColl*ncoll;
	  fMeanr5Coll += r5BColl*ncoll;
	  fMeanr2Cos2PhiColl += r2BColl*cos2PhiBColl*ncoll;
	  fMeanr2Sin2PhiColl += r2BColl*sin2PhiBColl*ncoll;
	  fMeanr2Cos3PhiColl += r2BColl*cos3PhiBColl*ncoll;
	  fMeanr2Sin3PhiColl += r2BColl*sin3PhiBColl*ncoll;
	  fMeanr2Cos4PhiColl += r2BColl*cos4PhiBColl*ncoll;
	  fMeanr2Sin4PhiColl += r2BColl*sin4PhiBColl*ncoll;
	  fMeanr2Cos5PhiColl += r2BColl*cos5PhiBColl*ncoll;
	  fMeanr2Sin5PhiColl += r2BColl*sin5PhiBColl*ncoll;
	  fMeanr3Cos3PhiColl += r3BColl*cos3PhiBColl*ncoll;
	  fMeanr3Sin3PhiColl += r3BColl*sin3PhiBColl*ncoll;
	  fMeanr4Cos4PhiColl += r4BColl*cos4PhiBColl*ncoll;
	  fMeanr4Sin4PhiColl += r4BColl*sin4PhiBColl*ncoll;
	  fMeanr5Cos5PhiColl += r5BColl*cos5PhiBColl*ncoll;
	  fMeanr5Sin5PhiColl += r5BColl*sin5PhiBColl*ncoll;
	  fMeanr2Com += r2BCom*((1-0.15)+0.15*ncoll);
	  fMeanr3Com += r3BCom*((1-0.15)+0.15*ncoll);
	  fMeanr4Com += r4BCom*((1-0.15)+0.15*ncoll);
	  fMeanr5Com += r5BCom*((1-0.15)+0.15*ncoll);
	  fMeanr2Cos2PhiCom += r2BCom*cos2PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr2Sin2PhiCom += r2BCom*sin2PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr2Cos3PhiCom += r2BCom*cos3PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr2Sin3PhiCom += r2BCom*sin3PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr2Cos4PhiCom += r2BCom*cos4PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr2Sin4PhiCom += r2BCom*sin4PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr2Cos5PhiCom += r2BCom*cos5PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr2Sin5PhiCom += r2BCom*sin5PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr3Cos3PhiCom += r3BCom*cos3PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr3Sin3PhiCom += r3BCom*sin3PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr4Cos4PhiCom += r4BCom*cos4PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr4Sin4PhiCom += r4BCom*sin4PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr5Cos5PhiCom += r5BCom*cos5PhiBCom*((1-0.15)+0.15*ncoll);
	  fMeanr5Sin5PhiCom += r5BCom*sin5PhiBCom*((1-0.15)+0.15*ncoll);
	}
    }
  
  if (fNpart>0)
    {
      fMeanXParts /= fNpart;
      fMeanYParts /= fNpart;
      fMeanX2Parts /= fNpart;
      fMeanY2Parts /= fNpart;
      fMeanXYParts /= fNpart;
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
      fMeanX2Parts = 0;
      fMeanY2Parts = 0;
      fMeanXYParts = 0;
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
  
 if (fNcom>0)
    {
      fMeanXCom /= fNcom;
      fMeanYCom /= fNcom;
      fMeanX2Com /= fNcom;
      fMeanY2Com /= fNcom;
      fMeanXYCom /= fNcom;
      fMeanr2Com /= fNcom;
      fMeanr3Com /= fNcom;
      fMeanr4Com /= fNcom;
      fMeanr5Com /= fNcom;
      fMeanr2Cos2PhiCom /= fNcom;
      fMeanr2Sin2PhiCom /= fNcom;
      fMeanr2Cos3PhiCom /= fNcom;
      fMeanr2Sin3PhiCom /= fNcom;
      fMeanr2Cos4PhiCom /= fNcom;
      fMeanr2Sin4PhiCom /= fNcom;
      fMeanr2Cos5PhiCom /= fNcom;
      fMeanr2Sin5PhiCom /= fNcom;
      fMeanr3Cos3PhiCom /= fNcom;
      fMeanr3Sin3PhiCom /= fNcom;
      fMeanr4Cos4PhiCom /= fNcom;
      fMeanr4Sin4PhiCom /= fNcom;
      fMeanr5Cos5PhiCom /= fNcom;
      fMeanr5Sin5PhiCom /= fNcom;
  }

 else
  {
    fMeanXCom = 0;
    fMeanYCom = 0;
    fMeanX2Com = 0;
    fMeanY2Com = 0;
    fMeanXYCom = 0;
    fMeanr2Cos2PhiCom =0;
    fMeanr2Sin2PhiCom =0;
    fMeanr2Cos3PhiCom =0;
    fMeanr2Sin3PhiCom =0;
    fMeanr2Cos4PhiCom =0;
    fMeanr2Sin4PhiCom =0;
    fMeanr2Cos5PhiCom =0;
    fMeanr2Sin5PhiCom =0;
    fMeanr3Cos3PhiCom =0;
    fMeanr3Sin3PhiCom =0;
    fMeanr4Cos4PhiCom =0;
    fMeanr4Sin4PhiCom =0;
    fMeanr5Cos5PhiCom =0;
    fMeanr5Sin5PhiCom =0;
  }

  if(fAN+fBN>0)
  {
    fMeanXSystem /= (fAN + fBN);
    fMeanYSystem /= (fAN + fBN);
    fMeanX2 /= (fAN + fBN);
    fMeanY2 /= (fAN + fBN);
    fMeanXY /= (fAN + fBN);
  }
  else
  {
    fMeanXSystem = 0;
    fMeanYSystem = 0;
    fMeanX2 = 0;
    fMeanY2 = 0;
    fMeanXY = 0; 
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
  


    
  //////////////////////////////////////////////////////////////////
  fSx2Parts=fMeanX2Parts-(fMeanXParts*fMeanXParts);
  fSy2Parts=fMeanY2Parts-(fMeanYParts*fMeanYParts);
  fSxyParts=fMeanXYParts-fMeanXParts*fMeanYParts;
  fSx2Coll=fMeanX2Coll-(fMeanXColl*fMeanXColl);
  fSy2Coll=fMeanY2Coll-(fMeanYColl*fMeanYColl);
  fSxyColl=fMeanXYColl-fMeanXColl*fMeanYColl;
  fSx2Com=fMeanX2Com-(fMeanXCom*fMeanXCom);
  fSy2Com=fMeanY2Com-(fMeanYCom*fMeanYCom);
  fSxyCom=fMeanXYCom-fMeanXCom*fMeanYCom;
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
  return (TMath::Sqrt((fSy2Parts-fSx2Parts)*(fSy2Parts-fSx2Parts)+4*fSxyParts*fSxyParts)/(fSy2Parts+fSx2Parts));
}

//_____________________________________________________________________________
Double_t AliGlauberMC::GetEccentricityPartColl() const
{
  //get participant eccentricity of binary collisions
  if (fNcoll<2) return 0.0;
  return (TMath::Sqrt((fSy2Coll-fSx2Coll)*(fSy2Coll-fSx2Coll)+4*fSxyColl*fSxyColl)/(fSy2Coll+fSx2Coll));
}

//_____________________________________________________________________________
Double_t AliGlauberMC::GetEccentricityPartCom() const
{
  //get participant eccentricity of binary collisions
  if (fNcom<2) return 0.0;
  return (TMath::Sqrt((fSy2Com-fSx2Com)*(fSy2Com-fSx2Com)+4*fSxyCom*fSxyCom)/(fSy2Com+fSx2Com));
}
//______________________________________________________________________________
Double_t AliGlauberMC::GetEccentricity() const
{
  //get standard eccentricity of participants
  if (fNpart<2) return 0.0;
  return ((fSy2Parts-fSx2Parts)/(fSy2Parts+fSx2Parts));
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEccentricityColl() const
{
  //get standard eccentricity of binary collisions
  if (fNcoll<2) return 0.0;
  return ((fSy2Coll-fSx2Coll)/(fSy2Coll+fSx2Coll));
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEccentricityCom() const
{
  //get standard eccentricity of combined weight 
  if (fNcom<2) return 0.0;
  //return ((fSy2Com-fSx2Com)/(fSy2Com+fSx2Com));
  return (((fMeanY2Com-(fMeanYCom*fMeanYCom)) - (fMeanX2Com-(fMeanXCom*fMeanXCom))) / ((fMeanY2Com-(fMeanYCom*fMeanYCom)) + (fMeanX2Com-(fMeanXCom*fMeanXCom))) );
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetStoa() const
{
  //get standard Transverse Overlap Area
  if (fNpart<2) return 0.0;
  return ( TMath::Pi()*(TMath::Sqrt(fSx2Parts))*(TMath::Sqrt(fSy2Parts)));
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
}
//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon3Part() const
{
  //get participant eccentricity of participants
  if (fNpart<2) return 0.0;
  //return (TMath::Sqrt(fMeanr2Cos3Phi*fMeanr2Cos3Phi+fMeanr2Sin3Phi*fMeanr2Sin3Phi)/fMeanr2);
  return (TMath::Sqrt(fMeanr3Cos3Phi*fMeanr3Cos3Phi+fMeanr3Sin3Phi*fMeanr3Sin3Phi)/fMeanr3);
}
//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon4Part() const
{
  //get participant eccentricity of participants
  if (fNpart<2) return 0.0;
  //return (TMath::Sqrt(fMeanr2Cos4Phi*fMeanr2Cos4Phi+fMeanr2Sin4Phi*fMeanr2Sin4Phi)/fMeanr2);
  return (TMath::Sqrt(fMeanr4Cos4Phi*fMeanr4Cos4Phi+fMeanr4Sin4Phi*fMeanr4Sin4Phi)/fMeanr4);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon5Part() const
{
  //get participant eccentricity of participants
  if (fNpart<2) return 0.0;
  //return (TMath::Sqrt(fMeanr2Cos5Phi*fMeanr2Cos5Phi+fMeanr2Sin5Phi*fMeanr2Sin5Phi)/fMeanr2);
  return (TMath::Sqrt(fMeanr5Cos5Phi*fMeanr5Cos5Phi+fMeanr5Sin5Phi*fMeanr5Sin5Phi)/fMeanr5);
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
  //return (TMath::Sqrt(fMeanr2Cos3PhiColl*fMeanr2Cos3PhiColl+fMeanr2Sin3PhiColl*fMeanr2Sin3PhiColl)/fMeanr2Coll);
  return (TMath::Sqrt(fMeanr3Cos3PhiColl*fMeanr3Cos3PhiColl+fMeanr3Sin3PhiColl*fMeanr3Sin3PhiColl)/fMeanr3Coll);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon4Coll() const
{
  //get epsilon4 of binary collisions
  if (fNcoll<2) return 0.0;
  //return (TMath::Sqrt(fMeanr2Cos4PhiColl*fMeanr2Cos4PhiColl+fMeanr2Sin4PhiColl*fMeanr2Sin4PhiColl)/fMeanr2Coll);
  return (TMath::Sqrt(fMeanr4Cos4PhiColl*fMeanr4Cos4PhiColl+fMeanr4Sin4PhiColl*fMeanr4Sin4PhiColl)/fMeanr4Coll);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon5Coll() const
{
  //get epsilon5 of binary collisions
  if (fNcoll<2) return 0.0;
  //return (TMath::Sqrt(fMeanr2Cos5PhiColl*fMeanr2Cos5PhiColl+fMeanr2Sin5PhiColl*fMeanr2Sin5PhiColl)/fMeanr2Coll);
  return (TMath::Sqrt(fMeanr5Cos5PhiColl*fMeanr5Cos5PhiColl+fMeanr5Sin5PhiColl*fMeanr5Sin5PhiColl)/fMeanr5Coll);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon2Com() const
{
  //get epsilon2 of binary collisions
  if (fNcom<2) return 0.0;
  return (TMath::Sqrt(fMeanr2Cos2PhiCom*fMeanr2Cos2PhiCom+fMeanr2Sin2PhiCom*fMeanr2Sin2PhiCom)/fMeanr2Com);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon3Com() const
{
  //get epsilon3 of binary collisions
  if (fNcom<2) return 0.0;
  //return (TMath::Sqrt(fMeanr2Cos3PhiCom*fMeanr2Cos3PhiCom+fMeanr2Sin3PhiCom*fMeanr2Sin3PhiCom)/fMeanr2Com);
  return (TMath::Sqrt(fMeanr3Cos3PhiCom*fMeanr3Cos3PhiCom+fMeanr3Sin3PhiCom*fMeanr3Sin3PhiCom)/fMeanr3Com);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon4Com() const
{
  //get epsilon4 of binary collisions
  if (fNcom<2) return 0.0;
  //return (TMath::Sqrt(fMeanr2Cos4PhiCom*fMeanr2Cos4PhiCom+fMeanr2Sin4PhiCom*fMeanr2Sin4PhiCom)/fMeanr2Com);
  return (TMath::Sqrt(fMeanr4Cos4PhiCom*fMeanr4Cos4PhiCom+fMeanr4Sin4PhiCom*fMeanr4Sin4PhiCom)/fMeanr4Com);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetEpsilon5Com() const
{
  //get epsilon5 of binary collisions
  if (fNcom<2) return 0.0;
  //return (TMath::Sqrt(fMeanr2Cos5PhiCom*fMeanr2Cos5PhiCom+fMeanr2Sin5PhiCom*fMeanr2Sin5PhiCom)/fMeanr2Com);
  return (TMath::Sqrt(fMeanr5Cos5PhiCom*fMeanr5Cos5PhiCom+fMeanr5Sin5PhiCom*fMeanr5Sin5PhiCom)/fMeanr5Com);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetPsi2() const
{
  return ((TMath::ATan2(fMeanr2Sin2Phi,fMeanr2Cos2Phi)+TMath::Pi())/2);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetPsi3() const
{
  //return ((TMath::ATan2(fMeanr2Sin3Phi,fMeanr2Cos3Phi)+TMath::Pi())/3);
  return ((TMath::ATan2(fMeanr3Sin3Phi,fMeanr3Cos3Phi)+TMath::Pi())/3);  
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetPsi4() const
{
  //return ((TMath::ATan2(fMeanr2Sin4Phi,fMeanr2Cos4Phi)+TMath::Pi())/4);
  return ((TMath::ATan2(fMeanr4Sin4Phi,fMeanr4Cos4Phi)+TMath::Pi())/4);
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetPsi5() const
{
  //return ((TMath::ATan2(fMeanr2Sin5Phi,fMeanr2Cos5Phi)+TMath::Pi())/5);
  return ((TMath::ATan2(fMeanr5Sin5Phi,fMeanr5Cos5Phi)+TMath::Pi())/5);  
}

/*_____________________________________________________________________________
Double_t AliGlauberMC::GetE43Part() const
{
  //get participant eccentricity of participants
  if (fNpart<2) return 0.0;
  //return ((TMath::Sqrt(fMeanr2Cos2Phi*fMeanr2Cos2Phi+fMeanr2Sin2Phi*fMeanr2Sin2Phi)/fMeanr2) * (TMath::Sqrt(fMeanr2Cos2Phi*fMeanr2Cos2Phi+fMeanr2Sin2Phi*fMeanr2Sin2Phi)/fMeanr2) * (TMath::Sqrt(fMeanr2Cos4Phi*fMeanr2Cos4Phi+fMeanr2Sin4Phi*fMeanr2Sin4Phi)/fMeanr2) * TMath::Cos (((TMath::ATan2(fMeanr2Sin2Phi,fMeanr2Cos2Phi)+TMath::Pi())/2) - ((TMath::ATan2(fMeanr2Sin4Phi,fMeanr2Cos4Phi)+TMath::Pi())/4)));
  return ((TMath::Sqrt(fMeanr2Cos2Phi*fMeanr2Cos2Phi+fMeanr2Sin2Phi*fMeanr2Sin2Phi)/fMeanr2) * (TMath::Sqrt(fMeanr2Cos2Phi*fMeanr2Cos2Phi+fMeanr2Sin2Phi*fMeanr2Sin2Phi)/fMeanr2) * (TMath::Sqrt(fMeanr4Cos4Phi*fMeanr4Cos4Phi+fMeanr4Sin4Phi*fMeanr4Sin4Phi)/fMeanr4) * TMath::Cos(((TMath::ATan2(fMeanr2Sin2Phi,fMeanr2Cos2Phi)+TMath::Pi())/2) - ((TMath::ATan2(fMeanr4Sin4Phi,fMeanr4Cos4Phi)+TMath::Pi())/4)));
}


//______________________________________________________________________________
Double_t AliGlauberMC::GetE43Coll() const
{
  //get epsilon5 of binary collisions
  if (fNcoll<2) return 0.0;
  //return ((TMath::Sqrt(fMeanr2Cos2PhiColl*fMeanr2Cos2PhiColl+fMeanr2Sin2PhiColl*fMeanr2Sin2PhiColl)/fMeanr2Coll) * (TMath::Sqrt(fMeanr2Cos2PhiColl*fMeanr2Cos2PhiColl+fMeanr2Sin2PhiColl*fMeanr2Sin2PhiColl)/fMeanr2Coll) * (TMath::Sqrt(fMeanr2Cos4PhiColl*fMeanr2Cos4PhiColl+fMeanr2Sin4PhiColl*fMeanr2Sin4PhiColl)/fMeanr2Coll) * TMath::Cos(((TMath::ATan2(fMeanr2Sin2Phi,fMeanr2Cos2Phi)+TMath::Pi())/2) - ((TMath::ATan2(fMeanr2Sin4Phi,fMeanr2Cos4Phi)+TMath::Pi())/4)));
  return ((TMath::Sqrt(fMeanr2Cos2PhiColl*fMeanr2Cos2PhiColl+fMeanr2Sin2PhiColl*fMeanr2Sin2PhiColl)/fMeanr2Coll) * (TMath::Sqrt(fMeanr2Cos2PhiColl*fMeanr2Cos2PhiColl+fMeanr2Sin2PhiColl*fMeanr2Sin2PhiColl)/fMeanr2Coll) * (TMath::Sqrt(fMeanr4Cos4PhiColl*fMeanr4Cos4PhiColl+fMeanr4Sin4PhiColl*fMeanr4Sin4PhiColl)/fMeanr4Coll) * TMath::Cos(((TMath::ATan2(fMeanr2Sin2Phi,fMeanr2Cos2Phi)+TMath::Pi())/2) - ((TMath::ATan2(fMeanr4Sin4Phi,fMeanr4Cos4Phi)+TMath::Pi())/4)));
}

//______________________________________________________________________________
Double_t AliGlauberMC::GetE43Com() const
{
  //get epsilon5 of binary collisions
  if (fNcom<2) return 0.0;
  //return ((TMath::Sqrt(fMeanr2Cos2PhiCom*fMeanr2Cos2PhiCom+fMeanr2Sin2PhiCom*fMeanr2Sin2PhiCom)/fMeanr2Com) * (TMath::Sqrt(fMeanr2Cos2PhiCom*fMeanr2Cos2PhiCom+fMeanr2Sin2PhiCom*fMeanr2Sin2PhiCom)/fMeanr2Com) * (TMath::Sqrt(fMeanr2Cos4PhiCom*fMeanr2Cos4PhiCom+fMeanr2Sin4PhiCom*fMeanr2Sin4PhiCom)/fMeanr2Com) * TMath::Cos(((TMath::ATan2(fMeanr2Sin2Phi,fMeanr2Cos2Phi)+TMath::Pi())/2) - ((TMath::ATan2(fMeanr2Sin4Phi,fMeanr2Cos4Phi)+TMath::Pi())/4)));
  return ((TMath::Sqrt(fMeanr2Cos2PhiCom*fMeanr2Cos2PhiCom+fMeanr2Sin2PhiCom*fMeanr2Sin2PhiCom)/fMeanr2Com) * (TMath::Sqrt(fMeanr2Cos2PhiCom*fMeanr2Cos2PhiCom+fMeanr2Sin2PhiCom*fMeanr2Sin2PhiCom)/fMeanr2Com) * (TMath::Sqrt(fMeanr4Cos4PhiCom*fMeanr4Cos4PhiCom+fMeanr4Sin4PhiCom*fMeanr4Sin4PhiCom)/fMeanr4Com) * TMath::Cos(((TMath::ATan2(fMeanr2Sin2Phi,fMeanr2Cos2Phi)+TMath::Pi())/2) - ((TMath::ATan2(fMeanr4Sin4Phi,fMeanr4Cos4Phi)+TMath::Pi())/4)));
}
//___________________________________________________________________________

Double_t AliGlauberMC::GetPsi4m2() const
{
  //return (TMath::Cos(4*((TMath::ATan2(fMeanr2Sin4Phi,fMeanr2Cos4Phi)+TMath::Pi())/4) - ((TMath::ATan2(fMeanr2Sin2Phi,fMeanr2Cos2Phi)+TMath::Pi())/2)));
  return (TMath::Cos(4*(((TMath::ATan2(fMeanr4Sin4Phi,fMeanr4Cos4Phi)+TMath::Pi())/4)-((TMath::ATan2(fMeanr2Sin2Phi,fMeanr2Cos2Phi)+TMath::Pi())/2))));
}
*/
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
                      "Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB:VarE:Stoa:VarEColl:VarECom:VarEPart:VarEPartColl:VarEPartCom:dNdEta:dNdEtaGBW:dNdEtaTwoNBD:xsect:tAA:Epsl2:Epsl3:Epsl4:Epsl5:E2Coll:E3Coll:E4Coll:E5Coll:E2Com:E3Com:E4Com:E5Com:Psi2:Psi3:Psi4:Psi5");
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
    Float_t v[45];
    v[0]  = GetNpart();
    v[1]  = GetNcoll();
    v[2]  = fBMC;
    v[3]  = fMeanXParts;
    v[4]  = fMeanYParts;
    v[5]  = fMeanX2Parts;
    v[6]  = fMeanY2Parts;
    v[7]  = fMeanXYParts;
    v[8]  = fSx2Parts;
    v[9]  = fSy2Parts;
    v[10] = fSxyParts;
    v[11] = fMeanXSystem;
    v[12] = fMeanYSystem;
    v[13] = fMeanXA;
    v[14] = fMeanYA;
    v[15] = fMeanXB;
    v[16] = fMeanYB;
    v[17] = GetEccentricity();
    v[18] = GetStoa();
    v[19] = GetEccentricityColl();
    v[20] = GetEccentricityCom();
    v[21] = GetEccentricityPart();
    v[22] = GetEccentricityPartColl();
    v[23] = GetEccentricityPartCom();
    if (fDoPartProd)
    {
      v[24] = GetdNdEta();
      v[25] = GetdNdEta();
      v[26] = v[24]+v[25];
    }
    else
    {
      v[24] = 0;
      v[25] = 0;
      v[26] = 0;
    }
    v[27]=fXSect;

    Float_t mytAA=-999;
    if (GetNcoll()>0) mytAA=GetNcoll()/fXSect;
    v[28]=mytAA;
    //_____________epsilon2,3,4,4_______
    v[29] = GetEpsilon2Part();
    v[30] = GetEpsilon3Part();
    v[31] = GetEpsilon4Part();
    v[32] = GetEpsilon5Part();
    v[33] = GetEpsilon2Coll();
    v[34] = GetEpsilon3Coll();
    v[35] = GetEpsilon4Coll();
    v[36] = GetEpsilon5Coll();
    v[37] = GetEpsilon2Com();
    v[38] = GetEpsilon3Com();
    v[39] = GetEpsilon4Com();
    v[40] = GetEpsilon5Com();
    v[41] = GetPsi2();
    v[42] = GetPsi3();
    v[43] = GetPsi4();
    v[44] = GetPsi5();
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
