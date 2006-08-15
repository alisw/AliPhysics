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

/* $Id$ */

//
// Realisation of AliLhcMonitor simulating an LHC interaction region.
// The interaction region is described by the two LHC beams,
// by the beta* and the crossing angle. 
// As a monitor it records the luminosity, average luminosity and beta*
// time evolution.
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//
#include "AliLhcIRegion.h"
#include "AliLhcBeam.h"
#include "AliLHC.h"

#include <TMath.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH1.h>

ClassImp(AliLhcIRegion)

AliLhcIRegion::AliLhcIRegion(AliLHC* lhc, const char* name, const char* title)
    :TNamed(name,title),
     fAccelerator(lhc),
     fBeam1(0),
     fBeam2(0),
     fLuminosity(0.),
     fLuminosity0(0.),
     fAverageLumi(0.),
     fBetaStar(0.),
     fBetaStar0(0.),
     fCrossingAngle(0.),
     fFrequency(0.),
     fLumiArray(0),
     fAverageLumiArray(0),
     fBetaStarArray(0)
{
// Constructor
}

AliLhcIRegion::AliLhcIRegion(const AliLhcIRegion& region): 
    TNamed(region), AliLhcMonitor(region),
    fAccelerator(0),
    fBeam1(0),
    fBeam2(0),
    fLuminosity(0.),
    fLuminosity0(0.),
    fAverageLumi(0.),
    fBetaStar(0.),
    fBetaStar0(0.),
    fCrossingAngle(0.),
    fFrequency(0.),
    fLumiArray(0),
    fAverageLumiArray(0),
    fBetaStarArray(0)
{
// copy constructor
}

AliLhcIRegion::~AliLhcIRegion()
{
// Destructor

  if (fLumiArray)        delete fLumiArray;
  if (fAverageLumiArray) delete fAverageLumiArray;
  if (fBetaStarArray)    delete fBetaStarArray;
}


AliLhcIRegion& AliLhcIRegion::operator=(const  AliLhcIRegion & /*rhs*/)
{
// Assignment operator
    return *this;
}

void AliLhcIRegion::Init()
{
  // Initialization
  printf("\n Initializing Interaction Region %s", GetName());
  printf("\n ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
  // Initial Luminosity
  fBeam1 = fAccelerator->Beam(0);
  fBeam2 = fAccelerator->Beam(1);
  fFrequency = 3.e10/(2.*TMath::Pi()*fAccelerator->Radius());

  Luminosity();

  fLuminosity0 = fLuminosity;  
  fAverageLumi=0.;
  fBetaStar0   = fBetaStar;
  printf("\n IR Beta*                 :%10.3e cm        ", fBetaStar);
  printf("\n IR Initial Luminosity    :%10.3e cm^-2s^-1 ", fLuminosity);
}

Float_t AliLhcIRegion::Luminosity()
{
  Float_t sigma1 = TMath::Sqrt(fBeam1->Emittance()*fBetaStar);
  Float_t sigma2 = TMath::Sqrt(fBeam2->Emittance()*fBetaStar);
  fLuminosity = fFrequency * 
    fBeam1->N()*fBeam2->N()/(2.*TMath::Pi()*(sigma1*sigma1+sigma2*sigma2));
  return fLuminosity;
}

void AliLhcIRegion::Update()
{
  Luminosity();
}

void AliLhcIRegion::SetMonitor(Int_t n)
{
  // Initialize monitors

  if (fLumiArray)        delete fLumiArray;
  if (fAverageLumiArray) delete fAverageLumiArray;
  if (fBetaStarArray)    delete fBetaStarArray;

  fLumiArray        = new Float_t[n];
  fAverageLumiArray = new Float_t[n];
  fBetaStarArray    = new Float_t[n];

  fAverageLumiArray[0] = 0.;
  fNmax = n;
}

void AliLhcIRegion::Record()
{
  // Record some time dependent quantities
  //
  Int_t n = fAccelerator->Nt();
  // Luminosity
  
  fLumiArray[n] = fLuminosity;
  
  // Average Luminosity respecting set-up and filling time
  if (fAccelerator->Time() > fAccelerator->SetUpTime())
    fAverageLumi+=fLuminosity*fAccelerator->TimeStep();
  
  fAverageLumiArray[n] = fAverageLumi/
    (Float_t(n+1)*fAccelerator->TimeStep()+fAccelerator->FillingTime())/
    fLuminosity0;
  
  // Beta* 
  fBetaStarArray[n] = fBetaStar;
}


void AliLhcIRegion::DrawPlots()
{
  //
  // Draw the monitor plots
  //
  Float_t* t = fAccelerator->TimeA();
  //
  char name1[20], name2[20], hname[20];
  sprintf(name1,"c%s",GetName());
  sprintf(name2,"b%s",GetName());
  char title[30];
  sprintf(title,"Luminosity Lifetime for %s",GetName());
  
  //
  sprintf(hname,"%s%d",name1,0);
  TH1 *g1 = new TH1F(hname,"Luminosity",fNmax,0,t[fNmax]);
  g1->SetMinimum(0);
  g1->SetMaximum(fLumiArray[0]*1.1);
  g1->SetStats(0);
  g1->GetXaxis()->SetTitle("t (h)");
  g1->GetYaxis()->SetTitle("L(t) (cm**-2 s**-1)");
  sprintf(hname,"%s%d",name1,1);
  TH1 *g2 = new TH1F(hname,"Luminosity",fNmax,0,t[fNmax]);
  g2->SetMinimum(0);
  g2->SetMaximum(1.1);
  g2->SetStats(0);
  g2->GetXaxis()->SetTitle("t (h)");
  g2->GetYaxis()->SetTitle("L(t)/L0");
  sprintf(hname,"%s%d",name1,3);

  TH1 *g3 = new TH1F(hname,"Average Luminosity",fNmax,0,t[fNmax]);
  g3->SetMinimum(0);
  g3->SetMaximum(1.1);
  g3->SetStats(0);
  
  g3->GetXaxis()->SetTitle("t (h)");
  g3->GetYaxis()->SetTitle("L(t)/L0");
  sprintf(hname,"%s%d",name1,3);

  TH1 *g4 = new TH1F(hname,"Beta*",fNmax,0,t[fNmax]);
  g4->SetMinimum(0);
  g4->SetMaximum(fBetaStarArray[0]*1.1);
  g4->SetStats(0);
  g4->GetXaxis()->SetTitle("t (h)");
  g4->GetYaxis()->SetTitle("beta* (cm)");

  TGraph* grLumi = new TGraph(fNmax, t, fLumiArray);
  grLumi->SetHistogram(g1);

  for (Int_t i=0; i<fNmax; i++) {
    fLumiArray[i]=fLumiArray[i]/fLuminosity0;
  }
  TGraph* grLumiN = new TGraph(fNmax, t, fLumiArray);
  grLumiN->SetHistogram(g2);

  TGraph* grLumiA = new TGraph(fNmax, t, fAverageLumiArray);
  grLumiA->SetHistogram(g3);
  grLumiA->SetLineStyle(2);

  TGraph* grLumiB = new TGraph(fNmax, t, fBetaStarArray);
  grLumiB->SetHistogram(g4);

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(grLumiN);
  mg->Add(grLumiA);

 

  TCanvas *c1 = new TCanvas(name1,title, 200, 10, 700, 500);
  
  c1->SetGrid();
  mg->Draw("AC");
  mg->GetXaxis()->SetTitle("t (h)");
  mg->GetYaxis()->SetTitle("L(t)/L0 and <L>/L0");
  mg->Draw("AC");

  TCanvas *c2 = new TCanvas(name2,title, 200, 10, 700, 500);
  c2->SetGrid();
  grLumiB->Draw("AC");

}




