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
// Class that holds all parameters about an LHC beam.
// The parameters can change with time.
// A monitor can be set that stores the time distribution of 
// emittance and number of particles per bunch.
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TMultiGraph.h>

#include "AliLhcBeam.h"
#include "AliLHC.h"

ClassImp(AliLhcBeam)

AliLhcBeam::AliLhcBeam(AliLHC* lhc):
    fAccelerator(lhc),
    fN(0),
    fN0(0),
    fNEmittance(0.),
    fEmittance(0.),
    fEmittance0(0.),
    fEmittanceL(0.),
    fEmittanceL0(0.),
    fEnergySpread(0.),
    fA(0),
    fZ(0),
    fEnergy(0.),
    fGamma(0.),
    fTimeArray(0),
    fEmittanceArray(0),
    fEmittanceLArray(0)
{
// Constructor
}

AliLhcBeam::AliLhcBeam(const AliLhcBeam& beam): 
    TNamed(beam), AliLhcMonitor(beam),
    fAccelerator(0),
    fN(0),
    fN0(0),
    fNEmittance(0.),
    fEmittance(0.),
    fEmittance0(0.),
    fEmittanceL(0.),
    fEmittanceL0(0.),
    fEnergySpread(0.),
    fA(0),
    fZ(0),
    fEnergy(0.),
    fGamma(0.),
    fTimeArray(0),
    fEmittanceArray(0),
    fEmittanceLArray(0)
{
// copy constructor
}

AliLhcBeam::~AliLhcBeam()
{
// Destructor

}

void AliLhcBeam::Init()
{
  // Initialization
  printf("\n Initializing Beam");
  printf("\n ^^^^^^^^^^^^^^^^^");
  // Scale energy with regidity 
  fEnergy     *= fZ/fA;
  fGamma       = fEnergy/0.938272;
  fEmittance   = fNEmittance/fGamma;
  fEmittance0  = fEmittance;
  fEmittanceL *= fZ;
  fEmittanceL0 = fEmittanceL; 
  fN0=fN;

  printf("\n Beam Energy                :%10.3e GeV", fEnergy);
  printf("\n Beam Normalized Emittance  :%10.3e cm ", fNEmittance);
  printf("\n Beam Particles per Bunch   :%10.3e    ", fN);
}

void AliLhcBeam::RemoveParticles(Float_t loss)
{
  fN-=loss;
}

void AliLhcBeam::IncreaseEmittance(Float_t de, Float_t del)
{
//
// Increase the emittance
  fEmittance    *= (1.+de);
  fEmittanceL   *= (1.+del);
  fEnergySpread *= (1.+del);
}

AliLhcBeam& AliLhcBeam::operator=(const  AliLhcBeam & /*rhs*/)
{
// Assignment operator
    return *this;
}
void AliLhcBeam::SetMonitor(Int_t n)
{
//
// Initialize a monitor with n time bins
  fNmax = n;
  if (fEmittanceArray)  delete fEmittanceArray;
  if (fEmittanceLArray) delete fEmittanceLArray;


  fEmittanceArray  = new Float_t[n];
  fEmittanceLArray = new Float_t[n];
}

void AliLhcBeam::Record()
{
    fEmittanceArray [fAccelerator->Nt()] = fEmittance/fEmittance0;
    fEmittanceLArray[fAccelerator->Nt()] = fEmittanceL/fEmittanceL0;
}


void AliLhcBeam::DrawPlots()
{
  // Draw monitor plots
  Float_t* t =  fAccelerator->TimeA();
  
  TH1 *e1 = new TH1F("e1","Hor. Emittance",fNmax,0,t[fNmax]);
  e1->SetMinimum(1);
  e1->SetMaximum(fEmittanceArray[fNmax]*1.1);
  e1->SetStats(0);
  e1->GetXaxis()->SetTitle("t (h)");
  e1->GetYaxis()->SetTitle("rel. Emittance (t)");

  TH1 *e2 = new TH1F("e2","Long. Emittance",fNmax,0,t[fNmax]);
  e2->SetMinimum(1);
  e2->SetMaximum(fEmittanceLArray[fNmax]*1.1);
  e2->SetStats(0);
  e2->GetXaxis()->SetTitle("t (h)");
  e2->GetYaxis()->SetTitle("rel. Emittance (t)");


  TGraph* grE  = new TGraph(fNmax, t, fEmittanceArray);
  grE->SetHistogram(e1);
  TGraph* grEl = new TGraph(fNmax, t, fEmittanceLArray);
  grEl->SetHistogram(e2);
  grEl->SetLineStyle(2);
 
  TMultiGraph* mg = new TMultiGraph();
  mg->Add(grE);
  mg->Add(grEl);

  TCanvas *c2 = new TCanvas("c2","Emittance Increase", 200, 10, 700, 500);
  c2->SetGrid();
  mg->Draw("AC");  
  mg->GetXaxis()->SetTitle("t (h)");
  mg->GetYaxis()->SetTitle("rel. Emittance(t)");
  mg->Draw("AC");
}




