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
// Realisation of AliLhcProcess for the fast simulation of the
// Intra Beam Scattering process
// in transverse and longitudinal direction.
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TMultiGraph.h>

#include "AliLhcProcessIBS.h"
#include "AliLHC.h"
#include "AliLhcIRegion.h"
#include "AliLhcBeam.h"

ClassImp(AliLhcProcessIBS)

  Double_t func(Double_t *x, Double_t *par);

AliLhcProcessIBS::AliLhcProcessIBS(AliLHC* lhc, const char* name, const char* title)
    :AliLhcProcess(lhc,name,title),
     fCrossSection(0.),
     fIRegions(0),
     fTaux(0.),
     fTaue(0.),
     fTauxArray(0),
     fTaueArray(0)
{
// Constructor
}

AliLhcProcessIBS::AliLhcProcessIBS(const AliLhcProcessIBS& ibs):
    AliLhcProcess(ibs),
    fCrossSection(0.),
    fIRegions(0),
    fTaux(0.),
    fTaue(0.),
    fTauxArray(0),
    fTaueArray(0)
{
// Copy Constructor
}


AliLhcProcessIBS::~AliLhcProcessIBS()
{
// Destructor

}

void AliLhcProcessIBS::Init()
{
  // Initialization
   const Float_t r0=1.535e-16;

   printf("\n Initializing Process %s", GetName());
   printf("\n ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");

   fIRegions = fAccelerator->IRegions();

   for (Int_t i = 0; i < 2; i++) {
     fBeam[i] = fAccelerator->Beam(i);
     fR[i]    = r0*fBeam[i]->Z()*fBeam[i]->Z()/fBeam[i]->A();
     fE[i]    = 0.938*fBeam[i]->A();
   }
}

void AliLhcProcessIBS::Evolve(Float_t dt)
{
//
// Evolve by one time step dt
  printf("\n Here process %s %f:", GetName(), dt);
   for (Int_t i=0; i<2; i++) {
     // Density
     Float_t sige    = fBeam[i]->EnergySpread();
     Float_t avbeta  = fAccelerator->AverageBeta();
     Float_t avd     = fAccelerator->AverageDisp();

     Float_t gamma   = fBeam[i]->Gamma();
     Float_t epsx    = fBeam[i]->Emittance()*gamma;
     Float_t epsy    = epsx;
     Float_t epse    = fBeam[i]->LongEmittance();
     Float_t sigxb   = TMath::Sqrt(epsx/gamma*avbeta);
     Float_t sigx    = TMath::Sqrt(sigxb*sigxb+(avd*avd*sige*sige));
     Float_t ssigx   = TMath::Sqrt(epsx/gamma/avbeta);
     Float_t sigy    = sigx;
     Float_t ssigy   = ssigx;

     Float_t asd = fBeam[i]->N()*fR[i]*fR[i]*fE[i]/
       (16.*TMath::Pi()*gamma*epsx*epsy*epse);
     
     // impact parameter
     Float_t d    = sige*avd/sigx;
     // 
     Float_t at = sige*TMath::Sqrt(1.-d*d)/(gamma*ssigx);
     Float_t bt = sige*TMath::Sqrt(1.-d*d)/(gamma*ssigy);
     Float_t ct = sige*TMath::Sqrt(1.-d*d)*TMath::Sqrt(4.*sigy/fR[i]);
     //
     Double_t par[3];
     par[0] = at;
     par[1] = bt;
     par[2] = ct;
     TF1 *fct = new TF1("func",func, 0., 1., 3);
     
     Float_t f = (Float_t) 8.0*TMath::Pi()*
	 fct->Integral(0., 1., par, 1.e-5);
     //
     fTaux =  1./(asd*f*(d*d-at*at/2.));
     fTaue =  1./(asd*f*(1.-d*d));
     //     printf("\n taux, taue %f %f", taux, taue);
     //     Float_t tauy = -2./at*at/asd/f;
     fBeam[i]->IncreaseEmittance(dt/fTaux, dt/fTaue);
   }

}

void AliLhcProcessIBS::SetMonitor(Int_t n)
{
  // Initialize Monitor
  if (fTauxArray) delete fTauxArray;
  if (fTaueArray) delete fTaueArray;
  fTauxArray = new Float_t[n];
  fTaueArray = new Float_t[n];
  fNmax = n;
}

void AliLhcProcessIBS::Record()
{
  // Record monitor quantities
    fTauxArray[fAccelerator->Nt()] = fTaux/3600.;
    fTaueArray[fAccelerator->Nt()] = fTaue/3600.;
}


void AliLhcProcessIBS::DrawPlots()
{
  // Draw monitor plots
  Float_t* t =  fAccelerator->TimeA();
  
  TH1 *t1 = new TH1F("t1","Hor. IBS growth time",fNmax,0,t[fNmax]);
  t1->SetMinimum(0);
  t1->SetMaximum(fTauxArray[fNmax]*1.1);
  t1->SetStats(0);
  t1->GetXaxis()->SetTitle("t (h)");
  t1->GetYaxis()->SetTitle("tau_x (t)");

  TH1 *t2 = new TH1F("t2","Long. IBS growth time",fNmax,0,t[fNmax]);
  t2->SetMinimum(0);
  t2->SetMaximum(fTaueArray[fNmax]*1.1);
  t2->SetStats(0);
  t2->GetXaxis()->SetTitle("t (h)");
  t2->GetYaxis()->SetTitle("tau_l (t)");

  TGraph* grTaux = new TGraph(fNmax, fAccelerator->TimeA(), fTauxArray);
  grTaux->SetHistogram(t1);

  TGraph* grTaue = new TGraph(fNmax, fAccelerator->TimeA(), fTaueArray);
  grTaue->SetHistogram(t2);
  grTaue->SetLineStyle(2);

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(grTaux);
  mg->Add(grTaue);

  TCanvas *c3 = new TCanvas("c3","IBS", 200, 10, 700, 500);
  c3->SetGrid();
  mg->Draw("AC");  
  mg->GetXaxis()->SetTitle("t (h)");
  mg->GetYaxis()->SetTitle("IBS Growth Time (h)");
  mg->Draw("AC");
}





AliLhcProcessIBS& AliLhcProcessIBS::operator=(const  AliLhcProcessIBS & /*rhs*/)
{
// Assignment operator
    return *this;
}

Double_t func(Double_t *x, Double_t *par)
{
  Double_t a  = par[0];
  Double_t b  = par[1];
  Double_t cc = par[2];

  const Double_t kbc = 0.5772;
  Double_t xx = x[0];
  Double_t x2=xx*xx;
  Double_t x1=1.0-x2;
  Double_t p=1.0/TMath::Sqrt(x2+a*a*x1);
  Double_t q=1.0/TMath::Sqrt(x2+b*b*x1);
  return (1.0-3.0*x2)*p*q*(2.0*TMath::Log(0.5*cc*(p+q))-kbc);
}

