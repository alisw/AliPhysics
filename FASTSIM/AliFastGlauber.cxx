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

// from AliRoot
#include "AliFastGlauber.h"
// from root
#include <TH1F.h>
#include <TF1.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TFile.h>

ClassImp(AliFastGlauber)

TF1*    AliFastGlauber::fWSb            = NULL;     
TF2*    AliFastGlauber::fWSbz           = NULL;    
TF1*    AliFastGlauber::fWSz            = NULL;     
TF1*    AliFastGlauber::fWSta           = NULL;    
TF2*    AliFastGlauber::fWStarfi        = NULL; 
TF2*    AliFastGlauber::fWAlmond        = NULL; 
TF1*    AliFastGlauber::fWStaa          = NULL;   
TF1*    AliFastGlauber::fWSgeo          = NULL;   
TF1*    AliFastGlauber::fWSbinary       = NULL;   
TF1*    AliFastGlauber::fWSN            = NULL;   
TF1*    AliFastGlauber::fWPathLength0   = NULL;   
TF1*    AliFastGlauber::fWPathLength    = NULL;
TF1*    AliFastGlauber::fWEnergyDensity = NULL;   
TF1*    AliFastGlauber::fWIntRadius     = NULL;   
Float_t AliFastGlauber::fbMax           = 0.;

AliFastGlauber::AliFastGlauber()
{
//  Default Constructor
//
//  Defaults for Pb
    SetWoodSaxonParameters(6.624, 0.549, 0.00, 7.69e-4);
    SetHardCrossSection();
    SetMaxImpact();
}

void AliFastGlauber::Init(Int_t mode)
{
// Initialisation
//
//  Wood-Saxon
//
    fWSb = new TF1("WSb", WSb, 0, fbMax, 4);
    fWSb->SetParameter(0, fWSr0);
    fWSb->SetParameter(1, fWSd);
    fWSb->SetParameter(2, fWSw);
    fWSb->SetParameter(3, fWSn);

    fWSbz = new TF2("WSbz", WSbz, 0, fbMax, 4);
    fWSbz->SetParameter(0, fWSr0);
    fWSbz->SetParameter(1, fWSd);
    fWSbz->SetParameter(2, fWSw);
    fWSbz->SetParameter(3, fWSn);

    fWSz = new TF1("WSz", WSz, 0, fbMax, 5);
    fWSz->SetParameter(0, fWSr0);
    fWSz->SetParameter(1, fWSd);
    fWSz->SetParameter(2, fWSw);
    fWSz->SetParameter(3, fWSn);

//
//  Thickness
//
    fWSta = new TF1("WSta", WSta, 0., fbMax, 0);
    
//
//  Overlap Kernel
//
    fWStarfi = new TF2("WStarfi", WStarfi, 0., fbMax, 0., TMath::Pi(), 1);
    fWStarfi->SetParameter(0, 0.);     
    fWStarfi->SetNpx(200);     
    fWStarfi->SetNpy(20);     
//
//  Almond shaped interaction region
//
    fWAlmond = new TF2("WAlmond", WAlmond, -fbMax, fbMax, -fbMax, fbMax, 1);
    fWAlmond->SetParameter(0, 0.);     
    fWAlmond->SetNpx(200);     
    fWAlmond->SetNpy(200);    
//
//  Path Length as a function of Phi
//    
    fWPathLength0 = new TF1("WPathLength0", WPathLength0, -TMath::Pi(), TMath::Pi(), 2);
    fWPathLength0->SetParameter(0, 0.);
//  Pathlength definition     
    fWPathLength0->SetParameter(1, 0.);     

    fWPathLength = new TF1("WPathLength", WPathLength, -TMath::Pi(), TMath::Pi(), 3);
//  Impact Parameter
    fWPathLength->SetParameter(0, 0.);    
//  Number of interactions used for average
    fWPathLength->SetParameter(1, 1000.);    
//  Pathlength definition
    fWPathLength->SetParameter(2, 0);    

    fWIntRadius = new TF1("WIntRadius", WIntRadius, 0., fbMax, 1);
    fWIntRadius->SetParameter(0, 0.);    


//
//  Overlap
//
    if (! mode) {
	fWStaa = new TF1("WStaa", WStaa, 0., fbMax, 0);
	fWStaa->SetNpx(100);
    } else {
	TFile* f = new TFile("$(ALICE_ROOT)/FASTSIM/data/glauberPbPb.root");
	fWStaa = (TF1*) f->Get("WStaa");
    }
    
//
    fWEnergyDensity = new TF1("WEnergyDensity", WEnergyDensity, 0., 2. * fWSr0, 1);
    fWEnergyDensity->SetParameter(0, fWSr0 + 1.);
    
//
//  Geometrical Cross-Section
//
    fWSgeo = new TF1("WSgeo", WSgeo, 0., fbMax, 0);
    fWSgeo->SetNpx(100);
//
//  Hard cross section (~ binary collisions)
//
    fWSbinary = new TF1("WSbinary", WSbinary, 0., fbMax, 1);
    fWSbinary->SetParameter(0, fSigmaHard); // mb
    fWSbinary->SetNpx(100);
//
// Hard collisions per event
//
    fWSN = new TF1("WSN", WSN, 0., fbMax, 1);
    fWSN->SetNpx(100);
}

void AliFastGlauber::DrawWSb()
{
//
//  Draw Wood-Saxon Nuclear Density Function
//
    TCanvas *c1 = new TCanvas("c1","Wood Saxon",400,10,600,700);
    c1->cd();
    fWSb->Draw();
}

void AliFastGlauber::DrawOverlap()
{
//
//  Draw Overlap Function
//
    TCanvas *c2 = new TCanvas("c2","Overlap",400,10,600,700);
    c2->cd();
    fWStaa->Draw();
}

void AliFastGlauber::DrawThickness()
{
//
//  Draw Thickness Function
//
    TCanvas *c3 = new TCanvas("c3","Thickness",400,10,600,700);
    c3->cd();
    fWSta->Draw();
}

void AliFastGlauber::DrawGeo()
{
//
//  Draw Geometrical Cross-Section
//
    TCanvas *c3 = new TCanvas("c3","Geometrical Cross-Section",400,10,600,700);
    c3->cd();
    fWSgeo->Draw();
}

void AliFastGlauber::DrawBinary()
{
//
//  Draw Wood-Saxon Nuclear Density Function
//
    TCanvas *c4 = new TCanvas("c4","Binary Cross-Section",400,10,600,700);
    c4->cd();
    fWSbinary->Draw();
}

void AliFastGlauber::DrawN()
{
//
//  Draw Binaries per event
//
    TCanvas *c5 = new TCanvas("c5","Binaries per event",400,10,600,700);
    c5->cd();
    fWSN->Draw();
}

void AliFastGlauber::DrawKernel(Double_t b)
{
//
//  Draw Kernel
//
    TCanvas *c6 = new TCanvas("c6","Kernel",400,10,600,700);
    c6->cd();
    fWStarfi->SetParameter(0, b);
    fWStarfi->Draw();
}

void AliFastGlauber::DrawAlmond(Double_t b)
{
//
//  Draw Interaction Almond
//
    TCanvas *c7 = new TCanvas("c7","Almond",400,10,600,700);
    c7->cd();
    fWAlmond->SetParameter(0, b);
    fWAlmond->Draw();
}

void AliFastGlauber::DrawPathLength0(Double_t b, Int_t iopt)
{
//
//  Draw Path Length
//
    TCanvas *c8 = new TCanvas("c8","Path Length",400,10,600,700);
    c8->cd();
    fWPathLength0->SetParameter(0, b);
    fWPathLength0->SetParameter(1, Double_t(iopt));
    fWPathLength0->SetMinimum(0.); 
    fWPathLength0->SetMaximum(10.); 
    fWPathLength0->Draw();
}

void AliFastGlauber::DrawPathLength(Double_t b , Int_t ni, Int_t iopt)
{
//
//  Draw Path Length
//
    TCanvas *c9 = new TCanvas("c9","Path Length",400,10,600,700);
    c9->cd();
    fWAlmond->SetParameter(0, b);

    fWPathLength->SetParameter(0, b);
    fWPathLength->SetParameter(1, Double_t (ni));
    fWPathLength->SetParameter(2, Double_t (iopt));
    fWPathLength->SetMinimum(0.); 
    fWPathLength->SetMaximum(10.); 
    fWPathLength->Draw();
}

void AliFastGlauber::DrawIntRadius(Double_t b)
{
//
//  Draw Interaction Radius
//
    TCanvas *c10 = new TCanvas("c10","Interaction Radius",400,10,600,700);
    c10->cd();
    fWIntRadius->SetParameter(0, b);
    fWIntRadius->SetMinimum(0.);
    fWIntRadius->Draw();
}

void AliFastGlauber::DrawEnergyDensity()
{
//
//  Draw energy density
//
    TCanvas *c11 = new TCanvas("c11","Energy Density",400, 10, 600, 700);
    c11->cd();
    fWEnergyDensity->SetMinimum(0.);
    fWEnergyDensity->Draw();
}

Double_t AliFastGlauber::WSb(Double_t* x, Double_t* par)
{
//
//  Woods-Saxon Parameterisation
//  as a function of radius
//
    Double_t xx  = x[0];
    Double_t r0  = par[0];
    Double_t d   = par[1];
    Double_t w   = par[2];
    Double_t n   = par[3];
    
    Double_t y  = n * (1.+w*(xx/r0)*(xx/r0))/(1.+TMath::Exp((xx-r0)/d));

    return y;
}

Double_t AliFastGlauber::WSbz(Double_t* x, Double_t* par)
{
//
//  Wood Saxon Parameterisation
//  as a function of z and  b
//
    Double_t bb  = x[0];
    Double_t zz  = x[1];
    Double_t r0  = par[0];
    Double_t d   = par[1];
    Double_t w   = par[2];
    Double_t n   = par[3];
    Double_t xx  = TMath::Sqrt(bb*bb+zz*zz);
    Double_t y  = n * (1.+w*(xx/r0)*(xx/r0))/(1.+TMath::Exp((xx-r0)/d));

    return y;
}

Double_t AliFastGlauber::WSz(Double_t* x, Double_t* par)
{
//
//  Wood Saxon Parameterisation
//  as a function of z for fixed b
//
    Double_t bb  = par[4];
    Double_t zz  = x[0];
    Double_t r0  = par[0];
    Double_t d   = par[1];
    Double_t w   = par[2];
    Double_t n   = par[3];
    Double_t xx  = TMath::Sqrt(bb*bb+zz*zz);
    Double_t y  = n * (1.+w*(xx/r0)*(xx/r0))/(1.+TMath::Exp((xx-r0)/d));

    return y;
}

Double_t AliFastGlauber::WSta(Double_t* x, Double_t* /*par*/)
{
//
//  Thickness function 
//
    Double_t b  = x[0];
    fWSz->SetParameter(4, b);
    Double_t y  = 2. * fWSz->Integral(0., fbMax);
    return y;
}



Double_t AliFastGlauber::WStarfi(Double_t* x, Double_t* par)
{
//
//  Kernel for overlap function
//
    Double_t b    = par[0];
    Double_t r1   = x[0];
    Double_t phi  = x[1];
    Double_t r2   = TMath::Sqrt(r1 * r1 + b * b - 2. * r1 * b * TMath::Cos(phi)); 
    Double_t y    = r1 * fWSta->Eval(r1) * fWSta->Eval(r2);
    return y;
}


Double_t AliFastGlauber::WAlmond(Double_t* x, Double_t* par)
{
//
//  Almond shaped interaction region
//
    Double_t b    = par[0];
    Double_t xx   = x[0] + b/2.;
    Double_t yy   = x[1];
    Double_t r1   = TMath::Sqrt(xx * xx + yy * yy);
    Double_t phi  = TMath::ATan2(yy,xx);
    
    Double_t r2   = TMath::Sqrt(r1 * r1 + b * b - 2. * r1 * b * TMath::Cos(phi)); 
//
//  Interaction probability calculated as product of thicknesses
//
    Double_t y    = fWSta->Eval(r1) * fWSta->Eval(r2);
    return y;
}

Double_t AliFastGlauber::WIntRadius(Double_t* x, Double_t* par)
{
//
//  Average radius at which interaction takes place
//
//  Radius in the Almond
    Double_t r    = x[0];
//  Impact parameter
    Double_t b    = par[0];
    fWAlmond->SetParameter(0, b);
//  Steps in phi
    Double_t dphi = 2. * TMath::Pi() / 100.;
//  Average over phi    
    Double_t phi  = 0.;
    Double_t y    = 0.;

    for (Int_t i = 0; i < 100; i++) {
	Double_t xx = r * TMath::Cos(phi);
	Double_t yy = r * TMath::Sin(phi);
	y   += fWAlmond->Eval(xx,yy);
	phi += dphi;
    } // phi loop
// Result multiplied by Jacobian (2 pi r)     
    return (2. * TMath::Pi() * y * r / 100.);
}

Double_t AliFastGlauber::WPathLength0(Double_t* x, Double_t* par)
{
//
//  Path Length as a function of phi for interaction point fixed at (0,0)
//
//
//  Steps in r 
    const Int_t    np  = 100;
    const Double_t dr  = fbMax/Double_t(np);
//  Impact parameter    
    Double_t b      = par[0];
//  Path Length definition
    Int_t    iopt   = Int_t(par[1]);
    
//  Phi direction in Almond
    Double_t phi0   = x[0];
    Double_t r  = 0.;
    Double_t rw = 0.;
    Double_t w  = 0.;
//  Step along radial direction phi   
    for (Int_t i = 0; i < np; i++) {
//
//  Transform into target frame
//
	Double_t xx   = r * TMath::Cos(phi0) + b / 2.;
	Double_t yy   = r * TMath::Sin(phi0);
	Double_t phi  = TMath::ATan2(yy, xx);
	
	Double_t r1   = TMath::Sqrt(xx * xx + yy * yy);
// Radius in projectile frame
	Double_t r2   = TMath::Sqrt(r1 * r1 + b * b - 2. * r1 * b * TMath::Cos(phi)); 
	Double_t y    = fWSta->Eval(r1) * fWSta->Eval(r2);

	rw += y * r;
	w  += y;
	r  += dr;
    } // radial steps
//
//  My length definition (is exact for hard disk)
    if (!iopt) {
	return (2. * rw / w);
    } else {
	return TMath::Sqrt(2. * rw * dr / fWSta->Eval(0.01) / fWSta->Eval(0.01));
    }
}

Double_t AliFastGlauber::WPathLength(Double_t* x, Double_t* par)
{
//
//  Path Length as a function of phi 
//  Interaction point from random distribution
//
//
//  r-steps
// 
    const Int_t    np   = 100;
    const Double_t dr  = fbMax/Double_t(np);
//  Number of interactions
    const Int_t    npi  = Int_t (par[1]);

//
//  Impact parameter    
    Double_t b      = par[0];
//  Path Length definition 
    Int_t    iopt   = Int_t(par[2]);
//  Phi direction
    Double_t phi0   = x[0];

    printf("phi0 %f \n", phi0);
    
//  Path length 
    Double_t l = 0.;
    
    for (Int_t in = 0; in < npi; in ++) {
	Double_t rw = 0.;
	Double_t w  = 0.;
	
	// Interaction point
	Double_t x0, y0;
	fWAlmond->GetRandom2(x0, y0);
// Initial radius
	Double_t r0  = TMath::Sqrt(x0 * x0 + y0 * y0);
	Int_t    nps = Int_t ((fbMax - r0)/dr) - 1;
	
	Double_t r  = 0.;
// Radial steps
	for (Int_t i = 0; (i < nps ); i++) {
	    
// Transform into target frame
	    Double_t xx   = x0 + r * TMath::Cos(phi0) + b / 2.;
	    Double_t yy   = y0 + r * TMath::Sin(phi0);
	    Double_t phi  = TMath::ATan2(yy, xx);
	    Double_t r1   = TMath::Sqrt(xx * xx + yy * yy);
// Radius in projectile frame
	    Double_t r2   = TMath::Sqrt(r1 * r1 + b * b - 2. * r1 * b * TMath::Cos(phi)); 
	    Double_t y    = fWSta->Eval(r1) * fWSta->Eval(r2);
	    
	    rw += y * r;
	    w  += y;
	    r  += dr;
	} // steps
// Average over interactions
	if (!iopt) {
	    l += (2. * rw / w);
	} else {
	    l+= TMath::Sqrt(2. * rw * dr / fWSta->Eval(0.01) / fWSta->Eval(0.01));
	}
    } // interactions
    return (l / Double_t(npi));
}

Double_t AliFastGlauber::WStaa(Double_t* x, Double_t* /*par*/)
{
//
//  Overlap function
//
    Double_t b    = x[0];
    fWStarfi->SetParameter(0, b);
/*
    Double_t al[2];
    Double_t bl[2];
    al[0] = 0.;
    al[1] = 0.;
    bl[0] = 6.6;
    bl[1] = TMath::Pi();
    Double_t err;
    
    Double_t y =  2. * fWStarfi->IntegralMultiple(2, al, bl, 0.001, err);
    printf("WStaa: %f %f %f\n", b, y, err);
*/
//
//  MC Integration
//
    Double_t y = 0;
    for (Int_t i = 0; i < 100000; i++)
    {
	Double_t phi = TMath::Pi() * gRandom->Rndm();
	Double_t b1  = fbMax       * gRandom->Rndm();	
	y += fWStarfi->Eval(b1, phi);
    }
    y *= 2. * 0.1 *  208. * 208. * TMath::Pi() * fbMax / 100000.;
    return y;
}

Double_t AliFastGlauber::WSgeo(Double_t* x, Double_t* /*par*/)
{
//
//  Geometrical Cross-Section
//
    Double_t b    = x[0];
    Double_t taa  = fWStaa->Eval(b);
    const Double_t sigma = 55.6; // mbarn
    
    Double_t y    = 2. * TMath::Pi() * b * (1. - TMath::Exp(- sigma * taa)); // fm
    return y;
}


Double_t AliFastGlauber::WSbinary(Double_t* x, Double_t* par)
{
//
//  Number of binary collisions
//
    Double_t b     = x[0];
    Double_t sigma = par[0];
    Double_t taa   = fWStaa->Eval(b);
    
    Double_t y    = 2. * TMath::Pi() * b * sigma * taa; // fm
    return y;
}

Double_t AliFastGlauber::WSN(Double_t* x, Double_t* /*par*/)
{
//
//  Number of hard processes per event
//
    Double_t b     = x[0];
    Double_t y     = fWSbinary->Eval(b)/fWSgeo->Eval(b);
    return y;
}

Double_t AliFastGlauber::WEnergyDensity(Double_t* x, Double_t* par)
{
//
//  Initial energy density as a function of the impact parameter
//
    Double_t b     = x[0];
    Double_t rA    = par[0];
//
//  Attention: area of transverse reaction zone in hard-sphere approximation !     
    Double_t saa   = (TMath::Pi() - 2. * TMath::ASin(b/ 2./ rA)) * rA * rA 
	- b * TMath::Sqrt(rA * rA - b * b/ 4.);
    Double_t taa   = fWStaa->Eval(b);
    
    return (taa/saa);
}

void AliFastGlauber::SimulateTrigger(Int_t n)
{
    //
    //  Simulates Trigger
    //
    TH1F* mbtH = new TH1F("mbtH", "MB Trigger b-Distribution",   100, 0., 20.);
    TH1F* hdtH = new TH1F("hdtH", "Hard Trigger b-Distribution", 100, 0., 20.);   
    TH1F* mbmH = new TH1F("mbmH", "MB Trigger Multiplicity Distribution",   100, 0., 8000.);
    TH1F* hdmH = new TH1F("hdmH", "Hard Trigger Multiplicity Distribution", 100, 0., 8000.);   

    mbtH->SetXTitle("b [fm]");
    hdtH->SetXTitle("b [fm]");    
    mbmH->SetXTitle("Multiplicity");
    hdmH->SetXTitle("Multiplicity");    

    TCanvas *c0 = new TCanvas("c0","Trigger Simulation",400,10,600,700);    
    c0->Divide(2,1);
    TCanvas *c1 = new TCanvas("c1","Trigger Simulation",400,10,600,700);    
    c1->Divide(1,2);

    //
    //
    Init(1);
    for (Int_t iev = 0; iev < n; iev++)
    {
	Float_t b, p, mult;
	GetRandom(b, p, mult);
	mbtH->Fill(b,1.);
	hdtH->Fill(b, p);
	mbmH->Fill(mult, 1.);
	hdmH->Fill(mult, p);

	c0->cd(1);
	mbtH->Draw();
	c0->cd(2);
	hdtH->Draw();	
	c0->Update();

	c1->cd(1);
	mbmH->Draw();
	c1->cd(2);
	hdmH->Draw();	
	c1->Update();
    }
}

void AliFastGlauber::GetRandom(Float_t& b, Float_t& p, Float_t& mult)
{
    //
    // Gives back a random impact parameter, hard trigger probability and multiplicity
    //
	b = fWSgeo->GetRandom();
	Float_t mu = fWSN->Eval(b);
	p = 1.-TMath::Exp(-mu);
	mult = 6000./fWSN->Eval(1.) * mu;
}

void AliFastGlauber::GetRandom(Int_t& bin, Bool_t& hard)
{
    //
    // Gives back a random impact parameter bin, and hard trigger decission
    //
	Float_t b  = fWSgeo->GetRandom();
	Float_t mu = fWSN->Eval(b) * fSigmaHard;
	Float_t p  = 1.-TMath::Exp(-mu);
	if (b < 5.) {
	    bin = 1;
	} else if (b <  8.6) {
	    bin = 2;
	} else if (b < 11.2) {
	    bin = 3;
	} else if (b < 13.2) {
	    bin = 4;
	} else if (b < 15.0) {
	    bin = 5;
	} else {
	    bin = 6;
	}
	
	hard = kFALSE;
	
	Float_t r = gRandom->Rndm();
	
	if (r < p) hard = kTRUE;
}


Float_t  AliFastGlauber::GetRandomImpactParameter(Float_t bmin, Float_t bmax)
{
    //
    // Gives back a random impact parameter in the range bmin .. bmax
    //

    Float_t b = -1.;
    while(b < bmin || b > bmax)
	b = fWSgeo->GetRandom();
    return b;
}

Double_t AliFastGlauber::CrossSection(Double_t b1, Double_t b2)
{
    //
    // Return cross-section integrated from b1 to b2 
    //
    
    return fWSgeo->Integral(b1, b2)/100.;
}

Double_t AliFastGlauber::FractionOfHardCrossSection(Double_t b1, Double_t b2)
{
    //
    // Return raction of hard cross-section integrated from b1 to b2 
    //
    
    return fWSbinary->Integral(b1, b2)/fWSbinary->Integral(0., 100.);
}


Double_t AliFastGlauber::Binaries(Double_t b)
{
    //
    // Return number of binary collisions normalized to 1 at b=0
    //
    
    return fWSN->Eval(b)/fWSN->Eval(0.001);
}
