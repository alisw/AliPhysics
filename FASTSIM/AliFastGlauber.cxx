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

/*
$Log$
Revision 1.3  2003/04/21 09:35:53  morsch
Protection against division by 0 in Binaries().

Revision 1.2  2003/04/14 14:23:44  morsch
Correction in Binaries().

Revision 1.1  2003/04/10 08:48:13  morsch
First commit.

*/

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

TF1*    AliFastGlauber::fWSb      = NULL;     
TF2*    AliFastGlauber::fWSbz     = NULL;    
TF1*    AliFastGlauber::fWSz      = NULL;     
TF1*    AliFastGlauber::fWSta     = NULL;    
TF2*    AliFastGlauber::fWStarfi  = NULL; 
TF1*    AliFastGlauber::fWStaa    = NULL;   
TF1*    AliFastGlauber::fWSgeo    = NULL;   
TF1*    AliFastGlauber::fWSbinary = NULL;   
TF1*    AliFastGlauber::fWSN      = NULL;   
Float_t AliFastGlauber::fbMax     = 0.;

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

void AliFastGlauber::DrawKernel()
{
//
//  Draw Kernel
//
    TCanvas *c6 = new TCanvas("c6","Kernel",400,10,600,700);
    c6->cd();
    fWStarfi->Draw();
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
//    if (y < 1.e-6) y = 0;
    return y;
}

Double_t AliFastGlauber::WSta(Double_t* x, Double_t* par)
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


Double_t AliFastGlauber::WStaa(Double_t* x, Double_t* par)
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

Double_t AliFastGlauber::WSgeo(Double_t* x, Double_t* par)
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

Double_t AliFastGlauber::WSN(Double_t* x, Double_t* par)
{
//
//  Number of hard processes per event
//
    Double_t b     = x[0];
    Double_t y;
    if (b != 0.) {
	y = fWSbinary->Eval(b)/fWSgeo->Eval(b);
    } else {
	y = fWSbinary->Eval(1.e-4)/fWSgeo->Eval(1.e-4);
    }
    return y;
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
	mult = 6000./fWSN->Eval(0.) * mu;
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
	} else {
	    bin = 5;
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
    
    return fWSN->Eval(b)/fWSN->Eval(0.);
}
