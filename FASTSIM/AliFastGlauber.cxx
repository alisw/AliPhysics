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
// Utility class to make simple Glauber type calculations 
//           for SYMMETRIC collision geometries (AA):
// Impact parameter, production points, reaction plane dependence
//
// The SimulateTrigger method can be used for simple MB and hard-process
// (binary scaling) trigger studies.
// 
// Some basic quantities can be visualized directly.
//
// The default set-up for PbPb or AUAu collisions can be read from a file 
// calling Init(1) or Init(2) if you want to read Almonds too.
//
// ***** If you change settings dont forget to call init afterwards, *****
// ***** in order to update the formulas with the new parameters.    *****
// 
// Author: andreas.morsch@cern.ch
//=================== Added by A. Dainese 11/02/04 ===========================
// Calculate path length for a parton with production point (x0,y0)
// and propagation direction (ux=cos(phi0),uy=sin(phi0)) 
// in a collision with impact parameter b and functions that make use
// of it.
//=================== Added by A. Dainese 05/03/04 ===========================
// Calculation of line integrals I0 and I1
//  integral0 =    \int_0^ellCut dl*(T_A*T_B)(x0+l*ux,y0+l*uy)
//  integral1 =    \int_0^ellCut dl*l*(T_A*T_B)(x0+l*ux,y0+l*uy)
// mostly for use in the Quenching class
//=================== Added by C. Loizdes 27/03/04 ===========================
// Handling of AuAu collisions
// More get/set functions
// Comments, units and clearing of code
//

// from AliRoot
#include "AliFastGlauber.h"
// from root
#include <Riostream.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TStyle.h>

ClassImp(AliFastGlauber)

Float_t AliFastGlauber::fgBMax           = 0.;
TF1*    AliFastGlauber::fgWSb            = NULL;     
TF2*    AliFastGlauber::fgWSbz           = NULL;    
TF1*    AliFastGlauber::fgWSz            = NULL;     
TF1*    AliFastGlauber::fgWSta           = NULL;    
TF2*    AliFastGlauber::fgWStarfi        = NULL; 
TF2*    AliFastGlauber::fgWAlmond        = NULL; 
TF1*    AliFastGlauber::fgWStaa          = NULL;   
TF1*    AliFastGlauber::fgWSgeo          = NULL;   
TF1*    AliFastGlauber::fgWSbinary       = NULL;   
TF1*    AliFastGlauber::fgWSN            = NULL;   
TF1*    AliFastGlauber::fgWPathLength0   = NULL;   
TF1*    AliFastGlauber::fgWPathLength    = NULL;
TF1*    AliFastGlauber::fgWEnergyDensity = NULL;   
TF1*    AliFastGlauber::fgWIntRadius     = NULL;   
TF2*    AliFastGlauber::fgWKParticipants = NULL; 
TF1*    AliFastGlauber::fgWParticipants  = NULL; 
TF2*    AliFastGlauber::fgWAlmondCurrent = NULL;    
TF2*    AliFastGlauber::fgWAlmondFixedB[40]; 
const Int_t AliFastGlauber::fgkMCInts = 100000;
Int_t AliFastGlauber::fgCounter = 0;       

AliFastGlauber::AliFastGlauber(): 
    fWSr0(0.),
    fWSd(0.), 
    fWSw(0.), 
    fWSn(0.), 
    fSigmaHard(0.),
    fSigmaNN(0.),  
    fA(0),         
    fBmin(0.),     
    fBmax(0.),     
    fEllDef(0),    
    fName()     
{
  //  Default Constructor 
  fgCounter++;
  if(fgCounter>1)
    Error("AliFastGlauber","More than one instance (%d) is not supported, check your code!",fgCounter);

  //  Defaults for Pb
  SetMaxImpact();
  SetLengthDefinition();
  SetPbPbLHC();
}

AliFastGlauber::AliFastGlauber(const AliFastGlauber & gl)
    :TObject(gl),
     fWSr0(0.),
     fWSd(0.), 
     fWSw(0.), 
     fWSn(0.), 
     fSigmaHard(0.),
     fSigmaNN(0.),  
     fA(0),         
     fBmin(0.),     
     fBmax(0.),     
     fEllDef(0),    
     fName()     
{
// Copy constructor
    gl.Copy(*this);
}

AliFastGlauber::~AliFastGlauber()
{
// Destructor
  fgCounter--;
  for(Int_t k=0; k<40; k++) delete fgWAlmondFixedB[k];
}

void AliFastGlauber::SetAuAuRhic()
{
  //Set all parameters for RHIC
  SetWoodSaxonParametersAu();
  SetHardCrossSection();
  SetNNCrossSection(42);
  SetNucleus(197);
  SetFileName("$(ALICE_ROOT)/FASTSIM/data/glauberAuAu.root");
}

void AliFastGlauber::SetPbPbLHC()
{
  //Set all parameters for LHC
  SetWoodSaxonParametersPb();
  SetHardCrossSection();
  SetNNCrossSection();
  SetNucleus();
  SetFileName();
}

void AliFastGlauber::Init(Int_t mode)
{
  // Initialisation
  //    mode = 0; all functions are calculated 
  //    mode = 1; overlap function is read from file (for Pb-Pb only)
  //    mode = 2; interaction almond functions are read from file 
  //              USE THIS FOR PATH LENGTH CALC.!
  //

  // 
  Reset(); 
  //

  //
  //  Wood-Saxon
  //
  fgWSb = new TF1("WSb", WSb, 0, fgBMax, 4);
  fgWSb->SetParameter(0, fWSr0);
  fgWSb->SetParameter(1, fWSd);
  fgWSb->SetParameter(2, fWSw);
  fgWSb->SetParameter(3, fWSn);

  fgWSbz = new TF2("WSbz", WSbz, 0, fgBMax, 0, fgBMax, 4);
  fgWSbz->SetParameter(0, fWSr0);
  fgWSbz->SetParameter(1, fWSd);
  fgWSbz->SetParameter(2, fWSw);
  fgWSbz->SetParameter(3, fWSn);

  fgWSz = new TF1("WSz", WSz, 0, fgBMax, 5);
  fgWSz->SetParameter(0, fWSr0);
  fgWSz->SetParameter(1, fWSd);
  fgWSz->SetParameter(2, fWSw);
  fgWSz->SetParameter(3, fWSn);

  //
  //  Thickness
  //
  fgWSta = new TF1("WSta", WSta, 0., fgBMax, 0);
    
  //
  //  Overlap Kernel
  //
  fgWStarfi = new TF2("WStarfi", WStarfi, 0., fgBMax, 0., TMath::Pi(), 1);
  fgWStarfi->SetParameter(0, 0.);     
  fgWStarfi->SetNpx(200);     
  fgWStarfi->SetNpy(20);    

  //
  //  Participants Kernel
  //
  fgWKParticipants = new TF2("WKParticipants", WKParticipants, 0., fgBMax, 0., TMath::Pi(), 3);
  fgWKParticipants->SetParameter(0, 0.);     
  fgWKParticipants->SetParameter(1, fSigmaNN);     
  fgWKParticipants->SetParameter(2, fA);     
  fgWKParticipants->SetNpx(200);     
  fgWKParticipants->SetNpy(20);      

  //
  //  Overlap and Participants
  //
  if (!mode) {
    fgWStaa = new TF1("WStaa", WStaa, 0., fgBMax, 1);
    fgWStaa->SetNpx(100);
    fgWStaa->SetParameter(0,fA);
    fgWStaa->SetNpx(100);
    fgWParticipants = new TF1("WParticipants", WParticipants, 0., fgBMax, 2);
    fgWParticipants->SetParameter(0, fSigmaNN);     
    fgWParticipants->SetParameter(1, fA);     
    fgWParticipants->SetNpx(100);
  } else {
    Info("Init","Reading overlap function from file %s",fName.Data()); 
    TFile* f = new TFile(fName.Data());
    if(!f->IsOpen()){
      Fatal("Init", "Could not open file %s",fName.Data());
    }
    fgWStaa = (TF1*) f->Get("WStaa");
    fgWParticipants = (TF1*) f->Get("WParticipants");
    delete f;
  }

  //
  //  Energy Density
  //
  fgWEnergyDensity = new TF1("WEnergyDensity", WEnergyDensity, 0., 2. * fWSr0, 1);
  fgWEnergyDensity->SetParameter(0, fWSr0 + 1.);
    
  //
  //  Geometrical Cross-Section
  //
  fgWSgeo = new TF1("WSgeo", WSgeo, 0., fgBMax, 1);
  fgWSgeo->SetParameter(0,fSigmaNN); //mbarn
  fgWSgeo->SetNpx(100);

  //
  //  Hard cross section (binary collisions)
  //
  fgWSbinary = new TF1("WSbinary", WSbinary, 0., fgBMax, 1);
  fgWSbinary->SetParameter(0, fSigmaHard); //mbarn
  fgWSbinary->SetNpx(100);

  //
  // Hard collisions per event
  //
  fgWSN = new TF1("WSN", WSN, 0., fgBMax, 1);
  fgWSN->SetNpx(100);

  //
  //  Almond shaped interaction region
  //
  fgWAlmond = new TF2("WAlmond", WAlmond, -fgBMax, fgBMax, -fgBMax, fgBMax, 1);
  fgWAlmond->SetParameter(0, 0.);     
  fgWAlmond->SetNpx(200);     
  fgWAlmond->SetNpy(200);  
  
  if(mode==2) {
    Info("Init","Reading interaction almonds from file: %s",fName.Data());
    Char_t almondName[100];
    TFile* ff = new TFile(fName.Data());
    for(Int_t k=0; k<40; k++) {
      sprintf(almondName,"WAlmondFixedB%d",k);
      fgWAlmondCurrent = (TF2*)ff->Get(almondName);
      fgWAlmondFixedB[k] = fgWAlmondCurrent;
    }
    delete ff;
  }

  fgWIntRadius = new TF1("WIntRadius", WIntRadius, 0., fgBMax, 1);
  fgWIntRadius->SetParameter(0, 0.);    

  //
  //  Path Length as a function of Phi
  //    
  fgWPathLength0 = new TF1("WPathLength0", WPathLength0, -TMath::Pi(), TMath::Pi(), 2);
  fgWPathLength0->SetParameter(0, 0.);
  fgWPathLength0->SetParameter(1, 0.); //Pathlength definition     

  fgWPathLength = new TF1("WPathLength", WPathLength, -TMath::Pi(), TMath::Pi(), 3);
  fgWPathLength->SetParameter(0, 0.);    //Impact Parameter
  fgWPathLength->SetParameter(1, 1000.); //Number of interactions used for average
  fgWPathLength->SetParameter(2, 0);     //Pathlength definition
}

void AliFastGlauber::Reset() const
{
  //
  // Reset dynamic allocated formulas
  // in case init is called twice

  if(fgWSb)            delete fgWSb;     
  if(fgWSbz)           delete fgWSbz;
  if(fgWSz)            delete fgWSz;
  if(fgWSta)           delete fgWSta;
  if(fgWStarfi)        delete fgWStarfi;
  if(fgWAlmond)        delete fgWAlmond;
  if(fgWStaa)          delete fgWStaa;
  if(fgWSgeo)          delete fgWSgeo;
  if(fgWSbinary)       delete fgWSbinary;
  if(fgWSN)            delete fgWSN;
  if(fgWPathLength0)   delete fgWPathLength0;
  if(fgWPathLength)    delete fgWPathLength;
  if(fgWEnergyDensity) delete fgWEnergyDensity;
  if(fgWIntRadius)     delete fgWIntRadius; 
  if(fgWKParticipants) delete fgWKParticipants; 
  if(fgWParticipants)  delete fgWParticipants;
}

void AliFastGlauber::DrawWSb() const
{
  //
  //  Draw Wood-Saxon Nuclear Density Function
  //
  TCanvas *c1 = new TCanvas("c1","Wood Saxon",400,10,600,700);
  c1->cd();
  Double_t max=fgWSb->GetMaximum(0,fgBMax)*1.01;
  TH2F *h2f=new TH2F("h2fwsb","Wood Saxon: #rho(r) = n (1-#omega(r/r_{0})^2)/(1+exp((r-r_{0})/d)) [fm^{-3}]",2,0,fgBMax,2,0,max);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("r [fm]");
  h2f->GetYaxis()->SetNoExponent(kTRUE);
  h2f->GetYaxis()->SetTitle("#rho [fm^{-3}]");
  h2f->Draw(); 
  fgWSb->Draw("same");
  TLegend *l1a = new TLegend(0.45,0.6,.90,0.8);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  sprintf(label,"r_{0} = %.2f fm",fWSr0);
  l1a->AddEntry(fgWSb,label,"");
  sprintf(label,"d = %.2f fm",fWSd);
  l1a->AddEntry(fgWSb,label,"");
  sprintf(label,"n = %.2e fm^{-3}",fWSn);
  l1a->AddEntry(fgWSb,label,"");
  sprintf(label,"#omega = %.2f",fWSw);
  l1a->AddEntry(fgWSb,label,"");
  l1a->Draw();
  c1->Update();
}

void AliFastGlauber::DrawOverlap() const
{
  //
  //  Draw Overlap Function
  //
  TCanvas *c2 = new TCanvas("c2","Overlap",400,10,600,700);
  c2->cd();
  Double_t max=fgWStaa->GetMaximum(0,fgBMax)*1.01;
  TH2F *h2f=new TH2F("h2ftaa","Overlap function: T_{AB} [mbarn^{-1}]",2,0,fgBMax,2,0, max);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("b [fm]");
  h2f->GetYaxis()->SetTitle("T_{AB} [mbarn^{-1}]");
  h2f->Draw(); 
  fgWStaa->Draw("same");
}

void AliFastGlauber::DrawParticipants() const
{
  //
  //  Draw Number of Participants Npart
  //
  TCanvas *c3 = new TCanvas("c3","Participants",400,10,600,700);
  c3->cd();
  Double_t max=fgWParticipants->GetMaximum(0,fgBMax)*1.01;
  TH2F *h2f=new TH2F("h2fpart","Number of Participants",2,0,fgBMax,2,0,max);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("b [fm]");
  h2f->GetYaxis()->SetTitle("N_{part}");
  h2f->Draw(); 
  fgWParticipants->Draw("same");
  TLegend *l1a = new TLegend(0.50,0.75,.90,0.9);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  sprintf(label,"#sigma^{inel.}_{NN} = %.1f mbarn",fSigmaNN);
  l1a->AddEntry(fgWParticipants,label,"");
  l1a->Draw();
  c3->Update();
}

void AliFastGlauber::DrawThickness() const
{
  //
  //  Draw Thickness Function
  //
  TCanvas *c4 = new TCanvas("c4","Thickness",400,10,600,700);
  c4->cd();
  Double_t max=fgWSta->GetMaximum(0,fgBMax)*1.01;
  TH2F *h2f=new TH2F("h2fta","Thickness function: T_{A} [fm^{-2}]",2,0,fgBMax,2,0,max);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("b [fm]");
  h2f->GetYaxis()->SetTitle("T_{A} [fm^{-2}]");
  h2f->Draw(); 
  fgWSta->Draw("same");
}

void AliFastGlauber::DrawGeo() const
{
  //
  //  Draw Geometrical Cross-Section
  //
  TCanvas *c5 = new TCanvas("c5","Geometrical Cross-Section",400,10,600,700);
  c5->cd();
  Double_t max=fgWSgeo->GetMaximum(0,fgBMax)*1.01;
  TH2F *h2f=new TH2F("h2fgeo","Differential Geometrical Cross-Section: d#sigma^{geo}_{AB}/db [fm]",2,0,fgBMax,2,0,max);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("b [fm]");
  h2f->GetYaxis()->SetTitle("d#sigma^{geo}_{AB}/db [fm]");
  h2f->Draw(); 
  fgWSgeo->Draw("same");
  TLegend *l1a = new TLegend(0.10,0.8,.40,0.9);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  sprintf(label,"#sigma_{NN}^{inel.} = %.1f mbarn",fSigmaNN);
  l1a->AddEntry(fgWSgeo,label,"");
  l1a->Draw();
  c5->Update();
}

void AliFastGlauber::DrawBinary() const
{
  //
  //  Draw Binary Cross-Section
  //
  TCanvas *c6 = new TCanvas("c6","Binary Cross-Section",400,10,600,700);
  c6->cd();
  Double_t max=fgWSbinary->GetMaximum(0,fgBMax)*1.01;
  TH2F *h2f=new TH2F("h2fbinary","Differential Binary Cross-Section: #sigma^{hard}_{NN} dT_{AB}/db [fm]",2,0,fgBMax,2,0,max);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("b [fm]");
  h2f->GetYaxis()->SetTitle("d#sigma^{hard}_{AB}/db [fm]");
  h2f->Draw(); 
  fgWSbinary->Draw("same");
  TLegend *l1a = new TLegend(0.50,0.8,.90,0.9);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  sprintf(label,"#sigma_{NN}^{hard} = %.1f mbarn",fSigmaHard);
  l1a->AddEntry(fgWSb,label,"");
  l1a->Draw();
  c6->Update();
}

void AliFastGlauber::DrawN() const
{
  //
  //  Draw Binaries per event (Ncoll)
  //
  TCanvas *c7 = new TCanvas("c7","Binaries per event",400,10,600,700);
  c7->cd();
  Double_t max=fgWSN->GetMaximum(0,fgBMax)*1.01;
  TH2F *h2f=new TH2F("h2fhardcols","Number of hard collisions: T_{AB} #sigma^{hard}_{NN}/#sigma_{AB}^{geo}",2,0,fgBMax,2,0,max);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("b [fm]");
  h2f->GetYaxis()->SetTitle("N_{coll}");
  h2f->Draw(); 
  fgWSN->Draw("same");
  TLegend *l1a = new TLegend(0.50,0.75,.90,0.9);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  sprintf(label,"#sigma^{hard}_{NN} = %.1f mbarn",fSigmaHard);
  l1a->AddEntry(fgWSN,label,"");
  sprintf(label,"#sigma^{inel.}_{NN} = %.1f mbarn",fSigmaNN);
  l1a->AddEntry(fgWSN,label,"");
  l1a->Draw();
  c7->Update();
}

void AliFastGlauber::DrawKernel(Double_t b) const
{
  //
  //  Draw Kernel
  //
  TCanvas *c8 = new TCanvas("c8","Kernel",400,10,600,700);
  c8->cd();
  fgWStarfi->SetParameter(0, b);
  TH2F *h2f=new TH2F("h2fkernel","Kernel of Overlap function: d^{2}T_{AB}/dr/d#phi [fm^{-3}]",2,0,fgBMax,2,0,TMath::Pi());
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("r [fm]");
  h2f->GetYaxis()->SetTitle("#phi [rad]");
  h2f->Draw(); 
  fgWStarfi->Draw("same");
  TLegend *l1a = new TLegend(0.65,0.8,.90,0.9);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  sprintf(label,"b = %.1f fm",b);
  l1a->AddEntry(fgWStarfi,label,"");
  l1a->Draw();
  c8->Update();
}

void AliFastGlauber::DrawAlmond(Double_t b) const
{
  //
  //  Draw Interaction Almond
  //
  TCanvas *c9 = new TCanvas("c9","Almond",400,10,600,700);
  c9->cd();
  fgWAlmond->SetParameter(0, b);
  TH2F *h2f=new TH2F("h2falmond","Interaction Almond [fm^{-4}]",2,0,fgBMax,2,0,fgBMax);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("x [fm]");
  h2f->GetYaxis()->SetTitle("y [fm]");
  h2f->Draw(); 
  fgWAlmond->Draw("same");
  TLegend *l1a = new TLegend(0.65,0.8,.90,0.9);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  sprintf(label,"b = %.1f fm",b);
  l1a->AddEntry(fgWAlmond,label,"");
  l1a->Draw();
  c9->Update();
}

void AliFastGlauber::DrawEnergyDensity() const
{
  //
  //  Draw energy density
  //
  TCanvas *c10 = new TCanvas("c10","Energy Density",400, 10, 600, 700);
  c10->cd();
  fgWEnergyDensity->SetMinimum(0.);
  Double_t max=fgWEnergyDensity->GetMaximum(0,fgWEnergyDensity->GetParameter(0))*1.01;
  TH2F *h2f=new TH2F("h2fenergydens","Energy density",2,0,fgBMax,2,0,max);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("b [fm]");
  h2f->GetYaxis()->SetTitle("fm^{-4}");
  h2f->Draw(); 
  fgWEnergyDensity->Draw("same");
  c10->Update();
}

void AliFastGlauber::DrawPathLength0(Double_t b, Int_t iopt) const
{
  //
  //  Draw Path Length
  //
  TCanvas *c11 = new TCanvas("c11","Path Length",400,10,600,700);
  c11->cd();
  fgWPathLength0->SetParameter(0, b);
  fgWPathLength0->SetParameter(1, Double_t(iopt));
  fgWPathLength0->SetMinimum(0.); 
  fgWPathLength0->SetMaximum(10.); 
  TH2F *h2f=new TH2F("h2fpathlength0","Path length",2,-TMath::Pi(), TMath::Pi(),2,0,10.);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("#phi [rad]");
  h2f->GetYaxis()->SetTitle("l [fm]");
  h2f->Draw(); 
  fgWPathLength0->Draw("same");
}

void AliFastGlauber::DrawPathLength(Double_t b , Int_t ni, Int_t iopt) const
{
  //
  //  Draw Path Length
  //
  TCanvas *c12 = new TCanvas("c12","Path Length",400,10,600,700);
  c12->cd();
  fgWAlmond->SetParameter(0, b);
  fgWPathLength->SetParameter(0, b);
  fgWPathLength->SetParameter(1, Double_t (ni));
  fgWPathLength->SetParameter(2, Double_t (iopt));
  fgWPathLength->SetMinimum(0.); 
  fgWPathLength->SetMaximum(10.); 
  TH2F *h2f=new TH2F("h2fpathlength","Path length",2,-TMath::Pi(), TMath::Pi(),2,0,10.);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("#phi [rad]");
  h2f->GetYaxis()->SetTitle("l [fm]");
  h2f->Draw(); 
  fgWPathLength->Draw("same");
}

void AliFastGlauber::DrawIntRadius(Double_t b) const
{
  //
  //  Draw Interaction Radius
  //
  TCanvas *c13 = new TCanvas("c13","Interaction Radius",400,10,600,700);
  c13->cd();
  fgWIntRadius->SetParameter(0, b);
  fgWIntRadius->SetMinimum(0);
  Double_t max=fgWIntRadius->GetMaximum(0,fgBMax)*1.01;
  TH2F *h2f=new TH2F("h2fintradius","Interaction Density",2,0.,fgBMax,2,0,max);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitle("r [fm]");
  h2f->GetYaxis()->SetTitle("[fm^{-3}]");
  h2f->Draw(); 
  fgWIntRadius->Draw("same");
}

Double_t AliFastGlauber::WSb(Double_t* x, Double_t* par)
{
  //
  //  Woods-Saxon Parameterisation
  //  as a function of radius (xx)
  //
  const Double_t kxx  = x[0];   //fm
  const Double_t kr0  = par[0]; //fm
  const Double_t kd   = par[1]; //fm   
  const Double_t kw   = par[2]; //no units
  const Double_t kn   = par[3]; //fm^-3 (used to normalize integral to one)
  Double_t y   = kn * (1.+kw*(kxx/kr0)*(kxx/kr0))/(1.+TMath::Exp((kxx-kr0)/kd));
  return y; //fm^-3
}

Double_t AliFastGlauber::WSbz(Double_t* x, Double_t* par)
{
  //
  //  Wood Saxon Parameterisation
  //  as a function of z and  b
  //
  const Double_t kbb  = x[0];   //fm
  const Double_t kzz  = x[1];   //fm
  const Double_t kr0  = par[0]; //fm
  const Double_t kd   = par[1]; //fm
  const Double_t kw   = par[2]; //no units
  const Double_t kn   = par[3]; //fm^-3 (used to normalize integral to one)
  const Double_t kxx  = TMath::Sqrt(kbb*kbb+kzz*kzz);
  Double_t y  = kn * (1.+kw*(kxx/kr0)*(kxx/kr0))/(1.+TMath::Exp((kxx-kr0)/kd));
  return y; //fm^-3
}

Double_t AliFastGlauber::WSz(Double_t* x, Double_t* par)
{
  //
  //  Wood Saxon Parameterisation
  //  as a function of z for fixed b
  //
  const Double_t kzz  = x[0];   //fm
  const Double_t kr0  = par[0]; //fm
  const Double_t kd   = par[1]; //fm
  const Double_t kw   = par[2]; //no units
  const Double_t kn   = par[3]; //fm^-3 (used to normalize integral to one)
  const Double_t kbb  = par[4]; //fm
  const Double_t kxx  = TMath::Sqrt(kbb*kbb+kzz*kzz);
  Double_t y  = kn * (1.+kw*(kxx/kr0)*(kxx/kr0))/(1.+TMath::Exp((kxx-kr0)/kd));
  return y; //fm^-3
}

Double_t AliFastGlauber::WSta(Double_t* x, Double_t* /*par*/)
{
  //
  //  Thickness function T_A
  //  as a function of b
  //
  const Double_t kb  = x[0];
  fgWSz->SetParameter(4, kb);
  Double_t y  = 2. * fgWSz->Integral(0., fgBMax);
  return y; //fm^-2
}

Double_t AliFastGlauber::WStarfi(Double_t* x, Double_t* par)
{
  //
  //  Kernel for overlap function: T_A(s)*T_A(s-b)
  //  as a function of r and phi
  const Double_t kr1  = x[0];
  const Double_t kphi = x[1];
  const Double_t kb   = par[0];
  const Double_t kr2  = TMath::Sqrt(kr1*kr1 + kb*kb - 2.*kr1*kb*TMath::Cos(kphi)); 
  Double_t y = kr1 * fgWSta->Eval(kr1) * fgWSta->Eval(kr2);
  return y; //fm^-3
}

Double_t AliFastGlauber::WStaa(Double_t* x, Double_t* par)
{
  //
  //  Overlap function 
  //  T_{AB}=Int d2s T_A(s)*T_B(s-b)
  //  as a function of b
  // (normalized to fA*fB)
  //
  const Double_t kb  = x[0];
  const Double_t ka = par[0];
  fgWStarfi->SetParameter(0, kb);

  // root integration seems to fail
  /* 
     Double_t al[2];
     Double_t bl[2];
     al[0] = 1e-6;
     al[1] = fgBMax;
     bl[0] = 0.;
     bl[1] = TMath::Pi();
     Double_t err;
     
     Double_t y =  2. *  208. * 208. * fgWStarfi->IntegralMultiple(2, al, bl, 0.001, err);
     printf("WStaa: %.5e %.5e %.5e\n", b, y, err);
  */

  //
  //  MC Integration
  //
  Double_t y = 0;
  

  for (Int_t i = 0; i < fgkMCInts; i++)
    {
	
      const Double_t kphi = TMath::Pi() * gRandom->Rndm();
      const Double_t kb1  = fgBMax      * gRandom->Rndm();	
      y += fgWStarfi->Eval(kb1, kphi);
    }
  y *= 2. * TMath::Pi() * fgBMax / fgkMCInts; //fm^-2
  y *= ka * ka * 0.1; //mbarn^-1
  return y;
}

Double_t AliFastGlauber::WKParticipants(Double_t* x, Double_t* par)
{
  //
  //  Kernel for number of participants
  //  as a function of r and phi
  //
  const Double_t kr1   = x[0];
  const Double_t kphi  = x[1];
  const Double_t kb    = par[0]; //fm
  const Double_t ksig  = par[1]; //mbarn
  const Double_t ka    = par[2]; //mass number
  const Double_t kr2   = TMath::Sqrt(kr1*kr1 +kb*kb - 2.*kr1*kb*TMath::Cos(kphi)); 
  const Double_t kxsi  = fgWSta->Eval(kr2) * ksig * 0.1; //no units
  /*
    Double_t y=(1-TMath::Power((1-xsi),aa))
   */
  Double_t a   = ka;
  Double_t sum = ka * kxsi;
  Double_t y   = sum;
  for (Int_t i = 1; i <= ka; i++)
    {
      a--;
      sum *= (-kxsi) * a / Float_t(i+1);
      y  += sum;
    }
  y *= kr1 * fgWSta->Eval(kr1);
  return y; //fm^-1
}

Double_t AliFastGlauber::WParticipants(Double_t* x, Double_t* par)
{
  //
  //  Number of Participants as 
  //  a function of b
  //
  const Double_t kb = x[0];
  const Double_t ksig  = par[0]; //mbarn
  const Double_t ka   = par[1];  //mass number
  fgWKParticipants->SetParameter(0, kb);
  fgWKParticipants->SetParameter(1, ksig);
  fgWKParticipants->SetParameter(2, ka);

  //
  //  MC Integration
  //
  Double_t y = 0;
  for (Int_t i = 0; i < fgkMCInts; i++)
    {
      const Double_t kphi = TMath::Pi() * gRandom->Rndm();
      const Double_t kb1  = fgBMax      * gRandom->Rndm();	
      y += fgWKParticipants->Eval(kb1, kphi);
    }
  y *= 2. *  ka * 2. * TMath::Pi() * fgBMax / fgkMCInts;
  return y; //no units
}

Double_t AliFastGlauber::WSgeo(Double_t* x, Double_t* par)
{
  //
  //  Geometrical Cross-Section
  //  as a function of b
  //
  const Double_t kb     = x[0];              //fm
  const Double_t ksigNN = par[0];            //mbarn
  const Double_t ktaa   = fgWStaa->Eval(kb); //mbarn^-1
  Double_t y     = 2. * TMath::Pi() * kb * (1. - TMath::Exp(- ksigNN * ktaa)); 
  return y; //fm
}

Double_t AliFastGlauber::WSbinary(Double_t* x, Double_t* par)
{
  //
  //  Number of binary hard collisions
  //  as a function of b
  //
  const Double_t kb     = x[0];              //fm
  const Double_t ksig   = par[0];            //mbarn
  const Double_t ktaa   = fgWStaa->Eval(kb); //mbarn^-1
  Double_t y = 2. * TMath::Pi() * kb * ksig * ktaa; 
  return y; //fm
}

Double_t AliFastGlauber::WSN(Double_t* x, Double_t* /*par*/)
{
  //
  //  Number of hard processes per event
  //  as a function of b
  const Double_t kb = x[0];
  Double_t y = fgWSbinary->Eval(kb)/fgWSgeo->Eval(kb);
  return y; //no units
}

Double_t AliFastGlauber::WEnergyDensity(Double_t* x, Double_t* par)
{
  //
  //  Initial energy density 
  //  as a function of the impact parameter
  //
  const Double_t kb     = x[0];
  const Double_t krA    = par[0];
  //
  //  Attention: area of transverse reaction zone in hard-sphere approximation !     
  const Double_t krA2=krA*krA;
  const Double_t kb2=kb*kb;  
  const Double_t ksaa = (TMath::Pi() - 2. * TMath::ASin(kb/ 2./ krA)) * krA2 
                      - kb * TMath::Sqrt(krA2 - kb2/ 4.); //fm^2
  const Double_t ktaa = fgWStaa->Eval(kb); //mbarn^-1
  Double_t y=ktaa/ksaa*10;
  return y; //fm^-4
}

Double_t AliFastGlauber::WAlmond(Double_t* x, Double_t* par)
{
  //
  //  Almond shaped interaction region
  //  as a function of cartesian x,y.
  //
  const Double_t kb    = par[0];
  const Double_t kxx   = x[0] + kb/2.;
  const Double_t kyy   = x[1];
  const Double_t kr1   = TMath::Sqrt(kxx*kxx + kyy*kyy);
  const Double_t kphi  = TMath::ATan2(kyy,kxx);
  const Double_t kr2   = TMath::Sqrt(kr1*kr1 + kb*kb - 2.*kr1*kb*TMath::Cos(kphi)); 
  //
  //  Interaction probability calculated as product of thicknesses
  //
  Double_t y    = fgWSta->Eval(kr1) * fgWSta->Eval(kr2);
  return y; //fm^-4
}

Double_t AliFastGlauber::WIntRadius(Double_t* x, Double_t* par)
{
  //
  //  Average interaction density over radius 
  //  at which interaction takes place
  //  as a function of radius
  //
  const Double_t kr    = x[0];
  const Double_t kb    = par[0];
  fgWAlmond->SetParameter(0, kb);
  //  Average over phi in small steps   
  const Double_t kdphi = 2. * TMath::Pi() / 100.;
  Double_t phi  = 0.;
  Double_t y    = 0.;
  for (Int_t i = 0; i < 100; i++) {
    const Double_t kxx = kr * TMath::Cos(phi);
    const Double_t kyy = kr * TMath::Sin(phi);
    y   += fgWAlmond->Eval(kxx,kyy);
    phi += kdphi;
  } // phi loop
  // Result multiplied by Jacobian (2 pi r)     
  y *= 2. * TMath::Pi() * kr / 100.;
  return y; //fm^-3
}

Double_t AliFastGlauber::WPathLength0(Double_t* x, Double_t* par)
{
  //
  //  Path Length as a function of phi 
  //  for interaction point fixed at (0,0)
  //  as a function of phi-direction
  //
  //  Phi direction in Almond
  const Double_t kphi0   = x[0];
  const Double_t kb      = par[0];
  //  Path Length definition
  const Int_t    kiopt   = Int_t(par[1]);

  //  Step along radial direction phi   
  const Int_t    kNp  = 100; // Steps in r 
  const Double_t kDr  = fgBMax/kNp;
  Double_t r  = 0.;
  Double_t rw = 0.;
  Double_t w  = 0.;
  for (Int_t i = 0; i < kNp; i++) {
    //
    //  Transform into target frame
    //
    const Double_t kxx   = r * TMath::Cos(kphi0) + kb / 2.;
    const Double_t kyy   = r * TMath::Sin(kphi0);
    const Double_t kphi  = TMath::ATan2(kyy, kxx);
    const Double_t kr1   = TMath::Sqrt(kxx*kxx + kyy*kyy);
    // Radius in projectile frame
    const Double_t kr2   = TMath::Sqrt(kr1*kr1 + kb*kb - 2.*kr1*kb*TMath::Cos(kphi)); 
    const Double_t ky    = fgWSta->Eval(kr1) * fgWSta->Eval(kr2);

    rw += ky * r;
    w  += ky;
    r  += kDr;
  } // radial steps

  Double_t y=0.;
  if (!kiopt)  // My length definition (is exact for hard disk)
    if(w) y= 2. * rw / w; 
  else {
    const Double_t knorm=fgWSta->Eval(1e-4);
    if(knorm) y =  TMath::Sqrt(2. * rw * kDr / knorm / knorm);
  }
  return y; //fm
}

Double_t AliFastGlauber::WPathLength(Double_t* x, Double_t* par)
{
  //
  //  Path Length as a function of phi 
  //  Interaction point from random distribution
  //  as a function of the phi-direction
  const Double_t kphi0   = x[0];
  const Double_t kb      = par[0];
  fgWAlmond->SetParameter(0, kb); 
  const Int_t    kNpi    = Int_t (par[1]); //Number of interactions
  const Int_t    kiopt   = Int_t(par[2]);  //Path Length definition 

  //
  //  r-steps
  // 
  const Int_t    kNp   = 100;
  const Double_t kDr  = fgBMax/Double_t(kNp);
  Double_t l = 0.;  //  Path length 
  for (Int_t in = 0; in < kNpi; in ++) {
    Double_t rw = 0.;
    Double_t w  = 0.;
    // Interaction point
    Double_t x0, y0;
    fgWAlmond->GetRandom2(x0, y0);
    // Initial radius
    const Double_t kr0  = TMath::Sqrt(x0*x0 + y0*y0);
    const Int_t    knps = Int_t ((fgBMax - kr0)/kDr) - 1;
	
    // Radial steps
    Double_t r  = 0.;
    for (Int_t i = 0; (i < knps ); i++) {
      // Transform into target frame
      const Double_t kxx   = x0 + r * TMath::Cos(kphi0) + kb / 2.;
      const Double_t kyy   = y0 + r * TMath::Sin(kphi0);
      const Double_t kphi  = TMath::ATan2(kyy, kxx);
      const Double_t kr1   = TMath::Sqrt(kxx*kxx + kyy*kyy);
      // Radius in projectile frame
      const Double_t kr2   = TMath::Sqrt(kr1*kr1 + kb*kb - 2.*kr1*kb*TMath::Cos(kphi)); 
      const Double_t ky    = fgWSta->Eval(kr1) * fgWSta->Eval(kr2);
	    
      rw += ky * r;
      w  += ky;
      r  += kDr;
    } // steps
    // Average over interactions
    if (!kiopt) {
      if(w) l += (2. * rw / w);
    } else {
      const Double_t knorm=fgWSta->Eval(1e-4);
      if(knorm) l+= 2. * rw * kDr / knorm / knorm;
    }
  } // interactions
  Double_t ret=0;
  if (!kiopt) 
    ret= l / kNpi;
  else 
    ret=TMath::Sqrt( l / kNpi);
  return ret; //fm
}

Double_t AliFastGlauber::CrossSection(Double_t b1, Double_t b2) const
{
  //
  // Return the geometrical cross-section integrated from b1 to b2 
  //
  return fgWSgeo->Integral(b1, b2)*10.; //mbarn
}

Double_t AliFastGlauber::HardCrossSection(Double_t b1, Double_t b2) const
{
  //
  // Return the hard cross-section integrated from b1 to b2 
  //
  return fgWSbinary->Integral(b1, b2)*10.; //mbarn
}

Double_t AliFastGlauber::FractionOfHardCrossSection(Double_t b1, Double_t b2) const
{
  //
  // Return fraction of hard cross-section integrated from b1 to b2 
  //
  return fgWSbinary->Integral(b1, b2)/fgWSbinary->Integral(0., 100.);
}

Double_t AliFastGlauber::NHard(Double_t b1, Double_t b2) const
{
  //
  //  Number of binary hard collisions 
  //  as a function of b (nucl/ex/0302016 eq. 19)
  //
  const Double_t kshard=HardCrossSection(b1,b2);
  const Double_t ksgeo=CrossSection(b1,b2); 
  if(ksgeo>0)
    return kshard/ksgeo;
  else return -1; 
}

Double_t AliFastGlauber::Binaries(Double_t b) const
{
  //
  // Return number of binary hard collisions normalized to 1 at b=0
  //
  if(b==0) b=1e-4;
  return fgWSN->Eval(b)/fgWSN->Eval(1e-4);
}

Double_t AliFastGlauber::MeanOverlap(Double_t b1, Double_t b2)
{
//
// Calculate the mean overlap for impact parameter range b1 .. b2
//
    Double_t sum  = 0.;
    Double_t sumc = 0.;
    Double_t b    = b1;
    
    while (b < b2-0.005) {
	Double_t  nc = GetNumberOfCollisions(b);
	sum  += 10. * fgWStaa->Eval(b) * fgWSgeo->Eval(b) * 0.01 / (1. - TMath::Exp(-nc));
	sumc += 10. * fgWSgeo->Eval(b) * 0.01;
	b += 0.01;
    }
    return (sum / CrossSection(b1, b2));
}


Double_t AliFastGlauber::MeanNumberOfCollisionsPerEvent(Double_t b1, Double_t b2)
{
//
// Calculate the mean number of collisions per event for impact parameter range b1 .. b2
//
    Double_t sum  = 0.;
    Double_t sumc = 0.;
    Double_t b    = b1;
    
    while (b < b2-0.005) {
	Double_t  nc = GetNumberOfCollisions(b);
	sum  += nc / (1. - TMath::Exp(-nc)) * 10. * fgWSgeo->Eval(b) * 0.01;
	sumc += 10. * fgWSgeo->Eval(b) * 0.01;
	b += 0.01;
    }
    return (sum / CrossSection(b1, b2));
}


Double_t AliFastGlauber::GetNumberOfBinaries(Double_t b) const
{
  //
  // Return number of binary hard collisions at b
  //
  if(b==0) b=1e-4;
  return fgWSN->Eval(b);
}

Double_t AliFastGlauber::Participants(Double_t  b) const
{
  //
  // Return the number of participants normalized to 1 at b=0
  //
  if(b==0) b=1e-4;
  return (fgWParticipants->Eval(b)/fgWParticipants->Eval(1e-4));
}

Double_t AliFastGlauber::GetNumberOfParticipants(Double_t  b) const
{
  //
  // Return the number of participants for impact parameter b
  //
  if(b==0) b=1e-4;
  return (fgWParticipants->Eval(b));
}

Double_t AliFastGlauber::GetNumberOfCollisions(Double_t  b) const
{
  //
  // Return the number of collisions for impact parameter b
  //
  if(b==0) b=1e-4;
  return (fgWStaa->Eval(b)*fSigmaNN);
}

Double_t AliFastGlauber::GetNumberOfCollisionsPerEvent(Double_t  b) const
{
  //
  // Return the number of collisions per event (at least one collision)
  // for impact parameter b
  //
    Double_t n = GetNumberOfCollisions(b);
    if (n > 0.) {
	return (n / (1. - TMath::Exp(- n)));
    } else {
	return (0.);
    }
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
  b = fgWSgeo->GetRandom();
  const Float_t kmu = fgWSN->Eval(b);
  p = 1.-TMath::Exp(-kmu);
  mult = 6000./fgWSN->Eval(1.) * kmu;
}

void AliFastGlauber::GetRandom(Int_t& bin, Bool_t& hard)
{
  //
  // Gives back a random impact parameter bin, and hard trigger decission
  //
  const Float_t kb  = fgWSgeo->GetRandom();
  const Float_t kmu = fgWSN->Eval(kb) * fSigmaHard;
  const Float_t kp  = 1.-TMath::Exp(-kmu);
  if (kb < 5.) {
    bin = 1;
  } else if (kb <  8.6) {
    bin = 2;
  } else if (kb < 11.2) {
    bin = 3;
  } else if (kb < 13.2) {
    bin = 4;
  } else if (kb < 15.0) {
    bin = 5;
  } else {
    bin = 6;
  }
  hard = kFALSE;
  const Float_t kr = gRandom->Rndm();
  if (kr < kp) hard = kTRUE;
}

Double_t  AliFastGlauber::GetRandomImpactParameter(Double_t bmin, Double_t bmax)
{
  //
  // Gives back a random impact parameter in the range bmin .. bmax
  //
  Float_t b = -1.;
  while(b < bmin || b > bmax)
    b = fgWSgeo->GetRandom();
  return b;
}

void AliFastGlauber::StoreFunctions() const
{
  //
  // Store in file functions
  //
  TFile* ff = new TFile(fName.Data(),"recreate");
  fgWStaa->Write("WStaa");
  fgWParticipants->Write("WParticipants");
  ff->Close();
  return;
}

//=================== Added by A. Dainese 11/02/04 ===========================

void AliFastGlauber::StoreAlmonds() const
{
  //
  // Store in file 
  // 40 almonds for b = (0.25+k*0.5) fm (k=0->39)
  //
  Char_t almondName[100];
  TFile* ff = new TFile(fName.Data(),"update");
  for(Int_t k=0; k<40; k++) {
    sprintf(almondName,"WAlmondFixedB%d",k);
    Double_t b = 0.25+k*0.5;
    Info("StoreAlmonds"," b = %f\n",b); 
    fgWAlmond->SetParameter(0,b);
    fgWAlmond->Write(almondName);
  }
  ff->Close();
  return;
}

void AliFastGlauber::SetCentralityClass(Double_t xsecFrLow,Double_t xsecFrUp)
{
  //
  // Set limits of centrality class as fractions 
  // of the geomtrical cross section  
  //
  if(xsecFrLow>1. || xsecFrUp>1. || xsecFrLow>xsecFrUp) {
    Error("SetCentralityClass", "Please set 0 <= xsecFrLow <= xsecFrUp <= 1\n");
    return;
  }

  Double_t bLow=0.,bUp=0.;
  Double_t xsecFr=0.;
  const Double_t knorm=fgWSgeo->Integral(0.,100.);
  while(xsecFr<xsecFrLow) {
    xsecFr = fgWSgeo->Integral(0.,bLow)/knorm;
    bLow += 0.1;
  }
  bUp = bLow;
  while(xsecFr<xsecFrUp) {
    xsecFr = fgWSgeo->Integral(0.,bUp)/knorm;
    bUp += 0.1;
  }

  Info("SetCentralityClass", "Centrality class: %4.2f-%4.2f; %4.1f < b < %4.1f fm",
	 xsecFrLow,xsecFrUp,bLow,bUp);
  fgWSbinary->SetRange(bLow,bUp);
  fBmin=bLow;
  fBmax=bUp;
  return;
}

void AliFastGlauber::GetRandomBHard(Double_t& b)
{
  //
  // Get random impact parameter according to distribution of 
  // hard (binary) cross-section, in the range defined by the centrality class
  //
  b = fgWSbinary->GetRandom();
  Int_t bin = 2*(Int_t)b;
  if( (b-(Int_t)b) > 0.5) bin++;
  fgWAlmondCurrent = fgWAlmondFixedB[bin]; 
  return;
}

void AliFastGlauber::GetRandomXY(Double_t& x,Double_t& y)
{
  //
  // Get random position of parton production point according to 
  // product of thickness functions
  //
  fgWAlmondCurrent->GetRandom2(x,y);
  return;
}

void AliFastGlauber::GetRandomPhi(Double_t& phi) 
{
  //
  // Get random parton azimuthal propagation direction
  //
  phi = 2.*TMath::Pi()*gRandom->Rndm();
  return;
}

Double_t AliFastGlauber::CalculateLength(Double_t b,Double_t x0,Double_t y0,Double_t phi0)
{
  // 
  // Calculate path length for a parton with production point (x0,y0)
  // and propagation direction (ux=cos(phi0),uy=sin(phi0)) 
  // in a collision with impact parameter b
  //

  // number of steps in l
  const Int_t    kNp  = 100;
  const Double_t kDl  = fgBMax/Double_t(kNp);

  if(fEllDef==1) {
    //
    // Definition 1:
    // 
    // ell = 2 * \int_0^\infty dl*l*(T_A*T_B)(x0+l*ux,y0+l*uy) /
    //           \int_0^\infty dl*(T_A*T_B)(x0+l*ux,y0+l*uy)
    //

    // Initial radius
    const Double_t kr0 = TMath::Sqrt(x0*x0 + y0*y0);
    const Int_t knps = Int_t ((fgBMax - kr0)/kDl) - 1;
    Double_t l  = 0.;
    Double_t integral1 = 0.;
    Double_t integral2 = 0.;
    // Radial steps
    for (Int_t i = 0; i < knps; i++) {
      
      // Transform into target frame
      const Double_t kxx   = x0 + l * TMath::Cos(phi0) + b / 2.;
      const Double_t kyy   = y0 + l * TMath::Sin(phi0);
      const Double_t kphi  = TMath::ATan2(kyy, kxx);
      const Double_t kr1   = TMath::Sqrt(kxx*kxx + kyy*kyy);
      // Radius in projectile frame
      const Double_t kr2   = TMath::Sqrt(kr1*kr1 + b*b - 2.*kr1*b*TMath::Cos(kphi)); 
      const Double_t kprodTATB = fgWSta->Eval(kr1) * fgWSta->Eval(kr2);
      
      integral1 += kprodTATB * l * kDl;
      integral2 += kprodTATB * kDl;
      l  += kDl;
    } // steps
    
    Double_t ell=0.;
    if(integral2)
      ell = (2. * integral1 / integral2);
    return ell;
  } else if(fEllDef==2) {
    //
    // Definition 2:
    // 
    // ell = \int_0^\infty dl*
    //       \Theta((T_A*T_B)(x0+l*ux,y0+l*uy)-0.5*(T_A*T_B)(0,0))
    //

    // Initial radius
    const Double_t kr0  = TMath::Sqrt(x0*x0 + y0*y0);
    const Int_t knps = Int_t ((fgBMax - kr0)/kDl) - 1;
    const Double_t kprodTATBHalfMax = 0.5*fgWAlmondCurrent->Eval(0.,0.);
    // Radial steps
    Double_t l  = 0.;
    Double_t integral = 0.;
    for (Int_t i = 0; i < knps; i++) {
      // Transform into target frame
      const Double_t kxx   = x0 + l * TMath::Cos(phi0) + b / 2.;
      const Double_t kyy   = y0 + l * TMath::Sin(phi0);
      const Double_t kphi  = TMath::ATan2(kyy, kxx);
      const Double_t kr1   = TMath::Sqrt(kxx*kxx + kyy*kyy);
      // Radius in projectile frame
      const Double_t kr2   = TMath::Sqrt(kr1*kr1 + b*b - 2.*kr1*b*TMath::Cos(kphi)); 
      const Double_t kprodTATB = fgWSta->Eval(kr1) * fgWSta->Eval(kr2);
      if(kprodTATB>kprodTATBHalfMax) integral += kDl;
      l  += kDl;
    } // steps
    Double_t ell = integral;
    return ell;
  } else {
    Error("CalculateLength","Wrong length definition setting: %d !\n",fEllDef);
    return -1.;
  }
}

void AliFastGlauber::GetLengthAndPhi(Double_t& ell,Double_t& phi,Double_t b)
{
  //
  // Return length from random b, x0, y0, phi0 
  // Return also phi0
  //
  Double_t x0,y0,phi0;
  if(b<0.) GetRandomBHard(b);
  GetRandomXY(x0,y0);
  GetRandomPhi(phi0);
  phi = phi0;
  ell = CalculateLength(b,x0,y0,phi0);
  return;
}

void AliFastGlauber::GetLength(Double_t& ell,Double_t b)
{
  //
  // Return length from random b, x0, y0, phi0 
  //
  Double_t phi;
  GetLengthAndPhi(ell,phi,b);
  return;
}

void AliFastGlauber::GetLengthsBackToBackAndPhi(Double_t& ell1,Double_t& ell2,Double_t &phi,Double_t b)
{
  //
  // Return 2 lengths back to back from random b, x0, y0, phi0 
  // Return also phi0 
 // 
  Double_t x0,y0,phi0;
  if(b<0.) GetRandomBHard(b);
  GetRandomXY(x0,y0);
  GetRandomPhi(phi0);
  const Double_t kphi0plusPi = phi0+TMath::Pi();
  phi = phi0;
  ell1 = CalculateLength(b,x0,y0,phi0);
  ell2 = CalculateLength(b,x0,y0,kphi0plusPi);
  return;
}

void AliFastGlauber::GetLengthsBackToBack(Double_t& ell1,Double_t& ell2,
					  Double_t b)
{
  //
  // Return 2 lengths back to back from random b, x0, y0, phi0 
  // 
  Double_t phi;
  GetLengthsBackToBackAndPhi(ell1,ell2,phi,b);
  return;
}

void AliFastGlauber::GetLengthsForPythia(Int_t n,Double_t* phi,Double_t* ell, Double_t b)
{
  //
  // Returns lenghts for n partons with azimuthal angles phi[n] 
  // from random b, x0, y0
  //
  Double_t x0, y0;
  if(b < 0.) GetRandomBHard(b);
  GetRandomXY(x0,y0);
  for(Int_t i = 0; i< n; i++) ell[i] = CalculateLength(b,x0,y0,phi[i]);
  return;
}

void AliFastGlauber::PlotBDistr(Int_t n)
{  
  // 
  // Plot distribution of n impact parameters
  //
  Double_t b;
  TH1F *hB = new TH1F("hB","dN/db",100,0,fgBMax); 
  hB->SetXTitle("b [fm]");
  hB->SetYTitle("dN/db [a.u.]");
  hB->SetFillColor(3);
  for(Int_t i=0; i<n; i++) {
    GetRandomBHard(b);
    hB->Fill(b);
  }
  TCanvas *cB = new TCanvas("cB","Impact parameter distribution",0,0,500,500);
  cB->cd();
  hB->Draw();
  return;
}

void AliFastGlauber::PlotLengthDistr(Int_t n,Bool_t save,const char *fname)
{
  //
  // Plot length distribution
  //
  Double_t ell;
  TH1F *hEll = new TH1F("hEll","Length distribution",64,-0.5,15); 
  hEll->SetXTitle("Transverse path length, L [fm]");
  hEll->SetYTitle("Probability");
  hEll->SetFillColor(2);
  for(Int_t i=0; i<n; i++) {
    GetLength(ell);
    hEll->Fill(ell);
  }
  hEll->Scale(1/(Double_t)n);
  TCanvas *cL = new TCanvas("cL","Length distribution",0,0,500,500);
  cL->cd();
  hEll->Draw();

  if(save) {
    TFile *f = new TFile(fname,"recreate");
    hEll->Write();
    f->Close();
  }
  return;
}

void AliFastGlauber::PlotLengthB2BDistr(Int_t n,Bool_t save,const char *fname)
{
  //
  // Plot lengths back-to-back distributions
  //
  Double_t ell1,ell2;
  TH2F *hElls = new TH2F("hElls","Lengths back-to-back",100,0,15,100,0,15); 
  hElls->SetXTitle("Transverse path length, L [fm]");
  hElls->SetYTitle("Transverse path length, L [fm]");
  for(Int_t i=0; i<n; i++) {
    GetLengthsBackToBack(ell1,ell2);
    hElls->Fill(ell1,ell2);
  }
  hElls->Scale(1/(Double_t)n);
  TCanvas *cLs = new TCanvas("cLs","Length back-to-back distribution",0,0,500,500);
  gStyle->SetPalette(1,0);
  cLs->cd();
  hElls->Draw("col,Z");
  if(save) {
    TFile *f = new TFile(fname,"recreate");
    hElls->Write();
    f->Close();
  }
  return;
}

void AliFastGlauber::PlotAlmonds() const
{
  //
  // Plot almonds for some impact parameters
  //
  TCanvas *c = new TCanvas("c","Almonds",0,0,500,500);
  gStyle->SetPalette(1,0);
  c->Divide(2,2);
  c->cd(1);
  fgWAlmondFixedB[0]->Draw("cont1");
  c->cd(2);
  fgWAlmondFixedB[10]->Draw("cont1");
  c->cd(3);
  fgWAlmondFixedB[20]->Draw("cont1");
  c->cd(4);
  fgWAlmondFixedB[30]->Draw("cont1");
  return;
}

//=================== Added by A. Dainese 05/03/04 ===========================

void AliFastGlauber::CalculateI0I1(Double_t& integral0,Double_t& integral1,
				   Double_t b,Double_t x0,Double_t y0,
                                   Double_t phi0,Double_t ellCut) const
{
  // 
  // Calculate integrals: 
  //  integral0 = \int_0^ellCut dl*(T_A*T_B)(x0+l*ux,y0+l*uy)
  //  integral1 = \int_0^ellCut dl*l*(T_A*T_B)(x0+l*ux,y0+l*uy)
  //
  // for a parton with production point (x0,y0)
  // and propagation direction (ux=cos(phi0),uy=sin(phi0)) 
  // in a collision with impact parameter b
  // 

  // number of steps in l
  const Int_t    kNp  = 100;
  const Double_t kDl  = fgBMax/Double_t(kNp);
    
  // Initial radius
  const Double_t kr0  = TMath::Sqrt(x0 * x0 + y0 * y0);
  const Int_t knps = Int_t ((fgBMax - kr0)/kDl) - 1;
    
  // Radial steps
  Double_t l  = 0.;
  integral0 = 0.;
  integral1 = 0.;
  Int_t i = 0;
  while((i < knps) && (l < ellCut)) {
    // Transform into target frame
    const Double_t kxx   = x0 + l * TMath::Cos(phi0) + b / 2.;
    const Double_t kyy   = y0 + l * TMath::Sin(phi0);
    const Double_t kphi  = TMath::ATan2(kyy, kxx);
    const Double_t kr1   = TMath::Sqrt(kxx*kxx + kyy*kyy);
    // Radius in projectile frame
    const Double_t kr2   = TMath::Sqrt(kr1*kr1 + b*b - 2.*kr1*b*TMath::Cos(kphi)); 
    const Double_t kprodTATB = fgWSta->Eval(kr1) * fgWSta->Eval(kr2);
    integral0 += kprodTATB * kDl;
    integral1 += kprodTATB * l * kDl;
    l  += kDl;
    i++;
  } // steps
  return;  
}

void AliFastGlauber::GetI0I1AndPhi(Double_t& integral0,Double_t& integral1,
				   Double_t& phi,
				   Double_t ellCut,Double_t b)
{
  //
  // Return I0 and I1 from random b, x0, y0, phi0 
  // Return also phi
  //
  Double_t x0,y0,phi0;
  if(b<0.) GetRandomBHard(b);
  GetRandomXY(x0,y0);
  GetRandomPhi(phi0);
  phi = phi0;
  CalculateI0I1(integral0,integral1,b,x0,y0,phi0,ellCut);
  return;
}

void AliFastGlauber::GetI0I1(Double_t& integral0,Double_t& integral1,
			     Double_t ellCut,Double_t b)
{
  //
  // Return I0 and I1 from random b, x0, y0, phi0 
  //
  Double_t phi;
  GetI0I1AndPhi(integral0,integral1,phi,ellCut,b);
  return;
}

void AliFastGlauber::GetI0I1BackToBackAndPhi(Double_t& integral01,Double_t& integral11,
					     Double_t& integral02,Double_t& integral12,
					     Double_t& phi,
					     Double_t ellCut,Double_t b)
{
  //
  // Return 2 pairs of I0 and I1 back to back from random b, x0, y0, phi0 
  // Return also phi0
  //
  Double_t x0,y0,phi0;
  if(b<0.) GetRandomBHard(b);
  GetRandomXY(x0,y0);
  GetRandomPhi(phi0);
  phi = phi0;
  const Double_t kphi0plusPi = phi0+TMath::Pi();
  CalculateI0I1(integral01,integral11,b,x0,y0,phi0,ellCut);
  CalculateI0I1(integral02,integral12,b,x0,y0,kphi0plusPi,ellCut);
  return;
}

void AliFastGlauber::GetI0I1BackToBackAndPhiAndXY(Double_t& integral01,Double_t& integral11,
					          Double_t& integral02,Double_t& integral12,
						  Double_t& phi,Double_t &x,Double_t &y,
						  Double_t ellCut,Double_t b)
{
  //
  // Return 2 pairs of I0 and I1 back to back from random b, x0, y0, phi0 
  // Return also phi0
  //
  Double_t x0,y0,phi0;
  if(b<0.) GetRandomBHard(b);
  GetRandomXY(x0,y0);
  GetRandomPhi(phi0);
  phi = phi0; x=x0; y=y0;
  const Double_t kphi0plusPi = phi0+TMath::Pi();
  CalculateI0I1(integral01,integral11,b,x0,y0,phi0,ellCut);
  CalculateI0I1(integral02,integral12,b,x0,y0,kphi0plusPi,ellCut);
  return;
}

void AliFastGlauber::GetI0I1BackToBack(Double_t& integral01,Double_t& integral11,
				       Double_t& integral02,Double_t& integral12,
				       Double_t ellCut,Double_t b)
{
  //
  // Return 2 pairs of I0 and I1 back to back from random b, x0, y0, phi0 
  //
  Double_t phi;
  GetI0I1BackToBackAndPhi(integral01,integral11,integral02,integral12,
			  phi,ellCut,b);
  return;
}

void AliFastGlauber::GetI0I1ForPythia(Int_t n,Double_t* phi,
				      Double_t* integral0,Double_t* integral1,
				      Double_t ellCut,Double_t b)
{
  //
  // Returns I0 and I1 pairs for n partons with azimuthal angles phi[n] 
  // from random b, x0, y0
  //
  Double_t x0,y0;
  if(b<0.) GetRandomBHard(b);
  GetRandomXY(x0,y0);
  for(Int_t i=0; i<n; i++) 
    CalculateI0I1(integral0[i],integral1[i],b,x0,y0,phi[i],ellCut);
  return;
}

void AliFastGlauber::GetI0I1ForPythiaAndXY(Int_t n,Double_t* phi,
				      Double_t* integral0,Double_t* integral1,
				      Double_t &x,Double_t& y,
				      Double_t ellCut,Double_t b)
{
  //
  // Returns I0 and I1 pairs for n partons with azimuthal angles phi[n] 
  // from random b, x0, y0 and return x0,y0
  //
  Double_t x0,y0;
  if(b<0.) GetRandomBHard(b);
  GetRandomXY(x0,y0);
  for(Int_t i=0; i<n; i++) 
    CalculateI0I1(integral0[i],integral1[i],b,x0,y0,phi[i],ellCut);
  x=x0;
  y=y0;
  return;
}

void AliFastGlauber::PlotI0I1Distr(Int_t n,Double_t ellCut,
				   Bool_t save,const char *fname)
{
  //
  // Plot I0-I1 distribution
  //
  Double_t i0,i1;
  TH2F *hI0I1s = new TH2F("hI0I1s","I_{0} versus I_{1}",1000,0,0.001,1000,0,0.01); 
  hI0I1s->SetXTitle("I_{0} [fm^{-3}]");
  hI0I1s->SetYTitle("I_{1} [fm^{-2}]");

  TH1F *hI0 = new TH1F("hI0","I_{0} = #hat{q}L / k",
		       1000,0,0.001); 
  hI0->SetXTitle("I_{0} [fm^{-3}]");
  hI0->SetYTitle("Probability");
  hI0->SetFillColor(3);
  TH1F *hI1 = new TH1F("hI1","I_{1} = #omega_{c} / k",
		       1000,0,0.01); 
  hI1->SetXTitle("I_{1} [fm^{-2}]");
  hI1->SetYTitle("Probability");
  hI1->SetFillColor(4);
  TH1F *h2 = new TH1F("h2","2 I_{1}^{2}/I_{0} = R / k",
		      100,0,0.02); 
  h2->SetXTitle("2 I_{1}^{2}/I_{0} [fm^{-1}]");
  h2->SetYTitle("Probability");
  h2->SetFillColor(2);
  TH1F *h3 = new TH1F("h3","2 I_{1}/I_{0} = L",
		      100,0,15); 
  h3->SetXTitle("2 I_{1}/I_{0} [fm]");
  h3->SetYTitle("Probability");
  h3->SetFillColor(5);
  TH1F *h4 = new TH1F("h4","I_{0}^{2}/(2 I_{1}) = #hat{q} / k",
		      100,0,0.00015); 
  h4->SetXTitle("I_{0}^{2}/(2 I_{1}) [fm^{-4}]");
  h4->SetYTitle("Probability");
  h4->SetFillColor(7);

  for(Int_t i=0; i<n; i++) {
    GetI0I1(i0,i1,ellCut);
    hI0I1s->Fill(i0,i1);
    hI0->Fill(i0);
    hI1->Fill(i1);
    h2->Fill(2.*i1*i1/i0);
    h3->Fill(2.*i1/i0);
    h4->Fill(i0*i0/2./i1);
  }
  hI0->Scale(1/(Double_t)n);
  hI1->Scale(1/(Double_t)n);
  h2->Scale(1/(Double_t)n);
  h3->Scale(1/(Double_t)n);
  h4->Scale(1/(Double_t)n);
  hI0I1s->Scale(1/(Double_t)n);

  TCanvas *cI0I1 = new TCanvas("cI0I1","I0 and I1",0,0,900,700);
  cI0I1->Divide(3,2);
  cI0I1->cd(1);
  hI0->Draw();
  cI0I1->cd(2);
  hI1->Draw();
  cI0I1->cd(3);
  h2->Draw();
  cI0I1->cd(4);
  h3->Draw();
  cI0I1->cd(5);
  h4->Draw();
  cI0I1->cd(6);
  gStyle->SetPalette(1,0);
  hI0I1s->Draw("col,Z");

  if(save) {
    TFile *f = new TFile(fname,"recreate");
    hI0I1s->Write();
    hI0->Write();
    hI1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    f->Close();
  }
  return;
}

void AliFastGlauber::PlotI0I1B2BDistr(Int_t n,Double_t ellCut,
				      Bool_t save,const char *fname)
{
  //
  // Plot I0-I1 back-to-back distributions
  //
  Double_t i01,i11,i02,i12;
  TH2F *hI0s = new TH2F("hI0s","I_{0}'s back-to-back",100,0,100,100,0,100); 
  hI0s->SetXTitle("I_{0} [fm^{-3}]");
  hI0s->SetYTitle("I_{0} [fm^{-3}]");
  TH2F *hI1s = new TH2F("hI1s","I_{1}'s back-to-back",100,0,100,100,0,100); 
  hI1s->SetXTitle("I_{1} [fm^{-2}]");
  hI1s->SetYTitle("I_{1} [fm^{-2}]");

  for(Int_t i=0; i<n; i++) {
    GetI0I1BackToBack(i01,i11,i02,i12,ellCut);
    hI0s->Fill(i01,i02);
    hI1s->Fill(i11,i12);
  }
  hI0s->Scale(1/(Double_t)n);
  hI1s->Scale(1/(Double_t)n);

  TCanvas *cI0I1s = new TCanvas("cI0I1s","I0 and I1 back-to-back distributions",0,0,800,400);
  gStyle->SetPalette(1,0);
  cI0I1s->Divide(2,1);
  cI0I1s->cd(1);
  hI0s->Draw("col,Z");
  cI0I1s->cd(2);
  hI1s->Draw("col,Z");

  if(save) {
    TFile *f = new TFile(fname,"recreate");
    hI0s->Write();
    hI1s->Write();
    f->Close();
  }
  return;
}

AliFastGlauber& AliFastGlauber::operator=(const  AliFastGlauber& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

void AliFastGlauber::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}

