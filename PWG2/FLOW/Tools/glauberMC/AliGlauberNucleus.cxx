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
//  AliGlauberNucleus implementation
//  support class for Glauber MC
//
//  origin: PHOBOS experiment
//  alification: Mikolaj Krzewicki, Nikhef, mikolaj.krzewicki@cern.ch
//
////////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TEllipse.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TF1.h>
#include <TRandom.h>
#include "AliGlauberNucleon.h"
#include "AliGlauberNucleus.h"

ClassImp(AliGlauberNucleus)

//______________________________________________________________________________
AliGlauberNucleus::AliGlauberNucleus(Option_t* iname, Int_t iN, Double_t iR, Double_t ia, Double_t iw, TF1* ifunc) : 
   fN(iN),fR(iR),fA(ia),fW(iw),fMinDist(-1),
   fF(0),fTrials(0),fFunction(ifunc),
   fNucleons(0)
{
   if (fN==0) {
      cout << "Setting up nucleus " << iname << endl;
      Lookup(iname);
   }
}

//______________________________________________________________________________
AliGlauberNucleus::~AliGlauberNucleus()
{
   if (fNucleons) {
      delete fNucleons;
   }
   delete fFunction;
}

//______________________________________________________________________________
void AliGlauberNucleus::Draw(Double_t xs, Int_t col)
{
   Double_t r = 0.5*sqrt(xs/TMath::Pi()/10.);
   TEllipse e;
   e.SetLineColor(col);
   e.SetFillColor(0);
   e.SetLineWidth(1);

   for (Int_t i = 0;i<fNucleons->GetEntries();++i) {
      AliGlauberNucleon* gn = (AliGlauberNucleon*) fNucleons->At(i);
      e.SetLineStyle(1);
      if (gn->IsSpectator()) e.SetLineStyle(3);
      e.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
   }
}

//______________________________________________________________________________
void AliGlauberNucleus::Lookup(Option_t* name)
{
   SetName(name);

   if      (TString(name) == "p")    {fN = 1;   fR = 0.6;   fA = 0;      fW =  0;      fF = 0;}
   else if (TString(name) == "d")    {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 1;}
   else if (TString(name) == "dh")   {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 3;}
   else if (TString(name) == "dhh")  {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 4;}
   else if (TString(name) == "O")    {fN = 16;  fR = 2.608; fA = 0.513;  fW = -0.051;  fF = 1;}
   else if (TString(name) == "Si")   {fN = 28;  fR = 3.34;  fA = 0.580;  fW = -0.233;  fF = 1;}
   else if (TString(name) == "S")    {fN = 32;  fR = 2.54;  fA = 2.191;  fW =  0.16;   fF = 2;}
   else if (TString(name) == "Ca")   {fN = 40;  fR = 3.766; fA = 0.586;  fW = -0.161;  fF = 1;}
   else if (TString(name) == "Ni")   {fN = 58;  fR = 4.309; fA = 0.517;  fW = -0.1308; fF = 1;}
   else if (TString(name) == "Cu")   {fN = 63;  fR = 4.2;   fA = 0.596;  fW =  0;      fF = 1;}
   else if (TString(name) == "W")    {fN = 186; fR = 6.58;  fA = 0.480;  fW =  0;      fF = 1;}
   else if (TString(name) == "Au")   {fN = 197; fR = 6.38;  fA = 0.535;  fW =  0;      fF = 1;}
   else if (TString(name) == "Pb")   {fN = 208; fR = 6.62;  fA = 0.546;  fW =  0;      fF = 1;}
   else if (TString(name) == "U")    {fN = 238; fR = 6.81;  fA = 0.6;    fW =  0;      fF = 1;}
   else {
      cout << "Could not find nucleus " << name << endl;
      return;
   }

   switch (fF)
   {
      case 0: // Proton
         fFunction = new TF1("prot","x*x*exp(-x/[0])",0,10);
         fFunction->SetParameter(0,fR);
         break;
      case 1: // 3pF
         fFunction = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,15);
         fFunction->SetParameters(fR,fA,fW);
         break;
      case 2: // 3pG
         fFunction = new TF1("3pg","x*x*(1+[2]*(x/[0])**2)/(1+exp((x**2-[0]**2)/[1]**2))",0,15);
         fFunction->SetParameters(fR,fA,fW);
         break;
      case 3: // Hulthen
         fFunction = new TF1("f3","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,10);
         fFunction->SetParameters(1/4.38,1/.85);
         break;
      case 4: // Hulthen HIJING
         fFunction = new TF1("f4","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,20);
         fFunction->SetParameters(2/4.38,2/.85);
         break;
      default:
         cerr << "Could not find function type " << fF << endl;
         return;
   }
}

//______________________________________________________________________________
void AliGlauberNucleus::SetR(Double_t ir)
{
   fR = ir;
   switch (fF)
   {
      case 0: // Proton
         fFunction->SetParameter(0,fR);
         break;
      case 1: // 3pF
         fFunction->SetParameter(0,fR);
         break;
      case 2: // 3pG
         fFunction->SetParameter(0,fR);
         break;
   }
}

//______________________________________________________________________________
void AliGlauberNucleus::SetA(Double_t ia)
{
   fA = ia;
   switch (fF)
   {
      case 0: // Proton
         break;
      case 1: // 3pF
         fFunction->SetParameter(1,fA);
         break;
      case 2: // 3pG
         fFunction->SetParameter(1,fA);
         break;
   }
}

//______________________________________________________________________________
void AliGlauberNucleus::SetW(Double_t iw)
{
   fW = iw;
   switch (fF)
   {
      case 0: // Proton
         break;
      case 1: // 3pF
         fFunction->SetParameter(2,fW);
         break;
      case 2: // 3pG
         fFunction->SetParameter(2,fW);
         break;
   }
}

//______________________________________________________________________________
void AliGlauberNucleus::ThrowNucleons(Double_t xshift)
{
   if (fNucleons==0) {
      fNucleons=new TObjArray(fN);
      fNucleons->SetOwner();
      for(Int_t i=0;i<fN;i++) {
	 AliGlauberNucleon *nucleon=new AliGlauberNucleon(); 
	 fNucleons->Add(nucleon); 
      }
   } 
   
   fTrials = 0;

   Double_t sumx=0;       
   Double_t sumy=0;       
   Double_t sumz=0;       

   Bool_t hulthen = (TString(GetName())=="dh");
   if (fN==2 && hulthen) { //special treatmeant for Hulten

      Double_t r = fFunction->GetRandom()/2;
      Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
      Double_t ctheta = 2*gRandom->Rndm() - 1 ;
      Double_t stheta = sqrt(1-ctheta*ctheta);
     
      AliGlauberNucleon *nucleon1=(AliGlauberNucleon*)(fNucleons->At(0));
      AliGlauberNucleon *nucleon2=(AliGlauberNucleon*)(fNucleons->At(1));
      nucleon1->Reset();
      nucleon1->SetXYZ(r * stheta * cos(phi) + xshift,
		       r * stheta * sin(phi),
		       r * ctheta);
      nucleon2->Reset();
      nucleon2->SetXYZ(-nucleon1->GetX() + 2*xshift,
		       -nucleon1->GetY(),
		       -nucleon1->GetZ());
      fTrials = 1;
      return;
   }

   for (Int_t i = 0; i<fN; i++) {
      AliGlauberNucleon *nucleon=(AliGlauberNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      while(1) {
         fTrials++;
         Double_t r = fFunction->GetRandom();
         Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
         Double_t ctheta = 2*gRandom->Rndm() - 1 ;
         Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
         Double_t x = r * stheta * cos(phi) + xshift;
         Double_t y = r * stheta * sin(phi);      
         Double_t z = r * ctheta;      
         nucleon->SetXYZ(x,y,z);
         if(fMinDist<0) break;
         Bool_t test=1;
         for (Int_t j = 0; j<i; j++) {
            AliGlauberNucleon *other=(AliGlauberNucleon*)fNucleons->At(j);
            Double_t xo=other->GetX();
            Double_t yo=other->GetY();
            Double_t zo=other->GetZ();
            Double_t dist = TMath::Sqrt((x-xo)*(x-xo)+
                                       (y-yo)*(y-yo)+
                                       (z-zo)*(z-zo));
	       
            if(dist<fMinDist) {
               test=0;
               break;
            }
         }
         if (test) break; //found nucleuon outside of mindist
      }
           
      sumx += nucleon->GetX();
      sumy += nucleon->GetY();
      sumz += nucleon->GetZ();
   }
      
   if(1) { // set the centre-of-mass to be at zero (+xshift)
      sumx = sumx/fN;  
      sumy = sumy/fN;  
      sumz = sumz/fN;  
      for (Int_t i = 0; i<fN; i++) {
         AliGlauberNucleon *nucleon=(AliGlauberNucleon*)(fNucleons->At(i));
         nucleon->SetXYZ(nucleon->GetX()-sumx-xshift,
                         nucleon->GetY()-sumy,
                         nucleon->GetZ()-sumz);
      }
   }
}

