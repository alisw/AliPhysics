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

#include "AliLog.h"
#include "AliHBTMonitorFunction.h"

//__________________________________________________________________
////////////////////////////////////////////////////////////////////
//
// class AliHBTMonitorFunction
//
// class AliHBTMonOneParticleFctn
// class AliHBTMonTwoParticleFctn
//
// class AliHBTMonOneParticleFctn1D
// class AliHBTMonOneParticleFctn2D
// class AliHBTMonOneParticleFctn3D
//
// class AliHBTMonTwoParticleFctn1D
// class AliHBTMonTwoParticleFctn2D
// class AliHBTMonTwoParticleFctn3D
//
// Base Classes for monitoring functions
// author: chajecki@if.pw.edu.pl
//
/******************************************************************/
/*
Base classes for monitor functions

          monitor function
               /    \
              /      \
             /        \
            /          \
           /            \
          /              \
         /                \
   one particle     two particle  
     /  |  \            /  |  \
    /   |   \          /   |   \
   1D  2D   3D        1D  2D   3D

Zbigniew.Chajecki@cern.ch

*/
///////////////////////////////////////////////////////////////////////


ClassImp( AliHBTMonitorFunction )

AliHBTMonitorFunction::AliHBTMonitorFunction():
 fParticleCut(new AliAODParticleEmptyCut())
{
  //ctor
}
/******************************************************************/
AliHBTMonitorFunction::AliHBTMonitorFunction(const char* name,const char* title):
 TNamed(name,title),
 fParticleCut(new AliAODParticleEmptyCut())
{
  //ctor
}
/******************************************************************/
AliHBTMonitorFunction::AliHBTMonitorFunction(const AliHBTMonitorFunction& /*in*/):
 TNamed(),
 fParticleCut(new AliAODParticleEmptyCut())
{
  //cpy ctor
  // We cannot copy because it is a mess with names (histogram and functions)
  MayNotUse("AliHBTMonitorFunction(const AliHBTMonitorFunction&");
}
/******************************************************************/
AliHBTMonitorFunction& AliHBTMonitorFunction::operator=(const AliHBTMonitorFunction& /*in*/)
{
  //assigment operator 
  // We cannot copy because it is a mess with names (histogram and functions)
  MayNotUse("operator=");
  return *this;
}
/******************************************************************/

AliHBTMonitorFunction::~AliHBTMonitorFunction()
 {
  //dtor
  delete fParticleCut;
 }
/******************************************************************/

void AliHBTMonitorFunction::Write()
 {
   //Writes an function to disk
   if (GetResult()) GetResult()->Write();
 }
/******************************************************************/

void AliHBTMonitorFunction::Init()
 {
   //Writes an function to disk
   AliDebug(1,"Entering");
   
   if (GetResult() == 0x0)
    {
      Warning("Init","Function has NULL result histogram!");
      return;
    }
   GetResult()->Reset();
   GetResult()->SetDirectory(0x0);
   AliDebug(1,"Done");
 }
/******************************************************************/

void AliHBTMonitorFunction::SetParticleCut(AliAODParticleCut* cut)
{
//Sets new Particle Cut. Old one is deleted
//Note that it is created new object instead of simple pointer set
//I do not want to have pointer 
//to object created somewhere else
//because in that case I could not believe that 
//it would always exist (sb could delete it)
//so we have always own copy

 if(!cut) 
   {
     Error("AliHBTMonitorFunction::SetParticleCut","argument is NULL");
     return;
   }
 delete fParticleCut;
 fParticleCut = (AliAODParticleCut*)cut->Clone();
 
}
/******************************************************************/

void AliHBTMonitorFunction::Rename(const Char_t * name)
 {
 //renames the function and histograms
  SetName(name);
  SetTitle(name);
  
  TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
  GetResult()->SetName(numstr.Data());
  GetResult()->SetTitle(numstr.Data());
 }
/******************************************************************/

void AliHBTMonitorFunction::Rename(const Char_t * name, const Char_t * title)
 {
 //renames and retitle the function and histograms
 
  SetName(name);
  SetTitle(title);
  
  TString numstrn = fName + " Result";  //name of the 
                                           //result histogram

  TString numstrt = fTitle + " Result";  //title of the 
                                           //result histogram
		   

  GetResult()->SetName(numstrn.Data());
  GetResult()->SetTitle(numstrt.Data());

 }

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTMonOneParticleFctn )  //z.ch.
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTMonTwoParticleFctn )  //z.ch.
/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTMonOneParticleFctn1D )
AliHBTMonOneParticleFctn1D::AliHBTMonOneParticleFctn1D():
 fResult(0x0)
 {
   //ctor
 }
/******************************************************************/

AliHBTMonOneParticleFctn1D::AliHBTMonOneParticleFctn1D(Int_t nbins, Double_t maxXval, Double_t minXval)
 {
   //ctor
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
   fResult   = new TH1D(numstr.Data(),numstr.Data(),nbins,minXval,maxXval);
   fResult->Sumw2();
 }

AliHBTMonOneParticleFctn1D::AliHBTMonOneParticleFctn1D(const Char_t *name, const Char_t *title,
                    Int_t nbins, Double_t maxXval, Double_t minXval)
	:AliHBTMonOneParticleFctn(name,title)
{
  //ctor
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
         
   fResult   = new TH1D(numstr.Data(),numstr.Data(),nbins,minXval,maxXval);
   fResult->Sumw2();
   fResult->SetDirectory(0x0);
   
}
/******************************************************************/
AliHBTMonOneParticleFctn1D::~AliHBTMonOneParticleFctn1D()
{
 //dtor
 delete fResult;
}
/******************************************************************/

void AliHBTMonOneParticleFctn1D::Process(AliVAODParticle* particle)
{
 //Fills the result
   particle = CheckParticle(particle);
   if(particle) fResult->Fill(GetValue(particle));
}
/******************************************************************/
/******************************************************************/
 
ClassImp( AliHBTMonOneParticleFctn2D )

AliHBTMonOneParticleFctn2D::AliHBTMonOneParticleFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval , 
                    Int_t nYbins, Double_t maxYval, Double_t minYval)

{
  //ctor
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
        
   fResult   = new TH2D(numstr.Data(),numstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval);
   fResult->Sumw2();
   fResult->SetDirectory(0x0);
}	  
/******************************************************************/

AliHBTMonOneParticleFctn2D::~AliHBTMonOneParticleFctn2D()
{
  //dtor
  delete fResult;
}
void AliHBTMonOneParticleFctn2D::Process(AliVAODParticle* particle)
{
  //fills the function for one particle
  particle = CheckParticle(particle);
  if(particle) 
   { 
     Double_t x,y;
     GetValues(particle,x,y);
     fResult->Fill(x,y);
   }
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTMonOneParticleFctn3D)

AliHBTMonOneParticleFctn3D::
AliHBTMonOneParticleFctn3D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                    Int_t nYbins, Double_t maxYval, Double_t minYval, 
                    Int_t nZbins, Double_t maxZval, Double_t minZval)

{
  //ctor
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
 
   fResult   = new TH3D(numstr.Data(),numstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval,
	       nZbins,minZval,maxZval);
   fResult->Sumw2();
   fResult->SetDirectory(0x0);

}	  
/******************************************************************/

AliHBTMonOneParticleFctn3D::~AliHBTMonOneParticleFctn3D()
{
  //dtor
  delete fResult;
}
/******************************************************************/


/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTMonTwoParticleFctn1D)

AliHBTMonTwoParticleFctn1D::
AliHBTMonTwoParticleFctn1D(Int_t nbins, Double_t maxval, Double_t minval)
 {
   //ctor
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
         
   fResult   = new TH1D(numstr.Data(),numstr.Data(),
                           nbins,minval,maxval);
   fResult->Sumw2();
   fResult->SetDirectory(0x0);
 }

AliHBTMonTwoParticleFctn1D::
AliHBTMonTwoParticleFctn1D(const Char_t* name, const Char_t* title,
                    Int_t nbins, Double_t maxval, Double_t minval)
	:AliHBTMonTwoParticleFctn(name,title)
 {
   //ctor
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram

   fResult   = new TH1D(numstr.Data(),numstr.Data(),
                           nbins,minval,maxval);
   fResult->Sumw2();
   fResult->SetDirectory(0x0);
 }


/******************************************************************/
AliHBTMonTwoParticleFctn1D::~AliHBTMonTwoParticleFctn1D()
{
  //dtor
  delete fResult;
}
/******************************************************************/
void AliHBTMonTwoParticleFctn1D::
Process(AliVAODParticle* trackparticle, AliVAODParticle* partparticle)
{
  //fills the function for one particle
  partparticle  = CheckParticle(partparticle);
  trackparticle = CheckParticle(trackparticle);
  if( partparticle && trackparticle) 
   { 
     Double_t x = GetValue(trackparticle,partparticle);
     fResult->Fill(x);
   }
}
/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTMonTwoParticleFctn2D)


AliHBTMonTwoParticleFctn2D::
AliHBTMonTwoParticleFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval , 
                    Int_t nYbins, Double_t maxYval, Double_t minYval)

{
  //ctor
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
         
   fResult   = new TH2D(numstr.Data(),numstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval);
   fResult->Sumw2();
   fResult->SetDirectory(0x0);
}	  
/******************************************************************/
AliHBTMonTwoParticleFctn2D::~AliHBTMonTwoParticleFctn2D()
{
  //dtor
  delete fResult;
}
/******************************************************************/
void AliHBTMonTwoParticleFctn2D::
Process(AliVAODParticle* trackparticle, AliVAODParticle* partparticle)
{
  //fills the function for one particle
  partparticle  = CheckParticle(partparticle);
  trackparticle = CheckParticle(trackparticle);
  if( partparticle && trackparticle) 
   { 
     Double_t x,y;
     GetValues(trackparticle,partparticle,x,y);
     fResult->Fill(x,y);
   }
}
/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp(AliHBTMonTwoParticleFctn3D)

void AliHBTMonTwoParticleFctn3D::
Process(AliVAODParticle* trackparticle, AliVAODParticle* partparticle)
{
  //fills the function for one particle
  partparticle  = CheckParticle(partparticle);
  trackparticle = CheckParticle(trackparticle);
  if( partparticle && trackparticle) 
   { 
     Double_t x,y,z;
     GetValues(trackparticle,partparticle,x,y,z);
     fResult->Fill(x,y,z);
   }
}
/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/

