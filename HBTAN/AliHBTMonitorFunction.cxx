#include "AliHBTMonitorFunction.h"
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
/******************************************************************/
/******************************************************************/

#include <Riostream.h>
ClassImp( AliHBTMonitorFunction )

AliHBTMonitorFunction::AliHBTMonitorFunction()
{
  fParticleCut = new AliHBTEmptyParticleCut(); //dummy cut
}
/******************************************************************/
AliHBTMonitorFunction::AliHBTMonitorFunction(const char* name,const char* title):TNamed(name,title)
{
  fParticleCut = new AliHBTEmptyParticleCut(); //dummy cut
}
/******************************************************************/

AliHBTMonitorFunction::~AliHBTMonitorFunction()
 {
  if (fParticleCut) delete fParticleCut;
 }
/******************************************************************/

void AliHBTMonitorFunction::
Write()
 {
   if (GetResult()) GetResult()->Write();
 }
/******************************************************************/

/******************************************************************/
void AliHBTMonitorFunction::SetParticleCut(AliHBTParticleCut* cut)
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
 fParticleCut = (AliHBTParticleCut*)cut->Clone();
 
}

/******************************************************************/

void AliHBTMonitorFunction::
Rename(const Char_t * name)
 {
 //renames the function and histograms
  SetName(name);
  SetTitle(name);
  
  TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
  GetResult()->SetName(numstr.Data());
  GetResult()->SetTitle(numstr.Data());
  
  
 }

void AliHBTMonitorFunction::
Rename(const Char_t * name, const Char_t * title)
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
AliHBTMonOneParticleFctn1D::AliHBTMonOneParticleFctn1D()
 {
   fResult = 0x0;
 }

AliHBTMonOneParticleFctn1D::
AliHBTMonOneParticleFctn1D(Int_t nbins, Double_t maxXval, Double_t minXval)
 {
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
         
   fResult   = new TH1D(numstr.Data(),numstr.Data(),nbins,minXval,maxXval);
 }

AliHBTMonOneParticleFctn1D::
AliHBTMonOneParticleFctn1D(const Char_t *name, const Char_t *title,
                    Int_t nbins, Double_t maxXval, Double_t minXval)
	:AliHBTMonOneParticleFctn(name,title)
{
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
         
   fResult   = new TH1D(numstr.Data(),numstr.Data(),nbins,minXval,maxXval);
}
/******************************************************************/
AliHBTMonOneParticleFctn1D::~AliHBTMonOneParticleFctn1D()
{
  delete fResult;
}
/******************************************************************/

void AliHBTMonOneParticleFctn1D::Process(AliHBTParticle* particle)
{
 //Fills the result
   particle = CheckParticle(particle);
   if(particle) fResult->Fill(GetValue(particle));
}
/******************************************************************/
/******************************************************************/
 
ClassImp( AliHBTMonOneParticleFctn2D )

AliHBTMonOneParticleFctn2D::
AliHBTMonOneParticleFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval , 
                    Int_t nYbins, Double_t maxYval, Double_t minYval)

{
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
        
   fResult   = new TH2D(numstr.Data(),numstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval);
	       
}	  

AliHBTMonOneParticleFctn2D::~AliHBTMonOneParticleFctn2D()
{
  delete fResult;
}
void AliHBTMonOneParticleFctn2D::Process(AliHBTParticle* particle)
{
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
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
 
   fResult   = new TH3D(numstr.Data(),numstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval,
	       nZbins,minZval,maxZval);

}	  
/******************************************************************/

AliHBTMonOneParticleFctn3D::~AliHBTMonOneParticleFctn3D()
{
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
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
         
   fResult   = new TH1D(numstr.Data(),numstr.Data(),
                           nbins,minval,maxval);
	       
 }

AliHBTMonTwoParticleFctn1D::
AliHBTMonTwoParticleFctn1D(const Char_t* name, const Char_t* title,
                    Int_t nbins, Double_t maxval, Double_t minval)
	:AliHBTMonTwoParticleFctn(name,title)
 {
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram

   fResult   = new TH1D(numstr.Data(),numstr.Data(),
                           nbins,minval,maxval);
	       
 }


/******************************************************************/
AliHBTMonTwoParticleFctn1D::~AliHBTMonTwoParticleFctn1D()
{
  delete fResult;
}
/******************************************************************/
void AliHBTMonTwoParticleFctn1D::
Process(AliHBTParticle* trackparticle, AliHBTParticle* partparticle)
{
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
   TString numstr = fName + " Result";  //title and name of the 
                                           //result histogram
         
   fResult   = new TH2D(numstr.Data(),numstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval);
	       
}	  
/******************************************************************/
AliHBTMonTwoParticleFctn2D::~AliHBTMonTwoParticleFctn2D()
{
  delete fResult;
}
/******************************************************************/
void AliHBTMonTwoParticleFctn2D::
Process(AliHBTParticle* trackparticle, AliHBTParticle* partparticle)
{
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
Process(AliHBTParticle* trackparticle, AliHBTParticle* partparticle)
{
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

