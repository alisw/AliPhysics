#include "AliHBTFunction.h"

/* $Id: */

//--------------------------------------------------------------------
//AliHBTFunction
//Author: Piotr Krzysztof Skowronski
//Piotr.Skowronski@cern.ch
//Base classes for HBT functions
/*
 OnePairFctn       Function                  TwoPairFctn
   | \    \        /  |   \                      /|\   
   |  \    \      /   |    \                    / | \    
   |   \    \   1D   2D    3D                  /  |  \    
   |    \    \  / \   |\   | \________________/__ |   \             
   |     \    \/   \  | \__|_____________    /   \|    \             
   |      \    \    \_|____|__           \  /     |     \            
   |       \ /  \___  |    |  \           \/      |\     \              
   |        /       \ |    |   \          /\      | \     \             
   |       / \       \|    |    \        /  \     |  \     \              
   |      /   \      /\    |     \      /    \    |   \     |              
   |     /     \    /  \   |      \    /      \   |    \    |               
   |    /       \  /    \  |       \  /        \  |     \   |
   |   /         \/      \ |        \/          \ |      \  |               
 OnePair1D   OnePair2D  OnePair3D  TwoPair1D  TwoPair2D TwoPair3D


 four particle functions are intendent to be resolution functions:
 it is mecessary to have simulated particle pair corresponding to given
 recontructed track pair in order to calculate function simualted value 
 and recontructed value, to be further histogrammed
*/
//--------------------------------------------------------------------- 

#include <TH2.h>
#include <TH3.h>

/******************************************************************/
/******************************************************************/

ClassImp( AliHBTFunction )

AliHBTFunction::AliHBTFunction():
  fPairCut(new AliHBTEmptyPairCut())   //dummy cut
{
//Default constructor
}
/******************************************************************/
AliHBTFunction::AliHBTFunction(const char* name,const char* title):
  TNamed(name,title),
  fPairCut(new AliHBTEmptyPairCut()) //dummy cut  
{
//Constructor  
}
/******************************************************************/

AliHBTFunction::~AliHBTFunction()
{
//destructor  
  delete fPairCut;
}
/******************************************************************/

void AliHBTFunction::WriteFunction()
{
//writes result of the function to file
   if (GetNumerator()) GetNumerator()->Write();
   if (GetDenominator()) GetDenominator()->Write();
   TH1* res = GetResult();
   if (res) res->Write();
}
/******************************************************************/

TH1* AliHBTFunction::GetRatio(Double_t normfactor)
 {
 //returns ratio of numerator and denominator
 //
   if (gDebug>0) Info("GetRatio","Norm. Factor is %f for %s",normfactor,GetName());
   
   if (normfactor == 0.0)
    {
      Error("GetRatio","Scaling Factor is 0. Null poiner returned");
      return 0x0;
    }
   TString str = fName + " ratio";
   TH1 *result = (TH1*)GetNumerator()->Clone(str.Data());
   
   result->SetTitle(str.Data());
   
   result->Divide(GetNumerator(),GetDenominator(),normfactor);
   
   return result;
   
 }
/******************************************************************/
void AliHBTFunction::SetPairCut(AliHBTPairCut* cut)
{
//Sets new Pair Cut. Old one is deleted
//Note that it is created new object instead of simple pointer set
//I do not want to have pointer 
//to object created somewhere else
//because in that case I could not believe that 
//it would always exist (sb could delete it)
//so we have always own copy

 if(!cut) 
   {
     Error("AliHBTFunction::SetPairCut","argument is NULL");
     return;
   }
 delete fPairCut;
 fPairCut = (AliHBTPairCut*)cut->Clone();
 
}

/******************************************************************/

void AliHBTFunction::Rename(const Char_t * name)
 {
 //renames the function and histograms
  SetName(name);
  SetTitle(name);
  
  TString numstr = fName + " Numerator";  //title and name of the 
                                           //numerator histogram
  TString denstr = fName + " Denominator";//title and name of the 
                                           //denominator histogram
 
  GetNumerator()->SetName(numstr.Data());
  GetNumerator()->SetTitle(numstr.Data());
  
  GetDenominator()->SetName(denstr.Data());
  GetDenominator()->SetTitle(denstr.Data());
  
 }

void AliHBTFunction::Rename(const Char_t * name, const Char_t * title)
 {
 //renames and retitle the function and histograms
 
  SetName(name);
  SetTitle(title);
  
  TString numstrn = fName + " Numerator";  //name of the 
                                           //numerator histogram

  TString numstrt = fTitle + " Numerator";  //title of the 
                                           //numerator histogram
		   
  TString denstrn = fName + " Denominator";//name of the 
                                           //denominator histogram

  TString denstrt = fTitle + " Denominator";//title of the 
                                           //denominator histogram
		   
 
  GetNumerator()->SetName(numstrn.Data());
  GetNumerator()->SetTitle(numstrt.Data());
  
  GetDenominator()->SetName(denstrn.Data());
  GetDenominator()->SetTitle(denstrt.Data());


 }
/******************************************************************/

void AliHBTFunction::InitFunction()
{
//Iniotializes fctn.: Resets histograms
//In case histograms are not created in ctor, builds with default parameters
  if ( !(GetNumerator()&&GetDenominator()) ) BuildHistos();
  GetNumerator()->Reset();
  GetDenominator()->Reset();
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTOnePairFctn )

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTTwoPairFctn)

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTFunction1D )

const Int_t AliHBTFunction1D::fgkDefaultNBins = 100;
const Float_t AliHBTFunction1D::fgkDefaultMin = 0.0;
const Float_t AliHBTFunction1D::fgkDefaultMax = 0.15;
const UInt_t AliHBTFunction1D::fgkDefaultNBinsToScale = 30;


AliHBTFunction1D::AliHBTFunction1D():
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScale(fgkDefaultNBinsToScale)
{//default constructor
}
/******************************************************************/

AliHBTFunction1D::AliHBTFunction1D(Int_t nbins, Float_t maxXval, Float_t minXval):
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScale(fgkDefaultNBinsToScale)
{
 //Constructor of Two Part One Dimentional Function 
 // nbins: number of bins in histograms - default 100
 // maxXval and minXval: range of histgram(s) default 0 - 0.15 (GeV)
 BuildHistos(nbins,maxXval,minXval);
}
/******************************************************************/
AliHBTFunction1D::AliHBTFunction1D(const Char_t *name, const Char_t *title):
 AliHBTFunction(name,title),
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScale(fgkDefaultNBinsToScale)
{//constructor
}
/******************************************************************/
AliHBTFunction1D::AliHBTFunction1D(const Char_t *name, const Char_t *title,
                                   Int_t nbins, Float_t maxXval, Float_t minXval):
 AliHBTFunction(name,title),
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScale(fgkDefaultNBinsToScale)
{
//constructor
  BuildHistos(nbins,maxXval,minXval);
}
/******************************************************************/

AliHBTFunction1D::~AliHBTFunction1D()
{
//destructor
  delete fNumerator;
  delete fDenominator;
}
/******************************************************************/
void AliHBTFunction1D::BuildHistos()
{
//builds histograms with default settings
 BuildHistos(fgkDefaultNBins,fgkDefaultMax,fgkDefaultMin);
}

/******************************************************************/

void AliHBTFunction1D::BuildHistos(Int_t nbins, Float_t max, Float_t min)
{
//builds numarator and denominator hitograms
  TString numstr = fName + " Numerator";  //title and name of the 
                                          //numerator histogram
  TString denstr = fName + " Denominator";//title and name of the 
                                          //denominator histogram
  fNumerator   = new TH1D(numstr.Data(),numstr.Data(),nbins,min,max);
  fDenominator = new TH1D(denstr.Data(),denstr.Data(),nbins,min,max);
   
  fNumerator->Sumw2();
  fDenominator->Sumw2();
}
/******************************************************************/

Double_t AliHBTFunction1D::Scale()
{
 //Calculates the factor that should be used to scale 
 //quatience of fNumerator and fDenominator to 1 at tail
 return Scale(fNumerator,fDenominator);
}
/******************************************************************/

Double_t AliHBTFunction1D::Scale(TH1D* num,TH1D* den)
{
 //Calculates the factor that should be used to scale 
 //quatience of num and den to 1 at tail
 
  if (gDebug>0) Info("Scale","Enetered Scale()");
  if(!num) 
   {
     Error("Scale","No numerator");
     return 0.0;
   }
  if(!den) 
   {
     Error("Scale","No denominator");
     return 0.0;
   }
  
  if(fNBinsToScale < 1) 
   {
    return 0.0;
    Error("Scale","Number of bins for scaling is smaller thnan 1");
   }
  UInt_t nbins = num->GetNbinsX();
  if (fNBinsToScale > nbins) 
   {
    Error("Scale","Number of bins for scaling is bigger thnan number of bins in histograms");
    return 0.0;
   }
  if (gDebug>0) Info("Scale","No errors detected");

  Double_t ratio;
  Double_t sum = 0;
  Int_t n = 0;
  
  Int_t offset = nbins - fNBinsToScale - 1; 

  for (UInt_t i = offset; i< nbins; i++)
   {
    if ( num->GetBinContent(i) > 0.0 )
     {
       ratio = den->GetBinContent(i)/num->GetBinContent(i);
       sum += ratio;
       n++;
     }
   }
  
  if(gDebug > 0) Info("Scale","sum=%f fNBinsToScale=%d n=%d",sum,fNBinsToScale,n);
  
  if (n == 0) return 0.0;
  Double_t ret = sum/((Double_t)n);

  if(gDebug > 0) Info("Scale","returning %f",ret);
  return ret;
} 

/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTFunction2D                                  //
//                                                   //
// Base Calss for 2-dimensinal Functions             //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

ClassImp( AliHBTFunction2D )

const Int_t AliHBTFunction2D::fgkDefaultNBinsX = 200;//default number of Bins in X axis in histograms
const Float_t AliHBTFunction2D::fgkDefaultMinX = 0.0;//Default min value of X axis in histograms
const Float_t AliHBTFunction2D::fgkDefaultMaxX = 1.5;//Default max value of X axis inhistograms

const Int_t AliHBTFunction2D::fgkDefaultNBinsY = 200;//default number of Bins in histograms
const Float_t AliHBTFunction2D::fgkDefaultMinY = -0.15;//Default min value of histograms
const Float_t AliHBTFunction2D::fgkDefaultMaxY =  0.15;//Default max value of histograms

const UInt_t AliHBTFunction2D::fgkDefaultNBinsToScaleX = 30;//Default number of X bins used for scaling to tale
const UInt_t AliHBTFunction2D::fgkDefaultNBinsToScaleY = 30;//Default number of bins used for scaling to tale

/******************************************************************/
AliHBTFunction2D::AliHBTFunction2D():
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScaleX(fgkDefaultNBinsToScaleX), 
 fNBinsToScaleY(fgkDefaultNBinsToScaleY)
{//default constructor
}
/******************************************************************/
AliHBTFunction2D::AliHBTFunction2D(const Char_t *name, const Char_t *title):
 AliHBTFunction(name,title),
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScaleX(fgkDefaultNBinsToScaleX), 
 fNBinsToScaleY(fgkDefaultNBinsToScaleY)
{//constructor
}
/******************************************************************/

AliHBTFunction2D::AliHBTFunction2D(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                   Int_t nYbins, Double_t maxYval, Double_t minYval):
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScaleX(fgkDefaultNBinsToScaleX), 
 fNBinsToScaleY(fgkDefaultNBinsToScaleY)
{
  BuildHistos(nXbins,maxXval,minXval,nYbins,maxYval,minYval);
}	  
/******************************************************************/

AliHBTFunction2D::AliHBTFunction2D(const Char_t *name, const Char_t *title,
                                   Int_t nXbins, Double_t maxXval, Double_t minXval,
                                   Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTFunction(name,title),
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScaleX(fgkDefaultNBinsToScaleX), 
 fNBinsToScaleY(fgkDefaultNBinsToScaleY)
{
  BuildHistos(nXbins,maxXval,minXval,nYbins,maxYval,minYval);
}	  
/******************************************************************/

AliHBTFunction2D::~AliHBTFunction2D()
{
  delete fNumerator;
  delete fDenominator;
}
/******************************************************************/

void AliHBTFunction2D::BuildHistos()
{
//Creates default histograms
  BuildHistos(fgkDefaultNBinsX,fgkDefaultMaxX,fgkDefaultMinX,
              fgkDefaultNBinsY,fgkDefaultMaxY,fgkDefaultMinY);
}
/******************************************************************/

void AliHBTFunction2D::BuildHistos(Int_t nxbins, Float_t xmax, Float_t xmin,
                                   Int_t nybins, Float_t ymax, Float_t ymin)
{
  //Builds numerator and denominator histograms (2d-case)
 TString numstr = fName + " Numerator";  //title and name of the 
                                           //numerator histogram
 TString denstr = fName + " Denominator";//title and name of the 
                                           //denominator histogram
         
 fNumerator   = new TH2D(numstr.Data(),numstr.Data(),
                         nxbins,xmin,xmax,nybins,ymin,ymax);
	       
 fDenominator = new TH2D(denstr.Data(),denstr.Data(),
                         nxbins,xmin,xmax,nybins,ymin,ymax);
 
 fNumerator->Sumw2();
 fDenominator->Sumw2();
}
/******************************************************************/

void AliHBTFunction2D::SetNumberOfBinsToScale(UInt_t xn, UInt_t yn)
{
//defines area used for scaling factor calculation
  fNBinsToScaleX = xn;
  fNBinsToScaleY = yn;
}
/******************************************************************/

Double_t AliHBTFunction2D::Scale()
{
  // Calculates the factor that should be used to scale 
  // quatience of fNumerator and fDenominator to 1 at 
  // given region
  if (gDebug>0) Info("Scale","Enetered Scale()");
  if(!fNumerator) 
   {
     Error("Scale","No numerator");
     return 0.0;
   }
  if(!fDenominator) 
   {
     Error("Scale","No denominator");
     return 0.0;
   }
  
  if( (fNBinsToScaleX < 1) || (fNBinsToScaleY < 1) ) 
   {
    return 0.0;
    Error("Scale","Number of bins for scaling is smaller thnan 1");
   }
  UInt_t nbinsX = fNumerator->GetNbinsX();
  if (fNBinsToScaleX > nbinsX) 
   {
    Error("Scale","Number of X bins for scaling is bigger thnan number of bins in histograms");
    return 0.0;
   }
   
  UInt_t nbinsY = fNumerator->GetNbinsX();
  if (fNBinsToScaleY > nbinsY) 
   {
    Error("Scale","Number of Y bins for scaling is bigger thnan number of bins in histograms");
    return 0.0;
   }

  if (gDebug>0) Info("Scale","No errors detected");

  Int_t offsetX = nbinsX - fNBinsToScaleX - 1; //bin that we start loop over bins in axis X
  Int_t offsetY = nbinsY - fNBinsToScaleY - 1; //bin that we start loop over bins in axis X

  Double_t ratio;
  Double_t sum = 0;
  Int_t n = 0;
  
  for (UInt_t j = offsetY; j< nbinsY; j++)
    for (UInt_t i = offsetX; i< nbinsX; i++)
     {
      if ( fNumerator->GetBinContent(i,j) > 0.0 )
       {
         ratio = fDenominator->GetBinContent(i,j)/fNumerator->GetBinContent(i,j);
         sum += ratio;
         n++;
       }
     }
  
  if(gDebug > 0) Info("Scale","sum=%f fNBinsToScaleX=%d fNBinsToScaleY=%d n=%d",sum,fNBinsToScaleX,fNBinsToScaleY,n);
  
  if (n == 0) return 0.0;
  Double_t ret = sum/((Double_t)n);

  if(gDebug > 0) Info("Scale","returning %f",ret);
  return ret;
} 

/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTFunction3D                                  //
//                                                   //
// Base Calss for 3-dimensinal Functions             //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

ClassImp( AliHBTFunction3D)

const Int_t AliHBTFunction3D::fgkDefaultNBinsX = 200;//default number of Bins in X axis in histograms
const Float_t AliHBTFunction3D::fgkDefaultMinX = 0.0;//Default min value of X axis in histograms
const Float_t AliHBTFunction3D::fgkDefaultMaxX = 1.5;//Default max value of X axis inhistograms

const Int_t AliHBTFunction3D::fgkDefaultNBinsY = 200;//default number of Bins in histograms
const Float_t AliHBTFunction3D::fgkDefaultMinY = -0.15;//Default min value of histograms
const Float_t AliHBTFunction3D::fgkDefaultMaxY =  0.15;//Default max value of histograms

const Int_t AliHBTFunction3D::fgkDefaultNBinsZ = 200;//default number of Bins in histograms
const Float_t AliHBTFunction3D::fgkDefaultMinZ = -0.15;//Default min value of histograms
const Float_t AliHBTFunction3D::fgkDefaultMaxZ =  0.15;//Default max value of histograms

const UInt_t AliHBTFunction3D::fgkDefaultNBinsToScaleX = 30;//Default number of X bins used for scaling to tale
const UInt_t AliHBTFunction3D::fgkDefaultNBinsToScaleY = 30;//Default number of bins used for scaling to tale
const UInt_t AliHBTFunction3D::fgkDefaultNBinsToScaleZ = 30;//Default number of bins used for scaling to tale

AliHBTFunction3D::AliHBTFunction3D():
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScaleX(fgkDefaultNBinsToScaleX), 
 fNBinsToScaleY(fgkDefaultNBinsToScaleY),
 fNBinsToScaleZ(fgkDefaultNBinsToScaleZ)
{
 //constructor
}
/******************************************************************/

AliHBTFunction3D::AliHBTFunction3D(const Char_t *name, const Char_t *title):
 AliHBTFunction(name,title),
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScaleX(fgkDefaultNBinsToScaleX), 
 fNBinsToScaleY(fgkDefaultNBinsToScaleY),
 fNBinsToScaleZ(fgkDefaultNBinsToScaleZ)
{
 //constructor
}  
/******************************************************************/

AliHBTFunction3D::AliHBTFunction3D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                  Int_t nYbins, Double_t maxYval, Double_t minYval, 
                  Int_t nZbins, Double_t maxZval, Double_t minZval):
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScaleX(fgkDefaultNBinsToScaleX), 
 fNBinsToScaleY(fgkDefaultNBinsToScaleY),
 fNBinsToScaleZ(fgkDefaultNBinsToScaleZ)
{
 //constructor
 BuildHistos( nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval);
}	  
/******************************************************************/

AliHBTFunction3D::AliHBTFunction3D(const Char_t *name, const Char_t *title,
                 Int_t nXbins, Double_t maxXval, Double_t minXval, 
                 Int_t nYbins, Double_t maxYval, Double_t minYval, 
                 Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTFunction(name,title),
 fNumerator(0x0),
 fDenominator(0x0),
 fNBinsToScaleX(fgkDefaultNBinsToScaleX), 
 fNBinsToScaleY(fgkDefaultNBinsToScaleY),
 fNBinsToScaleZ(fgkDefaultNBinsToScaleZ)
{
 //constructor
 BuildHistos( nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval);
}	  
/******************************************************************/


AliHBTFunction3D::~AliHBTFunction3D()
{
  delete fNumerator;
  delete fDenominator;
}
/******************************************************************/

void AliHBTFunction3D::BuildHistos()
{
//Creates default histograms
  BuildHistos(fgkDefaultNBinsX,fgkDefaultMaxX,fgkDefaultMinX,
              fgkDefaultNBinsY,fgkDefaultMaxY,fgkDefaultMinY,
              fgkDefaultNBinsZ,fgkDefaultMaxZ,fgkDefaultMinZ);
}
/******************************************************************/

void AliHBTFunction3D::BuildHistos(Int_t nxbins, Float_t xmax, Float_t xmin,
                                   Int_t nybins, Float_t ymax, Float_t ymin,
	               Int_t nzbins, Float_t zmax, Float_t zmin)
{
  //Builds numerator and denominator histograms (3d-case)
   TString numstr = fName + " Numerator";  //title and name of the 
                                           //numerator histogram
   TString denstr = fName + " Denominator";//title and name of the 
                                           //denominator histogram
         
   fNumerator   = new TH3D(numstr.Data(),numstr.Data(),
                           nxbins,xmin,xmin,nybins,ymin,ymin,nzbins,zmin,zmin);
	       
   fDenominator = new TH3D(denstr.Data(),denstr.Data(),
                           nxbins,xmin,xmin,nybins,ymin,ymin,nzbins,zmin,zmin);
   
   fNumerator->Sumw2();
   fDenominator->Sumw2();
}

/******************************************************************/

Double_t AliHBTFunction3D::Scale()
{
  // Calculates the factor that should be used to scale 
  // quatience of fNumerator and fDenominator to 1 at 
  // given volume
  if (gDebug>0) Info("Scale","Enetered Scale()");
  if(!fNumerator) 
   {
     Error("Scale","No numerator");
     return 0.0;
   }
  if(!fDenominator) 
   {
     Error("Scale","No denominator");
     return 0.0;
   }
  
  if( (fNBinsToScaleX < 1) || (fNBinsToScaleY < 1) || (fNBinsToScaleZ < 1)) 
   {
    return 0.0;
    Error("Scale","Number of bins for scaling is smaller thnan 1");
   }
  UInt_t nbinsX = fNumerator->GetNbinsX();
  if (fNBinsToScaleX > nbinsX) 
   {
    Error("Scale","Number of X bins for scaling is bigger thnan number of bins in histograms");
    return 0.0;
   }
   
  UInt_t nbinsY = fNumerator->GetNbinsX();
  if (fNBinsToScaleY > nbinsY) 
   {
    Error("Scale","Number of Y bins for scaling is bigger thnan number of bins in histograms");
    return 0.0;
   }

  UInt_t nbinsZ = fNumerator->GetNbinsZ();
  if (fNBinsToScaleZ > nbinsZ) 
   {
    Error("Scale","Number of Z bins for scaling is bigger thnan number of bins in histograms");
    return 0.0;
   }

  if (gDebug>0) Info("Scale","No errors detected");

  Int_t offsetX = nbinsX - fNBinsToScaleX - 1; //bin that we start loop over bins in axis X
  Int_t offsetY = nbinsY - fNBinsToScaleY - 1; //bin that we start loop over bins in axis Y
  Int_t offsetZ = nbinsZ - fNBinsToScaleZ - 1; //bin that we start loop over bins in axis Z

  Double_t ratio;
  Double_t sum = 0;
  Int_t n = 0;
  
  for (UInt_t k = offsetZ; k<nbinsZ; k++)
    for (UInt_t j = offsetY; j<nbinsY; j++)
      for (UInt_t i = offsetX; i<nbinsX; i++)
       {
        if ( fNumerator->GetBinContent(i,j,k) > 0.0 )
         {
           ratio = fDenominator->GetBinContent(i,j,k)/fNumerator->GetBinContent(i,j,k);
           sum += ratio;
           n++;
         }
       }
  
  if(gDebug > 0) 
    Info("Scale","sum=%f fNBinsToScaleX=%d fNBinsToScaleY=%d fNBinsToScaleZ=%d n=%d",
          sum,fNBinsToScaleX,fNBinsToScaleY,fNBinsToScaleZ,n);
  
  if (n == 0) return 0.0;
  Double_t ret = sum/((Double_t)n);

  if(gDebug > 0) Info("Scale","returning %f",ret);
  return ret;
} 
/******************************************************************/

void AliHBTFunction3D::SetNumberOfBinsToScale(UInt_t xn, UInt_t yn,UInt_t zn)
{
//sets up the volume to be used for scaling to tail
 fNBinsToScaleX = xn;
 fNBinsToScaleY = yn;
 fNBinsToScaleZ = zn;
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTOnePairFctn1D                               //
//                                                   //
// Base Calss for 1-dimensinal Functions that need   //
// one pair to fill function                         //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

ClassImp( AliHBTOnePairFctn1D )
/******************************************************************/

AliHBTOnePairFctn1D::AliHBTOnePairFctn1D(Int_t nbins, Float_t maxXval, Float_t minXval):
 AliHBTFunction1D(nbins,maxXval,minXval)
{
}
/******************************************************************/

AliHBTOnePairFctn1D::AliHBTOnePairFctn1D(const Char_t *name, const Char_t *title):
 AliHBTFunction1D(name,title)
{
}
/******************************************************************/

AliHBTOnePairFctn1D::AliHBTOnePairFctn1D(const Char_t *name, const Char_t *title,
                                         Int_t nbins, Float_t maxXval, Float_t minXval):
 AliHBTFunction1D(name,title,nbins,maxXval,minXval)
{
}
/******************************************************************/

void AliHBTOnePairFctn1D::ProcessSameEventParticles(AliHBTPair* pair)
{
 //Fills the numerator using pair from the same event
   pair = CheckPair(pair);
   if(pair) fNumerator->Fill(GetValue(pair));
}
/******************************************************************/
void AliHBTOnePairFctn1D::ProcessDiffEventParticles(AliHBTPair* pair)
 {
  //Fills the denumerator using mixed pairs
   pair = CheckPair(pair);
   if(pair) fDenominator->Fill(GetValue(pair));
  }

/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTOnePairFctn2D                               //
//                                                   //
// Base Calss for 2-dimensinal Functions that need   //
// one pair to fill function                         //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

ClassImp( AliHBTOnePairFctn2D )
/******************************************************************/

AliHBTOnePairFctn2D::AliHBTOnePairFctn2D(const Char_t *name, const Char_t *title):
 AliHBTFunction2D(name,title)
{
}
/******************************************************************/

AliHBTOnePairFctn2D::AliHBTOnePairFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTFunction2D(nXbins, maxXval, minXval, nYbins, maxYval, minYval)
{
}
/******************************************************************/

AliHBTOnePairFctn2D::AliHBTOnePairFctn2D(const Char_t *name, const Char_t *title,
                      Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTFunction2D(name,title, nXbins, maxXval, minXval, nYbins, maxYval, minYval)
{
}
/******************************************************************/

void AliHBTOnePairFctn2D::ProcessSameEventParticles(AliHBTPair* pair)
{
  // Fills the numerator using pairs from the same event
  pair = CheckPair(pair);
  if(pair) 
   { 
     Double_t x,y;
     GetValues(pair,x,y);
     fNumerator->Fill(x,y);
   }
}
/******************************************************************/

void AliHBTOnePairFctn2D::ProcessDiffEventParticles(AliHBTPair* pair)
{
  // Fills the denumerator using mixed pairs
  pair = CheckPair(pair);
  if(pair) 
   { 
     Double_t x,y;
     GetValues(pair,x,y);
     fDenominator->Fill(x,y);
   }
}
/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/
//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTOnePairFctn3D                               //
//                                                   //
// Base Calss for 3-dimensinal Functions that need   //
// one pair to fill function                         //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////
ClassImp( AliHBTOnePairFctn3D)

/******************************************************************/

AliHBTOnePairFctn3D::AliHBTOnePairFctn3D(const Char_t *name, const Char_t *title):
 AliHBTFunction3D(name,title)
{
}
/******************************************************************/

AliHBTOnePairFctn3D::AliHBTOnePairFctn3D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval, 
                      Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTFunction3D(nXbins, maxXval, minXval, nYbins, maxYval, minYval, nZbins, maxZval, minZval)
{
}
/******************************************************************/

AliHBTOnePairFctn3D::AliHBTOnePairFctn3D(const Char_t *name, const Char_t *title,
                      Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval, 
                      Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTFunction3D(name,title,nXbins, maxXval, minXval, nYbins, maxYval, minYval, nZbins, maxZval, minZval)
{
}
/******************************************************************/

void AliHBTOnePairFctn3D::ProcessSameEventParticles(AliHBTPair* pair)
{
//Reacts on pair coming from same event (real pairs)
//and fills numerator histogram
  pair  = CheckPair(pair);
  if( pair ) 
   { 
     Double_t x,y,z;
     GetValues(pair,x,y,z);
     fNumerator->Fill(x,y,z);
   }
}
/******************************************************************/

void AliHBTOnePairFctn3D::ProcessDiffEventParticles(AliHBTPair* pair)
{
//Reacts on pair coming from different events (mixed pairs)
//and fills denominator histogram
  pair  = CheckPair(pair);
  if( pair ) 
   { 
     Double_t x,y,z;
     GetValues(pair,x,y,z);
     fDenominator->Fill(x,y,z);
   }
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTTwoPairFctn1D                               //
//                                                   //
// Base Calss for 1-dimensinal Functions that need   //
// two pair (simulated and reconstructed)            //
// to fill function                                  //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////
ClassImp(AliHBTTwoPairFctn1D)
/******************************************************************/

AliHBTTwoPairFctn1D::AliHBTTwoPairFctn1D(Int_t nbins, Float_t maxXval, Float_t minXval):
 AliHBTFunction1D(nbins,maxXval,minXval)
{
}
/******************************************************************/

AliHBTTwoPairFctn1D::AliHBTTwoPairFctn1D(const Char_t *name, const Char_t *title):
 AliHBTFunction1D(name,title)
{
}
/******************************************************************/

AliHBTTwoPairFctn1D::AliHBTTwoPairFctn1D(const Char_t *name, const Char_t *title,
                                         Int_t nbins, Float_t maxXval, Float_t minXval):
 AliHBTFunction1D(name,title,nbins,maxXval,minXval)
{
}
/******************************************************************/

void AliHBTTwoPairFctn1D::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  // Fills the numerator using pairs from the same event
  partpair  = CheckPair(partpair);
  if( partpair ) 
   { 
     Double_t x = GetValue(trackpair,partpair);
     fNumerator->Fill(x);
   }
}
/******************************************************************/

void AliHBTTwoPairFctn1D::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  // Fills the denumerator usin mixed pairs
  partpair  = CheckPair(partpair);
  if( partpair )
   { 
     Double_t x = GetValue(trackpair,partpair);
     fDenominator->Fill(x);
   }
}
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTTwoPairFctn2D                               //
//                                                   //
// Base Calss for 2-dimensinal Functions that need   //
// two pair (simulated and reconstructed)            //
// to fill function                                  //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

ClassImp(AliHBTTwoPairFctn2D)

/******************************************************************/

AliHBTTwoPairFctn2D::AliHBTTwoPairFctn2D(const Char_t *name, const Char_t *title):
 AliHBTFunction2D(name,title)
{
}
/******************************************************************/

AliHBTTwoPairFctn2D::AliHBTTwoPairFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTFunction2D(nXbins, maxXval, minXval, nYbins, maxYval, minYval)
{
}
/******************************************************************/

AliHBTTwoPairFctn2D::AliHBTTwoPairFctn2D(const Char_t *name, const Char_t *title,
                      Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTFunction2D(name,title,nXbins, maxXval, minXval, nYbins, maxYval, minYval)
{
}
/******************************************************************/

void AliHBTTwoPairFctn2D::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//processes pair of particles coming from a same events (real pair)
  partpair  = CheckPair(partpair);  //check cuts
  if( partpair ) 
   { 
     Double_t x,y;
     GetValues(trackpair,partpair,x,y);
     fNumerator->Fill(x,y);
   }
}
/******************************************************************/

void AliHBTTwoPairFctn2D::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
//processes pair of particles coming from a different events (mixed pair)
  partpair  = CheckPair(partpair);
  if( partpair ) 
   { 
     Double_t x,y;
     GetValues(trackpair,partpair,x,y);
     fDenominator->Fill(x,y);
   }
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTTwoPairFctn3D                               //
//                                                   //
// Base Calss for 3-dimensinal Functions that need   //
// two pair (simulated and reconstructed)            //
// to fill function                                  //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

ClassImp(AliHBTTwoPairFctn3D)

/******************************************************************/

AliHBTTwoPairFctn3D::AliHBTTwoPairFctn3D(const Char_t *name, const Char_t *title):
 AliHBTFunction3D(name,title)
{
}
/******************************************************************/

AliHBTTwoPairFctn3D::AliHBTTwoPairFctn3D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval, 
                      Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTFunction3D(nXbins, maxXval, minXval, nYbins, maxYval, minYval, nZbins, maxZval, minZval)
{
}
/******************************************************************/

AliHBTTwoPairFctn3D::AliHBTTwoPairFctn3D(const Char_t *name, const Char_t *title,
                      Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval, 
                      Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTFunction3D(name,title,nXbins, maxXval, minXval, nYbins, maxYval, minYval, nZbins, maxZval, minZval)
{
}
/******************************************************************/

void AliHBTTwoPairFctn3D::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  // Fills th numerator using pairs from the same event
  partpair  = CheckPair(partpair);
  if( partpair ) 
   { 
     Double_t x,y,z;
     GetValues(trackpair,partpair,x,y,z);
     fNumerator->Fill(x,y,z);
   }
}
/******************************************************************/

void AliHBTTwoPairFctn3D::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  // Fills the denumerator using mixed pairs
  partpair  = CheckPair(partpair);
  if( partpair ) 
   { 
     Double_t x,y,z;
     GetValues(trackpair,partpair,x,y,z);
     fDenominator->Fill(x,y,z);
   }

}


/******************************************************************/
/******************************************************************/
/******************************************************************/

