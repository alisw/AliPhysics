#include "AliHBTFunction.h"
/******************************************************************/
/*
Piotr Krzysztof Skowronski
Piotr.Skowronski@cern.ch
Base classes for HBT functions

             function
              /    \
             /      \
            /        \
           /          \
          /            \
         /              \
        /                \
    two part          four part 
    /   |   \         /   |   \
   /    |    \       /    |    \
  1D   2D    3D     1D   2D    3D

 four particle functions are intendent to be resolution functions:
 it is mecessary to have simulated particle pair corresponding to given
 recontructed track pair in order to calculate function simualted value 
 and recontructed value to be further histogrammed
 
*/
/******************************************************************/
/******************************************************************/

#include <iostream.h>
ClassImp( AliHBTFunction )

AliHBTFunction::AliHBTFunction()
{
 
 fPairCut = new AliHBTEmptyPairCut(); //dummy cut
}

void AliHBTFunction::
Write()
 {
   if (GetNumerator()) GetNumerator()->Write();
   if (GetDenominator()) GetDenominator()->Write();
   TH1* res = GetResult();
   if (res) res->Write();
 }
/******************************************************************/

TH1* AliHBTFunction::
GetRatio(Double_t normfactor)
 {
   //if (gDebug>0) 
   cout<<"Mormfactor is "<<normfactor<<" for "<<fName<<endl;
   
   if (normfactor == 0.0)
    {
      Error("GetRatio","Scaling Factor is 0. Null poiner returned");
      return 0x0;
    }
   TString str = fName + " ratio";
   TH1 *result = (TH1*)GetNumerator()->Clone(str.Data());
   
   result->SetTitle(str.Data());
   //result->Sumw2();
   
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

void AliHBTFunction::
Rename(const Char_t * name)
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

void AliHBTFunction::
Rename(const Char_t * name, const Char_t * title)
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
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTTwoPartFctn )

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTFourPartFctn)

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTTwoPartFctn1D )

AliHBTTwoPartFctn1D::
AliHBTTwoPartFctn1D(Int_t nbins, Double_t maxXval, Double_t minXval)
 {
 //Constructor of Two Part One Dimentional Function 
 // nbins: number of bins in histograms - default 100
 // maxXval and minXval: range of histgram(s) default 0 - 0.15 (GeV)
 
 
   TString numstr = fName + " Numerator";  //title and name of the 
                                           //numerator histogram
   TString denstr = fName + " Denominator";//title and name of the 
                                           //denominator histogram
         
   fNumerator   = new TH1D(numstr.Data(),numstr.Data(),nbins,minXval,maxXval);
   fDenominator = new TH1D(denstr.Data(),denstr.Data(),nbins,minXval,maxXval);
   
   fNumerator->Sumw2();
   fDenominator->Sumw2();
   
   fNBinsToScale = 30;
   
 }
/******************************************************************/
AliHBTTwoPartFctn1D::~AliHBTTwoPartFctn1D()
{
  delete fNumerator;
  delete fDenominator;
}
/******************************************************************/

void AliHBTTwoPartFctn1D::ProcessSameEventParticles(AliHBTPair* pair)
{
 //Fills the numerator
   pair = CheckPair(pair);
   if(pair) fNumerator->Fill(GetValue(pair));
}
/******************************************************************/
void AliHBTTwoPartFctn1D::ProcessDiffEventParticles(AliHBTPair* pair)
 {
  //fills denumerator
   pair = CheckPair(pair);
   if(pair) fDenominator->Fill(GetValue(pair));

  }
/******************************************************************/
Double_t AliHBTTwoPartFctn1D::Scale()
{
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
  
  if(fNBinsToScale < 1) 
   {
    return 0.0;
    Error("Scale","Number of bins for scaling is smaller thnan 1");
   }
  Int_t nbins = fNumerator->GetNbinsX();
  if (fNBinsToScale > nbins) 
   {
    Error("Scale","Number of bins for scaling is bigger thnan number of bins in histograms");
    return 0.0;
   }
  Double_t ratios[fNBinsToScale];

  Int_t offset = nbins - fNBinsToScale - 1; 
  Int_t i;
  for ( i = offset; i< nbins; i++)
   {
    if ( fNumerator->GetBinContent(i) == 0.0 )
     {
       ratios[i - offset] = -1.0; //since we play with histograms negative is impossible 
                                  //so it is good flag
     }
    else
     {
       ratios[i - offset] = fDenominator->GetBinContent(i)/fNumerator->GetBinContent(i);
     }
   }
 
  Double_t sum = 0;
  Int_t skipped = 0;
  for (i = 0; i<fNBinsToScale; i++)
   {
    if (ratios[i] == -1.0) skipped++;
    else sum += ratios[i];
   }
  cout<<"sum="<<sum<<" fNBinsToScale="<<fNBinsToScale<<" skipped="<<skipped<<endl;

  return sum/(Double_t)(fNBinsToScale - skipped);
} 

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTTwoPartFctn2D )

AliHBTTwoPartFctn2D::
AliHBTTwoPartFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval , 
                    Int_t nYbins, Double_t maxYval, Double_t minYval)

{
   TString numstr = fName + " Numerator";  //title and name of the 
                                           //numerator histogram
   TString denstr = fName + " Denominator";//title and name of the 
                                           //denominator histogram
         
   fNumerator   = new TH2D(numstr.Data(),numstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval);
	       
   fDenominator = new TH2D(denstr.Data(),denstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval);
   
   fNumerator->Sumw2();
   fDenominator->Sumw2();

}	  
AliHBTTwoPartFctn2D::~AliHBTTwoPartFctn2D()
{
  delete fNumerator;
  delete fDenominator;
}
void AliHBTTwoPartFctn2D::ProcessSameEventParticles(AliHBTPair* pair)
{
  pair = CheckPair(pair);
  if(pair) 
   { 
     Double_t x,y;
     GetValues(pair,x,y);
     fNumerator->Fill(y,x);
   }
}

void AliHBTTwoPartFctn2D::ProcessDiffEventParticles(AliHBTPair* pair)
{
  pair = CheckPair(pair);
  if(pair) 
   { 
     Double_t x,y;
     GetValues(pair,x,y);
     fDenominator->Fill(y,x);
   }

}


/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTTwoPartFctn3D)

AliHBTTwoPartFctn3D::
AliHBTTwoPartFctn3D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                    Int_t nYbins, Double_t maxYval, Double_t minYval, 
                    Int_t nZbins, Double_t maxZval, Double_t minZval)

{
   TString numstr = fName + " Numerator";  //title and name of the 
                                           //numerator histogram
   TString denstr = fName + " Denominator";//title and name of the 
                                           //denominator histogram
         
   fNumerator   = new TH3D(numstr.Data(),numstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval,
	       nZbins,minZval,maxZval);
	       
   fDenominator = new TH3D(denstr.Data(),denstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval,
	       nZbins,minZval,maxZval);
   
   fNumerator->Sumw2();
   fDenominator->Sumw2();

}	  


AliHBTTwoPartFctn3D::~AliHBTTwoPartFctn3D()
{
  delete fNumerator;
  delete fDenominator;
}


/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTFourPartFctn2D)


AliHBTFourPartFctn2D::
AliHBTFourPartFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval , 
                    Int_t nYbins, Double_t maxYval, Double_t minYval)

{
   TString numstr = fName + " Numerator";  //title and name of the 
                                           //numerator histogram
   TString denstr = fName + " Denominator";//title and name of the 
                                           //denominator histogram
         
   fNumerator   = new TH2D(numstr.Data(),numstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval);
	       
   fDenominator = new TH2D(denstr.Data(),denstr.Data(),
                           nXbins,minXval,maxXval,
	       nYbins,minYval,maxYval);
   
   fNumerator->Sumw2();
   fDenominator->Sumw2();

}	  
AliHBTFourPartFctn2D::~AliHBTFourPartFctn2D()
{
  delete fNumerator;
  delete fDenominator;
}
void AliHBTFourPartFctn2D::
ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  partpair  = CheckPair(partpair);
  trackpair = CheckPair(trackpair);
  if( partpair && trackpair) 
   { 
     Double_t x,y;
     GetValues(trackpair,partpair,x,y);
     fNumerator->Fill(y,x);
   }
}

void AliHBTFourPartFctn2D::
ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  partpair  = CheckPair(partpair);
  trackpair = CheckPair(trackpair);
  if( partpair && trackpair)
   { 
     Double_t x,y;
     GetValues(trackpair,partpair,x,y);
     fDenominator->Fill(y,x);
   }

}

