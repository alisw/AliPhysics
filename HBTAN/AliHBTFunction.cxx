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
   if (res) GetResult()->Write();
 }
/******************************************************************/

TH1* AliHBTFunction::
GetRatio(Double_t normfactor)
 {
   TString str = fName + " ratio";
   TH1 *result = (TH1*)GetNumerator()->Clone(str.Data());
   
   result->SetTitle(str.Data());
   result->Sumw2();
   
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
   
 }
/******************************************************************/
AliHBTTwoPartFctn1D::~AliHBTTwoPartFctn1D()
{
  delete fNumerator;
  delete fDenominator;
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

