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
    one pair          two pair 
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
/******************************************************************/
AliHBTFunction::AliHBTFunction(const char* name,const char* title):TNamed(name,title)
{
  fPairCut = new AliHBTEmptyPairCut(); //dummy cut
}
/******************************************************************/

AliHBTFunction::~AliHBTFunction()
 {
  delete fPairCut;
 }
/******************************************************************/

void AliHBTFunction::Write()
 {
   if (GetNumerator()) GetNumerator()->Write();
   if (GetDenominator()) GetDenominator()->Write();
   TH1* res = GetResult();
   if (res) res->Write();
 }
/******************************************************************/

TH1* AliHBTFunction::GetRatio(Double_t normfactor)
 {
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

ClassImp( AliHBTOnePairFctn1D )
AliHBTOnePairFctn1D::AliHBTOnePairFctn1D():fNBinsToScale(30)
 {
   fNumerator = 0x0;
   fDenominator = 0x0;
 }

AliHBTOnePairFctn1D::
AliHBTOnePairFctn1D(Int_t nbins, Double_t maxXval, Double_t minXval)
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
AliHBTOnePairFctn1D::
AliHBTOnePairFctn1D(const Char_t *name, const Char_t *title,
                    Int_t nbins, Double_t maxXval, Double_t minXval)
	:AliHBTOnePairFctn(name,title)
{
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
AliHBTOnePairFctn1D::~AliHBTOnePairFctn1D()
{
  delete fNumerator;
  delete fDenominator;
}
/******************************************************************/

void AliHBTOnePairFctn1D::ProcessSameEventParticles(AliHBTPair* pair)
{
 //Fills the numerator
   pair = CheckPair(pair);
   if(pair) fNumerator->Fill(GetValue(pair));
}
/******************************************************************/
void AliHBTOnePairFctn1D::ProcessDiffEventParticles(AliHBTPair* pair)
 {
  //fills denumerator
   pair = CheckPair(pair);
   if(pair) fDenominator->Fill(GetValue(pair));

  }
/******************************************************************/
Double_t AliHBTOnePairFctn1D::Scale()
{
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
  if (gDebug>0) Info("Scale","No errors detected");

  Double_t ratio;
  Double_t sum = 0;
  Int_t N = 0;
  
  Int_t offset = nbins - fNBinsToScale - 1; 
  Int_t i;
  for ( i = offset; i< nbins; i++)
   {
    if ( fNumerator->GetBinContent(i) > 0.0 )
     {
       ratio = fDenominator->GetBinContent(i)/fNumerator->GetBinContent(i);
       sum += ratio;
       N++;
     }
   }
  
  if(gDebug > 0) Info("Scale","sum=%f fNBinsToScale=%d N=%d",sum,fNBinsToScale,N);
  
  if (N == 0) return 0.0;
  Double_t ret = sum/((Double_t)N);

  if(gDebug > 0) Info("Scale","returning %f",ret);
  return ret;
} 

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTOnePairFctn2D )

AliHBTOnePairFctn2D::
AliHBTOnePairFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval , 
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
AliHBTOnePairFctn2D::~AliHBTOnePairFctn2D()
{
  delete fNumerator;
  delete fDenominator;
}
void AliHBTOnePairFctn2D::ProcessSameEventParticles(AliHBTPair* pair)
{
  pair = CheckPair(pair);
  if(pair) 
   { 
     Double_t x,y;
     GetValues(pair,x,y);
     fNumerator->Fill(x,y);
   }
}

void AliHBTOnePairFctn2D::ProcessDiffEventParticles(AliHBTPair* pair)
{
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

ClassImp( AliHBTOnePairFctn3D)

AliHBTOnePairFctn3D::
AliHBTOnePairFctn3D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
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
/******************************************************************/

AliHBTOnePairFctn3D::~AliHBTOnePairFctn3D()
{
  delete fNumerator;
  delete fDenominator;
}
/******************************************************************/


/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTTwoPairFctn1D)

AliHBTTwoPairFctn1D::
AliHBTTwoPairFctn1D(Int_t nbins, Double_t maxval, Double_t minval)
 {
   TString numstr = fName + " Numerator";  //title and name of the 
                                           //numerator histogram
   TString denstr = fName + " Denominator";//title and name of the 
                                           //denominator histogram
         
   fNumerator   = new TH1D(numstr.Data(),numstr.Data(),
                           nbins,minval,maxval);
	       
   fDenominator = new TH1D(denstr.Data(),denstr.Data(),
                           nbins,minval,maxval);
   
   fNumerator->Sumw2();
   fDenominator->Sumw2();
   fNBinsToScale = 30;
 }

AliHBTTwoPairFctn1D::
AliHBTTwoPairFctn1D(const Char_t* name, const Char_t* title,
                    Int_t nbins, Double_t maxval, Double_t minval)
	:AliHBTTwoPairFctn(name,title)
 {
   TString numstr = fName + " Numerator";  //title and name of the 
                                           //numerator histogram
   TString denstr = fName + " Denominator";//title and name of the 
                                           //denominator histogram
         
   fNumerator   = new TH1D(numstr.Data(),numstr.Data(),
                           nbins,minval,maxval);
	       
   fDenominator = new TH1D(denstr.Data(),denstr.Data(),
                           nbins,minval,maxval);
   
   fNumerator->Sumw2();
   fDenominator->Sumw2();
   fNBinsToScale = 30;
 }


/******************************************************************/
AliHBTTwoPairFctn1D::~AliHBTTwoPairFctn1D()
{
  delete fNumerator;
  delete fDenominator;
}
void AliHBTTwoPairFctn1D::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
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
  partpair  = CheckPair(partpair);
  if( partpair )
   { 
     Double_t x = GetValue(trackpair,partpair);
     fDenominator->Fill(x);
   }
}
/******************************************************************/
Double_t AliHBTTwoPairFctn1D::Scale()
{
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
  if (gDebug>0) Info("Scale","No errors detected");

  Double_t ratio;
  Double_t sum = 0;
  Int_t N = 0;
  
  Int_t offset = nbins - fNBinsToScale - 1; 
  Int_t i;
  for ( i = offset; i< nbins; i++)
   {
    if ( fNumerator->GetBinContent(i) > 0.0 )
     {
       ratio = fDenominator->GetBinContent(i)/fNumerator->GetBinContent(i);
       sum += ratio;
       N++;
     }
   }
  
  if(gDebug > 0) Info("Scale","sum=%f fNBinsToScale=%d N=%d",sum,fNBinsToScale,N);
  
  if (N == 0) return 0.0;
  Double_t ret = sum/((Double_t)N);

  if(gDebug > 0) Info("Scale","returning %f",ret);
  return ret;
} 

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTTwoPairFctn2D)


AliHBTTwoPairFctn2D::
AliHBTTwoPairFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval , 
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
AliHBTTwoPairFctn2D::~AliHBTTwoPairFctn2D()
{
  delete fNumerator;
  delete fDenominator;
}
void AliHBTTwoPairFctn2D::
ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  partpair  = CheckPair(partpair);
  if( partpair ) 
   { 
     Double_t x,y;
     GetValues(trackpair,partpair,x,y);
     fNumerator->Fill(x,y);
   }
}

void AliHBTTwoPairFctn2D::
ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
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
ClassImp(AliHBTTwoPairFctn3D)

void AliHBTTwoPairFctn3D::
ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  partpair  = CheckPair(partpair);
  if( partpair ) 
   { 
     Double_t x,y,z;
     GetValues(trackpair,partpair,x,y,z);
     fNumerator->Fill(x,y,z);
   }
}

void AliHBTTwoPairFctn3D::
ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  partpair  = CheckPair(partpair);
  if( partpair ) 
   { 
     Double_t x,y,z;
     GetValues(trackpair,partpair,x,y,z);
     fDenominator->Fill(x,y,z);
   }

}

