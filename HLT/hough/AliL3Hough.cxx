#include <string.h>
#include <TH2.h>

#include "AliL3Hough.h"
#include "AliL3HoughTransformer.h"
#include "AliL3HoughMaxFinder.h"

ClassImp(AliL3Hough)

AliL3Hough::AliL3Hough()
{



}


AliL3Hough::AliL3Hough(Char_t *rootfile,TH2F *hist)
{
  
  fParamSpace = hist;
  strcpy(fInputFile,rootfile);
  
}

AliL3Hough::AliL3Hough(Char_t *rootfile,Int_t xbin,Double_t *xrange,Int_t ybin,Double_t *yrange)
{

  fParamSpace = new TH2F("fParamSpace","Parameter space",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  strcpy(fInputFile,rootfile);
}


AliL3Hough::~AliL3Hough()
{
  
  if(fHoughTransformer)
    delete fHoughTransformer;
}



void AliL3Hough::ProcessSlice(Int_t slice)
{
  

}

void AliL3Hough::ProcessPatch(Int_t patch)
{
    

}

void AliL3Hough::ProcessEtaSlice(Int_t patch,Double_t *eta)
{
  
  fHoughTransformer = new AliL3HoughTransformer(2,patch,eta);
  fHoughTransformer->GetPixels(fInputFile);
  fParamSpace->Reset();
  fHoughTransformer->InitTemplates(fParamSpace);
  fHoughTransformer->Transform2Circle(fParamSpace,0);
  
}
