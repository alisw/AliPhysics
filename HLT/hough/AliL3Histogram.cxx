//Author:        Anders Strand Vestbo
//Last Modified: 28.6.01

#include "AliL3Logging.h"
#include "AliL3Histogram.h"

//2D histogram class.

ClassImp(AliL3Histogram)

AliL3Histogram::AliL3Histogram()
{
  fNxbins = 0;
  fNybins = 0;
  fNcells = 0;
  fXmin = 0;
  fYmin = 0;
  fXmax = 0;
  fYmax = 0;
  fEntries = 0;
  fContent = 0;
}

  
AliL3Histogram::AliL3Histogram(Char_t *name,Char_t *id,Int_t nxbin,Double_t xmin,Double_t xmax,Int_t nybin,Double_t ymin,Double_t ymax) : TH2F(name,id,nxbin,xmin,xmax,nybin,ymin,ymax)
{
  
  strcpy(fName,name);
  fNxbins = nxbin;
  fNybins = nybin;
  fNcells = (nxbin+2)*(nybin+2);
  
  fXmin = xmin;
  fYmin = ymin;
  fXmax = xmax;
  fYmax = ymax;
  fEntries = 0;
  
  fContent = new Double_t[fNcells];
  Reset();
}

AliL3Histogram::~AliL3Histogram()
{
  //Destructor
  if(fContent)
    delete [] fContent;
}


void AliL3Histogram::Reset()
{
  
  for(Int_t i=0; i<fNcells; i++)
    fContent[i] = 0;
  fEntries=0;
}

void AliL3Histogram::Fill(Double_t x,Double_t y,Int_t weight)
{
  Int_t bin = FindBin(x,y);
  AddBinContent(bin,weight);

}

Int_t AliL3Histogram::FindBin(Double_t x,Double_t y)
{
  if(x < fXmin || x > fXmax)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::FindBin","array")<<AliL3Log::kDec<<
	"x-value out of range "<<x<<ENDLOG;
      return 0;
    }
  if(y < fYmin || y > fYmax)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::FindBin","array")<<AliL3Log::kDec<<
	"y-value out of range "<<y<<ENDLOG;
      return 0;
    }
  
  Int_t xbin = 1 + (Int_t)(fNxbins*(x-fXmin)/(fXmax-fXmin));
  Int_t ybin = 1 + (Int_t)(fNybins*(y-fYmin)/(fYmax-fYmin));

  return xbin + ybin*(fNxbins+2);
  
}

void AliL3Histogram::AddBinContent(Int_t xbin,Int_t ybin,Int_t weight)
{
  Int_t bin = xbin + ybin*(fNxbins+2);
  AddBinContent(bin,weight);

}

void AliL3Histogram::AddBinContent(Int_t bin,Int_t weight)
{
  if(bin < 0 || bin > fNcells)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::AddBinContent","array")<<AliL3Log::kDec<<
	"bin-value out of range "<<bin<<ENDLOG;
      return;
    }
  fEntries++;
  fContent[bin] += weight;
}

Double_t AliL3Histogram::GetXBinCenter(Int_t xbin)
{
  
  Double_t binwidth = (fXmax - fXmin) / fNxbins;
  return fXmin + (xbin-1) * binwidth + 0.5*binwidth;
  
}

Double_t AliL3Histogram::GetYBinCenter(Int_t ybin)
{
  
  Double_t binwidth = (fYmax - fYmin) / fNybins;
  return fYmin + (ybin-1) * binwidth + 0.5*binwidth;
  
}


void AliL3Histogram::Draw()
{
  TH2F *hist = new TH2F(fName,"",fNxbins,fXmin,fXmax,fNybins,fYmin,fYmax);
  for(Int_t bin=0; bin<fNcells; bin++)
    {
      hist->AddBinContent(bin,fContent[bin]);
      if(fContent[bin]!=0) printf("bin %d\n",bin);
    }
  //printf("ncells %d %d\n",(hist->GetNbinsX()+2)*(hist->GetNbinsY()+2),fNcells);
  printf("maxbin %d\n",hist->GetMaximumBin());
  
  hist->Draw();
}
