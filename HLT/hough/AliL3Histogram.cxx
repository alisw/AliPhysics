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
  fFirstXbin = 0;
  fLastXbin = 0;
  fFirstYbin = 0;
  fLastYbin = 0;
  fEntries = 0;
  fContent = 0;
}

  
AliL3Histogram::AliL3Histogram(Char_t *name,Char_t *id,Int_t nxbin,Double_t xmin,Double_t xmax,Int_t nybin,Double_t ymin,Double_t ymax) 
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
  fFirstXbin = 1;
  fFirstYbin = 1;
  fLastXbin = nxbin;
  fLastYbin = nybin;
  fRootHisto = 0;

  fContent = new Double_t[fNcells];
  Reset();
}

AliL3Histogram::~AliL3Histogram()
{
  //Destructor
  if(fContent)
    delete [] fContent;
  if(fRootHisto)
    delete fRootHisto;
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
  
  Int_t xbin = FindXbin(x);
  Int_t ybin = FindYbin(y);
  
  return GetBin(xbin,ybin);
}

Int_t AliL3Histogram::FindXbin(Double_t x)
{
  if(x < fXmin || x > fXmax)
    return 0;
  
  return 1 + (Int_t)(fNxbins*(x-fXmin)/(fXmax-fXmin));

}

Int_t AliL3Histogram::FindYbin(Double_t y)
{
  if(y < fYmin || y > fYmax)
    return 0;
  
  return 1 + (Int_t)(fNybins*(y-fYmin)/(fYmax-fYmin));

}

Int_t AliL3Histogram::GetBin(Int_t xbin,Int_t ybin)
{
  if(xbin < 0 || xbin > GetLastXbin())
    {
      LOG(AliL3Log::kError,"AliL3Histogram::GetBin","array")<<AliL3Log::kDec<<
	"xbin out of range "<<xbin<<ENDLOG;
      return 0;
    }
  if(ybin < 0 || ybin > GetLastYbin())
    {
      LOG(AliL3Log::kError,"AliL3Histogram::FindYbin","array")<<AliL3Log::kDec<<
	"ybin out of range "<<xbin<<ENDLOG;
      return 0;
    }
    
  return xbin + ybin*(fNxbins+2);
}

Double_t AliL3Histogram::GetBinContent(Int_t bin)
{
  if(bin >= fNcells)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::GetBinContent","array")<<AliL3Log::kDec<<
	"bin out of range "<<bin<<ENDLOG;
      return 0;
    }
  
  return fContent[bin];
}

void AliL3Histogram::SetBinContent(Int_t xbin,Int_t ybin,Int_t value)
{
  Int_t bin = GetBin(xbin,ybin);
  if(bin == 0) 
    return;
  SetBinContent(bin,value);
}

void AliL3Histogram::SetBinContent(Int_t bin,Int_t value)
{

  if(bin >= fNcells)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::SetBinContent","array")<<AliL3Log::kDec<<
	"bin out of range "<<bin<<ENDLOG;
      return;
    }
  if(bin == 0)
    return;
  fContent[bin]=value;
  
}

void AliL3Histogram::AddBinContent(Int_t xbin,Int_t ybin,Int_t weight)
{
  Int_t bin = GetBin(xbin,ybin);
  if(bin == 0)
    return;
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
  if(bin == 0)
    return;
  fEntries++;
  fContent[bin] += weight;
}

Double_t AliL3Histogram::GetBinCenterX(Int_t xbin)
{
  
  Double_t binwidth = (fXmax - fXmin) / fNxbins;
  return fXmin + (xbin-1) * binwidth + 0.5*binwidth;
  
}

Double_t AliL3Histogram::GetBinCenterY(Int_t ybin)
{
  
  Double_t binwidth = (fYmax - fYmin) / fNybins;
  return fYmin + (ybin-1) * binwidth + 0.5*binwidth;
  
}


void AliL3Histogram::Draw(Char_t *option)
{
  fRootHisto = new TH2F(fName,"",fNxbins,fXmin,fXmax,fNybins,fYmin,fYmax);
  for(Int_t bin=0; bin<fNcells; bin++)
    {
      fRootHisto->AddBinContent(bin,fContent[bin]);
    }
  //printf("ncells %d %d\n",(hist->GetNbinsX()+2)*(hist->GetNbinsY()+2),fNcells);
  
  fRootHisto->Draw(option);
  
}
