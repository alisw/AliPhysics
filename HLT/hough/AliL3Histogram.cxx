// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3Histogram.h"

#if __GNUC__ == 3
using namespace std;
#endif

/** \class AliL3Histogram
<pre>
//_____________________________________________________________
// AliL3Histogram
//
// 2D histogram class
//
</pre>
*/

//uncomment if you want overflow checks
//#define _IFON_

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
  fBinwidthX = 0;
  fBinwidthY = 0;  
  fFirstXbin = 0;
  fLastXbin = 0;
  fFirstYbin = 0;
  fLastYbin = 0;
  fEntries = 0;
  fContent = 0;
  fThreshold = 0;
#ifdef use_root
  fRootHisto = 0;
#endif
}

AliL3Histogram::AliL3Histogram(Char_t *name,Char_t */*id*/,
			       Int_t nxbin,Double_t xmin,Double_t xmax,
			       Int_t nybin,Double_t ymin,Double_t ymax) 
{
  strcpy(fName,name);

  fNxbins = nxbin;
  fNybins = nybin;
  fNcells = (nxbin+2)*(nybin+2);
  fXmin = xmin;
  fYmin = ymin;
  fXmax = xmax;
  fYmax = ymax;
  fBinwidthX = (fXmax - fXmin) / fNxbins;
  fBinwidthY = (fYmax - fYmin) / fNybins;
  
  fEntries = 0;
  fFirstXbin = 1;
  fFirstYbin = 1;
  fLastXbin = nxbin;
  fLastYbin = nybin;
#ifdef use_root
  fRootHisto = 0;
#endif
  fThreshold = 0;

  fContent = new Int_t[fNcells];
  Reset();
}

AliL3Histogram::~AliL3Histogram()
{
  //Destructor
  if(fContent)
    delete [] fContent;
#ifdef use_root
  if(fRootHisto)
    delete fRootHisto;
#endif
}

void AliL3Histogram::Reset()
{
  if(fContent)
    for(Int_t i=0; i<fNcells; i++) fContent[i] = 0;

  fEntries=0;
}

void AliL3Histogram::Fill(Double_t x,Double_t y,Int_t weight)
{
  Int_t bin = FindBin(x,y);
#ifdef _IFON_
  if(bin < 0)
    return;
#endif
  
  AddBinContent(bin,weight);
}

void AliL3Histogram::Fill(Double_t x,Int_t ybin,Int_t weight)
{
  Int_t xbin = FindXbin(x);
  Int_t bin = GetBin(xbin,ybin);
#ifdef _IFON_
  if(bin < 0)
    return;
#endif
  
  AddBinContent(bin,weight);
}

void AliL3Histogram::Fill(Int_t xbin,Double_t y,Int_t weight)
{
  Int_t ybin = FindYbin(y);
  Int_t bin = GetBin(xbin,ybin);
#ifdef _IFON_
  if(bin < 0)
    return;
#endif
  
  AddBinContent(bin,weight);
}

void AliL3Histogram::Fill(Int_t xbin,Int_t ybin,Int_t weight)
{
  Int_t bin = GetBin(xbin,ybin);
#ifdef _IFON_
  if(bin < 0)
    return;
#endif
  
  AddBinContent(bin,weight);
}

Int_t AliL3Histogram::FindBin(Double_t x,Double_t y) const
{
  Int_t xbin = FindXbin(x);
  Int_t ybin = FindYbin(y);
#ifdef _IFON_
  if(!xbin || !ybin)
    return -1;
#endif
  
  return GetBin(xbin,ybin);
}

Int_t AliL3Histogram::FindXbin(Double_t x) const
{
  if(x < fXmin || x > fXmax)
    return 0;
  
  return 1 + (Int_t)(fNxbins*(x-fXmin)/(fXmax-fXmin));
}

Int_t AliL3Histogram::FindYbin(Double_t y) const
{
  if(y < fYmin || y > fYmax)
    return 0;
  
  return 1 + (Int_t)(fNybins*(y-fYmin)/(fYmax-fYmin));
}

Int_t AliL3Histogram::GetBin(Int_t xbin,Int_t ybin) const
{
  if(xbin < fFirstXbin || xbin > fLastXbin)
    return 0;
  if(ybin < fFirstYbin || ybin > fLastYbin)
    return 0;
    
  return xbin + ybin*(fNxbins+2);
}

Int_t AliL3Histogram::GetBinContent(Int_t bin) const
{
  if(bin >= fNcells)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::GetBinContent","array")<<AliL3Log::kDec<<
	"bin out of range "<<bin<<ENDLOG;
      return 0;
    }
  
  if(fContent[bin] < fThreshold)
    return 0;
  return fContent[bin];
}

void AliL3Histogram::SetBinContent(Int_t xbin,Int_t ybin,Int_t value)
{
  Int_t bin = GetBin(xbin,ybin);
#ifdef _IFON_
  if(bin == 0) 
    return;
#endif

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
#ifdef _IFON_
  if(bin == 0)
    return;
#endif

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

void AliL3Histogram::Add(AliL3Histogram *h1,Double_t /*weight*/)
{
  //Adding two histograms. Should be identical.
  
  if(!h1)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::Add","Pointer")<<
	"Attempting to add a non-existing histogram"<<ENDLOG;
      return;
    }
  
  if(h1->GetNbinsX()!=fNxbins || h1->GetNbinsY()!=fNybins)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::Add","array")<<
	"Mismatch in the number of bins "<<ENDLOG;
      return;
    }

  if(h1->GetFirstXbin()!=fFirstXbin || h1->GetLastXbin()!=fLastXbin ||
     h1->GetFirstYbin()!=fFirstYbin || h1->GetLastYbin()!=fLastYbin)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::Add","array")<<
	"Mismatch in the bin numbering "<<ENDLOG;
      return;
    }
  
  for(Int_t bin=0; bin<fNcells; bin++)
    fContent[bin] += h1->GetBinContent(bin);
  
  fEntries += h1->GetNEntries();
}

Double_t AliL3Histogram::GetBinCenterX(Int_t xbin) const
{
  if(xbin < fFirstXbin || xbin > fLastXbin)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::GetBinCenterX","xbin")
	<<"Bin-value out of range "<<xbin<<ENDLOG;
      return -1;
    }

  return fXmin + (xbin-0.5) * fBinwidthX;
}

Double_t AliL3Histogram::GetBinCenterY(Int_t ybin) const
{
  if(ybin < fFirstYbin || ybin > fLastYbin)
    {
      LOG(AliL3Log::kError,"AliL3Histogram::GetBinCenterY","ybin")
	<<"Bin-value out of range "<<ybin<<ENDLOG;
      return -1;
    }

  return fYmin + (ybin-0.5) * fBinwidthY;
}

Double_t AliL3Histogram::GetPreciseBinCenterX(Float_t xbin) const
{
  if(xbin < (fFirstXbin-0.5) || xbin > (fLastXbin+0.5))
    {
      LOG(AliL3Log::kError,"AliL3Histogram::GetBinCenterX","xbin")
	<<"Bin-value out of range "<<xbin<<ENDLOG;
      return -1;
    }
  //  return fXmin + (xbin-1) * fBinwidthX + 0.5*fBinwidthX;
  return fXmin + (xbin-0.5) * fBinwidthX;
}

Double_t AliL3Histogram::GetPreciseBinCenterY(Float_t ybin) const
{
  if(ybin < (fFirstYbin-0.5) || ybin > (fLastYbin+0.5))
    {
      LOG(AliL3Log::kError,"AliL3Histogram::GetBinCenterY","ybin")
	<<"Bin-value out of range "<<ybin<<ENDLOG;
      return -1;
    }
  //  return fYmin + (ybin-1) * fBinwidthY + 0.5*fBinwidthY;
  return fYmin + (ybin-0.5) * fBinwidthY;
}

void AliL3Histogram::Draw(Char_t *option)
{
#ifdef use_root
  if(!fRootHisto)
    CreateRootHisto();
  
  for(Int_t xbin=GetFirstXbin(); xbin<=GetLastXbin(); xbin++)
    {
      for(Int_t ybin=GetFirstYbin(); ybin<=GetLastYbin(); ybin++)
	{
	  Int_t bin = GetBin(xbin,ybin);
	  fRootHisto->Fill(GetBinCenterX(xbin),GetBinCenterY(ybin),GetBinContent(bin));
	}
    }
  
  //fRootHisto->SetStats(kFALSE);
  fRootHisto->Draw(option);
  return;
#endif
  cerr<<"AliL3Histogram::Draw : You need to compile with ROOT in order to draw histogram"<<endl;
  
}

void AliL3Histogram::CreateRootHisto()
{
#ifdef use_root
  fRootHisto = new TH2F(fName,"",fNxbins,fXmin,fXmax,fNybins,fYmin,fYmax);
  return;
#endif
  cerr<<"AliL3Histogram::CreateRootHisto : You need to compile with ROOT in order to create ROOT histogram"<<endl;
}

ofstream& operator<<(ofstream &o, const AliL3Histogram &h)
{
  for(Int_t xbin=h.GetFirstXbin(); xbin<=h.GetLastXbin(); xbin++)
    {
      for(Int_t ybin=h.GetFirstYbin(); ybin<=h.GetLastYbin(); ybin++)
	{
	  Int_t bin = h.GetBin(xbin,ybin);
	  o << h.GetBinContent(bin) << " ";
	}
      o << endl;
    }
  return o;
}
