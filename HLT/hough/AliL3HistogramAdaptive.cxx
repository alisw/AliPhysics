// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include <string.h>
#include "AliL3StandardIncludes.h"
#include "AliL3Logging.h"
#include "AliL3HistogramAdaptive.h"
#include "AliL3Transform.h"
#include "AliL3Track.h"

//_____________________________________________________________
// AliL3HistogramAdaptive
//
// 2D histogram class

ClassImp(AliL3HistogramAdaptive)

AliL3HistogramAdaptive::AliL3HistogramAdaptive() : AliL3Histogram()
{
  
}

  
AliL3HistogramAdaptive::AliL3HistogramAdaptive(Char_t *name,Double_t minpt,Double_t maxpt,Double_t ptres,
					       Int_t nybins,Double_t ymin,Double_t ymax)
{
  strcpy(fName,name);
  
  fPtres = ptres;
  fXmin = -1*AliL3Transform::GetBFact()*AliL3Transform::GetBField()/minpt;
  fXmax = AliL3Transform::GetBFact()*AliL3Transform::GetBField()/minpt;

  fMinPt = minpt;
  fMaxPt = maxpt;
  fNxbins = InitPtBins();
  //cout<<"Setting "<<fNxbins<<" bins on x"<<endl;
  
  fNybins = nybins;
  fYmin = ymin;
  fYmax = ymax;
  fFirstXbin=1;
  fFirstYbin=1;
  fLastXbin = fNxbins;
  fLastYbin = fNybins;
  fNcells = (fNxbins+2)*(fNybins+2);
  
  fThreshold=0;
  fContent = new Int_t[fNcells];
  Reset();
}

AliL3HistogramAdaptive::~AliL3HistogramAdaptive()
{
  
}

Int_t AliL3HistogramAdaptive::InitPtBins()
{
  
  Double_t pt = fMinPt,delta_pt,local_pt;
  Int_t bin=0;

  while(pt < fMaxPt)
    {
      local_pt = pt;
            
      delta_pt = fPtres*local_pt;
      pt += delta_pt;
      bin++;
      //      cout<<"Setting "<<bin<<" at step "<<local_pt<<" "<<pt<<" interval "<<delta_pt<<endl;
    }

  return (bin+1)*2; //Both negative and positive kappa.
}


void AliL3HistogramAdaptive::Fill(Double_t x,Double_t y,Int_t weight)
{
  Int_t bin = FindBin(x,y);
  if(bin < 0)
    return;
  AddBinContent(bin,weight);

}

Int_t AliL3HistogramAdaptive::FindBin(Double_t x,Double_t y)
{
  
  Int_t xbin = FindXbin(x);
  Int_t ybin = FindYbin(y);
  
  if(xbin < 0) 
    return -1;
  return GetBin(xbin,ybin);
}

Int_t AliL3HistogramAdaptive::FindXbin(Double_t x)
{
  
  Double_t ptfind = fabs(AliL3Transform::GetBFact()*AliL3Transform::GetBField()/x);
  if(ptfind < fMinPt || ptfind > fMaxPt) return -1;
  //  cout<<"Looking for pt "<<ptfind<<endl;
  Double_t pt = fMinPt;
  Double_t delta_pt,local_pt;
  Int_t bin=0;
  while(pt < fMaxPt)
    {
      local_pt = pt;
      delta_pt = fPtres*local_pt;
      pt += delta_pt;
      
      if(ptfind >= local_pt && ptfind < pt)
	{
	  //	  cout<<"Found in range "<<local_pt<<" "<<pt<<endl;
	  break;
	}
      bin++;
    }
  if(bin >= fNxbins/2)
    cerr<<"AliL3HistogramAdaptive::FindXbin : Bin out of range : "<<bin<<endl;
  
  //  cout<<"Found xbin "<<bin<<" and x is "<<x<<endl;
  if(x < 0)
    {
      //        cout<<"returning xbin "<<bin<<endl;
      return bin;
    }
  else
    return fNxbins - 1 - bin;
}

Int_t AliL3HistogramAdaptive::FindYbin(Double_t y)
{
  if(y < fYmin || y > fYmax)
    return 0;
  
  return 1 + (Int_t)(fNybins*(y-fYmin)/(fYmax-fYmin));

}

Double_t AliL3HistogramAdaptive::GetBinCenterX(Int_t xbin)
{
  //    cout<<"Looking for bin "<<xbin<<endl;
  if(xbin < 0 || xbin > fNxbins)
    {
      cerr<<"AliL3HistogramAdaptive::GetBinCenterX : Xbin out of range "<<xbin<<endl;
      return 0;
    }
  Double_t pt = fMinPt;
  Double_t delta_pt=0,local_pt=0;
  Int_t bin=0;
  while(pt < fMaxPt)
    {
      local_pt = pt;
      delta_pt = fPtres*local_pt;
      pt += delta_pt;
      if(xbin == bin || xbin == fNxbins - 1 - bin)
	break;
      bin++;
    }
  //    cout<<"get center at ptinterval "<<local_pt<<" "<<pt;
  
  Double_t kappa = AliL3Transform::GetBFact()*AliL3Transform::GetBField()/(local_pt + 0.5*delta_pt);
  //    cout<<" found pt "<<local_pt+delta_pt*0.5<<" kappa "<<kappa<<" xbin "<<xbin<<" fNxbins/2-1 "<<fNxbins/2-1<<endl;
  if(xbin == bin)
    return -1.*kappa;
  else
    return kappa;
}

Double_t AliL3HistogramAdaptive::GetBinCenterY(Int_t ybin)
{
  
  Double_t binwidth = (fYmax - fYmin) / fNybins;
  return fYmin + (ybin-1) * binwidth + 0.5*binwidth;
  
}

void AliL3HistogramAdaptive::Draw(Char_t *option)
{
#ifdef use_root
  if(!fRootHisto)
    CreateRootHisto();
  
  Double_t kappa,psi;
  Int_t content,bin;
  for(Int_t i=0; i<fNxbins; i++)
    {
      kappa = GetBinCenterX(i);
      for(Int_t j=0; j<fNybins; j++)
	{
	  psi = GetBinCenterY(j);
	  bin = GetBin(i,j);
	  content = GetBinContent(bin);
	  fRootHisto->Fill(kappa,psi,content);
	}
    }
  fRootHisto->Draw(option);
  return;
#endif
  cerr<<"AliL3HistogramAdaptive::Draw : You need to compile with ROOT in order to draw histogram"<<endl;
}

void AliL3HistogramAdaptive::Print()
{
  cout<<"Printing content of histogram "<<fName<<endl;
  for(Int_t i=0; i<fNcells; i++)
    {
      if(GetBinContent(i)==0) continue;
      cout<<"Bin "<<i<<": "<<GetBinContent(i)<<endl;
    }

}
