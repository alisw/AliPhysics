//$Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV

#include <string.h>
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
  fKvalue=0;
  fPtstep=0;
}

  
AliL3HistogramAdaptive::AliL3HistogramAdaptive(Char_t *name,Int_t firstpatch,Int_t lastpatch,Double_t minpt,Double_t maxpt,
					       Int_t nybins,Double_t ymin,Double_t ymax)
{
  strcpy(fName,name);
  
  fXmin = -1*AliL3Transform::GetBFact()*AliL3Transform::GetBField()/minpt;
  fXmax = AliL3Transform::GetBFact()*AliL3Transform::GetBField()/minpt;

  fPtstep = 0.1;
  fMinPt = minpt;
  fMaxPt = maxpt;
  fNxbins = InitPtBins(firstpatch,lastpatch,0.1);
  
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
  if(fKvalue)
    delete [] fKvalue;
}

Int_t AliL3HistogramAdaptive::InitPtBins(Int_t firstpatch,Int_t lastpatch,Double_t xyresolution)
{
  
  Double_t pt = fMinPt;
  AliL3Track *track=0;
  Float_t xyz1[3],xyz2[3];
  Double_t angle1,angle2,diff_angle,s_tot,local_pt,pt_nplus1;
  Int_t nrows = AliL3Transform::GetLastRow(lastpatch)-AliL3Transform::GetFirstRow(firstpatch)+1;
  Int_t nptbins=0,index=-1;
  Int_t bounds = (Int_t)rint((fMaxPt-fMinPt)/fPtstep)+1;
  fKvalue = new Double_t[bounds];
  
  while(pt < fMaxPt + fPtstep)
    {
      track = new AliL3Track();
      track->SetPt(pt);
      track->SetPsi(0);
      track->SetCharge(1);
      track->SetFirstPoint(0,0,0);
      track->CalculateHelix();
      
      track->GetCrossingPoint(AliL3Transform::GetFirstRow(firstpatch),xyz1);
      track->GetCrossingPoint(AliL3Transform::GetLastRow(lastpatch),xyz2);
      
      angle1 = atan2((xyz2[1] - track->GetCenterY()),(xyz2[0] - track->GetCenterX()));
      if(angle1 < 0) angle1 += 2.*AliL3Transform::Pi();
      angle2 = atan2((xyz1[1] - track->GetCenterY()),(xyz1[0] - track->GetCenterX()));
      if(angle2 < 0) angle2 += 2.*AliL3Transform::Pi();
      diff_angle = angle1 - angle2;
      diff_angle = fmod(diff_angle,2*AliL3Transform::Pi());
      if((track->GetCharge()*diff_angle) > 0) diff_angle = diff_angle - track->GetCharge()*2.*AliL3Transform::Pi();
      
      s_tot = fabs(diff_angle)*track->GetRadius();
      
      index = (Int_t)rint((pt-fMinPt)/fPtstep);
      if(index < 0 || index > bounds)
	{
	  cerr<<"AliL3HistogramAdaptive::InitPtBins : Index out of range "<<index<<endl;
	  return 0;
	}
      //cout<<"Filling index "<<index<<" pt "<<pt<<" fMinPt "<<fMinPt<<" fPtstep "<<fPtstep<<endl;
      fKvalue[index] = xyresolution*sqrt(720/(nrows+4))/(AliL3Transform::GetBFact()*AliL3Transform::GetBField()*s_tot*s_tot);

      local_pt = pt;
      while(local_pt < pt + fPtstep)
	{
	  if(local_pt > fMaxPt) break;
	  pt_nplus1 = fKvalue[index]*local_pt*local_pt + local_pt;
	  //  cout<<"Pt "<<local_pt<<" binsize "<<(pt_nplus1-local_pt)<<endl;
	  local_pt = pt_nplus1;
	  nptbins++;
	}
      delete track;
      pt += fPtstep;
    }
  cout<<"Number of ptbins "<<nptbins*2<<endl;
  
  return nptbins*2; //Both negative and positive kappa.
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
  if(x < fXmin || x > fXmax)
    return 0;
  
  Double_t pt = fabs(AliL3Transform::GetBFact()*AliL3Transform::GetBField()/x);
  
  if(pt < fMinPt || pt > fMaxPt) 
    return -1;
  Double_t local_pt = fMinPt;
  Double_t pt_nplus1;
  Int_t index,bin=0;
  while(local_pt < fMaxPt)
    {
      index = (Int_t)rint((local_pt-fMinPt)/fPtstep);
      pt_nplus1 = fKvalue[index]*local_pt*local_pt + local_pt;
      //cout<<"Looking for pt "<<local_pt<<" pt_nplus1 "<<pt_nplus1<<" index "<<index<<endl;
      if(pt >= local_pt && pt < pt_nplus1)
	break;
      local_pt = pt_nplus1;
      bin++;
    }
  if(x < 0)
    return bin;
  else
    return fNxbins - bin;
}

Int_t AliL3HistogramAdaptive::FindYbin(Double_t y)
{
  if(y < fYmin || y > fYmax)
    return 0;
  
  return 1 + (Int_t)(fNybins*(y-fYmin)/(fYmax-fYmin));

}

Double_t AliL3HistogramAdaptive::GetBinCenterX(Int_t xbin)
{
  Double_t local_pt = fMinPt;
  Double_t pt_nplus1=0;
  Int_t index,bin=0;
  
  while(local_pt < fMaxPt)
    {
      index = (Int_t)rint((local_pt-fMinPt)/fPtstep);
      pt_nplus1 = fKvalue[index]*local_pt*local_pt + local_pt;
      if(xbin == bin || xbin == fNxbins - bin)
	break;
      local_pt = pt_nplus1;
      bin++;
    }
  
  Double_t binwidth = (pt_nplus1 - local_pt);
  Double_t kappa = AliL3Transform::GetBFact()*AliL3Transform::GetBField()/(local_pt + binwidth*0.5);
  
  cout<<"Returning xbin "<<xbin<<" out of total "<<fNxbins<<endl;
  if(xbin < fNxbins/2)
    return -1.*kappa;
  else
    return kappa;
      
  //Double_t binwidth = (fXmax - fXmin) / fNxbins;
  //return fXmin + (xbin-1) * binwidth + 0.5*binwidth;
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
