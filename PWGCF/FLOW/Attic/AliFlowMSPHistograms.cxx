/************************************************************************* 
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
*                                                                        * 
* Author: The ALICE Off-line Project.                                    * 
* Contributors are mentioned in the code where appropriate.              * 
*                                                                        * 
* Permission to use, copy, modify and distribute this software and its   * 
* documentation strictly for non-commercial purposes is hereby granted   * 
* without fee, provided that the above copyright notice appears in all   * 
* copies and that both the copyright notice and this permission notice   * 
* appear in the supporting documentation. The authors make no claims     * 
* about the suitability of this software for any purpose. It is          * 
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/
#include <iostream>

#include "AliFlowMSPHistograms.h"

#include <TObjString.h>
#include <TList.h>

AliFlowMSPHistograms::AliFlowMSPHistograms() 
   : TNamed(),fVarNames(),fPAvg(),fPProd(),fXName()
{
   // Default constructor for root IO
   //std::cerr << "Default constructor called" << std::endl;
}

AliFlowMSPHistograms::AliFlowMSPHistograms(const AliFlowMSPHistograms &x)
   : TNamed(),fVarNames(),fPAvg(x.fPAvg),fPProd(x.fPProd),fXName(x.fXName)
{
   // Copy constructor, clone all the TProfiles from the other object
   // The copy constructor made a shallow copy, i.e. all entries point to the originial profiles.
   // Therefore we need to Clone each entry.
   for( int i=0; i<fVarNames.GetEntries(); ++i ) {
      TObjString *org=(TObjString *)(fVarNames[i]);
      fVarNames[i]=new TObjString(*org);
   }
   fVarNames.SetOwner(kTRUE);

   for( int i=0; i<fPAvg.GetEntries(); ++i ) {
      TProfile *org=(TProfile *)(fPAvg[i]);
      fPAvg[i]=org->Clone();
   }
   fPAvg.SetOwner(kTRUE);

   for( int i=0; i<fPProd.GetEntries(); ++i ) {
      TProfile *org=(TProfile *)(fPProd[i]);
      fPProd[i]=org->Clone();
   }
   fPProd.SetOwner(kTRUE);
}

AliFlowMSPHistograms::AliFlowMSPHistograms(const int dim, const char *name, const int nBins, const double xmin, const double xmax)
   : TNamed(name,"MSP histograms"),fVarNames(),fPAvg(),fPProd(),fXName()
{
   Init(dim, name, nBins, xmin, xmax);
}

AliFlowMSPHistograms::AliFlowMSPHistograms(const int dim, const char *name, const int nBins, const double *xBins)
   : TNamed(name,"MSP histograms"),fVarNames(),fPAvg(),fPProd(),fXName()
{
   Init(dim, name, nBins, xBins);
}

AliFlowMSPHistograms::~AliFlowMSPHistograms()
{
   // Destructor. All histograms are owned by this object and will be deleted here.
   fPAvg.Clear();
   fPProd.Clear();
}

void AliFlowMSPHistograms::Init(const int dim, const char *name, const int nBins, const double *xBins)
{
   // Initialize all profiles after deleting the existing profiles
   // Variable binning is used for the xaxis
   
   fVarNames.SetOwner(kTRUE);
   fPAvg.SetOwner(kTRUE);
   fPProd.SetOwner(kTRUE);
   fVarNames.Clear();
   fPAvg.Clear();
   fPProd.Clear();

   fXName="";    // Name of x-axis 

   TString profileName;

   for(int i=0; i<dim; ++i) {
      TString vname;
      vname+="<";
      vname+=i;
      vname+=">";
      fVarNames.AddAtAndExpand(new TObjString(vname),i);
   }

   for(int i=0; i<dim; ++i) {
      profileName=name; profileName+=i;
      fPAvg[i] = new TProfile(profileName.Data(), "", nBins, xBins, "s");
      Avg(i)->SetXTitle(name);
      Avg(i)->SetDirectory(0);
      Avg(i)->Sumw2();
   }

   for(int i=0; i<dim; ++i) for(int j=0; j<=i; ++j) {
      profileName=name; profileName+=i; profileName+=j;
      fPProd[Idx(i,j)]=new TProfile(profileName.Data(), "", nBins, xBins, "s");
      Prod(i,j)->SetXTitle(name);
      Prod(i,j)->SetDirectory(0);
      Prod(i,j)->Sumw2();
   }

   UpdateTitles();
}

void AliFlowMSPHistograms::Init(const int dim, const char *name, const int nBins, const double xmin, const double xmax)
{
   // Initialize all profiles after deleting the existing profiles
   // Fixed size binning is used
   
   fVarNames.SetOwner(kTRUE);
   fPAvg.SetOwner(kTRUE);
   fPProd.SetOwner(kTRUE);
   fVarNames.Clear();
   fPAvg.Clear();
   fPProd.Clear();

   fXName="";    // Name of x-axis 

   TString profileName;

   for(int i=0; i<dim; ++i) {
      TString vname;
      vname+="<";
      vname+=i;
      vname+=">";
      fVarNames.AddAtAndExpand(new TObjString(vname),i);
   }

   for(int i=0; i<dim; ++i) {
      profileName=name; profileName+=i;
      fPAvg[i] = new TProfile(profileName.Data(), "", nBins, xmin, xmax, "s");
      Avg(i)->SetXTitle(name);
      Avg(i)->SetDirectory(0);
      Avg(i)->Sumw2();
   }

   for(int i=0; i<dim; ++i) for(int j=0; j<=i; ++j) {
      profileName=name; profileName+=i; profileName+=j;
      fPProd[Idx(i,j)]=new TProfile(profileName.Data(), "", nBins, xmin, xmax, "s");
      Prod(i,j)->SetXTitle(name);
      Prod(i,j)->SetDirectory(0);
      Prod(i,j)->Sumw2();
   }

   UpdateTitles();
}

void  AliFlowMSPHistograms::SetVarName(const char *name, const int i)
{
   if( i<0 || i>NVars() ) return;

   TObjString *old=(TObjString *)(fVarNames.At(i));
   if( old ) {
      delete old;
      fVarNames.RemoveAt(i);
   }
   fVarNames.AddAtAndExpand(new TObjString(name),i);
   UpdateTitles(i);
}

void  AliFlowMSPHistograms::SetVarName(const TString name, const int i)
{
   if( i<0 || i>=NVars() ) return;

   TObjString *old=(TObjString *)fVarNames.At(i);
   if( old ) {
      delete old;
      fVarNames.RemoveAt(i);
   }
   fVarNames.AddAtAndExpand(new TObjString(name),i);
   UpdateTitles(i);
}

const char * AliFlowMSPHistograms::VarName(const int i) const
{
   if( i<0 || i>=NVars() ) return 0;
   return ((TObjString *)fVarNames.At(i))->String().Data();
}

void AliFlowMSPHistograms::Fill(const double x, const double *par, const double *wpar)
{
   // Fill all profiles in the bin corresponding to x with the averages needed for the covariance matrix
   int dim=NVars();
   for(int i=0; i<dim; ++i) {
      Avg(i)->Fill(x,par[i],wpar[i]);
      for(int j=0; j<=i; ++j) {
         Prod(i,j)->Fill(x,par[i]*par[j],wpar[i]*wpar[j]);
      }
   }
}

double AliFlowMSPHistograms::Average(const int bin, const int i)const
{
   if( bin<=0 ) { // Average over all Xbins
      double stats[6];
      Avg(i)->GetStats(stats);
      double w=stats[0];            // sum of weights
      double y=stats[4];            // sum of w*y
      return (w>0?y/w:0);
   }
   return Avg(i)->GetBinContent(bin);
}

double AliFlowMSPHistograms::Covariance(const int bin, const int i, const int j)const
{
   double x,wx,y,wy,xy,wxy;

   if( bin<=0 ) {  // Special case: sum over all bins
      double stats[6];
      Avg(i)->GetStats(stats);
      wx=stats[0];                        // sum of weights from TProfile
      x=(wx>0?stats[4]/wx:0);             // weighted mean from TProfile

      Avg(j)->GetStats(stats);
      wy=stats[0];
      y=(wy>0?stats[4]/wy:0);
      
      Prod(i,j)->GetStats(stats);
      wxy=stats[0];
      xy=(wxy>0?stats[4]/wxy:0);
   }else{
      x=Avg(i)->GetBinContent(bin);       // weighted average per bin from TProfile
      wx=Avg(i)->GetBinEntries(bin);      // sum of weights per bin from TProfile
      y=Avg(j)->GetBinContent(bin);
      wy=Avg(j)->GetBinEntries(bin);
      xy=Prod(i,j)->GetBinContent(bin);
      wxy=Prod(i,j)->GetBinEntries(bin);
   }
   return Covariance(x, y, xy, wx, wy, wxy);
}

Long64_t AliFlowMSPHistograms::Merge(TCollection *c)
{
   // Merge all MSPHistograms in the collection with this
   // The operator+= is used for this 
   // No compatibility checks are done here, only the number of variables has to be compattible
   // The profiles are added, assuming the TProfile::Add function does more checks
   if( !c ) return 0;
   if( c->IsEmpty() ) return 0;

   // Loop over given collection of AliFlowMSPHistograms
   AliFlowMSPHistograms *toMerge=0;
   Long64_t count=0;
   TIter next(c);
   while ( (toMerge=dynamic_cast<AliFlowMSPHistograms *>(next())) ) {
      *this+=*toMerge;
      ++count;
   }
   return count;
}

AliFlowMSPHistograms &AliFlowMSPHistograms::operator+=(const AliFlowMSPHistograms &x)
{
   // Merge this with another AliFlowMSPHistograms object. 
   // The variable names are neither checked nor changed, only the number of variables is checked
   // The profiles are added using the TProfile::Add() function

   int nvars=NVars();
   if( x.NVars() != nvars ) {
      std::cout << "Trying to merge incompattible MSPHistograms. NVars: " << x.NVars() << "!=" << nvars << std::endl;
      return *this;
   }  

   for(int i=0; i<nvars; ++i) {  // Loop over averages histograms
      Avg(i)->Add(x.Avg(i));
      for(int j=0; j<=i; ++j) {  // Loop over products histograms
         Prod(i,j)->Add(x.Prod(i,j));
      }
   }
   return *this;
}

AliFlowMSPHistograms &AliFlowMSPHistograms::operator=(const AliFlowMSPHistograms &x)
{
   // Assignment, clone all the TProfiles from the other object
   // We need to Clone each entry and delete existing ones.
   for( int i=0; i<fVarNames.GetEntries(); ++i ) {
      TObjString *org=(TObjString *)(x.fVarNames[i]);
      delete fVarNames[i];
      fVarNames[i]=new TObjString(*org);
   }
   fVarNames.SetOwner(kTRUE);

   for( int i=0; i<fPAvg.GetEntries(); ++i ) {
      TProfile *org=(TProfile *)(x.fPAvg[i]);
      delete fPAvg[i];
      fPAvg[i]=org->Clone();
   }
   fPAvg.SetOwner(kTRUE);

   for( int i=0; i<fPProd.GetEntries(); ++i ) {
      TProfile *org=(TProfile *)(x.fPProd[i]);
      delete fPProd[i];
      fPProd[i]=org->Clone();
   }
   fPProd.SetOwner(kTRUE);
   return *this;
}

TObject *AliFlowMSPHistograms::Clone(const char *n) const
{
   AliFlowMSPHistograms *x=new AliFlowMSPHistograms(*this);
   if(n)x->SetName(n);
   return x;
}

// Private functions...................................................................................

double AliFlowMSPHistograms::Covariance(double x, double y, double xy, double wx, double wy, double wxy)const
{
   // Calculate the estimated covariance between two variables
   // inputs are the measured averages of x, y and xy and the corresponding sums of weights
   // These are taken from TProfile::GetBinContent(ibin) and TProfile::GetBinEntries(ibin)  
   // See also equations C12-C16 in Ante's thesis

   double enumerator=(xy-x*y)*wxy;
   double denominator= wx*wy-wxy;

   if(TMath::Abs(enumerator)<=1e-50*TMath::Abs(denominator)  ) { // TODO Get smallest Double_t from root ??
      return 0;
   }
   if( TMath::Abs(denominator)<=1e-50*TMath::Abs(enumerator) ) { // protect against infinity 
      return 1e50;
   }
   double result=enumerator/denominator;
   return result;
}

void AliFlowMSPHistograms::UpdateTitles(int k)
{
   // Update the titles of all (k<0) or the profiles for variable k using the variable and parameter names
   // Parameter and x variable names should be set before this function is called

   int dim=NVars();

   for(int i=0; i<dim; ++i) {
      if( i==k || k<0 ) {
         Avg(i)->SetXTitle(XName());
         TString profileTitle="<"; profileTitle+=VarName(i); profileTitle+=">"; 
         Avg(i)->SetYTitle(profileTitle.Data());
         profileTitle+=" vs "; profileTitle+=XName(); 
         Avg(i)->SetTitle(profileTitle.Data());
      }
   }

   for(int i=0; i<dim; ++i) for(int j=0; j<=i; ++j) {
      if( i==k || j==k || k<0 ) {
         Prod(i,j)->SetXTitle(XName());
         TString profileTitle="<"; profileTitle+=VarName(i); profileTitle+="*"; profileTitle+=VarName(j); profileTitle+=">";
         Prod(i,j)->SetYTitle(profileTitle.Data());
         profileTitle+=" vs "; profileTitle+=XName();
         Prod(i,j)->SetTitle(profileTitle.Data());
      }
   }
}

